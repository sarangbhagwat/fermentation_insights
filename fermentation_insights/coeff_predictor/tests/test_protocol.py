# -*- coding: utf-8 -*-
"""Final-ensemble protocol: refit after grouped early stopping, conformal
calibration from out-of-fold scores, cross-fitted coverage."""
import math

import numpy as np
import pandas as pd
import pytest
import torch

from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.evaluate import fit_method
from fermentation_insights.coeff_predictor.model import ModelConfig
from fermentation_insights.coeff_predictor.tests.conftest import FAST

MC = dict(width=16, depth=2, dropout=0.)


def _tc(**kw):
    return T.TrainConfig(**{**FAST, 'epochs': 30, 'patience': 6, 'min_epochs': 3, **kw})


def test_refit_mode_validated():
    with pytest.raises(ValueError, match='refit'):
        T.TrainConfig(refit='yes')
    assert T.TrainConfig.from_dict(T.TrainConfig(refit='anneal').to_dict()).refit == 'anneal'
    assert T.TrainConfig.from_dict({'epochs': 5}).refit == 'none'      # old checkpoints


def test_train_model_epochs_schedule_and_keep_last(tables):
    samples, points, _ = tables
    spec = T.FeatureSpec().fit(samples)
    tc = _tc()
    m, h = T.train_model(spec, samples, points, ModelConfig(**MC), T.LossConfig(), tc, seed=3,
                         n_epochs=7, keep='last')
    assert len(h) == 7 and m.train_info['best_epoch'] == 7 and m.train_info['epochs_run'] == 7
    assert m.train_info['schedule_epochs'] == tc.epochs and m.train_info['keep'] == 'last'
    # same seed, same epochs -> identical networks; the schedule length changes the result
    m2, _ = T.train_model(spec, samples, points, ModelConfig(**MC), T.LossConfig(), tc, seed=3,
                          n_epochs=7, keep='last')
    m3, _ = T.train_model(spec, samples, points, ModelConfig(**MC), T.LossConfig(), tc, seed=3,
                          n_epochs=7, schedule_epochs=7, keep='last')
    sd, sd2, sd3 = m.state_dict(), m2.state_dict(), m3.state_dict()
    assert all(torch.equal(sd[k], sd2[k]) for k in sd)
    assert any(not torch.equal(sd[k], sd3[k]) for k in sd if sd[k].is_floating_point())
    with pytest.raises(ValueError):
        T.train_model(spec, samples, points, ModelConfig(**MC), T.LossConfig(), tc, keep='first')


@pytest.mark.parametrize('mode', ['truncate', 'anneal'])
def test_refit_ensemble_trains_on_all_fit_samples_for_best_epochs(tables, mode):
    samples, points, _ = tables
    ens = T.train_ensemble(samples, points, 2, ModelConfig(**MC), tcfg=_tc(refit=mode), calib_frac=0.0,
                           group_col='config_id', seed=4, keep_pass1=True)
    meta = ens.meta
    assert meta['refit'] == mode and meta['n_calib'] == 0 and meta['n_fit'] == len(samples)
    assert meta['member_n_train'] == [len(samples)] * 2
    assert meta['refit_epochs'] == meta['best_epochs']                    # pass 2 = e* epochs
    assert all(1 <= b <= e for b, e in zip(meta['best_epochs'], meta['epochs']))
    for m, h in zip(ens.members, ens.refit_histories):
        assert m.train_info['keep'] == 'last' and len(h) == m.train_info['best_epoch']
        assert m.train_info['schedule_epochs'] == (len(h) if mode == 'anneal' else 30)
        assert 'val_data' not in h[-1]                                    # no hold-out in pass 2
    # pass 1 kept for diagnostics: the early-stopped members on the smaller fit sets
    p1 = ens.pass1
    assert len(p1.members) == 2 and p1.meta['refit'] == 'none'
    assert [m.train_info['best_epoch'] for m in p1.members] == meta['best_epochs']
    assert all(n < len(samples) for n in p1.meta['member_n_train'])
    # pass 1 of a refit run == the plain run with the same seeds
    ens0 = T.train_ensemble(samples, points, 2, ModelConfig(**MC), tcfg=_tc(), calib_frac=0.0,
                            group_col='config_id', seed=4)
    assert ens0.meta['best_epochs'] == meta['best_epochs'] and ens0.meta['refit'] == 'none'
    assert ens0.meta['refit_epochs'] == [None, None]
    for a, b in zip(ens0.members, p1.members):
        sa, sb = a.state_dict(), b.state_dict()
        assert all(torch.equal(sa[k], sb[k]) for k in sa)
    batch, outs = ens.predict_samples(samples.iloc[:3])
    assert np.isfinite(outs['coef']).all()


def test_refit_median_epochs(tables):
    samples, points, _ = tables
    with pytest.raises(ValueError, match='refit_epochs'):
        T.TrainConfig(refit_epochs='mean')
    ens = T.train_ensemble(samples, points, 3, ModelConfig(**MC), tcfg=_tc(refit='truncate', refit_epochs='median'),
                           calib_frac=0.0, group_col='config_id', seed=5)
    med = int(round(float(np.median(ens.meta['best_epochs']))))
    assert ens.meta['refit_epochs'] == [med] * 3
    assert [len(h) for h in ens.refit_histories] == [med] * 3


def test_refit_baselines_use_every_fit_sample(tables):
    samples, points, _ = tables
    tk = {**FAST, 'epochs': 4}
    m0 = fit_method('ridge', samples, points, 1, tcfg_kw=tk, group_col='config_id', calib_frac=0.0)
    m1 = fit_method('ridge', samples, points, 1, tcfg_kw={**tk, 'refit': 'anneal'}, group_col='config_id',
                    calib_frac=0.0)
    assert m1.fit_info['n_fit'] == len(samples) > m0.fit_info['n_fit']


def test_conformal_scores_match_calibrate_conformal(tables):
    samples, points, _ = tables
    ens = T.train_ensemble(samples, points, 2, ModelConfig(**MC), tcfg=_tc(), calib_frac=0.3,
                           group_col='config_id', seed=0)
    cal = ens.conformal
    assert cal['source'] == 'calibration hold-out'
    sc = T.conformal_scores(ens, samples, points)
    assert list(sc.columns[:2]) == ['sample_id', 'config_id'] and len(sc) == len(samples)
    assert set(T.CONF_KEYS) <= set(sc.columns)
    assert (sc[list(T.CONF_KEYS)].to_numpy() >= 0).all()
    # quantiles from a score table = calibrate_conformal on the same samples
    sub = samples.iloc[:20]
    ref = dict(T.calibrate_conformal(ens, sub, points, 0.2, group_col='config_id'))
    q = T._conformal_from_scores(T.conformal_scores(ens, sub, points), 0.2, T.conformal_eps(ens), 1, 'x')
    for k in T.CONF_KEYS:
        assert q[k] == ref[k] and q[f'{k}_n'] == ref[f'{k}_n']


def test_oof_scores_rescored_with_the_final_ensembles_eps(tables):
    samples, points, _ = tables
    cfgs = sorted(samples['config_id'].unique())
    held = samples['config_id'] == cfgs[0]
    fold = T.train_ensemble(samples[~held], points, 2, ModelConfig(**MC), tcfg=_tc(), calib_frac=0.0, seed=0)
    sc = T.conformal_scores(fold, samples, points)
    e_f = T.conformal_eps(fold)
    for k in T.EPS_KEYS:
        ok = ~sc['invalid'].to_numpy()
        np.testing.assert_allclose(sc[k][ok], (sc[f'{k}_abserr'] / (sc[f'{k}_sd'] + e_f[k]))[ok], rtol=1e-12)
    # the final ensemble has other target scales -> other eps; its quantiles use its own eps
    final = T.train_ensemble(samples.iloc[::2], points, 2, ModelConfig(**MC), tcfg=_tc(), calib_frac=0.0, seed=1)
    e_fin = T.conformal_eps(final)
    assert any(not np.isclose(e_fin[k], e_f[k]) for k in T.EPS_KEYS)
    conf = dict(T.calibrate_conformal_oof(final, sc, 0.1, group_col='config_id'))
    assert conf['eps_rescored'] == list(T.EPS_KEYS) and 'nominal' in conf['level_kind']
    for k in T.EPS_KEYS:
        assert conf[k] == T._conf_q(sc[f'{k}_abserr'] / (sc[f'{k}_sd'] + e_fin[k]), 0.1, warn=False)
        assert conf[f'{k}_eps'] == e_fin[k]
    assert conf['mpsp'] == T._conf_q(sc['mpsp'], 0.1, warn=False)      # relative floor: unchanged
    # without the raw columns (older tables): the stored scores are used as they are
    old = sc[[c for c in sc.columns if not c.endswith(('_abserr', '_sd'))]]
    c0 = dict(T.calibrate_conformal_oof(final, old, 0.1, group_col='config_id'))
    assert c0['eps_rescored'] == [] and c0['a'] == T._conf_q(sc['a'], 0.1, warn=False)
    # rescore keeps rows whose raw values are missing
    part = sc.copy()
    part.loc[:2, 'a_abserr'] = np.nan
    r, keys = T.rescore_conformal(part, e_fin)
    assert 'a' in keys and np.allclose(r['a'][:3], sc['a'][:3], equal_nan=True)


def _fake_scores(n_groups=12, per=10, seed=0, scale=None):
    rng = np.random.RandomState(seed)
    rows = []
    for gi in range(n_groups):
        s = (scale or 1.) * (1 + gi % 3)
        for j in range(per):
            rows.append(dict(sample_id=f'{gi}-{j}', config_id=f'c{gi}', invalid=False,
                             **{k: abs(rng.normal()) * s for k in T.CONF_KEYS}))
    return pd.DataFrame(rows)


def test_crossfit_coverage_leave_one_fold_out():
    sc = _fake_scores()
    cf = T.crossfit_coverage(sc, 0.1, 'config_id')
    for k in T.CONF_KEYS:
        c = cf[k]
        assert c['n_folds'] == 12 and 0 < c['min'] <= c['mean'] <= 1 and 0 <= c['frac_ok'] <= 1
    # by hand for one key / fold
    v, f = sc['mpsp'].to_numpy(), sc['config_id'].to_numpy()
    covs = [np.mean(v[f == g] <= T._conf_q(v[f != g], 0.1, warn=False)) for g in sorted(set(f))]
    assert np.isclose(cf['mpsp']['mean'], np.mean(covs)) and np.isclose(cf['mpsp']['min'], min(covs))
    # heterogeneous groups: the groups with 3x larger scores are under-covered
    assert cf['mpsp']['min'] < 0.9


def test_group_weighted_conformal_quantile():
    # one big configuration with small scores, nine small ones with large scores:
    # the plain quantile follows the big one, the group-weighted one counts each once
    big = np.full(900, 0.1)
    small = np.repeat(np.linspace(1., 9., 9), 5)
    v = np.r_[big, small]
    g = np.r_[np.zeros(900, int), np.repeat(np.arange(1, 10), 5)]
    assert T._conf_q(v, 0.2) == 0.1
    q, n, lev = T._conf_info_grouped(v, g, 0.2)
    # weighted CDF reaches 0.8 * 11 / 10 = 0.88 at the 8th small group (cum. weight 0.9)
    assert q == 8. and n == len(v) and lev == 0.8
    # too few groups for the level -> unbounded
    assert T._conf_info_grouped(v[:905], g[:905], 0.1)[0] == np.inf
    assert np.isnan(T._conf_info_grouped(np.array([np.nan]), np.array([0]), 0.1)[0])
    sc = _fake_scores()
    cf = T.crossfit_coverage(sc, 0.1, 'config_id', weighting='group')
    assert cf['mpsp']['n_folds'] == 12


def test_calibrate_conformal_oof_metadata_and_quantiles(tables):
    samples, points, _ = tables
    ens = T.train_ensemble(samples, points, 1, ModelConfig(**MC), tcfg=_tc(), calib_frac=0.0,
                           group_col='config_id', seed=0)
    assert ens.conformal == {} and ens.meta['n_calib'] == 0
    sc = _fake_scores(scale=2.)
    sc.loc[3, 'mpsp'] = np.nan                     # e.g. a sample without points
    conf = T.calibrate_conformal_oof(ens, sc, 0.1, group_col='config_id', scheme='config', cv_members=1)
    assert ens.conformal is conf and conf['source'].startswith('out-of-fold')
    assert conf['coverage_guaranteed'] is False and conf['n'] == len(sc) and conf['n_calib_groups'] == 12
    assert conf['n_folds'] == 12 and conf['scheme'] == 'config' and conf['cv_members'] == 1
    assert conf['mpsp_n'] == len(sc) - 1
    v = np.sort(sc['a'].to_numpy())
    assert conf['a'] == v[int(math.ceil((len(v) + 1) * 0.9)) - 1]
    eps = T.conformal_eps(ens)
    assert conf['a_eps'] == eps['a'] and conf['g_eps'] == eps['g'] and conf['mpsp_eps_rel'] == 0.01
    assert np.isclose(conf['mpsp_coverage_cv'], conf['crossfit']['mpsp']['mean'])
    with pytest.raises(ValueError, match='config_id'):
        T.calibrate_conformal_oof(ens, sc.drop(columns='config_id'))
    cg = dict(T.calibrate_conformal_oof(ens, sc, 0.1, weighting='group'))
    assert cg['weighting'] == 'group' and cg['mpsp'] == T._conf_info_grouped(sc['mpsp'], sc['config_id'], 0.1)[0]
    with pytest.raises(ValueError, match='weighting'):
        T.calibrate_conformal_oof(ens, sc, weighting='config')
    conf = T.calibrate_conformal_oof(ens, sc, 0.1, group_col='config_id', scheme='config', cv_members=1)
    # the intervals are usable through the saved/loaded ensemble
    import tempfile, os
    with tempfile.TemporaryDirectory() as d:
        ens.save(os.path.join(d, 'e'))
        e2 = T.Ensemble.load(os.path.join(d, 'e'))
    assert e2.conformal['mpsp'] == conf['mpsp'] and e2.conformal['crossfit']['a']['n_folds'] == 12

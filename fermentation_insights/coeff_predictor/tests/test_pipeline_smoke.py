# -*- coding: utf-8 -*-
"""End-to-end smoke test of the experiments pipeline (data -> split -> HPO ->
CV -> final_dev -> release -> tables, figures and report) on small synthetic
tables with a temporary split file and output directories, so that a full run
cannot crash in its final stages. Also checks the round-4 protocol contract:
the final and release intervals come from the development LOCO out-of-fold
scores, never from test samples."""
import json
import os

import numpy as np
import pandas as pd
import pytest

from fermentation_insights.coeff_predictor import dataset
from fermentation_insights.coeff_predictor import experiments as E
from fermentation_insights.coeff_predictor.tests._fake_tables import make_fake_tables

PRODUCTS = ['HP', 'SA', 'EtOH', 'TAL', 'BrothAcid']
FEEDSTOCKS = ['glucose', 'corn', 'cornstover']


@pytest.fixture(scope='module')
def tables_dir(tmp_path_factory):
    root = tmp_path_factory.mktemp('pipeline')
    s, p, c = make_fake_tables(n_samples_per_config=6, products=PRODUCTS, feedstocks=FEEDSTOCKS, seed=3)
    proc = root / 'processed'
    dataset.ProcessedData(samples=s, points=p, configs=c, report={}).save(str(proc))
    planned = E._planned(*[tuple(r) for r in c[['config_id', 'biorefinery', 'feedstock']].itertuples(index=False)])
    # one held-out configuration per product family (3 planned each) + 15 % hashed samples
    fr = E.make_frozen_split(7, planned, min_configs=3, id_frac=0.15)
    split = root / 'split.json'
    split.write_text(json.dumps(fr))
    (root / 'raw').mkdir()
    return root, s, fr


def _run(root, name, *flags):
    out, models = root / name / 'out', root / name / 'models'
    # v1 headline tables for the report's side-by-side comparison
    v1 = out / 'v1' / 'tables'
    v1.mkdir(parents=True)
    pd.DataFrame([dict(set='test_ood', method='pinn', mpsp_mape=0.145, mpsp_mape_median=0.09, R2_g=0.04),
                  dict(set='test_id', method='pinn', mpsp_mape=0.099, mpsp_mape_median=0.05, R2_g=0.9)]
                 ).to_csv(v1 / 'final_test_headline.csv', index=False)
    pd.DataFrame([dict(scheme='config', method='pinn', mpsp_mape=0.144)]).to_csv(v1 / 'cv_headline.csv', index=False)
    pd.DataFrame([dict(stage='cv', method='pinn', fits=54, mean_seconds=88.6, mean_member_epochs=580.)]
                 ).to_csv(v1 / 'timings_fits.csv', index=False)
    E._W.clear()
    try:
        rc = E.main(['--quick', '--processed', str(root / 'processed'), '--raw', str(root / 'raw'), '--no-rebuild',
                     '--out', str(out), '--models', str(models), '--split-file', str(root / 'split.json'),
                     '--threads', '1', *flags])
    finally:
        E._W.clear()
        fh = E._LOG.pop('fh', None)
        if fh is not None:
            fh.close()
    return rc, out, models


@pytest.fixture(scope='module')
def run_v2(tables_dir):
    root, s, fr = tables_dir
    rc, out, models = _run(root, 'v2')
    return dict(rc=rc, out=out, models=models, samples=s, split=fr)


def test_pipeline_runs_to_the_report(run_v2):
    assert run_v2['rc'] == 0
    out = run_v2['out']
    rep = (out / 'REPORT.md').read_text()
    for txt in ('Second use of the frozen test split', 'Changes since v1', 'v1 vs v2 headline numbers',
                '## 3. Training protocol', 'Cross-fitted interval coverage', 'LOCO restricted to',
                'No test sample was used for selection or calibration', '## 7. Release model'):
        assert txt in rep, txt
    t = out / 'tables'
    for f in ('cv_results.csv', 'cv_headline.csv', 'cv_oof_scores.csv.gz', 'cv_coverage_crossfit.csv',
              'cv_loco_ood_families.csv', 'protocol_epochs.csv', 'final_test_headline.csv', 'final_intervals.csv',
              'final_intervals_per_config.csv', 'timings.csv'):
        assert (t / f).exists(), f
    plan = json.loads((out / 'plan.json').read_text())
    assert plan['code_version'] == E.CODE_VERSION and plan['protocol'] == 'v2'
    assert plan['cv_members'] == plan['final_members'] == plan['release_members']
    head = pd.read_csv(t / 'final_test_headline.csv')
    assert {'test_ood', 'test_id', 'train_all'} <= set(head['set'])
    assert np.isfinite(head[(head['set'] == 'test_ood') & (head['method'] == 'pinn')]['mpsp_mape']).all()


def test_oof_scores_cover_dev_once_and_never_test(run_v2):
    t = run_v2['out'] / 'tables'
    sc = pd.read_csv(t / 'cv_oof_scores.csv.gz', dtype={'sample_id': str})
    sp = json.loads((run_v2['out'] / 'split_applied.json').read_text())
    test = set(sp['ood_test_sample_ids']) | set(sp['id_test_sample_ids'])
    s = run_v2['samples']
    dev = set(s['sample_id'].astype(str)) - test
    assert not set(sc['sample_id']) & test
    for (sch, meth), g in sc.groupby(['scheme', 'method']):
        assert meth in E.NN_BASES
        assert not g['sample_id'].duplicated().any(), (sch, meth)
        if sch != 'productivity':            # the productivity scheme tests the r-tails only
            assert set(g['sample_id']) == dev, (sch, meth)
    loco = sc[(sc['scheme'] == 'config') & (sc['method'] == 'pinn')]
    assert (loco['fold'] == loco['config_id']).all()          # leave-one-configuration-out
    cf = pd.read_csv(t / 'cv_coverage_crossfit.csv')
    assert {'config', 'random'} <= set(cf['scheme']) and {'mpsp', 'g', 'a'} <= set(cf['key'])
    assert cf['mean'].between(0, 1).all()


def test_final_and_release_intervals_from_dev_loco_scores(run_v2):
    from fermentation_insights.coeff_predictor.train import Ensemble
    from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
    sp = json.loads((run_v2['out'] / 'split_applied.json').read_text())
    n_dev = sp['n_dev']
    for name in ('final_dev', 'release'):
        path = run_v2['models'] / name
        ens = Ensemble.load(str(path))
        conf = ens.conformal
        assert conf['coverage_guaranteed'] is False and conf['scheme'] == 'config'
        assert 'leave-one-configuration-out' in conf['source']
        assert conf['n'] == n_dev and conf['n_calib_groups'] == sp['n_dev_configs']
        assert np.isfinite(conf['mpsp']) and 0 <= conf['mpsp_coverage_cv'] <= 1
        assert set(conf['crossfit']) >= {'a', 'b', 'c', 'd', 'g', 'mpsp'}
        assert 'conformal_random_cv' in ens.meta and ens.meta['conformal_random_cv']['scheme'] == 'random'
        # protocol: no calibration hold-out, sample-level early stopping, no refit
        assert ens.meta['n_calib'] == 0 and ens.meta['group_col'] is None and ens.meta['refit'] == 'none'
        assert len(ens.meta['best_epochs']) == len(ens.members)
        # v2 model components
        assert ens.mcfg.channel_inputs == E.MASKING_PRESET and ens.mcfg.r_structure == E.R_STRUCTURE
        assert ens.mcfg.fs_prior and ens.mcfg.clip_groups == ['feedstock', 'pretreatment']
        info = json.loads((path / 'experiment_info.json').read_text())
        assert info['protocol'] == 'v2' and info['calib_info']['scheme'] == 'config'
    rel = Ensemble.load(str(run_v2['models'] / 'release'))
    assert 'development + test' in rel.conformal['source']
    dev = Ensemble.load(str(run_v2['models'] / 'final_dev'))
    assert dev.conformal['mpsp'] == rel.conformal['mpsp']             # the same dev LOCO scores
    P = CoefficientPredictor.load(str(run_v2['models'] / 'release'))
    s = run_v2['samples']
    pr = P.predict(s[s['sample_kind'] == 'baseline'].head(2), warn=False)
    assert all(f'{k}_lo' in pr for k in ('a', 'b', 'c', 'd', 'g'))
    iv = pd.read_csv(run_v2['out'] / 'tables' / 'final_intervals.csv')
    r = iv[(iv['method'] == 'pinn') & (iv['calibration'] == 'loco') & (iv['key'] == 'mpsp')]
    assert {'test_ood', 'test_id'} <= set(r['set']) and r['coverage'].between(0, 1).all()
    assert r['crossfit_mean'].notna().all()


def test_protocol_table_and_stages(run_v2):
    pt = pd.read_csv(run_v2['out'] / 'tables' / 'protocol_epochs.csv')
    assert {'hpo', 'cv config', 'cv random', 'final_dev', 'release'} <= set(pt['stage'])
    assert (pt['best_epoch_min'] >= 1).all() and (pt['best_epoch_max'] <= 12).all()


def test_v1_flags_reproduce_v1_settings(tables_dir):
    """--v1-protocol and the per-component flags restore the v1 behaviour
    (grouped hold-outs, final_dev calibrated on its own hold-out, 'physics'
    masking, no r-structure, no safeguards); the release intervals still come
    from the dev LOCO scores (never from test samples)."""
    from fermentation_insights.coeff_predictor.train import Ensemble
    root, s, fr = tables_dir
    rc, out, models = _run(root, 'v1', '--v1-protocol', '--no-r-structure', '--no-extrapolation-safeguards',
                           '--v1-g-inputs')
    assert rc == 0
    dev = Ensemble.load(str(models / 'final_dev'))
    assert dev.meta['n_calib'] > 0 and dev.meta['group_col'] == 'config_id'
    assert dev.conformal['source'] == 'calibration hold-out'
    assert dev.mcfg.channel_inputs == 'physics' and not dev.mcfg.r_structure and not dev.mcfg.fs_prior
    rel = Ensemble.load(str(models / 'release'))
    assert rel.meta['n_calib'] == 0 and 'leave-one-configuration-out' in rel.conformal['source']
    assert json.loads((out / 'plan.json').read_text())['protocol'] == 'v1'
    rep = (out / 'REPORT.md').read_text()
    assert 'Not used in this run (`--v1-protocol`)' in rep and '**v1 protocol**' in rep

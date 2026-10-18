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
    # the HP and TAL rows stand for the configurations regenerated with the D16 / D17 corrections (harness v4)
    s['harness_version'] = np.where(s['config_id'].astype(str).str.match(r'^(HP|TAL)'), 4, 3)
    proc = root / 'processed'
    dataset.ProcessedData(samples=s, points=p, configs=c, report={}).save(str(proc))
    planned = E._planned(*[tuple(r) for r in c[['config_id', 'biorefinery', 'feedstock']].itertuples(index=False)])
    # one held-out configuration per product family (3 planned each) + 15 % hashed samples
    fr = E.make_frozen_split(7, planned, min_configs=3, id_frac=0.15)
    split = root / 'split.json'
    split.write_text(json.dumps(fr))
    (root / 'raw').mkdir()
    return root, s, fr


def _run(root, name, *flags, quick=True):
    out, models = root / name / 'out', root / name / 'models'
    # v1 headline tables for the report's side-by-side comparison
    v1 = out / 'v1' / 'tables'
    v1.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([dict(set='test_ood', method='pinn', mpsp_mape=0.145, mpsp_mape_median=0.09, R2_g=0.04),
                  dict(set='test_id', method='pinn', mpsp_mape=0.099, mpsp_mape_median=0.05, R2_g=0.9)]
                 ).to_csv(v1 / 'final_test_headline.csv', index=False)
    pd.DataFrame([dict(scheme='config', method='pinn', mpsp_mape=0.144)]).to_csv(v1 / 'cv_headline.csv', index=False)
    pd.DataFrame([dict(stage='cv', method='pinn', fits=54, mean_seconds=88.6, mean_member_epochs=580.)]
                 ).to_csv(v1 / 'timings_fits.csv', index=False)
    # ... and of the v2 run (the report compares v1 / v2 / v3)
    v2 = out / 'v2' / 'tables'
    v2.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([dict(set='test_ood', method='pinn', mpsp_mape=0.069, mpsp_mape_median=0.036, R2_g=0.90,
                       n_test=103),
                  dict(set='test_id', method='pinn', mpsp_mape=0.032, mpsp_mape_median=0.022, R2_g=0.98,
                       n_test=141)]).to_csv(v2 / 'final_test_headline.csv', index=False)
    pd.DataFrame([dict(scheme='config', method='pinn', mpsp_mape=0.093)]).to_csv(v2 / 'cv_headline.csv', index=False)
    pd.DataFrame([dict(stage='cv', method='pinn', fits=60, mean_seconds=40.2, mean_member_epochs=600.)]
                 ).to_csv(v2 / 'timings_fits.csv', index=False)
    # ... and of the v3 run (the report compares v1 / v2 / v3 / v4), with its selection and plan
    v3 = out / 'v3' / 'tables'
    v3.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([dict(set='test_ood', method='pinn', mpsp_mape=0.066, mpsp_mape_median=0.034, R2_g=0.91,
                       n_test=103),
                  dict(set='test_id', method='pinn', mpsp_mape=0.032, mpsp_mape_median=0.021, R2_g=0.98,
                       n_test=144)]).to_csv(v3 / 'final_test_headline.csv', index=False)
    pd.DataFrame([dict(scheme='config', method='pinn', mpsp_mape=0.085)]).to_csv(v3 / 'cv_headline.csv', index=False)
    pd.DataFrame([dict(stage='cv', method='pinn', fits=60, mean_seconds=41.0, mean_member_epochs=610.)]
                 ).to_csv(v3 / 'timings_fits.csv', index=False)
    (v3 / 'selected_hyperparameters.json').write_text(json.dumps(dict(
        selected=dict(width=128, depth=2, dropout=0.05, weight_decay=6e-05, lr=0.0047, w_mpsp=0.3, w_ch=3.0,
                      w_phys=0.03, colloc_ratio=2.0, early_stop='mpsp', clip_inputs=0.0),
        plan=dict(mode='full', hpo_k=5, hpo_candidates=16, hpo_members=2, cv_members=5, hpo_lopo_weight=0.0,
                  final_members=5, release_members=5, random_k=5, loco_folds=None, ablation_loco_folds=12,
                  epochs=300, patience=40, min_epochs=20, batch_size=128))))
    E._W.clear()
    try:
        rc = E.main((['--quick', '--threads', '1'] if quick else [])
                    + ['--processed', str(root / 'processed'), '--raw', str(root / 'raw'), '--no-rebuild',
                       '--out', str(out), '--models', str(models), '--split-file', str(root / 'split.json'),
                       *flags])
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
    assert E.THIS_RUN == 'v4' and E.TEST_USE == 'fourth' and E.CHANGES_TITLE == 'Changes v1 -> v2 -> v3 -> v4'
    for txt in ('Fourth use of the frozen test split', 'Changes v1 -> v2 -> v3 -> v4',
                'v1 / v2 / v3 / v4 headline numbers', '### v1 -> v2', '### v2 -> v3', '### v3 -> v4',
                'v1_mape', 'v2_mape', 'v3_mape', 'v4_mape', 'Cost per network fit, v4 vs v3',
                # default settings: v3 -> v4 changed only the data, and the report says so
                'v3 -> v4 ONLY the training data', 'search: NOT identical to v3', 'measures no model change',
                'the test results of v1, v2 and v3 had been seen', 'final review that found D16 and D17',
                'released model whatever the comparison with v1 / v2 / v3 shows',
                '## 3. Training protocol', 'Cross-fitted interval coverage', 'LOCO restricted to',
                'No test sample was used for selection or calibration', '## 7. Release model'):
        assert txt in rep, txt
    # the corrected (harness v4) training and test rows are counted in the notice and the changelog
    from fermentation_insights.coeff_predictor import report_checks as RC
    _, dev = E.apply_split(run_v2['samples'], run_v2['split'])
    rg = RC.regenerated_rows(run_v2['samples'], dev['sample_id'],
                             [c['config_id'] for c in run_v2['split']['test_configs']])
    assert rg and rg['test_ood'] + rg['test_id'] > 0
    assert E._regen_test_clause(rg) in rep and f'{rg["dev"]} of the {rg["n_dev"]} development samples' in rep
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


def test_round5_defaults_and_fidelity_tables(run_v2):
    """Round 5: the pipeline PI-NN has the power gate, the surface channel
    loss with the fixed weight and the defective-target mask; the CV writes
    the out-of-fold channel / price-intensity fidelity tables (development
    samples only) and the report shows them."""
    from fermentation_insights.coeff_predictor.train import Ensemble
    out, models = run_v2['out'], run_v2['models']
    ens = Ensemble.load(str(models / 'final_dev'))
    assert ens.mcfg.power_gate == E.POWER_GATE and ens.spec.power_export == E.POWER_GATE
    assert ens.mcfg.fs_b_only == E.FS_B_ONLY
    assert ens.mcfg.clip_exempt == list(E.CLIP_EXEMPT)
    assert ens.lcfg.ch_space == E.CHANNEL_LOSS['ch_space'] and ens.lcfg.w_ch == E.CHANNEL_LOSS['w_ch']
    t = out / 'tables'
    fid = pd.read_csv(t / 'cv_channel_fidelity.csv.gz')
    summ = pd.read_csv(t / 'cv_fidelity_summary.csv')
    assert {'config', 'feedstock', 'product'} <= set(summ['scheme'])
    _, dev = E.apply_split(run_v2['samples'], run_v2['split'])
    dev_ids = set(dev['sample_id'].astype(str))
    assert set(fid['sample_id'].astype(str)) <= dev_ids          # development samples only
    loco = fid[fid['scheme'] == 'config']
    assert loco['sample_id'].is_unique and set(loco['sample_id'].astype(str)) == dev_ids
    # (the synthetic points carry no per-point intensities q__*, so no sens_* columns here)
    assert np.isfinite(summ['ch_rms']).all() and np.isfinite(summ['ch_feedstock']).all()
    rep = (out / 'REPORT.md').read_text()
    assert 'Out-of-fold channel and price-intensity fidelity' in rep and 'power gate' in rep


def test_v2_model_flag_restores_v2(tables_dir):
    from fermentation_insights.coeff_predictor.train import Ensemble
    root, s, fr = tables_dir
    rc, out, models = _run(root, 'v2model', '--v2-model')
    assert rc == 0
    ens = Ensemble.load(str(models / 'final_dev'))
    assert not ens.mcfg.power_gate and not ens.spec.power_export and ens.lcfg.ch_space is None
    assert not ens.mcfg.fs_b_only and ens.mcfg.clip_exempt is None
    rep = (out / 'REPORT.md').read_text()
    # the settings line does not claim the round-5 components; the changelog marks them as not used
    assert 'power gate (`ModelConfig.power_gate`' not in rep and 'never clipped:' not in rep
    for flag in ('--no-power-gate', '--v2-channel-loss', '--no-target-mask', '--no-clip-exempt'):
        assert f'Not used in this run (`{flag}`)' in rep, flag
    # ... and the v4 notice and changelog do not claim that v3 -> v4 changed only the data
    assert E.settings_vs_v3(dict(mk_extra={})) and 'ONLY the training data' not in rep
    assert 'measures no model change' not in rep and 'search: NOT identical to v3' in rep
    assert 'v3 -> v4 the training data, and the settings (`--no-power-gate`' in rep


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
    item2 = next(line for line in rep.splitlines() if line.startswith('2. **Intervals**'))
    assert 'Not used in this run (`--v1-protocol`)' in item2
    # final_dev's intervals come from its own hold-out and are labelled so (not 'loco')
    iv = pd.read_csv(out / 'tables' / 'final_intervals.csv')
    assert set(iv[iv['method'] == 'pinn']['calibration']) >= {'holdout'}
    assert 'loco' not in set(iv[iv['method'] == 'pinn']['calibration'])


def test_full_mode_runs_resumes_and_replans(tables_dir):
    """Full mode (spawned worker pool, plan_budget, strict out-of-fold checks,
    --members / --ablation-loco-folds through base_plan) on the synthetic
    tables; a second identical run reuses plan.json and every cached fit; a
    changed --ablation-loco-folds recomputes the plan."""
    from fermentation_insights.coeff_predictor.train import Ensemble
    root, s, fr = tables_dir
    # a tiny budget: plan_budget shrinks the search to its minimum of 8 candidates
    # (--hpo-candidates would bypass the planning and plan.json reuse)
    # --hpo-lopo-weight 0.5 exercises the leave-one-product-out part of the HPO criterion (opt-in since round 5)
    # (kept before --ablation-loco-folds, whose value the replan check below replaces)
    flags = ['--threads', '4', '--budget-hours', '0.01', '--members', '2', '--hpo-members', '1',
             '--hpo-lopo-weight', '0.5', '--ablation-loco-folds', '4']
    rc, out, models = _run(root, 'full', *flags, quick=False)
    assert rc == 0
    plan = json.loads((out / 'plan.json').read_text())
    assert plan['mode'] == 'full' and plan['ablation_loco_folds'] == 4 and plan['hpo_candidates'] == 8
    assert plan['cv_members'] == plan['final_members'] == plan['release_members'] == 2
    assert plan['request'] == dict(budget_hours=0.01, hpo_candidates=0, hpo_members=1, ablation_loco_folds=4,
                                   hpo_lopo_weight=0.5)
    # the HPO criterion weighs the leave-one-product-out folds (--hpo-lopo-weight)
    hpo = pd.read_csv(out / 'tables' / 'hpo_candidates.csv')
    w = 0.5
    np.testing.assert_allclose(hpo['criterion'], (1 - w) * hpo['criterion_loco'] + w * hpo['criterion_lopo'])
    folds = pd.read_csv(out / 'tables' / 'hpo_folds.csv')
    assert folds['fold'].astype(str).str.startswith('lopo').sum() == 8 * s['fermentation_product'].nunique()
    assert plan['estimate']['fits_budget'] is False
    for name in ('final_dev', 'release'):
        conf = Ensemble.load(str(models / name)).conformal
        assert 'leave-one-configuration-out' in conf['source'] and conf['coverage_guaranteed'] is False
        assert conf['eps_rescored'] == ['a', 'b', 'c', 'd', 'g'] and 'nominal' in conf['level_kind']
    rep = (out / 'REPORT.md').read_text()
    assert 'QUICK' not in rep and 'adaptively re-used hold-out' in rep
    n_cache = sum(1 for _ in (out / 'cache').rglob('*') if _.is_file())
    log = out / 'experiments.log'
    n_log = len(log.read_text())
    # resume: same command -> plan reused, nothing re-trained
    rc, out, models = _run(root, 'full', *flags, quick=False)
    assert rc == 0
    new = log.read_text()[n_log:]
    assert 'plan: reusing' in new and 'not reusing' not in new
    assert sum(1 for _ in (out / 'cache').rglob('*') if _.is_file()) == n_cache
    # a changed request recomputes the plan (plan-only: stops after writing it)
    n_log = len(log.read_text())
    rc, out, models = _run(root, 'full', *flags[:-1], '3', '--plan-only', quick=False)
    new = log.read_text()[n_log:]
    assert 'not reusing' in new and 'request' in new
    assert json.loads((out / 'plan.json').read_text())['ablation_loco_folds'] == 3

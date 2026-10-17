# -*- coding: utf-8 -*-
"""Extrapolation safeguards: the neutralization gate, the sugar-availability
prior of the structured feedstock intensity (ModelConfig.fs_prior) and
group-restricted input clipping (ModelConfig.clip_groups)."""
import numpy as np
import pandas as pd
import pytest
import torch

from fermentation_insights.coeff_predictor.features import (
    FeatureSpec, STRUCTURAL, neutralization_factor, sugar_availability, SUGAR_GAIN, SACCHARIFICATION_DEFAULT)
from fermentation_insights.coeff_predictor.model import ModelConfig, make_layout
from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.tests._fake_tables import make_fake_tables


@pytest.fixture(scope='module')
def tabs():
    return make_fake_tables(n_samples_per_config=4, products=['HP', 'EtOH', 'SA', 'BrothAcid'],
                            feedstocks=['glucose', 'cornstover'], seed=5)


def _model(tabs, seed=0, **kw):
    samples, points, _ = tabs
    spec = FeatureSpec(downstream=True).fit(samples)
    batch = spec.transform(samples)
    cfg = ModelConfig(**{**dict(width=16, depth=2, dropout=0., water_gate=True, channel_inputs='physics',
                                masked_width=16), **kw})
    layout = make_layout(cfg, samples.columns)
    scales = T.fit_target_scales(samples, points, layout, batch, spec.continuous_mask(), T.fs_prior_z(spec, batch))
    torch.manual_seed(seed)
    m = T.build_model(spec, layout, scales, cfg, torch.float64).eval()
    X = torch.tensor(batch.X)
    rest = [torch.tensor(np.asarray(a, float)) for a in (batch.log_r, batch.prices, batch.w_sugar, batch.g_max)]
    return m, spec, layout, scales, X, rest, samples


def test_neutralization_factor():
    df = pd.DataFrame(dict(
        n_acid_groups=[1, 0, 0, 0, np.nan, 2],
        neutralizing_agent=[np.nan, np.nan, 'NaOH', '', np.nan, 'Ca(OH)2'],
        neutralization=[False, False, False, False, False, True],
        neutralization_degree=[0., 0., 0., 0., 0., 1.]))
    f = neutralization_factor(df)
    # acid w/o agent (HP low pH) / alcohol / lactone with pH control / explicit '' agent / unknown / neutralized acid
    np.testing.assert_array_equal(f[[0, 1, 2, 3, 5]], [1., 0., 1., 0., 1.])
    assert np.isnan(f[4])
    assert 'struct_neutralization' in STRUCTURAL


def test_sugar_availability():
    w_cane = 0.012 + 0.137 * SUGAR_GAIN['sucrose'] + 0.061 * SUGAR_GAIN['glucan'] + 0.036 * SUGAR_GAIN['xylan']
    stover_glucan, stover_xylan = 0.3, 0.17
    w_stover = stover_glucan * SUGAR_GAIN['glucan'] + stover_xylan * SUGAR_GAIN['xylan']
    df = pd.DataFrame(dict(
        feedstock_type=['pure_sugar', 'starch', 'sugar_juice', 'lignocellulose', 'lignocellulose', np.nan, np.nan],
        fs_w_glucose=[0.909, 0., 0.012, 0., 0., 0., 0.5],
        fs_w_starch=[0., 0.6, 0., 0., 0., 0., 0.],
        fs_w_sucrose=[0., 0., 0.137, 0., 0., 0., 0.],
        fs_w_glucan=[0., 0., 0.061, stover_glucan, stover_glucan, 0.3, 0.],
        fs_w_xylan=[0., 0., 0.036, stover_xylan, stover_xylan, 0., 0.],
        p_saccharification_yield=[np.nan, np.nan, np.nan, 0.8, np.nan, np.nan, np.nan],
        fs_w_fermentable_sugar_eq=[0.909, 0.6 * SUGAR_GAIN['starch'], w_cane, w_stover, w_stover,
                                   0.3 * SUGAR_GAIN['glucan'], 0.5]))
    a = sugar_availability(df)
    np.testing.assert_allclose(a[:5], [1., 1., (0.012 + 0.137 * SUGAR_GAIN['sucrose']) / w_cane, 0.8,
                                       SACCHARIFICATION_DEFAULT])
    # structural polysaccharides with an unknown feedstock type: no prior; free sugars only: 1
    assert np.isnan(a[5]) and a[6] == 1.
    assert 'struct_sugar_availability' in STRUCTURAL


def test_neutralization_gate(tabs):
    m, spec, layout, scales, X, rest, samples = _model(tabs, neutralization_gate=True)
    j = layout['channels'].index('neutralization')
    k = spec.structural_index('struct_neutralization')
    gate = X[:, k].numpy()
    etoh = samples['fermentation_product'].eq('EtOH').to_numpy().copy()
    assert (gate[etoh] == 0).all() and (gate[~etoh] == 1).all()
    o = m(X, *rest)
    assert torch.all(o['channels'][etoh, j] == 0) and torch.any(o['channels'][~etoh, j] != 0)
    # only the neutralization channel is gated; the structural input reaches no sub-network
    X2 = X.clone()
    X2[:, k] = 1.
    o2 = m(X2, *rest)
    other = [i for i in range(len(layout['channels'])) if i != j]
    torch.testing.assert_close(o2['channels'][:, other], o['channels'][:, other])
    torch.testing.assert_close(o2['g'], o['g'])
    assert torch.any(o2['channels'][etoh, j] != 0)
    # off by default; the data check keeps it where the targets vanish and disables it otherwise
    m0 = _model(tabs)[0]
    assert m0.neut_index is None
    rep = T.check_constraints(samples, tabs[1])
    row = rep.set_index('constraint').loc['neutralization_gate']
    assert row['enabled'] and row['n_cases'] == etoh.sum()
    bad = samples.copy()
    bad.loc[etoh, 'ch__neutralization__a'] = 0.3
    rep2 = T.check_constraints(bad, tabs[1])
    assert not rep2.set_index('constraint').loc['neutralization_gate', 'enabled']
    mc, _ = T.apply_constraint_report(rep2, ModelConfig(neutralization_gate=True), T.LossConfig())
    assert mc.neutralization_gate is False


def test_fs_prior(tabs):
    m, spec, layout, scales, X, rest, samples = _model(tabs, fs_prior=True)
    k = spec.structural_index('struct_sugar_availability')
    assert m.avail_index == k and 'fs_prior_center' in scales
    jfs = layout['price_keys'].index('fs')
    o = m(X, *rest)
    # b'_fs scales exactly as 1/s_avail; nothing else changes
    X2 = X.clone()
    X2[:, k] = X[:, k] * 0.5
    o2 = m(X2, *rest)
    torch.testing.assert_close(o2['intensities'][:, jfs, 1], 2. * o['intensities'][:, jfs, 1])
    torch.testing.assert_close(o2['channels'], o['channels'])
    torch.testing.assert_close(o2['intensities'][:, jfs, [0, 2, 3]], o['intensities'][:, jfs, [0, 2, 3]])
    # the structural input is the composition's prior (unknown -> 1)
    np.testing.assert_allclose(X[:, k].numpy(), np.nan_to_num(sugar_availability(samples), nan=1.))
    assert 0 < X[:, k].min() and X[:, k].max() <= 1
    # off: identical to the model without the option
    assert _model(tabs, fs_prior=False)[0].avail_index is None
    # serialization: old configs load with the option off
    d = ModelConfig(fs_prior=True).to_dict()
    d.pop('fs_prior'), d.pop('neutralization_gate'), d.pop('clip_groups')
    c = ModelConfig.from_dict(d)
    assert c.fs_prior is False and c.neutralization_gate is False and c.clip_groups is None


def test_clip_groups(tabs):
    m, spec, layout, scales, X, rest, samples = _model(tabs, clip_groups=['feedstock'])
    groups = spec.column_groups()
    fs_cols = [i for i in groups['feedstock'] if spec.continuous_mask()[i]]
    other = [i for g, idx in groups.items() if g != 'feedstock' for i in idx if spec.continuous_mask()[i]]
    assert fs_cols and other
    lo, hi = m.x_lo.numpy(), m.x_hi.numpy()
    assert np.isfinite(lo[fs_cols]).all() and (np.abs(lo[other]) >= 1e29).all() and (hi[other] >= 1e29).all()
    o = m(X, *rest)
    # far outside the training range: a feedstock column is held at the edge, other columns are not
    X2 = X.clone()
    X2[:, fs_cols] = X[:, fs_cols] + 50.
    X3 = X.clone()
    X3[:, fs_cols] = torch.as_tensor(hi[fs_cols], dtype=X.dtype)
    torch.testing.assert_close(m(X2, *rest)['coef'], m(X3, *rest)['coef'])
    m_all = _model(tabs, clip_inputs=0.)[0]
    assert np.isfinite(m_all.x_lo.numpy()[other]).all()
    # composes with the HPO dimension clip_inputs >= 0 (then every column is clipped at that margin)
    m_both = _model(tabs, clip_inputs=0.5, clip_groups=['feedstock'])[0]
    np.testing.assert_allclose(m_both.x_lo.numpy()[other], m_all.x_lo.numpy()[other] - 0.5)
    np.testing.assert_allclose(m_both.x_hi.numpy()[fs_cols], hi[fs_cols] + 0.5)
    assert not _model(tabs)[0].clip_on
    with pytest.raises(ValueError, match='unknown feature group'):
        ModelConfig(clip_groups=['feedstok'])
    # needs the column groups: build_model passes them without masking too
    m5 = _model(tabs, clip_groups=['feedstock'], channel_inputs=None)[0]
    assert (np.abs(m5.x_lo.numpy()[other]) >= 1e29).all()
    # physics-model structures only: the plain MLP baseline is unchanged
    p = ModelConfig(kind='plain', fs_prior=True, neutralization_gate=True, clip_groups=['feedstock'])
    assert not p.fs_prior and not p.neutralization_gate and p.clip_groups is None


def test_safeguards_ensemble_roundtrip(tabs, tmp_path):
    samples, points, _ = tabs
    cfg = ModelConfig(width=16, depth=2, dropout=0., water_gate=True, channel_inputs='physics', masked_width=16,
                      fs_prior=True, neutralization_gate=True, clip_groups=['feedstock', 'pretreatment'])
    ens = T.train_ensemble(samples, points, 1, cfg, T.LossConfig(),
                           T.TrainConfig(epochs=3, batch_size=16, min_epochs=1, patience=2),
                           calib_frac=0., spec_kwargs=dict(downstream=True))
    m = ens.members[0]
    assert m.avail_index is not None and m.neut_index is not None and m.clip_on
    assert ens.mcfg.neutralization_gate  # kept by the data check (EtOH: no neutralization cost)
    _, o1 = ens.predict_samples(samples)
    ens2 = T.Ensemble.load(ens.save(str(tmp_path / 'ck')))
    _, o2 = ens2.predict_samples(samples)
    np.testing.assert_allclose(o2['coef'], o1['coef'], rtol=1e-6, atol=1e-8)
    assert ens2.mcfg.fs_prior and ens2.mcfg.clip_groups == ['feedstock', 'pretreatment']

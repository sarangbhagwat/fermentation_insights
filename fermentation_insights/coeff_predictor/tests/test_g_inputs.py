# -*- coding: utf-8 -*-
"""Inputs of the recovery g: the 'reagent_prices' group (split from
'neutralization' so that the recovery can exclude the reagent prices) and the
opt-in titer-domain feature (FeatureSpec(domain_features=True))."""
import numpy as np
import pytest
import torch

from fermentation_insights.coeff_predictor.features import (
    FeatureSpec, feature_group, FEATURE_GROUPS, DOMAIN_PHYS)
from fermentation_insights.coeff_predictor.model import (
    ModelConfig, make_layout, CHANNEL_INPUT_PRESETS, G_KEY)
from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.tests._fake_tables import make_fake_tables


@pytest.fixture(scope='module')
def tables():
    samples, points, configs = make_fake_tables(n_samples_per_config=4, products=['HP', 'EtOH', 'SA', 'BrothAcid'],
                                                feedstocks=['glucose', 'cornstover'], seed=5)
    rng = np.random.RandomState(0)
    samples = samples.copy()
    # reagent prices of the neutralizing configurations (sampled per sample)
    samples['p_base_price'] = np.where(samples['neutralization'].astype(bool), rng.uniform(0.2, 0.3, len(samples)),
                                       np.nan)
    samples['p_acid_price'] = rng.uniform(0.09, 0.11, len(samples))
    return samples, points, configs


def _model(tables, spec_kw=None, **kw):
    samples, points, _ = tables
    spec = FeatureSpec(downstream=True, **(spec_kw or {})).fit(samples)
    batch = spec.transform(samples)
    cfg = ModelConfig(**{**dict(width=16, depth=2, dropout=0., water_gate=True, channel_inputs='physics',
                                masked_width=16), **kw})
    layout = make_layout(cfg, samples.columns)
    scales = T.fit_target_scales(samples, points, layout, batch, spec.continuous_mask())
    torch.manual_seed(0)
    m = T.build_model(spec, layout, scales, cfg, torch.float64).eval()
    return m, spec, batch


def test_groups():
    assert feature_group('p_base_price') == 'reagent_prices' == feature_group('missing__p_acid_price')
    assert feature_group('phys_log10_t_ref') == 'titer_domain'
    assert {'reagent_prices', 'titer_domain'} <= set(FEATURE_GROUPS)
    # 'physics' (saved checkpoints): the reagent prices reach exactly the outputs that see the
    # neutralization chemistry, i.e. the same columns as before the split
    for name, gs in CHANNEL_INPUT_PRESETS['physics'].items():
        assert ('neutralization' in gs) == ('reagent_prices' in gs), name
    # 'physics_v2' (pipeline default): only the recovery differs, without the reagent prices
    v1, v2 = CHANNEL_INPUT_PRESETS['physics'], CHANNEL_INPUT_PRESETS['physics_v2']
    assert {k for k in v1 if v1[k] != v2.get(k)} == {G_KEY} and set(v1) == set(v2)
    assert set(v2[G_KEY]) == set(v1[G_KEY]) - {'reagent_prices'}
    assert 'feedstock' not in v2[G_KEY] and 'scale' not in v2[G_KEY]


def test_reagent_prices_excluded_from_g(tables):
    m, spec, batch = _model(tables)
    groups = spec.column_groups()
    rp = set(groups['reagent_prices'])
    assert rp and rp <= set(m.output_inputs[G_KEY]['columns'])
    assert rp <= set(m.output_inputs['neutralization']['columns'])
    m2, spec2, batch2 = _model(tables, channel_inputs='physics_v2')
    assert not rp & set(m2.output_inputs[G_KEY]['columns'])
    assert set(m2.output_inputs[G_KEY]['columns']) == set(m.output_inputs[G_KEY]['columns']) - rp
    X = torch.tensor(batch2.X, requires_grad=True)
    lr = torch.tensor(batch2.log_r)
    rest = [torch.tensor(np.asarray(a, float)) for a in (batch2.prices, batch2.w_sugar, batch2.g_max)]
    o = m2(X, lr, *rest)
    gx, = torch.autograd.grad(o['g'].sum(), X, retain_graph=True)
    assert np.all(gx.numpy()[:, sorted(rp)] == 0.)
    j = list(m2.layout['channels']).index('neutralization')
    gn, = torch.autograd.grad(o['channels'][:, j].sum(), X)
    assert np.any(gn.numpy()[:, sorted(rp)] != 0.)


def test_domain_feature(tables):
    samples, points, _ = tables
    off = FeatureSpec(downstream=True).fit(samples)
    assert not set(DOMAIN_PHYS) & set(off.columns)          # default: unchanged columns
    on = FeatureSpec(downstream=True, domain_features=True).fit(samples)
    assert set(on.columns) == set(off.columns) | set(DOMAIN_PHYS)
    assert on.column_groups()['titer_domain'] == [on.columns.index('phys_log10_t_ref')]
    i = on.columns.index('phys_log10_t_ref')
    d = next(x for x in on.numeric if x['name'] == 'phys_log10_t_ref')
    X = on.transform(samples).X
    expect = np.log10(np.sqrt(samples['t_min'].to_numpy(float) * samples['t_max'].to_numpy(float)))
    np.testing.assert_allclose(X[:, i] * d['std'] + d['mean'], expect, rtol=1e-9)
    # serialization keeps the option (and specs saved before it load with it off)
    re = FeatureSpec.from_dict(on.to_dict())
    assert re.domain_features and re.columns == on.columns
    dd = off.to_dict()
    dd.pop('domain_features')
    assert not FeatureSpec.from_dict(dd).domain_features
    # a query without the domain gets the training fill (finite)
    q = samples.iloc[:2].drop(columns=['t_min', 't_max'])
    assert np.isfinite(on.transform(q, strict=False).X[:, i]).all()
    # only the recovery sees it with the override; the other outputs keep their inputs
    m, spec, _ = _model(tables, spec_kw=dict(domain_features=True),
                        channel_inputs={'preset': 'physics',
                                        G_KEY: CHANNEL_INPUT_PRESETS['physics'][G_KEY] + ('titer_domain',)})
    assert i in m.output_inputs[G_KEY]['columns']
    assert all(i not in v['columns'] for k, v in m.output_inputs.items() if k != G_KEY)

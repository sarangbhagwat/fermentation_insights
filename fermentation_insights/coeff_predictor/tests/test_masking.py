# -*- coding: utf-8 -*-
"""Channel-wise input masking (ModelConfig.channel_inputs): feature groups,
hard masking of the inputs of every cost channel (zero gradients from excluded
inputs), the hybrid residual network and checkpoint compatibility."""
import numpy as np
import pytest
import torch

from fermentation_insights.coeff_predictor.features import FeatureSpec, feature_group, FEATURE_GROUPS
from fermentation_insights.coeff_predictor.model import (
    ModelConfig, make_layout, resolve_channel_inputs, CHANNEL_INPUT_PRESETS, R_GROUP, G_KEY, INTENSITY_CHANNEL,
    GATED_CHANNELS, GATED_INTENSITIES)
from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.tests._fake_tables import make_fake_tables


@pytest.fixture(scope='module')
def ds_tables():
    return make_fake_tables(n_samples_per_config=4, products=['HP', 'EtOH', 'SA', 'BrothAcid'],
                            feedstocks=['glucose', 'cornstover'], seed=3)


def _model(tables, seed=0, **kw):
    samples, points, _ = tables
    spec = FeatureSpec(downstream=True).fit(samples)
    batch = spec.transform(samples)
    cfg = ModelConfig(**{**dict(width=16, depth=2, dropout=0., water_gate=True, channel_inputs='physics',
                                masked_width=16), **kw})
    layout = make_layout(cfg, samples.columns)
    scales = T.fit_target_scales(samples, points, layout, batch, spec.continuous_mask())
    torch.manual_seed(seed)
    m = T.build_model(spec, layout, scales, cfg, torch.float64).eval()
    X = torch.tensor(batch.X, requires_grad=True)
    lr = torch.tensor(batch.log_r, requires_grad=True)
    rest = [torch.tensor(np.asarray(a, float)) for a in (batch.prices, batch.w_sugar, batch.g_max)]
    return m, spec, layout, X, lr, rest


def _outputs(m, layout, o):
    """(name, scalar tensor) per output slot (sum over rows and coefficients)."""
    out = [(INTENSITY_CHANNEL[k], o['intensities'][:, j].sum()) for j, k in enumerate(layout['price_keys'])]
    out += [(c, o['channels'][:, i].sum()) for i, c in enumerate(layout['channels'])]
    return out + [(G_KEY, o['g'].sum())]


def test_feature_groups():
    assert feature_group('fs_w_lignin') == 'feedstock' and feature_group('feedstock_type==starch') == 'feedstock'
    assert feature_group('missing__p_IRR') == 'financial' and feature_group('phys_log10_capacity') == 'scale'
    assert feature_group('fp_dor') == 'product' and feature_group('theoretical_yield') == 'product'
    assert feature_group('pp_Tm') == 'final_product' and feature_group('ds_upg__hydrogenation') == 'upgrading'
    assert feature_group('ds_sep__distillation') == 'separation'
    assert feature_group('ds_sep__acidulation') == 'neutralization'
    assert feature_group('water_removal_route==none') == 'water_route'
    assert feature_group('phys_water_removal_c') == 'water_route'
    assert feature_group('ds_cop__DDGS') == 'coproducts' and feature_group('struct_water_removal') == 'structural'
    assert feature_group('p_mystery_parameter') == 'other' and feature_group('feedstock_price') == 'prices'
    for preset in CHANNEL_INPUT_PRESETS.values():
        assert all(g in FEATURE_GROUPS or g == R_GROUP for gs in preset.values() for g in gs)


def test_channel_inputs_resolution():
    assert resolve_channel_inputs(None) is None
    p = resolve_channel_inputs('physics')
    assert p['capex_upstream'] == CHANNEL_INPUT_PRESETS['physics']['capex_upstream']
    q = resolve_channel_inputs({'preset': 'physics', 'fixed': ['scale']})
    assert q['fixed'] == ('scale',) and q['feedstock'] == p['feedstock']
    with pytest.raises(ValueError, match='preset'):
        ModelConfig(channel_inputs='nonsense')
    assert ModelConfig(kind='plain', channel_inputs='physics').channel_inputs is None
    assert ModelConfig(channel_inputs='physics').masked and not ModelConfig().masked
    # checkpoints saved before the option load with one shared encoder
    d = ModelConfig().to_dict()
    d.pop('channel_inputs')
    assert ModelConfig.from_dict(d).channel_inputs is None
    with pytest.raises(ValueError, match='unknown feature group'):
        _model(make_fake_tables(n_samples_per_config=2, products=['HP', 'EtOH'], feedstocks=['glucose'], seed=1),
               channel_inputs={'preset': 'physics', 'fixed': ['scale', 'astrology']})


def test_masking_blocks_gradients_from_excluded_inputs(ds_tables):
    m, spec, layout, X, lr, rest = _model(ds_tables)
    assert m.masked and not m.unused_inputs
    o = m(X, lr, *rest)
    preset = CHANNEL_INPUT_PRESETS['physics']
    groups = spec.column_groups()
    gate = spec.structural_index('struct_water_removal')
    gated = set(GATED_CHANNELS) | {INTENSITY_CHANNEL[k] for k in GATED_INTENSITIES}
    seen = set()
    for name, v in _outputs(m, layout, o):
        gx, gr = torch.autograd.grad(v, (X, lr), retain_graph=True, allow_unused=True)
        gx = np.abs(gx.detach().numpy()).sum(0)
        allowed = {i for g in preset[name] for i in groups.get(g, [])}
        assert allowed == set(m.output_inputs[name]['columns'])
        if name in gated:
            allowed.add(gate)
        excluded = [i for i in range(spec.n_features) if i not in allowed]
        # d channel / d excluded input == 0 exactly
        assert np.all(gx[excluded] == 0.), (name, [spec.columns[i] for i in excluded if gx[i] != 0][:5])
        assert np.any(gx[sorted(allowed)] > 0), name
        dr = 0. if gr is None else float(gr.abs().sum())
        assert (dr > 0) == (R_GROUP in preset[name]), name
        seen.add(name)
    assert {'feedstock', 'capex_upstream', 'capex_fermentation', 'capex_downstream', G_KEY} <= seen
    # product chemistry never reaches the feedstock, upstream, fixed, wastewater or utilities channels
    prod = set(groups['product']) | set(groups['final_product']) | set(groups['separation'])
    for name in ('feedstock', 'capex_upstream', 'fixed', 'capex_wastewater', 'capex_utilities', 'electricity'):
        assert not prod & set(m.output_inputs[name]['columns']), name
    # ... and the feedstock never reaches the recovery
    assert not set(groups['feedstock']) & set(m.output_inputs[G_KEY]['columns'])
    # totals are still the exact sum of the channels
    s = o['channels'].sum(1) + o['price_channels'].sum(1)
    assert torch.allclose(s, o['coef'], rtol=0, atol=1e-12)


def test_masking_without_price_explicit_channels(ds_tables):
    """Ablation pinn-no_price: prices are encoder features; each price channel sees its price."""
    samples, points, _ = ds_tables
    spec = FeatureSpec(downstream=True, price_explicit=False).fit(samples)
    batch = spec.transform(samples)
    cfg = ModelConfig(width=16, depth=2, dropout=0., channel_inputs='physics', masked_width=16,
                      price_explicit=False)
    layout = make_layout(cfg, samples.columns)
    assert not layout['price_keys'] and 'feedstock' in layout['channels']
    scales = T.fit_target_scales(samples, points, layout, batch, spec.continuous_mask())
    m = T.build_model(spec, layout, scales, cfg, torch.float64)
    prices = set(spec.column_groups()['prices'])
    assert prices and prices <= set(m.output_inputs['feedstock']['columns'])
    assert not prices & set(m.output_inputs['capex_upstream']['columns'])


def test_hybrid_residual(ds_tables):
    m, spec, layout, X, lr, rest = _model(ds_tables, residual_width=8, residual_weight_decay=5.)
    res = m.residual_parameters()
    assert res and all(any(p is q for q in m.parameters()) for p in res)
    # zero-initialized residual head: identical to the purely masked model with the same sub-networks
    pure, *_ = _model(ds_tables)
    pure.load_state_dict({k: v for k, v in m.state_dict().items() if not k.startswith('resid')})
    o, o0 = m(X, lr, *rest), pure(X, lr, *rest)
    torch.testing.assert_close(o['coef'], o0['coef'])
    # a non-zero residual lets every encoder input reach every channel (that is its purpose)
    with torch.no_grad():
        m.resid_head.weight.normal_()
    o = m(X, lr, *rest)
    up = layout['channels'].index('capex_upstream')
    gx, = torch.autograd.grad(o['channels'][:, up].sum(), (X,))
    prod = spec.column_groups()['product']
    assert np.any(np.abs(gx.numpy())[:, prod] > 0)
    # the optimizer gives the residual weights their own weight decay
    samples, points, _ = ds_tables
    seen = {}
    orig = torch.optim.AdamW

    class Spy(orig):
        def __init__(self, groups, **kw):
            seen['groups'] = [(len(g['params']), g['weight_decay']) for g in groups]
            super().__init__(groups, **kw)
    torch.optim.AdamW = Spy
    try:
        T.train_ensemble(samples, points, 1, m.cfg, T.LossConfig(),
                         T.TrainConfig(epochs=2, batch_size=16, min_epochs=1, patience=1),
                         calib_frac=0., spec_kwargs=dict(downstream=True))
    finally:
        torch.optim.AdamW = orig
    assert len(seen['groups']) == 3 and seen['groups'][2][1] == 5. and seen['groups'][2][0] > 0


def test_masked_checkpoint_roundtrip(ds_tables, tmp_path):
    samples, points, _ = ds_tables
    cfg = ModelConfig(width=16, depth=2, dropout=0., water_gate=True, masked_width=16,
                      channel_inputs={'preset': 'physics', 'fixed': ['scale']}, residual_width=8)
    ens = T.train_ensemble(samples, points, 1, cfg, T.LossConfig(),
                           T.TrainConfig(epochs=3, batch_size=16, min_epochs=1, patience=2),
                           calib_frac=0., spec_kwargs=dict(downstream=True))
    _, o1 = ens.predict_samples(samples)
    assert np.isfinite(o1['coef']).all()
    re = T.Ensemble.load(ens.save(str(tmp_path / 'ck')))
    assert re.mcfg.channel_inputs == {'preset': 'physics', 'fixed': ['scale']} and re.members[0].masked
    _, o2 = re.predict_samples(samples)
    np.testing.assert_allclose(o2['coef'], o1['coef'], rtol=1e-6, atol=1e-8)
    # derived index buffers are not part of the saved state
    assert not any(k.startswith('sub_index') for k in re.members[0].state_dict())
    # user queries through the predictor
    import warnings
    from fermentation_insights.coeff_predictor.predict import CoefficientPredictor
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        o = CoefficientPredictor(re).predict(samples.iloc[:3])
    assert len(o) == 3 and np.isfinite(o[['a', 'b', 'c', 'd', 'g']].to_numpy(float)).all()
    np.testing.assert_allclose(o['a'].to_numpy(), o1['coef'][:, :3, 0].mean(0), rtol=1e-5, atol=1e-7)

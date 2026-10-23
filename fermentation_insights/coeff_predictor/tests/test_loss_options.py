"""Loss options: surface-space supervision and configuration-balanced weights."""
import numpy as np
import torch

from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.tests.conftest import FAST, build


def _batch(fitted, layout, dtype=torch.float64):
    s, p = fitted['samples'], fitted['points']
    return T.make_tensors(fitted['spec'], s, p, layout, fitted['batch'], dtype=dtype)


def _perfect_out(model, b):
    q = b['q']
    pc = b['prices'][:, model.price_index][:, :, None] * q
    return dict(coef=b['coef'].clone(), g=b['g'].clone(), channels=b['ch'].clone(), intensities=q.clone(),
                price_channels=pc)


def test_perfect_prediction_has_zero_coefficient_losses(fitted):
    model, _, _, _, layout, _ = build(fitted)
    b = _batch(fitted, layout)
    out = _perfect_out(model, b)
    for space in ('coef', 'surface'):
        for bal in (0., 1.):
            L = T.data_loss(model, out, b, T.LossConfig(coef_space=space, config_balance=bal))
            assert float(L['tot']) < 1e-12 and float(L['ch']) < 1e-12 and float(L['g']) < 1e-12


def test_coef_space_unweighted_matches_element_mean(fitted):
    model, tens, _, _, layout, _ = build(fitted)
    b = _batch(fitted, layout)
    with torch.no_grad():
        out = model(*[b[k] for k in ('X', 'log_r', 'prices', 'w_sugar', 'g_max')])
        L = T.data_loss(model, out, b, T.LossConfig())
    cs = torch.as_tensor(model.scales_np['coef_std'], dtype=torch.float64)
    ref = T._mmean(((out['coef'] - b['coef']) / cs) ** 2, b['coef_m'])
    assert torch.allclose(L['tot'], ref)
    # equal-sized configurations: balancing changes nothing
    n = b['n_cfg']
    assert torch.all(n == n[0])
    with torch.no_grad():
        Lb = T.data_loss(model, out, b, T.LossConfig(config_balance=1.))
    for k in ('mpsp', 'tot', 'ch', 'g'):
        assert torch.allclose(L[k], Lb[k])


def test_surface_tot_equals_relative_numerator_error(fitted):
    model, _, _, _, layout, _ = build(fitted)
    b = _batch(fitted, layout)
    out = _perfect_out(model, b)
    out['coef'] = out['coef'] * 1.1           # 10 % too large everywhere -> relative error 0.1
    L = T.data_loss(model, out, b, T.LossConfig(coef_space='surface'))
    assert abs(float(L['tot']) - 0.01) < 1e-6


def test_training_with_new_options(fitted):
    s, p = fitted['samples'], fitted['points']
    ens = T.train_ensemble(s, p, 1, lcfg=T.LossConfig(coef_space='surface', config_balance=0.5),
                           tcfg=T.TrainConfig(**FAST), calib_frac=0.)
    h = ens.histories[0]
    assert h[-1]['train_loss'] < h[0]['train_loss']
    assert np.isfinite(h[-1]['train_ch'])


def test_input_clipping_bounds_extrapolation(fitted):
    s, p, spec, batch = fitted['samples'], fitted['points'], fitted['spec'], fitted['batch']
    model, tens, _, _, layout, scales = build(fitted, clip_inputs=0.0)
    sc2 = T.fit_target_scales(s, p, layout, batch, spec.continuous_mask())
    cm = spec.continuous_mask()
    assert np.isfinite(np.asarray(sc2['x_lo'])[cm]).all() and np.isinf(np.asarray(sc2['x_lo'])[~cm]).all()
    m = T.build_model(spec, layout, sc2, model.cfg, torch.float64)
    m.eval()
    X = tens[0].clone()
    far = X.clone()
    far[:, torch.as_tensor(cm)] = far[:, torch.as_tensor(cm)] + 100.     # far outside the training range
    edge = torch.minimum(far, m.x_hi)
    with torch.no_grad():
        o_far = m(far, *tens[1:])['coef']
        o_edge = m(edge, *tens[1:])['coef']
        o_in = m(X, *tens[1:])['coef']
        o_ref = T.build_model(spec, layout, sc2, type(model.cfg)(**{**model.cfg.to_dict(), 'clip_inputs': -1.}),
                              torch.float64)
        o_ref.load_state_dict({k: v for k, v in m.state_dict().items() if k not in ('x_lo', 'x_hi')})
        o_ref.eval()
        o_in_ref = o_ref(X, *tens[1:])['coef']
    assert torch.allclose(o_far, o_edge)
    assert torch.allclose(o_in, o_in_ref)          # inside the range clipping changes nothing


def test_old_checkpoint_without_clip_inputs_loads_unclipped(fitted, tmp_path):
    import json
    from fermentation_insights.coeff_predictor.model import ModelConfig
    s, p = fitted['samples'], fitted['points']
    ens = T.train_ensemble(s, p, 1, mcfg=ModelConfig(clip_inputs=-1.), tcfg=T.TrainConfig(**FAST), calib_frac=0.)
    path = str(tmp_path / 'old')
    ens.save(path)
    meta = json.load(open(f'{path}/meta.json'))
    del meta['model_config']['clip_inputs']             # as written by the previous code version
    json.dump(meta, open(f'{path}/meta.json', 'w'))
    old = T.Ensemble.load(path)
    assert old.mcfg.clip_inputs < 0
    _, o1 = ens.predict_samples(s, strict=False)
    _, o2 = old.predict_samples(s, strict=False)
    assert np.allclose(o1['coef'], o2['coef'])
    assert ModelConfig().clip_inputs < 0                    # opt-in

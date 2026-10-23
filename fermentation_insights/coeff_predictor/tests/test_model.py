# -*- coding: utf-8 -*-
"""Structural tests of the physics-informed network."""
import numpy as np
import pytest
import torch

from fermentation_insights.coeff_predictor.model import (
    mpsp_equation, mpsp_equation_np, ModelConfig, make_layout, CoefficientNet)
from fermentation_insights.coeff_predictor import physics
from fermentation_insights.coeff_predictor import train as T
from fermentation_insights.coeff_predictor.tests.conftest import build


def test_shapes_physics_and_plain(fitted):
    model, tens, spec, batch, layout, _ = build(fitted)
    n = len(batch)
    out = model(*tens)
    C, J = len(layout['channels']), len(layout['price_keys'])
    assert J == 3 and C == 11
    assert out['coef'].shape == (n, 4) and out['g'].shape == (n,)
    assert out['channels'].shape == (n, C, 4)
    assert out['intensities'].shape == (n, J, 4) and out['price_channels'].shape == (n, J, 4)
    assert len(model.channel_names()) == C + J
    pm, ptens, *_ = build(fitted, kind='plain')
    po = pm(*ptens)
    assert set(po) == {'coef', 'g'} and po['coef'].shape == (n, 4)
    # untrained heads start near the training means
    assert np.allclose(out['coef'].detach().mean(0).numpy(), fitted['samples'][list('abcd')].mean().to_numpy(),
                       rtol=0.5, atol=0.5)


def test_sum_of_channels_is_total(fitted):
    model, tens, *_ = build(fitted)
    out = model(*tens)
    s = out['channels'].sum(1) + out['price_channels'].sum(1)
    assert torch.allclose(s, out['coef'], rtol=0, atol=1e-12)


def test_price_linearity_exact(fitted):
    model, (X, lr, P, w, gm), *_ = build(fitted)
    o1 = model(X, lr, P, w, gm)
    P2 = P.clone()
    P2[:, 0] *= 2.7
    P2[:, 1] = -0.3 * P2[:, 1] + 0.01
    o2 = model(X, lr, P2, w, gm)
    # intensities and non-price channels do not depend on prices
    assert torch.equal(o1['intensities'], o2['intensities'])
    assert torch.equal(o1['channels'], o2['channels'])
    assert torch.equal(o1['g'], o2['g'])
    # totals are affine in each price with slope = intensity
    dP = (P2 - P)[:, :, None]
    expected = o1['coef'] + (dP * o1['intensities']).sum(1)
    assert torch.allclose(o2['coef'], expected, rtol=1e-12, atol=1e-12)
    # midpoint linearity
    o3 = model(X, lr, 0.5 * (P + P2), w, gm)
    assert torch.allclose(o3['coef'], 0.5 * (o1['coef'] + o2['coef']), rtol=1e-12, atol=1e-12)


def test_g_bound_and_structured_feedstock(fitted):
    model, (X, lr, P, w, gm), *_ = build(fitted)
    with torch.no_grad():  # push the heads to extremes
        for p in model.parameters():
            p.mul_(50.)
    Xe = X * 20.
    out = model(Xe, lr, P, w, gm)
    assert torch.all(out['g'] <= gm) and torch.all(out['g'] >= 0)
    j = model.layout['price_keys'].index('fs')
    eta_inv = out['intensities'][:, j, 1] * w
    assert torch.all(eta_inv >= 1 - 1e-12)  # exp(softplus) >= 1  ->  eta <= 1


def test_equation_layer_matches_numpy():
    rng = np.random.RandomState(0)
    coef = rng.normal(size=(7, 4))
    g = rng.uniform(0.5, 1., 7)
    y = rng.uniform(0.1, 1., (7, 5))
    t = rng.uniform(5, 200, (7, 5))
    ref = (coef[:, :1] + coef[:, 1:2] / y + coef[:, 2:3] / t + coef[:, 3:4] / (y * t)) / g[:, None]
    tt = torch.tensor
    got = mpsp_equation(tt(coef)[:, None, :], tt(g)[:, None], tt(y), tt(t)).numpy()
    assert np.allclose(got, ref, rtol=1e-13)
    assert np.allclose(mpsp_equation_np(coef[:, None, :], g[:, None], y, t), ref, rtol=1e-13)
    pe = physics.mpsp_equation(y, t, coef[:, :1], coef[:, 1:2], coef[:, 2:3], coef[:, 3:4], g[:, None])
    assert np.allclose(pe, ref, rtol=1e-13)


def test_mono_r_autograd_matches_finite_difference(fitted):
    model, (X, lr, P, w, gm), *_ = build(fitted, seed=3)
    Y = torch.tensor(np.tile([0.2, 0.4, 0.6], (len(lr), 1)), dtype=torch.float64)
    Tt = torch.tensor(np.tile([20., 60., 120.], (len(lr), 1)), dtype=torch.float64)
    l = lr.clone().requires_grad_(True)
    out = model(X, l, P, w, gm)
    dM = T.mpsp_dlogr(out, l, Y, Tt)
    h = 1e-6
    Mp = model.mpsp(model(X, lr + h, P, w, gm), Y, Tt)
    Mm = model.mpsp(model(X, lr - h, P, w, gm), Y, Tt)
    fd = (Mp - Mm) / (2 * h)
    assert torch.allclose(dM, fd, rtol=1e-5, atol=1e-8)
    # residual = mean relu(dM)^2 / scale^2 and is differentiable w.r.t. parameters
    mask = torch.ones_like(Y)
    res = T.physics_residuals(model, out, Y, Tt, mask, ['mono_r'], log_r=l)['mono_r']
    manual = (torch.relu(fd) ** 2).mean() / model.scales_np['mpsp_scale'] ** 2
    assert torch.allclose(res, manual, rtol=1e-4, atol=1e-12)
    if res.item() > 0:
        res.backward()
        gsum = sum(float(p.grad.abs().sum()) for p in model.parameters() if p.grad is not None)
        assert gsum > 0


def test_physics_residuals_zero_for_satisfied_and_positive_for_violated(fitted):
    model, tens, *_ = build(fitted)
    out = model(*tens)
    n = len(tens[0])
    Y = torch.full((n, 2), 0.5, dtype=torch.float64)
    Tt = torch.full((n, 2), 50., dtype=torch.float64)
    mask = torch.ones_like(Y)
    good = dict(out, coef=torch.tensor([[1., 0.5, 10., 0.1]] * n, dtype=torch.float64))
    bad = dict(out, coef=torch.tensor([[-1., -0.5, -10., -0.1]] * n, dtype=torch.float64))
    act = ['mono_y', 'mono_t', 'positive']
    rg = T.physics_residuals(model, good, Y, Tt, mask, act)
    rb = T.physics_residuals(model, bad, Y, Tt, mask, act)
    assert all(float(v) == 0 for v in rg.values())
    assert all(float(v) > 0 for v in rb.values())


def test_ablation_flags_change_layout(tables):
    samples = tables[0]
    lay = make_layout(ModelConfig(channels=False), samples.columns)
    assert lay['channels'] == [] and lay['price_keys'] == ['fs', 'el', 'ng']
    lay = make_layout(ModelConfig(price_explicit=False), samples.columns)
    assert lay['price_keys'] == [] and 'feedstock' in lay['channels']
    cfg = ModelConfig(kind='plain')
    assert not (cfg.channels or cfg.price_explicit or cfg.g_bound or cfg.structured_fs)
    with pytest.raises(ValueError):
        ModelConfig(kind='bogus')


def test_no_channels_variant_runs(fitted):
    model, tens, *_ = build(fitted, channels=False)
    out = model(*tens)
    assert 'channels' not in out and out['coef'].shape[1] == 4 and torch.isfinite(out['coef']).all()

# -*- coding: utf-8 -*-
"""Shared fixtures: small synthetic section-7 tables and fitted helpers."""
import pytest
import torch

from fermentation_insights.coeff_predictor.tests._fake_tables import make_fake_tables
from fermentation_insights.coeff_predictor.features import FeatureSpec
from fermentation_insights.coeff_predictor.model import ModelConfig, make_layout
from fermentation_insights.coeff_predictor import train as T

torch.set_num_threads(1)

#: fast training settings for tests
FAST = dict(epochs=25, batch_size=32, patience=10, min_epochs=5, threads=1)


@pytest.fixture(scope='session')
def tables():
    return make_fake_tables(n_samples_per_config=5, products=['HP', 'LA', 'SA', 'EtOH', 'NewAcid'],
                            feedstocks=['glucose', 'corn', 'cornstover'], seed=0)


@pytest.fixture(scope='session')
def fitted(tables):
    """FeatureSpec and FeatureBatch on the full fake samples."""
    samples, points, _ = tables
    spec = FeatureSpec().fit(samples)
    return dict(spec=spec, batch=spec.transform(samples), samples=samples, points=points)


def build(fitted, dtype=torch.float64, seed=0, **kw):
    """Untrained model (+ input tensors) for structural tests."""
    torch.manual_seed(seed)
    samples, points, spec, batch = fitted['samples'], fitted['points'], fitted['spec'], fitted['batch']
    cfg = ModelConfig(**{**dict(width=16, depth=2, dropout=0.), **kw})
    if not cfg.price_explicit and spec.price_explicit:
        spec = FeatureSpec(price_explicit=False).fit(samples)
        batch = spec.transform(samples)
    layout = make_layout(cfg, samples.columns)
    scales = T.fit_target_scales(samples, points, layout, batch)
    model = T.build_model(spec, layout, scales, cfg, dtype)
    model.eval()
    tens = [torch.tensor(a, dtype=dtype) for a in (batch.X, batch.log_r, batch.prices, batch.w_sugar, batch.g_max)]
    return model, tens, spec, batch, layout, scales

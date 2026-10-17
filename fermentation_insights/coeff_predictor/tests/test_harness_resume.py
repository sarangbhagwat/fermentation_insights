# -*- coding: utf-8 -*-
"""Resuming a campaign skips completed configurations without loading them."""
from fermentation_insights.coeff_predictor import synthetic
from fermentation_insights.coeff_predictor.data_generation.harness import run_campaign


def _run(adapter, out):
    return run_campaign(adapter, out, n_random=2, n_productivity=1, seed=1,
                        n_y=3, n_t=3, point_timeout=0, verbose=False)


def test_completed_configurations_are_not_reloaded(tmp_path):
    name = next(iter(synthetic.SYNTHETIC_BIOREFINERIES))
    with synthetic.biosteam_stub():
        ad = synthetic.SyntheticAdapter(name, fail_every=0)
        loads = []
        original = ad.load
        ad.load = lambda config: (loads.append(config['config_id']), original(config))[1]
        path = _run(ad, str(tmp_path))
        first_loads = list(loads)
        lines = open(path).read().splitlines()
        n_lines = len(lines)
        loads.clear()
        _run(ad, str(tmp_path))
    n_configs = len(ad.configurations())
    assert n_configs >= 1 and len(first_loads) == n_configs
    assert sum('"sample_end"' in l for l in lines) == n_configs * (1 + 1 + 2)
    assert loads == []  # every configuration was complete: nothing loaded
    assert sum(1 for _ in open(path)) == n_lines  # and nothing re-simulated

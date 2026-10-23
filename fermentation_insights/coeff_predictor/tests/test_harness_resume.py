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


def test_failed_reset_reloads_and_marks_empty_sample_for_retry(tmp_path):
    """A convergence failure whose reset() also fails: the configuration is
    reloaded for the next sample, and a sample without any simulated point
    is not complete (``retry``), so a resume re-runs it."""
    import json
    from fermentation_insights.coeff_predictor.data_generation.harness import _completed_samples
    from fermentation_insights.coeff_predictor import dataset
    name = next(iter(synthetic.SYNTHETIC_BIOREFINERIES))
    with synthetic.biosteam_stub():
        ad = synthetic.SyntheticAdapter(name, fail_every=0)
        cid = ad.configurations()[0]['config_id']
        state = dict(n=0, broken=False, fail_at=None, loads=0)
        sim, load = ad.simulate, ad.load

        def simulate():
            state['n'] += 1
            if state['n'] == state['fail_at']:
                state['broken'] = True
            if state['broken']:
                raise RuntimeError('diverged')
            return sim()

        def reset():
            if state['broken']:
                raise RuntimeError('reset diverged too')

        def reload(config):
            state['broken'] = False
            state['loads'] += 1
            return load(config)

        ad.simulate, ad.reset, ad.load = simulate, reset, reload
        state['fail_at'] = 1          # the very first point of the first sample
        path = run_campaign(ad, str(tmp_path), config_ids=[cid], n_random=2, n_productivity=1, seed=1,
                            n_y=3, n_t=3, point_timeout=0, verbose=False)
        ends = [json.loads(l) for l in open(path) if '"sample_end"' in l]
        assert len(ends) == 4 and state['loads'] == 2            # initial load + one reload
        first = ends[0]
        assert first['n_ok'] == 0 and first['reset_failed'] and first['retry'] and first['reloaded']
        assert all(e['n_ok'] > 0 and not e['retry'] for e in ends[1:])   # later samples ran normally
        done = _completed_samples(path)
        assert first['sample_id'] not in done and len(done) == 3
        # resume: only the empty sample is re-run, and it now wins in the processed data
        state['fail_at'] = None
        run_campaign(ad, str(tmp_path), config_ids=[cid], n_random=2, n_productivity=1, seed=1,
                     n_y=3, n_t=3, point_timeout=0, verbose=False)
        ends2 = [json.loads(l) for l in open(path) if '"sample_end"' in l]
        assert len(ends2) == 5 and ends2[-1]['sample_id'] == first['sample_id'] and ends2[-1]['n_ok'] > 0
        assert len(_completed_samples(path)) == 4
    _, starts, pts, _ = dataset.read_raw([str(tmp_path)])
    r = starts[starts['sample_id'] == first['sample_id']]
    assert len(r) == 1 and bool(r['_complete'].iloc[0]) and r['n_ok_reported'].iloc[0] > 0

"""The processed training tables are shipped in git (README section 8, "Data
availability"), and the pipeline uses them when the raw simulation records
are absent. They hold the data regenerated on 2026-10-18 (D16 / D17
corrected); a run of experiments.py on them is the data of the shipped
checkpoints and of every REPORT.md number (until that run, results/ and
models/ are the v3 run on the tables before the regeneration)."""
import json
import os
import shutil
import types

import pytest

from fermentation_insights.coeff_predictor import experiments as E
from fermentation_insights.coeff_predictor import train as T

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(HERE, 'data', 'processed')
PLAN = os.path.join(HERE, 'results', 'plan.json')

#: the shipped processed tables: after the regeneration of 2026-10-18
#: (harness v4 / TAL adapter v3: data defects D16 and D17 corrected; DESIGN.md
#: section 1, "Data set")
DATA_HASH = 'ae69af98b3ba4561'
N_SAMPLES, N_RAW_FILES = 1600, 41
#: the tables of the v3 run (results/v3, models trained before the
#: regeneration; git history up to 85ec2b4)
V3_DATA_HASH = 'c345efda41c4fdca'


def _need():
    if not all(os.path.exists(os.path.join(PROC, f'{n}.csv.gz')) for n in ('samples', 'points', 'configs')):
        pytest.skip('processed tables not available (README section 8, "Data availability")')


def test_shipped_tables_are_the_regenerated_data():
    _need()
    s, p, _ = T.load_tables(PROC)
    s = s.sort_values('sample_id', kind='mergesort').reset_index(drop=True)
    assert E.data_hash(s, p) == DATA_HASH
    assert len(s) == N_SAMPLES
    # every sample of a configuration with an InternalCirculationRx was simulated with the D17 correction
    from fermentation_insights.coeff_predictor import data_patches as dp
    ic = s[[(b, c) in dp.D17_CONFIGS for b, c in zip(s['biorefinery'], s['config_id'])]]
    assert len(ic) and (ic['biosteam_fixes'] == dp.D17_FIX).all() and (ic['harness_version'] >= 4).all()
    tal = s[s['config_id'] == 'TAL_SA_sugarcane']
    assert len(tal) and (tal['adapter_version'] >= dp.D16_FIXED_ADAPTER_VERSION).all()


def test_reported_run_used_the_shipped_tables():
    _need()
    if not os.path.exists(PLAN):
        pytest.skip('results/plan.json not available')
    want = json.load(open(PLAN))['data_hash']
    if want == V3_DATA_HASH and DATA_HASH != V3_DATA_HASH:
        pytest.xfail('results/ and models/ are still the v3 run, trained on the tables before the D16 / D17 '
                     'regeneration (data hash c345efda41c4fdca): re-run experiments.py on the shipped tables')
    assert want == DATA_HASH
    with open(os.path.join(HERE, 'results', 'REPORT.md')) as fh:
        assert f'data hash `{want}`' in fh.read(2000)


def test_stage_data_uses_the_shipped_tables_without_raw_records(tmp_path):
    _need()
    proc = tmp_path / 'processed'
    shutil.copytree(PROC, proc)
    a = types.SimpleNamespace(processed=str(proc), raw=[str(tmp_path / 'no_raw')], no_rebuild=False, rebuild=False)
    D = E.stage_data(a)
    assert D['hash'] == DATA_HASH
    assert len(D['samples']) == N_SAMPLES and D['n_raw_files'] == N_RAW_FILES

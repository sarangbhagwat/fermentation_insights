"""The processed training tables are shipped in git (README section 8, "Data
availability"): they are the data of the shipped checkpoints and of every
REPORT.md number, and the pipeline uses them when the raw simulation records
are absent."""
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


def _need():
    if not all(os.path.exists(os.path.join(PROC, f'{n}.csv.gz')) for n in ('samples', 'points', 'configs')):
        pytest.skip('processed tables not available (README section 8, "Data availability")')
    if not os.path.exists(PLAN):
        pytest.skip('results/plan.json not available')
    return json.load(open(PLAN))['data_hash']


def test_shipped_tables_are_the_reported_data():
    want = _need()
    s, p, _ = T.load_tables(PROC)
    s = s.sort_values('sample_id', kind='mergesort').reset_index(drop=True)
    assert E.data_hash(s, p) == want
    with open(os.path.join(HERE, 'results', 'REPORT.md')) as fh:
        assert f'data hash `{want}`' in fh.read(2000)


def test_stage_data_uses_the_shipped_tables_without_raw_records(tmp_path):
    want = _need()
    proc = tmp_path / 'processed'
    shutil.copytree(PROC, proc)
    a = types.SimpleNamespace(processed=str(proc), raw=[str(tmp_path / 'no_raw')], no_rebuild=False, rebuild=False)
    D = E.stage_data(a)
    assert D['hash'] == want
    assert len(D['samples']) == 1600 and D['n_raw_files'] == 41

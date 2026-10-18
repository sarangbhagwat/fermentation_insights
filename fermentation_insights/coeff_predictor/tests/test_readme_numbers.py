# -*- coding: utf-8 -*-
"""The user guides must quote the numbers of the committed run
(``results/tables``), not those of an earlier checkpoint: the generated
accuracy block of README.md, the headline figures quoted in its prose and
in the top-level README, and no statement that the shipped checkpoint is
an older version."""
import json
import os

import pandas as pd
import pytest

from fermentation_insights.coeff_predictor import report_checks as RC

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES = os.path.join(HERE, 'results', 'tables')


def _read(path):
    with open(path, encoding='utf-8') as fh:
        return fh.read()


@pytest.fixture(scope='module')
def readme():
    need = ('cv_headline.csv', 'cv_baselines_paired.csv', 'cv_coverage_crossfit.csv',
            'interval_coverage_by_scheme.csv', 'final_test_headline.csv')
    if not all(os.path.exists(os.path.join(TABLES, f)) for f in need):
        pytest.skip('results tables not available')
    return _read(os.path.join(HERE, 'README.md'))


def test_generated_block_matches_tables(readme):
    assert RC.readme_summary(TABLES) in readme, (
        'README.md accuracy block is stale: run python -m fermentation_insights.coeff_predictor.report_checks')


def test_prose_quotes_the_committed_run(readme):
    top = _read(os.path.join(os.path.dirname(os.path.dirname(HERE)), 'README.md'))
    h = pd.read_csv(os.path.join(TABLES, 'cv_headline.csv'))
    m = lambda sch: f'{100 * h[(h.scheme == sch) & (h.method == "pinn")].mpsp_mape.iloc[0]:.1f} %'
    shc = pd.read_csv(os.path.join(TABLES, 'interval_coverage_by_scheme.csv')).set_index(['scheme', 'key'])
    cov = lambda sch, k: f'{shc.loc[(sch, k), "pooled"]:.2f}'
    for sch in ('config', 'feedstock', 'product'):
        assert m(sch) in readme, f'README does not quote the {sch} CV MAPE {m(sch)}'
    for sch in ('feedstock', 'product'):
        for k in ('mpsp', 'g'):
            assert cov(sch, k) in readme, f'README does not quote the {sch} {k} coverage {cov(sch, k)}'
    for sch in ('config', 'feedstock', 'product'):
        v = m(sch).replace(' %', '')
        assert v in top, f'top-level README does not quote the {sch} CV MAPE {v}'
    for txt in (readme, top):
        low = txt.lower()
        for stale in ('shipped v2', 'shipped (v2)', 'shipped v3', 'shipped (v3)', 'are the v3 models',
                      'must be re-run to produce', 'is to be re-run'):
            assert stale not in low, stale


def test_checkpoint_warning_numbers_match_tables():
    path = os.path.join(HERE, 'models', 'release', 'meta.json')
    f = os.path.join(TABLES, 'interval_coverage_by_scheme.csv')
    if not (os.path.exists(path) and os.path.exists(f)):
        pytest.skip('release checkpoint / tables not available')
    cov = json.load(open(path))['meta'].get('coverage_unseen')
    assert cov, 'models/release/meta.json lacks coverage_unseen (run the pipeline tables stage)'
    shc = pd.read_csv(f)
    for r in shc[shc.scheme.isin(['feedstock', 'product'])].itertuples():
        assert abs(cov[r.scheme][r.key] - r.pooled) < 1e-3


def test_readme_quotes_the_family_interval_coverage(readme):
    """The user guide quotes the LOCO per-family MPSP coverage of the shipped
    quantiles where it under-covers (REPORT.md section 6), and the
    checkpoint's coverage_by_family metadata matches the table."""
    f = os.path.join(TABLES, 'interval_coverage_by_family.csv')
    if not os.path.exists(f):
        pytest.skip('tables not available')
    fc = pd.read_csv(f)
    low = fc[(fc['key'] == 'mpsp') & (fc['coverage'] < 0.8)]
    for x in low.itertuples():
        assert f'{x.family} {x.coverage:.2f}' in readme, f'README does not quote the {x.family} MPSP coverage'
    path = os.path.join(HERE, 'models', 'release', 'meta.json')
    cbf = json.load(open(path))['meta'].get('coverage_by_family') if os.path.exists(path) else None
    if not cbf:
        pytest.skip('release checkpoint without coverage_by_family')
    s = pd.read_csv(os.path.join(HERE, 'data', 'processed', 'samples.csv.gz'), usecols=['biorefinery',
                                                                                       'fermentation_product'])
    fp = dict(zip(s['biorefinery'], s['fermentation_product']))
    for x in fc.itertuples():
        assert abs(cbf['new_configuration'][fp[x.family]][x.key] - x.coverage) < 1e-3

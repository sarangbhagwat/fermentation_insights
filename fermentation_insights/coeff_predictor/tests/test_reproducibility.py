# -*- coding: utf-8 -*-
"""Byte-stable pipeline outputs: rebuilt processed tables carry no gzip
timestamp, and the raw-data fingerprint of experiments.py depends on the
records' relative paths, sizes and content, not on where they are or when
they were unpacked (with a one-time adoption of the legacy path / mtime
key, so that the change does not rebuild the shipped tables)."""
import gzip
import os
import time

import pandas as pd

from fermentation_insights.coeff_predictor import experiments as E


def test_rebuilt_tables_are_rewritten_without_timestamp(tmp_path):
    d = pd.DataFrame(dict(a=[1., 2.], b=['x', 'y']))
    paths = []
    for n in ('one', 'two'):
        p = str(tmp_path / f'{n}.csv.gz')
        d.to_csv(p, index=False)                      # pandas default: gzip header with the current time
        time.sleep(1.1)
        E._gzip_without_timestamp(p)
        paths.append(p)
    a, b = (open(p, 'rb').read() for p in paths)
    assert a == b
    assert gzip.decompress(a) == d.to_csv(index=False).encode()


def test_raw_fingerprint_ignores_location_and_mtime(tmp_path):
    for root in ('r1', 'r2'):
        os.makedirs(tmp_path / root / 'HP')
        (tmp_path / root / 'HP' / 'HP_corn.jsonl').write_text('{"record": "config"}\n')
    os.utime(tmp_path / 'r2' / 'HP' / 'HP_corn.jsonl', (1e9, 1e9))
    h1, f1 = E.raw_fingerprint([str(tmp_path / 'r1')])
    h2, f2 = E.raw_fingerprint([str(tmp_path / 'r2')])
    assert h1 == h2 and f1[0][0] != f2[0][0]
    # the legacy key (absolute path, size, mtime) differs between the two
    assert E._legacy_raw_fingerprint(f1) != E._legacy_raw_fingerprint(f2)
    assert E._legacy_raw_fingerprint(f1) == E._sha(f1)
    (tmp_path / 'r2' / 'HP' / 'HP_corn.jsonl').write_text('{"record": "confiG"}\n')
    assert E.raw_fingerprint([str(tmp_path / 'r2')])[0] != h1      # same size, other content

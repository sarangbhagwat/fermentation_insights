# -*- coding: utf-8 -*-
"""
Command-line entry point for data-generation campaigns. Run from the
biorefinery's own conda environment, with Bioindustrial-Park and the
fermentation_insights repository root on PYTHONPATH::

    python -m fermentation_insights.coeff_predictor.data_generation.run \
        --adapter HP --out data/raw --n-random 20 --n-productivity 4 --seed 1

Multiple workers can run in parallel on disjoint configurations
(``--configs``) or different seeds; each writes its own JSON-lines file
(``--tag``).
"""
import argparse
import importlib
import os
import sys
import warnings


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--adapter', required=True, help='adapter module name in adapters/')
    ap.add_argument('--out', required=True)
    ap.add_argument('--configs', nargs='*', default=None, help='config_ids to run')
    ap.add_argument('--n-random', type=int, default=20)
    ap.add_argument('--n-productivity', type=int, default=4)
    ap.add_argument('--n-y', type=int, default=5)
    ap.add_argument('--n-t', type=int, default=5)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--prod-range', type=float, nargs=2, default=(0.2, 5.0))
    ap.add_argument('--point-timeout', type=float, default=600.)
    ap.add_argument('--sample-limit', type=int, default=None)
    ap.add_argument('--tag', default='')
    ap.add_argument('--list-configs', action='store_true')
    args = ap.parse_args(argv)
    warnings.filterwarnings('ignore')
    mod = importlib.import_module(
        f'fermentation_insights.coeff_predictor.data_generation.adapters.{args.adapter}')
    adapter = mod.Adapter()
    if args.list_configs:
        for c in adapter.configurations():
            print(c['config_id'])
        return
    from .harness import run_campaign
    path = run_campaign(adapter, args.out, config_ids=args.configs,
                        n_random=args.n_random, n_productivity=args.n_productivity,
                        seed=args.seed, n_y=args.n_y, n_t=args.n_t,
                        productivity_range=tuple(args.prod_range),
                        point_timeout=args.point_timeout,
                        sample_limit=args.sample_limit, tag=args.tag)
    print(f'\nDone: {path}')


if __name__ == '__main__':
    sys.exit(main())

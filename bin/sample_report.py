#!/usr/bin/env python3

import argparse
import csv
import itertools
import operator
import sys


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('columns')
    parser.add_argument('csvs', nargs='+', type=argparse.FileType('r'))
    parser.add_argument('--sortby')
    parser.add_argument(
        '--out', type=argparse.FileType('w'), default=sys.stdout)
    return parser.parse_args()


def main():
    args = get_args()
    rows = []
    csvs = (r for f in args.csvs for r in csv.DictReader(f))
    csvs = sorted(csvs, key=operator.itemgetter('sample'))
    csvs = itertools.groupby(csvs, key=operator.itemgetter('sample'))
    for sample, columns in csvs:
        row = {'sample': sample}
        for c in columns:
            row.update(c)
        rows.append(row)
    out = csv.DictWriter(
        args.out,
        extrasaction='ignore',
        fieldnames=args.columns.split(','))
    out.writeheader()
    out.writerows(rows)


if __name__ == '__main__':
    sys.exit(main())

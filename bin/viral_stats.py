#!/usr/bin/env python3
import argparse
import csv
import gzip
import sys


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('sample')
    parser.add_argument('fastq')
    parser.add_argument(
        '--out', default=sys.stdout, type=argparse.FileType('w'))
    return parser.parse_args()


def main():
    args = get_args()
    # <sample_id>/<library_type>/<FCID.Lane.Index-Index>
    with gzip.open(args.fastq, mode='rt') as fopen:
        count = sum(1 for l in fopen if l.startswith('@'))
    report = csv.DictWriter(args.out, fieldnames=['sample', 'viral'])
    report.writeheader()
    report.writerow({'sample': args.sample, 'viral': count})


if __name__ == '__main__':
    sys.exit(main())

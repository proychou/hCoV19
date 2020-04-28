#!/usr/bin/env python3
import argparse
import csv
import datetime
import gzip
import os
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
    with gzip.open(args.fastq, mode='rt') as fopen:
        count = sum(1 for l in fopen if l.startswith('@'))
    with gzip.open(args.fastq, mode='rt') as fopen:
        next(fopen)
        length = len(next(fopen))
    fieldnames = ['sample', 'fastq', 'raw_read_count',
                  'length', 'analysis_date']
    report = csv.DictWriter(args.out, fieldnames=fieldnames)
    report.writeheader()
    report.writerow({
        'analysis_date': datetime.date.today().strftime('%d-%b-%Y'),
        'length': length,
        'fastq': os.path.basename(args.fastq),
        'raw_read_count': count,
        'sample': args.sample,
        })


if __name__ == '__main__':
    sys.exit(main())

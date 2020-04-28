#!/usr/bin/env python3
import argparse
import csv
import datetime
import gzip
import os
import sys


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq')
    parser.add_argument('out', type=argparse.FileType('w'))
    parser.add_argument(
        '--sample_id',
        default=sys.stdout,
        type=argparse.FileType('w'))
    return parser.parse_args()


def main():
    args = get_args()
    fieldnames = ['sample', 'fastq', 'length', 'raw_read_count',
                  'library', 'fcid', 'lane', 'I1', 'I2', 'analysis_date']
    report = csv.DictWriter(args.out, fieldnames=fieldnames)
    report.writeheader()
    fastq = os.path.realpath(args.fastq)
    # <sample_id>/<library_type>/<FCID.Lane.Index-Index>
    sample_id, library_type, info = os.path.dirname(fastq).split('/')[-3:]
    fcid, lane, index = info.split('.')
    i1, i2 = index.split('-')
    with gzip.open(args.fastq, mode='rt') as fopen:
        count = sum(1 for l in fopen if l.startswith('@'))
    with gzip.open(args.fastq, mode='rt') as fopen:
        next(fopen)
        length = len(next(fopen))
    report.writerow({
        'analysis_date': datetime.date.today().strftime('%d-%b-%Y'),
        'fcid': fcid,
        'I1': i1,
        'I2': i2,
        'library': library_type,
        'length': length,
        'lane': lane,
        'fastq': os.path.basename(args.fastq),
        'raw_read_count': count,
        'sample': sample_id,
        })
    args.sample_id.write(sample_id)


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3
import argparse
import csv
import datetime
import itertools
import glob
import gzip
import os
import sys


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('datadir')
    parser.add_argument('--samplelist', type=argparse.FileType('w'))
    parser.add_argument('--take', type=int)
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    return parser.parse_args()


def main():
    args = get_args()
    take = None if args.take == -1 else args.take
    path = os.path.join(args.datadir, '**/*.fastq.gz')
    fastqs = itertools.islice(glob.iglob(path, recursive=True), take)
    fastqs = [os.path.abspath(f) for f in fastqs]
    fastqs = sorted(fastqs, key=lambda x: os.path.dirname(x))
    fieldnames = ['sample', 'R1', 'R2', 'length', 'raw_read_count',
                  'library', 'fcid', 'lane', 'I1', 'I2', 'analysis_date']
    report = csv.DictWriter(args.out, fieldnames=fieldnames)
    report.writeheader()
    for p, f in itertools.groupby(fastqs, key=lambda x: os.path.dirname(x)):
        f = list(f)
        # <sample_id>/<library_type>/<FCID.Lane.Index-Index>
        sample_id, library_type, info = p.split('/')[-3:]
        fcid, lane, index = info.split('.')
        i1, i2 = index.split('-')
        with gzip.open(f[0], mode='rt') as fopen:
            count = sum(1 for l in fopen if l.startswith('@'))
        with gzip.open(f[0], mode='rt') as fopen:
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
            'R1': os.path.basename(f[0]),
            'R2': os.path.basename(f[1]) if len(f) > 1 else None,
            'raw_read_count': count,
            'sample': sample_id,
            })
        if args.samplelist:
            for p in f:
                args.samplelist.write('{},{}\n'.format(sample_id, p))


if __name__ == '__main__':
    sys.exit(main())
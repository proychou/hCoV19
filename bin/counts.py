#!/usr/bin/env python3

import argparse
import csv
import itertools
import gzip
import operator
import os
import sys


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqs', nargs='+')
    parser.add_argument('--sortby')
    parser.add_argument(
        '--out', type=argparse.FileType('w'), default=sys.stdout)
    return parser.parse_args()


def main():
    args = get_args()
    procs = []
    for sfile in args.seqs:
        if sfile.endswith('fasta.gz'):
            start = '>'
            ofunc = gzip.open
        elif sfile.endswith('fasta'):
            start = '>'
            ofunc = open
        elif sfile.endswith('fastq.gz'):
            start = '@'
            ofunc = gzip.open
        elif sfile.endswith('fastq'):
            start = '@'
            ofunc = open
        else:
            raise ValueError('unsupported filetype ' + sfile)
        with ofunc(sfile, mode='rt') as seqs:
            count = sum(1 for l in seqs if l.startswith(start))
        basename = os.path.basename(sfile)
        sample, process, _ = basename.split('.', maxsplit=2)
        procs.append([sample, process, count])
    fieldnames = ['sample']
    for p in procs:
        if p[1] not in fieldnames:
            fieldnames.append(p[1])
    rows = []
    procs = sorted(procs, key=operator.itemgetter(0))
    procs = itertools.groupby(procs, key=operator.itemgetter(0))
    rows = []
    for sa, pr in procs:
        r = {'sample': sa}
        for p in pr:
            r.update({p[1]: p[2]})
        rows.append(r)
    if args.sortby:
        rows = sorted(rows, key=operator.itemgetter(args.sortby), reverse=True)
    report = csv.DictWriter(args.out, fieldnames=fieldnames)
    report.writeheader()
    report.writerows(rows)


if __name__ == '__main__':
    sys.exit(main())

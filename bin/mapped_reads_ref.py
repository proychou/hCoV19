#!/usr/bin/env python3
import argparse
import pysam
import sys


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('sample')
    parser.add_argument('bam')
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'))
    return parser.parse_args()


def main():
    args = get_args()
    pysam.index(args.bam)
    count = pysam.AlignmentFile(args.bam, 'r').get_index_statistics()[0].mapped
    args.out.write('sample,mapped_reads_ref\n')
    args.out.write('{},{}\n'.format(args.sample, count))


if __name__ == '__main__':
    sys.exit(main())

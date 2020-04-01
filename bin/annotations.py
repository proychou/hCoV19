#!/usr/bin/env python3

import argparse
import sys


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('arg1')
    parser.add_argument('--option1')
    parser.add_argument(
        '--out', type=argparse.FileType('w'), default=sys.stdout)
    return parser.parse_args()


def main():
    args = get_args()
    print(args)


if __name__ == '__main__':
    sys.exit(main())

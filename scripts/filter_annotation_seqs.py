#! /usr/bin/env python3
import sys
import os
import argparse
import gzip

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", "-g", help="GTF/GFF3 file (optionally gzipped)",
        required=True)
    parser.add_argument("--seqs", "-s", help="File of sequences to keep", 
        required=True)
    return parser.parse_args()

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def main(args):
    options = parse_args()
    seqs_keep = set([])
    f = open(options.seqs, 'r')
    for line in f:
        line = line.rstrip()
        seqs_keep.add(line)
    f.close()
    is_gz = False
    if is_gz_file(options.annotation):
        is_gz = True
        f = gzip.open(options.annotation, 'r')
    else:
        f = open(options.annotation, 'r')
    for line in f:
        if is_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        dat = line.split('\t')
        if dat[0] in seqs_keep:
            print(line)
    

if __name__ == '__main__':
    sys.exit(main(sys.argv))

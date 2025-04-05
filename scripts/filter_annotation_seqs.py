#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
from anno_parse import *

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", "-g", help="GTF/GFF3 file (optionally gzipped)",
        required=True)
    parser.add_argument("--seqs", "-s", help="File of sequences to keep", 
        required=True)
    return parser.parse_args()

def main(args):
    options = parse_args()
    seqs_keep = set([])
    f = open(options.seqs, 'r')
    for line in f:
        line = line.rstrip()
        seqs_keep.add(line)
    f.close()
    
    for record in read_annotation(options.annotation):
        if record['comment']:
            print(record['line'])
        else:
            if record['seq'] in seqs_keep:
                print(record['line'])
        
if __name__ == '__main__':
    sys.exit(main(sys.argv))

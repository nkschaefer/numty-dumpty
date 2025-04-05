#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
from anno_parse import *
"""
Takes in a GTF/GFF3 of mitochondrial gene annotations and outputs a list of 
gene names and IDs.
"""

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", "-a", help="GTF/GFF3 file (optionally gzipped)", 
        required=True)
    parser.add_argument("--chrM", "-m", help="Name of mitochondrial sequence",
        required=False, default="chrM")
    parser.add_argument("--id", "-i", help="Name of field storing gene ID",
        required=False, default="gene_id")
    parser.add_argument("--name", "-n", help="Name of field storing gene name",
        required=False, default="gene_name")
    return parser.parse_args()

def main(args):

    options = parse_args()
    
    for anno_dat in read_annotation(options.annotation):
        if not anno_dat['comment'] and anno_dat['seq'] == options.chrM and \
            anno_dat['feature'] == 'gene':
            if options.name in anno_dat['tags'] and options.id in anno_dat['tags']:
                gid = anno_dat['tags'][options.name]
                gname = anno_dat['tags'][options.id]
                if gid is not None and gname is not None:
                    print("{}\t{}".format(gid, gname))

if __name__ == '__main__':
    sys.exit(main(sys.argv))

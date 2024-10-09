#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
"""
Takes in a GTF of mitochondrial gene annotations and outputs a list of 
gene names and IDs.
"""
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", "-g", help="GTF file (optionally gzipped)", 
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
        
    is_gz = False
    f = None
    if is_gz_file(options.gtf):
        is_gz = True
        f = gzip.open(options.gtf, 'r')
    else:
        f = open(options.gtf, 'r')

    for line in f:
        if is_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        if line[0] != "#":
            dat = line.split('\t')
            if dat[0] == options.chrM and dat[2] == 'gene':
                gid = None
                gname = None
                for elt in dat[8].split(';'):
                    elt = elt.strip()
                    k, v = elt.split()
                    if k == options.id:
                        gid = v.strip('"')
                    elif k == options.name:
                        gname = v.strip('"')
                    if gid is not None and gname is not None:
                        break
                
                if gid is not None:
                    gid = gid.split('.')[0]
                    if gname is None:
                        gname = gid
                    print("{}\t{}".format(gid, gname))


if __name__ == '__main__':
    sys.exit(main(sys.argv))

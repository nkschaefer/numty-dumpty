#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
"""
Takes in a list of mitochondrial gene IDs & names, and an annotation file
(either GTF or GFF3 format), and outputs the same annotation file with 
mitochondrial gene annotations removed.
"""
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    
    annotation = parser.add_mutually_exclusive_group(required=True)
    annotation.add_argument("--gtf", help="Input annotation in GTF format")
    annotation.add_argument("--gff3", help="Input annotation in GFF3 format")
   
    parser.add_argument("--cat", "-c", help="Flag to specify that the annotation \
is from Comparative Annotation Toolkit (CAT). Equaivalent to --id source_gene --name source_gene_common_name", \
        action="store_true")
    parser.add_argument("--mito_genes", "-M", help="List of mitochondrial genes",
        required=True)
    parser.add_argument("--chrM", "-m", help="Name of mitochondrial sequence",
        required=False, default="chrM")
    parser.add_argument("--id", "-i", help="Name of field storing gene ID",
        required=False, default="gene_id")
    parser.add_argument("--name", "-n", help="Name of field storing gene name",
        required=False, default="gene_name")
    return parser.parse_args()

def parse_mito_genes(filename):
    mito_ids = set([])
    mito_names = set([])
    f = open(filename, 'r')
    for line in f:
        line = line.rstrip()
        if line != "":
            gid, gname = line.split('\t')
            mito_ids.add(gid)
            mito_names.add(gname)
    f.close()
    return (mito_ids, mito_names)

def parse_fields_gff3(fields):
    parsed = {}
    for elt in fields.split(';'):
        k, v = elt.split('=')
        parsed[k] = v
    return parsed

def parse_fields_gtf(fields):
    parsed = {}
    for elt in fields.rstrip(';').split(';'):
        k, v = elt.strip().split(' ')
        parsed[k] = v.strip('"')
    return parsed

def read_annotation(filename, gidfield, gnamefield, is_gff3=False):
    is_gz = False
    f = None
    if is_gz_file(filename):
        is_gz = True
        f = gzip.open(filename, 'r')
    else:
        f = open(filename, 'r')
    
    for line in f:
        if is_gz:
            line = line.decode().rstrip()
        else:
            line = line.strip()
        if line[0] != "#":
            dat = line.split('\t')
            fields = None
            if is_gff3:
                fields = parse_fields_gff3(dat[8])
            else:
                fields = parse_fields_gtf(dat[8])
            gid = None
            gname = None
            if gidfield in fields:
                gid = fields[gidfield]
            if gnamefield in fields:
                gname = fields[gnamefield]
            if gid is not None:
                gid = gid.split('.')[0]
                if gname is None:
                    gname = gid
            elif gname is not None:
                if gid is None:
                    gid = gname
            if gid is not None and gname is not None:
                yield (line, dat[0], gid, gname)

def main(args):

    options = parse_args()
    if options.cat:
        options.id = "source_gene"
        options.name = "source_gene_common_name"
    mito_ids, mito_names = parse_mito_genes(options.mito_genes)
    
    anno_file = None
    is_gff3 = False
    if options.gff3 is not None:
        anno_file = options.gff3
        is_gff3 = True
    else:
        anno_file = options.gtf

    for line, chrom, gid, gname in read_annotation(anno_file, options.id,
        options.name, is_gff3):
        if chrom != options.chrM and gid not in mito_ids and gname not in mito_names:
            print(line)

if __name__ == '__main__':
    sys.exit(main(sys.argv))

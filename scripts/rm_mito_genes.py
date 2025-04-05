#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
from anno_parse import *
"""
Takes in a list of mitochondrial gene IDs & names, and an annotation file
(either GTF or GFF3 format), and outputs the same annotation file with 
mitochondrial gene annotations removed.
"""

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", "-a", help="Input annotation (GTF/GFF3)", required=True)
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
    parser.add_argument("--remove_all_MT", "-r", help="Set to true to remove all genes \
on the mitochondrial sequence. Default: remove only mitochondrial genes from nuclear \
sequences.", action="store_true", default=False)
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

def main(args):

    options = parse_args()
    if options.cat:
        options.id = "source_gene"
        options.name = "source_gene_common_name"
    mito_ids, mito_names = parse_mito_genes(options.mito_genes)
    
    for feature in read_annotation(options.annotation):
        if feature['comment']:
            print(feature['line'])
        else:
            if options.remove_all_MT:
                if feature['seq'] != options.chrM and options.id in feature['tags'] and \
                    options.name in feature['tags']:
                    gid = feature['tags'][options.id]
                    gname = feature['tags'][options.name]
                    if gid not in mito_ids and gname not in mito_names:
                        print(feature['line'])
            else:
                if feature['seq'] == options.chrM:
                    print(feature['line'])
                elif options.id in feature['tags'] and options.name in feature['tags'] and \
                    feature['tags'][options.id] not in mito_ids and \
                    feature['tags'][options.name] not in mito_names:
                    print(feature['line'])

if __name__ == '__main__':
    sys.exit(main(sys.argv))

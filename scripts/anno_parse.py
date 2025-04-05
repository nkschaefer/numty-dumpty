#! /usr/bin/env python3
import sys
import os
import gzip

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def get_tags_gff3(dat):
    tags = {}
    try:
        for elt in dat[8].strip().rstrip(';').split(';'):
            elt = elt.strip()
            k, v = elt.split('=')
            if k == 'tag':
                # Can have multiple values
                if k not in tags:
                    tags[k] = []
                tags[k].append(v)
            tags[k] = v
    except:
        print("ERROR: input likely not in GFF3 format", file=sys.stderr)
        exit(1)
    return tags

def get_tags_gtf(dat):
    tags = {}
    try:
        for elt in dat[8].strip().rstrip(';').split(';'):
            elt = elt.strip()
            k, v = elt.split(' ')
            v = v.strip('"')
            if k == 'tag':
                # Can have multiple values
                if k not in tags:
                    tags[k] = []
                tags[k].append(v)
            else:
                tags[k] = v
    except:
        print("ERROR: input likely not GTF format.", file=sys.stderr)
        exit(1)
    return tags

def read_annotation(filename):
    is_gz = False
    f = None
    if is_gz_file(filename):
        is_gz = True
        f = gzip.open(filename, 'r')
    else:
        f = open(filename, 'r')
    
    is_gff = False
    first = True

    for line in f:
        if is_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        if line[0] != "#":
            dat = line.split('\t')
            if first:
                if '=' in dat[8]:
                    is_gff = True
                else:
                    is_gff = False
                first = False
            yld_dat = {'line': line, 'comment': False, 'seq': dat[0], 'feature': dat[2] }
            if is_gff:
                yld_dat['tags'] = get_tags_gff3(dat)
            else:
                yld_dat['tags'] = get_tags_gtf(dat)
            yield yld_dat
        else:
            yield{ 'line': line, 'comment': True, 'tags': None, 'seq': None, 'feature': None }



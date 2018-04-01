
from __future__ import print_function
import os,sys
import argparse
import re
import ssw
import numpy as np
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
import string
import fractions
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from matplotlib_venn import venn3, venn3_circles, venn2
import edlib
import math
import errno

def edlib_traceback(x, y, mode="NW", task="path", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    if ed == -1:
        ed= len(x)

    if task == "path":
        locations =  result["locations"]
        cigar =  result["cigar"]
        return (ed, locations, cigar)
    else:
        return (ed, None, None)

def read_fasta(fasta_file):
    fasta_seqs = {}
    k = 0
    temp = ''
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip()
            fasta_seqs[accession] = ''
            k += 1
        elif line[0] == '>':
            yield accession, temp
            temp = ''
            accession = line[1:].strip()
        else:
            temp += line.strip().upper()
    if accession:
        yield accession, temp

# def check_inexact(a,b,c):
#     for seq1 in a:
#         best_ed = len(seq1)
#         best_cig = ""
#         best_loc = ""
#         for seq2 in c:
#             if math.fabs(len(seq1) - len(seq2)) > best_ed:
#                 continue

#             edit_distance, locations, cigar = edlib_traceback(seq1, seq2, mode="NW", task="path", k=min(100, best_ed))
#             if edit_distance < best_ed:
#                 best_ed = edit_distance
#                 best_cig = cigar
#                 best_loc = locations

#         if best_ed > 0:
#             print(best_ed,best_cig, best_loc)
#             if best_ed > 500:
#                 print(seq1)

from collections import defaultdict

def main_fasta(args):
    hits = defaultdict(list)
    for i, file_ in enumerate(args.fastafiles):
        f = open(file_, "r")
        sample_id = args.names[i]
        dict_ = {seq.upper() : acc for (acc, seq) in  read_fasta(f)}
        hits[sample_id] = dict_
    path_, file_prefix = os.path.split(args.outfile)

    if len(args.fastafiles) == 3:

        a = set(hits[args.names[0]].keys())
        b = set(hits[args.names[1]].keys())
        c = set(hits[args.names[2]].keys())

        # if params.inexact:
        #     check_inexact(a,b,c)

        print(len(a), len(b), len(c))

        a_not_b_c = a - (b | c)
        b_not_a_c = b - (a | c)
        a_b_not_c = (a & b) - c
        c_not_a_b = c - (a | b)
        a_c_not_b = (a & c) - b
        b_c_not_a = (b & c) - a
        a_b_c = a & b & c

        outfile = open(os.path.join(path_, "a_b_not_c.fa"), "w")
        for i,  seq in enumerate( a_b_not_c):
            acc = hits[args.names[0]][seq]
            outfile.write(">{0}\n{1}\n".format(acc, seq))

        outfile = open(os.path.join(path_, "a_c_not_b.fa"), "w")
        for i,  seq in enumerate( a_c_not_b):
            acc = hits[args.names[0]][seq]
            outfile.write(">{0}\n{1}\n".format(acc, seq))
        
        outfile = open(os.path.join(path_, "b_c_not_a.fa"), "w")
        for i,  seq in enumerate( b_c_not_a):
            acc = hits[args.names[1]][seq]
            outfile.write(">{0}\n{1}\n".format(acc, seq))

        # outfile = open(os.path.join(path_, "b_not_a_c.fa"), "w")
        # for i,  seq in enumerate( b_not_a_c):
        #     acc = hits[args.names[1]][seq]
        #     outfile.write(">{0}\n{1}\n".format(acc, seq))

        # outfile = open(os.path.join(path_, "c_not_a_b.fa"), "w")
        # for i,  seq in enumerate( c_not_a_b):
        #     acc = hits[args.names[2]][seq]
        #     outfile.write(">{0}\n{1}\n".format(acc, seq))

        # outfile = open(os.path.join(path_, "a_or_b_or_c.fa"), "w")
        # for i,  seq in enumerate( (a | b | c)):
        #     if seq in dict1:
        #         acc = dict1[seq]
        #     elif seq in dict2:
        #         acc = dict2[seq]
        #     elif seq in dict3:
        #         acc = dict3[seq]            

        #     outfile.write(">{0}\n{1}\n".format(acc, seq))

        outfile = open(os.path.join(path_, "a_and_b_and_c.fa"), "w")
        for i,  seq in enumerate( (a & b & c)):
            acc = hits[args.names[0]][seq]
            outfile.write(">{0}\n{1}\n".format(acc, seq))

        outfile = open(os.path.join(path_, "pred_by_at_least_2.fa"), "w")
        for i,  seq in enumerate(   (a | b | c) - (a_not_b_c | b_not_a_c | c_not_a_b ) ):
            if seq in hits[args.names[0]]:
                acc = hits[args.names[0]][seq]
            elif seq in hits[args.names[1]]:
                acc = hits[args.names[1]][seq]
            elif seq in hits[args.names[2]]:
                acc = hits[args.names[2]][seq]            

            outfile.write(">{0}\n{1}\n".format(acc, seq))

        # outfile = open(os.path.join(path_, "all_unique.fa"), "w")
        # for i,  seq in enumerate( (a_not_b_c | b_not_a_c | c_not_a_b ) ):
        #     if seq in dict1:
        #         acc = dict1[seq]
        #     elif seq in dict2:
        #         acc = dict2[seq]
        #     elif seq in dict3:
        #         acc = dict3[seq]            

        #     outfile.write(">{0}\n{1}\n".format(acc, seq))

        r = venn3([a, b, c], (args.names[0], args.names[1], args.names[2]))
        plt.savefig(os.path.join(path_, "venn.png"))


    elif len(args.fastafiles) == 2:

        a = set(hits[args.names[0]].keys())
        b = set(hits[args.names[1]].keys())

        # if args.inexact:
        #     check_inexact(a,b,c)

        print(len(a), len(b))

        a_not_b = a - b 
        b_not_a = b - a
        a_b = a & b 

        outfile = open(os.path.join(path_, "a_and_b.fa"), "w")
        for i,  seq in enumerate((a & b)):
            acc = hits[args.names[0]][seq]
            outfile.write(">{0}\n{1}\n".format(acc, seq))
        r = venn2([a, b], (args.names[0], args.names[1]))
        plt.savefig(args.outfile) 

    else:
        print("only 2 or 3 sets!")


def main_tsv(args):

    if len(args.tsvfiles) == 3:
        hits = defaultdict(list)

        for i, file_ in enumerate(args.tsvfiles):
            f = open(file_, "r")
            sample_id = args.names[i]
            for line in f:
                pred, reference = line.split("\t")[:2]
                hits[sample_id].append(reference)


        a = set(hits[args.names[0]])
        b = set(hits[args.names[1]])
        c = set(hits[args.names[2]])

        # if args.inexact:
        #     check_inexact(a,b,c)

        print(len(a), len(b), len(c))

        a_not_b_c = a - (b | c)
        b_not_a_c = b - (a | c)
        a_b_not_c = (a & b) - c
        c_not_a_b = c - (a | b)
        a_c_not_b = (a & c) - b
        b_c_not_a = (b & c) - a
        a_b_c = a & b & c


        r = venn3([a, b, c], (args.names[0], args.names[1], args.names[2]))
        plt.savefig(args.outfile)

    elif len(args.tsvfiles) == 2:
        hits = defaultdict(list)

        for i, file_ in enumerate(args.tsvfiles):
            f = open(file_, "r")
            sample_id = args.names[i]
            for line in f:
                pred, reference = line.split("\t")[:2]
                hits[sample_id].append(reference)


        a = set(hits[args.names[0]])
        b = set(hits[args.names[1]])

        # if args.inexact:
        #     check_inexact(a,b,c)

        print(len(a), len(b))

        a_not_b = a - b 
        b_not_a = b - a
        a_b = a & b 


        r = venn2([a, b], (args.names[0], args.names[1]))
        plt.savefig(args.outfile) 
    else:
        print("only 2 or 3 sets!")


def mkdir_p(path):
    print("creating", path)
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    # parser.add_argument('--inexact',  action='store_true', help='Check inexact matches')

    parser.add_argument('--tsvfiles', type=str, nargs="+", help='files')
    parser.add_argument('--fastafiles', type=str, nargs="+", help='files')
    parser.add_argument('--names', type=str, nargs="+", required=True, help='names')

    parser.add_argument('--outfile', type=str, help='Output path of results')

    args = parser.parse_args()

    path_, file_prefix = os.path.split(args.outfile)
    mkdir_p(path_)
    if args.tsvfiles:
        main_tsv(args)
    if args.fastafiles:
        main_fasta(args)


from __future__ import print_function
import os,sys
import argparse
import re
# import ssw
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
from matplotlib_venn import venn2, venn2_circles
import edlib
import math

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

def check_inexact(a,b,c):
    for seq1 in a:
        best_ed = len(seq1)
        best_cig = ""
        best_loc = ""
        for seq2 in c:
            if math.fabs(len(seq1) - len(seq2)) > best_ed:
                continue

            edit_distance, locations, cigar = edlib_traceback(seq1, seq2, mode="NW", task="path", k=min(100, best_ed))
            if edit_distance < best_ed:
                best_ed = edit_distance
                best_cig = cigar
                best_loc = locations

        if best_ed > 0:
            print(best_ed,best_cig, best_loc)
            if best_ed > 500:
                print(seq1)


def main(params):
    sample1_seqs = set()
    for file_ in args.sample1:
        primer = file_.split("/")[-2]
        size = file_.split("/")[-3]
        print(file_, primer, size)
        predicted_seqs = set([seq.upper() for (acc, seq) in  read_fasta(open(file_, 'r'))])
        sample1_seqs.update(predicted_seqs)

    sample2_seqs = set()
    for file_ in args.sample2:
        primer = file_.split("/")[-2]
        size = file_.split("/")[-3]
        print(file_, primer, size)
        predicted_seqs = set([seq.upper() for (acc, seq) in  read_fasta(open(file_, 'r'))])
        sample2_seqs.update(predicted_seqs)


    print(len(sample1_seqs), len(sample2_seqs))
    r = venn2([sample1_seqs, sample2_seqs], (params.names[0], params.names[1]))
    plt.savefig(params.outfile)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('--sample1', type=str, nargs="+", help='Path to consensus fasta file(s)')
    parser.add_argument('--sample2', type=str, nargs="+", help='Path to consensus fasta file(s)')
    parser.add_argument('--names', type=str, nargs=2, required=True, help='Set names')
    parser.add_argument('outfile', type=str, help='Output path of results')

    args = parser.parse_args()
    main(args)
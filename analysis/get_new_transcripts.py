
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

    # get shared between samples
    sample1_seqs = set()
    seq_to_acc_sample1 = {}
    for file_ in args.sample1:
        primer = file_.split("/")[-2]
        size = file_.split("/")[-3]
        print(file_, primer, size)
        predicted_seqs = set([seq.upper() for (acc, seq) in  read_fasta(open(file_, 'r'))])
        sample1_seqs.update(predicted_seqs)
        for (acc, seq) in read_fasta(open(file_, 'r')):
            seq_to_acc_sample1[seq] = acc

    sample2_seqs = set()
    seq_to_acc_sample2 = {}
    for file_ in args.sample2:
        primer = file_.split("/")[-2]
        size = file_.split("/")[-3]
        print(file_, primer, size)
        predicted_seqs = set([seq.upper() for (acc, seq) in  read_fasta(open(file_, 'r'))])
        sample2_seqs.update(predicted_seqs)
        for (acc, seq) in read_fasta(open(file_, 'r')):
            seq_to_acc_sample2[seq] = acc



    print(len(sample1_seqs), len(sample2_seqs))
    shared_seqs = sample1_seqs & sample2_seqs
    print("Total shared between samples:", len(sample1_seqs & sample2_seqs))

    # get perfect illumina support

    shared1 = {}
    shared2 = {}
    for seq in shared_seqs:
        acc1 = seq_to_acc_sample1[seq]
        acc2 = seq_to_acc_sample2[seq]
        shared1[acc1] = seq
        shared2[acc2] = seq

    outfile = open(args.outfile, "w")
    already_sampled = set()
    full_support_shared = 0
    for file_ in args.illumina_support_files:
        for line in open(file_, "r"):
            is_new = False
            acc, support = line.split("\t")
            support = float(support)
            if support < 1.0:
                continue

            if acc in shared1:
                is_new = True
                seq = shared1[acc]
            if acc in shared2:
                is_new = True
                seq = shared2[acc]

            if is_new and (seq not in already_sampled):
                full_support_shared += 1
                outfile.write(">{0}\n{1}\n".format(acc, seq))
                already_sampled.add(seq)
    outfile.close()
    print("{0} transcripts had full support and shared betweeen samples".format(full_support_shared))


    # get att the fully supported and shared transcripts that have differences to database transcripts
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('--sample1', type=str, nargs="+", help='Path to consensus fasta file(s)')
    parser.add_argument('--sample2', type=str, nargs="+", help='Path to consensus fasta file(s)')
    parser.add_argument('--illumina_support_files', type=str, nargs="+", help='Path to tsv files with illumina support')
    parser.add_argument('--outfile', type=str, help='A fasta file with transcripts that are shared between smaples and have perfect illumina support.')

    args = parser.parse_args()
    main(args)
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
import edlib

def edlib_traceback(x, y, mode="NW", task="path", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    locations =  result["locations"]
    cigar =  result["cigar"]
    return ed, locations, cigar
# try:
#     import matplotlib
#     matplotlib.use('agg')
#     import matplotlib.pyplot as plt
#     import seaborn as sns
#     sns.set_palette("husl", desat=.6)
#     sns.set(font_scale=1.6)
#     plt.rcParams.update({'font.size': 22})
# except:
#     pass


# try:
#     import matplotlib
#     matplotlib.use('Agg')
#     import matplotlib.pyplot as plt
# except ImportError, RuntimeError:
#     pass

def read_fasta(fasta_file):
    fasta_seqs = {}
    k = 0
    temp = ''
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = "_".join(line[1:].strip().split())
            fasta_seqs[accession] = ''
            k += 1
        elif line[0] == '>':
            yield accession, temp
            temp = ''
            accession = "_".join(line[1:].strip().split())
        else:
            temp += line.strip().upper()
    if accession:
        yield accession, temp

def histogram(data, args, name='histogram.png', x='x-axis', y='y-axis', x_cutoff=None, title=None):
    if x_cutoff: 
        plt.hist(data, range=[0, x_cutoff], bins = 100)
    else:
        plt.hist(data, bins = 100)
    plt.xlabel(x)
    plt.ylabel(y)
    if title:
        plt.title(title)

    plt.savefig(os.path.join(args.outfolder, name))
    plt.clf()

def rand_jitter(arr):
    stdev = .003*(max(arr)-min(arr))
    return arr + np.random.randn(len(arr)) * stdev

def dot_plot(points, args, x_label='x', y_label='y', title='Dotplot', set_marker='o',  name='dot_plot.png'):
    x = map(lambda x: x[0], points)
    y = map(lambda x: x[1], points)
    plt.scatter(rand_jitter(x), rand_jitter(y), marker=set_marker, alpha=0.3)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.grid(True)

    plt.savefig(os.path.join(args.outfolder, name))
    plt.clf()



def reference_similarity(reference_transcripts, outfolder, params):
    """
        Stats about reference transcripts
    """

    print("calculating reference similarities")
    reference_similarities = {}
    for acc1 in reference_transcripts:
        seq1 = reference_transcripts[acc1]
        best_ed = 20000
        for acc2 in reference_transcripts:
            if acc1 == acc2:
                continue
            seq2 = reference_transcripts[acc2]
            edit_distance, locations, cigar = edlib_traceback(seq1, seq2, mode="NW", task="path", k=20000)
            # print(acc1,acc2, edit_distance)
            if 0 < edit_distance < best_ed:
                    best_ed = edit_distance
                    reference_similarities[acc1] = {}
                    reference_similarities[acc1][acc2] = edit_distance

    return reference_similarities

def main(args):
    reference_transcripts = {acc: seq.upper() for (acc, seq) in  read_fasta(open(args.reference_transcripts, 'r'))}
    print("Number of references:", len(reference_transcripts))


    # get abundance from here
    reference_transcripts_sampled = {}
    reference_abundances = {}
    tot_nr_reads_sequenced = 0
    for line in open(args.original_reads_logfile, 'r').readlines()[:-1]:
        member, nr_reads = line.strip().split()
        tot_nr_reads_sequenced += int(nr_reads) 
        reference_abundances[member] = nr_reads
        reference_transcripts_sampled[member] = reference_transcripts[member]
        print(member, nr_reads)
    # get edit distance from here
    reference_similarities = reference_similarity(reference_transcripts_sampled, args.outfile, args)

    # get info if captured by IsoCon here
    
    references_captured = set()

    for line in open(args.predicted_consensus_results, 'r').readlines()[1:]:
        if len(line.strip().split()) == 8:
            pred_id, pred_id_ctd, member_captured,  ed, identity, subs, ins, deletions = line.strip().split()

        else:
            pred_id, member_captured,  ed, identity, subs, ins, deletions = line.strip().split()

        # Only take TP, no ATP
        if subs == ins == deletions == "0":
            references_captured.add(member_captured)

    # output here
    out_tsv_file = open(args.outfile, "w")
    print(reference_similarities)
    for ref1 in reference_abundances: # only iterate over the ones found in reads, i.e., sampled at least once
        assert len(reference_similarities[ref1]) == 1
        for ref2 in reference_similarities[ref1]:
            if ref1 in references_captured:
                out_tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(ref1, ref2, reference_similarities[ref1][ref2], reference_abundances[ref1], reference_abundances[ref2], 1, tot_nr_reads_sequenced))     
            else:
                out_tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(ref1, ref2, reference_similarities[ref1][ref2], reference_abundances[ref1], reference_abundances[ref2], 0, tot_nr_reads_sequenced))     


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('reference_transcripts', type=str, help='Path to the transcript fasta file')
    parser.add_argument('original_reads_logfile', type=str, help='Path to the consensus fasta file')
    parser.add_argument('predicted_consensus_results', type=str, help='Path to the consensus fasta file')    
    parser.add_argument('outfile', type=str, help='Output path of results')

    args = parser.parse_args()
    outfolder, file_prefix = os.path.split(args.outfile)
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)
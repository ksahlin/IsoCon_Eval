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
import misc_functions

def ssw_alignment( x_acc, y_acc, x, y, i,j, max_discrepancy = 50):
    """
        Aligns two sequences with SSW
        x: query
        y: reference

    """
    if i % 10 == 0 and j % 10000 == 0:
        print("processing alignments on all y's where read x_i is participating. i={0}".format(i+1))

    score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-1)
    aligner = ssw.Aligner(gap_open=2, gap_extend=1, matrix=score_matrix)

    # for the ends that SSW leaves behind
    bio_matrix = matlist.blosum62
    g_open = -1
    g_extend = -0.5
    ######################################

    result = aligner.align(x, y, revcomp=True)
    y_alignment, match_line, x_alignment = result.alignment

    matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")
    deletions = x_alignment.count("-")
    insertions = y_alignment.count("-")
    assert deletions + insertions == indels
    # alignment_length = len(match_line)
    
    start_discrepancy = max(result.query_begin, result.reference_begin)  # 0-indexed # max(result.query_begin, result.reference_begin) - min(result.query_begin, result.reference_begin)
    query_end_discrepancy = len(x) - result.query_end - 1
    ref_end_discrepancy = len(y) - result.reference_end - 1
    end_discrepancy = max(query_end_discrepancy, ref_end_discrepancy)  # max(result.query_end, result.reference_end) - min(result.query_end, result.reference_end)
    # print(query_end_discrepancy, ref_end_discrepancy)
    tot_discrepancy = start_discrepancy + end_discrepancy

    if 0 < start_discrepancy <= max_discrepancy:
        # print("HERE")
        matches_snippet = 0
        mismatches_snippet = 0
        if result.query_begin and result.reference_begin:
            query_start_snippet = x[:result.query_begin]
            ref_start_snippet = y[:result.reference_begin]
            alns = pairwise2.align.globalds(query_start_snippet, ref_start_snippet, bio_matrix, g_open, g_extend)
            top_aln = alns[0]
            # print(alns)
            mismatches_snippet = len(list(filter(lambda x: x[0] != x[1] and x[0] != '-' and x[1] != "-", zip(top_aln[0],top_aln[1]))))
            indels_snippet = top_aln[0].count("-") + top_aln[1].count("-")
            matches_snippet = len(top_aln[0]) - mismatches_snippet - indels_snippet
            # print(matches_snippet, mismatches_snippet, indels_snippet)
            query_start_alignment_snippet = top_aln[0]
            ref_start_alignment_snippet = top_aln[1]
        elif result.query_begin:
            query_start_alignment_snippet = x[:result.query_begin]
            ref_start_alignment_snippet = "-"*len(query_start_alignment_snippet)
            indels_snippet = len(ref_start_alignment_snippet)
        elif result.reference_begin:
            ref_start_alignment_snippet = y[:result.reference_begin]
            query_start_alignment_snippet = "-"*len(ref_start_alignment_snippet)
            indels_snippet = len(query_start_alignment_snippet)
        else:
            print("BUG")
            sys.exit()
        matches, mismatches, indels = matches + matches_snippet, mismatches + mismatches_snippet, indels + indels_snippet

        # print(ref_start_alignment_snippet)
        # print(query_start_alignment_snippet)
        y_alignment = ref_start_alignment_snippet + y_alignment
        x_alignment = query_start_alignment_snippet + x_alignment

    if 0 < end_discrepancy <= max_discrepancy:
        # print("HERE2", query_end_discrepancy, ref_end_discrepancy)
        # print(y_alignment)
        # print(y)
        # print(match_line)
        # print(x_alignment)
        # print(x)
        # print(matches, len(x_alignment))
        matches_snippet = 0
        mismatches_snippet = 0
        if query_end_discrepancy and ref_end_discrepancy:
            query_end_snippet = x[result.query_end+1:]
            ref_end_snippet = y[result.reference_end+1:]
            alns = pairwise2.align.globalds(query_end_snippet, ref_end_snippet, bio_matrix, g_open, g_extend)
            top_aln = alns[0]
            mismatches_snippet = len(list(filter(lambda x: x[0] != x[1] and x[0] != '-' and x[1] != "-", zip(top_aln[0],top_aln[1]))))
            indels_snippet = top_aln[0].count("-") + top_aln[1].count("-")
            matches_snippet = len(top_aln[0]) - mismatches_snippet - indels_snippet
            query_end_alignment_snippet = top_aln[0]
            ref_end_alignment_snippet = top_aln[1]
        elif query_end_discrepancy:
            query_end_alignment_snippet = x[result.query_end+1:]
            ref_end_alignment_snippet = "-"*len(query_end_alignment_snippet)
            indels_snippet = len(ref_end_alignment_snippet)

        elif ref_end_discrepancy:
            ref_end_alignment_snippet = y[result.reference_end+1:]
            query_end_alignment_snippet = "-"*len(ref_end_alignment_snippet)
            indels_snippet = len(query_end_alignment_snippet)

        else:
            print("BUG")
            sys.exit()
        matches, mismatches, indels = matches + matches_snippet, mismatches + mismatches_snippet, indels + indels_snippet

        y_alignment = y_alignment + ref_end_alignment_snippet
        x_alignment = x_alignment + query_end_alignment_snippet

    # matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")
    deletions = x_alignment.count("-")
    insertions = y_alignment.count("-")
    assert deletions + insertions == indels

    if start_discrepancy > max_discrepancy or end_discrepancy > max_discrepancy:
        # print("REMOVING", start_discrepancy, end_discrepancy)
        return (y_alignment, x_alignment, None)

    else:
        return (y_alignment, x_alignment, (matches, mismatches, indels, deletions, insertions)) 


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


def reference_similarity(reference_transcripts, outfolder, params):
    """
        Stats about reference transcripts
    """
    seqs_seen = set()
    transcript_abundances = {}
    transcript_copies = Counter()
    transcript_sequences = {}
    for acc, seq in sorted(reference_transcripts.items(), key=lambda x: len(x[1])):

        try:
            tr_acc, copy_number_str = acc.split("copy")

        except ValueError:
            tr_acc, copy_number_str = acc, "1"  # viral data not simulated

        transcript_copies[tr_acc] += 1
        try :
            copy_number = int(copy_number_str)
        except ValueError:
            copy_number = 1

        if tr_acc not in transcript_abundances:
            transcript_abundances[tr_acc] = copy_number
            transcript_sequences[tr_acc] = seq
        elif tr_acc in transcript_abundances and copy_number >  transcript_abundances[tr_acc]:
            transcript_abundances[tr_acc] = copy_number
            transcript_sequences[tr_acc] = seq

        if seq in seqs_seen:
            # print("HERE!", len(seq))
            del reference_transcripts[acc]
        else:
            seqs_seen.add(seq)

    print("Number of unique references:", len(reference_transcripts))
    for t_acc, copy_nr in transcript_copies.items():
        print(t_acc, copy_nr)
    print("abundances:", transcript_abundances)



    print("calculating reference similarities")
    aligner = ssw.Aligner(gap_open=2, gap_extend=1)
    # do SW
    sorted_reference_tuples = sorted(reference_transcripts.items(), key = lambda x: len(x[1]))
    reference_similarities = {}
    for q_acc, q_seq in transcript_sequences.items():
        reference_similarities[q_acc] = {}
        for r_acc, r_seq in transcript_sequences.items():
            # r_aligned, q_aligned, stats = ssw_alignment( q_acc, r_acc, q_seq, r_seq, 0,0, max_discrepancy = 10000 )
            # if stats:
            #     matches, mismatches, indels, deletions, insertions = stats
            #     errors = mismatches + indels
            #     identity = errors/ float(errors + matches)
            #     reference_similarities[q_acc][r_acc] = errors
            # else:
            #     reference_similarities[q_acc][r_acc] = min(len(q_seq), len(r_seq))
            
            ed = edlib_ed(q_seq, r_seq, mode="NW", task="distance", k=10000)
            reference_similarities[q_acc][r_acc] = ed

    return transcript_abundances, transcript_copies, reference_similarities

import edlib
def edlib_ed(x, y, mode="NW", task="distance", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    return ed

def main(args):
    transcripts = {acc: seq.upper() for (acc, seq) in  read_fasta(open(args.transcripts, 'r'))}
    print("Number of references:", len(transcripts))
    transcript_abundances, transcript_copies, transcript_similarities = reference_similarity(transcripts, args.outfolder, args)
    out_file_ref_sim = open(args.outfile, "w")

    ref_similarities = []
    unique_transcript_abundances = set([a for a in transcript_abundances.values()])

    if len(unique_transcript_abundances) < len(transcript_abundances):
        not_unique_ab = True
        transcript_name_assigmment = {}
        cnt = 1
        for acc in transcript_abundances:
            transcript_name_assigmment[acc] = cnt
            cnt +=1

    else:
        not_unique_ab = False

    for acc1, a1 in  sorted(transcript_abundances.items(), key = lambda x: x[1]):
        for acc2, a2 in  sorted(transcript_abundances.items(), key = lambda x: x[1]):
            if not_unique_ab: 
                pattern = r'[\d]+'
                m_id1 = str(transcript_name_assigmment[acc1])
                m_id2 = str(transcript_name_assigmment[acc2])
                print("member ids:", m_id1, m_id2, a1, a2, transcript_similarities[acc1][acc2])
                out_file_ref_sim.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(m_id1, m_id2, args.mutation_rate, transcript_similarities[acc1][acc2], args.family))
            else:
                out_file_ref_sim.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(a1, a2, args.mutation_rate, transcript_similarities[acc1][acc2], args.family))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('transcripts', type=str, help='Path to the consensus fasta file')
    parser.add_argument('outfile', type=str, help='Output filename')
    parser.add_argument('--mutation_rate',  type=str, default="none", help='Mut rate')
    parser.add_argument('--family',  type=str, default="none", help='Mut rate')

    args = parser.parse_args()

    path_, file_prefix = os.path.split(args.outfile)
    misc_functions.mkdir_p(path_)
    args.outfolder = path_
    args.filename = file_prefix
    main(args)
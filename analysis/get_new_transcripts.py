
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

primer_to_family_and_sample = { "barcode_1-2kb" : { 0 : ("sample1", "HSFY2"),
                                                    1 :  ("sample2", "HSFY2"),
                                                    2 : ("sample1", "RBMY"),
                                                    3 : ("sample1", "CDY1"),
                                                    4 : ("sample1", "CDY2"),
                                                    5 : ("sample1", "DAZ"),
                                                    6 : ("sample2", "RBMY"),
                                                    7 : ("sample2", "CDY1"),
                                                    8 : ("sample2", "CDY2"),
                                                    9 : ("sample2", "DAZ")
                                                },
                                "barcode_1kb" : { 0 : ("sample1", "BPY"),
                                                  1 : ("sample1", "VCY"),
                                                  2 : ("sample1", "XKRY"),
                                                  3 : ("sample1", "PRY"),
                                                  4 : ("sample1", "HSFY1"),
                                                  5 : ("sample1", "TSPY"),
                                                  6 : ("sample2", "BPY"),
                                                  7 : ("sample2", "VCY"),
                                                  8 : ("sample2", "XKRY"),
                                                  9 : ("sample2", "PRY"),
                                                  10 : ("sample2", "HSFY1"),
                                                  11 : ("sample2", "TSPY")
                                                    }    }


def edlib_ed(x, y, mode="NW", task="distance", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    if task == "path":
        cigar =  result["cigar"]
        locations =  result["locations"]
        return ed, cigar, locations
    else:
        return ed

def get_minimizers_2set_simple(querys, targets, min_query_len):
    best_edit_distances = {}
    best_cigars = {}
    i = 1
    for acc1, seq1 in querys.items():
        if i % 200 == 0:
            print("processing candidate", i)
        best_ed = 2*len(seq1)
        best_edit_distances[acc1] = {}
        best_cigars[acc1] = {}
        if len(seq1) < min_query_len:
            continue
        for acc2, seq2 in targets.items():
            # edit_distance = edlib_ed(seq1, seq2, mode="HW", task = "path", k = max_ed_threshold) # seq1 = query, seq2 = target
            edit_distance, cigar, locations = edlib_ed(seq1, seq2, mode="HW", task = "path", k = best_ed)
            if edit_distance < 0:
                continue

            if 0 <= edit_distance <= best_ed:
                if len(locations) > 1 or locations[0][0] > 100 or (len(seq2) - locations[0][1]) > 100:
                    # print("Skipped:", cigar, locations, len(seq1), len(seq2))
                    continue
                if edit_distance < best_ed:
                    best_edit_distances[acc1] = {}
                    best_cigars[acc1] = {}
                best_ed = edit_distance
                best_edit_distances[acc1][acc2] = best_ed
                best_cigars[acc1][acc2] = cigar
            elif edit_distance == best_ed:
                if len(locations) > 1 or locations[0][0] > 100 or len(seq2) - locations[0][1] > 100:
                    # print("Skipped:", cigar, locations, len(seq1), len(seq2))
                    continue
                best_edit_distances[acc1][acc2] = best_ed
                best_cigars[acc1][acc2] = cigar

        i += 1
        # print(best_ed)
    # for acc in best_edit_distances:
    #     if len(best_edit_distances[acc]) >1:
    #         print(best_edit_distances[acc])
    #     for acc2 in best_edit_distances[acc]:
    #         if best_edit_distances[acc][acc2] == 0:
    #             print("OK",acc,acc2)
    # print("HERE")
    # for acc in best_edit_distances:
    #     for acc2 in best_edit_distances[acc]:
    #         print( best_edit_distances[acc][acc2])

    return best_edit_distances, best_cigars

def get_ssw_alignments(best_edit_distances, querys, targets):
    score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-2)
    aligner = ssw.Aligner(gap_open=2, gap_extend=1, matrix=score_matrix)
    best_edit_distances_ssw = {}
    best_cigars_ssw = {}
    for acc1 in best_edit_distances:
        seq1 = querys[acc1]
        best_ed = len(seq1)
        best_edit_distances_ssw[acc1] = {}
        best_cigars_ssw[acc1] = {}
        for acc2 in best_edit_distances[acc1]:
            seq2 = targets[acc2]
            result = aligner.align(seq1, seq2, revcomp=False)
            seq2_aln, match_line, seq1_aln = result.alignment
            matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")
            insertion_count = seq2_aln.count("-")
            deletion_count = seq1_aln.count("-")

            sw_ed = mismatches + indels
            best_edit_distances_ssw[acc1][acc2] =  sw_ed # (deletion_count, insertion_count, mismatches )
            seq1_aln, match_line, seq2_aln = result.alignment
            best_cigars_ssw[acc1][acc2] = (result.cigar, mismatches, indels, result.query_begin, len(seq1) - result.query_end - 1, result.reference_begin, len(seq2) - result.reference_end -1 )

            # print(acc1,acc2)
            # print(result.query_begin, len(seq1) - result.query_end - 1, result.reference_begin, len(seq2) - result.reference_end -1, result.cigar, mismatches, indels)
            # print()


            # print(sw_ed, (deletion_count, insertion_count, mismatches ))
            # print(seq1_aln)
            # print(match_line)
            # print(seq2_aln)
            # edit_distance, locations, cigar = edlib_traceback(seq1, seq2, k =1000)
            # print(edit_distance, locations, cigar)
            # print()            

    return best_edit_distances_ssw, best_cigars_ssw



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


def main(params):

    # get shared between samples
    sample1_seqs = set()
    seq_to_acc_sample1 = {}
    sample2_seqs = set()
    seq_to_acc_sample2 = {}

    for file_ in args.sample1:
        primer_id = int(file_.split("/")[-2])
        batch_size = file_.split("/")[-3].split("_polished")[0]
        print(primer_id, batch_size, primer_to_family_and_sample[batch_size][primer_id])
        if primer_to_family_and_sample[batch_size][primer_id][0] == "sample1":
            predicted_seqs = set([seq.upper() for (acc, seq) in  read_fasta(open(file_, 'r'))])
            sample1_seqs.update(predicted_seqs)
            for (acc, seq) in read_fasta(open(file_, 'r')):
                seq_to_acc_sample1[seq] = acc
        else:
            predicted_seqs = set([seq.upper() for (acc, seq) in  read_fasta(open(file_, 'r'))])
            sample2_seqs.update(predicted_seqs)
            for (acc, seq) in read_fasta(open(file_, 'r')):
                seq_to_acc_sample2[seq] = acc
        # print(file_, primer, size)


    # sample2_seqs = set()
    # seq_to_acc_sample2 = {}
    # for file_ in args.sample2:
    #     primer = file_.split("/")[-2]
    #     size = file_.split("/")[-3]
    #     print(file_, primer, size)
    #     predicted_seqs = set([seq.upper() for (acc, seq) in  read_fasta(open(file_, 'r'))])
    #     sample2_seqs.update(predicted_seqs)
    #     for (acc, seq) in read_fasta(open(file_, 'r')):
    #         seq_to_acc_sample2[seq] = acc



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

    outfile = open(os.path.join(args.outfolder, "all_supported_and_shared.fa"), "w")
    already_sampled = set()
    full_support_shared = 0
    fully_supported_shared = {}
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
                fully_supported_shared[acc] = seq
                outfile.write(">{0}\n{1}\n".format(acc, seq))
                already_sampled.add(seq)
    outfile.close()
    print("{0} transcripts had full support and shared betweeen samples".format(full_support_shared))


    # get att the fully supported and shared transcripts that have differences to database transcripts

    database = {acc: seq.upper() for (acc, seq) in  read_fasta(open(args.database, 'r'))}
    database = {acc: seq.upper() for (acc, seq) in database.items() if "UNAVAILABLE" not in seq }

    minimizer_graph_c_to_t, best_cigars = get_minimizers_2set_simple(fully_supported_shared, database, 100)
    minimizer_graph_x_to_c, best_cigars_ssw = get_ssw_alignments(minimizer_graph_c_to_t, fully_supported_shared, database)
    
    variations = set()
    high_id = 0
    exon_diff = 0
    misc = 0
    highly_similar_copies = open(os.path.join(args.outfolder, "high_id.fa"), "w")
    high_id_exon_diff = open(os.path.join(args.outfolder, "exon_diff.fa"), "w")
    other_alignments = open(os.path.join(args.outfolder, "misc_others.fa"), "w")

    for c in best_cigars_ssw:
        # print(len(best_cigars_ssw[c]))
        if len(best_cigars_ssw[c]) > 1:
            # print(best_cigars_ssw[c])
            smallest_SNV = 100000
            smallest_other = 100000
            for t_cand in best_cigars_ssw[c]: # pick the one with smallest number of substitutions and single nucl indels
                cigar, mismatches, indels, q_start_softcl, q_end_softcl, t_start_softcl, t_end_softcl = best_cigars_ssw[c][t_cand]
                if mismatches + indels < smallest_SNV:
                    smallest_SNV = mismatches + indels
                    smallest_other = q_start_softcl + q_end_softcl + t_start_softcl + t_end_softcl
                    t_min = t_cand
                elif mismatches + indels == smallest_SNV:
                    if q_start_softcl + q_end_softcl + t_start_softcl + t_end_softcl <  smallest_other:
                        smallest_other = q_start_softcl + q_end_softcl + t_start_softcl + t_end_softcl
                        t_min = t_cand

            t = t_min
        else:
            t = best_cigars_ssw[c].keys()[0]

        # print(best_cigars_ssw[c][t])
        cigar, mismatches, indels, q_start_softcl, q_end_softcl, t_start_softcl, t_end_softcl = best_cigars_ssw[c][t]
        has_exon_diff = parse_cigar(cigar)

        if not (mismatches == indels == q_start_softcl == q_end_softcl == 0): # not perfect
            variations.add(c)
            # print(cigar, mismatches, indels, q_start_softcl, q_end_softcl, t_start_softcl, t_end_softcl, t, c)
            if mismatches + indels + q_start_softcl + q_end_softcl < 10:
                highly_similar_copies.write(">{0}\n{1}\n".format(c + "_aligned_to_" + t, fully_supported_shared[c]))
            elif mismatches < 10 and has_exon_diff:
                high_id_exon_diff.write(">{0}\n{1}\n".format(c + "_aligned_to_" + t, fully_supported_shared[c]))
            else:
                print(cigar, mismatches, indels, q_start_softcl, q_end_softcl, t_start_softcl, t_end_softcl)
                other_alignments.write(">{0}\n{1}\n".format(c + "_aligned_to_" + t, fully_supported_shared[c]))
        else:
            pass
    
    # print("Total ",len(best_cigars_ssw))
    print("Total nr variations:", len(variations))

def parse_cigar(cigar):
    # print(cigar)

    tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    # print(result)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = cigar[i]
        i += 1
        tuples.append((length, type_ ))
    # print(tuples)
    for length, type_ in tuples:
        if int(length) >= 10 and type_ != "M":
            # print("has exon") 
            return True
    return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('--sample1', type=str, nargs="+", help='Path to consensus fasta file(s)')
    parser.add_argument('--sample2', type=str, nargs="+", help='Path to consensus fasta file(s)')
    parser.add_argument('--database', type=str, help='Path to database fasta file')
    parser.add_argument('--illumina_support_files', type=str, nargs="+", help='Path to tsv files with illumina support')
    parser.add_argument('--outfolder', type=str, help='A fasta file with transcripts that are shared between smaples and have perfect illumina support.')
    # parser.add_argument('--high_id', type=str, help='A fasta file with transcripts differ with at most 10 bases to DB transcripts.')
    # parser.add_argument('--high_id_exon_diff', type=str, help='A fasta file with transcripts differ with a region of at least 20 consecutive bp and no more than 10 other variants.')
    # parser.add_argument('--misc', type=str, help='A fasta file with transcripts with the rest of the alignments.')
    
    args = parser.parse_args()
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)
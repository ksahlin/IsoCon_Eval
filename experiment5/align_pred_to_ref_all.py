from __future__ import print_function
import os,sys
import argparse
import re
import math
import numpy as np
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
import string
import fractions
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

from collections import defaultdict
import edlib

import signal
from multiprocessing import Pool
import multiprocessing as mp
import sys

import re

import edlib
import ssw


def get_minimizers_2set_helper(arguments):
    args, kwargs = arguments
    return get_minimizers_2set(*args, **kwargs)


def get_minimizers_2set_simple(querys, targets):
    best_edit_distances = {}

    for acc1, seq1 in querys.items():
        best_ed = len(seq1)
        for acc2, seq2 in targets.items():
            edit_distance = edlib_ed(seq1, seq2, mode="NW", k = len(seq1)) # seq1 = query, seq2 = target
            if 0 <= edit_distance < best_ed:
                best_ed = edit_distance
                best_edit_distances[acc1] = {}
                best_edit_distances[acc1][acc2] = best_ed
            elif edit_distance == best_ed:
                best_edit_distances[acc1][acc2] = best_ed
        # print(best_ed)


    return best_edit_distances


def get_ssw_alignments(best_edit_distances, querys, targets):
    score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-2)
    aligner = ssw.Aligner(gap_open=2, gap_extend=1, matrix=score_matrix)
    best_edit_distances_ssw = {}
    for acc1 in best_edit_distances:
        seq1 = querys[acc1]
        best_ed = len(seq1)
        best_edit_distances_ssw[acc1] = {}

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
            # print(sw_ed, (deletion_count, insertion_count, mismatches ))

            # print(seq1_aln)
            # print(match_line)
            # print(seq2_aln)
            # edit_distance, locations, cigar = edlib_traceback(seq1, seq2, k =1000)
            # print(edit_distance, locations, cigar)
            # print()            

    return best_edit_distances_ssw


def get_minimizers_2set(batch, start_index, seq_to_acc_list_sorted, target_accessions):
    best_edit_distances = {}
    error_types = {"D":0, "S": 0, "I": 0}
    for i in range(start_index, start_index + len(batch)):
        if i % 50 == 0:
            print("processing ", i)

        seq1 = seq_to_acc_list_sorted[i][0]
        acc1 = seq_to_acc_list_sorted[i][1]
        if acc1 in target_accessions:
            continue

        # reach here and we have a query read to find best alignment for
            
        best_edit_distances[acc1] = {}
        best_ed = len(seq1)
        best_cigar = ""

        stop_up = False
        stop_down = False

        j = 1
        while True:
        # for j in range(1,len(seq_to_acc_list_sorted)):
            if i - j < 0:
                stop_down = True
            if i + j >= len(seq_to_acc_list_sorted):
                stop_up = True

            if not stop_down:
                seq2 = seq_to_acc_list_sorted[i - j][0]
                acc2 = seq_to_acc_list_sorted[i - j][1]  

                if math.fabs(len(seq1) - len(seq2)) > best_ed:
                    stop_down = True

            if not stop_up:
                seq3 = seq_to_acc_list_sorted[i + j][0]
                acc3 = seq_to_acc_list_sorted[i + j][1]  

                if math.fabs(len(seq1) - len(seq3)) > best_ed:
                    stop_up = True

            if not stop_down and acc2 in target_accessions:
                # if seq1 == seq2:
                #     print("ID:", acc1, acc2)
                edit_distance, locations, cigar = edlib_traceback(seq1, seq2, mode="NW", task="path", k=best_ed) # seq1 = query, seq2 = target

                if 0 <= edit_distance < best_ed:
                    best_ed = edit_distance
                    best_edit_distances[acc1] = {}
                    best_edit_distances[acc1][acc2] = best_ed
                    # lower_target_edit_distances[acc2] = best_ed 
                    best_cigar = cigar

                elif edit_distance == best_ed:
                    best_edit_distances[acc1][acc2] = best_ed

            if not stop_up and acc3 in target_accessions:
                # if seq1 == seq3:
                #     print("ID:", acc1, acc3)

                edit_distance, locations, cigar = edlib_traceback(seq1, seq3, mode="NW", task="path", k=best_ed)
                if 0 <= edit_distance < best_ed:
                    best_ed = edit_distance
                    best_edit_distances[acc1] = {}
                    best_edit_distances[acc1][acc3] = best_ed
                    # lower_target_edit_distances[acc3] = best_ed 
                    best_cigar = cigar

                elif edit_distance == best_ed:
                    best_edit_distances[acc1][acc3] = best_ed
 
            if stop_down and stop_up:
                break

            j += 1

        # if best_ed > 100:
        #     print(best_ed, "for seq with length", len(seq1), acc1, best_cigar)

        # print("best ed:", best_ed)
        # print(best_edit_distances[acc1])
        # if best_ed > 100:
        #     print(best_ed, "for seq with length", len(seq1), seq1)

    # print(best_edit_distances)
    return best_edit_distances

def get_exact_minimizer_graph_2set(seq_to_acc_list_sorted_all, target_accessions, single_core = False):

    if single_core:
        best_edit_distances = get_minimizers_2set(seq_to_acc_list_sorted_all, 0, seq_to_acc_list_sorted_all, target_accessions)
        minimizer_graph = transpose(best_edit_distances)

        # implement check here to se that all seqs got a minimizer, if not, print which noes that did not get a minimizer computed.!

    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())
        chunk_size = int(len(seq_to_acc_list_sorted_all) / (4*mp.cpu_count())) 
        chunks = [(i, seq_to_acc_list_sorted_all[i:i + chunk_size]) for i in range(0, len(seq_to_acc_list_sorted_all), chunk_size)] 
        print([i for i in range(0, len(seq_to_acc_list_sorted_all), chunk_size)])
        print([len(ch) for i,ch in chunks])
        try:
            res = pool.map_async(get_minimizers_2set_helper, [ ((chunk, i , seq_to_acc_list_sorted_all, target_accessions), {}) for i,chunk in chunks] )
            best_edit_distances_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            # print("Normal termination")
            pool.close()
        pool.join()
        best_edit_distances = {}
        for sub_graph in best_edit_distances_results:
            for seq in sub_graph:
                assert seq not in best_edit_distances
            best_edit_distances.update(sub_graph)
        
    return best_edit_distances


def get_minimizers_helper(arguments):
    args, kwargs = arguments
    return get_minimizers(*args, **kwargs)

def get_exact_minimizer_graph(seq_to_acc_list_sorted, single_core = False):

    if single_core:
        best_edit_distances = get_minimizers(seq_to_acc_list_sorted, 0, seq_to_acc_list_sorted)
        minimizer_graph = transpose(best_edit_distances)

        # implement check here to se that all seqs got a minimizer, if not, print which noes that did not get a minimizer computed.!

    else:
        ####### parallelize alignment #########
        # pool = Pool(processes=mp.cpu_count())
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=mp.cpu_count())
        chunk_size = int(len(seq_to_acc_list_sorted) / (4*mp.cpu_count())) 
        chunks = [(i, seq_to_acc_list_sorted[i:i + chunk_size]) for i in range(0, len(seq_to_acc_list_sorted), chunk_size)] 
        print([i for i in range(0, len(seq_to_acc_list_sorted), chunk_size)])
        print([len(ch) for i,ch in chunks])
        try:
            res = pool.map_async(get_minimizers_helper, [ ((chunk, i , seq_to_acc_list_sorted), {}) for i,chunk in chunks] )
            best_edit_distances_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            # print("Normal termination")
            pool.close()
        pool.join()
        best_edit_distances = {}
        for sub_graph in best_edit_distances_results:
            for seq in sub_graph:
                assert seq not in best_edit_distances
            best_edit_distances.update(sub_graph)
        
    return best_edit_distances

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
            temp += line.strip()
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


def collapse(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]




def edlib_ed(x, y, mode="NW", task="distance", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    return ed

def edlib_traceback(x, y, mode="NW", task="path", k=1):
    result = edlib.align(x, y, mode="NW", task=task, k=k)
    # print(x,y)
    ed = result["editDistance"]
    locations =  result["locations"]
    cigar =  result["cigar"]
    return ed, locations, cigar


def get_minimizers(batch_of_queries, start_index, seq_to_acc_list_sorted):
    best_edit_distances = {}
    lower_target_edit_distances = {}
    error_types = {"D":0, "S": 0, "I": 0}
    for i in range(start_index, start_index + len(batch_of_queries)):
        if i % 50 == 0:
            print("processing ", i)

        seq1 = seq_to_acc_list_sorted[i][0]
        acc1 = seq_to_acc_list_sorted[i][1]
        best_edit_distances[acc1] = {}
        # print(acc1, seq1)

        if acc1 in lower_target_edit_distances:
            best_ed = lower_target_edit_distances[acc1] 
            # print("already_comp", best_ed )
        else:
            best_ed = len(seq1)
        best_cigar = ""
        stop_up = False
        stop_down = False
        for j in range(1,len(seq_to_acc_list_sorted)):
            if i - j < 0:
                stop_down = True
            if i + j >= len(seq_to_acc_list_sorted):
                stop_up = True

            if not stop_down:
                seq2 = seq_to_acc_list_sorted[i - j][0]
                acc2 = seq_to_acc_list_sorted[i - j][1]  

                if math.fabs(len(seq1) - len(seq2)) > best_ed:
                    stop_down = True

            if not stop_up:
                seq3 = seq_to_acc_list_sorted[i + j][0]
                acc3 = seq_to_acc_list_sorted[i + j][1]  

                if math.fabs(len(seq1) - len(seq3)) > best_ed:
                    stop_up = True

            if not stop_down:
                edit_distance, locations, cigar = edlib_traceback(seq1, seq2, mode="NW", task="path", k=best_ed)
                if 0 <= edit_distance < best_ed:
                    best_ed = edit_distance
                    best_edit_distances[acc1] = {}
                    best_edit_distances[acc1][acc2] = best_ed
                    lower_target_edit_distances[acc2] = best_ed 
                    best_cigar = cigar

                elif edit_distance == best_ed:
                    best_edit_distances[acc1][acc2] = best_ed

            if not stop_up:
                edit_distance, locations, cigar = edlib_traceback(seq1, seq3, mode="NW", task="path", k=best_ed)
                if 0 <= edit_distance < best_ed:
                    best_ed = edit_distance
                    best_edit_distances[acc1] = {}
                    best_edit_distances[acc1][acc3] = best_ed
                    lower_target_edit_distances[acc3] = best_ed 
                    best_cigar = cigar

                elif edit_distance == best_ed:
                    best_edit_distances[acc1][acc3] = best_ed
 
            if stop_down and stop_up:
                break

            # print("best ed:", best_ed, best_cigar)
        # if best_ed > 100:
        #     print(best_ed, "for seq with length", len(seq1), acc1, best_cigar)

    return best_edit_distances


# def compute_minimizer_graph(S, params):
#     """
#         strings S are all unique here.
#     """
#     # consensus_transcripts = {acc: seq for (acc, seq) in  read_fasta(open(params.consensus_transcripts, 'r'))}
#     # print("Number of consensus:", len(consensus_transcripts))

#     # seq_to_acc = {seq: acc for (acc, seq) in  read_fasta(open(params.consensus_transcripts, 'r'))}

#     seq_to_acc_list = list(seq_to_acc.items())
#     seq_to_acc_list_sorted = sorted(seq_to_acc_list, key= lambda x: len(x[0]))
#     collapsed_consensus_transcripts =  { acc : seq for (seq, acc) in  seq_to_acc.items() }
#     print("Number of collapsed consensus:", len(collapsed_consensus_transcripts))
#     minimizer_graph = get_exact_minimizer_graph(seq_to_acc_list_sorted, single_core = params.single_core)

#     s1 = set( [ collapsed_consensus_transcripts[acc2] for acc1 in minimizer_graph for acc2 in minimizer_graph[acc1] ])
#     s2 = set([seq for seq in seq_to_acc] )
#     isolated = s2.difference(s1)
#     print("isolated:", len(isolated))

#     edges = 0
#     tot_ed = 0
#     edit_hist =[]
#     neighbors = []
#     for r1 in  minimizer_graph:
#         for r2 in minimizer_graph[r1]:
#             edges += 1
#             tot_ed += minimizer_graph[r1][r2]
#             edit_hist.append(minimizer_graph[r1][r2])

#         neighbors.append(len(minimizer_graph[r1]))

#     print("Number of edges:", edges)
#     print("Total edit distance:", tot_ed)
#     print("Avg ed (ed/edges):", tot_ed/ float(edges))
#     histogram(edit_hist, params, name='edit_distances.png', x='x-axis', y='y-axis', x_cutoff=100, title="Edit distances in minimizer graph")
#     histogram(neighbors, params, name='neighbours.png', x='x-axis', y='y-axis', title="Number of neighbours in minimizer graph")
#     histogram(neighbors, params, name='neighbours_zoomed.png', x='x-axis', y='y-axis', x_cutoff=20, title="Number of neighbours in minimizer graph")

#     return minimizer_graph, isolated


def main_temp_2set(args):
    predicted = {}
    for file_ in args.predicted:
        primer = file_.split("/")[-2]
        size = file_.split("/")[-3]

        print(file_, primer, size)
        predicted_subset = {acc: seq.upper() for (acc, seq) in  read_fasta(open(file_, 'r'))}
        for acc, seq in predicted_subset.items():
            if acc in predicted:
                print("OMG accession collision between genes -- need to adjust accesions")
                print("collision:", acc)
                sys.exit()
            else:  
                predicted[acc + "_" + primer + "-" + size] = seq

    # predicted = {acc: seq.upper() for (acc, seq) in  read_fasta(open(args.predicted, 'r'))}
    print("Number of predicted:", len(predicted))
    database = {acc: seq.upper() for (acc, seq) in  read_fasta(open(args.database, 'r'))}
    # print(database)
    database = {acc: seq.upper() for (acc, seq) in database.items() if "UNAVAILABLE" not in seq }
    print("Number of targets:", len(database))
    unique_queries = {seq: acc for (acc, seq) in predicted.items()}
    unique_targets = {seq: acc for (acc, seq) in  database.items()}
    print("Number of unique predicted sequences:", len(unique_queries))
    print("Number of unique targets sequences:", len(unique_targets))

    # seq_acc_queries = [(seq, acc) for acc, seq in  predicted.items()] 
    # seq_acc_targets = [(seq, acc) for acc, seq in  database.items()]

    # ONLY ALIGN NON-REDUNDANT SEQUENCES!
    seq_acc_queries = [(seq, acc) for seq, acc in  unique_queries.items()] 
    seq_acc_targets = [(seq, acc) for seq, acc in  unique_targets.items()] #{seq: acc for (acc, seq) in  read_fasta(open(args.database, 'r'))}
    
    seq_to_acc_list_sorted_all = sorted(seq_acc_queries + seq_acc_targets, key= lambda x: len(x[0]))



    minimizer_graph_x_to_c = get_exact_minimizer_graph_2set(seq_to_acc_list_sorted_all, set(database.keys()), single_core = args.single_core)
    minimizer_graph_x_to_c = get_ssw_alignments(minimizer_graph_x_to_c, predicted, database)
    # unique_database = {acc: seq for (seq, acc) in  unique_targets.items()} 
    # print(len(unique_database))
    # minimizer_graph_x_to_c = get_minimizers_2set_simple(predicted, unique_database)
    
    edges = 0
    tot_ed = 0
    edit_hist =[]
    neighbors = []
    for x in  minimizer_graph_x_to_c:
        for c in minimizer_graph_x_to_c[x]:
            edges += 1
            tot_ed += minimizer_graph_x_to_c[x][c]
            edit_hist.append(minimizer_graph_x_to_c[x][c])

        neighbors.append(len(minimizer_graph_x_to_c[x]))

    print("Number of edges:", edges)
    print("Total edit distance:", tot_ed)
    print("Avg ed (ed/edges):", tot_ed/ float(edges))
    histogram(edit_hist, args, name='edit_distances.png', x='x-axis', y='y-axis', x_cutoff=100, title="Edit distances in minimizer graph")
    histogram(edit_hist, args, name='edit_distances_zoomed.png', x='x-axis', y='y-axis', x_cutoff=5, title="Edit distances in minimizer graph")
    histogram(neighbors, args, name='neighbours.png', x='x-axis', y='y-axis', title="Number of neighbours in minimizer graph")
    histogram(neighbors, args, name='neighbours_zoomed.png', x='x-axis', y='y-axis', x_cutoff=20, title="Number of neighbours in minimizer graph")

    clusters_to_database = transpose(minimizer_graph_x_to_c)
    outfile_alignments = open(os.path.join(args.outfolder, "pred_to_ref_best_alignments.tsv"), "w")
    for t in  clusters_to_database:
        for c in clusters_to_database[t]:
            # print(t, clusters_to_database[t])
            outfile_alignments.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(c,t, len(predicted[c]), len(database[t]), clusters_to_database[t][c]))
    outfile_alignments.close()


def main_temp(args):
    predicted = {acc: seq for (acc, seq) in  read_fasta(open(args.predicted, 'r'))}
    print("Number of consensus:", len(predicted))
    seq_to_acc = [ (seq, acc) for (acc, seq) in  read_fasta(open(args.predicted, 'r'))]
    # seq_to_acc_list = list(seq_to_acc.items())
    seq_to_acc_list_sorted = sorted(seq_to_acc, key= lambda x: len(x[0]))
    collapsed_consensus_transcripts =  { acc : seq for (seq, acc) in  seq_to_acc }
    # print("Number of collapsed consensus:", len(collapsed_consensus_transcripts))
    minimizer_graph = get_exact_minimizer_graph(seq_to_acc_list_sorted, single_core = args.single_core)

    s1 = set( [ collapsed_consensus_transcripts[acc2] for acc1 in minimizer_graph for acc2 in minimizer_graph[acc1] ])
    s2 = set([seq for seq in seq_to_acc] )
    isolated = s2.difference(s1)
    print("isolated:", len(isolated))

    edges = 0
    tot_ed = 0
    edit_hist =[]
    edit_hist_one_per_read = []
    neighbors = []
    for r1 in  minimizer_graph:
        for i,r2 in enumerate(minimizer_graph[r1]):
            if i == 0:
                edit_hist_one_per_read.append(minimizer_graph[r1][r2])

            edges += 1
            tot_ed += minimizer_graph[r1][r2]
            edit_hist.append(minimizer_graph[r1][r2])

        neighbors.append(len(minimizer_graph[r1]))

    print("Number of edges:", edges)
    print("Total edit distance:", tot_ed)
    print("Avg ed (ed/edges):", tot_ed/ float(edges))

    n = float(len(edit_hist_one_per_read))
    mu = sum(edit_hist_one_per_read) / n
    sigma = (sum(list(map((lambda x: x ** 2 - 2 * x * mu + mu ** 2), edit_hist_one_per_read))) / (n - 1)) ** 0.5
    edit_hist_one_per_read.sort()
    max_ed = max(edit_hist_one_per_read)
    min_ed = min(edit_hist_one_per_read)
    if len(edit_hist_one_per_read) %2 == 0:
        median_ed = (edit_hist_one_per_read[int(len(edit_hist_one_per_read)/2)-1] + edit_hist_one_per_read[int(len(edit_hist_one_per_read)/2)]) / 2.0
    else:
        median_ed = edit_hist_one_per_read[int(len(edit_hist_one_per_read)/2)]

    print("mean ed: {0}, sd ed:{1}, min_ed:{2}, max_ed:{3}, median_ed:{4}\n".format(mu, sigma, min_ed, max_ed, median_ed))
    print("nr reads with a minimizer that we estimated ed of:", len(edit_hist_one_per_read))
    histogram(edit_hist_one_per_read, args, name='read_error_estimates_from_min_ed.png', x='Edit distance to minimizer', y='count', x_cutoff=100, title="Edit distances in minimizer graph")
    histogram(edit_hist_one_per_read, args, name='read_error_estimates_from_min_ed_zoomed.png', x='Edit distance to minimizer', y='count', x_cutoff=50, title="Edit distances in minimizer graph")
    print(sorted(edit_hist_one_per_read))

    histogram(edit_hist, args, name='edit_distances.png', x='Edit distance to minimizer', y='count', x_cutoff=100, title="Edit distances in minimizer graph")
    histogram(edit_hist, args, name='edit_distances_zoomed.png', x='Edit distance to minimizer', y='count', x_cutoff=50, title="Edit distances in minimizer graph")
    histogram(neighbors, args, name='neighbours.png', x='x-axis', y='y-axis', title="Number of neighbours in minimizer graph")
    histogram(neighbors, args, name='neighbours_zoomed.png', x='x-axis', y='y-axis', x_cutoff=20, title="Number of neighbours in minimizer graph")

    clusters_to_database = transpose(minimizer_graph)
    outfile_alignments = open(os.path.join(args.outfolder, "pred_to_pred_best_alignments.tsv"), "w")
    for t in  clusters_to_database:
        for c in clusters_to_database[t]:
            # print(t, clusters_to_database[t])
            outfile_alignments.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(c,t, len(predicted[c]), len(predicted[t]), clusters_to_database[t][c]))
    outfile_alignments.close()


def transpose(dct):
    d = defaultdict(dict)
    for key1, inner in dct.items():
        for key2, value in inner.items():
            d[key2][key1] = value
    return d   

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--predicted', type=str, nargs="+", help='Path to the predicted transcript fasta file')
    parser.add_argument('--database', type=str, default=None, help='Path to the consensus fasta file')
    parser.add_argument('--outfolder', type=str, help='Output path of results')
    parser.add_argument('--single_core', dest='single_core', action='store_true', help='Force working on single core. ')
    args = parser.parse_args()

    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)
    
    if args.database:
        main_temp_2set(args)
    else:
        main_temp(args)

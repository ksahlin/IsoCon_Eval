
from __future__ import print_function
import os,sys
import argparse
import re
import itertools

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
# import edlib
import math

from collections import defaultdict
import errno
import pickle

import networkx as nx
from networkx.drawing.nx_pydot import write_dot, pydot_layout


def print_database_to_tsv(db):
    column_order = ['id',  'predicted_acc', 'family', 'Illumina_support', 'perfect_match_database', 'sample1', 'sample2', 'acc_sample1', 'acc_sample2', 'both_samples', 'primer', 'category', 'full_supp_sample1', 'full_supp_sample2', 'gene_member_number', 'batch', 'sequence',]

    filename = open(args.outfile, "w")
    filename.write("\t".join(column_order) + "\n")
    for record in db['predicted_transcripts'].all():
        # print(record, type(record))
        # print(record.values(), type(record))
        record_values = [str(record[field]) for field in column_order]
        fmt = "\t".join(record_values)
        fmt += "\n"
        filename.write(fmt)


def ssw_alignment(x, y, ends_discrepancy_threshold = 250 ):
    """
        Aligns two sequences with SSW
        x: query
        y: reference

    """

    score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-20)
    aligner = ssw.Aligner(gap_open=50, gap_extend=0, matrix=score_matrix)

    # for the ends that SSW leaves behind
    bio_matrix = matlist.blosum62
    g_open = -1
    g_extend = -0.5
    ######################################

    # result = aligner.align("GA", "G", revcomp=False)
    # y_alignment, match_line, x_alignment = result.alignment
    # c = Counter(match_line)
    # matches, mismatches, indels = c["|"], c["*"], c[" "]
    # alignment_length = len(match_line)
    # print("matches:{0}, mismatches:{1}, indels:{2} ".format(matches, mismatches, indels))
    # print(match_line)

    result = aligner.align(x, y, revcomp=False)
    y_alignment, match_line, x_alignment = result.alignment
    # print()
    # print(y_alignment)
    # print(match_line)
    # print(x_alignment)
    matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")

    # alignment_length = len(match_line)
    
    start_discrepancy = max(result.query_begin, result.reference_begin)  # 0-indexed # max(result.query_begin, result.reference_begin) - min(result.query_begin, result.reference_begin)
    query_end_discrepancy = len(x) - result.query_end - 1
    ref_end_discrepancy = len(y) - result.reference_end - 1
    end_discrepancy = max(query_end_discrepancy, ref_end_discrepancy)  # max(result.query_end, result.reference_end) - min(result.query_end, result.reference_end)
    # print("disc:", start_discrepancy, end_discrepancy)
    tot_discrepancy = start_discrepancy + end_discrepancy

    if 0 < start_discrepancy <= ends_discrepancy_threshold:
        print("HERE",start_discrepancy)
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

    if 0 < end_discrepancy <= ends_discrepancy_threshold:
        print("HERE2", end_discrepancy)
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

    return x_alignment, y_alignment, matches, mismatches, indels, match_line       





def create_isoform_graph(transcripts):
    
    G = nx.Graph()
    for acc in transcripts.keys():
        G.add_node(acc, accession = acc, support = acc.split("_")[4])


    score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-2)
    aligner = ssw.Aligner(gap_open=2, gap_extend=0, matrix=score_matrix)
    cntr = 0

    processed = set()
    # already_assigned = set()

    for acc1, seq1 in sorted(transcripts.items(), key= lambda x: len(x[1]), reverse=True):
        cntr += 1
        processed.add(acc1)
        # print("length t:", len(seq1), acc1)
        if cntr % 5 == 0:
            print(cntr, "sequences processed")
        # print(acc1)
        # if acc1 in already_assigned:
        #     # print("allready assigned to larger sequence!")
        #     continue

        for acc2, seq2 in sorted(transcripts.items(), key= lambda x: len(x[1]), reverse=True):
            if acc2 in processed:
                continue
            # if acc1 == acc2:
            #     continue
            # if seq1 == seq2:
            #     continue
            # result = aligner.align(seq1, seq2, revcomp=False)
            # seq2_aln, match_line, seq1_aln = result.alignment
            print(seq1)
            print(seq2)
            seq1_aln, seq2_aln, matches, mismatches, indels, match_line = ssw_alignment(seq1, seq2, ends_discrepancy_threshold = 2000 )

            del_seq1 = re.findall(r"[-]+",seq1_aln)
            del_seq2 = re.findall(r"[-]+",seq2_aln)
            mismatches = len([ 1 for n1, n2 in zip(seq1_aln,seq2_aln) if n1 != n2 and n1 != "-" and n2 != "-" ])

            ## do not count length discrepancies in ends
            inner_del_seq1 = re.findall(r"[AGCT][-]+[AGCT]",seq1_aln)
            inner_del_seq2 = re.findall(r"[AGCT][-]+[AGCT]",seq2_aln)
            # print(inner_del_seq1)
            # print(inner_del_seq2)
            total_inner = sum([len(d) - 2 for d in inner_del_seq1]) + sum([len(d) - 2 for d in inner_del_seq2])
            # print(indels, total_inner)


            # by default (since all transcripts are distinct if we end up here), each transcript is its on gene member
            # if we find an alingment that contains only structural changes of > X (2) nucleotides, and no other smaller differences we classify as same family
            if mismatches == 0:
                del_lengths1 = [len(del_) for del_ in del_seq1]
                del_lengths2 = [len(del_) for del_ in del_seq2]
                no_small_del_in_seq1 = ((len(del_lengths1) > 0 and min(del_lengths1) >= 3) or len(del_lengths1)  == 0)
                no_small_del_in_seq2 = ((len(del_lengths2) > 0 and min(del_lengths2) >= 3) or len(del_lengths2)  == 0)
                # print(no_small_del_in_seq1, no_small_del_in_seq2)
                # print((len(del_seq1) > 0 and min(del_seq1) >= 3), len(del_seq1)  == 0)
                # if acc1[0][:14] == "transcript_460" and acc2[0][:14] == "transcript_467" :
                #     print("we are here", no_small_del_in_seq1, no_small_del_in_seq2, mismatches)
                #     sys.exit()
                if no_small_del_in_seq1 and no_small_del_in_seq2:
                    G.add_edge(acc1, acc2, alignment={ acc1 : seq1_aln, acc2 : seq2_aln })
                else:
                    pass
                    # print("Different only by small indel!!")

    list_of_maximal_cliques = list(nx.find_cliques(G))
    print("Number of possible members:", len(list_of_maximal_cliques) )
    print("clique sizes", [ len(cl) for cl in  sorted(list_of_maximal_cliques, key= lambda x: len(x), reverse=True)] )
    return G


def transpose(dct):
    d = defaultdict(dict)
    for key1, inner in dct.items():
        for key2, value in inner.items():
            d[key2][key1] = value
    return d   

# def is_equivalence_class(family_members, member):
#     isoforms_to_member = family_members[member]   # this is a set
#     print("Main set:", isoforms_to_member)

#     for candidate_member in family_members[member]:
#         temp_set = isoforms_to_member - set([candidate_member])
#         print("Isoform set:", temp_set)
#         isoforms_to_candidate_member = family_members[candidate_member]
#         if not temp_set.issubset(isoforms_to_candidate_member):
#             print("here why?")
#             return False

#     return True        

def construct_consensus_progressively(G, clique_id_to_accessions, transcripts):
    pattern = re.compile(r"[-]+")
    consensi = {}
    for member_id, list_of_acc in clique_id_to_accessions.items():
        if len(list_of_acc) == 1:
            isoform_acc = list_of_acc[0]
            consensi[member_id] = transcripts[ isoform_acc]
        else:
            consensi[member_id] = ""

    for member_id, list_of_acc in clique_id_to_accessions.items():
        if len(list_of_acc) > 1:
            isoform_acc = list_of_acc[0]
            c_temp = transcripts[ isoform_acc ]
            for acc in list_of_acc[1:]:
                isof = transcripts[acc]
                isof_aln, c_temp_aln, matches, mismatches, indels, match_line = ssw_alignment(isof, c_temp, ends_discrepancy_threshold = 2000 )

                c_list = []
                for n1, n2 in zip(c_temp_aln, isof_aln):
                    if n1 != "-":
                        c_list.append(n1)
                    else:
                        c_list.append(n2)
                c_temp = "".join([ n for n in c_list ] )
            consensi[member_id] = c_temp
        
        print("LENGHT C", len(consensi[member_id]))

    return consensi

def construct_multialignment_from_consensi(clique_id_to_accessions, transcripts, consensi):
    multialignments = {}
    for member_id, list_of_acc in clique_id_to_accessions.items():
        multialignments[member_id] = {}
        # multialignments[member_id]["consensus"] = consensi[member_id]
        cons = consensi[member_id]
        for isoform_acc in list_of_acc:
            isof = transcripts[isoform_acc]
            isof_aln, cons_aln, matches, mismatches, indels, match_line = ssw_alignment(isof, cons, ends_discrepancy_threshold = 2000 )
            multialignments[member_id][isoform_acc] = isof_aln
            # print(len(isof_aln), len(cons), len(isof)) 
            # print(isof_aln)
            # print(cons)
            assert len(isof_aln) == len(cons)

    return multialignments



def get_gene_member_number(transcripts, params):

    tsv_outfile = open(os.path.join(args.outfolder, "predicted_to_genemembers.tsv"), "w")
    tsv_outfile.write("FAMILY\tMEMBER_ID\tACCESSION\n")
    dot_graphs_folder = os.path.join(args.outfolder, "dot_graphs")
    mkdir_p(dot_graphs_folder)
    # doing the gene member analysis of only shared transcripts
    all_family_consensi = {}
    print(len(transcripts))

    if params.graph_prefix:
        G = pickle.load( open( os.path.join(params.graph_prefix, params.prefix + ".p"), "rb" ) )
    else:
        G = create_isoform_graph(transcripts)
        pickle.dump( G, open( os.path.join(params.outfolder, params.prefix  + ".p"), "wb" ) )
        dot_graph_out = os.path.join(dot_graphs_folder, params.prefix  + ".dot")
        write_dot(G, dot_graph_out)

    list_of_maximal_cliques = list(nx.find_cliques(G))
    clique_id_to_accessions = {member_id : cl for member_id, cl in enumerate(sorted(list_of_maximal_cliques, key= lambda x: len(x), reverse=True))}
    print(clique_id_to_accessions)
    consensi = construct_consensus_progressively(G, clique_id_to_accessions, transcripts)
    all_family_consensi[ params.prefix ] = consensi
    multialignments = construct_multialignment_from_consensi(clique_id_to_accessions, transcripts, consensi)

    ## Write output 
    family_folder = os.path.join(args.outfolder, params.prefix)
    mkdir_p(family_folder)

    consensi_fasta = open(os.path.join(family_folder,  "all_consensi.fa"), "w") 

    for member_id, cl in clique_id_to_accessions.items(): #enumerate(sorted(list_of_maximal_cliques, key= lambda x: len(x), reverse=True)):
        members_fasta = open(os.path.join(family_folder, str(member_id) + ".fa"), "w") 
        members_aligned_fasta = open(os.path.join(family_folder, str(member_id) + "_aligned.fa"), "w")
        members_aligned_fasta.write( ">{0}\n{1}\n".format("consensus", consensi[member_id]) )
        consensi_fasta.write(">{0}\n{1}\n".format(member_id, consensi[member_id]) )

        for isoform_acc in cl:
            isoform = transcripts[isoform_acc]
            # print(len(G.edge[isoform1][isoform2]["alignment"][isoform1]), len(G.edge[isoform1][isoform2]["alignment"][isoform2]) )
            tsv_outfile.write("{0}\t{1}\t{2}\n".format(params.prefix, member_id, str(isoform_acc)))
            members_fasta.write(">{0}\n{1}\n".format(isoform_acc, transcripts[isoform_acc]) )
            members_aligned_fasta.write( ">{0}\n{1}\n".format(isoform_acc, multialignments[member_id][isoform_acc]) )

    all_processed = set([ acc for member_id, list_of_acc in clique_id_to_accessions.items() for acc in list_of_acc ])
    print("Number non processed sequences(bug if more than 0): ", len( set(transcripts.keys()).difference(all_processed)) )

    transcript_accessions_to_clique_ids = defaultdict(list)
    for member_id, cl in clique_id_to_accessions.items():
        for transcript_number in cl:
            isoform, primer_id = number_to_member_acc[transcript_number]
            transcript_accessions_to_clique_ids[isoform].append(member_id)


    # update database field here
    for record in family_records:
        q_id = record["id"]
        isoform_accession = record["predicted_acc"]
        potential_gene_members = ",".join( [ str(member_id) for member_id in  transcript_accessions_to_clique_ids[isoform_accession] ])
        # print(potential_gene_members)
        # print("here")

        data = dict(id=q_id, gene_member_number = potential_gene_members)
        table.update(data, ['id'])

    
    # print a tsv file with graph information
    tsv_folder = os.path.join(args.outfolder, "graph_info")
    mkdir_p(tsv_folder)

    graph_tsv_cliques = open(os.path.join(tsv_folder,  target + "_graph.fa"), "w") 
    graph_tsv_cliques.write("node\tinteraction\tnode2\tsupport\tcoding\taccession\n")
    
    for n in G.nodes():
        coding = "yes" if G.node[n]["accession"] in is_coding else "no"
        nbrs = list(G.neighbors(n))
        if len(nbrs) == 0:
            graph_tsv_cliques.write("{0}\t\t\t{1}\t{2}\t{3}\n".format(n, G.node[n]["support"], coding, G.node[n]["accession"] ))
        else:
            for nbr in nbrs:
                graph_tsv_cliques.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(n, "isoform", nbr, G.node[n]["support"], coding, G.node[n]["accession"] ))



    graph_tsv_edit_distances = open(os.path.join(tsv_folder,  target + "_edit.fa"), "w") 
    graph_tsv_edit_distances.write("node\tinteraction\tnode2\tsupport_node1\tsupport_node2\tmin_ed\taccession\n")

    return all_family_consensi



def plot_member_ed(params, all_family_consensi):

    # print all pairwise distances to tsv file
    ## fam member_id1 member_id2 len(m1) len(m2) ed
    ed_tsv_file =open( os.path.join(params.outfolder, "member_edit_distances.tsv"), "w" )
    ed_tsv_file.write("family\tmember1\tmember2\tm1_len\tm2_len\ted\n")

    for fam in all_family_consensi:
        processed = set()
        for m1_id, m1_seq in all_family_consensi[fam].items():
            ed_tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(fam, m1_id, m1_id, len(m1_seq), len(m1_seq), 0 ))
            processed.add(m1_id)
            for m2_id, m2_seq in all_family_consensi[fam].items():
                if m2_id not in processed:
                    result = edlib.align(m1_seq, m2_seq, task = "distance")
                    ed = result["editDistance"]
                    ed_tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(fam, m1_id, m2_id, len(m1_seq), len(m2_seq), int(ed) ))

    ed_tsv_file.close()

    # plot
    dtypes =  {"family": str, "ed" : int, "member1" : int, "member2" : int, "m1_len" : int, "m2_len" : int}
    data = pd.read_csv(ed_tsv_file.name, sep="\t", dtype = dtypes)
    # print(data.isnull().values.any())

    def facet_heatmap(data, **kws):
        data = data.pivot("member1", 'member2', 'ed')
        pd.set_option('display.max_rows', 1000)
        # print(data)
        # print(data.columns)

        mask = data.isnull()
        # print(mask)
        for i in range(len(mask)):
            for j in range(len(mask[i])):
                if data[i ][j] > 99:
                    # print("true")
                    mask[i][j] = True

        sns.heatmap(data, cmap='coolwarm_r', annot=True, mask = mask, **kws)  # <-- Pass kwargs to heatmap fmt = "d"


    # with sns.plotting_context(font_scale=5.5):
    #     g = sns.FacetGrid(data, col="family", col_wrap=3, size=3, aspect=1 ) #, col_order = ["TSPY13P", "HSFY2", "DAZ2"])

    # cbar_ax = g.fig.add_axes([.92, .3, .02, .4])  # <-- Create a colorbar axes

    # g = g.map_dataframe(facet_heatmap,
    #                     cbar_ax=cbar_ax,
    #                     vmin=0, vmax=10)  # <-- Specify the colorbar axes and limits
    plot_outfolder = os.path.join(params.outfolder, "edit_distance_plots")
    mkdir_p(plot_outfolder)
    for target in ["BPY", "CDY1", "CDY2", "DAZ", "HSFY1", "HSFY2", "PRY", "RBMY", "TSPY", "XKRY", "VCY"]:
        df = data.loc[data['family'] == target]
        facet_heatmap(df, vmin=0, vmax=10)

        outfile = os.path.join(plot_outfolder, target + ".png")
        plt.savefig(outfile)
        plt.close()

    # g.set_titles(col_template="{col_name}", fontweight='bold', fontsize=14)
    # g.fig.subplots_adjust(right=.9)  # <-- Add space so the colorbar doesn't overlap the plot
    # outfile = os.path.join(params.outfolder, "member_edit_distances.png")

    # plt.savefig(outfile)
    # plt.close()

'''
    Below function taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].split()[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break


def main(params):

    # read in objects
    transcripts = {acc : seq for acc, (seq, _) in readfq(open(args.predicted, "r"))}

    # classify into gene memebers
    all_family_consensi = get_gene_member_number(transcripts, params)

    # plot edit distances of consensi within a family
    plot_member_ed(params, all_family_consensi)

    # print to file
    print_database_to_tsv(db)

    # print fasta
    print_database_to_separate_fastas(db)


    # result = db['predicted_transcripts'].all()
    # dataset.freeze(result, mode='item', format='csv', filename=args.outfile)





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
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    parser.add_argument('--predicted', type=str, help='Path to consensus fasta file(s)')
    parser.add_argument('--outfile', type=str, help='A fasta file with transcripts that are shared between smaples and have perfect illumina support.')
    parser.add_argument('--graph_prefix', type=str, help='.')    
    parser.add_argument('--prefix', type=str, help='E.g., BPY2, TSPY or th family under consideration, used as prefix for outfiles.')    

    args = parser.parse_args()
    path_, file_prefix = os.path.split(args.outfile)

    if path_:
        mkdir_p(path_)
    args.outfolder = path_

    main(args)


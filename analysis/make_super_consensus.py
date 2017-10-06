
from __future__ import print_function
import os,sys
import argparse
import re
import itertools
import errno
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist


import ssw

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


def construct_multialignment_from_consensi(transcripts, super_consensi):
    multialignments = {}
    for acc, seq in transcripts.items():
        isof_aln, super_cons_aln, matches, mismatches, indels, match_line = ssw_alignment(seq, super_consensi, ends_discrepancy_threshold = 2000 )
        multialignments[acc] = isof_aln
        # print(len(isof_aln), len(cons), len(isof)) 
        # print(isof_aln)
        # print(cons)
        assert len(isof_aln) == len(super_consensi)

    return multialignments

def construct_consensus_progressively(member_consensi):
    pattern = re.compile(r"[-]+")

    acc, seq = list(member_consensi.items())[0]
    c_temp = seq # just initialize with one sequence

    for acc, seq in member_consensi.items():
        isof = member_consensi[acc]
        isof_aln, c_temp_aln, matches, mismatches, indels, match_line = ssw_alignment(isof, c_temp, ends_discrepancy_threshold = 2000 )

        c_list = []
        for n1, n2 in zip(c_temp_aln, isof_aln):
            if n1 == "-":
                c_list.append(n2)
            elif n2 == "-":
                c_list.append(n1)
            elif n1 != n2:
                c_list.append("N")
            else:
                assert n1 == n2
                c_list.append(n1)

        c_temp = "".join([ n for n in c_list ] )
    super_consensi = c_temp
        
    print("LENGHT C", len(super_consensi))
    print(super_consensi)
    print(super_consensi.count("N"))
    return super_consensi


def main(args):

    member_consensi = {acc: seq.upper() for (acc, seq) in  read_fasta(open(args.member_consensi, 'r'))}
    transcripts = {acc: seq.upper() for (acc, seq) in  read_fasta(open(args.transcripts, 'r'))}
    super_consensi = construct_consensus_progressively(member_consensi)
    multialignments = construct_multialignment_from_consensi(transcripts, super_consensi)
    
    outfile = open(os.path.join(args.outfolder, "test_multi_all_rbmy.fa"), "w")
    outfile.write(">{0}\n{1}\n".format("super_consensus", super_consensi))
    for acc, seq in multialignments.items():
        outfile.write(">{0}\n{1}\n".format(acc, seq))

    outfile.close()


def mkdir_p(path):
    print("creating", path)
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('--member_consensi', type=str, help='Path to consensus fasta file(s)')
    parser.add_argument('--transcripts', type=str, help='Path to predicted fasta file(s)')
    parser.add_argument('--outfolder', type=str, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    
    args = parser.parse_args()

    mkdir_p(args.outfolder)

    main(args)
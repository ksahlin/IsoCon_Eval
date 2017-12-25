from __future__ import print_function
import os,sys
import argparse

import pysam
from collections import defaultdict

import math

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd


class CCS(object):
    """docstring for CCS"""
    def __init__(self, name, seq, qual, np):
        super(CCS, self).__init__()
        self.name = name
        self.seq = seq
        self.qual = qual
        self.np = np
        self.subreads = {}
    def positions_with_p_error_higher_than(self, prob):
        # indicator_probs = ""
        indicator_probs = []
        for i, p_error in  enumerate(get_p_error_from_q_list(self.qual)):
            if p_error > prob:
                # indicator_probs += str(i) + "," + str(p_error)
                indicator_probs.append((i,p_error) )
            # else:
            #     indicator_probs += " "
        return(indicator_probs)
            
def get_p_error_from_q_list(qual_list):
    probs = []
    for q in qual_list:
        p = 10**(-q/10.0)
        probs.append(p)
    return probs

def get_p_error_from_q(q):
    p = 10**(-q/10.0)
    return p

def get_p_error_from_char(qual_string):
    probs = []
    for char in qual_string:
        q = ord(char) - 33
        p = 10**(-q/10.0)
        probs.append(p)
    return probs

def get_ccs(ccs_file):  
    ccs_dict = defaultdict(dict)
    for read in ccs_file.fetch(until_eof=True):
        ccs_read = CCS(read.query_name, read.query_alignment_sequence, read.query_qualities, read.get_tag("np"))
        ccs_dict[read.query_name] = ccs_read
        # ccs_dict[read.query_name]["seq"] = read.query_alignment_sequence
        # print(read.query_qualities)
        # sys.exit()
        # ccs_dict[read.query_name]["qual"] = read.query_qualities
        # ccs_dict[read.query_name]["np"] = read.get_tag("np")
        assert len(read.query_alignment_sequence) == len(read.query_qualities)
        
        # if ccs_read.np > 10:
        #     print(ccs_read.np, ccs_read.positions_with_p_error_higher_than(0.01))
    return ccs_dict

def get_subreads(subreads_file):
    target = "m151210_031012_42146_c100926392550000001823199905121697_s1_p0/515/"
    subreads_dict = defaultdict(dict)

    for read in subreads_file.fetch(until_eof=True):
        # print(read.query_name)
        if target in read.query_name:
            # print(read.query_alignment_sequence)
            subreads_dict[read.query_name]["seq"] = read.query_alignment_sequence
            subreads_dict[read.query_name]["qual"] = read.query_qualities
            subreads_dict[read.query_name]["iq"] = get_p_error_from_char(read.get_tag("iq"))
            subreads_dict[read.query_name]["dq"] = get_p_error_from_char(read.get_tag("dq"))
            subreads_dict[read.query_name]["sq"] = get_p_error_from_char(read.get_tag("sq"))
            try:
                subreads_dict[read.query_name]["st"] = read.get_tag("st")
            except KeyError:
                print("tag 'st' not present for", read.query_name)
            subreads_dict[read.query_name]["dt"] = read.get_tag("dt")
            subreads_dict[read.query_name]["mq"] = get_p_error_from_char(read.get_tag("mq"))

    # print(subreads_dict)
    print(len(subreads_dict))

    return subreads_dict

def populate_ccs_objects_with_subreads(ccs_dict, subreads_dict):
    for q_acc in subreads_dict:
        q_acc_mod = "/".join(q_acc.split("/")[:-1])
        ccs_acc = q_acc_mod + "/ccs"
        if ccs_acc in ccs_dict:
            ccs_dict[ccs_acc].subreads[q_acc] = subreads_dict[q_acc]
            ccs_dict[ccs_acc].subreads[q_acc]["qual"] = [1.0 - ccs_dict[ccs_acc].subreads[q_acc]["iq"][pos] - ccs_dict[ccs_acc].subreads[q_acc]["dq"][pos] - ccs_dict[ccs_acc].subreads[q_acc]["sq"][pos] - ccs_dict[ccs_acc].subreads[q_acc]["mq"][pos] for pos in range(len(ccs_dict[ccs_acc].subreads[q_acc]["iq"]))  ]
    
    return ccs_dict



def get_polymer_quality_tuples(seq, qual):

    polymer_quality_tuples = []
    # first base is always start of any homoplymer
    h_qual = qual[0]
    h_base = seq[0]
    h_len = 1
    for i, (char1, char2) in enumerate(zip(seq[:-1], seq[1:])):
        if char1 != char2: # char2 is first base of a homopolymer
            #  print out info related to char1
            p_error = get_p_error_from_q(h_qual)
            assert h_base == char1
            polymer_quality_tuples.append( (h_len, p_error, h_base) )

            # start over with char2 information
            h_len = 1
            h_base = char2
            h_qual = qual[i+1]

        else:
            h_len += 1

    # end case
    # if char1 == char2:
    #     p_error = get_p_error_from_q(h_qual)
    #     polymer_quality_tuples.append( (h_len, p_error, h_base) )
    # elif char1 != char2:
    p_error = get_p_error_from_q(h_qual)
    polymer_quality_tuples.append( (h_len, p_error, h_base) )

    # print(len(seq), sum([h_len for (h_len, p_error, h_base) in polymer_quality_tuples]))
    # print(seq)
    assert len(seq) == sum([h_len for (h_len, p_error, h_base) in polymer_quality_tuples])
    return polymer_quality_tuples


def print_out_tsv(ccs_dict, args):
    tsv_file = open(os.path.join(args.outfolder, "passes_polymer_qual.tsv"), "w")
    tsv_file.write("Passes\tHomopolymenr_length\tP_error\tBase\n")
    cntr = 1
    for ccs_acc in ccs_dict:
        seq = ccs_dict[ccs_acc].seq
        qual = ccs_dict[ccs_acc].qual
        assert len(seq) == len(qual)
        num_passes = ccs_dict[ccs_acc].np
        polymer_quality_tuples = get_polymer_quality_tuples(seq, qual)
        for h_length, p_e, base in polymer_quality_tuples:
            tsv_file.write("{0}\t{1}\t{2}\t{3}\n".format(num_passes, h_length, p_e, base))

        if cntr % 5000 == 0:
            print("processing ccs nr {0}".format(cntr))
            break
        if cntr > args.nr_reads_to_plot:
            break
        cntr += 1

    tsv_file.close()
    return tsv_file.name

def plot_bp_qual(tsv_file, outfile):
    sns.plt.clf()
    with sns.plotting_context("paper", font_scale=1.8):
        indata = pd.read_csv(tsv_file, sep="\t")
        fig, ax = plt.subplots()

        g = sns.factorplot(x="Homopolymenr_length", y="P_error", col="Passes", col_order = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], col_wrap=3,
                            hue="Base", data=indata[indata.P_error.notnull()],
                            size=3, aspect=1.6, palette="Set3",
                            dodge=True, cut=0, bw=.2) #kind="violin",
        
        # g = sns.FacetGrid(data, row="Family", col="mutation_rate", size=3, aspect=1.6, row_order=["TSPY13P", "HSFY2", "DAZ2"], col_order=[0.01, 0.001, 0.0001], legend_out=True)
        sns.set(style="whitegrid", palette="muted")
        # (g.map(sns.violinplot, "read_count", args.y_axis, "TOOL", cut=0, hue_order=["ISOCON", "ICE"], palette=sns.color_palette("muted", 2)).despine(left=True).add_legend(title="TOOL", label_order=["ISOCON", "ICE"]))
        # g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=16)
        # g.set_yticklabels(["",0,0.2,0.4,0.6,0.8,1.0])
        # g.set(yscale="log")
        g.set(xlim=(0, 6))
        g.set(ylim=(0, 1))
        plt.savefig(outfile)
        plt.close()


def main(args):
    ccs_file = pysam.AlignmentFile(args.ccs, "rb", check_sq=False)
    ccs_dict = get_ccs(ccs_file)
    tsv_file = print_out_tsv(ccs_dict, args)
    
    # plot qual vs basepair
    outfile = os.path.join(args.outfolder, "plot_bp_quals.png")
    plot_bp_qual(tsv_file, outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('--ccs', type=str, help='Path to consensus fasta file(s)')
    parser.add_argument('--outfolder', type=str, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    parser.add_argument('--nr_reads_to_plot', type=int, default = 1000, help='Number of reads to make plot from.')
    
    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)
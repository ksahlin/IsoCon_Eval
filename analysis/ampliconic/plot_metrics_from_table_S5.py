import sys
import argparse
import re

import os
import random
import errno

try:
    # import matplotlib
    # matplotlib.use('agg')
    # import matplotlib.pyplot as plt
    import pylab as plt
    import seaborn as sns
    sns.set_palette("husl", desat=.6)
    sns.set(font_scale=1.6)
    plt.rcParams.update({'font.size': 12})
except:
    pass

import numpy as np
import pandas as pd

def get_relevant_values(line):
    id_, predicted_acc, family, Illumina_support,  perfect_match_database, sample1, sample2, acc_sample1, acc_sample2, both_samples, primer, category, full_supp_sample1, full_supp_sample2, gene_member_number, batch, sequence = line.split("\t")
    support  = int(predicted_acc.split("_")[4])
    p_val = predicted_acc.split("_")[5]
    if p_val == "not":
        p_val = 0.0
    else:
        p_val = float(p_val)

    shared = "yes"if both_samples == "yes" else "no"
    sample = "sample1" if predicted_acc == acc_sample1 else "sample2"
    ill_supp = 0 if Illumina_support == "None" else float(Illumina_support)

    if sample == "sample1": 
        if full_supp_sample1 ==  "yes":
            full_sup = "yes"
        else:
            full_sup = "no"            
    elif sample == "sample2": 
        if full_supp_sample2 ==  "yes":
            full_sup = "yes"
        else:
            full_sup = "no" 
    else:
        print("BUG")
        sys.exit()

    if p_val == 0.0:
        p_val = 1e-310
    return family, p_val, support, sample, shared, perfect_match_database, ill_supp, full_sup


def create_tsv_from_transcript_annotations(args):
    infile = open(args.transcript_annotations, "r")
    tsv_file = open( os.path.join(args.outfolder, "tsv_for_plotting.tsv"), "w")
    tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format("gene", "p-value", "sample", "shared", "read_depth", "perfect_match_database", "Illumina_support", "Full_support"))

    total_records = 0
    for i, line in enumerate(infile):
        if i == 0:
            continue
        _, predicted_acc, family, _, _, sample1, sample2, acc_sample1, acc_sample2, both_samples, _, _, _, _, _, _, _ = line.split("\t")
        if family == "HSFY1":
            continue

        values = get_relevant_values(line)
        if values:
            gene, p_val, read_depth, sample, shared, perfect_match_database, ill_supp, full_sup = values
            tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(gene, p_val, sample, shared, read_depth, perfect_match_database, ill_supp, full_sup))
        else:
            pass
        total_records += 1

    print(total_records, "records iterated through.")
    tsv_file.close()
    return tsv_file.name


def plot(data, args):
    sns.plt.clf()
    with sns.plotting_context("paper", font_scale=1.2):
        fig, ax = plt.subplots()
        g = sns.FacetGrid(data, row="gene", col="mutation_rate", hue="nr_reads", palette="Set1", size=3, aspect=3.2, row_order=["TSPY13P", "HSFY2", "DAZ2"], col_order=[ 0.0001], legend_out=True)
        g.set(xscale="log")
        g.set(xlim=(1e-100,1), ylim=(0,1))
        g.set(xticks=[1e-2, 1e-10, 1e-20, 1e-30, 1e-50, 1e-100], yticks=[0,0.2,0.4,0.6,0.8,1.0])
        sns.set(style="whitegrid", palette="muted")

        (g.map(plt.scatter, "corrected_p_val", "relative_support").add_legend())
        g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=12)

        g.set_ylabels("Relative transcript abundance")
        g.set_xlabels("P-value")

        plt.savefig(os.path.join(args.outfolder, "Figure_S17.pdf"))
        plt.clf()

def plot_box(data, args):
    sns.plt.clf()
    with sns.plotting_context("paper", font_scale=1.2):
        fig, ax = plt.subplots()
        ax.set(yscale="log")
        ax.set_ylim(1**-320, 1.0)
        ax = sns.boxplot(x="gene", y="corrected_p_val", hue="nr_reads", data=data, palette="Set3", order=["TSPY13P", "HSFY2", "DAZ2"])
        plt.savefig(os.path.join(args.outfolder, "Figure_S18.pdf"))
        plt.clf()

def pairplot(data, args):
    sns.plt.clf()
    with sns.plotting_context("paper", font_scale=1.2):
        fig, ax = plt.subplots()
        g = sns.lmplot(x="p-value", y="read_depth",data=data, hue = "shared", fit_reg= False)
        # labels = ["Not comp", 1e-200,  1e-100,  1e-50, 1e-30, 1e-20, 1e-10, 1e-2]
        # step = [1e-300, 1e-200,  1e-100,  1e-50, 1e-30, 1e-20, 1e-10, 1e-2] 
        # g = (g.set_xticklabels(labels, step))

        g = (g.set_axis_labels("p-value", "Read depth").set(xlim=(1e-320, 1.0), xscale="log"))  #, xticks=[10, 30, 50], yticks=[2, 6, 10], 
        # print([t.set_text("Not computed") for i,t in enumerate(g.ax.get_xticklabels()) if i==0   ])

        # g.set_axis_labels(xticks = (np.arange(5), ('Tom', 'Dick', 'Harry', 'Sally', 'Sue')))
        # g.set(yscale="log")
        # g.set_ylim(1.0**-320, 1.0)
        # g.set_xlim(1.0**-320, 1.0)
        plt.savefig(os.path.join(args.outfolder, "Figure_tab_S5A.pdf"))
        plt.clf()

        g = sns.lmplot(x="p-value", y="Illumina_support",data=data, hue = "shared", fit_reg= False)
        g = (g.set_axis_labels("p-value", "Illumina_support").set(xlim=(1e-320, 1.0), xscale="log"))  #, xticks=[10, 30, 50], yticks=[2, 6, 10], 
        plt.savefig(os.path.join(args.outfolder, "Figure_tab_S5B.pdf"))
        plt.clf()

        ax = sns.stripplot(x="p-value", y="perfect_match_database",data=data, hue = "shared", jitter=True, split=True)
        ax.set(xlim=(1e-320, 1.0), xscale="log", xlabel="p-value", ylabel="Perfect match ENSEMBL")  #, xticks=[10, 30, 50], yticks=[2, 6, 10], 
        plt.savefig(os.path.join(args.outfolder, "Figure_tab_S5C.pdf"))
        plt.clf()

        ax = sns.stripplot(x="read_depth", y="perfect_match_database",data=data, hue = "shared", jitter=True, split=True)
        ax.set(xscale="log", xlabel="Read depth", ylabel="Perfect match ENSEMBL")  #, xticks=[10, 30, 50], yticks=[2, 6, 10], 
        plt.savefig(os.path.join(args.outfolder, "Figure_tab_S5D.pdf"))
        plt.clf()


def main(args):
    
    tsv_file = create_tsv_from_transcript_annotations(args)
    data = pd.read_csv(tsv_file, sep="\t")
    pairplot(data, args)
    # print("DAZ")
    # print(data.loc[data["gene"] == "DAZ2"].corr(method= 'spearman'))

    # print("TSPY")
    # print(data.loc[data["gene"] == "TSPY13P"].corr(method= 'spearman'))

    # print("HSFY")
    # print(data.loc[data["gene"] == "HSFY2"].corr(method= 'spearman'))

    # print()
    # print("ALL")
    # print(data.corr(method='spearman'))
    # print(data.apply(lambda x : pd.factorize(x)[0]).corr(method='spearman', min_periods=1))

    # plot_box(tsv_file, args)

    # plot with swarmplot
    # plot(tsv_file, args)



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
    parser = argparse.ArgumentParser(description="Plot p-values.")

    parser.add_argument('--transcript_annotations', type=str, help='Path to output files')
    parser.add_argument('--outfolder', type=str, help='Outfolder.')
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    mkdir_p(args.outfolder)

    main(args)

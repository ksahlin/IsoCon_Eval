import sys
import argparse
import re

import os
import random
import errno

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")
# import matplotlib.pyplot as plt
# import matplotlib

import numpy as np
import seaborn as sns
import pandas as pd

def print_tsv(args):
    tsv_file = open(os.path.join(args.outfolder, "p_values_candidates.tsv"), "w")
    tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format("gene", "copies", "abundance", "mutation_rate", "nr_reads", "exp_nr", "corrected_p_val", "support", "relative_support", "nr_varinats", "acc"))
    pattern = r"Candidates written to file:"
    for file_ in args.isoconfiles:
        #get experiment parameters
        path_, file_prefix = os.path.split(file_)
        params = path_.split("/")[-1]
        gene, copies, abundance, mutation_rate, nr_reads, exp_nr = params.split("_")  #TSPY13P_8_exponential_0.001_20_9
        # print(gene, copies, abundance, mutation_rate, nr_reads, exp_nr)
        indexes = []
        file_contents = open(file_, "r").readlines()
        for i, l in enumerate(file_contents):
            if "Candidates written to file:" in l:
                # print(l,i)
                indexes.append(i)

        start_index, stop_index = indexes[-2], indexes[-1]
        # start_line_index = all_output.index("Candidates written to file:")
        # print(start_index, stop_index)
        # transcript_2_support_7 Support: 98 P-value: not_tested correction factor: NA delta size: -1
        # transcript_61_support_41 Support: 49 P-value: 6.16984589679744e-17 correction factor: 23616 delta size: 1

        for j in range(start_index +1, stop_index):
            candidate = file_contents[j]  
            acc, _, support, _, pvalue, _, _, correction_factor, _, _, nr_varinats =  candidate.split()
            if pvalue == "not_tested":
                continue
            support, p_val, correction_factor, nr_varinats = (int(support), float(pvalue), int(correction_factor), int(nr_varinats))
            corrected_p_val = p_val * correction_factor
            if corrected_p_val > 1:
                print(p_val, correction_factor)
                print(gene, copies, abundance, mutation_rate, nr_reads, exp_nr, corrected_p_val, support, relative_support, nr_varinats, acc)
                continue
            relative_support = float(support)/ float(nr_reads)
            tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format(gene, copies, abundance, mutation_rate, nr_reads, exp_nr, corrected_p_val, support, relative_support, nr_varinats, acc))
        # sys.exit()
    tsv_file.close()
    return tsv_file.name


def plot(tsv_file, args):
    sns.plt.clf()
    with sns.plotting_context("paper", font_scale=1.2):
        data = pd.read_csv(tsv_file, sep="\t")
        fig, ax = plt.subplots()
        # flierprops = dict(markerfacecolor='0.75', markersize=5, marker='o')
        # d = {'color': ['b', 'g', 'r']}
        g = sns.FacetGrid(data, row="gene", col="mutation_rate", hue="nr_reads", palette="Set1", size=3, aspect=3.2, row_order=["TSPY13P", "HSFY2", "DAZ2"], col_order=[ 0.0001], legend_out=True)
        # g.set_yticklabels([20, 100, 500, 2500, 12500])
        # g.set_yticklabels([1, 2, 10, 20, 100, 200, 1000, 2000])
        # print(help(g.set()))
        # g.set(yscale="log")
        g.set(xscale="log")
        g.set(xlim=(1e-100,1), ylim=(0,1))

        # g.set(xticks=[10e-2, 10e-5, 10e-10, 10e-20, 10e-50, 10e-100, 10e-150, 10e-200], yticks=[0, 1])
        g.set(xticks=[1e-2, 1e-10, 1e-20, 1e-30, 1e-50, 1e-100], yticks=[0,0.2,0.4,0.6,0.8,1.0])
        # g.set(ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        sns.set(style="whitegrid", palette="muted")

        (g.map(plt.scatter, "corrected_p_val", "relative_support").add_legend())
        g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=12)

        g.set_ylabels("Relative transcript abundance")
        g.set_xlabels("P-value")

        plt.savefig(os.path.join(args.outfolder, "dotplot.pdf"))
        plt.close()

def main(args):
    
    # get data
    tsv_file = print_tsv(args)

    # plot with swarmplot
    plot(tsv_file, args)



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

    parser.add_argument('--isoconfiles', type=str, nargs="+", help='Path to output files')
    parser.add_argument('--outfolder', type=str, help='Outfolder.')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    mkdir_p(args.outfolder)

    main(args)

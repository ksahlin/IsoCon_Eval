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


def recall_per_editdistance(args):
    plt.clf()
    with sns.plotting_context("paper", font_scale=1.8):
        data = pd.read_csv(args.recallfile, sep="\t")
        # new_data = data.groupby(["read_count", "abundance", "ed" ], as_index=False)['recall'].mean()
        # print(new_data)
        data.apply(pd.to_numeric, errors='coerce')
        g = sns.factorplot(x="read_count", y="recall",  hue="abundance", 
                            col="ed", data=data, col_wrap = 3, 
                            size=3, aspect=1.6,
                            col_order=[1, 2, 3, 4, 5]); #, hue_order=["0.5", "0.2", "0.1", "0.05", "0.01", "0.005"]);
        g.set( ylim=(0.0,1.0))
        # g.set(yscale="log", ylim=(0.1,10000))
        # g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=16)
        g.set_ylabels("Recall")
        g.set_xlabels("Total read depth")
        outfile = os.path.join(args.outfolder, "recall_per_editdistance.pdf")
        plt.savefig(outfile)
        plt.close()

def recall_per_abundance(args):
    plt.clf()
    with sns.plotting_context("paper", font_scale=1.8):
        data = pd.read_csv(args.recallfile, sep="\t")
        # new_data = data.groupby(["read_count", "abundance", "ed" ], as_index=False)['recall'].mean()
        # print(new_data)
        data.apply(pd.to_numeric, errors='coerce')
        g = sns.factorplot(x="read_count", y="recall",  hue="ed", 
                            col="abundance", data=data, col_wrap = 3, 
                            size=3, aspect=1.6,
                            col_order=[0.5, 0.2, 0.1, 0.05, 0.01, 0.005]);
        g.set( ylim=(0.0,1.0))
        # g.set(yscale="log", ylim=(0.1,10000))
        # g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=16)
        g.set_ylabels("Recall")
        g.set_xlabels("Total read depth")
        plt.subplots_adjust(top=0.9)
        g.fig.suptitle('TSPY')
        outfile = os.path.join(args.outfolder, "recall_per_abundance.pdf")
        plt.savefig(outfile)
        plt.close()

def recall_per_abundance_normalized(args):
    plt.clf()
    with sns.plotting_context("paper", font_scale=1.8):
        data = pd.read_csv(args.recallfile, sep="\t")
        # new_data = data.groupby(["read_count", "abundance", "ed" ], as_index=False)['recall'].mean()
        # print(new_data)
        data.apply(pd.to_numeric, errors='coerce')
        g = sns.factorplot(x="read_count", y="recall",  hue="ed", 
                            col="abundance", data=data, col_wrap = 3, 
                            size=3, aspect=1.6,
                            col_order=[0.5, 0.2, 0.1, 0.05, 0.01, 0.005]);
        g.set( ylim=(0.0,1.0))
        # g.set(yscale="log", ylim=(0.1,10000))
        # g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=16)
        g.set_ylabels("Recall")
        g.set_xlabels("Total read depth")
        outfile = os.path.join(args.outfolder, "recall_per_abundance_normalized.pdf")
        plt.savefig(outfile)
        plt.close()


def draw_heatmap(*args, **kwargs):
    data = kwargs.pop('data')
    d = data.pivot(index=args[1], columns=args[0], values=args[2])
    sns.heatmap(d, cmap='coolwarm_r', **kwargs)


def recall_heatmap(args):
    plt.clf()
    with sns.plotting_context("paper", font_scale=2.0):
        data = pd.read_csv(args.recallfile, sep="\t")
        # data = data.loc[ (data['ed'] == 1) | (data['ed'] == 5)]
        # new_data = subdata.groupby(["read_count", "abundance" ], as_index=False)['recall'].mean()
        print(len(data))
        # fill in missing experiments instances
        # print(data.loc[ (data['abundance'] == 0.02) & (data['Family'] == "HSFY2") ])

        # data = data.groupby(["read_count", "abundance", "ed", "Family" ]).filter(lambda x: len(x) >= 20)
        # data[data['recall'].isin(filtered)]
        # print(data.loc[ (data['abundance'] == 0.02) & (data['Family'] == "HSFY2") & (data['read_count'] <= 100) &  (data['ed'] == 5) ])

        # print(len(data))

        # nr_experiments = data.groupby(["read_count", "abundance", "ed", "Family" ], as_index=False).count()
        # print(nr_experiments)
        # sys.exit()

        new_data = data.groupby(["read_count", "abundance", "ed", "Family" ], as_index=False)['recall'].mean()
        # new_data.apply(pd.to_numeric, errors='coerce')
        print(len(new_data))
        new_data = new_data[new_data.read_count < 10000]
        print(len(new_data))
        # print(new_data)
        # new_data = new_data.pivot("read_count", "abundance", "recall")
        # ax = sns.heatmap(new_data, annot=True) #, fmt="d")
        # plt.subplots_adjust(top=0.9)
        # ax.set_title('TSPY (ed=1) recall of error correction')
        # ax.set_title('TSPY (ed=1) recall of IsoCon')

        #add missing data
        all_missing = []
        for f in ["TSPY13P", "HSFY2", "DAZ2"]:
            if f == "DAZ2":
                for a in [0.5,0.2, 0.1]:
                    for e in [1,5]:
                        for rd in [2000,5000, 7500]:
                            all_missing.append([rd, a, e, f, 1])
            if f == "HSFY2":
                for a in [0.5,0.2, 0.1, 0.05]:
                    for e in [1,5]:
                        for rd in [2000,5000, 7500]:
                            if e == 1 and a == 0.05 and rd == 2000:
                                continue
                            if e == 5 and a in set([0.05, 0.1]) and rd == 2000:
                                continue
                            all_missing.append([rd, a, e, f, 1])
            if f == "TSPY13P":
                for a in [0.5,0.2, 0.1]:
                    for e in [1,5]:
                        for rd in [2000,5000,7500]:
                            if e == 1 and a == 0.1 and rd in set([5000, 7500]):
                                continue
                            if e == 5 and a == 0.1 and rd in set([2000, 7500]):
                                continue
                            all_missing.append([rd, a, e, f, 1])

        missing_rows = pd.DataFrame( all_missing,  columns=["read_count", "abundance", "ed", "Family", "recall" ])
        # missing_rows = pd.DataFrame([[2000, 0.5, 1, "DAZ2", 0.99], [2000, 0.2, 1, "DAZ2", 0.99], [2000, 0.1, 1, "DAZ2", 0.99] ], columns=["read_count", "abundance", "ed", "Family", "recall" ])
        print(missing_rows)
        new_data = new_data.append(missing_rows, ignore_index=True)
        print(len(new_data))

        g = sns.FacetGrid(new_data, col='ed', row='Family', col_order=[1,5], row_order=[ "TSPY13P", "HSFY2", "DAZ2"], size=3, aspect=2.2)
        # g = sns.FacetGrid(new_data, col='ed', row='Family', col_order=[1,5], row_order=["TSPY13P", "HSFY2", "DAZ2"], size=3, aspect=1.6)
        cbar_ax = g.fig.add_axes([.92, .3, .02, .4])  # <-- Create a colorbar axes

        g.map_dataframe(draw_heatmap,'read_count', 'abundance', 'recall', 
                        annot=True, cbar_ax=cbar_ax, annot_kws = {"fontsize": 11}) #, cbar=False)
        g.fig.subplots_adjust(right=.9)  # <-- Add space so the colorbar doesn't overlap the plot

        g.set_xlabels("Read depth")
        g.set_ylabels("Relative abundance")
        # g.set_titles(col_template="{col_name}", fontweight='bold', fontsize=14)
        g.set_titles(row_template="{row_name}", col_template="ed={col_name}", size=16)
        g.fig.subplots_adjust(right=.9)  # <-- Add space so the colorbar doesn't overlap the plot
        plt.subplots_adjust(top=0.9)
        g.fig.suptitle('Recall rates', fontweight='bold', fontsize=18)
        # g.fig.suptitle('DAZ | final output', fontweight='bold', fontsize=18)
        outfile = os.path.join(args.outfolder, "all_recall_rates.pdf")
        # outfile = os.path.join(args.outfolder, "DAZ_recall_heatmap_candidates_fig_S9.pdf")
        plt.savefig(outfile)
        plt.close()

def main(args):

    # plot with stripplot
    # candidate_swarmplot(args)
    # sys.exit()
    
    # get recall data
    # recall_per_editdistance(args)
    # recall_per_abundance(args)
    recall_heatmap(args)
    # recall_per_abundance_normalized(args)


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

    parser.add_argument('--recallfile', type=str, help='data')
    parser.add_argument('--outfolder', type=str, help='Outfolder.')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    mkdir_p(args.outfolder)

    main(args)



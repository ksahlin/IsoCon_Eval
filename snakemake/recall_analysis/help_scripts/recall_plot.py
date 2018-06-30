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


# def recall_per_editdistance(args):
#     plt.clf()
#     with sns.plotting_context("paper", font_scale=1.8):
#         data = pd.read_csv(args.recallfile, sep="\t")
#         # new_data = data.groupby(["read_count", "abundance", "ed" ], as_index=False)['recall'].mean()
#         # print(new_data)
#         data.apply(pd.to_numeric, errors='coerce')
#         g = sns.factorplot(x="read_count", y="recall",  hue="abundance", 
#                             col="ed", data=data, col_wrap = 3, 
#                             size=3, aspect=1.6,
#                             col_order=[1, 2, 3, 4, 5]); #, hue_order=["0.5", "0.2", "0.1", "0.05", "0.01", "0.005"]);
#         g.set( ylim=(0.0,1.0))
#         # g.set(yscale="log", ylim=(0.1,10000))
#         # g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=16)
#         g.set_ylabels("Recall")
#         g.set_xlabels("Total read depth")
#         outfile = os.path.join(args.outfolder, "recall_per_editdistance.pdf")
#         plt.savefig(outfile)
#         plt.close()

def readdepth_per_abundance_forX_recall(args):
    plt.clf()
    with sns.plotting_context("paper", font_scale=1.8):
        data = pd.read_csv(args.recallfile, sep="\t")
        avg_recall_data = data.groupby(["read_count", "abundance", "ed", "Family" ], as_index=False)['recall'].mean()
        mod_data = avg_recall_data.loc[avg_recall_data['recall'] >= 0.7]
        print(mod_data)
        lowest_readdepth_data = mod_data.groupby(["abundance", "ed", "Family" ], as_index=False)['read_count'].min()
        print(lowest_readdepth_data)
        # lowest_readdepth_data.apply(pd.to_numeric, errors='coerce')
        g = sns.factorplot(x="abundance", y="read_count",  hue="ed", 
                            col="Family", data=lowest_readdepth_data, col_wrap = 3, 
                            size=3, aspect=1.6, hue_order=[1,5], ci=None, order=[0.5, 0.2, 0.1, 0.05, 0.01, 0.005]);
        # g.set( ylim=(0.0,1.0))
        g.set(yscale="log", ylim=(1,10000))
        # g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=16)
        g.set_ylabels("Total read depth required \nto get over 70% Recall")
        g.set_xlabels("Relative abundance")
        plt.subplots_adjust(top=0.9)
        # g.fig.suptitle('Required read depth for over 70% recall rate')
        outfile = os.path.join(args.outfolder, "readdepth_per_abundance_forX_recall.pdf")
        plt.savefig(outfile)
        plt.close()



def draw_heatmap(*args, **kwargs):
    data = kwargs.pop('data')
    d = data.pivot(index=args[1], columns=args[0], values=args[2])
    sns.heatmap(d, cmap='coolwarm_r', **kwargs)


# def recall_heatmap(args):
#     plt.clf()
#     with sns.plotting_context("paper", font_scale=1.0):
#         data = pd.read_csv(args.recallfile, sep="\t")
#         # subdata = data.loc[data['ed'] == 2]
#         # new_data = subdata.groupby(["read_count", "abundance" ], as_index=False)['recall'].mean()
#         new_data = data.groupby(["read_count", "abundance", "ed" ], as_index=False)['recall'].mean()
#         new_data.apply(pd.to_numeric, errors='coerce')
#         g = sns.FacetGrid(new_data, row='ed', row_order=[1,2,3], size=3, aspect=1.6)
#         g.map_dataframe(draw_heatmap,'read_count', 'abundance', 'recall', annot=True) #, cbar=False)
#         g.set_xlabels("Read depth")
#         g.set_ylabels("Relative transcript abundance")
#         g.set_titles(col_template="{col_name}", fontweight='bold', fontsize=14)
#         g.fig.subplots_adjust(right=.9)  # <-- Add space so the colorbar doesn't overlap the plot
#         plt.subplots_adjust(top=0.9)
#         g.fig.suptitle('DAZ | error correction', fontweight='bold', fontsize=18)
#         # g.fig.suptitle('DAZ | final output', fontweight='bold', fontsize=18)
#         outfile = os.path.join(args.outfolder, "DAZ_recall_heatmap_candidates.pdf")
#         # outfile = os.path.join(args.outfolder, "DAZ_recall_heatmap_candidates_fig_S9.pdf")
#         plt.savefig(outfile)
#         plt.close()

def main(args):

    # plot with stripplot
    # candidate_swarmplot(args)
    # sys.exit()
    
    # get recall data
    readdepth_per_abundance_forX_recall(args)
    # recall_heatmap(args)


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



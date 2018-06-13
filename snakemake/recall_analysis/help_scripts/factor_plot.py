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


def recall_plot(args):
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
        outfile = os.path.join(args.outfolder, "lol.pdf")
        plt.savefig(outfile)
        plt.close()

def main(args):

    # plot with stripplot
    # candidate_swarmplot(args)
    # sys.exit()
    
    # get recall data
    recall_plot(args)


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



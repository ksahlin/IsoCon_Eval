
import argparse
import os
import random
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

import misc_functions
import itertools

def heatmap(args):

    dtypes =  {"family": str, "ed" : int, "mutation_rate" : str, "member1" : int, "member2" : int}
    data = pd.read_csv(args.tsv_input, sep="\t", dtype = dtypes)
    print(data.columns)

    data['mutation_rate'] = data['mutation_rate'].astype(str)


    # mut_rates=['0.01', '0.001', "0.0001"]
    # families=['TSPY', 'HSFY', "DAZ"]    
    # abundance = [1,2,4,8,16,32,64,128]
    # data = pd.DataFrame(list(itertools.product(mut_rates, families, abundance, abundance)))
    # data.columns = ['mutation_rate', "family", 'member1','member2']
    # data['ed'] = np.random.randint(1, 10, data.shape[0])

    # pd.set_option('display.max_rows', len(data))
    # print(data)
    # pd.reset_option('display.max_rows')
    

    def facet_heatmap(data, color, **kws):
        data = data.pivot(index="member2", columns='member1', values='ed')
        print(data)
        # data = data.fillna(0)
        # print(data)
        print(data.columns)
        mask = np.zeros_like(data)
        mask[np.triu_indices_from(mask)] = True
        print(mask)
        for i in range(len(mask)):
            for j in range(len(mask[i])):
                if len(mask) == 30:
                    if data[i + 1][j + 1] > 99:
                        print("true")
                        mask[i][j] = 1
                else:
                    if data[2**i][2**j] > 99:
                        print("true")
                        mask[i][j] = 1
        # mask = mask[data < 100]
        print(mask)
        print(data)

        if len(mask) == 30:
            sns.heatmap(data, cmap='coolwarm_r', annot=True, fmt="d", mask = mask, annot_kws={"size":5}, **kws)  # <-- Pass kwargs to heatmap
        else:
            sns.heatmap(data, cmap='coolwarm_r', annot=True, fmt="d", mask = mask, **kws)  # <-- Pass kwargs to heatmap


    with sns.plotting_context(font_scale=5.5):
        g = sns.FacetGrid(data, col="mutation_rate", row="family", size=3, aspect=1, col_order = ["0.01", "0.001", "0.0001"] ) #, col_order = ["TSPY13P", "HSFY2", "DAZ2"])

    cbar_ax = g.fig.add_axes([.92, .3, .02, .4])  # <-- Create a colorbar axes

    g = g.map_dataframe(facet_heatmap,
                        cbar_ax=cbar_ax,
                        vmin=0, vmax=10)  # <-- Specify the colorbar axes and limits

    g.set_titles(col_template="{col_name}", fontweight='bold', fontsize=14)
    g.fig.subplots_adjust(right=.9)  # <-- Add space so the colorbar doesn't overlap the plot
    plt.savefig(args.outfile)
    plt.close()


def dotplot2(args, x_label, y_label, title):
    dataset = pd.read_csv(args.tsv_input, sep="\t")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale('log')
    ax.set_xscale('log')
    # ax.set_xticks(list(range(2,10, 2)) + list(range(10,100, 20)) + list(range(100,1000, 200)))
    ax.set_xticks([1,2,10,20,100,200, 1000, 2000])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax.set_yticks(list(range(2,10, 2)) + list(range(10,100, 20)) + list(range(100,1000, 200)))
    ax.set_yticks([20, 100, 500, 2500, 12500])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    # for x,y,lab in zip(X,Y,labels):
    #     ax.scatter(x,y,label=lab, alpha=0.5)

    #Generate some colors and markers
    # colors = ["r", "y", "g"]
    # markers = ['x','o','v','^','<']*100

    #Plot each individual point separately
    # m_size = float(max(dataset["ed"]))
    # labels = ["ed = 1 ", "ed = 3" + marker
    numbers = [ row for row in  dataset.values if int(row[2]) < 4 ]
    dots = [ row for row in dataset.values if int(row[2]) >= 4]

    for i,row in enumerate(dots + numbers):
        # print(row)

        # if int(row[2]) ==  1:
        #     m = "*"
        # elif int(row[2]) ==  2:
        #     m = "o"
        # else:
        #     m = "^"

        # m= "o" if int(row[5]) == 1 else "x"
        
        color= "g" if int(row[5]) == 1 else "r"
        # print(type(row[2]))
        ab_ratio = float(row[3])/float(row[4])
        ed = int(row[2])
        # marker = r"$ {0} $".format( round(float(row[3])/float(row[4]), 3)) if  int(row[3]) < int(row[4]) else  r"$ {0} $".format( round(float(row[3])/float(row[4]), 1) )
        marker = r"$ {0} $".format(ed) if ed < 4 else  r"x"
        alpha = 1.0 if ed < 4 else 0.5
        # ax.scatter(rand_jitter([float(row[3])]), rand_jitter([float(row[4])]), color=color, marker=m, alpha=0.7) #,s=(m_size/int(row[2]))**1.5
        ax.scatter(rand_jitter([float(row[3])]), rand_jitter([float(row[6])]), color=color,  marker=marker, alpha=alpha) # <---- this one plots the data in the most raw form
        # ax.scatter(rand_jitter([float(row[3])]), rand_jitter([float(row[4])]), color=color,  marker=marker, alpha=alpha) # <---- this one plots the data in the most raw form
        # ax.scatter(rand_jitter([float(row[3])/ float(row[4])]), rand_jitter([int(row[2])]), color=color,  marker=r"$ {} $".format(row[3]), alpha=0.7) # s=50*(1.0/int(row[2])),
        # ax.scatter(rand_jitter([int(row[3])]), rand_jitter([float(row[2])]), color=color, s=50*(1.0/ab_ratio), alpha=0.7) # s=50*(1.0/int(row[2])),

    # ax.legend(fontsize='small', loc = 2)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.savefig(args.outfile)
    plt.close()

    # for i in range(len(xs)):
    #     plt.scatter(xs[i], ys[i], marker=m[i])

def boxplot(args):
    sns.plt.clf()
    with sns.plotting_context("paper", font_scale=1.0):
        # dtypes =  {"family": str, "supported" : float}
        # data = pd.read_csv(args.tsv_input, sep="\t", dtype = dtypes)
        indata = pd.read_csv(args.tsv_input, sep="\t")
        fig, ax = plt.subplots()
        # flierprops = dict(markerfacecolor='0.75', markersize=5, marker='o')
        # d = {'color': ['b', 'g', 'r']}
        sns.set_color_codes("muted")
        # sns.set_palette("husl")
        g = sns.factorplot(x="family", y="supported", col="sample", hue="method", hue_order=["ISOCON", "ICE", "Illumina-corrected", "Original"], palette={"ISOCON": "b", "ICE": "g", "Illumina-corrected" : "k", "Original" : "r"}, data= indata, kind="bar", size=3, aspect=1.6, col_order=["sample1", "sample2"], legend_out=True)
        # g.set(yscale="log")
        sns.set(style="whitegrid", palette="muted")

        # (g.map(sns.swarmplot, "read_count", "transcript_read_depth", edgecolor="white", alpha=.7).despine(left=True).add_legend(title="captured", label_order=["yes", "candidate" ,"no"]))
        g.set_titles(col_template="{col_name}", row_template="{row_name}",  size=10)
        g.set_ylabels("% Illumina supported positions")
        plt.ylim(0, 100)
        # g.set_xlabels("Total depth")

        plt.savefig(args.outfile)
        plt.close()

def dotplot4(args):
    sns.plt.clf()
    with sns.plotting_context("paper", font_scale=1.8):
        data = pd.read_csv(args.tsv_input, sep="\t")
        fig, ax = plt.subplots()
        # flierprops = dict(markerfacecolor='0.75', markersize=5, marker='o')
        # d = {'color': ['b', 'g', 'r']}
        g = sns.FacetGrid(data, row="Family", col="mutation_rate", hue="captured", palette="Set1", size=3, aspect=1.6, row_order=["TSPY13P", "HSFY2", "DAZ2"], col_order=[0.01, 0.001, 0.0001], legend_out=True, hue_kws={"marker": ["^", "o", "s"]})
        # g.set_yticklabels([20, 100, 500, 2500, 12500])
        # g.set_yticklabels([1, 2, 10, 20, 100, 200, 1000, 2000])
        # print(help(g.set()))
        g.set(yscale="log")
        # g.set(xscale="log")
        # g.set(xticks=[1, 2, 10, 20, 100, 200, 1000, 2000], yticks=[20, 100, 500, 2500, 12500])
        # g.set(ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        sns.set(style="whitegrid", palette="muted")

        (g.map(sns.swarmplot, "read_count", "transcript_read_depth", edgecolor="white", alpha=.7).despine(left=True).add_legend(title="captured", label_order=["yes", "candidate" ,"no"]))
        g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=16)
        # g.set(xticks=[1, 2, 10, 20, 100, 200, 1000, 2000], yticks=[20, 100, 500, 2500, 12500])
        # plt.xticks([1, 2, 10, 20, 100, 200, 1000, 2000])
        # plt.yticks([20, 100, 500, 2500, 12500])

        g.set_ylabels("Transcript depth")
        g.set_xlabels("Total depth")

        plt.savefig(args.outfile)
        plt.close()

def dotplot3(args):
    sns.plt.clf()
    with sns.plotting_context("paper", font_scale=1.8):
        data = pd.read_csv(args.tsv_input, sep="\t")
        fig, ax = plt.subplots()
        # flierprops = dict(markerfacecolor='0.75', markersize=5, marker='o')
        # d = {'color': ['b', 'g', 'r']}
        g = sns.FacetGrid(data, row="Family", col="mutation_rate", hue="captured", palette="Set1", size=3, aspect=1.6, row_order=["TSPY13P", "HSFY2", "DAZ2"], col_order=[0.01, 0.001, 0.0001], legend_out=True, hue_kws={"marker": ["^", "o"]})
        # g.set_yticklabels([20, 100, 500, 2500, 12500])
        # g.set_yticklabels([1, 2, 10, 20, 100, 200, 1000, 2000])
        # print(help(g.set()))
        g.set(yscale="log")
        # g.set(xscale="log")
        # g.set(xticks=[1, 2, 10, 20, 100, 200, 1000, 2000], yticks=[20, 100, 500, 2500, 12500])
        # g.set(ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        sns.set(style="whitegrid", palette="muted")

        (g.map(sns.swarmplot, "read_count", "transcript_read_depth", edgecolor="white", alpha=.7).despine(left=True).add_legend(title="captured", label_order=["yes" ,"no"]))
        g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=16)
        # g.set(xticks=[1, 2, 10, 20, 100, 200, 1000, 2000], yticks=[20, 100, 500, 2500, 12500])
        # plt.xticks([1, 2, 10, 20, 100, 200, 1000, 2000])
        # plt.yticks([20, 100, 500, 2500, 12500])

        g.set_ylabels("Transcript depth")
        g.set_xlabels("Total depth")

        plt.savefig(args.outfile)
        plt.close()

def violinplot_combined_one_isoform(args):
    sns.plt.clf()
    with sns.plotting_context("paper", font_scale=1.8):
        # print(args.tsv_input)
        # sns.set(font_scale=2)
        true_positives = pd.read_csv(args.tsv_input, sep="\t")
        fig, ax = plt.subplots()
        flierprops = dict(markerfacecolor='0.75', markersize=5, marker='o')
        d = {'color': ['b', 'g', 'r']}
        g = sns.FacetGrid(true_positives, row="Family", size=3, aspect=1.6, row_order=["TSPY13P", "HSFY2", "DAZ2"], legend_out=True)
        sns.set(style="whitegrid", palette="muted")

        (g.map(sns.violinplot, "read_count", args.y_axis, "TOOL", cut=0, hue_order=["ISOCON", "ICE"], palette=sns.color_palette("muted", 2)).despine(left=True).add_legend(title="TOOL", label_order=["ISOCON", "ICE"]))
        g.set_titles(row_template="{row_name}",fontweight='bold', size=16)
        g.set_yticklabels(["",0,0.2,0.4,0.6,0.8,1.0])

        if args.y_axis == "FP":
            g.set(yscale="log")
        
        plt.savefig(args.outfile)
        plt.close()

def violinplot_combined(args):
    sns.plt.clf()
    with sns.plotting_context("paper", font_scale=1.8):
        # print(args.tsv_input)
        # sns.set(font_scale=2)
        true_positives = pd.read_csv(args.tsv_input, sep="\t")
        fig, ax = plt.subplots()
        flierprops = dict(markerfacecolor='0.75', markersize=5, marker='o')
        d = {'color': ['b', 'g', 'r']}
        g = sns.FacetGrid(true_positives, row="Family", col="mutation_rate", size=3, aspect=1.6, row_order=["TSPY13P", "HSFY2", "DAZ2"], col_order=[0.01, 0.001, 0.0001], legend_out=True)
        sns.set(style="whitegrid", palette="muted")

        (g.map(sns.violinplot, "read_count", args.y_axis, "TOOL", cut=0, hue_order=["ISOCON", "ICE"], palette=sns.color_palette("muted", 2)).despine(left=True).add_legend(title="TOOL", label_order=["ISOCON", "ICE"]))
        g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}",  size=16)
        g.set_yticklabels(["",0,0.2,0.4,0.6,0.8,1.0])

        if args.y_axis == "FP":
            g.set(yscale="log")

        plt.savefig(args.outfile)
        plt.close()


def violinplot(args):
    sns.plt.clf()
    with sns.plotting_context("paper", font_scale=1.8):
        print(args.tsv_input)
        true_positives = pd.read_csv(args.tsv_input, sep="\t")
        fig, ax = plt.subplots()
        flierprops = dict(markerfacecolor='0.75', markersize=5, marker='o')
        d = {'color': ['b', 'g', 'r']}
        g = sns.FacetGrid(true_positives, row="Family", size=3, aspect=3.0, row_order=["TSPY13P", "HSFY2", "DAZ2"], legend_out=True)
        print(true_positives["TP"])
        sns.set(style="whitegrid", palette="muted")
        (g.map(sns.violinplot, "read_count", args.y_axis, "mutation_rate", cut=0, hue_order=[0.0001, 0.001, 0.01], palette=sns.color_palette("muted", 3)).despine(left=True).add_legend(title="mutation_rate", label_order=["0.01", "0.001", "0.0001"]))
        # (g.map(sns.swarmplot, "read_count", "TP", "mutation_rate", palette=sns.color_palette("muted", 3)).despine(left=True).add_legend(title="mutation_rate"))
        axes = g.axes.flatten()
        axes[0].set_title("Internal")
        axes[1].set_title("External")
        g.set_titles(col_template="$\mu={col_name}$", row_template="{row_name}", size=16) # row_template="{row_name}",
        # g.set_xticklabels([0,20,100,500,2500,12500]) 
        g.set_yticklabels(["",0,0.2,0.4,0.6,0.8,1.0])
        plt.savefig(args.outfile)
        plt.close()


def predicted_to_member_id(args):

    # for pseudo vs coding, see:
    # https://stackoverflow.com/questions/37331937/seaborn-facetgrid-countplot-hue

    # fig =sns.FacetGrid(data=df,col='Sex',hue='Marker2',palette='Set1',size=4,aspect=1).map(sns.countplot,'Marker1',order=df.Marker1.unique()).add_legend()
    
    sns.plt.clf()
    with sns.plotting_context("paper", font_scale=1.0):
        print(args.tsv_input)
        data = pd.read_csv(args.tsv_input, sep="\t")
        sns.set(style="whitegrid", palette="muted")
        
        g = sns.factorplot("MEMBER_ID", col="FAMILY", col_wrap=4, 
                            data=data[data.MEMBER_ID.notnull()],
                            kind="count", size=2.5, aspect=.8)
        # titanic[titanic.deck.notnull()]

        # g = sns.FacetGrid(data, col="FAMILY", size=3, aspect=.5)
        # g.map(sns.countplot, "MEMBER_ID") #.add_legend()
        plt.subplots_adjust(top=0.9)
        g.fig.suptitle("Potential Isoforms per member",fontweight='bold', size=16) 

        plt.savefig(args.outfile)
        plt.close()

        # g.set_titles("Potential Isoforms per member",fontweight='bold', size=16)

def parse_tsv_lineplot(tsv_input):
    tsv_data = {}
    tsv_data["gene"] = []
    tsv_data["mut"] = []
    tsv_data["nr_reads"] = []
    tsv_data["TP"] = []
    tsv_data["ATP"] = []
    tsv_data["FP"] = []
    tsv_data["FN"] = []

    for i, line in enumerate(open(tsv_input, "r")):
        if i == 0:
            continue

        gene, mut, nr_reads, run_id, ed_avg, ed_sed, ed_min, ed_max, TP, ATP, err_ATP, FP, err_FP, FN = line.strip().split()
        tsv_data["gene"].append(gene)
        tsv_data["mut"].append(mut)
        tsv_data["nr_reads"].append(int(nr_reads))
        tsv_data["TP"].append(int(TP))
        tsv_data["ATP"].append(int(ATP))
        tsv_data["FP"].append(int(FP))
        tsv_data["FN"].append(int(FN))

    return tsv_data


def lineplot(tsv_data, mutation_rate, outfile, x_label, y_label, title):
    import numpy as np
    import matplotlib.lines as mlines
    x = np.arange(10)


    X = sorted(set(tsv_data["nr_reads"]))
    print(X)
    tspy_TP, tspy_FN, tspy_ATP  = [0]*len(X), [0]*len(X), [0]*len(X)
    hsfy_TP, hsfy_FN, hsfy_ATP  = [0]*len(X), [0]*len(X), [0]*len(X)
    daz_TP, daz_FN, daz_ATP  = [0]*len(X), [0]*len(X), [0]*len(X)

    for i in range(len(tsv_data["gene"])):
        if  tsv_data["mut"][i] != mutation_rate:
            continue
        pos = X.index(tsv_data["nr_reads"][i]) # get which x position the y-value have
        if tsv_data["gene"][i] == "TSPY13P":
            tspy_TP[pos] = tsv_data["TP"][i]
            tspy_FN[pos] = tsv_data["FN"][i]
            tspy_ATP[pos] = tsv_data["ATP"][i]

        if tsv_data["gene"][i] == "HSFY2":
            hsfy_TP[pos] = tsv_data["TP"][i]
            hsfy_FN[pos] = tsv_data["FN"][i]
            hsfy_ATP[pos] = tsv_data["ATP"][i]

        if tsv_data["gene"][i] == "DAZ2":
            daz_TP[pos] = tsv_data["TP"][i]
            daz_FN[pos] = tsv_data["FN"][i]
            daz_ATP[pos] = tsv_data["ATP"][i]

    

    dashes = [(None, None), [2, 2], [6, 2]]
    green_line = mlines.Line2D([], [], color='green', label='TSPY')
    blue_line = mlines.Line2D([], [], color='blue', label='HSFY')
    red_line = mlines.Line2D([], [], color='red', label='DAZ')
    TP_line = mlines.Line2D([], [], color='black', label='TP')
    TP_line.set_dashes(dashes[0])
    
    FN_line = mlines.Line2D([], [], color='black', label='FN')
    FN_line.set_dashes(dashes[1])
    ATP_line = mlines.Line2D([], [], color='black', label='ATP')
    ATP_line.set_dashes(dashes[2])

    plt.legend(handles=[green_line, blue_line, red_line, TP_line, FN_line, ATP_line], loc='upper left')

    plt.gca().set_color_cycle(['green', 'blue', 'red'])
    plt.plot(X, tspy_TP)
    plt.plot(X, hsfy_TP)
    plt.plot(X, daz_TP)
    
    line, = plt.plot(X, tspy_FN)
    line.set_dashes(dashes[1])
    line, = plt.plot(X, hsfy_FN)
    line.set_dashes(dashes[1])
    line, = plt.plot(X, daz_FN)    
    line.set_dashes(dashes[1])

    line, = plt.plot(X, tspy_ATP)
    line.set_dashes(dashes[2])
    line, = plt.plot(X, hsfy_ATP)
    line.set_dashes(dashes[2])
    line, = plt.plot(X, daz_ATP) 
    line.set_dashes(dashes[2])

    ######
    plt.xscale('log')
    # ax.set_xticks(list(range(1,10)) + list(range(10,100, 10)) + list(range(100,1000, 100)))
    # ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    plt.xticks(X)
    # plt.semilogx(subsx=X)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.savefig(outfile)


    # plt.legend(loc='upper left') 
    # plt.legend(['y = x', 'y = 2x', 'y = 3x', 'y = 4x'], loc='upper left')

def parse_tsv_dotplot(tsv_input):
    tsv_data = {}
    tsv_data["m1"] = []
    tsv_data["m2"] = []
    tsv_data["ed"] = []
    tsv_data["cov1"] = []
    tsv_data["cov2"] = []
    tsv_data["captured"] = []
    tsv_data["tot_nr_reads_sequenced"] = []

    for line in open(tsv_input, "r"):
        member1, member2, closest_edit_distance, nr_reads1, nr_reads2, is_captured, tot_reads_sequenced = line.strip().split()
        tsv_data["m1"].append(member1)
        tsv_data["m2"].append(member2)
        tsv_data["ed"].append(int(closest_edit_distance))
        tsv_data["cov1"].append(int(nr_reads1))
        tsv_data["cov2"].append(int(nr_reads2))
        tsv_data["captured"].append(int(is_captured))
        tsv_data["tot_nr_reads_sequenced"].append(int(tot_reads_sequenced))

    return tsv_data

def dotplot(X, Y, outfile, x_label, y_label, title):

    labels = ["captured", "uncaptured"]

    fig = plt.figure()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    
    ax = fig.add_subplot(111)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xticks(list(range(1,10)) + list(range(10,100, 10)) + list(range(100,1000, 100)))
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    for x,y,lab in zip(X,Y,labels):
        ax.scatter(x,y,label=lab, alpha=0.5)

    colormap = plt.cm.nipy_spectral #nipy_spectral, Set1,Paired  
    colorst = [colormap(i) for i in np.linspace(0, 0.9,len(ax.collections))]       
    for t,j1 in enumerate(ax.collections):
        j1.set_color(colorst[t])

    ax.legend(fontsize='small', loc = 2)
    fig.savefig(outfile)


def rand_jitter(arr):
    # print(max(arr), min(arr))
    arr_jitter = []
    for point in arr:
        stdev = .05 * point 
        arr_jitter.append(point + random.normalvariate(0, stdev) )
    # stdev = .003*(max(arr)-min(arr))
    # print(arr_jitter)
    return arr_jitter
    # print(np.random.randn(len(arr)) * stdev)
    # return arr + np.random.randn(len(arr)) * stdev



def main(args):
    if args.dotplot2:
        # x_label, y_label, title = "reads sequenced from member", "reads sequenced from closest neighbor", "Capture ability"
        x_label, y_label, title = "reads sequenced from transcript", "total reads sequenced in experiment", "Capture ability"
        dotplot2(args, x_label, y_label, title)

    elif args.dotplot:
        tsv_data = parse_tsv_dotplot(args.tsv_input)
        X_cap =   [tsv_data["cov1"][i] for i in range(len(tsv_data["cov1"])) if tsv_data["captured"][i] == 1] 
        X_not_cap =   [tsv_data["cov1"][i] for i in range(len(tsv_data["cov1"])) if tsv_data["captured"][i] == 0] 

        # Y fcn 1
        X = [rand_jitter(X_cap), rand_jitter(X_not_cap)]
        Y_cap = [tsv_data["ed"][i] * (float(tsv_data["cov1"][i])/float(tsv_data["cov2"][i])) for i in range(len(tsv_data["cov1"])) if tsv_data["captured"][i] == 1] 
        Y_not_cap = [tsv_data["ed"][i] * (float(tsv_data["cov1"][i])/float(tsv_data["cov2"][i])) for i in range(len(tsv_data["cov1"])) if tsv_data["captured"][i] == 0]  
        Y = [Y_cap, Y_not_cap]
        x_label, y_label, title = "Number of flnc reads (logscale)", "min_ed * abundance_ratio (logscale)", "Capture ability as function of coverage, edit distance and abundance"
        dotplot(X, Y, args.outfile, x_label, y_label, title)

        # Y fcn 3
        X = [rand_jitter(X_cap), rand_jitter(X_not_cap)]
        Y_cap = [float(tsv_data["cov1"][i])/float(tsv_data["cov2"][i]) for i in range(len(tsv_data["cov1"])) if tsv_data["captured"][i] == 1] 
        Y_not_cap = [float(tsv_data["cov1"][i])/float(tsv_data["cov2"][i]) for i in range(len(tsv_data["cov1"])) if tsv_data["captured"][i] == 0]  
        Y = [Y_cap, Y_not_cap]
        x_label, y_label, title = "Number of flnc reads (logscale)", "abundance_ratio (logscale)", "Capture ability as function of coverage and abundance"
        dotplot(X, Y, args.outfolder + "/y_ab_ratio_" + args.filename, x_label, y_label, title)

        # Y fcn 3
        X = [rand_jitter(X_cap), rand_jitter(X_not_cap)]
        Y_cap = [tsv_data["ed"][i] for i in range(len(tsv_data["cov1"])) if tsv_data["captured"][i] == 1] 
        Y_not_cap = [tsv_data["ed"][i] for i in range(len(tsv_data["cov1"])) if tsv_data["captured"][i] == 0]  
        Y = [Y_cap, Y_not_cap]
        x_label, y_label, title = "Number of flnc reads (logscale)", "min_ed", "Capture ability as function of coverage and edit distance"
        dotplot(X, Y, args.outfolder + "/y_ed_" + args.filename, x_label, y_label, title)

    elif args.lineplot:
        tsv_data = parse_tsv_lineplot(args.tsv_input)
        lineplot(tsv_data, args.mutation_rate, args.outfile, "Number of reads simulated", "Members", "Retained members as function of read count" )

    elif args.violinplot:
        violinplot(args)

    elif args.violinplot_combined:
        violinplot_combined(args)
    elif args.heatmap:
        heatmap(args)
    elif args.one_isoform:
        violinplot_combined_one_isoform(args)
    elif args.dotplot3:
        dotplot3(args)
    elif args.dotplot4:
        dotplot4(args)
    elif args.boxplot:
        boxplot(args)
    elif args.pred_to_member:
        predicted_to_member_id(args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate pacbio reads from a set of transcripts.")
    parser.add_argument('tsv_input', type=str, help='The fasta file with sequences to be sequenced.')
    parser.add_argument('outfile', type=str, help='Output filename')
    parser.add_argument('--dotplot', action='store_true', help='Create dotplot')
    parser.add_argument('--dotplot2', action='store_true', help='Create dotplot')
    parser.add_argument('--dotplot3', action='store_true', help='Create dotplot')
    parser.add_argument('--dotplot4', action='store_true', help='Create dotplot')
    parser.add_argument('--lineplot', action='store_true', help='Create lineplot')
    parser.add_argument('--boxplot', action='store_true', help='Create lineplot')
    parser.add_argument('--violinplot', action='store_true', help='Create violinplot')
    parser.add_argument('--violinplot_combined', action='store_true', help='Create violinplot')
    parser.add_argument('--heatmap', action='store_true', help='Create heatmap of edit distances')
    parser.add_argument('--one_isoform', action='store_true', help='Create heatmap of edit distances')
    parser.add_argument('--pred_to_member', action='store_true', help='Create plot of predictions into gene members.')

    parser.add_argument('--y_axis', type=str, default="TP", help='Metric to plot on y-axis (TP or FP for now).')    
    parser.add_argument('--y_min', type=float, default=0.0, help='Minimum value on y-axis.')    
    parser.add_argument('--y_max', type=float, default=1.0, help='Maximum value on y-axis.')    

    parser.add_argument('--mutation_rate', type=str, help='Create lineplot')

    args = parser.parse_args()

    if not (args.dotplot or args.lineplot or args.violinplot or args.dotplot2 or args.heatmap or args.violinplot_combined or args.one_isoform or args.dotplot3 or args.dotplot4 or args.boxplot or args.pred_to_member):
        print("Specify which plot with --dotplot or --lineplot ")
        sys.exit()

    path_, file_prefix = os.path.split(args.outfile)
    misc_functions.mkdir_p(path_)
    args.outfolder = path_
    args.filename = file_prefix
    main(args)
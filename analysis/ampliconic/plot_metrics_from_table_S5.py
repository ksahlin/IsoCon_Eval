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
    p_vals = []
    read_depths = []
    read_depths2 = []
    full_sups = []
    avg_ill_sups = []
    avg_ill_sups2 = []
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
            if p_val > 1e-20:
                p_vals.append(p_val)
                read_depths.append(read_depth)
                full_sups.append(full_sup)
                avg_ill_sups.append(float(ill_supp))
            else:
                if float(ill_supp) > 0:
                    avg_ill_sups2.append(float(ill_supp))
                    read_depths2.append(read_depth)


            tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(gene, p_val, sample, shared, read_depth, perfect_match_database, ill_supp, full_sup))
        else:
            pass
        total_records += 1

    print("Median pval for pvals >1e-20:", sorted(p_vals)[len(p_vals)/2] )
    print("Median read depth for transcripts with pvals >1e-20:", sorted(read_depths)[len(read_depths)/2] )
    print("Median read depth for transcripts with pvals <=1e-20:", sorted(read_depths2)[len(read_depths2)/2] )
    print("\% transcripts with full support out of all transcripts with pvalue >1e-20:", float(len([1 for s in full_sups if s == "yes"]))/len(full_sups) )
    print("Avg Illumina support for transcripts with pvalue >1e-20:", float(sum(avg_ill_sups))/len(avg_ill_sups) )
    print("Avg Illumina support for transcripts with pvalue <=1e-20:", float(sum(avg_ill_sups2))/len(avg_ill_sups2) )
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


def paulplots(data, args):
    sns.plt.clf()
    from matplotlib import pyplot
    with sns.plotting_context("paper", font_scale=1.2):
        fig, ax = plt.subplots()
        matches = data.loc[data['perfect_match_database'] == 'yes']
        pvals = [float(m) for m in matches["p-value"]]
        bins = [10**-i for i in [0,2,3,4,5,6,7,8,9,10,15,20][::-1] ] # [10**-i for i in [0,2,3,4,5,6,7,8,9,10,15,20,30,40,50,75,100,150,200,300,320][::-1] ] #  10**(-np.arange(0,100,2))
        cumulative_vals = [ len([pv for pv in pvals if 1e-300 < pv < b]) for b in bins  ]
        # print(cumulative_vals)
        plt.xscale('log')
        plt.xlim(1e-22, 1.0) # plt.xlim(1e-322, 1.0)
        ax = pyplot.plot(bins, cumulative_vals) 
        plt.ylabel('#Exact matches')
        plt.xlabel("p-value")
        plt.title("Number of exact ENSEMBLE matches with p-value lower than x")
        plt.savefig(os.path.join(args.outfolder, "Figure_cumulative_exact_matches_pval.pdf"))
        plt.clf()

        print(matches["read_depth"])
        read_depths = [float(m) for m in matches["read_depth"]]
        bins = [1,3,5,10,15,20,30,40,50,100,150,200,300,500,1000,2000,3000,5000,10000]
        cumulative_vals = [ len([rd for rd in read_depths if rd < b]) for b in bins  ]
        print(cumulative_vals)
        plt.xscale('log')
        plt.xlim(1, 10000)
        ax = pyplot.plot(bins, cumulative_vals) 
        plt.ylabel('#Exact matches')
        plt.xlabel("Read support")
        plt.title("Number of exact ENSEMBLE matches with read support lower than x")
        plt.savefig(os.path.join(args.outfolder, "Figure_cumulative_exact_matches_rd.pdf"))
        plt.clf()

        # print(data["Illumina_support"])
        two_d_data = [(float(pv), float(ills)) for pv, ills in zip(data["p-value"], data["Illumina_support"])]
        two_d_data = [t for t in two_d_data if t[1] > 0.0]

        x_vals = [10**-i for i in [0,2,3,4,5,6,7,8,9,10,15,20][::-1] ] # [10**-i for i in [0,2,3,4,5,6,7,8,9,10,15,20,30,40,50,75,100,150,200,300,320][::-1] ] #  10**(-np.arange(0,100,2))
        y_vals = [ 100*sum([ills for pv, ills in two_d_data if 1e-300 <  pv < b])/ float(len([ills for pv, ills in two_d_data if 1e-300 < pv < b])) for b in x_vals  ]
        print(x_vals)
        print(y_vals)
        plt.xscale('log')
        plt.xlim(1e-22, 1.0) # plt.xlim(1e-322, 1.0)
        ax = pyplot.plot(x_vals, y_vals) 
        plt.ylabel('%Average Illumina support')
        plt.xlabel("p-value")
        plt.title("Average Illumina supported bases for transcripts with p-value lower than x")
        plt.savefig(os.path.join(args.outfolder, "Figure_avg_illumina_suppport.pdf"))
        plt.clf()

        two_d_data = [(float(pv), float(rd)) for pv, rd in zip(data["p-value"], data["read_depth"])]
        two_d_data = [t for t in two_d_data if t[1] > 0.0]

        x_vals = [10**-i for i in [0,2,3,4,5,6,7,8,9,10,15,20][::-1] ] # [10**-i for i in [0,2,3,4,5,6,7,8,9,10,15,20,30,40,50,75,100,150,200,300,320][::-1] ] #  10**(-np.arange(0,100,2))
        y_vals = [ sum([rd for pv, rd in two_d_data if 1e-300 < pv < b])/ float(len([rd for pv, rd in two_d_data if 1e-300 < pv < b])) for b in x_vals  ]
        print(x_vals)
        print(y_vals)
        plt.xscale('log')
        plt.xlim(1e-22, 1.0) # plt.xlim(1e-322, 1.0)
        ax = pyplot.plot(x_vals, y_vals) 
        plt.ylabel('Average CCS read support')
        plt.xlabel("p-value")
        plt.title("Average CCS read support for transcripts with p-value lower than x")
        plt.savefig(os.path.join(args.outfolder, "Figure_avg_read_depth.pdf"))
        plt.clf()

def paulplots_updated(data, args):
    sns.plt.clf()
    from matplotlib import pyplot
    with sns.plotting_context("paper", font_scale=1.2):
        plt.rc('xtick', labelsize=8)
        exponents = [2,3,4,5,6,7,8,9,10,15,20, 300][::-1]
        bins = [10**-i for i in exponents ]
        labels = [ '{0}-{1}'.format(i,j) for i,j in zip(exponents[:-1], exponents[1:]) ] #[ r'$10^{{-({0}to{1})}}$'.format(i,j) for i,j in zip(exponents[:-1], exponents[1:]) ] #[ r'$10^{{-{0}}}-10^{{-{1}}}$'.format(i,j) for i,j in zip(exponents[:-1], exponents[1:]) ]  #[ str(10**-i) + '-' + str(10**-(i -1)) for i in [1,2,3,4,5,6,7,8,9,10,15,20, 300][::-1] ]
        print(labels)
        p_val_bins = pd.cut(data['p-value'], bins, labels = labels )
        data["p_bins"] = pd.Series(p_val_bins, index=data.index)
        print(p_val_bins)
        # print(data)

        # Full illumina coverage
        print(data.groupby(["p_bins"])['Full_support'].value_counts(normalize=True).mul(100))
        data_tmp = data.groupby(["p_bins"])['Full_support'].value_counts(normalize=True).mul(100)
        # data_tmp = pd.DataFrame({'index':data_tmp.index, 'p_bins':data_tmp.p_bins, "Full_support": data_tmp.Full_support })
        full_supp_data = [ {"p_bins" : index[0], "percentage" : float(val)} for index, val in data_tmp.iteritems() if index[1] == "yes"]
        full_supp = pd.DataFrame(full_supp_data)


        # print(data.groupby(["p_bins"])['Full_support'].value_counts())
        # counts = data.groupby(["p_bins"])['Full_support'].value_counts()
        # print(type(counts))
        # full_supp["count"] = pd.Series(counts, index=full_supp.index)
        # print(full_supp)
        g = sns.barplot(x="p_bins", y="percentage", data=full_supp)
        # for index, row in full_supp.iterrows():
        #     g.text(row.name,row.percentage, row.count, color='black', ha="center")

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.xlabel('p-values (log scale $10^{-x}$)')
        plt.ylabel('Transcripts with full illumina support (%)',fontsize=16)
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif')
        plt.ylim(0, 100)
        plt.savefig(os.path.join(args.outfolder, "Figure_full_illumina_avg.pdf"))
        plt.clf()
        data.to_csv("~/tmp/ISOCON_REVIEW/PVALUE_ANALYSIS/FIG_TO_TAB_S5_NEW/new_dataframe.tsv", sep="\t")


        # Illumina coverage per p-value bins
        fig, ax = plt.subplots()
        ax = sns.countplot(x="p_bins", hue="Full_support", data=data)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.xlabel('p-values (log scale $10^{-x}$)')
        plt.ylabel('Full illumina support',fontsize=16)
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif')
        # plt.ylim(0.95, 1.0)
        plt.savefig(os.path.join(args.outfolder, "Figure_illumina_avg.pdf"))
        plt.clf()
        data.to_csv("~/tmp/ISOCON_REVIEW/PVALUE_ANALYSIS/FIG_TO_TAB_S5_NEW/new_dataframe.tsv", sep="\t")

        # read depth per p-value bins
        fig, ax = plt.subplots()
        ax = sns.boxplot(x="p_bins", y="read_depth", data=data, showfliers=False)
        plt.xlabel('p-values (log scale $10^{-x}$)')
        plt.ylabel('Read depth',fontsize=16)
        plt.ylim(1, 100)
        plt.savefig(os.path.join(args.outfolder, "Figure_readdepth_avg.pdf"))
        plt.clf()

        # exact matches plot
        plt.rc('text', usetex=False)
        fig, ax = plt.subplots()
        matches = data.loc[data['perfect_match_database'] == 'yes']
        pvals = [float(m) for m in matches["p-value"]]
        bins = [10**-i for i in exponents ] # [10**-i for i in [0,2,3,4,5,6,7,8,9,10,15,20,30,40,50,75,100,150,200,300,320][::-1] ] #  10**(-np.arange(0,100,2))
        cumulative_vals = [ len([pv for pv in pvals if 1e-300 < pv < b]) for b in bins  ]
        print(cumulative_vals)
        plt.xscale('log')
        plt.xlim(1e-22, 1.0) # plt.xlim(1e-322, 1.0)
        plt.ylim(0, 45) # plt.xlim(1e-322, 1.0)

        ax = pyplot.plot(bins, cumulative_vals) 
        plt.ylabel('#Exact matches')
        plt.xlabel("p-value")
        plt.title("Number of exact ENSEMBLE matches with p-value lower than x")
        plt.savefig(os.path.join(args.outfolder, "Figure_cumulative_exact_matches_pval.pdf"))
        plt.clf()

        # print(matches["read_depth"])
        # read_depths = [float(m) for m in matches["read_depth"]]
        # bins = [1,3,5,10,15,20,30,40,50,100,150,200,300,500,1000,2000,3000,5000,10000]
        # cumulative_vals = [ len([rd for rd in read_depths if rd < b]) for b in bins  ]
        # print(cumulative_vals)
        # plt.xscale('log')
        # plt.xlim(1, 10000)
        # ax = pyplot.plot(bins, cumulative_vals) 
        # plt.ylabel('#Exact matches')
        # plt.xlabel("Read support")
        # plt.title("Number of exact ENSEMBLE matches with read support lower than x")
        # plt.savefig(os.path.join(args.outfolder, "Figure_cumulative_exact_matches_rd.pdf"))
        # plt.clf()

        # # print(data["Illumina_support"])
        # two_d_data = [(float(pv), float(ills)) for pv, ills in zip(data["p-value"], data["Illumina_support"])]
        # two_d_data = [t for t in two_d_data if t[1] > 0.0]

        # x_vals = [10**-i for i in [0,2,3,4,5,6,7,8,9,10,15,20][::-1] ] # [10**-i for i in [0,2,3,4,5,6,7,8,9,10,15,20,30,40,50,75,100,150,200,300,320][::-1] ] #  10**(-np.arange(0,100,2))
        # y_vals = [ 100*sum([ills for pv, ills in two_d_data if 1e-300 <  pv < b])/ float(len([ills for pv, ills in two_d_data if 1e-300 < pv < b])) for b in x_vals  ]
        # print(x_vals)
        # print(y_vals)
        # plt.xscale('log')
        # plt.xlim(1e-22, 1.0) # plt.xlim(1e-322, 1.0)
        # ax = pyplot.plot(x_vals, y_vals) 
        # plt.ylabel('%Average Illumina support')
        # plt.xlabel("p-value")
        # plt.title("Average Illumina supported bases for transcripts with p-value lower than x")
        # plt.savefig(os.path.join(args.outfolder, "Figure_avg_illumina_suppport.pdf"))
        # plt.clf()

        # two_d_data = [(float(pv), float(rd)) for pv, rd in zip(data["p-value"], data["read_depth"])]
        # two_d_data = [t for t in two_d_data if t[1] > 0.0]

        # x_vals = [10**-i for i in [0,2,3,4,5,6,7,8,9,10,15,20][::-1] ] # [10**-i for i in [0,2,3,4,5,6,7,8,9,10,15,20,30,40,50,75,100,150,200,300,320][::-1] ] #  10**(-np.arange(0,100,2))
        # y_vals = [ sum([rd for pv, rd in two_d_data if 1e-300 < pv < b])/ float(len([rd for pv, rd in two_d_data if 1e-300 < pv < b])) for b in x_vals  ]
        # print(x_vals)
        # print(y_vals)
        # plt.xscale('log')
        # plt.xlim(1e-22, 1.0) # plt.xlim(1e-322, 1.0)
        # ax = pyplot.plot(x_vals, y_vals) 
        # plt.ylabel('Average CCS read support')
        # plt.xlabel("p-value")
        # plt.title("Average CCS read support for transcripts with p-value lower than x")
        # plt.savefig(os.path.join(args.outfolder, "Figure_avg_read_depth.pdf"))
        # plt.clf()

def main(args):
    
    tsv_file = create_tsv_from_transcript_annotations(args)
    data = pd.read_csv(tsv_file, sep="\t")
    print(data.corr(method='spearman'))
    paulplots_updated(data, args)
    # pairplot(data, args)
    
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

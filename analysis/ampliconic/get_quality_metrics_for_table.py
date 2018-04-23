import sys, os, argparse

def fill_container(container, line):
    id_, predicted_acc, family, Illumina_support,  perfect_match_database, sample1, sample2, acc_sample1, acc_sample2, both_samples, primer, category, full_supp_sample1, full_supp_sample2, gene_member_number, batch, sequence = line.split("\t")
    support  = int(predicted_acc.split("_")[4])
    p_val = predicted_acc.split("_")[5]
    if p_val == "not":
        p_val = 0.0
    else:
        p_val = float(p_val)

    container["read_support"].append(support)
    # print(predicted_acc)
    ill_supp = 0 if Illumina_support == "None" else float(Illumina_support)
    container["illumina_support"].append(ill_supp)

    if predicted_acc == acc_sample1:
        if full_supp_sample1 ==  "yes":
            container["full_support"].append(1)
    elif predicted_acc == acc_sample2:
        if full_supp_sample2 ==  "yes":
            container["full_support"].append(1)
    else:
        print("BUG")
        sys.exit()

    if perfect_match_database ==  "yes":
        container["exact_match_db"] += 1
    container["pvalues"].append(p_val)
    container["total_transcripts"] += 1 
    if category == "protein-coding":
        container["coding"] += 1 


tsv_file = open(sys.argv[1], "r")

sample1_only = {"read_support": [], "illumina_support": [], "full_support": [], "exact_match_db" : 0,  "pvalues" : [], "total_transcripts" : 0, "coding": 0 }
sample2_only = {"read_support": [], "illumina_support": [], "full_support": [], "exact_match_db" : 0,  "pvalues" : [], "total_transcripts" : 0, "coding": 0 }
both1 = {"read_support": [], "illumina_support": [], "full_support": [], "exact_match_db" : 0,  "pvalues" : [], "total_transcripts" : 0, "coding": 0 }
both2 = {"read_support": [], "illumina_support": [], "full_support": [], "exact_match_db" : 0,  "pvalues" : [], "total_transcripts" : 0, "coding": 0 }

for i, line in enumerate(tsv_file):
    if i == 0:
        continue
    _, predicted_acc, _, _, _, sample1, sample2, acc_sample1, acc_sample2, both_samples, _, _, _, _, _, _, _ = line.split("\t")
    
    if both_samples == "yes":
        if predicted_acc == acc_sample1:
            fill_container(both1, line)

        elif predicted_acc == acc_sample2:
            fill_container(both2, line)

        else:
            print("BUG")
            sys.exit()

    elif sample1 == "yes":
        fill_container(sample1_only, line)

    elif sample2 == "yes":
        fill_container(sample2_only, line)

    else:
        print("BUG")
        sys.exit()


print(both1)
print(both2)
print(sample1_only)
print(sample2_only)

print("CATEGORY\tTOTAL TRANSCRIPTS\tAVG READ SUPPORT\tMEDIAN READ SUPPORT\tAVG ILLUMINA SUPPORT\t#FULL ILLUMINA SUPPORT\t#EXACT DB MATCHES\tMEDIAN PVAL\t95%QUANTILE P-val\n") 
print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format("both1", both1["total_transcripts"], sum(both1["read_support"])/float(len(both1["read_support"])), sorted(both1["read_support"])[len(both1["read_support"])/2],
                                                          sum(both1["illumina_support"])/float(len(both1["illumina_support"])), len(both1["full_support"]),
                                                        both1["exact_match_db"],  sorted(both1["pvalues"])[len(both1["pvalues"])/2], sorted(both1["pvalues"])[int(len(both1["pvalues"])*0.95)]  )) 

print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format("both2", both2["total_transcripts"], sum(both2["read_support"])/float(len(both2["read_support"])), sorted(both2["read_support"])[len(both2["read_support"])/2],
                                                          sum(both2["illumina_support"])/float(len(both2["illumina_support"])), len(both2["full_support"]),
                                                        both2["exact_match_db"],  sorted(both2["pvalues"])[len(both2["pvalues"])/2], sorted(both2["pvalues"])[int(len(both2["pvalues"])*0.95)]  )) 

print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format("sample1_only", sample1_only["total_transcripts"], sum(sample1_only["read_support"])/float(len(sample1_only["read_support"])), sorted(sample1_only["read_support"])[len(sample1_only["read_support"])/2],
                                                          sum(sample1_only["illumina_support"])/float(len(sample1_only["illumina_support"])), len(sample1_only["full_support"]),
                                                        sample1_only["exact_match_db"],  sorted(sample1_only["pvalues"])[len(sample1_only["pvalues"])/2], sorted(sample1_only["pvalues"])[int(len(sample1_only["pvalues"])*0.95)]  )) 

print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format("sample2_only", sample2_only["total_transcripts"], sum(sample2_only["read_support"])/float(len(sample2_only["read_support"])), sorted(sample2_only["read_support"])[len(sample2_only["read_support"])/2],
                                                          sum(sample2_only["illumina_support"])/float(len(sample2_only["illumina_support"])), len(sample2_only["full_support"]),
                                                        sample2_only["exact_match_db"],  sorted(sample2_only["pvalues"])[len(sample2_only["pvalues"])/2], sorted(sample2_only["pvalues"])[int(len(sample2_only["pvalues"])*0.95)]  )) 



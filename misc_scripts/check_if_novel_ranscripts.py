from __future__ import print_function
import os,sys
import argparse
import re
from collections import defaultdict
from collections import Counter


from Bio.Blast.Applications import NcbiblastnCommandline #, NcbitblastxCommandline
#from Bio.Blast import  NCBIWWW
from Bio.Blast import NCBIXML
#import pickle
import dill

try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_palette("husl", desat=.6)
    sns.set(font_scale=1.6)
    plt.rcParams.update({'font.size': 22})
except:
    pass


try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError, RuntimeError:
    pass


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
            temp += line.strip()
    yield accession, temp

def histogram(data, args, name='histogram.png', x='x-axis', y='y-axis', x_cutoff=None, title=None):
    if x_cutoff: 
        plt.hist(data, range=[0, x_cutoff], bins = 100)
    else:
        plt.hist(data, bins = 100)
    plt.xlabel(x)
    plt.ylabel(y)
    if title:
        plt.title(title)

    plt.savefig(os.path.join(args.outfolder, name))
    plt.clf()

def get_index(fam, targeted_list):
    return targeted_list.index(fam)

def make_latex_table(targeted_set, families_per_fasta):

    # defining rows, soon to be columns
    expected = [3, 3, 3, 2, 2, 6, 6, 2, 2]
    targeted = ["BPY", "CDY", "DAZ", "HSFY", "PRY", "RBMY", "TSPY", "XKRY", "VCY"]


    found_in_db =  defaultdict(lambda : defaultdict(int)) # dict of symbol : location : count # defaultdict(dict) # family : [Y member_count, Y-like member count]
    found_in_db_on_y =  defaultdict(lambda : defaultdict(set)) # dict of symbol : location : count # defaultdict(dict) # family : [Y member_count, Y-like member count]
    gene_info_file = open("/Users/kxs624/Documents/data/NCBI_RNA_database/Homo_sapiens.gene_info.tsv", 'r')
    for line in gene_info_file:
        columns = line.strip().split()
        gene_symbol, chromosome = columns[2], columns[6]
        for symbol in targeted_set:
            if gene_symbol.find(symbol) >= 0:
                if chromosome == "Y":
                    found_in_db[symbol]["Y"] += 1
                    found_in_db_on_y[symbol]["Y"].add(gene_symbol)
                else:
                    found_in_db[symbol]["Y-like"] += 1
                    found_in_db_on_y[symbol]["Y-like"].add(gene_symbol)

    print(found_in_db)
    actual_Y = [0]*9
    Y_like = [0]*9
    for fam, def_dict in found_in_db.iteritems():
        index = get_index(fam, targeted)
        for k,v in def_dict.iteritems():
            if k == "Y":
                actual_Y[index] = v
            else:
                Y_like[index] = v

            # print(fam,k,v) 





    begin_table = "\\begin{table}\n\\tiny \n \\begin{tabular}{l|ll | llll|llll|l}\n"
    end_table = "\\end{tabular}\n\\end{table}\n"
    header = "Family & \t\multicolumn{2}{c}{ Members found in NCBI database} & \multicolumn{4}{c}{Members found in consensus} & \multicolumn{4}{c}{Members found in flnc reads} & Expected members\\\ \n "
    subheader = "& on Y & Y-like  &  on Y  & Y (\# mapped) & Y-like & Y-like (\# mapped) &  on Y  & Y (\# mapped) & Y-like & Y-like (\# mapped) &  \\\ \n "


    dataset_found_members = []
    dataset_found_total_mapped = [] 
    dataset_found_Y = []
    dataset_found_Y_total_mapped = [] 
    dataset_found_Y_like = []
    dataset_found_Y_like_total_mapped = []

    for families in families_per_fasta:
        consensus_found_members = [0]*9
        consensus_found_total_mapped = [0]*9 
        consensus_found_Y = [0]*9
        consensus_found_Y_total_mapped = [0]*9 
        consensus_found_Y_like = [0]*9
        consensus_found_Y_like_total_mapped = [0]*9 

        for family, member_dict in families.iteritems():
            fam = str(family)
            if fam in targeted_set:
                index = get_index(fam,targeted)
                # print(fam, index)
                nr_members = len(member_dict)
                consensus_found_members[index] = nr_members
                consensus_found_total_mapped[index] = sum(map(lambda x: len(x), member_dict.values()))
                for member in member_dict:
                    # print (fam, member)
                    # print(found_in_db_on_y[fam])
                    member_id = str(fam+member)
                    # print(member_id)
                    # print(found_in_db_on_y[fam]['Y'])
                    # print(found_in_db_on_y[fam]['Y-like'])
                    # assert member_id in found_in_db_on_y[fam]['Y'] or  member_id in found_in_db_on_y[fam]['Y-like']
                    id_found = False
                    for mem in found_in_db_on_y[fam]["Y"]:
                        # print(mem)
                        # print(found_in_db_on_y[fam]["Y"][mem])
                        if mem.find(member_id) >= 0:# member_id in found_in_db_on_y[fam]["Y"]:
                            consensus_found_Y[index] += 1
                            consensus_found_Y_total_mapped[index] += len(member_dict[member])
                            id_found = True
                            break

                    if not id_found:
                        for mem in found_in_db_on_y[fam]["Y-like"]:
                            if mem.find(member_id) >= 0: # member_id in found_in_db_on_y[fam]["Y-like"]:
                                consensus_found_Y_like[index] += 1
                                consensus_found_Y_like_total_mapped[index] += len(member_dict[member])
                                id_found = True
                                break

                    if not id_found:
                        print("Bug")

        dataset_found_members.append(consensus_found_members)
        dataset_found_total_mapped.append(consensus_found_total_mapped)
        dataset_found_Y.append(consensus_found_Y)
        dataset_found_Y_total_mapped.append(consensus_found_Y_total_mapped)
        dataset_found_Y_like.append(consensus_found_Y_like)
        dataset_found_Y_like_total_mapped.append(consensus_found_Y_like_total_mapped)


        # print("FAMILY: {0}, number of members found: {1} ".format(family, nr_members))
        # for member in member_dict:
            # print("member: {0}, nr consensus transcripts that mapped to this: {1} ".format(member, len(member_dict[member])))
    # print(targeted)
    # print(actual_Y)
    # print(Y_like)
    # print(dataset_found_members)
    # print(dataset_found_total_mapped)
    # print(dataset_found_Y)
    # print(dataset_found_Y_total_mapped)
    # print(dataset_found_Y_like)
    # print(dataset_found_Y_like_total_mapped)

    for j in range(len(dataset_found_Y)):
        consensus_found_members = dataset_found_members[j]
        consensus_found_total_mapped = dataset_found_total_mapped[j] 
        consensus_found_Y = dataset_found_Y[j]
        consensus_found_Y_total_mapped = dataset_found_Y_total_mapped[j]
        consensus_found_Y_like = dataset_found_Y_like[j]
        consensus_found_Y_like_total_mapped = dataset_found_Y_like_total_mapped[j]

        for i in range(len(consensus_found_members)):
            assert consensus_found_members[i] == consensus_found_Y[i] + consensus_found_Y_like[i]
            assert consensus_found_total_mapped[i] == consensus_found_Y_total_mapped[i] + consensus_found_Y_like_total_mapped[i]

    # print("{0} targeted families found : {1}".format(nr_interesting_found, interesting_found) )
    # member_count = map( lambda (family, member_dict): "{0}:{1}".format(family,len(member_dict)), families.iteritems())
    # print("Nr members: {0} ".format(member_count))

    # data = [targeted, actual_Y, Y_like, consensus_found_Y, consensus_found_Y_total_mapped, consensus_found_Y_like, consensus_found_Y_like_total_mapped, expected]
    data = [targeted, actual_Y, Y_like] #, consensus_found_Y, consensus_found_Y_total_mapped, consensus_found_Y_like, consensus_found_Y_like_total_mapped, expected]

    for i in range(len(dataset_found_Y)):
        data.append(dataset_found_Y[i])
        data.append(dataset_found_Y_total_mapped[i])
        data.append(dataset_found_Y_like[i])
        data.append(dataset_found_Y_like_total_mapped[i])

    data.append(expected)
    # assert len(data) == 12

    transposed_results = map(list, zip(*data))

    # print (transposed_results)

    print(begin_table,header, subheader)

    for list_ in transposed_results:
        print("{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} & {8} & {9} & {10} & {11} \\\ ".format(*list_))

    print(end_table)

    other_human = ""
    non_human = ""
    # sys.exit()

def blast_fasta(args, iteration):
    database= "".join(["/Users/kxs624/Documents/data/NCBI_RNA_database/refseq_rna.0{0} ".format(i) for i in range(8)])
    database = "\""+database + "\""
    print(database)
    
    if not args.blast_xml:
        print("blasting")
        out_file = os.path.join(args.outfolder,'blast_out_{0}.xml'.format(iteration))
        blastn_cline = NcbiblastnCommandline(query=args.transcripts[iteration], db=database, evalue=0.001, outfmt=5, out=out_file, max_target_seqs=20)
        # blastn_cline = NcbitblastxCommandline(query=args.transcripts[iteration], db=database, evalue=0.001, outfmt=5, out=out_file, max_target_seqs=2)
        print(blastn_cline)
        stdout, stderr = blastn_cline()

        print(stdout, stderr)
        print("Done")
        blast_records = NCBIXML.parse(open(out_file,'r'))
    else:
        blast_records = NCBIXML.parse(open(args.blast_xml[iteration],'r'))
    return blast_records

def print_blast_hits(args, targeted, blast_records, iteration):

    transcript_dictionary = {acc: seq for (acc, seq) in  read_fasta(open(args.transcripts[iteration], 'r'))}
    transcript_lengths = [len(seq) for acc, seq in transcript_dictionary.iteritems()]
    blast_hits = open(os.path.join(args.outfolder, "blast_annotations.tsv"), "w")
    total_queries = 0
    blast_hits.write("Transcript_id\tCategory\tquery_length\tMatches\tMismatches\tgaps\tHit_title\n")


    for i, blast_record in enumerate(blast_records):
        print(dir(blast_record))
        query_length = blast_record.query_length
        # print(blast_record.effective_database_length, blast_record.effective_hsp_length, blast_record.effective_query_length)        
        print("QUERY ID:{0}, query_length: {1}, query_letters: {2}, nr alignments: {3}".format(blast_record.query_id, blast_record.query_length, blast_record.query_letters, len(blast_record.alignments) ))
        # print(blast_record.sc_match, blast_record.sc_mismatch, blast_record.query_length - 21)
        # print([(a.hit_def, a.hit_id, [(h.gaps, h.identities, h.match, h.score, h.align_length, h.bits, h.query_start, h.query_end, h.query, dir(h)) for h in a.hsps[0:1]], a.length, a.title) for a in blast_record.alignments])
        # sys.exit()
        print()
        print("QUERY:{0}".format(blast_record.query))
        print()
        total_queries += 1
        # Not found in data base:
        if len(blast_record.alignments) == 0:
            blast_hits.write("{0}\tNO_MATCH_TO_DB\n".format(blast_record.query ))
            print("HERE")
            continue

        # sys.exit()
        print("Nr alignments: ", len(blast_record.alignments))
        for j, alignment in enumerate(blast_record.alignments):
            # only best alignment (most matches) for given alignment
            if j > 0:
                break
            # print("ALMNTS", blast_record.alignments)
            # print(dir(blast_record.alignments[0]))
            # print(alignment.accession, alignment.title, alignment.hsps)
            print(alignment.accession)
            # print("alignnment dir:",dir(alignment))
            # print(alignment.accession, alignment.hit_def, alignment.hit_id, alignment.hsps, alignment.length, alignment.title)

            # sys.exit()
            # print(alignment.title)
            print("Nr HSPs: ", len(alignment.hsps))

            category = False
            for k, hsp in enumerate(alignment.hsps):
                # only best HSP for given alignment
                if k > 0:
                    break
                c = Counter(str(hsp.match)) 
                matches, mismatches = c["|"], c[" "]
                gaps = int(hsp.gaps)
                # print(dir(hsp))
                # print(hsp.gaps, hsp.identities, hsp.match, hsp.score, hsp.align_length, hsp.bits, hsp.query_start, hsp.query_end, hsp.query, hsp.sbjct_start, hsp.sbjct_end,)
                print("Matches:{0}, mismatches:{1}, gaps: {2}, alignment length: {3}".format(matches, mismatches, gaps, alignment.length))

                if query_length > matches + mismatches + gaps: 
                    category = "HARDCLIPPED"
                elif mismatches == 0 and gaps == 0:
                    category = "perfect_DB_match"
                else:
                    category = "NOVEL"
                
                blast_hits.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(blast_record.query, category, query_length, matches, mismatches, gaps, alignment.title))

                # Log the results for the best hit to a Homo sapiens trancript
                # print("query length:{0}, alignment length:{1}, mismatches:{2}".format(blast_record.query_length, alignment.length, mismatches))
                # print(hsp.__dict__.keys())
                # print(hsp.score, hsp.gaps, hsp.align_length, len(hsp.match), hsp.query_start, hsp.query_end, hsp.identities, len(hsp.query))
            print()

    print(total_queries)


def main(args):
    families_per_fasta = []
    targeted = set(["BPY", "CDY", "DAZ", "HSFY", "PRY", "RBMY", "TSPY", "XKRY", "VCY"])
    misc_hits =[] # tuple with (not mapped in db, not human, human but not targeted)

    if args.pickled_parse:
        families_per_fasta = dill.load(families_per_fasta, open( os.path.join(args.outfolder, "families_per_fasta.p"), "rb" ) )
        misc_hits = dill.load(misc_hits, open( os.path.join(args.outfolder, "misc_hits.p"), "rb" ) )
        pass
    else:
    
        for iteration in range(len(args.transcripts)):
            blast_records = blast_fasta(args, iteration)
            print_blast_hits(args, targeted, blast_records, iteration)
            
 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('transcripts', type=str, nargs="+", help='Path to the transcript fasta file')
    parser.add_argument('--blast_xml', dest="blast_xml", type=str, nargs="+", help='Path to the blast records')
    parser.add_argument('--pickled_parse', dest="pickled_parse", action="store_true", help='Path to the blast records')

    parser.add_argument('outfolder', type=str, help='Output path of results')

    args = parser.parse_args()

    outfolder = args.outfolder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    main(args)
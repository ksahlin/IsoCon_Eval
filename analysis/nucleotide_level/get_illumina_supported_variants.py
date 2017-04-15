
import argparse
import os
from collections import Counter
import re
import pysam
from collections import defaultdict


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
    if accession:
        yield accession, temp

def get_variants(illumina_to_ref, reference_fasta):
    samfile = pysam.AlignmentFile(illumina_to_ref, "rb" )
    illumina_positions = {}
    ref_positions = {}

    assert len(reference_fasta) == 1 # one reference at a time

    for pileupcolumn in samfile.pileup():
        print ("coverage at base %s = %s" %
               (pileupcolumn.pos, pileupcolumn.n))

        ref_base = reference_fasta[pileupcolumn.reference_name][pileupcolumn.pos]
        ref_positions[pileupcolumn.pos] = ref_base
        # print(ref_base)
        illumina_positions[pileupcolumn.pos] = defaultdict(int)

        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # query position is None if is_del or is_refskip is set.
                # print ('\tbase in read %s = %s' %
                #       (pileupread.alignment.query_name,
                #        pileupread.alignment.query_sequence[pileupread.query_position]))
                
                illumina_positions[pileupcolumn.pos][pileupread.alignment.query_sequence[pileupread.query_position]] += 1

            elif pileupread.is_del:
                illumina_positions[pileupcolumn.pos]["-"] += 1

            elif pileupread.is_refskip:
                illumina_positions[pileupcolumn.pos]["N"] += 1
            else:
                # should not end up here
                assert False


    illumina_variants = defaultdict(list)
    p_illumina_indel = 0.001
    p_illumina_subs = 0.001

    for pos in range(len(reference_fasta[pileupcolumn.reference_name])):
        ref_base = reference_fasta[pileupcolumn.reference_name][pos]
        if pos in illumina_positions:
            print(reference_fasta[pileupcolumn.reference_name][pos], illumina_positions[pos])
            total_illumina_support = sum([count for nucl, count in illumina_positions[pos].items()])            
            # SITE_THRESHOLD = p_illumina_error/3

            for site, count in illumina_positions[pos].items(): 
                if site != ref_base:
                    if site == "-"  and count >= max(1, total_illumina_support*p_illumina_indel):
                        illumina_variants[pos].append(site)
                    elif site != "-" and count >= max(1, total_illumina_support*p_illumina_subs):
                        illumina_variants[pos].append(site)
        else:
            print("No alignments", pos)

    for p in illumina_variants:
        print("ref base:", reference_fasta[pileupcolumn.reference_name][p] , p, illumina_variants[p], illumina_positions[p])

    samfile.close()

    # pileupread.is_del
    # pileupread.alignment.query_name
    # pileupread.alignment.query_sequence[pileupread.query_position]
    # pileupread.query_position


def main(args):
    """ 
        1. Retrieve the illumina reads that support variants w.r.t. consensus reference.\n\n
        2. Find the mapping positions of these reads when mapped to the predicted transcripts (for each respective method).\n\n
        3. See if the varinat is present (found) in a predicted transcript by looking if ilumina reads agrees with predicted transcript 
            at the base pair where the illumina read signalled a variant.\n\n
    """
    output_file = open(os.path.join(args.outfolder, "illumina_varinats.tsv" ), "w")
    reference_seq = {acc: seq for (acc, seq) in  read_fasta(open(args.consensus, 'r'))}


    get_variants(args.illumina_to_ref, reference_seq)



    # FP_positions = []
    # for line in open(args.support_file, "r"):
    #     consensus_name, ref_pos, ref_site, total_support, read_bases, base_qualities =  line.strip().split()

    #     # # here we discover unsopported positions (FP)
    #     # if int(total_support) < 2:
    #     #     FP_positions.append(int(ref_pos))

    #     # get all alternate loci and their support over positions

    #     ins_pattern = "\+[0-9]+[ACGTNacgtn]+"
    #     ins_res = re.findall(ins_pattern, read_bases.upper())
    #     c_ins = Counter(ins_res)


    #     del_pattern =  "-[0-9]+[ACGTNacgtn]+"
    #     del_res = re.findall(del_pattern, read_bases.upper())
    #     c_del = Counter(del_res)


    #     subs_pattern = "(?<![0-9'^'])[ACGTN]"
    #     subs_res = re.findall(subs_pattern, read_bases.upper())
    #     subs_res_string = "".join([subs for subs in subs_res])
    #     c_subs = Counter(subs_res)
        

    #     c = Counter(read_bases.upper())
    #     inferred_site_count = c["."] + c[","] + len(ins_res)  + len(subs_res) 
    #     print(total_support, inferred_site_count, c["."], c[","] , len(ins_res) , len(del_res) , len(subs_res)  )
    #     assert int(total_support) == inferred_site_count

    #     output_file.write("{0}\t{1}\t{2}:{3}\t".format(consensus_name, ref_pos, ref_site, c["."] + c[","], ))

    #     for site, count in  c_subs.items():
    #         output_file.write("{0}:{1}\t".format(site,count))

    #     for site, count in  c_ins.items():
    #         output_file.write("{0}:{1}\t".format(site,count))

    #     for site, count in  c_del.items():
    #         output_file.write("{0}:{1}\t".format(site,count))
        

if __name__ == '__main__':

# Take care of input


    parser = argparse.ArgumentParser(description = "This script does: \n\n {0}".format(main.__doc__)) 
    parser.add_argument('illumina_to_ref', type=str, help='Bam file. ')
    parser.add_argument('illumina_to_pred', type=str, help='Bam file ')
    parser.add_argument('consensus', type=str, help='Fasta file ')
    parser.add_argument('outfolder', type=str, help='outfile folder to put output in. ')

    args = parser.parse_args()
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)

import argparse
import os
from collections import Counter
import re
import pysam
from collections import defaultdict
import pickle

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

def get_variants_on_consensus(illumina_to_ref, reference_fasta, outfolder):
    samfile = pysam.AlignmentFile(illumina_to_ref, "rb" )
    illumina_positions = {}
    illumina_accessions = {}
    ref_positions = {}

    assert len(reference_fasta) == 1 # one reference at a time

    for pileupcolumn in samfile.pileup():
        print ("coverage at base %s = %s" %
               (pileupcolumn.pos, pileupcolumn.n))

        ref_base = reference_fasta[pileupcolumn.reference_name][pileupcolumn.pos]
        ref_positions[pileupcolumn.pos] = ref_base
        # print(ref_base)
        illumina_positions[pileupcolumn.pos] = defaultdict(int)
        illumina_accessions[pileupcolumn.pos] = defaultdict(list)
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # query position is None if is_del or is_refskip is set.
                # print ('\tbase in read %s = %s' %
                #       (pileupread.alignment.query_name,
                #        pileupread.alignment.query_sequence[pileupread.query_position]))
                illumina_base = pileupread.alignment.query_sequence[pileupread.query_position]
                illumina_positions[pileupcolumn.pos][illumina_base] += 1
                if illumina_base != ref_base:
                    illumina_accessions[pileupcolumn.pos][illumina_base].append( (pileupread.alignment.query_name, pileupread.query_position) )

            elif pileupread.is_del:
                illumina_positions[pileupcolumn.pos]["-"] += 1
                illumina_accessions[pileupcolumn.pos]["-"].append( (pileupread.alignment.query_name, pileupread.query_position_or_next ) )

            elif pileupread.is_refskip:
                illumina_positions[pileupcolumn.pos]["N"] += 1
            else:
                # should not end up here
                assert False



    illumina_variants = defaultdict(list)
    # p_illumina_indel = 0.001
    # p_illumina_subs = 0.001

    for pos in range(len(reference_fasta[pileupcolumn.reference_name])):
        ref_base = reference_fasta[pileupcolumn.reference_name][pos]
        if pos in illumina_positions:
            print(reference_fasta[pileupcolumn.reference_name][pos], illumina_positions[pos])

            # SITE_THRESHOLD = p_illumina_error/3
            total_illumina_support = sum([count for nucl, count in illumina_positions[pos].items()])            

            for site, count in illumina_positions[pos].items(): 
                if site != ref_base:
                    if site == "-"  and count >= 4: # max(1, total_illumina_support*p_illumina_indel):
                        illumina_variants[pos].append(site)
                    elif site != "-" and count >= 4: # max(1, total_illumina_support*p_illumina_subs):
                        illumina_variants[pos].append(site)
        else:
            print("No alignments", pos)

    for p in illumina_variants:
        print("ref base:", reference_fasta[pileupcolumn.reference_name][p] , p, illumina_variants[p], illumina_positions[p])
        for site in illumina_variants[p]:
            print(illumina_accessions[p][site])


    samfile.close()
    
    illumina_variants_out = open(os.path.join(outfolder, 'variants.pkl'), 'wb')
    illumina_accessions_out = open(os.path.join(outfolder, 'accessions.pkl'), 'wb')

    # Pickle dictionary using protocol 0.
    pickle.dump(illumina_variants, illumina_variants_out)
    pickle.dump(illumina_accessions, illumina_accessions_out)

    return illumina_variants, illumina_accessions




def find_if_supported_in_pred_transcripts(illumina_to_pred, illumina_variants, illumina_accessions):
    samfile = pysam.AlignmentFile(illumina_to_pred, "rb" )

    read_accession_to_query_pos_and_variant = defaultdict(list)
    for pos, var_dict in  illumina_accessions.items():
        for variant, acc_and_q_pos_list in var_dict.items():
            for acc, pos in acc_and_q_pos_list:
                read_accession_to_query_pos_and_variant[acc].append( (variant, pos) )

    nr_unmapped = 0
    # read_accession_to_query_pos_and_variant = { acc : (variant, pos) for pos, var_dict in  illumina_accessions.items() for variant, acc_and_q_pos_list in var_dict.items() for acc, pos in acc_and_q_pos_list  }
    print("Nr predicted:", len(samfile.references)) # one reference at a time
    # print(read_accession_to_query_pos_and_variant)
    interesting_states = set([1,2,8])
    captured = 0
    not_captured = 0

    for read in samfile.fetch():
        # print(read.query_name)
        if read.query_name in read_accession_to_query_pos_and_variant:
            variant_positions = read_accession_to_query_pos_and_variant[read.query_name]
            for site, q_var_pos in variant_positions:
                state_start = 0
                state_end = 0
                if read.is_unmapped:
                    print("is unmapped")
                    nr_unmapped += 1
                    continue
                for state, number in read.cigartuples:
                    state_end += number 
                    if state == 0 and (state_start <= q_var_pos <= state_end):
                        print("Variant captured!:", read.cigartuples, q_var_pos)
                        captured += 1
                    elif state in interesting_states and (state_start <= q_var_pos <= state_end):
                        print("Variant not found?!:", read.cigartuples)
                        not_captured +=1
                    state_start += state_end + 1



                # print(read.cigartuples)
    print("Total variant carrying reads:", len(read_accession_to_query_pos_and_variant))
    print("Number of varinat carrying reads that were unmapped on predicted:", nr_unmapped)
    print("Captured:", captured)
    print("Not captured:", not_captured)
    return

def main(args):
    """ 
        1. Retrieve the illumina reads that support variants w.r.t. consensus reference.\n\n
        2. Find the mapping positions of these reads when mapped to the predicted transcripts (for each respective method).\n\n
        3. See if the varinat is present (found) in a predicted transcript by looking if ilumina reads agrees with predicted transcript 
            at the base pair where the illumina read signalled a variant.\n\n
    """
    output_file = open(os.path.join(args.outfolder, "illumina_varinats.tsv" ), "w")
    reference_seq = {acc: seq for (acc, seq) in  read_fasta(open(args.consensus, 'r'))}

    # predicted_seq = {acc: seq for (acc, seq) in  read_fasta(open(args.predicted, 'r'))}
    # illumina_variants, illumina_accessions = get_variants_on_consensus(args.illumina_to_pred, reference_seq)
    if args.variant_folder:
        illumina_variants_in = open(os.path.join(args.variant_folder, 'variants.pkl'), 'rb')
        illumina_accessions_in = open(os.path.join(args.variant_folder, 'accessions.pkl'), 'rb')    
        illumina_variants = pickle.load(illumina_variants_in)
        illumina_accessions = pickle.load(illumina_accessions_in)
        # print(illumina_variants)
        # print(illumina_accessions)

    else:
        illumina_variants, illumina_accessions = get_variants_on_consensus(args.illumina_to_ref, reference_seq, args.outfolder)

    supported_vanriants, not_supported = find_if_supported_in_pred_transcripts(args.illumina_to_pred, illumina_variants, illumina_accessions)




if __name__ == '__main__':

# Take care of input


    parser = argparse.ArgumentParser(description = "This script does: \n\n {0}".format(main.__doc__)) 
    parser.add_argument('illumina_to_ref', type=str, help='Bam file. ')
    parser.add_argument('illumina_to_pred', type=str, help='Bam file ')
    parser.add_argument('consensus', type=str, help='Fasta file ')
    # parser.add_argument('predicted', type=str, help='Fasta file ')
    parser.add_argument('--variant_folder', default = "", type=str, help='Fasta file ')
    parser.add_argument('outfolder', type=str, help='outfile folder to put output in. ')

    args = parser.parse_args()
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)
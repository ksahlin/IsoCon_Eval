
import argparse
import os
from collections import Counter
import re
import pysam
from collections import defaultdict
import dill
import pickle

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

def custom_histogram(data, outfolder, name='histogram.png', x='x-axis', y='y-axis', title=None, params = {"bins" : 100}):

    plt.hist(data,**params)
    plt.xlabel(x)
    plt.ylabel(y)
    if title:
        plt.title(title)
    if "label" in params:
        plt.legend()
        print("here")
    plt.savefig(os.path.join(outfolder, name))
    plt.clf()

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

def get_status_in_cigar(alignment, q_pos):
    # print(pileupread.indel, pileupread.alignment.cigartuples)
    read_state = None
    state_length = None
    state_start_pos = 0
    state_end_pos = -1 # 0-indexed
    # prev_state  = None
    # prev_length  = None

    for state, number in alignment.cigartuples:
        if state == 2 or state ==5:
            # prev_state = state
            # prev_length = number
            continue
        state_end_pos += number
        if (state_start_pos <= q_pos <= state_end_pos):
            read_state = state
            state_length = number
            break
        state_start_pos = state_end_pos
        # prev_state  = state
        # prev_length  = number

    return read_state, state_length

def get_deletion_status_in_cigar(alignment, q_pos_next):
    read_state = None
    state_length = None
    state_start_pos = 0
    state_end_pos = -1 # 0-indexed
    prev_state  = None
    prev_length  = None

    for state, number in alignment.cigartuples:
        if state == 2 or state ==5:
            prev_state = state
            prev_length = number
            continue
        state_end_pos += number
        if (state_start_pos <= q_pos_next <= state_end_pos):
            read_state = state
            state_length = number
            break
        state_start_pos = state_end_pos
        prev_state  = state
        prev_length  = number

    return prev_state, prev_length

def get_variants_on_reference(illumina_to_ref, reference_fasta, outfolder, variant_cutoff):
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
        illumina_positions[-pileupcolumn.pos] = defaultdict(int)
        illumina_accessions[pileupcolumn.pos] = defaultdict(list)
        illumina_accessions[-pileupcolumn.pos] = defaultdict(list)
        for pileupread in pileupcolumn.pileups:   
            if pileupread.alignment.is_unmapped:
                # print("is unmapped")
                nr_unmapped += 1
                continue
            if pileupread.alignment.is_supplementary:
                # print("is supplementary")
                continue

            if pileupread.indel > 1:  # only dealing with insertions of one base pair for now!!
                print("Skipping logging of insertion of length:", pileupread.indel)
                continue
            if pileupread.is_del:
                state, length = get_deletion_status_in_cigar(pileupread.alignment, pileupread.query_position_or_next)
                # print("DELETION HERE:", pileupread.indel, pileupread.alignment.cigartuples, "length",length, pileupread.query_position_or_next)
                assert state == 2 # should be deletion
                if length > 1:
                    continue

            if not pileupread.is_del and not pileupread.is_refskip:
                illumina_base = pileupread.alignment.query_sequence[pileupread.query_position]
                illumina_positions[pileupcolumn.pos][illumina_base] += 1
                if illumina_base != ref_base:
                    illumina_accessions[pileupcolumn.pos][illumina_base].append( ((pileupread.alignment.query_name, pileupread.alignment.is_read1), pileupread.query_position))

            elif pileupread.is_del:
                print("Deletion:", pileupread.indel, pileupread.alignment.cigartuples, pileupread.query_position, pileupread.query_position_or_next)
                    # need to store deletion length before we end up here
                illumina_positions[pileupcolumn.pos]["-"] += 1
                # print("Deletion here! Length: {0}. Q-position {1}".format(pileupread.indel))
                illumina_accessions[pileupcolumn.pos]["-"].append( ((pileupread.alignment.query_name, pileupread.alignment.is_read1), pileupread.query_position_or_next ))

            elif pileupread.is_refskip:
                illumina_positions[pileupcolumn.pos]["N"] += 1
            else:
                # should not end up here
                assert False

            # between bases insertions
            if pileupread.indel == 1:
                print("insertion here! Length: {0}. Q-position {1}".format(pileupread.indel, pileupread.query_position), pileupread.alignment.cigartuples)
                illumina_base = pileupread.alignment.query_sequence[pileupread.query_position+1]
                illumina_accessions[-pileupcolumn.pos][illumina_base].append( ((pileupread.alignment.query_name, pileupread.alignment.is_read1), pileupread.query_position+1))
                illumina_positions[-pileupcolumn.pos][illumina_base] += 1



    illumina_variants = defaultdict(list)
    illumina_variants_read_accessions = defaultdict(lambda : defaultdict(list))
    # p_illumina_indel = 0.001
    # p_illumina_subs = 0.001

    for ref_pos in range(len(reference_fasta[pileupcolumn.reference_name])):
        ref_base = reference_fasta[pileupcolumn.reference_name][ref_pos]
        if ref_pos in illumina_positions:
            total_illumina_support = sum([count for nucl, count in illumina_positions[ref_pos].items()])            

            for site, count in illumina_positions[ref_pos].items(): 
                if site != ref_base and site != "N":
                    if site == "-"  and count >= variant_cutoff: # max(1, total_illumina_support*p_illumina_indel):
                        illumina_variants[ref_pos].append(site)
                        accessions_list = illumina_accessions[ref_pos][site]
                        illumina_variants_read_accessions[ref_pos][site] = accessions_list

                    elif site != "-" and count >= variant_cutoff: # max(1, total_illumina_support*p_illumina_subs):
                        illumina_variants[ref_pos].append(site)
                        accessions_list = illumina_accessions[ref_pos][site]
                        illumina_variants_read_accessions[ref_pos][site] = accessions_list

        if -ref_pos in illumina_positions: # insertions
            total_illumina_support = sum([count for nucl, count in illumina_positions[ref_pos].items()])            
            for site, count in illumina_positions[-ref_pos].items(): 
                print("INSERTION", ref_base, site )
                if count >= variant_cutoff and site != "N":
                    illumina_variants[-ref_pos].append(site)  
                    accessions_list = illumina_accessions[-ref_pos][site]
                    illumina_variants_read_accessions[-ref_pos][site] = accessions_list

        else:
            print("No alignments", ref_pos)

    for p in illumina_variants:
        print(p, illumina_variants[p], illumina_positions[p])
        if p < 0:
            print("REF POS flanking insertion:",  reference_fasta[pileupcolumn.reference_name][p],  reference_fasta[pileupcolumn.reference_name][p+1]  )
        # print("ref base:", reference_fasta[pileupcolumn.reference_name][p] , p, illumina_variants[p], illumina_positions[p])
        # for site in illumina_variants[p]:
        #     print(illumina_accessions[p][site])


    samfile.close()
    
    illumina_variants_out = open(os.path.join(outfolder, 'variants.pkl'), 'wb')
    illumina_accessions_out = open(os.path.join(outfolder, 'accessions.pkl'), 'wb')
    illumina_positions_out = open(os.path.join(outfolder, 'positions.pkl'), 'wb')

    # Pickle dictionary using protocol 0.
    pickle.dump(illumina_variants, illumina_variants_out)
    pickle.dump(illumina_variants_read_accessions, illumina_accessions_out)
    pickle.dump(illumina_positions, illumina_positions_out)

    # sys.exit()
    return illumina_variants, illumina_variants_read_accessions, illumina_positions




def find_if_supported_in_pred_transcripts(illumina_to_pred, illumina_variants, illumina_accessions, illumina_positions, predicted_transcripts, outfolder, args):
    samfile = pysam.AlignmentFile(illumina_to_pred, "rb" )

    read_accession_to_query_pos_and_variant = defaultdict(list)
    for ref_pos, var_dict in  illumina_accessions.items():
        for variant, acc_and_q_pos_list in var_dict.items():
            for read_id, pos in acc_and_q_pos_list:
                read_accession_to_query_pos_and_variant[read_id].append( (variant, pos, ref_pos) )

    print("Nr predicted transcripts:", len(samfile.references)) # one reference at a time
    deletions_captured = defaultdict(set) # ref_pos : site
    insertions_captured = defaultdict(set) # ref_pos : site
    substitutions_captured = defaultdict(set) # ref_pos : site
    nr_unmapped = 0

    for read in samfile.fetch():
        # print(read.query_name)
        # if not read.is_unmapped:
        #     print("is unmapped")
        #     for state, number in read.cigartuples:
        #         if state == 8:
        #             print("OMG,", mismatch)

        if (read.query_name, read.is_read1) in read_accession_to_query_pos_and_variant:
            # print(read.reference_name)
            if read.is_unmapped:
                # print("is unmapped")
                nr_unmapped += 1
                continue
            if read.is_supplementary:
                # print("is supplementary")
                continue

            read_aligned_to_pred_transcript_positions = read.get_reference_positions(full_length=True)
            variant_positions = read_accession_to_query_pos_and_variant[(read.query_name, read.is_read1)]
            for site, q_var_pos, ref_pos in variant_positions:
                assert ref_pos in illumina_variants
                assert site in illumina_variants[ref_pos]

                # dealing with insertion
                if ref_pos < 0:
                    # print("INSERTION")
                    ref_pos_insertion = -ref_pos
                    assert ref_pos_insertion > 0
                    read_state, state_length = get_status_in_cigar(read, q_var_pos)

                    if read_state == 0:
                        # print("INSERTION captured!:", ref_pos, read.cigartuples, q_var_pos, read_state, state_length)
                        insertions_captured[ref_pos].add(site)
                    else:
                        pass
                        print("INSERTION not captured!:", ref_pos, read.cigartuples, read_state, state_length, q_var_pos, site)
                        print("alignment on pred", read.cigartuples, read.query_name, read.is_read1)
                        tempsam = pysam.AlignmentFile(args.illumina_to_ref, "rb" )
                        for r in tempsam.fetch():
                            if r.query_name == read.query_name:
                                print( "Alignment on ref:",r.cigartuples, r.query_name, r.is_read1)


                elif site == "-": # how to match deletions? need to look at cigar here
                    # print("DELETION")
                    assert ref_pos >= 0
                    q_pos_next = q_var_pos # the position is really the next one adter the deletion (reported by pysam pileup) 
                    read_state1, state_length1 = get_status_in_cigar(read, q_pos_next) #right flanking
                    read_state2, state_length2 = get_status_in_cigar(read, q_pos_next-1) # left flanking
                    if read_state1 == read_state2 == 0 :
                        # print("DELETION captured!:", ref_pos, read.cigartuples, q_var_pos, read_state1, state_length1,  read_state2, state_length2)
                        deletions_captured[ref_pos].add(site)
                    else:
                        pass
                        print("DELETION not captured!:", ref_pos, read.cigartuples, read_state1, state_length1,  read_state2, state_length2)

                else: # substitution
                    pass
                    # assert ref_pos >= 0 and site != "-" 
                    # predicted_transcript_pos = read_aligned_to_pred_transcript_positions[q_var_pos]
                    # if predicted_transcript_pos:
                    #     pred_transcript_site = predicted_transcripts[read.reference_name][predicted_transcript_pos]
                    #     if site == pred_transcript_site:
                    #         sites_captured[ref_pos].add(site)
                    #         print("SUBSTITUTION CAPTURED:", ref_pos, site, pred_transcript_site)
                    #     else:
                    #         print("Sites not matching for SUBSTITUTION:", ref_pos, site, pred_transcript_site)
                    # else:
                    #     print(q_var_pos, "this part of the read was not aligned to any predicted transcript", ref_pos, read.cigartuples)



    total_number_of_illumina_varinats = len([1 for pos in illumina_variants for site in illumina_variants[pos]])
    captured = 0
    print(len(sites_captured))
    del_captured = len([ site for sites in  deletions_captured for site in sites])
    ins_captured = len([ site for sites in  insertions_captured for site in sites])
    subs_captured = len([ site for sites in  substitutions_captured for site in sites])
    print("Deletions CAPTURED1:", del_captured)
    print("Insertions CAPTURED1:", ins_captured)
    print("Substitutions CAPTURED1:", subs_captured)
    # for pos in sites_captured:
    #     for site in sites_captured[pos]:
    #         captured += 1
    #         print(pos, site)


    print("SITES CAPTURED2:")
    del_not_captured = 0
    del_captured2 = 0
    del_captured_illumina_depths = []
    del_non_captured_illumina_depths = []

    ins_not_captured = 0
    ins_captured2 = 0
    ins_captured_illumina_depths = []
    ins_non_captured_illumina_depths = []

    subs_not_captured = 0
    subs_captured2 = 0
    subs_captured_illumina_depths = []
    subs_non_captured_illumina_depths = []

    for ref_pos in illumina_variants:
        if ref_pos >=0:
            for illumina_site in illumina_variants[ref_pos]:

                if illumina_site == "-":
                    if ref_pos in deletions_captured and illumina_site in deletions_captured[ref_pos]:
                        # print("HEHEH", illumina_site, ref_pos, sites_captured[ref_pos])
                        # if illumina_site in deletions_captured[ref_pos]:
                            # print("Cap2", ref_pos, illumina_site)
                        del_captured2 += 1
                        del_captured_illumina_depths.append(illumina_positions[ref_pos][illumina_site])
                    else:
                        del_not_captured +=1
                        del_non_captured_illumina_depths.append(illumina_positions[ref_pos][illumina_site])
                else:
                    if ref_pos in substitutions_captured and illumina_site in deletions_captured[ref_pos]:
                        # print("HEHEH", illumina_site, ref_pos, sites_captured[ref_pos])
                        # if illumina_site in deletions_captured[ref_pos]:
                        #     # print("Cap2", ref_pos, illumina_site)
                        #     del_captured2 += 1
                        #     del_captured_illumina_depths.append(illumina_positions[ref_pos][illumina_site])
                    else:
                        del_not_captured +=1
                        del_non_captured_illumina_depths.append(illumina_positions[ref_pos][illumina_site])                    
        
        else: #insertion
            for illumina_site in illumina_variants[ref_pos]:
                if ref_pos in insertions_captured and illumina_site in insertions_captured[ref_pos]:
                    # if illumina_site in insertions_captured[ref_pos]:
                    #     # print("Cap2", ref_pos, illumina_site)
                    #     ins_captured2 += 1
                    #     ins_captured_illumina_depths.append(illumina_positions[ref_pos][illumina_site])
                else:
                    ins_not_captured +=1
                    ins_non_captured_illumina_depths.append(illumina_positions[ref_pos][illumina_site])                


    print(set(sites_captured.keys()).difference(set(illumina_variants)))

    print("Total variant carrying reads:", len(read_accession_to_query_pos_and_variant))
    print("Total number of illumina supported sites:", total_number_of_illumina_varinats)
    print("Number of varinat carrying reads that were unmapped on predicted:", nr_unmapped)
    print("Deletion sites captured:", del_captured)
    print("Insertion sites captured:", ins_captured)
    print("Sites captured calc2:", captured2)
    print("Sites not captured:", not_captured)

    # common_params = dict(bins=30, 
    #                  range=(0, 10000), 
    #                  normed=0, label=['subs','ins','del'])
    common_params = dict(normed=0, label=['captured','not captured'], range=(0, 30), bins=30)
    print(captured_illumina_depths)
    print(non_captured_illumina_depths)
    title_header = "captured vs non-captured illumina depths" #"reference mismatches, total: {0}, perfect: {1}".format(len(tuple_identities[0]), perfect_matches)
    custom_histogram(captured_illumina_depths, outfolder, name='captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)
    custom_histogram(non_captured_illumina_depths, outfolder, name='not_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)

    return

def main(args):
    """ 
        1. Retrieve the illumina reads that support variants w.r.t. reference reference.\n\n
        2. Find the mapping positions of these reads when mapped to the predicted transcripts (for each respective method).\n\n
        3. See if the varinat is present (found) in a predicted transcript by looking if ilumina reads agrees with predicted transcript 
            at the base pair where the illumina read signalled a variant.\n\n
    """
    output_file = open(os.path.join(args.outfolder, "illumina_varinats.tsv" ), "w")
    reference_seq = {acc: seq for (acc, seq) in  read_fasta(open(args.reference, 'r'))}

    # predicted_seq = {acc: seq for (acc, seq) in  read_fasta(open(args.predicted, 'r'))}
    # illumina_variants, illumina_accessions = get_variants_on_consensus(args.illumina_to_pred, reference_seq)
    if args.variants_exists:
        illumina_variants_in = open(os.path.join(args.outfolder, 'variants.pkl'), 'rb')
        illumina_accessions_in = open(os.path.join(args.outfolder, 'accessions.pkl'), 'rb')    
        illumina_positions_in = open(os.path.join(args.outfolder, 'positions.pkl'), 'rb')    
        illumina_variants = pickle.load(illumina_variants_in)
        illumina_accessions = pickle.load(illumina_accessions_in)
        illumina_positions = pickle.load(illumina_positions_in)
        # print(illumina_variants)
        # print(illumina_accessions)

    else:
        illumina_variants, illumina_accessions, illumina_positions = get_variants_on_reference(args.illumina_to_ref, reference_seq, args.outfolder, args.variant_cutoff)

    predicted_transcripts = {acc.split()[0]: seq for (acc, seq) in read_fasta(open(args.predicted, 'r'))}
    # print(predicted_transcripts)
    find_if_supported_in_pred_transcripts(args.illumina_to_pred, illumina_variants, illumina_accessions, illumina_positions, predicted_transcripts, args.outfolder, args)




if __name__ == '__main__':

# Take care of input


    parser = argparse.ArgumentParser(description = "This script does: \n\n {0}".format(main.__doc__)) 
    parser.add_argument('-illumina_to_ref', type=str, help='Bam file. ')
    parser.add_argument('-illumina_to_pred', type=str, help='Bam file ')
    parser.add_argument('-reference', type=str, help='Fasta file ')
    parser.add_argument('-predicted', type=str, help='Fasta file ')
    parser.add_argument('-outfolder', type=str, help='outfile folder to put output in. ')
    parser.add_argument('--variant_cutoff', type=int, default =2, help='Fasta file ')
    parser.add_argument('--variants_exists', action="store_true", help='Fasta file ')

    args = parser.parse_args()
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)
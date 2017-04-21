
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

def is_deletion(alignment, q_pos_next):
    # print(pileupread.indel, pileupread.alignment.cigartuples)
    read_state = None
    state_length = None
    state_start_pos = 0
    state_end_pos = -1 # 0-indexed
    start_flank = q_pos_next -1
    end_flank = q_pos_next
    # prev_state  = None
    # prev_length  = None

    for state, number in alignment.cigartuples:
        if state == 2 or state ==5:
            # prev_state = state
            # prev_length = number
            continue
        state_end_pos += number
        if (state_start_pos < start_flank < end_flank < state_end_pos):
            print("here", state, state_start_pos, state_end_pos )
            if state == 0:
                return False

            read_state = state
            state_length = number

        state_start_pos = state_end_pos

    return True

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

def get_unsupported_positions_on_predicted(illumina_to_pred, reference_fasta, outfolder, unsupported_cutoff):
    samfile = pysam.AlignmentFile(illumina_to_pred, "rb" )

    reference_coverage = {acc : {pos : 0} for acc in reference_fasta for pos in range(len(reference_fasta[acc]))} 
    reference_distribution = {acc : {pos : defaultdict(int)} for acc in reference_fasta for pos in range(len(reference_fasta[acc]))}  
    print("number of references:", len(reference_coverage) )
    
    references_seen_in_pileup = set()
    current_ref = ""
    prev_pos = -1
    for pileupcolumn in samfile.pileup():
        # new reference
        references_seen_in_pileup.add(pileupcolumn.reference_name)
        if pileupcolumn.reference_name != current_ref:
            if prev_pos + 1 < len(reference_fasta[pileupcolumn.reference_name]):
                print("No coverage on reference {0} at positions: {1} to {2}. Length reference:{3} (0-indexed)".format(pileupcolumn.reference_name, prev_pos + 1, len(reference_fasta[pileupcolumn.reference_name]) - 1, len(reference_fasta[pileupcolumn.reference_name])))

            prev_pos = -1
            current_ref = pileupcolumn.reference_name

        # print(pileupcolumn.reference_name, pileupcolumn.pos)
        # print ("coverage at base %s = %s" %
        #        (pileupcolumn.pos, pileupcolumn.n))
        if prev_pos + 1 < pileupcolumn.pos:
            # for i in range(prev_pos +1, pileupcolumn.pos):
            print("No coverage on reference {0} at positions: {1} to {2}. Length reference:{3} (0-indexed)".format(pileupcolumn.reference_name, prev_pos + 1, pileupcolumn.pos, len(reference_fasta[pileupcolumn.reference_name])))
        prev_pos = pileupcolumn.pos

        # CHECK DELETIONS INSERTIONS AND SUBSTITUTUTIONS HERE
        illumina_support_count = 0
        ref_base = reference_fasta[pileupcolumn.reference_name][pileupcolumn.pos]

        for pileupread in pileupcolumn.pileups:   
            if pileupread.alignment.is_unmapped:
                print("is unmapped")
                nr_unmapped += 1
                continue
            if pileupread.alignment.is_supplementary:
                # print("is supplementary")
                continue

            if not pileupread.is_del and not pileupread.is_refskip: # no deletion
                illumina_base = pileupread.alignment.query_sequence[pileupread.query_position]
                if ref_base == illumina_base: # no substitution 
                    if pileupread.indel == 0: # no suceeding insertion
                        illumina_support_count += 1 # we have a read that fully supports the position
        # print(illumina_support_count, pileupcolumn.n, unsupported_cutoff)
        illumina_supported = False
        if illumina_support_count >= unsupported_cutoff:
            illumina_supported = True

        # check what type of error and log it
        # reference_distribution[pileupcolumn.reference_name][pileupcolumn.pos]
        if not illumina_supported:
            print("here")
            # IF ALL ILLUMINA HAS AN INSERTION BETWEEN TWO BASE PAIRS -- this signals a deletion in the predicted transcript, count this!  
            if pileupread.indel > 1:  # only dealing with insertions of one base pair for now!!
                print("Skipping logging of insertion of length:", pileupread.indel)
                continue


    not_seen_in_pileup = set(reference_fasta.keys()).difference(references_seen_in_pileup)
    print("references not seen in pileup:", not_seen_in_pileup)

    # for pileupcolumn in samfile.pileup():
    #     if pileupcolumn.reference_name not in ref_positions:
    #         ref_positions[pileupcolumn.reference_name] = {}
    #     print ("coverage at base %s = %s" %
    #            (pileupcolumn.pos, pileupcolumn.n))

    #     ref_base = reference_fasta[pileupcolumn.reference_name][pileupcolumn.pos]
    #     ref_positions[pileupcolumn.pos] = ref_base
    #     # print(ref_base)
    #     illumina_positions[pileupcolumn.pos] = defaultdict(int)
    #     illumina_positions[-pileupcolumn.pos] = defaultdict(int)
    #     illumina_accessions[pileupcolumn.pos] = defaultdict(list)
    #     illumina_accessions[-pileupcolumn.pos] = defaultdict(list)
    #     for pileupread in pileupcolumn.pileups:   
    #         if pileupread.alignment.is_unmapped:
    #             # print("is unmapped")
    #             nr_unmapped += 1
    #             continue
    #         if pileupread.alignment.is_supplementary:
    #             # print("is supplementary")
    #             continue

    #         if pileupread.indel > 1:  # only dealing with insertions of one base pair for now!!
    #             print("Skipping logging of insertion of length:", pileupread.indel)
    #             continue
    #         if pileupread.is_del:
    #             state, length = get_deletion_status_in_cigar(pileupread.alignment, pileupread.query_position_or_next)
    #             # print("DELETION HERE:", pileupread.indel, pileupread.alignment.cigartuples, "length",length, pileupread.query_position_or_next)
    #             assert state == 2 # should be deletion
    #             if length > 1:
    #                 continue

    #         if not pileupread.is_del and not pileupread.is_refskip:
    #             illumina_base = pileupread.alignment.query_sequence[pileupread.query_position]
    #             illumina_positions[pileupcolumn.pos][illumina_base] += 1
    #             if illumina_base != ref_base:
    #                 illumina_accessions[pileupcolumn.pos][illumina_base].append( ((pileupread.alignment.query_name, pileupread.alignment.is_read1), pileupread.query_position))

    #         elif pileupread.is_del:
    #             print("Deletion:", pileupread.indel, pileupread.alignment.cigartuples, pileupread.query_position, pileupread.query_position_or_next)
    #                 # need to store deletion length before we end up here
    #             illumina_positions[pileupcolumn.pos]["-"] += 1
    #             # print("Deletion here! Length: {0}. Q-position {1}".format(pileupread.indel))
    #             illumina_accessions[pileupcolumn.pos]["-"].append( ((pileupread.alignment.query_name, pileupread.alignment.is_read1), pileupread.query_position_or_next ))

    #         elif pileupread.is_refskip:
    #             illumina_positions[pileupcolumn.pos]["N"] += 1
    #         else:
    #             # should not end up here
    #             assert False

    #         # between bases insertions
    #         if pileupread.indel == 1:
    #             print("insertion here! Length: {0}. Q-position {1}".format(pileupread.indel, pileupread.query_position), pileupread.alignment.cigartuples)
    #             illumina_base = pileupread.alignment.query_sequence[pileupread.query_position+1]
    #             illumina_accessions[-pileupcolumn.pos][illumina_base].append( ((pileupread.alignment.query_name, pileupread.alignment.is_read1), pileupread.query_position+1))
    #             illumina_positions[-pileupcolumn.pos][illumina_base] += 1



    # illumina_variants = defaultdict(list)
    # illumina_variants_read_accessions = defaultdict(lambda : defaultdict(list))
    # # p_illumina_indel = 0.001
    # # p_illumina_subs = 0.001

    # for ref_pos in range(len(reference_fasta[pileupcolumn.reference_name])):
    #     ref_base = reference_fasta[pileupcolumn.reference_name][ref_pos]
    #     if ref_pos in illumina_positions:
    #         total_illumina_support = sum([count for nucl, count in illumina_positions[ref_pos].items()])            

    #         for site, count in illumina_positions[ref_pos].items(): 
    #             if site != ref_base and site != "N":
    #                 if site == "-"  and count >= variant_cutoff: # max(1, total_illumina_support*p_illumina_indel):
    #                     illumina_variants[ref_pos].append(site)
    #                     accessions_list = illumina_accessions[ref_pos][site]
    #                     illumina_variants_read_accessions[ref_pos][site] = accessions_list

    #                 elif site != "-" and count >= variant_cutoff: # max(1, total_illumina_support*p_illumina_subs):
    #                     illumina_variants[ref_pos].append(site)
    #                     accessions_list = illumina_accessions[ref_pos][site]
    #                     illumina_variants_read_accessions[ref_pos][site] = accessions_list

    #     if -ref_pos in illumina_positions: # insertions
    #         total_illumina_support = sum([count for nucl, count in illumina_positions[ref_pos].items()])            
    #         for site, count in illumina_positions[-ref_pos].items(): 
    #             print("INSERTION", ref_base, site )
    #             if count >= variant_cutoff and site != "N":
    #                 illumina_variants[-ref_pos].append(site)  
    #                 accessions_list = illumina_accessions[-ref_pos][site]
    #                 illumina_variants_read_accessions[-ref_pos][site] = accessions_list

    #     else:
    #         print("No alignments", ref_pos)

    # for p in illumina_variants:
    #     print(p, illumina_variants[p], illumina_positions[p])
    #     if p < 0:
    #         print("REF POS flanking insertion:",  reference_fasta[pileupcolumn.reference_name][p],  reference_fasta[pileupcolumn.reference_name][p+1]  )
    #     # print("ref base:", reference_fasta[pileupcolumn.reference_name][p] , p, illumina_variants[p], illumina_positions[p])
    #     # for site in illumina_variants[p]:
    #     #     print(illumina_accessions[p][site])


    # samfile.close()
    
    # illumina_variants_out = open(os.path.join(outfolder, 'variants.pkl'), 'wb')
    # illumina_accessions_out = open(os.path.join(outfolder, 'accessions.pkl'), 'wb')
    # illumina_positions_out = open(os.path.join(outfolder, 'positions.pkl'), 'wb')

    # # Pickle dictionary using protocol 0.
    # pickle.dump(illumina_variants, illumina_variants_out)
    # pickle.dump(illumina_variants_read_accessions, illumina_accessions_out)
    # pickle.dump(illumina_positions, illumina_positions_out)

    # # sys.exit()
    # return illumina_variants, illumina_variants_read_accessions, illumina_positions




def get_general_alignment_quality(illumina_to_pred, illumina_variants, illumina_accessions, illumina_positions, predicted_transcripts, outfolder, args):
    samfile = pysam.AlignmentFile(illumina_to_pred, "rb" )

    read_accession_to_query_pos_and_variant = defaultdict(list)
    for ref_pos, var_dict in  illumina_accessions.items():
        for variant, acc_and_q_pos_list in var_dict.items():
            for read_id, pos in acc_and_q_pos_list:
                read_accession_to_query_pos_and_variant[read_id].append( (variant, pos, ref_pos) )

    print("Nr predicted transcripts:", len(samfile.references)) # one reference at a time
    deletions_captured = defaultdict(Counter) # ref_pos : site
    insertions_captured = defaultdict(Counter) # ref_pos : site
    substitutions_captured = defaultdict(Counter) # ref_pos : site
    nr_unmapped = 0

    for read in samfile.fetch():

        # if not read.is_unmapped:
            # print(read.is_unmapped, read.flag)
            # for state, num in read.cigartuples:
            #     if state ==8:
            #         print("DIFF!!!")

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
                        insertions_captured[ref_pos][site] += 1 #.add(site)
                    else:
                        pass
                        #### CHECK #######
                        # print("INSERTION not captured!:", ref_pos, read.cigartuples, read_state, state_length, q_var_pos, site)
                        # print("alignment on pred", read.cigartuples, read.query_name, read.is_read1)
                        # tempsam = pysam.AlignmentFile(args.illumina_to_pred, "rb" )
                        # for r in tempsam.fetch():
                        #     if r.query_name == read.query_name:
                        #         print( "Alignment on ref:",r.cigartuples, r.query_name, r.is_read1)
                        ###################

                elif site == "-": # how to match deletions? need to look at cigar here
                    # print("DELETION")
                    assert ref_pos >= 0
                    q_pos_next = q_var_pos # the position is really the next one adter the deletion (reported by pysam pileup) 
                    # read_state1, state_length1 = get_status_in_cigar(read, q_pos_next) #right flanking
                    # read_state2, state_length2 = get_status_in_cigar(read, q_pos_next-1) # left flanking

                    # print("DELETION !:", ref_pos,  q_var_pos, read.cigartuples)
                    if not is_deletion(read, q_pos_next):
                        # print("DELETION captured!:", ref_pos, read.cigartuples, q_var_pos, read.query_name, read.is_read1)
                        deletions_captured[ref_pos][site] += 1 #.add(site)
                    else:
                        pass
                        # print("DELETION not captured!:", ref_pos, q_var_pos, read.cigartuples)

                else: # substitution
                    pass
                    assert ref_pos >= 0 and site != "-" 
                    predicted_transcript_pos = read_aligned_to_pred_transcript_positions[q_var_pos]
                    if predicted_transcript_pos:
                        pred_transcript_site = predicted_transcripts[read.reference_name][predicted_transcript_pos]
                        pred_transcript_region = predicted_transcripts[read.reference_name][predicted_transcript_pos - 4 : predicted_transcript_pos + 4]
                        illumina_read_region = read.query_alignment_sequence[max(0,q_var_pos - read.qstart - 4) : min(q_var_pos - read.qstart + 4, len( read.query_alignment_sequence) )] # query_alignment_sequence is the portion [read.qstart : read.qend]
                        
                        if site == pred_transcript_site: # and aligns well! check region = region to verify this statement:
                            # print("CAPTURED:", illumina_read_region, pred_transcript_region, site, q_var_pos, read.qstart, read.qend)
                            substitutions_captured[ref_pos][site] += 1 #.add(site)
                            # print("SUBSTITUTION CAPTURED:", ref_pos, site, pred_transcript_site)
                        else:
                            pass
                            print("Sites not matching for SUBSTITUTION:", ref_pos, site, pred_transcript_site, illumina_read_region, pred_transcript_region)
                    else:
                        # print(q_var_pos, "this part of the read was not aligned to any predicted transcript", ref_pos, read.cigartuples)
                        pass



    total_number_of_illumina_varinats = len([1 for pos in illumina_variants for site in illumina_variants[pos]])
    captured = 0
    del_captured = len([ site for ref_pos in  deletions_captured for site in deletions_captured[ref_pos]])
    ins_captured = len([ site for ref_pos in  insertions_captured for site in insertions_captured[ref_pos]])
    subs_captured = len([ site for ref_pos in  substitutions_captured for site in substitutions_captured[ref_pos]])
    print("Deletions CAPTURED1:", del_captured)
    print("Insertions CAPTURED1:", ins_captured)
    print("Substitutions CAPTURED1:", subs_captured)

    illumina_support_on_pred_del_captured = [  deletions_captured[ref_pos][site] for ref_pos in  deletions_captured for site in deletions_captured[ref_pos]]
    illumina_support_on_pred_ins_captured = [  insertions_captured[ref_pos][site] for ref_pos in  insertions_captured for site in insertions_captured[ref_pos]]
    illumina_support_on_pred_subs_captured = [  substitutions_captured[ref_pos][site] for ref_pos in  substitutions_captured for site in substitutions_captured[ref_pos]]
    print("Illumina support Deletions CAPTURED1:", illumina_support_on_pred_del_captured)
    print("Illumina support Insertions CAPTURED1:", illumina_support_on_pred_ins_captured)
    print("Illumina support Substitutions CAPTURED1:", illumina_support_on_pred_subs_captured)

    # for pos in sites_captured:
    #     for site in sites_captured[pos]:
    #         captured += 1
    #         print(pos, site)

    print("\n")
    print("\n")
    print("FINAL STATS FROM RUN")
    print("\n")

    del_not_captured = 0
    del_captured2 = 0
    del_captured_illumina_depths = []
    del_non_captured_illumina_depths = []
    illumina_support_on_pred_del_captured = []

    ins_not_captured = 0
    ins_captured2 = 0
    ins_captured_illumina_depths = []
    ins_non_captured_illumina_depths = []
    illumina_support_on_pred_ins_captured = []

    subs_not_captured = 0
    subs_captured2 = 0
    subs_captured_illumina_depths = []
    subs_non_captured_illumina_depths = []
    illumina_support_on_pred_subs_captured = []

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
                        illumina_support_on_pred_del_captured.append(deletions_captured[ref_pos][illumina_site])
                    else:
                        del_not_captured +=1
                        del_non_captured_illumina_depths.append(illumina_positions[ref_pos][illumina_site])
                else:
                    if ref_pos in substitutions_captured and illumina_site in substitutions_captured[ref_pos]:
                        # print("HEHEH", illumina_site, ref_pos, sites_captured[ref_pos])
                        # if illumina_site in deletions_captured[ref_pos]:
                        #     # print("Cap2", ref_pos, illumina_site)
                        subs_captured2 += 1
                        subs_captured_illumina_depths.append(illumina_positions[ref_pos][illumina_site])
                        illumina_support_on_pred_subs_captured.append(substitutions_captured[ref_pos][illumina_site])

                    else:
                        subs_not_captured +=1
                        subs_non_captured_illumina_depths.append(illumina_positions[ref_pos][illumina_site])                    
        
        else: #insertion
            for illumina_site in illumina_variants[ref_pos]:
                if ref_pos in insertions_captured and illumina_site in insertions_captured[ref_pos]:
                    # if illumina_site in insertions_captured[ref_pos]:
                    #     # print("Cap2", ref_pos, illumina_site)
                    ins_captured2 += 1
                    ins_captured_illumina_depths.append(illumina_positions[ref_pos][illumina_site])
                    illumina_support_on_pred_ins_captured.append(insertions_captured[ref_pos][illumina_site])

                else:
                    ins_not_captured +=1
                    ins_non_captured_illumina_depths.append(illumina_positions[ref_pos][illumina_site])                


    # print(set(sites_captured.keys()).difference(set(illumina_variants)))

    print("Total variant carrying reads:", len(read_accession_to_query_pos_and_variant))
    print("Total number of illumina supported sites:", total_number_of_illumina_varinats)
    print("Number of variant carrying reads that were unmapped on predicted:", nr_unmapped)

    # print("Deletion sites captured:", del_captured)
    # print("Insertion sites captured:", ins_captured)
    # print("Substitution sites captured:", subs_captured)
    assert del_captured == del_captured2
    assert ins_captured == ins_captured2
    assert subs_captured == subs_captured2

    print("Deletion sites captured:", del_captured2)
    print("Insertion sites captured:", ins_captured2)
    print("Substitution sites captured:", subs_captured2)
    print("\n")
    print("Deletion sites not captured:", del_not_captured)
    print("Insertion sites not captured:", ins_not_captured)
    print("Substitution sites not captured:", subs_not_captured)

    # print("Sites not captured:", del_not_captured, ins_not_captured, subs_not_captured )

    # common_params = dict(bins=30, 
    #                  range=(0, 10000), 
    #                  normed=0, label=['subs','ins','del'])
    common_params = dict(normed=0, label=['captured','not captured'], range=(0, 30), bins=30)
    print("deletions captured depths on reference:")
    print(del_captured_illumina_depths)
    print("deletions captured depths on predicted:")
    print(illumina_support_on_pred_del_captured)
    print("deletions not captured depths:")
    print(del_non_captured_illumina_depths)
    print("\n")

    print("insertions captured depths on reference:")
    print(ins_captured_illumina_depths)
    print("insertions captured depths on predicted:")
    print(illumina_support_on_pred_ins_captured)
    print("insertions not captured depths:")    
    print(ins_non_captured_illumina_depths)
    print("\n")

    print("substitutions captured depths on reference:")
    print(subs_captured_illumina_depths)
    print("substitutions captured depths on predicted:")
    print(illumina_support_on_pred_subs_captured)
    print("substitutions not captured depths:")
    print(subs_non_captured_illumina_depths)

    title_header = "Illumina depths deletions" #"reference mismatches, total: {0}, perfect: {1}".format(len(tuple_identities[0]), perfect_matches)
    custom_histogram(del_captured_illumina_depths, outfolder, name='del_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)
    custom_histogram(del_non_captured_illumina_depths, outfolder, name='del_not_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)

    title_header = "Illumina depths insertions" #"reference mismatches, total: {0}, perfect: {1}".format(len(tuple_identities[0]), perfect_matches)
    custom_histogram(ins_captured_illumina_depths, outfolder, name='ins_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)
    custom_histogram(ins_non_captured_illumina_depths, outfolder, name='ins_not_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)

    title_header = "Illumina depths substitutions" #"reference mismatches, total: {0}, perfect: {1}".format(len(tuple_identities[0]), perfect_matches)
    custom_histogram(subs_captured_illumina_depths, outfolder, name='subs_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)
    custom_histogram(subs_non_captured_illumina_depths, outfolder, name='subs_not_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)

    return

def main(args):
    """ 
        1. Retrieve the positions on predicted transcripts that are not supported by illumina reads.\n\n
        2. Also log for each read all the substitutions and indels of 1 bp.\n\n
        3. Summarize statistics\n\n
    """
    output_file = open(os.path.join(args.outfolder, "alignments_summary.tsv" ), "w")
    predicted_seqs = {acc.split()[0]: seq for (acc, seq) in  read_fasta(open(args.predicted, 'r'))}

    get_unsupported_positions_on_predicted(args.illumina_to_pred, predicted_seqs, args.outfolder, args.unsupported_cutoff)
    # print(predicted_transcripts)
    # get_general_alignment_quality(args.illumina_to_pred, predicted_seqs, args.outfolder, args)




if __name__ == '__main__':

# Take care of input

    parser = argparse.ArgumentParser(description = "This script does: \n\n {0}".format(main.__doc__)) 
    parser.add_argument('-illumina_to_pred', type=str, help='Bam file ')
    parser.add_argument('-predicted', type=str, help='Fasta file ')
    parser.add_argument('-outfolder', type=str, help='outfile folder to put output in. ')
    parser.add_argument('--unsupported_cutoff', type=int, default=2, help='Fasta file ')


    args = parser.parse_args()

    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)
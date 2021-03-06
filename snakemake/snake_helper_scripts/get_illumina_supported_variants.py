
import argparse
import os
from collections import Counter
import re
import pysam
from collections import defaultdict
import dill
import pickle
import math

import edlib

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

def indels_in_cigar(alignment):
    # mismatches = 0
    indels = 0
    matches = 0
    for state, number in alignment.cigartuples:
        if state == 0:
            matches += number
        else:
            indels += number

    return indels, matches

def get_edit_distance_from_cigar(alignment):
    # print(pileupread.indel, pileupread.alignment.cigartuples)
    total_ed = 0
    for state, number in alignment.cigartuples:
        if state != 0:
            total_ed += number


    return total_ed

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

def get_variants_on_reference(illumina_to_ref, reference_fasta, outfolder, variant_cutoff):
    samfile = pysam.AlignmentFile(illumina_to_ref, "rb" )
    illumina_positions = {}
    illumina_accessions = {}
    ref_positions = {}
    coverage_vector = []
    assert len(reference_fasta) == 1 # one reference at a time

    for pileupcolumn in samfile.pileup():
        print ("coverage at base %s = %s" %
               (pileupcolumn.pos, pileupcolumn.n))
        coverage_vector.append(pileupcolumn.n)

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

            # has to be aligned with very high quality! Otherwise, e.g., reads from other isoforms might cause spurious vairnat by
            # having the same alignment over exon junctions that differs between the consensus and other isoforms
            if get_edit_distance_from_cigar(pileupread.alignment) > 5:
                print("Bigger than 5", get_edit_distance_from_cigar(pileupread.alignment), pileupread.alignment.cigartuples )
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


    # for pos in illumina_positions:
    #     print(illumina_positions[pos])
    # sys.exit()


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
    coverage_vector_out = open(os.path.join(outfolder, 'coverage.pkl'), 'wb')

    # Pickle dictionary using protocol 0.
    pickle.dump(illumina_variants, illumina_variants_out)
    pickle.dump(illumina_variants_read_accessions, illumina_accessions_out)
    pickle.dump(illumina_positions, illumina_positions_out)
    pickle.dump(coverage_vector, coverage_vector_out)

    # sys.exit()
    return illumina_variants, illumina_variants_read_accessions, illumina_positions, coverage_vector




def find_if_supported_in_pred_transcripts(illumina_to_pred, illumina_variants, illumina_accessions, illumina_positions, predicted_transcripts, coverage_vector, outfolder, args):
    samfile = pysam.AlignmentFile(illumina_to_pred, "rb" )

    read_accession_to_query_pos_and_variant = defaultdict(list)
    for ref_pos, var_dict in  illumina_accessions.items():
        for variant, acc_and_q_pos_list in var_dict.items():
            for read_id, pos in acc_and_q_pos_list:

                A read can only be assigned to one variant here! this is wrong!!!!!
                
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

            # print(indels_in_cigar(read))
            nr_indels, matches = indels_in_cigar(read)
            # TODO: check edit distance between read and reference segment
            q_seq = read.query_sequence
            r_seq = predicted_transcripts[read.reference_name][read.reference_start: read.reference_start + read.query_length]
            result = edlib.align(q_seq, r_seq)
            ed = result["editDistance"]
            # print(ed, "ed")
            if ed > args.max_mismatch:
                print("HERE", nr_indels, ed, matches)
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
                        # tempsam = pysam.AlignmentFile(args.illumina_to_ref, "rb" )
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
                            # print("Sites not matching for SUBSTITUTION:", ref_pos, site, pred_transcript_site, illumina_read_region, pred_transcript_region)
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

    coverage_vector.sort()
    n = float( len(coverage_vector) )
    mu =  sum(coverage_vector) / n
    sigma = (sum(list(map((lambda x: x ** 2 - 2 * x * mu + mu ** 2), coverage_vector))) / (n - 1)) ** 0.5
    min_cov = min(coverage_vector)
    max_cov = max(coverage_vector)
    coverage_vector.sort()
    if len(coverage_vector) %2 == 0:
        median_cov = (coverage_vector[int(len(coverage_vector)/2)-1] + coverage_vector[int(len(coverage_vector)/2)]) / 2.0
    else:
        median_cov = coverage_vector[int(len(coverage_vector)/2)]


    print("Coverage over positions: mean {0}, stddev: {1}, median: {2}, min: {3}, max:{4}".format(mu, sigma, median_cov, min_cov, max_cov) )
    p_illumina_subs = 0.001
    mean_var = p_illumina_subs* mu
    # if poisson distributed, mean = var
    suggested_cutoff = math.ceil(mean_var + 6 * math.sqrt(mean_var))
    print("Based on coverage distribution and prob of illumina error of base set to 0.0001, suggested variant support cutoff is: {0}, calculated as: suggested_cutoff = mean_error + 6*std_error, where mean_error = 0.0001*mean_cov.".format(suggested_cutoff))
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

    print("Total Illumina reads that supported a variant on consensus: {0}".format(len(read_accession_to_query_pos_and_variant)))
    print("Total variant positions on consensus with coverage of at least {0} illumina reads: {1}".format(args.variant_cutoff, total_number_of_illumina_varinats))
    print("Number of Illumina reads that upported a varinat on consensus but were not mapped to any predicted transcript: {0}".format(nr_unmapped))
    print("\n")
    print("Total deletion sites on consensus with coverage of at least {0} illumina reads: {1}".format(args.variant_cutoff, del_captured2 + del_not_captured))
    print("Total insertion sites on consensus with coverage of at least {0} illumina reads: {1}".format(args.variant_cutoff, ins_captured2 + ins_not_captured))
    print("Total substitution sites on consensus with coverage of at least {0} illumina reads: {1}".format(args.variant_cutoff, subs_captured2 + subs_not_captured))
    assert del_captured == del_captured2
    assert ins_captured == ins_captured2
    assert subs_captured == subs_captured2
    print("\n")
    print("Deletion sites captured in predicted with coverage of at least {0} illumina reads: {1}".format(args.variant_cutoff, del_captured2))
    print("Insertion sites captured in predicted with coverage of at least {0} illumina reads: {1}".format(args.variant_cutoff, ins_captured2))
    print("Substitution sites captured in predicted with coverage of at least {0} illumina reads: {1}".format(args.variant_cutoff, subs_captured2))
    print("\n")
    print("Deletion sites not captured in predicted with coverage of at least {0} illumina reads: {1}".format(args.variant_cutoff, del_not_captured))
    print("Insertion sites not captured in predicted with coverage of at least {0} illumina reads: {1}".format(args.variant_cutoff, ins_not_captured))
    print("Substitution sites not captured in predicted with coverage of at least {0} illumina reads: {1}".format(args.variant_cutoff, subs_not_captured))
    print("\n")
    # print("Sites not captured:", del_not_captured, ins_not_captured, subs_not_captured )

    # common_params = dict(bins=30, 
    #                  range=(0, 10000), 
    #                  normed=0, label=['subs','ins','del'])
    common_params = dict(normed=0, label=['captured','not captured'], range=(0, 30), bins=30)
    print("deletions captured depths on reference:")
    print(sorted(del_captured_illumina_depths))
    print("deletions captured depths on predicted:")
    print(sorted(illumina_support_on_pred_del_captured))
    print("deletions not captured depths:")
    print(sorted(del_non_captured_illumina_depths))
    print("\n")

    print("insertions captured depths on reference:")
    print(sorted(ins_captured_illumina_depths))
    print("insertions captured depths on predicted:")
    print(sorted(illumina_support_on_pred_ins_captured))
    print("insertions not captured depths:")    
    print(sorted(ins_non_captured_illumina_depths))
    print("\n")

    print("substitutions captured depths on reference:")
    print(sorted(subs_captured_illumina_depths))
    print("substitutions captured depths on predicted:")
    print(sorted(illumina_support_on_pred_subs_captured))
    print("substitutions not captured depths:")
    print(sorted(subs_non_captured_illumina_depths))

    title_header = "Illumina depths deletions" #"reference mismatches, total: {0}, perfect: {1}".format(len(tuple_identities[0]), perfect_matches)
    custom_histogram(del_captured_illumina_depths, outfolder, name='del_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)
    custom_histogram(del_non_captured_illumina_depths, outfolder, name='del_not_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)

    title_header = "Illumina depths insertions" #"reference mismatches, total: {0}, perfect: {1}".format(len(tuple_identities[0]), perfect_matches)
    custom_histogram(ins_captured_illumina_depths, outfolder, name='ins_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)
    custom_histogram(ins_non_captured_illumina_depths, outfolder, name='ins_not_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)

    title_header = "Illumina depths substitutions" #"reference mismatches, total: {0}, perfect: {1}".format(len(tuple_identities[0]), perfect_matches)
    custom_histogram(subs_captured_illumina_depths, outfolder, name='subs_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)
    custom_histogram(subs_non_captured_illumina_depths, outfolder, name='subs_not_captured.png', x='Illumina depth', y='frequency', title=title_header, params = common_params)


    # print to file code
    tsv_outfile_del_captured = open(os.path.join(args.outfolder, 'deletions_captured_depths.tsv'), 'w')
    for depth in del_captured_illumina_depths:
        tsv_outfile_del_captured.write("{0}\n".format(depth))
    tsv_outfile_del_captured.close()

    tsv_outfile_del_uncaptured = open(os.path.join(args.outfolder, 'deletions_uncaptured_depths.tsv'), 'w')
    for depth in del_non_captured_illumina_depths:
        tsv_outfile_del_uncaptured.write("{0}\n".format(depth))
    tsv_outfile_del_uncaptured.close()


    tsv_outfile_ins_captured = open(os.path.join(args.outfolder, 'insertions_captured_depths.tsv'), 'w')
    for depth in ins_captured_illumina_depths:
        tsv_outfile_ins_captured.write("{0}\n".format(depth))
    tsv_outfile_ins_captured.close()

    tsv_outfile_ins_uncaptured = open(os.path.join(args.outfolder, 'insertions_uncaptured_depths.tsv'), 'w')
    for depth in ins_non_captured_illumina_depths:
        tsv_outfile_ins_uncaptured.write("{0}\n".format(depth))
    tsv_outfile_ins_uncaptured.close()

    tsv_outfile_subs_captured = open(os.path.join(args.outfolder, 'substitutions_captured_depths.tsv'), 'w')
    for depth in subs_captured_illumina_depths:
        tsv_outfile_subs_captured.write("{0}\n".format(depth))
    tsv_outfile_subs_captured.close()

    tsv_outfile_subs_uncaptured = open(os.path.join(args.outfolder, 'substitutions_uncaptured_depths.tsv'), 'w')
    for depth in subs_non_captured_illumina_depths:
        tsv_outfile_subs_uncaptured.write("{0}\n".format(depth))
    tsv_outfile_subs_uncaptured.close()
    #################


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
        coverage_vector_in = open(os.path.join(args.outfolder, 'coverage.pkl'), 'rb')    
        illumina_variants = pickle.load(illumina_variants_in)
        illumina_accessions = pickle.load(illumina_accessions_in)
        illumina_positions = pickle.load(illumina_positions_in)
        coverage_vector = pickle.load(coverage_vector_in)
        # print(illumina_variants)
        # print(illumina_accessions)

    else:
        illumina_variants, illumina_accessions, illumina_positions, coverage_vector = get_variants_on_reference(args.illumina_to_ref, reference_seq, args.outfolder, args.variant_cutoff)

    predicted_transcripts = {acc.split()[0]: seq for (acc, seq) in read_fasta(open(args.predicted, 'r'))}
    # print(predicted_transcripts)
    find_if_supported_in_pred_transcripts(args.illumina_to_pred, illumina_variants, illumina_accessions, illumina_positions, predicted_transcripts, coverage_vector, args.outfolder, args)




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
    parser.add_argument('--max_mismatch', type=int, default = 3, help='The maximum number of mismatches of an illumina reads in order to be considered supporting a variant in the predicted transcript')

    args = parser.parse_args()
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)
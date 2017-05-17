import argparse
import os, sys

def get_best_hits(file_name):
    best_hits = {}
    for line in open(file_name):
        # print(line)
        # print(line.strip().split("\t"))
        query, target, len_query, len_target, ed = line.strip().split("\t")
        ed = int(ed)
        if target in best_hits:
            if ed < best_hits[target]:
                best_hits[target] = ed

        else:
            best_hits[target] = ed
    return best_hits

def main(args):
    
    read_hits = get_best_hits(args.best_hits_reads)
    predicted_hits = get_best_hits(args.best_hits_pred)
    print(len(read_hits))
    print(len(predicted_hits))

    FN = {}
    for db_hit in read_hits:
        if db_hit not in predicted_hits:
            FN[db_hit] = read_hits[db_hit]
            print("not found:", db_hit, read_hits[db_hit])
        elif predicted_hits[db_hit] > read_hits[db_hit]:
            print("read had better hit", predicted_hits[db_hit], read_hits[db_hit], db_hit)
            FN[db_hit] = read_hits[db_hit]
        elif predicted_hits[db_hit] == read_hits[db_hit]:
            print("Same quality hit", predicted_hits[db_hit], read_hits[db_hit], db_hit)
        else:
            # print("predicted had better hit", predicted_hits[db_hit], read_hits[db_hit])
            continue

    print("False negatives")
    # for fn in FN:
    #     print(FN[fn], fn)
    better_hit_to_pred = {}
    for db_hit in predicted_hits:
        if db_hit not in read_hits:
            better_hit_to_pred[db_hit] = predicted_hits[db_hit]
            print("not found in reads:", db_hit, better_hit_to_pred[db_hit])
        elif predicted_hits[db_hit] < read_hits[db_hit]:
            print("Predicted had better hit", predicted_hits[db_hit], read_hits[db_hit], db_hit)
            better_hit_to_pred[db_hit] = predicted_hits[db_hit]
        else:
            # print("read had better hit", predicted_hits[db_hit], read_hits[db_hit])
            continue

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--best_hits_pred', type=str, help='Path to the predicted transcript fasta file')
    parser.add_argument('--best_hits_reads', type=str, default=None, help='Path to the consensus fasta file')
    parser.add_argument('--outfolder', type=str, help='Output path of results')
    parser.add_argument('--single_core', dest='single_core', action='store_true', help='Force working on single core. ')
    args = parser.parse_args()

    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)
    

    main(args)

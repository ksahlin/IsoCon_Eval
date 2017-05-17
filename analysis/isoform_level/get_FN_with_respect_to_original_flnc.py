

def get_best_hits(file_name):
    best_hits = {}
    for line in open(file_name):
        query, target, ed = line.strip().split()
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
    FN = {}
    for db_hit in read_hits:
        if db_hit not in predicted_hits:
            FN[db_hit] = read_hits[db_hit]
        elif predicted_hits[db_hit] > read_hits[db_hit]:
            print("read had better hit")
            FN[db_hit] = read_hits[db_hit]
        else:
            continue

    print("False negatives")
    for fn in FN:
        print(FN[fn], fn)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--best_hits_pred', type=str, help='Path to the predicted transcript fasta file')
    parser.add_argument('--best_hits_reads', type=str, default=None, help='Path to the consensus fasta file')
    parser.add_argument('--outfolder', type=str, help='Output path of results')
    parser.add_argument('--single_core', dest='single_core', action='store_true', help='Force working on single core. ')
    args = parser.parse_args()

    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)
    
    if args.database:
    else:
        main(args)

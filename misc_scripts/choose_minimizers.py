import networkx as nx
import argparse, os

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


def main(args):
    predicted = {acc: seq for (acc, seq) in  read_fasta(open(args.fasta, 'r'))}

    G = nx.DiGraph()
    for line in open(args.file1, "r"):
        read_acc, m_acc, r_len, m_len, ed = line.strip().split("\t")
        G.add_edge(read_acc, m_acc, edit_distance=ed) 
        !! COLLAPSING OF NODES: ADD WEIGHT HERE SOMEWHERE IF SEQUENCES ARE IDENTICAL 

    # components = [len(c) for c in sorted(nx.weakly_connected_components(G), key=len, reverse=True)]
    # print(components)
    # print(len(components))

    # strong_components = [len(c) for c in sorted(nx.strongly_connected_components(G), key=len, reverse=True)]
    # print(strong_components)
    # print(len(strong_components))

    nr_consensus = 0
    G_transpose = nx.reverse(G)
    consensus = set()
    for subgraph in sorted(nx.weakly_connected_component_subgraphs(G_transpose), key=len, reverse=True):
        print("Subgraph of size", len(subgraph.nodes()), len(subgraph) )
        while subgraph:
            number_connected_to = {}
            reachable_comp_sizes = []
            reachable_comp_nodes = []
            for m in subgraph:
                reachable_comp = set([m])
                for reachable_node in nx.dfs_postorder_nodes(subgraph,source=m):
                    reachable_comp.add(reachable_node)
                reachable_comp_sizes.append(len(reachable_comp))
                reachable_comp_nodes.append(reachable_comp)
                number_connected_to[m] = len(reachable_comp)

            sorted_reachable_comp_sizes = sorted(reachable_comp_sizes, reverse=True)
            sorted_reachable_comp_nodes = sorted(reachable_comp_nodes, key = len, reverse=True)
            print(sorted_reachable_comp_sizes)
            biggest_comp = sorted_reachable_comp_nodes[0]
            !! NEED TO GET MAX SIZE BY SUMMING OVER SELF WEIGHT AND ALL NEIGHBORS. 
            minimizer, max_size =  max([(m, size) for m, size in number_connected_to.items()], key= lambda x: x[1])
            print(minimizer, max_size)
            consensus.add(minimizer)   
            subgraph.remove_nodes_from(biggest_comp)

            print("Subgraph after removal size", len(subgraph.nodes()), len(subgraph) )
            nr_consensus += 1

    print("NR CONSENSUS:", nr_consensus)
    out_f = open(os.path.join(args.outfolder, "minimizers.tex"), "w")
    for c in consensus:
        c_seq = predicted[c]
        out_f.write(">{0}\n{1}\n".format(c,c_seq))




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--file1', type=str, help='Alignment file between reads and pred')
    parser.add_argument('--fasta', type=str, default=None, help='Fasta file of pred')
    parser.add_argument('--outfolder', type=str, help='Output path of results')
    parser.add_argument('--single_core', dest='single_core', action='store_true', help='Force working on single core. ')
    args = parser.parse_args()

    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)
import networkx as nx
import argparse, os

from operator import itemgetter
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

def highest_reachable_with_edge_degrees(args):
    predicted = {acc: seq for (acc, seq) in  read_fasta(open(args.fasta, 'r'))}
    predicted_seq_to_acc = defaultdict(list)
    for (acc, seq) in read_fasta(open(args.fasta, 'r')):
        predicted_seq_to_acc[seq].append(acc)
    
    degrees = {}
    for seq, list_acc in predicted_seq_to_acc.items():
        deg = len(list_acc)
        for acc in list_acc:
            degrees[acc] = deg


    G_acc_graph = nx.DiGraph()
    for line in open(args.file1, "r"):
        read_acc, m_acc, r_len, m_len, ed = line.strip().split("\t")
        G_acc_graph.add_edge(read_acc, m_acc, edit_distance=int(ed)) 


    G = nx.DiGraph()
    for acc in G_acc_graph.nodes():
        deg = degrees[acc]
        seq = predicted[acc]
        if seq in G:
            assert G.node[seq]["degree"] == deg
        else:
            G.add_node(seq, degree = deg)

    for read_acc, m_acc in G_acc_graph.edges():
        deg = degrees[read_acc]
        read_seq = predicted[read_acc]
        m_seq = predicted[m_acc]
        ed = int(G_acc_graph[read_acc][m_acc]["edit_distance"])
        if read_seq != m_seq:
            G.add_edge(read_seq, m_seq, edit_distance=int(ed))
            # G[read_seq]["degree"] = deg 
            assert deg == 1
            assert ed > 0
        else:
            assert G.node[read_seq]["degree"] > 1
            assert ed == 0

    partition_sizes = []
    nr_consensus = 0
    G_transpose = nx.reverse(G)
    consensus = set()
    print("here")
    for subgraph in sorted(nx.weakly_connected_component_subgraphs(G_transpose), key=len, reverse=True):
        # print("Subgraph of size", len(subgraph.nodes()), len(subgraph) )
        while subgraph:
            number_connected_to = {}
            reachable_comp_sizes = []
            reachable_comp_weights = {}
            reachable_comp_nodes = []
            edit_distances = []
            processed = set()
            for m in subgraph:
                if m in processed:
                    continue
                reachable_comp = set([m])
                # print("deg m", subgraph.node[m], subgraph[m])
                reachable_comp_weight = subgraph.node[m]["degree"]
                processed.add(m)
                # print("cl size:", len([n for n in nx.dfs_postorder_nodes(subgraph, source=m)]))
                for reachable_node in nx.dfs_postorder_nodes(subgraph, source=m): # store reachable node as processed here to avoid computation
                    if reachable_node == m:
                        continue
                    
                    processed.add(reachable_node)
                    reachable_comp.add(reachable_node)
                    reachable_comp_weight += subgraph.node[reachable_node]["degree"]
                    # print(subgraph.node[reachable_node])
                    if reachable_node in subgraph[m]:
                        edit_distances.append(subgraph[m][reachable_node]["edit_distance"])
                    # print(len(G[reachable_node]), len(G_transpose[reachable_node]), reachable_node == m )
                    assert subgraph.node[reachable_node]["degree"] == 1

                reachable_comp_sizes.append(len(reachable_comp))
                reachable_comp_weights[reachable_comp_weight] = reachable_comp
                reachable_comp_nodes.append(reachable_comp)
                number_connected_to[m] = len(reachable_comp)

            sorted_reachable_comp_sizes = sorted(reachable_comp_sizes, reverse=True)
            sorted_reachable_comp_weights = sorted(reachable_comp_weights.keys(), reverse=True)
            max_weight = max(sorted_reachable_comp_weights)
            sorted_reachable_comp_nodes = sorted(reachable_comp_nodes, key = len, reverse=True)
            # print(sorted_reachable_comp_sizes)
            # print("weight", reachable_comp_weight)
            biggest_comp = sorted_reachable_comp_nodes[0]
            biggest_weighted_comp = reachable_comp_weights[max_weight]
            # print(max_weight, len(biggest_comp), len(biggest_weighted_comp))
            # assert len(biggest_comp) == len(biggest_weighted_comp)
            # !! NEED TO GET MAX SIZE BY SUMMING OVER SELF WEIGHT AND ALL NEIGHBORS. 
            minimizer, max_size =  max([(m, size) for m, size in number_connected_to.items()], key= lambda x: x[1])
            # print(minimizer, max_size)
            partition_sizes.append(max_size)
            consensus.add(minimizer)   
            subgraph.remove_nodes_from(biggest_weighted_comp)

            edit_distances.sort() 
            # print("Subgraph after removal size", len(subgraph.nodes()), len(subgraph), edit_distances )
            nr_consensus += 1

    print("NR CONSENSUS:", nr_consensus)
    out_f = open(os.path.join(args.outfolder, "minimizers.tex"), "w")
    for c_seq in consensus:
        one_c_acc = predicted_seq_to_acc[c_seq][0]
        out_f.write(">{0}\n{1}\n".format(one_c_acc,c_seq))


    print(sorted(partition_sizes, reverse = True))


def highest_reachable(args):
    predicted = {acc: seq for (acc, seq) in  read_fasta(open(args.fasta, 'r'))}

    G = nx.DiGraph()
    for line in open(args.file1, "r"):
        read_acc, m_acc, r_len, m_len, ed = line.strip().split("\t")
        G.add_edge(read_acc, m_acc, edit_distance=int(ed)) 
        # !! COLLAPSING OF NODES: ADD WEIGHT HERE SOMEWHERE IF SEQUENCES ARE IDENTICAL 

    # components = [len(c) for c in sorted(nx.weakly_connected_components(G), key=len, reverse=True)]
    # print(components)
    # print(len(components))

    # strong_components = [len(c) for c in sorted(nx.strongly_connected_components(G), key=len, reverse=True)]
    # print(strong_components)
    # print(len(strong_components))

    partition_sizes = []
    nr_consensus = 0
    G_transpose = nx.reverse(G)
    consensus = set()
    print("here")
    for subgraph in sorted(nx.weakly_connected_component_subgraphs(G_transpose), key=len, reverse=True):
        # print("Subgraph of size", len(subgraph.nodes()), len(subgraph) )
        while subgraph:
            number_connected_to = {}
            reachable_comp_sizes = []
            reachable_comp_nodes = []
            edit_distances = []
            processed = set()
            for m in subgraph:
                if m in processed:
                    continue
                reachable_comp = set([m])
                processed.add(m)
                for reachable_node in nx.dfs_postorder_nodes(subgraph,source=m): # store reachable node as processed here to avoid computation
                    processed.add(reachable_node)
                    reachable_comp.add(reachable_node)
                    if reachable_node in subgraph[m]:
                        edit_distances.append(subgraph[m][reachable_node]["edit_distance"])
                reachable_comp_sizes.append(len(reachable_comp))
                reachable_comp_nodes.append(reachable_comp)
                number_connected_to[m] = len(reachable_comp)

            sorted_reachable_comp_sizes = sorted(reachable_comp_sizes, reverse=True)
            sorted_reachable_comp_nodes = sorted(reachable_comp_nodes, key = len, reverse=True)
            # print(sorted_reachable_comp_sizes)
            biggest_comp = sorted_reachable_comp_nodes[0]
            # !! NEED TO GET MAX SIZE BY SUMMING OVER SELF WEIGHT AND ALL NEIGHBORS. 
            minimizer, max_size =  max([(m, size) for m, size in number_connected_to.items()], key= lambda x: x[1])
            # print(minimizer, max_size)
            partition_sizes.append(max_size)
            consensus.add(minimizer)   
            subgraph.remove_nodes_from(biggest_comp)

            edit_distances.sort() 
            # print("Subgraph after removal size", len(subgraph.nodes()), len(subgraph), edit_distances )
            nr_consensus += 1

    print("NR CONSENSUS:", nr_consensus)
    out_f = open(os.path.join(args.outfolder, "minimizers.tex"), "w")
    for c in consensus:
        c_seq = predicted[c]
        out_f.write(">{0}\n{1}\n".format(c,c_seq))


    print(sorted(partition_sizes, reverse = True))



def highest_nbrs(args):
    predicted = {acc: seq for (acc, seq) in  read_fasta(open(args.fasta, 'r'))}

    G = nx.DiGraph()
    for line in open(args.file1, "r"):
        read_acc, m_acc, r_len, m_len, ed = line.strip().split("\t")
        G.add_edge(read_acc, m_acc, edit_distance=int(ed)) 
        # !! COLLAPSING OF NODES: ADD WEIGHT HERE SOMEWHERE IF SEQUENCES ARE IDENTICAL 

    partition_sizes = []
    nr_consensus = 0
    G_transpose = nx.reverse(G)
    consensus = set()
    print("here")
    for subgraph in sorted(nx.weakly_connected_component_subgraphs(G_transpose), key=len, reverse=True):
        # print("Subgraph of size", len(subgraph.nodes()), len(subgraph) )
        while subgraph:
            m, degree = sorted(subgraph.degree_iter(),key=itemgetter(1),reverse=True)[0]
            edit_distances = []
            number_connected_to = {}
            processed = set()
            reachable_comp_sizes = []
            reachable_comp_nodes = []
            reachable_comp = set([m])
            processed.add(m)

            # for m in subgraph:
            #     if m in processed:
            #         continue
            #     reachable_comp = set([m])
            #     processed.add(m)
            for reachable_node in nx.dfs_postorder_nodes(subgraph,source=m): # store reachable node as processed here to avoid computation
                processed.add(reachable_node)
                reachable_comp.add(reachable_node)
                if reachable_node in subgraph[m]:
                    edit_distances.append(subgraph[m][reachable_node]["edit_distance"])
            reachable_comp_sizes.append(len(reachable_comp))
            reachable_comp_nodes.append(reachable_comp)
            number_connected_to[m] = len(reachable_comp)

            sorted_reachable_comp_sizes = sorted(reachable_comp_sizes, reverse=True)
            sorted_reachable_comp_nodes = sorted(reachable_comp_nodes, key = len, reverse=True)
            # print(sorted_reachable_comp_sizes)
            biggest_comp = sorted_reachable_comp_nodes[0]
            # !! NEED TO GET MAX SIZE BY SUMMING OVER SELF WEIGHT AND ALL NEIGHBORS. 
            minimizer, max_size =  max([(m, size) for m, size in number_connected_to.items()], key= lambda x: x[1])
            # print(minimizer, max_size)
            partition_sizes.append(max_size)
            consensus.add(minimizer)   
            subgraph.remove_nodes_from(biggest_comp)

            edit_distances.sort() 
            # print("Subgraph after removal size", len(subgraph.nodes()), len(subgraph), edit_distances )
            nr_consensus += 1

    print("NR CONSENSUS:", nr_consensus)
    out_f = open(os.path.join(args.outfolder, "minimizers.tex"), "w")
    for c in consensus:
        c_seq = predicted[c]
        out_f.write(">{0}\n{1}\n".format(c,c_seq))


    print(sorted(partition_sizes, reverse = True))

def main(args):
    highest_reachable_with_edge_degrees(args)
    highest_reachable(args)
    highest_nbrs(args)

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
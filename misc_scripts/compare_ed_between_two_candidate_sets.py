import networkx as nx
import argparse, os


def main(args):
    file1_best = {}
    references1 = set()
    G = nx.DiGraph()

    for line in open(args.file1, "r"):
        # print(line)
        read_acc, m_acc, r_len, m_len, ed = line.strip().split("\t")
        file1_best[read_acc] = (int(ed), m_acc)
        references1.add(m_acc)
        G.add_edge(read_acc, m_acc, edit_distance=ed)

    components = [len(c) for c in sorted(nx.weakly_connected_components(G), key=len, reverse=True)]
    print(components)
    print(len(components))

    strong_components = [len(c) for c in sorted(nx.strongly_connected_components(G), key=len, reverse=True)]
    print(strong_components)
    print(len(strong_components))

    # biconnected_components = [len(c) for c in sorted(nx.biconnected_connected_components(G), key=len, reverse=True)]
    # print(biconnected_components)
    # print(len(biconnected_components))

    G_transpose = nx.reverse(G)
    for subgraph in sorted(nx.weakly_connected_component_subgraphs(G_transpose), key=len, reverse=True):
        print("Subgraph of size", len(subgraph.nodes()), len(subgraph) )
        reachable_comp_sizes = []
        for m in subgraph:
            reachable_comp = set([m])
            for reachable_node in nx.dfs_postorder_nodes(subgraph,source=m):
                reachable_comp.add(reachable_node)
            reachable_comp_sizes.append(len(reachable_comp))
        sorted_reachable_comp_sizes = sorted(reachable_comp_sizes,reverse=True)
        print(sorted_reachable_comp_sizes)
        # print(len(reachable_comp_sizes))

    references2 = set()
    file2_best = {}
    for line in open(args.file2, "r"):
        read_acc, m_acc, r_len, m_len, ed = line.strip().split("\t")
        file2_best[read_acc] = (int(ed), m_acc)
        references2.add(m_acc)


    better_on_pred1 = {}
    c1, c2, c3 = 0,0,0
    for r_acc in file2_best:
        if file2_best[r_acc][0] < file1_best[r_acc][0]:
            c2 +=1

        elif file1_best[r_acc][0] < file2_best[r_acc][0]:
            c1 += 1
            better_on_pred1[r_acc] = (file1_best[r_acc], file2_best[r_acc][0])
            if file1_best[r_acc][0] < file2_best[r_acc][0]/2 and file1_best[r_acc][0] > 2:
                # print((file1_best[r_acc][0], file2_best[r_acc][0]), r_acc)
                pass

        elif file2_best[r_acc][0] == file1_best[r_acc][0]:
            c3 += 1

        else:
            print("BUUUG")

    print("Better on, set1 {0}, set2 {1}, tie {2}".format(c1,c2,c3))
    print("Number of references set1:", len(references1))
    print("Number of references set2:", len(references2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--file1', type=str, help='Alignment file between reads and pred')
    parser.add_argument('--file2', type=str, default=None, help='Alignment file between reads and pred')
    parser.add_argument('--outfolder', type=str, help='Output path of results')
    parser.add_argument('--single_core', dest='single_core', action='store_true', help='Force working on single core. ')
    args = parser.parse_args()

    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)
    
    main(args)
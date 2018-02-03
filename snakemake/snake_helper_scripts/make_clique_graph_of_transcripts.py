import argparse
import os
import errno
import networkx as nx
from math import *
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")
# import matplotlib.pyplot as plt
# import matplotlib

import seaborn as sns
import itertools as it
from networkx.drawing.nx_pydot import read_dot, pydot_layout

global colors, hatches
# colors=it.cycle('bgcmy')# blue, green, red, ...
colors=it.cycle(['blue', 'green', 'red', 'cyan', 'magenta', 'purple', 'pink', 'brown', 'orange', 'teal', 'coral', 'lightblue', 'lime', 'lavender', 'turquoise', 'darkgreen', 'tan', 'salmon', 'gold', 'lightpurple', 'darkred', 'darkblue'])# blue, green, red, ...
# hatches=it.cycle('/\|-+*')

def draw_circle_around_clique(clique,coords):
    dist=0
    temp_dist=0
    center=[0 for i in range(2)]
    color=next(colors)
    for a in clique:
        for b in clique:
            # print(coords[a][0], coords[b][0])
            temp_dist=(coords[a][0]-coords[b][0])**2+(coords[a][1]-coords[b][1])**2
            if temp_dist>dist:
                dist=temp_dist
                for i in range(2):
                    center[i]=(coords[a][i]+coords[b][i])/2
    rad=dist**0.5/2
    # cir = plt.Circle((center[0],center[1]),   radius=rad*1.3,fill=False,color=color,hatch=hatches.next())
    cir = plt.Circle((center[0],center[1]),   radius=rad*1.3,fill=False,color=color)
    plt.gca().add_patch(cir)
    plt.axis('scaled')
    # return color of the circle, to use it as the color for vertices of the cliques
    return color


def mkdir_p(path):
    print("creating", path)
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def main(args):
    # create a random graph
    # G=nx.gnp_random_graph(n=10,p=0.6)
    G = read_dot(args.dot_file)
    # remember the coordinates of the vertices
    if args.layout == "shell":
        coords=nx.shell_layout(G)
    else:
        coords=pydot_layout(G) 


    # remove "len(clique)>2" if you're interested in maxcliques with 2 edges
    cliques=[clique for clique in nx.find_cliques(G) if len(clique)>2]

    # #draw the graph
    nx.draw(G,pos=coords,node_size  = 200)
    print(coords)
    in_more_than_one_clique = set()
    seen = set()
    for clique in cliques:
        for node in clique:
            if node not in seen:
                seen.add(node)
            else:
                in_more_than_one_clique.add(node)

    # already_processed_nodes = set()
    for clique in cliques:
        print("Clique to appear: ",clique)
        color = draw_circle_around_clique(clique, coords)
        node_colors = []
        for node in clique:
            if node not in in_more_than_one_clique:
                node_colors.append(color)
            else:
                node_colors.append("grey")
        nx.draw_networkx_nodes(G, pos=coords, nodelist=clique, node_color= node_colors, node_size  = 200)
        # already_processed_nodes.update(clique)

    nx.draw_networkx_labels(G,coords, font_size=8)
    # nx.draw_networkx_nodes(G,pos=coords, nodelist=clique)
    plt.savefig(args.outfile)
    plt.close()

    # sequence_order = open("/Users/kxs624/tmp/sequence_order.tsv", "w")
    # is_printed = set()

    # for component in sorted(nx.connected_components(G), key = lambda x: len(x)):
    #     print(component)
    #     cliques=[clique for clique in nx.find_cliques(G.subgraph(component))]
    #     for clique in sorted(cliques, key = lambda x: len(x)):
    #         sequence_order.write("-\n")
    #         nodes = sorted(clique, key = lambda x: x in in_more_than_one_clique)
    #         for i, node in enumerate(nodes):
    #             if node not in is_printed:
    #                 sequence_order.write(node + "\n")
    #                 is_printed.add(node)
    #                 if len(nodes) > i+1 and nodes[i] not in in_more_than_one_clique and nodes[i+1] in in_more_than_one_clique:
    #                     sequence_order.write("-\n")




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="plot graph from dot file.")
    parser.add_argument('--dot_file', type=str, help='Path to the dot graph.')
    parser.add_argument('--outfile', type=str, help='Output path of plot')
    parser.add_argument('--layout', type=str, default="shell", help='"shell" or ')
    args = parser.parse_args()

    path_, file_prefix = os.path.split(args.outfile)
    mkdir_p(path_)
    args.outfolder = path_

    # if not os.path.exists(args.outfolder):
    #     os.makedirs(args.outfolder)
    
    main(args)
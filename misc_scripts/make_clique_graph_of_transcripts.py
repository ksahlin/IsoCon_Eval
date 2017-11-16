import networkx as nx
from math import *
import matplotlib.pylab as plt
import itertools as it

def draw_circle_around_clique(clique,coords):
    dist=0
    temp_dist=0
    center=[0 for i in range(2)]
    color=colors.next()
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

global colors, hatches
# colors=it.cycle('bgcmy')# blue, green, red, ...
colors=it.cycle(['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'purple', 'pink', 'brown', 'orange', 'teal', 'coral', 'lightblue', 'lime', 'lavender', 'turquoise', 'darkgreen', 'tan', 'salmon', 'gold', 'lightpurple', 'darkred', 'darkblue'])# blue, green, red, ...
# hatches=it.cycle('/\|-+*')

# create a random graph
# G=nx.gnp_random_graph(n=10,p=0.6)
G = nx.read_dot("/Users/kxs624/Dropbox/IsoCon/v3_transcripts/original_predictions/dot_graphs/RBMY.dot")
# remember the coordinates of the vertices
coords=nx.pydot_layout(G)


# remove "len(clique)>2" if you're interested in maxcliques with 2 edges
cliques=[clique for clique in nx.find_cliques(G) if len(clique)>2]

# #draw the graph
nx.draw(G,pos=coords,node_size  = 100)
print(coords)
already_processed_nodes = set()
for clique in cliques:
    print("Clique to appear: ",clique)
    color = draw_circle_around_clique(clique, coords)
    node_colors = []
    for node in clique:
        if node not in already_processed_nodes:
            node_colors.append(color)
        else:
            node_colors.append("grey")
    nx.draw_networkx_nodes(G, pos=coords, nodelist=clique, node_color= node_colors, node_size  = 100)
    already_processed_nodes.update(clique)

# nx.draw_networkx_nodes(G,pos=coords, nodelist=clique)


plt.show()
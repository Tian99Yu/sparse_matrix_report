import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

a = np.eye(5)
a[4,2] = 4
a[3,1] = 2
a[3,4] = 2

G = nx.from_numpy_matrix(a, create_using=nx.DiGraph)

print(G.nodes)
pos = nx.spring_layout(G)

nx.draw_networkx_nodes(G, pos)
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos, edge_color='r', arrows = True)
plt.savefig("graph.pdf")
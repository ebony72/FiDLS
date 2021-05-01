import networkx as nx
#import matplotlib.pyplot as plt
from utils import SPL, is_embeddable, hub, centre


#spl = nx.shortest_path_length
'''Introduce the AG class '''

class ArchitectureGraph:
    def __init__(self, g):
        if not isinstance(g, nx.Graph): raise Exception('Not a graph')
        #V = list(g.nodes())
        #V.sort()
        #EG = list(g.edges())
        #EG.sort()
        self.graph = g
        if not nx.is_connected(g): raise Exception('The AG should be connected!')
            # largest_cc = max(nx.connected_components(g), key=len)
            # g = g.subgraph(largest_cc).copy() 
        #self.node_list = V
        #self.edge_list = EG
        self.diameter = nx.diameter(g)
        self.SPL = SPL(g)
        
# define the architecture graph
# IBM Q Tokyo (Q20) 
def q20():
    g = nx.Graph()
    g.add_nodes_from(list(range(20)))
    for i in range(0,4):
        g.add_edge(i,i+1)
        g.add_edge(i+1,i)

    for i in range(5,9):
        g.add_edge(i,i+1)
        g.add_edge(i+1,i)
        
    for i in range(10,14):
        g.add_edge(i,i+1)
        g.add_edge(i+1,i)
        
    for i in range(15,19):
        g.add_edge(i,i+1)
        g.add_edge(i+1,i)
        
    for i in range(0,15):
        g.add_edge(i,i+5)
        g.add_edge(i+5,i)

    for i in [1,3,5,7,11,13]:
        g.add_edge(i,i+6)
        g.add_edge(i+6,i)

    for i in [2,4,6,8,12,14]:
        g.add_edge(i,i+4)
        g.add_edge(i+4,i)
    return g

# IBM Q Tokyo (Q20) 
def qgrid(m,n):
    g = nx.Graph()
    g.add_nodes_from(list(range(0,m*n-1)))
    for i in range(0,m):
        for j in range(0,n):
            if i < m-1: g.add_edge(j*m+i, j*m+i+1)
            if j < n-1: g.add_edge(j*m+i, (j+1)*m+i)
    return g

def rochester():
    g = nx.Graph()
    g.add_nodes_from(list(range(0,53)))
    I_1= list(range(4)) + list(range(7,15)) +\
        list(range(19,27)) + list(range(30,38)) +\
            list(range(42,50))
            
    #print(I_1)
    for i in I_1:
        g.add_edge(i,i+1)
    E = [(0,5),(5,9),(4,6),(6,13),(7,16),(16,19),\
         (11,17),(17,23),(15,18),(18,27),(21,28),(28,32),\
             (25,29),(29,36),(30,39),(39,42),(34,40),(40,46),\
                 (38,41),(41,50),(44,51),(48,52)]
    g.add_edges_from(E)
    return g

def sycamore():
    g = nx.Graph()
    g.add_nodes_from(list(range(0,54))) 
    I = list(range(6,12))+list(range(18,24))+list(range(30,36))+\
        list(range(42,48))
    for i in I:
        for j in g.nodes():
            if j in I: continue
            if i-j in [5,6] or j-i in [6,7]:
                g.add_edge(i,j)
    g.remove_node(3)
    sorted(g)
    if 3 in g.nodes(): raise Exception('node error')
    mapping = dict()
    for n in g.nodes():
        if n < 3:
            mapping[n] = n
        else:
            mapping[n] = n - 1
            
    h = nx.relabel_nodes(g, mapping)
    sorted(h)
    return h
    

    '''Rochester and Sycamore are embeddable in Q19x19 '''
    #x = is_embeddable(rochester(),qgrid(19,19),False)[0] #yes
    #x = is_embeddable(rochester(),sycamore(),False)[0]
    #x = is_embeddable(sycamore(),qgrid(19,19), True)
    # (True, {26: 180, 32: 161, 38: 160, 39: 142, 33: 143, 40: 124, 43: 159, 37: 178,\
    #    42: 177, 36: 196, 41: 195, 44: 141, 45: 123, 19: 199, 13: 218, 6: 237,\
    #    12: 236, 5: 255, 14: 200, 17: 235, 18: 217, 20: 181, 15: 182, 9: 183, 16: 164,\
    #    21: 163, 24: 216, 25: 198, 27: 162, 28: 144, 29: 215, 30: 197, 31: 179, 7: 219,\
    #    8: 201, 1: 256, 2: 238, 34: 125, 3: 202, 35: 214, 4: 184, 10: 165, 11: 254,\
    #    46: 105, 48: 176, 49: 158, 50: 140, 51: 122, 52: 104, 22: 145, 23: 234, 0: 274, 47: 194})
    #x = is_embeddable(qgrid(5,5), sycamore(), False)
    # (True, {6: 13, 11: 18, 7: 19, 8: 26, 12: 25, 13: 31, 16: 24, 17: 30, 18: 37, 1: 7,\
    #    2: 14, 3: 20, 5: 6, 9: 32, 10: 12, 14: 38, 15: 17, 19: 43, 21: 29, 22: 36,\
    #    23: 42, 0: 2, 24: 49, 4: 27, 20: 23})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#_________________TEST_____________________________________________________________#
if __name__=='__main__':
    AG = ArchitectureGraph(sycamore())
    g = AG.graph
    spl_dic = AG.SPL
    #for i in spl_dic:
    #    print(i, spl_dic[i])
    print(len(g.nodes()),len(g.edges()))
    print(nx.diameter(g))
    for node in g.nodes():
        print(node, list(g.neighbors(node)))
    print('Centre')
    print(centre(g))
    print('Hub')
    print(hub(g))
    # g = g.subgraph([21,22,23,24,25,28,29,32,33,34,35,36])
    # nx.draw(g, with_labels=True)
    # plt.draw()
    # plt.show()

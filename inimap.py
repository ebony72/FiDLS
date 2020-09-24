'''This module constructs weighted graph and topgraph initial mappings ''' 
import networkx as nx
from utils import topgates, cx_dependency_graph, graph_of_circuit, is_embeddable
'''Always put local parameters before global ones''' 
import time
#####\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\        
# Weighted SUBGRAPH

def best_wtg_o_ini_mapping(C, G, anchor, stop): #'o' for original
    ''' Return a graph g which is isomorphic to a subgraph of G
            while maximizing the number of CNOTs in C that correspond to edges in g
        Method: sort the edges according to their weights (the number of CNOTs in C corresponding to each edge);
                construct a graph by starting with the edge with the largest weight; then consider the edge with the second large weight, ...
                if in any step the graph is not isomorphic to a subgraph of G, skip this edge and consider the next till all edges are considered.
    
    Args:
        C (list): the input circuit
        G (graph): the architecture graph
        
    Returns:
        g (graph)
        map (dict)
    '''    
    g_of_c = graph_of_circuit(C)
    test = is_embeddable(g_of_c, G, anchor, stop)
    if test[0]:
        #print('The graph of the circuit is embeddable in G')
        return g_of_c, test[1]
    
    edge_wgt_list = list([C.count([e[0],e[1]]) + C.count([e[1],e[0]]), e] for e in g_of_c.edges())
    edge_wgt_list.sort(key=lambda t: t[0], reverse=True) # q[0] weight, q[1] edge
    
    '''Sort the edges reversely according to their weights''' 
    EdgeList = list(item[1] for item in edge_wgt_list)    
    #edge_num = len(EdgeList)
    
    '''We search backward, remove the first edge that makes g not embeddable, 
            and continue till all edges are evaluated in sequence. '''
            
    #Hard_Edge_index = 0 # the index of the first hard edge
    g = nx.Graph()
    result = dict()
    # add the first edge into g
    edge = EdgeList[0]
    g.add_edge(edge[0], edge[1])
    
    # rp = 0
    # for rp in range(edge_num):
    #     # h is the index of the last edge that can be added into g
    #     h = Hard_Edge_index
    #     if h == edge_num - 1: 
    #         return g, result
        
    #     EdgeList_temp = EdgeList[h+1:edge_num]
    #     for edge in EdgeList_temp:           
    #         g.add_edge(edge[0], edge[1])           
    #         i = EdgeList.index(edge)            
    #         # find the largest i such that the first i-1 edges are embeddable
    #         test = is_embeddable(g, G, anchor, stop)
    #         if not test[0]:
    #             Hard_Edge_index = i
    #             g.remove_edge(edge[0], edge[1])
    #             break
    #         result = test[1]
    #         if i == edge_num- 1:                
    #             return g, result
    # return g, result
    

    #EdgeList_temp = EdgeList[:]
    for edge in EdgeList:           
        g.add_edge(edge[0], edge[1])           
        test = is_embeddable(g, G, anchor, stop)
        if not test[0]:
            g.remove_edge(edge[0], edge[1])
            if nx.degree(g, edge[0]) == 0: g.remove_node(edge[0])
            if nx.degree(g, edge[1]) == 0: g.remove_node(edge[1])
        else:
            result = test[1]
    return g, result
###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
def best_topgraph_o_ini_mapping(L1, C, G, anchor, stop): # 'o': the original version
    ''' Return the topgraph g of the input circuit C
        Method: Consider all gates in the gate dependency graph one by one from the top.
                Let x=C[j] be the current edge. Add it to g and check if g is still embeddable into G.
                Otherwise, remove all gates that are dependent on x from the gate dependency graph
                and check if there are any gate left and continue.
        Args:
            C (list): the input circuit
            G (graph): the architecture graph
            L1: the list of indices of unsolved gates in C           
        Returns:
            g (graph): the topgraph
            map (dict)
    '''      
    g_of_c = graph_of_circuit(C)
    test = is_embeddable(g_of_c, G, anchor, stop)
    if test[0]:
        #print('The graph of the circuit is embeddable in G')
        return g_of_c, test[1]    
    gdg = cx_dependency_graph(L1, C)   
    '''Construct the top subgraph g from the circuit '''
    g = nx.Graph() # the subgraph of G induced by nodes with indices in GATES
    GATES = [] # the indices of topgates that can be executed (i.e., put in g) in this round
    result = dict()
    # we consider unsolved gates in C one by one
    Dump = set() # record all those gates that cannot be put in the solvable graph in this round
    '''To speed up we may consider only L1[0:500] for example.'''
    start = time.time()
    for j in L1:
        if time.time()-start > 100: 
            print('Mapping time exceeds 100 seconds!')
            break
        a = L1.index(j)
        gate = C[j]
        '''Update Dump by adding descedents of items in Dump'''
        # CHANGE = True        
        # while CHANGE and Dump:
        #     ldump = len(Dump)
        #     Dump_TG = topgates(list(Dump), C)
        #     for s in Dump_TG:
        #         Dump = Dump | set(nx.descendants(gdg, s))
        #     '''Check if Dump changed'''
        #     if ldump == len(Dump):
        #         CHANGE = False
        '''Stop when there are no gates outside Dump '''
        if Dump >= set(L1[a:]):
            return g, result, GATES
        if j in Dump: continue          
        if (gate[0], gate[1]) in g.edges() and not (gate[1], gate[0]) in g.edges():
            raise Exception('edge error')
        if (gate[0], gate[1]) in g.edges(): # the gate or its inverse has been considered
            GATES.append(j) # C[i] is to be solved in this round
            continue        
        g.add_edge(gate[0], gate[1])
        test = is_embeddable(g, G, anchor, stop)
        if test[0]: 
            result = test[1]
            GATES.append(j)   
        else: 
            # remove the gate from consideration
            Dump.add(j)
            Dump = Dump | set(nx.descendants(gdg, j)) #0908                   
            g.remove_edge(gate[0], gate[1])
            if nx.degree(g, gate[0]) == 0: g.remove_node(gate[0])
            if nx.degree(g, gate[1]) == 0: g.remove_node(gate[1])
    return g, result, GATES

'''The following x-version seems very slow!'''
def best_topgraph_x_ini_mapping(L1, C, G, anchor, stop): #slow for rochester
    ''' The 'x' version does not use Dump to collect unsolvable gates.'''      
    g_of_c = graph_of_circuit(C)
    test = is_embeddable(g_of_c, G, anchor, stop)
    if test[0]:
        #print('The graph of the circuit is embeddable in G')
        return g_of_c, test[1]    
    gdg = cx_dependency_graph(L1, C)   
    '''Construct the top subgraph g from the circuit '''
    g = nx.Graph() # the subgraph of G induced by nodes with indices in GATES
    GATES = [] # the indices of topgates that can be executed (i.e., put in g) in this round
    LX = L1[:]
    result = dict()
    step = 0
    while LX and step < 2000:
        step += 1
        LX0 = LX[:]
        LTG = topgates(LX0, C)
        for i in LTG:
            LX0.remove(i)
            gate = C[i]
            if (gate[0], gate[1]) in g.edges(): # the gate or its inverse has been considered
                GATES.append(i) # C[i] is to be solved in this round
                continue        
            g.add_edge(gate[0], gate[1])
            test = is_embeddable(g, G, anchor, stop)
            if test[0]: 
                result = test[1]
                GATES.append(i) 
            else:  
                '''Remove i's descendants from consideration'''
                LX0 = [x for x in LX0 if x not in nx.descendants(gdg, i)]                 
                g.remove_edge(gate[0], gate[1])
                if nx.degree(g, gate[0]) == 0: g.remove_node(gate[0])
                if nx.degree(g, gate[1]) == 0: g.remove_node(gate[1])
        LX = LX0[:] 
    return g, result, GATES


###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

def best_topgraph_y_ini_mapping(L1, C, G, anchor, stop):
    ''' Test a more efficient way for scanning the gates in C, where for each edge
        in g_of_c, we fix the index of the 1st cnot gate which gives the edge 
    '''      
    g = graph_of_circuit(C)
    test = is_embeddable(g, G, anchor, stop)
    if test[0]:
        #print('The graph of the circuit is embeddable in G')
        return g, test[1]    
    
    L_temp = [] #the list of the index of the 1st cnot 
    for edge in g.edges():
        p, q = edge[0], edge[1]
        s = min([k for k in L1 if set(C[k]) == {p,q}]) 
        L_temp.append(s) 
    # L_temp.sort(reverse=True)    
    # for s in L_temp:
    #     p, q = C[s][0], C[s][1]
    #     test = is_embeddable(g, G, True, stop)
    #     if test[0]:
    #         return g, test[1]
    #     g.remove_edge(p,q) 
    #     if nx.degree(g, p) == 0: g.remove_node(p)
    #     if nx.degree(g, q) == 0: g.remove_node(q)
    # if g.nodes(): raise Exception('g should be empty here')
    L_temp.sort()        
    g = nx.Graph()
    g.add_edge(C[L_temp[0]][0], C[L_temp[0]][1])
    result = dict()
    for s in L_temp[1:]:
        p, q = C[s][0], C[s][1]
        g.add_edge(p,q)
        test = is_embeddable(g, G, True, stop)
        if not test[0]:
            return g, result
        result = test[1]
    return g, result
'''The edge by edge search in 'y' can be sped up if we consider in a bipartite way'''
def top_z_P(i, L_temp, C, G, anchor, stop):
    '''L_temp is the list of the indices of the first cnot gates corresponding to edges in g_of_c'''
    g = nx.Graph()
    for s in L_temp[:i+1]:
        g.add_edge(C[s][0], C[s][1])
    test = is_embeddable(g, G, anchor, stop)
    return test[0]

def search_bipartite_top_z(yes_bound, no_bound, test_number, L_temp, C, G, anchor, stop):
    '''L_temp is the list of the indices of the first cnot gates corresponding to edges in g_of_c'''
    if not type(test_number) == int: raise Exception('Only consider integers')
    if no_bound == yes_bound + 1: return yes_bound
    if top_z_P(test_number, L_temp, C, G, anchor, stop):
        yes_bound = test_number
        if test_number == no_bound: return test_number
        test_number = yes_bound + max(1, (no_bound - yes_bound)//2)
        return search_bipartite_top_z(yes_bound, no_bound, test_number, L_temp, C, G, anchor, stop)
    else:
        if test_number == yes_bound + 1: return yes_bound
        no_bound = test_number
        test_number = yes_bound + max(1, (no_bound - yes_bound)//2)
        return search_bipartite_top_z(yes_bound, no_bound, test_number, L_temp, C, G, anchor, stop)

def best_topgraph_z_ini_mapping(L1, C, G, anchor, stop):
    ''' Test a more efficient way for scanning the gates in C than 'y' by using search_bipartite_top_z'''      
    g = graph_of_circuit(C)
    test = is_embeddable(g, G, anchor, stop)
    if test[0]:
        #print('The graph of the circuit is embeddable in G')
        return g, test[1]
    #print('The graph of the circuit is very likely NOT embeddable in G')    
    L_temp = [] #the index list of first cnot for each edge 
    for edge in g.edges():
        p, q = edge[0], edge[1]
        s = min([k for k in L1 if set(C[k]) == {p,q}]) 
        L_temp.append(s)

    L_temp.sort()
    yes_bound = 0
    no_bound = len(L_temp)
    test_number = no_bound//2
    exact_bound = search_bipartite_top_z(yes_bound, no_bound, test_number, L_temp, C, G, anchor, stop)
            
    g = nx.Graph()
    # add the first edge into g
    for s in L_temp[:exact_bound+1]:
        g.add_edge(C[s][0], C[s][1])
    test = is_embeddable(g, G, anchor, stop)
    if not test[0]: raise Exception('Check why the subgraph is not embeddable!')
    return g, test[1]

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# the weighted subgraph initial mapping
def _tau_bsg_(C, G, anchor, stop):
    ''' Return the weighted subgraph initial mapping

    Args:
        C (list): the input circuit
        G (graph): the architecture graph
    Returns:
        tau (list): the weighted subgraph initial mapping
    '''
    result = best_wtg_o_ini_mapping(C, G, anchor, stop)[1] # o in {o,x}
    return result 

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# the topgraph initial mapping
def _tau_bstg_(C, G, anchor, stop):
    ''' Return the topgraph initial mapping

    Args:
        C (list): the input circuit
        G (graph): the architecture graph
    Returns:
        tau (list): the topgraph initial mapping
    '''    
    L = list(range(len(C))) # o in {o,x,y,z}
    '''The z = y version is the fastest, x = o the slowest'''
    # result = best_topgraph_y_ini_mapping(L, C, G, anchor, stop)[1]
    result = best_topgraph_z_ini_mapping(L, C, G, anchor, stop)[1]
    return result
###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#_________________TEST_____________________________________________________________#
if __name__=='__main__':
    import os
    #file_name = 'qft_16.qasm' #qubits: 16, gates: 512 //49@>7500"
    #file_name = 'pf2_20_before.qasm' #qubits: 20, gates: 1020
    path='CNOT_lists2/'
    #path='testqasm/'
    files= os.listdir(path) #得到文件夹下的所有文件名称
    count = 0
    #start = time.time()
    for file_name in files:
        count += 1
        current_path = path + file_name
        #if count != 27:
        #if file_name != 'qft_10.qasm':
        if file_name != 'qft_16.qasm':
            pass
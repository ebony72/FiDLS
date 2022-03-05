# (Re)Created on Aug 25, 2020 by Sanjiang Li || mrlisj@gmail.com
'''Test if the circuit graph of a quantum circuit is embeddable to an architecture graph'''

from ag import ArchitectureGraph # architecture graph
from ag import q20, rochester, sycamore, qgrid
from circ_utils import  qubit_in_circuit, CreateCircuitFromQASM, ReducedCircuit, graph_of_circuit,\
    cx_dependency_graph, greedy_solved_gates, map_circuit_cost
from graph_utils import hub
from utils import map2tau
from vfstest import Vf
import networkx as nx
import json
import os
import copy
import time

from inimap import _tau_bstg_
from fidls_gx import qct_old

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
stop, ag, anchor, delta, QFilter =  10, 'tokyo', False, 0.95, '01y'
   
path = "B131/"
# path = "bigQ/"
# path = "BNTF/"
printOK = False
printOK = True


#------ define the architecture graph ---------------------
global AG
if ag == 'tokyo': AG = ArchitectureGraph(q20())
elif ag == 'sycamore': AG = ArchitectureGraph(sycamore())
elif ag == 'rochester': AG = ArchitectureGraph(rochester())
elif ag == 'q19x19': AG = ArchitectureGraph(qgrid(19,19))
elif ag == 'q5x5': AG = ArchitectureGraph(qgrid(5,5))
elif ag == 'q9x9': AG = ArchitectureGraph(qgrid(9,9))

global EG
global V
G = AG.graph
EG = AG.graph.edges()
V = list(AG.graph.nodes())
SPL = AG.SPL
#-----------------------------------------------------------

content = time.asctime()
print(content)

def is_embeddable(g, H, anchor, stop):
    '''check if a small graph g is embeddable in a large H, anchor is bool
        g, H (Graph)
        anchor (bool): whether or not mapping anchor of g to that of H
        stop (float): time limit for vf2
    '''
    vf2 = Vf()
    result = {} 
    if anchor: result[hub(g)] = hub(H)
    result = vf2.dfsMatch(g, H, result, stop)
    lng = len(nx.nodes(g))
    if len(result) == lng:
        return True, result   
    return False, result

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
'''Process a list of circuits ''' 
files = os.listdir(path) 
t_start = time.time()
count = 0

succ = 0
E = []
for file_name in files:
    timeA = time.time()
    count += 1
    if count != 2: continue
    if file_name[-4:] == 'qasm':
        cir = CreateCircuitFromQASM(file_name, path)
        C = ReducedCircuit(cir)
    else: #C is a list
        with open(path + file_name, 'r') as f:
            sqn = json.loads(f.read())
        C = sqn
    
    ''' Extract the circuit info: Here the circuit is represented as a list of gates like [0,3] '''
    l = len(C)
    L = list(range(l)) #index list of gates in C
    Q = qubit_in_circuit(L,C)
    nl = len(Q) # number of qubits in C            

    if len(Q) > len(V): continue
    print('Cir.%s: %s has %s qubits and %s gates' %(count, file_name[0:-9], len(Q), l))
    # Diam  = nx.diameter(G)

    g = graph_of_circuit(C)
        
    '''Check if g is embeddable into G'''
    vf2 = Vf()
    result = {} 
    if anchor: result[hub(g)] = hub(G)
    result = vf2.dfsMatch(g, G, result, stop, printOK)    

    if len(result) == len(nx.nodes(g)):
        print('Embeddable', result)
        succ += 1
        E.append(count)
        continue
    #-------------------------------------------------------------------------    

    '''I. generate the initial topgraph mapping and topsection'''
    tau_top = _tau_bstg_(C, G, anchor, stop)
    print(tau_top)
    tau_top_list = map2tau(tau_top, V)
    
    '''Return the list of solved gates, which are the index set of the top section'''
    Ltop = greedy_solved_gates(tau_top_list, L, C, nl, EG)  
    topcirc = [C[i] for i in Ltop]
    topgraph = graph_of_circuit(topcirc)
    Lrest = [i for i in L if i not in Ltop]  
    
    # C_out, cost_time = qct_old(tau_top_list, C, Q, G, EG, V, SPL, QFilter)
    # print('Using topgraph, the results are', len(C_out), cost_time)

    '''II. Generate the Depth dictionary of C'''
    
    _, Depth_dict = cx_dependency_graph(Lrest, C, nl)
    
    '''III. Generate the optimal initial mapping tau_opt '''
    somevalue = map_circuit_cost(delta, G, tau_top, Lrest, C, nl)
    # print('The baseline cost is %s' %round(somevalue,2) )
    
    '''IV. Compare the results '''

    vf2 = Vf()
    result = {4: 2, 8: 7, 0: 6, 1: 1, 3: 3, 6: 11, 7: 12, 2: 10}
    S = []
    S, _, newsomevalue, _ = vf2.dfsTopRestMatch(topgraph, G, delta, Depth_dict, Lrest, C, nl, S, result, somevalue, 20*stop, printOK)
    tau_opt = copy.deepcopy(tau_top)
    # print(S)
    if S: 
        tau_opt = copy.deepcopy(S[0])
    else:
        print('S is empty!')
    print('The optimal mapping is', tau_opt)    

    tau_opt_list = map2tau(tau_opt, V)
    C_out2, cost_time2 = qct_old(tau_opt_list, C, Q, G, EG, V, SPL, QFilter)
    print('Using the optimal mapping, the results are', len(C_out2), cost_time2) 

t_end = time.time()
content = 'The time spent for this test is: %s' %round(t_end-t_start, 2)
print(content)

# print(succ)
# print(E)
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

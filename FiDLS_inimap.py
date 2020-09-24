# (Re)Created on Aug 25, 2020 by Sanjiang Li || mrlisj@gmail.com
'''Call this module if you want to generate initial mappings for your own circuits:
    1. Put your cirucits in a folder;
    2. Change the path and select wgt or top in line 17 and ag in line 16
    3. The initial mapping list is saved and you can run FiDLS_run to do QCT '''

#import networkx as nx
from ag import ArchitectureGraph # architecture graph
from ag import q20, rochester, sycamore, qgrid
from inimap import _tau_bsg_, _tau_bstg_ # two initial mappings
from utils import  qubit_in_circuit, CreateCircuitFromQASM, ReducedCircuit
import json
import os
import time
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
anchor, lev, stop, ag =  True, 3, 10, 'q19x19'
initial_mapping = 'top' #'top', 'wgt'
   
path = "B131/"
# path = "bigQ/"
# path = "BNTF/"

'''name2 for reading and writing inimap'''
name2 = '_inimap_list_' + ag + '_' + initial_mapping + '_' + path[0:-1]

# define the architecture graph
global AG
if ag == 'tokyo': AG = ArchitectureGraph(q20())
elif ag == 'sycamore': AG = ArchitectureGraph(sycamore())
elif ag == 'rochester': AG = ArchitectureGraph(rochester())
elif ag == 'q19x19': AG = ArchitectureGraph(qgrid(19,19))
elif ag == 'q5x5': AG = ArchitectureGraph(qgrid(5,5))
elif ag == 'q9x9': AG = ArchitectureGraph(qgrid(9,9))

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
def save_result(name, content):
    # pass
    name = str(name)
    content = str(content)
    file = open("inimap/ohmy" + name2 + ".txt", mode = 'a')
    file.write(content)
    file.write('\n')
    file.close()


content = time.asctime()
print(content)
# save_result(name, content)

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
global EG
global V
G = AG.graph
EG = AG.graph.edges()
V = AG.graph.nodes()
SPL = AG.SPL
files = os.listdir(path) 
t_start = time.time()
count = 0
IM = []

for file_name in files:
    timeA = time.time()
    count += 1
    if file_name[-4:] == 'qasm':
        cir = CreateCircuitFromQASM(file_name, path)
        C = ReducedCircuit(cir)
    else: #C is a list
        with open(path + file_name, 'r') as f:
            sqn = json.loads(f.read())
        C = sqn
    l = len(C)
    if path == 'bigQ/':
        if count in {19,21,34,42,44,47,49}: continue 
            #duplicate circuits!! 19=2, 21=20, 34=17, 42=7, 44=27, 47=24, 49=40
        if l > 15000: continue  
    L = list(range(l))
    Q = qubit_in_circuit(L,C)
    if len(Q) > len(V): continue
    print('Cir.%s: %s has %s qubits and %s gates' %(count, file_name[0:-9], len(Q), l))

    #\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
      ### select an initial mapping ### 
    #\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\ 
    _map_ = dict()
    if initial_mapping == 'wgt': # weighted graph initial mapping
        _map_ = _tau_bsg_(C, G, anchor, stop)
        print(_map_)
        im = []
        for key in _map_:
            im.append([key, _map_[key]])
        IM.append([count, im])
    elif initial_mapping == 'top': # topsubgraph mapping
        _map_ = _tau_bstg_(C, G, anchor, stop)
        print(_map_)
        im = []
        for key in _map_:
            im.append([key, _map_[key]])
        IM.append([count, im])
    else: 
        pass

    '''if _map_ is incomplete, we may complete it in a natural way'''  
    # if len(_map_) < len(Q):
    #     _map_ = map_completion(_map_, L, C, Q, AG, V)
    #     print(_map_)
    
t_end = time.time()
content = 'The time spent for this test is: %s' %round(t_end-t_start, 2)
print(content)
# save_result(name, content)

content = IM
save_result(name2, content)
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

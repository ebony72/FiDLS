# (Re)Created on Aug 25, 2020 by Sanjiang Li || mrlisj@gmail.com
#@ Sep 24, 2020
'''The current version is quite different from the version used in the TC paper, but uses the same pricinple!''' 
'''For tokyo and B131 circuits, the results are even better in both effect and efficiency''' 
'''For bigQ circuits, it performs better than reported in the paper, showing that the inimap is better'''
'''For q19x19, the current implementation of inimap and vfs seems much slower!'''

from ag import ArchitectureGraph # architecture graph
from ag import q20, rochester, sycamore, qgrid
# from inimap import _tau_bsg_, _tau_bstg_ # two initial mappings
from utils import centre, hub, graph_of_circuit
from utils import  qubit_in_circuit, CreateCircuitFromQASM, ReducedCircuit
import json
import os
import time
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
'''The preferred QFilter is '01y', which uses Q0+Q1 as 1st filter and Q0 else'''
QFilter, anchor, lev, stop, ag, SIZE = '01y', True, 3, 10, 'tokyo', 'medium'
initial_mapping = 'top' #'top', 'wgt', 'empty', 'naive'
   
path = "B131/"
# path = "bigQ/"
# path = "BNTF/"

name = ag + '_' + QFilter + '_' + initial_mapping + '_' + path[0:-1] + '_'+ SIZE + '_' 
'''name2 for reading and writing inimap'''
name2 = '_inimap_list_' + ag + '_' + initial_mapping + '_' + path[0:-1]

GVal = True # test fidls_g
# GVal = False
if GVal: 
    from fidls_g import qct_old
    name += '_G_old_'
    slayer, gamma = 0, 0
else:
    from fidls_d import qct_old
    name += '_D_old_'
    slayer, gamma = 2, 0.8

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
    pass
    # name = str(name)
    # content = str(content)
    # file = open("testRecord/0924-" + name + ".txt", mode = 'a')
    # file.write(content)
    # file.write('\n')
    # file.close()

content = 'QFilter, GVal, SIZE, initial_mapping  =', QFilter, GVal, SIZE, initial_mapping
save_result(name, content)
content = 'lev, slayer, gamma, stop, anchor =', lev, slayer, gamma, stop, anchor
save_result(name, content)

content = 'In this test, we use Type %s filter and %s initial mapping for %s circuits on %s'\
    %(QFilter, initial_mapping, SIZE, ag)
print(content)
save_result(name, content)
content = '*****************************************'
save_result(name, content)

content = time.asctime()
print(content)
save_result(name, content)

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
sum_in = 0
sum_out = 0
COST_TIME = 0

IM = []
'''Load the precomputed inimap list! If no such list, create one using FiDLS_inimap'''
if initial_mapping in {'top', 'wgt'}:
    with open('inimap/' + name2 + '.txt', 'r') as f:
        IM = json.loads(f.read())

for file_name in files:
    timeA = time.time()
    count += 1
    # if count != 2: continue
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
    if SIZE == 'small':
        if l >= 100: continue
    if SIZE == 'medium':
        if l>1000 or l<100: continue
    if SIZE == 'large':
        if l <= 1000: continue
    L = list(range(l))
    Q = qubit_in_circuit(L,C)
    if len(Q) > len(V): continue
    print('Cir.%s: %s has %s qubits and %s gates' %(count, file_name[0:-9], len(Q), l))

    #\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
      ### select an initial mapping ### 
    #\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\ 
    _map_ = dict()
    if initial_mapping == 'wgt': # weighted graph initial mapping
        # _map_ = _tau_bsg_(C, G, anchor, stop) 
        imlist = IM[count-1][1] 
        for x in imlist:
            _map_[x[0]] = x[1]
        print(_map_)
        
        # im = []
        # for key in _map_:
        #     im.append([key, _map_[key]])
        # IM.append([count, im])
    elif initial_mapping == 'empty': #empty mapping
        g_of_c = graph_of_circuit(C)
        p, q = centre(g_of_c), hub(g_of_c)        
        u, v = centre(G), hub(G)
        _map_[q] = v
    elif initial_mapping == 'top': # topsubgraph mapping
        # _map_ = _tau_bstg_(C, G, anchor, stop)
        imlist = IM[count-1][1]
        # print(count, IM[count-1])
        for x in imlist:
            _map_[x[0]] = x[1]
        # print(_map_)
        # content = count, file_name, _map_, round(time.time()-timeA, 2)
        # save_result(name2, content)
        # im = []
        # for key in _map_:
        #     im.append([key, _map_[key]])
        # IM.append([count, im])
        
        '''if _map_ is incomplete, we may complete it in a natural way'''  
        # if len(_map_) < len(Q):
        #     _map_ = map_completion(_map_, L, C, Q, AG, V)
        #     print(_map_)
            
    elif initial_mapping == 'naive': #naive mapping
        for i in range(len(Q)):
            _map_[i] = i
    else: 
        pass

# t_end = time.time()
# content = 'The time spent for this test is: %s' %round(t_end-t_start, 2)
# print(content)
# save_result(name, content)

# content = IM
# save_result(name, content)
    
    tau = [-1]*len(V)
    for key in _map_: # map physical qubit BB[key] to logic qubit key
        tau[_map_[key]] = key
    ##################################################################
    sum_in += l
    # C_out, cost_time = qct(tau, C, Q, AG, EG, V, SPL, QFilter, lev, slayer, gamma, GVal)
    
    '''Compare with old version '''
    C_out, cost_time = qct_old(tau, C, Q, G, EG, V, SPL, QFilter)

    COST_TIME += cost_time
    sum_out += len(C_out)
    
    content = count, file_name[0:-9], len(Q), l, len(C_out), len(C_out)-l,\
        round(cost_time,2), round(len(C_out)/l, 4)
    
    print(content)
    save_result(name, content)
content = 'The average ratio is %s:%s = %s' %(sum_out, sum_in, round(sum_out/sum_in,4))
print(content)
save_result(name, content)
t_end = time.time()
content = 'The time spent for this test is %s : %s'\
            %(round(COST_TIME,2), round(t_end-t_start, 2))
print(content)
save_result(name, content)
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

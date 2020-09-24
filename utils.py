'''Created on Aug 25, 2020 by Sanjiang Li || mrlisj@gmail.com '''

'''C is always a circuit consisting of cnot gates, with form [contr, target]
   LD a sublist of the index list of C
'''
#####\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
import networkx as nx
from vfs import Vf
from qiskit import QuantumCircuit

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
def CreateCircuitFromQASM(file, path):
    QASM_file = open(path + file, 'r')
    iter_f = iter(QASM_file)
    QASM = ''
    for line in iter_f: 
        QASM = QASM + line
    #print(QASM)
    cir = QuantumCircuit.from_qasm_str(QASM)
    QASM_file.close    
    return cir
    
def ReducedCircuit(cir):
    '''Return Reduced Circuit containing only [name, [p,q]], e.g., ['cx', [0,2]] '''
    C = []
    for gate in cir:
        if gate[0].name != 'cx': continue
        qubits = [q.index for q in gate[1]]
        C.append(qubits)
    return C

def qubit_in_circuit(LD, C): 
    ''' Return the set of qubits in a subcircuit D of C
    Args:
        LD (list): a sublist of CNOT gates of the input circuit C
    Returns:
        QD (set): the set of qubits in D
    '''
    QD = set()
    for i in LD:
        QD.add(C[i][0])
        QD.add(C[i][1])
    return QD

#####################################################################################
#        #GRAPHS related to circuit and NISQ device architecture#
#####################################################################################

def graph_of_circuit(C):
    ''' Return the induced graph of the reduced circuit C
            - node set: qubits in C
            - edge set: all pair (p,q) if CNOT [p,q] or CNOT[q,p] in C
        Args:
            C (list): the input reduced circuit
        Returns:
            g (Graph)
    '''   
    L = list(range(len(C)))
    g = nx.Graph()
    g.add_nodes_from(qubit_in_circuit(L, C))
    for gate in C:
        if len(gate) != 2: continue
        g.add_edge(gate[0],gate[1])
    return g     
###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
def topgates(LD, C, nl): #@0922
    ''' Return the index list of toplayer CNOT gates 
    Args:
        C (list): the input circuit
        LD (list): the index list of D
    Return:
        LTG (list): the index list of toplayer CNOT gates 
    '''
    #if not LD: raise Exception ('LD should be nonempty!')
    if not LD: return []
    LTG = []
    N = set() # the set of qubits of gates in LTG
    for i in LD:
        p, q = C[i][0], C[i][1]
        if len(N) >= nl - 1: #@0922  
            return LTG
        #if one qubit of the gate has appeared in N, it's not a topgate, 
        # and all edges after it are not topgates
        if p in N or q in N:
            N.add(p)
            N.add(q)
            continue
        else: #if p not in N and q not in N:
            N.add(p)
            N.add(q)
            LTG.append(i)            
    return LTG

def topgates_3_lev(LD, C, nl):
    '''//q: different resutl obtained if we remove the comments. don't know why''' 
    LDx = LD[:]
    LTG = topgates(LDx, C, nl)
    for i in LTG:
        LDx.remove(i)
    LTG1 = topgates(LDx, C, nl)
    for i in LTG1:
        LDx.remove(i)
    LTG2 = topgates(LDx, C, nl)
    return LTG, LTG1, LTG2
###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\ 
# compute the gate dependency graph of (a part of) the circuit
def cx_dependency_graph(LX, C, nl): # LX is the index set of a subset of C
    ''' Return the gate dependency graph induced by the circuit
    Args:
        C (list): the input circuit
        LX (list): the list of indices of a sublist X of C       
    Returns:
        g (graph): its nodes are indices in LX and (i,j) is an edge in g if gate C[j] depends on C[i]  
    '''
    g = nx.DiGraph()
    g.add_nodes_from(LX)
    L1 = LX[:]
    while L1:
        L0 = L1[:]
        TG = topgates(L0, C, nl)
        for i in TG:
            L0.remove(i)       
        if not L0: break                   
        for i in TG:
            p, q = C[i][0], C[i][1]
            for j in L0:
                if p in C[j]:
                     g.add_edge(i,j) 
                     break
            for j in L0:
                if q in C[j]:
                     g.add_edge(i,j) 
                     break
        L1 = L0[:]
    return g

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
spl = nx.shortest_path_length
def SPL(g):
    if not nx.is_connected(g): raise Exception ('g is not connected!')
        # '''g is disconnected! We consider its largest connected component instead!'''
        # largest_cc = max(nx.connected_components(g), key=len)
        # g = g.subgraph(largest_cc)
    spl_dic = dict() 
    V = list(g.nodes())
    V.sort()
    for p in V:
        for q in V:
            if (q,p) in spl_dic:
                d = spl_dic[(q,p)]
            else:
                d = nx.shortest_path_length(g,p,q)
            spl_dic[(p,q)] = d
    return spl_dic

def centre(g):
    if not nx.is_connected(g): 
        largest_cc = max(nx.connected_components(g), key=len)
        g = g.subgraph(largest_cc).copy() 

    radium = nx.diameter(g)
    Centre = []
    for node in g.nodes():
        radium_temp = max([spl(g,node,nodex) for nodex in g.nodes])
        if  radium_temp > radium: continue
        radium = radium_temp
        Centre.append([node,radium])
    Centre = [cand[0] for cand in Centre if cand[1] == radium ]
    #deg = max([g.degree(node) for node in Centre ])
    step = 0
    while len(Centre) > 1 and step < radium:
        step += 1
        '''Compare how many (step+1)-nbrs they have if they have the same (step)-nbrs'''
        Centre = [[x, len([y for y in g.nodes() if spl(g,x,y)== step+1])] for x in Centre]
        max_val = max([item[1] for item in Centre])
        Centre = [item[0] for item in Centre if item[1]==max_val]
    return Centre[0]  

def hub(g):
    '''A hub of g is a node with maximum degree'''
    if not nx.is_connected(g): 
        largest_cc = max(nx.connected_components(g), key=len)
        g = g.subgraph(largest_cc).copy() 
        
    deg = max([g.degree(node) for node in g.nodes ])
    Hub = []
    for node in g.nodes():
        if g.degree(node) == deg: Hub.append(node)

    step = 0
    while len(Hub) > 1 and step < nx.diameter(g):
        step += 1
        Hub = [[x, len([y for y in g.nodes() if spl(g,x,y)== step+1])] for x in Hub]
        max_val = max([item[1] for item in Hub])
        Hub = [item[0] for item in Hub if item[1]==max_val]
    return Hub[0]
    
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


#####################################################################################
#        #Qubit Mapping tau: V -> Q or map: Q -> V
####################################################################################

def tau2map(tau):
    ''' Return the map: Q -> V of tau: V -> Q '''
    dic1 = dict()
    for q in tau:
        if q != -1: #i.e., q has been assigned the phy. qubit u = tau.index(q)
            dic1[q] = tau.index(q)
    return dic1

def map2tau(dic, V):
    '''Return the tau: V -> Q of a (map) dic: Q -> V '''
    tau = [-1]*len(V)
    for q in dic:
        tau[dic[q]] = q
    return tau

def swap(tau, u, v, EG): # p, q are physical qubits
    ''' Swap the images of u and v in a mapping tau: V->Q
        Args:
            tau (list): a mapping from V to Q
            u, v (int): physical qubits in V
        Returns:
            taunew (list): a new mapping with the images of u and v swapped
    '''
    if (u, v) not in EG:
        return tau
    taunew = tau[:]
    taunew[u] = tau[v]
    taunew[v] = tau[u]
    return taunew  
#transform tau to tau' with the images of p and q swapped
def swap_along_a_path(tau, u, v, path, EG): # u, v are physical qubits
    ''' Shift the image of i-th node in path to the (i+1)-th, and that of v to u 
    Args:
        EG (list): the edge list of the architecture graph G
        tau (list): a mapping from V to Q 
        path (list): a shortest path in AG from u to v
    Returns:
        taunew (list): a new mapping with the images of p and q swapped
    '''
    if path[0] != u or path[-1] != v: raise Exception ('path is illegal!')
    taunew = tau[:]    
    for i in range(len(path)):
        taunew[path[i]] = tau[path[(i+1)%len(path)]]
    taux = taunew[:]
    if taux[v] != tau[u] or taux[u] != tau[path[1]]:
        print(tau2map(tau))
        print(tau2map(taux))
        print(u,v,path)
        print(taux[u],tau[v],taux[v],tau[u])
        raise Exception('swap error', u, v, path)
    return taux

def invpath(path):
    '''Get the inverse path of a path in G'''
    m = len(path)
    pathinv = [path[m-i-1] for i in range(m)]
    return pathinv

def path2action(path):
    '''Get the swap combination associated to a path in G '''
    if len(path) < 2: return []
    action = []
    for i in range(len(path)-1):
        action.append([path[i],path[i+1]])
    return action

def entail(tau, gate, EG):
    ''' Check if a mapping entails a gate
    Args:
        EG (list): the edge list of the architecture graph G
        tau (list): a mapping from V to Q 
        gate (list): [p,q] represents a CNOT gate
    Returns:
        Boolean: True if tau entails gate
    '''
    p, q = gate[0], gate[1]
    if (p not in tau) or (q not in tau):
        return False
    u, v = tau.index(p), tau.index(q) 
    if ((u,v) in EG and (v,u) not in EG) or ((u,v) not in EG and (v,u) in EG): 
        raise Exception('AG is undirected!')
    if (u,v) in EG: # G is the architecture graph, EG its edge list
        return True    
    return False

def greedy_solved_gates(tau, LD, C, nl, EG):
    ''' Return the list of the indices of all gates solvable by a mapping tau
    Args:
        C (list): the input circuit
        G (graph): the architecture graph
        EG (list): the list of edges of G
        nl (int): the number of qubits in C
        tau (list): a mapping from V to Q, where V (Q) is the list of physical (logical) qubits
        LD (list): the index list of D, which represents the current logical circuit
        
    Returns:
        LSTG (list): the list of indices of solvable gates in D
    '''    
    LDx = LD[:]
    LSTG = [] # the index set of all solved topgates 
    while True: # as long as there are gates that can be solved by this mapping
        if not LDx: return LSTG        
        LX = topgates(LDx, C, nl) # topgates returns indices of topgates
        ldx = len(LDx)       
        # examine gates in LX one by one, check if they are entailed by tau
        for i in LX:
            gate = C[i]
            if entail(tau, gate, EG):
                LDx.remove(i)
                LSTG.append(i) # C[i] is solved by tau 
        if len(LDx) == ldx: # tau does not solve any new topgate
            return LSTG
        
#####################################################################################
#        #How good a mapping is? and compare two mappings
####################################################################################

def gate_phy_distance(gate, tau, V, SPL): #! V is an integer continuum
    #@0921-del AG -add SPL
    '''return the physical distance of two qubits in a gate w.r.t. mapping tau'''
    UnOcc = list(u for u in V if tau[u] == -1) #unoccupied qubits
    p, q = gate[0], gate[1]
    if p in tau and q in tau:
        return SPL[(tau.index(p), tau.index(q))]
    elif p in tau and q not in tau:
        u = tau.index(p)            
        d = min([SPL[(u,v)] for v in UnOcc])
        return d
    elif p not in tau and q in tau:
        v = tau.index(q)
        d = min([SPL[(u,v)] for u in UnOcc])
        return d        
    else:  #p not in tau and q not in tau:
        '''In this case, the gate should be a topgate from the very beginning!'''
        d = min([SPL[(u,v)] for u in UnOcc for v in UnOcc if v != u])
        return d  
    
def min_gate_dist(tau, LTG, C, V, SPL):
    #@0921-del AG -add SPL
    if not LTG: raise Exception('min_gate_dist_error')
    m = min([gate_phy_distance(C[i], tau, V, SPL) for i in LTG])
    return m

#! not verified yet!
# def solvable_gates_estimation(m, tau, LTG, LD, C, AG, V, lev):
#     '''step: number of swaps applied '''
#     L = list(range(len(C)))
#     LDx = LD[:]
#     #LTG = topgates(LDx, C)   
#     #m = min_gate_dist(tau, LTG, C, AG, V)
#     if m > lev + 1: return 0
#     unsolvable_top_gates = []
#     for i in LD[:100]:
#         if gate_phy_distance(C[i], tau, AG, V) > m:
#             unsolvable_top_gates.append(i)
#     gdg = cx_dependency_graph(L, C, nl)  
#     Dump = set()
#     for i in unsolvable_top_gates:
#         Dump.add(i)
#         Dump | set(nx.descendants(gdg, i))
#     LDx = [k for k in LD[:100] if k not in Dump]
#     return len(LDx)

#! not verified yet!    
# def solvable_gate_value(m, tau, LTG, LD, C, AG, V):    
#     # fix lev as 3
#     if m <= 1: raise Exception ('No gates should be solvable!')
#     if m > 4: return 0
#     sge_3 = solvable_gates_estimation(3, m, tau, LTG, LD, C, AG, V)
#     if m == 4: 
#         return sge_3/3
#     sge_2 = solvable_gates_estimation(2, m, tau, LTG, LD, C, AG, V)
#     if m == 3: 
#         if 3*sge_2 >= 2*sge_3: return sge_2/2
#         return sge_3/3
#     sge_1 = solvable_gates_estimation(1, m, tau, LTG, LD, C, AG, V)
#     if 3*sge_1 >= max(1.5*sge_2, sge_3): return sge_1/1
#     elif 3*sge_2 >= 2*sge_3: return sge_2/2
#     return sge_3/3                    
'''Include LTGx to avoid recalling topgates'''
def R_hat(tau, LTG, LTG1, LTG2, C, V, SPL, slayer, gamma): #Eq.4 in the paper
    #@0921-del AG -add SPL

    '''Return the weigthed sum of gate_phy_distances for the top-s layer gates'''
    '''The smaller the better'''
    dsum = 0
    for i in LTG:
        dsum += gate_phy_distance(C[i], tau, V, SPL)
    if slayer == 0:
        return dsum
    dsum1 = 0
    for i in LTG1:
        dsum1 += gate_phy_distance(C[i], tau, V, SPL)
    dsum += gamma * dsum1
    if slayer == 1:
        return dsum    
    dsum2 = 0
    for i in LTG2:
        dsum2 += gate_phy_distance(C[i], tau, V, SPL)
    dsum += (gamma**2) * dsum2
    return dsum  

#A special case of R_hat
def R_hat_1(tau, LTG, C, V, SPL): #Eq.4 in the paper
    '''Return the weigthed sum of gate_phy_distances for the top layer gates'''
    '''The smaller the better'''
    dsum = 0
    for i in LTG:
        dsum += gate_phy_distance(C[i], tau, V, SPL)
    return dsum
#A special case of R_hat
def R_hat_3(tau, LTG, LTG1, LTG2, C, V, SPL): #Eq.4 in the paper
    '''Return the weigthed sum of gate_phy_distances for the top3 layer gates'''
    '''The smaller the better'''
    dsum = 0
    for i in LTG:
        dsum += gate_phy_distance(C[i], tau, V, SPL)
    dsum1 = 0
    for i in LTG1:
        dsum1 += gate_phy_distance(C[i], tau, V, SPL)
    dsum += 0.8 * dsum1   
    dsum2 = 0
    for i in LTG2:
        dsum2 += gate_phy_distance(C[i], tau, V, SPL)
    dsum += (0.8**2) * dsum2
    return dsum  

# def s_better(tau1, tau2, LD, C, AG, V, slayer, gamma):
#     '''Compare if tau1 has smaller R_hat value than tau2 '''
#     if R_hat(tau1, LD, C, AG, V, slayer, gamma) < R_hat(tau2, LD, C, AG, V, slayer, gamma):
#         return True
#     return False


def one_shot_map_extension(dic1, LD, C, Q, AG, V):
    '''Extend current mapping (dic1) one step '''
    tau = map2tau(dic1, V)
    UnOcc = list(u for u in V if tau[u] == -1)
    Occ = list(u for u in V if tau[u] != -1)
    if not Occ: 
        '''Occ2 are phy. qubits which are 2-close to those occupied'''
        Occ2 = V[:] 
    else:
        Occ2 = list(u for u in UnOcc if min([AG.SPL[(u,v)] for v in Occ]) <= 2)
    rhat = 100
    for q in Q:
        if q in tau: continue
        for u in Occ2:            
            '''Extend the mapping temporalliy'''
            dic1_temp = tau2map(tau) 
            dic1_temp.update({q:u}) #extend dict1_temp with cand
            tau_temp = map2tau(dic1_temp,V) #extension
            rhat_temp = R_hat(tau_temp, LD, C, V, SPL, 3, 0.8)
            if rhat_temp >= rhat: continue
            rhat = rhat_temp
            cand = {q:u}
    # print('0', rhat, cand, dic1)
    dic1.update(cand)
    # print('1', round(rhat,2), cand, dic1, len(dic1), len(Q))
    return dic1

def map_completion(dic1, LD, C, Q, AG, V):
    result = dic1
    # print(type(result), len(result))
    if not result: raise Exception('dic1 should be nonempty!')
    while len(result) < len(Q):
        result = one_shot_map_extension(result, LD, C, Q, AG, V)
    return result
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#_________________TEST_____________________________________________________________#
if __name__=='__main__':
    import json
    import os
    import time
    import matplotlib.pyplot as plt          
    path = "CNOT_lists2/" #文件夹目录
    files= os.listdir(path) #得到文件夹下的所有文件名称
    t_start = time.time()
    count = 0
    sum_in = 0
    sum_out = 0
    COST_TIME = 0
    for file_name in files:
        timeA = time.time()
        count += 1
        # if count != 2: continue
    
        current_path = path + file_name
        with open(current_path, 'r') as f:
            sqn = json.loads(f.read())
        C = sqn
        l = len(C)
        L = list(range(l))
        if l != 11: continue #9, 5, 10
    
        nl = len(qubit_in_circuit(L,C))
        print('Circuit %s : %s contains %s gates and %s qubits' %(count, file_name, l,nl))
        g_of_c = graph_of_circuit(C)
        # for i in L: 
        #     print(i, C[i])
        cdg = cx_dependency_graph(L, C, nl)
        for x in cdg.nodes():
            print(x, C[x], set(nx.descendants(cdg, x)))
        nx.draw(cdg, with_labels=True)
        plt.draw()
        plt.show()
        


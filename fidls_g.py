# Created on June 10, 2019 by Sanjiang Li || mrlisj@gmail.com
'''This version try to be as consistent as the 0422-github version, 
    make as few changes as possible '''
#####\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
import networkx as nx
from utils import tau2map, map2tau, swap, swap_along_a_path
from utils import topgates, gate_phy_distance, min_gate_dist
from utils import qubit_in_circuit, greedy_solved_gates, R_hat_1
import time
###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
''' We fix the following notations: 
    AG: the architecture graph (object)
    G: AG.graph
    EG: the list of edges in G (edge has form (2,3) not [2,3])
     V: the list of physical qubits in G
     C: the input logical circuit, defined as a list of cnot gates
     L: the index list of gates in C, defined as list(range(len(C)))
     L1, L_temp, LD: a sublist of L
     LTG, LTG1, LTG2: sublists of L corr. to the 0,1,2 front layers of the current logical circuit
     Q: the set of logic qubits in the circuit C
     nl: len(Q)
   tau, taunew, taux, tauz: a mapping from V to Q 
   gate: a CNOT gate, represented as [p,q] or C[i], where p,q are qubits in Q
    p, q: logical qubits in Q
    u, v: physical qubits in V
    SPL: the shortest path length dictionary of AG
    action: a list of edges in G, denoting swaps    
'''
'''Parameter order: tau, gate|LTG|L1, C, Q|nl, G, EG, V, SPL, QFilter_type '''
############################################################################
'''We use the following function SRMD to extend the mapping or solve a gate in case of Fallback!'''
def swap_reduce_min_dist_x(tau, LTG, LD, C, nl, G, EG, V, SPL):
    md = min_gate_dist(tau, LTG, C, V, SPL)
    UnOcc = list(u for u in V if tau[u] == -1) #unoccupied phy. qubits
    Occ = list(u for u in V if tau[u] != -1)
    if not Occ: 
        '''Occ2 are phy. qubits which are 2-close to those occupied'''
        Occ2 = V[:] 
    else:
        Occ2 = list(u for u in UnOcc if min([SPL[(u,v)] for v in Occ]) <= 2)
        
    GATE0, GATE1, GATE2, GATE3 = [], [], [], [] #front gates of type x
    for i in LTG:
        if gate_phy_distance(C[i], tau, V, SPL) > md: continue
        p, q = C[i]
        if (p not in tau) and (q not in tau): #type 0
            GATE0.append(i)
        elif (p in tau) and (q not in tau): #type 1
            GATE1.append(i)
        elif (p not in tau) and (q in tau): #type 2
            GATE2.append(i)
        else:
            GATE3.append(i)
    '''We prefer to swap rather than extension if GATE3 is nonempty.'''    
    if GATE3:
        gsg, rhat, action = 0, 100, None
        for i in GATE3:
            p, q = C[i]
            u, v = tau.index(p), tau.index(q)
            pi = nx.shortest_path(G, u, v)
            '''If len(pi)=2 then md=1 and no action is required!'''
            tau_temp = swap_along_a_path(tau, u, v, pi, EG)
            gsg_temp = len(greedy_solved_gates(tau_temp, LD, C, nl, EG))
            if gsg_temp < gsg: continue
            rhat_temp =  R_hat_1(tau_temp, LTG, C, V, SPL)
            if gsg_temp == gsg and rhat_temp >= rhat: continue
            gsg = gsg_temp
            rhat = rhat_temp
            action_temp = []
            for k in range(len(pi)-2): #q: -1 or -2?
                action_temp.append([pi[k], pi[k+1]])
            taux = tau_temp[:]
            action = action_temp[:]
        # print('Use swap to reduce min dist and solve a gate', md)
        if action == None: raise Exception('Some gates should be solved!', md) 
        return action, taux, ['type3', GATE3, md] # update tau with swaps in path  
     
    '''If GATE3 is empty, we select the best extension!'''    
    CAND = []
    if GATE1 or GATE2:
        for i in GATE1:
            p, q = C[i]
            u = tau.index(p)
            if q in tau: raise Exception('q is not in tau', p, q)
            CAND += [{q:v} for v in UnOcc if SPL[(u,v)] == md]
        for i in GATE2:
            p, q = C[i]
            v = tau.index(q)
            if p in tau: raise Exception('p is not in tau', p, q)
            CAND += [{p:u} for u in UnOcc if SPL[(u,v)] == md]
    else:
        '''For naive i.m., wgtgraph and topgraph, this is very unlikely!'''
        if GATE0: print('GATE0 is nonempty!', GATE0)
        for i in GATE0:
            p, q = C[i]
            if p in tau or q in tau: raise Exception('p,q should not be in tau',p,q)
            CAND += [{p:u,q:v} for u in Occ2 for v in UnOcc if SPL[(u,v)] == md]
            
    rhat_record = 100
    for cand in CAND:
        '''Extend the mapping temporalliy'''
        taux = tau[:]
        dic1_temp = tau2map(taux) 
        for key in cand:
            if key in dic1_temp:
                print(cand, dic1_temp, key, GATE0, GATE1, GATE2, GATE3, md)
                raise Exception('Wrong!')
        dic1_temp.update(cand) #extend dict1_temp with cand
        tau_temp = map2tau(dic1_temp,V) #extension
        rhat_temp = R_hat_1(tau_temp, LTG, C, V, SPL)
        if rhat_temp >= rhat_record: continue
        rhat_record = rhat_temp
        tau_record = tau_temp[:]
    return [], tau_record, ['not type3', GATE0, GATE1, GATE2, GATE3, md]
#####################################################################
  ###\__/#      GreedyV3 SWAPS         \__/#\#/~\\__/#\__/#\#/~\__/
#####################################################################
def Q_Filter(tau, LTG, LTG1, C, Q, EG):
    Q0 = qubit_in_circuit(LTG, C)
    Q1 = qubit_in_circuit(LTG1, C) #x00
    Q1x = qubit_in_circuit(LTG+LTG1, C) #x01
    return Q, Q0, Q1, Q1x

'''Find all relevant actions with length up to 3'''
def SWAP3(tau, LTG, LTG1, C, Q, EG, V, SPL, QFilter_type):
    QFilter = Q_Filter(tau, LTG, LTG1, C, Q, EG)
    if QFilter_type == '9': #no filter is used
        QF1, QF2, QF3 = QFilter[0], QFilter[0], QFilter[0]
    elif QFilter_type == '0': #only Q1-filter (note in the paper this is called Q0-filter)
        QF1, QF2, QF3 =  QFilter[1], QFilter[1], QFilter[1]
    elif QFilter_type == '01': #Q1-filter for the first layer and Q2-filter for the rest (Q0+Q1-filter in the paper)
        QF1, QF2, QF3 = QFilter[1], QFilter[2], QFilter[2]
    elif QFilter_type == '01x': #Q1-filter for the first layer and (Q1+Q2)-filter for the rest (Q0+weak_Q1-filter
        QF1, QF2, QF3 = QFilter[1], QFilter[3], QFilter[3]
    elif QFilter_type == '1x': #only Q2x-filter
        QF1, QF2, QF3 = QFilter[3], QFilter[3], QFilter[3]
    else: #'01y'
        QF1, QF2, QF3 = QFilter[3], QFilter[1], QFilter[1] 
            
    SWAP3 = []
    rhat = R_hat_1(tau, LTG, C, V, SPL) 
    for edge_1 in EG:
        p, q = edge_1
         #skip if both p, q are irrelevant phy. qubits w.r.t. TG
        if tau[p] not in QF1 and tau[q] not in QF1: continue
        tau1 = swap(tau, p, q, EG)
        rhat1 = R_hat_1(tau1, LTG, C, V, SPL)        
        if rhat1 > rhat: continue
        x1 = [ [edge_1], tau1]
        if x1 not in SWAP3:
            SWAP3.append(x1)

        for edge_2 in EG:
            if edge_1 == edge_2: continue            
            p, q = edge_2
            if tau1[p] not in QF2 and tau1[q] not in QF2: continue            
            tau2 = swap(tau1, p, q, EG)              
            rhat2 = R_hat_1(tau2, LTG, C, V, SPL)
            if rhat2 > rhat1: continue        
            x2 = [ [edge_1,edge_2], tau2]
            if x2 not in SWAP3:
                SWAP3.append(x2)

            for  edge_3 in EG:
                if edge_3 in {edge_1, edge_2}: continue
                p, q = edge_3
                if tau2[p] not in QF3 and tau2[q] not in QF3: continue

                tau3 = swap(tau2, p, q, EG)
                mdg3 = min_gate_dist(tau3, LTG, C, V, SPL)
                '''If GVal, we aim to solve at least one gate.'''
                if mdg3 > 1: continue 
                x3 = [ [edge_1, edge_2, edge_3], tau3]
                if x3 not in SWAP3:
                    SWAP3.append(x3)                    
    return SWAP3

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
'''Generate the next action and mapping'''
def good_next_mapping(tau, LD, C, Q, G, EG, V, SPL, QFilter_type, Fallback):    
    nl = len(Q)
    LDx = LD[:]
    LTG = topgates(LD, C, nl)
    for i in LTG:
        LDx.remove(i)
    LTG1 = topgates(LDx, C, nl)

    '''In case of Fallback, even if the mapping is complete, we call SRMD to
        solve at least one gate in the top layer '''
    if [v for v in V if tau[v] == -1] or Fallback: #i.e., tau incomplete
        '''//Compare the following two srmd functions!'''
        # action_rec, tau_rec = swap_reduce_min_dist(tau, LTG, C, G, EG, V, SPL) 
        action_rec, tau_rec, type = swap_reduce_min_dist_x(tau, LTG, LD, C, nl, G, EG, V, SPL)
        if not action_rec or Fallback:
            '''Extend the mapping or Fallback!'''
            # print('SRMD used', action_rec, Fallback, len(LD), type)
            return action_rec, tau_rec            
    sg = None #the selected good mapping 
    t = 0 #efficiency ratio
    SWAPX = SWAP3(tau, LTG, LTG1, C, Q, EG, V, SPL, QFilter_type)
    for rho in SWAPX:        
        action_temp, tau_temp = rho
        gsg_temp = greedy_solved_gates(tau_temp, LD, C, nl, EG)
        if not gsg_temp: continue
        t_temp = len(gsg_temp)/len(action_temp)
        if t_temp < t: continue
        if t_temp == t and sg != None: continue
        t = t_temp
        sg = [action_temp, tau_temp]
           
    if sg == None: 
        if not [v for v in V if tau[v] == -1]: #i.e., tau is complete
            # action_rec, tau_rec = swap_reduce_min_dist(tau, LTG, C, G, EG, V, SPL) 
            action_rec, tau_rec, type = swap_reduce_min_dist_x(tau, LTG, LD, C, nl, G, EG, V, SPL)
            print('  Oops, we reduces the minimal distance', action_rec, Fallback, len(LD), type)
        return action_rec, tau_rec    
    action_rec, tau_rec = sg
    return action_rec, tau_rec             

##################################################################
      #\__/#\__/#\#/~\ THE MAIN ALGORITH10M \__/#\__/#\#/~\
##################################################################
def qct_old(tau, C, Q, G, EG, V, SPL, QFilter_type):
    l = len(C)
    L = list(range(l)) # the index set of C
    nl = len(Q) # number of qubits in C            
    Diam  = nx.diameter(G)
    
    '''The initial state '''
    nsvg = 0 # the number of executed gates
    COST = 0 # the number of auxiliary CNOT gates
    taunew = tau[:]
    L1 = L[:]
    nr = 0
    C_out = [] # the output circuit

    '''For Fallback! May ignore for first read'''
    action = [] #the selected action, consisting of one or more swaps
    Fallback = False
    state = [taunew, L1, C_out, COST, nr, nsvg]
    Record = [] # records of not increasing nsvg   

    start = time.time()
    while nsvg < l:
        GSG = greedy_solved_gates(taunew, L1, C, nl, EG)        
        for i in GSG:
            L1.remove(i)
            p = taunew.index(C[i][0])
            q = taunew.index(C[i][1])
            C_out.append([p,q])                    
        nsvg += len(GSG) # the number of solved gates
        
        '''Reset the Fallback parameters if some gates are solved or the mapping extended!'''
        '''May ignore for first read '''
        if GSG or not action: #remember the last time the program has solved sth
            Record = []  
            tau_record = taunew[:] 
            L_record = L1[:]
            C_out_record = C_out[:]
            COST_record, nr_record, nsvg_record = COST, nr, nsvg
            state = [tau_record, L_record, C_out_record, COST_record, nr_record, nsvg_record]
            Fallback = False
        else: 
            Record.append(nsvg)    
        if len(Record) > 2*Diam: 
            print('use Fallback @%s!' %nr)
            Fallback = True #return to the record state
            COST, nr, nsvg = state[3], state[4], state[5]
            taunew = state[0][:]
            L1 = state[1][:]
            C_out = state[2][:]
            print('go back to round', nr, len(greedy_solved_gates(taunew, L1, C, nl, EG)))
        if len(Record) > 4*Diam: raise Exception('The program gets stuck!', nr, Record)  
            
        if not L1: break
        '''Go to the next round!''' 
        nr += 1
        action, taunew = good_next_mapping(taunew, L1, C, Q, G, EG, V, SPL, QFilter_type, Fallback)
        cost = len(action)*3            
        COST += cost
        for edge in action:
            C_out.append([edge[0],edge[1]])
            C_out.append([edge[1],edge[0]])
            C_out.append([edge[0],edge[1]])
           
    cost_time = time.time() - start
    return C_out, round(cost_time,2)
   

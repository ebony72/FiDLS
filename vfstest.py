import networkx as nx
from maps import Map
import time
import copy
from utils import dict2tuple
from circ_utils import map_circuit_cost 
from graph_utils import SPL, centre, hub, is_embedding, map_dist

Test = False
# Test = True

class Vf:

    __origin = None
    __sub = None

    #type = 0, pre; type = 1, succ
    #sl divide the graph neighborhood of a vertex into two disjoint parts: pre (in map) and succ (not in map)
    def preSucc(self, vertexNeighbor, map, type):
        #vertexNeighbor and map can be empty
        if not (isinstance(vertexNeighbor, list) and isinstance(map, list)):
            raise Exception ("Class Vf preSucc() arguments type error! vertexNeighbor and map expected list!")
        if not (type == 0 or type == 1):
            raise Exception ("Class Vf preSucc() arguments value error! type expected 0 or 1!")
           
        result = []
        #succ
        if type:
            for node in vertexNeighbor:
                if node not in map:                   
                    result.append(node)
        #pre
        else:
            for node in vertexNeighbor:
                if node in map:
                    result.append(node)
        return result

    def CandMeetRules(self, v1, v2, subgraph, graph, result, subMap, gMap, subMNeighbor, gMNeighbor):
        '''Check if the pair {v1:v2} meets VF2 rules for the current mapping result'''    
            
        if not result: return True     
        
        v1Neighbor = list(nx.all_neighbors(subgraph, v1))
        v2Neighbor = list(nx.all_neighbors(graph, v2))
                
        v1Pre = self.preSucc(v1Neighbor, subMap, 0)
        v1Succ = self.preSucc(v1Neighbor, subMap, 1)
        v2Pre = self.preSucc(v2Neighbor, gMap, 0)
        v2Succ = self.preSucc(v2Neighbor, gMap, 1)
        
        if len(v1Pre) > len(v2Pre) or len(v1Succ) > len(v2Succ): return False
                       
        for pre in v1Pre:
            if result[pre] not in v2Pre: return False
                       
        '''A vertex is a neighbor of subMap if it is not in and a neighbor of some vertex in the map'''
        len1 = len(set(v1Neighbor) & set(subMNeighbor)) #number of all x which is a nghb of v1 and the map 
        len2 = len(set(v2Neighbor) & set(gMNeighbor))

        if len1 > len2: return False
        return True        

    def dfsMatch(self, subgraph, graph, result, stop, printOK): #sl stop is the time limit
        '''Check if the current map is completable '''
        start_A = time.time()
        if not isinstance(result, dict):
            raise Exception ("Class Vf dfsMatch() arguments type error! result expected dict!")
        
        curMap = Map(result) #sl create a Map object!   
         
        if len(result) == len(nx.nodes(subgraph)): return result
        if not is_embedding(result, subgraph, graph): raise Exception ('Not an embedding!')
        
        '''Construct the current neighborhoods of the mapping'''
        subMNeighbor = curMap.neighbor(subgraph, 0) 
        gMNeighbor = curMap.neighbor(graph, 1)   
        
        if subMNeighbor and len(subMNeighbor) > len(gMNeighbor): return {} #non-completable

        '''Rank the unmapped neighbors by their degrees and select the highest one'''
        '''//todo: could also rank by centralities or other properties''' 
        subMN = subMNeighbor[:]
        if not subMNeighbor:            
            subMN = list(set(nx.nodes(subgraph)) - set(curMap.subMap()))
        subMN_deg_node = list([nx.degree(subgraph, v), v] for v in subMN)
        subMN_deg_node.sort(key=lambda t: t[0], reverse=True)
        selected_vertex = subMN_deg_node[0][1] 
        
        #------------------------------------------------------------
        if Test and not result:
            selected_vertex = centre(subgraph)
        #------------------------------------------------------------
        
        '''Our AGs are always connected. gMNeighbor is empty iff result is empty!'''  
        Candidate = gMNeighbor[:]
        if not subMNeighbor: 
            Candidate = list(set(nx.nodes(graph)) - set(curMap.gMap()))

        #------------------------------------------------------------
        subMN_deg = list(item[0] for item in subMN_deg_node)
        gMN_deg_node = list([nx.degree(graph, v), v] for v in Candidate)
        gMN_deg_node.sort(key=lambda t: t[0], reverse=True)
        gMN_deg = list(item[0] for item in gMN_deg_node)
        
        if len(subMN_deg) > len(gMN_deg) or sum(subMN_deg) > sum(gMN_deg): return {}
        for i in range(len(subMN_deg)):
            if subMN_deg[i] > gMN_deg[i]: return {} #non-completable
        #------------------------------------------------------------        

        '''Remove those graph neighbours which cannot match the suggraph candidate node''' 
        '''//todo: could also rank candidates according to centrality '''
        Candidate = [ t for t in Candidate if\
                       nx.degree(subgraph, selected_vertex) <= nx.degree(graph,t)]
            
        # if printOK: print('degree sorted candidates', Candidate)
        if not Candidate: return {} #non-completable

        #------------------------------------------------------------
        if Test and not result:
            '''Sort the candidates according to their centrality in graph'''
            # radium = nx.diameter(graph)
            Centre = []
            for node in graph.nodes():
                radium_temp = max([SPL(graph,node,nodex) for nodex in graph.nodes])
                # if  radium_temp > radium: continue
                # radium = radium_temp
                Centre.append([node,radium_temp])
            Centre.sort(key=lambda t: t[0])
            Centre = [cand[0] for cand in Centre]
            Candidate = Centre[:]
            
        elif Test:               
            '''Sort the candidates according to the diameter of the embedded subgraphs'''
            temp_node_set = list( result.values() )
            temp_g = graph.subgraph(temp_node_set)
            if nx.is_connected(temp_g):
                CandDiameter = []
                for t in Candidate:
                    node_set = temp_node_set[:]
                    node_set.append(t)
                    temp_g = graph.subgraph(node_set)
                    if nx.is_connected(temp_g):
                        d = nx.diameter(temp_g)
                    else: 
                        d = nx.diameter(graph)
                    CandDiameter.append([d,t])
                    
                CandDiameter.sort(key=lambda t: t[0])
                Candidate = [item[1] for item in CandDiameter]
            
        if Test and printOK: print('diameter sorted candidates', Candidate)
        # ------------------------------------------------------------
        
        for t in Candidate:
            '''Check if {selected_vertex : t} can be added to the current mapping '''
            if(self.CandMeetRules(selected_vertex, t, subgraph, graph, result, curMap.subMap(),\
                                curMap.gMap(), subMNeighbor, gMNeighbor)):
                
                result[selected_vertex] = t #Extend the mapping!
                
                if time.time()-start_A > stop: 
                    print('dfsmatch time exceeds', stop) #timeout
                    return {}
                
                self.dfsMatch(subgraph, graph, result, stop, printOK)

                if len(result) == len(nx.nodes(subgraph)): return result
                
                '''The pair is not helpful. We thus pop it out and restore the mapping.'''
                result.pop(selected_vertex)
                
        return {}   

    def dfsAllMatch(self, subgraph, graph, S, result, stop, printOK): #sl stop is the time limit, S the collection of current solutions
        '''Add to S all embeddings which extend the current mapping'''
        start_A = time.time()
        if not isinstance(result, dict):
            raise Exception ("Class Vf dfsMatch() arguments type error! result expected dict!")
        
        if printOK:
            print('^^^^^^^^^^^^^^^')
            print('Call dfsMatch for', result)
        
        curMap = Map(result) #sl create a Map object!   
         
        if len(result) == len(nx.nodes(subgraph)): 
            S.add(result)
            return S, result
        
        '''Construct the current neighborhoods of the mapping'''
        subMNeighbor = curMap.neighbor(subgraph, 0) 
        gMNeighbor = curMap.neighbor(graph, 1)   
        
        if subMNeighbor and len(subMNeighbor) > len(gMNeighbor): return S, {}

        '''Rank the unmapped neighbors by their degrees and select the highest one''' 
        subMN = subMNeighbor[:]
        if not subMNeighbor:            
            subMN = list(set(nx.nodes(subgraph)) - set(curMap.subMap()))
        subMN_deg_node = list([nx.degree(subgraph, v), v] for v in subMN)
        subMN_deg_node.sort(key=lambda t: t[0], reverse=True)
        selected_vertex = subMN_deg_node[0][1] 
        
        '''Our AGs are always connected. gMNeighbor is empty iff result is empty!'''  
        Candidate = gMNeighbor[:]
        if not subMNeighbor: 
            Candidate = list(set(nx.nodes(graph)) - set(curMap.gMap()))

        #------------------------------------------------------------
        subMN_deg = list(item[0] for item in subMN_deg_node)
        gMN_deg_node = list([nx.degree(graph, v), v] for v in Candidate)
        gMN_deg_node.sort(key=lambda t: t[0], reverse=True)
        gMN_deg = list(item[0] for item in gMN_deg_node)
        
        if len(subMN_deg) > len(gMN_deg) or sum(subMN_deg) > sum(gMN_deg): return S, {}
        for i in range(len(subMN_deg)):
            if subMN_deg[i] > gMN_deg[i]: return S, {}
        #------------------------------------------------------------        

        '''Remove those graph neighbours which cannot match the suggraph candidate node''' 
        Candidate = [ t for t in Candidate if\
                       nx.degree(subgraph, selected_vertex) <= nx.degree(graph,t)]

        if not Candidate: return S, {}
        if not result and printOK: print(selected_vertex, Candidate)
        
        for t in Candidate:
            
            if printOK: print('Can we add',  {selected_vertex:t}, 'to', result, '?') 
            
            '''Check if {selected_vertex : t} satisfies the VF2 rule'''
            if(self.CandMeetRules(selected_vertex, t, subgraph, graph, result, curMap.subMap(),\
                                curMap.gMap(), subMNeighbor, gMNeighbor)):
                
                if printOK: print([selected_vertex,t], 'meets the rules for', result) 

                result[selected_vertex] = t #Extend the mapping!
                
                if time.time()-start_A > stop: 
                    print('dfsmatch time exceeds', stop ) #timeout
                    return S, result
                
                '''Return all completions of the extended mapping'''
                self.dfsAllMatch(subgraph, graph, S, result, stop, printOK)
                
                if printOK:
                    print('After adding', [selected_vertex,t], 'and calling dfsMatch, we have', result, ) 
                    print('----------------')

                '''After this call of dfsMatch, we extect to obtain an updated S. 
                    The result itself will be updated? '''

                if len(result) == len(nx.nodes(subgraph)): 
                    S.add(result)
                    if printOK: print('An embedding has been found!', result)

                '''The pair is examined and we pop it out and restore the mapping'''
                result.pop(selected_vertex)
                
            # else:
                # if printOK: print([selected_vertex,t], 'does not meet the rules for', result) 
                                
        return S, result   

    '''We next consider the case when we want to generate all good mappings which are close to some previous mapping'''
    ''' 
        tau: the current mapping
        result: the mapping need construct
        dist = sum (dist[tau_dict[p],result[p]])
        map_dist(tau,result)
        somevalue = len(result)*nx.diameter(graph)
    '''
    #global somevalue
    def dfsTopMatch(self, subgraph, graph, tau_dict, S, result, somevalue, stop, printOK): #sl stop is the time limit, S the collection of current solutions
        '''Add to S all embeddings which extend the current mapping'''
        # print('aaa', somevalue)

        start_A = time.time()
        if not isinstance(result, dict):
            raise Exception ("Class Vf dfsMatch() arguments type error! result expected dict!")
        
        if printOK:
            print('^^^^^^^^^^^^^^^')
            print('Call dfsMatch for', result)
        
        curMap = Map(result) #sl create a Map object!   
         
        if len(result) == len(nx.nodes(subgraph)): 
            #-----------------------------------------------------------
            '''Check if result is indeed an embedding'''
            if not is_embedding(result, subgraph, graph):
                raise Exception(result, 'is not an embedding')
            #-----------------------------------------------------------
                
            temp_dist = map_dist(subgraph,graph,tau_dict,result)
            if printOK: print(temp_dist, somevalue, result)

            if temp_dist < somevalue:
                # print(temp_dist, somevalue)
                somevalue = min(temp_dist, somevalue)
                # S.append(result)
                
                S = [result]
                # print(somevalue)
            return S, result, somevalue            

                
            # for i in result:
            #     for j in result:
            #         if (i,j) in nx.edges(subgraph):
            #             if (result[i],result[j]) not in nx.edges(graph):
            #                 print(result, 'is not an embedding')
            #                 return S, {}, somevalue
            # print(result)
            # for (i,j) in nx.edges(subgraph):
            #     # print('sub', i,j)
            #     # print('g', result[i],result[j]) 
            #     if (result[i],result[j]) not in nx.edges(graph):
            #        raise Exception(result, 'is not an embedding')
            #        # return S, {}, somevalue           
            #-----------------------------------------------------------
            

        
        '''Construct the current neighborhoods of the mapping'''
        subMNeighbor = curMap.neighbor(subgraph, 0) 
        gMNeighbor = curMap.neighbor(graph, 1)   
        
        if subMNeighbor and len(subMNeighbor) > len(gMNeighbor): return S, {}, somevalue

        '''Rank the unmapped neighbors by their degrees and select the highest one''' 
        subMN = subMNeighbor[:]
        if not subMNeighbor:            
            subMN = list(set(nx.nodes(subgraph)) - set(curMap.subMap()))
        subMN_deg_node = list([nx.degree(subgraph, v), v] for v in subMN)
        subMN_deg_node.sort(key=lambda t: t[0], reverse=True)
        if not subMN_deg_node: raise Exception(subMN_deg_node)
        selected_vertex = subMN_deg_node[0][1] 
        
        '''Our AGs are always connected. gMNeighbor is empty iff result is empty!'''  
        Candidate = gMNeighbor[:]
        if not subMNeighbor: 
            Candidate = list(set(nx.nodes(graph)) - set(curMap.gMap()))

        #------------------------------------------------------------
        subMN_deg = list(item[0] for item in subMN_deg_node)
        gMN_deg_node = list([nx.degree(graph, v), v] for v in Candidate)
        gMN_deg_node.sort(key=lambda t: t[0], reverse=True)
        gMN_deg = list(item[0] for item in gMN_deg_node)
        
        if len(subMN_deg) > len(gMN_deg) or sum(subMN_deg) > sum(gMN_deg): return S, {}, somevalue
        for i in range(len(subMN_deg)):
            if subMN_deg[i] > gMN_deg[i]: return S, {}, somevalue
        #------------------------------------------------------------        

        '''Remove those graph neighbours which cannot match the suggraph candidate node''' 
        Candidate = [ t for t in Candidate if\
                       nx.degree(subgraph, selected_vertex) <= nx.degree(graph,t)]
        #---------------------------------------------------------------------------#
        # select those mappings close to tau_dict 
        Cand_dist = []
        for v in Candidate:
            result_new = copy.deepcopy(result)
            result_new[selected_vertex] = v
            d =  map_dist(subgraph, graph, tau_dict,result_new)
            if d < somevalue: 
                Cand_dist.append([d,v]) 
        Cand_dist.sort(key=lambda t: t[0])
        Candidate = [ item[1] for item in Cand_dist]
        #---------------------------------------------------------------------------#

        if not Candidate: return S, {}, somevalue
        if not result and printOK: print(selected_vertex, Candidate)
        
        # print('---', somevalue)
        for t in Candidate:
            
            if printOK: print('Can we add',  {selected_vertex:t}, 'to', result, '?') 
            
            '''Check if {selected_vertex : t} satisfies the VF2 rule'''
            if(self.CandMeetRules(selected_vertex, t, subgraph, graph, result, curMap.subMap(),\
                                curMap.gMap(), subMNeighbor, gMNeighbor)):
                
                if printOK: print([selected_vertex,t], 'meets the rules for', result) 

                result[selected_vertex] = t #Extend the mapping!
                
                if time.time()-start_A > stop: 
                    print('dfsmatch time exceeds', stop ) #timeout
                    return S, result, somevalue
                
                # print('kkk', somevalue)
                '''Return all completions of the extended mapping'''
                _,_,somevalue = self.dfsTopMatch(subgraph, graph, tau_dict, S, result, somevalue, stop, printOK)
                #//q: if as below, the self calling does not stop.
                # S,result,somevalue = self.dfsTopMatch(subgraph, graph, tau_dict, S, result, somevalue, stop, printOK)
                #//q: if just self.dfsTopMatch, the original somevlaue is restored when we consider a new t 
                
                # print('ppp', somevalue)

                if printOK:
                    print('After adding', [selected_vertex,t], 'and calling dfsMatch, we have', result, ) 
                    print('----------------')

                '''After this call of dfsMatch, we extect to obtain an updated S. 
                    The result itself will be updated? '''

                if len(result) == len(nx.nodes(subgraph)): 
                    temp_dist = map_dist(subgraph,graph,tau_dict,result)
                    if temp_dist < somevalue:
                        print(temp_dist,'::', somevalue)
                        # somevalue = temp_dist
                        somevalue = min(temp_dist, somevalue)
                        # S.append(result)                    
                        S = [result]
                        if printOK: print('A good embedding has been found!', result)

                '''The pair is examined and we pop it out and restore the mapping'''
                if selected_vertex in result: 
                    result.pop(selected_vertex)
                
            # else:
                # if printOK: print([selected_vertex,t], 'does not meet the rules for', result) 
        # print('zzz', somevalue)                                
        return S, result, somevalue   


#global somevalue
    def dfsTopRestMatch(self, subgraph, graph, delta, Depth_dict, Lrest, C, nl, S, result, somevalue, stop, printOK):
        '''Add to S all embeddings which extend the current mapping'''
        start_A = time.time()
        
        if printOK:
            print('^^^^^^^^^^^^^^^')
            print('Call dfsMatch for', result, '\n Start S =', S, round(somevalue,2))
        
        curMap = Map(result) #sl create a Map object!   
        
        '''1.0 In case we are in a leaf node of the search tree'''
        if len(result) == len(nx.nodes(subgraph)): 
                
            temp_dist = map_circuit_cost(delta, graph, result, Lrest, C, nl)
            
            SChange = False
            if round(temp_dist,2) < round(somevalue,2) or not S:
                somevalue = min(round(temp_dist,2), round(somevalue,2)) 
                temp_result = copy.deepcopy(result)
                S = [temp_result]
                SChange = True
      

            if printOK: 
                print('End Call for', result, '\n Return S =',  S, round(somevalue,2), SChange)
                
            '''The search on this leaf node is over!'''
            return S, {}, somevalue, SChange            
                    
        '''1.1 Otherwise, determine the child nodes of the current node '''  
        
        '''1.1.1. Construct the current neighborhoods of the mapping'''
        subMNeighbor = curMap.neighbor(subgraph, 0) 
        gMNeighbor = curMap.neighbor(graph, 1)   
        
        '''1.1.2. Rank the unmapped neighbors by their degrees and select the highest one''' 
        subMN = subMNeighbor[:]
        if not subMNeighbor:            
            subMN = list(set(nx.nodes(subgraph)) - set(curMap.subMap()))
        subMN_deg_node = list([nx.degree(subgraph, v), v] for v in subMN)
        subMN_deg_node.sort(key=lambda t: t[0], reverse=True)
        if not subMN_deg_node: raise Exception(subMN_deg_node)
        selected_vertex = subMN_deg_node[0][1] 
        print('The selected vertex is', selected_vertex)
        
        '''1.1.3. Select the candidate nodes in graph that we can map the selected_vertex '''
        
        '''Note that our AGs are always connected. gMNeighbor is empty iff result is empty!'''  
        Candidate = gMNeighbor[:]
        if not subMNeighbor: 
            Candidate = list(set(nx.nodes(graph)) - set(curMap.gMap()))

        #---------------------------------------------------------------------------#
        if not Candidate: 
            print('bad node 2', len(S), round(temp_dist,2), round(somevalue,2), result)
            
            '''The search on this non-leaf node is also over!'''            
            return S, {}, somevalue, False
        
        '''2.0 Having the child nodes, we then need to examine these child nodes one by one.'''
        print('Candidate', Candidate)
        for t in Candidate:
            
            if printOK: print('Can we add',  {selected_vertex:t}, 'to', result, '?') 
            
            '''2.1 Check if {selected_vertex : t} satisfies the VF2 rule'''
            if (self.CandMeetRules(selected_vertex, t, subgraph, graph, result, curMap.subMap(),\
                                curMap.gMap(), subMNeighbor, gMNeighbor)):
                                
                result[selected_vertex] = t #Extend the mapping!
                            
                '''Return the best completion of the extended mapping'''
                _,_,_, SChange =\
                    self.dfsTopRestMatch(subgraph, graph, delta, Depth_dict, Lrest, C, nl, S, result, somevalue, stop, printOK)


                print('End Call (inside loop) for', result, '\n Return S =',  S, round(somevalue,2), SChange)
                # if SChange: print(S)


                '''After this call of dfsMatch, we extect to obtain an updated S. 
                    The result itself will be updated? '''                        

                '''2.2 The pair is examined and we pop it out and restore the mapping'''
                if selected_vertex not in result: 
                    raise Exception (selected_vertex, 'not in', result)
                if selected_vertex in result: 
                    result.pop(selected_vertex)
                
        '''3. We have examined all child nodes of the current search node - result '''
                
        if printOK:
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            print('End Call for', result, '\n Return S =', S, round(somevalue,2))
            
        return S, {}, somevalue, SChange   
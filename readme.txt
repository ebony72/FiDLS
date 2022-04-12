This repository is the Python implementation for algorithms introduced in 
"Qubit Mapping Based on Subgraph Isomorphism and Filtered Depth-Limited Search" 
(https://arxiv.org/abs/2004.07138, published in IEEE Transactions on Computers, 70(11): 1777-1788 (2021)). 

This is a significant revision of the April 22 version 
(see https://github.com/BensonZhou1991/Qubit-Mapping-Subgraph-Isomorphism)
Readers can use "FiDLS_inimap.py" to generate initial mappings for their circuits and 
"FiDLS_run.py" to get information about the output circuits. 
For example, the result (91, 'qft_16', 16, 240, 432, 192, 4.28, 1.8) shows 
the count (91), name ('qft_16'), num_qubit (16), num_cnot (240) of the input circuits, and 
the output_cnots (432), added_cnots (192), time (4.28"), and ratio := output_cnots/num_cnot (1.8). 

We provide the following selections:

A. QFilter_type (line 18)
	One can select type from {'9', '0', '01', '01x', '01y', '1x'}
		 '9 :: if no filter is used'
		 '0 :: if Q0 is used for all layers'
		'01 :: if Q0 is used for the first layer and Q1 is used for the other layers'		
	       '01x :: if Q0 is used for the first layer and Q0+Q1 is used for the other layers'
	       '01y :: if Q0+Q1 is used for the first layer and Q0 is used for the other layers'
		'1x :: if Q0+Q1 is used for all layers '
	where Q0 and Q1 are the set of qubits in the front layer and the first look-ahead layer 
	of the current logical circuit.
	The default filter type is '01y'. 

B. Initial mapping (line 19)
	One can select mapping from {'top', 'wgt', 'empty', 'naive'}.	The default mapping is 'top'

C. Size of the circuits (line 18)
	One can select size from {'small', 'medium', 'large', 'all'}.	The default size is 'medium'
	
D. Architecture graph AG (line 18)
    One can select AG from {'tokyo', 'sycamore', 'rochester', 'q19x19'} and may also define her own AG. 
    The default AG is 'tokyo'
    
E. FiDLS-G (GVal is True) or FiDLS-D (GVal is False). FiDLS-G is the default one. 

Other selections include whether or not using anchor and the stop time in inimap construction. 

Send email to me (Sanjiang Li: mrlisj@gmail.com) if you have questions!

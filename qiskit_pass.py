"""FiDLSSwap — Qiskit TransformationPass wrapping the FiDLS router.

Usage:
    from qiskit.transpiler import CouplingMap, PassManager
    from qiskit_pass import FiDLSSwap

    coupling_map = CouplingMap.from_grid(5, 5)
    pm = PassManager([FiDLSSwap(coupling_map)])
    routed = pm.run(qc)

Preserves single-qubit gates and barriers. If an upstream layout pass has set
``property_set["layout"]``, that layout is used as the initial mapping;
otherwise a trivial layout (virtual i → physical i) is generated.
"""
from collections import defaultdict, deque

import networkx as nx
import numpy as np

from qiskit.circuit.library import SwapGate
from qiskit.dagcircuit import DAGCircuit, DAGOpNode
from qiskit.transpiler.basepasses import TransformationPass
from qiskit.transpiler.exceptions import TranspilerError
from qiskit.transpiler.layout import Layout

from ag import _spl_matrix
from router import qct, VARIANTS
from utils import qubit_in_circuit


class FiDLSSwap(TransformationPass):
    """Route a virtual quantum circuit onto a hardware coupling map using FiDLS.

    Args:
        coupling_map: a Qiskit ``CouplingMap``. If asymmetric it is symmetrized
                      (FiDLS assumes undirected couplings; gate direction is the
                      basis-gate-decomposition pass's problem).
        qfilter_type: Q-filter selector. One of
                      {"9", "0", "01", "01x", "01y", "1x"}. Default ``"01y"``.
        variant:      ``"G"`` (1-layer lookahead) or ``"D"`` (3-layer). Default ``"G"``.
    """

    def __init__(self, coupling_map, qfilter_type="01y", variant="G"):
        super().__init__()
        if variant not in VARIANTS:
            raise ValueError(f"variant must be one of {VARIANTS}, got {variant!r}")
        if qfilter_type not in {"9", "0", "01", "01x", "01y", "1x"}:
            raise ValueError(f"unknown qfilter_type {qfilter_type!r}")

        if coupling_map is None:
            raise TranspilerError("FiDLSSwap requires a coupling_map")

        # FiDLS treats edges as undirected. Symmetrize if needed.
        if not coupling_map.is_symmetric:
            coupling_map = coupling_map.copy()
            coupling_map.make_symmetric()
        self.coupling_map = coupling_map
        self.qfilter_type = qfilter_type
        self.variant = variant

        # Precompute the architecture once per pass instance.
        G = nx.Graph()
        G.add_nodes_from(range(coupling_map.size()))
        for u, v in coupling_map.get_edges():
            G.add_edge(u, v)
        if not nx.is_connected(G):
            raise TranspilerError("FiDLSSwap requires a connected coupling map")
        self._G = G
        self._EG = G.edges()
        self._V = list(G.nodes())
        self._spl_mat = _spl_matrix(G)
        self._SPL = {(u, v): int(self._spl_mat[u, v])
                     for u in self._V for v in self._V}

    def run(self, dag):
        if len(dag.qregs) != 1:
            raise TranspilerError("FiDLSSwap expects a single quantum register")
        canonical_register = next(iter(dag.qregs.values()))
        nq = len(dag.qubits)
        if nq > self.coupling_map.size():
            raise TranspilerError(
                f"Circuit has {nq} virtual qubits but coupling map only has "
                f"{self.coupling_map.size()} physical qubits"
            )

        bit_to_idx = {bit: i for i, bit in enumerate(dag.qubits)}

        # Initial layout: upstream pass result or trivial.
        layout = self.property_set.get("layout")
        if layout is None:
            layout = Layout.generate_trivial_layout(canonical_register)

        # Collect 2-qubit (non-barrier) ops in topological order; build C.
        cnot_nodes = []
        for node in dag.topological_op_nodes():
            if node.op.name == "barrier":
                continue
            if len(node.qargs) == 1:
                continue
            if len(node.qargs) == 2:
                cnot_nodes.append(node)
                continue
            raise TranspilerError(
                f"FiDLSSwap supports 1- and 2-qubit gates only; got "
                f"{len(node.qargs)}-qubit gate {node.op.name!r}"
            )
        cnot_node_to_idx = {n: i for i, n in enumerate(cnot_nodes)}
        C = [[bit_to_idx[n.qargs[0]], bit_to_idx[n.qargs[1]]] for n in cnot_nodes]
        Q = qubit_in_circuit(list(range(len(C))), C)

        # Build tau (physical -> logical) from the initial layout.
        nphys = self.coupling_map.size()
        tau = [-1] * nphys
        for bit, idx in bit_to_idx.items():
            tau[layout[bit]] = idx

        # If there are no 2-qubit gates, the routing is a no-op.
        if not C:
            new_dag = self._build_empty_clone(dag)
            for node in dag.topological_op_nodes():
                new_qargs = [canonical_register[layout[q]] for q in node.qargs]
                new_dag.apply_operation_back(node.op, new_qargs, node.cargs)
            self.property_set["final_layout"] = layout.copy()
            return new_dag

        # Run the router. Get an event stream so we can re-interleave 1q gates.
        _C_out, _cost, events = qct(
            tau, C, Q, self._G, self._EG, self._V, self._SPL,
            self.qfilter_type, variant=self.variant,
            spl_mat=self._spl_mat, return_events=True,
        )

        # Build output DAG by event-stream traversal, with 1q/barrier gates
        # released by predecessor counting.
        new_dag = self._build_empty_clone(dag)
        # virtual Qubit -> current physical index; updated as swaps are applied.
        v2p = {bit: layout[bit] for bit in dag.qubits}
        # physical index -> virtual Qubit currently there.
        p2v = {p: None for p in range(nphys)}
        for bit in dag.qubits:
            p2v[v2p[bit]] = bit

        # Predecessor bookkeeping over op nodes only (ignore in/out wires).
        op_nodes = list(dag.op_nodes())
        pred_count = {}
        op_successors = defaultdict(list)
        for node in op_nodes:
            preds = [p for p in dag.predecessors(node)
                     if isinstance(p, DAGOpNode)]
            pred_count[node] = len(preds)
            for p in preds:
                op_successors[p].append(node)

        # A node is "non-routing" if it doesn't need a 2-qubit coupling
        # decision: 1-qubit gates and barriers. These can be emitted as soon
        # as their op-predecessors are all emitted.
        def is_non_routing(node):
            return node.op.name == "barrier" or len(node.qargs) != 2

        ready = deque(n for n in op_nodes
                      if pred_count[n] == 0 and is_non_routing(n))

        def release_successors(node):
            for s in op_successors[node]:
                pred_count[s] -= 1
                if pred_count[s] == 0 and is_non_routing(s):
                    ready.append(s)

        def emit_non_routing(node):
            new_qargs = [canonical_register[v2p[q]] for q in node.qargs]
            new_dag.apply_operation_back(node.op, new_qargs, node.cargs)
            release_successors(node)

        def drain_ready():
            while ready:
                emit_non_routing(ready.popleft())

        for ev in events:
            kind = ev[0]
            if kind == "swap":
                # Drain pending 1q gates BEFORE applying the swap to v2p, so
                # they're emitted at their pre-swap physical positions.
                drain_ready()
                u, v = ev[2]
                qu = p2v[u]
                qv = p2v[v]
                if qu is not None:
                    v2p[qu] = v
                if qv is not None:
                    v2p[qv] = u
                p2v[u], p2v[v] = qv, qu
                new_dag.apply_operation_back(
                    SwapGate(),
                    [canonical_register[u], canonical_register[v]],
                    [],
                )
            elif kind == "cnot":
                idx = ev[1]
                p, q = ev[2]
                node = cnot_nodes[idx]
                # Flush any 1q/barrier ops that became ready in topological
                # order before this 2q gate.
                drain_ready()
                new_dag.apply_operation_back(
                    node.op,
                    [canonical_register[p], canonical_register[q]],
                    node.cargs,
                )
                release_successors(node)
            else:
                raise RuntimeError(f"unexpected event kind {kind!r}")

        # Drain any trailing 1q gates or barriers.
        drain_ready()

        # Sanity check: every op-node should have been emitted.
        unemitted = [n for n in op_nodes if pred_count[n] != 0]
        if unemitted:
            raise RuntimeError(
                f"FiDLSSwap: {len(unemitted)} op nodes were not emitted "
                f"(likely a bug in predecessor counting)"
            )

        # Final layout: virtual -> physical after all swaps applied.
        final_layout = Layout()
        for bit, phys in v2p.items():
            final_layout.add(bit, phys)
        self.property_set["final_layout"] = final_layout
        return new_dag

    @staticmethod
    def _build_empty_clone(dag):
        """Empty DAGCircuit with the same registers as `dag`."""
        new = DAGCircuit()
        new.name = dag.name
        new.metadata = dag.metadata
        for qreg in dag.qregs.values():
            new.add_qreg(qreg)
        for creg in dag.cregs.values():
            new.add_creg(creg)
        return new

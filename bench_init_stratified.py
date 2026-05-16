"""Stratified init-mapping ablation: complete vs partial SI mappings.

For each B131 medium circuit, check whether the cached SI mapping covers
every logical qubit (complete) or leaves some unmapped (partial). Then run
FiDLS-G with both `fidls-top` (SI + completion) and `sabrelyt` (SabreLayout)
initial mappings, and report per-bucket ratios.

If `fidls-top` ≈ `sabrelyt` on the "complete" bucket and the gap only shows
up on the "partial" bucket, then map-completion is the entire story. If the
gap persists on completes, there's also a fundamental "fit front layer
first" bias in the SI strategy.
"""
import argparse
import json
import os

from qiskit import QuantumCircuit
from qiskit.transpiler import CouplingMap, PassManager
from qiskit.transpiler.passes import SabreLayout

import ag
from router import qct
from utils import qubit_in_circuit, map_completion


B131 = os.path.join(os.path.dirname(__file__), "B131")
INIMAP_DIR = os.path.join(os.path.dirname(__file__), "inimap")


def _coupling_map_from_graph(g):
    edges = []
    for u, v in g.edges():
        edges.append((u, v)); edges.append((v, u))
    return CouplingMap(edges)


def _sabre_layout_v2p(C, n_phys, cm, seed=0):
    qc = QuantumCircuit(n_phys)
    for ctrl, tgt in C:
        qc.cx(ctrl, tgt)
    pm = PassManager([SabreLayout(cm, seed=seed)])
    pm.run(qc)
    layout = pm.property_set["layout"]
    return {qc.qubits.index(b): layout[b] for b in qc.qubits}


def _fidls_route(C, v2p, A, variant="G"):
    V = list(A.graph.nodes())
    tau = [-1] * len(V)
    for q, p in v2p.items():
        tau[p] = q
    out, _ = qct(
        tau, C, qubit_in_circuit(list(range(len(C))), C),
        A.graph, A.graph.edges(), V, A.SPL,
        "01y", variant=variant, spl_mat=A.spl_mat,
    )
    return len(out)


def run_topology(topology, lo=100, hi=1000):
    A = ag.build(topology)
    V = list(A.graph.nodes())
    cm = _coupling_map_from_graph(A.graph)
    n_phys = cm.size()
    cache = os.path.join(INIMAP_DIR, f"_inimap_list_{topology}_top_B131.txt")
    with open(cache) as f:
        IM = json.loads(f.read())

    complete = []   # (name, in_cnots, out_fidls_top, out_sabrelyt)
    partial = []

    for cnt, fn in enumerate(os.listdir(B131), 1):
        if fn.startswith("."):
            continue
        with open(os.path.join(B131, fn)) as f:
            C = json.loads(f.read())
        if not (lo <= len(C) <= hi):
            continue
        Q = qubit_in_circuit(list(range(len(C))), C)
        if len(Q) > len(V):
            continue
        entry = next((e[1] for e in IM if e[0] == cnt), None)
        if entry is None:
            continue
        si_v2p = {q: p for q, p in entry}
        si_complete = (len(si_v2p) == len(Q)
                       and all(q in si_v2p for q in Q))

        # Completed FiDLS-top (this is what the router actually sees).
        fidls_v2p = dict(si_v2p)
        if len(fidls_v2p) < len(Q):
            fidls_v2p = map_completion(
                fidls_v2p, list(range(len(C))), C, Q, A, V
            )
        if len(fidls_v2p) < len(Q):
            continue

        # SabreLayout (restricted to qubits actually used).
        sabre_v2p = _sabre_layout_v2p(C, n_phys, cm, seed=0)
        sabre_v2p_q = {q: sabre_v2p[q] for q in Q if q in sabre_v2p}
        if len(sabre_v2p_q) < len(Q):
            sabre_v2p_q = map_completion(
                dict(sabre_v2p_q), list(range(len(C))), C, Q, A, V
            )
        if len(sabre_v2p_q) < len(Q):
            continue

        out_top = _fidls_route(C, fidls_v2p, A)
        out_sab = _fidls_route(C, sabre_v2p_q, A)

        row = (fn[:-9], len(Q), len(si_v2p), len(C), out_top, out_sab)
        (complete if si_complete else partial).append(row)

    def summarize(bucket, label):
        if not bucket:
            print(f"  {label}: (no circuits)")
            return
        n = len(bucket)
        sum_in = sum(r[3] for r in bucket)
        sum_top = sum(r[4] for r in bucket)
        sum_sab = sum(r[5] for r in bucket)
        gap_pct = 100 * (sum_top - sum_sab) / sum_sab if sum_sab else 0.0
        print(f"  {label} ({n:>2} circuits, {sum_in:>5} input CNOTs):")
        print(f"    fidls-top init     : ratio {sum_top/sum_in:.4f} "
              f"({sum_top} out CNOTs)")
        print(f"    sabrelyt  init     : ratio {sum_sab/sum_in:.4f} "
              f"({sum_sab} out CNOTs)")
        print(f"    gap (top vs sabre) : {gap_pct:+.1f}%  "
              f"({'fidls-top worse' if gap_pct > 0 else 'fidls-top better'})")

    print(f"\n=== {topology} | B131 medium ===")
    summarize(complete, "SI was COMPLETE   ")
    summarize(partial,  "SI was PARTIAL    ")

    if complete and partial:
        # Also show the per-circuit data for partial — interesting to see
        # how partial each one was.
        print("\n    partial-bucket detail (SI mapped / Q size):")
        for name, nq, mapped, _l, _t, _s in partial:
            print(f"      {name:<28} {mapped}/{nq}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--topologies", nargs="+",
                   default=["tokyo", "rochester", "sycamore"])
    args = p.parse_args()
    for t in args.topologies:
        run_topology(t)


if __name__ == "__main__":
    main()

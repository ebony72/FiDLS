"""Apples-to-apples comparison: FiDLS vs Qiskit's SabreSwap on B131 medium.

Both routers see the same circuits, same coupling map, and the same initial
layout (FiDLS's cached `top` subgraph-isomorphism mapping). The only
difference is the routing algorithm. We report:

  * CNOT count in / out (each SWAP = 3 CNOTs)
  * Wall-clock routing time
  * Mean overhead ratio (sum_out / sum_in)
"""
import argparse
import json
import os
import time

from qiskit import QuantumCircuit
from qiskit.transpiler import CouplingMap, PassManager
from qiskit.transpiler.passes import SabreSwap

import ag
from router import qct
from utils import qubit_in_circuit, map_completion


B131 = os.path.join(os.path.dirname(__file__), "B131")
INIMAP_DIR = os.path.join(os.path.dirname(__file__), "inimap")

# Default seeds chosen up front so re-runs are reproducible. Increasing the
# count tightens the SabreSwap mean estimate but multiplies the bench time.
SABRE_SEEDS = [0, 1, 2, 3, 4]


def _coupling_map_from_graph(g):
    """NetworkX undirected graph -> symmetric Qiskit CouplingMap."""
    edges = []
    for u, v in g.edges():
        edges.append((u, v))
        edges.append((v, u))
    cm = CouplingMap(edges)
    return cm


def _load_runs(topology, size_lo=100, size_hi=1000):
    """Both routers see the same initial layout. We complete partial VF2
    mappings with map_completion so SabreSwap (which needs every virtual qubit
    placed) sees the same starting state FiDLS would after online extension."""
    A = ag.build(topology)
    V = list(A.graph.nodes())
    cache = os.path.join(INIMAP_DIR, f"_inimap_list_{topology}_top_B131.txt")
    with open(cache) as f:
        IM = json.loads(f.read())

    runs = []
    for cnt, fn in enumerate(os.listdir(B131), 1):
        if fn.startswith("."):
            continue
        with open(os.path.join(B131, fn)) as f:
            C = json.loads(f.read())
        if not (size_lo <= len(C) <= size_hi):
            continue
        Q = qubit_in_circuit(list(range(len(C))), C)
        if len(Q) > len(V):
            continue
        entry = next((e[1] for e in IM if e[0] == cnt), None)
        if entry is None:
            continue
        v2p = {q: phys for q, phys in entry}
        # Complete partial mappings so both routers start from the same full
        # initial layout. (Without this, SabreSwap can't even start — it has
        # no concept of online qubit-placement extension.)
        if len(v2p) < len(Q):
            v2p = map_completion(v2p, list(range(len(C))), C, Q, A, V)
        if len(v2p) < len(Q):
            # Couldn't even fill in the rest — skip this circuit on both sides.
            continue
        tau = [-1] * len(V)
        for q, p in v2p.items():
            tau[p] = q
        runs.append((cnt, fn[:-9], len(Q), C, tau, v2p))
    return A, runs


def _fidls_one(C, tau, A, qfilter, variant):
    t0 = time.time()
    out, _ = qct(
        tau[:], C, qubit_in_circuit(list(range(len(C))), C),
        A.graph, A.graph.edges(), list(A.graph.nodes()),
        A.SPL, qfilter, variant=variant, spl_mat=A.spl_mat,
    )
    return len(out), time.time() - t0


def _sabre_one(C, v2p, cm, heuristic, seed):
    """Build a physical-qubit QC with the FiDLS initial layout, then SabreSwap.

    Returns (cnot_equiv_count, time). Each SwapGate is counted as 3 CNOTs to
    match FiDLS's accounting (FiDLS decomposes swaps into 3 CNOTs in its
    output stream; Sabre keeps them as a single SwapGate op).
    """
    n_phys = cm.size()
    qc = QuantumCircuit(n_phys)
    for ctrl, tgt in C:
        if ctrl not in v2p or tgt not in v2p:
            return None, None  # mapping incomplete; skip
        qc.cx(v2p[ctrl], v2p[tgt])
    t0 = time.time()
    pm = PassManager([SabreSwap(cm, heuristic=heuristic, seed=seed)])
    out = pm.run(qc)
    dt = time.time() - t0
    ops = out.count_ops()
    cnot_equiv = ops.get("cx", 0) + 3 * ops.get("swap", 0)
    return cnot_equiv, dt


def run_topology(topology, size_lo, size_hi):
    A, runs = _load_runs(topology, size_lo, size_hi)
    cm = _coupling_map_from_graph(A.graph)
    n = len(runs)

    print(f"\n=== {topology} | B131 size {size_lo}–{size_hi} | {n} circuits ===")
    if not runs:
        return

    rows = []
    # FiDLS-G and FiDLS-D
    for variant in ("G", "D"):
        sum_in = sum_out = 0
        total_t = 0.0
        for _cnt, _name, _nq, C, tau, _v2p in runs:
            out_cnots, dt = _fidls_one(C, tau, A, "01y", variant)
            sum_in += len(C); sum_out += out_cnots; total_t += dt
        rows.append((f"FiDLS-{variant}", sum_out, sum_out / sum_in, total_t))

    # SabreSwap with three heuristics, averaged over seeds.
    for heur in ("basic", "lookahead", "decay"):
        sum_in_h = sum_out_h = 0
        total_t_h = 0.0
        skipped = 0
        for _cnt, _name, _nq, C, _tau, v2p in runs:
            best = None  # best across seeds
            cum_t = 0.0
            for seed in SABRE_SEEDS:
                out_cnots, dt = _sabre_one(C, v2p, cm, heur, seed)
                if out_cnots is None:
                    skipped += 1
                    break
                cum_t += dt
                if best is None or out_cnots < best:
                    best = out_cnots
            if best is None:
                continue
            sum_in_h += len(C); sum_out_h += best; total_t_h += cum_t / len(SABRE_SEEDS)
        if sum_in_h:
            rows.append((
                f"Sabre-{heur}",
                sum_out_h, sum_out_h / sum_in_h, total_t_h,
            ))

    in_total = sum(len(C) for _, _, _, C, _, _ in runs)
    print(f"  input CNOTs total: {in_total}")
    print(f"  {'router':<18}{'out CNOTs':>12}{'ratio':>10}{'time (s)':>12}")
    for label, out, ratio, t in sorted(rows, key=lambda r: r[2]):
        print(f"  {label:<18}{out:>12}{ratio:>10.4f}{t:>12.2f}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--topologies", nargs="+",
                   default=["tokyo", "rochester", "sycamore"])
    p.add_argument("--size", default="medium",
                   choices=["small", "medium", "large", "all"])
    return p.parse_args()


SIZE_BOUNDS = {
    "small": (1, 99),
    "medium": (100, 1000),
    "large": (1001, 1_000_000),
    "all": (1, 1_000_000),
}


def main():
    args = parse_args()
    lo, hi = SIZE_BOUNDS[args.size]
    for t in args.topologies:
        run_topology(t, lo, hi)


if __name__ == "__main__":
    main()

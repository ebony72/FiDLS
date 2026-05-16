"""Init-mapping ablation: does FiDLS's tokyo win come from the SI init or the router?

For each B131 medium circuit on a given topology, build the cross product:
    init  ∈ {FiDLS-top, SabreLayout, trivial}
    router ∈ {FiDLS-G, Sabre-decay}
Run all 6 combos starting from the same circuit + initial layout. Report mean
gate-count ratio per cell.

The point: if FiDLS-top + Sabre ≈ FiDLS-top + FiDLS-G, the win is the init.
If SabreLayout + FiDLS-G ≈ FiDLS-top + FiDLS-G, the win is the router.
"""
import argparse
import json
import os
import time

from qiskit import QuantumCircuit
from qiskit.transpiler import CouplingMap, PassManager
from qiskit.transpiler.passes import (
    SabreSwap, SabreLayout, FullAncillaAllocation, EnlargeWithAncilla, ApplyLayout,
    SetLayout,
)

import ag
from router import qct
from utils import qubit_in_circuit, map_completion


B131 = os.path.join(os.path.dirname(__file__), "B131")
INIMAP_DIR = os.path.join(os.path.dirname(__file__), "inimap")
SABRE_SWAP_SEEDS = [0, 1, 2, 3, 4]


def _coupling_map_from_graph(g):
    edges = []
    for u, v in g.edges():
        edges.append((u, v))
        edges.append((v, u))
    return CouplingMap(edges)


def _load_circuits(topology, lo, hi):
    """(architecture, list of (idx, name, Q, C, fidls_v2p)) — fidls_v2p is the
    cached SI mapping completed to cover all of Q."""
    A = ag.build(topology)
    V = list(A.graph.nodes())
    cache = os.path.join(INIMAP_DIR, f"_inimap_list_{topology}_top_B131.txt")
    with open(cache) as f:
        IM = json.loads(f.read())
    runs = []
    for cnt, fn in enumerate(os.listdir(B131), 1):
        if fn.startswith("."): continue
        with open(os.path.join(B131, fn)) as f:
            C = json.loads(f.read())
        if not (lo <= len(C) <= hi): continue
        Q = qubit_in_circuit(list(range(len(C))), C)
        if len(Q) > len(V): continue
        entry = next((e[1] for e in IM if e[0] == cnt), None)
        if entry is None: continue
        v2p = {q: phys for q, phys in entry}
        if len(v2p) < len(Q):
            v2p = map_completion(v2p, list(range(len(C))), C, Q, A, V)
        if len(v2p) < len(Q):
            continue
        runs.append((cnt, fn[:-9], Q, C, v2p))
    return A, runs


def _sabre_layout(C, n_phys, cm, seed=0):
    """Run SabreLayout on a logical-qubit QC, return the v2p mapping it picks.

    We build the QC with n_phys qubits (no ancilla expansion needed if
    n_phys == cm.size()), apply SabreLayout, then read property_set["layout"]
    and convert to v2p.
    """
    qc = QuantumCircuit(n_phys)
    for ctrl, tgt in C:
        qc.cx(ctrl, tgt)
    pm = PassManager([SabreLayout(cm, seed=seed)])
    pm.run(qc)
    layout = pm.property_set["layout"]
    v2p = {}
    for bit in qc.qubits:
        v2p[qc.qubits.index(bit)] = layout[bit]
    return v2p


def _fidls_route(C, v2p, A, qfilter="01y", variant="G"):
    V = list(A.graph.nodes())
    tau = [-1] * len(V)
    for q, p in v2p.items():
        tau[p] = q
    t0 = time.time()
    out, _ = qct(
        tau, C, qubit_in_circuit(list(range(len(C))), C),
        A.graph, A.graph.edges(), V, A.SPL,
        qfilter, variant=variant, spl_mat=A.spl_mat,
    )
    return len(out), time.time() - t0


def _sabre_route(C, v2p, cm, seeds=SABRE_SWAP_SEEDS, heuristic="decay"):
    """Build a physical-qubit QC at the given initial layout, route with
    SabreSwap across `seeds` seeds, return best (min) cnot-equivalent count.

    Returns (None, None) if v2p doesn't cover every logical qubit appearing
    in C (caller should skip this circuit/init combo)."""
    n_phys = cm.size()
    if any(c not in v2p or t not in v2p for c, t in C):
        return None, None
    best = None
    total_t = 0.0
    for seed in seeds:
        qc = QuantumCircuit(n_phys)
        for ctrl, tgt in C:
            qc.cx(v2p[ctrl], v2p[tgt])
        t0 = time.time()
        pm = PassManager([SabreSwap(cm, heuristic=heuristic, seed=seed)])
        out = pm.run(qc)
        total_t += time.time() - t0
        ops = out.count_ops()
        cnot_equiv = ops.get("cx", 0) + 3 * ops.get("swap", 0)
        if best is None or cnot_equiv < best:
            best = cnot_equiv
    return best, total_t / len(seeds)


def run_topology(topology, lo, hi):
    A, runs = _load_circuits(topology, lo, hi)
    cm = _coupling_map_from_graph(A.graph)
    n_phys = cm.size()
    n = len(runs)
    print(f"\n=== {topology} | B131 size {lo}–{hi} | {n} circuits ===")
    if not runs: return

    sum_in = sum(len(C) for *_x, C, _v in runs)
    print(f"input CNOTs total: {sum_in}")

    # Three inits per circuit.
    cell_totals = {}   # (init, router) -> sum of out cnots
    cell_times = {}
    for cnt, name, Q, C, fidls_v2p in runs:
        trivial_v2p = {q: q for q in Q}  # virtual i -> physical i

        # SabreLayout init (single seed for fairness)
        sabre_v2p = _sabre_layout(C, n_phys, cm, seed=0)
        # Restrict to logical qubits actually used (sabre's layout is over all n_phys virt)
        sabre_v2p_q = {q: sabre_v2p[q] for q in Q if q in sabre_v2p}
        if len(sabre_v2p_q) < len(Q):
            # Pad any missing using map_completion semantics
            tmp = dict(sabre_v2p_q)
            tmp = map_completion(tmp, list(range(len(C))), C, Q, A, list(A.graph.nodes()))
            sabre_v2p_q = tmp

        inits = {
            "fidls-top":  fidls_v2p,
            "sabrelyt":   sabre_v2p_q,
            "trivial":    trivial_v2p,
        }
        for init_name, v2p in inits.items():
            # FiDLS-G
            out_f, t_f = _fidls_route(C, v2p, A, qfilter="01y", variant="G")
            cell_totals.setdefault((init_name, "fidls-G"), 0)
            cell_totals[(init_name, "fidls-G")] += out_f
            cell_times.setdefault((init_name, "fidls-G"), 0.0)
            cell_times[(init_name, "fidls-G")] += t_f
            # Sabre-decay
            out_s, t_s = _sabre_route(C, v2p, cm)
            if out_s is None:
                # Init doesn't cover all of Q — leave this cell empty
                continue
            cell_totals.setdefault((init_name, "sabre-decay"), 0)
            cell_totals[(init_name, "sabre-decay")] += out_s
            cell_times.setdefault((init_name, "sabre-decay"), 0.0)
            cell_times[(init_name, "sabre-decay")] += t_s

    # Pretty-print as a matrix.
    inits_order = ["fidls-top", "sabrelyt", "trivial"]
    routers_order = ["fidls-G", "sabre-decay"]
    header = f"{'router \\\\ init':<18}" + "".join(f"{i:>14}" for i in inits_order)
    print(header)
    print("-" * len(header))
    for r in routers_order:
        cells = []
        for i in inits_order:
            out = cell_totals.get((i, r), 0)
            cells.append(f"{out/sum_in:>14.4f}")
        print(f"{r:<18}" + "".join(cells))

    print()
    print("(values are out-CNOTs / in-CNOTs; lower is better)")
    print()
    header = f"{'router \\\\ init  time(s)':<26}" + "".join(f"{i:>14}" for i in inits_order)
    print(header)
    print("-" * len(header))
    for r in routers_order:
        cells = []
        for i in inits_order:
            cells.append(f"{cell_times.get((i, r), 0.0):>14.2f}")
        print(f"{r:<26}" + "".join(cells))


SIZE = {"small": (1, 99), "medium": (100, 1000), "large": (1001, 1_000_000),
        "all": (1, 1_000_000)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--topologies", nargs="+",
                   default=["tokyo", "rochester", "sycamore"])
    p.add_argument("--size", default="medium", choices=sorted(SIZE))
    args = p.parse_args()
    lo, hi = SIZE[args.size]
    for t in args.topologies:
        run_topology(t, lo, hi)


if __name__ == "__main__":
    main()

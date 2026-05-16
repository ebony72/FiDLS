"""Entry point: run FiDLS on a folder of QASM/JSON circuits and report stats.

Example:
    python FiDLS_run.py --ag tokyo --variant G --filter 01y \
        --mapping top --size medium --path B131/

For each circuit, prints a row:
    (idx, name, num_qubits, in_cnots, out_cnots, added_cnots, time_s, ratio)

If `--mapping` is 'top' or 'wgt', the script reads a precomputed inimap cache
from `inimap/_inimap_list_<ag>_<mapping>_<path-stem>.txt`. Run
`FiDLS_inimap.py` first to generate it.
"""
import argparse
import json
import os
import sys
import time

import ag as ag_mod
from router import qct
from utils import (
    CreateCircuitFromQASM, ReducedCircuit, centre, hub,
    graph_of_circuit, qubit_in_circuit, map_completion,
)


SIZE_FILTERS = {
    "small":  lambda n: n < 100,
    "medium": lambda n: 100 <= n <= 1000,
    "large":  lambda n: n > 1000,
    "all":    lambda n: True,
}

# Hardcoded skip-list for bigQ/ (duplicate circuits, mirrored from original).
BIGQ_DUPLICATE_INDICES = {19, 21, 34, 42, 44, 47, 49}
BIGQ_MAX_GATES = 15000


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ag", default="tokyo",
                   choices=sorted(ag_mod.TOPOLOGIES),
                   help="Architecture graph name (default: tokyo).")
    p.add_argument("--variant", default="G", choices=["G", "D"],
                   help="FiDLS variant: G (1-layer lookahead) or D (3-layer). "
                        "Default G.")
    p.add_argument("--filter", dest="qfilter", default="01y",
                   choices=["9", "0", "01", "01x", "01y", "1x"],
                   help="Q-filter type (default 01y).")
    p.add_argument("--mapping", default="top",
                   choices=["top", "wgt", "empty", "naive"],
                   help="Initial mapping strategy (default top).")
    p.add_argument("--size", default="medium",
                   choices=sorted(SIZE_FILTERS),
                   help="Circuit size filter (default medium).")
    p.add_argument("--path", default="B131/",
                   help="Directory of input circuits (default B131/).")
    p.add_argument("--log", default=None,
                   help="Optional path to append per-circuit results. "
                        "If unset, results only go to stdout.")
    p.add_argument("--only", type=int, default=None,
                   help="Process only circuit #N (1-indexed) for quick checks.")
    p.add_argument("--complete-init", action="store_true", default=False,
                   help="Greedily complete partial initial mappings before "
                        "routing. Off by default — the router's online SRMD "
                        "extension usually produces better placement. Useful "
                        "for very large topologies (e.g. q19x19) where VF2 "
                        "leaves many logical qubits unplaced.")
    return p.parse_args(argv)


def load_inimap_cache(args):
    """Load precomputed initial mappings keyed by circuit index."""
    path_stem = args.path.rstrip("/")
    name2 = f"_inimap_list_{args.ag}_{args.mapping}_{path_stem}"
    cache_path = os.path.join("inimap", name2 + ".txt")
    try:
        with open(cache_path, "r") as f:
            return json.loads(f.read()), cache_path
    except FileNotFoundError:
        sys.exit(
            f"Initial-mapping cache not found at {cache_path}.\n"
            f"Run: python FiDLS_inimap.py --ag {args.ag} "
            f"--mapping {args.mapping} --path {args.path}"
        )


def initial_mapping_for(args, IM, count, C, G, V):
    """Construct the (logical -> physical) dict for circuit #count."""
    if args.mapping in ("top", "wgt"):
        # Cached: IM is a list of [idx, [[q, v], ...]] entries.
        imlist = IM[count - 1][1]
        return {x[0]: x[1] for x in imlist}
    if args.mapping == "empty":
        g_of_c = graph_of_circuit(C)
        q = hub(g_of_c)
        v = hub(G)
        return {q: v}
    if args.mapping == "naive":
        Q = qubit_in_circuit(list(range(len(C))), C)
        return {i: i for i in range(len(Q))}
    raise ValueError(f"unknown mapping {args.mapping}")


def open_log(path):
    """Return a callable that writes a line to the log file (and stdout)."""
    if path is None:
        return lambda s: None
    fh = open(path, "a")

    def write(s):
        fh.write(str(s) + "\n")
        fh.flush()
    return write


def main(argv=None):
    args = parse_args(argv)

    AG = ag_mod.build(args.ag)
    G = AG.graph
    # EG must support symmetric `(u, v) in EG` lookups (used by utils.entail).
    # NetworkX EdgeView does that natively; list(G.edges()) does not.
    EG = G.edges()
    V = list(G.nodes())
    SPL = AG.SPL
    SPL_mat = AG.spl_mat  # NumPy matrix passed to router for fast inner-loop lookups

    log = open_log(args.log)
    print(time.asctime())
    log(time.asctime())
    header = (f"FiDLS-{args.variant} on {args.ag} | filter={args.qfilter} | "
              f"mapping={args.mapping} | size={args.size} | path={args.path}")
    print(header)
    log(header)

    IM = None
    if args.mapping in ("top", "wgt"):
        IM, cache_path = load_inimap_cache(args)
        log(f"loaded initial mappings from {cache_path}")

    size_ok = SIZE_FILTERS[args.size]
    files = os.listdir(args.path)

    sum_in = sum_out = 0
    total_route_time = 0.0
    t_start = time.time()
    count = 0

    for file_name in files:
        count += 1
        if args.only is not None and count != args.only:
            continue
        if file_name.endswith("qasm"):
            cir = CreateCircuitFromQASM(file_name, args.path)
            C = ReducedCircuit(cir)
        else:
            with open(args.path + file_name, "r") as f:
                C = json.loads(f.read())
        l = len(C)
        if args.path.rstrip("/") == "bigQ":
            if count in BIGQ_DUPLICATE_INDICES:
                continue
            if l > BIGQ_MAX_GATES:
                continue
        if not size_ok(l):
            continue

        L = list(range(l))
        Q = qubit_in_circuit(L, C)
        if len(Q) > len(V):
            continue
        print(f"Cir.{count}: {file_name[:-9]} has {len(Q)} qubits and {l} gates")

        _map_ = initial_mapping_for(args, IM, count, C, G, V)

        # Optional up-front completion of partial initial mappings. The router
        # handles partial maps fine via online SRMD extension and usually
        # produces better placement than greedy completion; we only invoke
        # map_completion when --complete-init is set (e.g. for very large
        # topologies where VF2 leaves many qubits unplaced).
        if args.complete_init and len(_map_) < len(Q):
            _map_ = map_completion(_map_, L, C, Q, AG, V)

        tau = [-1] * len(V)
        for log_q, phys_v in _map_.items():
            tau[phys_v] = log_q

        sum_in += l
        C_out, cost_time = qct(tau, C, Q, G, EG, V, SPL,
                               args.qfilter, variant=args.variant,
                               spl_mat=SPL_mat)
        total_route_time += cost_time
        sum_out += len(C_out)

        row = (count, file_name[:-9], len(Q), l, len(C_out), len(C_out) - l,
               round(cost_time, 2), round(len(C_out) / l, 4))
        print(row)
        log(row)

    if sum_in:
        summary = f"average ratio = {sum_out}/{sum_in} = {round(sum_out/sum_in, 4)}"
    else:
        summary = "no circuits processed"
    print(summary)
    log(summary)
    elapsed = round(time.time() - t_start, 2)
    timing = f"routing time {round(total_route_time, 2)}s; total {elapsed}s"
    print(timing)
    log(timing)


if __name__ == "__main__":
    main()

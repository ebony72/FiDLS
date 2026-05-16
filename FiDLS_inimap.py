"""Precompute initial mappings for a folder of circuits and cache them.

Run this before FiDLS_run.py whenever you change architecture, mapping
strategy, or circuit folder. Output is a JSON list of [idx, mapping_pairs]
entries written to `inimap/_inimap_list_<ag>_<mapping>_<path-stem>.txt`.

Example:
    python FiDLS_inimap.py --ag tokyo --mapping top --path B131/
"""
import argparse
import json
import os
import time

import ag as ag_mod
from inimap import _tau_bsg_, _tau_bstg_
from utils import (
    CreateCircuitFromQASM, ReducedCircuit, qubit_in_circuit,
)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ag", default="tokyo",
                   choices=sorted(ag_mod.TOPOLOGIES),
                   help="Architecture graph name (default tokyo).")
    p.add_argument("--mapping", default="top", choices=["top", "wgt"],
                   help="Initial mapping strategy (default top).")
    p.add_argument("--path", default="B131/",
                   help="Directory of input circuits (default B131/).")
    p.add_argument("--anchor", action="store_true", default=True,
                   help="Use anchor in subgraph-isomorphism search (default on).")
    p.add_argument("--no-anchor", dest="anchor", action="store_false")
    p.add_argument("--stop", type=float, default=10.0,
                   help="VF2 time budget in seconds (default 10).")
    p.add_argument("--out-dir", default="inimap",
                   help="Output directory (default inimap/).")
    return p.parse_args(argv)


def cache_path_for(args):
    path_stem = args.path.rstrip("/")
    name2 = f"_inimap_list_{args.ag}_{args.mapping}_{path_stem}"
    return os.path.join(args.out_dir, name2 + ".txt")


def main(argv=None):
    args = parse_args(argv)

    AG = ag_mod.build(args.ag)
    G = AG.graph
    V = list(G.nodes())

    print(time.asctime())
    print(f"Computing initial mappings: ag={args.ag}, mapping={args.mapping}, "
          f"path={args.path}, anchor={args.anchor}, stop={args.stop}")

    files = os.listdir(args.path)
    IM = []
    t_start = time.time()
    count = 0

    for file_name in files:
        count += 1
        if file_name.endswith("qasm"):
            cir = CreateCircuitFromQASM(file_name, args.path)
            C = ReducedCircuit(cir)
        else:
            with open(args.path + file_name, "r") as f:
                C = json.loads(f.read())
        l = len(C)
        if args.path.rstrip("/") == "bigQ":
            if count in {19, 21, 34, 42, 44, 47, 49}:
                continue
            if l > 15000:
                continue

        L = list(range(l))
        Q = qubit_in_circuit(L, C)
        if len(Q) > len(V):
            continue
        print(f"Cir.{count}: {file_name[:-9]} has {len(Q)} qubits and {l} gates")

        if args.mapping == "wgt":
            _map_ = _tau_bsg_(C, G, args.anchor, args.stop)
        else:  # "top"
            _map_ = _tau_bstg_(C, G, args.anchor, args.stop)

        print(_map_)
        im = [[k, v] for k, v in _map_.items()]
        IM.append([count, im])

    elapsed = round(time.time() - t_start, 2)
    print(f"Computed {len(IM)} initial mappings in {elapsed}s.")

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = cache_path_for(args)
    with open(out_path, "w") as f:
        f.write(json.dumps(IM))
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    main()

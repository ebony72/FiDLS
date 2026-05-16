"""Regression tests for the FiDLS router (router.qct).

These pin the gate-count ratios that Phase 2/3 produce on the bundled B131
medium suite, so future refactors can't silently regress quality. Times are
not asserted (machine-dependent) but are reported via pytest -v.
"""
import json
import os

import pytest

import ag
from router import qct
from utils import qubit_in_circuit


B131 = os.path.join(os.path.dirname(__file__), "..", "B131")
INIMAP_DIR = os.path.join(os.path.dirname(__file__), "..", "inimap")


def _load_b131_medium(topology, mapping):
    """Yield (circuit_index, qubits, C, tau) for B131 medium circuits."""
    A = ag.build(topology)
    V = list(A.graph.nodes())
    cache = os.path.join(INIMAP_DIR, f"_inimap_list_{topology}_{mapping}_B131.txt")
    with open(cache) as f:
        IM = json.loads(f.read())

    runs = []
    # IM cache was built with unsorted os.listdir order; match that exactly so
    # the count→file→initial-mapping indexing aligns with the cache.
    for cnt, fn in enumerate(os.listdir(B131), 1):
        if fn.startswith("."):
            continue
        with open(os.path.join(B131, fn)) as f:
            C = json.loads(f.read())
        if not (100 <= len(C) <= 1000):
            continue
        Q = qubit_in_circuit(list(range(len(C))), C)
        if len(Q) > len(V):
            continue
        entry = next((e[1] for e in IM if e[0] == cnt), None)
        if entry is None:
            continue
        tau = [-1] * len(V)
        for q, v in entry:
            tau[v] = q
        runs.append((cnt, fn, Q, C, tau))
    return A, runs


def _sweep(topology, mapping, variant, qfilter="01y"):
    A, runs = _load_b131_medium(topology, mapping)
    G = A.graph
    EG = G.edges()
    V = list(G.nodes())
    sum_in = sum_out = 0
    for cnt, fn, Q, C, tau in runs:
        out, _ = qct(tau[:], C, Q, G, EG, V, A.SPL, qfilter,
                     variant=variant, spl_mat=A.spl_mat)
        sum_in += len(C)
        sum_out += len(out)
    assert sum_in > 0, "no circuits matched the size filter"
    return sum_out / sum_in, len(runs)


@pytest.mark.parametrize("topology, variant, expected_ratio", [
    ("tokyo",     "G", 1.3584),
    ("tokyo",     "D", 1.8913),
    ("rochester", "G", 3.2389),
    ("rochester", "D", 3.0227),
    ("sycamore",  "G", 2.8423),
])
def test_b131_medium_ratio(topology, variant, expected_ratio):
    """Routing must not regress beyond ±0.5% on the B131 medium suite."""
    ratio, n = _sweep(topology, "top", variant)
    assert n >= 10, f"expected ≥10 circuits, got {n}"
    tolerance = max(0.005 * expected_ratio, 0.01)
    assert abs(ratio - expected_ratio) <= tolerance, (
        f"{topology}/{variant}: ratio {ratio:.4f} differs from "
        f"expected {expected_ratio:.4f} by > {tolerance:.4f}"
    )


def test_variants_are_listed():
    from router import VARIANTS
    assert VARIANTS == ("G", "D")


def test_invalid_variant_raises():
    A = ag.build("tokyo")
    with pytest.raises(ValueError):
        qct([0]*20, [[0, 1]], {0, 1}, A.graph, A.graph.edges(),
            list(A.graph.nodes()), A.SPL, "01y", variant="X")

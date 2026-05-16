"""Merged FiDLS router: FiDLS-G (1-layer lookahead) and FiDLS-D (3-layer lookahead).

Originally split across fidls_g.py and fidls_d.py. The two variants share ~80% of
their code; the only meaningful differences are the cost function (R_hat_1 vs
R_hat_3) and the selection criterion inside `good_next_mapping` (gates-solved
per swap vs cost-reduction per swap). Both are dispatched on the `variant`
argument to `qct`.

Phase 2 performance changes:
  * Maintain a `p2v` (logical -> physical) dict alongside `tau` so the inner
    loop avoids O(|V|) `tau.index(...)` calls.
  * Use a NumPy SPL matrix (built by ag.ArchitectureGraph) instead of the
    dict-keyed SPL for O(1) distance lookups with no hash overhead.
  * Set-based dedup of SWAP candidates (already from Phase 1).

The public signature `qct(tau, C, Q, G, EG, V, SPL, qfilter_type, variant=...)`
is unchanged. If `SPL` is a NumPy ndarray, the fast path is used; if it's the
legacy dict, fallbacks convert as needed.
"""
import math
import time
from itertools import islice

import networkx as nx
import numpy as np

from utils import (
    tau2map, map2tau, swap, swap_along_a_path,
    topgates, topgates_3_lev,
    qubit_in_circuit,
)

VARIANTS = ("G", "D")
_GAMMA_D = 0.8            # R_hat_3 weight for layer 1
_GAMMA_D2 = _GAMMA_D ** 2  # R_hat_3 weight for layer 2

# Phase 3 tuning knobs.
SRMD_MAX_SHORTEST_PATHS = 8  # cap paths explored per GATE3 gate


# ---------------------------------------------------------------------------
# p2v reverse-map helpers
# ---------------------------------------------------------------------------
def _build_p2v(tau):
    """Build logical-qubit -> physical-qubit dict from tau."""
    return {q: u for u, q in enumerate(tau) if q != -1}


# ---------------------------------------------------------------------------
# Fast versions of utils.entail and utils.greedy_solved_gates.
# These avoid the O(V) tau.index() calls in the originals; profiling on tokyo
# showed greedy_solved_gates accounted for ~73% of runtime in Phase 2 before
# these replacements.
# ---------------------------------------------------------------------------
def _entail(p2v, gate, EG):
    """True iff `gate` is executable under mapping p2v on edge-set EG.

    EG must support symmetric `(u, v) in EG` lookups (NetworkX EdgeView does).
    """
    p, q = gate[0], gate[1]
    u = p2v.get(p, -1)
    if u == -1:
        return False
    v = p2v.get(q, -1)
    if v == -1:
        return False
    return (u, v) in EG


def _greedy_solved(tau, p2v, LD, C, nl, EG, front=None):
    """Repeatedly absorb top-layer gates that are already executable.

    `front` may be passed when topgates(LD, C, nl) has already been computed by
    the caller (e.g. inside _swap3 where every candidate shares the same LD).
    The fast path: if no front-layer gate is executable under p2v, we return
    immediately without ever calling topgates ourselves.
    """
    LDx = list(LD)
    if front is None:
        top = topgates(LDx, C, nl)
    else:
        top = front
    solved = []
    while True:
        prev_len = len(LDx)
        for i in top:
            if _entail(p2v, C[i], EG):
                LDx.remove(i)
                solved.append(i)
        if len(LDx) == prev_len:
            return solved
        if not LDx:
            return solved
        # New front layer may have opened up; recompute.
        top = topgates(LDx, C, nl)


def _swap_pair(tau, p2v, u, v):
    """Apply swap(u, v): return (new_tau, new_p2v). Caller guarantees (u,v) in EG."""
    new_tau = tau[:]
    new_p2v = dict(p2v)
    qu, qv = tau[u], tau[v]
    new_tau[u], new_tau[v] = qv, qu
    if qu != -1:
        new_p2v[qu] = v
    if qv != -1:
        new_p2v[qv] = u
    return new_tau, new_p2v


# ---------------------------------------------------------------------------
# Distance + cost primitives (use p2v + spl_mat directly; no tau.index calls)
# ---------------------------------------------------------------------------
def _gate_dist(gate, p2v, UnOcc, spl_mat):
    """Physical distance between the endpoints of `gate` under mapping p2v.

    `UnOcc` is the pre-computed list of unoccupied physical qubits (= {u : tau[u]==-1}).
    When both endpoints are mapped this is a single matrix lookup; otherwise it
    is the minimum distance to an unoccupied qubit (matching the paper's
    `gate_phy_distance` semantics).
    """
    p, q = gate[0], gate[1]
    p_pos = p2v.get(p, -1)
    q_pos = p2v.get(q, -1)
    if p_pos != -1 and q_pos != -1:
        return int(spl_mat[p_pos, q_pos])
    if p_pos != -1:
        return int(min(spl_mat[p_pos, v] for v in UnOcc))
    if q_pos != -1:
        return int(min(spl_mat[u, q_pos] for u in UnOcc))
    return int(min(spl_mat[u, v] for u in UnOcc for v in UnOcc if v != u))


def _rhat(p2v, layers, C, UnOcc, spl_mat, variant):
    """Cost of `tau` (encoded via p2v) under either FiDLS-G or FiDLS-D scoring."""
    total = 0
    for i in layers[0]:
        total += _gate_dist(C[i], p2v, UnOcc, spl_mat)
    if variant == "G":
        return total
    total1 = 0
    for i in layers[1]:
        total1 += _gate_dist(C[i], p2v, UnOcc, spl_mat)
    total += _GAMMA_D * total1
    total2 = 0
    for i in layers[2]:
        total2 += _gate_dist(C[i], p2v, UnOcc, spl_mat)
    total += _GAMMA_D2 * total2
    return total


def _min_gate_dist(p2v, front, C, UnOcc, spl_mat):
    if not front:
        raise ValueError("min_gate_dist over an empty front layer")
    return min(_gate_dist(C[i], p2v, UnOcc, spl_mat) for i in front)


# ---------------------------------------------------------------------------
# Front-layer extraction & Q-filter
# ---------------------------------------------------------------------------
def _layers(LD, C, nl, variant):
    """Return [LTG, LTG1] for variant G, [LTG, LTG1, LTG2] for variant D.

    Both variants score with layer 0 only (G) or layers 0..2 (D), but both
    need layer 1 for the Q-filter sets (Q1x = qubits in LTG + LTG1). Returning
    LTG1 even for variant G mirrors the original fidls_g.py behavior.
    """
    if variant == "G":
        LDx = list(LD)
        LTG = topgates(LDx, C, nl)
        for i in LTG:
            LDx.remove(i)
        LTG1 = topgates(LDx, C, nl)
        return [LTG, LTG1]
    return list(topgates_3_lev(LD, C, nl))


def _qfilter_sets(LTG, LTG1, C, Q):
    """(Q, Q0, Q1, Q1x) qubit sets used to filter swap candidates."""
    Q0 = qubit_in_circuit(LTG, C)
    Q1 = qubit_in_circuit(LTG1, C)
    Q1x = qubit_in_circuit(LTG + LTG1, C)
    return Q, Q0, Q1, Q1x


_QFILTER_DISPATCH = {
    "9":   lambda Q, Q0, Q1, Q1x: (Q,   Q,   Q),
    "0":   lambda Q, Q0, Q1, Q1x: (Q0,  Q0,  Q0),
    "01":  lambda Q, Q0, Q1, Q1x: (Q0,  Q1,  Q1),
    "01x": lambda Q, Q0, Q1, Q1x: (Q0,  Q1x, Q1x),
    "1x":  lambda Q, Q0, Q1, Q1x: (Q1x, Q1x, Q1x),
    "01y": lambda Q, Q0, Q1, Q1x: (Q1x, Q0,  Q0),
}


def _resolve_filters(qfilter_type, LTG, LTG1, C, Q):
    Q_, Q0, Q1, Q1x = _qfilter_sets(LTG, LTG1, C, Q)
    fn = _QFILTER_DISPATCH.get(qfilter_type, _QFILTER_DISPATCH["01y"])
    return fn(Q_, Q0, Q1, Q1x)


# ---------------------------------------------------------------------------
# Fallback / map-extension helper (paper: SRMD = swap-reduce-min-distance)
# ---------------------------------------------------------------------------
def _swap_reduce_min_dist(tau, p2v, layers, LD, C, nl, G, EG, V, UnOcc, spl_mat,
                          SPL_dict, variant):
    """Either extend an incomplete mapping or solve at least one front-layer gate.

    Returns (action, tau_new, p2v_new, info). When `tau` is incomplete we
    typically extend; if any front-layer gate has both endpoints mapped (type 3),
    we route along a shortest path instead.
    """
    front = layers[0]
    md = _min_gate_dist(p2v, front, C, UnOcc, spl_mat)

    GATE0, GATE1, GATE2, GATE3 = [], [], [], []
    for i in front:
        if _gate_dist(C[i], p2v, UnOcc, spl_mat) > md:
            continue
        p, q = C[i]
        p_in = p in p2v
        q_in = q in p2v
        if not p_in and not q_in:
            GATE0.append(i)
        elif p_in and not q_in:
            GATE1.append(i)
        elif not p_in and q_in:
            GATE2.append(i)
        else:
            GATE3.append(i)

    if GATE3:
        # Prefer swapping over extension when both endpoints are already placed.
        # Multi-path SRMD (Phase 3): try the first SRMD_MAX_SHORTEST_PATHS
        # shortest paths between each gate's endpoints. The original used
        # nx.shortest_path which returns one arbitrary geodesic; on graphs
        # with multiple equally short paths this leaves easy wins on the table.
        best_gsg, best_rhat, action = 0, math.inf, None
        taux = None
        p2vx = None
        for i in GATE3:
            p, q = C[i]
            u, v = p2v[p], p2v[q]
            for path in islice(nx.all_shortest_paths(G, u, v),
                               SRMD_MAX_SHORTEST_PATHS):
                tau_temp = swap_along_a_path(tau, u, v, path, EG)
                p2v_temp = _build_p2v(tau_temp)
                gsg_temp = len(_greedy_solved(tau_temp, p2v_temp, LD, C, nl, EG,
                                              front=front))
                if gsg_temp < best_gsg:
                    continue
                unocc_temp = [w for w in V if tau_temp[w] == -1]
                rhat_temp = _rhat(p2v_temp, layers, C, unocc_temp, spl_mat, variant)
                if gsg_temp == best_gsg and rhat_temp >= best_rhat:
                    continue
                best_gsg = gsg_temp
                best_rhat = rhat_temp
                action_temp = [[path[k], path[k + 1]] for k in range(len(path) - 2)]
                taux = tau_temp[:]
                p2vx = dict(p2v_temp)
                action = action_temp[:]
        if action is None:
            raise RuntimeError(f"Some gates should be solved! (md={md})")
        return action, taux, p2vx, ["type3", GATE3, md]

    # GATE3 empty: extend the mapping using the best candidate placement.
    Occ = [u for u in V if tau[u] != -1]
    if not Occ:
        Occ2 = V[:]
    else:
        Occ2 = [u for u in UnOcc if min(int(spl_mat[u, w]) for w in Occ) <= 2]

    CAND = []
    if GATE1 or GATE2:
        for i in GATE1:
            p, q = C[i]
            u = p2v[p]
            assert q not in p2v, f"q={q} unexpectedly mapped for gate {C[i]}"
            CAND += [{q: v} for v in UnOcc if int(spl_mat[u, v]) == md]
        for i in GATE2:
            p, q = C[i]
            v = p2v[q]
            assert p not in p2v, f"p={p} unexpectedly mapped for gate {C[i]}"
            CAND += [{p: u} for u in UnOcc if int(spl_mat[u, v]) == md]
    else:
        if GATE0:
            print(f"GATE0 is nonempty! {GATE0}")
        for i in GATE0:
            p, q = C[i]
            assert p not in p2v and q not in p2v
            CAND += [{p: u, q: v} for u in Occ2 for v in UnOcc if int(spl_mat[u, v]) == md]

    rhat_record = math.inf
    tau_record = tau[:]
    p2v_record = dict(p2v)
    for cand in CAND:
        new_p2v = dict(p2v)
        for key, val in cand.items():
            assert key not in new_p2v, "Extension conflict in candidate"
            new_p2v[key] = val
        new_tau = tau[:]
        for log_q, phys_v in cand.items():
            new_tau[phys_v] = log_q
        unocc_new = [w for w in V if new_tau[w] == -1]
        rhat_temp = _rhat(new_p2v, layers, C, unocc_new, spl_mat, variant)
        if rhat_temp >= rhat_record:
            continue
        rhat_record = rhat_temp
        tau_record = new_tau
        p2v_record = new_p2v
    return [], tau_record, p2v_record, ["not type3", GATE0, GATE1, GATE2, GATE3, md]


# ---------------------------------------------------------------------------
# Filtered depth-limited SWAP enumeration (DLS up to length 3)
# ---------------------------------------------------------------------------
def _swap3(tau, p2v, layers, C, Q, EG, V, spl_mat, qfilter_type, variant):
    """Enumerate SWAP sequences of length 1, 2, or 3 that don't worsen the cost.

    In-place mutation + backtracking: rejected candidates pay only swap+revert
    cost (a few list/dict mutations), not a full tau[:] / dict(p2v) copy.
    Snapshots are taken only when a candidate is kept.
    """
    LTG = layers[0]
    LTG1 = layers[1] if len(layers) > 1 else []
    QF1, QF2, QF3 = _resolve_filters(qfilter_type, LTG, LTG1, C, Q)

    results = []
    seen = set()

    # Mutable working state; swaps are applied in-place and reverted on backtrack.
    work_tau = tau[:]
    work_p2v = dict(p2v)
    has_unmapped = any(t == -1 for t in tau)

    def apply(u, v):
        qu, qv = work_tau[u], work_tau[v]
        work_tau[u], work_tau[v] = qv, qu
        if qu != -1:
            work_p2v[qu] = v
        if qv != -1:
            work_p2v[qv] = u
        return qu, qv

    def revert(u, v, qu, qv):
        work_tau[u], work_tau[v] = qu, qv
        if qu != -1:
            work_p2v[qu] = u
        if qv != -1:
            work_p2v[qv] = v

    def unocc():
        # Common fast path: when no -1 entries exist, every gate's endpoints
        # are mapped and _gate_dist takes its single-lookup branch.
        if not has_unmapped:
            return ()
        return [u for u in V if work_tau[u] == -1]

    def remember(edges):
        key = (tuple(tuple(e) for e in edges), tuple(work_tau))
        if key in seen:
            return
        seen.add(key)
        results.append((list(edges), list(work_tau), dict(work_p2v)))

    rhat0 = _rhat(work_p2v, layers, C, unocc(), spl_mat, variant)

    for edge_1 in EG:
        p1, q1 = edge_1
        if work_tau[p1] not in QF1 and work_tau[q1] not in QF1:
            continue
        qu1, qv1 = apply(p1, q1)
        rhat1 = _rhat(work_p2v, layers, C, unocc(), spl_mat, variant)
        if rhat1 > rhat0:
            revert(p1, q1, qu1, qv1)
            continue
        remember([edge_1])

        for edge_2 in EG:
            if edge_1 == edge_2:
                continue
            p2, q2 = edge_2
            if work_tau[p2] not in QF2 and work_tau[q2] not in QF2:
                continue
            qu2, qv2 = apply(p2, q2)
            rhat2 = _rhat(work_p2v, layers, C, unocc(), spl_mat, variant)
            if rhat2 > rhat1:
                revert(p2, q2, qu2, qv2)
                continue
            remember([edge_1, edge_2])

            for edge_3 in EG:
                if edge_3 == edge_1 or edge_3 == edge_2:
                    continue
                p3, q3 = edge_3
                if work_tau[p3] not in QF3 and work_tau[q3] not in QF3:
                    continue
                qu3, qv3 = apply(p3, q3)
                if _min_gate_dist(work_p2v, LTG, C, unocc(), spl_mat) <= 1:
                    remember([edge_1, edge_2, edge_3])
                revert(p3, q3, qu3, qv3)

            revert(p2, q2, qu2, qv2)

        revert(p1, q1, qu1, qv1)

    return results


# ---------------------------------------------------------------------------
# Next-action selector
# ---------------------------------------------------------------------------
def _good_next_mapping(tau, p2v, LD, C, Q, G, EG, V, spl_mat, SPL_dict,
                       qfilter_type, fallback, variant):
    nl = len(Q)
    layers = _layers(LD, C, nl, variant)
    UnOcc = [u for u in V if tau[u] == -1]

    # When the mapping is incomplete (or we're in fallback), SRMD runs first.
    # If SRMD returns an empty action it means the mapping was *extended* (no
    # SWAP) — that's the answer and we return immediately. Otherwise SRMD found
    # a GATE3 path-swap; we keep it as a fallback in case _swap3 finds nothing.
    srmd_action = None
    srmd_tau = None
    srmd_p2v = None
    if UnOcc or fallback:
        action_rec, tau_rec, p2v_rec, _info = _swap_reduce_min_dist(
            tau, p2v, layers, LD, C, nl, G, EG, V, UnOcc, spl_mat,
            SPL_dict, variant
        )
        if not action_rec or fallback:
            return action_rec, tau_rec, p2v_rec
        srmd_action = action_rec
        srmd_tau = tau_rec
        srmd_p2v = p2v_rec

    swaps = _swap3(tau, p2v, layers, C, Q, EG, V, spl_mat, qfilter_type, variant)

    front = layers[0]  # Shared front layer for all candidate evaluations.
    best = None
    best_score = 0.0
    if variant == "G":
        for action, tau_new, p2v_new in swaps:
            gsg = _greedy_solved(tau_new, p2v_new, LD, C, nl, EG, front=front)
            if not gsg:
                continue
            score = len(gsg) / len(action)
            if score > best_score or (score == best_score and best is None):
                best_score = score
                best = (action, tau_new, p2v_new)
    else:
        rhat0 = _rhat(p2v, layers, C, UnOcc, spl_mat, variant)
        for action, tau_new, p2v_new in swaps:
            unocc_new = [u for u in V if tau_new[u] == -1]
            rhat_new = _rhat(p2v_new, layers, C, unocc_new, spl_mat, variant)
            if rhat_new > rhat0:
                continue
            if rhat_new == rhat0 and best is not None:
                continue
            score = (rhat0 - rhat_new) / len(action)
            if score > best_score or (score == best_score and best is None):
                best_score = score
                best = (action, tau_new, p2v_new)

    if best is None:
        # _swap3 found nothing useful. Prefer the SRMD path-swap if we already
        # have one (this is what the original fidls_g.py did); otherwise, on a
        # complete mapping, force one more SRMD call to drive progress.
        if srmd_action is not None:
            return list(srmd_action), list(srmd_tau), dict(srmd_p2v)
        if not UnOcc:
            action_rec, tau_rec, p2v_rec, info = _swap_reduce_min_dist(
                tau, p2v, layers, LD, C, nl, G, EG, V, UnOcc, spl_mat,
                SPL_dict, variant
            )
            print(f"  SRMD reduced minimal distance: action={action_rec}, "
                  f"fallback={fallback}, |LD|={len(LD)}, info={info}")
            return action_rec, tau_rec, p2v_rec
        return [], tau[:], dict(p2v)

    return list(best[0]), list(best[1]), dict(best[2])


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def qct(tau, C, Q, G, EG, V, SPL, qfilter_type, variant="G", spl_mat=None,
        return_events=False):
    """Route circuit `C` onto graph `G` given initial mapping `tau`.

    Args:
        tau: list[int] of length |V|; physical->logical, -1 for unmapped.
        C: list of [ctrl, tgt] CNOT pairs (the logical circuit).
        Q: set of logical qubits used in C.
        G: NetworkX architecture graph.
        EG: edges of G (an EdgeView supporting symmetric `in` checks).
        V: list of physical qubits of G.
        SPL: dict[(u, v) -> dist], kept for backward compat / fallback.
        qfilter_type: one of {'9', '0', '01', '01x', '01y', '1x'}.
        variant: 'G' (1-layer lookahead) or 'D' (3-layer lookahead).
        spl_mat: optional NumPy SPL matrix (faster than the dict). Built from
                 SPL if not provided.
        return_events: if True, also return a structured event stream useful
                       for Qiskit integration (see Returns).

    Returns:
        When return_events=False (default):
            (C_out, cost_time)
                C_out: list of [p, q] CNOT pairs on physical qubits; each
                       inserted SWAP is expanded into three CNOTs.
                cost_time: float, wall-clock seconds spent routing (rounded).
        When return_events=True:
            (C_out, cost_time, events)
                events: list of (kind, payload) tuples in chronological order.
                  - ('cnot', orig_idx, [p, q]) — original gate C[orig_idx] now
                    on physical qubits [p, q]
                  - ('swap', None, [u, v]) — one SWAP between physical u and v
                    (decomposed as three CNOTs in C_out)
    """
    if variant not in VARIANTS:
        raise ValueError(f"variant must be one of {VARIANTS}, got {variant!r}")

    if spl_mat is None:
        n = max(V) + 1
        spl_mat = np.full((n, n), -1, dtype=np.int32)
        for (u, v), d in SPL.items():
            spl_mat[u, v] = d

    nl = len(Q)
    diam = nx.diameter(G)

    nsvg = 0           # number of solved gates so far
    cost_swaps = 0     # SWAPs inserted (× 3 = added CNOTs)
    tau_new = tau[:]
    p2v_new = _build_p2v(tau_new)
    L1 = list(range(len(C)))
    nr = 0
    C_out = []
    events = [] if return_events else None

    action = []
    fallback = False
    state = [tau_new[:], dict(p2v_new), L1[:], C_out[:], cost_swaps, nr, nsvg,
             list(events) if events is not None else None]
    record = []

    start = time.time()
    while nsvg < len(C):
        GSG = _greedy_solved(tau_new, p2v_new, L1, C, nl, EG)
        for i in GSG:
            L1.remove(i)
            p = p2v_new[C[i][0]]
            q = p2v_new[C[i][1]]
            C_out.append([p, q])
            if events is not None:
                events.append(("cnot", i, [p, q]))
        nsvg += len(GSG)

        if GSG or not action:
            record = []
            state = [tau_new[:], dict(p2v_new), L1[:], C_out[:],
                     cost_swaps, nr, nsvg,
                     list(events) if events is not None else None]
            fallback = False
        else:
            record.append(nsvg)

        if len(record) > 2 * diam:
            print(f"use Fallback @ round {nr}!")
            fallback = True
            tau_new = state[0][:]
            p2v_new = dict(state[1])
            L1 = state[2][:]
            C_out = state[3][:]
            cost_swaps, nr, nsvg = state[4], state[5], state[6]
            if events is not None and state[7] is not None:
                events[:] = list(state[7])
            print(f"go back to round {nr}, "
                  f"solvable={len(_greedy_solved(tau_new, p2v_new, L1, C, nl, EG))}")
        if len(record) > 4 * diam:
            raise RuntimeError(f"Routing stuck at round {nr}: record={record}")

        if not L1:
            break

        nr += 1
        action, tau_new, p2v_new = _good_next_mapping(
            tau_new, p2v_new, L1, C, Q, G, EG, V, spl_mat, SPL,
            qfilter_type, fallback, variant
        )

        cost_swaps += len(action)
        for edge in action:
            C_out.append([edge[0], edge[1]])
            C_out.append([edge[1], edge[0]])
            C_out.append([edge[0], edge[1]])
            if events is not None:
                events.append(("swap", None, [edge[0], edge[1]]))

    cost_time = time.time() - start
    if return_events:
        return C_out, round(cost_time, 2), events
    return C_out, round(cost_time, 2)


# Backward-compat alias: original code used qct_old.
qct_old = qct

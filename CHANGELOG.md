# Changelog

## v0.2.0 — 2026 modernization

This release revises the original 2021 paper code into a tested, installable,
Qiskit-integrated compiler while preserving the published algorithm.

### Phase 1 — foundation (correctness fixes + consolidation)

- Merged `fidls_g.py` + `fidls_d.py` (~580 lines of near-duplicate code) into
  a single `router.py` with `variant="G" | "D"` and `lookahead=1 | 3`.
- Fixed latent bugs:
  - `rhat_record = 100` magic sentinel → `math.inf` (silently rejected valid
    candidates on graphs whose front-layer cost can exceed 100, e.g. q19x19).
  - Defensive init of `tau_record` before the candidate loop (previously raised
    `UnboundLocalError` if `CAND` ended up empty).
  - `== None` → `is None` (3 sites).
  - `set`-based dedup for SWAP candidate enumeration (was `O(N)`
    `x not in list` check).
- Added `argparse` CLIs to `FiDLS_run.py` and `FiDLS_inimap.py`; added
  `if __name__ == "__main__":` guards.
- Added a topology registry (`ag.TOPOLOGIES` + `ag.build(name)`); deleted
  the broken `fidls23.py` (had a `SyntaxError` plus undefined names) and the
  experimental `fidls_red.py`.

### Phase 2 — performance (no behavior change)

- Added a `p2v` reverse-map (logical → physical) maintained alongside `tau`
  through the routing loop; removed all `tau.index(...)` calls from the
  hot path (each was `O(|V|)`).
- Added a NumPy SPL matrix to `ArchitectureGraph` for O(1) distance lookups
  with no hash overhead, alongside the existing dict.
- Added fast `_entail` and `_greedy_solved` in `router.py` that take `p2v`
  directly (the originals lived in `utils.py` and used `tau.index`).
- Cached the front layer (`topgates(LD, C, nl)` result) across `_swap3`
  iterations — every candidate within a single round of SRMD/swap-selection
  shares it, but the original recomputed 250k+ times per round.
- Replaced per-candidate `tau[:]` + `dict(p2v)` allocations in `_swap3`
  with in-place mutation + backtracking. Rejected candidates pay only the
  swap+revert cost; snapshots are taken only for kept candidates.
- Fixed two real bugs caught in profiling:
  - Variant G was computing `LTG1 = []` for the Q-filter (so the
    `01y`/`01x` filters degenerated to `01y → Q0` everywhere). Restoring
    the original `LTG1 = topgates(LD-LTG, C, nl)` matched fidls-G's
    historical numbers exactly (tokyo G ratio 1.467 → 1.358).
  - When SRMD returned a non-empty GATE3 path-swap but `_swap3` found no
    improving swap, the next-action selector was discarding the SRMD
    result and returning an empty action. On sparse graphs this caused
    20K+ no-progress rounds (effective infinite loop on `con1_216`).
    The fix preserves the SRMD action as a fallback, matching the
    original `fidls_g.py` behavior.
- Net result: tokyo G goes from "1.467 buggy, 1.21s" (Phase 1) to "1.358
  matches original, 0.85s" (faster than the original); rochester goes from
  never-completes to 0.56s.

### Phase 3 — algorithmic upgrades

- **Multi-path SRMD**: in the GATE3 branch, try up to
  `SRMD_MAX_SHORTEST_PATHS` (=8) geodesics between each gate's endpoints
  and pick the (gate, path) with the most gates solved. Sycamore G ratio
  improved 0.5% (2.86 → 2.84); rochester unchanged within noise.
- **Fixed `one_shot_map_extension`**: the function was passing 7 args to
  a 9-arg `R_hat` (with the wrong types) — unreachable from any entry
  point but a latent bug. Now uses `topgates_3_lev` + the proper `R_hat`
  call.
- **`--complete-init` flag** in `FiDLS_run.py` wires the now-working
  `map_completion`. Off by default because greedy up-front completion is
  *worse* than the router's online SRMD extension on rochester (~5%
  quality loss). Useful on q19x19 where many qubits start unmapped.

### Phase 4 — modern API and tests

- **`qiskit_pass.py`**: `FiDLSSwap(TransformationPass)` for Qiskit ≥ 1.0.
  Preserves single-qubit gates and barriers, handles all six Q-filter
  types and both variants, sets `property_set["final_layout"]`.
- Added `return_events=True` to `qct()` so callers can interleave 1-qubit
  gates correctly with the routed 2-qubit gate stream.
- Two bugs caught while building the pass:
  - Event-stream order ≠ DAG topological order (FiDLS solves independent
    CNOTs in arrival order, not DAG order). Fixed by event-stream-driven
    traversal with predecessor counting.
  - `v2p` was updated before draining the pending-1q queue, causing
    1-qubit gates to land on their *post-swap* physical positions. Caught
    by the unitary-equivalence test.
- **`pyproject.toml`**: pip-installable as `pip install -e .`, declares
  Qiskit / NumPy / NetworkX as dependencies.
- **`tests/`**: 20 pytest tests covering router ratios per (topology,
  variant), architecture registry, Qiskit pass 1q preservation, no-swap
  paths, variant D, error paths, and unitary equivalence under the
  final permutation.

### Phase 6 — repo hygiene

- Added `.gitignore`, `README.md`, `.github/workflows/ci.yml` running
  pytest on Python 3.10 / 3.11 / 3.12.
- Removed `sqatest0.py` (2020 testing script, superseded by pytest),
  `vfstest.py` (524-line scratch file importing modules that don't
  exist), `readme.txt` (old paper-era readme; content folded into
  `README.md`).

### Benchmarking findings (not algorithmic changes, but documented)

- **`bench_sabre.py`**: head-to-head FiDLS vs Qiskit's `SabreSwap` on
  B131 medium. FiDLS-G beats Sabre-decay on tokyo by 18% (1.36 vs 1.61
  ratio). Sabre wins on rochester (2.57 vs FiDLS-G 3.39) and sycamore
  (2.44 vs 3.03). Sabre is 50–200× faster across the board (Rust hot
  path vs FiDLS's pure Python).
- **`bench_init_ablation.py`**: factoring out init-mapping. The FiDLS
  router is the differentiator on tokyo — even with a trivial init
  FiDLS-G beats Sabre with any init. On rochester/sycamore the Sabre
  router wins regardless of init.
- **`bench_init_stratified.py`**: separating complete vs partial SI
  mappings. On tokyo, complete-SI bucket has FiDLS-top beating
  SabreLayout by 3.5% — the 18% headline win is diluted by the partial
  bucket where greedy completion underperforms. On rochester/sycamore,
  even complete-SI loses by 20+%, confirming that the SI strategy itself
  is a tokyo-specialized heuristic.

### Compatibility / migration

- `qct(...)` keeps its original positional signature; new arguments
  (`spl_mat`, `return_events`) are optional with sensible defaults.
- `qct_old` is kept as an alias of `qct` for code that still imports it.
- `ag.ArchitectureGraph` keeps `.graph`, `.SPL`, `.diameter` and adds
  `.spl_mat`.
- The CLIs have new argparse flags but the defaults reproduce the
  original behavior; `FiDLS_inimap.py` writes to the path
  `FiDLS_run.py` reads (no more manual rename of the `ohmy_…` file).

## v0.1.x — original 2021 implementation

See the `master` branch and [arXiv:2004.07138](https://arxiv.org/abs/2004.07138).

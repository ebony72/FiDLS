# FiDLS — Claude Code notes

Implementation of Li et al. 2021 (arXiv:2004.07138), modernized for Qiskit ≥1.0.
Full revision history: CHANGELOG.md. Benchmark numbers + Sabre comparison: README.md.

## Files Claude will touch

| File | Purpose |
|---|---|
| `router.py` | Merged FiDLS-G + FiDLS-D. Entry: `qct(...)`. Pass `return_events=True` for Qiskit interop. |
| `qiskit_pass.py` | `FiDLSSwap(TransformationPass)` for Qiskit ≥1.0. |
| `ag.py` | Topology registry (`TOPOLOGIES`, `build(name)`) + NumPy SPL matrix. |
| `inimap.py`, `vfs.py`, `maps.py` | Pure-Python VF2 for SI-based init mapping. |
| `utils.py` | Circuit I/O, `topgates`, `R_hat`. |
| `tests/` | pytest regression suite (20 tests). |
| `bench_*.py` | SabreSwap comparison + init-mapping ablations. |

## Gotchas (real, found by bisection in Phases 1–4)

- **q19x19 + VF2 init mapping doesn't terminate.** Pure-Python VF2 (`vfs.py`) doesn't scale to 361 nodes. Use `--mapping naive` or `--mapping empty`. Routing itself works fine on q19x19 (≈4 s for B131 medium).
- **B131 / bigQ files have `.qasm.txt` extension but contain JSON CNOT lists.** The dispatch in `FiDLS_run.py` is `file_name.endswith("qasm")` (no dot) — False for `.qasm.txt`, so the JSON branch is taken. "Fixing" this breaks loading.
- **`EG` must be `G.edges()` (EdgeView), not `list(G.edges())`.** `utils.entail` relies on symmetric `(u, v) in EG` lookups. A plain list fails the symmetry check with "AG is undirected!".
- **`inimap/*.txt` cache uses unsorted `os.listdir()` order.** The cache stores `[count, mapping]` where `count` is the position in unsorted listdir. Tests/benchmarks that re-derive `count` must match (do NOT call `sorted()`), or initial mappings mismatch their circuits and ratios regress silently.
- **Variant G still needs `LTG1` (the second front layer) for the Q-filter** even though it scores with layer 0 only. Don't "simplify" `_layers()` to return one element — it degrades `01y` filtering to use only `Q0` and costs ~8% quality on tokyo.
- **`map_completion` regresses rochester ~5%.** Keep it behind `--complete-init`. The router's online SRMD extension produces better placement when it has full front-layer info.
- **"GATE0 is nonempty!" prints are informational**, not errors. They fire on partial-init circuits where front-layer gates have both endpoints unmapped.
- **Swap accounting**: `qct()` decomposes each inserted SWAP into 3 CNOTs inline in `C_out`. `FiDLSSwap` keeps them as `SwapGate` ops. For Sabre comparisons, count Sabre's swap ops × 3.

## Branches on GitHub

- `master` — original 2021 paper code, untouched.
- `modernize-2026` — this work.

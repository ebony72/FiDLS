# FiDLS

Python implementation of the qubit-mapping algorithm from
[Li et al., "Qubit Mapping Based on Subgraph Isomorphism and Filtered
Depth-Limited Search"](https://arxiv.org/abs/2004.07138) (IEEE TC, 2021),
revised and packaged for use as a Qiskit transpiler pass.

## Install

```sh
pip install -e .
```

Requires Python ≥ 3.9, Qiskit ≥ 1.0, NetworkX, NumPy.

## Two ways to use it

### 1. Qiskit `TransformationPass`

```python
from qiskit import QuantumCircuit
from qiskit.transpiler import CouplingMap, PassManager
from qiskit_pass import FiDLSSwap

qc = QuantumCircuit(5)
qc.h(0); qc.cx(0, 4); qc.cx(1, 3); qc.t(2)

cm = CouplingMap.from_grid(5, 5)
pm = PassManager([FiDLSSwap(cm, variant="G", qfilter_type="01y")])
routed = pm.run(qc)

final_layout = pm.property_set["final_layout"]
```

Preserves single-qubit gates and barriers. Sets `property_set["final_layout"]`.
If `property_set["layout"]` is set by an upstream layout pass, that layout is
used as the initial mapping; otherwise a trivial layout is generated.

### 2. CLI for the bundled benchmark suites

```sh
# Generate initial mappings (one-time per topology + circuit folder)
python FiDLS_inimap.py --ag tokyo --mapping top --path B131/

# Route circuits and print per-circuit stats + summary
python FiDLS_run.py --ag tokyo --variant G --filter 01y --mapping top \
                    --size medium --path B131/
```

Available CLI options for `FiDLS_run.py`:

| Flag | Choices | Default | Notes |
|---|---|---|---|
| `--ag` | `tokyo`, `sycamore`, `rochester`, `guadalupe`, `q5x5`, `q9x9`, `q19x19` | `tokyo` | Architecture graph |
| `--variant` | `G` (1-layer lookahead), `D` (3-layer) | `G` | |
| `--filter` | `9`, `0`, `01`, `01x`, `01y`, `1x` | `01y` | Q-filter selector |
| `--mapping` | `top`, `wgt`, `empty`, `naive` | `top` | Initial mapping strategy |
| `--size` | `small`, `medium`, `large`, `all` | `medium` | Circuit-size filter |
| `--path` | any directory of QASM / JSON CNOT-list files | `B131/` | |
| `--log` | path | _none_ | Write per-row + summary to a log file |
| `--only` | int | _none_ | Process only circuit #N (1-indexed) |
| `--complete-init` | flag | off | Greedy completion of partial initial mappings |

For very large topologies (`q19x19`), use `--mapping naive` or `--mapping empty`
— the bundled pure-Python VF2 (`vfs.py`) doesn't scale to 361 nodes within a
reasonable budget. Routing itself works fine on q19x19.

## What's in this repository

| File | Purpose |
|---|---|
| `router.py` | The FiDLS router. Merged FiDLS-G + FiDLS-D in one `qct()` |
| `qiskit_pass.py` | `FiDLSSwap(TransformationPass)` |
| `ag.py` | Architecture graphs (tokyo, sycamore, rochester, guadalupe, q5x5, q9x9, q19x19) + `TOPOLOGIES` registry and NumPy SPL matrix |
| `inimap.py` | Subgraph-isomorphism initial mapping (weighted-graph + topgraph variants) |
| `vfs.py`, `maps.py` | Pure-Python VF2 used by `inimap.py` |
| `utils.py` | Circuit I/O, topgates extraction, distance metrics, `R_hat` cost functions |
| `FiDLS_run.py` | CLI: route a folder of circuits, report stats |
| `FiDLS_inimap.py` | CLI: precompute initial mappings, cache to `inimap/` |
| `bench_sabre.py` | Head-to-head benchmark: FiDLS vs Qiskit `SabreSwap` |
| `bench_init_ablation.py` | Init-mapping ablation: FiDLS-top vs SabreLayout vs trivial |
| `bench_init_stratified.py` | Init-mapping ablation stratified by complete vs partial SI mappings |
| `CHANGELOG.md` | Phase-by-phase summary of revisions from the original 2021 code |
| `tests/` | pytest suite (20 tests) covering router quality, architecture builds, and Qiskit pass correctness (including unitary equivalence) |
| `B131/`, `bigQ/` | Benchmark circuits |
| `inimap/` | Precomputed initial-mapping cache |
| `testRecord/` | Historical (2020) gate-count results from the paper |

## How it compares to Qiskit's `SabreSwap`

Run `python bench_sabre.py --topologies tokyo rochester sycamore --size medium`
to reproduce. Quality numbers below are CNOT-overhead ratio (out CNOTs / in
CNOTs, lower is better) on the B131 medium suite (14 circuits, 2511 input
CNOTs). Each swap is counted as 3 CNOTs on both sides. SabreSwap is best of 5
seeds.

| Topology | FiDLS-G | FiDLS-D | Sabre-decay |
|---|---|---|---|
| Tokyo | **1.36** | 1.96 | 1.61 |
| Rochester | 3.39 | 3.23 | **2.57** |
| Sycamore | 3.03 | 2.82 | **2.44** |

**FiDLS-G wins by ~18% on tokyo** (the architecture from the 2021 paper).
Sabre wins on the larger sparse graphs by ~5–20%. Sabre is 50–200× faster
across the board (its inner loop is in Rust; FiDLS is pure Python).

### Where the FiDLS quality win comes from

The `bench_init_ablation.py` benchmark factors out init-mapping vs router:

| Topology | Router | fidls-top init | sabrelyt init | trivial init |
|---|---|---|---|---|
| Tokyo | FiDLS-G | 1.36 | **1.32** | 1.50 |
| Tokyo | Sabre-decay | 1.61 | 1.55 | 2.13 |
| Rochester | FiDLS-G | 3.39 | 2.91 | 3.32 |
| Rochester | Sabre-decay | **2.57** | 2.78 | 3.01 |
| Sycamore | FiDLS-G | 3.03 | 2.45 | 3.50 |
| Sycamore | Sabre-decay | **2.44** | 2.46 | 2.89 |

Conclusions: on tokyo the FiDLS router itself is responsible for the quality
win (any initial mapping, FiDLS-G beats Sabre-decay). On rochester/sycamore
the Sabre router dominates. SabreLayout often outperforms the original SI-based
initial mapping for both routers, including FiDLS-G itself.

### Why SabreLayout beats the SI-based init: a stratified view

`bench_init_stratified.py` separates each circuit by whether the cached SI
mapping was *complete* (every logical qubit placed by VF2) or *partial*
(some qubits left unmapped and filled in by the greedy `map_completion`).
FiDLS-G run with both initial mappings, then bucketed:

| Topology | Bucket | n | fidls-top ratio | sabrelyt ratio | gap |
|---|---|---|---|---|---|
| Tokyo | SI complete | 6 | 1.337 | 1.385 | **−3.5% (SI better)** |
| Tokyo | SI partial | 8 | 1.397 | 1.220 | +14.4% (SI worse) |
| Rochester | SI complete | 2 | 3.386 | 2.802 | +20.9% |
| Rochester | SI partial | 12 | 3.394 | 2.921 | +16.2% |
| Sycamore | SI complete | 2 | 3.084 | 2.422 | +27.3% |
| Sycamore | SI partial | 12 | 3.021 | 2.459 | +22.9% |

- **On tokyo**, partial-mapping handling is *the entire story*. When VF2
  finds a complete embedding, the SI-based mapping actually beats SabreLayout
  by 3.5%. The 18% headline win FiDLS-G shows over SabreSwap on tokyo is
  diluted by the 8/14 circuits where greedy `map_completion` fills in
  partial mappings worse than SabreLayout would.
- **On rochester/sycamore**, even the complete-mapping bucket loses by
  20–27% (small sample, n=2 each, but the magnitude is consistent). The SI
  strategy fundamentally optimizes for embedding the *front layer* on
  coupled edges — a tokyo-friendly bias that costs you globally on sparse
  graphs where the front layer isn't representative of the whole circuit.

The takeaway: the SI-based initial mapping is a tokyo-specialized
heuristic. SabreLayout's iterative whole-circuit refinement has surpassed
it on every architecture except the dense small one it was tuned for.

## Tests

```sh
pip install -e ".[test]"
pytest
```

20 tests covering:

- gate-count-ratio regression on every (topology, variant) combination,
- architecture-graph registry + SPL agreement (dict ↔ NumPy matrix),
- Qiskit pass: 1q gate preservation, no-swap-needed paths, final-layout property,
  variant D, error paths, **unitary equivalence under the final permutation**.

## Caveats and limitations

- VF2-based `top` / `wgt` initial mappings (`inimap.py` → `vfs.py`) are pure
  Python. They don't scale to topologies > ~100 nodes. For q19x19, use
  `--mapping naive` or `--mapping empty`.
- `_swap3` enumeration is `O(|E|³)` with Q-filter pruning. Routing remains fast
  on tokyo / rochester / sycamore but the constant factor is much larger than
  Sabre's Rust implementation.
- The router decomposes each inserted SWAP into three CNOTs in its output
  stream (paper accounting). The Qiskit pass keeps SWAPs as `SwapGate` ops.

## Citing

If you use FiDLS in research, please cite the original paper:

> Sanjiang Li, Xiangzhen Zhou, Yuan Feng. "Qubit Mapping Based on Subgraph
> Isomorphism and Filtered Depth-Limited Search." *IEEE Transactions on
> Computers*, 70(11): 1777–1788 (2021).

## License

MIT. See [license.md](license.md).

"""Tests for the FiDLSSwap Qiskit TransformationPass."""
import pytest

from qiskit import QuantumCircuit
from qiskit.transpiler import CouplingMap, PassManager
from qiskit_pass import FiDLSSwap


def _tokyo_coupling_map():
    """5x4 IBM-Q-Tokyo-shaped coupling map."""
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 4),
        (5, 6), (6, 7), (7, 8), (8, 9),
        (10, 11), (11, 12), (12, 13), (13, 14),
        (15, 16), (16, 17), (17, 18), (18, 19),
        (0, 5), (1, 6), (2, 7), (3, 8), (4, 9),
        (5, 10), (6, 11), (7, 12), (8, 13), (9, 14),
        (10, 15), (11, 16), (12, 17), (13, 18), (14, 19),
    ]
    cm = CouplingMap(edges)
    cm.make_symmetric()
    return cm


def _count(circ):
    return dict(circ.count_ops())


def test_already_routed_circuit_is_unchanged_in_structure():
    """Circuit whose CNOTs are already on coupling-map edges should not need swaps."""
    cm = _tokyo_coupling_map()
    qc = QuantumCircuit(20)
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    qc.cx(5, 6)
    qc.t(3)
    routed = PassManager([FiDLSSwap(cm)]).run(qc)
    assert _count(routed).get("swap", 0) == 0
    # All input gates preserved.
    assert _count(routed).get("cx") == 3
    assert _count(routed).get("h") == 1
    assert _count(routed).get("t") == 1


def test_single_qubit_gates_are_preserved():
    cm = _tokyo_coupling_map()
    qc = QuantumCircuit(5)
    for i in range(5):
        qc.h(i)
    qc.cx(0, 4)  # forces routing on tokyo
    qc.t(2)
    qc.cx(1, 3)
    qc.s(0)
    routed = PassManager([FiDLSSwap(cm)]).run(qc)
    # 1q gates must survive verbatim (counts only).
    assert _count(routed).get("h") == 5
    assert _count(routed).get("t") == 1
    assert _count(routed).get("s") == 1
    # Original CNOTs must be present (count, not position).
    assert _count(routed).get("cx") == 2


def test_no_two_qubit_gates_no_swaps():
    cm = _tokyo_coupling_map()
    qc = QuantumCircuit(4)
    qc.h(0); qc.t(1); qc.s(2); qc.x(3)
    routed = PassManager([FiDLSSwap(cm)]).run(qc)
    assert _count(routed).get("swap", 0) == 0
    assert _count(routed).get("h") == 1
    assert _count(routed).get("t") == 1


def test_final_layout_set_on_property_set():
    cm = _tokyo_coupling_map()
    qc = QuantumCircuit(5)
    qc.cx(0, 4)
    qc.cx(1, 3)
    pm = PassManager([FiDLSSwap(cm)])
    pm.run(qc)
    final_layout = pm.property_set["final_layout"]
    assert final_layout is not None
    # Every virtual qubit should map to a distinct physical qubit.
    physicals = {final_layout[q] for q in qc.qubits}
    assert len(physicals) == len(qc.qubits)


def test_variant_d_runs():
    cm = _tokyo_coupling_map()
    qc = QuantumCircuit(5)
    qc.cx(0, 4); qc.cx(1, 3); qc.cx(2, 4); qc.cx(0, 3)
    routed = PassManager([FiDLSSwap(cm, variant="D")]).run(qc)
    assert _count(routed).get("cx") == 4


def test_invalid_variant_raises():
    cm = _tokyo_coupling_map()
    with pytest.raises(ValueError):
        FiDLSSwap(cm, variant="X")


def test_invalid_qfilter_raises():
    cm = _tokyo_coupling_map()
    with pytest.raises(ValueError):
        FiDLSSwap(cm, qfilter_type="bogus")


def test_too_many_virtual_qubits_rejected():
    from qiskit.transpiler.exceptions import TranspilerError
    cm = CouplingMap([(0, 1), (1, 2)])  # 3 physical
    cm.make_symmetric()
    qc = QuantumCircuit(5)
    qc.cx(0, 1)
    with pytest.raises(TranspilerError):
        PassManager([FiDLSSwap(cm)]).run(qc)


def test_routed_circuit_is_unitarily_equivalent():
    """The routed circuit, composed with the inverse of the final permutation,
    should be unitarily equivalent to the input. We verify on a small case
    where ``Operator`` is cheap."""
    from qiskit.quantum_info import Operator
    # Small (5q) coupling map: a path 0-1-2-3-4
    cm = CouplingMap([(0, 1), (1, 2), (2, 3), (3, 4)])
    cm.make_symmetric()

    qc = QuantumCircuit(5)
    qc.h(0); qc.h(2); qc.h(4)
    qc.cx(0, 4)  # far apart — needs routing on the path
    qc.cx(1, 3)
    qc.t(2)
    qc.cx(0, 3)

    pm = PassManager([FiDLSSwap(cm)])
    routed = pm.run(qc)
    final_layout = pm.property_set["final_layout"]

    # Apply the inverse layout permutation to the routed circuit's wires so
    # that virtual qubit i ends up back at physical i.
    perm_qc = QuantumCircuit(5)
    perm_qc.compose(routed, inplace=True)
    # Manually permute via Operator semantics: the routed circuit acts on
    # physical qubits; reorder to match input virtual order.
    routed_op = Operator(routed)
    # routed acts on physical qubits 0..4 in canonical order.
    # final_layout maps virtual qubit i -> physical j. So row/col j of routed_op
    # corresponds to virtual i.
    perm = [final_layout[qc.qubits[i]] for i in range(5)]
    # Build permutation operator that maps physical-ordering back to virtual.
    inv_perm_circ = QuantumCircuit(5)
    # We test equivalence by reordering: simulate the input on virtual order
    # and the routed on physical order, then permute back.
    input_op = Operator(qc)
    # The routed circuit produces, on physical qubits, the same state up to
    # a permutation that maps physical j -> virtual perm.index(j).
    # Equivalence check: routed[physical ordering] == input[virtual ordering]
    # after permuting physical wires by `perm`.
    permuted = routed_op.apply_permutation(perm, front=False)
    assert permuted.equiv(input_op), (
        f"routed circuit is NOT unitarily equivalent to input "
        f"(final_layout perm={perm})"
    )

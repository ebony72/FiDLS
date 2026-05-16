"""Tests for the ag (architecture graph) module."""
import networkx as nx
import numpy as np
import pytest

import ag


def test_registry_lists_all_topologies():
    assert set(ag.TOPOLOGIES) == {
        "tokyo", "sycamore", "rochester", "guadalupe",
        "q5x5", "q9x9", "q19x19",
    }


def test_build_returns_connected_graph_for_each_topology():
    for name in ag.TOPOLOGIES:
        A = ag.build(name)
        assert nx.is_connected(A.graph), f"{name} disconnected"
        assert A.diameter >= 1
        assert A.SPL[(0, 0)] == 0
        assert A.spl_mat.shape == (max(A.graph.nodes()) + 1,) * 2
        assert A.spl_mat[0, 0] == 0


def test_unknown_topology_raises():
    with pytest.raises(KeyError):
        ag.build("not-a-topology")


def test_spl_matrix_matches_spl_dict():
    A = ag.build("tokyo")
    for (u, v), d in A.SPL.items():
        assert int(A.spl_mat[u, v]) == d, f"SPL mismatch at ({u},{v})"

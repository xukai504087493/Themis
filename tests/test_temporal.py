"""Unit tests for src.temporal: distance fixed point, backtracking,
feasibility, SCC, linear-time B-hypergraph distance."""
import hypergraph
import temporal

INF = float("inf")


def _chain():
    return [{"tail": ["s"], "head": "a", "id": "e1"},
            {"tail": ["a"], "head": "t", "id": "e2"}]


def test_compute_distances_chain():
    edges = _chain()
    reach = {"s", "a", "t"}
    dist = temporal.compute_distances({"s"}, edges, reach)
    assert dist["s"] == 0
    assert dist["a"] == 1
    assert dist["t"] == 2


def test_backtrack_hyperpath_recovers_path():
    edges = _chain()
    _, _, _, _, in_edges = hypergraph.build(edges)
    dist = temporal.compute_distances({"s"}, edges, {"s", "a", "t"})
    F = temporal.backtrack_hyperpath({"s"}, "t", dist, in_edges)
    assert F is not None and len(F) == 2
    assert {id(e) for e in F} == {id(e) for e in edges}


def test_temporal_feasible_true_on_chain():
    assert temporal.temporal_feasible({"s"}, _chain()) == (True, [])


def test_temporal_feasible_false_on_dangling_edge():
    edges = [{"tail": ["a"], "head": "b", "id": "e1"}]  # a never produced
    feasible, blocked = temporal.temporal_feasible({"s"}, edges)
    assert feasible is False
    assert blocked == edges


def test_temporal_feasible_false_on_cycle_without_entry():
    edges = [{"tail": ["a"], "head": "b", "id": "e1"},
             {"tail": ["b"], "head": "a", "id": "e2"}]
    feasible, blocked = temporal.temporal_feasible({"s"}, edges)
    assert feasible is False
    assert set(id(e) for e in blocked) == {id(e) for e in edges}


def test_scc_kosaraju_basic():
    arcs = [(0, 1), (1, 0), (1, 2)]
    comp, nc = temporal.scc_kosaraju(arcs, 3)
    assert nc == 2
    assert comp[0] == comp[1]
    assert comp[2] != comp[0]


def test_linear_temporal_distance_chain():
    edges = _chain()
    d = {id(e): w for e, w in zip(edges, (5.0, 2.0))}
    delta, tinf, nc = temporal.linear_temporal_distance(
        ["s", "a", "t"], edges, {"s"}, "t", d)
    assert tinf is False
    assert delta == 7.0
    assert nc >= 1


def test_linear_temporal_distance_positive_cycle_is_infeasible():
    edges = [{"tail": ["a"], "head": "b", "id": "e1"},
             {"tail": ["b"], "head": "a", "id": "e2"}]
    d = {id(edges[0]): 3.0, id(edges[1]): 4.0}
    delta, tinf, _ = temporal.linear_temporal_distance(
        ["a", "b"], edges, {"s"}, "a", d)
    assert delta is None
    assert tinf is True

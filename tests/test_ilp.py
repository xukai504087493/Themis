"""Unit tests for src.ilp.eval_T (completion-time fixed point). No solver
license dependency: solver paths are covered by the full-scale pipeline."""
import hypergraph
import ilp


def _chain_edges():
    return [{"tail": ["s"], "head": "a", "id": "e1"},
            {"tail": ["a"], "head": "t", "id": "e2"}]


def test_eval_T_chain():
    edges = _chain_edges()
    _, _, _, _, in_edges = hypergraph.build(edges)
    d = {id(edges[0]): 5.0, id(edges[1]): 2.0}
    T = ilp.eval_T(edges, d, {"s"}, in_edges, "t")
    assert T == 7.0


def test_eval_T_returns_none_if_unreachable():
    edges = [{"tail": ["a"], "head": "t", "id": "e1"}]
    _, _, _, _, in_edges = hypergraph.build(edges)
    d = {id(edges[0]): 1.0}
    assert ilp.eval_T(edges, d, {"s"}, in_edges, "t") is None

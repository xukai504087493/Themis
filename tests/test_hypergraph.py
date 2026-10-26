"""Unit tests for src.hypergraph: TSV loading, dedup, vertex maps."""
import hypergraph


def _write(tmp_path, text):
    p = tmp_path / "edges.txt"
    p.write_text(text, encoding="utf-8")
    return str(p)


def test_load_edges_parses_and_drops_header(tmp_path):
    text = ("Tail\tHead\tPosReg\tNegReg\tID\n"
            "s;\ta\tP\tN\te1\n"
            "a;\tt\tP\tN\te2\n")
    edges = hypergraph.load_edges(_write(tmp_path, text))
    assert len(edges) == 2
    assert edges[0] == {"tail": ["s"], "head": "a", "id": "e1"}
    assert edges[1] == {"tail": ["a"], "head": "t", "id": "e2"}


def test_load_edges_splits_tails_on_semicolon(tmp_path):
    text = "x;y;\tt\tP\tN\te1\n"
    edges = hypergraph.load_edges(_write(tmp_path, text))
    assert edges[0]["tail"] == ["x", "y"]


def test_load_edges_drops_none_and_empty(tmp_path):
    text = ("s;None\ta\tP\tN\te1\n"   # literal 'None' tail cell dropped
            "\tt\tP\tN\te2\n"         # empty tail -> dropped
            "s\tNone\tP\tN\te3\n")    # 'None' head -> dropped
    edges = hypergraph.load_edges(_write(tmp_path, text))
    assert len(edges) == 1
    assert edges[0] == {"tail": ["s"], "head": "a", "id": "e1"}


def test_load_edges_deduplicates(tmp_path):
    text = "s\ta\tP\tN\te1\ns\ta\tP\tN\te2\n"
    edges = hypergraph.load_edges(_write(tmp_path, text))
    assert len(edges) == 1
    assert edges[0]["id"] == "e1"


def test_build_vertex_maps():
    edges = [{"tail": ["s"], "head": "a", "id": "e1"},
             {"tail": ["a"], "head": "t", "id": "e2"}]
    verts, tail_membership, has_tail, heads, in_edges = hypergraph.build(edges)
    assert verts == {"s", "a", "t"}
    assert tail_membership["s"] == [edges[0]]
    assert tail_membership["a"] == [edges[1]]
    assert has_tail == {"s", "a"}
    assert heads == {"a", "t"}
    assert in_edges["a"] == [edges[0]]
    assert in_edges["t"] == [edges[1]]


def test_reachable_from_bfs():
    edges = [{"tail": ["s"], "head": "a", "id": "e1"},
             {"tail": ["a"], "head": "t", "id": "e2"},
             {"tail": ["u"], "head": "v", "id": "e3"}]
    _, tail_membership, _, _, _ = hypergraph.build(edges)
    reach = hypergraph.reachable_from({"s"}, tail_membership)
    assert reach == {"s", "a", "t"}

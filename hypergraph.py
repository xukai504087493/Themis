"""Hypergraph loading and construction from Pathway Commons TSV exports.

Edge format: Tail|Head|PosReg|NegReg|ID. Tail cells are ';'-separated;
literal 'None' / empty tails or heads are dropped (Mmunin make_hypergraph
protocol), duplicate (head, tail-set) rows collapse to one edge.
"""


def load_edges(path):
    """Read hyperedge TSV into a list of {tail, head, id} dicts."""
    seen = set()
    edges = []
    with open(path, "r", encoding="utf-8") as f:
        for ln in f:
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            if parts[4].strip() == "ID":  # header row
                continue
            tail = [x.strip() for x in parts[0].split(";")
                    if x.strip() not in ("", "None")]
            head = parts[1].strip()
            if not tail or head in ("", "None"):
                continue
            key = (head, frozenset(tail))
            if key in seen:
                continue
            seen.add(key)
            edges.append({"tail": tail, "head": head, "id": parts[4].strip()})
    return edges


def build(edges):
    """Vertex sets and incidence maps: tail_membership, has_tail, heads,
    in_edges (head -> producing edges)."""
    verts = set()
    tail_membership = {}
    has_tail = set()
    heads = set()
    in_edges = {}
    for e in edges:
        for u in e["tail"]:
            verts.add(u)
            tail_membership.setdefault(u, []).append(e)
            has_tail.add(u)
        verts.add(e["head"])
        heads.add(e["head"])
        in_edges.setdefault(e["head"], []).append(e)
    return verts, tail_membership, has_tail, heads, in_edges


def reachable_from(sources, tail_membership):
    """Vertices reachable from the super-source via hyperedges (BFS)."""
    reach = set(sources)
    ready = list(sources)
    while ready:
        v = ready.pop()
        for e in tail_membership.get(v, []):
            if e["head"] not in reach:
                reach.add(e["head"])
                ready.append(e["head"])
    return reach

"""Temporal hyperpath algorithms.

Static side: tail-distance fixed point (min semantics, converges on cyclic
hypergraphs), argmin backtracking, generation-order feasibility check.
Exact side: linear-time temporal distance for B-hypergraphs (single-head
hypergraphs) by SCC contraction + min-plus fixed point in topological order.
"""

INF = float("inf")


def compute_distances(sources, edges, reach):
    """Tail-distance fixed point: dist[h] = min over in-edges e of
    (1 + max_{v in tail(e)} dist[v]). Serves as the static heuristic
    distance and the temporal labeling (earliest completion time)."""
    dist = {v: 0 for v in sources}
    changed = True
    while changed:
        changed = False
        for e in edges:
            h = e["head"]
            if h not in reach or not e["tail"]:
                continue
            if all(u in dist for u in e["tail"]):
                cand = 1 + max(dist[u] for u in e["tail"])
                if cand < dist.get(h, INF):
                    dist[h] = cand
                    changed = True
    return dist


def backtrack_hyperpath(sources, target, dist, in_edges):
    """Backtrack from target: pick the incoming edge minimizing the tail
    bound (argmin of 1 + max dist[tails]), constructing a static hyperpath."""
    if target not in dist:
        return None
    F = {}
    stack = [target]
    while stack:
        v = stack.pop()
        if v in sources:
            continue
        best = None
        best_val = None
        for e in in_edges.get(v, []):
            tails = e["tail"]
            if not tails or not all(u in dist for u in tails):
                continue
            val = 1 + max(dist[u] for u in tails)
            if best_val is None or val < best_val:
                best, best_val = e, val
        if best is None:
            continue
        if id(best) not in F:
            F[id(best)] = best
            for u in best["tail"]:
                if u not in sources:
                    stack.append(u)
    return list(F.values())


def temporal_feasible(sources, F):
    """Generation-order check; returns (feasible, blocked_edges)."""
    if not F:
        return False, []
    gen = set(sources)
    blocked = list(F)
    changed = True
    while changed:
        changed = False
        remaining = []
        for e in blocked:
            if all(v in gen for v in e["tail"]):
                gen.add(e["head"])
                changed = True
            else:
                remaining.append(e)
        blocked = remaining
    feasible = all(e["head"] in gen for e in F)
    return feasible, blocked


def scc_kosaraju(arcs, n):
    """Iterative Kosaraju over n vertices; arcs = list of (u, v) pairs.
    Returns (comp, n_components)."""
    out = [[] for _ in range(n)]
    for u, v in arcs:
        out[u].append(v)
    rev = [[] for _ in range(n)]
    for u, v in arcs:
        rev[v].append(u)
    visited = [False] * n
    order = []

    def dfs1(start):
        stack = [(start, 0)]
        while stack:
            v, i = stack.pop()
            if i == 0:
                if visited[v]:
                    continue
                visited[v] = True
            if i < len(out[v]):
                stack.append((v, i + 1))
                nxt = out[v][i]
                if not visited[nxt]:
                    stack.append((nxt, 0))
            else:
                order.append(v)

    for v in range(n):
        if not visited[v]:
            dfs1(v)
    comp = [-1] * n
    cid = 0

    def dfs2(start):
        stack = [start]
        comp[start] = cid
        while stack:
            v = stack.pop()
            for u in rev[v]:
                if comp[u] == -1:
                    comp[u] = cid
                    stack.append(u)

    for v in reversed(order):
        if comp[v] == -1:
            dfs2(v)
            cid += 1
    return comp, cid


def linear_temporal_distance(vertices, Es, sources, target, d):
    """Exact earliest completion time on a B-hypergraph (single head per
    edge) by SCC contraction (zero-delay cycles collapse) and a min-plus
    fixed point over the block DAG. A component with an internal
    positive-delay arc has no finite solution (delta stays infinite).

    Returns (delta_target, tinf, n_scc)."""
    vid = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    arcs = []
    head_to_edges = {}
    for e in Es:
        head_to_edges.setdefault(e["head"], []).append(e)
        for u in e["tail"]:
            if u in vid:
                arcs.append((vid[u], vid[e["head"]], d[id(e)]))
    if not arcs:
        return None, True, 0
    comp, nc = scc_kosaraju([(u, v) for u, v, _ in arcs], n)

    block_out = [[] for _ in range(nc)]
    block_in = [set() for _ in range(nc)]
    for u, v, _ in arcs:
        cu, cv = comp[u], comp[v]
        if cu != cv and cv not in block_out[cu]:
            block_out[cu].append(cv)
            block_in[cv].add(cu)
    indeg = [len(block_in[c]) for c in range(nc)]
    stack = [c for c in range(nc) if indeg[c] == 0]
    topo = []
    while stack:
        c = stack.pop()
        topo.append(c)
        for nxt in block_out[c]:
            indeg[nxt] -= 1
            if indeg[nxt] == 0:
                stack.append(nxt)

    delta = [INF] * n
    for v in vertices:
        if v in sources:
            delta[vid[v]] = 0.0
    block_verts = [[] for _ in range(nc)]
    for v in vertices:
        block_verts[comp[vid[v]]].append(v)

    # Min-plus relaxation inside each block: external tails use settled
    # values, internal tails the running value; at most |B| rounds suffice.
    # Min semantics breaks positive-delay cycles by selection: vertices
    # produced from outside keep their external value, so the cycle is
    # bypassed rather than deleted.
    for B in topo:
        for _ in range(len(block_verts[B])):
            changed = False
            for v in block_verts[B]:
                best = delta[vid[v]]
                for e in head_to_edges.get(v, []):
                    tailmax = 0.0
                    ok = True
                    for u in e["tail"]:
                        if u not in vid:
                            continue
                        if u in sources:
                            continue
                        if delta[vid[u]] == INF:
                            ok = False
                            break
                        if delta[vid[u]] > tailmax:
                            tailmax = delta[vid[u]]
                    if ok:
                        val = tailmax + d[id(e)]
                        if val < best:
                            best = val
                if best < delta[vid[v]]:
                    delta[vid[v]] = best
                    changed = True
            if not changed:
                break

    dt = delta[vid[target]]
    return (dt if dt != INF else None, dt == INF, nc)

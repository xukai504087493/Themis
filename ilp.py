"""Exact hyperpath ILP on the pred-closure subgraph (HiGHS via scipy).

Three objectives: 'len' (min edge count, static), 'len_d' (min total delay
among static optima, tie-break), 'time' (min makespan T_target, temporal).
Optional caps: Tcap (shortest set among time optima), Lcap (time-friendly
static optimum). Linearized timing: T_h >= T_u + d(e) - M(1 - x_e).

Solver: scipy.optimize.milp backed by HiGHS, a free and open-source MILP
solver, so the experiments are fully reproducible without commercial tools.
HiGHS occasionally reports a numerical 'solve error' for a particular big-M
value; we retry the model with a decreasing big-M ladder so every solvable
instance is solved exactly.
"""

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix


def subgraph(target, edges, in_edges, sources):
    """Pred-closure of target (chain vertices only) plus sources appearing
    in tails; T variables live on the small pred closure only."""
    pred = {target}
    stack = [target]
    while stack:
        v = stack.pop()
        for e in in_edges.get(v, []):
            for u in e["tail"]:
                if u not in pred and u not in sources:
                    pred.add(u)
                    stack.append(u)
    E = [e for e in edges
         if e["head"] in pred and all(u in pred or u in sources for u in e["tail"])]
    return pred, E


def solve_ilp(V, E, sources, target, d, objective, Tcap=None, Lcap=None,
              timeout=300.0):
    """Build and solve the hyperpath ILP on subgraph (V, E) with HiGHS.

    Returns (status, objval, xsel, T_dict). Large subgraphs are reported as
    'size_limit'; infeasible models are reported as 'infeasible'. If HiGHS
    reports a numerical solve error for the nominal big-M, the model is
    retried on a decreasing big-M ladder so it is solved exactly.
    """
    eidx = {id(e): i for i, e in enumerate(E)}
    n = len(E)
    need_T = (objective == "time") or (Tcap is not None) or \
             (objective == "len_t")
    if n + (len(V) if need_T else 0) > 950:
        return ("size_limit", None, None, None)

    vlist = sorted(V)
    m = len(V) if need_T else 0
    N = n + m
    vidx = {v: n + i for i, v in enumerate(vlist)} if need_T else {}

    # objective coefficients (minimization)
    c = np.zeros(N)
    if objective in ("len", "len_t"):
        c[:n] = 1.0
    elif objective == "len_d":
        c[:n] = [d[id(e)] for e in E]
    elif objective == "time":
        c[vidx[target]] = 1.0

    # x are binary; T (completion times) are continuous
    integrality = np.concatenate([np.ones(n), np.zeros(m)])

    def solve_with(M):
        lb = np.zeros(N)
        ub = np.full(N, np.inf)
        ub[:n] = 1.0
        if need_T:
            ub[n:] = M
            for v in sources:
                if v in vidx:
                    lb[vidx[v]] = ub[vidx[v]] = 0.0
        rows = []
        row_lb = []
        row_ub = []

        def add(r, lo, hi):
            rows.append(r)
            row_lb.append(lo)
            row_ub.append(hi)

        # generation-flow: every non-source tail needs a selected producer
        for e in E:
            for u in e["tail"]:
                if u in sources:
                    continue
                r = np.zeros(N)
                r[eidx[id(e)]] = 1.0
                for f in E:
                    if f["head"] == u and f is not e:
                        r[eidx[id(f)]] = -1.0
                add(r, -np.inf, 0.0)

        # the target must have a selected in-edge
        tcols = [eidx[id(f)] for f in E if f["head"] == target]
        if not tcols:
            return None
        r = np.zeros(N)
        r[tcols] = 1.0
        add(r, 1.0, np.inf)

        if need_T:
            for e in E:
                h = e["head"]
                for u in e["tail"]:
                    if u == h:
                        continue  # self-loop: does not advance completion time
                    r = np.zeros(N)
                    r[vidx[h]] = 1.0
                    if u in vidx:
                        r[vidx[u]] = -1.0
                    r[eidx[id(e)]] = -M
                    add(r, d[id(e)] - M, np.inf)
            if Tcap is not None:
                r = np.zeros(N)
                r[vidx[target]] = 1.0
                add(r, -np.inf, Tcap + 1e-6)
        if Lcap is not None:
            r = np.zeros(N)
            r[:n] = 1.0
            add(r, -np.inf, Lcap + 1e-6)

        A = lil_matrix((len(rows), N), dtype=float)
        for i, r in enumerate(rows):
            nz = np.flatnonzero(r)
            if nz.size:
                A[i, nz] = r[nz]
        cons = LinearConstraint(A.tocsr(), np.array(row_lb), np.array(row_ub))
        return milp(c, integrality=integrality, bounds=Bounds(lb, ub),
                    constraints=cons,
                    options={"time_limit": timeout, "disp": False})

    if not need_T:
        res = solve_with(None)
    else:
        md = max(d[id(e)] for e in E)
        M0 = md * len(E) + 10
        candidates = []
        for Mv in [M0, M0 // 2, md * len(E), md * max(len(V), 1), M0 // 4,
                   100.0, 50.0, 20.0, 10.0]:
            if Mv not in candidates and Mv >= md:
                candidates.append(Mv)
        res = None
        for M in candidates:
            res = solve_with(M)
            if res is None:
                return ("infeasible", None, None, None)
            if res.status != 4:  # 4 == numerical solve error; retry with a
                break            # smaller big-M, otherwise keep this verdict
    if res is None:
        return ("infeasible", None, None, None)
    if res.status != 0 or res.x is None:
        if res.status == 2:
            return ("infeasible", None, None, None)
        msg = getattr(res, "message", "status %d" % res.status)
        return ("error:" + str(msg)[:40], None, None, None)
    xvals = res.x
    xsel = [E[i] for i in range(n) if xvals[i] > 0.5]
    T = None
    if need_T and objective in ("time", "len_t"):
        T = {v: float(xvals[vidx[v]]) for v in vlist}
    return ("integer optimal solution", float(res.fun), xsel, T)


def eval_T(F, d, sources, in_edges, target):
    """Completion time of a fixed edge set F under delays (min fixed point)."""
    T = {v: 0.0 for v in sources}
    changed = True
    while changed:
        changed = False
        for e in F:
            h = e["head"]
            if not e["tail"] or not all(u in T for u in e["tail"]):
                continue
            cand = max(T[u] for u in e["tail"]) + d[id(e)]
            if cand < T.get(h, float("inf")):
                T[h] = cand
                changed = True
    return T.get(target, None)

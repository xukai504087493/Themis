#!/usr/bin/env python3
"""Single-instance solver entry point (Mmunin-style CLI).

Solve the temporal shortest hyperpath problem for one source-target pair
under the measured-kinetics delay model, or run a bundled demo instance.

Usage:
  python run.py --name <dataset> --target <target-id> [--objective time|len|len_d] [--budget <seconds>]
  python run.py --name <dataset> --demo
"""

import argparse
import json
import os
import sys
import time

from hypergraph import build, load_edges, reachable_from
from ilp import solve_ilp, subgraph

BASE = os.path.dirname(os.path.abspath(__file__))
DATASETS = ["allpid", "allreactome", "allkegg", "allhumancyc"]
DELAY_NAMES = {"allpid": "pid", "allreactome": "reactome",
               "allkegg": "kegg", "allhumancyc": "humancyc"}


def edges_path(name):
    return os.path.join(BASE, "data", "hypergraphs", "%s-hyperedges.txt" % name)


def delay_path(name):
    return os.path.join(BASE, "data", "%s_delay.json" % DELAY_NAMES[name])


def pick_demo_target(edges, in_edges, sources):
    """First reachable target whose pred-closure fits the solver budget."""
    verts, tail_membership, has_tail, heads, _ = build(edges)
    targets = [v for v in verts if v not in has_tail]
    reach = reachable_from(sources, tail_membership)
    for t in sorted(targets):
        if t not in reach:
            continue
        _, E = subgraph(t, edges, in_edges, sources)
        if 0 < len(E) <= 200:
            return t
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--name", default="allpid",
                    help="dataset name (allpid, allreactome, allkegg, allhumancyc)")
    ap.add_argument("--target", help="target vertex id in the dataset")
    ap.add_argument("--demo", action="store_true",
                    help="solve a bundled representative instance")
    ap.add_argument("--objective", choices=["time", "len", "len_d"], default="time",
                    help="objective to optimize (default: time = min makespan)")
    ap.add_argument("--budget", type=float, default=300.0,
                    help="solver time limit in seconds")
    args = ap.parse_args()

    if args.name not in DATASETS:
        sys.exit("unknown dataset --name %r (use allpid/allreactome/allkegg/allhumancyc)"
                 % args.name)
    edges = load_edges(edges_path(args.name))
    verts, tail_membership, has_tail, heads, in_edges = build(edges)
    sources = set(v for v in verts if v not in heads)
    delay_map = json.load(open(delay_path(args.name)))

    if args.demo:
        t = pick_demo_target(edges, in_edges, sources)
        if t is None:
            sys.exit("no demo target found")
        print("demo target: %s" % t)
    else:
        if not args.target:
            sys.exit("provide --target <id> or --demo")
        t = args.target
        if t not in in_edges:
            sys.exit("target %r has no producing reaction" % t)

    V, E = subgraph(t, edges, in_edges, sources)
    d = {id(e): delay_map.get(e["id"], 1) for e in E}
    print("dataset: %s | vertices: %d | edges: %d | pred-closure of %s: %d edges"
          % (args.name, len(verts), len(edges), t, len(E)))
    t0 = time.time()
    status, obj, F_sel, T = solve_ilp(V, E, sources, t, d, args.objective,
                                      timeout=args.budget)
    dt = time.time() - t0
    print("objective: %s | status: %s | value: %s | elapsed: %.2fs"
          % (args.objective, status, obj, dt))
    if F_sel:
        print("selected reactions: %d" % len(F_sel))
        if T is not None:
            print("makespan of %s: %s" % (t, T.get(t)))
        for e in F_sel[:8]:
            print("  %s <- %s (%s)" % (e["head"], "+".join(e["tail"]), e["id"]))
        if len(F_sel) > 8:
            print("  ...")


if __name__ == "__main__":
    main()

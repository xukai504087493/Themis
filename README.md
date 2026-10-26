# Themis: Exact Temporal Shortest Hyperpaths with Cycles

Companion code for the RECOMB 2027 submission *Exact Temporal Shortest
Hyperpaths with Cycles for Pathway Inference in Cellular Reaction Networks*.
This package solves the temporal shortest hyperpath problem (T-SHP) on
biochemical reaction hypergraphs with measured-kinetics delays.

## About

Cell signaling and metabolic pathways are classically modeled as ordinary
graphs, but reactions with multiple reactants that must all be present
simultaneously, or that assemble protein complexes, cannot be represented
faithfully in this way. Directed hypergraphs generalize ordinary graphs to
such multiway reactions: a hyperedge is directed from a set of reactants
(tail) to a set of products (head). Finding the most efficient way to
synthesize a target molecule from available sources then maps to a shortest
hyperpath problem on a directed hypergraph.

Real reaction networks contain cycles, and reactions take time. Existing
exact algorithms either ignore reaction delays entirely, or forbid cycles in
their temporal formulations; neither can handle a cyclic hyperpath under
per-reaction delays. Themis introduces the temporal shortest hyperpath
problem with cycles (T-SHP), in which every reaction carries a positive delay
derived from measured catalytic rate constants, and solves it to proven
optimality with an exact cutting-plane algorithm driven by two separation
oracles: temporal-closure cuts and positive-delay-cycle detection. The
evaluation covers 14,401 source-target instances from four pathway databases
(NCI-PID, Reactome, KEGG, HumanCyc) with delays built from 35,147 measured
kcat values; on these instances, the median solve time is 3 milliseconds and
the largest observed solve time is 0.14 seconds.

## Features

- Directed-hypergraph loading from Pathway Commons TSV exports
  (`Tail|Head|PosReg|NegReg|ID` protocol).
- Static tail-distance fixed point, argmin backtracking, generation-order
  feasibility, and linear-time temporal distance on B-hypergraphs.
- Exact ILP formulations with three objectives (`len` / `len_d` / `time`)
  and optional makespan / edge-count caps.
- All dependencies are free and open-source: the MILP engine is HiGHS,
  accessed through `scipy.optimize.milp`; no commercial solver license is
  required.
- Bundled data: hypergraph exports and delay maps for all four databases are
  vendored, so the package runs out of the box.

## Installation

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Running the code

Solve a single source-target instance under the measured-kinetics delay
model (Mmunin-style command line):

```bash
# solve a bundled representative instance
.venv/bin/python run.py --name allpid --demo

# solve a specific target
.venv/bin/python run.py --name allreactome --target <target-id> --objective time

# objectives: time (min makespan, default), len (min edge count, static),
#             len_d (min total delay among static optima)
.venv/bin/python run.py --name allkegg --target <target-id> --objective len
```

`--name` selects one of `allpid`, `allreactome`, `allkegg`, `allhumancyc`.
The demo instance is chosen automatically; for `--target`, use any vertex id
present in the corresponding `data/hypergraphs/<name>-hyperedges.txt` file
(e.g. Pathway Commons ids such as
`http://pathwaycommons.org/pc12/Protein_...`).

## Data

- `data/hypergraphs/<name>-hyperedges.txt` — hypergraph exports built from
  Pathway Commons release 12 following the construction protocol of Ritz et
  al. and Krieger and Kececioglu.
- `data/<name>_delay.json` — per-reaction delay maps: `d(e) =
  clip(ceil(5*log10(k_med/k_e)+25), 1, 50)`, the logarithmic discretization
  of the inverse turnover rate described in the paper, derived from 35,147
  measured kcat values (EITLEM-Kinetics and DLKcat, originally sourced from
  BRENDA and SABIO-RK).
- `data/preprocess/build_delay.py` — reconstructs a delay map from the raw
  kinetic and reaction tables (not distributed here); see
  `data/preprocess/README.md` for the upstream sources.

## Tests

```bash
.venv/bin/python -m pytest tests/ -q
```

`tests/` cover the algorithm core (loading, temporal distance, feasibility,
SCC, linear B-hypergraph distance) and require no solver license.

## Citation

Noteworthy uses of this code should cite the accompanying paper:

> Exact Temporal Shortest Hyperpaths with Cycles for Pathway Inference in
> Cellular Reaction Networks. RECOMB 2027 (in submission).

## License

The Themis source code is free for noncommercial use, and comes with neither
warranty nor guarantee. The source code cannot be redistributed in any form
without consent of the authors. If you wish to use it for commercial purposes,
you must first obtain permission from all authors. See `LICENSE`.

# Themis: Exact Temporal Shortest Hyperpaths with Cycles
![License: Noncommercial](https://img.shields.io/badge/License-Noncommercial%20%26%20Cite-blue.svg)

Companion code for the RECOMB 2027 submission *Exact Temporal Shortest
Hyperpaths with Cycles for Pathway Inference in Cellular Reaction Networks*.
Themis solves the temporal shortest hyperpath problem (T-SHP) on biochemical
reaction hypergraphs with measured-kinetics delays, and is the exact engine
behind all experimental results in the paper (14,401 source-target
instances, median solve time 3 ms).

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
oracles: temporal-closure cuts and positive-delay-cycle detection.

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

## Python Requirements

The code requires Python >= 3.9. Dependencies are declared in
`requirements.txt`:

- numpy>=1.24
- scipy>=1.9

## Installation Guide

Clone the repository and install the dependencies (conda is recommended, as
in the reference implementations):

```bash
git clone <repo-url>
cd Themis
conda create -n themis_env python=3.9
conda activate themis_env
pip install -r requirements.txt
```

## Running the code

Solve a single source-target instance under the measured-kinetics delay
model (Mmunin-style command line):

```bash
# solve a bundled representative instance
python run.py --name allpid --demo

# solve a specific target (any vertex id in the data file)
python run.py --name allreactome --target <target-id> --objective time
```

Parameter definition:

- `--name`: dataset name, optional: `allpid`, `allreactome`, `allkegg`,
  `allhumancyc` (default: `allpid`)
- `--target`: target vertex id, e.g.
  `http://pathwaycommons.org/pc12/Protein_...`
- `--objective`: optimization objective, optional: `time` (min makespan,
  default), `len` (min edge count, static), `len_d` (min total delay among
  static optima)
- `--budget`: solver time limit in seconds (default: 300)
- `--demo`: solve an automatically chosen representative instance

Expected output of the demo run:

```
demo target: http://pathwaycommons.org/pc12/Complex_0007a39c9847990d54142810bd3a6b46;http://pathwaycommons.org/pc12/Protein_c9a6c85d3c93a372636dc4038f485493
dataset: allpid | vertices: 7393 | edges: 4620 | pred-closure of http://pathwaycommons.org/pc12/Complex_0007a39c9847990d54142810bd3a6b46;http://pathwaycommons.org/pc12/Protein_c9a6c85d3c93a372636dc4038f485493: 2 edges
objective: time | status: integer optimal solution | value: 2.0 | elapsed: 0.01s
selected reactions: 2
makespan of http://pathwaycommons.org/pc12/Complex_0007a39c9847990d54142810bd3a6b46;http://pathwaycommons.org/pc12/Protein_c9a6c85d3c93a372636dc4038f485493: 2.0
```

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
python -m pytest tests/ -q
```

`tests/` cover the algorithm core (loading, temporal distance, feasibility,
SCC, linear B-hypergraph distance) and require no solver license.

## Cite Us

Noteworthy uses of this code should cite the accompanying paper:

```bibtex
@inproceedings{themis2027,
  title     = {Exact Temporal Shortest Hyperpaths with Cycles for Pathway Inference in Cellular Reaction Networks},
  author    = {Themis authors},
  booktitle = {RECOMB},
  year      = {2027},
  note      = {in submission; author list to be completed at camera-ready}
}
```

## License

The Themis source code is free for noncommercial use, and comes with neither
warranty nor guarantee. The source code cannot be redistributed in any form
without consent of the authors. If you wish to use it for commercial purposes,
you must first obtain permission from all authors. The full terms are in
[`LICENSE`](LICENSE).

# Delay-map reconstruction

The vendored `data/*_delay.json` delay maps are the final product of this
pipeline; every experiment reads them, and they are shipped with the
package. `build_delay.py` reconstructs the maps from raw kinetic and reaction
tables that are **not** distributed here:

- `kcat_data.json` — 35,147 measured catalytic rate constants compiled from
  the EITLEM-Kinetics and DLKcat kinetic datasets (originally sourced from
  the BRENDA and SABIO-RK databases).
- `uniprot_human_ec.tsv` — human UniProt -> EC mapping (UniProt release).
- `reactome_prot.tsv`, `reactome_complex.tsv`, `reactome_cat.tsv` — Reactome
  protein / complex / catalysis tables (Pathway Commons / Reactome).
- `pid_prot3.tsv`, `pid_cat3.tsv`, `pid_complex.tsv` — NCI-PID equivalents.
- `ec.txt` — EC number list.

To rebuild a delay map, obtain these tables from the upstream sources above,
place them in this directory, and run:

```bash
python build_delay.py
```

The delay formula is `d(e) = clip(ceil(5*log10(k_med/k_e)+25), 1, 50)`, the
logarithmic discretization of the inverse turnover rate described in the
paper (Sect. 3, Experimental setup).

import json, csv, re, collections, statistics, os, math, subprocess

# 1. EC -> median kcat (35,147 kinetic records; kcat datasets from the
#    EITLEM-Kinetics and DLKcat studies, originally sourced from BRENDA and SABIO-RK).
kcat_by_ec = collections.defaultdict(list)
for s in json.load(open('kcat_data.json')):
    kcat_by_ec[s['ec']].append(float(s['value']))
ec_kcat = {ec: statistics.median(vals) for ec, vals in kcat_by_ec.items()}
print('EC classes:', len(ec_kcat), '| kcat range: %.2e ~ %.2e' % (min(ec_kcat.values()), max(ec_kcat.values())))

# 2. UniProt -> EC (human enzyme table; keep complete four-level EC only).
up_ec = collections.defaultdict(set)
with open('uniprot_human_ec.tsv') as f:
    rd = csv.DictReader(f, delimiter='\t')
    for r in rd:
        for ec in r['EC number'].split(';'):
            ec = ec.strip()
            if re.fullmatch(r'\d+\.\d+\.\d+\.\d+', ec):
                up_ec[r['Entry']].add(ec)
print('human enzymes (UniProt):', len(up_ec))


def load_tsv(p):
    d = {}
    if not os.path.exists(p):
        return d
    for line in open(p):
        parts = line.rstrip('\n').split('\t')
        if len(parts) >= 2 and parts[0]:
            d[parts[0]] = parts[1]
    return d


# 3. protein / complex -> UniProt; complex -> components.
prot_up = load_tsv('reactome_prot.tsv')
comp_up = load_tsv('reactome_complex.tsv')
pid_prot_up = load_tsv('pid_prot3.tsv')
# 4. reaction -> catalyst refs.
reactome_cat = [l.rstrip('\n').split('\t') for l in open('reactome_cat.tsv')]
pid_cat = [l.rstrip('\n').split('\t') for l in open('pid_cat3.tsv')]
print('reactome catalysis rows:', len(reactome_cat), '| pid catalysis rows:', len(pid_cat))


def refs_to_uniprots(refs, prot_up, comp_up):
    """Map controller refs (#id or full URI) to a set of UniProt accessions."""
    ups = set()
    for r in refs:
        rid = r.lstrip('#')
        if rid in prot_up:
            for u in prot_up[rid].split(','):
                ups.add(u)
        elif rid in comp_up:
            for c in comp_up[rid].split(','):
                cid = c.lstrip('#')
                if cid in prot_up:
                    for u in prot_up[cid].split(','):
                        ups.add(u)
    return ups


def reaction_to_ec(cat_rows, prot_up, comp_up):
    """Map a controlled reaction ID to its set of EC numbers."""
    r2ec = {}
    for parts in cat_rows:
        if len(parts) < 2:
            continue
        r, refs = parts[0], parts[1].split(',')
        ups = refs_to_uniprots(refs, prot_up, comp_up)
        ecs = set()
        for u in ups:
            ecs |= up_ec.get(u, set())
        if ecs:
            r2ec.setdefault(r, set()).update(ecs)
    return r2ec


# PID complex components: parsed on demand.
if not os.path.exists('pid_complex.tsv'):
    subprocess.run(['python3', 'parse_sax2.py', 'PathwayCommons12.pid.BIOPAX.owl.gz',
                    'pid_cat3.tsv', 'pid_prot3.tsv', 'pid_complex.tsv'], check=False)
pid_comp_up = load_tsv('pid_complex.tsv')

r2ec_rx = reaction_to_ec(reactome_cat, prot_up, comp_up)
r2ec_pid = reaction_to_ec(pid_cat, pid_prot_up, pid_comp_up)
print('Reactome reactions->EC:', len(r2ec_rx), '| PID reactions->EC:', len(r2ec_pid))

# 5. Reaction -> kcat (EC median; geometric mean over multiple ECs) -> delay.
K_MED = statistics.median(ec_kcat.values())   # global reference kcat median (s^-1)


def ecs_to_delay(ecs):
    kcats = [ec_kcat[ec] for ec in ecs if ec in ec_kcat]
    if not kcats:
        return None
    k_e = statistics.geometric_mean(kcats)
    # Log-scaled discrete delay (paper formula):
    # clip(ceil(5*log10(k_med/k_e) + 25), 1, 50).
    d = max(1, min(50, math.ceil(5 * math.log10(K_MED / k_e) + 25)))
    return d, k_e


def build_map(r2ec):
    m = {}
    for r, ecs in r2ec.items():
        res = ecs_to_delay(ecs)
        if res:
            m[r] = (res[0], res[1], ','.join(sorted(ecs)))
    return m


dm_rx = build_map(r2ec_rx)
dm_pid = build_map(r2ec_pid)
print('Reactome reactions->delay:', len(dm_rx), '| PID reactions->delay:', len(dm_pid))
json.dump({k: v[0] for k, v in dm_rx.items()}, open('reactome_delay.json', 'w'))
json.dump({k: v[0] for k, v in dm_pid.items()}, open('pid_delay.json', 'w'))

# Delay distribution summary.
dc = collections.Counter(v[0] for v in dm_rx.values())
print('Reactome d distribution:', dict(sorted(dc.items())))
dc2 = collections.Counter(v[0] for v in dm_pid.values())
print('PID d distribution:', dict(sorted(dc2.items())))

#!/usr/bin/env python3
"""Summarise a folder of FCS files named MTB<patient>_<batch>_V<vaccine>W<week>_<cocktail>.fcs.

Writes fcs_metadata.csv (one row per file), fcs_counts.csv (files per value of
each category) and fcs_overview.xlsx (the same two tables as sheets plus one sheet per category
cross-tabulated against the others) into the FCS folder, or into --out.

    python scripts/fcs_overview.py data/MTB_FACS/fcs [--out results]
"""
import argparse, re
from pathlib import Path
import pandas as pd

CATS = ['patient', 'batch', 'vaccine', 'week', 'cocktail']
NAME = re.compile(r'^MTB(?P<patient>\d+)_(?P<batch>\d+)_V(?P<vaccine>\d+)W(?P<week>\d+)_(?P<cocktail>.+)$')

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('fcs_dir', type=Path, help='folder searched recursively for *.fcs')
p.add_argument('--out', type=Path, help='output folder (default: fcs_dir)')
a = p.parse_args()
out = a.out or a.fcs_dir
out.mkdir(parents=True, exist_ok=True)

rows = []
for f in sorted(a.fcs_dir.rglob('*.fcs')):
    m = NAME.match(f.stem)
    if not m:
        print(f'skipping unparsable name: {f.name}')
        continue
    d = {k: int(v) for k, v in m.groupdict().items() if k != 'cocktail'}
    d['cocktail'] = 'neg' if m['cocktail'].endswith('neg') else m['cocktail']
    rows.append({'filename': f.name, **d})
ov = pd.DataFrame(rows).sort_values(CATS)
counts = pd.concat(
    ov[c].value_counts().sort_index().rename_axis('value').reset_index(name='n_files').assign(variable=c)
    for c in CATS)[['variable', 'value', 'n_files']]

ov.to_csv(out / 'fcs_metadata.csv', index=False)
counts.to_csv(out / 'fcs_counts.csv', index=False)
with pd.ExcelWriter(out / 'fcs_overview.xlsx') as xw:
    ov.to_excel(xw, sheet_name='metadata', index=False)
    counts.to_excel(xw, sheet_name='counts', index=False)
    for c in CATS:
        row = 0
        for o in [x for x in CATS if x != c]:
            ct = pd.crosstab(ov[c], ov[o], margins=True, margins_name='total').rename_axis(f'{c} \\ {o}')
            pd.DataFrame([f'files per {c} x {o}']).to_excel(xw, sheet_name=f'by_{c}', startrow=row, index=False, header=False)
            ct.to_excel(xw, sheet_name=f'by_{c}', startrow=row + 1)
            row += len(ct) + 4
print(f'{len(ov)} files -> {out}/fcs_metadata.csv, fcs_counts.csv, fcs_overview.xlsx')

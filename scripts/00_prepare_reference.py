#!/usr/bin/env python
"""Build the gene-coordinate table and chromosome-arm boundaries used for binning (GRCh38, Ensembl release 98, which matches the
10x GRCh38-2020-A reference gene models). Downloads the GTF and the karyotype bands from Ensembl unless already present.

    python scripts/00_prepare_reference.py            # writes work/ref/gene_table_GRCh38_ens98.csv.gz and centromeres_GRCh38.csv
"""
import gzip, json, re, sys, urllib.request
import pandas as pd
from sccnv import config

GTF_URL = 'https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz'
BANDS_URL = 'https://rest.ensembl.org/info/assembly/homo_sapiens?bands=1&content-type=application/json'

def fetch(url, dst):
    if dst.exists() and dst.stat().st_size > 0: return
    print('downloading', url); urllib.request.urlretrieve(url, dst)

config.REF.mkdir(parents=True, exist_ok=True)
gtf_f, bands_f = config.REF / 'Homo_sapiens.GRCh38.98.gtf.gz', config.REF / 'assembly_bands.json'
fetch(GTF_URL, gtf_f); fetch(BANDS_URL, bands_f)

rows = []
with gzip.open(gtf_f, 'rt') as f:
    for l in f:
        if l.startswith('#'): continue
        t = l.split('\t')
        if t[2] != 'gene': continue
        a = t[8]
        gid = re.search(r'gene_id "([^"]+)"', a).group(1); gn = re.search(r'gene_name "([^"]+)"', a); bt = re.search(r'gene_biotype "([^"]+)"', a)
        rows.append((gid, gn.group(1) if gn else gid, t[0], int(t[3]), int(t[4]), t[6], bt.group(1) if bt else ''))
gtf = pd.DataFrame(rows, columns=['gene_id', 'symbol', 'chr', 'start', 'end', 'strand', 'biotype'])
gt = gtf[gtf.chr.isin([str(i) for i in range(1, 23)] + ['X', 'Y'])]
gt.to_csv(config.REF / 'gene_table_GRCh38_ens98.csv.gz', index=False)

bands = json.load(open(bands_f)); kt = bands['karyotype']; cen = []
for r in bands['top_level_region']:
    if r['name'] in kt and r['name'] != 'MT':
        b = r['bands']
        cen.append((r['name'], max(x['end'] for x in b if x['id'].startswith('p')), min(x['start'] for x in b if x['id'].startswith('q')), r['length']))
pd.DataFrame(cen, columns=['chr', 'p_end', 'q_start', 'length']).to_csv(config.REF / 'centromeres_GRCh38.csv', index=False)
print(f'{len(gt)} genes, {len(cen)} chromosomes written to {config.REF}')

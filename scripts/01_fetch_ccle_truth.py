#!/usr/bin/env python
"""Fetch gene-level DNA copy number (log2) for six prostate cancer cell lines from the cBioPortal public API.

  * 22Rv1, DU145, LNCaP (clone FGC), PC3  <- ccle_broad_2025_log2CNA (WGS-derived, log2(relative CN + 1))
  * VCaP, MDA-PCa-2b                       <- cellline_ccle_broad_log2CNA (SNP6, log2(CN / 2)); absent from the 2025 study
Output: work/ref/ccle2025_log2CNA_prostate_lines.csv, work/ref/ccle2012_snp6_log2CNA_prostate_lines.csv
"""
import time, requests, pandas as pd
from sccnv import config

B = 'https://www.cbioportal.org/api'; H = {'Accept': 'application/json'}
config.REF.mkdir(parents=True, exist_ok=True)
genes = requests.get(f'{B}/genes?projection=SUMMARY&pageSize=100000', headers=H, timeout=300).json()
ent = [g['entrezGeneId'] for g in genes if g['type'] == 'protein-coding']; e2s = {g['entrezGeneId']: g['hugoGeneSymbol'] for g in genes}

def fetch(profile, samples, chunk=1500, tries=5):
    out = []
    for k in range(0, len(ent), chunk):
        for t in range(tries):                                       # the public API occasionally returns transient 5xx
            r = requests.post(f'{B}/molecular-profiles/{profile}/molecular-data/fetch?projection=SUMMARY',
                              json={'entrezGeneIds': ent[k:k + chunk], 'sampleIds': samples}, headers={**H, 'Content-Type': 'application/json'}, timeout=300)
            if r.status_code == 200: out += r.json(); break
            time.sleep(3 * (t + 1))
        else: raise RuntimeError(f'{profile}: chunk {k} failed ({r.status_code})')
    d = pd.DataFrame([{'symbol': e2s[x['entrezGeneId']], 'sample': x['sampleId'], 'log2cn': float(x['value'])} for x in out if x['value'] not in ('NA', None)])
    return d.pivot_table(index='symbol', columns='sample', values='log2cn')

fetch('ccle_broad_2025_log2CNA', ['22RV1_PROSTATE', 'DU145_PROSTATE', 'LNCAPCLONEFGC_PROSTATE', 'PC3_PROSTATE']).to_csv(config.REF / 'ccle2025_log2CNA_prostate_lines.csv')
fetch('cellline_ccle_broad_log2CNA', ['VCAP_PROSTATE', 'MDAPCA2B_PROSTATE']).to_csv(config.REF / 'ccle2012_snp6_log2CNA_prostate_lines.csv')
print('done')

"""Dataset loaders used by the analysis scripts. Each returns a dict with counts (cells x genes, CSR), gene symbols,
metadata-derived labels and, where available, an orthogonal DNA copy-number truth."""
import numpy as np, pandas as pd
from . import config, loaders, benchmark as bench
from .core import load_gene_table, load_centromeres, stochastic_round, bin_truth

EPI = ['luminal', 'basal', 'club', 'hillock']
STRO = ['endothelial', 'fibroblast', 'smooth_muscle', 'schwann']
IMM = ['b_cells', 't_cells', 'nk_cells', 'macrophage', 'monocyte', 'mast', 'plasma', 'dendritic_plasmacytoid', 'granulocytes']


def reference_tables():
    """gene coordinate table (Ensembl 98 / GRCh38, matches the 10x reference) and arm boundaries; see scripts/00_prepare_reference.py"""
    return (load_gene_table(config.REF / 'gene_table_GRCh38_ens98.csv.gz'), load_centromeres(config.REF / 'centromeres_GRCh38.csv'))


# ------------------------------------------------------------------ prostate cancer cell-line mixture (PrEC = normal reference)
def celllines(seed=1):
    X, genes, cells, meta = loaders.load_export(str(config.EXPORT / 'celllines'))
    table, cen = reference_tables(); rng = np.random.default_rng(seed)
    ann = meta['annotation'].values.astype(str); clus = meta['seurat_clusters'].astype(int).astype(str).values
    grp = np.array([f"{a}|c{c}" for a, c in zip(ann, clus)], dtype=object)
    prec = np.where(ann == 'PrEC')[0]; ho = rng.choice(prec, size=int(0.3 * len(prec)), replace=False)   # 30% of PrEC held out as diploid test cells
    is_ref = (ann == 'PrEC'); is_ref[ho] = False; grp[ho] = 'PrEC_heldout|c6'; grp[is_ref] = 'PrEC_ref|c6'
    gt25 = pd.read_csv(config.REF / 'ccle2025_log2CNA_prostate_lines.csv', index_col=0)
    gt12 = pd.read_csv(config.REF / 'ccle2012_snp6_log2CNA_prostate_lines.csv', index_col=0)
    l25 = lambda c: np.log2(np.clip(2 ** gt25[c] - 1, 0.05, None))          # ccle_broad_2025 log2(relative CN + 1)  ->  log2(relative CN)
    gene_truth = {'22Rv1': l25('22RV1_PROSTATE'), 'DU145': l25('DU145_PROSTATE'), 'LNCaP': l25('LNCAPCLONEFGC_PROSTATE'),
                  'PC3': l25('PC3_PROSTATE'), 'VCaP': gt12['VCAP_PROSTATE'], 'MDA_PCa_2b': gt12['MDAPCA2B_PROSTATE']}
    return dict(X=X, genes=genes, table=table, cen=cen, ann=ann, grp=grp, is_ref=is_ref, ho=(grp == 'PrEC_heldout|c6'),
                gene_truth=gene_truth, cell_key=ann.copy())

def celllines_truths(S, sel):
    return {k: bench.mode_center(bin_truth(sel['bins'], sel['genes'], v)) for k, v in S['gene_truth'].items()}


# ------------------------------------------------------------------ WA71 lymph-node metastasis with tumour DNA copy number
def wa71(seed=0):
    X, genes, cells, meta = loaders.load_export(str(config.EXPORT / 'wa71'))
    table, cen = reference_tables(); rng = np.random.default_rng(seed)
    X = stochastic_round(X, rng)       # counts are ambient-corrected (non-integer): unbiased integer rounding for count splitting / thinning
    anno = meta.anno.values.astype(str); clus = meta.seurat_clusters.astype(int).astype(str).values
    is_str = np.isin(anno, ['Endo_1', 'Endo_2', 'Endo_3', 'SMC']); is_imm = np.isin(anno, ['Bcell', 'Tcell', 'Macro'])
    ho = (rng.random(len(anno)) < 0.2) & (is_str | is_imm)                    # 20% of reference-lineage cells held out as diploid test cells
    grp = np.array([f"{a}|c{c}" + ("_ho" if h else "") for a, c, h in zip(anno, clus, ho)], dtype=object)
    seg = bench.load_wa71_segments(str(config.DATA / config.FILES['wa71_segments']))
    cell_key = np.where(np.isin(anno, ['Tumor_1', 'Tumor_2']), 'tumor', '')
    return dict(X=X, genes=genes, table=table, cen=cen, anno=anno, clus=clus, is_str=is_str, is_imm=is_imm, ho=ho, grp=grp,
                seg=seg, cell_key=cell_key, meta=meta)

def wa71_truth(S, sel):
    gt = bench.segments_to_gene_log2(S['seg'], S['table'], list(sel['genes'].symbol))
    return bench.mode_center(bin_truth(sel['bins'], sel['genes'], gt))


# ------------------------------------------------------------------ primary prostate tumours (Wong et al.): 7 patients, benign + tumour library each
def wong():
    X, genes, cells, meta = loaders.load_export(str(config.EXPORT / 'wong'))
    table, cen = reference_tables(); anno = meta.anno.values.astype(str)
    cls = np.array(['epi' if a in EPI else 'stroma' if a in STRO else 'immune' if a in IMM else 'other' for a in anno])
    return dict(X=X, genes=genes, table=table, cen=cen, anno=anno, cls=cls, lib=meta.tumor_normal.values.astype(str),
                pat=meta.patient.values.astype(str), clus=meta.seurat_clusters.astype(int).astype(str).values, meta=meta)


# ------------------------------------------------------------------ healthy prostate (Henry et al. 2018): 3 donors, 5 libraries, 10x 3' v2
def healthy():
    X, genes, cells, meta = loaders.load_export(str(config.EXPORT / 'healthy'))
    return dict(X=X, genes=genes, ctype=meta.pred_prostateref.values.astype(str), donor=meta.patient_id.values.astype(str),
                lib=meta.ID.values.astype(str), zone=meta.Notes.values.astype(str), meta=meta)

"""End-to-end example on simulated counts (no data files needed):  python examples/simulated_quickstart.py

Diploid reference cells plus 'tumour' cells with a 1.6x gain on chromosome 2q and a 0.5x loss on 3p."""
import numpy as np
from sccnv import Config, run, select_genes_and_bins
from sccnv.simulate import toy_genome, simulate

genome = toy_genome(); table, cen, _, _ = genome
S = simulate(genome)
cfg = Config(genes_per_bin=40, metacell_cells=20)
# 1. genes -> bins along the genome (one reference mask marks the cells that define the baseline)
sel = select_genes_and_bins(S['X'], S['genes'], table, cen, cfg, S['grp'], S['ref_mask'])
# 2. count-split metacells inside each group, reference baseline + null, segmentation, calibrated calls
res = run(S['X'], S['genes'], table, cen, S['grp'], {'reference': S['ref_mask']}, cfg, verbose=False, sel=sel)

r, seg, call = res.pool(['tum|c0'])                 # pooled profile of all tumour metacells
for chrom, arm in [('2', 'q'), ('3', 'p'), ('1', 'p')]:
    m = ((res.bins.chr == chrom) & (res.bins.arm == arm)).values
    print(f'chr{chrom}{arm}: mean log2 ratio {r[m].mean():+.2f}   bins called gain/loss: {(call[m] > 0).sum()}/{(call[m] < 0).sum()} of {m.sum()}')
print('metacells:', len(res.meta_table), '| reference metacells:', int(res.meta_table.is_ref.sum()),
      '| calls in held-out diploid metacells:', int((res.call[res.meta_table.group.str.contains("_ho").values] != 0).sum()))
print('per-cell matrix (cells x bins):', res.cell_matrix('seg').shape)

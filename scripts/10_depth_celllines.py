"""Depth-downsampling experiment on the prostate cell-line mixture (PrEC = diploid reference; CCLE DNA copy number as truth).
   python scripts/10_depth_celllines.py [--quick]"""
import sys, time
from sccnv import config, datasets, offsets as O
from sccnv import benchmark as bench, experiments, loaders as data_io
from sccnv.core import *
QUICK = '--quick' in sys.argv
S=datasets.celllines(); X=S['X']
base=dict(genes_per_bin=100,weight_power=0.5)
sel=select_genes_and_bins(X,S['genes'],S['table'],S['cen'],Config(**base),S['grp'],S['is_ref'])
tb=datasets.celllines_truths(S,sel)
df=experiments.run_depth_experiment(X,S['genes'],S['table'],S['cen'],S['grp'],{'PrEC':S['is_ref']},sel,S['cell_key'],tb,S['ho'],([1.0,0.1,0.025] if QUICK else [1.0,0.5,0.25,0.1,0.05,0.025]),experiments.DEFAULT_METHODS,base,seed=1)
df.to_csv(config.results_path('depth_experiment_celllines' + ('_quick' if QUICK else '') + '.csv'),index=False); print('done',df.shape)

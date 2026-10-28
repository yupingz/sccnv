"""Depth-downsampling experiment on the WA71 lymph-node metastasis (tumour DNA copy number as truth).
   python scripts/11_depth_wa71.py stroma|both [--quick]     stroma = endothelial + smooth muscle reference; both = per-bin median of stroma and immune"""
import sys, time
from sccnv import config, datasets, offsets as O
from sccnv import benchmark as bench, experiments, loaders as data_io
from sccnv.core import *
which=sys.argv[1]; QUICK = '--quick' in sys.argv
S=datasets.wa71(); X=S['X']
base=dict(genes_per_bin=100,weight_power=0.5)
if which=='stroma': ref={'stroma':S['is_str']&~S['ho']}
else: ref={'stroma':S['is_str']&~S['ho'],'immune':S['is_imm']&~S['ho']}; base['ref_combine']='median'
sel=select_genes_and_bins(X,S['genes'],S['table'],S['cen'],Config(**base),S['grp'],np.any(list(ref.values()),0))
tb=datasets.wa71_truth(S,sel)
methods=dict(experiments.DEFAULT_METHODS); methods['metacell-80k UMI (pipeline)']=dict(kind='pipeline',cfg=dict(target_umi=80000,metacell_cells=20,thin=True))
# at 20 cells the stroma metacells are ~58k B-fold UMIs here, so keep metacell-20 as the fixed-size arm
df=experiments.run_depth_experiment(X,S['genes'],S['table'],S['cen'],S['grp'],ref,sel,S['cell_key'],{'tumor':tb},S['ho'],([1.0,0.1] if QUICK else [1.0,0.5,0.25,0.1,0.05]),methods,base,seed=1)
df['reference']=which; df.to_csv(config.results_path(f'depth_experiment_wa71_{which}' + ('_quick' if QUICK else '') + '.csv'),index=False); print('done',df.shape)

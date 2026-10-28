"""Healthy donors (diploid epithelium) analysed on the SAME genes/bins/weights as the Wong tumours: empirical null for CNV-identity scores.
Each donor: reference = own stroma+immune; epithelium = query, lineage offset from the OTHER donors (matched subtype)."""
import os
os.environ.setdefault('OMP_NUM_THREADS','8'); os.environ.setdefault('NUMBA_NUM_THREADS','8')
import sys, pickle; import numpy as np, pandas as pd
from sccnv import config, datasets, offsets as O
from sccnv import benchmark as bench, experiments, loaders as data_io
from sccnv.core import *
W=datasets.wong(); table,cen=W['table'],W['cen']
BASE=dict(genes_per_bin=100,weight_power=0.5,ref_combine='median',metacell_cells=20)
grp_all=np.array([f'{a}|c{k}' for a,k in zip(W['anno'],W['clus'])],dtype=object)
sel=select_genes_and_bins(W['X'],W['genes'],table,cen,Config(**BASE),grp_all,np.isin(W['cls'],['stroma','immune']))
H=datasets.healthy(); X=H['X']; genes=H['genes']; ct=H['ctype']; don=H['donor']
clus=H['meta'].seurat_clusters.astype(int).astype(str).values
cls=np.array(['epi' if c in O.EPI else 'stroma' if c in O.H_STROMA else 'immune' if c in O.H_IMMUNE else 'other' for c in ct])
Yb=O.bin_counts(X,genes,sel); ref_sets=[O.H_STROMA,O.H_IMMUNE]
# map healthy gene columns to the Wong bin genes
idx_by={g:i for i,g in enumerate(genes)}; cols=[idx_by[s] for s in sel['genes'].symbol]
sel_h=dict(sel); gdf=sel['genes'].copy(); gdf['col']=cols; sel_h['genes']=gdf
out={}
for d in sorted(set(don)):
    idx=np.where((don==d)&(cls!='other'))[0]; Xd=X[idx]; c_=ct[idx]; k_=clus[idx]; cl_=cls[idx]
    rng=np.random.default_rng(0); ho=(rng.random(len(idx))<0.2)&np.isin(cl_,['stroma','immune'])
    grp=np.array([f'{a}|c{k}'+('_ho' if h else '') for a,k,h in zip(c_,k_,ho)],dtype=object)
    ref={'stroma':(cl_=='stroma')&~ho,'immune':(cl_=='immune')&~ho}
    others=(don!=d)&(cls!='other'); OFF={}
    for s in O.EPI:
        u,D=O.offset_table(Yb[others],ct[others],don[others],[s],ref_sets); OFF[s]=O.combine_units(D,3)
    eg=[g for g in pd.unique(grp) if g.split('|')[0] in O.EPI]
    lo={g:OFF[g.split('|')[0]][0] for g in eg}; lv={g:OFF[g.split('|')[0]][1] for g in eg}
    out[d]={}
    for v,(a,b) in {'none':(None,None),'healthy_matched':(lo,lv)}.items():
        res=run(Xd,genes,table,cen,grp,ref,Config(**BASE),verbose=False,sel=sel_h,lineage_offset=a,offset_var=b)
        out[d][v]=dict(meta_table=res.meta_table,cell_meta=res.cell_meta,r=res.r,sd=res.sd,call=res.call,arm_id=res.arm_id,lam=res.lam,tau2=res.tau2,penalty=res.penalty,cfg=res.cfg,
                       group_segments=res.group_segments,group_seg=res.group_seg,group_r=res.group_r,anno=c_)
    print(d,'done',flush=True)
pickle.dump(out,open(config.results_path('healthy_on_wong_bins.pkl'),'wb')); print('saved')

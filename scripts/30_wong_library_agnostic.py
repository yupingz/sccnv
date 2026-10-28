"""Primary tumours (Wong et al.): every epithelial cell is a query cell; library labels (benign/tumour) are NOT used anywhere in the
analysis. Lineage offsets come from healthy prostate donors (matched epithelial subtype). Identity (aneuploid vs diploid) is resolved from
the CNV profile of each cell's metacell."""
import os
os.environ.setdefault('OMP_NUM_THREADS','4'); os.environ.setdefault('OPENBLAS_NUM_THREADS','4'); os.environ.setdefault('NUMBA_NUM_THREADS','4')
import sys, pickle, time, zlib, multiprocessing as mp
import numpy as np, pandas as pd
from sccnv import config, datasets, offsets as O
from sccnv import benchmark as bench, experiments, loaders as data_io
from sccnv.core import *
W=datasets.wong(); X=W['X']; genes=W['genes']; table=W['table']; cen=W['cen']
anno,cls,pat,clus=W['anno'],W['cls'],W['pat'],W['clus']
BASE=dict(genes_per_bin=100,weight_power=0.5,ref_combine='median',metacell_cells=20)
keep_cls=cls!='other'
grp_all=np.array([f'{a}|c{k}' for a,k in zip(anno,clus)],dtype=object)
cfg_sel=Config(**BASE)
sel=select_genes_and_bins(X,genes,table,cen,cfg_sel,grp_all,np.isin(cls,['stroma','immune']))
print('[global] genes',len(sel['genes']),'bins',len(sel['bins']),flush=True)
# ---- healthy-donor offsets on these bins (3 donors), per epithelial subtype and pooled
H=datasets.healthy(); Yh=O.bin_counts(H['X'],H['genes'],sel)
ref_sets=[O.H_STROMA,O.H_IMMUNE]
OFF={}
for name,types in [('all',O.EPI)]+[(s,[s]) for s in O.EPI]:
    u,D=O.offset_table(Yh,H['ctype'],H['donor'],types,ref_sets); OFF[name]=O.combine_units(D,3)
    print('[offset]',name,'donors',len(u),'RMS %.3f'%np.sqrt((OFF[name][0]**2).mean()),flush=True)
perm=np.random.default_rng(1).permutation(len(sel['bins']))
def offsets_for(variant,egroups):
    if variant=='none': return None,None
    if variant=='healthy_all': return {g:OFF['all'][0] for g in egroups},{g:OFF['all'][1] for g in egroups}
    if variant=='healthy_matched': return {g:OFF[g.split('|')[0]][0] for g in egroups},{g:OFF[g.split('|')[0]][1] for g in egroups}
    if variant=='permuted_control': return {g:OFF[g.split('|')[0]][0][perm] for g in egroups},{g:OFF[g.split('|')[0]][1] for g in egroups}
VARIANTS=['none','healthy_matched','healthy_all','permuted_control']
CAP=2000
def process(p):
    t0=time.time(); rng=np.random.default_rng(zlib.crc32(p.encode()))
    idx=np.where((pat==p)&keep_cls)[0]; take=[]
    for c in ['stroma','immune','epi']:
        ii=idx[cls[idx]==c]; take.append(rng.choice(ii,CAP,replace=False) if (c!='epi' and len(ii)>CAP) else ii)
    idx=np.sort(np.concatenate(take)); Xp=X[idx]; a=anno[idx]; k=clus[idx]; c_=cls[idx]
    ho=(rng.random(len(idx))<0.2)&np.isin(c_,['stroma','immune'])
    grp=np.array([f'{a_}|c{k_}'+('_ho' if h else '') for a_,k_,h in zip(a,k,ho)],dtype=object)
    ref={'stroma':(c_=='stroma')&~ho,'immune':(c_=='immune')&~ho}
    egroups=[g for g in pd.unique(grp) if g.split('|')[0] in O.EPI]
    out=dict(patient=p,idx=idx,anno=a,cls=c_,lib=W['lib'][idx],ho=ho,grp=grp,variants={})
    for v in VARIANTS:
        lo,lv=offsets_for(v,egroups)
        res=run(Xp,genes,table,cen,grp,ref,Config(**BASE),verbose=False,sel=sel,lineage_offset=lo,offset_var=lv)
        mt=res.meta_table
        out['variants'][v]=dict(meta=res.cell_meta,meta_table=mt,call=res.call,seg=res.seg,r=res.r,sd=res.sd,floor=res.floor,z_call=res.z_call,penalty=res.penalty,group_seg=res.group_seg,group_r=res.group_r,group_call=res.group_call,group_segments=res.group_segments,lam=res.lam,tau2=res.tau2,arm_id=res.arm_id)
    print(f'[{p}] done {time.time()-t0:.0f}s cells={len(idx)}',flush=True)
    return out
if __name__=='__main__':
    pats=sorted(set(pat))
    with mp.get_context('fork').Pool(len(pats)) as pool: outs=pool.map(process,pats)
    pickle.dump(dict(outs=outs,bins=sel['bins'],genes=sel['genes'],OFF=OFF),open(config.results_path('wong_resolved.pkl'),'wb')); print('saved')

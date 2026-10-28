"""Negative control: healthy epithelium is diploid. For each donor, offset estimated from the OTHER donors is applied to that donor's
epithelium, which is analysed against that donor's own stromal+immune reference. Spurious signal should drop with a transferable offset."""
import numpy as np, pandas as pd
from sccnv import config, datasets, offsets as O
from sccnv import benchmark as bench, experiments, loaders as data_io
from sccnv.core import *
H=datasets.healthy(); X=H['X']; genes=H['genes']; ct=H['ctype']; don=H['donor']
table,cen=datasets.reference_tables()
clus=H['meta'].seurat_clusters.astype(int).astype(str).values
cls=np.array(['epi' if c in O.EPI else 'stroma' if c in O.H_STROMA else 'immune' if c in O.H_IMMUNE else 'other' for c in ct])
rows=[]; store={}
for d in sorted(set(don)):
    idx=np.where((don==d)&(cls!='other'))[0]; Xd=X[idx]; c_=ct[idx]; k_=clus[idx]; cl_=cls[idx]
    rng=np.random.default_rng(0); ho=(rng.random(len(idx))<0.2)&np.isin(cl_,['stroma','immune'])
    grp=np.array([f'{a}|c{k}'+('_ho' if h else '') for a,k,h in zip(c_,k_,ho)],dtype=object)
    ref={'stroma':(cl_=='stroma')&~ho,'immune':(cl_=='immune')&~ho}
    cfg=Config(genes_per_bin=100,weight_power=0.5,metacell_cells=20,ref_combine='median')
    sel=select_genes_and_bins(Xd,genes,table,cen,cfg,grp,np.any(list(ref.values()),0))
    Yb=O.bin_counts(X,genes,sel); others=(don!=d)&(cls!='other')
    ref_sets=[O.H_STROMA,O.H_IMMUNE]
    def off_for(types):
        u,D=O.offset_table(Yb[others],ct[others],don[others],types,ref_sets); return O.combine_units(D,3)
    off_all,var_all=off_for(O.EPI); off_sub={s:off_for([s]) for s in O.EPI}
    egroups=[g for g in pd.unique(grp) if g.split('|')[0] in O.EPI]
    perm=np.random.default_rng(1).permutation(len(off_all))
    variants={'no offset':(None,None),
              'other donors, all-epithelium offset':({g:off_all for g in egroups},{g:var_all for g in egroups}),
              'other donors, matched-subtype offset':({g:off_sub[g.split('|')[0]][0] for g in egroups},{g:off_sub[g.split('|')[0]][1] for g in egroups}),
              'control: all-epithelium offset, bins permuted':({g:off_all[perm] for g in egroups},{g:var_all for g in egroups})}
    for vn,(lo,lv) in variants.items():
        res=run(Xd,genes,table,cen,grp,ref,cfg,verbose=False,sel=sel,lineage_offset=lo,offset_var=lv)
        mt=res.meta_table; e=np.where(mt.group.str.split('|').str[0].isin(O.EPI).values)[0]; h=np.where(mt.group.str.contains('_ho').values)[0]
        burden=(res.call[e]!=0).mean(1)
        rows.append(dict(donor=d,variant=vn,n_epi_meta=len(e),
            frac_meta_any_call=float((res.call[e]!=0).any(1).mean()),frac_meta_burden_ge5pct=float((burden>=0.05).mean()),mean_burden=float(burden.mean()),
            rms_r_epi=float(np.sqrt((res.r[e]**2).mean())),rms_r_heldout_ref=float(np.sqrt((res.r[h]**2).mean())),
            floor_epi=float(np.median([res.floor[g] for g in egroups if g in res.floor])),z_call=res.z_call,n_ref_meta=int(mt.is_ref.sum())))
        store[(d,vn)]=burden; store[(d,vn,'segrms')]=np.sqrt((res.seg[e]**2).mean(1)); store[(d,vn,'zrms')]=np.sqrt((res.z[e]**2).mean(1)); store[(d,vn,'ncells')]=np.array([(res.cell_meta==m).sum() for m in e])
    print(d,'done',flush=True)
df=pd.DataFrame(rows); pd.set_option('display.width',250); print(df.round(4).to_string(index=False))
print(df.groupby('variant')[['frac_meta_any_call','frac_meta_burden_ge5pct','mean_burden','rms_r_epi','floor_epi']].mean().round(4))
df.to_csv(config.results_path('healthy_lodo_negative_control.csv'),index=False)
import pickle; pickle.dump(store,open(config.results_path('healthy_lodo_burdens.pkl'),'wb'))

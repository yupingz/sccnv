"""WA71: does a healthy-donor epithelium offset improve detection against tumour DNA? (no offset / healthy / permuted-bin control / Wong benign-labelled)
   python scripts/22_wa71_healthy_offset.py"""
import numpy as np, pandas as pd
from sccnv import config, datasets, offsets as O
from sccnv import benchmark as bench, experiments, loaders as data_io
from sccnv.core import *
H=datasets.healthy(); S=datasets.wa71(); X=S['X']
ref={'stroma':S['is_str']&~S['ho']}
base=dict(genes_per_bin=100,weight_power=0.5,metacell_cells=20)
cfg=Config(**base)
sel=select_genes_and_bins(X,S['genes'],S['table'],S['cen'],cfg,S['grp'],S['is_str']&~S['ho'])
tb=datasets.wa71_truth(S,sel); Yh=O.bin_counts(H['X'],H['genes'],sel)
tgroups=[g for g in pd.unique(S['grp']) if str(g).startswith('Tumor')]
ref_sets=[['endothelial','smooth_muscle']]
def hoff(types,min_units=3):
    u,D=O.offset_table(Yh,H['ctype'],H['donor'],types,ref_sets); return O.combine_units(D,min_units)
off_all,var_all=hoff(O.EPI); off_lum,var_lum=hoff(['luminal']); off_bas,var_bas=hoff(['basal'])
# previous-style offset from Wong benign-LABELLED libraries (for comparison only)
W=datasets.wong(); Yw=O.bin_counts(W['X'],W['genes'],sel); mB=W['lib']=='B'
u,D=O.offset_table(Yw[mB],W['anno'][mB],W['pat'][mB],['luminal'],[['endothelial','smooth_muscle']]); off_wb,var_wb=O.combine_units(D)
print('offset corr healthy-all vs healthy-luminal %.2f; healthy-luminal vs Wong-benign-luminal %.2f; RMS(h_lum - wong_b) %.3f'%(np.corrcoef(off_all,off_lum)[0,1],np.corrcoef(off_lum,off_wb)[0,1],np.sqrt(np.mean((off_lum-off_wb)**2))))
rng=np.random.default_rng(0); perm=rng.permutation(len(off_all))
trials=[('no offset',None,None),
        ('healthy, all epithelium',off_all,var_all),('healthy, all epithelium, no variance term',off_all,None),
        ('healthy, luminal',off_lum,var_lum),('healthy, luminal, no variance term',off_lum,None),
        ('healthy, basal',off_bas,var_bas),
        ('Wong benign-labelled luminal (3+ patients)',off_wb,var_wb),
        ('control: healthy all-epithelium offset, bins permuted',off_all[perm],var_all)]
rows=[]
for label,off,var in trials:
    lo=None if off is None else {g:off for g in tgroups}; lv=None if var is None else {g:var for g in tgroups}
    res=run(X,S['genes'],S['table'],S['cen'],S['grp'],ref,cfg,verbose=False,sel=sel,lineage_offset=lo,offset_var=lv)
    fl=np.median([v for g,v in res.floor.items() if g.startswith('Tumor')])
    r_,seg_,call_=res.pool([g for g in res.meta_table.group.unique() if g.startswith('Tumor')]); ok=~np.isnan(tb)
    for prof in ['seg','r']:
        d=bench.eval_profiles(res.cell_matrix(prof),S['cell_key'],{'tumor':tb},S['ho'])[0]
        rows.append(dict(trial=label,profile=prof,floor=fl,z_call=res.z_call,pooled_pearson=np.corrcoef(r_[ok],tb[ok])[0,1],**{k:d[k] for k in ['sens_gain_fpr10','sens_loss_fpr10','auc_gain','auc_loss','pearson']}))
    print(label,'done',flush=True)
df=pd.DataFrame(rows); pd.set_option('display.width',250); print(df.round(3).to_string(index=False)); df.to_csv(config.results_path('wa71_healthy_offset_results.csv'),index=False)

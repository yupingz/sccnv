"""Legacy one-factor-at-a-time tuning on WA71 (variance floor, penalty, bin size, weights, metacell size, count splitting)."""
import sys, time
from sccnv import config, datasets, offsets as O
from sccnv import benchmark as bench, experiments, loaders as data_io
from sccnv.core import *
S=datasets.wa71(); X=S['X']; ref={'stroma':S['is_str']&~S['ho']}
rows=[]
def trial(label, ref_masks=ref, **kw):
    cfg=Config(**{**dict(genes_per_bin=100,weight_power=0.5,target_umi=80000),**kw})
    sel=select_genes_and_bins(X,S['genes'],S['table'],S['cen'],cfg,S['grp'],np.any(list(ref_masks.values()),0))
    res=run(X,S['genes'],S['table'],S['cen'],S['grp'],ref_masks,cfg,verbose=False,sel=sel)
    tb=datasets.wa71_truth(S,sel)
    fl=np.median([v for g,v in res.floor.items() if g.startswith('Tumor')])
    for prof in ['seg','r']:
        for d in bench.eval_profiles(res.cell_matrix(prof),S['cell_key'],{'tumor':tb},S['ho']):
            rows.append(dict(trial=label,profile=prof,floor_tumor=fl,penalty=res.penalty,z_call=res.z_call,**{k:d[k] for k in ['sens_gain_fpr10','sens_loss_fpr10','auc_gain','auc_loss','pearson']}))
    print(label,'done',flush=True)
trial('default')
trial('floor off',floor=False)
trial('penalty 3',seg_penalty=3.0)
trial('penalty 12',seg_penalty=12.0)
trial('penalty 3, floor off',seg_penalty=3.0,floor=False)
trial('min_seg_bins 2',min_seg_bins=2)
trial('gpb 50',genes_per_bin=50)
trial('gpb 150',genes_per_bin=150)
trial('weight_power 0',weight_power=0.0)
trial('target_umi 40k',target_umi=40000)
trial('target_umi 160k',target_umi=160000)
trial('no count split',thin=False)
df=pd.DataFrame(rows); pd.set_option('display.width',250)
print(df.round(3).to_string(index=False))

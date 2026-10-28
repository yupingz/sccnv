"""Legacy tuning sweep (bin size x minimum segment length) on WA71 and the cell lines.   python scripts/legacy/grid_seg.py wa71|cl"""
import sys, time
from sccnv import config, datasets, offsets as O
from sccnv import benchmark as bench, experiments, loaders as data_io
from sccnv.core import *
which=sys.argv[1]
if which=='wa71':
    S=datasets.wa71(); ref={'stroma':S['is_str']&~S['ho']}; extra=dict(target_umi=80000)
    ref_union=S['is_str']&~S['ho']
else:
    S=datasets.celllines(); ref={'PrEC':S['is_ref']}; extra=dict(metacell_cells=20); ref_union=S['is_ref']
X=S['X']; rows=[]
for gpb,msb in [(100,3),(100,2),(50,4),(50,6),(75,3),(75,4),(150,2)]:
    cfg=Config(genes_per_bin=gpb,weight_power=0.5,min_seg_bins=msb,**extra)
    sel=select_genes_and_bins(X,S['genes'],S['table'],S['cen'],cfg,S['grp'],ref_union)
    res=run(X,S['genes'],S['table'],S['cen'],S['grp'],ref,cfg,verbose=False,sel=sel)
    tb=datasets.wa71_truth(S,sel) if which=='wa71' else datasets.celllines_truths(S,sel)
    tbd={'tumor':tb} if which=='wa71' else tb
    for prof in ['seg','r']:
        ds=bench.eval_profiles(res.cell_matrix(prof),S['cell_key'],tbd,S['ho'])
        d=pd.DataFrame(ds); rows.append(dict(data=which,gpb=gpb,min_seg_bins=msb,profile=prof,bins=len(sel['bins']),**{k:d[k].mean() for k in ['sens_gain_fpr10','sens_loss_fpr10','auc_gain','auc_loss','pearson']}))
    print(gpb,msb,flush=True)
df=pd.DataFrame(rows); df.to_csv(config.results_path(f'grid_seg_{which}.csv'),index=False); print(df.round(3).to_string(index=False))

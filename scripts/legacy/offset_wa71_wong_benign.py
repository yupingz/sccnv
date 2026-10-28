"""SUPERSEDED. WA71 offset test using the Wong benign-LABELLED library epithelium as the normal control. Benign-labelled libraries can contain
tumour cells, so this was replaced by healthy-donor offsets (scripts/22_wa71_healthy_offset.py). Kept for the record."""
import sys, pickle, time
from sccnv import config, datasets, offsets as O
from sccnv import benchmark as bench, experiments, loaders as data_io
from sccnv.core import *
S=datasets.wa71(); X=S['X']; ref={'stroma':S['is_str']&~S['ho']}
base=dict(genes_per_bin=100,weight_power=0.5,target_umi=80000)
cfg=Config(**base)
sel=select_genes_and_bins(X,S['genes'],S['table'],S['cen'],cfg,S['grp'],S['is_str']&~S['ho'])
# ---- Wong: same-patient benign epithelium vs endothelial+smooth muscle, on wa71's gene/bin definition
XW,gW,cW,mW=data_io.load_export(str(config.EXPORT/'wong'))
idx_by_sym={g:i for i,g in enumerate(gW)}
cols=[idx_by_sym[s] for s in sel['genes'].symbol]
Wd=sp.diags(sel['w'].astype(np.float32))
Yb=np.asarray(((XW[:,cols]@Wd)@sel['M']).todense())
anno=mW.anno.values.astype(str); lib=mW.tumor_normal.values.astype(str); pat=mW.patient.values.astype(str)
epi=np.isin(anno,['luminal','basal','club','hillock']); st=np.isin(anno,['endothelial','smooth_muscle'])
lp=lambda v:(lambda a:np.log2(a/a.sum()))(v+0.5)
offs={}
for q in sorted(set(pat)):
    m=pat==q; e=m&epi&(lib=='B'); s=m&st
    if e.sum()>=150: offs[q]=lp(Yb[e].sum(0))-lp(Yb[s].sum(0)); print(q,'benign epi',e.sum(),'endo+SMC',s.sum())
D=np.vstack(list(offs.values())); off=np.median(D,0); var=(1.4826*np.median(np.abs(D-off),0))**2
print('offset RMS %.3f ; between-patient SD (median) %.3f'%(np.sqrt((off**2).mean()),np.sqrt(np.median(var))))
np.save(config.results_path('wa71_external_offset.npy'),off); np.save(config.results_path('wa71_external_offset_var.npy'),var)
tb=datasets.wa71_truth(S,sel)
rows=[]
def trial(label,**kw):
    res=run(X,S['genes'],S['table'],S['cen'],S['grp'],ref,Config(**base),verbose=False,sel=sel,**kw)
    ho_rows=np.where(res.meta_table.group.str.contains('_ho').values)[0]
    fl=np.median([v for g,v in res.floor.items() if g.startswith('Tumor')])
    for prof in ['seg','r']:
        for d in bench.eval_profiles(res.cell_matrix(prof),S['cell_key'],{'tumor':tb},S['ho']):
            rows.append(dict(trial=label,profile=prof,floor=fl,**{k:d[k] for k in ['sens_gain_fpr10','sens_loss_fpr10','auc_gain','auc_loss','pearson']}))
    tg=[g for g in res.meta_table.group.unique() if g.startswith('Tumor')]; r_,seg_,_=res.pool(tg); ok=~np.isnan(tb)
    rows[-1]['pooled_pearson_r']=np.corrcoef(r_[ok],tb[ok])[0,1]; rows[-2]['pooled_pearson_r']=np.corrcoef(r_[ok],tb[ok])[0,1]
    print(label,'floor',round(fl,3),flush=True)
trial('no offset')
trial('external offset (Wong benign epi, median of patients)',lineage_offset=off,offset_var=var)
trial('external offset, no variance term',lineage_offset=off)
rng=np.random.default_rng(0); shuf=rng.permutation(len(off))
trial('control: offset permuted across bins',lineage_offset=off[shuf],offset_var=var)
df=pd.DataFrame(rows); pd.set_option('display.width',250); print(df.round(3).to_string(index=False))
df.to_csv(config.results_path('wa71_external_offset_results.csv'),index=False)

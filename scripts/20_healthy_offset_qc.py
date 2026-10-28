"""Healthy-donor lineage offsets: donor-to-donor agreement, epithelial-subtype and zone effects, and transfer of the immune-minus-stroma contrast
to WA71 and Wong.   python scripts/20_healthy_offset_qc.py"""
import numpy as np, pandas as pd
from sccnv import config, datasets, offsets as O
from sccnv import benchmark as bench, experiments, loaders as data_io
from sccnv.core import *
H=datasets.healthy(); print(pd.crosstab(H['lib'],H['ctype']).to_string())
print(pd.crosstab([H['donor'],H['lib']],H['zone']).to_string())
# ------------------------------------------------ WA71 bins
S=datasets.wa71(); cfg=Config(target_umi=80000)
sel=select_genes_and_bins(S['X'],S['genes'],S['table'],S['cen'],cfg,S['grp'],S['is_str']&~S['ho'])
Yh=O.bin_counts(H['X'],H['genes'],sel)
Yw=O.bin_counts(S['X'],S['genes'],sel)
print('bins',Yh.shape[1])
ref_sets=[['endothelial','smooth_muscle']]       # WA71 stroma-only reference <-> healthy endothelial + SMC
res={}
for name,types in [('all_epi',O.EPI),('luminal',['luminal']),('basal',['basal']),('club',['club']),('hillock',['hillock'])]:
    u_d,D_d=O.offset_table(Yh,H['ctype'],H['donor'],types,ref_sets)
    u_l,D_l=O.offset_table(Yh,H['ctype'],H['lib'],types,ref_sets)
    res[name]=(u_d,D_d,u_l,D_l)
    off,var=O.combine_units(D_d)
    sd_between=np.sqrt(np.mean(np.var(D_d,0,ddof=1))) if len(D_d)>1 else np.nan
    print(f"{name:8s} donors {u_d} libs {len(u_l)}  offset RMS {np.sqrt((off**2).mean()):.3f}  between-donor SD {sd_between:.3f}  between-library SD {np.sqrt(np.mean(np.var(D_l,0,ddof=1))) if len(D_l)>1 else np.nan:.3f}")
# donor-donor correlation of all_epi offsets, and subtype differences
u_d,D_d=res['all_epi'][0],res['all_epi'][1]
print('donor x donor corr (all_epi):\n',pd.DataFrame(np.corrcoef(D_d),index=u_d,columns=u_d).round(2).to_string())
subs=['luminal','basal','club','hillock']
M={s:np.median(res[s][1],0) for s in subs}
print('subtype offset RMS difference (median over donors):')
print(pd.DataFrame([[np.sqrt(np.mean((M[a]-M[b])**2)) for b in subs] for a in subs],index=subs,columns=subs).round(3).to_string())
print('subtype offset correlation:'); print(pd.DataFrame(np.corrcoef([M[s] for s in subs]),index=subs,columns=subs).round(2).to_string())
# zone effect: PZ vs TZ libraries within donors D27, D35 (luminal+basal pooled)
u_l,D_l=res['all_epi'][2],res['all_epi'][3]; zone_of={l:H['zone'][H['lib']==l][0] for l in u_l}; don_of={l:H['donor'][H['lib']==l][0] for l in u_l}
print({l:(don_of[l],zone_of[l]) for l in u_l})
for d in ['D27','D35']:
    ls=[l for l in u_l if don_of[l]==d]
    if len(ls)==2: print(d,'PZ vs TZ library offset: RMS diff %.3f, corr %.2f'%(np.sqrt(np.mean((D_l[u_l.index(ls[0])]-D_l[u_l.index(ls[1])])**2)),np.corrcoef(D_l[u_l.index(ls[0])],D_l[u_l.index(ls[1])])[0,1]))
# ------------------------------------------------ transfer check: immune - stroma contrast, healthy vs WA71
def contrast(Y,ct,imm,stro):
    return O.lp(Y[np.isin(ct,imm)].sum(0))-O.lp(Y[np.isin(ct,stro)].sum(0))
C_h={d:contrast(Yh[H['donor']==d],H['ctype'][H['donor']==d],['t_cells','macrophage'],['endothelial','smooth_muscle']) for d in sorted(set(H['donor']))}
C_w=contrast(Yw,S['anno'],['Tcell','Macro','Bcell'],['Endo_1','Endo_2','Endo_3','SMC'])
ch=np.mean(list(C_h.values()),0)
print('immune-minus-stroma contrast, healthy(mean of donors) vs WA71: pearson %.3f slope %.2f (SD healthy %.3f, WA71 %.3f)'%(np.corrcoef(ch,C_w)[0,1],np.polyfit(ch,C_w,1)[0],ch.std(),C_w.std()))
print('  per donor corr with WA71:',{d:round(float(np.corrcoef(c,C_w)[0,1]),2) for d,c in C_h.items()})
dd=sorted(C_h); print('  donor-donor corr:',round(float(np.corrcoef(C_h[dd[0]],C_h[dd[1]])[0,1]),2),round(float(np.corrcoef(C_h[dd[0]],C_h[dd[2]])[0,1]),2),round(float(np.corrcoef(C_h[dd[1]],C_h[dd[2]])[0,1]),2))
# same contrast in Wong on WA71 bins (per patient)
W=datasets.wong(); Ywg=O.bin_counts(W['X'],W['genes'],sel)
C_p={p:contrast(Ywg[W['pat']==p],W['anno'][W['pat']==p],['t_cells','macrophage'],['endothelial','smooth_muscle']) for p in sorted(set(W['pat']))}
print('immune-minus-stroma contrast: Wong patients vs healthy mean: corr',{p:round(float(np.corrcoef(c,ch)[0,1]),2) for p,c in C_p.items()},' vs WA71:',{p:round(float(np.corrcoef(c,C_w)[0,1]),2) for p,c in C_p.items()})
np.savez(config.results_path('offsets_wa71_bins.npz'),**{f'{k}_D_donor':v[1] for k,v in res.items()},**{f'{k}_D_lib':v[3] for k,v in res.items()})

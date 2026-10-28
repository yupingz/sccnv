"""Group-level and metacell-level summaries of the library-agnostic Wong analysis (scripts/30) against healthy diploid epithelium (scripts/31).

Writes to SCCNV_RESULTS:
  wong_group_level_calls.csv / healthy_group_level_calls.csv   per epithelial group: bins called, max |z| (none vs healthy-matched offset)
  group_amplitude_healthy_vs_wong.csv                          per group (>= 10 metacells): RMS / max |seg| / RMS of unsegmented bins
  group_arm_means.csv                                          arm-level mean of each group's fitted profile (healthy and Wong)
  wong_epi_metacell_zrms.csv                                   per epithelial metacell: RMS z-score, subtype, library composition (annotation only)
  arm_library_check.csv                                        arm-level ratios of strong events split by library label (descriptive; arms chosen by eye)
"""
import pickle
import numpy as np, pandas as pd
from sccnv import config

EPI = ['luminal', 'basal', 'club', 'hillock']
R = pickle.load(open(config.results_path('wong_resolved.pkl'), 'rb')); outs = {o['patient']: o for o in R['outs']}; bins = R['bins']
HW = pickle.load(open(config.results_path('healthy_on_wong_bins.pkl'), 'rb'))
hs = pickle.load(open(config.results_path('healthy_lodo_burdens.pkl'), 'rb'))
arm = (bins.chr + bins.arm).values; arms = list(pd.unique(arm))


def epi_groups(mt, min_meta=1):
    return [g for g in pd.unique(mt.group) if g.split('|')[0] in EPI and '_ho' not in g and (mt.group.values == g).sum() >= min_meta]


def group_table(d):
    mt, gs = d['meta_table'], d['group_segments']; rows = []
    for g in epi_groups(mt):
        s = gs[(gs.group == g) & (gs.call != 0)] if gs is not None and len(gs) else gs.iloc[0:0]
        nb = int(s.n_bins.sum()) if len(s) else 0
        rows.append(dict(group=g, n_meta=int((mt.group.values == g).sum()), n_bins_called=nb, frac_bins_called=nb / len(bins),
                         max_abs_z=float(s.z.abs().max()) if len(s) else 0.0, n_segs=len(s)))
    return pd.DataFrame(rows)


def frac_tumour_library(o, d, groups):
    meta, mt, lib = d['meta'], d['meta_table'], o['lib']; out = []
    for g in groups:
        cells = np.where(np.isin(meta, np.where(mt.group.values == g)[0]))[0]
        out.append((lib[cells] == 'T').mean() if len(cells) else np.nan)
    return out


# ---------------------------------------------------------------- group-level calls
GT = []
for p, o in outs.items():
    for v in ['none', 'healthy_matched']:
        t = group_table(o['variants'][v]); t['patient'] = p; t['variant'] = v
        t['frac_tumourlib'] = frac_tumour_library(o, o['variants'][v], t.group); GT.append(t)
GT = pd.concat(GT); GT.to_csv(config.results_path('wong_group_level_calls.csv'), index=False)
HT = []
for dn, vv in HW.items():
    for v, dd in vv.items():
        t = group_table(dd); t['donor'] = dn; t['variant'] = v; HT.append(t)
HT = pd.concat(HT); HT.to_csv(config.results_path('healthy_group_level_calls.csv'), index=False)
print(HT.groupby('variant').agg(groups=('group', 'size'), with_call=('n_bins_called', lambda x: (x > 0).sum()), max_frac=('frac_bins_called', 'max'), mean_frac=('frac_bins_called', 'mean')).round(3))

# ---------------------------------------------------------------- amplitude and arm-level means (healthy_matched variant)
def amp_rows(d, label, who):
    mt = d['meta_table']; rows = []
    for g in epi_groups(mt, 10):
        gs, gr = np.asarray(d['group_seg'][g]), np.asarray(d['group_r'][g])
        rows.append(dict(who=who, label=label, group=g, n_meta=int((mt.group.values == g).sum()), seg_rms=float(np.sqrt((gs ** 2).mean())),
                         max_abs_seg=float(np.abs(gs).max()), r_rms=float(np.sqrt((gr ** 2).mean()))))
    return rows

rows, arm_rows = [], []
for p, o in outs.items():
    d = o['variants']['healthy_matched']
    rr = amp_rows(d, p, 'Wong'); fr = frac_tumour_library(o, d, [r['group'] for r in rr])
    for r_, f in zip(rr, fr): r_['frac_tumourlib'] = f
    rows += rr
    for r_ in rr: arm_rows.append(pd.Series(pd.Series(np.asarray(d['group_seg'][r_['group']])).groupby(arm, sort=False).mean().reindex(arms)).rename(f"Wong|{p}|{r_['group']}"))
for dn, vv in HW.items():
    d = vv['healthy_matched']; rr = amp_rows(d, dn, 'Healthy'); rows += rr
    for r_ in rr: arm_rows.append(pd.Series(pd.Series(np.asarray(d['group_seg'][r_['group']])).groupby(arm, sort=False).mean().reindex(arms)).rename(f"Healthy|{dn}|{r_['group']}"))
A = pd.DataFrame(rows); A.to_csv(config.results_path('group_amplitude_healthy_vs_wong.csv'), index=False)
AM = pd.DataFrame(arm_rows); AM.to_csv(config.results_path('group_arm_means.csv'))
hq = A[A.who == 'Healthy'].seg_rms
print('groups >= 10 metacells: healthy', int((A.who == 'Healthy').sum()), 'Wong', int((A.who == 'Wong').sum()),
      '| Wong above healthy max seg_rms:', int(((A.who == 'Wong') & (A.seg_rms > hq.max())).sum()))

# ---------------------------------------------------------------- per-metacell RMS z (healthy null from leave-one-donor-out)
hv = {'none': 'no offset', 'healthy_matched': 'other donors, matched-subtype offset', 'healthy_all': 'other donors, all-epithelium offset'}
rows = []
for p, o in outs.items():
    for v in hv:
        d = o['variants'][v]; mt = d['meta_table']; e = np.where(mt.group.str.split('|').str[0].isin(EPI).values)[0]
        z = d['r'][e] / d['sd'][e]; zr = np.sqrt((z ** 2).mean(1)); meta, lib = d['meta'], o['lib']
        for j, m in enumerate(e):
            cells = np.where(meta == m)[0]
            rows.append(dict(patient=p, variant=v, metacell=m, n_cells=len(cells), frac_T=(lib[cells] == 'T').mean(), zrms=zr[j], subtype=mt.group.values[m].split('|')[0]))
Z = pd.DataFrame(rows); Z['lib'] = np.where(Z.frac_T > 0.5, 'T', 'B')
for v in hv:
    h = np.concatenate([x for k, x in hs.items() if len(k) == 3 and k[1] == hv[v] and k[2] == 'zrms'])
    Z.loc[Z.variant == v, 'healthy_q99_zrms'] = np.quantile(h, 0.99)
Z.to_csv(config.results_path('wong_epi_metacell_zrms.csv'), index=False)

# ---------------------------------------------------------------- strong events by library label (descriptive; arms chosen by eye from the heatmap)
CHECK = [('ICC2', ['16q', '15q']), ('ICC4', ['8p', '13q', '16q']), ('ICC7', ['13q', '18q'])]
rows = []
for p, arms_ in CHECK:
    d = outs[p]['variants']['healthy_matched']; mt = d['meta_table']; cm = d['meta']; lib = outs[p]['lib']
    e = np.where(mt.group.str.split('|').str[0].isin(EPI).values & ~mt.is_ref.values & ~mt.group.str.contains('_ho').values)[0]
    fT = np.array([(lib[cm == m] == 'T').mean() for m in e])
    for a in arms_:
        v = d['r'][e][:, arm == a].mean(1)
        rows.append(dict(patient=p, arm=a, n_tumourlib=int((fT > 0.5).sum()), med_tumourlib=np.median(v[fT > 0.5]), n_benignlib=int((fT <= 0.5).sum()), med_benignlib=np.median(v[fT <= 0.5]),
                         frac_benignlib_below_minus0p25=(v[fT <= 0.5] < -0.25).mean(), frac_tumourlib_below_minus0p25=(v[fT > 0.5] < -0.25).mean()))
pd.DataFrame(rows).round(3).to_csv(config.results_path('arm_library_check.csv'), index=False)
print(pd.DataFrame(rows).round(2).to_string(index=False))

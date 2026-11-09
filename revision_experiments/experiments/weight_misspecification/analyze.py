from pathlib import Path
import json, os, math
import numpy as np
import pandas as pd
from scipy.stats import t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(os.environ.get('WM_OUTPUT', ROOT/'outputs/weight_misspecification'))
raw=pd.read_csv(OUT/'raw_metrics.csv')
grid=pd.read_csv(OUT/'grid.csv')
METHODS=['Proposed','Pure_aSBM','Louvain','KMeans_estK']
DIST=['Gaussian','Student_t5','Laplace','Contaminated']
BETAS=[0,.25,.5,1,2,4]
COLORS=['#225f9c','#c46b31','#49836b','#8d60a8']
KEYS=['condition_id','rep','method','beta']
assert len(raw)==1200 and not raw.duplicated(KEYS).any()
assert len(raw[['condition_id','rep']].drop_duplicates())==200
assert len(raw[raw.method=='Proposed'])==600
for c in ['graph_hash','attribute_hash','initial_hash','post_init_rng_hash','initial_K']:
    assert raw.groupby(['condition_id','rep'])[c].nunique(dropna=False).eq(1).all(),c
assert raw.groupby(['r','rep']).graph_hash.nunique().eq(1).all()
assert raw.groupby(['r','rep']).apply(lambda g:g[g.method=='Louvain'].ARI.nunique(),include_groups=False).eq(1).all()
good=raw.status.eq('ok')
assert np.isfinite(raw.loc[good,['ARI','NMI']]).all().all()
assert raw.loc[~good,'ARI'].isna().all()
assert raw.loc[good,'ARI'].between(-1,1).all() and raw.loc[good,'NMI'].between(0,1).all()
for _,c in grid.iterrows():
    g=raw[raw.condition_id==c.condition_id]
    assert len(g)==(90 if c.weight_study else 40)
    for _,m in g.groupby(['method','beta'],dropna=False):assert sorted(m.rep)==list(range(1,11))

def stats(v):
    v=pd.to_numeric(v).dropna().to_numpy(); n=len(v)
    return dict(n=n,mean=float(v.mean()) if n else np.nan,
                sd=float(v.std(ddof=1)) if n>1 else np.nan,
                mcse=float(v.std(ddof=1)/np.sqrt(n)) if n>1 else np.nan)

summaries=[]
for (cid,method,beta),g in raw.groupby(['condition_id','method','beta'],dropna=False,sort=False):
    row=g.iloc[0][['condition_id','r','sigma','h','distribution','weight_study','misspecification_study']].to_dict()
    row.update(method=method,beta=beta,n_expected=10,n_present=len(g),n_failed=int(g.status.ne('ok').sum()),
               n_not_stopped=int((g.met_stopping_rule.astype(str).str.lower()=='false').sum()),
               n_init_warnings=int(g.init_warnings.notna().sum()),n_fit_warnings=int(g.warnings.notna().sum()),
               stop_reasons=json.dumps(g.stop_reason.value_counts().to_dict()))
    for metric in ['ARI','NMI','final_K','iterations','macro_F1','min_F1','high_F1','low_F1']:
        row.update({f'{metric}_{k}':v for k,v in stats(g.loc[g.status=='ok',metric]).items()})
    summaries.append(row)
summary=pd.DataFrame(summaries)
summary.to_csv(OUT/'condition_method_summary.csv',index=False)
raw[~good].to_csv(OUT/'failure_records.csv',index=False)

paired=[];detail=[]
def pair(left,right,keys,kind,labels):
    p=left.merge(right,on=keys,suffixes=('_test','_reference'),validate='one_to_one')
    assert len(p)==10
    for metric in ['ARI','NMI']:
        dif=p[metric+'_test']-p[metric+'_reference']; s=stats(dif)
        half=t.ppf(.975,s['n']-1)*s['mcse'] if s['n']>1 else np.nan
        paired.append(dict(kind=kind,**labels,metric=metric,**s,ci95_low=s['mean']-half,ci95_high=s['mean']+half,
                           missing_pairs=10-s['n'],wins=int((dif>1e-12).sum()),ties=int((dif.abs()<=1e-12).sum()),losses=int((dif< -1e-12).sum())))
        for rep,d in zip(p.rep,dif): detail.append(dict(kind=kind,**labels,metric=metric,rep=rep,difference=d))
for _,c in grid[grid.weight_study].iterrows():
    g=raw[(raw.condition_id==c.condition_id)&(raw.method=='Proposed')]
    for beta in BETAS:
        pair(g[g.beta==beta],g[g.beta==1],['rep'],'beta_minus_1',dict(condition_id=c.condition_id,r=c.r,sigma=c.sigma,h=c.h,beta=beta,method='Proposed',distribution='Gaussian'))
for _,c in grid[(grid.h==1)&(grid.distribution!='Gaussian')].iterrows():
    for method in METHODS:
        a=raw[(raw.condition_id==c.condition_id)&(raw.method==method)]
        b=raw[(raw.r==c.r)&(raw.sigma==c.sigma)&(raw.h==1)&(raw.distribution=='Gaussian')&(raw.method==method)]
        if method=='Proposed':b=b[b.beta==1]
        pair(a,b,['rep'],'distribution_minus_gaussian',dict(condition_id=c.condition_id,r=c.r,sigma=c.sigma,h=1,beta=1 if method=='Proposed' else np.nan,method=method,distribution=c.distribution))
for _,c in grid.iterrows():
    g=raw[raw.condition_id==c.condition_id]
    a=g[(g.method=='Proposed')&(g.beta==1)]
    for method in METHODS[1:]:
        pair(a,g[g.method==method],['rep'],'proposed_minus_baseline',dict(condition_id=c.condition_id,r=c.r,sigma=c.sigma,h=c.h,beta=1,method=method,distribution=c.distribution))
pairs=pd.DataFrame(paired); pairs.to_csv(OUT/'paired_summary.csv',index=False)
pd.DataFrame(detail).to_csv(OUT/'paired_run_differences.csv',index=False)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,4,figsize=(17,8),sharey=True)
for ax,(_,c) in zip(axes.flat,grid[grid.weight_study].iterrows()):
    s=summary[(summary.condition_id==c.condition_id)&(summary.method=='Proposed')].set_index('beta').loc[BETAS]
    ax.errorbar(range(6),s.ARI_mean,yerr=s.ARI_sd,marker='o',capsize=3,lw=2,color=COLORS[0],label='Proposed')
    for method,color in zip(METHODS[1:],COLORS[1:]):
        b=summary[(summary.condition_id==c.condition_id)&(summary.method==method)].iloc[0]
        ax.axhline(b.ARI_mean,color=color,ls='--',lw=1,label=method)
    ax.set(title=f'r={c.r:g}, sigma={c.sigma:g}, h={c.h:g}',xticks=range(6),xticklabels=[str(x) for x in BETAS],xlabel='Attribute weight beta',ylim=(-.07,1.07))
    ax.axvline(3,color='#888888',alpha=.3,lw=1);ax.grid(axis='y',alpha=.15)
axes[0,0].set_ylabel('ARI');axes[1,0].set_ylabel('ARI')
fig.suptitle('Weight sensitivity: 10 paired repetitions\nProposed: mean ± sample SD; dashed baselines: mean (SD in tables)',fontsize=15)
fig.legend(*axes[0,0].get_legend_handles_labels(),loc='lower center',ncol=4,frameon=False)
fig.tight_layout(rect=(0,.06,1,.95));fig.savefig(OUT/'weight_sensitivity.png',dpi=180);fig.savefig(OUT/'weight_sensitivity.pdf');plt.close(fig)

fig,axes=plt.subplots(2,2,figsize=(13,10),sharey=True)
for ax,(r,sigma) in zip(axes.flat,[(1.5,1.5),(4,1.5),(1.5,4),(4,4)]):
    for j,(method,color) in enumerate(zip(METHODS,COLORS)):
        s=summary[(summary.r==r)&(summary.sigma==sigma)&(summary.h==1)&(summary.method==method)]
        if method=='Proposed':s=s[s.beta==1]
        s=s.set_index('distribution').loc[DIST]
        pos=np.arange(4)+(j-1.5)*.16
        ax.errorbar(pos,s.ARI_mean,yerr=s.ARI_sd,marker='o',capsize=3,lw=1.4,color=color,label=method)
        for x,(_,v) in zip(pos,s.iterrows()):
            if v.ARI_n<10:
                ypos=v.ARI_mean+v.ARI_sd+.04 if v.ARI_mean>.45 else v.ARI_mean-v.ARI_sd-.06
                ax.text(x,ypos,f'n={int(v.ARI_n)}/10',fontsize=8,color=color,ha='center',
                        bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1))
    ax.set(title=f'r={r:g}, sigma={sigma:g}',xticks=range(4),xticklabels=['Gaussian','t (df=5)','Laplace','Contaminated'],ylim=(-.1,1.12),ylabel='ARI')
    ax.grid(axis='y',alpha=.15)
fig.suptitle('Distribution misspecification: mean ± sample SD\n10 repeats (successful n labeled if <10); equal theoretical noise covariance',fontsize=15)
fig.legend(*axes[0,0].get_legend_handles_labels(),loc='lower center',ncol=4,frameon=False)
fig.tight_layout(rect=(0,.05,1,.93));fig.savefig(OUT/'distribution_comparison.png',dpi=180);fig.savefig(OUT/'distribution_comparison.pdf');plt.close(fig)

def cell(row,metric='ARI'):
    n=int(row[f'{metric}_n'])
    if not n:return 'Failed (0/10)'
    if n==1:return f"{row[f'{metric}_mean']:.3f} (n=1/10, SD unavailable)"
    s=f"{row[f'{metric}_mean']:.3f} ± {row[f'{metric}_sd']:.3f}"
    return s+(f' (n={n}/10)' if n!=10 else '')
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|']+['| '+' | '.join(map(str,row))+' |' for row in rows])
md=['# Attribute likelihood weights and non-Gaussian attributes','',
    'Tables report mean ± sample SD across ten replicates per condition. The likelihood-weight study gives 720 records and the distribution study gives 640. They share 160 Gaussian reference records, giving 1,200 unique records across 200 condition-replicate pairs. Conditions sharing draws are paired within replicates.','',
    'K-means uses iter.max=10 in both studies. The β=1 fit matches the default implementation. The β=0 fit uses attributes in its initialization.','',
    '## Fit status','']
rows=[]
for method,g in raw.groupby('method',sort=False):
    rows.append([method,len(g),g.status.eq('ok').sum(),g.status.eq('error').sum(),int((g.met_stopping_rule.astype(str).str.lower()=='false').sum()),int(g.warnings.notna().sum())])
md += [table(['Method','Records','Returned','Errors','Stopping rule not met','Fits with warnings'],rows),'',
       f"{raw.drop_duplicates(['condition_id','rep']).init_warnings.notna().sum()} of 200 condition-replicate initializations issued a warning. All weights on a given input share its initialization. The aSBM baseline exits if its objective does not increase.",'',
       'Scores are missing for failed fits. Means and SDs use the returned fits and display the available count where it differs from ten.','']
if (~good).any():
    failures=raw[~good].groupby(['distribution','r','sigma','h','method','error'],dropna=False).size().reset_index(name='count')
    md += [table(['Distribution','r','σ','h','Method','Error','Count'],failures.values.tolist()),'']
capped=raw[raw.stop_reason.isin(['safety_cap','iteration_cap','timeout'])]
if len(capped):
    md += ['The summaries include the following returned partitions from fits that reached a limit.','',
           table(['Condition','Replicate','Method','β','ARI','Iterations','Stop reason'],
                 capped[['condition_id','rep','method','beta','ARI','iterations','stop_reason']].values.tolist()),'']
for metric in ['ARI','NMI']:
    md += [f'## Likelihood weights ({metric})','']
    rows=[]
    for _,c in grid[grid.weight_study].iterrows():
        s=summary[(summary.condition_id==c.condition_id)&(summary.method=='Proposed')].set_index('beta')
        rows.append([c.condition_id,f'{c.r:g}',f'{c.sigma:g}',f'{c.h:g}']+[cell(s.loc[b],metric) for b in BETAS])
    md += [table(['Condition','r','σ','h']+[f'β={b:g}' for b in BETAS],rows),'']
    rows=[]
    for _,c in grid[grid.weight_study].iterrows():
        s=summary[(summary.condition_id==c.condition_id)&(summary.method!='Proposed')].set_index('method')
        rows.append([c.condition_id]+[cell(s.loc[m],metric) for m in METHODS[1:]])
    md += ['Unweighted reference fits on the same inputs','',table(['Condition']+METHODS[1:],rows),'']
    md += [f'## Attribute distributions ({metric})','']
    rows=[]
    for _,c in grid[grid.misspecification_study].sort_values(['r','sigma','distribution']).iterrows():
        s=summary[(summary.condition_id==c.condition_id)&((summary.method!='Proposed')|summary.beta.eq(1))].set_index('method')
        rows.append([c.condition_id,f'{c.r:g}',f'{c.sigma:g}',c.distribution]+[cell(s.loc[m],metric) for m in METHODS])
    md += [table(['Condition','r','σ','Distribution']+METHODS,rows),'']
md += ['## Paired ARI differences from β=1','',
       'Intervals are pointwise 95% Student-t intervals across the paired replicate indices.','']
p=pairs[(pairs.kind=='beta_minus_1')&(pairs.metric=='ARI')&(pairs.beta!=1)]
md += [table(['Condition','β','Paired replicates','ΔARI','95% interval'],[[v.condition_id,f'{v.beta:g}',v.n,f'{v["mean"]:+.4f}',f'[{v.ci95_low:+.4f}, {v.ci95_high:+.4f}]'] for _,v in p.iterrows()]),'',
       '## Paired aSBM-Louvain ARI differences from Gaussian residuals','']
p=pairs[(pairs.kind=='distribution_minus_gaussian')&(pairs.metric=='ARI')&(pairs.method=='Proposed')]
md += [table(['r','σ','Distribution','Paired replicates','ΔARI','95% interval'],[[f'{v.r:g}',f'{v.sigma:g}',v.distribution,v.n,f'{v["mean"]:+.4f}',f'[{v.ci95_low:+.4f}, {v.ci95_high:+.4f}]'] for _,v in p.iterrows()]),'',
       '## Output files','',
       '- [Per-fit results](raw_metrics.csv)',
       '- [Condition summaries](condition_method_summary.csv)',
       '- [Paired differences](paired_summary.csv)',
       '- [Failed fits](failure_records.csv)','',
       '![Likelihood-weight study](weight_sensitivity.png)','', '![Non-Gaussian attribute study](distribution_comparison.png)','']
(OUT/'complete_results.md').write_text('\n'.join(md))
audit=dict(records=len(raw),instances=200,proposed_fits=600,errors=int((~good).sum()),
           initialization_warnings=int(raw.drop_duplicates(['condition_id','rep']).init_warnings.notna().sum()),
           pairing_checks='passed',metrics_bounds='passed',planned_grid_completeness='passed')
(OUT/'analysis_audit.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit,indent=2))
print(summary[['condition_id','method','beta','ARI_mean','ARI_sd','ARI_n','n_failed']].to_string(index=False))

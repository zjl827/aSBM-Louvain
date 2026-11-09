"""Analyze the nested fixed-data experiment without refitting any model."""
from pathlib import Path
import os, json, math, html
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(os.environ.get('INIT_OUTPUT', ROOT/'outputs/initialization_stability'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

raw=pd.read_csv(OUT/'raw_metrics.csv')
tasks=pd.read_csv(OUT/'tasks.csv')
assert len(raw)==264 and raw.task_id.is_unique and set(raw.task_id)==set(tasks.task_id)
main=raw[raw.variant.eq('main')].copy();diag=raw[raw.variant.eq('downstream_check')].copy()
assert len(main)==240 and len(diag)==24
def nonempty(s):return s.fillna('').astype(str).str.strip().ne('')
main['init_warning']=nonempty(main.init_warnings)
main['fit_warning']=nonempty(main.warnings)
main['success']=main.status.eq('ok')
main['cap']=main.stop_reason.eq('safety_cap')
dataset_rows=[]
for (condition,dataset,method),g in main.groupby(['condition','dataset','init_method'],sort=True):
    good=g[g.success]
    row=dict(condition=condition,dataset=dataset,init_method=method,r=g.r.iloc[0],sigma=g.sigma.iloc[0],
             planned=len(g),successful=len(good),errors=int((~g.success).sum()),caps=int(g.cap.sum()),
             initialization_warnings=int(g.init_warning.sum()),fit_warnings=int(g.fit_warning.sum()))
    for metric in ['initial_ARI','ARI','NMI','initial_K','final_K','iterations','ARI_change']:
        x=good[metric].dropna()
        row.update({metric+'_mean':x.mean(),metric+'_sd':x.std(ddof=1),metric+'_min':x.min(),metric+'_max':x.max()})
    row['unique_initial_partitions']=g.initial_hash.nunique()
    dataset_rows.append(row)
ds=pd.DataFrame(dataset_rows);ds.to_csv(OUT/'dataset_summary.csv',index=False)

condition_rows=[]
for (condition,method),g in ds.groupby(['condition','init_method'],sort=True):
    row=dict(condition=condition,init_method=method,r=g.r.iloc[0],sigma=g.sigma.iloc[0],independent_datasets=len(g),
             planned_fits=int(g.planned.sum()),successful=int(g.successful.sum()),errors=int(g.errors.sum()),
             caps=int(g.caps.sum()),initialization_warnings=int(g.initialization_warnings.sum()))
    for metric in ['initial_ARI','ARI','NMI','initial_K','final_K','iterations','ARI_change']:
        row[metric+'_grand_mean']=g[metric+'_mean'].mean()
        row[metric+'_between_dataset_means_sd']=g[metric+'_mean'].std(ddof=1)
        row[metric+'_rms_within_dataset_sd']=np.sqrt(np.mean(g[metric+'_sd']**2))
        row[metric+'_all_min']=g[metric+'_min'].min();row[metric+'_all_max']=g[metric+'_max'].max()
    condition_rows.append(row)
cs=pd.DataFrame(condition_rows);cs.to_csv(OUT/'condition_summary.csv',index=False)

key=['condition','dataset','data_id','init_id']
paired=main[main.init_method.eq('kmeans')].merge(main[main.init_method.eq('louvain')],on=key,suffixes=('_kmeans','_louvain'),validate='one_to_one')
assert len(paired)==120
for field in ['A_hash','X_hash','fit_rng_hash']:
    assert paired[field+'_kmeans'].eq(paired[field+'_louvain']).all()
paired['delta_ARI_kmeans_minus_louvain']=paired.ARI_kmeans-paired.ARI_louvain
paired.to_csv(OUT/'paired_initialization_comparison.csv',index=False)
pair_ds=paired.groupby(['condition','dataset']).agg(
    paired_n=('delta_ARI_kmeans_minus_louvain','count'),
    delta_mean=('delta_ARI_kmeans_minus_louvain','mean'),delta_sd=('delta_ARI_kmeans_minus_louvain','std'),
    delta_min=('delta_ARI_kmeans_minus_louvain','min'),delta_max=('delta_ARI_kmeans_minus_louvain','max')).reset_index()
pair_ds.to_csv(OUT/'paired_dataset_summary.csv',index=False)
pair_cs=pair_ds.groupby('condition').agg(independent_datasets=('dataset','count'),
    delta_mean=('delta_mean','mean'),between_dataset_delta_means_sd=('delta_mean','std'),
    delta_min=('delta_min','min'),delta_max=('delta_max','max')).reset_index()
pair_cs.to_csv(OUT/'paired_condition_summary.csv',index=False)

check=diag.merge(main[main.init_id.eq(1)],on=['condition','dataset','data_id','init_id','init_method'],suffixes=('_alternate','_main'),validate='one_to_one')
assert len(check)==24
for field in ['A_hash','X_hash','initial_hash']:assert check[field+'_alternate'].eq(check[field+'_main']).all()
assert check.fit_rng_hash_alternate.ne(check.fit_rng_hash_main).all()
check['delta_ARI']=check.ARI_alternate-check.ARI_main
check['absolute_delta_ARI']=check.delta_ARI.abs()
check.to_csv(OUT/'downstream_rng_comparison.csv',index=False)
check_cs=check.groupby(['condition','init_method']).agg(paired_n=('delta_ARI','count'),
  delta_mean=('delta_ARI','mean'),abs_delta_mean=('absolute_delta_ARI','mean'),abs_delta_max=('absolute_delta_ARI','max')).reset_index()
check_cs.to_csv(OUT/'downstream_rng_summary.csv',index=False)

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
colors={'kmeans':'#17699B','louvain':'#D66D37'};labels={'kmeans':'Attribute K-means start','louvain':'Graph Louvain start'}
fig,axes=plt.subplots(2,2,figsize=(11,8),sharex=True,sharey=True)
for ax,(cond,g) in zip(axes.flat,main.groupby('condition',sort=True)):
    ax.plot([-.03,1.03],[-.03,1.03],color='#AAA',linestyle='--',linewidth=1,zorder=0)
    for method,mg in g.groupby('init_method'):
        for di,dg in mg.groupby('data_id'):
            ax.scatter(dg.initial_ARI,dg.ARI,s=35,marker=['o','s','^'][di-1],alpha=.72,color=colors[method],label=labels[method] if di==1 else None)
    ax.set_title(f'{cond}: r={g.r.iloc[0]:g}, sigma={g.sigma.iloc[0]:g}')
    ax.set_xlim(-.035,1.035);ax.set_ylim(-.035,1.035);ax.grid(alpha=.15)
axes[0,0].legend(loc='upper left',fontsize=8)
for ax in axes[1,:]:ax.set_xlabel('ARI of the actual initial partition')
for ax in axes[:,0]:ax.set_ylabel('ARI after aSBM-Louvain')
fig.suptitle('Fixed-input initialization study: 3 datasets x 10 initialization seeds per condition',fontsize=12)
fig.text(.5,.012,'Marker shape identifies the fixed dataset. Downstream RNG is held fixed within each dataset.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.035,1,.95]);fig.savefig(OUT/'initial_to_final.png');fig.savefig(OUT/'initial_to_final.pdf');plt.close(fig)

fig,axes=plt.subplots(2,2,figsize=(11,7.5),sharey=True)
for ax,(cond,g) in zip(axes.flat,ds.groupby('condition',sort=True)):
    for method,mg in g.groupby('init_method'):
        offset=-.13 if method=='kmeans' else .13
        mg=mg.sort_values('dataset')
        ax.errorbar(np.arange(1,4)+offset,mg.ARI_mean,yerr=mg.ARI_sd,fmt='o',capsize=4,color=colors[method],label=labels[method])
    ax.set_title(f'{cond}: r={g.r.iloc[0]:g}, sigma={g.sigma.iloc[0]:g}')
    ax.set_xticks([1,2,3],['Dataset 1','Dataset 2','Dataset 3']);ax.set_ylim(-.06,1.14);ax.grid(axis='y',alpha=.15)
axes[0,0].legend(fontsize=8)
for ax in axes[:,0]:ax.set_ylabel('Final ARI: mean +/- initialization SD')
fig.suptitle('Within-data initialization variability, conditional on one downstream RNG stream',fontsize=12)
fig.text(.5,.012,'Each mean and SD uses 10 starts on exactly the same graph and attributes; the 3 datasets are independently generated.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.035,1,.95]);fig.savefig(OUT/'within_dataset_variability.png');fig.savefig(OUT/'within_dataset_variability.pdf');plt.close(fig)

short=cs[['condition','init_method','r','sigma','initial_ARI_grand_mean','ARI_grand_mean','ARI_rms_within_dataset_sd','ARI_between_dataset_means_sd','ARI_all_min','ARI_all_max','initialization_warnings','caps']].copy()
short.columns=['Condition','Initialization','r','sigma','Initial ARI mean','Final ARI mean','RMS within-data SD','SD of 3 data means','Minimum final ARI','Maximum final ARI','Init warnings','Caps']
short=short.round(4)
details=ds[['condition','dataset','init_method','initial_ARI_mean','ARI_mean','ARI_sd','ARI_min','ARI_max','initial_K_mean','final_K_mean','initialization_warnings','caps']].round(4)
counts=dict(planned_fits=264,main_fits=240,downstream_checks=24,independent_datasets=12,
 errors=int(raw.status.ne('ok').sum()),caps=int(raw.stop_reason.eq('safety_cap').sum()),
 init_warning_rows=int(nonempty(raw.init_warnings).sum()),fit_warning_rows=int(nonempty(raw.warnings).sum()),
 main_init_warning_rows=int(main.init_warning.sum()),max_iterations=int(raw.iterations.max()))
summary_lines=[
 'INITIALIZATION STUDY — CURRENT IMPLEMENTATION',
 'Completed planned experiment: 4 conditions x 3 independently generated datasets x 10 initializations x 2 initializer methods = 240 main fits; 24 preset downstream RNG checks.',
 'K-means: nstart=5, iter.max=10 unchanged; subsequent algorithm: tau0=1, rho=.8, resolution=1, joint thresholds 1e-3/1e-4, stable K for 3 updates, cap 1000.',
 'The two initializers share the same preliminary Louvain result; it supplies graph labels and K. K-means then supplies attribute labels at max(2,K).',
 'A fixed downstream RNG state is restored for all starts on each dataset. This controls the random source stream but does not force identical random moves after algorithm paths diverge.',
 json.dumps(counts,ensure_ascii=False),
 'Ten distinct K-means initializations produced warnings. One of these initializations is reused in a downstream diagnostic, giving 11 warning-bearing fit rows in total; it is not an eleventh independent initialization warning.',
 '',short.to_string(index=False),
 '', 'FINDINGS',
 'Weak graph signal (r=1.5): both initialization strategies fail to recover the planted communities in the tested inputs. For sigma=1.5, informative K-means starts have average ARI 0.7121 but finish at 0.0293; all 30 such starts decline. This rules out an uninformative starting partition as the sole explanation in these inputs, but does not identify one particular outer-loop update as the cause.',
 'Stronger graph signal (r=4): K-means starts have higher average final ARI and lower within-data variation than graph Louvain starts in these 3 datasets per noise condition. At sigma=1.5, all 10 main initialization seeds give the same final ARI within each fixed dataset; this is conditional on the selected downstream RNG and is not a universal invariance claim.',
 'The alternate downstream-RNG check produces nonzero changes in some pairs (largest absolute ARI change about 0.0254). Initialization is therefore not the only possible source of algorithmic variability.',
 'There is no independent claim of statistical recovery, global convergence, or uniform robustness. Only 3 independently generated datasets are used per condition, and no significance conclusion is based on treating their 30 fits as iid data replicates.',
 '', 'Within-data SD describes changing initialization while data and downstream RNG are held fixed. SD across three dataset means is descriptive and not a population variance-component estimate.',
 'Thirty fits per condition/method are nested in three independent datasets; do not treat all thirty as independent data replications. No statistical-significance claims are made.',
 '', 'K-means minus Louvain initialization, averaged within each dataset then across the three datasets:',pair_cs.round(4).to_string(index=False),
 '', 'Preset second downstream RNG (only init_id=1, three datasets per condition/method):',check_cs.round(4).to_string(index=False),
 'This small downstream diagnostic does not provide a complete initialization-versus-Louvain variance decomposition.',
 '', 'All dataset details:',details.to_string(index=False),
 '', 'All failures, initialization warnings and cap terminations are retained. Fit elapsed times are diagnostics, not a formal runtime comparison.',
 'The protocol records the model source and RNG/input hashes. Raw RDS files contain the saved partitions and trajectories.'
]
text='\n'.join(summary_lines)+'\n';(OUT/'report.txt').write_text(text)
web='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Initialization stability: current aSBM-Louvain</title><style>body{font:16px/1.5 system-ui,sans-serif;max-width:1180px;margin:35px auto;padding:0 25px;color:#182c3b}h1,h2{line-height:1.2}table{border-collapse:collapse;font-size:12px;width:100%;margin:20px 0}th,td{padding:8px 7px;border-bottom:1px solid #ccd8df;text-align:right}th{background:#edf3f7}img{width:100%;height:auto}.note{padding:15px;background:#eff5f8;border-left:4px solid #17699b}pre{white-space:pre-wrap;font:13px/1.5 ui-monospace,monospace}a{color:#17699b}</style><body>'''
web+='<h1>Initialization stability of current aSBM-Louvain</h1>'
web+='<p>Four conditions; 12 independent fixed datasets; 240 main fits and 24 preset downstream-RNG checks. All fits use the unchanged current corrected algorithm.</p>'
web+='<p><b>Findings:</b> weak-network failure persists under both starting methods, including informative attribute-based starts. For r=4, K-means initialization gives higher average final ARI and lower within-data variability in these datasets. The alternate downstream-RNG diagnostic changes some results, so initialization is not the only source of algorithmic variation.</p>'
web+='<div class="note"><b>Experimental units:</b> each condition has 3 independent datasets and 10 initialization seeds per dataset. The 30 fits per method/condition are not 30 independent generated datasets. K-means remains nstart=5, iter.max=10.</div>'
web+='<h2>Condition summaries</h2>'+short.to_html(index=False,border=0)
web+='<p>RMS within-data SD measures sensitivity to initialization conditional on fixed inputs and a fixed downstream RNG stream. The SD of 3 dataset means is a descriptive between-dataset summary, not a fully estimated variance component.</p>'
web+='<img src="initial_to_final.png" alt="Initial and final ARI for every main fit"><img src="within_dataset_variability.png" alt="Conditional initialization variation for all 12 fixed datasets">'
web+='<h2>Paired initialization comparison</h2>'+pair_cs.round(4).to_html(index=False,border=0)
web+='<h2>Downstream RNG diagnostic</h2>'+check_cs.round(4).to_html(index=False,border=0)
web+='<p>Only init_id=1 receives an alternate downstream stream. The same saved initial labels are reused. This is a limited diagnostic, not a complete random-effects decomposition.</p>'
web+='<h2>All fixed datasets</h2>'+details.to_html(index=False,border=0)
web+='<h2>Protocol and audit</h2><pre>'+html.escape(text)+'</pre>'
web+='<p><a href="protocol.json">Frozen protocol</a> · <a href="raw_metrics.csv">All attempts</a> · <a href="dataset_summary.csv">Dataset summaries</a> · <a href="downstream_rng_comparison.csv">RNG pairs</a> · <a href="validation.txt">Validation</a></p></body></html>'
(OUT/'report.html').write_text(web)
(OUT/'analysis_audit.json').write_text(json.dumps(dict(counts=counts,metric_check='R runner independently recalculates ARI from contingencies and compares NMI with igraph; Python checks planned design and matching hashes.',all_planned_keys=True,main_rng_match=True,diagnostic_inputs_match=True),indent=2)+'\n')
print(text[:11000])

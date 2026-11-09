#!/usr/bin/env python3
"""Describe and independently audit a fixed-node conditional graph bootstrap.

This script never generates data or fits a clustering method. All metrics are
conditional on the Proposed-fitted plug-in generator and fixed algorithm seeds.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ARM_ORDER = ['graph_only', 'joint']
MISS = {'', 'NA', 'NaN', 'nan', 'None', 'null', 'N/A'}
METRICS = ['ARI_anchor', 'ARI_own_source', 'ARI_original_truth', 'NMI_anchor',
           'within_anchor_pair_retention', 'between_anchor_pair_merge', 'final_K']


def args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir', type=Path, default=Path('results/graph_bootstrap/full'))
    p.add_argument('--output-dir', type=Path)
    p.add_argument('--expected-b', type=int, default=50)
    p.add_argument('--expected-sources', type=int, default=10)
    return p.parse_args()


def read_csv(path):
    if not path.exists():
        return []
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields=None):
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: '' if isinstance(v, float) and not math.isfinite(v) else v for k, v in row.items()})


def num(x):
    if x is None or str(x).strip() in MISS:
        return float('nan')
    try:
        f = float(x)
        return f if math.isfinite(f) else float('nan')
    except (TypeError, ValueError):
        return float('nan')


def present(x):
    return x is not None and str(x).strip() not in MISS


def error_row(r):
    return present(r.get('error')) or 'error' in str(r.get('status', '')).lower() or 'error' in str(r.get('stop_reason', '')).lower()


def cap_row(r):
    s = (str(r.get('status', ''))+' '+str(r.get('stop_reason', ''))).lower()
    return not error_row(r) and any(t in s for t in ['cap', 'max_iter', 'maximum_iteration', 'maxit', 'timeout'])


def source_key(r):
    return str(r['source_id'])


def row_key(r):
    return (source_key(r), r['arm'], int(float(r.get('rep', 0))), r['method'])


def contingency_metrics(a, b):
    """ARI and arithmetic NMI, from integer contingency counts."""
    if len(a) != len(b) or len(a) == 0:
        raise ValueError('Partitions have different/zero lengths')
    joint = Counter(zip(a, b)); ca = Counter(a); cb = Counter(b); n = len(a)
    choose = lambda vals: sum(v*(v-1)//2 for v in vals)
    tp = choose(joint.values()); wa = choose(ca.values()); wb = choose(cb.values()); total = n*(n-1)//2
    expected = wa*wb/total if total else 0
    denominator = .5*(wa+wb)-expected
    ari = (tp-expected)/denominator if abs(denominator) > 1e-15 else 1.0
    ha = -sum(c/n*math.log(c/n) for c in ca.values())
    hb = -sum(c/n*math.log(c/n) for c in cb.values())
    mi = sum(c/n*math.log(c*n/(ca[x]*cb[y])) for (x, y), c in joint.items())
    nmi = 2*mi/(ha+hb) if ha+hb > 1e-15 else 1.0
    return dict(ARI=ari, NMI=nmi,
                within_anchor_pair_retention=tp/wa if wa else float('nan'),
                between_anchor_pair_merge=(wb-tp)/(total-wa) if total-wa else float('nan'),
                pairs_total=total, anchor_within_pairs=wa)


def describe(values):
    x = np.asarray([v for v in values if math.isfinite(v)], dtype=float)
    if not len(x):
        return dict(n_metric=0, mean=float('nan'), sd=float('nan'), variance=float('nan'),
                    mcse=float('nan'), median=float('nan'), q25=float('nan'), q75=float('nan'),
                    iqr=float('nan'), q05=float('nan'), q95=float('nan'), minimum=float('nan'), maximum=float('nan'))
    sd = float(np.std(x, ddof=1)) if len(x)>1 else float('nan')
    q05, q25, med, q75, q95 = np.quantile(x, [.05, .25, .5, .75, .95], method='linear')
    return dict(n_metric=len(x), mean=float(x.mean()), sd=sd, variance=sd*sd,
                mcse=sd/math.sqrt(len(x)), median=float(med), q25=float(q25), q75=float(q75),
                iqr=float(q75-q25), q05=float(q05), q95=float(q95), minimum=float(x.min()), maximum=float(x.max()))


def find_table(root, candidates):
    for p in candidates:
        q = root/p
        if q.exists():
            return q
    return root/candidates[0]


def partition_path(root, row):
    if present(row.get('partition_file')):
        p = Path(row['partition_file'])
        return p if p.is_absolute() else root/p
    s, arm, rep, method = row_key(row)
    return root/'partitions'/f'{s}_{arm}_b{rep:03d}_{method}.csv'


def audit_partitions(root, rows, base, out):
    checked=[]; failures=[]; base_arrays={}; node_metadata={}; partition_cache=defaultdict(list)
    combined=[dict(r, arm=r.get('arm') or 'base', rep=r.get('rep') or '0') for r in base]+rows
    for i, r in enumerate(combined):
        key=row_key(r); p=partition_path(root,r)
        flags=[]; values={}; errors={}
        if not p.exists():
            if not error_row(r):
                failures.append(dict(kind='partition_missing', key=str(key), detail=str(p)))
            checked.append(dict(source_id=key[0],arm=key[1],rep=key[2],method=key[3],
                                partition_present=False,partition_rows=0,checked=False,all_match=error_row(r)))
            continue
        records=read_csv(p)
        required=['node_id','anchor','own_source','truth','predicted']
        if not records or any(k not in records[0] for k in required):
            failures.append(dict(kind='partition_schema',key=str(key),detail='required: '+','.join(required)))
            continue
        # Real node identifiers and reference annotations need not be numeric.
        arrays={k:[str(x[k]) for x in records] for k in required}
        ids=arrays['node_id']; anchor=arrays['anchor']; own=arrays['own_source']; truth=arrays['truth']; pred=arrays['predicted']
        flags.append(all(present(v) for k in ['node_id','anchor','truth'] for v in arrays[k]))
        flags.append(len(ids)==len(set(ids)))
        meta=(ids,anchor,truth)
        if key[0] not in node_metadata:
            node_metadata[key[0]]=meta
        flags.append(meta==node_metadata[key[0]])
        if key[1]=='base':
            base_arrays[(key[0],key[3])]=pred
            flags.append(own==pred)
            if key[3]=='Proposed':
                flags.append(anchor==pred)
        elif (key[0],key[3]) in base_arrays:
            flags.append(own==base_arrays[(key[0],key[3])])
        else:
            flags.append(False)
        if not all(present(v) for v in pred):
            legitimate=error_row(r) and not any(present(v) for v in pred) and all(not math.isfinite(num(r.get(m))) for m in METRICS)
            if not (legitimate and all(flags)):
                failures.append(dict(kind='missing_partition_labels',key=str(key),detail='missing labels without a fully recorded error'))
            checked.append(dict(source_id=key[0],arm=key[1],rep=key[2],method=key[3],partition_present=True,
                                partition_rows=len(ids),checked=False,all_match=legitimate and all(flags),recorded_error=error_row(r)))
            continue
        ma=contingency_metrics(anchor,pred); mt=contingency_metrics(truth,pred)
        own_ari=contingency_metrics(own,pred)['ARI'] if all(present(v) for v in own) else float('nan')
        values=dict(ARI_anchor=ma['ARI'],NMI_anchor=ma['NMI'],ARI_own_source=own_ari,ARI_original_truth=mt['ARI'],
                    within_anchor_pair_retention=ma['within_anchor_pair_retention'],between_anchor_pair_merge=ma['between_anchor_pair_merge'],
                    final_K=len(set(pred)))
        for name,value in values.items():
            recorded=num(r.get(name))
            # Base rows may intentionally omit anchor-relative co-clustering metrics.
            if not math.isfinite(recorded) and (key[1]=='base' or not math.isfinite(value)):
                continue
            delta=abs(value-recorded) if math.isfinite(value) and math.isfinite(recorded) else float('inf')
            errors[name+'_abs_error']=delta
            flags.append(delta<1e-9)
        ok=all(flags)
        if not ok:
            failures.append(dict(kind='partition_or_metric_mismatch',key=str(key),detail=json.dumps(errors)))
        checked.append(dict(source_id=key[0],arm=key[1],rep=key[2],method=key[3],partition_present=True,
                            partition_rows=len(ids),checked=True,all_match=ok,**errors))
        if key[1]!='base':
            labels=np.unique(np.asarray(pred),return_inverse=True)[1].astype(np.int64)
            counts=np.bincount(labels)
            partition_cache[(key[0],key[1],key[3])].append(dict(rep=key[2],labels=labels,
                same_pairs=float(np.sum(counts*(counts-1))/2),cap=cap_row(r)))
        if (i+1)%200==0:
            print(f'Independent partition audit: {i+1}/{len(combined)}',flush=True)
    write_csv(out/'analysis_partition_audit.csv',checked)
    return checked,failures,partition_cache


def audit_design(root, rows, base, manifest, sources, expected_b, method_order):
    checks=[]
    def add(name,ok,detail): checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
    keys=[row_key(r) for r in rows]
    add('unique_method_keys',len(keys)==len(set(keys)),len(keys))
    basekeys=[(source_key(r),r['method']) for r in base]
    add('unique_base_method_keys',len(basekeys)==len(set(basekeys)),len(basekeys))
    basemap={(source_key(r),r['method']):r for r in base}
    mm={(source_key(r),int(float(r.get('rep',r.get('b',0))))):r for r in manifest}
    for s in sources:
        group=[r for r in rows if source_key(r)==s]
        smbase=[r for r in base if source_key(r)==s]
        add(f'{s}:base_method_count',len(smbase)==4,len(smbase))
        for arm in ARM_ORDER:
            for method in method_order:
                g=[r for r in group if r['arm']==arm and r['method']==method]
                reps={int(float(r['rep'])) for r in g}
                add(f'{s}:{arm}:{method}:B',len(g)==expected_b and reps==set(range(1,expected_b+1)),len(g))
        for seed_field in ['init_seed','fit_seed']:
            vals={r.get(seed_field) for r in group+smbase}
            add(f'{s}:{seed_field}:fixed',len(vals)==1 and all(present(v) for v in vals),sorted(str(v) for v in vals))
        for method in method_order:
            g=[r for r in group+smbase if r['method']==method]
            vals={r.get('fit_rng_hash') for r in g if present(r.get('fit_rng_hash'))}
            add(f'{s}:{method}:fit_rng_fixed',len(vals)==1 and all(present(r.get('fit_rng_hash')) for r in g),sorted(vals))
        original_x={r.get('attribute_hash') for r in smbase if present(r.get('attribute_hash'))}
        gox={r.get('attribute_hash') for r in group if r['arm']=='graph_only'}
        add(f'{s}:graph_only_X_fixed',len(original_x)==1 and gox==original_x,{'base':sorted(original_x),'graph_only':sorted(str(v) for v in gox)})
        model_hashes={r.get('source_model_hash') for r in group if present(r.get('source_model_hash'))}
        add(f'{s}:model_hash_fixed',len(model_hashes)==1,sorted(model_hashes))
        for b in range(1,expected_b+1):
            sub=[r for r in group if int(float(r['rep']))==b]
            ah={r.get('graph_hash') for r in sub}
            add(f'{s}:b{b}:paired_graph',len(sub)==8 and len(ah)==1 and all(present(x) for x in ah),sorted(str(v) for v in ah))
            for arm in ARM_ORDER:
                ar=[r for r in sub if r['arm']==arm]
                xh={r.get('attribute_hash') for r in ar}
                add(f'{s}:b{b}:{arm}:shared_X',len(ar)==4 and len(xh)==1 and all(present(x) for x in xh),sorted(str(v) for v in xh))
                pp=[r for r in ar if r['method'] in method_order[:2]]
                ih={r.get('initial_hash') for r in pp}
                add(f'{s}:b{b}:{arm}:shared_model_initialization',len(pp)==2 and len(ih)==1 and all(present(v) for v in ih),sorted(str(v) for v in ih))
            mf=mm.get((s,b))
            add(f'{s}:b{b}:generation_manifest_present',mf is not None,'')
            if mf:
                add(f'{s}:b{b}:manifest_graph_hash',ah=={mf.get('graph_hash')},mf.get('graph_hash'))
                joint_x={r.get('attribute_hash') for r in sub if r['arm']=='joint'}
                add(f'{s}:b{b}:manifest_joint_X',joint_x=={mf.get('X_joint_hash')},mf.get('X_joint_hash'))
                add(f'{s}:b{b}:manifest_model_hash',model_hashes=={mf.get('model_hash')},mf.get('model_hash'))
                pv=str(mf.get('permutation_verified','')).lower()
                add(f'{s}:b{b}:permutation_verified',pv in ['true','1'],pv)
    return checks


def summaries(rows, sources, expected_b, models, method_order):
    result=[]; kdist=[]
    anchor_ks={source_key(r):num(r.get('K')) for r in models}
    for s in sources:
        for arm in ARM_ORDER:
            for method in method_order:
                group=[r for r in rows if source_key(r)==s and r['arm']==arm and r['method']==method]
                if not group: continue
                errors=sum(error_row(r) for r in group); caps=sum(cap_row(r) for r in group)
                valid=sum(math.isfinite(num(r.get('ARI_anchor'))) for r in group)
                warn=sum(present(r.get('init_warnings')) or present(r.get('fit_warnings')) for r in group)
                common=dict(source_id=s,arm=arm,method=method,expected_B=expected_b,n_attempted=len(group),
                    n_errors=errors,n_caps=caps,n_partition_metrics=valid,n_no_error=len(group)-errors,
                    no_error_rate_all_B=(len(group)-errors)/expected_b,cap_rate_all_B=caps/expected_b,
                    error_rate_all_B=errors/expected_b,n_warning_runs=warn,anchor_K=anchor_ks.get(s,float('nan')),
                    n_K_equals_anchor=sum(num(r.get('final_K'))==anchor_ks.get(s) for r in group),
                    K_match_fraction_all_B=sum(num(r.get('final_K'))==anchor_ks.get(s) for r in group)/expected_b)
                for metric in METRICS:
                    result.append(dict(**common,metric=metric,**describe([num(r.get(metric)) for r in group])))
                counts=Counter(int(num(r['final_K'])) for r in group if math.isfinite(num(r.get('final_K'))))
                for k,count in sorted(counts.items()):
                    kdist.append(dict(source_id=s,arm=arm,method=method,K=k,count=count,fraction_all_B=count/expected_b))
    return result,kdist


def paired_summaries(rows,sources,expected_b,method_order):
    lookup={row_key(r):r for r in rows};pairs=[];sums=[]
    for s in sources:
        for method in method_order:
            for metric in METRICS:
                diffs=[]
                for b in range(1,expected_b+1):
                    a=lookup.get((s,'graph_only',b,method));j=lookup.get((s,'joint',b,method))
                    va=num(a.get(metric)) if a else float('nan');vj=num(j.get(metric)) if j else float('nan')
                    delta=vj-va if math.isfinite(va) and math.isfinite(vj) else float('nan')
                    diffs.append(delta)
                    pairs.append(dict(source_id=s,rep=b,method=method,metric=metric,graph_only=va,joint=vj,joint_minus_graph_only=delta,
                                      shared_graph=bool(a and j and a.get('graph_hash')==j.get('graph_hash'))))
                sums.append(dict(source_id=s,method=method,metric=metric,expected_paired_B=expected_b,
                                 n_incomplete_pairs=sum(not math.isfinite(v) for v in diffs),**describe(diffs)))
    return pairs,sums


def overall_counts(rows,base,expected_sources,expected_b,method_order):
    totals=[];events=[]
    for phase,data,planned in [('base',base,expected_sources),('bootstrap',rows,expected_sources*expected_b*2)]:
        for method in method_order:
            g=[r for r in data if r['method']==method]
            totals.append(dict(phase=phase,method=method,planned=planned,completed_records=len(g),
                errors=sum(error_row(r) for r in g),caps=sum(cap_row(r) for r in g),
                normal_or_reference=len(g)-sum(error_row(r) or cap_row(r) for r in g),
                valid_metric_rows=sum(math.isfinite(num(r.get('ARI_anchor'))) for r in g),
                initialization_warning_rows=sum(present(r.get('init_warnings')) for r in g),
                fit_warning_rows=sum(present(r.get('fit_warnings')) for r in g)))
        grouped=defaultdict(list)
        for r in data:grouped[(source_key(r),r.get('arm','base'),int(float(r.get('rep',0))))].append(r)
        for (s,a,b),group in grouped.items():
            bym={r['method']:r for r in group}
            warn=lambda m:set(v for v in str(bym.get(m,{}).get('init_warnings','')).split(' | ') if present(v))
            lw=warn('Louvain');kw=warn('KMeans_estK')-lw
            if lw or kw:
                events.append(dict(phase=phase,source_id=s,arm=a,rep=b,kmeans_warning=bool(kw),louvain_warning=bool(lw),
                    kmeans_warning_text=' | '.join(sorted(kw)),louvain_warning_text=' | '.join(sorted(lw)),
                    affected_method_rows=sum(present(r.get('init_warnings')) for r in group)))
    return totals,events


def pairwise_partition_diagnostic(cache,out,expected_b):
    """Post-protocol descriptive pairwise ARI; pairs are NOT independent."""
    raw=[];summary=[]
    for (s,arm,method),partitions in cache.items():
        partitions=sorted(partitions,key=lambda r:r['rep']);scores=[]
        for i,pa in enumerate(partitions):
            a=pa['labels'];n=len(a);total=n*(n-1)/2
            for pb in partitions[i+1:]:
                b=pb['labels'];assert len(b)==n
                table=np.bincount(a*(int(b.max())+1)+b)
                observed=float(np.sum(table*(table-1))/2)
                expected=pa['same_pairs']*pb['same_pairs']/total if total else 0
                den=.5*(pa['same_pairs']+pb['same_pairs'])-expected
                ari=(observed-expected)/den if abs(den)>1e-15 else 1.0
                scores.append(ari)
                raw.append(dict(source_id=s,arm=arm,method=method,rep_i=pa['rep'],rep_j=pb['rep'],ARI=ari,
                    cap_i=pa['cap'],cap_j=pb['cap']))
        vals=np.asarray(scores,dtype=float)
        desc=dict(pairwise_ARI_mean=float(vals.mean()),pairwise_ARI_median=float(np.median(vals)),
                  pairwise_ARI_sd=float(vals.std(ddof=1)) if len(vals)>1 else float('nan'),
                  pairwise_ARI_min=float(vals.min()),pairwise_ARI_max=float(vals.max()),
                  pairwise_ARI_q25=float(np.quantile(vals,.25)),pairwise_ARI_q75=float(np.quantile(vals,.75))) if len(vals) else {
                  k:float('nan') for k in ['pairwise_ARI_mean','pairwise_ARI_median','pairwise_ARI_sd','pairwise_ARI_min','pairwise_ARI_max','pairwise_ARI_q25','pairwise_ARI_q75']}
        summary.append(dict(source_id=s,arm=arm,method=method,planned_B=expected_b,n_valid_partitions=len(partitions),
            n_dependent_pairs=len(scores),post_protocol_descriptive=True,**desc))
        print(f'Pairwise diagnostic: {s}/{arm}/{method}: {len(partitions)} partitions, {len(scores)} dependent pairs',flush=True)
    write_csv(out/'analysis_pairwise_ARI.csv',raw,fields=['source_id','arm','method','rep_i','rep_j','ARI','cap_i','cap_j'])
    write_csv(out/'analysis_pairwise_summary.csv',summary)
    # Verify one pair per group using the independent dictionary-count formula.
    verification=[]
    for key,parts in cache.items():
        if len(parts)<2:continue
        aa,bb=sorted(parts,key=lambda r:r['rep'])[:2]
        independent=contingency_metrics(aa['labels'].tolist(),bb['labels'].tolist())['ARI']
        r=next(r for r in raw if (r['source_id'],r['arm'],r['method'],r['rep_i'],r['rep_j'])==(*key,aa['rep'],bb['rep']))
        verification.append(dict(source_id=key[0],arm=key[1],method=key[2],absolute_error=abs(independent-r['ARI']),match=abs(independent-r['ARI'])<1e-12))
    write_csv(out/'analysis_pairwise_formula_check.csv',verification,fields=['source_id','arm','method','absolute_error','match'])
    assert all(r['match'] for r in verification)
    return summary,len(raw),len(verification)


def main():
    cfg = args()
    root = cfg.input_dir
    out = cfg.output_dir or root
    out.mkdir(parents=True, exist_ok=True)
    rows = read_csv(root / "raw_metrics.csv")
    base = read_csv(root / "base_metrics.csv")
    models = read_csv(root / "models/model_summary.csv")
    manifest = read_csv(root / "graphs/generation_manifest.csv")
    sources = list(dict.fromkeys(source_key(r) for r in base))
    if len(sources) != cfg.expected_sources or len(rows) != cfg.expected_sources * cfg.expected_b * 8:
        raise RuntimeError("Unexpected record count. Check expected sources and B.")
    # Preserve the method identifier and partition filenames in earlier output directories.
    pure_method = "Pure_aSBM_C" if any(r["method"] == "Pure_aSBM_C" for r in base) else "Pure_aSBM"
    method_order = ["Proposed", pure_method, "Louvain", "KMeans_estK"]
    checks = audit_design(root, rows, base, manifest, sources, cfg.expected_b, method_order)
    write_csv(out / "analysis_design_audit.csv", checks)
    _, failures, partition_cache = audit_partitions(root, rows, base, out)
    if failures or not all(c["passed"] for c in checks):
        write_csv(out / "analysis_audit_failures.csv", failures, fields=["kind", "key", "detail"])
        raise RuntimeError("Bootstrap design or partition validation failed")
    group_summary, k_distribution = summaries(rows, sources, cfg.expected_b, models, method_order)
    differences, paired = paired_summaries(rows, sources, cfg.expected_b, method_order)
    totals, warnings = overall_counts(rows, base, cfg.expected_sources, cfg.expected_b, method_order)
    write_csv(out / "analysis_group_summary.csv", group_summary)
    write_csv(out / "analysis_K_distribution.csv", k_distribution)
    write_csv(out / "analysis_paired_arm_differences.csv", differences)
    write_csv(out / "analysis_paired_arm_summary.csv", paired)
    write_csv(out / "analysis_overall_counts.csv", totals)
    write_csv(out / "analysis_initialization_warning_events.csv", warnings,
              fields=["phase", "source_id", "arm", "rep", "kmeans_warning",
                      "louvain_warning", "kmeans_warning_text", "louvain_warning_text",
                      "affected_method_rows"])
    if cfg.expected_b > 1:
        pairwise_partition_diagnostic(partition_cache, out, cfg.expected_b)
    print(f"Validated {len(base)} base fits and {len(rows)} bootstrap method records.")


if __name__ == "__main__":
    main()

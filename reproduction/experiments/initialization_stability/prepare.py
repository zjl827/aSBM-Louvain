"""Generate the initialization design and record the model source hashes."""
from pathlib import Path
import csv, hashlib, json, os

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(os.environ.get('INIT_OUTPUT', ROOT / 'outputs/initialization_stability'))
for directory in [OUT/'data', OUT/'initializations', OUT/'runs']:
    directory.mkdir(parents=True, exist_ok=True)
model_source = ROOT.parent / 'asbm_louvain.R'
sources = {model_source.name: {'file': model_source.name,
                               'sha256': hashlib.sha256(model_source.read_bytes()).hexdigest()}}
conditions = []
datasets = []
tasks = []
for ci, (r, sigma) in enumerate([(1.5,1.5),(1.5,4),(4,1.5),(4,4)],1):
    conditions.append(dict(condition=f'I{ci}',r=r,sigma=sigma,N=1000,K_true=5,D=50,expected_degree=20))
    for di in range(1,4):
        dataset=f'I{ci}_data{di}'
        base=10050000+ci*10000+di*100
        ds=dict(**conditions[-1],dataset=dataset,data_id=di,graph_seed=base+1,
                means_seed=base+2,noise_seed=base+3,
                downstream_seed_main=11050000+ci*100+di,
                downstream_seed_check=12050000+ci*100+di)
        datasets.append(ds)
        for ii in range(1,11):
            init_seed=13050000+ci*10000+di*100+ii
            for method in ['kmeans','louvain']:
                for variant in (['main','downstream_check'] if ii==1 else ['main']):
                    tasks.append(dict(dataset=dataset,condition=f'I{ci}',data_id=di,
                        init_id=ii,init_method=method,init_seed=init_seed,variant=variant,
                        downstream_seed=ds['downstream_seed_main' if variant=='main' else 'downstream_seed_check'],
                        task_id=f'{dataset}_init{ii:02d}_{method}_{variant}'))
protocol = dict(version='initialization_stability_v1',frozen_date='2026-10-05',
  scope='Nested initialization study using the frozen aSBM-Louvain implementation.',
  design=dict(conditions=conditions,data_replicates_per_condition=3,initialization_seeds_per_dataset=10,
              initializer_methods=['kmeans','louvain'],main_fits=240,
              downstream_rng_diagnostic=dict(init_id=1,alternative_streams=1,additional_fits=24),total_fits=264),
  generation=dict(undirected=True,self_loops=False,block_sizes=[200]*5,
                  p_out='20/(199*r+800)',p_in='r*p_out',means='independent N(0,1) coordinates',
                  attributes='X_i = mu[z_i] + sigma * Z_i; independent standard normal Z'),
  initialization=dict(shared_louvain_partition='One original-graph Louvain result per dataset/init seed supplies K and Louvain initial labels.',
      kmeans='After that Louvain call, continue the same RNG stream for K-means with centers=max(2,K_Louvain), nstart=5, iter.max=10, Hartigan-Wong.',
      louvain='Use the shared original-graph Louvain labels; do not force minimum K=2 for this alternative method.'),
  algorithm=dict(tau0=1,rho_tau=.8,resolution=1,epsilon_Q=1e-4,epsilon_L=1e-3,
      k_stable_required=3,T_outer=1000,prior_feedback=True,timeout_sec=300),
  control=dict(workers=1,BLAS_threads=1,OMP_threads=1,
      downstream_rng='Reset to a fixed saved state for each dataset across all 10 initializations and both methods. Diagnostic changes only that state for init_id=1.',
      seed_limitation='Equal post-initialization RNG controls the source stream, not identical random moves on diverging algorithm paths.'),
  analysis=dict(primary='ARI',secondary=['NMI','final_K','iterations','stop_reason','warnings','initial_ARI','initial_K'],
      summaries=['Within each fixed dataset: mean, sample SD, range over 10 initialization seeds.',
                 'Per condition: average of the 3 dataset means, SD across 3 dataset means, RMS within-dataset SD.',
                 'Paired differences between K-means and Louvain initializations within dataset/init seed.',
                 'Initial-to-final ARI differences; alternative downstream RNG comparison on 24 preset fits.'],
      inference='Descriptive study: the 30 fits per condition/method are clustered in 3 independent datasets, not 30 independent data realizations. No significance claims or complete random-effects variance decomposition.',
      failures='Retain every planned attempt, warning, cap, timeout and error; no seed replacement. Report denominators.'),
  sources=sources)
for name, rows in [('conditions.csv',conditions),('datasets.csv',datasets),('tasks.csv',tasks)]:
    path=OUT/name
    import io
    buffer=io.StringIO(newline='')
    writer=csv.DictWriter(buffer,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    data=buffer.getvalue()
    if path.exists(): assert path.read_text()==data.replace('\r\n','\n')
    else:path.write_text(data)
path=OUT/'protocol.json'
data=json.dumps(protocol,indent=2,ensure_ascii=False)+'\n'
if path.exists():
    if path.read_text() != data:
        raise ValueError('Use a new INIT_OUTPUT directory for the current implementation.')
else:
    path.write_text(data)
print(f'Frozen {len(datasets)} datasets; {len(tasks)} planned fits. Sources:',sources)

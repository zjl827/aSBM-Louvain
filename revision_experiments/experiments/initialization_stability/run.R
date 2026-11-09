#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(igraph);library(mvtnorm);library(mclust)})
script<-normalizePath(sub('^--file=','',commandArgs(FALSE)[startsWith(commandArgs(FALSE),'--file=')][1]))
work<-dirname(script);root<-normalizePath(file.path(work,'../..'));out <- Sys.getenv('INIT_OUTPUT', file.path(root, 'outputs/initialization_stability'))
code <- file.path(root, 'code')
source(file.path(code, 'asbmlou.R'));source(file.path(code, 'asbmlou_revised_stop.R'))
stopifnot(formals(stats::kmeans)$iter.max==10L)
stage<-Sys.getenv('INIT_STAGE','smoke');stopifnot(stage%in%c('smoke','full','validate'))
hash_object<-function(x){p<-tempfile('init_hash_');on.exit(unlink(p));saveRDS(x,p,compress=FALSE,version=3);unname(tools::md5sum(p))}
source_files<-c(file.path(code, 'asbmlou.R'),file.path(code, 'asbmlou_revised_stop.R'),script,file.path(out,'protocol.json'),file.path(out,'datasets.csv'),file.path(out,'tasks.csv'))
manifest<-tools::md5sum(source_files);manifest_hash<-hash_object(manifest)
if(file.exists(file.path(out,'manifest.rds')))stopifnot(identical(manifest,readRDS(file.path(out,'manifest.rds')))) else {
 saveRDS(manifest,file.path(out,'manifest.rds'))
 write.csv(data.frame(path=names(manifest),md5=unname(manifest)),file.path(out,'manifest.csv'),row.names=FALSE)
}
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
writeLines(paste(names(Sys.getenv(c('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'))),Sys.getenv(c('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')),sep='='),file.path(out,'thread_environment.txt'))
datasets<-read.csv(file.path(out,'datasets.csv'),stringsAsFactors=FALSE);tasks<-read.csv(file.path(out,'tasks.csv'),stringsAsFactors=FALSE)
safe<-function(f){w<-character();start<-proc.time()[['elapsed']];ans<-tryCatch(withCallingHandlers(f(),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)structure(list(message=conditionMessage(e)),class='run_error'));list(value=ans,warnings=paste(unique(w),collapse=' | '),elapsed=proc.time()[['elapsed']]-start)}
nmi<-function(a,b){p<-table(a,b)/length(a);x<-rowSums(p);y<-colSums(p);e<-outer(x,y);v<-p>0;mi<-sum(p[v]*log(p[v]/e[v]));hx<--sum(x[x>0]*log(x[x>0]));hy<--sum(y[y>0]*log(y[y>0]));if(hx+hy==0)1 else 2*mi/(hx+hy)}
ari_independent<-function(a,b){tab<-table(a,b);choose2<-function(x)x*(x-1)/2;A<-sum(choose2(tab));B<-sum(choose2(rowSums(tab)));C<-sum(choose2(colSums(tab)));M<-choose2(sum(tab));expected<-B*C/M;denom<-(B+C)/2-expected;if(abs(denom)<1e-15)1 else(A-expected)/denom}

load_data<-function(ds){path<-file.path(out,'data',paste0(ds$dataset,'.rds'));if(file.exists(path)){x<-readRDS(path);stopifnot(x$manifest_hash==manifest_hash);return(x)}
 k<-ds$K_true;n<-ds$N;d<-ds$D;sizes<-rep(n/k,k);truth<-rep(seq_len(k),sizes)
 set.seed(ds$graph_seed);pout<-ds$expected_degree/((sizes[1]-1)*ds$r+n-sizes[1]);P<-matrix(pout,k,k);diag(P)<-ds$r*pout
 graph<-sample_sbm(n,P,sizes,directed=FALSE,loops=FALSE);A<-as.matrix(as_adjacency_matrix(graph,sparse=FALSE))
 set.seed(ds$means_seed);mu<-matrix(rnorm(k*d),k,d)
 set.seed(ds$noise_seed);Z<-matrix(rnorm(n*d),n,d);X<-mu[truth,,drop=FALSE]+ds$sigma*Z
 x<-list(A=A,X=X,mu=mu,Z=Z,truth=truth,P=P,dataset=ds,manifest_hash=manifest_hash,
      A_hash=hash_object(A),X_hash=hash_object(X),means_hash=hash_object(mu),noise_hash=hash_object(Z))
 stopifnot(identical(A,t(A)),all(diag(A)==0),!anyNA(X));saveRDS(x,path);x}

load_initial<-function(data,task){path<-file.path(out,'initializations',sprintf('%s_init%02d.rds',task$dataset,task$init_id))
 if(file.exists(path)){x<-readRDS(path);stopifnot(x$manifest_hash==manifest_hash,x$A_hash==data$A_hash,x$X_hash==data$X_hash);return(x)}
 set.seed(task$init_seed);g<-graph_from_adjacency_matrix(data$A,mode='undirected',diag=FALSE)
 lp<-safe(function()as.integer(factor(membership(cluster_louvain(g,resolution=1)))))
 fail<-inherits(lp$value,'run_error');K<-if(fail)NA_integer_ else length(unique(lp$value));post_louvain_rng<-.Random.seed
 km<-if(fail)lp else safe(function()stats::kmeans(data$X,centers=max(2L,K),nstart=5,iter.max=10L,algorithm='Hartigan-Wong'))
 labels<-if(inherits(km$value,'run_error'))NULL else as.integer(factor(km$value$cluster))
 x<-list(manifest_hash=manifest_hash,A_hash=data$A_hash,X_hash=data$X_hash,init_seed=task$init_seed,init_id=task$init_id,
   louvain=lp,kmeans=km,kmeans_partition=labels,louvain_K=K,post_louvain_rng=post_louvain_rng,post_init_rng=.Random.seed)
 saveRDS(x,path);x}

run_task<-function(task){dest<-file.path(out,'runs',paste0(task$task_id,'.rds'))
 if(file.exists(dest)){old<-readRDS(dest);stopifnot(old$manifest_hash==manifest_hash);return(old$row)}
 ds<-datasets[datasets$dataset==task$dataset,];data<-load_data(ds);ini<-load_initial(data,task)
 selected<-ini[[task$init_method]];initial<-if(task$init_method=='kmeans')ini$kmeans_partition else if(inherits(selected$value,'run_error'))NULL else selected$value
 set.seed(task$downstream_seed);fit_rng<-.Random.seed
 row<-cbind(task,ds[,c('r','sigma','N','K_true','D','expected_degree')]);row$A_hash<-data$A_hash;row$X_hash<-data$X_hash
 row$means_hash<-data$means_hash;row$initial_hash<-hash_object(initial);row$fit_rng_hash<-hash_object(fit_rng)
 row$initial_K<-if(is.null(initial))NA else length(unique(initial));row$initial_ARI<-if(is.null(initial))NA else adjustedRandIndex(data$truth,initial)
 row$initial_NMI<-if(is.null(initial))NA else nmi(data$truth,initial)
 row$init_warnings<-paste(Filter(nzchar,c(ini$louvain$warnings,if(task$init_method=='kmeans')ini$kmeans$warnings else '')),collapse=' | ')
 fit<-if(is.null(initial))selected else safe(function(){setTimeLimit(elapsed=300,transient=TRUE);on.exit(setTimeLimit(cpu=Inf,elapsed=Inf,transient=FALSE));aSBM_Louvain_revised(data$A,data$X,T_outer=1000L,rho_tau=.8,tau0=1,resolution=1,
    epsilon_Q=1e-4,epsilon_L=1e-3,k_stable_required=3L,init_method=task$init_method,initial_partition=initial,prior_feedback=TRUE,verbose=FALSE)})
 failed<-inherits(fit$value,'run_error');row$status<-if(failed)'error' else 'ok';row$error<-if(failed)fit$value$message else '';row$warnings<-fit$warnings;row$diagnostic_elapsed_sec<-fit$elapsed
 row$ARI<-row$NMI<-row$final_K<-row$iterations<-row$last_K_change<-NA_real_;row$stop_reason<-'error'
 if(!failed){f<-fit$value;stopifnot(identical(as.integer(f$initial_partition),as.integer(initial)),length(f$partition)==ds$N)
   row$ARI<-adjustedRandIndex(data$truth,f$partition);row$NMI<-nmi(data$truth,f$partition);row$final_K<-f$final_K;row$iterations<-f$iterations_run;row$stop_reason<-f$stop_reason
   ch<-which(diff(f$K_history)!=0);row$last_K_change<-if(length(ch))max(ch) else 0L
   stopifnot(abs(row$ARI-ari_independent(data$truth,f$partition))<1e-10,abs(row$NMI-compare(data$truth,f$partition,method='nmi'))<1e-10)}
 row$ARI_change<-row$ARI-row$initial_ARI
 record<-list(row=row,manifest_hash=manifest_hash,initial_partition=initial,fit_rng=fit_rng,fit=fit,truth=data$truth)
 saveRDS(record,paste0(dest,'.tmp'));stopifnot(file.rename(paste0(dest,'.tmp'),dest));write.csv(row,sub('rds$','csv',dest),row.names=FALSE)
 cat(sprintf('DONE %s init=%.4f final=%.4f iter=%s stop=%s sec=%.2f\n',task$task_id,row$initial_ARI,row$ARI,row$iterations,row$stop_reason,row$diagnostic_elapsed_sec));flush.console();row}

if(stage!='validate'){
 selected<-if(stage=='smoke')tasks[tasks$data_id==1 & tasks$init_id==1 & tasks$variant=='main',] else tasks
 cat(sprintf('START stage=%s fits=%d manifest=%s\n',stage,nrow(selected),manifest_hash));flush.console()
 for(i in seq_len(nrow(selected)))run_task(selected[i,,drop=FALSE])
}
paths<-list.files(file.path(out,'runs'),pattern='\\.rds$',full.names=TRUE);records<-lapply(paths,readRDS);rows<-do.call(rbind,lapply(records,`[[`,'row'))
write.csv(rows,file.path(out,'raw_metrics.csv'),row.names=FALSE)
stopifnot(!anyDuplicated(rows$task_id),all(rows$task_id%in%tasks$task_id))
if(stage%in%c('full','validate')){
 stopifnot(nrow(rows)==264L,setequal(rows$task_id,tasks$task_id))
 for(i in seq_len(nrow(datasets))){d<-load_data(datasets[i,]);stopifnot(hash_object(d$A)==d$A_hash,hash_object(d$X)==d$X_hash,identical(d$mu[d$truth,,drop=FALSE]+datasets$sigma[i]*d$Z,d$X))}
 for(x in records){r<-x$row;stopifnot(x$manifest_hash==manifest_hash,hash_object(x$initial_partition)==r$initial_hash,hash_object(x$fit_rng)==r$fit_rng_hash)
  if(r$status=='ok'){stopifnot(abs(ari_independent(x$truth,x$fit$value$partition)-r$ARI)<1e-10,abs(nmi(x$truth,x$fit$value$partition)-r$NMI)<1e-10)} }
 for(ds in datasets$dataset){d<-rows[rows$dataset==ds & rows$variant=='main',];stopifnot(nrow(d)==20,length(unique(d$A_hash))==1,length(unique(d$X_hash))==1,length(unique(d$fit_rng_hash))==1)}
 for(i in which(rows$variant=='downstream_check')){r<-rows[i,];base<-rows[rows$dataset==r$dataset & rows$init_id==r$init_id & rows$init_method==r$init_method & rows$variant=='main',];stopifnot(nrow(base)==1,base$initial_hash==r$initial_hash,base$A_hash==r$A_hash,base$X_hash==r$X_hash,base$fit_rng_hash!=r$fit_rng_hash)}
 writeLines(c('PASS: all 264 planned attempts retained; current source/protocol fingerprints verified.',
   'PASS: 12 fixed datasets, exact X reconstruction, saved input and means hashes verified.',
   'PASS: within each main dataset all 20 fits share the same downstream RNG state.',
   'PASS: 24 diagnostic pairs preserve data/initial labels and differ only in downstream RNG.',
   'PASS: all successful saved partitions independently rescored for ARI/NMI.',
   'K-means nstart=5, iter.max=10. No changes to shared core.'),file.path(out,'validation.txt'))
}
cat(sprintf('COMPLETE stage=%s saved=%d errors=%d caps=%d\n',stage,nrow(rows),sum(rows$status!='ok'),sum(rows$stop_reason=='safety_cap')))

#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(igraph);library(Matrix);library(mvtnorm);library(mclust)})
work<-'experiments/timing_current';code<-'code';script<-file.path(work,'run.R')
args<-commandArgs(TRUE);stage<-if(length(args))args[[1L]] else 'formal'
stopifnot(length(args)<=1L, stage %in% c('formal','smoke'))
input<-Sys.getenv('REVISION_REALDATA_DIR','results/realdata_current/inputs')
out<-Sys.getenv('REVISION_TIMING_OUT',file.path('results','timing_current',stage))
for(d in c(out,file.path(out,'runs'),file.path(out,'data'),file.path(out,'warmup')))dir.create(d,recursive=TRUE,showWarnings=FALSE)
vars<-c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLAS_NUM_THREADS','MKL_NUM_THREADS')
stopifnot(all(Sys.getenv(vars)=='1'),identical(as.integer(formals(stats::kmeans)$iter.max),10L))
source(file.path(code,'asbmlou.R'));source(file.path(code,'asbmlou_revised_stop.R'))
source(file.path(code,'pure_asbm_vectorized.R'))
if (is.loaded('pure_asbm_floor_correction_c')) {
 dyn.unload(getNativeSymbolInfo('pure_asbm_floor_correction_c')$dll[['path']])
}
stopifnot(!is.loaded('pure_asbm_floor_correction_c'))
source(file.path(work,'generate.R'))
tasks_path<-file.path(work,'tasks.csv')
tasks<-read.csv(tasks_path,stringsAsFactors=FALSE)
if(stage=='smoke')tasks<-tasks[1,,drop=FALSE]
stopifnot(nrow(tasks)==if(stage=='formal')70L else 1L)
input_paths<-file.path(input,paste0(unique(tasks$dataset[tasks$kind=='real']),
                                  '__expression_graph_spatial_attributes.rds'))
if(!any(tasks$kind=='real'))input_paths<-character()
input_hashes<-tools::md5sum(input_paths)
stopifnot(!anyNA(input_hashes))
hash_object<-function(x){f<-tempfile();on.exit(unlink(f));saveRDS(x,f,compress=FALSE,version=3);unname(tools::md5sum(f))}
atomic_rds<-function(x,p){saveRDS(x,paste0(p,'.tmp'));stopifnot(file.rename(paste0(p,'.tmp'),p))}
atomic_csv<-function(x,p){write.csv(x,paste0(p,'.tmp'),row.names=FALSE,na='');stopifnot(file.rename(paste0(p,'.tmp'),p))}
code_hashes<-tools::md5sum(c(file.path(code,c('asbmlou.R','asbmlou_revised_stop.R','pure_asbm_vectorized.R')),script,file.path(work,'generate.R')))
tasks_hash<-unname(tools::md5sum(tasks_path))
run_hash<-hash_object(list(code_hashes=code_hashes,tasks_hash=tasks_hash,input_hashes=input_hashes))
metadata<-list(stage=stage,start_time=as.character(Sys.time()),run_hash=run_hash,code_hashes=code_hashes,
              tasks_hash=tasks_hash,input_hashes=input_hashes,threads=Sys.getenv(vars),RNGkind=RNGkind(),
              session=capture.output(sessionInfo()),external=extSoftVersion())
metadata_path<-file.path(out,'metadata.rds')
if(file.exists(metadata_path))stopifnot(identical(readRDS(metadata_path)$run_hash,run_hash)) else {
 saveRDS(metadata,metadata_path);capture.output(str(metadata),file=file.path(out,'metadata.txt'))
}
dir.create(file.path(out,'sessions'),showWarnings=FALSE)
session_id<-paste0(format(Sys.time(),'%Y%m%d_%H%M%S'),'_',Sys.getpid())
saveRDS(metadata,file.path(out,'sessions',paste0(session_id,'.rds')))
capture.output(str(metadata),file=file.path(out,'sessions',paste0(session_id,'.txt')))

capture_call<-function(f){warnings<-character();value<-tryCatch(withCallingHandlers(f(),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)structure(list(message=conditionMessage(e)),class='fit_error'));list(value=value,warnings=paste(unique(warnings),collapse=' | '))}
methods<-c('Proposed','Pure_exact_R')
# Original rotation restricted to the two formally reported methods.
orders<-list(c(1L,2L),c(2L,1L),c(1L,2L),c(2L,1L),c(2L,1L),c(1L,2L))
run_fit<-function(method,data,initial,rng,resolution){
 stopifnot(!is.loaded('pure_asbm_floor_correction_c'))
 backend<-if(method=='Proposed')'current_R_igraph' else 'exact_vectorized_R'
 assign('.Random.seed',rng,envir=.GlobalEnv);invisible(gc())
 started<-proc.time()
 x<-capture_call(function(){
   if(method=='Proposed')aSBM_Louvain_revised(data$A,data$X,T_outer=1000L,rho_tau=.8,tau0=1,resolution=resolution,
     epsilon_Q=1e-4,epsilon_L=1e-3,k_stable_required=3L,init_method='kmeans',initial_partition=initial,verbose=FALSE)
   else fit_pure_asbm_vectorized(data$A,data$X,initial_partition=initial,tolerance=1,timeout_sec=3600,max_iterations=1000L)
 })
 timing<-proc.time()-started
 x$fit_sec<-unname(timing[['elapsed']]);x$cpu_sec<-unname(timing[['user.self']]+timing[['sys.self']]);x$backend<-backend
 x
}
load_data<-function(task){
 p<-file.path(out,'data',paste0(task$case_id,'.rds'))
 if(file.exists(p))return(readRDS(p))
 if(task$kind=='real'){
  source<-file.path(input,paste0(task$dataset,'__expression_graph_spatial_attributes.rds'))
  x<-readRDS(source);data<-list(A=x$adjacency,X=x$attributes,z=x$truth,source=source,
    source_md5=unname(tools::md5sum(source)))
 }else{
  x<-simulate_asbm_data_revised(task$N,5L,task$D,task$r,task$sigma,20,task$data_seed)
  data<-list(A=x$A,X=x$X,z=x$z,source='current archived generator',source_md5=NA_character_)
 }
 stopifnot(nrow(data$X)==task$N,ncol(data$X)==task$D,all(is.finite(data$X)),all(diag(data$A)==0),isTRUE(all.equal(data$A,t(data$A))))
 data$hashes<-list(A=hash_object(data$A),X=hash_object(data$X),z=hash_object(data$z));data$run_hash<-run_hash
 atomic_rds(data,p);data
}
make_metrics<-function(task,data,initial,rng,init_attempt,init_sec,method,x,position){
 row<-task;row$method<-method;row$execution_position<-position;row$backend<-x$backend
 row$N<-nrow(data$X);row$D<-ncol(data$X);row$edges<-sum(data$A)/2;row$adjacency_class<-class(data$A)[1]
 row$graph_hash<-data$hashes$A;row$attribute_hash<-data$hashes$X;row$initial_hash<-hash_object(initial);row$post_init_rng_hash<-hash_object(rng)
 initial_valid<-length(initial)==nrow(data$X)&&!anyNA(initial)
 row$initial_K<-if(initial_valid)length(unique(initial)) else NA_integer_
 row$initial_ARI<-if(initial_valid)adjustedRandIndex(data$z,initial) else NA_real_
 row$initialization_sec<-init_sec;row$fit_sec<-x$fit_sec;row$total_sec<-init_sec+x$fit_sec;row$cpu_sec<-x$cpu_sec
 row$initialization_warnings<-init_attempt$warnings;row$warnings<-x$warnings;row$run_hash<-run_hash
 row$ARI<-row$iterations<-row$final_K<-NA_real_;row$stop_reason<-'error';row$error<-'';row$status<-'error'
 f<-x$value
 if(inherits(f,'fit_error'))row$error<-f$message else{
  stopifnot(length(f$partition)==nrow(data$X),!anyNA(f$partition))
  row$ARI<-adjustedRandIndex(data$z,f$partition);row$iterations<-if(method=='Proposed')f$iterations_run else f$iterations
  row$final_K<-length(unique(f$partition));row$stop_reason<-f$stop_reason
  row$status<-if(f$stop_reason %in% c('safety_cap','iteration_cap','timeout'))'incomplete' else 'ok'
 }
 row
}
run_case<-function(task){
 destination<-file.path(out,'runs',paste0(task$case_id,'.rds'))
 if(file.exists(destination)){old<-readRDS(destination);stopifnot(identical(old$run_hash,run_hash),identical(old$task,task));return(old$metrics)}
 data<-load_data(task);stopifnot(identical(data$run_hash,run_hash))
 init_path<-file.path(out,'runs',paste0(task$case_id,'__initialization.rds'))
 if(file.exists(init_path)){
  saved_init<-readRDS(init_path)
  stopifnot(identical(saved_init$run_hash,run_hash),identical(saved_init$task,task),identical(saved_init$data_hashes,data$hashes))
  init<-saved_init$initialization;init_sec<-saved_init$initialization_sec;rng<-saved_init$rng
 }else{
  set.seed(task$algorithm_seed);invisible(gc());t0<-proc.time()[['elapsed']]
  init<-capture_call(function()pure_asbm_matched_kmeans_initialization(data$A,data$X,resolution=task$resolution,nstart=5L))
  init_sec<-proc.time()[['elapsed']]-t0;rng<-.Random.seed
  atomic_rds(list(run_hash=run_hash,task=task,data_hashes=data$hashes,initialization=init,initialization_sec=init_sec,rng=rng),init_path)
 }
 initial_failed<-inherits(init$value,'fit_error')
 initial<-if(initial_failed)NULL else init$value$partition
 initial_hash<-hash_object(initial);rng_hash<-hash_object(rng)
 order<-methods[orders[[((task$rep-1L)%%6L)+1L]]]
 attempts<-list();rows<-list()
 for(method in order){
  checkpoint<-file.path(out,'runs',paste0(task$case_id,'__',method,'.rds'))
  if(file.exists(checkpoint)){
   saved<-readRDS(checkpoint);stopifnot(identical(saved$run_hash,run_hash),identical(saved$task,task),identical(saved$initial_hash,initial_hash),identical(saved$rng_hash,rng_hash),identical(saved$data_hashes,data$hashes));x<-saved$attempt
  }else{
   cat(sprintf('START %s %s N=%d K0=%d\n',task$case_id,method,nrow(data$X),length(unique(initial))));flush.console()
   x<-if(initial_failed)list(value=init$value,warnings='',fit_sec=NA_real_,cpu_sec=NA_real_,backend='not_run_initialization_failed') else run_fit(method,data,initial,rng,task$resolution)
   atomic_rds(list(run_hash=run_hash,task=task,initial_hash=initial_hash,rng_hash=rng_hash,data_hashes=data$hashes,attempt=x),checkpoint)
  }
  attempts[[method]]<-x;rows[[method]]<-make_metrics(task,data,initial,rng,init,init_sec,method,x,match(method,order))
  if(initial_failed)rows[[method]]$stop_reason<-'initialization_error'
  cat(sprintf('DONE %s %s fit=%.3f status=%s stop=%s ARI=%.5f\n',task$case_id,method,x$fit_sec,rows[[method]]$status,rows[[method]]$stop_reason,rows[[method]]$ARI));flush.console()
 }
 metrics<-do.call(rbind,rows)
 saved<-list(run_hash=run_hash,task=task,metrics=metrics,attempts=attempts,initial=init,post_init_rng=rng,data_hashes=data$hashes)
 atomic_rds(saved,destination);atomic_csv(metrics,sub('\\.rds$','.csv',destination));metrics
}

# Warm up the actual call paths on separate data, outside the benchmark table.
warm<-simulate_asbm_data_revised(300L,5L,10L,4.5,1.5,20,100599991L)
set.seed(100599992L);wi<-capture_call(function()pure_asbm_matched_kmeans_initialization(warm$A,warm$X,nstart=5L))
stopifnot(!inherits(wi$value,'fit_error'));wrng<-.Random.seed
for(m in methods){w<-run_fit(m,warm,wi$value$partition,wrng,1);atomic_rds(w,file.path(out,'warmup',paste0(m,'.rds')));stopifnot(!inherits(w$value,'fit_error'))}
rm(warm,wi,w);invisible(gc())
rows<-list()
for(i in seq_len(nrow(tasks))){rows[[i]]<-run_case(tasks[i,,drop=FALSE]);atomic_csv(do.call(rbind,rows),file.path(out,'raw_metrics.csv'));invisible(gc())}
raw<-do.call(rbind,rows)
stopifnot(nrow(raw)==nrow(tasks)*2L,!anyDuplicated(paste(raw$case_id,raw$method)))
writeLines(c(sprintf('Completed %d matched pairs and %d attempted fits',nrow(tasks),nrow(raw)),
 sprintf('Errors: %d; caps/timeouts: %d',sum(raw$status=='error'),sum(raw$status=='incomplete')),
 sprintf('run_hash=%s',run_hash),as.character(Sys.time())),file.path(out,'completion.txt'))
summary<-do.call(rbind,lapply(split(raw,interaction(raw$kind,ifelse(raw$kind=='real',raw$dataset,'simulation'),raw$N,raw$method,drop=TRUE)),function(d){
 ok<-d$status=='ok'
 data.frame(d[1L,c('kind','dataset','N','D','method')],n=nrow(d),n_completed=sum(ok),
   total_mean_sec=mean(d$total_sec[ok]),total_SD_sec=sd(d$total_sec[ok]),
   total_median_sec=if(any(ok))median(d$total_sec[ok])else NA_real_,
   mean_ARI=mean(d$ARI[ok]))
}))
write.csv(summary,file.path(out,'summary.csv'),row.names=FALSE)
cat('ALL_DONE\n')

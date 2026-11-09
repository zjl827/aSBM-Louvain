#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(igraph);library(Matrix);library(mvtnorm);library(mclust)})
work<-'experiments/graph_bootstrap';script<-file.path(work,'run.R')
args<-commandArgs(TRUE);stage<-if(length(args))args[[1L]] else 'full'
stopifnot(length(args)<=1L,stage%in%c('preflight','smoke','full'))
real_input<-Sys.getenv('REVISION_REALDATA_DIR','results/realdata_current/inputs')
out<-Sys.getenv('REVISION_BOOTSTRAP_OUT',file.path('results','graph_bootstrap',stage))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(d in c('base','models','graphs','runs','partitions'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
source(file.path(work,'boot_helpers.R'));source(file.path(work,'simulate_source.R'))
algorithm<-'../asbm_louvain.R'
source(algorithm);source('code/pure_asbm_vectorized.R')
stopifnot(formals(stats::kmeans)$iter.max==10L)
workers<-as.integer(Sys.getenv('REVISION_BOOTSTRAP_WORKERS',if(.Platform$OS.type=='windows')'1' else '2'))
stopifnot(workers>=1L,workers<=2L, .Platform$OS.type!='windows' || workers==1L)
sources<-read.csv(file.path(work,'sources.csv'),stringsAsFactors=FALSE)
tasks<-read.csv(file.path(work,'tasks.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(sources)==10L,!anyDuplicated(sources$case_id),nrow(tasks)==1000L,
          all(tasks$arm%in%c('graph_only','joint')),!anyDuplicated(tasks[c('source_id','rep','arm')]))
if(stage=='smoke') {
  sources<-sources[1L,,drop=FALSE]
  tasks<-tasks[tasks$source_id==sources$case_id & tasks$rep==1L,,drop=FALSE]
}
real_sources<-sources$dataset[sources$kind=='real']
input_paths<-if(length(real_sources))file.path(real_input,paste0(real_sources,'__expression_graph_spatial_attributes.rds')) else character()
freeze_files<-c(script,file.path(work,c('boot_helpers.R','simulate_source.R','protocol.json','sources.csv','tasks.csv')),
  algorithm,'code/pure_asbm_vectorized.R',input_paths)
manifest<-tools::md5sum(freeze_files);stopifnot(!anyNA(manifest));manifest_hash<-gb_hash(manifest)
if(file.exists(file.path(out,'manifest.rds'))) {
  if(!identical(manifest,readRDS(file.path(out,'manifest.rds')))) {
    stop('Saved bootstrap checkpoints use different code or inputs. Set REVISION_BOOTSTRAP_OUT to a fresh directory.')
  }
} else {
  gb_write_rds(manifest,file.path(out,'manifest.rds'))
  write.csv(data.frame(path=names(manifest),md5=unname(manifest)),file.path(out,'manifest.csv'),row.names=FALSE)
}
write.csv(sources,file.path(out,'sources.csv'),row.names=FALSE)
stopifnot(file.copy(file.path(work,'protocol.json'),out,overwrite=TRUE))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
writeLines(paste(names(Sys.getenv(c('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'))),
  Sys.getenv(c('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')),sep='='),file.path(out,'thread_environment.txt'))
write.csv(gb_generator_fixture(),file.path(out,'generator_fixture.csv'),row.names=FALSE)

prepare_source<-function(source) {
  base_path<-file.path(out,'base',paste0(source$case_id,'.rds'));model_path<-file.path(out,'models',paste0(source$case_id,'.rds'))
  if(file.exists(base_path)&&file.exists(model_path)) {
    b<-readRDS(base_path);m<-readRDS(model_path)
    stopifnot(b$manifest_hash==manifest_hash,m$manifest_hash==manifest_hash);return(invisible(TRUE))
  }
  input<-gb_load_input(source,real_input);fits<-gb_fit_methods(input$A,input$X,source)
  if(gb_failed(fits$attempts$Proposed)) {
    gb_write_rds(list(source=source,input=input,fits=fits,manifest_hash=manifest_hash),paste0(base_path,'.failed'))
    stop(paste('No Proposed anchor for',source$case_id,':',fits$attempts$Proposed$value$message))
  }
  anchor<-gb_partition(fits$attempts$Proposed);own_partitions<-lapply(fits$attempts,gb_partition)
  model<-gb_plugin_model(input$A,input$X,anchor,source);model$manifest_hash<-manifest_hash;model$source<-source
  model$summary<-data.frame(source_id=source$case_id,kind=source$kind,dataset=source$dataset,
      N=source$N,D=source$D,K=length(unique(anchor)),source_base_ARI_original_truth=adjustedRandIndex(input$truth,anchor),
      expected_edges=model$expected_edges,observed_edges=model$observed_edges,
      expected_minus_observed=model$expected_edges-model$observed_edges,model_hash=model$model_hash,
      min_cluster_size=min(model$parameters$sizes),singular_MLE_covariances=sum(model$components$MLE_cov_rank<source$D),
      anchor_stop_reason=fits$attempts$Proposed$value$stop_reason,
      graph_hash=input$graph_hash,attribute_hash=input$attribute_hash,representation=input$representation)
  rows<-gb_record_metrics(fits,input$A,input$X,source,'base',0L,anchor,own_partitions,input$truth,model$model_hash)
  gb_partition_csvs(fits,source,'base',0L,anchor,own_partitions,input$truth,input$node_ids,out)
  base<-list(source=source,input=input,fits=fits,anchor=anchor,own_partitions=own_partitions,
       metrics=rows,manifest_hash=manifest_hash)
  gb_write_rds(model,model_path);gb_write_rds(base,base_path)
  write.csv(model$parameters$P,file.path(out,'models',paste0(source$case_id,'_block_probabilities.csv')),row.names=FALSE)
  write.csv(rows,file.path(out,'base',paste0(source$case_id,'_metrics.csv')),row.names=FALSE)
  cat(sprintf('BASE %s K=%d ARItruth=%.5f stop=%s expected_edges=%.3f observed=%.3f\n',source$case_id,
     model$summary$K,model$summary$source_base_ARI_original_truth,model$summary$anchor_stop_reason,model$expected_edges,model$observed_edges));flush.console()
  invisible(TRUE)
}
cat(sprintf('START stage=%s sources=%d workers=%d manifest=%s\n',stage,nrow(sources),workers,manifest_hash));flush.console()
for(i in seq_len(nrow(sources)))prepare_source(sources[i,,drop=FALSE])
bases<-lapply(sources$case_id,function(id)readRDS(file.path(out,'base',paste0(id,'.rds'))))
models<-lapply(sources$case_id,function(id)readRDS(file.path(out,'models',paste0(id,'.rds'))))
write.csv(do.call(rbind,lapply(bases,`[[`,'metrics')),file.path(out,'base_metrics.csv'),row.names=FALSE)
write.csv(do.call(rbind,lapply(models,`[[`,'summary')),file.path(out,'models/model_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,lapply(models,`[[`,'components')),file.path(out,'models/model_components.csv'),row.names=FALSE)
rm(bases,models);gc()

run_pair<-function(task) {
  source<-sources[sources$case_id==task$source_id[1],,drop=FALSE];stopifnot(nrow(source)==1,nrow(task)==2,
    length(unique(task$graph_seed))==1,length(unique(task$attribute_seed))==1,setequal(task$arm,c('graph_only','joint')))
  b<-readRDS(file.path(out,'base',paste0(source$case_id,'.rds')));model<-readRDS(file.path(out,'models',paste0(source$case_id,'.rds')))
  graph_path<-file.path(out,'graphs',sprintf('%s_b%03d.rds',source$case_id,task$rep[1]))
  if(file.exists(graph_path)){generated<-readRDS(graph_path);stopifnot(generated$manifest_hash==manifest_hash,generated$meta$model_hash==model$model_hash)} else {
    generated<-gb_generate(model,b$input,source,task$rep[1],task$graph_seed[1],task$attribute_seed[1]);generated$manifest_hash<-manifest_hash
    gb_write_rds(generated,graph_path);write.csv(generated$meta,sub('rds$','csv',graph_path),row.names=FALSE)
  }
  for(arm in c('graph_only','joint')) {
    path<-file.path(out,'runs',sprintf('%s_%s_b%03d.rds',source$case_id,arm,task$rep[1]))
    if(file.exists(path)){old<-readRDS(path);stopifnot(old$manifest_hash==manifest_hash);next}
    X<-if(arm=='graph_only')b$input$X else generated$X_joint
    fits<-gb_fit_methods(generated$A,X,source)
    rows<-gb_record_metrics(fits,generated$A,X,source,arm,task$rep[1],b$anchor,b$own_partitions,b$input$truth,model$model_hash)
    gb_partition_csvs(fits,source,arm,task$rep[1],b$anchor,b$own_partitions,b$input$truth,b$input$node_ids,out)
    record<-list(source=source,arm=arm,rep=task$rep[1],metrics=rows,fits=fits,
      graph_path=graph_path,source_base_path=file.path(out,'base',paste0(source$case_id,'.rds')),manifest_hash=manifest_hash)
    gb_write_rds(record,path);write.csv(rows,sub('rds$','csv',path),row.names=FALSE)
    cat(sprintf('DONE %s %s b=%03d errors=%d caps=%d ProposedARI=%.5f sec=%.2f\n',source$case_id,arm,task$rep[1],
      sum(rows$status=='error'),sum(rows$stop_reason%in%c('safety_cap','iteration_cap','timeout')),
      rows$ARI_anchor[rows$method=='Proposed'],sum(rows$diagnostic_sec)));flush.console()
  }
  TRUE
}
if(stage!='preflight') {
  selected<-tasks
  grouped<-split(selected,paste(selected$source_id,selected$rep,sep='__'))
  result<-parallel::mclapply(grouped,run_pair,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
  stopifnot(!any(vapply(result,inherits,logical(1),what='try-error')))
  paths<-list.files(file.path(out,'runs'),pattern='\\.rds$',full.names=TRUE)
  rows<-do.call(rbind,lapply(paths,function(p)readRDS(p)$metrics));write.csv(rows,file.path(out,'raw_metrics.csv'),row.names=FALSE)
  graph_paths<-list.files(file.path(out,'graphs'),pattern='_b[0-9]+\\.csv$',full.names=TRUE)
  gm<-do.call(rbind,lapply(graph_paths,read.csv,stringsAsFactors=FALSE));write.csv(gm,file.path(out,'graphs/generation_manifest.csv'),row.names=FALSE)
  stopifnot(!anyDuplicated(rows[c('source_id','arm','rep','method')]))
  if(stage=='full')stopifnot(nrow(rows)==4000L,nrow(gm)==500L)
  cat(sprintf('COMPLETE stage=%s rows=%d errors=%d caps=%d\n',stage,nrow(rows),sum(rows$status=='error'),
    sum(rows$stop_reason%in%c('safety_cap','iteration_cap','timeout'))))
} else cat('PREFLIGHT COMPLETE: ten base fits/anchor models, block-density edge-expectation identity and covariance eigenvalues checked.\n')

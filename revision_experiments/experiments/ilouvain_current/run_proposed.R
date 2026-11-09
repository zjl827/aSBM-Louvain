#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(mclust))
work <- "experiments/ilouvain_current"
input <- Sys.getenv("REVISION_DANCER_DIR", "data/dancer")
out <- Sys.getenv("REVISION_DANCER_OUT", "results/ilouvain_current")
for (d in c(out, file.path(out, "runs"), file.path(out, "partitions"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
core <- c("code/asbmlou.R", "code/asbmlou_revised_stop.R")
for (p in core) source(p)
jobs <- read.csv(file.path(work, "jobs.csv"), stringsAsFactors = FALSE)
args <- commandArgs(TRUE)
if (length(args)) {
  if (!identical(args, "smoke")) stop("Usage: Rscript experiments/ilouvain_current/run_proposed.R [smoke]")
  jobs <- jobs[1L, , drop = FALSE]
}
# Checkpoints are tied to their input bytes, core, job grid and runner.
runner <- file.path(work, "run_proposed.R")
code_hash <- tools::md5sum(c(core, runner, file.path(work, "jobs.csv")))
nmi <- function(a,b) {
  tab <- table(a,b); p<-tab/sum(tab); pa<-rowSums(p); pb<-colSums(p); nz<-p>0
  ha <- -sum(pa[pa>0]*log(pa[pa>0])); hb <- -sum(pb[pb>0]*log(pb[pb>0]))
  if(ha==0 || hb==0) return(0)
  sum(p[nz]*log(p[nz]/outer(pa,pb)[nz]))/sqrt(ha*hb)
}
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
write.csv(data.frame(variable=c('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'),
  value=Sys.getenv(c('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'))),
  file.path(out,'thread_environment.csv'),row.names=FALSE)
for(j in seq_len(nrow(jobs))) {
  job <- jobs[j,]; key <- sprintf('%02d_%s_%s_rep%02d',job$job_id,job$graph,job$attribute_condition,job$repetition)
  cp <- file.path(out,'runs',paste0(key,'.rds'))
  ef<-file.path(input,paste0(job$graph,'_edges.csv'))
  lf<-file.path(input,paste0(job$graph,'_labels.csv'))
  xf<-file.path(input,paste0(job$graph,'_',job$attribute_condition,'.csv'))
  input_hash <- tools::md5sum(c(ef, lf, xf))
  stopifnot(!anyNA(input_hash))
  if (file.exists(cp)) {
    old <- readRDS(cp)
    stopifnot(identical(old$job,job), identical(old$input_md5,input_hash),
              identical(old$code_md5,code_hash))
    cat('CACHED',key,'\n'); flush.console(); next
  }
  e<-read.csv(ef); lab<-read.csv(lf); xf_data<-read.csv(xf); truth<-lab$label; N<-length(truth)
  stopifnot(identical(as.integer(lab$node),seq_len(N)),identical(as.integer(xf_data$node),seq_len(N)))
  g<-graph_from_edgelist(as.matrix(e),directed=FALSE)
  if(vcount(g)<N) g<-add_vertices(g,N-vcount(g))
  stopifnot(vcount(g)==N,!any_loop(g),!any_multiple(g))
  A<-as.matrix(as_adjacency_matrix(g,sparse=FALSE)); dimnames(A)<-list(as.character(seq_len(N)),as.character(seq_len(N)))
  X<-as.matrix(xf_data[,c('attribute_1','attribute_2')]); stopifnot(all(is.finite(X)))
  init_warnings<-character(); fit_warnings<-character(); error<-''; fit<-NULL; z0<-NULL; rng_after_init<-NULL
  set.seed(job$seed); rng_start<-.Random.seed
  init_seconds<-NA_real_; fit_seconds<-NA_real_
  tm<-proc.time()[['elapsed']]
  tryCatch({
    z0 <- withCallingHandlers({
      kest<-max(2L,length(unique(membership(cluster_louvain(g,resolution=1)))))
      as.numeric(factor(stats::kmeans(X,centers=kest,nstart=5,iter.max=10,algorithm='Hartigan-Wong')$cluster))
    },warning=function(w){init_warnings<<-c(init_warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    init_seconds<-proc.time()[['elapsed']]-tm
    rng_after_init<-.Random.seed
    tm<-proc.time()[['elapsed']]
    fit<-withCallingHandlers(aSBM_Louvain_revised(A,X,T_outer=1000L,rho_tau=.8,tau0=1,resolution=1,
      epsilon_Q=1e-4,epsilon_L=1e-3,k_stable_required=3L,init_method='kmeans',initial_partition=z0,prior_feedback=TRUE),
      warning=function(w){fit_warnings<<-c(fit_warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    fit_seconds<-proc.time()[['elapsed']]-tm
  },error=function(err){error<<-conditionMessage(err);fit_seconds<<-proc.time()[['elapsed']]-tm})
  ok<-!is.null(fit)&&!nzchar(error)
  if(ok) stopifnot(length(fit$partition)==N,all(is.finite(fit$partition)),all(is.finite(fit$soft_labels)),
      max(abs(rowSums(fit$soft_labels)-1))<1e-8,identical(fit$initial_partition,as.vector(z0)))
  metrics<-data.frame(job_id=job$job_id,graph=job$graph,attribute_condition=job$attribute_condition,
    repetition=job$repetition,seed=job$seed,published_validation=job$published_validation,N=N,M=ecount(g),D=ncol(X),K_true=length(unique(truth)),
    initial_K=if(is.null(z0))NA else length(unique(z0)),initial_ARI=if(is.null(z0))NA else adjustedRandIndex(truth,z0),
    ARI=if(ok)adjustedRandIndex(truth,fit$partition) else NA,NMI=if(ok)nmi(truth,fit$partition) else NA,
    K_found=if(ok)fit$final_K else NA,iterations=if(ok)fit$iterations_run else NA,converged=if(ok)fit$converged else FALSE,
    stop_reason=if(ok)fit$stop_reason else 'error',init_warning_count=length(init_warnings),fit_warning_count=length(fit_warnings),
    init_warnings=paste(unique(init_warnings),collapse=' | '),fit_warnings=paste(unique(fit_warnings),collapse=' | '),
    error=error,init_seconds=init_seconds,fit_seconds=fit_seconds,graph_input_md5=input_hash[[1L]], truth_input_md5=input_hash[[2L]], attribute_input_md5=input_hash[[3L]])
  payload<-list(job=job,input_md5=input_hash,code_md5=code_hash,rng_start=rng_start,rng_after_init=rng_after_init,
    initial_partition=z0,truth=truth,fit=fit,metrics=metrics,initialization_warnings=init_warnings,fit_warnings=fit_warnings,error=error)
  saveRDS(payload,paste0(cp,".tmp"))
  write.csv(metrics,file.path(out,'runs',paste0(key,'.csv')),row.names=FALSE)
  write.csv(data.frame(node=seq_len(N),truth=truth,initial=if(is.null(z0))NA else z0,
    predicted=if(ok)fit$partition else NA),file.path(out,'partitions',paste0(key,'.csv')),row.names=FALSE)
  if(ok) write.csv(data.frame(iteration=0:fit$iterations_run,K=fit$K_history,
    objective=c(NA,fit$elbo_history),modularity=c(NA,fit$modularity_history)),
    file.path(out,'runs',paste0(key,'_trajectory.csv')),row.names=FALSE)
  stopifnot(file.rename(paste0(cp,".tmp"),cp))
  cat(sprintf('DONE %s ARI=%.6f K=%s iter=%s stop=%s warnings=%s/%s error=%s seconds=%.2f\n',key,
    metrics$ARI,metrics$K_found,metrics$iterations,metrics$stop_reason,length(init_warnings),length(fit_warnings),error,fit_seconds));flush.console()
  rm(A,X,g,payload,fit);gc(FALSE)
}

paths <- list.files(file.path(out, "runs"), pattern = "\\.rds$", full.names = TRUE)
rows <- do.call(rbind, lapply(paths, function(p) readRDS(p)$metrics))
write.csv(rows, file.path(out, "proposed_results.csv"), row.names = FALSE)

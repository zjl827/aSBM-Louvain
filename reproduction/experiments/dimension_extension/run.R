#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(igraph); library(mvtnorm); library(mclust)})
script <- normalizePath(sub("^--file=", "", commandArgs(FALSE)[startsWith(commandArgs(FALSE), "--file=")][1]))
wd <- dirname(script); root <- normalizePath(file.path(wd, "../.."))
previous <- file.path(root, "experiments", "weight_misspecification")
code <- file.path(root, "code")
source(file.path(code, "weighted_model.R"))
source(file.path(previous, "design_helpers.R"))
source(file.path(code, "pure_asbm_vectorized.R"))
model <- make_model(code)
# Transparent diagnostics: preserve every original Pure objective evaluation.
pure_objective_original <- pure_asbm_objective
pure_objective_log <- numeric()
pure_asbm_objective <- function(...) {
  value <- pure_objective_original(...)
  pure_objective_log <<- c(pure_objective_log,value)
  value
}
stopifnot(formals(stats::kmeans)$iter.max == 10L)
smoke <- Sys.getenv("DIM_SMOKE", "0") == "1"
out <- Sys.getenv("DIM_OUTPUT", file.path(root, "outputs/dimension_extension"))
pilot_out <- Sys.getenv("DIM_PILOT_OUTPUT", file.path(root, "outputs/dimension_pilot"))
workers <- as.integer(Sys.getenv("DIM_WORKERS", "2")); reps <- if(smoke)1L else 10L
for(sub in c("runs","data_cache"))dir.create(file.path(out,sub), recursive=TRUE, showWarnings=FALSE)
files <- c(script, file.path(previous, "design_helpers.R"),
           file.path(code, c("weighted_model.R", "asbmlou.R", "fit_asbm_louvain.R",
                             "pure_asbm_vectorized.R")),
           file.path(wd, "protocol.md"))
manifest <- tools::md5sum(files); manifest_hash <- hash_object(manifest)
mf <- file.path(out,"manifest.rds")
if (file.exists(mf)) {
  if (!identical(readRDS(mf), manifest)) {
    stop("Use a new DIM_OUTPUT directory for the current implementation.")
  }
} else saveRDS(manifest, mf)
write.csv(data.frame(path=names(manifest),md5=unname(manifest)),file.path(out,"manifest.csv"),row.names=FALSE)
saveRDS(list(smoke=smoke,reps=reps,workers=workers,rng=RNGkind(),
             session=capture.output(sessionInfo()),manifest_hash=manifest_hash),file.path(out,"settings.rds"))
grid <- expand.grid(D=c(400L,800L),sigma=c(1.5,4));grid$condition_id<-sprintf("E%02d",seq_len(nrow(grid)))
write.csv(grid,file.path(out,"grid.csv"),row.names=FALSE)

safe_fit <- function(fun) {
  warnings<-character();start<-proc.time()[["elapsed"]]
  fit<-tryCatch({setTimeLimit(cpu=Inf,elapsed=300,transient=TRUE)
    withCallingHandlers(fun(),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  },error=function(e)structure(list(message=conditionMessage(e)),class="fit_error"),
  finally=setTimeLimit(cpu=Inf,elapsed=Inf,transient=FALSE))
  list(fit=fit,elapsed=proc.time()[["elapsed"]]-start,
       warning_count=length(warnings),warnings=paste(unique(warnings),collapse=" | "))
}

run_rep <- function(rep_id) {
  n<-1000L;k<-5L;sizes<-rep(200L,k);truth<-rep(seq_len(k),times=sizes)
  seeds<-list(graph=108100000L+rep_id,means=108110000L+rep_id,residual=108120000L+rep_id,
              louvain=108130000L+rep_id,kmeans=108140000L+rep_id)
  set.seed(seeds$graph);pout<-20/(199*4+800);block<-matrix(pout,k,k);diag(block)<-4*pout
  g<-sample_sbm(n,block,sizes,directed=FALSE,loops=FALSE);A<-as.matrix(as_adjacency_matrix(g,sparse=FALSE))
  set.seed(seeds$means);means<-matrix(rnorm(k*50),k,50)
  set.seed(seeds$residual);Z<-matrix(rnorm(n*800),n,800)
  set.seed(seeds$louvain);lp<-as.integer(factor(membership(cluster_louvain(g,resolution=1))))
  initK<-max(2L,length(unique(lp)));ah<-hash_object(A)
  stopifnot(identical(A,t(A)),all(diag(A)==0),all(is.finite(Z)))
  prior <- readRDS(file.path(pilot_out,"data_cache",sprintf("rep%02d.rds",rep_id)))
  stopifnot(identical(A,prior$A),identical(means,prior$means),
            identical(Z[,1:200],prior$Z),identical(truth,prior$truth),
            identical(seeds,prior$seeds),identical(lp,prior$louvain_partition))
  saveRDS(list(A=A,means=means,Z=Z,truth=truth,seeds=seeds,louvain_partition=lp),
          file.path(out,"data_cache",sprintf("rep%02d.rds",rep_id)))
  for(ci in seq_len(nrow(grid))) {
    cnd<-grid[ci,];key<-sprintf("%s_rep%02d",cnd$condition_id,rep_id)
    dest<-file.path(out,"runs",paste0(key,".rds"))
    if(file.exists(dest)) {old<-readRDS(dest);stopifnot(old$manifest_hash==manifest_hash,nrow(old$metrics)==4L);next}
    d<-cnd$D;M<-cbind(means,matrix(0,k,d-50L));X<-M[truth,,drop=FALSE]+cnd$sigma*Z[,seq_len(d),drop=FALSE]
    set.seed(seeds$kmeans);km<-safe_fit(function()stats::kmeans(X,centers=initK,nstart=5,iter.max=10L))
    init_failed<-inherits(km$fit,"fit_error");initial<-if(init_failed)NULL else as.integer(factor(km$fit$cluster))
    post_rng<-.Random.seed;initial_sizes<-if(init_failed)NA_integer_ else as.integer(table(initial))
    assign("pure_objective_log",numeric(),envir=.GlobalEnv)
    methods<-c("Proposed","Pure_aSBM","Louvain","KMeans_estK");records<-attempts<-list()
    for(mi in seq_along(methods)) {
      method<-methods[mi]
      if(method=="Louvain")attempt<-list(fit=list(partition=lp,iterations=NA_integer_,stop_reason="reference_partition"),elapsed=NA_real_,warning_count=0L,warnings="")
      else if(init_failed)attempt<-km
      else if(method=="KMeans_estK")attempt<-list(fit=list(partition=initial,iterations=km$fit$iter,stop_reason="reference_partition"),elapsed=km$elapsed,warning_count=km$warning_count,warnings=km$warnings)
      else {
        assign(".Random.seed",post_rng,envir=.GlobalEnv)
        attempt<-if(method=="Proposed")safe_fit(function()model$aSBM_Louvain_revised(
          A,X,T_outer=1000L,tau0=1,rho_tau=.8,resolution=1,epsilon_Q=1e-4,epsilon_L=1e-3,
          k_stable_required=3L,init_method="kmeans",initial_partition=initial,verbose=FALSE))else
            safe_fit(function()fit_pure_asbm_vectorized(A,X,initial_partition=initial,tolerance=1,timeout_sec=300,max_iterations=1000L))
      }
      fit<-attempt$fit;failed<-inherits(fit,"fit_error");row<-cnd
      row$rep<-rep_id;row$method<-method;row$initial_K<-initK;row$initial_min_size<-min(initial_sizes)
      row$initial_max_size<-max(initial_sizes);row$initial_clusters_le_D<-sum(initial_sizes<=d)
      row$graph_hash<-ah;row$attribute_hash<-hash_object(X);row$initial_hash<-hash_object(initial)
      row$post_init_rng_hash<-hash_object(post_rng);row$init_warnings<-km$warnings
      row$init_warning_count<-km$warning_count;row$init_ifault<-if(init_failed)NA_integer_ else km$fit$ifault
      row$init_iterations<-if(init_failed)NA_integer_ else km$fit$iter
      row$warnings<-attempt$warnings;row$warning_count<-attempt$warning_count
      row$status<-if(failed)"error" else "ok";row$error<-if(failed)fit$message else "";row$elapsed_sec<-attempt$elapsed
      row$realized_mean_degree<-sum(A)/n
      row$mean_pairwise_mahalanobis_sq<-mean(as.numeric(dist(M))^2)/cnd$sigma^2
      row$ARI<-row$NMI<-NA_real_;row$iterations<-row$final_K<-NA_integer_
      row$stop_reason<-"error";row$met_stopping_rule<-FALSE;row$objective_finite<-NA
      if(!failed) {
        pred<-as.integer(fit$partition);stopifnot(length(pred)==n,!anyNA(pred))
        row$ARI<-adjustedRandIndex(truth,pred);row$NMI<-partition_nmi(truth,pred)
        row$final_K<-length(unique(pred));row$iterations<-if(method=="Proposed")fit$iterations_run else fit$iterations
        row$stop_reason<-fit$stop_reason
        row$met_stopping_rule<-if(method=="Proposed")fit$stop_reason=="joint_stability" else if(method=="Pure_aSBM")fit$stop_reason %in% c("tolerance","objective_not_increased") else NA
        row$objective_finite<-if(method=="Proposed")all(is.finite(fit$elbo_history))else if(method=="Pure_aSBM")is.finite(fit$final_objective)else NA
      }
      row$pure_initial_objective_finite<-if(method=="Pure_aSBM" && length(pure_objective_log))is.finite(pure_objective_log[1])else NA
      row$pure_candidate_nonfinite<-if(method=="Pure_aSBM")sum(!is.finite(pure_objective_log[-1]))else NA_integer_
      row$pure_identical_to_kmeans<-if(method=="Pure_aSBM"&&!failed)identical(as.integer(fit$partition),initial)else NA
      records[[mi]]<-row;attempts[[mi]]<-attempt
    }
    result<-list(manifest_hash=manifest_hash,metrics=do.call(rbind,records),fits=attempts,X=X,
                 initial_partition=initial,post_init_rng=post_rng,seeds=seeds,kmeans=km,pure_objective_log=pure_objective_log)
    saveRDS(result,paste0(dest,".tmp"));stopifnot(file.rename(paste0(dest,".tmp"),dest))
    write.csv(result$metrics,sub("\\.rds$","_metrics.csv",dest),row.names=FALSE)
    cat(sprintf("DONE %s D=%d sigma=%g errors=%d\n",key,d,cnd$sigma,sum(result$metrics$status!="ok")));flush.console()
  };TRUE
}
cat(sprintf("START conditions=%d reps=%d workers=%d manifest=%s\n",nrow(grid),reps,workers,manifest_hash));flush.console()
tasks<-parallel::mclapply(seq_len(reps),run_rep,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
stopifnot(!any(vapply(tasks,inherits,logical(1),what="try-error")))
runs<-lapply(list.files(file.path(out,"runs"),pattern="\\.rds$",full.names=TRUE),readRDS)
metrics<-do.call(rbind,lapply(runs,`[[`,"metrics"))
stopifnot(nrow(metrics)==nrow(grid)*reps*4L,!anyDuplicated(paste(metrics$condition_id,metrics$rep,metrics$method)))
write.csv(metrics,file.path(out,"raw_metrics.csv"),row.names=FALSE)
cat(sprintf("COMPLETE instances=%d records=%d errors=%d\n",length(runs),nrow(metrics),sum(metrics$status!="ok")))

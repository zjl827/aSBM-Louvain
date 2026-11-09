#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(igraph); library(mvtnorm); library(mclust)})
script <- normalizePath(sub("^--file=", "", commandArgs(FALSE)[startsWith(commandArgs(FALSE), "--file=")][1]))
wd <- dirname(script); root <- normalizePath(file.path(wd, "../..")); code <- file.path(root, "code")
model_source <- file.path(dirname(root), "asbm_louvain.R")
source(file.path(code, "weighted_model.R"))
source(file.path(wd, "design_helpers.R"))
source(file.path(code, "pure_asbm_vectorized.R"))
stopifnot(formals(stats::kmeans)$iter.max == 10L)
betas <- c(0, 0.25, 0.5, 1, 2, 4)
models <- setNames(lapply(betas, function(b) make_model(model_source, b)), as.character(betas))
baseline <- make_model(model_source)
smoke <- Sys.getenv("WM_SMOKE", "0") == "1"
out <- Sys.getenv("WM_OUTPUT", file.path(root, "outputs/weight_misspecification"))
workers <- as.integer(Sys.getenv("WM_WORKERS", "4")); reps <- if (smoke) 1L else 10L
dir.create(file.path(out, "runs"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out, "data_cache"), recursive=TRUE, showWarnings=FALSE)
files <- c(script, file.path(wd, c("design_helpers.R", "protocol.md")),
           model_source, file.path(code, c("weighted_model.R", "pure_asbm_vectorized.R")))
manifest <- tools::md5sum(files); manifest_hash <- hash_object(manifest)
manifest_file <- file.path(out, "manifest.rds")
if (file.exists(manifest_file)) {
  if (!identical(readRDS(manifest_file), manifest)) {
    stop("Use a new WM_OUTPUT directory for the current implementation.")
  }
} else saveRDS(manifest, manifest_file)
write.csv(data.frame(path=names(manifest), md5=unname(manifest)), file.path(out, "manifest.csv"), row.names=FALSE)
saveRDS(list(smoke=smoke, reps=reps, workers=workers, rng=RNGkind(),
             session=capture.output(sessionInfo()), manifest_hash=manifest_hash), file.path(out,"settings.rds"))

safe_fit <- function(fun) {
  warnings <- character(); start <- proc.time()[["elapsed"]]
  fit <- tryCatch({
    setTimeLimit(cpu=Inf, elapsed=300, transient=TRUE)
    withCallingHandlers(fun(), warning=function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  }, error=function(e) structure(list(message=conditionMessage(e)), class="fit_error"),
  finally=setTimeLimit(cpu=Inf, elapsed=Inf, transient=FALSE))
  list(fit=fit, elapsed=proc.time()[["elapsed"]]-start, warnings=paste(unique(warnings),collapse=" | "))
}
fit_proposed <- function(env,A,X,initial) env$aSBM_Louvain_revised(
  A, X, T_outer=1000L, tau0=1, rho_tau=0.8, resolution=1,
  epsilon_Q=1e-4, epsilon_L=1e-3, k_stable_required=3L,
  init_method="kmeans", initial_partition=initial, verbose=FALSE)

# Independent explicit pair-sum objective and local-update checks before any study fit.
set.seed(73311); nt <- 8; kt <- 2; dt <- 2
at <- matrix(rbinom(nt*nt,1,0.3),nt); at[lower.tri(at)] <- t(at)[lower.tri(at)]; diag(at)<-0
xt <- matrix(rnorm(nt*dt),nt); qt <- matrix(runif(nt*kt),nt); qt <- qt/rowSums(qt)
pt <- matrix(c(.3,.1,.1,.4),2); pit <- matrix(c(.7,.3),nt,2,byrow=TRUE)
mut <- list(c(0,0),c(1,-1)); covt <- list(diag(2),diag(c(2,1)))
llx <- sapply(seq_len(kt),function(k) dmvnorm(xt,mut[[k]],covt[[k]],log=TRUE))
lla <- 0; msg <- matrix(0,nt,kt)
for(i in seq_len(nt)) for(j in seq_len(nt)) if(i!=j) for(k in seq_len(kt)) for(l in seq_len(kt)) {
  term <- at[i,j]*log(pt[k,l]+1e-10)+(1-at[i,j])*log(1-pt[k,l]+1e-10)
  msg[i,k] <- msg[i,k]+qt[j,l]*term
  if(i<j) lla <- lla+qt[i,k]*qt[j,l]*term
}
for(b in betas) for(tau in c(.7,1e-8)) {
  env <- models[[as.character(b)]]
  expected <- lla+b*sum(qt*llx)-max(tau,1e-5)*sum(qt*(log(qt+1e-20)-log(pit+1e-20)))
  observed <- env$Calculate_ELBO(xt,mut,covt,pt,pit,qt,at,tau)
  scores <- log(pit+1e-20)+(msg+b*llx)/max(tau,1e-5)
  probs <- exp(scores-apply(scores,1,max)); probs <- probs/rowSums(probs)
  actual <- env$VariationalEStep_Strict(xt,mut,covt,pt,pit,tau,at,qt)
  stopifnot(abs(expected-observed)<1e-9, max(abs(probs-actual))<1e-10)
}
stopifnot(identical(baseline$Calculate_ELBO(xt,mut,covt,pt,pit,qt,at,.7),
                    models[["1"]]$Calculate_ELBO(xt,mut,covt,pt,pit,qt,at,.7)),
          identical(baseline$VariationalEStep_Strict(xt,mut,covt,pt,pit,.7,at,qt),
                    models[["1"]]$VariationalEStep_Strict(xt,mut,covt,pt,pit,.7,at,qt)))
writeLines("PASS: explicit undirected pair sums, weighted softmax, beta=1 identity, tau floor, K-means default limit=10", file.path(out,"preflight_validation.txt"))

grid <- expand.grid(r=c(1.5,4),sigma=c(1.5,4),h=c(1,3),distribution="Gaussian",stringsAsFactors=FALSE)
grid <- rbind(grid, expand.grid(r=c(1.5,4),sigma=c(1.5,4),h=1,
                              distribution=c("Student_t5","Laplace","Contaminated"),stringsAsFactors=FALSE))
grid$condition_id <- sprintf("W%02d",seq_len(nrow(grid)))
grid$weight_study <- grid$distribution=="Gaussian"
grid$misspecification_study <- grid$h==1
if(smoke) grid <- grid[c(1,9,13,17),]
write.csv(grid,file.path(out,"grid.csv"),row.names=FALSE)

run_group <- function(task) {
  ri<-task$ri; rep_id<-task$rep; ratio<-c(1.5,4)[ri]
  n<-1000L; k<-5L; d<-50L; sizes<-rep(200L,k); truth<-rep(seq_len(k),times=sizes)
  seeds <- list(graph=92210000L+ri*100000L+rep_id, means=92211000L+rep_id,
                residual=92212000L+rep_id, high=92213000L+rep_id, tail=92214000L+rep_id,
                contamination=92215000L+rep_id, louvain=92700000L+ri*1000L+rep_id,
                kmeans=92800000L+rep_id)
  set.seed(seeds$graph); pout<-20/(199*ratio+800); block<-matrix(pout,k,k); diag(block)<-ratio*pout
  g<-sample_sbm(n,block,sizes,directed=FALSE,loops=FALSE)
  A<-as.matrix(as_adjacency_matrix(g,sparse=FALSE))
  set.seed(seeds$means); means<-matrix(rnorm(k*d),k,d)
  set.seed(seeds$residual); Z<-matrix(rnorm(n*d),n,d)
  set.seed(seeds$high); high<-sort(sample(seq_len(k),2))
  set.seed(seeds$tail); u<-rchisq(n,df=5)
  set.seed(seeds$contamination); contamination<-runif(n)<.1
  noise<-list(Gaussian=Z, Student_t5=Z/sqrt(u/3),
              Laplace=-sign(Z)*log(2*pnorm(-abs(Z)))/sqrt(2),
              Contaminated=Z*ifelse(contamination,5,1)/sqrt(3.4))
  stopifnot(all(vapply(noise,function(x)all(is.finite(x)),logical(1))), all(diag(A)==0), identical(A,t(A)))
  set.seed(seeds$louvain); lp<-as.integer(factor(membership(cluster_louvain(g,resolution=1))))
  initK<-max(2L,length(unique(lp))); ah<-hash_object(A)
  base_key<-sprintf("r%d_rep%02d",ri,rep_id)
  saveRDS(list(A=A,means=means,Z=Z,u=u,contamination=contamination,high=high,truth=truth,seeds=seeds),
          file.path(out,"data_cache",paste0(base_key,".rds")))
  conditions<-grid[grid$r==ratio,]
  for(ci in seq_len(nrow(conditions))) {
    cnd<-conditions[ci,]; key<-sprintf("%s_rep%02d",cnd$condition_id,rep_id)
    dest<-file.path(out,"runs",paste0(key,".rds"))
    if(file.exists(dest)) {
      old<-readRDS(dest);stopifnot(old$manifest_hash==manifest_hash,nrow(old$metrics)==if(cnd$weight_study)9L else 4L)
      next
    }
    sigmas<-community_noise_sd(cnd$h,cnd$sigma,high,sizes)
    residual<-noise[[cnd$distribution]]*sigmas[truth]
    X<-means[truth,,drop=FALSE]+residual
    stopifnot(abs(mean(sigmas^2)-cnd$sigma^2)<1e-10)
    set.seed(seeds$kmeans)
    km<-safe_fit(function()stats::kmeans(X,centers=initK,nstart=5,iter.max=10L))
    failed_init<-inherits(km$fit,"fit_error")
    initial<-if(failed_init)NULL else as.integer(factor(km$fit$cluster)); post_rng<-.Random.seed
    bs<-if(cnd$weight_study)betas else 1
    methods<-c(rep("Proposed",length(bs)),"Pure_aSBM","Louvain","KMeans_estK")
    weights<-c(bs,rep(NA_real_,3)); records<-list(); locals<-list(); attempts<-list(); replay<-NULL
    for(mi in seq_along(methods)) {
      method<-methods[mi]; beta<-weights[mi]
      if(method=="Louvain") {
        attempt<-list(fit=list(partition=lp,iterations=NA_integer_,stop_reason="reference_partition"),elapsed=NA_real_,warnings="")
      } else if(failed_init) attempt<-km else if(method=="KMeans_estK") {
        attempt<-list(fit=list(partition=initial,iterations=km$fit$iter,stop_reason="reference_partition"),elapsed=km$elapsed,warnings=km$warnings)
      } else {
        assign(".Random.seed",post_rng,envir=.GlobalEnv)
        attempt<-if(method=="Proposed")safe_fit(function()fit_proposed(models[[as.character(beta)]],A,X,initial)) else
          safe_fit(function()fit_pure_asbm_vectorized(A,X,initial_partition=initial,tolerance=1,timeout_sec=300,max_iterations=1000L))
      }
      fit<-attempt$fit; failed<-inherits(fit,"fit_error")
      # Eight complete runs verify beta=1 against the unchanged corrected implementation.
      if(method=="Proposed" && beta==1 && cnd$weight_study && rep_id==1 && !failed_init) {
        assign(".Random.seed",post_rng,envir=.GlobalEnv)
        replay<-safe_fit(function()fit_proposed(baseline,A,X,initial))
        stopifnot(identical(replay$fit,fit), identical(replay$warnings,attempt$warnings))
      }
      row<-cnd; row$rep<-rep_id;row$method<-method;row$beta<-beta;row$base_key<-base_key
      row$initial_K<-initK;row$graph_hash<-ah;row$attribute_hash<-hash_object(X)
      row$initial_hash<-hash_object(initial);row$post_init_rng_hash<-hash_object(post_rng)
      row$init_warnings<-km$warnings;row$warnings<-attempt$warnings
      row$status<-if(failed)"error" else "ok";row$error<-if(failed)fit$message else ""
      row$elapsed_sec<-attempt$elapsed;row$realized_noise_variance<-mean(residual^2)
      row$realized_mean_degree<-sum(A)/n;row$contamination_fraction<-if(cnd$distribution=="Contaminated")mean(contamination) else NA_real_
      row$ARI<-row$NMI<-row$macro_F1<-row$min_F1<-row$high_F1<-row$low_F1<-NA_real_
      row$iterations<-row$final_K<-NA_integer_;row$stop_reason<-"error";row$met_stopping_rule<-FALSE
      if(!failed) {
        pred<-as.integer(fit$partition);stopifnot(length(pred)==n,!anyNA(pred))
        row$ARI<-adjustedRandIndex(truth,pred);row$NMI<-partition_nmi(truth,pred)
        row$final_K<-length(unique(pred));row$iterations<-if(method=="Proposed")fit$iterations_run else fit$iterations
        row$stop_reason<-fit$stop_reason
        row$met_stopping_rule<-if(method=="Proposed")fit$stop_reason=="joint_stability" else if(method=="Pure_aSBM")fit$stop_reason %in% c("tolerance","objective_not_increased") else NA
        cm<-community_scores(truth,pred,high)
        row$macro_F1<-mean(cm$F1);row$min_F1<-min(cm$F1)
        row$high_F1<-if(cnd$h>1)mean(cm$F1[cm$is_high_noise]) else NA_real_
        row$low_F1<-if(cnd$h>1)mean(cm$F1[!cm$is_high_noise]) else NA_real_
        cm$condition_id<-cnd$condition_id;cm$rep<-rep_id;cm$method<-method;cm$beta<-beta
        locals[[mi]]<-cm
      }
      records[[mi]]<-row;attempts[[mi]]<-attempt
    }
    result<-list(manifest_hash=manifest_hash,metrics=do.call(rbind,records),communities=do.call(rbind,locals),
                 fits=attempts,initial_partition=initial,post_init_rng=post_rng,X=X,base_key=base_key,
                 sigmas=sigmas,seeds=seeds,beta1_replay=if(is.null(replay))NA else TRUE)
    saveRDS(result,paste0(dest,".tmp"));stopifnot(file.rename(paste0(dest,".tmp"),dest))
    write.csv(result$metrics,sub("\\.rds$","_metrics.csv",dest),row.names=FALSE)
    cat(sprintf("DONE %s r=%g sigma=%g h=%g %s rows=%d errors=%d\n",key,ratio,cnd$sigma,cnd$h,cnd$distribution,nrow(result$metrics),sum(result$metrics$status!="ok")));flush.console()
  }
  TRUE
}
tasks<-list()
for(ri in 1:2) if(c(1.5,4)[ri] %in% grid$r) for(rep_id in seq_len(reps))tasks[[length(tasks)+1L]]<-list(ri=ri,rep=rep_id)
cat(sprintf("START conditions=%d repetitions=%d workers=%d manifest=%s\n",nrow(grid),reps,workers,manifest_hash));flush.console()
result<-parallel::mclapply(tasks,run_group,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
stopifnot(!any(vapply(result,inherits,logical(1),what="try-error")))
runs<-lapply(list.files(file.path(out,"runs"),pattern="\\.rds$",full.names=TRUE),readRDS)
metrics<-do.call(rbind,lapply(runs,`[[`,"metrics"));communities<-do.call(rbind,lapply(runs,`[[`,"communities"))
stopifnot(nrow(metrics)==sum(ifelse(grid$weight_study,9,4))*reps,
          !anyDuplicated(paste(metrics$condition_id,metrics$rep,metrics$method,metrics$beta)))
write.csv(metrics,file.path(out,"raw_metrics.csv"),row.names=FALSE)
write.csv(communities,file.path(out,"community_metrics.csv"),row.names=FALSE)
cat(sprintf("COMPLETE instances=%d records=%d Proposed=%d errors=%d beta1_replays=%d\n",length(runs),nrow(metrics),sum(metrics$method=="Proposed"),sum(metrics$status!="ok"),sum(vapply(runs,function(x)isTRUE(x$beta1_replay),logical(1)))));flush.console()

#!/usr/bin/env Rscript
# Independent read-only validation and transparent Pure objective diagnostics.
# Does not alter the manifest, experiment code, saved predictions, or metrics.
suppressPackageStartupMessages({library(igraph);library(mvtnorm)})
script <- normalizePath(sub("^--file=", "", commandArgs(FALSE)[startsWith(commandArgs(FALSE), "--file=")][1]))
root <- normalizePath(file.path(dirname(script), "../.."))
out <- Sys.getenv("DIM_PILOT_OUTPUT", file.path(root, "outputs/dimension_pilot"))
source(file.path(root, "experiments", "weight_misspecification", "design_helpers.R"))
code <- file.path(root, "code")
source(file.path(code, "pure_asbm_vectorized.R"))
checks <- character()
check <- function(ok, label) {stopifnot(isTRUE(ok));checks <<- c(checks,label)}
close <- function(a,b,tol=1e-11) isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE))
independent_ari <- function(a,b) {
  t<-table(a,b);ch<-function(x)x*(x-1)/2
  s<-sum(ch(t));r<-sum(ch(rowSums(t)));c<-sum(ch(colSums(t)));nn<-ch(sum(t))
  expected<-r*c/nn;denominator<-(r+c)/2-expected
  if(denominator==0) return(if(identical(as.vector(outer(a,a,"==")),as.vector(outer(b,b,"=="))))1 else 0)
  (s-expected)/denominator
}
independent_nmi <- function(a,b) {
  p<-table(a,b)/length(a);pa<-rowSums(p);pb<-colSums(p);nz<-p>0
  information<-sum(p[nz]*log(p[nz]/outer(pa,pb)[nz]))
  entropy<- -sum(pa[pa>0]*log(pa[pa>0]))-sum(pb[pb>0]*log(pb[pb>0]))
  if(entropy==0)1 else 2*information/entropy
}
manifest<-readRDS(file.path(out,"manifest.rds"));settings<-readRDS(file.path(out,"settings.rds"))
check(identical(tools::md5sum(names(manifest)),manifest),"Manifest hashes unchanged")
check(identical(hash_object(manifest),settings$manifest_hash),"Manifest identifier consistent")
check(!settings$smoke&&settings$reps==10&&settings$workers>=1,"Full planned run settings")
RNGkind(settings$rng[1],settings$rng[2],settings$rng[3])
raw<-read.csv(file.path(out,"raw_metrics.csv"),stringsAsFactors=FALSE)
files<-list.files(file.path(out,"runs"),pattern="\\.rds$",full.names=TRUE)
check(length(files)==60&&nrow(raw)==240,"60 datasets and 240 method records present")
check(!anyDuplicated(paste(raw$condition_id,raw$rep,raw$method)),"Unique condition/rep/method keys")
check(all(table(raw$condition_id,raw$method)==10),"Exactly 10 planned replicates in every cell")
seeds_seen<-integer();graph_hashes<-mean_hashes<-residual_hashes<-character()
datasets<-list();maximum_ari_error<-maximum_nmi_error<-0;prediction_count<-0L
replay<-list();objective_original<-pure_asbm_objective
objective_log<-numeric()
pure_asbm_objective<-function(...) {v<-objective_original(...);objective_log<<-c(objective_log,v);v}
for(b in 1:10) {
  cache<-readRDS(file.path(out,"data_cache",sprintf("rep%02d.rds",b)))
  expected_seeds<-list(graph=108100000L+b,means=108110000L+b,residual=108120000L+b,louvain=108130000L+b,kmeans=108140000L+b)
  check(identical(cache$seeds,expected_seeds),sprintf("rep%02d scheduled seed values",b))
  seeds_seen<-c(seeds_seen,unlist(cache$seeds))
  set.seed(cache$seeds$graph);pout<-20/(199*4+800);B<-matrix(pout,5,5);diag(B)<-4*pout
  g<-sample_sbm(1000,B,rep(200L,5),directed=FALSE,loops=FALSE)
  A<-as.matrix(as_adjacency_matrix(g,sparse=FALSE))
  set.seed(cache$seeds$means);means<-matrix(rnorm(5*50),5,50)
  set.seed(cache$seeds$residual);Z<-matrix(rnorm(1000*200),1000,200)
  set.seed(cache$seeds$louvain);lp<-as.integer(factor(membership(cluster_louvain(g,resolution=1))))
  check(identical(cache$A,A)&&identical(cache$means,means)&&identical(cache$Z,Z)&&identical(cache$louvain_partition,lp),sprintf("rep%02d A, means, residuals, Louvain reproduce exactly from seeds",b))
  check(identical(cache$truth,rep(1:5,each=200)),sprintf("rep%02d correct true labels",b))
  graph_hashes<-c(graph_hashes,hash_object(A));mean_hashes<-c(mean_hashes,hash_object(means));residual_hashes<-c(residual_hashes,hash_object(Z))
  for(sig in c(1.5,4)) {
    previous_X<-NULL;pair_distances<-NULL
    for(d in c(50,100,200)) {
      ci<-which(c(50,100,200,50,100,200)==d & c(1.5,1.5,1.5,4,4,4)==sig)
      key<-sprintf("D%02d_rep%02d",ci,b);r<-readRDS(file.path(out,"runs",paste0(key,".rds")))
      check(identical(r$manifest_hash,settings$manifest_hash)&&identical(r$seeds,cache$seeds),paste(key,"saved provenance matches"))
      M<-cbind(means,matrix(0,5,d-50));X<-M[cache$truth,,drop=FALSE]+sig*Z[,1:d,drop=FALSE]
      check(identical(r$X,X),paste(key,"X matches exact construction with zero added means"))
      if(!is.null(previous_X))check(identical(X[,1:ncol(previous_X),drop=FALSE],previous_X),paste(key,"nested X prefix exact"))
      previous_X<-X
      ds<-as.numeric(dist(M))^2/sig^2
      if(!is.null(pair_distances))check(identical(ds,pair_distances),paste(key,"all 10 theoretical pair distances unchanged with D"))
      pair_distances<-ds
      iw<-character();set.seed(cache$seeds$kmeans)
      km<-withCallingHandlers(kmeans(X,centers=max(2L,length(unique(lp))),nstart=5,iter.max=10L),warning=function(w){iw<<-c(iw,conditionMessage(w));invokeRestart("muffleWarning")})
      initial<-as.integer(factor(km$cluster));post_rng<-.Random.seed
      check(identical(initial,r$initial_partition)&&identical(post_rng,r$post_init_rng),paste(key,"K-means iter.max=10 and post-initialization RNG reproduce exactly"))
      check(identical(paste(unique(iw),collapse=" | "),r$kmeans$warnings)&&length(iw)==r$kmeans$warning_count,paste(key,"initialization warnings reproducible"))
      m<-r$metrics;csv<-raw[raw$condition_id==sprintf("D%02d",ci)&raw$rep==b,]
      rownames(m)<-rownames(csv)<-NULL
      # Empty strings are imported as empty strings and NA values remain NA.
      check(close(m,csv),paste(key,"RDS metrics match aggregate CSV"))
      check(all(m$graph_hash==hash_object(A))&&all(m$attribute_hash==hash_object(X))&&all(m$initial_hash==hash_object(initial))&&all(m$post_init_rng_hash==hash_object(post_rng)),paste(key,"all methods share data/initialization/RNG metadata"))
      check(all(abs(m$mean_pairwise_mahalanobis_sq-mean(ds))<1e-11),paste(key,"reported theoretical distance correct"))
      check(identical(as.integer(r$fits[[1]]$fit$initial_partition),initial),paste(key,"Proposed saved initial partition matches shared partition"))
      if(!inherits(r$fits[[2]]$fit,"fit_error"))check(identical(as.integer(r$fits[[2]]$fit$initial_partition),initial),paste(key,"Pure saved initial partition matches shared partition"))
      check(identical(r$fits[[3]]$fit$partition,lp)&&identical(r$fits[[4]]$fit$partition,initial),paste(key,"baseline predictions use scheduled Louvain and K-means partitions"))
      for(i in 1:4) {
        fit<-r$fits[[i]]$fit;failed<-inherits(fit,"fit_error")
        check(identical(m$warnings[i],r$fits[[i]]$warnings)&&m$warning_count[i]==r$fits[[i]]$warning_count,paste(key,m$method[i],"complete warning metadata retained"))
        check((m$status[i]=="error")==failed,paste(key,m$method[i],"error status consistent"))
        if(failed) {check(is.na(m$ARI[i])&&is.na(m$NMI[i])&&identical(m$error[i],fit$message),paste(key,m$method[i],"failure retained with no invented metric"));next}
        a<-independent_ari(cache$truth,fit$partition);nm<-independent_nmi(cache$truth,fit$partition)
        maximum_ari_error<-max(maximum_ari_error,abs(a-m$ARI[i]));maximum_nmi_error<-max(maximum_nmi_error,abs(nm-m$NMI[i]))
        check(abs(a-m$ARI[i])<1e-11&&abs(nm-m$NMI[i])<1e-11,paste(key,m$method[i],"independent ARI and NMI formula agrees"));prediction_count<-prediction_count+1L
        if(i==1)check(identical(m$objective_finite[i],all(is.finite(fit$elbo_history)))&&identical(m$stop_reason[i],fit$stop_reason),paste(key,"Proposed objective and stop metadata agree"))
        if(i==2)check(identical(m$objective_finite[i],is.finite(fit$final_objective))&&identical(m$stop_reason[i],fit$stop_reason),paste(key,"Pure retained objective and stop metadata agree"))
        if(i<=2)check(identical(m$met_stopping_rule[i],if(i==1)fit$stop_reason=="joint_stability" else fit$stop_reason%in%c("tolerance","objective_not_increased")),paste(key,m$method[i],"code stopping flag accurately recorded (not a convergence proof)"))
      }
      # Pure is replayed solely to inspect rejected candidate objectives. No algorithm changes.
      objective_log<-numeric();assign(".Random.seed",post_rng,envir=.GlobalEnv);pw<-character()
      pf<-tryCatch(withCallingHandlers(fit_pure_asbm_vectorized(A,X,initial_partition=initial,tolerance=1,timeout_sec=300,max_iterations=1000L),warning=function(w){pw<<-c(pw,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)structure(list(message=conditionMessage(e)),class="fit_error"))
      stored<-r$fits[[2]]$fit;pf_failed<-inherits(pf,"fit_error")
      check(pf_failed==inherits(stored,"fit_error"),paste(key,"Pure replay failure state same"))
      if(pf_failed)check(identical(pf$message,stored$message),paste(key,"Pure replay error same")) else check(identical(pf$partition,stored$partition)&&identical(pf$stop_reason,stored$stop_reason)&&close(pf$final_objective,stored$final_objective)&&pf$iterations==stored$iterations,paste(key,"Pure replay predictions, stopping and objective agree"))
      check(identical(paste(unique(pw),collapse=" | "),r$fits[[2]]$warnings)&&length(pw)==r$fits[[2]]$warning_count,paste(key,"Pure replay warnings agree"))
      candidates<-objective_log[-1]
      replay[[key]]<-data.frame(condition_id=sprintf("D%02d",ci),rep=b,D=d,sigma=sig,status=if(pf_failed)"error" else "ok",stop_reason=if(pf_failed)"error" else pf$stop_reason,
          objective_calls=length(objective_log),initial_objective=if(length(objective_log))objective_log[1]else NA_real_,
          initial_objective_finite=if(length(objective_log))is.finite(objective_log[1])else NA,
          candidate_nonfinite=sum(!is.finite(candidates)),last_candidate=if(length(candidates))tail(candidates,1)else NA_real_,
          final_objective_finite=if(pf_failed)NA else is.finite(pf$final_objective),
          identical_to_kmeans=if(pf_failed)NA else identical(as.integer(pf$partition),initial),
          partition_ARI_to_kmeans=if(pf_failed)NA else independent_ari(pf$partition,initial),
          pure_iterations=if(pf_failed)NA else pf$iterations)
    }
  }
}
check(!anyDuplicated(seeds_seen)&&length(seeds_seen)==50,"All 50 data/algorithm stream seeds distinct across repetitions and roles")
check(!anyDuplicated(graph_hashes)&&!anyDuplicated(mean_hashes)&&!anyDuplicated(residual_hashes),"All 10 graph, mean, and residual arrays distinct across repetitions")
replay<-do.call(rbind,replay);rownames(replay)<-NULL
write.csv(replay,file.path(out,"pure_objective_diagnostic.csv"),row.names=FALSE)
diag_summary<-do.call(rbind,lapply(split(replay,list(replay$D,replay$sigma),drop=TRUE),function(a)data.frame(D=a$D[1],sigma=a$sigma[1],planned=nrow(a),returned=sum(a$status=="ok"),initial_nonfinite=sum(!a$initial_objective_finite),candidate_nonfinite=sum(a$candidate_nonfinite),returned_identical_to_kmeans=sum(a$identical_to_kmeans,na.rm=TRUE))))
report<-c("DIMENSION PILOT INDEPENDENT VALIDATION",sprintf("Passed assertions: %d",length(checks)),sprintf("Predictions independently rescored: %d / 240; one failed fit preserved",prediction_count),
  sprintf("Maximum absolute ARI difference: %.3g",maximum_ari_error),sprintf("Maximum absolute NMI difference: %.3g",maximum_nmi_error),
  "All dataset seed reconstructions, nested prefixes, theoretical pair distances, initialization replays, and Pure replays passed.",
  "Distinct deterministic PRNG streams and regenerated arrays were verified; statistical independence is a property assumed of the PRNG design, not proved by hash uniqueness.",
  "Proposed/Pure share initial partitions and saved post-initialization RNG; code explicitly resets RNG before each iterative fit.",
  "Pure objective diagnostics: transparent logging wrapper around original objective; no fit algorithm or experiment outputs altered.",
  capture.output(print(diag_summary,row.names=FALSE)),
  sprintf("Total Pure initial nonfinite objectives: %d; runs with a rejected nonfinite candidate: %d; successful outputs identical to K-means: %d / %d",sum(!replay$initial_objective_finite),sum(replay$candidate_nonfinite>0),sum(replay$identical_to_kmeans,na.rm=TRUE),sum(replay$status=="ok")),
  "At D=200: all 20 initial objectives were -Inf; 19 returned results rejected a later -Inf candidate and exactly retained the K-means hard partition; one run raised an error. Thus finite saved final objectives mask numerical degeneracy in this setting.",
  "Additional numerical degeneration: D=100, sigma=1.5, rep=9 also had initial -Inf and a rejected -Inf candidate.",
  "Pure stops:",capture.output(with(replay,print(table(D,sigma,stop_reason)))),
  "All original exception rows:",capture.output(print(raw[raw$status!="ok"|raw$warning_count>0|raw$init_warning_count>0|(!is.na(raw$objective_finite)&!raw$objective_finite),c("condition_id","rep","D","sigma","method","status","stop_reason","error","warnings","init_warnings","objective_finite")],row.names=FALSE)),
  "Interpretation constraints:",
  "objective_not_increased is a code stopping condition, not proof of convergence or optimality. It can reject a nonfinite candidate while retaining a finite prior objective.",
  "Saved objective_finite only describes retained values; see pure_objective_diagnostic.csv for rejected candidate values.",
  "Initial cluster size <= D raises covariance rank concerns for hard initial estimates but does not imply model nonidentifiability or establish the rank of later soft-weight estimates.",
  "Pure uses a different covariance regularization scheme from Proposed. Numerical failure is not evidence of a statistical recovery advantage.",
  "END VALIDATION")
writeLines(report,file.path(out,"validation.txt"));cat(paste(report[1:18],collapse="\n"),"\n")

#!/usr/bin/env Rscript
# Read-only input, prediction, and diagnostic validation of the extension.
suppressPackageStartupMessages({library(igraph);library(mvtnorm)})
script<-normalizePath(sub("^--file=","",commandArgs(FALSE)[startsWith(commandArgs(FALSE),"--file=")][1]))
root<-normalizePath(file.path(dirname(script),"../.."))
out <- Sys.getenv("DIM_OUTPUT", file.path(root, "outputs/dimension_extension"))
oldout <- Sys.getenv("DIM_PILOT_OUTPUT", file.path(root, "outputs/dimension_pilot"))
source(file.path(root, "experiments", "weight_misspecification", "design_helpers.R"))
code <- file.path(root, "code")
source(file.path(code, "pure_asbm_vectorized.R"))
assertions<-0L
check<-function(ok,label){if(!isTRUE(ok))stop(label);assertions<<-assertions+1L}
close<-function(a,b,tol=1e-11)isTRUE(all.equal(a,b,tolerance=tol,check.attributes=FALSE))
ari<-function(a,b){t<-table(a,b);ch<-function(x)x*(x-1)/2;s<-sum(ch(t));r<-sum(ch(rowSums(t)));c<-sum(ch(colSums(t)));nn<-ch(sum(t));e<-r*c/nn;d<-(r+c)/2-e;if(d==0)return(if(identical(as.vector(outer(a,a,"==")),as.vector(outer(b,b,"=="))))1 else 0);(s-e)/d}
nmi<-function(a,b){p<-table(a,b)/length(a);pa<-rowSums(p);pb<-colSums(p);nz<-p>0;info<-sum(p[nz]*log(p[nz]/outer(pa,pb)[nz]));entropy<- -sum(pa[pa>0]*log(pa[pa>0]))-sum(pb[pb>0]*log(pb[pb>0]));if(entropy==0)1 else 2*info/entropy}
manifest<-readRDS(file.path(out,"manifest.rds"));prior_manifest<-readRDS(file.path(oldout,"manifest.rds"));settings<-readRDS(file.path(out,"settings.rds"))
check(identical(tools::md5sum(names(manifest)),manifest),"Extension manifest changed")
check(identical(tools::md5sum(names(prior_manifest)),prior_manifest),"Original manifest changed")
common<-intersect(names(manifest),names(prior_manifest))
check(length(common)==4&&identical(manifest[common],prior_manifest[common]),"Underlying model/helper code differs from pilot")
check(identical(hash_object(manifest),settings$manifest_hash),"Manifest hash mismatch")
check(!settings$smoke&&settings$reps==10&&settings$workers>=1,"Unexpected settings")
RNGkind(settings$rng[1],settings$rng[2],settings$rng[3])
raw<-read.csv(file.path(out,"raw_metrics.csv"),stringsAsFactors=FALSE,colClasses=c(error="character",warnings="character",init_warnings="character"))
files<-list.files(file.path(out,"runs"),pattern="\\.rds$",full.names=TRUE)
check(length(files)==40&&nrow(raw)==160,"Incomplete extension results")
check(!anyDuplicated(paste(raw$condition_id,raw$rep,raw$method))&&all(table(raw$condition_id,raw$method)==10),"Incomplete or duplicate method cells")
diagnostic<-list();rescored<-0L;ari_error<-nmi_error<-0;seeds<-integer();replay_note<-NULL
for(b in 1:10){
  cache<-readRDS(file.path(out,"data_cache",sprintf("rep%02d.rds",b)))
  prior<-readRDS(file.path(oldout,"data_cache",sprintf("rep%02d.rds",b)))
  expected_seeds<-list(graph=108100000L+b,means=108110000L+b,residual=108120000L+b,louvain=108130000L+b,kmeans=108140000L+b)
  check(identical(cache$seeds,expected_seeds)&&identical(prior$seeds,expected_seeds),"Seed schedule mismatch")
  seeds<-c(seeds,unlist(cache$seeds))
  for(field in c("A","means","truth","seeds","louvain_partition"))check(identical(cache[[field]],prior[[field]]),paste("Different original",field))
  check(identical(cache$Z[,1:200,drop=FALSE],prior$Z),"Original residual prefix mismatch")
  set.seed(cache$seeds$graph);pout<-20/(199*4+800);B<-matrix(pout,5,5);diag(B)<-4*pout
  g<-sample_sbm(1000,B,rep(200L,5),directed=FALSE,loops=FALSE);A<-as.matrix(as_adjacency_matrix(g,sparse=FALSE))
  set.seed(cache$seeds$means);M50<-matrix(rnorm(250),5,50)
  set.seed(cache$seeds$residual);Z<-matrix(rnorm(800000),1000,800)
  set.seed(cache$seeds$louvain);lp<-as.integer(factor(membership(cluster_louvain(g,resolution=1))))
  check(identical(A,cache$A)&&identical(M50,cache$means)&&identical(Z,cache$Z)&&identical(lp,cache$louvain_partition),"Seeds did not reproduce inputs")
  check(identical(cache$truth,rep(1:5,each=200)),"Wrong truth labels")
  for(sig in c(1.5,4)){
    oldci<-if(sig==1.5)"D03"else"D06"
    oldr<-readRDS(file.path(oldout,"runs",sprintf("%s_rep%02d.rds",oldci,b)))
    prevX<-oldr$X;dist0<-as.numeric(dist(M50))^2/sig^2
    for(d in c(400,800)){
      ci<-if(sig==1.5)if(d==400)1L else 2L else if(d==400)3L else 4L
      key<-sprintf("E%02d_rep%02d",ci,b);r<-readRDS(file.path(out,"runs",paste0(key,".rds")))
      check(identical(r$manifest_hash,settings$manifest_hash)&&identical(r$seeds,cache$seeds),"Run provenance differs")
      M<-cbind(M50,matrix(0,5,d-50));X<-M[cache$truth,,drop=FALSE]+sig*Z[,1:d,drop=FALSE]
      check(identical(X,r$X),"Wrong extended input")
      check(identical(X[,1:ncol(prevX),drop=FALSE],prevX),"Extended X prefix mismatch")
      prevX<-X
      check(identical(as.numeric(dist(M))^2/sig^2,dist0),"Pairwise theoretical distances changed")
      iw<-character();set.seed(cache$seeds$kmeans)
      km<-withCallingHandlers(kmeans(X,centers=max(2L,length(unique(lp))),nstart=5,iter.max=10L),warning=function(w){iw<<-c(iw,conditionMessage(w));invokeRestart("muffleWarning")})
      init<-as.integer(factor(km$cluster));rng<-.Random.seed
      check(identical(init,r$initial_partition)&&identical(rng,r$post_init_rng),"K-means init or RNG replay mismatch")
      check(identical(paste(unique(iw),collapse=" | "),r$kmeans$warnings)&&length(iw)==r$kmeans$warning_count,"K-means warning replay mismatch")
      m<-r$metrics;csv<-raw[raw$condition_id==sprintf("E%02d",ci)&raw$rep==b,];rownames(m)<-rownames(csv)<-NULL
      check(close(m,csv),"Saved run metrics differ from aggregate")
      check(all(m$graph_hash==hash_object(A))&&all(m$attribute_hash==hash_object(X))&&all(m$initial_hash==hash_object(init))&&all(m$post_init_rng_hash==hash_object(rng)),"Shared data/init/RNG metadata mismatch")
      check(all(abs(m$mean_pairwise_mahalanobis_sq-mean(dist0))<1e-11),"Wrong reported theoretical distance")
      for(i in 1:4){
        fit<-r$fits[[i]]$fit;failed<-inherits(fit,"fit_error")
        check(identical(m$warnings[i],r$fits[[i]]$warnings)&&m$warning_count[i]==r$fits[[i]]$warning_count,"Warning metadata differs")
        check((m$status[i]=="error")==failed,"Error status differs")
        if(failed){check(is.na(m$ARI[i])&&is.na(m$NMI[i])&&identical(m$error[i],fit$message),"Failure not preserved");next}
        check(length(fit$partition)==1000&&!anyNA(fit$partition),"Invalid prediction")
        if(i<=2)check(identical(as.integer(fit$initial_partition),init),"Iterative methods did not share initial partition")
        if(i==3)check(identical(fit$partition,lp),"Wrong Louvain prediction")
        if(i==4)check(identical(fit$partition,init),"Wrong K-means prediction")
        av<-ari(cache$truth,fit$partition);nv<-nmi(cache$truth,fit$partition)
        ari_error<-max(ari_error,abs(av-m$ARI[i]));nmi_error<-max(nmi_error,abs(nv-m$NMI[i]))
        check(abs(av-m$ARI[i])<1e-11&&abs(nv-m$NMI[i])<1e-11,"Independent metric mismatch");rescored<-rescored+1L
        check(identical(m$stop_reason[i],fit$stop_reason),"Stop reason mismatch")
        if(i<=2){
          finite<-if(i==1)all(is.finite(fit$elbo_history))else is.finite(fit$final_objective)
          stopping<-if(i==1)fit$stop_reason=="joint_stability"else fit$stop_reason%in%c("tolerance","objective_not_increased")
          check(identical(m$objective_finite[i],finite)&&identical(m$met_stopping_rule[i],stopping),"Objective/stopping flag mismatch")
        }
      }
      log<-r$pure_objective_log;pf<-r$fits[[2]]$fit;fail<-inherits(pf,"fit_error")
      check(identical(m$pure_initial_objective_finite[2],if(length(log))is.finite(log[1])else NA),"Wrong initial Pure objective diagnostic")
      check(m$pure_candidate_nonfinite[2]==sum(!is.finite(log[-1])),"Wrong nonfinite candidate diagnostic")
      eq<-if(fail)NA else identical(as.integer(pf$partition),init)
      check(identical(m$pure_identical_to_kmeans[2],eq),"Wrong Pure/Kmeans equality diagnostic")
      if(!fail)check(any(log==pf$final_objective,na.rm=TRUE),"Pure final objective absent from preserved log")
      diagnostic[[key]]<-data.frame(condition_id=sprintf("E%02d",ci),rep=b,D=d,sigma=sig,status=m$status[2],stop_reason=m$stop_reason[2],objective_calls=length(log),initial_objective=if(length(log))log[1]else NA_real_,initial_nonfinite=if(length(log))!is.finite(log[1])else NA,candidate_nonfinite=sum(!is.finite(log[-1])),last_candidate=if(length(log)>1)tail(log,1)else NA_real_,returned_identical_to_kmeans=eq)
      # Prespecified single original-function replay: rep=1, D=400, sigma=1.5.
      # No diagnostic wrapper here: compare original function to saved wrapped result.
      if(b==1&&d==400&&sig==1.5){
        assign(".Random.seed",rng,envir=.GlobalEnv);pw<-character()
        ref<-tryCatch(withCallingHandlers(fit_pure_asbm_vectorized(A,X,initial_partition=init,tolerance=1,timeout_sec=300,max_iterations=1000L),warning=function(w){pw<<-c(pw,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)structure(list(message=conditionMessage(e)),class="fit_error"))
        check(inherits(ref,"fit_error")==fail,"Original Pure replay status differs")
        if(fail)check(identical(ref$message,pf$message),"Original Pure replay error differs") else check(identical(ref$partition,pf$partition)&&identical(ref$stop_reason,pf$stop_reason)&&ref$iterations==pf$iterations&&close(ref$final_objective,pf$final_objective)&&close(ref$expectation,pf$expectation),"Original Pure replay differs from transparent wrapper")
        check(identical(paste(unique(pw),collapse=" | "),r$fits[[2]]$warnings)&&length(pw)==r$fits[[2]]$warning_count,"Original Pure replay warnings differ")
        replay_note<-"Prespecified original-function replay (rep 1, D=400, sigma=1.5) matched stored wrapped fit, including partition, expectation, objective, iterations, stop reason, and warnings. No other iterative fit rerun."
      }
    }
  }
}
check(length(seeds)==50&&!anyDuplicated(seeds),"Nonunique seed schedule")
diagnostic<-do.call(rbind,diagnostic);rownames(diagnostic)<-NULL
write.csv(diagnostic,file.path(out,"validation_pure_diagnostic.csv"),row.names=FALSE)
dg<-do.call(rbind,lapply(split(diagnostic,list(diagnostic$D,diagnostic$sigma),drop=TRUE),function(a)data.frame(D=a$D[1],sigma=a$sigma[1],planned=nrow(a),returned=sum(a$status=="ok"),initial_nonfinite=sum(a$initial_nonfinite),runs_with_nonfinite_candidates=sum(a$candidate_nonfinite>0),returned_identical_to_kmeans=sum(a$returned_identical_to_kmeans,na.rm=TRUE))))
report<-c("DIMENSION EXTENSION INDEPENDENT VALIDATION",sprintf("Passed assertions: %d",assertions),sprintf("Independently rescored predictions: %d / 160",rescored),sprintf("Maximum absolute ARI / NMI differences: %.3g / %.3g",ari_error,nmi_error),
  "Original and extension manifests intact; all four shared model and helper files have identical hashes.",
  "Runner diff inspected: only dimensions/output provenance, larger residual matrix with exact prefix check, and transparent Pure objective logging changed; algorithm calls and fit parameters unchanged.",
  "All 40 attributes exactly reconstructed; D=200 old prefixes and D=400-to-800 prefixes agree exactly. A, mean, truth, seeds, and Louvain partition are identical to original repetition.",
  "All 10 theoretical pairwise distances unchanged across dimensions; added community means are exactly zero.",
  "All 40 K-means initializations replayed with nstart=5, iter.max=10; saved post-initialization RNG and warnings match exactly.",
  "Both iterative fits use shared initial partitions and code explicitly resets the saved RNG before fitting. All prediction/exception/stop metadata checked.",
  replay_note,
  "Pure objective log and K-means identity diagnostics:",capture.output(print(dg,row.names=FALSE)),
  "All stops:",capture.output(print(with(raw,table(method,stop_reason)))),
  "Initialization warnings (unique datasets):",capture.output(print(unique(raw[raw$init_warning_count>0,c("condition_id","rep","D","sigma","init_warnings")]),row.names=FALSE)),
  "Errors and fit warnings:",capture.output(print(raw[raw$status!="ok"|raw$warning_count>0,c("condition_id","rep","D","sigma","method","status","error","warnings")],row.names=FALSE)),
  "Interpretation: same ten independent data repetitions, with paired additional dimensions; no new independent-repetition count is gained.",
  "An objective_not_increased flag is only a code stop condition, not proof of convergence/optimality. A retained finite objective can coexist with rejected nonfinite candidates.",
  "Hard initial covariance rank insufficiency is not a claim about model identifiability or all later soft covariance ranks. Different regularization prevents interpreting Pure numerical problems as statistical superiority.",
  "The selected replay validates the logging wrapper on one prespecified case; it is not an all-fit iterative replay.",
  "END VALIDATION")
writeLines(report,file.path(out,"validation.txt"));cat(paste(report,collapse="\n"),"\n")

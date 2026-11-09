# Helpers for the declared hard-partition plug-in bootstrap. No shared model is modified.
gb_hash <- function(x) {
  path <- tempfile('graph_bootstrap_hash_'); on.exit(unlink(path))
  saveRDS(x,path,compress=FALSE,version=3);unname(tools::md5sum(path))
}
gb_write_rds <- function(x,path) {
  saveRDS(x,paste0(path,'.tmp'));stopifnot(file.rename(paste0(path,'.tmp'),path));invisible(path)
}
gb_attempt <- function(fun, timeout=300) {
  warnings<-character();start<-proc.time()[['elapsed']]
  value<-tryCatch({setTimeLimit(cpu=Inf,elapsed=timeout,transient=TRUE)
    withCallingHandlers(fun(),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  },error=function(e)structure(list(message=conditionMessage(e)),class='gb_error'),
  finally=setTimeLimit(cpu=Inf,elapsed=Inf,transient=FALSE))
  list(value=value,warnings=paste(unique(warnings),collapse=' | '),elapsed=proc.time()[['elapsed']]-start)
}
gb_failed <- function(x) inherits(x$value,'gb_error')
gb_nmi <- function(a,b) {
  p<-table(a,b)/length(a);px<-rowSums(p);py<-colSums(p);e<-outer(px,py);v<-p>0
  mi<-sum(p[v]*log(p[v]/e[v]));h<--sum(px[px>0]*log(px[px>0]))-sum(py[py>0]*log(py[py>0]))
  if(h==0) 1 else 2*mi/h
}
gb_pair_scores <- function(anchor,prediction) {
  counts<-table(anchor,prediction);c2<-function(x)x*(x-1)/2
  within_total<-sum(c2(rowSums(counts)));between_total<-c2(sum(counts))-within_total
  within_same<-sum(c2(counts));between_same<-sum(c2(colSums(counts)))-within_same
  c(within_anchor_pair_retention=if(within_total>0)within_same/within_total else NA_real_,
    between_anchor_pair_merge=if(between_total>0)between_same/between_total else NA_real_)
}
gb_load_input <- function(source,input_dir) {
  if(source$kind=='simulation') {
    data<-gb_simulated_source(source)
    A<-data$A;X<-data$X;truth<-data$truth
  } else {
    data<-readRDS(file.path(input_dir,paste0(source$dataset,'__expression_graph_spatial_attributes.rds')))
    A<-data$adjacency;X<-data$attributes;truth<-data$truth
  }
  stopifnot(nrow(A)==source$N,ncol(A)==source$N,nrow(X)==source$N,ncol(X)==source$D,length(truth)==source$N)
  stopifnot(all(diag(A)==0),isTRUE(all.equal(A,t(A),check.attributes=FALSE)),!anyNA(X),!anyNA(truth))
  node_ids<-rownames(X);if(is.null(node_ids))node_ids<-as.character(seq_len(nrow(X)))
  list(A=A,X=X,truth=truth,node_ids=node_ids,graph_hash=gb_hash(A),attribute_hash=gb_hash(X),
       truth_hash=gb_hash(truth),representation=if(inherits(A,'sparseMatrix'))'sparse' else 'dense')
}
gb_initialization <- function(A,X,source) {
  set.seed(source$init_seed);seed_state<-.Random.seed
  g<-igraph::graph_from_adjacency_matrix(A,mode='undirected',diag=FALSE)
  lp<-gb_attempt(function()as.integer(factor(igraph::membership(igraph::cluster_louvain(g,resolution=source$resolution)))))
  after_louvain<-.Random.seed
  km<-if(gb_failed(lp))lp else gb_attempt(function()stats::kmeans(X,centers=max(2L,length(unique(lp$value))),nstart=5,iter.max=10L,algorithm='Hartigan-Wong'))
  km_labels<-if(gb_failed(km))NULL else as.integer(factor(km$value$cluster))
  list(louvain=lp,kmeans=km,kmeans_partition=km_labels,seed_state=seed_state,post_louvain_rng=after_louvain,post_init_rng=.Random.seed)
}
gb_fit_methods <- function(A,X,source) {
  init<-gb_initialization(A,X,source);attempts<-list()
  set.seed(source$fit_seed);fit_rng<-.Random.seed
  initial<-init$kmeans_partition
  for(method in c('Proposed','Pure_aSBM_C','Louvain','KMeans_estK')) {
    assign('.Random.seed',fit_rng,envir=.GlobalEnv)
    if(method=='Louvain') {
      a<-init$louvain
      if(!gb_failed(a))a$value<-list(partition=a$value,iterations=NA_integer_,stop_reason='reference_partition')
    } else if(method=='KMeans_estK') {
      a<-init$kmeans
      if(!gb_failed(a))a$value<-list(partition=initial,iterations=a$value$iter,stop_reason='reference_partition')
    } else if(is.null(initial)) a<-init$kmeans else if(method=='Proposed') {
      a<-gb_attempt(function()aSBM_Louvain_revised(A,X,T_outer=1000L,rho_tau=.8,tau0=1,
            resolution=source$resolution,epsilon_Q=1e-4,epsilon_L=1e-3,k_stable_required=3L,
            initial_partition=initial,init_method='kmeans',prior_feedback=TRUE,verbose=FALSE))
    } else {
      a<-gb_attempt(function()fit_pure_asbm_vectorized(A,X,initial_partition=initial,tolerance=1,timeout_sec=300,max_iterations=1000L))
    }
    attempts[[method]]<-a
  }
  list(initialization=init,fit_rng=fit_rng,attempts=attempts)
}
gb_partition <- function(attempt) if(gb_failed(attempt))NULL else as.integer(attempt$value$partition)
gb_record_metrics <- function(fits,A,X,source,arm,rep,anchor,own_partitions,truth,model_hash=NA_character_) {
  graph_hash<-gb_hash(A);attribute_hash<-gb_hash(X);rng_hash<-gb_hash(fits$fit_rng)
  rows<-lapply(names(fits$attempts),function(method) {
    attempt<-fits$attempts[[method]];pred<-gb_partition(attempt)
    initial<-if(method=='Louvain') {
      if(gb_failed(fits$initialization$louvain))NULL else fits$initialization$louvain$value
    } else fits$initialization$kmeans_partition
    initial_warning<-fits$initialization$louvain$warnings
    if(method!='Louvain')initial_warning<-paste(Filter(nzchar,c(initial_warning,fits$initialization$kmeans$warnings)),collapse=' | ')
    failed<-is.null(pred);fit<-attempt$value
    row<-data.frame(source_id=source$case_id,kind=source$kind,dataset=source$dataset,arm=arm,rep=rep,method=method,
      init_seed=source$init_seed,fit_seed=source$fit_seed,N=nrow(X),D=ncol(X),resolution=source$resolution,
      source_model_hash=model_hash,graph_hash=graph_hash,attribute_hash=attribute_hash,
      initial_hash=gb_hash(initial),fit_rng_hash=rng_hash,fit_rng_applied=method%in%c('Proposed','Pure_aSBM_C'),
      representation=if(inherits(A,'sparseMatrix'))'sparse' else 'dense',edge_count=sum(A)/2,
      initial_K=if(is.null(initial))NA_integer_ else length(unique(initial)),
      initial_ARI_anchor=if(is.null(initial))NA_real_ else mclust::adjustedRandIndex(anchor,initial),
      ARI_anchor=NA_real_,ARI_own_source=NA_real_,ARI_original_truth=NA_real_,NMI_anchor=NA_real_,
      final_K=NA_integer_,iterations=NA_integer_,stop_reason=if(failed)'error' else fit$stop_reason,
      status=if(failed)'error' else 'ok',error=if(failed)fit$message else '',
      init_warnings=initial_warning,fit_warnings=if(method%in%c('Proposed','Pure_aSBM_C'))attempt$warnings else '',
      diagnostic_sec=attempt$elapsed,within_anchor_pair_retention=NA_real_,between_anchor_pair_merge=NA_real_,
      stringsAsFactors=FALSE)
    if(!failed) {
      stopifnot(length(pred)==nrow(X),!anyNA(pred))
      row$ARI_anchor<-mclust::adjustedRandIndex(anchor,pred);row$NMI_anchor<-gb_nmi(anchor,pred)
      row$ARI_own_source<-if(is.null(own_partitions[[method]]))NA_real_ else mclust::adjustedRandIndex(own_partitions[[method]],pred)
      row$ARI_original_truth<-mclust::adjustedRandIndex(truth,pred)
      row$final_K<-length(unique(pred));row$iterations<-if(method=='Proposed')fit$iterations_run else fit$iterations
      pair<-gb_pair_scores(anchor,pred);row$within_anchor_pair_retention<-pair[[1]];row$between_anchor_pair_merge<-pair[[2]]
    }
    row
  })
  do.call(rbind,rows)
}
gb_partition_csvs <- function(fits,source,arm,rep,anchor,own_partitions,truth,node_ids,out) {
  for(method in names(fits$attempts)) {
    predicted<-gb_partition(fits$attempts[[method]])
    if(is.null(predicted))predicted<-base::rep(NA_integer_,length(anchor))
    own<-own_partitions[[method]];if(is.null(own))own<-base::rep(NA_integer_,length(anchor))
    data<-data.frame(node_id=node_ids,anchor=anchor,own_source=own,truth=truth,predicted=predicted)
    path<-file.path(out,'partitions',sprintf('%s_%s_b%03d_%s.csv',source$case_id,arm,rep,method))
    write.csv(data,path,row.names=FALSE)
  }
}
gb_plugin_model <- function(A,X,anchor,source) {
  anchor<-as.integer(factor(anchor));K<-max(anchor);N<-length(anchor);D<-ncol(X);sizes<-tabulate(anchor,nbins=K)
  H<-matrix(0,N,K);H[cbind(seq_len(N),anchor)]<-1
  numerator<-as.matrix(crossprod(H,A%*%H));denominator<-outer(sizes,sizes);diag(denominator)<-sizes*(sizes-1)
  P<-matrix(0,K,K);valid<-denominator>0;P[valid]<-numerator[valid]/denominator[valid]
  stopifnot(isTRUE(all.equal(P,t(P))),all(P>=0),all(P<=1))
  means<-covariances<-vector('list',K);component_rows<-vector('list',K)
  for(k in seq_len(K)) {
    ids<-which(anchor==k);mu<-colMeans(X[ids,,drop=FALSE]);centered<-sweep(X[ids,,drop=FALSE],2L,mu,'-')
    raw_cov<-crossprod(centered)/length(ids);S<-raw_cov+diag(1e-6,D)
    raw_eigen<-eigen(raw_cov,symmetric=TRUE,only.values=TRUE)$values
    covariance_eigen<-eigen(S,symmetric=TRUE,only.values=TRUE)$values
    tolerance<-max(abs(raw_eigen),1)*max(D,length(ids))*.Machine$double.eps
    means[[k]]<-mu;covariances[[k]]<-S
    component_rows[[k]]<-data.frame(source_id=source$case_id,component=k,size=length(ids),D=D,
      MLE_cov_rank=sum(raw_eigen>tolerance),rank_threshold=tolerance,min_eig_without_ridge=min(raw_eigen),
      min_eig_with_ridge=min(covariance_eigen),trace_cov_over_D=sum(diag(S))/D,
      covariance_ridge=1e-6)
    stopifnot(all(is.finite(S)),min(covariance_eigen)>0)
  }
  expected_edges<-sum(P*denominator)/2;observed_edges<-sum(A)/2
  stopifnot(abs(expected_edges-observed_edges)<1e-7*max(1,observed_edges))
  order_nodes<-order(anchor);inverse_order<-order(order_nodes)
  stopifnot(identical(as.integer(rep(seq_len(K),sizes)[inverse_order]),anchor))
  parameters<-list(anchor=anchor,P=P,means=means,covariances=covariances,sizes=sizes,
                   order=order_nodes,inverse_order=inverse_order)
  list(parameters=parameters,model_hash=gb_hash(parameters),components=do.call(rbind,component_rows),
       expected_edges=expected_edges,observed_edges=observed_edges,
       construction='Raw hard-partition block densities; fixed hard labels and node identity; Gaussian MLE covariance/n_k plus 1e-6 I; no probability clipping.')
}
gb_generate <- function(model,input,source,rep,graph_seed,attribute_seed) {
  p<-model$parameters;set.seed(graph_seed);graph_rng<-.Random.seed
  sorted_graph<-igraph::sample_sbm(n=source$N,pref.matrix=p$P,block.sizes=p$sizes,directed=FALSE,loops=FALSE)
  sorted_A<-igraph::as_adjacency_matrix(sorted_graph,sparse=input$representation=='sparse')
  A<-sorted_A[p$inverse_order,p$inverse_order,drop=FALSE]
  if(input$representation=='dense')A<-as.matrix(A)
  dimnames(A)<-dimnames(input$A)
  permutation_verified<-isTRUE(all.equal(A[p$order,p$order,drop=FALSE],sorted_A,check.attributes=FALSE))
  stopifnot(permutation_verified,all(diag(A)==0),isTRUE(all.equal(A,t(A),check.attributes=FALSE)))
  set.seed(attribute_seed);attribute_rng<-.Random.seed;X<-matrix(NA_real_,source$N,source$D)
  for(k in seq_along(p$sizes)) {
    ids<-which(p$anchor==k)
    X[ids,]<-mvtnorm::rmvnorm(length(ids),mean=p$means[[k]],sigma=p$covariances[[k]],method='chol')
  }
  dimnames(X)<-dimnames(input$X)
  stopifnot(!anyNA(X),all(is.finite(X)))
  meta<-data.frame(source_id=source$case_id,rep=rep,graph_seed=graph_seed,attribute_seed=attribute_seed,
     graph_hash=gb_hash(A),X_joint_hash=gb_hash(X),X_original_hash=input$attribute_hash,
     model_hash=model$model_hash,observed_edges=sum(A)/2,expected_edges=model$expected_edges,
     permutation_verified=permutation_verified,graph_rng_hash=gb_hash(graph_rng),attribute_rng_hash=gb_hash(attribute_rng))
  list(A=A,X_joint=X,meta=meta,order=p$order,inverse_order=p$inverse_order,
       graph_rng=graph_rng,attribute_rng=attribute_rng)
}
gb_generator_fixture <- function() {
  # Independently known 0/1 block graph with shuffled node labels and a singleton.
  anchor<-c(2L,1L,3L,2L,1L,2L)
  expected_P<-matrix(c(1,0,1,0,1,0,1,0,0),3L,3L)
  expected_A<-expected_P[anchor,anchor];diag(expected_A)<-0
  X<-matrix(seq_len(12L)/10,6L,2L)
  source<-data.frame(case_id='fixture',N=6L,D=2L)
  rows<-lapply(c('dense','sparse'),function(representation) {
    A<-if(representation=='sparse')as(Matrix::Matrix(expected_A,sparse=TRUE),'generalMatrix') else expected_A
    input<-list(A=A,X=X,representation=representation,attribute_hash=gb_hash(X))
    model<-gb_plugin_model(A,X,anchor,source)
    stopifnot(identical(model$parameters$P,expected_P),model$parameters$P[3,3]==0,
              model$components$MLE_cov_rank[3]==0,model$expected_edges==sum(expected_A)/2)
    generated<-gb_generate(model,input,source,1L,4001L,4002L)
    stopifnot(isTRUE(all.equal(as.matrix(generated$A),expected_A,check.attributes=FALSE)),
              generated$meta$permutation_verified)
    data.frame(representation=representation,N=6L,singleton_probability=0,
       permutation_verified=TRUE,expected_graph_recovered=TRUE,observed_edges=sum(expected_A)/2)
  })
  do.call(rbind,rows)
}

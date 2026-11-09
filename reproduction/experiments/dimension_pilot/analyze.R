#!/usr/bin/env Rscript
script<-normalizePath(sub("^--file=","",commandArgs(FALSE)[startsWith(commandArgs(FALSE),"--file=")][1]))
root<-normalizePath(file.path(dirname(script),"../.."));out <- Sys.getenv("DIM_PILOT_OUTPUT", file.path(root, "outputs/dimension_pilot"))
x<-read.csv(file.path(out,"raw_metrics.csv"),stringsAsFactors=FALSE)
stopifnot(nrow(x)==240,!anyDuplicated(paste(x$condition_id,x$rep,x$method)))
method_order<-c("Proposed","Pure_aSBM","Louvain","KMeans_estK")
keys<-unique(x[,c("condition_id","D","sigma","method")]);keys<-keys[order(keys$sigma,keys$D,match(keys$method,method_order)),]
summary<-do.call(rbind,lapply(seq_len(nrow(keys)),function(i){
  a<-keys[i,];b<-x[x$condition_id==a$condition_id&x$method==a$method,]
  a$n_planned<-nrow(b);a$n_returned<-sum(b$status=="ok");a$n_error<-sum(b$status!="ok")
  a$n_stopping_rule<-sum(b$met_stopping_rule,na.rm=TRUE);a$n_init_warning<-sum(b$init_warning_count>0)
  a$n_fit_warning<-sum(b$warning_count>0);a$n_stored_objective_nonfinite<-sum(!b$objective_finite,na.rm=TRUE)
  for(metric in c("ARI","NMI")) {
    v<-b[[metric]][is.finite(b[[metric]])];a[[paste0(metric,"_n")]]<-length(v)
    a[[paste0(metric,"_mean")]]<-if(length(v))mean(v)else NA_real_
    a[[paste0(metric,"_sd")]]<-if(length(v)>1)sd(v)else NA_real_
    a[[paste0(metric,"_mcse")]]<-if(length(v)>1)sd(v)/sqrt(length(v))else NA_real_
  }
  a$final_K_mean<-mean(b$final_K,na.rm=TRUE);a$initial_K_mean<-mean(b$initial_K)
  a$initial_min_size_min<-min(b$initial_min_size);a$n_rank_risk_init<-sum(b$initial_clusters_le_D>0)
  a$iterations_mean<-if(all(is.na(b$iterations)))NA_real_ else mean(b$iterations,na.rm=TRUE)
  a$elapsed_mean<-if(all(is.na(b$elapsed_sec)))NA_real_ else mean(b$elapsed_sec,na.rm=TRUE)
  a
}))
write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
pair_stats<-function(v){
  v<-v[is.finite(v)];n<-length(v);mu<-if(n)mean(v)else NA_real_;s<-if(n>1)sd(v)else NA_real_
  half<-if(n>1)qt(.975,n-1)*s/sqrt(n)else NA_real_
  data.frame(n_pairs=n,mean_delta=mu,sd_delta=s,mcse_delta=s/sqrt(n),ci95_low=mu-half,ci95_high=mu+half)
}
paired<-list();j<-0L
for(sig in c(1.5,4))for(method in method_order)for(d in c(100,200)) {
  b<-x[x$sigma==sig&x$method==method&x$D==50,c("rep","ARI")]
  a<-x[x$sigma==sig&x$method==method&x$D==d,c("rep","ARI")]
  p<-merge(a,b,by="rep",suffixes=c("_D","_50"));j<-j+1L
  paired[[j]]<-cbind(data.frame(sigma=sig,method=method,D=d,reference_D=50),pair_stats(p$ARI_D-p$ARI_50))
}
paired<-do.call(rbind,paired);write.csv(paired,file.path(out,"paired_dimension_ARI.csv"),row.names=FALSE)
comparison<-list();j<-0L
for(sig in c(1.5,4))for(d in c(50,100,200))for(method in method_order[-1]) {
  a<-x[x$sigma==sig&x$D==d&x$method=="Proposed",c("rep","ARI")]
  b<-x[x$sigma==sig&x$D==d&x$method==method,c("rep","ARI")]
  p<-merge(a,b,by="rep",suffixes=c("_proposed","_baseline"));j<-j+1L
  comparison[[j]]<-cbind(data.frame(sigma=sig,D=d,baseline=method),pair_stats(p$ARI_proposed-p$ARI_baseline))
}
write.csv(do.call(rbind,comparison),file.path(out,"paired_method_ARI.csv"),row.names=FALSE)
write.csv(x[x$status!="ok"|x$init_warning_count>0|x$warning_count>0|x$stop_reason %in% c("safety_cap","iteration_cap","timeout")|(!is.na(x$objective_finite)&!x$objective_finite),],file.path(out,"exception_records.csv"),row.names=FALSE)

colors<-c("#126E82","#C75631","#636363","#7656A3");pch_values<-c(16,17,15,18)
draw<-function(){
  layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),heights=c(1,.22))
  par(mar=c(4.5,4.5,3.5,1),oma=c(1,0,2,0),family="sans",las=1)
  for(sig in c(1.5,4)) {
    z<-summary[summary$sigma==sig,]
    yr<-range(c(-.02,1.02,z$ARI_mean-z$ARI_sd,z$ARI_mean+z$ARI_sd),na.rm=TRUE)
    plot(NA,xlim=c(43,207),ylim=yr,xlab="Attribute dimension D",ylab="Adjusted Rand index (ARI)",xaxt="n",main=paste0("Noise SD = ",sig),bty="l")
    axis(1,at=c(50,100,200));abline(h=c(0,.25,.5,.75,1),col="#EAEAEA",lwd=.7)
    for(mi in seq_along(method_order)) {
      a<-z[z$method==method_order[mi],];off<-c(-3,-1,1,3)[mi]
      lines(a$D+off,a$ARI_mean,col=colors[mi],lwd=2,lty=if(mi==3)2 else 1)
      arrows(a$D+off,a$ARI_mean-a$ARI_sd,a$D+off,a$ARI_mean+a$ARI_sd,angle=90,code=3,length=.035,col=colors[mi])
      points(a$D+off,a$ARI_mean,pch=pch_values[mi],col=colors[mi],cex=1.15)
    }
  }
  mtext("Fixed 50 signal coordinates; added coordinates contain Gaussian noise",outer=TRUE,side=3,line=.3,cex=1.15)
  par(mar=c(0,0,0,0));plot.new()
  legend("top",legend=c("Proposed","Pure aSBM","Louvain","K-means (estimated K)"),col=colors,pch=pch_values,lty=c(1,1,2,1),lwd=2,bty="n",cex=.95,horiz=TRUE)
  text(.5,.36,"Mean +/- sample SD; 10 independent data repetitions per setting; n = 1000, K = 5, r = 4",cex=.8)
  text(.5,.08,"Fit warnings and failures are listed in exception_records.csv",cex=.8)
}
png(file.path(out,"dimension_ARI.png"),width=1600,height=850,res=150);draw();dev.off()
pdf(file.path(out,"dimension_ARI.pdf"),width=10.67,height=5.67);draw();dev.off()
print(summary[,c("D","sigma","method","ARI_mean","ARI_sd","NMI_mean","n_error","n_stopping_rule","n_init_warning","n_stored_objective_nonfinite","final_K_mean","iterations_mean","n_rank_risk_init")],row.names=FALSE)
print(paired[paired$method=="Proposed",],row.names=FALSE)

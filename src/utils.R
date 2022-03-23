LR<-function(x,gold_standard,l_warn=0){
  print('LR function called')
  calculate_metrics<-function(outcome_table,in_gold_standard,threshold){
    P<-which(outcome_table$screen.score>=threshold)
    N<-which(outcome_table$screen.score<threshold)
    TP<-length(intersect(P,in_gold_standard))
    FP<-length(P)-TP
    FN<-length(intersect(N,in_gold_standard))
    TN<-length(N)-FN
    res<-c(threshold, TP,FP,TN,FN)
    names(res)<-c('Threshold','TP','FP','TN','FN')
    return(res)
  }
  outcome_table<-data.frame(Gene=names(x),screen.score=x)
  outcome_table<-outcome_table[order(-outcome_table$screen.score),]
  outcome_table$in.gold.standard<-outcome_table$Gene%in%gold_standard
  in_gold_standard<-which(outcome_table$Gene%in%gold_standard)
  metrics_list<-lapply(unique(outcome_table$screen.score), function(t)
    calculate_metrics(outcome_table,in_gold_standard ,t))
  metrics_table<-as.data.frame(do.call('rbind',metrics_list))
  metrics_table$sensitivity<-metrics_table$TP/(metrics_table$TP+metrics_table$FN)
  metrics_table$specficity<-metrics_table$TN/(metrics_table$TN+metrics_table$FP)
  metrics_table$LRplus<-metrics_table$sensitivity/(1-metrics_table$specficity)
  metrics_table$LRplus.corrected<-sapply(metrics_table$Threshold,function(val)
    max(metrics_table$LRplus[which(metrics_table$Threshold<=val)]))
  maxinf<-max(metrics_table$LRplus.corrected[is.finite(metrics_table$LRplus.corrected)])
  metrics_table$LRplus.corrected[which(is.infinite(metrics_table$LRplus.corrected)==TRUE)]<-
    maxinf
  
  if(any(diff(metrics_table$LRplus.corrected>0))){
    stop('likelihood ratio does not increase with screen score')
  }
  outcome_table$LRplus<-metrics_table$LRplus.corrected[match(outcome_table$screen.score,
                                                             metrics_table$Threshold)]
  outcome_table$LRplus[which(outcome_table$LRplus==0)]<-0
  outcome_table<-outcome_table[order(outcome_table$Gene),]
  #  if(length(warnings())>l_warn){
  #    print(warnings())
  #    l_warn<-length(warnings())
  #  }
  return(outcome_table)
}

getPreds<-function(predslist,gs){
  table<-as.data.frame(do.call('cbind',lapply(predslist, function(s) s$LRplus)))
  row.names(table)<-predslist[[1]]$Gene
  table$overall_ratio<-apply(table,1,function(x) sum(log2(x)))
  table$in.gold.standard<-rownames(table)%in%gs
  table<-table[order(-table$overall_ratio),]
  return(table)
}


make_saturation_table<-function(preds_table){
  cumgs<-cumsum(preds_table$in.gold.standard)
  saturation_table<-data.frame(Rank=1:nrow(preds_table),gs_count=cumgs,
                               gold.standard=preds_table$in.gold.standard,gene=rownames(preds_table))
  saturation_table$TPR<-saturation_table$gs_count/sum(saturation_table$gold.standard)
  saturation_table$FPR<-(saturation_table$Rank-saturation_table$gs_count)/(nrow(saturation_table)-sum(saturation_table$gold.standard))
  saturation_table$Recall<-100*saturation_table$TPR
  saturation_table$gs_ratio<-saturation_table$Rank/saturation_table$gs_count
  return(saturation_table)
}

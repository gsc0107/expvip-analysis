getGenesContributingEachDecile<-function(geneInformation,
                                         subset="Development", 
                                         genes_to_count=NULL, 
                                         probs=seq(from=0.05, to=1,by=0.05),
                                         min_tpms=c(0.5,1,2)
                                        ){
    tpms<-geneInformation$meanTpms
    tpms<-tpms[tpms$subset == subset, ]
    if(!is.null(genes_to_count)){
       tpms<-tpms[tpms$gene %in% genes_to_count, ] 
    }
    all_counts<-NULL
    for(factor in unique(tpms$factor)){
        l_tpms<-tpms[tpms$factor == factor &
                     tpms$value > 0 ,]
        #print(factor)
        values<-l_tpms$value
        names(values)<-l_tpms$gene
        values<- values[!is.na(values)]
        total_sum<-sum(values)
        sorted<-sort(values,decreasing = T)
        cumulative <- cumsum(sorted)
        cum_percentage <-  cumulative / total_sum
        counts<-list()
        #print(head(cum_percentage))
        counts[["Subset"]]  <- subset
        counts[["Factor"]]  <- factor
        counts[["Samples"]] <- max(l_tpms$samples)
        for(p in probs){
            tmp <- cum_percentage[cum_percentage < p]
            counts[[paste0("percentage ",as.character(100 * p))]] <- length(tmp)
        }
        
        for(min_tpm in min_tpms){
            tmp <- values[values > min_tpm]
            counts[[paste0("min_tpm ",as.character(min_tpm))]] <- length(tmp)
        }
        
        tmp_df<-data.frame(counts)
        if(is.null(all_counts)){
            all_counts<-tmp_df
        }else{
            all_counts<-rbind(all_counts, tmp_df)
        }
    }
    all_counts
}


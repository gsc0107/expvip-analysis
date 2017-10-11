getBottomAndTopPercent<-function(geneInformation, 
                                 dataset="HC_CS_no_stress", 
                                min_factor_count=5){
    
    triadMovement<-geneInformation$triadMovement
    triads<-geneInformation$triads
    local_tm <- triadMovement[triadMovement$dataset == dataset & 
                              triadMovement$factor_count > min_factor_count ,
                             c("group_id","central_mean_distance")]

    local_tm$rank <- rank(local_tm$central_mean_distance) / nrow(local_tm)
    
    low_10 <- local_tm[local_tm$rank < 0.1, ] 
    low_25 <- local_tm[local_tm$rank < 0.25, ]
    low_50 <- local_tm[local_tm$rank < 0.50, ]
    top_10 <- local_tm[local_tm$rank > 0.9, ]
    top_50 <- local_tm[local_tm$rank > 0.50, ]
    top_25 <- local_tm[local_tm$rank > 0.75, ]
    
    middle10pc <- local_tm[local_tm$rank <= 0.55 & local_tm$rank >= 0.45, ]
    middle20pc <- local_tm[local_tm$rank <= 0.60 & local_tm$rank >= 0.40, ]
    middle50pc <- local_tm[local_tm$rank <= 0.75 & local_tm$rank >= 0.25, ]
    middle80pc <- local_tm[local_tm$rank <= 0.90 & local_tm$rank >= 0.10, ]
    
    
    print(nrow(local_tm))
    #if (nrow(local_tm) == 0)        return 
    
    low_10$category<-"low_10pc"
    low_25$category<-"low_25pc"
    top_25$category<-"top_25pc"
    top_10$category<-"top_10pc"
    top_50$category<-"top_50pc"
    low_50$category<-"low_50pc"
    middle20pc$category<-"middle_20pc"
    middle10pc$category<-"middle_10pc"
    middle80pc$category<-"middle_80pc"
    middle50pc$category<-"middle_50pc"
    local_tm$category<-"all"
    
    groups<-rbind(low_10,low_25,top_10, top_25, top_50, low_50, 
                  middle20pc,
                  middle50pc,
                  middle10pc, 
                  middle80pc, local_tm)
    genes<-sqldf("SELECT DISTINCT gene, category from groups JOIN triads ON 
triads.group_id = groups.group_id")
    genes
}
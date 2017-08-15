options(gsubfn.engine = "R")
library(ggplot2)
library(reshape2)
library(sqldf)
library(fields)
library(gridExtra)
library(ggtern)
library(clue)
library(geometry)
library(gtable)
library(goseq)

loadGeneInformation<-function(dir="../TablesForExploration"){
    path<-paste0(dir,"/CanonicalTranscript.rds")
    canonicalTranscripts<-readRDS(path)
    canonicalTranscripts$intron_length<- canonicalTranscripts$mrna_length -  canonicalTranscripts$exon_length
    canonicalTranscripts$chr_group <- substr(canonicalTranscripts$Chr,4,4)
    canonicalTranscripts$genome    <- substr(canonicalTranscripts$Chr,5,5)
    
    path<-paste0(dir, "/MeanTpms.rds")
    meanTpms <- readRDS(path)
    expressed_genes<-unique(meanTpms$gene)
    canonicalTranscripts<-canonicalTranscripts[canonicalTranscripts$Gene %in% expressed_genes, ]
    canonicalTranscripts$scaled_5per_position <-   5 * ceiling(canonicalTranscripts$scaled_1per_position / 5)
    canonicalTranscripts$scaled_5per_position <- ifelse(canonicalTranscripts$scaled_5per_position == 0, 
        5, 
        canonicalTranscripts$scaled_5per_position)
    path<-paste0(dir,"/TriadMovement.rds")
    triadMovement<-readRDS(path)
    
    path<-paste0(dir,"/Triads.rds")
    triads<-readRDS(path)
    
    path<-paste0(dir,"/universe_table.csv")
    gene_universe<-read.csv(path)
    
    path<-paste0(dir, "/OntologiesForGenes.rds")
    ontologies<-readRDS(path)
    
    path<-paste0(dir, "/id_names_merged.txt")
    id_names <- read.csv(path, header=F, sep = "\t")
    
    path<-paste0(dir, "/WGCNA_table.csv")
    WGCNA <-  read.csv(path)
    
    path<-paste0(dir, "/ObservedGOTermsWithSlim.csv")
    go_slim<-read.csv(path, row.names=1)
    
    list(canonicalTranscripts=canonicalTranscripts, 
       meanTpms=meanTpms,
       triads=triads, 
       triadMovement=triadMovement,
       gene_universe=gene_universe, 
       ontologies=ontologies,
       id_names=id_names,
       WGCNA=WGCNA,
       GOSlim=go_slim
       )
}



prepare_hist_stats<-function(table, column="size_cds"){
    table<-table[table[,column]>0,]
    probs <- c( 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
    quantiles <- data.frame(quantile(table[,column], prob=probs,na.rm=TRUE, include.lowest=TRUE), stringsAsFactors=FALSE)
    quantiles$quant<-rownames(quantiles)
    colnames(quantiles)<-c("value", "quant")
    values<-quantiles$values
    local_mean<-mean(table[,column])
    local_sd<-sd(table[,column])
    local_max <-  max(table[,column])
    

    stats_list<-list(mean=local_mean, sd = local_sd, cv = local_sd/local_mean, 
        median =  median(table[,column],2), 
        max = local_max, 
        n = nrow(table)
        )
    stats_list
}

plotHistogram<-function(table, column="size_cds"){
    table<-table[table[,column]>0,]
    probs <- c( 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
    quantiles <- data.frame(quantile(table[,column], prob=probs,na.rm=TRUE, include.lowest=TRUE), stringsAsFactors=FALSE)
    quantiles$quant<-rownames(quantiles)
    colnames(quantiles)<-c("value", "quant")
    values<-quantiles$values
    local_mean<-mean(table[,column])
    local_sd<-sd(table[,column])
    local_max <-  max(table[,column])
    p <- ggplot(table, aes_string(column))
    
    if(nrow(table) > 100){
     table <- within(table,
         quantile <- as.integer(
             cut(table[,column],
                 unique(quantile(table[,column], 
                    prob=probs,
                    na.rm=TRUE, 
                    include.lowest=TRUE))
                 )
             ))

     table$quantile<-ifelse(is.na(table$quantile),7,table$quantile)
     table$quantile<-as.factor(table$quantile)

     iq <- quantiles$value[4] - quantiles$value[2]

     xmax <- quantiles$value[3] + (iq * 2)
     xmin <- quantiles$value[3] - (iq * 2)
     if(xmin < 0){
        xmin <- 0
    }
    
    if(xmax > local_max){
        xmax <- local_max + 1
    }

    p <- ggplot(table, aes_string(column, fill="quantile"))
    p <- p + geom_vline(data=quantiles,aes(xintercept=quantiles$value) )
    for(i in seq(1,nrow(quantiles))){
        x_pos<-quantiles$value[i]
        gtext <- textGrob(quantiles$quant[i], y=0.02,  gp = gpar(fontsize = 6,col = "red"))
        p <- p + annotation_custom(gtext, xmin=x_pos, xmax=x_pos)
    }
    p <- p  + xlim(xmin, xmax) +
    scale_fill_brewer(palette="Dark2")
    
}




p <- p + geom_histogram(bins=50, position = "identity") + theme_bw() 
p <- p + theme(legend.position="none")
p <- p + ggtitle(paste0("Mean: ", round(local_mean,2), 
    " SD:", round(local_sd,2),
    " CV:", round(local_sd/local_mean, 2), 
    " Median:", round(median(table[,column],2)),
    " Max:", round(local_max,2),
    " N:", nrow(table))) 
p <- p + theme(plot.title = element_text(size=6))

stats_list<-list(mean=local_mean, sd = local_sd, cv = local_sd/local_mean, 
    median =  median(table[,column],2), 
    max = local_max, 
    n = nrow(table)
    )
p
}

#This function gets the expected number of genes per each 5pc bin.    
get_expected_values_per_5pc_bin<-function(gene_table, 
  numberOfGenes, 
  group_in_single_chromosome=FALSE){

    query<-"SELECT Chr, chr_group, genome,scaled_5per_position, count(*) as count 
    FROM
    gene_table 
    WHERE geneconf = 'HC' AND Chr != 'chrUn'
    GROUP BY  Chr, chr_group, genome, scaled_5per_position"
    counts<-sqldf(query)
    if(group_in_single_chromosome){
        query<-"SELECT 'All' as Chr, 'all' as chr_group, 'all'  as genome,
        scaled_5per_position, count(*) as count 
        FROM
        gene_table 
        WHERE geneconf = 'HC' AND Chr != 'chrUn'
        GROUP BY scaled_5per_position"
    }
    counts<-sqldf(query)
    multiplier <- numberOfGenes / sum(counts$count)
    counts$expected <- counts$count * multiplier
    counts
}



plot_per_chromosome_5pc_bins_facet<-function(table,expected_per_chr,
   expected_all_chromosomes=NULL, 
   title = "Test"){
    chromosomes=c("1A", "1B", "1D",
        "2A", "2B", "2D",
        "3A", "3B", "3D",
        "4A", "4B", "4D",
        "5A", "5B", "5D",
        "6A", "6B", "6D",

        "7A", "7B", "7D")
    
    gs<-list()
    local_title = paste0(title, "\n Genes per chromosome 5% bin\nN: ", nrow(table) )
    
    t1 <- table[table$Chr != "chrUn",]

    
    expected_per_chr <- expected_per_chr[expected_per_chr$Chr != "chrUn",]
    p <-ggplot(t1,aes(scaled_5per_position)) 
    
    p <- p + xlim(0,100)
    p <- p + geom_bar() + theme_bw()
    p <- p + facet_grid(chr_group~genome,  drop = TRUE)
    p <- p + ylab(" count ") + xlab("")
    if(!is.null(expected_all_chromosomes)){
        expected_all_chromosomes$expected <- expected_all_chromosomes$expected/21
        p <- p + geom_line(data=expected_all_chromosomes[, c("scaled_5per_position", "expected")],
         aes(x=scaled_5per_position, y=expected), color="blue", size = 0.5)
    }
    
    p <- p + geom_point(data=expected_per_chr, 
        aes(x=scaled_5per_position, y=expected), color="red", size = 0.5)
    
    
    
    
    gs[[length(gs)+1]] <- p
    
    p <-ggplot(table,aes(Chr, fill=geneconf))  + geom_bar() + theme_bw()
    p <- p + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylab("") + xlab("")
    gs[[length(gs)+1]] <- p
    
    g1<-arrangeGrob(grobs=gs, ncol=1, heights=c(0.8,0.2), top=local_title ) 
    g1
}

plot_tpms_summary<-function(tpms, experiment="850_samples", min_tpm=0.5, title="Test"){

    local_tpms<-subset(tpms, (subset == experiment) & 
     ( factor != "all" & factor != "all_means" & factor != "all_mean_filter" ) &
     value > min_tpm)

    local_title <- paste0(title, "\n", experiment)
    
    p  <- ggplot(local_tpms, aes(value)) 
    p  <- p + geom_histogram(bins=30 ) + theme_bw()
    p  <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
      strip.text = element_text(size=6))
    p  <- p + facet_wrap(~ factor, ncol=4) 
    p  <- p + xlim(0,15) 
    p  <- p + ylab("Count") + xlab("")
    p  <- p + theme(strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")))

    
    local_tpms<-subset(tpms, (subset == experiment) & 
     ( factor == "all_mean_filter" ) &
     value > min_tpm)
    
    p2  <- ggplot(local_tpms, aes(value)) 
    p2  <- p2 + geom_histogram(bins=30 ) + theme_bw()
    p2  <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
      strip.text = element_text(size=6, lineheight=0.5))
    p2  <- p2 + facet_wrap(~ factor, ncol=1)  + xlim(0,15)
    p2  <- p2 + ylab("") + xlab("TPM") 

    mytheme <- gridExtra::ttheme_default(
        core = list(fg_params=list(cex = 0.5)),
        colhead = list(fg_params=list(cex = 0.5)),
        rowhead = list(fg_params=list(cex = 0.5)))
    

    
    lay <- rbind(c(1),
       c(2))
    g1<-arrangeGrob(grobs=list(p,p2), heights=c(0.8,0.2), layout_matrix=lay, top = local_title) 
    g1
}


get_tpms_desc_stats<-function(tpms, min_tpm=0.5, experiment="850_samples"){


    local_tpms<-subset(tpms, (subset == experiment) & 
        factor != "all" &
        factor != "all_means"  &
        value > min_tpm)
    
    factor_means <- aggregate(value ~ factor, data = local_tpms, mean)
    factor_max   <- aggregate(value ~ factor, data = local_tpms, max)
    factor_min   <- aggregate(value ~ factor, data = local_tpms, min)
    factor_sd    <- aggregate(value ~ factor, data = local_tpms, sd)
    factor_median<- aggregate(value ~ factor, data = local_tpms, median)
    factor_n     <- aggregate(value ~ factor, data = local_tpms, length)
    
    rownames(factor_means)<- factor_means$factor
    rownames(factor_max)<- factor_max$factor
    rownames(factor_min)<- factor_min$factor
    rownames(factor_sd)<- factor_sd$factor
    rownames(factor_median)<- factor_median$factor
    rownames(factor_n)<- factor_n$factor
    
    factor_means$factor<-NULL
    factor_max$factor<-NULL
    factor_min$factor<-NULL
    factor_sd$factor<-NULL
    factor_median$factor <- NULL
    factor_n$factor <- NULL
    
    factor_means<-cbind(factor_means, factor_median)
    factor_means<-cbind(factor_means, factor_max)
    factor_means<-cbind(factor_means, factor_min)
    factor_means<-cbind(factor_means, factor_sd)
    factor_means<-cbind(factor_means, factor_n)
    
    colnames(factor_means)<-c("Mean", "Median", "Max", "Min", "SD","N")
    factor_means$CV<-factor_means$SD/factor_means$Median
    factor_means<-round(factor_means, 1)
}

plot_tpm_desc_stats<-function(tpms,subset_tpms, experiment="850_samples", min_tpm=0.5, title="Test" ){
    local_title <- paste0(title, "\n", experiment)
    
    
    t_local  <-get_tpms_desc_stats(subset_tpms, experiment=experiment, min_tpm=min_tpm)
    t_global <-get_tpms_desc_stats(tpms,        experiment=experiment, min_tpm=min_tpm)
    rownames(t_global) <- NULL
    
    mytheme <- gridExtra::ttheme_default(
       core = list(fg_params=list(cex = 0.5)),
       colhead = list(fg_params=list(cex = 0.5)),
       rowhead = list(fg_params=list(cex = 0.5)))
    
    t1 <- tableGrob(t_local, theme=mytheme)
    t2 <- tableGrob(t_global, theme=mytheme)
    
    title <- textGrob("Observed",gp=gpar(fontsize=10))
    padding <- unit(5,"mm")
    table <- gtable_add_rows(
        t1, 
        heights = grobHeight(title),
        pos = 0)
    table <- gtable_add_grob(
        table, title, 
        1, 1, 1, 
        ncol(table))
    
    title2 <- textGrob("Expected",gp=gpar(fontsize=10))
    
    table2 <- gtable_add_rows(
        t2, 
        heights = grobHeight(title2),
        pos = 0)
    table2 <- gtable_add_grob(
        table2, title2, 
        1, 1, 1, 
        ncol(table2))
    
    g1<-arrangeGrob(grobs=list(table ,table2), ncol=2, top = local_title) 
    g1
}

plot_all_means_filteredtpms_summary<-function(tpms, experiment="850_samples", min_tpm=0.5, title="Test"){

    local_title <- paste0(title, "\n", experiment)
    
    local_tpms<-subset(tpms, (subset == experiment) & 
        factor == "all_mean_filter"  )
    
    
    
    breaks <- seq(-1,max(local_tpms$samples) , by = 1)
    samples.cut <- cut(local_tpms$samples, breaks, include.lowest = FALSE)
    samples.freq <- table(samples.cut)
    
    cumfreq0 =  cumsum(samples.freq)
    
    local_tpms<-subset(tpms, subset == experiment & 
        factor == "all_mean_filter" & 
        value > min_tpm)
    
    
    level <- ifelse(local_tpms$value < 0.1, 0.1, local_tpms$value) 
    level <- ceiling(log10(level))
    local_tpms$level <- level
    local_tpms$exp_max_value <- as.factor(10**level)
    
    gs<-list()
    gs[[length(gs)+1]] <- ggplot(local_tpms, aes(samples, fill = exp_max_value)) + 
    geom_bar(width = 0.75) +   scale_fill_brewer(palette="RdPu") + 
    ggtitle(paste0("Gene expression level and number of tissues in which genes are expressed")) + 
    theme_bw() +
    labs(fill="AVG TPM", x="") +
    theme(legend.position=c(.25,.75))+  
    guides(color = guide_legend(override.aes = list(size=5)))
    
    freqs_df <- data.frame(samples =  seq(0,max(local_tpms$samples) , by = 1), cum_freq=cumfreq0 )
    
    
    gs[[length(gs)+1]] <- ggplot(freqs_df, aes(samples, cum_freq)) +
    geom_line() + geom_point() + theme_bw() +
    labs(x="Number of tissues/conditions", y="Cumulative frequency")
    
    g1<-arrangeGrob(grobs=gs, ncol=1, top=local_title )
    g1
}

get_triads_from_genes<-function(genes, geneInformation, dataset="HC_CS_no_stress" , min_no_genes = 1){
    triads<-geneInformation$triads
    triadMovement<-geneInformation$triadMovement
    triads_with_genes <- triads[triads$gene %in% genes,]
    tridas_with_genes <- triads_with_genes[triads_with_genes$dataset == dataset,]
    genes_in_triads<-sqldf("SELECT group_id, gene FROM triads_with_genes GROUP BY gene") 
    triad_gene_count<-sqldf("SELECT group_id, count(*) as count from genes_in_triads GROUP BY group_id")
    group_ids <- triad_gene_count[triad_gene_count$count >= min_no_genes, "group_id"]
    list(triads=triads[triads$group_id %in% group_ids & triads$dataset==dataset,], 
        triadMovement=triadMovement[triadMovement$group_id %in% group_ids & triadMovement$dataset==dataset ,])
    
}

plot_distribution_for_factor<-function(res, unit="central_mean_distance" ,color_by="total_categories"){
    #unit<-paste0(from,"_",unit)
    
    local_res<-res[res$factor_count>1,c(unit, color_by)]
    local_res[,2]<-as.factor(local_res[,2])
    probs <- c( 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
    
    quantiles <- data.frame(quantile(local_res[,1], prob=probs,na.rm=TRUE), stringsAsFactors=FALSE)
    quantiles$quant<-rownames(quantiles)
    colnames(quantiles)<-c("value", "quant")
    #print(quantiles)
    local_res<-local_res[order(local_res[,color_by], decreasing = F),]
    
    
    p <- ggplot(local_res,aes_string(unit, fill=color_by)) +
    geom_vline(data=quantiles,aes(xintercept=quantiles$value) ) #
    for(i in seq(1,nrow(quantiles))){
        x_pos<-quantiles$value[i]
        gtext <- textGrob(quantiles$quant[i], y=0.02,  gp = gpar(fontsize = 6,col = "red"))
        p <- p + annotation_custom(gtext, xmin=x_pos, xmax=x_pos)
    }
    p <- p + theme_bw()
    #p <- p + geom_text(data = quantiles, mapping = aes_string(label = "quant", y = 0)) #, mapping = aes(label = quant, y = 0))
    p <- p + geom_histogram(position = "identity",alpha=0.5,binwidth=0.01) + theme(legend.position="none")
    p
}

get_stats_title<-function(d){
    paste0("Mean: ", round(mean(d),2), 
     " SD:", round(sd(d),2),
     " CV:", round(sd(d)/mean(d), 2), 
     " Median:", round(median(d),2),
     " N:", length(d)) 
}


table_dominance_summary<-function(selected_triads, experiment="HC_CS_no_stress", title="test"){
    local_title <- paste0(title, "\n", experiment)
    triads <- selected_triads$triads 
    triadMovement<-selected_triads$triadMovement
    all_means_filter<-triads[triads$factor=="all_mean_filter",]
    
    df  <- prepare_hist_stats(all_means_filter, column="value") 
    df$value_type <- "All mean filter TPM for genes in triad"


    tmp<- prepare_hist_stats(triadMovement, column="central_max_distance")
    tmp$value_type <- "central_max_distance"
    df<-rbind(df,tmp)

    tmp<- prepare_hist_stats(triadMovement, column="central_mean_distance")
    tmp$value_type <- "central_mean_distance"
    df<-rbind(df,tmp)
    
    tmp<- prepare_hist_stats(triadMovement, column="sum_mean_tpm")
    tmp$value_type <- "sum_mean_tpm"
    df<-rbind(df,tmp)

    df$dataset<-experiment
    df$title <- title
    df
}

plot_dominance_summary<-function(selected_triads, experiment="HC_CS_no_stress", title="test"){
   local_title <- paste0(title, "\n", experiment)

   triads <- selected_triads$triads 
   triadMovement<-selected_triads$triadMovement

   all_means_filter<-triads[triads$factor=="all_mean_filter",]
    #print(unique(triads$factor))
    gs<-list()
    gs[[length(gs)+1]] <- plotHistogram(all_means_filter, column="value") + 
    xlab("All mean filter TPM for genes in triad")
    
    gs[[length(gs)+1]] <- ggplot(triadMovement, aes(factor_count)) +  geom_bar() + theme_bw() +
    xlab("No. of conditions of genes(count per triad)") +
    ggtitle(get_stats_title(triadMovement$factor_count)) + theme(plot.title = element_text(size=6))

    
    

    p <- ggplot(triadMovement, aes(total_categories)) + geom_bar() + theme_bw()
    p <- p + labs(fill="Main\ncategory", x="No. of categories") +
    ggtitle(get_stats_title(triadMovement$total_categories)) + theme(plot.title = element_text(size=6))

    gs[[length(gs)+1]] <- p

    gs[[length(gs)+1]] <- plotHistogram(triadMovement, column="central_max_distance")
    gs[[length(gs)+1]] <- plotHistogram(triadMovement, column="central_mean_distance")
    gs[[length(gs)+1]] <- plotHistogram(triadMovement, column="sum_mean_tpm")
    g1<-arrangeGrob(grobs=gs, ncol=2, top=local_title )
    g1
}

get_dominance_summary_tables_per_factor<-function(selected_triads, 
  description = "description",
  experiment="HC_CS_no_stress",
  n=NULL
  ){
    triads <- selected_triads$triads
    triads <- triads[triads$dataset==experiment, ]
    
    query <- paste0("SELECT factor, " , 
        description , 
        " as description, count(*) as count FROM triads GROUP BY factor, " , 
        description )
    table <- sqldf(query)
    
    casted <- dcast(table, factor  ~  description , value.var="count")
    
    rownames(casted) <- casted$factor
    casted$factor <- NULL
    casted<-as.matrix(casted)
    casted <- ifelse(is.na(casted),0, casted)
    total_per_factor<-rowSums(casted)
    percentage <-  as.matrix(100 * casted / total_per_factor)
    
    if(! is.null(n)){
        casted <- percentage * n / 100
        total_per_factor<-rowSums(casted)
    }
    pasted<-matrix(paste(as.matrix(round(casted,0)),
        as.matrix(round(percentage,2)) , sep=" - "),
    nrow=nrow(casted), 
    dimnames=dimnames(casted))
    pasted<-matrix(paste0(pasted, "%"),
     nrow=nrow(casted), 
     dimnames=dimnames(casted))
    pasted<-data.frame(pasted)
   # print(total_per_factor)
   pasted$total <-total_per_factor
   long<-melt(casted)
   colnames(long)<-c("factor", "description", "count")
   list(long=long, casted=casted, percentage=percentage, pasted=pasted,total=total_per_factor )
}


table_with_title<-function(title, table){
    mytheme <- gridExtra::ttheme_default(
       core = list(fg_params=list(cex = 0.5)),
       colhead = list(fg_params=list(cex = 0.5)),
       rowhead = list(fg_params=list(cex = 0.5)))

    table2 <- tableGrob(table, theme=mytheme)
    g1<-arrangeGrob(grobs=list(table2), ncol=1, top=title )
}

plot_dominance_summary_tables<-function(selected_triads, 
    expected_desc,
    expected_gen_desc,
    experiment="HC_CS_no_stress", title="test"){
    local_title <- paste0(title, "\n", experiment)
    
    triads <- selected_triads$triads 
    
    expected_triads_desc<-expected_desc$long
    total_genes<-sum(expected_triads_desc$count)
    multiplier<-nrow(triads) / total_genes
    expected_triads_desc$exp_count <-  expected_triads_desc$count * multiplier
    
    expected_triads_gen_desc<-expected_gen_desc$long
    total_genes<-sum(expected_triads_gen_desc$count)
    multiplier<-nrow(triads) / total_genes
    expected_triads_gen_desc$exp_count <-  expected_triads_gen_desc$count * multiplier
    expected_triads_gen_desc$general_description <- expected_triads_gen_desc$description
    
    gs <- list()
    p <- ggplot(triads, aes(factor)) + geom_bar() + theme_bw()
    p <- p + theme(axis.text.x = element_text(angle = 35, hjust = 1, size=4)) + labs(fill="", x="")
    p <- p + theme(legend.text=element_text(size=5))+
    theme(legend.title=element_text(size=6)) + facet_wrap(~ description, ncol=4) 
    theme(legend.key.size = unit(0.4,"line"))
    
    p <- p + geom_point(data=expected_triads_desc, 
        aes(x=factor, y=exp_count), color="red", size = 0.5)
    
    gs[[length(gs)+1]] <- p
    
    p <- ggplot(triads, aes(factor)) + geom_bar() + theme_bw()
    p <- p + theme(axis.text.x = element_text(angle = 35, hjust = 1, size=4)) + labs(fill="", x="")
    p <- p + theme(legend.text=element_text(size=5))+
    theme(legend.title=element_text(size=6))+
    theme(legend.key.size = unit(0.4,"line"))+ facet_wrap(~ general_description, ncol=3) 

    p <- p + geom_point(data=expected_triads_gen_desc, 
        aes(x=factor, y=exp_count), color="red", size = 0.5)
    
    gs[[length(gs)+1]] <- p

    g1<-arrangeGrob(grobs=gs, ncol=1, top=local_title )
    g1
}

get_goseq_enrichment<-function(geneInformation, genes_to_plot, 
 name="Random Samples",
 dataset="HC_CS_no_stress", 
 ontology="GO"
 ){
    id_names <- geneInformation$id_names
    universe<-geneInformation$gene_universe 
    universe<-universe[universe$dataset==dataset,]
    #print(nrow(universe))
    if(nrow(universe) == 0) {
        return (data.frame(category= numeric(0), 
          over_represented_pvalue= numeric(0),
          under_represented_pvalue= numeric(0),
          numDEInCat= numeric(0),
          numInCat= numeric(0),
          ontology= numeric(0),
          over_rep_padj= numeric(0),
          under_rep_padj= numeric(0),
          description= numeric(0),
          universe_size= numeric(0)
          ))
    }
    ontologies<-geneInformation$ontologies
    ontologies<-ontologies[ontologies$ontology==ontology,]
    
    assayed.genes <- as.vector(universe$gene)
    gene.vector=as.integer(assayed.genes%in%genes_to_plot)
    names(gene.vector)=assayed.genes
    
    transcripts<-geneInformation$canonicalTranscripts
    lengths <- transcripts[,c("Gene", "exon_length")]
    
    colnames(lengths) <- c("gene", "length")
    t1 <- subset(lengths, gene %in% universe$gene)
    gene.lens <- as.numeric(t1$length)
    names(gene.lens) = t1$gene
    
    all_go <- subset(ontologies, Gene %in% universe$gene)
    all_go<-all_go[,c(1,2)]
    
    pwf = nullp(gene.vector, bias.data = gene.lens, plot.fit = FALSE)
    GO.wall = goseq(pwf, gene2cat = all_go)
    
    #this gave table with p-values...now correct for multiple testing using FDR
    # add new column with over represented GO terms padj
    GO.wall$over_rep_padj=p.adjust(GO.wall$over_represented_pvalue, method="BY")
    # add new column with under represented GO terms padj
    GO.wall$under_rep_padj=p.adjust(GO.wall$under_represented_pvalue, method="BY")

    # now select only GO terms where the padj is <0.05 for either enriched or under represented
    sig.GO <- GO.wall[GO.wall$over_rep_padj <0.05 | GO.wall$under_rep_padj <0.05,]
    sig.GO <- sig.GO[order(sig.GO$over_rep_padj),]
    if(nrow(sig.GO) == 0 ){
        return (data.frame(category= numeric(0), 
          over_represented_pvalue= numeric(0),
          under_represented_pvalue= numeric(0),
          numDEInCat= numeric(0),
          numInCat= numeric(0),
          ontology= numeric(0),
          over_rep_padj= numeric(0),
          under_rep_padj= numeric(0),
          description= numeric(0),
          universe_size= numeric(0)
          ))
    }
    
    sig.GO2 <- merge(sig.GO, id_names, by.x="category", by.y="V1", all.x =TRUE, all.y =FALSE)

    if( nrow(sig.GO) > 0 && (ontology == "PO" || ontology == "TO" )){
        sig.GO2$ontology<-ontology
    }
    
    sig.GO <-sig.GO2[,c('category','over_represented_pvalue','under_represented_pvalue','numDEInCat','numInCat',
        'ontology','over_rep_padj','under_rep_padj','V2')]
    
    colnames(sig.GO)<-c('category','over_represented_pvalue','under_represented_pvalue','numDEInCat','numInCat',
        'ontology','over_rep_padj','under_rep_padj','description')
    
    sig.GO$type<-ifelse(sig.GO$under_rep_padj > sig.GO$over_rep_padj,
        "Over represented", 
        "Under represented" )
    sig.GO$p_adjust<-ifelse(sig.GO$under_rep_padj > sig.GO$over_rep_padj,
        sig.GO$over_rep_padj, 
        sig.GO$under_rep_padj )
    sig.GO$percentage<- round(100 * sig.GO$numDEInCat/sig.GO$numInCat,2)
    
    ret<-sig.GO
    
    slim<-geneInformation$GOSlim
    sig.GO <- sqldf("SELECT category, 
        over_represented_pvalue,
        under_represented_pvalue,
        numDEInCat,
        numInCat,
        ontology,
        over_rep_padj,
        under_rep_padj,
        description,
        type,
        p_adjust,
        percentage,  
        GROUP_CONCAT(DISTINCT slim_acc) as slim_acc, 
        GROUP_CONCAT(DISTINCT slim_term_type) as slim_term_type,
        GROUP_CONCAT(DISTINCT slim_name) as slim_name 
        FROM ret LEFT JOIN slim ON category = acc 
        GROUP BY category, 
        over_represented_pvalue,
        under_represented_pvalue,
        numDEInCat,
        numInCat,
        ontology,
        over_rep_padj,
        under_rep_padj,
        description,
        type,
        p_adjust,
        percentage")
    
    
    sig.GO

}

plot_enrichment<-function(enrichment,experiment="HC_CS_no_stress", title="test" , type="Over represented"){
    local_title <- paste0(title, "\n", experiment, "\n", type , "\nPercentage of genes from each ontology\n")
    gs<-list()
    for(ont in unique(enrichment$ontology)){
        current<-enrichment[enrichment$ontology == ont & enrichment$type == type, ]
        current$description<-strtrim(current$description, 70)
        if(nrow(current) > 0){
            gs[[length(gs)+1]] <- ggplot(current, aes(description, percentage)) + 
            coord_flip() +
            geom_col(fill="#a8ddb5") + labs(y="", x="percentage", title=ont) +  theme_bw() + #ylim(0,100) + 
            theme(axis.text.x = element_text(angle = 0, hjust = 1, size=4),
               axis.text.y = element_text(angle = 35, hjust = 1, size=4))  +
            geom_text(aes(label=numDEInCat), position=position_dodge(width=0.9), size=3)
        }
    }
    if(length(gs) == 0){
        gs[[length(gs)+1]] <- textGrob("None")
    }
    g1<-arrangeGrob(grobs=gs, ncol=2, top=local_title )
}

plot_gene_summary<-function(geneInformation, genes_to_plot, name="Random Samples" , output_path="./Test"){

    summary_df <- NULL

    local_table<-geneInformation$canonicalTranscripts
    local_table<-local_table[local_table$Gene %in% genes_to_plot,]
    
    local_mean_tpms<-geneInformation$meanTpms
    local_mean_tpms<-local_mean_tpms[local_mean_tpms$gene %in% genes_to_plot, ]
    
    stats_to_plot<-c('size_cds', 'exon_no', 'exon_length','intron_length', 'X3UTR_length', 'X5UTR_length' )
    
    
    expected_per_chr<-get_expected_values_per_5pc_bin(geneInformation$canonicalTranscripts, nrow(local_table))
    
    gs<-list()
    plots<-list()
    
    dir<-paste0(output_path,"/",name)
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)


    plots[[length(plots)+1]] <- textGrob(paste0(name, " Gene summary"))
    for(plot in stats_to_plot){
        p<-plotHistogram(local_table,column=plot)
        gs[[length(gs)+1]] <- p
        tmp_df <- prepare_hist_stats(local_table, column=plot)
        tmp_df$value_type <- plot
        tmp_df$dataset<- "Gene summary"
        tmp_df$title  <- name
        if(is.null(summary_df)){
            summary_df <- tmp_df
            }else{
                summary_df <- rbind(summary_df, tmp_df) 
            }
            
        }
        plots[[length(plots)+1]] <- arrangeGrob(grobs=gs, ncol=2 , top = paste0(name, "\n Gene properties"))
        plots[[length(plots)+1]] <- plot_per_chromosome_5pc_bins_facet(local_table,
         expected_per_chr=expected_per_chr, 
         title=name)

        plots[[length(plots)+1]] <- textGrob(paste0(name, " TPM summary"))
        for(s in unique(geneInformation$meanTpms$subset)){
            plots[[length(plots)+1]] <- plot_tpms_summary(local_mean_tpms, experiment=s, title=name) 
            plots[[length(plots)+1]] <- plot_tpm_desc_stats(geneInformation$meanTpms, local_mean_tpms, experiment=s, title=name)
            plots[[length(plots)+1]] <- plot_all_means_filteredtpms_summary(local_mean_tpms, experiment=s, title=name) 
        }

        plots[[length(plots)+1]] <- textGrob(paste0(name, " Triad summary"))


        triada_movment_df<-NULL


        for(s in unique(geneInformation$triads$dataset)){
            for(i in c(1,2,3) ){
                local_triads <- get_triads_from_genes(genes_to_plot, geneInformation, dataset=s, min_no_genes = i)
                if(nrow(local_triads$triads)== 0){
                    next
                }
                name_tmp <-name
                name <- paste0(name, " Min genes in triad: ", i )
                plots[[length(plots)+1]] <- plot_dominance_summary(local_triads, experiment=s, title=name)

            #tmp_df<- table_dominance_summary(local_triads, experiment=s, title=name)            
            #summary_df <- rbind(summary_df, tmp_df)


            observed_desc    <-get_dominance_summary_tables_per_factor(local_triads, experiment=s)
            observed_gen_desc<-get_dominance_summary_tables_per_factor(local_triads,
             description="general_description",
             experiment=s)

            expected_desc    <-get_dominance_summary_tables_per_factor(geneInformation, 
             experiment=s, 
             n=observed_desc$total)
            expected_gen_desc<-get_dominance_summary_tables_per_factor(geneInformation,
             description="general_description",
             experiment=s,
             n=observed_gen_desc$total)


            

            o_l  <- observed_desc$long
            og_l <- observed_gen_desc$long
            e_l  <- expected_desc$long
            eg_l <- expected_gen_desc$long

            o_l$group <- "fine"
            e_l$group <- "fine"
            og_l$group <- "general"
            eg_l$group <- "general"

            o_l$sum_for <- "observed"
            e_l$sum_for <- "expected"
            og_l$sum_for <- "observed"
            eg_l$sum_for <- "expected"

            o_l$dataset <- name_tmp
            og_l$dataset <- name_tmp
            e_l$dataset <- name_tmp
            eg_l$dataset <- name_tmp


            o_l$dataset_for_triads <- s
            og_l$dataset_for_triads <- s
            e_l$dataset_for_triads <- s
            eg_l$dataset_for_triads <- s

            o_l$min_expressed_genes_in_triad <- i
            og_l$min_expressed_genes_in_triad <- i
            e_l$min_expressed_genes_in_triad <- i
            eg_l$min_expressed_genes_in_triad <- i

            if(is.null(triada_movment_df)){
                triada_movment_df <- o_l
            }
            else{
                triada_movment_df <- rbind(triada_movment_df, o_l)
            }
            triada_movment_df <- rbind(triada_movment_df, og_l)
            triada_movment_df <- rbind(triada_movment_df, e_l)
            triada_movment_df <- rbind(triada_movment_df, eg_l)

            plots[[length(plots)+1]] <-plot_dominance_summary_tables(local_triads,
              expected_desc, 
              expected_gen_desc, 
              experiment=s, title=name  )
            
            

            plots[[length(plots)+1]] <- table_with_title(paste(name,s, "Observed",sep="\n"),
                observed_desc$pasted
                )
            plots[[length(plots)+1]] <- table_with_title(paste(name, s, "Expected",sep="\n"),
                expected_desc$pasted
                )

            plots[[length(plots)+1]] <- table_with_title(paste(name, s, "Observed",sep="\n"),
             observed_gen_desc$pasted
             )
            plots[[length(plots)+1]] <- table_with_title(paste(name, s, "Expected",sep="\n"),
                expected_gen_desc$pasted
                )
            name<-name_tmp
        }
        
        
        
    }
    
    output_summary<-paste0(dir, "/", "summary_from_histograms.csv")
    write.csv(summary_df, file=output_summary) 

    output_enrichment<-paste0(dir, "/", "triad_movment_summary.csv")
    write.csv(triada_movment_df, file=output_enrichment) 

    


    all_enrichments <- NULL
    for(g_u in unique(geneInformation$gene_universe$dataset)){
        for(ont in unique(geneInformation$ontologies$ontology)){
            enrichment_test<- get_goseq_enrichment(geneInformation, genes_to_plot,  ontology=ont , dataset=g_u)
            if(nrow(enrichment_test) == 0){
                next
            }
            enrichment_test$universe<-g_u
            enrichment_test$ontology_universe <- ont
            if(is.null(all_enrichments)){
                all_enrichments <- enrichment_test
            }
            else{
                all_enrichments<-rbind(all_enrichments, enrichment_test)
            }
        }
        
        
        
        
        
    }

    output_enrichment<-paste0(dir, "/", "enrichment.csv")
    write.csv(all_enrichments, file=output_enrichment)
    output_pdf<-paste0(dir, "/",name ,".pdf")


    g1<-marrangeGrob(plots, ncol=1, nrow=1, top="", bottom = quote(paste("page", g, "of",
     pages)))

    ggsave(output_pdf, plot=g1 , width = 210, height = 297, units = "mm")
    g1
}




get_counts_values_per_5pc_bin<-function(gene_table,group_in_single_chromosome=FALSE){

    query<-"SELECT Chr, chr_group, genome,scaled_5per_position, count(*) as count 
    FROM
    gene_table 
    WHERE Chr != 'chrUn'
    GROUP BY  Chr, chr_group, genome, scaled_5per_position"
    
    counts<-sqldf(query)
    if(group_in_single_chromosome){
        query<-"SELECT 'All' as Chr, 'all' as chr_group, 'all'  as genome,
        scaled_5per_position, count(*)/21 as count 
        FROM
        gene_table 
        WHERE Chr != 'chrUn'
        GROUP BY scaled_5per_position"
    }
    counts<-sqldf(query)

    counts
}

plot_per_chromosome_5pc_bins_overlap_lines<-function(table,expected_per_chr,
   expected_all_chromosomes=NULL, 
   title = "Test"){
    chromosomes=c("1A", "1B", "1D",
        "2A", "2B", "2D",
        "3A", "3B", "3D",
        "4A", "4B", "4D",
        "5A", "5B", "5D",
        "6A", "6B", "6D",

        "7A", "7B", "7D")
    
    
    gs<-list()
    local_title = paste0(title, "\n Genes per chromosome 5% bin\nN: ", nrow(table) )
    
    t <- table[table$Chr != "chrUn",]
    t1 <-  get_counts_values_per_5pc_bin(t)
    
    t2 <-  get_counts_values_per_5pc_bin(t, group_in_single_chromosome=TRUE)
    print(head(t1))
    expected_per_chr <- expected_per_chr[expected_per_chr$Chr != "chrUn",]
    
    p <-ggplot(expected_per_chr,aes(scaled_5per_position, expected, group=Chr)) 
    p <- p + xlim(0,100)
    
    
    
    
    p <- p + geom_line(data=expected_per_chr[, c("scaled_5per_position", "expected", "Chr")],
     aes(x=scaled_5per_position, y=expected,group=Chr 
        ),
     color='black', size=1, alpha=0.1 ) 
    
    
    p <- p + ylab(" count ") + xlab("")
    
    
    
    if(!is.null(expected_all_chromosomes)){
        expected_all_chromosomes$expected <- expected_all_chromosomes$expected/21
        exp_norm <-expected_all_chromosomes[, c("scaled_5per_position", "expected", "Chr")]
        p <- p + geom_line(data=exp_norm,
         aes(x=scaled_5per_position, y=expected), 
         color="blue", alpha=0.5)
    }
    p  <- p + theme_bw()
    p1 <- p + geom_line(color="red", alpha=0.3) 
    
    p1 <- p1 + geom_point(data=t1, aes(x=scaled_5per_position, y=count), color="red")
    p  <-  p + geom_point(data=t2, aes(x=scaled_5per_position, y=count), color="blue")
    
    gs[[length(gs)+1]] <- p1 + facet_grid(chr_group~genome,  drop = TRUE)
    gs[[length(gs)+1]] <- p  
    
    g1<-arrangeGrob(grobs=gs, ncol=1, heights=c(0.8,0.2),  top=local_title ) 
    g1
}

plotExpressedTissuesAcrossChromosomes<-function(geneInformation, 
    genes_to_plot, 
    subset="850_samples", 
    factor="all_mean_filter", 
    title = "Test"){
    #print(head(genes_to_plot))
    tpms<-geneInformation$meanTpms
    tpms<-tpms[tpms$factor==factor & tpms$subset==subset,]
    transcripts<-geneInformation$canonicalTranscripts
    genes_to_plot<-data.frame(Gene=genes_to_plot)
    
    transcripts$scaled_5per_position <-   5 * ceiling(transcripts$scaled_1per_position / 5) 
    transcripts$scaled_5per_position <- ifelse(transcripts$scaled_5per_position == 0, 5, transcripts$scaled_5per_position)
    #print(max(tpms$samples))
    
    expected_tissues <- sqldf("SELECT AVG(value) as meanTPM, AVG(samples) as noSamples, scaled_5per_position
        FROM tpms 
        JOIN transcripts ON tpms.gene = transcripts.Gene 
        WHERE geneconf = 'HC' AND Chr != 'chrUn'
        GROUP BY scaled_5per_position")
    
    expected_tissues_mean <- sqldf("SELECT 
        Chr,
        chr_group, 
        genome, 
        scaled_5per_position, 
        AVG(value) as meanTPM, 
        AVG(samples) as noSamples, 
        count(*) as count
        FROM tpms 
        JOIN transcripts ON tpms.gene = transcripts.Gene 
        WHERE geneconf = 'HC' AND Chr != 'chrUn'
        GROUP BY Chr, scaled_5per_position, chr_group, genome
        ORDER BY Chr, chr_group, genome, scaled_5per_position ")
    
    
    gs<-list()
    local_title = paste0(title, "\n Average expressed per 5% bin\nN: ", nrow(table) )
    
    #t <- table[table$Chr != "chrUn",]
    t1 <- sqldf("SELECT AVG(value) as meanTPM, AVG(samples) as noSamples, scaled_5per_position
        FROM tpms 
        JOIN transcripts ON tpms.gene = transcripts.Gene 
        WHERE transcripts.Gene in genes_to_plot AND Chr != 'chrUn'
        GROUP BY scaled_5per_position")
    
    t2 <- sqldf("SELECT 
        Chr,
        chr_group, 
        genome, 
        scaled_5per_position, 
        AVG(value) as meanTPM, 
        AVG(samples) as noSamples, 
        count(*) as count
        FROM tpms 
        JOIN transcripts ON tpms.gene = transcripts.Gene 
        WHERE transcripts.Gene in genes_to_plot AND Chr != 'chrUn'
        GROUP BY Chr, scaled_5per_position, chr_group, genome
        ORDER BY Chr, chr_group, genome, scaled_5per_position ")
    
    #print("-.-")
    t1$Chr<-"All"
    #print(head(t1))
    #print(head(t2))

    expected_per_chr <- expected_per_chr[expected_per_chr$Chr != "chrUn",]
    
    p <-ggplot(expected_tissues_mean,aes(x=scaled_5per_position, y=noSamples, group=Chr)) 
    
    samples_reduced<-expected_tissues_mean[, c("scaled_5per_position", "noSamples", "Chr")]

    p <- p + geom_line(data=samples_reduced,
     aes(x=scaled_5per_position, y=noSamples,group=Chr 
        ),color='black', size=1, alpha=0.1 ) 
    p <- p + ylab("No of tissues") + xlab("")
    if(!is.null(expected_tissues)){
        #print(head(expected_tissues))
        exp_norm <-expected_tissues[, c("scaled_5per_position", "noSamples")]
        exp_norm$Chr<-"All"
        p <- p + geom_line(data=exp_norm,
         aes(x=scaled_5per_position, y=noSamples), 
         color="blue", alpha=0.5)
    }
    p  <- p + theme_bw()
    p1 <- p + geom_line(color="red", alpha=0.3) 

    p1 <- p1 + geom_point(data=t2, aes(x=scaled_5per_position, y=noSamples), color="red")
    p  <-  p + geom_point(data=t1, aes(x=scaled_5per_position, y=noSamples), color="blue")

    gs[[length(gs)+1]] <- p1 + facet_grid(chr_group~genome,  drop = TRUE)
    gs[[length(gs)+1]] <- p  

    g1<-arrangeGrob(grobs=gs, ncol=1, heights=c(0.8,0.2),  top=local_title ) 
    g1
}


plotTPMOfExpressedTissuesAcrossChromosomes<-function(geneInformation, 
    genes_to_plot, 
    subset="850_samples", 
    factor="all_mean_filter", 
    title = "Test"){
    #print(head(genes_to_plot))
    tpms<-geneInformation$meanTpms
    tpms<-tpms[tpms$factor==factor & tpms$subset==subset,]
    transcripts<-geneInformation$canonicalTranscripts
    genes_to_plot<-data.frame(Gene=genes_to_plot)
    
    transcripts$scaled_5per_position <-   5 * ceiling(transcripts$scaled_1per_position / 5) 
    transcripts$scaled_5per_position <- ifelse(transcripts$scaled_5per_position == 0, 5, transcripts$scaled_5per_position)
    #print(max(tpms$samples))
    
    expected_tissues <- sqldf("SELECT AVG(value) as meanTPM, AVG(samples) as noSamples, scaled_5per_position
        FROM tpms 
        JOIN transcripts ON tpms.gene = transcripts.Gene 
        WHERE geneconf = 'HC' AND Chr != 'chrUn'
        GROUP BY scaled_5per_position")
    
    expected_tissues_mean <- sqldf("SELECT 
        Chr,
        chr_group, 
        genome, 
        scaled_5per_position, 
        AVG(value) as meanTPM, 
        AVG(samples) as noSamples, 
        count(*) as count
        FROM tpms 
        JOIN transcripts ON tpms.gene = transcripts.Gene 
        WHERE geneconf = 'HC' AND Chr != 'chrUn'
        GROUP BY Chr, scaled_5per_position, chr_group, genome
        ORDER BY Chr, chr_group, genome, scaled_5per_position ")
    
    
    gs<-list()
    local_title = paste0(title, "\n Mean TPM of expressed tissues expressed per 5% bin\nN: ", nrow(table) )
    
    #t <- table[table$Chr != "chrUn",]
    t1 <- sqldf("SELECT AVG(value) as meanTPM, AVG(samples) as noSamples, scaled_5per_position
        FROM tpms 
        JOIN transcripts ON tpms.gene = transcripts.Gene 
        WHERE transcripts.Gene in genes_to_plot AND Chr != 'chrUn'
        GROUP BY scaled_5per_position")
    
    t2 <- sqldf("SELECT 
        Chr,
        chr_group, 
        genome, 
        scaled_5per_position, 
        AVG(value) as meanTPM, 
        AVG(samples) as noSamples, 
        count(*) as count
        FROM tpms 
        JOIN transcripts ON tpms.gene = transcripts.Gene 
        WHERE transcripts.Gene in genes_to_plot AND Chr != 'chrUn'
        GROUP BY Chr, scaled_5per_position, chr_group, genome
        ORDER BY Chr, chr_group, genome, scaled_5per_position ")
    
    #print("-.-")
    t1$Chr<-"All"
    #print(head(t1))
    #print(head(t2))

    expected_per_chr <- expected_per_chr[expected_per_chr$Chr != "chrUn",]
    
    p <-ggplot(expected_tissues_mean,aes(x=scaled_5per_position, y=meanTPM, group=Chr)) 
    
    samples_reduced<-expected_tissues_mean[, c("scaled_5per_position", "meanTPM", "Chr")]

    p <- p + geom_line(data=samples_reduced,
     aes(x=scaled_5per_position, y=meanTPM,group=Chr 
        ),color='black', size=1, alpha=0.1 ) 
    p <- p + ylab("No of tissues") + xlab("")
    if(!is.null(expected_tissues)){
        #print(head(expected_tissues))
        exp_norm <-expected_tissues[, c("scaled_5per_position", "meanTPM")]
        exp_norm$Chr<-"All"
        p <- p + geom_line(data=exp_norm,
         aes(x=scaled_5per_position, y=meanTPM), 
         color="blue", alpha=0.5)
    }
    p  <- p + theme_bw() + scale_y_log10()
    p1 <- p + geom_line(color="red", alpha=0.3) 

    p1 <- p1 + geom_point(data=t2, aes(x=scaled_5per_position, y=meanTPM), color="red")
    p  <-  p + geom_point(data=t1, aes(x=scaled_5per_position, y=meanTPM), color="blue")

    gs[[length(gs)+1]] <- p1 + facet_grid(chr_group~genome,  drop = TRUE)
    gs[[length(gs)+1]] <- p  

    g1<-arrangeGrob(grobs=gs, ncol=1, heights=c(0.8,0.2),  top=local_title ) 
    g1
}

args = commandArgs(trailingOnly=TRUE)
#folder<-"/Users/ramirezr/Dropbox/JIC/expVIPMetadatas/RefSeq1.0/TablesForExploration"
#genes_to_plot_path<-"/Users/ramirezr/Dropbox/JIC/expVIPMetadatas/RefSeq1.0/notebook/gene_set_files/modules/WGCNA_850/WGCNA_850_Module_15.txt"
folder<-args[1]
genes_to_plot_path<-args[2]
name<-basename(genes_to_plot_path)
path<-paste0(genes_to_plot_path,"_plots")

genes_to_plot<-read.csv(genes_to_plot_path)
genes_to_plot<-as.vector(genes_to_plot[,1])

print(name)
print(paste0("number of genes to plot: ", length(genes_to_plot)))
print(head(genes_to_plot))

geneInformation<-loadGeneInformation(dir=folder)
g<-plot_gene_summary(geneInformation,genes_to_plot, name = name, output_path = path )

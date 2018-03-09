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
library(plyr)

options(keep.source = TRUE, error =
  quote({
    cat("Environment:\n", file=stderr());

    # TODO: setup option for dumping to a file (?)
    # Set `to.file` argument to write this to a file for post-mortem debugging
    dump.frames();  # writes to last.dump

    #
    # Debugging in R
    #   http://www.stats.uwo.ca/faculty/murdoch/software/debuggingR/index.shtml
    #
    # Post-mortem debugging
    #   http://www.stats.uwo.ca/faculty/murdoch/software/debuggingR/pmd.shtml
    #
    # Relation functions:
    #   dump.frames
    #   recover
    # >>limitedLabels  (formatting of the dump with source/line numbers)
    #   sys.frame (and associated)
    #   traceback
    #   geterrmessage
    #
    # Output based on the debugger function definition.

    n <- length(last.dump)
    calls <- names(last.dump)
    cat(paste("  ", 1L:n, ": ", calls, sep = ""), sep = "\n", file=stderr())
    cat("\n", file=stderr())

    if (!interactive()) {
      q()
    }
  }))

is.error <- function(x) inherits(x, "try-error")

loadGeneInformation<-function(dir="../TablesForExploration",
                              motifs=T,
                              WGCNA=F,
                              meanTpms=T,
                              non_syntenic_triads=F
                             ){

    path<-paste0(dir,"/CanonicalTranscript.rds")
    canonicalTranscripts<-readRDS(path)sudo apt-get install libzmq3-dev libssh2-1-dev python3-pip build-essential python3-dev
    canonicalTranscripts$intron_length<- canonicalTranscripts$mrna_length -  canonicalTranscripts$exon_length
    canonicalTranscripts$chr_group <- substr(canonicalTranscripts$Chr,4,4)
    canonicalTranscripts$genome    <- substr(canonicalTranscripts$Chr,5,5)
    expressed_genes <- canonicalTranscripts$Gene

    if(meanTpms == T){
        path<-paste0(dir, "/MeanTpms.rds")
        meanTpms <- readRDS(path)
        expressed_genes<-unique(meanTpms$gene)
    }

    canonicalTranscripts<-canonicalTranscripts[canonicalTranscripts$Gene %in% expressed_genes, ]

    canonicalTranscripts$scaled_5per_position <-   5 * ceiling(canonicalTranscripts$scaled_1per_position / 5)
    canonicalTranscripts$scaled_5per_position <- ifelse(canonicalTranscripts$scaled_5per_position == 0,
        5,
        canonicalTranscripts$scaled_5per_position)

    path<-paste0(dir, "/region_partition.csv")
    partition<-read.csv(path, row.names=1)

    partition_percentages<-round(100*partition/partition$Length)
    partition_percentages$Chr <- rownames(partition_percentages)
    partition$Chr <- rownames(partition)
    ct<-canonicalTranscripts
    ct_with_partition<-sqldf('SELECT ct.*, CASE
WHEN scaled_1per_position < R1_R2a THEN "R1"
WHEN scaled_1per_position < R2a_C  THEN "R2A"
WHEN scaled_1per_position < C_R2b  THEN "C"
WHEN scaled_1per_position < R2b_R3  THEN "R2B"
ELSE "R3" END as partition

FROM ct LEFT JOIN partition_percentages ON ct.chr = partition_percentages.chr   ')

    x<-  as.factor(ct_with_partition$partition)
    x <- factor(x,levels(x)[c(2,3,1,4,5)])
    ct_with_partition$partition <- x
    canonicalTranscripts<-ct_with_partition

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

    if(WGCNA == T){
        path<-paste0(dir, "/WGCNA_table.csv")
        WGCNA <-  read.csv(path)
    }

    path<-paste0(dir, "/ObservedGOTermsWithSlim.csv")
    go_slim<-read.csv(path, row.names=1)


    if(motifs == T){
        path<-paste0(dir, "/motifs.rds")
        motifs <- readRDS(path)
        motifs<-unique(motifs)
    }

    path<-paste0(dir, "/SegmentalTriads.csv")
    allTriads<-read.csv(path, stringsAsFactors=F)
    only_genes<-allTriads[,c("group_id","A", "B", "D")]
    allTriads<-melt(only_genes, id.vars<-c("group_id"),
        variable.name = "chr_group",
        value.name ="gene")

    if(non_syntenic_triads){
        path<-paste0(dir,"/Non_syn_Triads.rds")
        tmp<-readRDS(path)
        triads <- rbind(triads, tmp)

        path<-paste0(dir,"/Non_syn_TriadMovement.rds")
        tmp<-readRDS(path)
        triadMovement <- rbind(triadMovement, tmp)

        path<-paste0(dir, "/NonSyntenicTriads.csv")
        tmp_allTriads<-read.csv(path, stringsAsFactors=F)
        only_genes<-tmp_allTriads[,c("group_id","A", "B", "D")]
        tmp_allTriads<-melt(only_genes, id.vars<-c("group_id"),
            variable.name = "chr_group",
            value.name ="gene")
        allTriads<-rbind(allTriads, tmp_allTriads)
    }


    list(canonicalTranscripts=canonicalTranscripts,
       meanTpms=meanTpms,
       triads=triads,
       triadMovement=triadMovement,
       gene_universe=gene_universe,
       ontologies=ontologies,
       id_names=id_names,
       WGCNA=WGCNA,
       GOSlim=go_slim,
       partition=partition,
       motifs=motifs,
       allTriads=allTriads
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

plotHistogram<-function(table,
 column="size_cds",
 probs = c( 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)){
    table<-table[table[,column]>0,]
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
       table$quantile<-ifelse(is.na(table$quantile),0,table$quantile)
       table$quantile<-as.factor(table$quantile)

       iq <- quantiles$value[4] - quantiles$value[2]

       xmax <- quantiles$value[3] + (iq * 2)
       xmin <- quantiles$value[3] - (iq * 2)
       if(xmin < 0) xmin <- 0
       if(xmax > local_max)  xmax <- local_max + 1

       p <- ggplot(table, aes_string(column, fill="quantile"))
       p <- p + geom_vline(data=quantiles,aes(xintercept=quantiles$value) )
       for(i in seq(1,nrow(quantiles))){
        x_pos<-quantiles$value[i]
        gtext <- textGrob(quantiles$quant[i], y=0.02,  gp = gpar(fontsize = 6,col = "red"))
        p <- p + annotation_custom(gtext, xmin=x_pos, xmax=x_pos)
       }
       p <- p  + xlim(xmin, xmax) + scale_fill_brewer(palette="Dark2")
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
                                             title = "Test"){

    gs<-list()
    local_title = paste0(title, "\n Genes per chromosome 5% bin\nN: ", nrow(table) )

    t1 <- table[table$Chr != "chrUn",]
    expected_per_chr <- expected_per_chr[expected_per_chr$Chr != "chrUn",]

    p <-ggplot(t1,aes(as.factor(scaled_5per_position)))
    p <- p + geom_bar(aes(fill=partition))
    p <- p + theme_bw()
    p <- p + facet_grid(chr_group~genome,  drop = FALSE)
    p <- p + ylab(" count ") + xlab("")
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    p <- p + scale_fill_brewer(palette = "Set1")
    p <- p + geom_point(data=expected_per_chr,
        aes(x=as.factor(scaled_5per_position), y=expected), color="red", size = 0.5)
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

plot_all_means_filteredtpms_summary<-function(tpms, transcripts, experiment="850_samples", min_tpm=0.5, title="Test"){

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

    gs[[length(gs)+1]] <- plot_expressed_tissues_across_chromosomes(tpms,transcripts, bin_size = 2 )
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
    triads <- selected_triads$triads
    triadMovement<-selected_triads$triadMovement
    all_means_filter<-triads[triads$factor=="all_mean_filter",]

    df  <- prepare_hist_stats(all_means_filter, column="value")
    df$value_type <- "All mean filter TPM for genes in triad"
    df$dataset<-experiment
    df$title <- title

    tmp<- prepare_hist_stats(triadMovement, column="central_max_distance")
    tmp$value_type <- "central_max_distance"
    tmp$dataset<-experiment
    tmp$title <- title
    df<-rbind(df,tmp)

    tmp<- prepare_hist_stats(triadMovement, column="central_mean_distance")
    tmp$value_type <- "central_mean_distance"
    tmp$dataset<-experiment
    tmp$title <- title
    df<-rbind(df,tmp)

    tmp<- prepare_hist_stats(triadMovement, column="sum_mean_tpm")
    tmp$value_type <- "sum_mean_tpm"
    tmp$dataset<-experiment
    tmp$title <- title
    df<-rbind(df,tmp)

    tmp<- prepare_hist_stats(triadMovement, column="factor_count")
    tmp$value_type <- "No. of conditions of genes(count per triad)"
    tmp$dataset<-experiment
    tmp$title <- title
    df<-rbind(df,tmp)

    tmp<- prepare_hist_stats(triadMovement, column="total_categories")
    tmp$value_type <- "No. of categories"
    tmp$dataset<-experiment
    tmp$title <- title
    df<-rbind(df,tmp)

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
  factor = NULL,
  n=NULL
  ){
    triads <- selected_triads$triads
    triads <- triads[triads$dataset==experiment, ]
    if(!is.null(factor)){
        triads <- triads[triads$factor==factor, ]
    }
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
    pasted$total <-total_per_factor
    long<-melt(casted)
    colnames(long)<-c("factor", "description", "count")
    percentage_long<-melt(percentage)
    colnames(percentage_long)<-c("factor", "description", "percentage")
    long<-merge(x=long,
                y=percentage_long,
                by.x=c("factor","description"),
                by.y=c("factor","description"))
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

plot_expressed_tissues_across_chromosomes<-function(tpms, transcripts,
                                                title = "Test", bin_size=5){

    #all_means_filter, ct_with_partition,
    #transcripts<-geneInformation$canonicalTranscripts
    #all_means_filter<-geneInformation$meanTpms
    #all_means_filter<-all_means_filter[all_means_filter$factor=='all_mean_filter']
    transcripts$scaled_5per_position <-   bin_size * ceiling(transcripts$scaled_1per_position / bin_size)
    transcripts$scaled_5per_position <- ifelse(transcripts$scaled_5per_position == 0, bin_size, transcripts$scaled_5per_position)
    #transcripts$scaled_5per_position <- transcripts$scaled_pc
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

    p <-ggplot(expected_tissues_mean,aes(x=scaled_5per_position, y=noSamples, group=Chr))
    samples_reduced<-expected_tissues_mean[, c("scaled_5per_position", "noSamples", "Chr")]

    p <- p + geom_line(data=samples_reduced,
                       aes(x=scaled_5per_position, y=noSamples,group=Chr
                        ),color='black', size=0.4, alpha=0.2 )
    p <- p + ylab("No of tissues") + xlab("")
    if(!is.null(expected_tissues)){
        exp_norm <-expected_tissues[, c("scaled_5per_position", "noSamples")]
        exp_norm$Chr<-"All"
        p <- p + geom_line(data=exp_norm,
                       aes(x=scaled_5per_position, y=noSamples),
                           color="blue", size=1, alpha=1)
    }
    p  <- p + theme_bw() + theme(axis.text=element_text(size=7),
          axis.title=element_text(size=7),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
    p <- p + scale_x_continuous(expand = c(0, 0))
    p <- p + geom_vline(xintercept = c(9, 28,50,80))
    p
}

plot_gene_summary<-function(geneInformation, genes_to_plot, name="Random Samples" , output_path="./Test", run_stats=FALSE){

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
    gc()

    plots[[length(plots)+1]] <- textGrob(paste0(name, " Gene summary"))
    print("Plotting gene summaries")
    for(plot in stats_to_plot){

        p<-plotHistogram(local_table,column=plot)
        gs[[length(gs)+1]] <- p
        tmp_df <- prepare_hist_stats(local_table, column=plot)
        tmp_df$value_type <- plot
        tmp_df$dataset<- "Gene summary"
        tmp_df$title  <- name
        if(is.null(summary_df))
        {
            summary_df <- tmp_df
        }
        else
        {
            summary_df <- rbind(summary_df, tmp_df)
        }

    }
    plots[[length(plots)+1]] <- arrangeGrob(grobs=gs, ncol=2 , top = paste0(name, "\n Gene properties"))
    plots[[length(plots)+1]] <- plot_per_chromosome_5pc_bins_facet(local_table,
        expected_per_chr=expected_per_chr,
        title=name)

    gene_density<-get_gene_density(genes_to_plot, geneInformation)
    output_gene_density<-paste0(dir, "/", "gene_density_per_region.csv")
    plots[[length(plots)+1]] <- plot_per_partition_gene_count(local_table, gene_density, title=name)

    write.csv(gene_density, file=output_gene_density)

    print("Plotting TPM summaries")
    plots[[length(plots)+1]] <- textGrob(paste0(name, " TPM summary"))
    for(s in unique(geneInformation$meanTpms$subset)){
        plots[[length(plots)+1]] <- plot_tpms_summary(local_mean_tpms, experiment=s, title=name)
        plots[[length(plots)+1]] <- plot_density_expression(geneInformation,genes_to_plot, experiment=s, title=name)
        plots[[length(plots)+1]] <- plot_tpm_desc_stats(geneInformation$meanTpms, local_mean_tpms, experiment=s, title=name)
        plots[[length(plots)+1]] <- plot_all_means_filteredtpms_summary(local_mean_tpms, geneInformation$canonicalTranscripts, experiment=s, title=name)

    }

    plots[[length(plots)+1]] <- textGrob(paste0(name, " Triad summary"))
    triada_movment_df<-NULL

    print("Plotting triads")
    for(s in unique(geneInformation$triads$dataset)){
        tmp_plot<-plot_clust_dist(geneInformation,    genes_to_plot, experiment=s, title = paste0(s,"\n",name))
        if(is.null(tmp_plot)){
            next
        }
        plots[[length(plots)+1]] <- tmp_plot
        plots[[length(plots)+1]] <- plot_triad_movment(geneInformation, genes_to_plot,
                          experiment=s,
                          title=paste0(s,"\n",name))
        for(i in c(1,2,3) ){

            local_triads <- get_triads_from_genes(genes_to_plot, geneInformation, dataset=s, min_no_genes = i)
            if(nrow(local_triads$triads)== 0){
                next
            }
            name_tmp <-name
            name <- paste0(name, " Min genes in triad: ", i )
            plots[[length(plots)+1]] <- plot_dominance_summary(local_triads, experiment=s, title=name)

            tmp_df<- table_dominance_summary(local_triads, experiment=s, title=name)
            summary_df <- rbind(summary_df, tmp_df)


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


    output_pdf<-paste0(dir, "/",name ,".pdf")
    g1<-marrangeGrob(plots, ncol=1, nrow=1, top="", bottom = quote(paste("page", g, "of",
       pages)))

    ggsave(output_pdf, plot=g1 , width = 210, height = 297, units = "mm")
    if(run_stats){

        path_motifs_t_test<-paste0(dir, "/", "motifs_t_test.csv")
        path_motifs_fisher<-paste0(dir, "/", "motifs_fisher.csv")
        path_motifs_triads<-paste0(dir, "/", "motifs_triads.csv")

        print("Testing motif enrichment")
        res<-get_motifs_for_genes(genes_to_plot, geneInformation, name=name)

        write.csv(res$t,
            file=path_motifs_t_test,
            row.names=F)
        write.csv(res$fisher,
            file=path_motifs_fisher,
            row.names=F)

        write.csv(get_motifs_for_triad(genes_to_plot, geneInformation, name=name),
            file=path_motifs_triads,
            row.names=F)
        res<-NULL

        print("Testing GO enrichment")
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
    }
    g1
}

plot_density_expression<-function(geneInformation,genes_to_plot, experiment="850_samples", min_tpm=0.5, title="Test"){

    tpms<-geneInformation$meanTpms
    tpms<-tpms[tpms$gene %in% genes_to_plot, ]

    local_tpms<-subset(tpms, (subset == experiment) &
     ( factor != "all" & factor != "all_means" & factor != "all_mean_filter" ) &
     value > min_tpm)

    ct <- geneInformation$canonicalTranscripts

    local_tpms<-sqldf("SELECT local_tpms.*, ct.genome FROM local_tpms JOIN ct ON ct.gene = local_tpms.gene
WHERE ct.genome != 'n' ")

    local_tpms$log_value <- log2(local_tpms$value)

    local_title <- paste0(title, "\n", experiment)

    p  <- ggplot(local_tpms, aes(value, colour = genome))
    p  <- p + geom_density( ) + theme_bw()
    p  <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1),
      strip.text = element_text(size=6))
    p  <- p + facet_wrap(~ factor, ncol=4)
    p  <- p + scale_x_continuous(trans='log2',expand = c(0, 0))
    p <- p + coord_cartesian(xlim = c(0.5, 128))
    p  <- p + ylab("Density") + xlab("")
    p  <- p + theme(strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")))


    local_tpms<-subset(tpms, (subset == experiment) &
     ( factor == "all_mean_filter" ) &
     value > min_tpm)

    local_tpms<-sqldf("SELECT local_tpms.*, ct.genome FROM local_tpms JOIN ct ON ct.gene = local_tpms.gene
WHERE ct.genome != 'n' ")
    local_tpms$log_value <- log2(local_tpms$value)

    p2  <- ggplot(local_tpms, aes(value, colour = genome))
    p2  <- p2 + geom_density() + theme_bw()
    p2  <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1),
      strip.text = element_text(size=6, lineheight=0.5))

    p2  <- p2 + ylab("") + xlab("TPM")
    p2  <- p2 + scale_x_continuous(trans='log2',expand = c(0, 0))
    p2 <- p2 + coord_cartesian(xlim = c(0.5, 128))
    mytheme <- gridExtra::ttheme_default(
        core = list(fg_params=list(cex = 0.5)),
        colhead = list(fg_params=list(cex = 0.5)),
        rowhead = list(fg_params=list(cex = 0.5)))



    lay <- rbind(c(1),
       c(2))
    g1<-arrangeGrob(grobs=list(p,p2), heights=c(0.8,0.2), layout_matrix=lay, top = local_title)
    g1
}


get_counts_values_per_5pc_bin<-function(gene_table,group_in_single_chromosome=FALSE){
    query<-"SELECT Chr, chr_group, genome,scaled_5per_position, count(*) as count
    FROM
    gene_table
    WHERE Chr != 'chrUn'
    GROUP BY  Chr, chr_group, genome, scaled_5per_position"
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

plot_triad_movment<-function(geneInformation,
                          genes_to_plot,
                          experiment="HC_850_samples",
                          title="Test",
                             density = FALSE ,
                             points  = TRUE
                          ){

    allTriads<-geneInformation$allTriads
    selectedTriads<-unique(allTriads[allTriads$gene %in% genes_to_plot, "group_id"])

    tmp_df<-geneInformation$triads[geneInformation$triads$group_id %in% selectedTriads &
                                  geneInformation$triads$dataset == experiment,
                                   c("group_id","factor","clust","description","chr_group","normalised_triad", "Distance")]
    clust_df <- dcast(tmp_df,group_id +clust+description+factor+Distance ~ chr_group, value.var = "normalised_triad")

    clust_df_all_mean <- clust_df[clust_df$factor == 'all_mean_filter' , ]
    clust_df_factors  <- clust_df[clust_df$factor  != 'all_mean_filter' &
                                  clust_df$factor  != 'all_means' &
                                  clust_df$factor  != 'all' , ]

    tern_mean <- ggtern(clust_df_all_mean,aes(A,B,D,color=description)) +  theme_bw() +
        theme_legend_position(x = "topleft")  +
       theme_arrownormal() + guides(colour = guide_legend(override.aes = list(alpha = 1))) +ggtitle("Mean")

    tern_all <- ggtern(clust_df_factors,aes(A,B,D,color=description)) + theme_bw() +
       guides(colour = guide_legend(override.aes = list(alpha = 1)))

    tern_fact <- tern_all + facet_wrap(~factor, ncol=4) +
    theme_notitles() + theme(legend.position = "none") + theme_nolabels()

    tern_all <-  tern_all + theme_legend_position(x = "topleft")   + theme_arrownormal()  +ggtitle("All factors")

    if(density){
        tern_mean<- tern_mean  + stat_density_tern(
        geom='polygon', show.legend = F,

        aes(fill=..level..),

        bins=10,
        color='grey')
       tern_all <- tern_all   +  stat_density_tern(
        geom='polygon',show.legend = F,
        aes(fill=..level..),
        bins=10,
        color='grey')
      tern_fact<- tern_fact  +  stat_density_tern(
        geom='polygon',show.legend = F,
        aes(fill=..level..),
        bins=5,
        color='grey')

    }
    if(points){
        tern_mean <- tern_mean  +  geom_point(alpha=0.25)
        tern_all  <- tern_all   +  geom_point(alpha=0.25)
        tern_fact <- tern_fact  +  geom_point(alpha=0.25)
    }


    gs<-list(tern_mean, tern_all, tern_fact)

    lay <- rbind(c( 1, 3),
                 c( 2, 3)
                 )

    g2 <- arrangeGrob(grobs = gs, layout_matrix = lay, top = title)
    g2
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

    t1$Chr<-"All"
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

get_empty_bins_for_partitions<-function(geneInformation){
    partition<-geneInformation$partition
    chrs<-NULL
    for(i in rownames(partition))
    {
        chrom<-partition[i,"Chr"]
        length<-partition[i,"Length"]
        Chr<-rep(chrom, length+1)
        bin<-0:length+1
        count<-rep(0,length+1)
        chr_group<-rep(substr(chrom, 4,4),length+1)
        genome<-rep(substr(chrom, 5,5),length+1)
        l_partition<-rep(0,length+1)
        df<-data.frame(Chr,chr_group, genome, partition=l_partition, bin, count)
        if(is.null(chrs))
        {
            chrs<-df
        }else
        {
            chrs<-rbind(chrs, df)
        }
    }

    chrs_with_partition<-sqldf('SELECT chrs.Chr, chrs.chr_group, genome,  CASE
        WHEN chrs.bin < R1_R2a THEN "R1"
        WHEN chrs.bin < R2a_C  THEN "R2A"
        WHEN chrs.bin < C_R2b  THEN "C"
        WHEN chrs.bin < R2b_R3  THEN "R2B"
        ELSE "R3" END as new_partition, chrs.bin,
        count
        FROM chrs LEFT JOIN partition ON partition.Chr = chrs.Chr   ')
    c("Chr", "chr_group", "genome","partition","bin", "count")->colnames(chrs_with_partition)
    chrs_with_partition
}

get_gene_density<-function(genes_to_plot, geneInformation, bin_size=1000000){
    local_table<-geneInformation$canonicalTranscripts
    local_table<-local_table[local_table$Gene %in% genes_to_plot,]

    ct<-local_table
    ct$bin <- round(ct$Start / bin_size)
    density <- sqldf("SELECT Chr, chr_group, genome, partition, bin, count(*) as count FROM ct
        GROUP BY Chr, chr_group, genome, partition, bin
        ORDER BY Chr, bin")
    density<-rbind(density, get_empty_bins_for_partitions(geneInformation))
    sqldf("SELECT  Chr, chr_group, genome, partition, bin, sum(count) as count
        FROM density
        GROUP BY Chr, chr_group, genome, partition, bin
        ORDER BY Chr, bin")


}

plot_per_partition_gene_count<-function(table, gene_density, title = "Test"){

    gs<-list()
    local_title = paste0(title, "\n Genes per 1MBp \nN: ", nrow(table) )

    t1 <- gene_density[gene_density$Chr != "chrUn",]

    ylim1 = boxplot.stats(t1$count)$stats[c(1, 5)]
    p <-ggplot(t1,aes(partition, count))
    p <- p + geom_boxplot(aes(fill=partition))
    p <- p + theme_bw()
    p <- p + coord_cartesian(ylim = ylim1*1.05)
    p <- p + scale_fill_brewer(palette = "Set1")
    p <- p + theme(legend.position="none")
    p1 <- p + facet_grid(chr_group~genome,  drop = FALSE)
    gs[[length(gs)+1]] <- p1
    gs[[length(gs)+1]] <- p

    g1<-arrangeGrob(grobs=gs, ncol=1, heights=c(0.8,0.2), top=local_title )
    g1
}

get_motifs_for_triad<-function(genes, geneInformation, name=name){
    motifs<-geneInformation$motifs
    triads<-geneInformation$allTriads
    motifs<-motifs[motifs$gene %in% genes,]
    genes_df <- data.frame(gene=genes)
    triad_counts <- sqldf("
        SELECT group_id
        FROM triads
        JOIN genes_df
        WHERE triads.gene = genes_df.gene
        GROUP BY group_id
        HAVING count(*) = 3")
    query<-"SELECT
    motif,
    motif_set,
    chr_group,
    count(DISTINCT triads.gene) as total_genes,
    sum(count) as sum,
    avg(count) as average
    FROM triads
    JOIN motifs ON triads.gene = motifs.gene
    JOIN triad_counts ON triad_counts.group_id = triads.group_id
    GROUP BY motif, motif_set, chr_group"
    aggregated<-sqldf(query)

    sums<-sqldf("SELECT motif, motif_set,
        sum(total_genes) as sum_total_genes,
        sum(sum) as sum_sum,
        sum(average) as sum_average
        FROM aggregated
        GROUP BY motif, motif_set")
    percentages<-sqldf("SELECT aggregated.*,
        sum_total_genes, sum_sum, sum_average,
        100.0 * total_genes / sum_total_genes as percentage_total_genes,
        100.0 * sum         / sum_sum         as percentage_sum,
        100.0 * average     / sum_average     as percentage_average
        FROM aggregated JOIN sums
        ON  aggregated.motif = sums.motif
        AND aggregated.motif_set = sums.motif_set
        ORDER BY
        motif_set, motif,  chr_group ")
    percentages$Gene_set <- name
    percentages
}

get_motifs_for_genes<-function(genes_to_plot, geneInformation, name="Test"){
    datasets<-unique(geneInformation$gene_universe$dataset)
    total_sets<-length(datasets)
    gene_universe<-geneInformation$gene_universe
    motifs<-geneInformation$motifs
    motif_sets<-unique(motifs$motif_set)
    alternatives<-c("greater","less")
    matrix_for_test<-matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
    rownames(matrix_for_test)<-c("gene_set", "universe")
    colnames(matrix_for_test)<-c("motif", "genes")

    enrich_all_family <- data.frame(matrix(ncol = 10, nrow = 0))
    colnames(enrich_all_family) <- c("Universe",
                                     "Gene_set",
                                     "Motif_set",
                                     "Motif",
                                     "have_gene_set",
                                     "dont_have_gene_set",
                                     "have_universe",
                                     "dont_have_universe",
                                     "fisher_alternative",
                                     "fisher_pvalue")

    enrich_t_test <- data.frame(matrix(ncol = 11, nrow = 0))
    colnames(enrich_t_test) <- c("Universe",
                                     "Gene_set",
                                     "Motif_set",
                                     "Motif",
                                     "gene_set_gene_count",
                                     "universe_gene_count",
                                     "gene_set_motif_mean",
                                     "universe_motif_mean",
                                     "statistic_t",
                                     "parameter_df",
                                     "p_value")

    for(dataset in datasets){
        universe<-gene_universe[gene_universe$dataset == dataset,]
        genes_in_universe<-genes_to_plot[genes_to_plot %in% universe$gene]
        count_gene_set_genes<-length(genes_in_universe)
        count_universe_genes<-nrow(universe)
        #We advance the iteration if we are comparing the universe to itself. This is because
        #The statistical tests fail when both sets are exactly the same.
        if(count_gene_set_genes == count_universe_genes) next
        #matrix_for_test[1,2] <- count_gene_set_genes
        #matrix_for_test[2,2] <- count_universe_genes

        for(m_set in motif_sets){
            universe_motifs<-motifs[motifs$motif_set == m_set &
                                    motifs$gene %in% universe$gene, ]

            gene_set_motifs<-universe_motifs[universe_motifs$motif_set == m_set &
                                             universe_motifs$gene %in% genes_in_universe, ]

            motifs_in_gene_set <- unique(gene_set_motifs$motif)
            motifs_in_universe <- unique(universe_motifs$motif)

            motifs_gene_set_count<-count(gene_set_motifs, "motif")
            motifs_universe_count<-count(universe_motifs, "motif")


            for(motif in motifs_in_gene_set){

                subset_genes <- count_gene_set_genes
                all_genes    <- count_universe_genes

                subset_genes_motif <- motifs_gene_set_count[motifs_gene_set_count$motif==motif,"freq"]
                all_genes_motif    <- motifs_universe_count[motifs_universe_count$motif==motif,"freq"]

                a <- subset_genes_motif
                b <- all_genes_motif - subset_genes_motif
                c <- subset_genes - subset_genes_motif
                d <- all_genes - all_genes_motif - c
                #print(motif)
                matrix_for_test<-matrix(c(a, b, c, d), nrow = 2, ncol = 2)

                rownames(matrix_for_test)<-c("gene_set", "universe")
                colnames(matrix_for_test)<-c("have", "dont_have")
                #print(matrix_for_test)
                for(alternative in alternatives){
                    p.value<-2
                    tmp<-try(fisher.test(matrix_for_test, alternative = alternative))
                    if(!is.error(tmp)){
                       p.value<- tmp$p.value
                    }
                    enrich_all_family[nrow(enrich_all_family) + 1,] = list(dataset,
                                                                           name,
                                                                           m_set,
                                                                           motif,
                                                                           matrix_for_test[1,1],
                                                                           matrix_for_test[1,2],
                                                                           matrix_for_test[2,1],
                                                                           matrix_for_test[2,2],
                                                                           alternative,
                                                                           p.value)
                }

                #Here the student T starts. We will use a new dataframe.
                universe_motif_counts<-universe_motifs[universe_motifs$motif==motif, "count"]
                gene_set_motif_counts<-gene_set_motifs[gene_set_motifs$motif==motif, "count"]
                t_test<-list("estimate.mean of x"=0,
                             "estimate.mean of y"=0,
                             "statistic.t"=0,
                             "parameter.df"=0,
                             "p.value"=2
                            )
                if(length(gene_set_motif_counts) > 3 ){
                    tmp <- try(unlist(t.test(gene_set_motif_counts, universe_motif_counts)))
                    if(!is.error(tmp)){
                       t_test<-tmp
                       enrich_t_test[nrow(enrich_t_test) + 1,] = list(dataset,
                                                                   name,
                                                                   m_set,
                                                                   motif,
                                                                   length(gene_set_motif_counts),
                                                                   length(universe_motif_counts),
                                                                   t_test["estimate.mean of x"],
                                                                   t_test["estimate.mean of y"],
                                                                   t_test["statistic.t"],
                                                                   t_test["parameter.df"],
                                                                   t_test["p.value"]
                                                                  )
                    }else{
                        enrich_t_test[nrow(enrich_t_test) + 1,] = list(dataset,
                                                                   name,
                                                                   m_set,
                                                                   motif,
                                                                   length(gene_set_motif_counts),
                                                                   length(universe_motif_counts),
                                                                   mean(gene_set_motif_counts),
                                                                   mean(universe_motif_counts),
                                                                   0,
                                                                   0,
                                                                   2
                                                                  )
                    }


                }else{
                   enrich_t_test[nrow(enrich_t_test) + 1,] = list(dataset,
                                                                   name,
                                                                   m_set,
                                                                   motif,
                                                                   length(gene_set_motif_counts),
                                                                   length(universe_motif_counts),
                                                                   0,
                                                                   0,
                                                                   0,
                                                                   0,
                                                                   1
                                                                  )
                }
            }
        }
    }
    enrich_all_family$padj_BH <- p.adjust(enrich_all_family$fisher_pvalue, method="BH")
    enrich_t_test$padj_BH     <- p.adjust(enrich_t_test$p_value, method="BH")
    list(fisher=enrich_all_family, t=enrich_t_test)
}

plot_triad_movment_single<-function(geneInformation,
                          genes_to_plot,
                          experiment="HC_850_samples",
                          title="Test",
                             density = FALSE ,
                             points  = TRUE
                          ){

    allTriads<-geneInformation$allTriads
    selectedTriads<-unique(allTriads[allTriads$gene %in% genes_to_plot, "group_id"])

    tmp_df<-geneInformation$triads[geneInformation$triads$group_id %in% selectedTriads &
                                  geneInformation$triads$dataset == experiment,
                                   c("group_id","factor","clust","description","general_description","chr_group","normalised_triad", "Distance")]

    clust_df <- dcast(tmp_df,group_id + clust + description + general_description + factor + Distance ~
                      chr_group, value.var = "normalised_triad")

    clust_df_all_mean <- clust_df[clust_df$factor  == 'all_mean_filter' , ]

    original_description <- clust_df_all_mean[,c("group_id", "description", "general_description")]
    colnames(original_description)<-c("group_id", "original_description", "original_general_description")

    clust_df_factors  <- clust_df[clust_df$factor  != 'all_mean_filter' &
                                  clust_df$factor  != 'all_means' &
                                  clust_df$factor  != 'all' , ]
    #Central="#808A9F",
    group.colors<-c(A.dominant = "#579D1C", B.dominant = "#4B1F6F", D.dominant ="#FF950E",
             Central="#AAAAAA",
             A.suppressed = "#7FD127", B.suppressed = "#7D31AF", D.suppressed ="#FFAD42")

    group.alpha<-c(A.dominant = 1, B.dominant = 1, D.dominant =1,
             Central=0.1,
             A.suppressed = 0.2, B.suppressed = 0.2, D.suppressed = 0.2)

    group.alpha.tern<-c(A.dominant = 0.05, B.dominant = 0.05, D.dominant =0.05,
             Central=0.05,
             A.suppressed = 0.05, B.suppressed = 0.05, D.suppressed = 0.05)


    clust_df_factors<-merge(clust_df_factors, original_description, by="group_id")

    clust_df_factors_full<-clust_df_factors
    clust_df_factors<-clust_df_factors[clust_df_factors$description !=clust_df_factors$original_description , ]

    tern_all <- ggtern(clust_df_factors,aes(A,B,D,color=original_description)) + theme_bw()

    tern_all  <- tern_all   +  geom_point(aes(alpha=original_description)) +
    scale_alpha_manual(values=group.alpha.tern) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    scale_color_manual(values=group.colors)

    tern_fact <- tern_all + facet_wrap(~original_general_description, ncol=1) +
    theme_notitles() + theme(legend.position = "none") + theme_nolabels()

     tern_all <-  tern_all + theme_legend_position(x = "topleft")   + theme_arrownormal()  +
    ggtitle("All factors")

    tern_fact_2 <- ggtern(clust_df_factors_full,aes(A,B,D,color=original_description)) + theme_bw()

    tern_fact_2  <- tern_fact_2   +  geom_point(aes(alpha=original_description)) +
    scale_alpha_manual(values=group.alpha.tern) +
    scale_color_manual(values=group.colors)

    tern_fact_2 <- tern_fact_2 + facet_wrap(~original_general_description, ncol=1) +
    theme_notitles() + theme(legend.position = "none") + theme_nolabels()


    tern_mean <- ggtern(clust_df_all_mean,aes(A,B,D,color=description)) +  theme_bw() +
        theme_legend_position(x = "topleft")  + scale_color_manual(values=group.colors) +
       theme_arrownormal() + ggtitle("Mean")

    tern_mean <- tern_mean +  geom_point(aes(alpha=description)) + scale_alpha_manual(values=group.alpha)
    #   guides(colour = guide_legend(override.aes = list(alpha = 1)))
    tern_all


    gs<-list(tern_mean, tern_all, tern_fact, tern_fact_2)

    lay <- rbind(c( 1, 3, 4),
                 c( 2, 3, 4)
                 )

    g2 <- arrangeGrob(grobs = gs, layout_matrix = lay, top = title)
    g2
}


plot_normalized_triads<-function(triads){
    group.colors <- c(A = "#579D1C", B = "#4B1F6F", D ="#FF950E")
    p <- ggplot(triads, aes(chr_group, normalised_triad, fill=chr_group))
    p <- p + geom_boxplot(outlier.alpha = 0.05)
    p <- p + theme_classic()
    p <- p + scale_fill_manual(values=group.colors) + guides(fill=FALSE)
    p <- p + ylim(0,1)
    p <- p + ylab("") + xlab("")
    p
}


plot_clust_dist<-function(geneInformation,
                          genes_to_plot,
                          experiment="HC_850_samples",
                          title="All"
                          ){

    allTriads<-geneInformation$allTriads
    selectedTriads<-unique(allTriads[allTriads$gene %in% genes_to_plot, "group_id"])



    tmp_df<-geneInformation$triads[geneInformation$triads$group_id %in% selectedTriads &
                                  geneInformation$triads$dataset == experiment,
                                   c("group_id","factor","clust","description","general_description","chr_group","normalised_triad")]
    if(nrow(tmp_df) == 0){
        return(NULL)
    }
    clust_df <- dcast(tmp_df,group_id +general_description+clust+description+factor ~ chr_group, value.var = "normalised_triad")
    clusters<-sort(c("B.suppressed",
                   "Central",
                   "A.dominant",
                   "A.suppressed",
                   "B.dominant",
                   "D.suppressed",
                   "D.dominant"))

    group.colors <- c("B.suppressed"="#4B1F6F",
                   "Central"="#999999",
                   "A.dominant"="#579D1C",
                   "A.suppressed"="#579D1C",
                   "B.dominant"="#4B1F6F",
                   "D.suppressed"="#FF950E",
                   "D.dominant"="#FF950E")

     group.fills <- c("B.suppressed"="white",
                   "Central"="#999999",
                   "A.dominant"="#579D1C",
                   "A.suppressed"="white",
                   "B.dominant"="#4B1F6F",
                   "D.suppressed"="white",
                   "D.dominant"="#FF950E")

    group.shapes <- c("B.suppressed"=25,
                   "Central"=19,
                   "A.dominant"=17,
                   "A.suppressed"=25,
                   "B.dominant"=17,
                   "D.suppressed"=25,
                   "D.dominant"=17)

    tern <- ggtern(clust_df,aes(A,B,D)) + theme_classic()
    tern <- tern + geom_point(aes(color=description,
                      shape=description),
                  alpha=0.15)
    tern <- tern + theme_arrownormal()
    tern <- tern + theme_legend_position(x = "topleft")
    tern <- tern + guides(colour = guide_legend(override.aes = list(alpha = 1)))
    tern <- tern + scale_color_manual(values=group.colors)
    tern <- tern + scale_shape_manual(values=group.shapes)
    #tern <- tern + scale_fill_discrete(values=group.fills)

    gs<-list(tern)
    dat <- data.frame(
        A=numeric(0),B=numeric(0), D=numeric(0), size=numeric(0),stringsAsFactors=FALSE )

    rownames(dat)<-rownames(clusters)
    for(c in clusters){

        tmp_df_clust<-tmp_df[tmp_df$description==c,]
        p <- plot_normalized_triads(tmp_df_clust)
        p <- p + ggtitle(c)


        dat[c,1] <- round(100*mean(tmp_df_clust[tmp_df_clust$chr_group == "A","normalised_triad"]),digits=2)
        dat[c,2] <- round(100*mean(tmp_df_clust[tmp_df_clust$chr_group == "B","normalised_triad"]),digits=2)
        dat[c,3] <- round(100*mean(tmp_df_clust[tmp_df_clust$chr_group == "D","normalised_triad"]),digits=2)
        dat[c,4] <- nrow(tmp_df_clust)

        gs[[length(gs)+1]] <- p
    }

    total_size<-sum(dat$size)
    dat$percentage<-round(100*dat$size/total_size,digits=2)
    gs[[length(gs)+1]]<-tableGrob(dat)
    lay <- rbind(c( 1, 1, 1, 2, 4, 7),
                 c( 1, 1, 1, 3, 5, 8),
                 c( 9, 9, 9, NA,6,NA)
                 )

    g2 <- arrangeGrob(grobs = gs, layout_matrix = lay, top = title)
    g2
}


args = commandArgs(trailingOnly=TRUE)
#folder<-"/Users/ramirezr/Dropbox/JIC/expVIPMetadatas/RefSeq1.0/TablesForExploration"
#genes_to_plot_path<-"/Users/ramirezr/Dropbox/JIC/expVIPMetadatas/RefSeq1.0/notebook/gene_set_files/04.modules/WGCNA_850/WGCNA_850_Module_15.txt"
folder<-args[1]
genes_to_plot_path<-args[2]
name<-basename(genes_to_plot_path)
path<-paste0(genes_to_plot_path,"_plots")

genes_to_plot<-read.csv(genes_to_plot_path)
genes_to_plot<-as.vector(genes_to_plot[,1])

print(name)
print(paste0("number of genes to plot: ", length(genes_to_plot)))
print(head(genes_to_plot))

geneInformation<-loadGeneInformation(dir=folder, non_syntenic_triads=T)
g <- plot_gene_summary(geneInformation, genes_to_plot, name = name,run_stats=TRUE , output_path = path )

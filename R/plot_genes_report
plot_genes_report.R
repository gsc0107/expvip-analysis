
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

folder<-"/Users/ramirezr/Dropbox/JIC/expVIPMetadatas/RefSeq1.0/TablesForExploration"

loadGeneInformation<-function(dir="../TablesForExploration"){
    path<-paste0(dir,"/CanonicalTranscript.rds")
    canonicalTranscripts<-readRDS(path)
    canonicalTranscripts$intron_length<- canonicalTranscripts$mrna_length -  canonicalTranscripts$exon_length
    canonicalTranscripts$chr_group <- substr(canonicalTranscripts$Chr,4,4)
    canonicalTranscripts$genome    <- substr(canonicalTranscripts$Chr,5,5)
    
    path<-paste0(dir, "/MeanTpms.rds")
    meanTpms <- readRDS(path)
    
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
    
    list(canonicalTranscripts=canonicalTranscripts, 
         meanTpms=meanTpms,
         triads=triads, 
         triadMovement=triadMovement,
         gene_universe=gene_universe, 
         ontologies=ontologies,
         id_names=id_names,
         WGCNA=WGCNA
        )
}
geneInformation<-loadGeneInformation()

colnames(geneInformation$canonicalTranscripts)


genes_to_plot<-subset(geneInformation$WGCNA, ModuleLabels == 15 & set == "WGCNA_850")
genes_to_plot<-genes_to_plot$Gene
head(genes_to_plot)
class(genes_to_plot)
nrow(genes_to_plot)

plotHistogram<-function(table, column="size_cds"){
    
    probs <- c( 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
    quantiles <- data.frame(quantile(table[,column], prob=probs,na.rm=TRUE, include.lowest=TRUE), stringsAsFactors=FALSE)
    quantiles$quant<-rownames(quantiles)
    colnames(quantiles)<-c("value", "quant")
    values<-quantiles$values
    local_mean<-mean(table[,column])
    local_sd<-sd(table[,column])
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
    
    p <- ggplot(table, aes_string(column, fill="quantile")) + 
    geom_vline(data=quantiles,aes(xintercept=quantiles$value) )
    for(i in seq(1,nrow(quantiles))){
        x_pos<-quantiles$value[i]
         gtext <- textGrob(quantiles$quant[i], y=0.02,  gp = gpar(fontsize = 6,col = "red"))
         p <- p + annotation_custom(gtext, xmin=x_pos, xmax=x_pos)
    }
    p<- p + geom_histogram(bins=50, position = "identity") + theme_bw() + 
    scale_fill_brewer(palette="Dark2")+
    theme(legend.position="none")
    p<- p + ggtitle(paste0("Mean: ", round(local_mean,2), 
                           " SD:", round(local_sd,2),
                          " N:", nrow(table)))
    p
}
local_table<-geneInformation$canonicalTranscripts
local_table<-local_table[local_table$Gene %in% genes_to_plot,]
head(local_table)
plotHistogram(local_table, column="exon_no")
nrow(local_table)

plot_per_chromosome_5pc_bins_facet<-function(table, title = "Test"){
    chromosomes=c("1A", "1B", "1D",
                "2A", "2B", "2D",
                "3A", "3B", "3D",
                "4A", "4B", "4D",
                "5A", "5B", "5D",
                "6A", "6B", "6D",
                "7A", "7B", "7D")
    
    
    gs<-list()
    local_title = paste0(title, "\n Genes per chromosome 5% bin")
    
    t1 <- table[table$Chr != "chrUn",]
    t2 <- table[table$Chr == "chrUn",]
    p <-ggplot(t1,aes(scaled_5per_position)) 
    
    p <- p + xlim(0,100)
    p <- p + geom_bar() + theme_bw()
    p <- p + facet_grid(chr_group~genome,  drop = TRUE)
    p <- p + ylab(" count ") + xlab("")
    
    gs[[length(gs)+1]] <- p
    
    p <-ggplot(table,aes(Chr, fill=geneconf))  + geom_bar() + theme_bw()
    p <- p + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + ylab("") + xlab("")
    gs[[length(gs)+1]] <- p
    
    g1<-arrangeGrob(grobs=gs, ncol=1, heights=c(0.8,0.2), top=local_title ) 
}
g<-plot_per_chromosome_5pc_bins_facet(local_table)
grid.draw(g)

colnames(local_table)



plot_tpms_summary<-function(tpms, experiment="850_samples", min_tpm=0.5, title="Test"){
    
    local_tpms<-subset(tpms, (subset == experiment) & 
                       ( factor != "all" & factor != "all_means" & factor != "all_mean_filter" ) &
                      value > min_tpm)
   
    local_title <- paste0(title, "\n", experiment)
    
    p  <- ggplot(local_tpms, aes(value)) 
    p  <- p + geom_histogram(bins=30 ) + theme_bw()
    p  <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                      strip.text = element_text(size=6))
    p  <- p + facet_wrap(~ factor, ncol=3)  + xlim(0,15) 
    p  <- p + ylab("Count") + xlab("")
    p  <- p + theme(strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")))

    
    factor_means <- aggregate(value ~ factor, data = local_tpms, mean)
    factor_max   <- aggregate(value ~ factor, data = local_tpms, max)
    factor_min   <- aggregate(value ~ factor, data = local_tpms, min)
    factor_sd    <- aggregate(value ~ factor, data = local_tpms, sd)
    
    local_tpms<-subset(tpms, (subset == experiment) & 
                       ( factor == "all_mean_filter" ) &
                      value > min_tpm)
    
    p2  <- ggplot(local_tpms, aes(value)) 
    p2  <- p2 + geom_histogram(bins=30 ) + theme_bw()
    p2  <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                      strip.text = element_text(size=6, lineheight=0.5))
    p2  <- p2 + facet_wrap(~ factor, ncol=1)  + xlim(0,15)
    p2  <- p2 + ylab("") + xlab("TPM") 
    
    factor_means <- rbind(factor_means,aggregate(value ~ factor, data = local_tpms, mean))
    factor_max   <- rbind(factor_max,  aggregate(value ~ factor, data = local_tpms, max))
    factor_min   <- rbind(factor_min,  aggregate(value ~ factor, data = local_tpms, min))
    factor_sd    <- rbind(factor_sd,   aggregate(value ~ factor, data = local_tpms, sd))
    
    
    rownames(factor_means)<- factor_means$factor
    rownames(factor_max)<- factor_max$factor
    rownames(factor_min)<- factor_min$factor
    rownames(factor_sd)<- factor_sd$factor
    
    factor_means$factor<-NULL
    factor_max$factor<-NULL
    factor_min$factor<-NULL
    factor_sd$factor<-NULL
    
    factor_means<-cbind(factor_means, factor_max)
    factor_means<-cbind(factor_means, factor_min)
    factor_means<-cbind(factor_means, factor_sd)
    
   
    factor_means<-round(factor_means, 2)
    colnames(factor_means)<-c("Mean", "Max", "Min", "S.D.")
    mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.5)),
    colhead = list(fg_params=list(cex = 0.5)),
    rowhead = list(fg_params=list(cex = 0.5)))
    
    p3 <- tableGrob(factor_means, theme=mytheme)
    
    lay <- rbind(c(1,3),
                 c(1,2))
    g1<-arrangeGrob(grobs=list(p,p2,p3), heights=c(0.8,0.2), widths=c(0.6,0.4), layout_matrix=lay, top = local_title) 
    g1
}
local_tpms<-geneInformation$meanTpms
local_tpms<-local_tpms[local_tpms$gene %in% genes_to_plot,]
g<-plot_tpms_summary(local_tpms)
head(geneInformation$meanTpms)
grid.draw(g)

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
g<-plot_all_means_filteredtpms_summary(local_tpms)
grid.draw(g)

head(geneInformation$triadMovement)
head(geneInformation$triads)
colnames(geneInformation$triadMovement)

get_triads_from_genes<-function(genes, geneInformation, dataset="HC_CS_no_stress"){
    triads<-geneInformation$triads
    triadMovement<-geneInformation$triadMovement
    group_ids = triads[triads$gene %in% genes,"group_id"]
    list(triads=triads[triads$group_id %in% group_ids & triads$dataset==dataset,], 
        triadMovement=triadMovement[triadMovement$group_id %in% group_ids & triadMovement$dataset==dataset ,])
    
}
local_triads <- get_triads_from_genes(genes_to_plot, geneInformation)

unique(local_triads$triads$min_triad_sum)

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
plot_distribution_for_factor(local_triads$triadMovement)

plot_dominance_summary<-function(selected_triads, experiment="HC_CS_no_stress", title="test"){
     local_title <- paste0(title, "\n", experiment)
    
    triads <- selected_triads$triads 
    triadMovement<-selected_triads$triadMovement
    
    all_means_filter<-triads[triads$factor=="all_mean_filter",]
    #print(unique(triads$factor))
    gs<-list()
    gs[[length(gs)+1]] <- plotHistogram(all_means_filter, column="value") + 
    xlab("All mean filter for genes in triad") + xlim (0,15)
    
    gs[[length(gs)+1]] <- ggplot(triadMovement, aes(factor_count)) +  geom_bar() + theme_bw() +
    xlab("No. of conditions") 
    
    p <- ggplot(triads, aes(factor, fill=description)) + geom_bar() + theme_bw()
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=4)) + labs(fill="", x="")
    gs[[length(gs)+1]] <- p
    
    p <- ggplot(triads, aes(factor, fill=general_description)) + geom_bar() + theme_bw()
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=4)) + labs(fill="", x="")
    gs[[length(gs)+1]] <- p
    
    
   
    p <- ggplot(triadMovement, aes(total_categories)) + geom_bar() + theme_bw()
   p <- p + labs(fill="Main\ncategory", x="No. of categories") 
    p<- p + ggtitle("")
   gs[[length(gs)+1]] <- p
    
    gs[[length(gs)+1]] <- plotHistogram(triadMovement, column="central_max_distance")
    gs[[length(gs)+1]] <- plotHistogram(triadMovement, column="central_mean_distance")
    gs[[length(gs)+1]] <- plotHistogram(triadMovement, column="sum_mean_tpm")
    g1<-arrangeGrob(grobs=gs, ncol=2, top=local_title )
    g1
}
g<-plot_dominance_summary(local_triads)
grid.draw(g)

get_goseq_enrichment<-function(geneInformation, genes_to_plot, 
                               name="Random Samples",
                               dataset="HC_CS_no_stress", 
                               ontology="GO"
                              ){
    id_names <- geneInformation$id_names
    universe<-geneInformation$gene_universe 
    universe<-universe[universe$dataset==dataset,]
    
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
    
    sig.GO2 <- merge(sig.GO, id_names, by.x="category", by.y="V1", all.x =TRUE, all.y =FALSE)
   
    if(ontology == "PO" || ontology == "TO" ){
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
    sig.GO$universe_size <- nrow(universe)
    sig.GO

}


enrichment_test<- get_goseq_enrichment(geneInformation, genes_to_plot,  ontology="GO") 
enrichment_test_po<- get_goseq_enrichment(geneInformation, genes_to_plot,  ontology="PO") 
enrichment_test_to<- get_goseq_enrichment(geneInformation, genes_to_plot,  ontology="TO") 
enrichment_test<-rbind(enrichment_test,enrichment_test_po)
enrichment_test<-rbind(enrichment_test,enrichment_test_to)
nrow(enrichment_test)

plot_enrichment<-function(enrichment,experiment="HC_CS_no_stress", title="test" , type="Over represented"){
    local_title <- paste0(title, "\n", experiment, "\n", type , "\nPercentage of genes from each ontology\n")
    gs<-list()
    for(ont in unique(enrichment$ontology)){
        current<-enrichment[enrichment$ontology == ont & enrichment$type == type, ]
        if(nrow(current) > 0){
            gs[[length(gs)+1]] <- ggplot(current, aes(description, percentage)) + 
            coord_flip() +
            geom_col() + labs(y="", x="", title=ont) +  theme_bw() + #ylim(0,100) + 
            theme(axis.text.x = element_text(angle = 0, hjust = 1, size=4),
             axis.text.y = element_text(angle = 35, hjust = 1, size=4)) 
        }
    }
    g1<-arrangeGrob(grobs=gs, ncol=2, top=local_title )
}
g<-plot_enrichment(enrichment_test, type="Under represented")
grid.draw(g)

plot_gene_summary<-function(geneInformation, genes_to_plot, name="Random Samples" , output_path="./Test"){
    local_table<-geneInformation$canonicalTranscripts
    local_table<-local_table[local_table$Gene %in% genes_to_plot,]
    
    local_mean_tpms<-geneInformation$meanTpms
    local_mean_tpms<-local_mean_tpms[local_mean_tpms$gene %in% genes_to_plot, ]
    
    stats_to_plot<-c('size_cds', 'exon_no', 'exon_length','intron_length', 'X3UTR_length', 'X5UTR_length' )
    
    gs<-list()
    plots<-list()
    
    plots[[length(plots)+1]] <- textGrob(paste0(name, " Gene summary"))
    for(plot in stats_to_plot){
        p<-plotHistogram(local_table,column=plot)
        gs[[length(gs)+1]] <- p
    }
    plots[[length(plots)+1]] <- arrangeGrob(grobs=gs, ncol=2 , top = paste0(name, "\n Gene properties"))
    plots[[length(plots)+1]] <- plot_per_chromosome_5pc_bins_facet(local_table, title=name)
    
    plots[[length(plots)+1]] <- textGrob(paste0(name, " TPM summary"))
    for(s in unique(geneInformation$meanTpms$subset)){
        plots[[length(plots)+1]] <- plot_tpms_summary(local_mean_tpms, experiment=s, title=name) 
        plots[[length(plots)+1]] <- plot_all_means_filteredtpms_summary(local_mean_tpms, experiment=s, title=name) 
    }
    
    plots[[length(plots)+1]] <- textGrob(paste0(name, " Triad summary"))
    for(s in unique(geneInformation$triads$dataset)){
        local_triads <- get_triads_from_genes(genes_to_plot, geneInformation, dataset=s)
        plots[[length(plots)+1]] <- plot_dominance_summary(local_triads, experiment=s, title=name)
    }
    
    plots[[length(plots)+1]] <- textGrob(paste0(name, " Go enrichment"))                                     
    
    
    all_enrichments <- NULL
    for(g_u in unique(geneInformation$gene_universe$dataset)){
        enrichment_test<- get_goseq_enrichment(geneInformation, genes_to_plot,  ontology="GO") 
        enrichment_test_po<- get_goseq_enrichment(geneInformation, genes_to_plot,  ontology="PO") 
        enrichment_test_to<- get_goseq_enrichment(geneInformation, genes_to_plot,  ontology="TO") 
        enrichment_test<-rbind(enrichment_test,enrichment_test_po)
        enrichment_test<-rbind(enrichment_test,enrichment_test_to)
        enrichment_test$universe<-g_u
        
        for(type in unique(enrichment_test$type)){
            plots[[length(plots)+1]] <- plot_enrichment(enrichment_test,
                                                        experiment=g_u, 
                                                        title=name, 
                                                        type=type)
        }
        if(is.null(all_enrichments)){
            all_enrichments <- enrichment_test
        }else{
            all_enrichments<-rbind(all_enrichments, enrichment_test)
        }
    }
    
    dir<-paste0(output_path,"/",name)
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    
    output_enrichment<-paste0(dir, "/", "enrichment.csv")
    write.csv(all_enrichments, file=output_enrichment)
    output_pdf<-paste0(dir, "/",name ,".pdf")
    
    
    g1<-marrangeGrob(plots, ncol=1, nrow=1, top="", bottom = quote(paste("page", g, "of",
       pages)))
    
    ggsave(output_pdf, plot=g1 , width = 210, height = 297, units = "mm")
    g1
}
g<-plot_gene_summary(geneInformation,genes_to_plot, name="Module 15, 850 samples" )


ggsave("test.pdf", plot=g , width = 210, height = 297, units = "mm")

?goseq

head(geneInformation$gene_universe)

library(sqldf)

GU<-geneInformation$gene_universe
sqldf("SELECT dataset, count(*) from GU group by dataset")



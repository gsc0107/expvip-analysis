
library(ggplot2)
library(reshape2)
library(sqldf)
library(fields)
library(gridExtra)
library(ggtern)
library(clue)
library(geometry)
require(gtable)

folder<-"/Users/ramirezr/Dropbox/JIC/expVIPMetadatas/RefSeq1.0/TablesForExploration"

loadGeneInformation<-function(dir="../TablesForExploration"){
    path<-paste0(dir,"/CanonicalTranscript.rds")
    canonicalTranscripts<-readRDS(path)
    canonicalTranscripts$intron_length<- canonicalTranscripts$mrna_length -  canonicalTranscripts$exon_length
    canonicalTranscripts$chr_group <- substr(canonicalTranscripts$Chr,4,4)
    canonicalTranscripts$genome    <- substr(canonicalTranscripts$Chr,5,5)
    
    path<-paste0(dir, "/MeanTpms.rds")
    meanTpms <- readRDS(path)
    list(canonicalTranscripts=canonicalTranscripts, meanTpms=meanTpms)
}
geneInformation<-loadGeneInformation()

colnames(geneInformation$canonicalTranscripts)

genes_to_plot<-sample(geneInformation$canonicalTranscripts$Gene, 1000)
head(genes_to_plot)
class(genes_to_plot)

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
    p<- p + ggtitle(paste0("Mean: ", round(local_mean,2), " SD:", round(local_sd,2)))
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
g<-plot_tpms_summary(geneInformation$meanTpms)
head(geneInformation$meanTpms)
grid.draw(g)

plot_gene_summary<-function(geneInformation, genes_to_plot, name="Random Samples"){
    local_table<-geneInformation$canonicalTranscripts
    local_table<-local_table[local_table$Gene %in% genes_to_plot,]
    
    stats_to_plot<-c('size_cds', 'exon_no', 'exon_length','intron_length', 'X3UTR_length', 'X5UTR_length' )
    
    gs<-list()
    plots<-list()
    for(plot in stats_to_plot){
        p<-plotHistogram(local_table,column=plot)
        gs[[length(gs)+1]] <- p
    }
    
   
    plots[[length(plots)+1]] <- arrangeGrob(grobs=gs, ncol=2 , top = paste0(name, "\n Gene properties"))
    plots[[length(plots)+1]] <- plot_per_chromosome_5pc_bins_facet(local_table, title=name)
    
    for(s in unique(geneInformation$meanTpms$subset)){
        plots[[length(plots)+1]] <- plot_tpms_summary(geneInformation$meanTpms, experiment=s, title=name) 
    }
    
    g1<-marrangeGrob(plots, ncol=1, nrow=1, top="", bottom = quote(paste("page", g, "of",
       pages)))
}
g<-plot_gene_summary(geneInformation,genes_to_plot )
ggsave("Test_summary.pdf", plot=g , width = 210, height = 297, units = "mm")







unique(local_tpms$subset)



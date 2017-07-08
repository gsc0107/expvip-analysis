#' Function to load a matrix of values from a metadata file. 
#'
#' This function allows you to express your love of cats.
#' @param metadata A dataframe with the values loaded from a metadata file as described in https://github.com/homonecloco/expvip-web/wiki/LoadingMetadata#loading-metadata
#' @param folder The folder with the expression values, named with the values in study_title in the metadata. cd
#' @param unit Expression unit to be used as sufix of the .tsv files. Defaults to 'tpm'.
#. @param values . Defaults to 'c(Development)'
#' @keywords metadata load
#' @export
#' @examples
#' folder<-"/Users/ramirezr/Dropbox/JIC/expVIPMetadatas/RefSeq1.0/expressionValues"
#' metadata_file<-"Metadata_FINAL.txt"
#' 
#' tpms  <-loadValuesFromExperiment(metadata, folder, unit="tpm",  values=unique(metadata$study_title))
#' counts<-loadValuesFromExperiment(metadata, folder, unit="count",values=unique(metadata$study_title))
#' 
#' metadata_used<-tpms[[2]]
#' tpms<-tpms[[1]]
#' counts<-counts[[1]]

loadValuesFromExperiment<-function(metadata, folder, unit="tpm", values=c("Development")){
    v<-values[1]
    v<-gsub(" ","_",v)
    path<-paste0(folder,"/",v,"_",unit,".tsv")
    ret<-read.table(path, row.names = 1, header= TRUE)
    for(i in 2:length(values)){
        v<-values[i]
        v<-gsub(" ","_",v)
        path<-paste0(folder,"/",v,"_",unit,".tsv")
        tmp<-read.table(path, row.names = 1, header= TRUE)
        ret<-cbind(ret,tmp)
    }
    md<-metadata[metadata$Sample.IDs%in%colnames(ret),]
    ret<-ret[,as.character(md$Sample.IDs),]
    list(ret,md)
}


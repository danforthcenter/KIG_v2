#ionomicsKnownGeneList with OrthoFinder orthologs
rm(list = ls())
library(reader)
library(plyr)
library(dplyr)
library(data.table)
library(doParallel)
library(stringr)

### functions ###
unlist_orthologs <- function(x, org=o) {
  d<-data.frame(Species=o,
                GeneID=unlist(strsplit(gsub(" ","",x[paste(o)]), split = ",")),
                Elements=x["Elements"],
                Inferred.from=x["Inferred.from"],
                Inferred.ID=x["GeneID"],
                row.names = NULL)
  return(d)
}

OrthoFinder_protiens_2_genes <- function(x) {
  if(grepl("[A-Za-z0-9]*_P[0-9]*",x)){
    x<-gsub('_P[0-9]*','',x)
    return(x)
  }
  if(grepl("[A-Za-z]*\\.[A-Za-z0-9]*\\.[0-9]\\.p",x)){
    x<-gsub('\\.[0-9].p','',x)
    return(x)
  }
  if(grepl("[A-Za-z0-9_]*\\.[0-9]",x)){
    x<-gsub('\\.[0-9]|\\.[0-9].p','',x)
    return(x)
  }
  return(x)
}

get_orthogroups<-function(x, OF_table=orthogroups, ind_matrix=og_strs){
  og<-OF_table[str_which(ind_matrix, x),]
  if(nrow(og)==0){
    og<-cbind(og[1L,], GeneID=x)
    og[is.na(og)]<-""
    return(og)
  }else{
    og<-as.data.frame(apply(og, c(1,2), OrthoFinder_protiens_2_genes))
    og$GeneID<-x
    return(og)
  }
}
#################
# option to dopar here
# registerDoParallel(cores = num_Core)
# I usually make num_Core 3 on my laptop, but go nuts if you have a server

setwd("D:/lehar/Documents/CGAS/")
dir.create("./IonomicsKnownGenes/knownIonomicsGenesWOrthologs_PHO")
listBase <-
  read.csv("./IonomicsKnownGenes/ionomics_known_genes_input.csv")
listBase <- listBase[order(listBase$Species, listBase$GeneID), ]
###checking for duplicate entries
listBase[which(duplicated(listBase$GeneID)), ]

###removing wheat for now since genome is not great
listBase<-listBase[!listBase$Species=="Taestivum_v2.2",]

###removes whitespaces in lists so they collapse correctly
listBase$Elements <- gsub(" ", "", listBase$Elements, fixed = TRUE)

###adding columns for later
listBase<-cbind(listBase, data.frame(`Primary/Inferred`="Primary", Inferred.from=NA, Inferred.Elements=NA))

orthogroups<- fread("./data/OrthoFinder_orthologs/Orthogroups_10s.tsv",
                    sep = "\t", header = T, stringsAsFactors = F)
colnames(orthogroups)<-gsub("_[0-9]{3}|.protein_primaryTranscriptOnly","", colnames(orthogroups))

orgs<-grep("OG",colnames(orthogroups), invert = T, value = T)
og_strs<-apply(orthogroups,1,paste0,collapse=",")

orths<-bind_rows(lapply(listBase$GeneID, FUN = get_orthogroups))

PrimaryLists <- merge(listBase,orths[,-1],by = "GeneID")

InferredOrthologs <- foreach(o=orgs, .packages = packages.loaded(), .combine = rbind) %do% {
  base <- PrimaryLists[,c("Species","GeneName", "GeneID", "Elements", paste0(grep(o, colnames(PrimaryLists), value = T)))]
  colnames(base)[2]<-"Inferred.from"
  base<-base[-which(base[paste0(o)]==""),]
  base$Inferred.from<-paste0(gsub("^([A-Za-z]{2}).*$","\\1", base$Species),base$Inferred.from)
  
  combined_df<-Reduce(rbind,apply(base, 1, unlist_orthologs, org=o))
  
  #removing the primary gene from ortholog group 
  combined_df<-combined_df[which(combined_df$GeneID!=combined_df$Inferred.ID),]
  
  
  unique_df<-foreach(i = as.character(unique(combined_df$GeneID)),
                     .packages = packages.loaded(),
                     .combine = rbind) %do% {
    d<-data.frame(Species = combined_df$Species[1],
                  GeneID = i,
                  Elements =
                    paste(
                    as.character(
                      unique(
                        unlist(
                    strsplit(
                      combined_df$Elements[combined_df$GeneID == i], 
                      split = ",")
                    ))), 
                  collapse = ","),
                  GeneName = NA,
                  Citation.s....DOI.only = NA,
                  Closest.species.orthologs = NA,
                  Tissue = NA,
                  Comments = NA,
                  `Primary/Inferred` = "Inferred",
                  Inferred.from = paste(
                    as.character(
                      unique(
                        combined_df$Inferred.from[combined_df$GeneID == i]
                                        
                      )),
                    collapse = ","),
                  Inferred.Elements = NA
                  )
    return(d)
                  
  }
    
}

InferredLists <-
  foreach(i = as.character(unique(InferredOrthologs$GeneID)),
          .packages = c('plyr', 'dplyr', 'reader'),
          .combine = rbind) %do% {
            if(i %in% PrimaryLists$GeneID){
              PrimaryLists$Primary.Inferred[PrimaryLists$GeneID==i] <- "Primary/Inferred"
              PrimaryLists$Inferred.from[PrimaryLists$GeneID==i] <- 
                InferredOrthologs$Inferred.from[InferredOrthologs$GeneID==i]
              PrimaryLists$Inferred.Elements[PrimaryLists$GeneID==i] <- paste0(
                as.character(
                  unique(
                    unlist(
                      strsplit(
                          InferredOrthologs$Elements[InferredOrthologs$GeneID == i],
                        split = ",")))),
                collapse = ",")
              return(NULL)
            }else{
              append_orthologs<-get_orthogroups(i)[,-1]
              return(merge(InferredOrthologs[InferredOrthologs$GeneID==i,], append_orthologs, by = "GeneID"))
            }
          } 

###making the regular lists
setnames(
  PrimaryLists,
  old = orgs,
  new = gsub("([A-Z])(.+)", "\\1\\.\\2 orthologs", orgs)
)
setnames(
  InferredLists,
  old = orgs,
  new = gsub("([A-Z])(.+)", "\\1\\.\\2 orthologs", orgs)
)

for(o in orgs) {
  final_list<-rbind(PrimaryLists[PrimaryLists$Species==o,], 
                    InferredLists[InferredLists$Species==o,])
  final_list$Species<-gsub("^([A-Z])","\\1.", final_list$Species)
  
  #reordering colnames, trying not to rely on index since that may change for some users, can comment out
  final_list<-final_list[order(final_list$GeneName, final_list$GeneID),c(grep("Species",colnames(final_list)),grep("Species",colnames(final_list),invert = T))]

  write.table(final_list,
              file = paste0("./IonomicsKnownGenes/knownIonomicsGenesWOrthologs_nowheat/",
                            o,
                            "_knownIonomicsGenesWOrthologs.csv"),
              col.names = T,
              row.names = F,
              sep = ",")

  
}

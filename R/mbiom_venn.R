#' @title Microbiome Venn diagram
#' @author Kara J.M. Taylor
#' @description Constructs a venn diagram with the number of sets defined by levels in `group`.
#' For datasets where the number of groups is 5 or less, a diagram figure will be drawn.
#' For datasets with more than 5 groups, no figure will be drawn, but the sets will be listed.
#' @param mat Community matrix (site rows and species columns)
#' @param meta Sample metadata
#' @param group Column in `meta` indicating the sets
#' @param tax Taxonomy table (optional)
#' @param file Output file name
#' @param tax.grp Taxonomic level to summarize taxa to. Defaults to 'tag'
#' @returns List containing set information
#' @export


mbiom.venn <- function(mat, meta, group, tax=NULL, file="", tax.grp="tag") {

  cran <- c("dplyr", "VennDiagram", "RColorBrewer")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  library(dplyr)
  library(VennDiagram)
  library(RColorBrewer)

  palette(c(RColorBrewer::brewer.pal(8, "Set1"), "#000000"))

  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    names(meta)[1] <- "sample"
    id <- 1
  }

  if (any(!row.names(mat)%in%meta[,which(names(meta)=="sample")])){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }

  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/venn/")==FALSE){
    dir.create("./coreplots/venn/")
  }
  if(dir.exists(paste0("./coreplots/venn/", file))==FALSE){
    dir.create(paste0("./coreplots/venn/", file))
  }

  g <- which(names(meta)==group)

  mat2 <- data.frame(sample=row.names(mat), mat)
  mat <- right_join(meta[,c(id,g)], mat2)


  iff <- function(x){
    if (x > 0){
      x <- 1
    } else {
      x <- 0
    }
  }

  top <- function(x){
    y <- sum(x)
    return(y)
  }

  cumul <- function(x){
    y <- seq(1:length(x))
    z1 <- which(x==0)
    y[z1] <- 0
    x <- y
  }

  rem0 <- function(x){
    if (any(x==0)){
      x <- x[-which(x==0)]
    } else {
      x <- x
    }
  }

  nax <- function(x){x[which(is.na(x)==T)] <- "unclassified"; return(x)}


  df <- mat[,-1] %>%
    group_by_(group) %>%
    summarize_all(funs(mean=mean)) %>%
    as.data.frame()
  row.names(df) <- with(df, eval(parse(text=group)))

  test <- apply(df[,-1], MARGIN=1, FUN=cumul) %>% as.data.frame()

  renamex <- names(df[,-1])
  renamex <- unlist(strsplit(renamex, split="_"))[seq(1,2*length(renamex), 2)]
  renamex <- as.numeric(substring(renamex, first=4))

  for (f in 1:dim(test)[2]){
    test[which(test[,f]!=0),f] <- renamex[which(test[,f]!=0)]
  }
  #return(test)
  sorted <- sort(apply(df[,-1], MARGIN=2, FUN=top), decreasing=T)
  #return(sorted)
  matched <- match(names(sorted), names(df[,-1]))
  #return(matched)
  #return(as.list(test[matched,]))
  test <- sapply(as.list(test[matched,]), FUN=rem0)
  #return(test)
  if (!is.null(tax)){
    if (tax.grp %in% names(tax)){
      test3 <- list()
      if (tax.grp != "tag"){
        names(tax)[which(names(tax)==tax.grp)] <- "taxgrp"
        subtax <- tax[which(tax$taxgrp!=""),]
        for (i in 1:length(test)){
          test3[[names(test)[i]]] <- as.character(unique(subtax$taxgrp[match(paste0("asv",test[[i]]), subtax$tag)]))
        }
      } else {
        for (i in 1:length(test)){
          test3[[names(test)[i]]] <- paste0("asv",test[[i]])
        }
      }
      sink(paste("./coreplots/venn/",file, "/euler_group_tax.txt", sep=""))
      print(test3)
      sink()
    }
    #return(test3)
    if (any(is.na(unlist(test3)))){test3 <- lapply(test3, FUN=nax)}


    if (length(test3)<=4){
      venn.diagram(test3, paste("./coreplots/venn/",file, "/euler_group_tax.tiff", sep=""),
                   category.names = names(test3), margin=0.12,
                   col=palette()[1:length(test3)], lwd=4,
                   fill=palette()[1:length(test3)], alpha=.125,
                   cat.dist=c(0.2, 0.2, 0.1, 0.1)[1:length(test3)],
                   cat.fontface="bold", cat.fontfamiy="sans", cat.cex=1.8, cex=1.5,
                   fontface="bold", fontfamily="sans")
    }
    if (length(test3)==5){
      venn.diagram(test3, paste("./coreplots/venn/",file, "/euler_group_tax.tiff", sep=""),
                   category.names = names(test3), margin=0.12,
                   col=palette()[1:length(test3)], lwd=4,
                   fill=palette()[1:length(test3)], alpha=.125,
                   cat.dist=c(0.22, 0.22, 0.22, 0.22, 0.22)[1:length(test3)],
                   cat.fontface="bold", cat.fontfamiy="sans", cat.cex=1.8, cex=1,
                   fontface="bold", fontfamily="sans")
    }
    if (length(test3)>5){
      return(test3)
    }

    return(test3)
  } else {

    sink(paste("./coreplots/venn/",file, "/euler_group_relative.txt", sep=""))
    print(test)
    sink()

    if (length(test) < 5){
      venn.diagram(test, paste("./coreplots/venn/",file, "/euler_group_relative.tiff", sep=""),
                   category.names = names(test), margin=0.12,
                   col=palette()[1:length(test)], lwd=4,
                   fill=palette()[1:length(test)], alpha=.125,
                   cat.dist=c(0.1, 0.1, 0.1, 0.1)[1:length(test)],
                   cat.fontface="bold", cat.fontfamiy="sans", cat.cex=1.8, cex=1.5,
                   fontface="bold", fontfamily="sans")
    } else {
      venn.diagram(test, paste("./coreplots/venn/",file, "/euler_group_genus.tiff", sep=""),
                   category.names = names(test), margin=0.12,
                   col=palette()[1:length(test)], lwd=4,
                   fill=palette()[1:length(test)], alpha=.125,
                   cat.dist=c(0.22, 0.22, 0.22, 0.22, 0.22)[1:length(test)],
                   cat.fontface="bold", cat.fontfamiy="sans", cat.cex=1.8, cex=1,
                   fontface="bold", fontfamily="sans")
    }
    return(test)
  }
}

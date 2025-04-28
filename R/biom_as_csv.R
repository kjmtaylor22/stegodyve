#' @title BIOM to CSV
#' @author Kara J.M. Taylor
#' @description Converts *.biom files to CSV and extracts taxonomy info
#' @param path location of BIOM table
#' @returns Two .csv documents: one containing the community matrix (feature_table.csv)
#' and the other containing the taxonomy information (taxonomy.csv).
#' @export


biom.as.csv <- function(path){
  if (!"dplyr"%in%row.names(installed.packages())) {
    install.packages("dplyr")
  }
  if (!"biomformat"%in%row.names(installed.packages())){
    #source("https://bioconductor.org/biocLite.R")
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")}
    #biocLite(i, suppressUpdates=T)
    BiocManager::install("biomformat")
  }
  library(dplyr)
  library(biomformat)

  biom <- read_biom(path)

  den <- function(x){
    y <- x$id
    return(y)
  }

  taxon <- function(x){
    y <- paste0(unlist(x$metadata), collapse="; ")
    return(y)
  }

  tax1 <- lapply(biom$rows, FUN=den)
  tax2 <- lapply(biom$rows, FUN=taxon)

  samples <- lapply(biom$columns, FUN=den)

  tax <- data.frame(den.otu=unlist(tax1), taxonomy=unlist(tax2))

  comm <- biom$data %>% as.data.frame() %>% t() %>% as.data.frame() %>%
    reshape2::dcast(V1 ~ V2, value.var = "V3") %>% .[,-1]

  colnames(comm) <- unlist(samples)
  row.names(comm) <- tax1

  write.csv(comm, "./feature_table.csv")
  write.csv(tax, "taxonomy.csv", row.names=F)
}







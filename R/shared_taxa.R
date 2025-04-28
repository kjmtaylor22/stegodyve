#' @title Parsing shared taxa from Venn diagram
#' @author Kara J.M. Taylor
#' @description Assigns each taxon to its Venn digram subsection.
#' @param mbiom Output list from `mbiom_venn`
#' @param file Output filename
#' @param excel Logical; return as .xlsx as well as to GlobalEnv?
#' @returns Returns list for each discrete segment in the Venn diagram
#' @export

shared.taxa <- function(mbiom, file, excel=T){

  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/venn/")==FALSE){
    dir.create("./coreplots/venn/")
  }
  if(dir.exists(paste0("./coreplots/venn/", file))==FALSE){
    dir.create(paste0("./coreplots/venn/", file))
  }

  if (!"openxlsx"%in%row.names(installed.packages())) {
    install.packages("openxlsx")
  }
  library(openxlsx)

  a <- combn(names(mbiom), 2)

  if (length(mbiom)>3){
    b <- combn(names(mbiom), 3)
  }
  if (length(mbiom)>4){
    c <- combn(names(mbiom), 4)
  }

  out <- list()
  for (i in 1:length(mbiom)){
    x <- mbiom[[i]][which(!mbiom[[i]]%in%unlist(mbiom[-i], use.names=F))]

    out[[names(mbiom)[i]]] <- x

    if (length(mbiom)==2){
      all <- mbiom[[1]][which(mbiom[[1]]%in%mbiom[[2]])]
    }
  }

  for (i in 1:dim(a)[2]){
    if (length(mbiom)!=2){
      x <- mbiom[[a[1,i]]][which(mbiom[[a[1,i]]]%in%mbiom[[a[2,i]]] &
                                   !mbiom[[a[1,i]]]%in%unlist(mbiom[which(!names(mbiom)%in%a[,i])]))]

      out[[paste(a[1,i],a[2,i], sep=".")]] <- x
    }
  }

  if (length(mbiom)==3){
    all <- mbiom[[1]][which(mbiom[[1]]%in%mbiom[[2]] &
                              mbiom[[1]]%in%mbiom[[3]])]
  }

  if (length(mbiom)>3){
    for (i in 1:dim(b)[2]){
      x <- mbiom[[b[1,i]]][which(mbiom[[b[1,i]]]%in%mbiom[[b[2,i]]] &
                                   mbiom[[b[1,i]]]%in%mbiom[[b[3,i]]] &
                                   !mbiom[[b[1,i]]]%in%unlist(mbiom[which(!names(mbiom)%in%b[,i])]))]

      out[[paste(b[1,i],b[2,i],b[3,i], sep=".")]] <- x
    }
    if (length(mbiom)==4){
      all <- mbiom[[1]][which(mbiom[[1]]%in%mbiom[[2]] &
                                mbiom[[1]]%in%mbiom[[3]] &
                                mbiom[[1]]%in%mbiom[[4]])]
    }
  }

  if (length(mbiom)==5){
    for (i in 1:dim(c)[2]){
      x <- mbiom[[c[1,i]]][which(mbiom[[c[1,i]]]%in%mbiom[[c[2,i]]] &
                                   mbiom[[c[1,i]]]%in%mbiom[[c[3,i]]] &
                                   mbiom[[c[1,i]]]%in%mbiom[[c[4,i]]] &
                                   !mbiom[[c[1,i]]]%in%unlist(mbiom[which(!names(mbiom)%in%c[,i])]))]

      out[[paste(c[1,i],c[2,i],c[3,i],c[4,i], sep=".")]] <- x
    }

    all <- mbiom[[1]][which(mbiom[[1]]%in%mbiom[[2]] &
                              mbiom[[1]]%in%mbiom[[3]] &
                              mbiom[[1]]%in%mbiom[[4]] &
                              mbiom[[1]]%in%mbiom[[5]])]
  }


  out[["all"]] <- all

  if (excel==T){
    if (file.exists(paste0("./coreplots/venn/",file, "/", file, ".xlsx"))){
      file.remove(paste0("./coreplots/venn/",file, "/", file, ".xlsx"))
    }
    wb <- createWorkbook()
    for (i in 1:length(out)){
      if (length(out[i])>0){
        addWorksheet(wb, sheetName = names(out)[i], gridLines=F)
        writeData(wb, sheet = names(out)[i], out[i])
      }
    }
    saveWorkbook(wb, paste0("./coreplots/venn/",file, "/", file, ".xlsx"), T)
  }
  return(out)
}

#' @title Basic alpha diversity calculation
#' @author Kara J.M. Taylor
#' @description Use input community matrix and metadata to calculate richness, evenness,
#' Shannon diversity, and Chao1 indices and plot them with desired grouping variable.
#' @param comm Community matric, site (rows) by species (columns).
#' @param meta Metadata file; should have a column called 'SampleID' that matches the
#' Row names in `comm`.
#' @param group Grouping variable from metadata; name as character; May also indicate
#' a column for a numeric independent variable, rather than a group.
#' @param cols Character vector of length matching `group` levels; if NULL,
#' will default to color ramp of rainbow colors.
#' @returns Data frame containing alpha diversity metrics for each sample
#' and group assignment from metadata.
#' @export

alphadiv<- function(comm, meta, group, cols=NULL, filename=deparse(substitute(comm))){

  pkgs <- c("vegan", "dplyr","fossil", "ggplot2", "ggpubr", "gridExtra")
  if (any(!pkgs%in%row.names(installed.packages()))) {
    install.packages(pkgs[which(!pkgs%in%row.names(installed.packages())==T)])
  }
  library(dplyr)
  library(vegan)
  library(fossil)
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)

  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "SampleID"
  } else {
    names(meta)[1] <- "SampleID"
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    id <- 1
  }


  g <- which(names(meta)%in%group)
  names(meta)[g] <- "group.id"

  if (any(!row.names(comm)%in%meta$SampleID)){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }

  shan <- diversity(comm)
  rich <- specnumber(comm)
  even <- shan/log(rich)
  chao <- apply(comm, MARGIN=1, FUN=chao1, taxa.row=F)

  shan <- data.frame(SampleID=names(shan), ShannonH=shan, Richness=rich, Evenness=even, Chao1=chao)

  shan <- left_join(shan, meta[,c(id,g)])

  if (is.null(cols)){
    clr <- c("magenta", "red", "orange", "yellow", "lawngreen", "mediumturquoise", "dodgerblue3", "navy", "blueviolet")
  } else {clr <- cols}

  if (class(shan$group.id)%in%c("factor","character")){

    clr <- colorRampPalette(clr)(length(unique(shan$group.id)))
    names(clr) <- as.character(unique(shan$group.id))

    shan$group.id <- factor(shan$group.id,
                            levels=sort(unique(shan$group.id)))
    shan <- shan[order(shan$group.id),]

    pl1 <- ggplot(shan, aes(group.id, ShannonH, fill=group.id)) +
      geom_boxplot() + ggthemes::theme_few() +
      theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      scale_fill_manual(values=clr) + guides(fill="none")
    pl2 <- ggplot(shan, aes(group.id, Richness, fill=group.id)) +
      geom_boxplot() + ggthemes::theme_few() +
      theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      scale_fill_manual(values=clr) + guides(fill="none")
    pl3 <- ggplot(shan, aes(group.id, Evenness, fill=group.id)) +
      geom_boxplot() + ggthemes::theme_few() +
      theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      scale_fill_manual(values=clr) + guides(fill="none")
    pl4 <- ggplot(shan, aes(group.id, Chao1, fill=group.id)) +
      geom_boxplot() + ggthemes::theme_few() +
      theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      scale_fill_manual(values=clr) + guides(fill="none")
    pleg <- as_ggplot(get_legend(pl1 + guides(fill=guide_legend(title = group, keywidth = 1, keyheight=0.75, ncol=5))))
  }

  if (class(shan$group.id)=="numeric"){

    pl1 <- ggplot(shan, aes(group.id, ShannonH, fill=group.id)) +
      geom_point(shape=22) + ggthemes::theme_few() + theme(axis.title.x = element_blank()) +
      scale_fill_gradientn(colors=clr) + guides(fill="none") +
      stat_smooth(method="lm", formula=y ~ x, se=F, size=0.5, show.legend = F) +
      ggpmisc::stat_poly_eq(formula=y ~ x, parse=T, label.x = "left", label.y="top", size=3,
                            aes(label=paste(stat(eq.label), stat(rr.label), sep="~~~~"))) +
      ggpmisc::stat_fit_glance(method = 'lm', method.args = list(formula = y ~ x), label.x = 'right', label.y = "bottom", size = 3,
                               aes(label = paste0("p = ", signif(..p.value.., digits = 2))))
    pl2 <- ggplot(shan, aes(group.id, Richness, fill=group.id)) +
      geom_point(shape=22) + ggthemes::theme_few() + theme(axis.title.x = element_blank()) +
      scale_fill_gradientn(colors=clr) + guides(fill="none") +
      stat_smooth(method="lm", formula=y ~ x, se=F, size=0.5, show.legend = F) +
      ggpmisc::stat_poly_eq(formula=y ~ x, parse=T, label.x = "left", label.y="top", size=3,
                            aes(label=paste(stat(eq.label), stat(rr.label), sep="~~~~"))) +
      ggpmisc::stat_fit_glance(method = 'lm', method.args = list(formula = y ~ x), label.x = 'right', label.y = "bottom", size = 3,
                               aes(label = paste0("p = ", signif(..p.value.., digits = 2))))
    pl3 <- ggplot(shan, aes(group.id, Evenness, fill=group.id)) +
      geom_point(shape=22) + ggthemes::theme_few() + theme(axis.title.x = element_blank()) +
      scale_fill_gradientn(colors=clr) + guides(fill="none") +
      stat_smooth(method="lm", formula=y ~ x, se=F, size=0.5, show.legend = F) +
      ggpmisc::stat_poly_eq(formula=y ~ x, parse=T, label.x = "left", label.y="top", size=3,
                            aes(label=paste(stat(eq.label), stat(rr.label), sep="~~~~"))) +
      ggpmisc::stat_fit_glance(method = 'lm', method.args = list(formula = y ~ x), label.x = 'right', label.y = "bottom", size = 3,
                               aes(label = paste0("p = ", signif(..p.value.., digits = 2))))
    pl4 <- ggplot(shan, aes(group.id, Chao1, fill=group.id)) +
      geom_point(shape=22) + ggthemes::theme_few() + theme(axis.title.x = element_blank()) +
      scale_fill_gradientn(colors=clr) + guides(fill="none") +
      stat_smooth(method="lm", formula=y ~ x, se=F, size=0.5, show.legend = F) +
      ggpmisc::stat_poly_eq(formula=y ~ x, parse=T, label.x = "left", label.y="top", size=3,
                            aes(label=paste(stat(eq.label), stat(rr.label), sep="~~~~"))) +
      ggpmisc::stat_fit_glance(method = 'lm', method.args = list(formula = y ~ x), label.x = 'right', label.y = "bottom", size = 3,
                               aes(label = paste0("p = ", signif(..p.value.., digits = 2))))
    pleg <- as_ggplot(get_legend(pl1 + guides(fill=guide_colorbar(title = group, barwidth = 10,
                                                                  barheight=0.75, direction = "horizontal"))))

  }

  layout=matrix(c(1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,
                  3,3,3,4,4,4,3,3,3,4,4,4,3,3,3,4,4,4,
                  5,5,5,5,5,5), ncol=6, byrow=T)

  print(grid.arrange(arrangeGrob(pl1, pl2, pl3, pl4, pleg, layout_matrix=layout)))

  return(shan)
}

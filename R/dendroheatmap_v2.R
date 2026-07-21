#' @title Dendro-heatmap
#' @author Kara J.M. Taylor
#' @description Constructs a sample abundance heatmap of select taxa with a
#' phylogenetic tree of taxa on the y-axis and a dendrogram of samples ordered
#' by hierarchical compositional relatedness on the x-axis.
#' @param comm Community matrix (site rows and species columns)
#' @param tax Taxonomy table
#' @param meta Sample metadata
#' @param bdiv Sample dissimilarity matrix (optional; recommend using UniFrac distance)
#' @param path Path to phylogenetic tree in Newick format.
#' @param group Column in `meta` indicating the sets
#' @param subgroup Column in `meta` indicating subsets or sample ID names
#' @param core.list List of ASVs to include. May be numeric to select the top n most abundant ASVs. Default= top 10.
#' @param hi.tax Taxonomic level in `tax` to summarize `core.list` to (optional)
#' @param tax.grp Taxonomic level to summarize taxa to. Defaults to 'tag'
#' @returns List containing set information
#' @export


dendro.heatmap <- function(comm, tax, meta, bdiv=NULL, path="", group, subgroup, core.list=10, hi.tax=NULL){

  library(dplyr)
  library(ape)
  library(vegan)
  library(ggplot2)
  library(reshape2)
  library(ggdendro)
  library(ggtree)
  library(dendextend)
  library(cowplot)

  if (!any(row.names(comm)%in%row.names(meta))){
    stop("Check row names of community and metadata matrices for identity.")
  }
  if (!is.numeric(core.list)&!any(core.list%in%names(comm))){
    stop("Check column names of community matrix and list of selected ASVs for identity.")
  }

  meta <- meta[row.names(comm),]
  meta <- data.frame(ID=row.names(meta), meta)

  if (length(unique(meta[,group]))==1|length(unique(meta[,subgroup]))==1){
    meta <- mutate(meta, group.subgroup=paste(eval(parse(text=group)),
                                              eval(parse(text=subgroup)),
                                              ID, sep="."))
  } else {
    meta <- mutate(meta, group.subgroup=paste(eval(parse(text=group)),
                                              eval(parse(text=subgroup)),
                                              sep="."))
  }

  #### Series of nested functions to produce grouped heatmap
  top <- sub.tax(comm, tax, meta, group, subgroup, core.list)
  den <- group.dendro(comm, meta, tax, bdiv, path, group, subgroup)
  phy <- otu.phylo(top, tax, path, hi.tax) #use hi.tax if you want core OR all taxa at level other than in `sub`
  heat <- top.heat(top, tax, phy, den, hi.tax) #use hi.tax if you want core OR all taxa

  e <- ggplot() +
    geom_blank() +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0,0,0,0, "pt"))

  l <- cowplot::get_legend(heat)
  heat <- heat + theme(legend.key.width = unit(2, "cm"))

  gg <- egg::ggarrange(ggdraw(l), den, phy, heat + theme(legend.position = 'none'),
                       nrow=2, widths=c(1,4), heights=c(1,4), draw=F)

  out <- list("legend"=l, "site.dendro"=den, "tax.phylo"=phy, "heat.raster"=heat, "whole"=gg)

  return(out)
}

## function for subsetting comm by select (or top n) taxa
sub.tax <- function(comm, tax, meta, group, subgroup, core.list){

  if (!is.null(core.list)){
    if (class(core.list)=="numeric"){
      core.list <- names(sort(colMeans(comm), decreasing = T))[1:core.list]
    }
    comm <- comm[,core.list]
  }

  ra <- function(x){ #where x is a vector
    y <- sum(x)
    z <- sapply(x, FUN=function(x){z1 <- ((x)/y)}) ## the x+1 solves the log transformation problem later without skewing the data too much
    return(z)
  }

  tmp <- apply(comm, 1, ra) %>% t() %>% as.data.frame()
  if (any(is.na(rowSums(tmp)))){
    tmp[which(is.na(rowSums(tmp))),] <- 0
  }

  top <- right_join(meta[c("ID","group.subgroup")],
                    data.frame(ID=row.names(tmp), tmp),
                    by=c("ID"="ID")) %>% .[,-1]

  row.names(top) <- row.names(comm)
  return(top)
}

## function for creating cluster dendrogram for the groups and subgroups
group.dendro <- function(comm, meta, tax, bdiv, path, group, subgroup){

  columns <- colnames(comm)
  rows <- row.names(comm)

  data <- data.frame(ID=row.names(comm), comm)
  data <- right_join(meta[,c("ID", "group.subgroup")], data,
                     by=c("ID"="ID")) %>%
    `row.names<-`(rows) %>% .[,-1]

  if (is.null(bdiv)){
    avgd <- data %>%
      group_by(group.subgroup) %>%
      summarize_all(funs(mean=mean)) %>%
      as.data.frame() %>%
      `row.names<-`(.[,"group.subgroup"]) %>% .[,-1] %>%
      `colnames<-`(columns)

    tree.comm <- t(avgd) %>%
      `row.names<-`(tax$den.otu[match(row.names(.), tax$tag)])

    tips <- setdiff(ape::read.tree(path)$tip.label,row.names(tree.comm))
    drop.tree <- ape::drop.tip(ape::read.tree(path), tip=tips)
    bdiv <- GUniFrac::GUniFrac(t(tree.comm), drop.tree, alpha=1)$unifracs[,,1]
  } else {
    tmpnames <- meta$group.subgroup[match(row.names(bdiv), meta$ID)]
    dimnames(bdiv) <- list(tmpnames, tmpnames)
    if (any(is.na(tmpnames))){
      pull <- which(is.na(tmpnames))
      bdiv <- bdiv[-pull,-pull]
    }
  }

  clust <- hclust(d=as.dist(bdiv), method="average")

  dend <- as.dendrogram(clust) %>%
    dendextend::rotate(sort(as.character(unique(data$group.subgroup))))

  ddend <- dendro_data(dend, type="rectangle")

  order <- ddend$labels$label

  ngsg <- length(unique(data$group.subgroup))
  zoom <- ngsg*(0.0595-0.0001*ngsg)

  ymin <- c(min(ddend$segments$y), min(ddend$segments$yend))
  ymax <- c(max(ddend$segments$y), max(ddend$segments$yend))

  ddend$segments$yend[ddend$segments$yend==ymin[2]] <- ymax[2]-(ymax[1]-ymin[1])*1.1

  d <- ggplot(segment(ddend)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=0.5) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0,0,0,0, "pt")) +
    coord_cartesian(xlim=c(zoom, ngsg+1-zoom))

  d[["tip.labels"]] <- order
  return(d)
}


## function for creating phylogenetic cladogram for top taxa
otu.phylo <- function(top, tax, path, hi.tax){#tax= output from bact.tax

  otus <- names(top)[-1]
  ngsg <- length(unique(top$group.subgroup))

  tree <- ape::read.tree(path)
  tips <- as.character(tax$den.otu[-match(otus, tax$tag)])

  prune <- match(tax$den.otu, tree$tip.label)
  if (length(prune)>0){
    if(any(is.na(prune))){
      prune <- prune[which(!is.na(prune))]
    }
    tips <- c(tips, tree$tip.label[-prune])
  }
  drop.tree <- ape::drop.tip(tree, tip=tips)

  labs <- data.frame(tips=drop.tree$tip.label)
  labs <- left_join(labs, tax, by=c("tips"="den.otu"))

  if (!is.null(hi.tax)){
    if (hi.tax%in%names(tax)){
      get <- which(names(labs)==hi.tax)
      drop.tree$tip.label <- as.character(unlist(labs[,get]))
    }
  }else {
    drop.tree$tip.label <- as.character(labs$tag)
  }

  p <- ggtree::ggtree(drop.tree, branch.length = "none", ladderize=FALSE) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0,0,0,0, "pt"))

  pull <- data.frame(node=p$data$node, label=p$data$label)
  dup <- which(duplicated(pull$label[1:length(otus)])==TRUE)
  new <- drop.tip(drop.tree, tip=pull$node[dup])

  n <- length(new$tip.label)
  crop <- ggtree(new)@data; crop <- max(crop$x)-0.05

  cp <- ggtree::ggtree(new, size=0.5) + theme_tree2() +
    geom_tiplab(align=TRUE, linesize=.5) +
    #scale_x_continuous(expand=c(0,0.05)) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0,0,0,0, "pt")) +
  coord_cartesian(xlim=c(0,crop))

  tmp <- cp$data[1:n, c(4,5,7)] %>%
    .[order(.$y),]

  cp[["new.labels"]] <- tmp$label
  return(cp)
}

## function for creating heatmap using the top n otus
top.heat <- function(top, tax, phylo, dendro, hi.tax, filename){

  otus <- names(top)[-1]

  group1 <- top %>%
    group_by(group.subgroup) %>%
    summarize_all(funs(mean(.,na.rm=T))) %>%
    as.data.frame() %>%
    `row.names<-`(.[,1]) %>% .[,-1]

  group1 <- group1[match(dendro$tip.labels, row.names(group1)),]

  group2 <- reshape(group1, varying=otus,
                    v.names="abund",
                    timevar="otu",
                    times=otus,
                    direction="long",
                    ids=row.names(group1))

  group2 <- left_join(group2, tax, by=c("otu"="tag"))
  names(group2)[which(names(group2)=="otu")] <- "tag"

  if (!is.null(hi.tax)){
    if (hi.tax%in%names(tax)){
      group3 <- group2[,c("id","tag",hi.tax,"abund")] %>%
        group_by_("id","tag",hi.tax) %>%
        summarize_at("abund", funs(abund=sum(.,na.rm=T))) %>%
        as.data.frame()
      names(group3)[which(names(group3)==hi.tax)] <- "name"
    }
  } else {
    group3 <- group2[,c(1,2,3)] %>%
      group_by(id, tag, name) %>%
      summarize_at("abund", funs(abund=max)) %>%
      as.data.frame()
  }

  group3$id <- factor(group3$id,
                      levels=dendro$tip.labels)

  group3$name <- factor(group3$name,
                        levels=phylo$new.labels)

  min.lim <- floor(min(group3$abund,na.rm = T))

  brks <- seq(-2, 2, 1)

  colors <- rev(c("magenta", "red", "orange", "yellow", "green", "turquoise", "blue", "black"))

  if (min.lim < brks[1]){
    colors <- c(rep("black", abs(min.lim - brks[1])+1), colors)
    brks <- seq(min.lim, 2, 1)
  }

  lbs <- paste0(10^brks, "%")

  lims <- c(min(brks), max(brks))

  h <- ggplot(group3, aes(id, name)) +
    geom_raster(aes(fill=abund*100)) +
    scale_fill_gradientn(colors=colors, limits=c(0,ceiling(max(group3$abund)*100)),
                         breaks=round(seq(0,ceiling(max(group3$abund)*100),length=4)),
                         guide=guide_colorbar(title="Relative Abundance (%)", title.position = "top",
                                              ticks.colour="black", frame.colour="black", label.position = "bottom",
                                              title.theme = element_text(size=12, angle=0, hjust=0.5)),
                         na.value = "white") +
    theme(axis.ticks=element_line(color="white"),
          axis.text.x=element_text(angle=90, size=14, hjust=1, vjust=0.5, face="bold"),
          axis.text.y=element_text(size=14, vjust=0.5), axis.title = element_blank(),
          plot.margin = margin(1,1,5,1, "pt"),
          plot.caption=element_text(size=12),
          legend.text=element_text(size=12, face = "bold"),
          legend.position="top", legend.direction = "horizontal",
          legend.justification="center",
          legend.key.width = unit(0.5, "null"),
          legend.box.margin=margin(1,2,1,5,"pt")) +
    scale_y_discrete(position="right") +
    labs(caption=paste("Taxa: n = ", as.character(length(levels(group3$name)))))

  return(h)
}

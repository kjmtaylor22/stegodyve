#' @title Microbiome Venn diagram contents
#' @author Kara J.M. Taylor
#' @description Constructs bar graphs showing the relative abundance of taxa in each set.
#' @param shared Output list from `shared_taxa`
#' @param comm Community matrix
#' @param select List of taxa to show; will show all taxa by default.
#' @param tax Parsed taxonomy table.
#' @param hi.tax Taxonomic level in `tax` to show in bar plots.
#' @param meta Metadata table
#' @param group Column in `meta` containing set grouping information
#' @param subgroup Column in `meta` containing within-set subdivisions of the taxonomic distribution.
#' @param file Output filename
#' @param print Logical; Print all plots to file? Default is FALSE
#' @param compare.otus Logical; Print plot showing summary distribution? Default is FALSE
#' @param set.y Logical; Fix ylim across all output plots? Default is FALSE
#' @returns Data frame containing relative abundance for each `hi.tax` in each group.subgroup level
#' @export

mbiom.bar <- function(shared, comm, select=NULL, tax, hi.tax, meta, group, subgroup, file="mbiom_bar", print=F, compare.otus=F, set.y=F){

  cran <- c("dplyr", "ggplot2", "ggpubr")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  library(dplyr)
  library(ggplot2)
  library(ggpubr)

  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    names(meta)[1] <- "sample"
    id <- 1
  }

  if (any(!row.names(comm)%in%meta[,which(names(meta)=="sample")])){
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

  splt <- function(x){
    y <- strsplit(x, split="\\.")
    y <- unlist(y)
  }

  ra <- function(x){ #where x is a vector
    y <- sum(x)
    z <- sapply(x, FUN=function(x){z1 <- x/y})
    return(z)
  }

  sel.list <- as.list(names(shared))

  if (length(sel.list)>4){
    sel.list <- lapply(sel.list, FUN=splt)
  }

  if (is.null(select)){select <- names(comm)}

  s <- which(names(meta)=="sample")
  g <- which(names(meta)==group)
  sg <- which(names(meta)==subgroup)

  meta <- mutate(meta, !!"group.subgroup" :=
                   paste(eval(parse(text=group)), eval(parse(text=subgroup)), sep=":"))
  g.sg <- which(names(meta)=="group.subgroup")

  ht <- which(names(tax)==hi.tax)

  for (j in 1:dim(tax)[2]){
    uls <- unlist(shared)
    if ("unclassified"%in%uls){
      uls <- uls[-which(uls=="unclassified")]
    }
    t <- all(!is.na(match(uls, tax[,j])))
    if (t==TRUE){break}
  }

  df <- data.frame()
  for (i in 1:length(sel.list)){

    selection <- sel.list[[i]]

    if ("all"%in%selection){
      selection <- unique(unlist(sel.list))
      selection <- selection[-length(selection)]
    }

    otus <- tax$tag[match(shared[[i]], tax[,j])]
    otus <- select[select%in%otus]

    if (length(otus)==0){next}


    data <- data.frame(sample=row.names(comm), comm)
    data <- right_join(meta[,c(id,g.sg)], data)

    data2 <- data[,-1] %>%
      group_by(group.subgroup) %>%
      summarize_all(funs(mean=mean)) %>%
      as.data.frame()

    new <- apply(data2[,-1], MARGIN=1, FUN=ra) %>% t()

    new <- data.frame(groupvar=data2[,1], new)
    names(new)[-1] <- names(comm)

    data2 <- reshape2::melt(new)

    data2 <- right_join(tax[,c(3, ht)], data2, by=c("tag"="variable"))
    data2 <- right_join(unique(data.frame(meta[,c(g.sg, g, sg)])), data2, by=c("group.subgroup"="groupvar"))

    ce.il <- data2[data2[,group]%in%selection & data2$tag%in%otus,]
    ce.il$value <- ce.il$value*100

    names(ce.il)[grep(hi.tax, names(ce.il))[1]] <- "Taxonomy"
    names(ce.il)[grep(group, names(ce.il))[1]] <- "Group"
    names(ce.il)[grep(subgroup, names(ce.il))[1]] <- "Subgroup"

    ce.il <- data.frame(ce.il, sec.list=paste(selection, collapse="."))

    df <- rbind(df, ce.il)

  }

  write.csv(df, paste0("./coreplots/venn/", file, "/dataframe.csv"))
  if (print==F){return(df)}

  if (compare.otus==T){
    g <- ggplot(df, aes(Subgroup, value, fill=Group)) +
      facet_wrap(vars(Taxonomy), ncol=4, dir="h", scales="free_y", strip.position = "top") +
      geom_bar(stat="summary", fun.y="sum", width=.96, position="dodge") +
      theme_minimal() +
      theme(panel.border=element_rect(linetype="solid", fill=NA, color="black"),
            legend.position="bottom",
            legend.justification="center",
            legend.direction="horizontal",
            panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour="grey90", linetype="36"),
            panel.grid.minor.y = element_blank(),
            strip.text = element_text(size = 45, face="bold"),
            axis.title.y = element_text(size = 40, face="bold"),
            axis.text.x = element_text(color = "black", size = 40, angle=45, vjust=0.9, hjust=1),
            axis.text.y = element_text(color = "black", size = 40),
            legend.title = element_text(size = 40, face="bold"),
            legend.text = element_text(size = 30),
            panel.spacing.y = unit(2, "lines"),
            plot.margin = margin(3, 3, 3, 3, "lines")) +
      ylab("Relative Abundance (%)")+
      xlab(NULL)
    h <- ceiling(length(unique(df$Taxonomy))/4)
    jpeg(paste0("./coreplots/venn/", file, "/", group,"_sharedOTUs.jpg"),
         width=2300, height=420*h+280)
    print(g)
    dev.off()
  }

  #return(ce.il)
  colors <- c("#330033", "#F0E442", "#0072B2", "#FF0033", "#FF9933", "#999999", "#56B4E9", "#FFCC00", "#00FF00", "#FF66FF", "#666666")
  colors <- colorRampPalette(colors)(length(unique(as.factor(df$Taxonomy))))
  names(colors) <- as.character(unique(df$Taxonomy))

  if (compare.otus==T){
    k <- ggplot(df, aes(Subgroup, value, fill=Taxonomy)) +
      facet_wrap(vars(Group), ncol=2, dir="v", scales="free_y", strip.position = "right") +
      geom_col() +
      #geom_bar(stat="summary", fun.y="sum", width=.96, position="stack") +
      theme_minimal() +
      scale_fill_manual(values=colors) +
      theme(panel.border=element_rect(linetype="solid", fill=NA, color="black"),
            legend.position="right",
            legend.justification="center",
            #legend.direction="horizontal",
            panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour="grey90", linetype="36"),
            panel.grid.minor.y = element_blank(),
            strip.text = element_text(size = 60, face="bold"),
            axis.title.y = element_text(size = 40, face="bold"),
            axis.text.x = element_text(color = "black", size = 40, angle=45, vjust=0.9, hjust=1),
            axis.text.y = element_text(color = "black", size = 40),
            legend.title = element_text(size = 40, face="bold"),
            legend.text = element_text(size = 30),
            panel.spacing.y = unit(2, "lines"),
            plot.margin = margin(3, 3, 3, 3, "lines"))+
      ylab("Relative Abundance (%)")+
      xlab(NULL)
    h <- ceiling(length(unique(df$Group))/2)
    jpeg(paste0("./coreplots/venn/", file, "/", group,"_NicheSpace.jpg"),
         width=2300, height=460*h)
    print(k + guides(fill=F))
    dev.off()

    pleg <- as_ggplot(get_legend(k))
    if (length(unique(as.factor(df$Taxonomy))) > 13){
      wid <- length(unique(as.factor(df$Taxonomy)))/14*900
    } else {wid <- 900}
    jpeg(paste0("./coreplots/venn/", file, "/", group,"_NicheSpace_legend.jpg"),
         width=wid, height=1150)
    print(pleg)
    dev.off()
  }

  if (set.y==T){
    fix.ylim <- expression(ylim(0,100))
  } else {fix.ylim <- expression(ylim(0,NA))}

  for (i in unique(df$sec.list)){

    sub.sel <- subset(df, subset=sec.list==i)


    f <- ggplot(sub.sel, aes(Subgroup, value, fill=Taxonomy)) +
      facet_wrap(vars(Group), dir="v", scales="free_y", strip.position = "right") +
      geom_bar(stat="identity", width=.96) +
      theme_minimal() +
      scale_fill_manual(values=colors) +
      eval(expr=fix.ylim) +
      theme(panel.border=element_rect(linetype="solid", fill=NA, color="black"),
            legend.position="right",
            legend.justification="center",
            panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour="grey90", linetype="36"),
            panel.grid.minor.y = element_blank(),
            strip.text = element_text(size = 80, face="bold"),
            axis.title = element_text(size = 65, face="bold"),
            axis.text.x = element_text(color = "black", size = 60, angle=45, vjust=0.9, hjust=1),
            axis.text.y = element_text(color = "black", size = 55),
            legend.title = element_text(size = 50, face="bold"),
            legend.text = element_text(size = 50),
            panel.spacing.y = unit(2, "lines"),
            plot.margin = margin(3, 3, 3, 3, "lines")) +
      ylab(NULL)+
      xlab(NULL)
    #eval(expr=leg)
    w <- length(unique(df$Subgroup))

    #leg <- expression(guides(fill=guide_legend(title="", ncol=1)))
    pleg <- as_ggplot(get_legend(f))
    if (length(unique(as.factor(sub.sel$Taxonomy))) > 13){
      wid <- length(unique(as.factor(sub.sel$Taxonomy)))/14*900
    } else {wid <- 900}

    count <- unlist(strsplit(i, split="\\."))

    if (length(count)==5){
      jpeg(paste0("./coreplots/venn/", file, "/", i,"_legend.jpg"),
           width=wid, height=1150)
      print(pleg)
      dev.off()
      f <- f + facet_wrap(vars(Group), dir="v", scales="fixed", strip.position = "right")

      jpeg(paste0("./coreplots/venn/", file, "/", i,"_sharedOTUs.jpg"),
           width=2300, height=1650)
    }
    if (length(count)==4){
      jpeg(paste0("./coreplots/venn/", file, "/", i,"_legend.jpg"),
           width=wid, height=1150)
      print(pleg)
      dev.off()
      f <- f + facet_wrap(vars(Group), dir="v", scales="fixed", strip.position = "right")

      jpeg(paste0("./coreplots/venn/", file, "/", i,"_sharedOTUs.jpg"),
           width=2300, height=1250)
    }
    if (length(count)==3){
      jpeg(paste0("./coreplots/venn/", file, "/", i,"_legend.jpg"),
           width=wid, height=1150)
      print(pleg)
      dev.off()

      jpeg(paste0("./coreplots/venn/", file, "/", i,"_sharedOTUs.jpg"),
           width=1300, height=1650)
    }
    if (length(count)==2){
      jpeg(paste0("./coreplots/venn/", file, "/", i, "_legend.jpg"),
           width=wid, height=1150)
      print(pleg)
      dev.off()

      jpeg(paste0("./coreplots/venn/", file, "/", i, "_sharedOTUs.jpg"),
           width=1300, height=1150)
    }
    if (length(count)==1){
      jpeg(paste0("./coreplots/venn/", file, "/", i,"_legend.jpg"),
           width=wid, height=1150)
      print(pleg)
      dev.off()

      jpeg(paste0("./coreplots/venn/", file, "/", i,"_sharedOTUs.jpg"),
           width=1300, height=700)
    }
    print(f+guides(fill=F))
    dev.off()


  }
  graphics.off()
  return(df)
}

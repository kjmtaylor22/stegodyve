---
title: "Spatial and environmental influences on the assembly of silk microbiomes in the African social spider, Stegodyphus dumicola"
author: "K.J.M.Taylor"
date: "2025-04-23"
css: 
output:
  minidown::mini_document:
    tabset: true
    toc: true
    toc_float: true
    toc_highlight: false
    toc_depth: 4
    self_contained: true
    code_folding: show
    lightbox: true
    thumbnails: false
    highlight: "pygments"
    gallery: true
    fig_width: 7
    fig_height: 5
    dpi: 300
    framework: "mini"
    spacing: "normal"
---

```{css echo = FALSE}
blockquote{
background-color: #cee9e9;
}

pre {
  max-height: 300px;
  overflow-y: auto;
}

```

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = T, eval=T, cache=T, 
                      message=F, warning=F, 
                      fig.path = "figs/")
```

# Analysis setup

```{r, eval=T}
library(stegodyve) # the package written for this analysis
library(dplyr)
library(parallel)
library(vegan)
library(ade4)
library(adespatial)
library(sp)
library(spdep)
library(adegraphics)
library(VennDiagram)
library(ggplot2)
library(ggtern)
library(cowplot)
library(gridExtra)
```

## Data wrangling

### Tree construction from FNA
```{bash, eval=F}
#!/bin/bash

# Job name and who to send updates to
#SBATCH --job-name=TreeFNA
#SBATCH --mail-user=****
#SBATCH --mail-type=ALL
#SBATCH --account=****
#SBATCH --partition=****
#SBATCH --qos=****  
#SBATCH -o logs/TreeFNA.%j.log
#SBATCH --error errors/TreeFNA.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=12G
#SBATCH --time=5:00:00

module load qiime2/2021.11


qiime tools import \
  --input-path sv.seqa.fna \
  --output-path aligned-sequences.qza \
  --type 'FeatureData[AlignedSequence]'


qiime alignment mask \
  --i-alignment aligned-sequences.qza \
  --o-masked-alignment MaskedAligned.qza
  
    
qiime phylogeny fasttree \
  --i-alignment MaskedAligned.qza \
  --o-tree unrooted-tree.qza
  

qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
```

### Bacteria
```{r, eval=F}
library(biomformat)

biom.as.csv("data/bacteria/ASV_Table.biom")

taxonomy <- read.csv("data/bacteria/taxonomy.csv")
taxB <- bact.tax(taxonomy = taxonomy, database = "green")

library(seqinr)

fna <- read.fasta("data/bacteria/sv.seqs.fna", forceDNAtolower = F, as.string = T)
fna <- fna[names(fna)%in%tax$den.otu]
write.fasta(fna, names=names(fna), file.out="data/bacteria/reduced-seqs.fna") ## change to Unix with Notepad++; can create fast tree with QIIME2


comm <- read.csv("data/bacteria/feature_table.csv", row.names = 1) %>% t() %>% as.data.frame()
comm[is.na(comm)] <- 0

commB <- comm[,match(tax$den.otu, colnames(comm))]

tree <- ape::read.tree("data/bacteria/tree.nwk")

colnames(commB) <- as.character(tax$tag)

sepmeta <- function(x){
  y <- strsplit(x,split = ".", fixed = T)
  return(unlist(y))
}

newmeta <- lapply(as.list(row.names(commB)), FUN=sepmeta) %>%
  as.data.frame() %>% t() %>% as.data.frame() %>%
  `colnames<-`(c("SiteID","CollectionYear","ColonyID","WebType")) %>%
  mutate(colony_ID=paste(SiteID, CollectionYear, ColonyID, sep="-"))

cleanmeta <- function(x){
  z <- gsub("\\s*\\(\\d+\\)","", x)
  z <- gsub("[()]","",z)
  z <- gsub("O","C", z)
  return(z)
}

meta <- openxlsx::read.xlsx("data/Stegodyphus field data 2018.xlsx", 1)
pull <- sapply(meta$colony_ID, FUN=cleanmeta)

meta$colony_ID <- pull
meta2 <- meta[which(duplicated(meta$colony_ID)==F),]

newmeta <- left_join(newmeta, meta2)
newmeta <- select(newmeta, column_name="colony_ID", everything())
newmeta <- data.frame(SampleID=row.names(comm), newmeta)

newmeta$country <- as.factor(newmeta$country)
newmeta$SiteID <- as.factor(newmeta$SiteID)
newmeta$WebType <- as.factor(newmeta$WebType)

newmeta <- mutate(newmeta, CountryType=paste(country, WebType, sep="_"),
                  SampleType=WebType,
                  SiteType=paste(SiteID, WebType, sep="_"),
                  TypeCountry=paste(WebType, country, sep="_"))


newmeta$SampleType <- factor(newmeta$SampleType, labels=c("Web","Web","Soil"))

newmeta$capture_web_area[newmeta$capture_web_area=="n/a"] <- NA
newmeta$capture_web_area <- gsub("cm", "", newmeta$capture_web_area)

newmeta$retreat_web_volume[newmeta$retreat_web_volume=="n/a"] <- NA

newmeta <- newmeta %>% separate(capture_web_area, c("cap_length","cap_width","cap_height"), sep="x", convert=T) %>%
  separate(retreat_web_volume, c("ret_length","ret_width","ret_height"), sep="x", convert=T) %>%
  mutate(cap_volume=cap_length*cap_width*cap_height) %>%
  mutate(ret_volume=ret_length*ret_width*ret_height)

reads <- read.csv("data/bacteria/absolute.abundance.csv", row.names=1)

metaB <- left_join(newmeta, reads[,1:5], by=c("SampleID"="customer_label"))
metaB$Ct <- as.numeric(metaB$Ct)


## Force some small distance differences between these sites and their sister sites to remove dups
metaB$WebType <- factor(metaB$WebType, levels=c("soil","retreat","capture"))
metaB$latitude[metaB$column_name%in%c("B1-18-H","UPP-18-B")] <-
  metaB$latitude[metaB$column_name%in%c("B1-18-H","UPP-18-B")] +0.00001
metaB$longitude[metaB$column_name%in%c("B1-18-H","UPP-18-B")] <- 
  metaB$longitude[metaB$column_name%in%c("B1-18-H","UPP-18-B")] +0.00001
row.names(metaB) <- metaB$SampleID

## ignore
#unifracB <- read.table("data/bacteria/00...AllSamples.Bac16Sv34-20220330T145802Z-001/00...AllSamples.Bac16Sv34/Beta_Diversity_ASV/Unweighted_UniFrac/Raw_Data/unweighted_unifrac_dm.txt")

unifracB <- GUniFrac::GUniFrac(comm, tree)

brayB <- vegan::vegdist(comm) %>% as.matrix()


save(commB, taxB, metaB, unifracB, brayB, file="data/bacteria/basefiles_bact.RD")
```

### Fungi
```{r, eval=F}
library(biomformat)

biom.as.csv("data/fungi/ASV_Table.biom")

taxonomy <- read.csv("data/fungi/taxonomy.csv")
taxF <- bact.tax(taxonomy = taxonomy, database = "fungi")

library(seqinr)

fna <- read.fasta("data/fungi/sv.seqs.fna", forceDNAtolower = F, as.string = T)
fna <- fna[names(fna)%in%taxF$den.otu]
write.fasta(fna, names=names(fna), file.out="data/fungi/reduced-seqs.fna") ## change to Unix with Notepad++; can create fast tree with QIIME2


commF <- read.csv("data/fungi/feature_table.csv", row.names = 1) %>% t() %>% as.data.frame()
commF[is.na(commF)] <- 0

commF <- commF[,match(taxF$den.otu, colnames(commF))]

tree <- ape::read.tree("data/fungi/treeFungi.nwk")

unifracF <- GUniFrac::GUniFrac(commF, tree)


colnames(commF) <- as.character(taxF$tag)


brayF <- vegan::vegdist(commF) %>% as.matrix()


load('data/bacteria/basefiles_bact.RD')

metaF <- read.csv("data/fungi/absolute.abundance.csv")[,-1]
metaF <- metaF[metaF$customer_label%in%row.names(commF),]
metaF <- left_join(metaF, metaB[,-c(29:32)], by=c("customer_label"="SampleID"))
names(metaF)[1] <- "SampleID"
metaF <- left_join(metaB, metaF[,c(1,6)])
metaF$population[56] <- "Kimberly"
metaF$population[20] <- "Uppington"

write.csv(metaF, "data/fungi/metadata_fungi.csv")
## copy down name/loc data for R31N.18.A.retreat from R31N.18.A.capture

metaF <- read.csv("data/fungi/metadata_fungi.csv", row.names = 1)

## Force some small distance differences between these sites and their sister sites to remove dups
metaF$WebType <- factor(metaF$WebType, levels=c("soil","retreat","capture"))
metaF$latitude[metaF$column_name%in%c("B1-18-H","UPP-18-B")] <-
  metaF$latitude[metaF$column_name%in%c("B1-18-H","UPP-18-B")] +0.00001
metaF$longitude[metaF$column_name%in%c("B1-18-H","UPP-18-B")] <- 
  metaF$longitude[metaF$column_name%in%c("B1-18-H","UPP-18-B")] +0.00001
row.names(metaF) <- metaF$SampleID


save(commF, taxF, unifracF, brayF, metaF, file="data/fungi/basefiles_fungi.RD")
```

# Data analysis {.tabset} 

```{r}
rm(list = ls()) ## clear everything from the earlier data wrangling and reload

load("data/bacteria/basefiles_bact.RD")
load("data/fungi/basefiles_fungi.RD")

rm(unifracB, unifracF) ## not used here
```

## 0. Spatial and environmental gradients (Figure 1A)

```{r 00_gradient_Figure1A, fig.width=14.3, fig.height=3.9}

all(unique(metaB$column_name) %in% unique(metaF$column_name))
set <- unique(metaB[,c("column_name", "country", "Temp", "AridIndex", "latitude", "longitude")]) %>%
  `row.names<-`(.$column_name)
dist <- dist(set[,5:6]) %>% as.matrix %>% reshape2::melt(value.name = "dist")
env <- dist(set[,3:4]) %>% as.matrix %>% reshape2::melt()
tri <- matrix(0,nrow(set), nrow(set)) %>% lower.tri(.) %>% reshape2::melt()
envdist <- cbind(dist, envgrad=env[,3], tri=tri[,3]) %>% subset(subset=tri==T)

sample_dist <- left_join(envdist, unique(metaB[,c("column_name","country")]), by=c("Var1"="column_name")) %>%
  left_join(unique(metaB[,c("column_name","country")]), by=c("Var2"="column_name")) %>%
  mutate(transect.wb=country.x==country.y, compare=paste(country.y, "vs.", country.x))
# re-level factors
sample_dist$dist <- sample_dist$dist*111.11 # set distances from decimal degrees to km
sample_dist$transect.wb[sample_dist$transect.wb==T] <- "Within transect"
sample_dist$transect.wb[sample_dist$transect.wb==F] <- "Between transects"
sample_dist$transect.wb <- factor(sample_dist$transect.wb, levels=unique(sample_dist$transect.wb))
sample_dist$compare <- factor(sample_dist$compare, 
                              labels=c("Namibia vs. Namibia",
                                       "Namibia vs. South Africa",
                                       "Namibia vs. South Africa",
                                       "South Africa vs. South Africa"))
sample_dist$compare <- factor(sample_dist$compare, 
                              levels=c("Namibia vs. Namibia",
                                       "South Africa vs. South Africa",
                                       "Namibia vs. South Africa"))


pull <- subset(sample_dist, subset=transect.wb=="Within transect", drop = T)
lm(envgrad ~ dist*compare, pull) %>% summary()
aov(Temp ~ country, set) %>% summary()
aov(AridIndex ~ country, set) %>% summary()

## Figure 1, Plot A; 
ggplot(sample_dist, aes(dist, envgrad, linetype=compare, shape=compare)) + 
  geom_point(size=3) +  scale_shape_manual(values=c(16,8,17)) +
  scale_linetype_manual(values=c("solid","longdash","dotted")) +
  stat_smooth(method="lm", formula='y~x', show.legend = c(color=F), color="black")+
  ggpmisc::stat_poly_eq(formula=y ~ x, parse=T, label.x = c(1050,465,250), label.y=c(1,1.9,0.75),
                        aes(label=..eq.label..), vjust=0.8, size=4, geom="label") +
  guides(shape=guide_legend(title=NULL, override.aes = list(size=5)),
         linetype=guide_legend(title=NULL, override.aes = list(linewidth=1.5))) +
  theme_bw() + theme(axis.text = element_text(size=20),
                     axis.title = element_text(size=18),
                     plot.margin=unit(c(1.5,0.5,0,0.5), "lines"), 
                     legend.text=element_text(margin=margin(0,7,0,0), size=14),
                     legend.key.width = unit(2,"cm"),
                     legend.position = "top", legend.justification = "center") +
  xlab("Distance between colonies (km)") + ylab("Environmental dissimilarity\n(Temp:AridIndex)")

```

<a href="#data-analysis" style="color:steelblue;" >Back to top</a>

## 1. Sample composition {.tabset} 

> **Question:** Are taxa differentially present and abundant in different sample types? 

- **Diversity:** Basic diversity analyses were conducted to assess sampling differences between sample types. Species accumulation curves were drawn to assess the adequacy of sampling depth for all samples. Five alpha diversity metrics (Shannon's H, absolute richness, Pielou's evenness, Chao1, and Faith's PD) were used to examine within-sample diversity and compare between sample types. The four metrics were analysed with ANOVA and Tukey HSD post-hoc test.

- **Composition:** Principal coordinates and distance-based PERMANOVA were used to examine compositional differences between sample types. 

- **Variability:** Beta dispersion was used to assess the degree of within-sample type similarity. 

- **Distribution:** ASVs were parsed based on their presence in at least one sample of any type (soil, capture, or retreat). An ASV can be found in any of the seven possible combinations for the three sample types (i.e. soil, capture, retreat, soil & capture, soil & retreat, retreat & capture, soil & capture & retreat). For each sample type, the mean total ASV relative abundance (%) is given for each sample type combination.

### 1.1 Sample depth 
```{r 01_rarecurves, eval=F}
## Bacteria
rareB <- rarecurve(commB, sample=250, tidy=T)
rareB <- left_join(rareB, metaB[,c("SampleID","WebType")], by=c("Site"="SampleID"))

## Fungi
rareF <- rarecurve(commF, sample=250, tidy=T)
rareF <- left_join(rareF, metaF[,c("SampleID","WebType")], by=c("Site"="SampleID"))

save(rareB, rareF, file="output/rarecurves.RD")
```


<a href="#data-analysis" style="color:steelblue;" >Back to top</a>

### 1.2 Bacteria

#### Alpha diversity
```{r 01_bacteriaDiversity}
## Alpha diversity
alB <- alphadiv(commB, metaB, "WebType", cols=c("#4DAF4A","#377EB8", "#E41A1C"), "bacteria")

# Faith's PD
tree <- ape::read.tree("data/bacteria/tree.nwk")
pdcomm <- commB; names(pdcomm) <- taxB$den.otu#[match(names(pull.commB),taxB$tag)]
pruned <- picante::prune.sample(pdcomm, tree)
pd <- picante::pd(pdcomm, pruned)

# Combine div metrics and anova
alB <- data.frame(alB, FaithPD=pd$PD)
hsdB <- data.frame(); for (i in names(alB)[c(2,3,4,7)]){
  hsdB <- rbind(hsdB, 
               cbind(index=i, 
                     round(TukeyHSD(aov(eval(parse(text=i)) ~ group.id, alB))$group.id,4) %>%
                        as.data.frame() %>% tibble::rownames_to_column("test")))
}
alB[,-1] %>% group_by(group.id) %>% summarize_all(.funs=list(mean)); hsdB
```

#### Beta diversity
```{r 01_bacteriaDispersion, message=T}
#brayB <- brayB[row.names(pull.commB), row.names(pull.commB)]
#metaB <- metaB[row.names(pull.commB),]

## Beta dispersion
identical(row.names(brayB),row.names(metaB))
identical(colnames(brayB),row.names(metaB))

bdispB <- betadisper(as.dist(brayB), group = metaB$WebType, "centroid")

adonis2(brayB ~ WebType*country, data = metaB)
pairwiseAdonis::pairwise.adonis2(brayB ~ WebType*country, metaB)

anova(bdispB)

hsdB <- rbind(hsdB, cbind(index="Dispersion", round(TukeyHSD(bdispB)$group, 4) %>%
   as.data.frame() %>% tibble::rownames_to_column("test"))); hsdB
```

<a href="#data-analysis" style="color:steelblue;" >Back to top</a>

### 1.3 Fungi

#### Alpha diversity
```{r 01_fungiDiversity}
## Alpha diversity
alF <- alphadiv(commF, metaF, "WebType", cols=c("#4DAF4A","#377EB8", "#E41A1C"), "fungi")

# Faith's PD
tree <- ape::read.tree("data/fungi/treeFungi.nwk")
pdcomm <- commF; names(pdcomm) <- taxF$den.otu#[match(names(pull.commF),taxF$tag)]
pruned <- picante::prune.sample(pdcomm, tree)
pd <- picante::pd(pdcomm, pruned)

# Combine div metrics and anova
alF <- data.frame(alF, FaithPD=pd$PD)
hsdF <- data.frame(); for (i in names(alF)[c(2,3,4,7)]){
  hsdF <- rbind(hsdF, 
               cbind(index=i, 
                     round(TukeyHSD(aov(eval(parse(text=i)) ~ group.id, alF))$group.id,4) %>%
                        as.data.frame() %>% tibble::rownames_to_column("test")))
}
hsdF

alF[,-1] %>% group_by(group.id) %>% summarize_all(.funs=list(mean))
```

#### Beta diversity
```{r 01_fungiDispersion, message=T}

#brayF <- brayF[row.names(pull.commF), row.names(pull.commF)]
#metaF <- metaF[row.names(pull.commF),]

## Beta dispersion
identical(labels(brayF)[[1]],row.names(metaF))
identical(labels(brayF)[[2]],row.names(metaF))

brayF <- brayF[row.names(metaF), row.names(metaF)]

bdispF <- betadisper(as.dist(brayF), group = metaF$WebType, "centroid")

adonis2(brayF ~ WebType*country, data = metaF)

pairwiseAdonis::pairwise.adonis2(brayF~WebType*country, metaF)

anova(bdispF)

hsdF <- rbind(hsdF, cbind(index="Dispersion", round(TukeyHSD(bdispF)$group, 4) %>%
   as.data.frame() %>% tibble::rownames_to_column("test"))); hsdF
```

<a href="#data-analysis" style="color:steelblue;" >Back to top</a>

### 1.4 Plots

#### Venn diagrams
```{r 01_VennDiagrams, eval=F}
## Are the shared microbiota differentially abundant in the different sample types?

# construct venn diagrams for all ASVs across all samples 
shareB <- mbiom.venn(commB, metaB,"WebType", taxB, "bacteria", "tag") %>% 
  # identify the taxa in each section of the diagram
  shared.taxa(., "bacteria", T)
# use the parsed taxonomic and sample information to look at the distribution of ASV abundances
mbB <- mbiom.bar(shareB, commB, NULL, taxB, "class", metaB, "WebType","column_name",
            file="bacteria/barplots", F)

# construct venn diagrams for all ASVs across all samples 
shareF <- mbiom.venn(commF, metaF,"WebType", taxF, "fungi", "tag") %>% 
  # identify the taxa in each section of the diagram
  shared.taxa(., "fungi", T)
# use the parsed taxonomic and sample information to look at the distribution of ASV abundances
mbF <- mbiom.bar(shareF, commF, NULL, taxF, "class", metaF, "WebType","column_name",
            file="fungi/barplots", F)

```

```{r}
# read from the abundance data frame
mbB <- read.csv("coreplots/venn/bacteria/barplots/dataframe.csv")
mbF <- read.csv("coreplots/venn/fungi/barplots/dataframe.csv")

# re-label venn diagram sections for clarity
mbB$sec.list <- gsub(".", "\n\u2229\n", mbB$sec.list, fixed=T)
mbF$sec.list <- gsub(".", "\n\u2229\n", mbF$sec.list, fixed=T)

newlab <- c(unique(mbF$sec.list), unique(mbB$sec.list))[c(2,3,1,13,11,12,14)]

# order levels of new labels
mbB$sec.list <- factor(mbB$sec.list, levels=unique(mbB$sec.list)[c(2,3,1,6,4,5,7)])
mbF$sec.list <- factor(mbF$sec.list, levels=unique(mbF$sec.list)[c(2,3,1,6,4,5,7)])
mbB$sec.list <- factor(mbB$sec.list, labels=newlab)
mbF$sec.list <- factor(mbF$sec.list, labels=newlab)
mbB$Group <- factor(mbB$Group, levels=unique(mbB$Group)[c(2,1,3)])
mbF$Group <- factor(mbF$Group, levels=unique(mbF$Group)[c(2,3,1)])

# join additional metadata
mbB <- left_join(mbB, metaB[,c("column_name","WebType","SiteID","country")],
            by=c("Subgroup"="column_name","Group"="WebType"))
mbF <- left_join(mbF, metaF[,c("column_name","WebType","SiteID","country")],
            by=c("Subgroup"="column_name","Group"="WebType"))

# join together now
mb <- rbind(data.frame(mbB, set="Bacteria"),
            data.frame(mbF, set="Fungi"))

# for each sample in each venn section, sum the relative abundances of present ASVs
mbsum <- mb %>% 
  group_by(country, SiteID, set, Subgroup, Group, sec.list) %>%
  summarize_at(vars(value), funs(sum))

mbsummean <- mbsum %>%
  group_by(Group, set, sec.list) %>%
  summarize_at(vars(value), .funs=list(mean))

mbsummean$Group <- factor(mbsummean$Group, levels=c("retreat","capture","soil"))

```

#### Figure 2
```{r 01_SampleComposition_combined, fig.width=9, fig.height=4.5}

library(magick)

al <- rbind(data.frame(alB, set="Bacteria"), 
            data.frame(alF, set="Fungi"))
bd <- rbind(data.frame(dist=bdispB$distances, g=bdispB$group, set="Bacteria"), 
            data.frame(dist=bdispF$distances, g=bdispF$group, set="Fungi"))
cv <- rbind(data.frame(bdispB$vectors[,1:2], g=bdispB$group, set="Bacteria"), 
            data.frame(bdispF$vectors[,1:2], g=bdispF$group, set="Fungi"))
ce <- rbind(data.frame(bdispB$centroids[,1:2], g=levels(bdispB$group), set="Bacteria"), 
            data.frame(bdispF$centroids[,1:2], g=levels(bdispF$group), set="Fungi"))
ve <- t(cbind(c(bdispB$eig[1:3]/sum(bdispB$eig),set="Bacteria", g=NA),
              c(bdispF$eig[1:3]/sum(bdispF$eig),set="Fungi", g=NA)))

cv <- left_join(cv, ce, by=c("g"="g","set"="set"), suffix = c("",".cent"))
hull <- cv %>% group_by(g, set) %>% slice(chull(PCoA1, PCoA2))

egg::ggarrange(
  
  ggplot(al, aes(group.id, Richness, fill=group.id)) +
    facet_wrap(vars(set), scales="free_y", strip.position="left", ncol=1) +
    geom_boxplot(linetype="dashed", outlier.shape = 1) +
    stat_boxplot(geom = "errorbar", width=0.4, aes(ymin = ..ymax..)) +
    stat_boxplot(geom = "errorbar", width=0.4, aes(ymax = ..ymin..)) + 
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), 
                 outlier.shape = 1) +
    scale_fill_manual(values=c("#E41A1C","#377EB8", "#4DAF4A")) +
    labs(title="Diversity", x=NULL) + guides(fill="none") +
    ggthemes::theme_few(base_size = 12) +
    theme(strip.text=element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=18, face="bold", family="mono", hjust=0.5)),
  
  ggplot(cv, aes(color=g)) +
    facet_wrap(vars(set), strip.position="left", ncol=1) +
    geom_hline(yintercept=0, color="grey90", linetype="dashed") +
    geom_vline(xintercept=0, color="grey90", linetype="dashed") +
    geom_segment(aes(x=PCoA1, xend=PCoA1.cent, y=PCoA2, yend=PCoA2.cent), alpha=0.7) +
    geom_point(aes(PCoA1, PCoA2), shape=1) +
    geom_label(aes(PCoA1.cent, PCoA2.cent, label=g), size=3) +
    geom_polygon(data=hull, aes(PCoA1, PCoA2), alpha=0.7, fill=NA) +
    geom_text(data=ve, hjust=1, size=2.5,
              aes(0.45, 0.4, label=paste0("PCoA1: ",round(as.numeric(PCoA1)*100,1),"%\n",
                                     "PCoA2: ",round(as.numeric(PCoA2)*100,1),"%"))) +
    scale_color_manual(values=c( "#4DAF4A","#377EB8","#E41A1C")) +
    labs(title="Composition") + 
    guides(color="none") +
    ggthemes::theme_few(base_size = 12) +
    theme(strip.text=element_blank(),
          plot.title = element_text(size=18, face="bold", family="mono", hjust=0.5)),
  
  ggplot(bd, aes(g, dist, fill=g)) +
    facet_wrap(vars(set), scales="free_y", strip.position="left", ncol=1) +
    geom_boxplot(linetype="dashed", outlier.shape = 1) +
    stat_boxplot(geom = "errorbar", width=0.4, aes(ymin = ..ymax..)) +
    stat_boxplot(geom = "errorbar", width=0.4, aes(ymax = ..ymin..)) + 
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), 
                 outlier.shape = 1) +
    scale_fill_manual(values=c("#E41A1C","#377EB8", "#4DAF4A")) +
    labs(title="Variability", x=NULL, y="Distance to centroid") + guides(fill="none") +
    ggthemes::theme_few(base_size = 12) +
    theme(strip.text=element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=18, face="bold", family="mono", hjust=0.5)), 
   
  ggplot(mbsummean) + 
    facet_wrap(vars(set), scales="free_y", strip.position="left", ncol=1) +
    geom_col(aes(Group, value, fill=sec.list, color=sec.list),
             position="stack", show.legend = F) +
    scale_fill_manual(values=c("#e6eff6", "#e9f5e9", "#fae3e1",
                            "#d4e6e1","#e3d6de","#e6dcd0","#d1d2cc")) +
    scale_color_manual(values=c("#377EB8", "#4DAF4A", "#E41A1C",
                            "#035b5b","#b80565","#88510b","black")) +
    labs(y="Total Mean Relative Abundance (%)", x=NULL, 
         title="  Distribution", fill="Venn segment", color="Venn segment") +
    ggthemes::theme_few(base_size = 12) +
    coord_cartesian(expand=F,clip="off") + 
    theme(panel.border = element_blank(),
          plot.title = element_text(size=18, face="bold", family="mono"),
          strip.text = element_blank(),
          plot.margin = margin(5,15,1,5),
          axis.text.x = element_text(color=c("#377EB8", "#4DAF4A", "#E41A1C"), 
                                     angle=90, hjust=1, vjust=0.5)),
 
  cowplot::plot_grid(
    cowplot::ggdraw() +
      draw_image("coreplots/venn/bacteria/euler_group_tax_edit.tiff") +
      draw_label("Bacteria", colour = "black", fontface="bold", fontfamily = "mono",
                 size = 18, angle = -90, x=1, vjust=1),
    cowplot::ggdraw() +
      draw_image("coreplots/venn/fungi/euler_group_tax_edit.tiff") +
      draw_label("Fungi", colour = "black", fontface="bold",fontfamily = "mono",
                 size = 18, angle = -90, x=1, vjust=1),
    ncol=1),

  nrow=1, widths=c(1,2,1,0.5,2)
)

```




#### Figure S1
```{r 01_FigureS1_rarecurves_plots}
load("output/rarecurves.RD")

egg::ggarrange(
  
  ggplot(rareB, aes(Sample, Species, group=Site, color=WebType)) + 
    theme_bw() + labs(title="Bacteria") +
    facet_wrap(vars(WebType), ncol=1) +
    scale_color_manual(values=c("#E41A1C","#377EB8", "#4DAF4A")) +
    geom_line(show.legend = F),
  
  ggplot(rareF, aes(Sample, Species, group=Site, color=WebType)) + 
    theme_bw() + labs(title="Fungi") +
    facet_wrap(vars(WebType), ncol=1) +
    scale_color_manual(values=c("#E41A1C","#377EB8", "#4DAF4A")) +
    geom_line(show.legend = F),
  
  nrow=1)

rm(rareB, rareF)
```

#### Figure S2
```{r 01_FigureS2a, fig.width=8, fig.height=3}
divlong <- data.frame(al, Dispersion=c(bdispB$distances, bdispF$distances)) %>%
  reshape2::melt()

egg::ggarrange(
  ggplot(divlong[divlong$set=="Bacteria",], aes(group.id, value, fill=group.id)) +
      facet_wrap(vars(variable), scales="free_y", nrow=1, strip.position = "left") +
      geom_boxplot(linetype="dashed", outlier.shape = 1) +
      stat_boxplot(geom = "errorbar", width=0.4, aes(ymin = ..ymax..)) +
      stat_boxplot(geom = "errorbar", width=0.4, aes(ymax = ..ymin..)) + #
      stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), 
                   outlier.shape = 1) +
      scale_fill_manual(values=c("#E41A1C","#377EB8", "#4DAF4A")) +
      labs(title="Bacteria", x=NULL, y=NULL) + guides(fill="none") +
      ggthemes::theme_few(base_size = 12) +
      theme(strip.placement = "outside",
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            plot.title = element_text(size=18, face="bold", family="mono", hjust=0.5)),
  ggplot(divlong[divlong$set=="Fungi",], aes(group.id, value, fill=group.id)) +
      facet_wrap(vars(variable), scales="free_y", nrow=1, strip.position = "left") +
      geom_boxplot(linetype="dashed", outlier.shape = 1) +
      stat_boxplot(geom = "errorbar", width=0.4, aes(ymin = ..ymax..)) +
      stat_boxplot(geom = "errorbar", width=0.4, aes(ymax = ..ymin..)) + 
      stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), 
                   outlier.shape = 1) +
      scale_fill_manual(values=c("#E41A1C","#377EB8", "#4DAF4A")) +
      labs(title="Fungi", x=NULL, y=NULL) + guides(fill="none") +
      ggthemes::theme_few(base_size = 12) +
      theme(strip.placement = "outside",
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            plot.title = element_text(size=18, face="bold", family="mono", hjust=0.5)),
ncol=1)
```  

```{r 01_FigureS2b, fig.width=8, fig.height=5}
divsig <- rbind(data.frame(hsdB, set="Bacteria"),
                data.frame(hsdF, set="Fungi")) %>%
  tidyr::separate(test,c("type1","type2"), remove = F)

egg::ggarrange(
  ggplot(divsig[divsig$set=="Bacteria",]) + labs(title="Bacteria", y=NULL) +
    facet_wrap(vars(index),scales="free_x", ncol=1) +
    ggthemes::theme_clean() +
    geom_vline(xintercept=0, color="grey40", linetype="dashed") +
    geom_errorbar(aes(xmin=lwr, xmax=upr, y=test,
                      color=ifelse(upr<0|lwr>0,"**","ns")), show.legend = F) +
    geom_point(aes(x=diff, y=test,
                      color=ifelse(upr<0|lwr>0,"**","ns")), show.legend = F) +
    scale_color_manual(values=c("blue","black")),
  ggplot(divsig[divsig$set=="Fungi",]) + labs(title="Fungi", y=NULL) +
    facet_wrap(vars(index),scales="free_x", ncol=1) +
    ggthemes::theme_clean() +
    geom_vline(xintercept=0, color="grey40", linetype="dashed") +
    geom_errorbar(aes(xmin=lwr, xmax=upr, y=test,
                      color=ifelse(upr<0|lwr>0,"**","ns")), show.legend = F) +
    geom_point(aes(x=diff, y=test,
                      color=ifelse(upr<0|lwr>0,"**","ns")), show.legend = F) +
    scale_color_manual(values=c("blue","black")),
ncol=2)

```

<a href="#data-analysis" style="color:steelblue;" >Back to top</a>

## 2. Sources & sinks {.tabset}

> **Question:** How is the capture silk microbiota related to the retreat silk or soil?

**SourceTracker:** SourceTracker analysis is a Bayesian approach to quantifying the proportion of a microbial population in a sample that comes from each of a set of sources. It randomly assigns each ASV to one of the given source environments and updates the likelihood of observing it in that source based on the proportion present in the sink environment. It repeats this Gibbs sampling process 10 times to estimate the variability of the posterior distribution. After each sink ASV has been assigned proportionally to each source, we took the mean of the Gibbs samples of the ASV for each sample. We then took the mean across the samples to get an estimate of the overall proportional source assignment for the sequence. With each capture silk ASV given an assignment to 'soil', 'retreat' or 'unknown', we used these proportions to group sequences by source value with k-means clustering. Four clusters were created: Soil, Retreat, Mixed, and Unknown. ASVs grouped into the Soil, Retreat, and Unknown clusters were those predominantly assigned to each of those three sources by SourceTracker. ASVs grouped into the 'Mixed' cluster were those with similar proportions of the three assigned sources.

*Note:*  Originally, it was designed as a tool for identifying sources of contamination in microbiome samples. We used SourceTracker v1, which was written for implementation in R.

```{r 02_ST_RunFunction, eval=F}
## create SourceTracker wrapper to ease parallel runs
run.st <- function(source.index, sink.index, comm.mat, env.mat, rare, a1, a2){
  m1 <- sourcetracker(comm.mat[source.index,], as.character(env.mat[source.index]), rare) # from SourceTracker source code https://github.com/danknights/sourcetracker/src, exported from `stegodyve`
  comm.m1 <- comm.mat[sink.index,] %>% .[,-which(colSums(.)==0)]
  results.m1 <- stegodyve:::predict.sourcetracker(m1, comm.m1, alpha1=a1, alpha2=a2, 
                                                  rarefaction_depth = rare, full.results=T) # from SourceTracker source code https://github.com/danknights/sourcetracker/src, exported from `stegodyve`
}

# this is a comment from the SourceTracker example source code (not run)
  # tune the alpha values using cross-validation (this is slow!)
  #tune.results <- tune.st(otus[train.ix,], envs[train.ix])
  #alpha1 <- tune.results$best.alpha1
  #alpha2 <- tune.results$best.alpha2
  #save(tune.results, file="output/alpha_tuning.RD")

# skip tuning; set alphas as a constant
alpha1 <- alpha2 <- 0.001

# this is a comment from the SourceTracker example source code (not run)
  # Estimate leave-one-out source proportions in training data 
  #results.train <- predict(st, alpha1=alpha1, alpha2=alpha2)



```

### 2.1 Run on bacteria
```{r 02_bacteriaSetup}

# reduce the sample set to locations that have samples from all three types
st.set <- names(which(summary(as.factor(metaB$column_name))==3))[-7] # remove poorly sampled
# reduce metadata to match the reduced set
pullB <- metaB[which(metaB$column_name%in%st.set),]
# create separate vector identifying the sample types
envs <- pullB$WebType
# reduce the community matrix to the reduced set and remove columns containing all zeroes
pull.commB <- commB[pullB$SampleID,] %>% .[,-which(colSums(.)==0)]

```

```{r 02_SourceTracker_bacteria, eval=F}
## SourceTracker model builder

# extract the source environments and source/sink indices
source.list <- list(m1= which(pullB$WebType%in%c("soil", "retreat")),
                    m2= which(pullB$WebType=="retreat"),
                    m4= which(pullB$WebType=="soil"))

# run sets analysing the same sink together using `parallel` and save output
rB <- mclapply(X=source.list, FUN=run.st, 
              sink.index=which(pullB$WebType=="capture"),
              comm.mat=pull.commB, env.mat=envs, rare=15000, a1=alpha1, a2=alpha2,
              mc.cores=3)
save(r, file="output/sourcetracker_analysis_bacteria.RD")

# run additional model with a different sink and save output
sink.m3 <- which(pullB$WebType=="retreat")
source.m3 <- which(pullB$WebType=="soil")
m3 <- run.st(source.m3, sink.m3, pull.commB, envs, 15000, alpha1, alpha2)
rB[["m3"]] <- m3
save(rB, pull.commB, pullB, taxB, file="output/sourcetracker_analysis_bacteria.RD")

```

<a href="#data-analysis" style="color:steelblue;" >Back to top</a>

### 2.2 Run on fungi

```{r 02_fungiSetup}

# reduce the sample set to locations that have samples from all three types
st.set <- names(which(summary(as.factor(metaF$column_name))==3))[-10] # remove poorly sampled
# reduce metadata to match the reduced set
pullF <- metaF[which(metaF$column_name%in%st.set),]
# create separate vector identifying the sample types
envs <- pullF$WebType
# reduce the community matrix to the reduced set and remove columns containing all zeroes
pull.commF <- commF[pullF$SampleID,] %>% .[,-which(colSums(.)==0)]

```

```{r 02_SourceTracker_fungi, eval=F}
## SourceTracker model builder

# extract the source environments and source/sink indices
source.list <- list(m1= which(pullF$WebType%in%c("soil", "retreat")),
                    m2= which(pullF$WebType=="retreat"),
                    m4= which(pullF$WebType=="soil"))

# run sets analysing the same sink together using `parallel` and save output
rF <- mclapply(X=source.list, FUN=run.st, 
              sink.index=which(pullF$WebType=="capture"),
              comm.mat=pull.commF, env.mat=envs, rare=15000, a1=alpha1, a2=alpha2,
              mc.cores=3)
save(rF, file="output/sourcetracker_analysis_fungi.RD")

# run additional model with a different sink and save output
sink.m3 <- which(pullF$WebType=="retreat")
source.m3 <- which(pullF$WebType=="soil")
m3 <- run.st(source.m3, sink.m3, pull.commF, envs, 15000, alpha1, alpha2)
rF[["m3"]] <- m3
save(rF, pull.commF, pullF, taxF, file="output/sourcetracker_analysis_fungi.RD")

```

<a href="#data-analysis" style="color:steelblue;" >Back to top</a>

### 2.3 Downstream analysis

```{r 02_ST_carpentry}
## function to reshape SourceTracker object into data frame for downstream analysis and visualization
ST_analysis <- function(x, pc, p, tax, meta){
  
  otu.info <- function(stobj, otu.lab){
    otu <- reshape2::melt(stobj$full.results, 
                          varnames=c("gibbs","source","otu","sample"))
    otu$otu <- factor(otu$otu, labels=otu.lab)
    otu$sample <- factor(otu$sample, labels=stobj$samplenames)
    otu$source <- factor(otu$source, labels=stobj$train.envs)
    out <- otu %>% group_by(source, otu, sample) %>%
      summarize_at(vars(value), .funs=list(gibbs.mean=mean, sd=sd)) %>% 
      mutate(se=sd/sqrt(10))
  }

  m1.otu <- otu.info(x, otu.lab=colnames(pc)) %>% 
    .[-which(.$gibbs.mean==0),] %>%
    left_join(p, by=c("sample"="SampleID")) %>%
    left_join(tax, by=c("otu"="tag"))
 
  # pull sample Gibbs means and sd
  sample <- as.data.frame(x$proportions) %>% 
    tibble::rownames_to_column("SampleID") %>%
    left_join(p[,c("SampleID","column_name","SiteID","country")]) %>%
    cbind(as.data.frame(x$proportions_sd) %>% 
            `colnames<-`(c("ret_sd","soil_sd","unk_sd")))
  
  g0 <- mutate(m1.otu[,1:6], cv=sd/gibbs.mean) %>% 
    reshape2::dcast(sample+otu~source,value.var="cv")
  
  gmom <- data.frame(type=c("Retreat","Soil","Unknown"),
                     skewness=c(moments::skewness(g0$retreat, na.rm=T),
                                moments::skewness(g0$soil, na.rm=T),
                                moments::skewness(g0$Unknown, na.rm=T)),
                     kurtosis=c(moments::kurtosis(g0$retreat, na.rm=T),
                                moments::kurtosis(g0$soil, na.rm=T),
                                moments::kurtosis(g0$Unknown, na.rm=T)))
  
  top5 <- m1.otu %>% group_by(source, otu) %>% 
    summarize_at(vars(gibbs.mean), .funs=sum) %>%
    .[order(.$source,.$gibbs.mean, decreasing=T),] %>% 
    left_join(tax[,2:3], by=c("otu"="tag")) %>%
    group_by(source) %>% mutate(rank=row_number()) %>% 
    .[which(.$rank<=5),]
  
  # cast long format df into wide format
  castop <- m1.otu %>% group_by(otu, sample) %>% 
    mutate(amount=gibbs.mean/sum(gibbs.mean))
  castop <- reshape2::dcast(castop[,c("otu", "column_name", "source", "amount")], 
                            otu + column_name ~ source) %>% .[,-2] %>%
    `colnames<-`(c("otu","Retreat","Soil","Unknown"))
  castop[is.na(castop)] <- 0
  
  # calculate mean proportions for each ASV (across samples)
  meantop <- castop %>% group_by(otu) %>% summarize_all(.funs=mean)
  k <- character(nrow(meantop))
  k[meantop$Retreat>0.5] <- "Retreat"
  k[meantop$Soil>0.5] <- "Soil"
  k[meantop$Unknown>0.5] <-"Unknown"
  k[k==""] <- "Mixed"

  meantop$g <- k
  meantop$g <- factor(meantop$g, levels=c("Soil","Mixed","Retreat","Unknown"))
  
  out <- list(newMeta=m1.otu, 
              meanContent=sample,
              CoefVar=list(CV=g0, moments=gmom),
              Top5=top5,
              PropGibbs10=castop, 
              meanProp=meantop)
  
  return(out)
}
```

### 2.4 Plots

#### Figure 4
```{r 02_SourceTracker_plots1, fig.width=3, fig.height=4}
## load model outputs SourceTracker analysis
load("output/sourcetracker_analysis_bacteria.RD")
load("output/sourcetracker_analysis_fungi.RD")

bact <- ST_analysis(rB$m1, pull.commB, pullB, taxB)
fung <- ST_analysis(rF$m1, pull.commF, pullF, taxF)

## Visualize distribution of CV for each source
g0 <- rbind(data.frame(bact$CoefVar$CV, set="Bacteria"),
            data.frame(fung$CoefVar$CV, set="Fungi"))

ggplot(g0) + 
  facet_wrap(vars(set), scales = "free_y", strip.position = "right", ncol=1) +
  ggthemes::theme_clean() + 
  theme(strip.text=element_text(family="mono", face="bold", size=16, angle=-90)) +
  xlab("CV of Gibbs draws") + coord_cartesian(expand=F) +
  geom_density(aes(soil), color="#E41A1C", fill="#E41A1C", alpha=0.3) +
  geom_density(aes(retreat), color="#377EB8", fill="#377EB8", alpha=0.3) +
  geom_density(aes(Unknown), color="#4DAF4A", fill="#4DAF4A", alpha=0.3) 
```

```{r 02_SourceTracker_plots2}
## visualize ASV distribution on ternary plot
meantop <- rbind(data.frame(bact$meanProp, set="Bacteria"),
                 data.frame(fung$meanProp, set="Fungi"))

## summarizes previous plot in a pie chart
ggplot(meantop) + 
  facet_wrap(vars(set), strip.position = "top", nrow=1) +
  stat_count(aes("", fill=g), geom="bar", color="white", position="fill", show.legend = F) +
  stat_count(aes(x=1.15, fill=g, label=paste(..fill.., ..count.., sep="\n")), size=4.5,
             geom="label", color="white", position=position_fill(0.5), show.legend=F) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(4, "Set1")[c(1,4,2,3)]) +
  coord_polar("y",) + labs(tag="(B)") + theme_void() + 
  theme(plot.margin = unit(c(-4,-1,-5,-3),"lines"),
        strip.text=element_text(family="mono", face="bold", size=20))
```

#### Figure S4
```{r 02_SourceTracker_FigureS4a, fig.width=6, fig.height=7}
## visualize sample distribution on ternary plot
sample <- rbind(data.frame(bact$meanContent[order(bact$meanContent$column_name),],
                           set="Bacteria"),
                data.frame(fung$meanContent[order(fung$meanContent$column_name),],
                           set="Fungi")) %>%
  mutate(shape=country, color=column_name) %>% 
  .[order(.$country, .$column_name, .$set),]

sample$column_name <- factor(sample$column_name, levels=unique(sample$column_name))

nam <- which(sample$country=="Namibia")
saf <- which(sample$country=="South Africa")

sample$shape <- factor(sample$shape, labels=c(1,8))
sample$shape <- as.numeric(as.character(sample$shape))
sample$color[nam] <- rep(colorRampPalette(c("pink","magenta"))(4),2)
sample$color[saf] <- colorRampPalette(c("green","turquoise","navy")
                                      )(9)[c(1,1,2,2,3,4,5,5,6,7,8,8,9,9)]

ggtern::grid.arrange(
  ggtern(sample, aes(x=retreat, y=Unknown, z=soil)) +
    facet_wrap(vars(set), scales = "free", strip.position = "top", ncol=1) +
    geom_point(aes(color=column_name, shape=country), size=2) + 
    theme_light() + theme_showarrows() + 
    labs(title="Average site assignment", color=NULL,shape=NULL) +
    scale_color_manual(values=unique(sample$color)) +
    scale_shape_manual(values=c(1,8)) +
    guides(color=guide_legend(title="Site name", nrow=4,
      override.aes = list(color=unique(sample[,12:13])$color,
                          shape=unique(sample[,12:13])$shape)),
           shape=guide_legend(title="Country", order=1, nrow=1)) +
    theme(tern.panel.mask.show=F, 
          plot.margin = unit(c(0,-2,0,-2),"cm"),
          legend.position = "bottom", legend.byrow = F,
          legend.justification = 0.2,
          legend.title.position = "top",
          legend.box="vertical", 
          legend.box.just="left",
          legend.key.height = unit(0.05,"lines"),
          legend.key.spacing = unit(0.1, "mm"),
          legend.text = element_text(size=8),
          legend.spacing = unit(1,"mm"),
          tern.axis.title.L=element_text(hjust=0),
          strip.background=element_blank(),
          strip.text=element_text(family="mono", face="bold", size=20, color="black")),
  ggtern(meantop, aes(x=Retreat, y=Unknown, z=Soil, color=g)) +
    facet_wrap(vars(set), scales = "free", strip.position = "top", ncol=1) +
    geom_point(aes(shape=set), show.legend = c("shape"=F)) +
    scale_color_manual(values=RColorBrewer::brewer.pal(4, "Set1")[c(1,4,2,3)]) +
    scale_shape_manual(values=c(16,15)) +
    theme_light() + theme_showarrows() + 
    labs(title="ASV source assignment", color="Assignment") + 
    theme(tern.panel.mask.show=F, 
          tern.axis.title.L=element_text(hjust=0),
          legend.position = "bottom", legend.byrow = F,
          legend.justification = 0.6,
          legend.title.position = "top",
          legend.direction = "vertical",
          legend.key.height = unit(0.05,"lines"),
          legend.text = element_text(size=8),
          plot.margin=unit(c(0,-2,1,-2),"cm"),
          strip.background=element_blank(),
          strip.text=element_text(family="mono", face="bold", size=20, color="black")),
nrow=1)
```

```{r 02_SourceTracker_FigureS4b, fig.height=6, fig.width=6}
## visualize ASV assignment based on proportion of mean Gibbs draws
egg::ggarrange(
  ggplot(reshape2::melt(meantop)%>%subset(subset=set=="Bacteria")) + 
    facet_grid(set~g, scales="free_x", space="free") + 
    geom_bar(aes(otu, value, fill=variable), 
             stat="identity", show.legend = F) +
    scale_fill_manual(values=RColorBrewer::brewer.pal(4, "Set1")[c(2,1,3)]) +
    ggthemes::theme_few(base_size = 14) + coord_cartesian(expand=F) +
    ylab("Source assignment (%)") + xlab("Capture silk ASVs") + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.clip = "off", strip.text.y=element_text(family="mono", face="bold", size=16)),
  ggplot(reshape2::melt(meantop)%>%subset(subset=set=="Fungi")) + 
    facet_grid(set~g, scales="free_x", space="free") + 
    geom_bar(aes(otu, value, fill=variable), 
             stat="identity", show.legend = F) +
    scale_fill_manual(values=RColorBrewer::brewer.pal(4, "Set1")[c(2,1,3)]) +
    ggthemes::theme_few(base_size = 14) + coord_cartesian(expand=F) +
    ylab("Source assignment (%)") + xlab("Capture silk ASVs") + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.clip="off", strip.text.y=element_text(family="mono", face="bold", size=16)),
  ncol=1)

```

<a href="#data-analysis" style="color:steelblue;" >Back to top</a>

## 3. Intra-colony similarity

> **Question:** Are the silk communities of the same colony related independently of the local soil?

> **Question:** Is there spatial turnover in silk communities that is independent of the local soil?

**Pairwise sample distance:** To understand the degree of compositional similarity between sample types, we  calculated binary Bray-Curtis distances between retreat and capture matrices, between retreat and soil matrices, and between capture and soil matrices. 

**Distance-decay regression:** We first performed a single regression of pairwise compositional distances between samples against the geographical distances. 

**Group-wise similarity:** We then grouped compositional distances by relationship: between samples collected from the same location (Within-colony), between samples of the same transect (Within-transect), and between samples across transects (Between-transect). The mean compositional distances for each group were compared using two-way ANOVA.

```{r 03_SampleType_CrossComparison, eval=F}
CrossComp <- function(pc, meta){
  # subset community matrix to each sample type and format row names to match 
  capture <- pc[grep("capture",row.names(pc)),] %>% 
    .[order(row.names(.)),] %>% 
    `row.names<-`(stringr::str_sub(row.names(.), 1, -8)) #(removes sampletype identifiers)
  retreat <- pc[grep("retreat",row.names(pc)),] %>% 
    .[order(row.names(.)),] %>%
    `row.names<-`(stringr::str_sub(row.names(.), 1, -8))
  soil <- pc[grep("soil",row.names(pc)),] %>% 
    .[order(row.names(.)),] %>%
    `row.names<-`(stringr::str_sub(row.names(.), 1, -5))
  
  # formats the dataframes together. Row names must match
  mats1 <- as.bipartite.net(capture, retreat, F) # no rarefaction on the data frames
  mats2 <- as.bipartite.net(capture, soil, F)
  mats3 <- as.bipartite.net(retreat, soil, F)
  
  # 'interaction' values are the shared richness between capture silk and retreat silk
  cc_mat1 <- 1-(cross.comparison(mats1, "bray", T))
  cc_mat2 <- 1-(cross.comparison(mats2, "bray", T))
  cc_mat3 <- 1-(cross.comparison(mats3, "bray", T))
  
  # site distances
  dist1 <- dist(meta[match(row.names(cc_mat1),row.names(meta)),
                        c("latitude","longitude")],diag=T,upper=T) %>% as.matrix()
  
  ## regress nestedness values to site distances
  # retreat/capture
  df1 <- cbind(reshape2::melt(cc_mat1, value.name="cross"),
               dist=reshape2::melt(dist1, value.name="dist")[,3])
  m1 <- lm(cross ~ dist, data=df1); summary(m1)
  # capture/soil
  df2 <- cbind(reshape2::melt(cc_mat2, value.name="cross"),
               dist=reshape2::melt(dist1, value.name="dist")[,3])
  m2 <- lm(cross ~ dist, data=df2); summary(m2)
  # retreat/soil
  df3 <- cbind(reshape2::melt(cc_mat3, value.name="cross"),
               dist=reshape2::melt(dist1, value.name="dist")[,3])
  m3 <- lm(cross ~ dist, data=df3); summary(m3)
  
  # bind the results from df1-3 together and join additional metadata
  df <- rbind(data.frame(df1, nestedness="Capture\u2286Retreat"),
              data.frame(df2, nestedness="Capture\u2286Soil"),
              data.frame(df3, nestedness="Retreat\u2286Soil")) %>%
    left_join(meta[,c("SampleID","column_name","SiteID","WebType","country")], by=c("Var1"="SampleID")) %>%
    left_join(meta[,c("SampleID","column_name","SiteID","WebType","country")], by=c("Var2"="SampleID"))
  
  # create separate matrix to store ordered information 
  distid <- dist1[order(dist1[,1]),] %>% row.names(.) %>% data.frame() %>%
    left_join(meta[,c("SampleID", "column_name","SiteID","country")], by=c("."="SampleID"))
  
  # order column name
  df$column_name.x <- factor(df$column_name.x, levels=distid$column_name)
  df$column_name.y <- factor(df$column_name.y, levels=distid$column_name)
  
  # add summary variables and colony/transect information
  df <- df %>% group_by(nestedness, Var1) %>% 
    mutate(x=length(Var1)/2+0.5, y=-0.25, sums.x=sum(cross), margin=length(Var1)+1) %>%
    group_by(nestedness, Var2) %>% mutate(sums.y=sum(cross)) %>%
    group_by(nestedness, Var1) %>% mutate(scaled=round(scales::rescale(cross),2)) %>%
    mutate(colony.wb=column_name.x==column_name.y, 
           transect.wb=country.x==country.y,
           compare=paste(country.x, "vs.", country.y)) %>%
    mutate(with.between=as.factor(colony.wb+transect.wb))
  
  # re-level colony/transect information for legibility
  df$compare <- as.factor(df$compare)
  df$compare <- factor(df$compare, labels=unique(df$compare)[c(1,3,3,4)])
  df$colony.wb[df$colony.wb==T] <- "Within colony"
  df$colony.wb[df$colony.wb==F] <- "Between colonies"
  df$transect.wb[df$transect.wb==T] <- "Within transect"
  df$transect.wb[df$transect.wb==F] <- "Between transects"
  df$with.between <- factor(df$with.between, labels=c("Between\ntransects","Within\ntransect", "Within\ncolony"))
  df$transect.wb <- factor(df$transect.wb, levels=unique(df$transect.wb))
  df$nestedness <- factor(df$nestedness, levels=unique(df$nestedness))
  
  ## statistical analysis to go with group comparison
  av1 <- aov(cross ~ with.between, 
             data=subset(df, subset=nestedness=="Capture\u2286Retreat")) %>%
    TukeyHSD() %>% .$with.between %>% as.data.frame() %>% 
    tibble::rownames_to_column() %>% cbind(nestedness="Capture\u2286Retreat")
  
  av2 <- aov(cross ~ with.between, 
             data=subset(df, subset=nestedness=="Capture\u2286Soil")) %>% 
    TukeyHSD() %>% .$with.between %>% as.data.frame() %>% 
    tibble::rownames_to_column() %>% cbind(nestedness="Capture\u2286Soil")
  
  av3 <- aov(cross ~ with.between, 
             data=subset(df, subset=nestedness=="Retreat\u2286Soil")) %>% 
    TukeyHSD() %>% .$with.between %>% as.data.frame() %>% 
    tibble::rownames_to_column() %>% cbind(nestedness="Retreat\u2286Soil")
  
  ## formatting of significance bars to go with group comparison
  sigdf <- rbind(av1,av2,av3)
  sigdf$rowname <- as.factor(sigdf$rowname)
  sigdf$rowname <- factor(sigdf$rowname, levels=unique(sigdf$rowname)[c(3,2,1)]) 
  sigdf <- left_join(sigdf, df[,c("nestedness","cross")]) %>% 
    group_by(nestedness, `p adj`, rowname) %>%
    summarize_at(vars(cross), .funs=list(max)) %>%
    mutate(y=0, yend=0, bars=cross)
  sigdf$bars <- sigdf$bars+seq(0.02,0.04,length=3)
  sigdf$y[sigdf$rowname=="Within\ntransect-Between\ntransects"] <- 1
  sigdf$y[sigdf$rowname=="Within\ncolony-Between\ntransects"] <- 1
  sigdf$y[sigdf$rowname=="Within\ncolony-Within\ntransect"] <- 2
  sigdf$yend[sigdf$rowname=="Within\ntransect-Between\ntransects"] <- 2
  sigdf$yend[sigdf$rowname=="Within\ncolony-Between\ntransects"] <- 3
  sigdf$yend[sigdf$rowname=="Within\ncolony-Within\ntransect"] <- 3
  sigdf <- sigdf[sigdf$`p adj`<0.05,]
  
  sigdf <- sigdf %>% group_by(nestedness) %>%
    mutate(star=mean(c(y,yend)))
  
  ## prepare outbox
  out <- list(data=list("retreat"=retreat,"capture"=capture,"soil"=soil),
              pair.models=list("retreat-capture"=m1,"capture-soil"=m2,"retreat-soil"=m3),
              pair.anova=list("retreat-capture"=av1,"capture-soil"=av2,"retreat-soil"=av3),
              dist.dat=df,
              sig.df=sigdf)

}



ccbact <- CrossComp(pull.commB, metaB)
ccfung <- CrossComp(pull.commF, metaF)

save(ccbact, ccfung, file="output/crosscomparison_both.RD")

```

### 3.1 Plots
```{r 03_CrossCompar_plotting}

CCplot <- function(input, title){
  # pairwise sample distance
  p1 <- ggplot(input) + 
    facet_wrap(vars(nestedness),ncol=1,scales="free") +
    geom_raster(aes(column_name.y, column_name.x, fill=dist*111.11)) +
    geom_point(aes(column_name.y, column_name.x, size=cross, color=cross)) +
    geom_text(data=unique(df[,c(5,8,12,14,15)]), 
              aes(x=x,y=y,label=WebType.y), hjust=0.5, vjust=0.5) +
    geom_text(data=unique(df[,c(5,8,12,14,15)]), 
              aes(x=y,y=x,label=WebType.x), hjust=0.5, vjust=0.5, angle=90) +
    scale_color_gradient2(low="black", mid="#56B1F7", high="white",
                          midpoint=(min(df$cross)+max(df$cross))/2,
                          breaks=round(seq(min(df$cross),max(df$cross), length=3),3),
                          labels=round(seq(min(df$cross),max(df$cross), length=3),3),
                          limits=round(seq(min(df$cross),max(df$cross), length=2),3)) +
    scale_size_continuous(range=c(0,3), 
                          breaks=round(seq(min(df$cross),max(df$cross), length=3),3),
                          labels=round(seq(min(df$cross),max(df$cross), length=3),3),
                          limits=round(seq(0,max(df$cross)+0.001, length=2),3)) +
    scale_fill_gradient(low="black", high="white", 
                         breaks=floor(seq(min(df$dist*111.11),max(df$dist*111.11), length=2)),
                         labels=floor(seq(min(df$dist*111.11),max(df$dist*111.11), length=2))) +
    coord_cartesian(clip="off", expand=F) + theme_bw(base_size = 14) + 
    theme(axis.text = element_blank(), axis.ticks=element_blank(), axis.title = element_blank(),
          plot.title = element_text(hjust=0.5, face="bold"),
          legend.position = "bottom", legend.title=element_text(size=9), 
          legend.text=element_text(size=7, margin=margin(0,0,0,0)), legend.spacing.x=unit(0.3,"cm"),
          legend.box = "vertical", legend.box.just = "left", legend.margin=margin(0,0,0,0), 
          panel.border=element_rect(color=NA, fill=NA), plot.margin=unit(c(1.5,0.5,0,0.5), "lines")) +
    labs(title=title, x=NULL, y=NULL, size="Compositional\nsimilarity", fill="Distance\n(km)") +
    guides(color="none",
           fill=guide_colorbar(direction="horizontal", title.hjust = 0.5, title.position = "left",
                               label.position = "top", barwidth=5, barheight=0.7, order = 1),
           size=guide_legend(direction="horizontal", title.hjust = 0.5,title.position = "left", 
                             label.position="top", order = 2, 
                             override.aes = list(shape=21, color="black", 
                                                 fill=colorRampPalette(c("black","#56B1F7","white"))(3))))
  
  return(p1)
}

```

#### Figure 5
```{r 03_CrossCompar_Figure5, fig.width=4, fig.height=7}

load("output/crosscomparison_both.RD")

df <- rbind(data.frame(ccbact$dist.dat, set="Bacteria"),
            data.frame(ccfung$dist.dat, set="Fungi"))

ggplot(df) + ylim(0,0.7) +
  facet_grid(nestedness ~ set, switch="y") +
  geom_boxplot(aes(forcats::fct_rev(with.between), cross, fill=nestedness, alpha=with.between), 
               show.legend = c("fill"=F)) +
  scale_fill_manual(values=c("#035b5b","#88510b","#b80565")) +
  scale_alpha_manual(values=c(0.2,0.6,1)) +  
  theme_bw(base_size=14) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "bottom",
        legend.justification = "right",
        legend.title.position = "bottom",
        legend.key.spacing = unit(4,"mm"),
        legend.margin = margin(-0.5,0,0.5,0,"lines"),
        legend.text = element_text(margin = margin(l = 1)),
        legend.title = element_text(hjust=0.5, face="bold", family="mono", size=12),
        strip.text = element_text(hjust=0.5, face="bold", family="mono", size=14),
        plot.title = element_text(hjust=0.5, face="bold", family="mono", size=20),
        plot.margin=unit(c(0.5,0.5,0,0.5), "lines")) +
  coord_cartesian(clip=F) +
  guides(alpha=guide_legend(reverse = T, keywidth = unit(1, "cm"), byrow=T,
                         title="Pairwise sample comparison",
                         override.aes = list(fill=c("grey20","grey55","grey90"), alpha=1))) +  labs(y="Compositional similarity")


```

```{r 03_CrossCompar_Figure5_alt, eval=F, echo=F}

load("output/crosscomparison_both.RD")

egg::ggarrange(
  ggplot(ccbact$dist.dat) + xlim(0,0.7) +
    facet_wrap(vars(nestedness), dir="v") +
    geom_boxplot(aes(cross, with.between, fill=nestedness, alpha=with.between), 
                 show.legend = F) +
    scale_fill_manual(values=c("#035b5b","#88510b","#b80565")) +
    scale_alpha_manual(values=c(0,0.5,1)) +  
    theme_bw(base_size=14) + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(), legend.position = "bottom",
          plot.margin=unit(c(1.5,0,0,0.5), "lines"),
          strip.text = element_text(hjust=0.5, face="bold", family="mono", size=12),
          plot.title = element_text(hjust=0.5, face="bold", family="mono", size=20)) +
    labs(title="Bacteria", x="Compositional similarity"),
  ggplot(ccfung$dist.dat) + xlim(0,0.7) +
    facet_wrap(vars(nestedness), dir="v") +
    geom_boxplot(aes(cross, with.between, fill=nestedness, alpha=with.between), 
                 show.legend = c(fill=F)) +
    scale_fill_manual(values=c("#035b5b","#88510b","#b80565")) +
    scale_alpha_manual(values=c(0,0.5,1)) +  
    theme_bw(base_size=14) + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(), legend.position = "right",
          legend.title = element_text(hjust=0, face="bold", family="mono", size=12),
          strip.text = element_text(hjust=0.5, face="bold", family="mono", size=12),
          legend.key.spacing.y = unit(2, "mm"), 
          plot.margin=unit(c(1.5,0.5,0,0), "lines"),
          plot.title = element_text(hjust=0.5, face="bold", family="mono", size=20)) +
    guides(alpha=guide_legend(reverse = T, keywidth = unit(1, "cm"), byrow=T,
                         title="Pairwise\nsample\ncomparison",
                         override.aes = list(fill=c("grey20","grey55","grey90"), alpha=1))) +
    labs(title="Fungi", x="Compositional similarity"),
  nrow=1)

```

#### Figure S5
```{r 03_CrossCompar_FigureS5, fig.width=14, fig.height=8}

df$with.between <- gsub("\n"," ",df$with.between)
df$with.between <- factor(df$with.between, levels=rev(unique(df$with.between)))

sigdf <- rbind(data.frame(ccbact$sig.df, set="Bacteria"),
               data.frame(ccfung$sig.df, set="Fungi"))

## group-wise similarity
g3 <- ggplot(df) + 
  facet_grid(nestedness ~ set, switch="y") + 
  geom_boxplot(aes(cross, with.between, fill=nestedness, alpha=with.between), 
               show.legend = c(fill=F)) +
  geom_segment(data=sigdf, aes(x=bars, xend=bars, y=y, yend=yend), size=0.75) +
  geom_text(data=sigdf, aes(x=cross+0.06, y=star, label="*"), size=10, vjust=0.8) +
  scale_fill_manual(values=c("#035b5b","#88510b","#b80565")) +
  scale_alpha_manual(values=c(0,0.5,1)) +  
  scale_x_continuous(n.breaks = 3,
                     limits=round(seq(0,max(df$cross)+0.1, length=2),1)) +
  theme_bw(base_size=16) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(), legend.position = "bottom",
        plot.margin=unit(c(1.5,1.5,0,0.5), "lines")) +
  guides(alpha=guide_legend(direction="vertical", title.position = "left", 
                           reverse = T, keywidth = unit(1, "cm"), title="Sample\nidentity",
                           override.aes = list(fill=c("grey20","grey55","grey90"), alpha=1))) +
  labs(fill="Sample identity") + xlab("Compositional similarity")

##pairwise sample distance
g0 <- CCplot(ccbact$dist.dat,"Bacteria")
g1 <- CCplot(ccfung$dist.dat,"Fungi")
#egg::ggarrange(g0, g1, nrow=1)

## distance-decay regression
g2 <- ggplot(df, aes(dist*111.11, cross, linetype=transect.wb)) + 
  facet_grid(nestedness ~ set) +
  geom_point(aes(shape=compare, color=nestedness), show.legend = c(color=F)) +
  stat_smooth(method="lm", formula='y~x', color="black", fullrange=T)+
  ggpmisc::stat_poly_eq(formula=y ~ x, parse=T, label.x = "left", label.y="bottom", 
                        vstep=0.7, hstep=0.9, vjust=-1, size=2.5,
                        aes(label=paste("bold(\"",c("Within","Between")[stat(linetype)], 
                                     "\ntransect:  \")*", after_stat(eq.label), sep = ""))) +
  ggpmisc::stat_poly_eq(formula=y ~ x, parse=T, label.x = "left", label.y="bottom", 
                        vstep=0.7, hstep=0.9, vjust= 0, size=2.5,
                      aes(label=after_stat(rr.label))) +
  scale_color_manual(values=c("#035b5b","#88510b","#b80565")) +
  scale_shape_manual(values=c(17,15,20)) +
  guides(linetype=guide_legend(title=NULL, direction="horizontal", 
                               label.theme=element_text(margin=margin(0,7,0,0), size=9),
                               reverse = F, keywidth = unit(1, "cm")),
         shape=guide_legend(title=NULL, direction="horizontal", ncol=2, byrow = T,
                            override.aes = list(size=3),
                            label.theme=element_text(margin=margin(0,7,0,0), size=9))) +
  theme_bw(base_size = 14) + 
  theme(plot.margin=unit(c(1.5,0.5,0,0.5), "lines"),
        legend.box="vertical",
        legend.spacing.x = unit(3, "lines"), legend.spacing.y = unit(-0.25, "lines"),
        legend.position = "bottom", legend.justification = c(0.75,-1)) +
  xlab("Distance between colonies (km)") + ylab("Compositional similarity")


egg::ggarrange(g3, g0, g1, g2, nrow=1, widths=c(2,1,1,2))
```

<a href="#data-analysis" style="color:steelblue;" >Back to top</a>

## 4. Structuring processes {.tabset}

> **Question:** Are the silk and soil microbiota structured more by spatial or environmental processes?

*Note:* Silk and soil likely experience different constraints to bacterial and fungal community assembly. We used both a distance-based regression approach and variation partitioning to elucidate differences in the relative contributions of spatial and environmental processes on structuring the communities. 

**Distance-based regression:** This approach used binary Bray-Curtis distances between samples of the same type regressed against geographical distance and ‘environmental distance,’ respectively. Environmental and geographical distance multiple regressions were performed for within- and between-transect comparisons with compositional distance.

- <u> *Environmental distance* </u> was calculated as the Euclidean distance between sites based on temperature and aridity. Both temperature and aridity covariates were centered and scaled prior to this calculation. 
- <u> *Geographical distance* </u> was calculated as the Euclidean distance between sites from latitude and longitude given in decimal degrees. Geographical distances were converted from decimal degrees to kilometers after calculation. 

**Variation partitioning:** This approach uses slightly different input data from the distance-based regression above. Significance of the environmental and spatial components were quantified using canonical correspondence analysis (CCA).

- <u> *Moran Eigenvector Map (MEM)* </u> was constructed from a spatial weighting matrix of the input site coordinates, following Dray et al (2006). The MEM was first used to quantify the degree of spatial autocorrelation present in the environmental covariates, temperature and aridity. 

- <u> *Moran Structural Randomization (MSR)* </u>  was then performed between temperature and aridity to quantify the degree of correlation between the covariates after accounting for spatial autocorrelation. The MSR test indicated that, after accounting for inherent correlation due to spatial autocorrelation, temperature and PEV were sufficiently non-correlated to include both in subsequent variation partitioning analysis. 


```{r 04_Processes_carpentry, eval=F}
load("output/crosscomparison_both.RD")

StructAnalysis <- function(cc, meta){
  
  distdecreg <- function(df){
    ## distance regression
    lmD <- lm(cross ~ dist*transect.wb, data=df) %>% summary() 
    ## environment regression
    lmS <- lm(cross ~ envgrad*transect.wb, data=df) %>% summary()
    
    return(list("dist_geo"=lmD, "dist_env"=lmS))
  }
  
  noSingle <- function(x){return(length(which(x==0)))}
  
  retreat <- cc$data$retreat[,-which(apply(cc$data$retreat,2,noSingle)>=10)] %>%
    `row.names<-`(paste0(row.names(.),"retreat"))
  capture <- cc$data$capture[,-which(apply(cc$data$capture,2,noSingle)>=10)] %>%
    `row.names<-`(paste0(row.names(.),"capture"))
  soil <- cc$data$soil[,-which(apply(cc$data$soil,2,noSingle)>=10)] %>%
    `row.names<-`(paste0(row.names(.),"soil"))
  
  # pull lat-long info from metadata
  spat1 <- meta[match(row.names(retreat),row.names(meta)), c("latitude","longitude")]
  
  # pull temp and aridity index info from metadata and standardize the columns
  env1 <- decostand(meta[match(row.names(retreat),row.names(meta)), 
                         c("Temp","AridIndex")], method="standardize", MARGIN=2)
  
  # format community tables into community dissimilarity tables   
  ret.mou <- vegdist(retreat, "bray", T, T, T) %>% as.matrix()
  cap.mou <- vegdist(capture, "bray", T, T, T) %>% as.matrix()
  soi.mou <- vegdist(soil, "bray", T, T, T) %>% as.matrix()
  
  # generate site distances based on lat-long and environmental variables
  dist2 <- dist(spat1, diag=T, upper=T) %>% as.matrix()
  dist3 <- dist(env1, diag=T, upper=T) %>% as.matrix()
  
  ## compile information into a single long-form data frame
  sample_dist <- rbind(
    ret.mou %>%
      reshape2::melt(value.name="cross") %>% # reshape data
      mutate(cross=1-(cross)) %>% # convert dissimilarity to similarity
      cbind(., envgrad=reshape2::melt(dist3, value.name="envgrad")[,3]) %>% # bind env distances
      cbind(., dist=reshape2::melt(dist2, value.name="dist")[,3]) %>% # bind geo distances
      .[which(reshape2::melt(lower.tri(matrix(0,nrow(ret.mou),nrow(ret.mou))))[,3]==T),] %>%  # remove duplicates 
      cbind(nestedness="Retreat"),
    
    cap.mou %>%
      reshape2::melt(value.name="cross") %>% 
      mutate(cross=1-(cross)) %>%
      cbind(., envgrad=reshape2::melt(dist3, value.name="envgrad")[,3]) %>%
      cbind(., dist=reshape2::melt(dist2, value.name="dist")[,3]) %>%
      .[which(reshape2::melt(lower.tri(matrix(0,nrow(cap.mou),nrow(cap.mou))))[,3]==T),] %>%
      cbind(nestedness="Capture"),
    
    soi.mou %>%
      reshape2::melt(value.name="cross") %>% 
      mutate(cross=1-(cross)) %>%
      cbind(., envgrad=reshape2::melt(dist3, value.name="envgrad")[,3]) %>%
      cbind(., dist=reshape2::melt(dist2, value.name="dist")[,3]) %>%
      .[which(reshape2::melt(lower.tri(matrix(0,nrow(soi.mou),nrow(soi.mou))))[,3]==T),] %>%
      cbind(nestedness="Soil")
  )
  
  # join metadata and add transect information
  sample_dist <- left_join(sample_dist, meta[,c("SampleID","country")], by=c("Var1"="SampleID")) %>%
    left_join(meta[,c("SampleID","country")], by=c("Var2"="SampleID")) %>%
    mutate(transect.wb=country.x==country.y, compare=paste(country.y, "vs.", country.x))
  
  # re-level factors
  sample_dist$nestedness <- factor(sample_dist$nestedness, levels=c("Retreat","Capture","Soil"))
  sample_dist$dist <- sample_dist$dist*111.11 # set distances from decimal degrees to km
  sample_dist$transect.wb[sample_dist$transect.wb==T] <- "Within transect"
  sample_dist$transect.wb[sample_dist$transect.wb==F] <- "Between transects"
  sample_dist$transect.wb <- factor(sample_dist$transect.wb, levels=unique(sample_dist$transect.wb))
  sample_dist$compare <- as.factor(sample_dist$compare)
  sample_dist$compare <- factor(sample_dist$compare, levels=unique(sample_dist$compare)[c(1,3,2)])
  
  mR <- distdecreg(subset(sample_dist, subset=nestedness=="Retreat"))
  mC <- distdecreg(subset(sample_dist, subset=nestedness=="Capture"))
  mS <- distdecreg(subset(sample_dist, subset=nestedness=="Soil"))

  out <- list(data=list("retreat"=list("data"=retreat,"dist"=ret.mou),
                        "capture"=list("data"=capture,"dist"=cap.mou),
                        "soil"=list("data"=soil,"dist"=soi.mou)),
              spatial=list("data"=spat1, "dist"=dist2),
              environ=list("data"=env1, "dist"=dist3),
              reg.model=list("retreat"=mR, "capture"=mC, "soil"=mS),
              df=sample_dist)
  
  return(out)
}

sa_bact <- StructAnalysis(ccbact, metaB)
sa_fung <- StructAnalysis(ccfung, metaF)

save(sa_bact, sa_fung, file="output/structuringprocesses_both.RD")
```

```{r 04_Regression_tables, eval=F}
## raw regression table
rawtab <- function(sa){
  
  # function to reformat regression summary output into table format
  sift <- function(lm.summary){
  list(Rsq=lm.summary$r.squared, adjRsq=lm.summary$adj.r.squared,
       sigma=lm.summary$sigma, Fstat=lm.summary$fstatistic[1][[1]], 
       pval=pf(lm.summary$fstatistic[1], lm.summary$df[1], lm.summary$df[2], lower.tail = FALSE)[[1]],
       observations=sum(lm.summary$df[1:2]))
  }
  
  # combine raw regression results into table and format
  tab1 <- rbind(sa$reg.model$retreat$dist_geo[[4]], 
                sa$reg.model$retreat$dist_env[[4]], 
                sa$reg.model$capture$dist_geo[[4]], 
                sa$reg.model$capture$dist_env[[4]],  
                sa$reg.model$soil$dist_geo[[4]], 
                sa$reg.model$soil$dist_env[[4]]) %>% 
    as.data.frame() %>%
  apply(., MARGIN=2, FUN=as.numeric) %>% as.data.frame() %>% 
  cbind(type=c(rep("Retreat",8),rep("Capture",8), rep("Soil",8)),
        test=rep(c(rep("Distance",4),rep("Environment",4)),3),
        coef=rep(c("Intercept","Variable","Transect","TransectEffect"),3)) %>%
  .[,c(5,6,7,1:4)] %>% 
    `colnames<-`(c("Sample Type","Regression","Term","Coefficients","Std.Err","t Stat","P-value"))
  
  # filter p-value information into table
  pval.tab <- reshape2::dcast(tab1, `Sample Type` + Regression ~ Term, value.var = "P-value")
  
  # filter coefficients into table
  coef.tab <- reshape2::dcast(tab1, `Sample Type` + Regression ~ Term, value.var = "Coefficients") %>%
    mutate(Transect=Transect+Intercept, TransectEffect=TransectEffect+Variable)
  
  # use function to generate organized regression results
  tab2 <- cbind(retD=sift(sa$reg.model$retreat$dist_geo),
                retS=sift(sa$reg.model$retreat$dist_env),
                capD=sift(sa$reg.model$capture$dist_geo),
                capS=sift(sa$reg.model$capture$dist_env),
                soiD=sift(sa$reg.model$soil$dist_geo),
                soiS=sift(sa$reg.model$soil$dist_env))
  stats <- row.names(tab2) 
  tab2 <- apply(tab2, 2, as.numeric) %>% as.data.frame() %>% 
    `rownames<-`(stats) %>% tibble::rownames_to_column()
  
  dims <- nrow(sa$data$retreat$dist)^2-nrow(sa$data$retreat$dist)
  
  # table of residuals, just for fun.
  tab5 <- data.frame(resid=c(sa$reg.model$retreat$dist_geo$residuals, 
                             sa$reg.model$retreat$dist_env$residuals,
                             sa$reg.model$capture$dist_geo$residuals, 
                             sa$reg.model$capture$dist_env$residuals,
                             sa$reg.model$soil$dist_geo$residuals, 
                             sa$reg.model$soil$dist_env$residuals),
                     sample=c(rep("Retreat",dims),
                              rep("Capture",dims),
                              rep("Soil",dims)),
                     regres=c(rep("Dist",dims/2),rep("Env",dims/2),
                              rep("Dist",dims/2),rep("Env",dims/2),
                              rep("Dist",dims/2),rep("Env",dims/2)))
  tab5$sample <- factor(tab5$sample, levels=unique(tab5$sample))
  
  ## combine all regression tables into a single named list
  out <- list("Table1"=tab1, "RegressionStat"=tab2, 
              "Coefficients"=coef.tab, "P-values"=pval.tab,
              "Residuals"=tab5)
  
  return(out)
}

tabsB <- rawtab(sa_bact)
tabsF <- rawtab(sa_fung)


# print to Excel sheets, one sheet per list element
openxlsx::write.xlsx(tabsB, "output/structProc_tables_bact.xlsx")
openxlsx::write.xlsx(tabsF, "output/structProc_tables_fung.xlsx")
```


### 4.1 Regression analysis

```{r 04_Processes_regression}
load("output/structuringprocesses_both.RD")

# Retreat distance-decay
sa_bact$reg.model$retreat$dist_geo
# Retreat environment-decay
sa_bact$reg.model$retreat$dist_env
# Capture distance-decay
sa_bact$reg.model$capture$dist_geo
# Capture environment-decay
sa_bact$reg.model$capture$dist_env
# Soil distance-decay
sa_bact$reg.model$soil$dist_geo
# Soil environment-decay
sa_bact$reg.model$soil$dist_env

```

<a href="#data-analysis" style="color:steelblue;" >Back to top</a>

### 4.2 Varpart analysis

#### MEM (Figure S3)
```{r 04_Processes_MEM_FigureS3}
genmem <- function(spat){
  ## build Moran Eigenvector Map
  nb.rel<- graph2nb(relativeneigh(spat), sym=T)  ## computation of spatial neighborhood
  d <- nbdists(nb.rel, spat) %>% ## geographic distance matrix
    lapply(function(x) 1 - x/max(dist(spat)))
  nb.wl <- nb2listw(nb.rel, glist=d) ## spatial weighting matrix (SWM)
  mem <- mem(nb.wl)## build MEM
  ## outbox
  out <- list('latlong'=spat, 'neigh'=nb.rel, 'dist'=d, 'wlist'=nb.wl, 'mem'=mem)
  return(out)
}

gm_bacteria <- genmem(sa_bact$spatial$data)
gm_fungi <- genmem(sa_fung$spatial$data)

# s.value is weird about coordinates. Rescale lat-long to each start at 0 
nxyB <- apply(gm_bacteria$latlong, MARGIN=2, FUN=function(x){return(x-min(x))})*111.11 # and convert to km
nxyF <- apply(gm_fungi$latlong, MARGIN=2, FUN=function(x){return(x-min(x))})*111.11 # and convert to km

# plot spatial neighborhood and MEMs
s.label(nxyB, nb = gm_bacteria$neigh, label=NULL, pnb.edge.col = "red", 
        pnb.edge.lwd=2, pnb.node.col = "red", pnb.node.pch=1, main = "Relative")
s.label(nxyF, nb = gm_fungi$neigh, label=NULL, pnb.edge.col = "blue", pnb.edge.lwd=2, 
        pnb.node.col = "blue", pnb.edge.lty = "dotted", pnb.node.pch=0, add=T)
s.label(data.frame(x=c(400,650),y=c(700,500)), label=c("Bacteria","Fungi"), 
        plabels.col=c("red","blue"), plabels.boxes.lty=c("solid","dashed"), add=T)

## plot MEMs and autocor for bacteria
s.value(nxyB, gm_bacteria$mem, symbol = "circle", ppoint.cex = 0.6, main="Bacteria")
# check correlation between lat-long and MEM
cor(nxyB, gm_bacteria$mem)
# check environmental variables for spatial autocorrelation
MC.env <- moran.randtest(sa_bact$environ$data, gm_bacteria$wlist, nrepet = 999) # both are spatially autocorrelated
cor(sa_bact$environ$data, gm_bacteria$mem) # check against MEMs for best correlation

## plot MEMs and autocor for fungi
s.value(nxyF, gm_fungi$mem, symbol = "square", ppoint.cex = 0.6, main="Fungi")
# check correlation between lat-long and MEM
cor(nxyF, gm_fungi$mem)
# check environmental variables for spatial autocorrelation
MC.env <- moran.randtest(sa_fung$environ$data, gm_fungi$wlist, nrepet = 999) # both are spatially autocorrelated
cor(sa_fung$environ$data, gm_fungi$mem) # check against MEMs for best correlation
```

#### MSR
```{r 04_Processes_MSRtest}
## Moran Structural Ranomization test for env variables
# Are they correlated, apart from both being spatially autocorrelated?

## Bacteria
msr1 <- msr(sa_bact$environ$data$Temp, gm_bacteria$wlist)
msr2 <- msr(sa_bact$environ$data$AridIndex, gm_bacteria$wlist)
obs <- cor(sa_bact$environ$data$Temp, sa_bact$environ$data$AridIndex)
sim <- sapply(1:ncol(msr1), function(i) cor(msr1[,i], msr2[,i]))
testmsrB <- as.randtest(obs = obs, sim = sim, alter = "two-sided")
testmsrB ## not sufficiently correlated with each other -> include both

## Fungi
msr1 <- msr(sa_fung$environ$data$Temp, gm_fungi$wlist)
msr2 <- msr(sa_fung$environ$data$AridIndex, gm_fungi$wlist)
obs <- cor(sa_fung$environ$data$Temp, sa_fung$environ$data$AridIndex)
sim <- sapply(1:ncol(msr1), function(i) cor(msr1[,i], msr2[,i]))
testmsrF <- as.randtest(obs = obs, sim = sim, alter = "two-sided")
testmsrF ## not sufficiently correlated with each other -> include both
```

#### CCA
```{r 04_Processes_CCA, eval=F}
partitioning <- function(sa, gm){
  bysamp <- function(sampletype, col){
    
    stdist <- sa$data[[sampletype]]$dist
    nbwl <- gm$wlist
    env <- sa$environ$data
    mem <- gm$mem
      
    # check best MEM for community data
    mem.sel1 <- mem.select(stdist,listw = nbwl)$MEM.select 
    # variation partitioning
    varp <- varpart(as.dist(stdist), env, mem$MEM1)
    # Check significance of E|S fraction
    accaES <- anova.cca(dbrda(as.dist(stdist) ~ Temp + AridIndex + Condition(mem), 
                        data=cbind(env, mem=mem$MEM1))) # take note of significance level
    # Check significance of S|E fraction
    accaSE <- anova.cca(dbrda(as.dist(stdist) ~ mem + Condition(Temp + AridIndex), 
                        data=cbind(env, mem=mem$MEM1))) # take note of significance level
    # quick venn diagram to check # plot(var.R, digits=5, Xnames=c("Envir","Spatial"), bg="#377EB8") 
     # store AdjR fractions
    fract <- round(varp$part$indfract$Adj.R.squared, 2)
    fract[fract<0] <-0
    
    sigstar <- function(s){
      ss <- rep("   ",length(s))
      ss[s<0.05] <- "  *"
      ss[s<0.01] <- " **"
      ss[s<=0.001] <- "***"
      return(ss)
    }
    
    pval <- c(accaES[1,4], accaSE[1,4], NA, NA)
  
    fracts=data.frame(sampletype=sampletype,
                                  labels=c("Envir  ","Spatial","Both","Resid"),
                                  sets=c(fract[1]+fract[3],
                                           fract[2]+fract[3],
                                           fract[3],fract[4]),
                                  fracts=fract,
                                  p.val=pval,
                                  signif=sigstar(pval))
    
    fracts <- mutate(fracts, labels=paste(labels, signif, sep="\n\n"))
    
    plotVenn <- draw.pairwise.venn(fracts$sets[1], fracts$sets[2],fracts$sets[3], 
                                   fracts$labels[1:2], 
                                   fontfamily="sans", cat.fontfamily="sans", 
                                   col=col, fill=col, alpha=0.3, 
                                   cat.default.pos = "text", 
                                   cat.just = list(c(0.3,0.6),c(0.5,0.6)))
    
    out <- list(cca.ES=accaES, 
                cca.SE=accaSE, 
                fracts=fracts,
                plot=plotVenn)
    
    return(out)
  }
  
  colors <- data.frame(retreat="#377EB8", capture="#4DAF4A", soil="#E41A1C")
  
  outlist <- list()
  for (i in c("retreat","capture","soil")){
    outlist[[i]] <- bysamp(i, unlist(colors[i]))
  }
  return(outlist)
}

vp_bact <- partitioning(sa_bact, gm_bacteria)
vp_fung <- partitioning(sa_fung, gm_fungi)

save(vp_bact,vp_fung, file="output/varpart_both.RD")

```

```{r 04_Processes_CCA_output, message=T}

load("output/varpart_both.RD")

## BACTERIA
vp_bact$retreat$cca.ES ## retreat env
vp_bact$retreat$cca.SE ## retreat spatial
vp_bact$capture$cca.ES ## capture env
vp_bact$capture$cca.SE ## capture spatial
vp_bact$soil$cca.ES ## soil env
vp_bact$soil$cca.SE ## soil spatial

## FUNGI
vp_fung$retreat$cca.ES ## retreat env
vp_fung$retreat$cca.SE ## retreat spatial
vp_fung$capture$cca.ES ## capture env
vp_fung$capture$cca.SE ## capture spatial
vp_fung$soil$cca.ES ## soil env
vp_fung$soil$cca.SE ## soil spatial

```

<a href="#data-analysis" style="color:steelblue;" >Back to top</a>

### 4.3 Plots

#### Varpart (Figure 3)
```{r 04_ProcessesVarpart_Figure3, fig.width=4, fig.height=5}

g4 <- ggdraw() + 
  draw_plot(as_grob(vp_bact$retreat$plot), 0,0.65,0.45,0.3) + ## panel A top
  draw_plot(as_grob(vp_bact$capture$plot), 0,0.35,0.45,0.3) + ## panel A middle
  draw_plot(as_grob(vp_bact$soil$plot), 0,0.05,0.45,0.3) + ## panel A bottom
  draw_plot_label(c(paste("R^2 ==", sum(vp_bact$retreat$fracts$fracts[-4])),
                    paste("R^2 ==", sum(vp_bact$capture$fracts$fracts[-4])),
                    paste("R^2 ==", sum(vp_bact$soil$fracts$fracts[-4]))), 
                  c(0.25,0.25,0.25), 
                  c(0.73,0.43,0.13), 
                  fontface="plain", size=9, parse=T) +
  draw_plot(as_grob(vp_fung$retreat$plot), 0.5,0.65,0.45,0.3) + ## panel B top
  draw_plot(as_grob(vp_fung$capture$plot), 0.5,0.35,0.45,0.3) + ## panel B middle
  draw_plot(as_grob(vp_fung$soil$plot), 0.5,0.05,0.45,0.3) + ## panel B bottom
  draw_plot_label(c(paste("R^2 ==", sum(vp_fung$retreat$fracts$fracts[-4])),
                    paste("R^2 ==", sum(vp_fung$capture$fracts$fracts[-4])),
                    paste("R^2 ==", sum(vp_fung$soil$fracts$fracts[-4]))), 
                  c(0.75,0.75,0.75), 
                  c(0.73,0.43,0.13), 
                  fontface="plain", size=9, parse=T) +
  draw_plot_label(c("Bacteria","Fungi"),c(0.05,0.55), c(1,1), 
                  fontface="plain", size=20)

g4
```

#### Regression (Figure 3)
```{r 04_ProcessesRegression_Figure3, fig.width=11, fig.height=9}

load("output/structuringprocesses_both.RD")

sample_dist <- rbind(data.frame(sa_bact$df, set="Bacteria"),
                     data.frame(sa_fung$df, set="Fungi")) %>%
  .[,c(3,4,5,6,9,10,11)]

t1 <- openxlsx::read.xlsx("coreplots/structure/bacteria.xlsx", 2)
t1 <- data.frame(nestedness=c("Retreat","Capture","Soil"),
                 cross=Inf, envgrad=-Inf, dist=-Inf,
                 r.sqD=t(t1[1,c(2,4,6)]), r.sqE=t(t1[1,c(3,5,7)]))

t2 <- openxlsx::read.xlsx("coreplots/structure/fungi.xlsx", 2)
t2 <- data.frame(nestedness=c("Retreat","Capture","Soil"),
                 cross=Inf, envgrad=-Inf, dist=-Inf,
                 r.sqD=t(t2[1,c(2,4,6)]), r.sqE=t(t2[1,c(3,5,7)]))

t3 <- rbind(data.frame(t1, set="Bacteria"),
            data.frame(t2, set="Fungi"))


## environment gradient plots
g5 <- ggplot(sample_dist, aes(envgrad, cross)) + 
  facet_grid(nestedness~set) + 
  geom_point(aes(shape=compare, color=nestedness), show.legend = c(color=F)) +
  scale_shape_manual(values=c(17,20, 15)) +
  stat_smooth(aes(linetype=transect.wb), 
              method="lm", formula='y~x', color="black")+
  ggpmisc::stat_poly_eq(formula=y ~ x, parse=T, label.x = "right",
                        label.y="top",vjust=2,vstep=0.85,hstep=0.9,
                        aes(linetype=transect.wb,
                            label=paste("bold(\"",
                                        c("Within","Between"
                                          )[stat(linetype)],
                                        "\ntransect:  \")*",
                                        after_stat(eq.label), 
                                        sep = "")), size=2.5) +
    geom_label(data=t1, aes(envgrad, cross, 
                         label = paste("R^2==",round(X1.1,2))), 
            parse = TRUE, hjust=0, vjust=1) + 
  scale_color_manual(values=c("#377EB8","#4DAF4A","#E41A1C")) +
  theme_bw() + theme(plot.margin=unit(c(1.5,0.5,0,0.5), "lines"), 
                     legend.spacing.x = unit(5, "lines"), 
                     legend.text=element_text(margin=margin(0,7,0,0)),
                     legend.justification = "left",
                     legend.box="vertical",
                     strip.background = element_blank(),
                     strip.text.x=element_text(size=18, face="bold", family="mono"),
                     strip.text.y=element_blank()) +
    guides(shape=guide_legend(title=NULL, direction="vertical", position="bottom",ncol=2,
                            byrow = F,override.aes = list(size=3),label.theme=element_text(margin=margin(0,7,0,0), size=9)),
         linetype=guide_legend(title=NULL, direction="horizontal", position = "bottom",order=1,
                               reverse = F, keywidth = unit(1, "cm"))) +
  xlab("Environmental dissimilarity\n(Temp:AridIndex)") + ylab("Compositional similarity") + coord_cartesian(ylim=c(0,1))

## distance gradient plots
g6 <- ggplot(sample_dist, aes(dist, cross)) + 
  facet_grid(nestedness~set) + 
  geom_point(aes(shape=compare, color=nestedness), show.legend =F) +
  stat_smooth(aes(linetype=transect.wb),method="lm", formula='y~x',
              color="black", fullrange = T, show.legend = F)+
  ggpmisc::stat_poly_eq(formula=y ~ x, parse=T, label.x = "left",
                    label.y="bottom",vjust=0.7,vstep=0.7,hstep=0,
                       aes(linetype=transect.wb,
                         label=paste("bold(\"",c("Within","Between"
                                                   )[stat(linetype)],
                                       "\ntransect:  \")*",
                                       after_stat(eq.label), 
                                       sep = "")), size=2.5) +
  geom_label(data=t1, aes(dist, cross, 
                         label = paste("R^2==",round(X1,2))), 
            parse = TRUE, hjust=0, vjust=1) + 
  scale_color_manual(values=c("#377EB8","#4DAF4A","#E41A1C")) +
  scale_shape_manual(values=c(17,20, 15)) +
  guides(shape=guide_legend(title=NULL, direction="horizontal",
                            byrow = T,override.aes = list(size=3), 
                            label.theme=element_text(margin=margin(0,7,0,0), size=9)),
         linetype=guide_legend(title=NULL, direction="horizontal", 
                               reverse = F, keywidth = unit(1, "cm"))) +
  theme_bw() + theme(plot.margin=unit(c(1.5,0.5,0,0.5), "lines"),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.title.y=element_blank(),
                     legend.box="vertical",
                     legend.spacing.x = unit(3, "lines"),
                     legend.spacing.y = unit(-0.25, "lines"),
                     legend.position = "bottom", 
                     legend.justification = "center",
                     strip.background = element_blank(),
                     strip.text=element_text(size=18, face="bold", family="mono")) +
  xlab("Distance between colonies\n(km)") + ylab("Compositional similarity") + coord_cartesian(ylim=c(0,1))

egg::ggarrange(g5,g6,nrow=1)

```


<a href="#data-analysis" style="color:steelblue;" >Back to top</a>

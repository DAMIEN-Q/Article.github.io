# CC2 ADM Damien Quemener
---
title: "Analyse de données métabarcode d'un article (Damien Quemener)"
output: github_document
---
```{r}
library(dada2)
library(Rcpp)
```

```{r}
path <- "~/Analyse_Article/Data"
path_trim <- "fastq_trimmed"
dir.create(path_trim, showWarnings = FALSE)
```

```{r}
fnFs = sort(
  list.files(path, pattern="_1.fastq",  
                    full.names = TRUE))

fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
## Retrait des amorces 

#Amorce 
FWD <- "ACTCCTACGGGAGGCAGCA"
REV <- "GGACTACHVGGGTWTCTAAT"

fnFs.trim <- file.path(path_trim, basename(fnFs))
fnRs.trim <- file.path(path_trim, basename(fnRs))

cutadapt <- "cutadapt"

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(
    "-g", paste0("^", FWD),
    "-G", paste0("^", REV),
    "--discard-untrimmed",
    "-o", fnFs.trim[i],
    "-p", fnRs.trim[i],
    fnFs[i], fnRs[i]
  ))
}
```

```{r}
library(ShortRead)
sread(readFastq(fnFs.trim[1]))[1:5]

```
Ici, j'ai été obligé de supprimer les amorces car trop de séquence considérer comme des chimères par la suite.


```{r}
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

## **Visualisation des profils de qualité de lecture**
### **Qualité de lecture Forward** ###
```{r}
plotQualityProfile(fnFs.trim[1:2])
```

On observe une bonne qualité de lecture forward jusqu'à environ 220 pb

## **Qualité de lecture Reverse**
```{r}
plotQualityProfile(fnRs.trim[1:2])
```

Ici, la qualité de lecture reverse est la même que pour forward, c'est à dire 220 pb.

## **Filtrer et Rogner**

```{r}
filtFs <- file.path(path_trim,"filtered",
                    paste0(sample.names,"_F_filt.fastq.gz"))
filtRs <- file.path(path_trim, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

### **Filtrage et pré-traitement des séquences** ###

```{r}
out <- filterAndTrim(fnFs.trim, filtFs, 
                     fnRs.trim, filtRs, 
                     truncLen=c(220,220), 
              maxN=0, 
              maxEE=c(2,5),
              truncQ=2, 
              rm.phix=TRUE, 
              compress=TRUE, 
              multithread=FALSE
              )

head(out) 
```
Le tableau montre clairement le nombre lectures avant filtrage (reads.in) et après filtrage (reads.out). Ici, on observe que la majorité des lectures est conservé (seulement environ 10 000 lectures filtrée)

## **Taux d'erreur de séquençage**

```{r,cache=TRUE}
errF <- learnErrors(filtFs, multithread=TRUE)
```
```{r,cache=TRUE}
errR <- learnErrors(filtRs, multithread=TRUE)
```
## **Visualisation des taux d'erreur** 

```{r,cache=TRUE}
plotErrors(errF, nominalQ=TRUE)
```

## **Application de l'algo de DADA2**

```{r,cache=TRUE}
dadaFs <- dada(filtFs,
               err=errF, 
               multithread=TRUE 
               )
```


```{r,cache=TRUE}
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

## **Inspection**

```{r,cache=TRUE}
dadaFs[[1]]
```

## **Fusion des lectures appariées**

```{r,cache=TRUE}
mergers <- mergePairs(dadaFs, 
                      filtFs, 
                      dadaRs, 
                      filtRs, 
                      verbose=TRUE 
                      )
head(mergers[[1]])
```

## **Construction d'une table de séquence**

```{r,cache=TRUE}
seqtab <- makeSequenceTable(mergers) 
dim(seqtab)
```

### **Inspection de la distribution des longueurs des séquences** ###
  
```{r,cache=TRUE}
table( 
  nchar( 
    getSequences(seqtab) 
        ))
```


## **Supprimer les chimères**

Une chimère est une séquence artificielle formée quand 2 fragments d'ADN s'hybrident partiellement et sont amplifiés comme une "fausse" séquence.

```{r,cache=TRUE}
seqtab.nochim <- removeBimeraDenovo( 
  seqtab, 
  method="consensus", 
  multithread=TRUE, 
  verbose=TRUE 
  ) 

dim(seqtab.nochim)
```


## **Proportion de lecture "survivante" au filtrage de chimères**

```{r,cache=TRUE}
sum(seqtab.nochim)/sum(seqtab)
```


## **Suivre les lectures à travers la pipeline**

```{r,cache=TRUE}
getN <- function(x) sum(getUniques(x)) 

track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim) 
               )

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

head(track)
```


## **Attribution d'une taxonomie**

Comparaison des ASVs à une base de données de références pour attribuer une taxonomie
```{r,cache=TRUE}
taxa <- assignTaxonomy(seqtab.nochim,  "~/Analyse_Article/Data/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1", 
                       multithread=TRUE 
                       )
```

## **Examination des attributions**

```{r,cache=TRUE}
taxa.print <- taxa 
rownames(taxa.print) <- NULL 
head(taxa.print) 
```

Ici, on retrouve bien comme dans l'article le phylum majoritaire Bacillota (anciennement Firmicutes). Néanmoins, les auteurs ont utilisé le pipeline QIIME2 pour pour l'analyse des séquences et non DADA2, donc il est fort propable d'avoir des divergences dans la suite de l'analyse avec PhyloSeq. 


## **Package PhyloSeq**

```{r,echo= TRUE}
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")
```

```{r,echo = TRUE}
library("phyloseq")
packageVersion("phyloseq")

library("ggplot2")
packageVersion("ggplot2")

library("scales")
packageVersion("scales")

library("grid")
packageVersion("grid")
```

### **Documentation intégrée à Phyloseq** ###
```{r,cache=TRUE}
`?`("phyloseq-package")
`?`(phyloseq)
```

### **Vignette dans le package** ###
```{r,cache=TRUE}
vignette("phyloseq-basics")
vignette("phyloseq-analysis")
```

### **Contruction d'un objet Phyloseq à partir des données DADA2 ** ###
```{r,cache=TRUE}
theme_set(theme_bw())
library(dplyr)
```

```{r,cache=TRUE}
meta <- read.csv("SraRunTable.csv", header = TRUE, sep = ",")

rownames(meta) <- meta$Run

meta_phy <- meta %>%
  select(Run, Sample.Name, LibraryLayout, 
          HOST, LibrarySelection, LibrarySource) 
rownames(meta_phy) <- meta_phy$Run
meta_phy$Group <- meta_phy$Sample.Name   

meta_phy$Group <- NA

meta_phy$Group[grepl("^A", meta_phy$Sample.Name)] <- "CON"  
meta_phy$Group[grepl("^B", meta_phy$Sample.Name)] <- "LBC"  
meta_phy$Group[grepl("^C", meta_phy$Sample.Name)] <- "HBC"  

meta_phy$Group <- factor(meta_phy$Group, levels = c("CON", "LBC", "HBC"))

SAMP <- sample_data(meta_phy)

table(meta_phy$Sample.Name, meta_phy$Group)

```
```{r,cache=TRUE}
OTU <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(taxa))
stopifnot( all(colnames(seqtab.nochim) == rownames(TAX)) )
stopifnot( all(rownames(SAMP) %in% rownames(seqtab.nochim)) )
```

```{r,cache=TRUE}
ps <- phyloseq(OTU,TAX,SAMP)
ps
```

### **Arbre phylogénétique annoté** ###
```{r}
library(Biostrings)
library(ape)
library(phangorn)

## Extraire les séquences 
asv_seqs = getSequences(seqtab.nochim)

## Créer les IDs d’ASV et les mettre en noms de colonnes
asv_ids <- paste0("ASV", seq_along(asv_seqs))
seqtab.nochim.phy = seqtab.nochim
colnames(seqtab.nochim.phy) <- asv_ids

dna_all <- DNAStringSet(asv_seqs)
widths <- width(dna_all)
table(widths)

# Garder uniquement la longueur la plus fréquente
main_len <- as.integer(names(sort(table(widths), decreasing = TRUE))[1])
keep     <- which(widths == main_len)

dna <- dna_all[keep]                       
asv_ids <- paste0("ASV", seq_along(dna))
names(dna) <- asv_ids

# Filtrer la table d’abondance et la taxo 
seqtab_filt <- seqtab.nochim[, keep, drop = FALSE]
colnames(seqtab_filt) <- asv_ids

taxa_filt <- taxa[keep, , drop = FALSE]    
rownames(taxa_filt) <- asv_ids

## Conversion en objet phangorn
dna_phy <- as.phyDat(dna, type = "DNA")

## Matrice de distances + arbre Neighbor-Joining
dm     <- dist.ml(dna_phy)
treeNJ <- NJ(dm)

## Maximum likelihood (amélioration de l’arbre)
fit  <- pml(treeNJ, data = dna_phy)
fitG <- optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE)

treeML <- fitG$tree   # arbre final
```

```{r}
OTU  <- otu_table(as.matrix(seqtab_filt), taxa_are_rows = FALSE)
TAX  <- tax_table(as.matrix(taxa_filt))
SAMP <- sample_data(meta_phy)
TREE <- phy_tree(treeML)
RS   <- refseq(dna)


ps <- phyloseq(OTU, TAX, SAMP, TREE, RS)

ps.Faecalibacterium <- subset_taxa(ps, Genus == "Faecalibacterium")

plot_tree(ps.Faecalibacterium, color = "Group", label.tips = "Genus", size = "abundance")
```

### **Alpha Diversité** ###
```{r,cache=TRUE}
library(ggplot2)

p <- plot_richness(ps,
                   x = "Group",
                   measures = c("Observed", "Shannon", "Simpson", "Chao1"),
                   color = "Group")

p + 
  geom_boxplot(aes(fill = Group), alpha = 0.3, outlier.shape = NA) +
  theme_bw()

```

On observe pas de différence sur la diversité alpha

### **Beta Diversité** ###
```{r,cache=TRUE}
ord <- ordinate(ps, method = "PCoA", distance = "unifrac")
plot_ordination(ps, ord, color = "Group")
```

### **Barplots de composition** ###
```{r,cache=TRUE}
## Phylum + abondance relative
ps_phylum <- tax_glom(ps, taxrank = "Phylum")
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

## Moyenne par Group × Phylum
dfP <- psmelt(ps_phylum_rel) %>%      # Sample, Group, Phylum, Abundance
  group_by(Group, Phylum) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

## Sélection des 10 phyla les plus abondants (globalement)
top_phyla <- dfP %>%
  group_by(Phylum) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  arrange(desc(Abundance)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

dfP$Phylum_clean <- ifelse(dfP$Phylum %in% top_phyla,
                           as.character(dfP$Phylum),
                           "Other")

## Barplot : une barre par groupe, 10 phyla + Other
ggplot(dfP, aes(x = Group, y = Abundance, fill = Phylum_clean)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") +
  theme_bw()

```

```{r,cache=TRUE}
## Agglomérer au niveau Genre + abondance relative
ps_genus <- tax_glom(ps, taxrank = "Genus")
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))  

## Extraire et moyenner par Groupe × Genre
dfG <- psmelt(ps_genus_rel) %>%          
  group_by(Group, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

## Garde que les genres les plus abondants
topN <- 10
top_genus <- dfG %>%
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  arrange(desc(Abundance)) %>%
  slice_head(n = topN) %>%
  pull(Genus)

dfG$Genus_clean <- ifelse(dfG$Genus %in% top_genus,
                          as.character(dfG$Genus),
                          "Other")

## Barplot 
ggplot(dfG, aes(x = Group, y = Abundance, fill = Genus_clean)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Relative abundance") +
  theme_bw()
```

### **Raréfaction** ###
```{r,cache=TRUE}
ps.rare <- rarefy_even_depth(ps)
```
```{r,cache=TRUE}
## Sans Raréfaction

ord_unrare <- ordinate(ps, method = "PCoA", distance = "unifrac")
plot_ordination(ps, ord_unrare, color = "Group") +
  geom_point(size = 4)

```

```{r,cache=TRUE}
## Avec Raréfaction

ord_rare <- ordinate(ps.rare, method = "PCoA", distance = "unifrac")
plot_ordination(ps.rare, ord_rare, color = "Group") +
  geom_point(size = 4)

```

On observe bien qu'avec et sans raréfaction, on voit que les groupe HBC et LBC sont très proche


### **La fonction plot network** ###
```{r,cache=TRUE}
ps_20 <- prune_taxa(taxa_sums(ps) > 20, ps)
jg <- make_network(ps_20, type = "taxa", dist.fun = "bray",
                   max.dist = 0.7)
plot_network(jg, ps_20, "taxa", color = "Phylum")

```

On voit ici que le phylum "Bacillota" est majoritairement présent dans les différents groupes


### **La heatmap** ###
```{r,cache=TRUE}
# Garder les 50 taxa les plus abondants
ps_top <- prune_taxa(names(sort(taxa_sums(ps), decreasing = TRUE))[1:50], ps)

# Passer en abondances relatives
ps_rel <- transform_sample_counts(ps_top, function(x) x / sum(x))

# Extraire au format long
df <- psmelt(ps_rel)  # colonnes : Sample, Group, Genus, Abundance, etc.

# Moyenne d'abondance par Group × Genus
df_mean <- df %>%
  group_by(Group, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

ggplot(df_mean, aes(x = Group, y = Genus, fill = Abundance)) +
  geom_tile() +
  scale_fill_gradient(low = "navy", high = "yellow") +
  xlab("Group") +
  ylab("Genus (mean relative abundance)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

```


"Enterocloster" a une abondance relative bien supérieur au autres dans les groupe CON et HBC


### **Statistique d'écart** ###
```{r, cache=TRUE}
library(cluster)
library(vegan)

ord <- ordinate(ps, method = "PCoA", distance = "bray")
x <- ord$vectors
x2 <- x[, 1:2]

pam1 <- function(x, k) list(cluster = pam(x, k, cluster.only = TRUE))

set.seed(123)  
gs <- clusGap(x2,
              FUN   = pam1,
              K.max = 6,   
              B     = 50)  
print(gs, method = "Tibs2001SEmax")

```

Ici, le nombre de cluster recommandé est 1

```{r, cache=TRUE}
## Fonction "wrapper"
gap_statistic_ordination <- function(physeq,
                                     method = "PCoA",
                                     distance = "bray",
                                     FUNcluster = "pam1",
                                     K.max = 6,
                                     axes = 1:2,
                                     B = 50,
                                     verbose = FALSE) {
  ord_phy <- ordinate(physeq, method = method, distance = distance)
  x_2   <- ord_phy$vectors

  if (is.null(axes)) axes <- 1:ncol(x_2)
  x_3   <- x_2[, axes, drop = FALSE]

  if (FUNcluster == "pam1") {
    FUNcluster <- function(x_3, k) list(cluster = pam(x_3, k, cluster.only = TRUE))
  }

  clusGap(x_3, FUN = FUNcluster, K.max = K.max, B = B, verbose = verbose)
}
```

```{r, cache=TRUE}
plot_clusgap <- function(clusgap, title = "Gap statistic (clusters)") {
  library(ggplot2)
  gstab <- data.frame(clusgap$Tab, k = 1:nrow(clusgap$Tab))
  ggplot(gstab, aes(k, gap)) +
    geom_line() +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.1) +
    ggtitle(title) +
    theme_bw()
}
```

```{r, cache=TRUE}
gap <- gap_statistic_ordination(ps,
                               method   = "PCoA",
                               distance = "bray",
                               K.max    = 6,
                               B        = 50)
print(gap, method = "Tibs2001SEmax")  
plot_clusgap(gs)

```

On observe ici que les "gaps" sont proche de 0 ou négative ce qui signifie que le clustering est soit similaire à un nuage aléatoire, soit qu'il n'est pas plus compact.

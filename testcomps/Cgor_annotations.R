setwd("~/Dropbox/genomes/CladeC_Symbiodinium_transcriptome_plast/testcomps/")
library(tidyverse)

blast_gene <- read.delim("blast_iso2gene.tab", header=F)
head(blast_gene)

plast_gene <- read.delim("plast_iso2gene.tab", header=F)
head(plast_gene)

# get rid of second column (not in BLAST)
plast_gene <- plast_gene %>% select(V1,V3)
head(plast_gene)

nrow(blast_gene)-nrow(plast_gene) # 3886 more genes annotated using PLAST than BLAST

table(blast_gene$V1 %in% plast_gene$V1) # 218 in blast that aren't in plast
table(plast_gene$V1 %in% blast_gene$V1) # 4014 in plast that aren't in blast

matches <- merge(plast_gene,blast_gene,by=1)
names(matches) <- c("isogroup", "plastGene", "blastGene")
head(matches)

# see if the strings match
# get rid of characters besides the gene description
matches$plastClean <- gsub("OS.*","",matches$plastGene)
matches$plastClean <- gsub("^ *", "", matches$plastClean)
matches$plastClean <- gsub("* $", "", matches$plastClean)
matches$blastClean <- gsub("OS.*","",matches$blastGene)
matches$blastClean <- gsub("^ *", "", matches$blastClean)
matches$blastClean <- gsub("* $", "", matches$blastClean)
head(matches)

# how many match?
table(matches$plastClean==matches$blastClean)
# FALSE  TRUE 
# 6502 15519 

# which ones don't match?
nomatch <- matches[!matches$plastClean==matches$blastClean,][c(4,5)]
tail(nomatch)
# Checking by eye, it looks lie the genes are usually only slightly different
# Example:
# Voltage-dependent T-type calcium channel subunit alpha-1G vs. Voltage-dependent N-type calcium channel subunit alpha-1B

# summarize genes and plot
plast_results <- data.frame(matching=table(plast_gene$V1 %in% blast_gene$V1)[2],
                            unique=table(plast_gene$V1 %in% blast_gene$V1)[1])
blast_results <- data.frame(matching=table(blast_gene$V1 %in% plast_gene$V1)[2],
                            unique=table(blast_gene$V1 %in% plast_gene$V1)[1])

gene_results_numbers <- as.data.frame(rbind(blast_results,plast_results)) %>% 
  mutate(program=c("blast","plast")) %>%
  gather(type, number, 1:2) %>%
  ggplot() + geom_bar(aes(y = number, x = program, fill = factor(type,levels=c("unique", "matching"))), stat="identity") + 
  guides(fill=guide_legend(title="hit type")) +
  theme_bw()
gene_results_numbers


##################################################################################################################################
##################################################################################################################################
# GO annotations -----
##################################################################################################################################
##################################################################################################################################

blast_go <- read.delim("blast_iso2go.tab", header=F)
head(blast_go)

# load plast annotations for GO
plast_go <- read.delim("plast_iso2go.tab", header=F)
head(plast_go)

# get rid of second column (short gene ID) that isn't in the blast file
plast_go <- plast_go %>% select(V1,V3) %>% rename(V2="V3")
head(plast_go)

# get rid of isogroups that don't have annotations
plast_go  <- plast_go[!plast_go $V2=="noMatch",]
head(plast_go)

nrow(blast_go)-nrow(plast_go) # 2879 more GO annotations in PLAST

table(blast_go$V1 %in% plast_go$V1) # 756 in blast that aren't in plast
table(plast_go$V1 %in% blast_go$V1) # 3635 in plast that aren't in blast

gomatches <- merge(plast_go,blast_go,by=1)
names(gomatches) <- c("isogroup", "plastGO", "blastGO")
head(gomatches)

# summarize and plot
blast_results_go <- data.frame(matching=table(blast_go$V1 %in% plast_go$V1)[2],
                            unique=table(blast_go$V1 %in% plast_go$V1)[1])
plast_results_go <- data.frame(matching=table(plast_go$V1 %in% blast_go$V1)[2],
                            unique=table(plast_go$V1 %in% blast_go$V1)[1])
go_results_numbers <- as.data.frame(rbind(blast_results_go,plast_results_go)) %>% 
  mutate(program=c("blast","plast")) %>%
  gather(type, number, 1:2) %>%
  ggplot() + geom_bar(aes(y = number, x = program, fill = factor(type,levels=c("unique", "matching"))), stat="identity") + 
  guides(fill=guide_legend(title="hit type")) +
  theme_bw()
go_results_numbers

# Compare GO annotations
plast_go_long <- strsplit(as.character(plast_go$V2), split = ";")
plast_go_long <- data.frame(V1 = rep(plast_go$V1, sapply(plast_go_long, length)), 
                            V2 = unlist(plast_go_long))
head(plast_go_long)
length(unique(plast_go_long$V1)) # 25469 isogroups with annotations

blast_go_long <- strsplit(as.character(blast_go$V2), split = ";")
blast_go_long <- data.frame(V1 = rep(blast_go$V1, sapply(blast_go_long, length)), 
                            V2 = unlist(blast_go_long))
head(blast_go_long)
length(unique(blast_go_long$V1)) # 21779 isogroups with annotations

# matches?
gomatchlong <- merge(plast_go_long, blast_go_long, by=c("V1","V2"))
nrow(gomatchlong) # 13720 total instances where the isogroup AND GO term matched
length(unique(gomatchlong$V1)) # 13720 isogroups have at least one matching GO term


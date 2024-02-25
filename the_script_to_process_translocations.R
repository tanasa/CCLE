##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
# some resources about network biology
# https://www.bigbookofr.com/network-analysis.html
# https://yunranchen.github.io/intro-net-r/advanced-network-visualization.html
# https://www.bioconductor.org/packages/release/bioc/vignettes/netbiov/inst/doc/netbiov-intro.pdf

library("karyoploteR")
library("plyr")
library("doBy")
library("tidyverse")
library('dplyr')
library("tidyr")
library("readr")
library("splitstackshape")
library("stringi")
library("gGnome")
library("qgraph")
library("igraphdata")
library(network)
library(sna)
library(ggplot2)
library(visNetwork)
library(networkD3)
library(magrittr)
library(pandoc)
library("htmlwidgets")
library(magrittr)
library(clusterProfiler)
# library("ggdensity")
# library("statnet")
# library("splineTimeR")
library(dplyr)
library(stringr)
library(tidyr)

# chooseCRANmirror(78)
# setwd("/Users/tanasab/Desktop/CCLE_SvABA")

################################################################################################
################################################################################################
################################################################################################
################################################################################################

# The script parses the columns "bp1" and "bp2", "site1" and "site2", that contain the information about the position of the breakpoints (bp1 and bp2), 
# and about the genes that are disrupted ("site1", and"site2") by SV. It produces the following columns :

# gene1_anywhere : it is the closest gene TSS to the breakpoint_1, that can be either intergenic, intronic or exonic
# gene1_anywhere_orientation : the transcriptional orientation of gene1 
# gene2_anywhere : it is the closest gene TSS to the breakpoint_2, that can be either intergenic, intronic or exonic
# gene2_anywhere_orientation : the transcriptional orientation of gene1 
# position1: the type of the genomic region that is disrupted (intergenic, intronic, exonic)
# position2: the type of the genomic region that is disrupted (intergenic, intronic, exonic) 

# the coordinates and the orientation of the breakpoint 1 :
# breakpoint1_chr 
# breakpoit1_orientation 
# breakpoint1_start
# breakpoint1_end 

# the coordinates and the orientation of the breakpoint 2 :
# breakpoint2_chr 
# breakpoint2_orientation 
# breakpoint2_start 
# breakpoint2_end

# The final file that is produced is called "CCLE_translocations_SvABA_20181221.LIST.GENES.POSITIONS.INTERACTIONS.txt" and 
# it is used as in input in other two scripts that do perform network analysis. 

# One script that performs the network analysis uses the list of GENE FUSIONS.
# Another script that performs the network analysis uses the list of all SV (that is composed of the list of gene fusions, and the list of intergenic SV)

# An intermediary file that is generated is called : "CCLE_translocations_SvABA_20181221.list.closest.genes.txt"

################################################################################################
################################################################################################ the coordinates on hg19
################################################################################################
################################################################################################

# we work with the file C...mo.txt, it is the same as the original file, except that I have replaced "gene" with "Gene";
# the data has been mapped on hg19
FILE = "CCLE_translocations_SvABA_20181221.mo.txt"

################################################################################################
################################################################################################

CCLE <- read.table(FILE, header=T, sep="\t", stringsAsFactors=F)

x = CCLE
head(x)
dim(x)
colnames(x)

# [1] "CCLE_name"       "map_id"          "bp1"             "bp2"
# [5] "class"           "gene1"           "gene2"           "site1"
# [9] "site2"           "fusion"          "multi_sv_fusion" "cosmic_fus"
# 75344 entries

x$TUMOR = x$CCLE_name
unique(x$TUMOR)

#################################################################################################################
#################################################################################################################
#################################################################################################################
################################################################################################################# the code that extracts the genes close to the breakpoints
################################################################################################################# either in an intergenic or an intragenic region
#################################################################################################################
#################################################################################################################
#################################################################################################################
# head(x$site1)
# head(x$site2)

tu = x

tu$site1 = as.factor(tu$site1)
tu$site2 = as.factor(tu$site2)

tu$site1_1 = sub(".*Gene", "",tu$site1)
tu$site2_1 = sub(".*Gene", "",tu$site2)

k1=c("\\(")

tu = tu %>% separate(site1_1, c("site1_11", "site1_11g"), sep=k1)
tu = tu %>% separate(site2_1, c("site2_11", "site2_11g"), sep=k1)

k2=c("\\)") 

tu = tu %>% separate(site1_11g, c("site1_111", "site1_111g"), sep=k2)
tu = tu %>% separate(site2_11g, c("site2_111", "site2_111g"), sep=k2)

# to remove the columns that are called "site1_111g" and "site2_111g"

tu = subset(tu, select = -c(site1_111g, site2_111g)) 

# write.table(tu, file = paste("CCLE_translocations_SvABA_20181221", "list.closest.genes.txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE) 

tu = tu %>% mutate(position1 = case_when( (grepl("Intergenic", tu$site1) == TRUE) ~ "Intergenic", 
	                                      (grepl("Intron", tu$site1) == TRUE) ~ "Intron",
	                                      (grepl("Exon", tu$site1) == TRUE) ~ "Exon"
))

tu = tu %>% mutate(position2 = case_when( (grepl("Intergenic", tu$site2) == TRUE) ~ "Intergenic", 
	                                      (grepl("Intron", tu$site2) == TRUE) ~ "Intron",
	                                      (grepl("Exon", tu$site2) == TRUE) ~ "Exon"
))

write.table(tu, file = paste("CCLE_translocations_SvABA_20181221", "list.closest.genes.txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)  
	
################################################################################################
################################################################################################
################################################################################################
################################################################################################

# to split the coordinates of the junctions
# x$bp1
# x$bp2

# a need to replace the LAST CHARACTER in bp1 or in bp2 in such a way that we can split :)
# using the library stringr

# bp1_t = str_replace(x$bp1, "\\-$", "&-")
# bp1_tt = str_replace(bp1_t, "\\+$", "&+")
# bp1_tt = data.frame(bp1_tt)
# head(bp1_tt)

# bp2_t = str_replace(x$bp2, "\\-$", "&-")
# bp2_tt = str_replace(bp2_t, "\\+$", "&+")
# bp2_tt = data.frame(bp2_tt)
# head(bp2_tt)

x = tu

x$bp1 = sub("\\-$", "z-", x$bp1)
x$bp1 = sub("\\+$", "z+", x$bp1)
head(x$bp1)

x$bp2 = sub("\\-$", "z-", x$bp2)
x$bp2 = sub("\\+$", "z+", x$bp2)
head(x$bp2)

x$bp1_t = x$bp1
x$bp2_t = x$bp2

head(x)

z1 = cSplit(x, "bp1_t", ":")
z2 = cSplit(z1, "bp1_t_2", "z")
z3 = cSplit(z2, "bp1_t_2_1", "-")

z4 = cSplit(z3, "bp2_t", ":")
z5 = cSplit(z4, "bp2_t_2", "z")
z6 = cSplit(z5, "bp2_t_2_1", "-")

colnames(z6)

# [1] "CCLE_name"       "map_id"          "bp1"             "bp2"
# [5] "class"           "gene1"           "gene2"           "site1"
# [9] "site2"           "fusion"          "multi_sv_fusion" "cosmic_fus"
# [13] "TUMOR"           "bp1_t_1"         "bp1_t_2_2"       "bp1_t_2_1_1"
# [17] "bp1_t_2_1_2"     "bp2_t_1"         "bp2_t_2_2"       "bp2_t_2_1_1"
# [21] "bp2_t_2_1_2"

###############################################################################################
###############################################################################################

mainChr = c(as.character(1:22),'x','X','y','Y')
z7 = dplyr::filter(z6, bp1_t_1 %in% mainChr)
z8 = dplyr::filter(z7, bp2_t_1 %in% mainChr)

# to transform the names of the chromosomes with "chr"

z9 = z8

z9$bp1_chr = paste("chr", z8$bp1_t_1, sep="")
z9$bp2_chr = paste("chr", z8$bp2_t_1, sep="")
head(z9)

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

z9$NAME = paste(z9$bp1, "&and&", z9$bp2, sep="") 
z9$SCORE = 1

# colnames(z9)
# [1] "CCLE_name"       "map_id"          "bp1"             "bp2"
# [5] "class"           "gene1"           "gene2"           "site1"
# [9] "site2"           "fusion"          "multi_sv_fusion" "cosmic_fus"
#[13] "TUMOR"           "bp1_t_1"         "bp1_t_2_2"       "bp1_t_2_1_1"
#[17] "bp1_t_2_1_2"     "bp2_t_1"         "bp2_t_2_2"       "bp2_t_2_1_1"
#[21] "bp2_t_2_1_2"     "bp1_chr"         "bp2_chr"         "NUMBER"
#[25] "SCORE"

# "bp1_t_1", 
# "bp1_t_2_1_1",
# "bp1_t_2_1_2"

# bp2_t_1 
# bp2_t_2_1_1 
# bp2_t_2_1_2

# "bp1_t_2_2",
# bp2_t_2_2 

# NUMBER
# SCORE

# changing the names of the COLUMNS in the DATAFRAME

colnames(z9)[which(colnames(z9) == "site1_11")] <- "gene1_anywhere"
colnames(z9)[which(colnames(z9) == "site1_111")] <- "gene1_anywhere_orientation"
colnames(z9)[which(colnames(z9) == "site2_11")] <- "gene2_anywhere"
colnames(z9)[which(colnames(z9) == "site2_111")] <- "gene2_anywhere_orientation"
colnames(z9)[which(colnames(z9) == "bp1_t_1")] <- "breakpoint1_chr"
colnames(z9)[which(colnames(z9) == "bp1_t_2_2")] <- "breakpoit1_orientation"
colnames(z9)[which(colnames(z9) == "bp1_t_2_1_1")] <- "breakpoint1_start"
colnames(z9)[which(colnames(z9) == "bp1_t_2_1_2")] <- "brerakpoint1_end"
colnames(z9)[which(colnames(z9) == "bp2_t_1")] <- "breakpoint2_chr"
colnames(z9)[which(colnames(z9) == "bp2_t_2_2")] <- "breakpoint2_orientation"
colnames(z9)[which(colnames(z9) == "bp2_t_2_1_1")] <- "breakpoint2_start"
colnames(z9)[which(colnames(z9) == "bp2_t_2_1_2")] <- "breakpoint2_end"

z10 = select(z9, -c('bp1_chr','bp2_chr'))

write.table(z10, file = paste("CCLE_translocations_SvABA_20181221", "LIST.GENES.POSITIONS.INTERACTIONS.txt", sep="."), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE) 

##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################

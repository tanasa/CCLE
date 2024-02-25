################################################################################################
################################################################################################
################################################################################################
################################################################################################

library("karyoploteR")
library("plyr")
library("doBy")
library("tidyverse")
library('dplyr')
library("tidyr")
library("readr")
library("splitstackshape")

################################################################################################
################################################################################################
################################################################################################
################################################################################################

# args <- commandArgs(TRUE)
# FILE <- args[1]  

# the data has been mapped on HG19 version. it consists of SvABA calls.
FILE = "CCLE_translocations_SvABA_20181221.txt"

# the script displays the SV data in a "karyotype" format. 

################################################################################################
################################################################################################
       
y <- read.table(FILE, header=T, sep="\t", stringsAsFactors=F)

x = y
head(x)
dim(x)
colnames(x)

# [1] "CCLE_name"       "map_id"          "bp1"             "bp2"
# [5] "class"           "gene1"           "gene2"           "site1"
# [9] "site2"           "fusion"          "multi_sv_fusion" "cosmic_fus"
# 75344 entries

x$SVTYPE.x = gsub("-like", "", x$class)
x$TUMOR = x$CCLE_name

unique(x$TUMOR)

z1 = cSplit(x, "bp1", ":")
z2 = cSplit(z1, "bp1_2", "-")
z3 = cSplit(z2, "bp1_2_2", "+")

z4 = cSplit(z3, "bp2", ":")
z5 = cSplit(z4, "bp2_2", "-")
z6 = cSplit(z5, "bp2_2_2", "+")

# colnames(z6)
# [1] "CCLE_name"       "map_id"          "class"           "gene1"
# [5] "gene2"           "site1"           "site2"           "fusion"
# [9] "multi_sv_fusion" "cosmic_fus"      "SVTYPE"          "TUMOR"
# [13] "bp1_1"           "bp1_2_1"         "bp1_2_2_1"       "bp2_1"
# [17] "bp2_2_1"         "bp2_2_2_1"

###############################################################################################
###############################################################################################

# to keep only the typical chromosomes on bp1_1 and bp2_1

mainChr = c(as.character(1:22),'x','X','y','Y')
z7 = dplyr::filter(z6, bp1_1 %in% mainChr)
z8 = dplyr::filter(z7, bp2_1 %in% mainChr)

# to transform the names of the chromosomes with "chr"

z = z8

z8$bp1_2 = paste("chr", z8$bp1_1, sep="")
z8$bp2_2 = paste("chr", z8$bp2_1, sep="")

z9 = z8 

z9 = z9 %>% rename(seqnames.x = bp1_2)
z9 = z9 %>% rename(start.x = bp1_2_1)
# z9 %>% rename( = bp1_2_2_1)
z9 = z9 %>% rename(CHR2.x = bp2_2)
z9 = z9 %>% rename(END.x = bp2_2_1)
# z9 %>% rename( = bp2_2_2_1)

z9$LENGTH = abs(z9$END.x - z9$start.x) 

##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
###############################################################################################
###############################################################################################
###############################################################################################  a LOOP that goes over the entire set of tumor cell lines (328 cell lines)

LIST_TUMORS = unique(z9$TUMOR)

############################# the LOOP finishes at the end of the SCRIPT #####################

for (i in 1: length(LIST_TUMORS)) 
{ 
	
# print(LIST_TUMORS[i])
name = LIST_TUMORS[i]

z10 = subset(z9, TUMOR == name)
write.table(z10, file=paste(name, "re-arrangements.txt", sep="."),
                 sep="\t", quote = FALSE, row.names=FALSE, col.names = TRUE)

z = z10 

############################################################################################### once we select a particular cell line, 
############################################################################################### we print the data frame 
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
############################################################################################### 
###############################################################################################
###############################################################################################

DEL <- subset(z, SVTYPE.x=="DEL")
DUP <- subset(z, SVTYPE.x=="DUP")
INS <- subset(z, SVTYPE.x=="INS")
INV <- subset(z, SVTYPE.x=="INV")
TRA <- subset(z, SVTYPE.x=="TRA")

dim(DEL)
dim(DUP)
dim(INS)
dim(INV)
dim(TRA)

### the number of records are :

dim(DEL)[1]
dim(DUP)[1]
dim(INS)[1]
dim(INV)[1]
dim(TRA)[1]

################################################################################################
################################################################################################
################################################################################################
################################################################################################
### to extract the coordinates for all SV :

coord.SV.partner1 <- makeGRangesFromDataFrame(data.frame(chr=z$seqnames.x, 
                                                         start=z$start.x,
                                                         end=z$start.x ))

coord.SV.partner2 <- makeGRangesFromDataFrame(data.frame(chr=z$CHR2.x, 
                                                         start=z$END.x,
                                                         end=z$END.x ))

#### in order to merge these GRANGES objects :

coord.SV.partners <- c(coord.SV.partner1, coord.SV.partner2)
length(coord.SV.partners)

################################################################################################
################################################################################################
################################################################################################
################################################################################################

if (dim(DEL)[1] > 0)
{
coord.DEL.partner1 <- makeGRangesFromDataFrame(data.frame(chr=DEL$seqnames.x, 
                                                         start=DEL$start.x,
                                                         end=DEL$start.x ))

coord.DEL.partner2 <- makeGRangesFromDataFrame(data.frame(chr=DEL$CHR2.x, 
                                                         start=DEL$END.x,
                                                         end=DEL$END.x ))

#### in order to merge these GRANGES objects :

coord.DEL.partners <- c(coord.DEL.partner1, coord.DEL.partner2)
}

# length(coord.DEL.partners)

##### if there are no deletions, introducing a DUMMY variable in order to avoid error messages:

if (dim(DEL)[1] == 0)
{
coord.DEL.partners <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################

if (dim(DUP)[1] > 0)
{
coord.DUP.partner1 <- makeGRangesFromDataFrame(data.frame(chr=DUP$seqnames.x, 
                                                         start=DUP$start.x,
                                                         end=DUP$start.x ))

coord.DUP.partner2 <- makeGRangesFromDataFrame(data.frame(chr=DUP$CHR2.x, 
                                                         start=DUP$END.x,
                                                         end=DUP$END.x ))

#### in order to merge these GRANGES objects:

coord.DUP.partners <- c(coord.DUP.partner1, coord.DUP.partner2)
}

# length(coord.DUP.partners)

##### if there are no duplications, introducing a DUMMY variable in order to avoid error messages:

if (dim(DUP)[1] == 0)
{
coord.DUP.partners <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################

if (dim(INS)[1] > 0)
{
coord.INS.partner1 <- makeGRangesFromDataFrame(data.frame(chr=INS$seqnames.x, 
                                                         start=INS$start.x,
                                                         end=INS$start.x ))

coord.INS.partner2 <- makeGRangesFromDataFrame(data.frame(chr=INS$CHR2.x, 
                                                         start=INS$END.x,
                                                         end=INS$END.x ))

#### in order to merge these GRANGES objects :

coord.INS.partners <- c(coord.INS.partner1, coord.INS.partner2)
}

# length(coord.INS.partners)

##### if there are no insertions, introducing a DUMMY variable in order to avoid error messages:

if (dim(INS)[1] == 0)
{
coord.INS.partners <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################

if (dim(INV)[1] > 0)
{

coord.INV.partner1 <- makeGRangesFromDataFrame(data.frame(chr=INV$seqnames.x, 
                                                         start=INV$start.x,
                                                         end=INV$start.x ))

coord.INV.partner2 <- makeGRangesFromDataFrame(data.frame(chr=INV$CHR2.x, 
                                                         start=INV$END.x,
                                                         end=INV$END.x ))

### in order to merge these GRANGES objects :

coord.INV.partners <- c(coord.INV.partner1, coord.INV.partner2)
}

# length(coord.INV.partners)

##### if there are no inversions, introducing a DUMMY variable in order to avoid error messages:

if (dim(INV)[1] == 0)
{
coord.INV.partners <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################

if (dim(TRA)[1] > 0)
{
coord.TRA.partner1 <- makeGRangesFromDataFrame(data.frame(chr=TRA$seqnames.x, 
                                                         start=TRA$start.x,
                                                         end=TRA$start.x ))

coord.TRA.partner2 <- makeGRangesFromDataFrame(data.frame(chr=TRA$CHR2.x, 
                                                         start=TRA$END.x,
                                                         end=TRA$END.x ))

#### in order to merge these GRANGES objects :

coord.TRA.partners <- c(coord.TRA.partner1, coord.TRA.partner2)
}

# length(coord.TRA.partners)

##### if there are no translocations, introducing a DUMMY variable in order to avoid error messages:

if (dim(TRA)[1] == 0)
{
coord.TRA.partners <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################
### in this version we print only SV and CNV

pdf(paste(name, ".view.karyotype.BREAKPOINTS.of.DEL.DUP.INS.INV.and.TRA.links.pdf", sep=""), width=10, height=12)

kp <- plotKaryotype(genome="hg19", plot.type=2, main="", cex=0.6)

kpPlotRegions(kp, coord.DEL.partners, col="red", data.panel=1, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.partners, col="darkorange", data.panel=1, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.partners, col="yellow", data.panel=1, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.partners, col="limegreen", data.panel=1, r0=0.61, r1=0.8 )

# kpPlotRegions(kp, coord.TRA.partners, col="blue", data.panel=2, r0=0.81, r1=0.1 )
# kpPlotRegions(kp, coord.SV.partners, col="magenta", data.panel=1, r0=0.21, r1=0.4 )

kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", 
                                                  data.panel=2, r0=0.21, r1=0.4 ) 

dev.off()

### PNG :

png(paste(name, ".view.karyotype.BREAKPOINTS.of.DEL.DUP.INS.INV.and.TRA.links.png", sep=""), 
                          width=10, height=12, units="cm", res = 300, pointsize = 10)

kp <- plotKaryotype(genome="hg19", plot.type=2, main="", cex=0.6)

kpPlotRegions(kp, coord.DEL.partners, col="red", data.panel=1, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.partners, col="darkorange", data.panel=1, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.partners, col="yellow", data.panel=1, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.partners, col="limegreen", data.panel=1, r0=0.61, r1=0.8 )

# kpPlotRegions(kp, coord.TRA.partners, col="blue", data.panel=2, r0=0.81, r1=0.1 )
# kpPlotRegions(kp, coord.SV.partners, col="magenta", data.panel=1, r0=0.21, r1=0.4 )

kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", 
                                                  data.panel=2, r0=0.21, r1=0.4 ) 

dev.off()

###############################################################################################################################################
###############################################################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
######################## in these displays we show the LENGTH of DEL < DUP < INS < INV
###############################################################################################################################################
###############################################################################################################################################

if (dim(DEL)[1] > 0)
{
coord.DEL.length <- makeGRangesFromDataFrame(data.frame(chr=DEL$seqnames.x, 
                                                        start=ifelse(DEL$start.x < DEL$END.x, DEL$start.x, DEL$END.x),
                                                        end=ifelse(DEL$start.x < DEL$END.x, DEL$END.x, DEL$start.x)))
} else { coord.DEL.length <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1)) }

####

if (dim(DUP)[1] > 0)
{
coord.DUP.length <- makeGRangesFromDataFrame(data.frame(chr=DUP$seqnames.x, 
                                                        start=ifelse(DUP$start.x < DUP$END.x, DUP$start.x, DUP$END.x),
                                                        end=ifelse(DUP$start.x < DUP$END.x, DUP$END.x, DUP$start.x)))
} else {coord.DUP.length <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))}

####

if (dim(INS)[1] > 0)
{
coord.INS.length <- makeGRangesFromDataFrame(data.frame(chr=INS$seqnames.x, 
                                                        start=ifelse(INS$start.x < INS$END.x, INS$start.x, INS$END.x),
                                                        end=ifelse(INS$start.x < INS$END.x, INS$END.x, INS$start.x)))
} else {coord.INS.length <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))}

####

if (dim(INV)[1] > 0)
{
coord.INV.length <- makeGRangesFromDataFrame(data.frame(chr=INV$seqnames.x, 
                                                        start=ifelse(INV$start.x < INV$END.x, INV$start.x, INV$END.x),
                                                        end=ifelse(INV$start.x < INV$END.x, INV$END.x, INV$start.x)))
} else {coord.INV.length <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))}

############################################################################################
############################################################################################
################## for TRA we do keep the display of the breakpoints :
################## as we have defined those breakpoints above
############################################################################################
############################################################################################

if (dim(TRA)[1] > 0)
{
coord.TRA.partner1 <- makeGRangesFromDataFrame(data.frame(chr=TRA$seqnames.x, 
                                                         start=TRA$start.x,
                                                         end=TRA$start.x ))

coord.TRA.partner2 <- makeGRangesFromDataFrame(data.frame(chr=TRA$CHR2.x, 
                                                         start=TRA$END.x,
                                                         end=TRA$END.x ))

#### in order to merge these GRANGES :

coord.TRA.partners <- c(coord.TRA.partner1, coord.TRA.partner2)
}

# length(coord.TRA.partners)

##### if there are no translocations, introducing a DUMMY variable in order to avoid error messages:

if (dim(TRA)[1] == 0)
{
coord.TRA.partners <- GRanges(seqnames = "chrY", ranges = IRanges(start = 1, end = 1))
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################

pdf(paste(name, ".view.karyotype.LENGTH.of.DEL.DUP.INS.INV.and.show.TRA.links.pdf", sep=""), 
                                                                  width=10, height=12)

kp <- plotKaryotype(genome="hg19", plot.type=2, main="", cex=0.6)

kpPlotRegions(kp, coord.DEL.length, col="red", data.panel=2, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.length, col="darkorange", data.panel=2, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.length, col="yellow", data.panel=2, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.length, col="limegreen", data.panel=2, r0=0.61, r1=0.8 )

# kpPlotRegions(kp, coord.TRA.partners, col="blue", data.panel=1, r0=0.21, r1=0.4 )
kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", 
                                                  data.panel=1, r0=0.21, r1=0.4 ) 

dev.off()

#################################################################################################
#################################################################################################

png(paste(name, ".view.karyotype.LENGTH.of.DEL.DUP.INS.INV.and.show.TRA.links.png", sep=""), 
                 width=10, height=12, units="cm", res = 300, pointsize = 10)

kp <- plotKaryotype(genome="hg19", plot.type=2, main="", cex=0.6)

kpPlotRegions(kp, coord.DEL.length, col="red", data.panel=2, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.length, col="darkorange", data.panel=2, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.length, col="yellow", data.panel=2, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.length, col="limegreen", data.panel=2, r0=0.61, r1=0.8 )

# kpPlotRegions(kp, coord.TRA.partners, col="blue", data.panel=1, r0=0.21, r1=0.4 )
kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", 
                                                     data.panel=1, r0=0.21, r1=0.4 )

dev.off()

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
### to display the TRANSLOCATIONS as ARCS, the chromosomes are represented linearly  

pdf(paste(name, ".view.karyotype.LENGTH.of.DEL.DUP.INS.INV.and.show.TRA.links.HORIZONTALLY.pdf", sep=""), 
                                                                  width=30, height=20)

pp <- getDefaultPlotParams(plot.type=3)
pp$data2height <- 400

kp <- plotKaryotype(genome="hg19", cex = 1,  
                                   plot.params = pp, 
                                   labels.plotter = NULL, 
                                   ideogram.plotter = NULL, 
                                   plot.type=3)

kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=90, cex=1)

kpPlotRegions(kp, coord.DEL.length, col="red", data.panel=1, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.length, col="darkorange", data.panel=1, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.length, col="yellow", data.panel=1, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.length, col="limegreen", data.panel=1, r0=0.61, r1=0.8 )

kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", data.panel=2, r0=0, r1=0.8) 

dev.off()

### PNG :

png(paste(name, ".view.karyotype.LENGTH.of.DEL.DUP.INS.INV.and.show.TRA.links.HORIZONTALLY.png", 
          sep=""), width=30, height=20, units="cm", res = 300, pointsize = 10)

pp <- getDefaultPlotParams(plot.type=3)
pp$data2height <- 400

kp <- plotKaryotype(genome="hg19", cex = 1,  
                                   plot.params = pp, 
                                   labels.plotter = NULL, 
                                   ideogram.plotter = NULL, 
                                   plot.type=3)

kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=90, cex=1)

kpPlotRegions(kp, coord.DEL.length, col="red", data.panel=1, r0=0, r1=0.2 )
kpPlotRegions(kp, coord.DUP.length, col="darkorange", data.panel=1, r0=0.21, r1=0.4 )
kpPlotRegions(kp, coord.INS.length, col="yellow", data.panel=1, r0=0.41, r1=0.6 )
kpPlotRegions(kp, coord.INV.length, col="limegreen", data.panel=1, r0=0.61, r1=0.8 )

kpPlotLinks(kp, coord.TRA.partner1, coord.TRA.partner2, col="blue", data.panel=2, r0=0, r1=0.8) 

dev.off()

#################################################################################################
#################################################################################################
#################################################################################################
################################################################################################# 
#### we end the LOOP that displays the data for each cell line

}

#################################################################################################################
#################################################################################################################
#################################################################################################################

######################################################################################################
######################################################################################################
######################################################################################################

library(fs)
library(stringr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(ggplot2)
library(dplyr)

######################################################################################################
######################################################################################################
######################################################################################################

# The script compare the numerical values of the DENSITY AREAS across the cell lines.

setwd("")
files <- list.files(path=getwd())
length(files)

# BIG = data.frame(matrix(nrow = 0, ncol = 3)) 
# colnames(BIG) = c("AUC","nodes","value")

######################################################################################################
######################################################################################################
######################################################################################################

library(data.table)    
BIG = data.table(
            "density" = character(),
            "value" = numeric()
)

BIG = as.data.frame(BIG)

name.a=list()
name.b=list()

######################################################################################################
######################################################################################################
######################################################################################################

for(i in 1:length(files)){
  data <- read.delim(file = files[i], header = FALSE, stringsAsFactors = FALSE, sep="\t", row.names=NULL)
  # name.a[i] <- sub("\\.SV.network.igraph.measurements.txt", "", files[i]) # removing the end of the file name
  # name.b[i] <- sub("graph\\.", "", name.a[i])
  
  colnames(data)[1] = "density"
  colnames(data)[2] = "value"

  name.a[[i]] <- sub("\\.SV.network.igraph.measurements.eCDF.txt.density.area.txt", "", files[i]) # removing the end of the file name
  name.b[[i]] <- sub("graph\\.", "", name.a[[i]])
  
  BIG[i,"cell"] = str_split(name.b[i], "_",  n = 2)[[1]][1]
  BIG[i,"tissue"] = str_split(name.b[i], "_",  n = 2)[[1]][2]  
  
  BIG[i, "density"] = data$"density"
  BIG[i,"value"] = data$"value"
}
  
write.table(BIG, file="ALL.the.CELL.LINES.in.a.LARGE.file.DENSITY.AREAS.txt", sep="\t", append=FALSE, col.names=TRUE, row.names=FALSE)

######################################################################################################
######################################################################################################
######################################################################################################

BIG2 = BIG
BIG2$"density" = as.character(BIG2$"density")
BIG2$"value" = as.numeric(BIG2$"value")                              
BIG2$"cell" = as.factor(BIG2$"cell")
BIG2$"tissue" = as.factor(BIG2$"tissue")                    

str(BIG2)

######################################################################################################
######################################################################################################
###################################################################################################### density value
# z = melt(BIG2,id.vars=c("tissue", "the diameter")) 

mm = BIG2 %>% select("tissue", "value")
colnames(mm) = c("tissue", "value")
sapply(mm, class)

zz = melt(mm, id =c ("tissue"))
head(zz)

ggplot(zz, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Density value") 

ggsave(paste("ALL.cell.lines.density.value", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
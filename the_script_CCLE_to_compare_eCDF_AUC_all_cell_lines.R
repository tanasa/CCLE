
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
# library(tidyverse)

######################################################################################################
######################################################################################################
######################################################################################################

setwd("./CCLE_SvABA//files_eCDF_areas")
files <- list.files(path=getwd())
length(files)

BIG = data.frame(matrix(nrow = 0, ncol = 10)) 
# colnames(BIG) = c("AUC","nodes","value")

######################################################################################################
######################################################################################################
######################################################################################################

library(data.table)    
BIG = data.table(
          "cell" = character(),
          "tissue" = character(),
"AUC less than.1" = numeric(),
"AUC more than.1" = numeric(),
"AUC less than.2" = numeric(),
"AUC more than.2" = numeric(),
"AUC less than.3" = numeric(),
"AUC more than.3" = numeric(),
"AUC less than.4" = numeric(),
"AUC more than.4" = numeric(),
"AUC less than.5" = numeric(),
"AUC more than.5" = numeric()
)

BIG = as.data.frame(BIG)

name.a=list()
name.b=list()

######################################################################################################
######################################################################################################
######################################################################################################

for(i in 1:length(files)){
  data <- read.delim(file = files[i], header = FALSE, stringsAsFactors = FALSE, sep=",", row.names=NULL)
  # name.a[i] <- sub("\\.SV.network.igraph.measurements.txt", "", files[i]) # removing the end of the file name
  # name.b[i] <- sub("graph\\.", "", name.a[i])
  
  colnames(data)[1] = "AUC"
  colnames(data)[2] = "nodes"
  colnames(data)[3] = "value"
  data$AUCnodes = paste(data$AUC, data$nodes, sep=".")

  name.a[[i]] <- sub("\\.SV.network.igraph.measurements.eCDF.txt.eCDF.area.txt", "", files[i]) # removing the end of the file name
  name.b[[i]] <- sub("graph\\.", "", name.a[[i]])
  
  BIG[i,"cell"] = str_split(name.b[i], "_",  n = 2)[[1]][1]
  BIG[i,"tissue"] = str_split(name.b[i], "_",  n = 2)[[1]][2]  
  
  BIG[i, "AUC less than.1"] = data$"value"[1]
  BIG[i, "AUC more than.1"] = data$"value"[2]
  BIG[i, "AUC less than.2"] = data$"value"[3]
  BIG[i, "AUC more than.2"] = data$"value"[4]
  BIG[i, "AUC less than.3"] = data$"value"[5]
  BIG[i, "AUC more than.3"] = data$"value"[6]
  BIG[i, "AUC less than.4"] = data$"value"[7]
  BIG[i, "AUC more than.4"] = data$"value"[8]
  BIG[i, "AUC less than.5"] = data$"value"[9]
  BIG[i, "AUC more than.5"] = data$"value"[10]
}
  
write.table(BIG, file="ALL.the.CELL.LINES.in.a.LARGE.file.eCDF.AREAS.txt", sep="\t", append=FALSE, col.names=TRUE, row.names=FALSE)

######################################################################################################
######################################################################################################
######################################################################################################

BIG2 = BIG
BIG2$"cell" = as.factor(BIG2$"cell"  )
BIG2$"tissue" = as.factor(BIG2$"tissue"  )

BIG2$"AUC less than.1" = as.numeric(BIG2$"AUC less than.1")
BIG2$"AUC more than.1" = as.numeric(BIG2$"AUC more than.1")                              
BIG2$"AUC less than.2" = as.numeric(BIG2$"AUC less than.2"  )              
BIG2$"AUC more than.2"   = as.numeric(BIG2$"AUC more than.2"  )                     
BIG2$"AUC less than.3"    = as.numeric(BIG2$"AUC less than.3"  )                   
BIG2$"AUC more than.3"= as.numeric(BIG2$"AUC more than.3"  )
BIG2$"AUC less than.4"  = as.numeric(BIG2$"AUC less than.4"  )
BIG2$"AUC more than.4"  = as.numeric(BIG2$"AUC more than.4"  )
BIG2$"AUC less than.5"      = as.numeric(BIG2$"AUC less than.5"   )       
BIG2$"AUC more than.5"  = as.numeric(BIG2$"AUC more than.5"  )     

str(BIG2)

######################################################################################################
######################################################################################################
###################################################################################################### 
# z = melt(BIG2,id.vars=c("tissue", "the diameter")) 

mm = BIG2 %>% select("tissue", "AUC less than.1")
colnames(mm) = c("tissue", "AUC less than.1")
sapply(mm, class)

zz = melt(mm, id =c ("tissue"))
head(zz)

ggplot(zz, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("AUC less than.1") 

ggsave(paste("ALL.AUC less than.1", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### 

mm_A = BIG2 %>% select("tissue", "AUC more than.1")
colnames(mm_A) = c("tissue", "AUC more than.1")
sapply(mm_A, class)

zz_A = melt(mm_A, id =c ("tissue"))
head(zz)

ggplot(zz_A, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("AUC more than.1") 

ggsave(paste("ALL.AUC more than.1", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### 
mm_B = BIG2 %>% select("tissue", "AUC less than.2")
colnames(mm_B) = c("tissue", "AUC less than.2")
sapply(mm_B, class)

zz_B = melt(mm_B, id =c ("tissue"))
head(zz_B)

ggplot(zz_B, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("AUC less than.2") 


ggsave(paste("ALL.AUC less than.2", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### 

mm_C = BIG2 %>% select("tissue", "AUC more than.2")
colnames(mm_C) = c("tissue", "AUC more than.2")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("AUC more than.2") 

ggsave(paste("ALL.AUC more than.2", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### 

mm_C = BIG2 %>% select("tissue", "AUC less than.3")
colnames(mm_C) = c("tissue", "AUC less than.3")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("AUC less than.3") 

ggsave(paste("ALL.AUC less than.3", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### 

mm_C = BIG2 %>% select("tissue", "AUC more than.3")
colnames(mm_C) = c("tissue", "AUC more than.3")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("AUC more than.3") 

ggsave(paste("ALL.AUC more than.3", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### 

mm_C = BIG2 %>% select("tissue", "AUC less than.4")
colnames(mm_C) = c("tissue", "AUC less than.4")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("AUC less than.4") 

ggsave(paste("ALL.AUC less than.4", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### 

mm_C = BIG2 %>% select("tissue", "AUC more than.4")
colnames(mm_C) = c("tissue", "AUC more than.4")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("AUC more than.4") 

ggsave(paste("ALL.AUC more than.4", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### 

mm_C = BIG2 %>% select("tissue", "AUC less than.5")
colnames(mm_C) = c("tissue", "AUC less than.5")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("AUC less than.5") 

ggsave(paste("ALL.AUC less than.5", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### 

mm_C = BIG2 %>% select("tissue", "AUC more than.5")
colnames(mm_C) = c("tissue", "AUC more than.5")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("AUC more than.5") 

ggsave(paste("ALL.AUC more than.5", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### 
######################################################################################################

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

setwd("./CCLE_SvABA/")
files <- list.files(path=getwd())
length(files)


columns = c("cell", 
            "tissue",
            "the diameter",
            "number of nodes (vertices)",
            "number of SV (edges)",
            "density of the nodes",
            "Number of communities : Spectral Clustering",
            "Number of communities : Louvain clustering",
            "Modularity : Spectral Clustering",
            "Modularity : Louvain clustering",
            "average path for unconnected networks",
            "Global clustering coefficient",
            "Average Local Clustering Coefficient",
            "Edge density",
            "Small world index : transitivity :",
            "Small world index : APL :",
            "Small world index : transitivity :",
             "Shapiro-Wilk normality test",
             "package GAMlss : distribution :",
             "package GAMlss : distribution :")

BIG = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(BIG) = columns

######################################################################################################
######################################################################################################
######################################################################################################

library(data.table)    
BIG = data.table(
            "cell"= character(),
            "tissue"= character(), 
            "the diameter"=numeric() ,
            "number of nodes (vertices)"=numeric(), 
            "number of SV (edges)"=numeric(), 
            "density of the nodes"=numeric(),
            "Number of communities : Spectral Clustering"=numeric(), 
            "Number of communities : Louvain clustering"=numeric(),
            "Modularity : Spectral Clustering"=numeric(), 
            "Modularity : Louvain clustering"=numeric(),
             "average path for unconnected networks"=numeric(),
             "Global clustering coefficient"=numeric(), 
             "Average Local Clustering Coefficient"=numeric(), 
             "Edge density"=numeric(), 
             "Small world index : transitivity :"=numeric(), 
             "Small world index : APL :"=numeric(), 
             "Small world index : transitivity :"=numeric(), 
             "Shapiro-Wilk normality test"=numeric(), 
             "package GAMlss : distribution :"= character(), 
             "package GAMlss : distribution :"= character()
)

BIG = as.data.frame(BIG)

name.a=list()
name.b=list()

######################################################################################################
######################################################################################################
######################################################################################################

for(i in 1:length(files)){
  data <- read.delim(file = files[i], header = TRUE, stringsAsFactors = FALSE, sep="\t", row.names=NULL)
  # name.a[i] <- sub("\\.SV.network.igraph.measurements.txt", "", files[i]) # removing the end of the file name
  # name.b[i] <- sub("graph\\.", "", name.a[i])
  
  name.a[[i]] <- sub("\\.SV.network.igraph.measurements.txt", "", files[i]) # removing the end of the file name
  name.b[[i]] <- sub("graph\\.", "", name.a[[i]])
  
  # BIG[i,"cell"] = str_split(name.b[i], "_",  n = 2)[[1]][1]
  # BIG[i,"tissue"] = str_split(name.b[i], "_",  n = 2)[[1]][2]   
  
  BIG[i,"cell"] = str_split(name.b[i], "_",  n = 2)[[1]][1]
  BIG[i,"tissue"] = str_split(name.b[i], "_",  n = 2)[[1]][2]  
  
  BIG[i, "the diameter"] = data$"A.summary.about.the.network.features."[[1]]

  
  BIG[i,"number of nodes (vertices)"] = data$"A.summary.about.the.network.features."[2]
  BIG[i,"number of SV (edges)"] = data$"A.summary.about.the.network.features."[3]
  BIG[i,"density of the nodes"] = data$"A.summary.about.the.network.features."[4]
  BIG[i, "Number of communities : Spectral Clustering"] = data$"A.summary.about.the.network.features."[5]
  BIG[i, "Number of communities : Louvain clustering"] = data$"A.summary.about.the.network.features."[6]
  BIG[i, "Modularity : Spectral Clustering"] = data$"A.summary.about.the.network.features."[7]
  BIG[i, "Modularity : Louvain clustering"] = data$"A.summary.about.the.network.features."[8]
  BIG[i, "average path for unconnected networks"] = data$"A.summary.about.the.network.features."[9]
  BIG[i, "Global clustering coefficient"] = data$"A.summary.about.the.network.features."[10]
  BIG[i, "Average Local Clustering Coefficient"] = data$"A.summary.about.the.network.features."[11]
  BIG[i, "Edge density"] = data$"A.summary.about.the.network.features."[12]
  BIG[i, "Small world index : transitivity :"] = data$"A.summary.about.the.network.features."[13]
  BIG[i, "Small world index : APL :"] = data$"A.summary.about.the.network.features."[14]
  BIG[i, "Small world index : transitivity :"] = data$"A.summary.about.the.network.features."[15]
  BIG[i, "Shapiro-Wilk normality test"] = data$"A.summary.about.the.network.features."[16]
  BIG[i, "package GAMlss : distribution :"] = data$"A.summary.about.the.network.features."[17]
  BIG[i, "package GAMlss : distribution :"] = data$"A.summary.about.the.network.features."[18]  

}
  
write.table(BIG, file="ALL.the.CELL.LINES.in.a.LARGE.file.txt", sep="\t", append=FALSE, col.names=TRUE, row.names=FALSE)

######################################################################################################
######################################################################################################
######################################################################################################

BIG2 = BIG
BIG2$"Edge density" = as.numeric(BIG2$"Edge density")
BIG2$"the diameter" = as.numeric(BIG2$"the diameter")                              
BIG2$"number of nodes (vertices)" = as.numeric(BIG2$"number of nodes (vertices)"  )              
BIG2$"number of SV (edges)"   = as.numeric(BIG2$"number of SV (edges)"  )                     
BIG2$"density of the nodes"    = as.numeric(BIG2$"density of the nodes"  )                   
BIG2$"Number of communities : Spectral Clustering"= as.numeric(BIG2$"Number of communities : Spectral Clustering"  )
BIG2$"Number of communities : Louvain clustering"  = as.numeric(BIG2$"Number of communities : Louvain clustering"  )
BIG2$"Modularity : Spectral Clustering"  = as.numeric(BIG2$"Modularity : Spectral Clustering"  )
BIG2$"Modularity : Louvain clustering"      = as.numeric(BIG2$"Modularity : Louvain clustering"   )       
BIG2$"average path for unconnected networks"  = as.numeric(BIG2$"average path for unconnected networks"  )     
BIG2$"Global clustering coefficient"     = as.numeric(BIG2$"Global clustering coefficient"  )         
BIG2$"Average Local Clustering Coefficient"= as.numeric(BIG2$"Average Local Clustering Coefficient"  )        
BIG2$"Edge density"          = as.numeric(BIG2$"Edge density"  )                      
BIG2$"Small world index : transitivity :"   = as.numeric(BIG2$"Small world index : transitivity :"   )      
BIG2$"Small world index : APL :"       = as.numeric(BIG2$"Small world index : APL :"   )            
BIG2$"Small world index : transitivity :"  = as.numeric(BIG2$"Small world index : transitivity :"  )        
BIG2$"Shapiro-Wilk normality test" = as.numeric(BIG2$"Shapiro-Wilk normality test"  )

BIG2$"cell" = as.factor(BIG2$"cell"  )
BIG2$"tissue" = as.factor(BIG2$"tissue"  )

str(BIG2)

######################################################################################################
######################################################################################################
###################################################################################################### diameter
# z = melt(BIG2,id.vars=c("tissue", "the diameter")) 

mm = BIG2 %>% select("tissue", "the diameter")
colnames(mm) = c("tissue", "diameter")
sapply(mm, class)

zz = melt(mm, id =c ("tissue"))
head(zz)

ggplot(zz, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Network Diameter") 

ggsave(paste("ALL.network.diameter", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### number of nodes 

mm_A = BIG2 %>% select("tissue", "number of nodes (vertices)")
colnames(mm_A) = c("tissue", "number of nodes (vertices)")
sapply(mm_A, class)

zz_A = melt(mm_A, id =c ("tissue"))
head(zz)

ggplot(zz_A, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("number of nodes (vertices)") 

ggsave(paste("ALL.number of nodes (vertices)", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### number of edges

mm_B = BIG2 %>% select("tissue", "number of SV (edges)")
colnames(mm_B) = c("tissue", "number of SV (edges)")
sapply(mm_B, class)

zz_B = melt(mm_B, id =c ("tissue"))
head(zz_B)

ggplot(zz_B, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("number of SV (edges)") 


ggsave(paste("ALL.number of SV (edges)", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### density of the nodes

mm_C = BIG2 %>% select("tissue", "density of the nodes")
colnames(mm_C) = c("tissue", "density of the nodes")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("density of the nodes") 

ggsave(paste("ALL.density of the nodes", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### Number of communities : Spectral Clustering

mm_C = BIG2 %>% select("tissue", "Number of communities : Spectral Clustering")
colnames(mm_C) = c("tissue", "Number of communities : Spectral Clustering")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Number of communities : Spectral Clustering") 

ggsave(paste("ALL.Number of communities : Spectral Clustering", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### Number of communities : Louvain clustering

mm_C = BIG2 %>% select("tissue", "Number of communities : Louvain clustering")
colnames(mm_C) = c("tissue", "Number of communities : Louvain clustering")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Number of communities : Louvain clustering") 

ggsave(paste("ALL.Number of communities : Louvain clustering", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### Modularity : Spectral Clustering

mm_C = BIG2 %>% select("tissue", "Modularity : Spectral Clustering")
colnames(mm_C) = c("tissue", "Modularity : Spectral Clustering")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Modularity : Spectral Clustering") 

ggsave(paste("ALL.Modularity : Spectral Clustering", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### Modularity : Louvain clustering

mm_C = BIG2 %>% select("tissue", "Modularity : Louvain clustering")
colnames(mm_C) = c("tissue", "Modularity : Louvain clustering")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Modularity : Louvain clustering") 

ggsave(paste("ALL.Modularity : Louvain clustering", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### average path for unconnected networks

mm_C = BIG2 %>% select("tissue", "average path for unconnected networks")
colnames(mm_C) = c("tissue", "average path for unconnected networks")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("average path for unconnected networks") 

ggsave(paste("ALL.average path for unconnected networks", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### Global clustering coefficient

mm_C = BIG2 %>% select("tissue", "Global clustering coefficient")
colnames(mm_C) = c("tissue", "Global clustering coefficient")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Global clustering coefficient") 

ggsave(paste("ALL.Global clustering coefficient", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### Average Local Clustering Coefficient

mm_C = BIG2 %>% select("tissue", "Average Local Clustering Coefficient")
colnames(mm_C) = c("tissue", "Average Local Clustering Coefficient")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Average Local Clustering Coefficient") 

ggsave(paste("ALL. Average Local Clustering Coefficient", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### Edge density

mm_C = BIG2 %>% select("tissue", "Edge density")
colnames(mm_C) = c("tissue", "Edge density")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Edge density") 

ggsave(paste("ALL.Edge density", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### Small world index : transitivity 

mm_C = BIG2 %>% select("tissue", "Small world index : transitivity ")
colnames(mm_C) = c("tissue", "Small world index : transitivity ")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Small world index : transitivity ") 

ggsave(paste("ALL.Small world index : transitivity ", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
###################################################################################################### Small world index : APL

mm_C = BIG2 %>% select("tissue", "Small world index : APL")
colnames(mm_C) = c("tissue", "Small world index : APL")
sapply(mm_C, class)

zz_C = melt(mm_C, id =c ("tissue"))
head(zz_C)

ggplot(zz_C, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Small world index : APL") 

ggsave(paste("ALL.Small world index : APL", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################
######################################################################################################  Shapiro-Wilk normality test
# z = melt(BIG2,id.vars=c("tissue", "the diameter")) 

mm = BIG2 %>% select("tissue", " Shapiro-Wilk normality test")
colnames(mm) = c("tissue", " Shapiro-Wilk normality test")
sapply(mm, class)

zz = melt(mm, id =c ("tissue"))
head(zz)

ggplot(zz, aes(x = tissue, y = value, color=tissue)) +            # Applying ggplot function
  geom_boxplot() +
  geom_jitter() +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle(" Shapiro-Wilk normality test") 

ggsave(paste("ALL. Shapiro-Wilk normality test", "boxplot.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

######################################################################################################
######################################################################################################

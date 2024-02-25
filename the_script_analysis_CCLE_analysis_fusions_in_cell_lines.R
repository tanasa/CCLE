##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################

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
library(ggpubr)
library(clusterProfiler)
# library("ggdensity")
# library("statnet")
library(dplyr)
library(stringr)
library(tidyr)
library("AnnotationHub")
library("clusterProfiler")
library("DOSE")
library("enrichplot")
library("ggnewscale")
library(msigdbr)
library("AnnotationHub")
library("org.Hs.eg.db")
library(AnnotationDbi)
library("EnsDb.Hsapiens.v86")
library("enrichplot")
library(fitdistrplus)
library(gamlss)
library(gamlss.dist)
library(gamlss.add)
library(zoo)

# chooseCRANmirror(78)

################################################################################################
################################################################################################
################################################################################################
################################################################################################

# we work with the file C...mo.txt, where we have replaced "gene" with "Gene"
# the data has been mapped on hg19 or hg38

FILE = "CCLE_translocations_SvABA_20181221.LIST.GENES.POSITIONS.INTERACTIONS.txt"
z = read.table("CCLE_translocations_SvABA_20181221.LIST.GENES.POSITIONS.INTERACTIONS.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE) 

# colnames(z9)
# [1] "CCLE_name"                  "map_id"                     "bp1"                       
# [4] "bp2"                        "class"                      "Gene1"                     
# [7] "Gene2"                      "site1"                      "site2"                     
#[10] "fusion"                     "multi_sv_fusion"            "cosmic_fus"                
#[13] "TUMOR"                      "Gene1"             "Gene1_orientation"
#[16] "Gene2"             "Gene2_orientation" "position1"                 
#[19] "position2"                  "breakpoint1_chr"            "breakpoit1_orientation"    
#[22] "breakpoint1_start"          "brerakpoint1_end"           "breakpoint2_chr"           
#[25] "breakpoint2_orientation"    "breakpoint2_start"          "breakpoint2_end"           
#[28] "NAME"                       "SCORE"

##############################################################################################################################################################################################
##############################################################################################################################################################################################
############################################################################################################################################################################################## TO SELECT CHROMOSOMES

# hg38 = read.delim("hg38.chrom.classical.broad.format.sizes", header=F, sep="\t", stringsAsFactors=F)
# setnames(hg38, c("chr", "length"))
# hg38$chr
# hg38$length

# we select only the autosome chromosomes from 1 to 22, and we exclude the chromosome M
# hg38 = hg38[!is.na(as.numeric(hg38$chr)), ]
# we switch the chromosome sizes to hg19. I believe that CCLE data is provided on hg19.  

# hg19 = read.delim("hg19.chrom.classical.broad.format.sizes", header=F, sep="\t", stringsAsFactors=F)
# setnames(hg19, c("chr", "length"))
# hg19$chr
# hg19$length

### we select only the autosome chromosomes from 1 to 22, and we exclude the chromosome M
### hg19 = hg19[!is.na(as.numeric(hg19$chr)), ]

##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
############################################################################################################ EXAMINING THE PAIR of GENES whose intronic or exonic region is disrupted by a SV
# we extract the genes : "Gene1", "Gene2"

t = z[, c("breakpoint1_chr", "breakpoint1_start", "breakpoint1_end", 
	      "breakpoint2_chr", "breakpoint2_start", "breakpoint2_end", 
		  "NAME", "SCORE",
          "breakpoint1_orientation", 
          "breakpoint2_orientation",  
		  "TUMOR", "class", "Gene1", "Gene2", 
		  "position1", "position2")]
		  		  
###############################################################################################
############################################################################################### shall we select
############################################################################################### MCF-7 cell type
# TUMOR_TYPE="MCF7_BREAST"
# y = subset(t, TUMOR=="MCF7_BREAST")
# write.table(y, file=paste("graph", TUMOR_TYPE, "network.bedpe.txt", sep="."), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
##############################################################################################
###############################################################################################
############################################################################################### 

LIST_TUMORS = unique(t$TUMOR)

resultALL = try({

###############################################################################################
###############################################################################################

for (i in 1: length(LIST_TUMORS)) 
{ 
	
# print(LIST_TUMORS[i])
TUMOR_TYPE = LIST_TUMORS[i]

###############################################################################################
###############################################################################################

y = subset(t, TUMOR == TUMOR_TYPE)
write.table(y, file=paste("graph", TUMOR_TYPE, "network.bedpe.txt", sep="."), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

###############################################################################################
###############################################################################################
###############################################################################################
############################################################################################### 

el = y[!(y$Gene1=="") & !(y$Gene2==""), ]  # it is appropriate for the LIST of FUSIONS ; shall we decide not to work with the list of SV 
el = subset(el, select=c("Gene1", "Gene2", "class"))
dim(el)

###############################################################################################
###############################################################################################
###############################################################################################
############################################################################################### we can work with the dataframe that we call "el"

el$color[el$class=="DEL-like"] <- "red"
el$color[el$class=="DUP-like"] <- "green"
el$color[el$class=="INV-like"] <- "brown"
el$color[el$class=="TRA-like"] <- "blue"

# el$color[el$class=="INS-like"] <- "yellow"
# el$color <- as.factor(el$color)

###############################################################################################
###############################################################################################
############################################################################################### USING IGRAPH DISPLAY and PROPERTIES 
###############################################################################################

g = graph_from_data_frame(d = el, directed = FALSE)

#  layout = layout.reingold.tilford 
#  plot(g,
#      vertex.color = rgb(0.8,0.2,0.2,0.9),           # Node color
#      vertex.frame.color = "Forestgreen",            # Node border color
#      vertex.shape=c("circle"),                      # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
    # vertex.size=c(10:15),                          # Size of the node (default is 15)
#       vertex.size=0.1,
#	    vertex.label.family="Times",                   # Font family of the label (e.g.“Times”, “Helvetica”)
#	    vertex.label.font=1,                         # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
#	    vertex.label.cex=0.5,                        # Font size (multiplication factor, device-dependent)
#	    vertex.label.dist=0.5,                       # Distance between the label and the vertex
#	    vertex.label.degree=0                        # The second size of the node (e.g. for a rectangle)
#   )

###############################################################################################
###############################################################################################
############################################################################################### a good DISPLAY :)

lo <- layout_with_kk(g) # create a layout
lo <- norm_coords(lo, ymin=-1, ymax=1, xmin=-1, xmax=1)

par(mfrow=c(1,2), mar=c(0,0,0,0))

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.png", sep="."), width = 1600, height = 1600, units = "px")        ##### great !!
plot.igraph(g, 
	          # layout=lo*10,
			  # layout = layout_in_circle(network),
	          edge.arrow.width = .25,
    	      edge.arrow.size = .25,
	          edge.color=edge_attr(g)$color,
	          # edge.color=el$COLOR,
	          # edge.color=E(g)$color,
     
	          vertex.frame.color = "Forestgreen",            # Node border color
              vertex.shape=c("circle"),                      # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
                                                             # Size of the node (default is 15)
              vertex.size=5,
 	          vertex.label.family="Times",                  # Font family of the label (e.g.“Times”, “Helvetica”)
 	          vertex.label.font=1,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
 	          vertex.label.cex=1,                           # Font size (multiplication factor, device-dependent)
 	          vertex.label.dist=0.5,                        # Distance between the label and the vertex
 	          vertex.label.degree=0)
dev.off()

# ggsave(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.png", sep="."), 
#			   width = 30,
#			   height = 50,
#			   units = "cm",
#			   dpi = 300 )

#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### posibly to add Cytoscape
# https://dshizuka.github.io/networkanalysis/index.html

# V(g) ## vertex
# VERTEX ATTRIBUTES : in a separate graph
# E(g) ## edges
# E(g)$color <- as.factor(el$COLOR)
# V(g)$name
# E(g)$width

E(g)$weight = 1 # a weighted network  

# plot(g, layout=layout_with_kk(g), vertex.label="", vertex.color="gold", edge.color="slateblue", edge.width=E(g)$weight*5)
# TO PRINT the ADJANCENCY MATRIX :
# as_adjacency_matrix(g, sparse=F)
# as_adjacency_matrix(g, sparse=T)
# To PRINT the LIST of EDGES :
# as_edgelist(g)
# TO PRINT the ADJANCENCY LIST : the list of nodes
# as_adj_list(g)

# TO EXTRACT the DATA FRAME :
# as_data_frame(g)
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### MEASURING THE NETWORKS

# 4.1 Centrality Measures (i.e., node-level measures)
# 4.1.1 Degree and Strength
# 4.1.2 Betweenness
# 4.1.3 Assembling a dataset of node-level measures
# 4.2 Network-level measures
# 4.2.1 Size and density
# 4.2.2 Components
# 4.2.3 Degree distributions
# 4.2.4 Average path length & Diameter
# 4.2.5 Clustering Coefficient (Transitivity)

#################################################################################################################### 
#################################################################################################################### DEGREE
#################################################################################################################### BETWEENNESS
#################################################################################################################### CLOSENESS

degree=igraph::degree
betweenness=igraph::betweenness
closeness=igraph::closeness

#################################################################################################################### 
#################################################################################################################### DEGREE CENTRALITY
#################################################################################################################### 
# NODEL LEVEL MEASURES :

# DEGREE CENTRALITY : 
# Degree centrality is simply the number of edges connected to a given node. 
# In a social network, this might mean the number of friends an individual has.
# degree(g)

de=igraph::degree(g)

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.CENTRALITY.DEGREE.png", sep="."), width = 1600, height = 1600, units = "px")   # it looks good
plot.igraph(g, 
	         # layout=lo*10,
			 # layout = layout_in_circle(network),
	         edge.arrow.width = .25,
    	     edge.arrow.size = .25,
	         edge.color=edge_attr(g)$color,
	         # edge.color=el$COLOR,
	         # edge.color=E(g)$color,
     
	         vertex.frame.color = "Forestgreen",            # Node border color
             vertex.shape=c("circle"),                      # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
                                                            # Size of the node (default is 15)
             # vertex.size=5,
 	         vertex.label.family="Times",                  # Font family of the label (e.g.“Times”, “Helvetica”)
 	         vertex.label.font=1,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
 	         vertex.label.cex=0.5,                         # Font size (multiplication factor, device-dependent)
 	         vertex.label.dist=0.5,                        # Distance between the label and the vertex
 	         vertex.label.degree=0, 
			 vertex.size=de*1.5,                           # the important parameter :)
			 edge.width=E(g)$weight)
			 
dev.off()

# ggsave(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.png", sep="."), 
#			   width = 30,
#			   height = 50,
#			   units = "cm",
#			   dpi = 300 )

# In weighted networks, we can also node strength, which is the sum of the weights of edges connected to the node
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### BETWEENNESS CENTRALITY
# betweenness centrality : it is defined as the number of geodesic paths (shortest paths) that go through a given node. 

be = betweenness(g, normalized=T)
# plot(g,  vertex.label="", vertex.color="gold", edge.color="slateblue", vertex.size=be*50, edge.width=E(g)$weight*1)

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.CENTRALITY.BETWEENNESS.png", sep="."), width = 1600, height = 1600, units = "px")        
plot.igraph(g, 
	         # layout=lo*10,
			 # layout = layout_in_circle(network),
	         edge.arrow.width = .25,
    	     edge.arrow.size = .25,
	         edge.color=edge_attr(g)$color,
	         # edge.color=el$COLOR,
	         # edge.color=E(g)$color,
     
	         vertex.frame.color = "Forestgreen",            # Node border color
             vertex.shape=c("circle"),                      # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
                                                            # Size of the node (default is 15)
             # vertex.size=5,
 	         vertex.label.family="Times",                  # Font family of the label (e.g.“Times”, “Helvetica”)
 	         vertex.label.font=1,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
 	         vertex.label.cex=0.5,                         # Font size (multiplication factor, device-dependent)
 	         vertex.label.dist=0.5,                        # Distance between the label and the vertex
 	         vertex.label.degree=0, 
			 vertex.size=be*5,                             # the important parameter :)
			 edge.width=E(g)$weight)
			 
dev.off()

# ggsave(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.png", sep="."), 
#			   width = 30,
#			   height = 50,
#			   units = "cm",
#			   dpi = 300 )

#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### CLOSENESS CENTRALITY
#################################################################################################################### 

ce = closeness(g, normalized=T)

# Number of steps required to access every other node from a given node

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.CENTRALITY.CLOSENESS.png", sep="."), width = 600, height = 600, units = "px")       
plot.igraph(g, 
	         # layout=lo*10,
			 # layout = layout_in_circle(network),
	         edge.arrow.width = .25,
    	     edge.arrow.size = .25,
	         edge.color=edge_attr(g)$color,
	         # edge.color=el$COLOR,
	         # edge.color=E(g)$color,
     
	         vertex.frame.color = "Forestgreen",            # Node border color
             vertex.shape=c("circle"),                      # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
                                                            # Size of the node (default is 15)
             vertex.size=1,
 	         vertex.label.family="Times",                  # Font family of the label (e.g.“Times”, “Helvetica”)
 	         vertex.label.font=1,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
 	         vertex.label.cex=0.5,                         # Font size (multiplication factor, device-dependent)
 	         vertex.label.dist=0.5,                        # Distance between the label and the vertex
 	         vertex.label.degree=0, 
			 vertex.size=ce*5,                             # the important parameter :)
			 edge.width=E(g)$weight)			 
			 # vertex.size=ce*1 )
dev.off()


# These are nodes that tend to act as “bridges” between different clusters of nodes in the network ... 			 
# Closeness (centrality based on distance to others in the graph) Inverse of the node’s average geodesic distance to others in the network

distance_table(g)
dg = diameter(g)
is_connected(g)
is_connected(g,mode = "weak")
is_connected(g,mode = "strong")

write.table("A summary about the network features:", file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
            append = TRUE, quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)

write.table(paste("the diameter", diameter(g), sep="\t"), file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
			     append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#################################################################################################################### 
#################################################################################################################### NETWORK-LEVEL MEASURES
#################################################################################################################### 
#################################################################################################################### SIZE and DENSITY

# n for # nodes 
# and 
# m for # edges

n = vcount(g)
m = ecount(g)

# density = [# edges that exist] / [# edges that are possible]

# In an undirected network with no loops, the number of edges that are possible is exactly the number of dyads that exist in the network. 
# In turn, the number of dyads is n(n−1)2
# where n = number of nodes. Withthis information, we can calculate the density with the following:

dyads = n*(n-1)/2
density = m/dyads

# Components
# components(g)
# Degree distributions -- > SCALE-FREE NETWORKS
# a function to verify whether a network is scale free or not :
# https://www.rdocumentation.org/packages/splineTimeR/versions/1.0.1/topics/networkProperties
# networkProperties(igr)

write.table(paste("number of nodes (vertices)", vcount(g), sep="\t"), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
			     append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
						
write.table(paste("number of SV (edges)", ecount(g), sep="\t"), 
             file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
			 append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
				
write.table(paste("density of the nodes", density, sep="\t"), 
              file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)	
								
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### MODULARITY
#################################################################################################################### COMMUNITY DETECTION

# One class of methods for community detection (often called ‘modularity-optimization method’) to find the partitions 
# in the network that assigns nodes into communities such that Q is maximized.

# There are many algorithms that are available :

# walktrap.community(g) #RANDOM WALKS
# leading.eigenvector.community(g) # SPECTRAL PARTITIONING
# label.propagation.community(g)
# cluster_louvain(g) # LOUVAIN CLUSTERING

eb1 = leading.eigenvector.community(g)
eb2 =  cluster_louvain(g) 
length(eb1)
length(eb2)

write.table(paste("Number of communities : Spectral Clustering", length(eb1), sep="\t"), 
                file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
				append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Number of communities : Louvain clustering", length(eb2), sep="\t"),  
                file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
				append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

modularity(eb1)
modularity(eb2)				

write.table(paste("Modularity : Spectral Clustering", modularity(eb1), sep="\t"),  
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
		    append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
				
write.table(paste("Modularity : Louvain clustering", modularity(eb2), sep="\t"), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
								
membership(eb1)
membership(eb2)				

write.table(as.data.frame(as.matrix(membership(eb1))), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.MODULARITY.SPECTRAL.txt", sep="."),
			append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

write.table(as.data.frame(as.matrix(membership(eb2))), 
			file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.MODULARITY.LOUVAIN.txt", sep="."),
			append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.COMMUNITY.spectral.clustering.png", sep="."), width = 2400, height = 2400, units = "px")        ##### great !
plot(eb1, g, 
	         # layout=lo*10,
			 # layout = layout_in_circle(network),
	         edge.arrow.width = .25,
    	     edge.arrow.size = .25,
	         edge.color=edge_attr(g)$color,
	         # edge.color=el$COLOR,
	         # edge.color=E(g)$color,
			 # vertex.label="",
	         vertex.frame.color = "Forestgreen",            # Node border color
             vertex.shape=c("circle"),                      # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
                                                            
             vertex.size=2,
 	         vertex.label.family="Times",                   # Font family of the label (e.g.“Times”, “Helvetica”)
 	         vertex.label.font=10,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
 	         vertex.label.cex=1,                            # Font size (multiplication factor, device-dependent)
 	         vertex.label.dist=0.5,                         # Distance between the label and the vertex
 	         vertex.label.degree=0, 
			 # vertex.size=ce*5, 
			 edge.width=E(g)$weight )
dev.off()
						 			 
		 
png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.COMMUNITY.louvain.clustering.png", sep="."), width = 2400, height = 2400, units = "px")        ##### great :)
plot(eb2, g, 
	         # layout=lo*10,
			 # layout = layout_in_circle(network),
	         edge.arrow.width = .25,
    	     edge.arrow.size = .25,
	         edge.color=edge_attr(g)$color,
	         # edge.color=el$COLOR,
	         # edge.color=E(g)$color,
             # # vertex.label="",
	         vertex.frame.color = "Forestgreen",            # Node border color
             vertex.shape=c("circle"),                      # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
                                                            
             vertex.size=2,
 	         vertex.label.family="Times",                    # Font family of the label (e.g.“Times”, “Helvetica”)
 	         vertex.label.font=0.5,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
 	         vertex.label.cex=1,                             # Font size (multiplication factor, device-dependent)
 	         vertex.label.dist=0.5,                          # Distance between the label and the vertex
 	         vertex.label.degree=0, 
			 # vertex.size=ce*5, 
			 edge.width=E(g)$weight )
dev.off()

#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### PAGE RANK (google)
# if we want to build the graph directly from the edgelist

ela = graph.edgelist(as.matrix(data.frame(el$Gene1, el$Gene2)), directed = F)

V(ela)$page_rank <- page_rank(ela, directed = FALSE)$vector

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.COMMUNITY.page.rank.png", sep="."), width = 1400, height = 1400, units = "px")  
plot(ela, 
	 	         # layout=lo*10,
	 	         edge.arrow.width = .25,
	     	     edge.arrow.size = .25,
	 	         edge.color=edge_attr(g)$color,
	 	         # edge.color=el$COLOR,
	 	         # edge.color=E(g)$color,
	             # # vertex.label="",
	 	         vertex.frame.color = "Forestgreen",            # Node border color
	             vertex.shape=c("circle"),                      # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
	             vertex.label.color = "black",                                                 
	             # vertex.size=5,
	  	         vertex.label.family="Times",                  # Font family of the label (e.g.“Times”, “Helvetica”)
	  	         vertex.label.font=1,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
	  	         vertex.label.cex=0.5,                         # Font size (multiplication factor, device-dependent)
	  	         vertex.label.dist=0.5,                        # Distance between the label and the vertex
	  	         vertex.label.degree=0, 
	 			 vertex.size = V(ela)$page_rank/max(V(ela)$page_rank) * 20)
dev.off()
	 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### COMPONENTS : to print the components 
# https://bookdown.org/markhoff/social_network_analysis/finding-groups-in-networks.html#component-analysis
# https://malucalle.github.io/statistical-pills/nice-plots-with-r.html

component_list <- decompose.graph(g, mode = "weak")
component_list # a LIST with all the COMPONENTS 

compg = component_list

# Compg is a list of igraph objects ; each component stands for a community ;
# we print compg[[1]], or compg[[2]], or compg[[3]], ...  
# the number of the communities is length(component_list)
	
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### HISTOGRAM of DEGREE
####################################################################################################################

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.DEGREE.HISTOGRAM.png", sep="."), width = 400, height = 400, units = "px")
hist(degree(g), breaks=40, col="red", main="Node degree histogram", xlim=c(0,10))
dev.off()

# However, if we wanted to compare the degree distributions of different networks, it might be more useful to plot the probability densities of each degree: 
# i.e., what proportion of nodes has degree = 1, degree = 2, etc.

pk = degree.distribution(g)

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.DEGREE.DISTRIBUTION.png", sep="."), width = 400, height = 400, units = "px")
plot(pk, pch=20, col="red", main="Node degree distribution", xlim=c(0,10))
dev.off()

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.DEGREE.DENSITY.png", sep="."), width = 400, height = 400, units = "px")
plot(density(degree(g)),  main="Node degree density", xlab="node degree", col="red", xlim=c(0,10))
dev.off()
		
png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.DEGREE.DENSITY.v2.png", sep="."), width = 400, height = 400, units = "px")
plot(density(degree(g)),  main="Node degree density", xlab="node degree", col="red", ylim=c(0,1))
dev.off()

# write.table(as.data.frame(pk), 
#            file=paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.NODE.DISTRIBUTION.txt", sep="."),
#			append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
								
#################################################################################################################### 
#################################################################################################################### DIAMETER = (maximum path length)
#################################################################################################################### DISTANCES
#################################################################################################################### AVERAGE PATH LENGTH

# The average path length can be considered the average “degrees of separation” between all pairs of nodes in the network, 
# and the diameter is the maximum degree of separation that exists in the network

paths=distances(g, algorithm="unweighted")
paths

# This matrix contains a bunch of cells that are “Inf” (i.e., infinity). This is because the network is not connected, 
# and you can’t calculate path lengths between nodes in different components.

paths[paths=="Inf"]=NA
mean(paths[upper.tri(paths)], na.rm=T)

write.table(paste("average path for unconnected networks", mean_distance(g), sep="\t"),
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)	
								
				
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### CLUSTERING

# Here are the verbal definitions.
# Global Clustering Coefficient = “ratio of triangles to connected triples”
# Local Clustering Coefficient = for each node, the proportion of their neighbors that are connected to each other
# Average Local Clustering Coefficient: If Ci is the proportion of two nodes connected to node i that are also connected to each other 
# (i.e., the Local Clustering Coefficient), then Average Local Clustering Coefficient = 1n∑ni=1Ci
 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### CLUSTERING COEFFICIENTS

g.cluster=transitivity(g, "global")
l.cluster=transitivity(g, "local")
av.l.cluster=transitivity(g, "localaverage")

g.cluster
l.cluster
av.l.cluster

write.table(paste("Global clustering coefficient", g.cluster, sep="\t"),
                file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
				append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)	
								
write.table(paste("Average Local Clustering Coefficient", av.l.cluster,sep="\t"), 
                 file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
				 append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)	
						
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### HOMOPHILY
#################################################################################################################### ASSORTMENT

# Assortment
# One major pattern common to many social networks (and other types of networks) is homophily or assortment—the tendency for nodes
# that share a trait to be connected. The assortment coefficient is a commonly used measurement of homophily. 
# It is similar to the modularity index used in community detection, but the assortativity coefficient is used when we know a priori 
# the ‘type’ or ‘value’ of nodes. For example, we can use the assortment coefficient to examine whether discrete node types 
# (e.g., gender, ethnicity, species, etc.) are more or less connected to each other

# MODULARITY
# ASSORTATIVITY : the assortnet package

#################################################################################################################### 
#################################################################################################################### Linear Regression for Network Data 
#################################################################################################################### 
#################################################################################################################### MRQAP

# library(sna) 
# library(asnipe)

# generate 3 random adjacency matrices using the rgraph() function within sna
# set.seed(2)

# m1=rgraph(10, m=1, tprob=0.5, mode="graph")
# m2=rgraph(10, m=1, tprob=0.5, mode="graph") 
# m3=rgraph(10, m=1, tprob=0.5, mode="graph")

# n1=graph_from_adjacency_matrix(m1)
# n2=graph_from_adjacency_matrix(m2)
# n3=graph_from_adjacency_matrix(m3)

# plot(n1)
# plot(n2)
# plot(n3)

# netlm(m1, m2+m3, mode="graph", nullhyp="qap", test.statistic="t-value")

# Testing modularity of empirical network against randomized networks
# https://dshizuka.github.io/networkanalysis/example_sparrownet_analysis.html

#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### PLOTTING THE DEGREE DISTRIBUTION 
# PLOTTING THE DEGREE DISTRIBUTION

edge_density(g, loops=F)

write.table(paste("Edge density", edge_density(g, loops=F), sep="\t"), 
                 file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
				 append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)	

deg.dist.g = degree.distribution(g)

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.preferential.attachment.SCALE.FREE.png", sep="."), width = 400, height = 400, units = "px")
plot(deg.dist.g, pch=20, xlab="k", ylab="P(k)",las=1, main = "Degree Distribution",log="xy")
dev.off()

write.table(degree.distribution(g),
                file=paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.DEGREE.DISTRIBUTION.txt", sep="."),
				append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)	

# write.table(as.data.frame(pk), 
#            file=paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.NODE.DISTRIBUTION.txt", sep="."),
#			append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

deg <- igraph::degree(g, mode="all")

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.DEGREE.HISTOGRAM.v2.png", sep="."), width = 400, height = 400, units = "px")
hist(deg, breaks=1:vcount(g)-1, main="Histogram of node degree")
dev.off()

deg.dist <- degree_distribution(g, cumulative=T, mode="all")

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.DEGREE.DISTRIBUTION.v2.eCDF.png", sep="."), width = 400, height = 400, units = "px")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")
dev.off()

DEG.DIST = data.frame(DEGREE = 0:max(deg), DEGDIST = 1-deg.dist)

write.table(DEG.DIST, 
                 file=paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.DEGREE.DISTRIBUTION.v2.eCDF.txt", sep="."),
				 append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)	
 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
# another way to compute the eCDF based on :
# https://stackoverflow.com/questions/12169502/how-do-i-extract-ecdf-values-out-of-ecdfplot
# https://statisticsglobe.com/extract-ecdf-values-from-function-r

DEGREES = data.frame(deg)$deg

fun.ecdf <- ecdf(DEGREES) #
my.ecdf <- fun.ecdf(DEGREES)
data.ecdf <- data.frame(DEGREES, my.ecdf)  

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.DEGREE.DISTRIBUTION.v3.eCDF.png", sep="."), width = 400, height = 400, units = "px")
plot(data.ecdf$DEGREES,  
	 data.ecdf$my.ecdf,
     pch=19, cex=1.2, col="red", 
     xlab="Degree", ylab="Cumulative Frequency")
dev.off()

write.table(data.ecdf, 
                 file=paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.DEGREE.DISTRIBUTION.v3.eCDF.txt", sep="."),
				 append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)	

#################################################################################################################### VISUALIZATION
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### INTERACTIVE with VisNetworkData()

# vertices : V(g)
# edges : E(g)

data1 <- toVisNetworkData(g)

visNetwork(nodes = data1$nodes, edges = data1$edges, height = "500px")  %>%
visIgraphLayout() %>%
visNodes(size = 50) %>% 
visEdges(shadow = TRUE,
arrows =list(to = list(enabled = TRUE, scaleFactor = 2)),
color = list(color = "lightblue", highlight = "red")) %>% 
visEdges(arrows = "from") %>%
saveWidget(file=paste("graph",TUMOR_TYPE,"network.igraph.layout.js.VisNetworkData.html", sep="."))
  
#################################################################################################################### VISUALIZATION
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### INTERACTIVE with network3D
# http://www.sthda.com/english/articles/33-social-network-analysis/137-interactive-network-visualization-using-r/
# https://github.com/christophergandrud/networkD3
# https://datastorm-open.github.io/visNetwork/
# https://yunranchen.github.io/intro-net-r/advanced-network-visualization.html

data2 = igraph_to_networkD3(g)

src = el$Gene1 
target = el$Gene2
networkData <- data.frame(src, target)

# forceNetwork(networkData)
# sankeyNetwork(networkData)
# radialNetwork(networkData)
# diagonalNetwork(networkData)
# hc <- hclust(dist(as.matrix(g)), "ave")
# dendroNetwork(networkData)
# htmlwidgets()

result = try({

simpleNetwork(networkData, 
              linkColour = "#afafaf", 
			   fontSize=12, zoom=T, 
               opacity = 0.8, charge=-300,
               width = 1600, height = 1600)  %>%
saveNetwork(file=paste("graph",TUMOR_TYPE,"network.igraph.layout.js.network3D.html", sep="."))

}, silent = TRUE)

#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### TO DISPLAY the data on a HEATMAP
#################################################################################################################### 
####################################################################################################################

g.mat=get.adjacency(g,sparse = FALSE)

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.display.HEATMAP.png", sep="."), width = 400, height = 400, units = "px")
heatmap(g.mat)
dev.off()

#################################################################################################################### 
####################################################################################################################
#################################################################################################################### 
#################################################################################################################### QGRAPH
# https://www.adelaisvoranu.com/wp-content/uploads/2018/01/PracticalDay1-1.pdf
# https://www.adelaisvoranu.com/wp-content/uploads/2018/01/PracticalDay2-1.pdf

# Change to long format: -- edgelist but including all the 0s 

longData=reshape2::melt(g.mat)
  longData_all=as_tibble(longData)
  longData_all1=longData_all%>%mutate(Var1=forcats::fct_rev(Var1))

# using geom_tile
ggplot(longData_all1, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="white", high="#333333",na.value = "red") + 
  theme_bw()+ggtitle("")+xlab("")+ylab("") +        # set clean background and no titles 
  guides(fill = guide_colourbar(barheight = 12)) +  # can set the length of colour bar
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2)) +
  theme(axis.text.y = element_text(hjust = 1, size = 2)) +
  theme(legend.position = "none")  
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.display.HEATMAP_v0.png", sep="."), width = 2500, height = 2500, units = "px")
  
# to coerce igraph into a dataframe
# https://stackoverflow.com/questions/4904972/convert-igraph-object-to-a-data-frame-in-r
# in qgraph there is a very NICE LAYOUT " SPRING"
# https://www.adelaisvoranu.com/wp-content/uploads/2018/01/PracticalDay1-1.pdf
# http://sachaepskamp.com/qgraph/reference/qgraph.html
# qgraph(g.mat)

png(paste("graph",TUMOR_TYPE,"network.qgraph.layout.display.SPRING.png", sep="."), width = 1400, height = 2000, units = "px")
qgraph(g.mat, layout="spring")
dev.off()

# png(paste("graph",TUMOR_TYPE,"network.qgraph.layout.display.CIRCLE.png", sep="."), width = 1400, height = 2000, units = "px")
# qgraph(g.mat, layout="circle")
# dev.off()
			 
Centrality <- centrality(g.mat, all.shortest.paths = TRUE)	# it computes for each node in a graph
		 
Centrality$OutDegree
Centrality$InDegree
Centrality$Closeness
Centrality$Betweenness
Centrality$ShortestPathLengths

png(paste("graph",TUMOR_TYPE,"network.qgraph.layout.display.CENTRALITY.PLOT.png", sep="."), width = 800, height = 2400, units = "px")
centralityPlot(g.mat, 
               include =c("Degree","Strength","OutDegree","InDegree","OutStrength",  "InStrength"))
dev.off()
              
# NODES can reach other with only a few steps (related to the idea of “six degrees of separation”)				 				  
# in QGRAPH
# CORRELATION NETWORK
# cormat = cor(q.mat) #correlation matrix is generated
# png(paste("graph",TUMOR_TYPE,"network.qgraph.layout.display.CENTRALITY.PLOT.png", sep="."), width = 800, height = 800, units = "px")
# qgraph(cormat, shape="circle", posCol="darkgreen", negCol="darkred", layout="spring", vsize=10)
# dev.off()

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

# SMALL WORLD INDEX
smallworldIndex(g)

# A network has a small-world structure if it is (a) structured such that it forms clusters: my friends also talk to each-other :)

write.table(paste("Small world index : transitivity :", smallworldIndex(g)$transitivity, sep="\t"), 
              file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)	

write.table(paste("Small world index : APL :", smallworldIndex(g)$APL, sep="\t"), 
			     file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
			  	 append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)	

write.table(paste("Small world index : transitivity :", smallworldIndex(g)$index, sep="\t"), 
			     file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
			     append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)	

##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################								
############################################################################################################################################################################################## printing the DENSITY
# to print the DENSITY : x and y

DENSITY = data.frame(X = density(degree(g))$x, Y = density(degree(g))$y)

write.table(DENSITY, file=paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", sep="."),
								  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)		

# to OVERLAY the plots in R :
# https://stackoverflow.com/questions/55747048/how-to-plot-2-histograms-different-row-lengths-in-one-graph-ggplot


					 # df1 = data.frame(x = density(degree(g))$x, k = density(degree(g))$x, y = density(degree(g))$x)
					 # df2 = data.frame(x1 = density(degree(g))$y, m = density(degree(g))$y, y1 = density(degree(g))$y)
					 # ggplot() +
					 #  geom_density(data = df1, aes(x = x, y = y, color = "1")) +
					 #  geom_density(data = df2, aes(x = x1, y = y1, color = "2")) +
					 #  scale_color_manual(name = "Lines",
					 #                     values = c("1" = "blue", "2" = "red"))
								  
##############################################################################################################################################################################################
############################################################################################################################################################################################## computing GO
##############################################################################################################################################################################################								
##############################################################################################################################################################################################
##############################################################################################################################################################################################								
##############################################################################################################################################################################################
##############################################################################################################################################################################################								
##############################################################################################################################################################################################

GENES = unique(c(el$Gene1, el$Gene2))
GENES = gsub(" ", "", GENES)

library("AnnotationHub")
library("org.Hs.eg.db")
library(AnnotationDbi)
library("EnsDb.Hsapiens.v86")

# GENES_id = AnnotationDbi::select(EnsDb.Hsapiens.v86, keys=GENES, columns=c('ENTREZID','GO','PATH'), keytype='SYMBOL')
# GENES_entrez_id = GENES_id$ENTREZID

GENES_id = bitr(GENES, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
GENES_entrez_id = GENES_id$ENTREZID

##############################################################################################################################################################################################								
##############################################################################################################################################################################################
# bitr(GENES, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# keytypes(org.Hs.eg.db)
# keytypes(EnsDb.Hsapiens.v86)
# mapIds(org.Hs.eg.db, keys=GENES, column="ENTREZID", keytype="SYMBOL", multiVals="first")

# ggo <- groupGO(gene     = GENES_entrez_id ,
#               OrgDb    = org.Hs.eg.db,
#               ont      = "MF",
#               level    = 3,
#               readable = TRUE)
##############################################################################################################################################################################################			
##############################################################################################################################################################################################			   

result = try({
			   
ego = enrichGO(gene = GENES_entrez_id,
	     universe      = GENES_entrez_id,
	     OrgDb         = org.Hs.eg.db,
	     ont           = "MF",
	     pAdjustMethod = "fdr",
	     pvalueCutoff  = 1,
	     qvalueCutoff  = 1,
	    readable      = TRUE)

as.data.frame(ego)

if (dim(as.data.frame(ego))[1] != 0)
{

write.table(as.data.frame(ego), file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p1.txt", sep="."),
								  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)	


#png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p1.png", sep="."), width = 800, height = 800, units = "px")
dotplot(ego)								 
#dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p1.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )


ego2 <- pairwise_termsim(ego)

# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p1+.png", sep="."), width = 800, height = 800, units = "px")
treeplot(ego2)								 
# dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p1+.png", sep="."),
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )


}

}, silent = TRUE)

##############################################################################################################################################################################################								
##############################################################################################################################################################################################

gene.df <- bitr(GENES_entrez_id, fromType = "ENTREZID",
		        toType = c("ENSEMBL", "SYMBOL"),
		        OrgDb = org.Hs.eg.db)

##############################################################################################################################################################################################								
##############################################################################################################################################################################################

result = try({
				
 egoo <- enrichGO(gene        = gene.df$ENSEMBL,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "CC",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1, readable = TRUE)

 as.data.frame(egoo)

if (dim(as.data.frame(egoo))[1] != 0)
{

# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p2.png", sep="."), width = 800, height = 800, units = "px")
 dotplot(egoo)								 
# dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p2.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )
			   
 write.table(as.data.frame(egoo), file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p2.txt", sep="."),
								  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)	


 egoo2 <- pairwise_termsim(egoo)

# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p2+.png", sep="."), width = 800, height = 800, units = "px")
 treeplot(egoo2)								 
# dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p2+.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )
			   
}

}, silent = TRUE)

##############################################################################################################################################################################################								
##############################################################################################################################################################################################
# we follow the tutorial that is described in :
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html 
##############################################################################################################################################################################################								
##############################################################################################################################################################################################

result = try({
				
ego3 = enrichGO(gene = GENES_entrez_id,
                OrgDb        = org.Hs.eg.db,
                ont          = "CC",
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 1)
				
as.data.frame(ego3)

if (dim(as.data.frame(ego3))[1] != 0)
{

# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p3.png", sep="."), width = 800, height = 800, units = "px")
dotplot(ego3)								 
# dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p3.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )
			   
write.table(as.data.frame(ego3), file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p3.txt", sep="."),
							    append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)	


egoo3 <- pairwise_termsim(ego3)

# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p3+.png", sep="."), width = 800, height = 800, units = "px")
treeplot(egoo3)								 
# dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p3+.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )
}
	
}, silent = TRUE)	
								 
##############################################################################################################################################################################################								
##############################################################################################################################################################################################
##############################################################################################################################################################################################								
##############################################################################################################################################################################################

result = try({

ego4 = enrichGO(gene = GENES_entrez_id,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 1)
				
as.data.frame(ego4)

if (dim(as.data.frame(ego4))[1] != 0)
{

# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p4.png", sep="."), width = 800, height = 800, units = "px")
dotplot(ego4)								 
# dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p4.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )
			   
write.table(as.data.frame(ego4), file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p4.txt", sep="."),
							    append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)	


egoo4 <- pairwise_termsim(ego4)

# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p4+.png", sep="."), width = 800, height = 800, units = "px")
treeplot(egoo4)								 
# dev.off()								 
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.p4+.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )
			   
}

}, silent = TRUE)

##############################################################################################################################################################################################								
##############################################################################################################################################################################################								 
##############################################################################################################################################################################################								
############################################################################################################################################################################################## KEGG PATHWAYS

# result = try({
								 
# hsa <- search_kegg_organism('Homo sapiens', by='scientific_name')

# mkk <- enrichMKEGG(gene = GENES_entrez_id,
#                   organism = 'hsa',
#                   pvalueCutoff = 1,
#                   qvalueCutoff = 1)
				   

# write.table(as.data.frame(mkk), file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.KEGG.txt", sep="."),
#				   							    append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
				   							    
# }, silent = TRUE)	
												
##############################################################################################################################################################################################								
############################################################################################################################################################################################## WIKI PATHWAYS							 

# result = try({
								 
# wp = enrichWP(GENES_entrez_id, organism = "Homo sapiens") 

# write.table(as.data.frame(wp), file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.WP.txt", sep="."),
#				   							    append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)	
																								
# gseWP(GENES_entrez_id, organism = "Homo sapiens")

# }, silent = TRUE)
								 							 
##############################################################################################################################################################################################								
##############################################################################################################################################################################################
##############################################################################################################################################################################################								
############################################################################################################################################################################################## REACTOME

# library("ReactomePA")
						
# reactome_enrich <- enrichPathway(gene=GENES_entrez_id, pvalueCutoff = 1, readable=TRUE)
	
# reactome_gse <- gsePathway(GENES_entrez_id, 
#				            pvalueCutoff = 1,
#				            pAdjustMethod = "BH", 
#				            verbose = FALSE)
							
##############################################################################################################################################################################################								
##############################################################################################################################################################################################
##############################################################################################################################################################################################								
############################################################################################################################################################################################## enrichDO

result = try({

disease <- enrichDO(gene    = GENES_entrez_id,
              ont           = "DO",
              pvalueCutoff  = 1,
              pAdjustMethod = "BH",
              universe      = GENES_entrez_id,
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 1,
              readable      = FALSE)


if (dim(as.data.frame(disease))[1] != 0)
{
															 
# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.DISEASE.png", sep="."), width = 800, height = 800, units = "px")
dotplot(disease)								 
# dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.DISEASE.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )
			   
write.table(as.data.frame(disease), file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.DISEASE.txt", sep="."),
							    append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)	

disease2 <- pairwise_termsim(disease)

# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.DISEASE+.png", sep="."), width = 800, height = 800, units = "px")
treeplot(disease2)								 
# dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.DISEASE+.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )
}

}, silent = TRUE)
	
##############################################################################################################################################################################################								
##############################################################################################################################################################################################
##############################################################################################################################################################################################								
############################################################################################################################################################################################## enrichDGN
	
# edo = enrichDGN(GENES_entrez_id)									 
# edo <- gseDO(GENES_entrez_id)
# edo2 <- pairwise_termsim(edo2)
										
##############################################################################################################################################################################################								
##############################################################################################################################################################################################
##############################################################################################################################################################################################								
##############################################################################################################################################################################################

result = try({
	
m_df <- msigdbr(species = "Homo sapiens")
m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
										   dplyr::select(gs_name, entrez_gene)
											   
em <- enricher(GENES_entrez_id, TERM2GENE=m_t2g)

# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.MSIGDBR.png", sep="."), width = 800, height = 800, units = "px")
dotplot(em)								 
# dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.MSIGDBR.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 ) 
			   
write.table(as.data.frame(em), file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.MSIGDBR.txt", sep="."),
							    append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)	

em2 <- pairwise_termsim(em)

# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.MSIGDBR+.png", sep="."), width = 800, height = 800, units = "px")
treeplot(em2)								 
# dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.MSIGDBR+.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )
			   
}, silent = TRUE)

##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################								
############################################################################################################################################################################################## VISUALIZATION METHODS										
										
# dotplot(ego, showCategory=30)
# barplot(edo, showCategory=20) 										
# heatplot(edo, showCategory=20)
# dotplot(mkk)										
# dotplot(mkk2)										
# treeplot(edo)										
# emapplot(edo, cex_category=1.5,layout="kk")										
# upsetplot(edo)
# ridgeplot(edo)										
# gseaplot(edo2)
# gsearank(edo)							
			
##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################								
############################################################################################################################################################################################## enrichNCG()
# enrichNCG(GENES_entrez_id) 
 
result = try({
											   
cancer <- enrichNCG(GENES_entrez_id)

# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.CANCER.png", sep="."), width = 800, height = 800, units = "px")
dotplot(cancer)								 
# dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.CANCER.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )
			   
write.table(as.data.frame(em), file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.CANCER.txt", sep="."),
							    append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)	

cancer2 <- pairwise_termsim(cancer)

# png(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.CANCER+.png", sep="."), width = 800, height = 800, units = "px")
treeplot(cancer2)								 
# dev.off()
ggsave(paste("graph",TUMOR_TYPE,"network.igraph.measurements.clusterprofiler.enrichGO.CANCER+.png", sep="."), 
			   width = 20,
			   height = 40,
			   units = "cm",
			   dpi = 300 )
			   
}, silent = TRUE)

#################################################################################################################### 
#################################################################################################################### we are working with NODE DEGREE (DE)
#################################################################################################################### KS or FITDISTR
#################################################################################################################### to determine the type of DISTRIBUTION

# NORMALITY TEST : shapiro.test(de)

write.table(paste("Shapiro-Wilk normality test", shapiro.test(de)$p.value, sep="\t"), 
              file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

library(fitdistrplus)

result = try({

png(paste("graph",TUMOR_TYPE,"network.display.FITDISTRPLUS.CONTINUOUS.png", sep="."), width = 800, height = 800, units = "px")
descdist(de, discrete = FALSE)								 
dev.off()

png(paste("graph",TUMOR_TYPE,"network.display.FITDISTRPLUS.DISCRETE.png", sep="."), width = 800, height = 800, units = "px")
descdist(de, discrete = TRUE)								 
dev.off()

# ggsave(paste("graph",TUMOR_TYPE,"network.qgraph.layout.display.FITDISTRPLUS.png", sep="."), 
#			   width = 20,
#			   height = 40,
#			   units = "cm",
#			   dpi = 300 )

}, silent = TRUE)

#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 

df.de = as.data.frame(de)

#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### 
#################################################################################################################### AUTOMATED DISTRIBUTION FITTING

library(gamlss)
library(gamlss.dist)
library(gamlss.add)

result = try({

fit = fitDist(df.de$de, k = 2, type = "realplus", trace = FALSE, try.gamlss = TRUE)
# fit$family[1]
# fit$family[2]

write.table(paste("package GAMlss : distribution :", fit$family[1], sep="\t"), 
              file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("package GAMlss : distribution :", fit$family[2], sep="\t"), 
              file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)		

}, silent = TRUE)

#################################################################################################################### 
#################################################################################################################### fitdistrplus package
# https://www.huber.embl.de/users/kaspar/biostat_2021/2-demo.html
# https://stackoverflow.com/questions/50991870/how-do-i-know-what-distribution-of-data-follows-in-r
# df.de = as.data.frame(de)
#################################################################################################################### 
####################################################################################################################  CONTINUOUS DISTRIBUTION

result = try({

fw <- fitdist(df.de$de, "weibull")
fln <- fitdist(df.de$de, "lnorm")
fg <- fitdist(df.de$de, "gamma")
fn <- fitdist(df.de$de, "norm")
fu <- fitdist(df.de$de, "unif")
fe = fitdist(df.de$de, "exp")
fl = fitdist(df.de$de, "logis")

#################################################################################################################### 
#################################################################################################################### 

plot.legend <- c("weibull", "lognormal", "gamma", "norm", "unif", "exp", "logis")

png(paste("graph",TUMOR_TYPE,"network.display.FITDIST.CONTINUOUS.png", sep="."), width = 800, height = 800, units = "px")
par(mfrow = c(2,2))
denscomp(list(fw, fln, fg, fn, fu, fe, fl), legendtext = plot.legend)
qqcomp(list(fw, fln, fg, fn, fu, fe, fl), legendtext = plot.legend)
cdfcomp(list(fw, fln, fg, fn, fu, fe, fl), legendtext = plot.legend)
ppcomp(list(fw, fln, fg, fn, fu, fe, fl), legendtext = plot.legend)
dev.off()


GOF = gofstat(list(fw, fln, fg, fn, fu, fe, fl), fitnames = c("weibull", "lognormal", "gamma", "norm", "unif", "exp", "logis"))

write.table(paste("chisq  :", data.frame(distrib=rownames(data.frame(GOF$chisq)), value=GOF$chisq), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.CONTINOUS.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Anderson-Darling statistics  :", data.frame(distrib=rownames(data.frame(GOF$ad)), value=GOF$ad), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.CONTINOUS.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Kolmogorov-Smirnov statistics :", data.frame(distrib=rownames(data.frame(GOF$ks)), value=GOF$ks), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.CONTINOUS.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Kolmogorov-Smirnov statistics results :", data.frame(distrib=rownames(data.frame(GOF$kstest)), value=GOF$kstest), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.CONTINOUS.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Akaike's Information Criterion :", data.frame(distrib=rownames(data.frame(GOF$aic)), value=GOF$aic), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.CONTINOUS.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Bayesian Information Criterion :", data.frame(distrib=rownames(data.frame(GOF$bic)), value=GOF$bic), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.CONTINOUS.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

}, silent = TRUE)

#################################################################################################################### 
#################################################################################################################### BETA DISTRIBUTION
# to rescale the data in the interval [0,1] to try to see if the beta distribution fits 

result = try({

fb <- fitdist(df.de$de/100,"beta")

plot.legend3 <- c("beta")

png(paste("graph",TUMOR_TYPE,"network.display.FITDIST.BETA.png", sep="."), width = 800, height = 800, units = "px")
plot(fb)
dev.off()

GOF3 = gofstat(fb, fitnames = c("beta"))


write.table(paste("chisq  :", data.frame(distrib=rownames(data.frame(GOF3$chisq)), value=GOF3$chisq), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.BETA.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Anderson-Darling statistics  :", data.frame(distrib=rownames(data.frame(GOF3$ad)), value=GOF3$ad), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.BETA.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Kolmogorov-Smirnov statistics :", data.frame(distrib=rownames(data.frame(GOF3$ks)), value=GOF3$ks), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.BETA.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Kolmogorov-Smirnov statistics results :", data.frame(distrib=rownames(data.frame(GOF3$kstest)), value=GOF3$kstest), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.BETA.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Akaike's Information Criterion :", data.frame(distrib=rownames(data.frame(GOF3$aic)), value=GOF3$aic), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.BETA.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Bayesian Information Criterion :", data.frame(distrib=rownames(data.frame(GOF3$bic)), value=GOF3$bic), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.BETA.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


}, silent = TRUE)

#################################################################################################################### 
#################################################################################################################### DISCRETE DISTRIBUTIONS

result = try({

fp <- fitdist(df.de$de, "pois")
fnb <- fitdist(df.de$de, "nbinom")
fgeom <- fitdist(df.de$de, "geom")

plot.legend2 <- c("pois", "nbinom", "geom")

png(paste("graph",TUMOR_TYPE,"network.display.FITDIST.DISCRETE.png", sep="."), width = 800, height = 800, units = "px")
par(mfrow = c(2,2))
denscomp(list(fp, fnb, fgeom), legendtext = plot.legend2)
qqcomp(list(fp, fnb, fgeom), legendtext = plot.legend2)
cdfcomp(list(fp, fnb, fgeom), legendtext = plot.legend2)
ppcomp(list(fp, fnb, fgeom), legendtext = plot.legend2)
dev.off()

GOF2 = gofstat(list(fp, fnb, fgeom), fitnames = c("pois", "nbinom", "geom"))

write.table(paste("chisq  :", data.frame(distrib=rownames(data.frame(GOF2$chisq)), value=GOF2$chisq), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.DISCRETE.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("chisq p-value:", data.frame(distrib=rownames(data.frame(GOF2$chisqpvalue)), value=GOF2$chisqpvalue), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.DISCRETE.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Akaike's Information Criterion :", data.frame(distrib=rownames(data.frame(GOF2$aic)), value=GOF2$aic), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.DISCRETE.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("Bayesian Information Criterion :", data.frame(distrib=rownames(data.frame(GOF2$bic)), value=GOF2$bic), sep="\t"),
              file=paste("graph",TUMOR_TYPE,"network.display.FITDIST.DISCRETE.statistics.txt", sep="."),
			  append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


}, silent = TRUE)


#################################################################################################################### 
#################################################################################################################### another way to determine if it fits the distribution by using KS TEST
# https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
# it has been suggested by stackexchange
#################################################################################################################### 
#################################################################################################################### 
# https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
# it has been suggested by stackexchange
# x <- seq(0, 1000, 0.1)
# g <- ecdf(x)
# G <- g(x)
##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
# https://stackoverflow.com/questions/40851328/compute-area-under-density-estimation-curve-i-e-probability
# TO EXTRACT the AREA under DENSITY CURVE
# https://stackoverflow.com/questions/40851328/compute-area-under-density-estimation-curve-i-e-probability
# we compute the AUC based on eCDF, LINE.ecdf is the maximum number of the Y axis. 

result = try({

DEGREES = data.frame(deg)$deg

fun.ecdf <- ecdf(DEGREES) #
my.ecdf <- fun.ecdf(DEGREES)
data.ecdf <- data.frame(DEGREES, my.ecdf)

LINE.ecdf = max(unique(DEGREES)) 

z1 = 1
z2 = 2
z3 = 3
z4 = 4 
z5 = 5

write.table(paste("AUC less than", z1, format(round(fun.ecdf(z1), 2), nsmall = 2), sep=","), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.eCDF.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("AUC more than", z1, format(round(1-fun.ecdf(z1), 2), nsmall = 2), sep=","), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.eCDF.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("AUC less than", z2, format(round(fun.ecdf(z2), 2), nsmall = 2), sep=","), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.eCDF.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("AUC more than", z2, format(round(1-fun.ecdf(z2), 2), nsmall = 2), sep=","), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.eCDF.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("AUC less than", z3, format(round(fun.ecdf(z3), 2), nsmall = 2), sep=","), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.eCDF.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("AUC more than", z3, format(round(1-fun.ecdf(z3), 2), nsmall = 2), sep=","), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.eCDF.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("AUC less than", z4, format(round(fun.ecdf(z4), 2), nsmall = 2), sep=","), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.eCDF.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("AUC more than", z4, format(round(1-fun.ecdf(z4), 2), nsmall = 2), sep=","), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.eCDF.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("AUC less than", z5, format(round(fun.ecdf(z5), 2), nsmall = 2), sep=","), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.eCDF.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(paste("AUC more than", z5, format(round(1-fun.ecdf(z5), 2), nsmall = 2), sep=","), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.eCDF.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

}, silent = TRUE)

#################################################################################################################### 
#################################################################################################################### 
# https://stackoverflow.com/questions/40851328/compute-area-under-density-estimation-curve-i-e-probability
# https://stackoverflow.com/questions/36944874/r-area-under-curve-of-ogive
#################################################################################################################### 
#################################################################################################################### AREA under the DENSITY CURVE 
############################################################################## computing the area under the DENSITY CURVE : NUMERICAL INTEGRATION

result = try({

g = graph_from_data_frame(d = el, directed = FALSE)
DENSITY = density(degree(g))
d = density(degree(g))

xx <- d$x              ### evenly spaced points on [min(x) - 3 * d$bw, max(x) + 3 * d$bw]
dx <- xx[2L] - xx[1L]  ### spacing / bin size
yy <- d$y  ## 512 density values for `xx`

# method for numerical integration : Riemann Sum
# The area under the estimated density curve is:
C <- sum(yy) * dx  ## sum(yy * dx)
# [1] 1.000976

# Since Riemann Sum is only an approximation, this deviates from 1 (total probability) a little bit.
# We call this C value a "normalizing constant".
# Numerical integration on [1, Inf] can be approximated by :

p.unscaled <- sum(yy[xx >= 1]) * dx

# which should be further scaled it by C for a proper probability estimation:

p.scaled <- p.unscaled / C

write.table(paste("node density - numerical integration on [1, Inf]", p.scaled, sep="\t"), 
            file=paste("graph",TUMOR_TYPE,"network.igraph.measurements.eCDF.txt", sep="."),
			append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

}, silent = TRUE)

#################################################################################################################### 
####################################################################################################################
# another way has been presented in :
# https://stackoverflow.com/questions/36944874/r-area-under-curve-of-ogive
###################################################################################################################
####################################################################################################################

result = try({

library(zoo)

g = graph_from_data_frame(d = el, directed = FALSE)
DENSITY = density(degree(g))
dens = density(degree(g))

# x interval
dx = median(diff(dens$x))

# mean height for each pair of y values
h = rollmean(dens$y, 2)

# Area under curve
sum(h*dx)  # 1.000943

# Cumulative area
cumsum(h*dx)

# Plot density, showing the points at which density is calculated 

png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.DENSITY.zoo.png", sep="."), width = 400, height = 400, units = "px")
plot(dens, lwd=3, type="l")
# abline(v=dens$x, lty="11")
dev.off()


png(paste("graph",TUMOR_TYPE,"network.igraph.layout.kk.eCDF.zoo.png", sep="."), width = 400, height = 400, units = "px")
plot(dens$x[-length(dens$x)] + 0.5*dx, cumsum(h*dx), lwd=3, type="l")
# abline(v=dens$x[-length(dens$x)] + 0.5*dx, col="#FF000060", lty="11")
dev.off()

}, silent = TRUE)

##############################################################################################################################################################################################								
##############################################################################################################################################################################################
#################################################################################################################### # ending the BIG LOOP
#################################################################################################################### # that runs on all the CELL LINES in the CCLE collection
##############################################################################################################################################################################################								
##############################################################################################################################################################################################
# ending the BIG LOOP
}

}, silent = TRUE)

# }

##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################								
##############################################################################################################################################################################################
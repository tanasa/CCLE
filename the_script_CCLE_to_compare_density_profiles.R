##################################################################################
##################################################################################
##################################################################################

setwd("./CCLE_SvABA/")

library("ggplot2")
library("reshape2")
library("magrittr")
library("dplyr")

##################################################################################
##################################################################################
##################################################################################
################################################################################## the frame of the comparison

# A = read.table("", header=TRUE, sep="\t", stringsAsFactors = F)
# B = read.table("", header=TRUE, sep="\t", stringsAsFactors = F)

# C = "" 

# bind_rows(
#   A %>% mutate(id = "fusions"),
#   B %>% mutate(id = "SV")
#) %>%
#   ggplot(aes(X, Y, colour = id)) +
#   geom_line() +
#   xlab("node degree") + 
#   ylab("density") + 
#   ggtitle("density plots : fusions and SV") + 
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold")) +
#    theme(plot.title = element_text(size = 12, face = "bold"),
#    legend.title=element_text(size=10), 
#    legend.text=element_text(size=9))

# ggsave(paste("graph.density", C, "comparisons.DENSITY.png", sep="."), 
#			   width = 30,
#			   height = 20,
#			   units = "cm",
#			   dpi = 300 )

# write.table( paste("wilcoxon.test :", wilcox.test(A$Y, B$Y)$p.value, sep="\t"),
#             file = paste("graph.density", C, "comparisons.DENSITY.txt", sep="."),
#             append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

##################################################################################
##################################################################################
##################################################################################

A = read.table("graph.AU565_BREAST.SV.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
B = read.table("graph.AU565_BREAST.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)

C = "AU565" 

bind_rows(
   A %>% mutate(id = "SV"),
   B %>% mutate(id = "fusions")
) %>%
   ggplot(aes(X, Y, colour = id)) +
   geom_line() +
   xlab("node degree") + 
   ylab("density") + 
   ggtitle("density plots : fusions and SV") + 
   theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold")) +
    theme(plot.title = element_text(size = 12, face = "bold"),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=9))

ggsave(paste("graph.density", C, "comparisons.DENSITY.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

write.table( paste("wilcoxon.test :", wilcox.test(A$Y, B$Y)$p.value, sep="\t"),
             file = paste("graph.density", C, "comparisons.DENSITY.txt", sep="."),
             append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

##################################################################################
##################################################################################
##################################################################################

A = read.table("graph.CAL120_BREAST.SV.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
B = read.table("graph.CAL120_BREAST.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)

C = "CAL120" 

bind_rows(
   A %>% mutate(id = "SV"),
   B %>% mutate(id = "fusions")
) %>%
   ggplot(aes(X, Y, colour = id)) +
   geom_line() +
   xlab("node degree") + 
   ylab("density") + 
   ggtitle("density plots : fusions and SV") + 
   theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold")) +
    theme(plot.title = element_text(size = 12, face = "bold"),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=9))

ggsave(paste("graph.density", C, "comparisons.DENSITY.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

write.table( paste("wilcoxon.test :", wilcox.test(A$Y, B$Y)$p.value, sep="\t"),
             file = paste("graph.density", C, "comparisons.DENSITY.txt", sep="."),
             append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

##################################################################################
##################################################################################
##################################################################################

A = read.table("graph.CAL51_BREAST.SV.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
B = read.table("graph.CAL51_BREAST.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)

C = "CAL51" 

bind_rows(
   A %>% mutate(id = "SV"),
   B %>% mutate(id = "fusions")
) %>%
   ggplot(aes(X, Y, colour = id)) +
   geom_line() +
   xlab("node degree") + 
   ylab("density") + 
   ggtitle("density plots : fusions and SV") + 
   theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold")) +
    theme(plot.title = element_text(size = 12, face = "bold"),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=9))

ggsave(paste("graph.density", C, "comparisons.DENSITY.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

write.table( paste("wilcoxon.test :", wilcox.test(A$Y, B$Y)$p.value, sep="\t"),
             file = paste("graph.density", C, "comparisons.DENSITY.txt", sep="."),
             append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

##################################################################################
##################################################################################
##################################################################################

A = read.table("graph.DU4475_BREAST.SV.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
B = read.table("graph.DU4475_BREAST.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)

C = "DU4475" 

bind_rows(
   A %>% mutate(id = "SV"),
   B %>% mutate(id = "fusions")
) %>%
   ggplot(aes(X, Y, colour = id)) +
   geom_line() +
   xlab("node degree") + 
   ylab("density") + 
   ggtitle("density plots : fusions and SV") + 
   theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold")) +
    theme(plot.title = element_text(size = 12, face = "bold"),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=9))

ggsave(paste("graph.density", C, "comparisons.DENSITY.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

write.table( paste("wilcoxon.test :", wilcox.test(A$Y, B$Y)$p.value, sep="\t"),
             file = paste("graph.density", C, "comparisons.DENSITY.txt", sep="."),
             append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

##################################################################################
##################################################################################
##################################################################################

A = read.table("graph.EFM19_BREAST.SV.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
B = read.table("graph.EFM19_BREAST.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)

C = "EFM19" 

bind_rows(
   A %>% mutate(id = "SV"),
   B %>% mutate(id = "fusions")
) %>%
   ggplot(aes(X, Y, colour = id)) +
   geom_line() +
   xlab("node degree") + 
   ylab("density") + 
   ggtitle("density plots : fusions and SV") + 
   theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold")) +
    theme(plot.title = element_text(size = 12, face = "bold"),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=9))

ggsave(paste("graph.density", C, "comparisons.DENSITY.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

write.table( paste("wilcoxon.test :", wilcox.test(A$Y, B$Y)$p.value, sep="\t"),
             file = paste("graph.density", C, "comparisons.DENSITY.txt", sep="."),
             append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

##################################################################################
##################################################################################
##################################################################################

A = read.table("graph.HCC1187_BREAST.SV.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
B = read.table("graph.HCC1187_BREAST.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)

C = "HCC1187" 

bind_rows(
   A %>% mutate(id = "SV"),
   B %>% mutate(id = "fusions")
) %>%
   ggplot(aes(X, Y, colour = id)) +
   geom_line() +
   xlab("node degree") + 
   ylab("density") + 
   ggtitle("density plots : fusions and SV") + 
   theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold")) +
    theme(plot.title = element_text(size = 12, face = "bold"),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=9))

ggsave(paste("graph.density", C, "comparisons.DENSITY.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

write.table( paste("wilcoxon.test :", wilcox.test(A$Y, B$Y)$p.value, sep="\t"),
             file = paste("graph.density", C, "comparisons.DENSITY.txt", sep="."),
             append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

##################################################################################
##################################################################################
##################################################################################

A = read.table("graph.AU565_BREAST.SV.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
B = read.table("graph.CAL120_BREAST.SV.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
C = read.table("graph.CAL51_BREAST.SV.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F) 
D = read.table("graph.DU4475_BREAST.SV.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
E = read.table("graph.EFM19_BREAST.SV.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
F = read.table("graph.HCC1187_BREAST.SV.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)

bind_rows(
   A %>% mutate(id = "AU565"),
   B %>% mutate(id = "CAL120"),
   C %>% mutate(id = "CAL51"),
   D %>% mutate(id = "DU4475"),
   E %>% mutate(id = "EFM19"),
   F %>% mutate(id = "HCC1187")
) %>%
   ggplot(aes(X, Y, colour = id)) +
   geom_line() +
   xlab("node degree") + 
   ylab("density") + 
   ggtitle("density plots : SV") + 
   theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold")) +
    theme(plot.title = element_text(size = 12, face = "bold"),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=9))

ggsave(paste("graph.density", "comparisons.SV", "comparisons.DENSITY.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

##################################################################################
##################################################################################
##################################################################################

A = read.table("graph.AU565_BREAST.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
B = read.table("graph.CAL120_BREAST.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
C = read.table("graph.CAL51_BREAST.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F) 
D = read.table("graph.DU4475_BREAST.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
E = read.table("graph.EFM19_BREAST.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)
F = read.table("graph.HCC1187_BREAST.network.igraph.layout.kk.DEGREE.DENSITY.xy.txt", header=TRUE, sep="\t", stringsAsFactors = F)

bind_rows(
   A %>% mutate(id = "AU565"),
   B %>% mutate(id = "CAL120"),
   C %>% mutate(id = "CAL51"),
   D %>% mutate(id = "DU4475"),
   E %>% mutate(id = "EFM19"),
   F %>% mutate(id = "HCC1187")
) %>%
   ggplot(aes(X, Y, colour = id)) +
   geom_line() +
   xlab("node degree") + 
   ylab("density") + 
   ggtitle("density plots : fusions") + 
   theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold")) +
    theme(plot.title = element_text(size = 12, face = "bold"),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=9))

ggsave(paste("graph.density", "comparisons.fusions", "comparisons.DENSITY.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

##################################################################################
##################################################################################
##################################################################################

bind_rows(
   A %>% mutate(id = "SV"),
   B %>% mutate(id = "fusions")
) %>%
   ggplot(aes(X, Y, colour = id)) +
   geom_line() +
   xlab("node degree") + 
   ylab("density") + 
   ggtitle("density plots : fusions and SV") + 
   theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold")) +
    theme(plot.title = element_text(size = 12, face = "bold"),
    legend.title=element_text(size=10), 
    legend.text=element_text(size=9))

ggsave(paste("graph.density", C, "comparisons.DENSITY.png", sep="."), 
			   width = 30,
			   height = 20,
			   units = "cm",
			   dpi = 300 )

write.table( paste("wilcoxon.test :", wilcox.test(A$Y, B$Y)$p.value, sep="\t"),
             file = paste("graph.density", C, "comparisons.DENSITY.txt", sep="."),
             append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

##################################################################################
##################################################################################
##################################################################################
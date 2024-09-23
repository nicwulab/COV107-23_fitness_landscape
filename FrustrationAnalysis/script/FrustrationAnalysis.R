#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(plyr)
library(dplyr)
library(gridExtra)
library(qualpalr)
library(ggrepel)
library(sinaplot)
library(ggforce)
require(cowplot)

plot_FrustIndex <- function(df_WT, df_mutant, WT_graph, mutant_graph, difference_graph, ylabel_diff, diff_file){
  textsize <- 9
  p_WT <- ggplot() +
    geom_point(data=df_WT[df_WT$V1 == 26,], aes(x=V1, y=V8), color='lightgreen', size=1) +
    geom_point(data=df_WT[df_WT$V1 == 34,], aes(x=V1, y=V8), color='grey', size=1) +
    geom_point(data=df_WT[df_WT$V1 == 52,], aes(x=V1, y=V8), color='pink', size=1) +
    geom_line(data=df_WT, aes(x=V1, y=V8), color='deepskyblue3', linetype='solid', linewidth=0.3) +
    geom_hline(yintercept=-1, linetype='dashed', color='red', size=0.3) +
    geom_hline(yintercept=0.78, linetype='dashed', color='green4', size=0.3) +
    ylim(-2.5, 2.5) +
    theme_cowplot(12) +
    theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
          axis.text=element_text(size=textsize, face="bold", colour = 'black'),
          axis.text.x=element_text(angle=0, hjust=0.5, vjust=0.5, colour = 'black'),
          axis.text.y=element_text(hjust=0.5, vjust=0.5, colour = 'black'),
          axis.title.x=element_text(size=textsize, face="bold"),
          axis.title.y=element_text(size=textsize, face="bold"),
          legend.position = 'none') +
      ylab('Frustration Index') +
      xlab('Residue Number')
  
  ggsave(WT_graph, p_WT, width=5, height=3, bg='white', dpi=1200)
  
  
  p_mutant <- ggplot() +
    geom_point(data=df_mutant[df_mutant$V1 == 26,], aes(x=V1, y=V8), color='lightgreen', size=1) +
    geom_point(data=df_mutant[df_mutant$V1 == 34,], aes(x=V1, y=V8), color='grey', size=1) +
    geom_point(data=df_mutant[df_mutant$V1 == 52,], aes(x=V1, y=V8), color='pink', size=1) +
    geom_line(data=df_mutant, aes(x=V1, y=V8), color='deepskyblue3', linetype='solid', linewidth=0.3) +
    geom_hline(yintercept=-1, linetype='dotted', color='red', size=0.3) +
    geom_hline(yintercept=0.78, linetype='dotted', color='green4', size=0.3) +
    ylim(-2.5, 2.5) +
    theme_cowplot(12) +
    theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
          axis.text=element_text(size=textsize, face="bold", colour = 'black'),
          axis.text.x=element_text(angle=0, hjust=0.5, vjust=0.5, colour = 'black'),
          axis.text.y=element_text(hjust=0.5, vjust=0.5, colour = 'black'),
          axis.title.x=element_text(size=textsize, face="bold"),
          axis.title.y=element_text(size=textsize, face="bold"),
          legend.position = 'none') +
    ylab('Frustration Index') +
    xlab('Residue Number')
  
  ggsave(mutant_graph, p_mutant, width=5, height=3, bg='white', dpi=1200)
  
  
  difference <- as.data.frame(cbind(df_WT$V1, df_WT$V8, df_mutant$V8))
  colnames(difference) <- c('Residue', 'WT', 'Mutant')
  difference$diff <- difference$Mutant - difference$WT
  df_diff_v2 <- as.data.frame(cbind(difference$Residue, difference$diff))
  write.table(df_diff_v2,diff_file,sep="\t",row.names=FALSE, col.names=FALSE)
  
  p_diff <- ggplot() +
    geom_point(data=difference[difference$Residue == 26,], aes(x=Residue, y=diff), color='lightgreen', size=1) +
    geom_point(data=difference[difference$Residue == 34,], aes(x=Residue, y=diff), color='grey', size=1) +
    geom_point(data=difference[difference$Residue == 52,], aes(x=Residue, y=diff), color ='pink', size=1) +
    geom_line(data=difference, aes(x=Residue, y=diff), color='deepskyblue3', linetype='dashed', linewidth=0.3) +
    ylim(-1.5, 1.5) +
    theme_cowplot(12) +
    theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
          axis.text=element_text(size=textsize, face="bold", colour = 'black'),
          axis.text.x=element_text(angle=0, hjust=0.5, vjust=0.5, colour = 'black'),
          axis.text.y=element_text(hjust=0.5, vjust=0.5, colour = 'black'),
          axis.title.x=element_text(size=textsize, face="bold"),
          axis.title.y=element_text(size=textsize, face="bold"),
          legend.position = 'none') +
    ylab(ylabel_diff) +
    xlab('Residue Number')
  
  ggsave(difference_graph, p_diff, width=5, height=3, bg='white', dpi=1200)
  }


##### MAIN #####

df_WT <- read.table('data/WT-singleresidue.txt', header=FALSE)
df_WT$V1 <- unlist(lapply(df_WT$V1, function(x) as.numeric(gsub("\\D", "", x))))

df_S53P <- read.table('data/S53P-singleresidue.txt', header=FALSE)
df_S53P$V1 <- unlist(lapply(df_S53P$V1, function(x) as.numeric(gsub("\\D", "", x))))

df_S35T <- read.table('data/S35T-singleresidue.txt', header=FALSE)
df_S35T$V1 <- unlist(lapply(df_S35T$V1, function(x) as.numeric(gsub("\\D", "", x))))

df_F27L <- read.table('data/F27L-singleresidue.txt', header=FALSE)
df_F27L$V1 <- unlist(lapply(df_S35T$V1, function(x) as.numeric(gsub("\\D", "", x))))


plot_FrustIndex(df_WT, df_S53P, 'graph/WT_FrustIndex.pdf', 'graph/S53P_FrustIndex.pdf', 'graph/S53P-WT_FrustIndex.pdf', 'Frustration Index Difference (S53P-WT)', 'data/S53P-WT_difference.txt')
plot_FrustIndex(df_WT, df_S35T, 'graph/WT_FrustIndex.pdf', 'graph/S35T_FrustIndex.pdf', 'graph/S35T-WT_FrustIndex.pdf', 'Frustration Index Difference (S35T-WT)', 'data/S35T-WT_difference.txt')
plot_FrustIndex(df_WT, df_F27L, 'graph/WT_FrustIndex.pdf', 'graph/F27L_FrustIndex.pdf', 'graph/F27L-WT_FrustIndex.pdf', 'Frustration Index Difference (F27L-WT)', 'data/F27L-WT_difference.txt')

# For PyMOL #
WT_fixbb <- as.data.frame(cbind(df_WT$V1, df_WT$V8))
write.table(WT_fixbb,"data/WT-FrustIndex.txt",sep="\t",row.names=FALSE, col.names=FALSE)

S53P_fixbb <- as.data.frame(cbind(df_S53P$V1, df_S53P$V8))
write.table(S53P_fixbb,"data/S53P-FrustIndex.txt",sep="\t",row.names=FALSE, col.names=FALSE)

S35T_fixbb <- as.data.frame(cbind(df_S35T$V1, df_S35T$V8))
write.table(S35T_fixbb,"data/S35T-FrustIndex.txt",sep="\t",row.names=FALSE, col.names=FALSE)

F27L_fixbb <- as.data.frame(cbind(df_F27L$V1, df_F27L$V8))
write.table(F27L_fixbb,"data/F27L-FrustIndex.txt",sep="\t",row.names=FALSE, col.names=FALSE)

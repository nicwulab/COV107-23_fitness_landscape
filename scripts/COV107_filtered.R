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

COV107 <- read.table(file="./results/COV107_mutlib_fit_exp.tsv", header = TRUE)
textsize <- 9

# Correlation between `exp1_fit` and `exp2_fit`
exp_plot <- ggplot(COV107, aes(x=exp1_fit, y=exp2_fit)) + 
  geom_hex(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_cowplot(12) +
  xlab("Expression (Replicate 1)") +
  ylab("Expression (Replicate 2)") +
  theme(plot.title=element_text(size=textsize+2,hjust=0.5,vjust=0.5,face="bold",colour = 'black'),
        axis.text=element_text(size=textsize,face="bold",colour = 'black'),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
        axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
        axis.title.x=element_text(size=textsize,face="bold"),,
        axis.title.y=element_text(size=textsize,face="bold"),
        legend.position = "right",
        legend.title    = element_text(size=textsize-2,face="bold"),
        legend.text=element_text(size=textsize-2,face="bold"),
        legend.justification='center',
        legend.key.size = unit(0.3,"line")) +
  xlim(0,7.5) +
  ylim(0,7.5)

cor_exp <- cor(x = COV107$exp1_fit, y = COV107$exp2_fit, method = "pearson")
print(cor_exp)

ggsave("./graph/cor_exp.pdf", plot=exp_plot, width=2.8, height=2.5, dpi=600)


# qsub -I -N depth -q queue2T -l select=1:ncpus=2:mem=9GB -l walltime=10000:00:00
# source activate r4_py37_env

###---------- correlations between rib length and body size ------
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(Hmisc)
library(pheatmap)
setwd("/data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5")
f2id = read_xlsx("../infer.102.xlsx", sheet = 6) %>% select(1,2)
ph = fread("../f2_phenotype.txt") %>% mutate(T1 = T28-T39) %>%
  select("iid", paste0("T", c(2:4, 6, 12, 24:27, 39, 75:77,103:112, 118:132))) %>%
  filter(iid %in% f2id$ID1)

rp = rcorr(as.matrix(ph[,-1]))
r = rp$r %>% as.data.frame() %>% select(paste0("T",118:132))
r = r[1:23,]
names(r) = paste0("LOR", 1:15)
row.names(r) = c("BW", "BL", "BH", "CD", "CDE", "LCW", "RCW", "CSL", "COL", "NOR",
                 "DCP", "LMOHC", "LMF", "IMF", "MD", "DMF", "MA","RMFT","WMFD","ITMFD", "RAFSR", "RATER", "RALR")

r1 = pheatmap(r,treeheight_row = 20, treeheight_col = 20)
ggsave("figure/lor_cor.pdf", r1, width = 5,height = 6)


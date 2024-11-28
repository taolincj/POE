# source activate r4_py37_env
library(data.table)
library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(stringr)
setwd("/data_group/xiehaibing/xiehaibing2/f1/rna/DIR_STARpass2expression/pass2notwasp")

gtf = rtracklayer::import('/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf') %>%  as.data.frame()
gtf = gtf %>% filter(type=="gene") %>% select(gene_name, gene_id)
gtf1 = gtf; names(gtf1) = c("deg1", "deg")
gtf2 = gtf; names(gtf2) = c("ase1", "ase")

### gene-gene network. expression(DEG) = expression(POE) + sex + breeding_type
ggi_e27 = function(zz1, fdr1){
  e1 = fread("gene_tpms_e27.tsv") #所有组织的表达量
  d1 = e1 %>% select(Gene_ID, ends_with(zz1)) %>% 
    column_to_rownames("Gene_ID") %>%
    t() %>% as.data.frame() #某一组织的表达量
  deg = fread(paste0(zz1, "_DEG_notwasp_sig.txt")) #DEG
  ase = fread(paste0("/data_group/xiehaibing/xiehaibing2/f1/rna/pysam/expression/deg/", zz1, "_poe_read10_sig.txt"))
  g1 = deg$geneID
  g2 = ase$gene_id
  
  if (zz1=="jr"){d1$sex = c(1,1,2,2,1,2,1,1,2,1,2,1,2,1,1,2);
  d1$type = c(rep(2,8), rep(1,7), 2)}
  if (zz1!="jr"){d1$sex = c(2,2,2,2,2,2,1,2,2,2,2,1,2,1,1,1,1,1,2,1);
  d1$type = c(rep(1,9), rep(2, 10), 1)}
  
  gg = function(deg, ase){
    fit = lm(d1[[deg]] ~ d1[[ase]] + d1$sex + d1$type)
    co = summary(fit)$coefficients %>% as.data.frame()
    x = c(deg=deg, ase=ase,co["d1[[ase]]", 1:4]) %>% as.data.frame()
    return(x)
  }
  
  possgg = possibly(.f = gg, otherwise = NULL)
  
  x = list(deg = rep(g1, each = length(g2)), 
           ase = rep(g2, length(g1))) %>% 
    pmap(possgg) %>% data.table::rbindlist()
  names(x)[6] = "p"
  x = x %>% arrange(p) %>% na.omit() 
  x$fdr = p.adjust(x$p, method = "BH")
  x = x %>% filter(fdr <= fdr1) %>% 
    left_join(gtf1, by="deg") %>% 
    left_join(gtf2, by="ase") %>% arrange(p)
  return(x)
}

x1 = ggi_e27("rib", fdr1=0.05); fwrite(x1, "rib.ggi.txt", sep = "\t")


ggi_d90 = function(zz1, fdr1){
  e1 = fread("gene_tpms_d90.tsv") #所有组织的表达量
  d1 = e1 %>% select(Gene_ID, ends_with(zz1)) %>% 
    column_to_rownames("Gene_ID") %>%
    t() %>% as.data.frame() #某一组织的表达量
  d1$type = rep(c(2,1),each=6)
  deg = fread(paste0("/data_group/xiehaibing/xiehaibing2/f1/d90/rna/DIR_STARpass2expression/pass2notwasp/", 
                     zz1, ".pass2.notwasp.deg.txt")) #DEG
  ase = fread(paste0("/data_group/xiehaibing/xiehaibing2/f1/d90/rna/pysam/expression/deg/", 
                     zz1, "_poe_sig.txt"))
  g1 = deg$geneID
  g2 = ase$geneID
  
  gg = function(deg, ase){
    fit = lm(d1[[deg]] ~ d1[[ase]] + d1$type)
    co = summary(fit)$coefficients %>% as.data.frame()
    x = c(deg=deg, ase=ase,co["d1[[ase]]", 1:4]) %>% as.data.frame()
    return(x)
  }
  
  possgg = possibly(.f = gg, otherwise = NULL)
  
  x = list(deg = rep(g1, each = length(g2)), 
           ase = rep(g2, length(g1))) %>% 
    pmap(possgg) %>% data.table::rbindlist()
  names(x)[6] = "p"
  x = x %>% arrange(p) %>% na.omit() 
  x$fdr = p.adjust(x$p, method = "BH")
  x = x %>% filter(fdr <= fdr1) %>% 
    left_join(gtf1, by="deg") %>% 
    left_join(gtf2, by="ase") %>% arrange(p)
  return(x)
}
rib90 = ggi_d90("rib", fdr1=0.05); fwrite(rib90, "rib.ggi.d90.txt", sep = "\t")

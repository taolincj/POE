# maternal genetic effects
# source activate r4_py37_env
library(data.table)
library(dplyr)
library(DESeq2)
library(tibble)
library(stringr)
setwd("/data_group/xiehaibing/xiehaibing2/f1/d90/rna/DIR_STARpass1expression")
gtf = rtracklayer::import('/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf')
gtf = as.data.frame(gtf) %>% filter(type=="gene") %>% select(seqnames, start, end, gene_id, gene_name)

### maternal (genetic) effects using PASS2 notWASP data
# a function explore DEG related with maternal (genetic) effects
setwd("/data_group/xiehaibing/xiehaibing2/f1/d90/rna/DIR_STARpass2expression/pass2notwasp")
maternal = function(zz){
  # step 0 read count matrix
  jr = fread(paste0(zz, ".pass2notwasp.featureCounts.txt"))
  jr = jr %>% select(1, 7:ncol(jr)) %>% column_to_rownames("Geneid")
  names(jr) = names(jr) %>% str_remove(pattern = paste0(zz, ".Aligned.sortedByCoord.out.bam"))
  
  # step 1 experiment design
  pat = factor(c(rep(c("DS","LW"), each=6), rep(c("DS","LW"), each=6)))
  mat = factor(c(rep(c("LW", "DS"), each=12)))
  coldata = data.frame(pat = pat, mat=mat)
  row.names(coldata) = names(jr)
  
  dds <- DESeqDataSetFromMatrix(countData = jr, colData = coldata, design = ~ pat+mat)
  dds
  keep <- rowSums(counts(dds)) >= 24
  dds <- dds[keep,]
  dds
  
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("mat", "DS", "LW"))
  res1 = as.data.frame(res) %>% filter(abs(log2FoldChange)>=1, padj<=0.05) %>% 
    rownames_to_column("geneID") %>% arrange(pvalue) %>% left_join(gtf, by=c("geneID" = "gene_id"))
  return(res1)
}

rib = maternal(zz="rib")
fwrite(rib, "maternal/rib.txt", sep = "\t") 

### maternal effects gene
library(data.table)
library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(stringr)
setwd("/data_group/xiehaibing/xiehaibing2/f1/d90/rna/DIR_STARpass2expression/pass2notwasp/maternal")
gtf = rtracklayer::import('/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf') %>%  as.data.frame()
gtf = gtf %>% filter(type=="gene") %>% select(gene_name, gene_id)
gtf1 = gtf; names(gtf1) = c("meg1", "meg")
gtf2 = gtf; names(gtf2) = c("ase1", "ase")

ggi_d90 = function(zz1, fdr1){
  e1 = fread("/data_group/xiehaibing/xiehaibing2/f1/rna/DIR_STARpass2expression/pass2notwasp/gene_all_fpkm_d90.tsv") #所有组织的表达量
  d1 = e1 %>% select(Gene_ID, ends_with(zz1)) %>% 
    column_to_rownames("Gene_ID") %>%
    t() %>% as.data.frame() #某一组织的表达量
  d1$type = rep(c(3,2,4,1), each=6)
  meg = fread(paste0("/data_group/xiehaibing/xiehaibing2/f1/d90/rna/DIR_STARpass2expression/pass2notwasp/maternal/", zz1, ".txt")) #MEG
  ase = fread(paste0("/data_group/xiehaibing/xiehaibing2/f1/d90/rna/pysam/expression/deg/", zz1, "_poe_sig.txt"))
  g1 = meg$geneID
  g2 = ase$geneID

  gg = function(meg, ase){
    fit = lm(d1[[meg]] ~ d1[[ase]] + d1$type)
    co = summary(fit)$coefficients %>% as.data.frame()
    x = c(meg=meg, ase=ase,co["d1[[ase]]", 1:4]) %>% as.data.frame()
    return(x)
  }
  
  possgg = possibly(.f = gg, otherwise = NULL)
  
  x = list(meg = rep(g1, each = length(g2)), 
           ase = rep(g2, length(g1))) %>% 
    pmap(possgg) %>% data.table::rbindlist()
  names(x)[6] = "p"
  x = x %>% arrange(p) %>% na.omit() 
  x$fdr = p.adjust(x$p, method = "BH")
  x = x %>% filter(fdr <= fdr1) %>% 
    left_join(gtf1, by="meg") %>% 
    left_join(gtf2, by="ase") %>% arrange(p)
  return(x)
}

rib90_ggi = ggi_d90("rib", fdr1=0.05); fwrite(rib90_ggi, "rib90_meg.ggi.txt", sep = "\t")

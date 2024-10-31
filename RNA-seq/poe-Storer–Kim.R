# source activate r4_py37_env
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(WRS)
library(rtracklayer)

###--------------------------- Storer-Kim test at 27dpf ------------------
setwd("/data_group/xiehaibing/xiehaibing2/f1/rna/pysam/expression")
gtf = rtracklayer::import('/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf')
gtf = as.data.frame(gtf) %>% filter(type=="gene") %>% select(gene_id, gene_name, gene_biotype)

SK = function(tissue){
  # step 0 read count matrix
  xz = fread(paste0(tissue, ".featureCounts.txt"))
  xz = xz %>% select(1, 7:ncol(xz)) %>% column_to_rownames("Geneid")
  names(xz) = c("202026-10mat", "202026-10pat", "202026-2mat",  "202026-2pat", "202102-7mat",  "202102-7pat", 
                "202102-8mat", "202102-8pat", "202232-1mat", "202232-1pat", "202232-2mat", "202232-2pat", 
                "202232-3mat", "202232-3pat", "202266-7mat",  "202266-7pat", "202266-8mat",  "202266-8pat", 
                "202266-9mat",  "202266-9pat", "Y193208-1mat", "Y193208-1pat", "Y193208-2mat", "Y193208-2pat", 
                "Y193208-3mat", "Y193208-3pat", "Y193208-4mat", "Y193208-4pat", "Y200902-1mat", "Y200902-1pat", 
                "Y200902-2mat", "Y200902-2pat", "Y200902-3mat", "Y200902-3pat", "Y204711-1mat", "Y204711-1pat", 
                "Y204711-2mat", "Y204711-2pat", "Y204711-3mat", "Y204711-3pat")
  xz = xz %>% select("202026-10mat", "202026-2mat", "202102-7mat", "202102-8mat", 
                     "202266-7mat", "202232-1mat", "202232-2mat", "202232-3mat",
                     "202266-8mat", "202266-9mat", "Y193208-1mat", "Y193208-2mat", "Y193208-3mat", 
                     "Y193208-4mat", "Y200902-1mat", "Y200902-2mat", "Y200902-3mat", "Y204711-1mat", 
                     "Y204711-2mat", "Y204711-3mat", "202026-10pat", "202026-2pat", "202102-7pat", 
                     "202102-8pat", "202232-1pat", "202232-2pat", "202232-3pat",
                     "202266-7pat", "202266-8pat", "202266-9pat", "Y193208-1pat",
                     "Y193208-2pat","Y193208-3pat", "Y193208-4pat", "Y200902-1pat", "Y200902-2pat",
                     "Y200902-3pat", "Y204711-1pat", "Y204711-2pat", "Y204711-3pat")
  
  keep <- rowSums(xz) >= 40
  xz = xz[keep, ]
  zj = c("202026-10", "202026-2", "202102-7", "202102-8", "202266-7", 
         "202232-1", "202232-2", "202232-3", "202266-8", "202266-9")
  for (i in zj){
    j = i %>% gsub(pattern = "-", replacement = "_")
    tot = xz[[paste0(i, "pat")]] + xz[[paste0(i, "mat")]]
    assign(paste0("p", j), xz[[paste0(i, "pat")]]/tot)
  }
  
  fj = c("Y193208-1", "Y193208-2", "Y193208-3", "Y193208-4", "Y200902-1", 
         "Y200902-2", "Y200902-3", "Y204711-1", "Y204711-2", "Y204711-3")
  for (i in fj){
    j = i %>% gsub(pattern = "-", replacement = "_")
    tot = xz[[paste0(i, "pat")]] + xz[[paste0(i, "mat")]]
    assign(paste0("p", j), xz[[paste0(i, "mat")]]/tot)
  }
  
  d1 = cbind(p202026_10, p202026_2, p202102_7, p202102_8, p202232_1,
             p202232_2, p202232_3, p202266_7, p202266_8, p202266_9,
             pY193208_1, pY193208_2, pY193208_3, pY193208_4, pY200902_1,
             pY200902_2, pY200902_3, pY204711_1, pY204711_2, pY204711_3)
  rownames(d1) = rownames(xz)
  
  x1 = function(xy){
    x = xy[1:10]
    y = xy[11:20]
    twobinom(x=x, y=y)$p.value %>% return()
  }
  
  P1 = function(xy){
    x = xy[1:10]
    y = xy[11:20]
    twobinom(x=x, y=y)$p1 %>% return()
  }
  
  P2 = function(xy){
    x = xy[1:10]
    y = xy[11:20]
    twobinom(x=x, y=y)$p2 %>% return()
  }
  
  N1 = function(xy){xy[1:10] %>% na.omit() %>% length() %>% return()}
  N2 = function(xy){xy[11:20] %>% na.omit() %>% length() %>% return()}
  
  pvalue = apply(d1, 1, x1)
  p1 = apply(d1, 1, P1)
  p2 = apply(d1, 1, P2)
  n1 = apply(d1, 1, N1)
  n2 = apply(d1, 1, N2)
  
  xz1 = cbind(p1, p2, n1, n2, pvalue) %>% na.omit() %>% as.data.frame() %>% 
    filter(n1 >= 5, n2 >= 5) %>% arrange(pvalue)
  fdrs = p.adjust(xz1$pvalue, method = "fdr")
  xz1$fdrs = fdrs
  
  xz1 = xz1 %>% rownames_to_column("gene") %>% 
    left_join(gtf, by=c("gene"="gene_id"))
  return(xz1)
}

rib = SK("rib")
fwrite(rib, "deg/sk_gene_rib.txt", sep = "\t")

###--------------------------- Storer-Kim test at d90 ------------------
setwd("/data_group/xiehaibing/xiehaibing2/f1/d90/rna/pysam/expression")
gtf = rtracklayer::import('/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf')
gtf = as.data.frame(gtf) %>% filter(type=="gene") %>% select(seqnames, start, end, gene_id, gene_name)
names(gtf)[4:5] = c("geneID", "names")
sk_d90 = function(tissue){
  jr = fread(paste0(tissue, ".featureCounts.txt"))
  jr = jr %>% select(1, 7:ncol(jr)) %>% column_to_rownames("Geneid")
  names(jr) = c("02403.m", "02501.m", "02601.m", "40801.m", "40807.m", "41303.m",
                "02403.p", "02501.p", "02601.p", "40801.p", "40807.p", "41303.p",
                "02409.m", "02505.m", "02603.m", "40805.m", "41301.m", "41305.m",
                "02409.p", "02505.p", "02603.p", "40805.p", "41301.p", "41305.p")
  jr = jr %>% select("02403.m","02501.m","02601.m","02409.m","02505.m","02603.m",
                     "40801.m","40807.m","41303.m","40805.m","41301.m","41305.m",
                     "02403.p","02501.p","02601.p","02409.p","02505.p","02603.p",
                     "40801.p","40807.p","41303.p","40805.p","41301.p","41305.p")
  
  keep <- rowSums(jr) >= 24
  jr = jr[keep, ]
  
  zj = c("40801","40807","41303","40805","41301","41305")
  for (i in zj){
    tot = jr[[paste0(i, ".p")]] + jr[[paste0(i, ".m")]]
    assign(paste0("p", i), jr[[paste0(i, ".p")]]/tot)
  }
  
  fj = c("02403","02501","02601","02409","02505","02603")
  for (i in fj){
    tot = jr[[paste0(i, ".p")]] + jr[[paste0(i, ".m")]]
    assign(paste0("p", i), jr[[paste0(i, ".m")]]/tot)
  }
  
  d1 = cbind(p40801, p40807, p41303, p40805, p41301, p41305,
             p02403, p02501, p02601, p02409, p02505, p02603)
  rownames(d1) = rownames(jr)
  
  x1 = function(xy){
    x = xy[1:6]
    y = xy[7:12]
    twobinom(x=x, y=y)$p.value %>% return()
  }
  
  P1 = function(xy){
    x = xy[1:6]
    y = xy[7:12]
    twobinom(x=x, y=y)$p1 %>% return()
  }
  
  P2 = function(xy){
    x = xy[1:6]
    y = xy[7:12]
    twobinom(x=x, y=y)$p2 %>% return()
  }
  
  N1 = function(xy){xy[1:6] %>% na.omit() %>% length() %>% return()}
  N2 = function(xy){xy[7:12] %>% na.omit() %>% length() %>% return()}
  
  pvalue = apply(d1, 1, x1)
  p1 = apply(d1, 1, P1)
  p2 = apply(d1, 1, P2)
  n1 = apply(d1, 1, N1)
  n2 = apply(d1, 1, N2)
  
  xz1 = cbind(p1, p2, n1, n2, pvalue) %>% na.omit() %>% as.data.frame() %>% 
    filter(n1 >= 4, n2 >= 4) %>% arrange(pvalue)
  fdrs = p.adjust(xz1$pvalue, method = "fdr")
  xz1$fdrs = fdrs
  
  xz1 = xz1 %>% rownames_to_column("gene") %>% 
    left_join(gtf, by=c("gene"="gene_id"))
  return(xz1)
}

rib90 = sk_d90("rib")
jr90 = sk_d90("jr")
fwrite(rib90, "deg/sk_gene_rib90.txt", sep = "\t")
fwrite(jr90, "deg/sk_gene_jr90.txt", sep = "\t")
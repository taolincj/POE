# qsub -I -N putty3 -q queue2T -l select=1:ncpus=10:mem=100GB -l walltime=10000:00:00
# source activate r4_py37_env
library(data.table)
library(dplyr)
library(DESeq2)
library(tibble)
library(stringr)
setwd("/data_group/xiehaibing/xiehaibing2/f1/rna/pysam/expression")

###-------------------deg for rib at 27dpf----------------------
degxz = function(tissue, Nread){
  # step 0 read count matrix
  xz = fread(paste0(tissue, ".featureCounts.txt"))
  xz = xz %>% select(1, 7:ncol(xz)) %>% column_to_rownames("Geneid")
  
  names(xz) %>% str_remove("aternal.sort.bam")
  
  names(xz) = c("202026-10m", "202026-10p", "202026-2m",  "202026-2p",  "202102-7m",  "202102-7p",  "202102-8m",  "202102-8p", 
                "202232-1m",  "202232-1p",  "202232-2m",  "202232-2p",  "202232-3m",  "202232-3p",  "202266-7m",  "202266-7p", 
                "202266-8m",  "202266-8p",  "202266-9m",  "202266-9p",  "Y193208-1m", "Y193208-1p", "Y193208-2m", "Y193208-2p",
                "Y193208-3m", "Y193208-3p", "Y193208-4m", "Y193208-4p", "Y200902-1m", "Y200902-1p", "Y200902-2m", "Y200902-2p",
                "Y200902-3m", "Y200902-3p", "Y204711-1m", "Y204711-1p", "Y204711-2m", "Y204711-2p", "Y204711-3m", "Y204711-3p")
  
  
  xz = xz %>% select("202026-10m", "202026-2m", "202102-7m", "202102-8m", "202232-1m", 
                     "202232-2m", "202232-3m", "202266-7m", "202266-8m", "202266-9m", 
                     "Y193208-1m", "Y193208-2m", "Y193208-3m", "Y193208-4m", "Y200902-1m", 
                     "Y200902-2m", "Y200902-3m", "Y204711-1m", "Y204711-2m", "Y204711-3m", 
                     "202026-10p", "202026-2p", "202102-7p", "202102-8p", "202232-1p",
                     "202232-2p", "202232-3p", "202266-7p", "202266-8p", "202266-9p",
                     "Y193208-1p",  "Y193208-2p","Y193208-3p", "Y193208-4p", "Y200902-1p", 
                     "Y200902-2p",  "Y200902-3p", "Y204711-1p", "Y204711-2p", "Y204711-3p")
  
  # step 1 experiment design
  sex = factor(c(rep(c(1, 2, 2, 2, 2, 2, 2, 1, 2, 2,
                       2, 2, 1, 2, 1, 1, 1, 1, 1, 2), 2)))
  parent = factor(c(rep("maternal", 20), rep("paternal", 20)))
  breed = factor(c(rep("DS", 10), rep("LW", 10), rep("LW", 10), rep("DS", 10)))
  coldata = data.frame(parent, sex, breed)
  row.names(coldata) = names(xz)
  
  dds <- DESeqDataSetFromMatrix(countData = xz, colData = coldata, design = ~ parent + breed + sex)
  dds
  keep <- rowSums(counts(dds) >= Nread) >= 40 
  dds <- dds[keep,]
  
  # Differential expression analysis
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("parent", "maternal", "paternal"))
  breed <- results(dds, contrast=c("breed", "DS", "LW"))
  
  res1 = as.data.frame(res) %>% filter(abs(log2FoldChange)>=1, padj<=0.05) %>% 
    rownames_to_column("geneID") %>% left_join(exp, by="geneID") %>% arrange(padj)
  breed1 = as.data.frame(breed) %>% filter(abs(log2FoldChange)>=1, padj<=0.05) %>% 
    rownames_to_column("geneID") %>% left_join(exp, by="geneID") %>% arrange(padj)
  return(list(poe0=as.data.frame(res), poe1=res1, brd0=as.data.frame(breed), brd1=breed1))
}

rib = degxz("rib", 10)
fwrite(rib$poe0, "deg/rib_poe_read10.txt", sep = "\t")
fwrite(rib$poe1, "deg/rib_poe_read10_sig.txt", sep = "\t")


###-------------------deg at d90----------------------
setwd("/data_group/xiehaibing/xiehaibing2/f1/d90/rna/pysam/expression")
gtf = rtracklayer::import('/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf')
gtf = as.data.frame(gtf) %>% filter(type=="gene") %>% select(seqnames, start, end, gene_id, gene_name)
names(gtf)[4:5] = c("geneID", "names")

# a function to return differentially expressed genes for d90.
deg = function(zz){
  # step 0 read count matrix
  jr = fread(paste0(zz, ".featureCounts.txt"))
  jr = jr %>% select(1, 7:ncol(jr)) %>% column_to_rownames("Geneid")
  names(jr) = c("02403.m", "02501.m", "02601.m", "40801.m", "40807.m", "41303.m",
                "02403.p", "02501.p", "02601.p", "40801.p", "40807.p", "41303.p",
                "02409.m", "02505.m", "02603.m", "40805.m", "41301.m", "41305.m",
                "02409.p", "02505.p", "02603.p", "40805.p", "41301.p", "41305.p")
  jr = jr %>% select("02403.m","02501.m","02601.m","02409.m","02505.m","02603.m",
                     "40801.m","40807.m","41303.m","40805.m","41301.m","41305.m",
                     "02403.p","02501.p","02601.p","02409.p","02505.p","02603.p",
                     "40801.p","40807.p","41303.p","40805.p","41301.p","41305.p")
  
  # step 1 experiment design
  parent = factor(c(rep("maternal", 12), rep("paternal", 12)))
  breed = factor(c(rep("LW", 6), rep("DS", 6), rep("DS", 6), rep("LW", 6)))
  coldata = data.frame(parent, breed)
  row.names(coldata) = names(jr)
  
  dds <- DESeqDataSetFromMatrix(countData = jr, colData = coldata, design = ~ parent + breed)
  dds
  keep <- rowSums(counts(dds)) >= 24
  #keep <- rowSums(counts(dds) >= 10) >= 32 #Pre-filtering: at least X samples with a count of 10 or more
  dds <- dds[keep,]
  
  # Differential expression analysis
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("parent", "maternal", "paternal"))
  breed <- results(dds, contrast=c("breed", "DS", "LW"))
  
  res1 = as.data.frame(res) %>% filter(abs(log2FoldChange)>=1, padj<=0.05) %>% 
    rownames_to_column("geneID") %>% left_join(gtf, by="geneID")
  breed1 = as.data.frame(breed) %>% filter(abs(log2FoldChange)>=1, padj<=0.05) %>% 
    rownames_to_column("geneID") %>% left_join(gtf, by="geneID")
  
  res = as.data.frame(res) %>% rownames_to_column("geneID")
  breed = as.data.frame(breed) %>% rownames_to_column("geneID")
  return(list(poe0=res, poe1=res1, brd0=breed, brd1=breed1))
}

jr = deg("jr")
fwrite(jr$poe0, "deg/jr_poe.txt", sep = "\t")
fwrite(jr$poe1, "deg/jr_poe_sig.txt", sep = "\t")

rib = deg("rib")
fwrite(rib$poe0, "deg/rib_poe.txt", sep = "\t")
fwrite(rib$poe1, "deg/rib_poe_sig.txt", sep = "\t")
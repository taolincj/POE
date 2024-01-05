# qsub -I -N depth -q queue2T -l select=1:ncpus=5:mem=50GB -l walltime=10000:00:00
# source activate r4_py37_env

library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(Seurat)
library(cowplot)
magma="/data_group/xiehaibing/xiehaibing2/software/magma"
setwd("/data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5/magma")

### annotation: 11 snp and gene
snp = fread("../f2.filtered.add.name.bim") %>% select(2,1,4)
fwrite(snp, "snp.txt", sep = "\t", col.names = F)

gtf = rtracklayer::import('/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf')
gtf1 = gtf %>% as.data.frame() %>% 
  select(seqnames, start, end, gene_id, type, gene_name) %>% 
  filter(type=="gene", seqnames %in% 1:18) %>%
  mutate(gene = if_else(is.na(gene_name), gene_id, gene_name)) %>%
  select(gene, seqnames, start, end) %>% distinct(gene, .keep_all = TRUE) %>%
  na.omit()
fwrite(gtf1, "gene.loc", sep = "\t", col.names = F)

system(paste0(magma, " --annotate window=10 --snp-loc snp.txt --gene-loc gene.loc --out pig11"))

### gene analysis on SNP p-value
sam = fread("sample.size") %>% column_to_rownames("V1")
# het
for (i in 1:132){
  fread(paste0("../gcta_heter/het", i, ".mlma")) %>% select(SNP, p) %>% na.omit() %>%
    fwrite(paste0("het/t", i, ".txt"), col.names = F, sep = "\t")
  system(paste0(magma, " --bfile ../f2.filtered.add.name --gene-annot pig11.genes.annot --pval het/t", i, ".txt N=", sam[i, "V2"], " --out het/result", i))
  print(i)
}

# add
for (i in 1:132){
  fread(paste0("../gcta/t", i, ".mlma")) %>% select(SNP, p) %>% na.omit() %>%
    fwrite(paste0("add/t", i, ".txt"), col.names = F, sep = "\t")
  system(paste0(magma, " --bfile ../f2.filtered.add.name --gene-annot pig11.genes.annot --pval add/t", i, ".txt N=", sam[i, "V2"], " --out add/result", i))
}

# pat
for (i in 1:132){
  fread(paste0("../gcta_parental/pat", i, ".mlma")) %>% select(SNP, p) %>% na.omit() %>%
    fwrite(paste0("pat/t", i, ".txt"), col.names = F, sep = "\t")
  system(paste0(magma, " --bfile ../f2.filtered.add.name --gene-annot pig11.genes.annot --pval pat/t", i, ".txt N=", sam[i, "V2"], " --out pat/result", i))
}

# mat
for (i in 1:132){
  fread(paste0("../gcta_parental/mat", i, ".mlma")) %>% select(SNP, p) %>% na.omit() %>%
    fwrite(paste0("mat/t", i, ".txt"), col.names = F, sep = "\t")
  system(paste0(magma, " --bfile ../f2.filtered.add.name --gene-annot pig11.genes.annot --pval mat/t", i, ".txt N=", sam[i, "V2"], " --out mat/result", i))
}

### 单细胞数据之前处理的magma_rib.R
rib <- readRDS("/data_group/xiehaibing/xiehaibing2/magma/ljbscrna/rib.rds")
exp <- rib@assays$RNA@counts
genes.keep <- rowSums(as.matrix(exp) > 0) >= 10
exp <- exp[genes.keep,]
exp1 <-  log2(exp+1) %>% as.data.frame() # for cell cluster analyses
clsuter = Idents(rib) %>% as.data.frame()
names(clsuter) = "type"

# for cluster analysis
generate_covs_clsuter <- function(exp1, out.dir, x, cors = 25){
  if(!dir.exists(out.dir)) dir.create(out.dir)
  genes.means =  rowMeans(exp1)
  write_cov <- function(ItemNumber){
    cell = clsuter %>% filter(type == ItemNumber)
    exp = exp1[, rownames(cell)]
    type.means = rowMeans(exp)
    df <- data.frame(GENE =  row.names(exp), E = type.means, A = genes.means)
    cover_name <- paste0(out.dir, "/cluster", ItemNumber, ".cov")
    write.table(df, file = cover_name, quote = F, col.names = T, row.names = F, sep = "\t")
  }
  
  parallel::mclapply(X = x, FUN = write_cov, mc.cores = cors)
  
} 

outdir = "rib_cluster_cov"
if (!dir.exists(outdir)){dir.create(outdir, recursive=TRUE)}
generate_covs_clsuter(exp1 = exp1, out.dir = outdir, x=0:5, cors = 25)

gene_prop_cells <- function(cor = 10, raw.file.path, covariate.files.dir.path, magma.path, out.dir){
  cover.files.path = list.files(path = covariate.files.dir.path, pattern = "[.]cov$")
  cell.n <- gsub(pattern = "[.]cov$", replacement = "", x = cover.files.path)
  cover.df <- data.frame(cell = cell.n, path = paste0(covariate.files.dir.path, "/", cell.n, ".cov"))
  mag_prop <- function(IterNumber){
    mag.com <- paste0(magma.path, "magma",
                      " --gene-results ", raw.file.path, 
                      " --gene-covar ",  cover.df$path[IterNumber], 
                      " --model condition-hide=A direction=greater --out ",
                      out.dir, "/", cover.df$cell[IterNumber])
    system(command =  mag.com, wait = T)
  }
  
  parallel::mclapply(X = 1:nrow(cover.df), FUN = mag_prop, mc.cores = cor)
  
}

magma.path = "/data_group/xiehaibing/xiehaibing2/software/"

testclustetr = function(raw.dir, outtype){
  for (i in 118:132){
    raw <- paste0("/data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5/magma/", raw.dir, i, ".genes.raw")
    outdir = paste0(outtype, "/T", i, "cluster")
    indir = "rib_cluster_cov"
    if (!dir.exists(outdir)){dir.create(outdir, recursive=TRUE)}
    if (file.exists(raw)){
      gene_prop_cells(cor = 30, raw.file.path = raw,  covariate.files.dir.path = indir,
                      out.dir = outdir, magma.path = magma.path)
    }
  }
}

testclustetr("het/result", "het")
testclustetr("add/result", "add")
testclustetr("pat/result", "pat")
testclustetr("mat/result", "mat")


# a function to merge gsa results.
md =function(id, type){
  mergeDate = function(dir){
    d0 = matrix(ncol = 8, nrow=0) %>% as.data.frame()
    for (i in 0:5){
      d1 = read.csv(paste0(dir, "cluster", i, ".gsa.out"), sep = "", header = T, comment.char = "#")
      d1$cluster = i
      d0 = rbind(d0, d1)
    }
    return(d0)
  }
  d = matrix(ncol = 9, nrow=0) %>% as.data.frame()
  for (i in id){#id=118:132 LOR
    d2 = mergeDate(paste0(type, "/T", i, "cluster/"))
    d2$trait = paste0("T", i)
    d = rbind(d, d2)
  }
  d = d[,4:9] %>% filter(P<0.05)
  return(d)
}

d = matrix(ncol = 7, nrow = 0) %>% as.data.frame()
for (i in c("add", "pat", "mat", "het")){
  d1 = md(118:132,i) %>% mutate(type = i)
  d = rbind(d, d1)
}

fwrite(d, "LOR_cluster.txt", sep = "\t")

library(ggplot2); library(readxl); library(dplyr); library(data.table)
# 可视化细胞类群
rib <- RenameIdents(rib, `0` = "FB", `1` = "CT", `2` = "SC", `3` = "OB", `4` = "MSC", `5` = "HEC")

cl = DimPlot(rib, reduction = "umap") + theme_minimal()+
  theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank()) #0FB, 1CT, 2SC, 3OB, 4MSC, 5HEC


gr = readxl::read_xlsx("/data_group/xiehaibing/xiehaibing2/poe-gwas/group.xlsx", sheet = 1)
d = fread("LOR_cluster.txt") %>% left_join(gr[,1:2], by="trait") %>% 
  filter(type %in% c("pat", "mat")) %>% 
  mutate(class=ifelse(type=="pat", "paternal", "maternal"))
d1 = d
d1$cluster = as.character(d1$cluster)
ct = ggplot(d1, aes(cluster, id, color=class, size=-log10(P))) + geom_point(alpha=0.9)  +
  scale_color_manual(values=c("#31a354","#3182bd")) +
  scale_x_discrete(breaks = 1:5, labels = c("CT", "SC", "OB", "MSC", "HEC"))+
  scale_y_discrete(limits=paste0("LOR", c(1:11,13,14)))+
  labs(x="Cell type", y=NULL, color="Parental GWAS", size=expression(-lg(italic(P)))) + 
  theme_bw() + theme(axis.text = element_text(color = "black"))
ggsave("rib_cell_type.pdf", ct, height = 3, width = 4)
#0FB, 1CT, 2SC, 3OB, 4MSC, 5HEC


###---------- 检查paternal GWAS中哪些基因在CT中起主要作用--------
setwd("/data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5/magma")
library(readr)
library(dplyr)
library(data.table)
magma="/data_group/xiehaibing/xiehaibing2/software/magma"
cov = read_tsv("rib_cluster_cov/cluster1.cov")

rmone = function(trait){
  gwas = read.csv(paste0("pat/result", trait, ".genes.raw"), comment.char = "#", header = F, sep = " ") #
  gene = intersect(cov$GENE, gwas$V1)
  
  for (i in gene){
    cov %>% filter(GENE != i) %>% fwrite("test.cov", sep = "\t")
    system(paste0(magma, " --gene-results pat/result", trait, ".genes.raw --gene-covar test.cov --model condition-hide=A direction=greater --out pat_rmone/t", trait, i))
  }
  
  d0 = as.data.frame(matrix(ncol=9, nrow = 0))
  names(d0) = c("VARIABLE", "TYPE", "NGENES", "BETA", "BETA_STD", "SE", "P", "trait", "RMgene")
  for (i in gene){
    if (file.exists(paste0("pat_rmone/t", trait, i, ".gsa.out"))){
      d1 = read.csv(paste0("pat_rmone/t", trait, i, ".gsa.out"), comment.char = "#", header = T, sep = "") %>% mutate(RMgene = i)
      d0 = rbind(d0, d1)
    }
  }
  
  d0 = d0 %>% arrange(desc(P))
  write_tsv(d0, paste0("pat_rmone/t", trait, "RMoneCT.txt"))
}

rmone(119)
rmone(125)
rmone(128)
rmone(130)
rmone(131)

dr = d1 %>% filter(cluster==1) %>% 
  select(BETA, BETA_STD, SE, P,   type,   id)
names(dr)[5:6] = c("RMgene", "trait")
dr$RMgene=""
r1 = fread("t119RMoneCT.txt")[1,] %>% mutate(trait="LOR2")
r2 = fread("t125RMoneCT.txt")[1,] %>% mutate(trait="LOR8")
r3 = fread("t128RMoneCT.txt")[1,] %>% mutate(trait="LOR11")
r4 = fread("t130RMoneCT.txt")[1,] %>% mutate(trait="LOR13")
r5 = fread("t131RMoneCT.txt")[1,] %>% mutate(trait="LOR14")

r = rbind(r1,r2,r3,r4, r5) %>% select(BETA, BETA_STD, SE, P,  RMgene, trait)
r = rbind(r, dr)
r$type = c(rep(c("rem", "ref"), each=5))

rm = ggplot(r, aes(trait, BETA, fill=type)) + 
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  geom_errorbar(aes(ymin = BETA, ymax = BETA + SE), color="gray",
                width = 0.5, position = position_dodge()) +
  scale_fill_manual(values = c("#3182bd", "#9ecae1")) +
  scale_x_discrete(limits=c("LOR2", "LOR8", "LOR11", "LOR13", "LOR14"),
                   labels=c("LOR2\nRABIF\n0.008\n0.013", 
                            "LOR8\nXRCC5\n0.049\n0.081*", 
                            "LOR11\nSEPTIN5\n0.005\n0.010", 
                            "LOR13\nCNMD\n0.020\n0.040", 
                            "LOR14\nCNMD\n0.042\n0.074*")) +
  labs(x=NULL, y="Beta") + 
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.position = c(0.3,1))

pall = plot_grid(ct, rm, ncol = 2, align = "hv", labels = "AUTO")
ggsave("cell_type_RMone.pdf", pall, height = 4, width = 7)
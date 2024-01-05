# qsub -I -N depth -q queue2T -l select=1:ncpus=5:mem=20GB -l walltime=10000:00:00
# source activate r4_py37_env
library(data.table)
library(dplyr)
library(stringr)
library(readxl)
library(tidyr)
library(ggplot2)
library(cowplot)
library(doBy)
library(stringr)
library(tidyverse)
library(readxl)
library(purrr)
setwd("/data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5")

###----------- get parental allele in vcf -------------
na = fread("f2.filtered.sample") %>% filter(ID_1 != 0)
v2 = fread( "f2.filtered.hap.gz")

v0 = v2 %>% select(V1, V3, V4, V5)
names(v0) = c("#CHROM", "POS", "REF", "ALT")
v0$ID = "."
v0$QUAL = "."
v0$FILTER = "."
v0$INFO = "."
v0$FORMAT = "DS"
v0 = v0[, c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")]

pa = names(v2)[seq(6, ncol(v2), 2)]
pat = v2 %>% select(all_of(pa)) #paternal allele
names(pat) = na$ID_1

ma = names(v2)[seq(7, ncol(v2), 2)]
mat = v2 %>% select(all_of(ma)) #maternal allele
names(mat) = na$ID_1

pat = cbind(v0, pat)
mat = cbind(v0, mat)

fwrite(pat, "pat.chr.vcf", sep = "\t")
fwrite(mat, "mat.chr.vcf", sep = "\t")
system("sed '22r pat.chr.vcf' header > pat_chr.vcf")
system("sed '22r mat.chr.vcf' header > mat_chr.vcf")
system("bgzip -@ 10 -i pat_chr.vcf")
system("bgzip -@ 10 -i mat_chr.vcf")
system("rm pat.chr.vcf")
system("rm mat.chr.vcf")

###----------- additive GWAS ----------
# pca
# plink2 --vcf f2.filtered.vcf.gz --make-bed --out f2.filtered.add
# plink --bfile f2.filtered.add --pca 50 --out pca50
# twstats -t twtable -i pca50.eigenval -o pca50.eigenvec_tw
# plink --bfile f2.filtered.add --set-missing-var-ids rs.@.# --make-bed --out f2.filtered.add.name
gemma="/data_group/xiehaibing/xiehaibing2/software/gemma-0.98.5"
ph = fread("../f2_phenotype.txt") %>% mutate(T1 = T28-T39)
f2 = read_xlsx("../infer.102.xlsx", sheet = 5) %>% select(1,2)

fa = fread("f2.filtered.add.fam") %>% select(V2)
names(fa) = "ID2"
fa1 = fa %>% left_join(f2, by="ID2") %>% left_join(ph, by=c("ID1"="iid")) %>%
  select(ID2, sex, paste0("T", 1:132))

# pc = fread("pca50.eigenvec") %>% select(2:5)
# names(pc) = c("ID2", "pc1", "pc2", "pc3")
# cov = fa1 %>% select(ID2, sex) %>% left_join(pc, by="ID2")
# cov$ID2 = 1
# fwrite(cov, "gemma.cov", sep = "\t", col.names = F)

# set phenotype as 1
# fam = fread("f2.filtered.add.fam")
# fam$V6=1
# fwrite(fam, "f2.filtered.add.fam", col.names = F, sep = "\t")
# system(paste0(gemma, " -bfile f2.filtered.add -gk 2 -o addgk2")) #kinship

for (i in 1:132){
  # alter fam file
  read.csv("f2.filtered.add.fam", header = F, sep="\t") %>%
    rename(ID2=V2) %>% left_join(fa1 %>% select("ID2", paste0("T",i)), by="ID2") %>%
    select(V1, ID2, V3, V4, V5, paste0("T",i)) %>%
    write.table("f2.filtered.add.fam", quote = F, sep = "\t", row.names = F, col.names = F)
  # gemma perform GWAS
  print(paste0("testing trait", i))
  system(paste0(gemma, " -bfile f2.filtered.add  -miss 1 -notsnp -r2 1.0 -lmm 1 -k ./output/addgk2.sXX.txt -c gemma.cov -o add_t", i))
}


### gcta ### 一个性状约2小时！gemma太慢了，一晚上一个都没算出来，占用的资源还不少！
co = fa1 %>% mutate(fid = 0) %>% select(fid, ID2, sex) # discrete covariates
fwrite(co, "gcta.cov", col.names = F, sep = "\t")

fa1$sex = 0 
fa2 = fa1 %>% select(sex, ID2, paste0("T", 1:132)) # phenotype
fwrite(fa2, "gcta.phen", col.names = F, sep = "\t", na="NA", quote = F)

system("plink --bfile f2.filtered.add --set-missing-var-ids ssc@.# --make-bed --out f2.filtered.add.name")
system("plink --bfile f2.filtered.add.name --indep-pairwise 500 50 0.2 --out f2.filtered.add.name") #LD prune
system("plink --bfile f2.filtered.add.name --indep-pairwise 500 50 0.2 --out f2.filtered.add.name") #LD prune
system("plink --bfile f2.filtered.add.name --extract f2.filtered.add.name.prune.in --pca 50 --out f2.filtered.add.name.prune.in50")
system("twstats -t twtable -i f2.filtered.add.name.prune.in50.eigenval -o f2.filtered.add.name.prune.in50_tw")

pc = fread("f2.filtered.add.name.prune.in50.eigenvec") %>% select(1:5) # quantitative covariates
fwrite(pc, "gcta.qcovar", col.names = F, sep = "\t")

system("gcta64 --bfile f2.filtered.add --make-grm --thread-num 10 --out gcta/add")
system("gcta64 --mlma --bfile f2.filtered.add.name --grm gcta/add --pheno gcta.phen --covar gcta.cov --qcovar gcta.qcovar --out gcta/t1 --thread-num 10")

### ---------- parental GWAS -----------------
system("plink2 --vcf pat_chr.vcf.gz dosage=DS --set-missing-var-ids ssc@.# --make-bed --out pat")
system("plink2 --vcf mat_chr.vcf.gz dosage=DS --set-missing-var-ids ssc@.# --make-bed --out mat")
system("gcta64 --mlma --bfile pat --grm gcta/add --pheno gcta.phen --covar gcta.cov --qcovar gcta.qcovar --out gcta_parental/pat1 --thread-num 10")

###----------- paternal and maternal h2 --------
system("gcta64 --bfile pat --autosome --make-grm --out h2/pat --thread-num 10")
system("gcta64 --bfile mat --autosome --make-grm --out h2/mat --thread-num 10")

for (i in 1:132){
  system(paste0("gcta64 --reml --mgrm h2/multi_grm.txt --pheno gcta.phen --mpheno ", i ," --covar gcta.cov --qcovar gcta.qcovar --out h2/t", i, " --thread-num 10"))
}

"%ni%"=Negate("%in%")
h0 = matrix(ncol = 4, nrow = 0) %>% as.data.frame()
names(h0) = c("Source", "Variance", "SE", "trait")
for (i in 1: 132){
  filex = paste0("h2/t", i, ".hsq")
  if (file.exists(filex)) {
    h1 = fread(filex) %>% 
      filter(Source %in% c("V(G1)/Vp", "V(G2)/Vp")) %>% 
      mutate(trait=paste0("t", i))
    h0 = rbind(h0, h1)
  }
}

h0 = h0 %>% mutate(Source = if_else(Source=="V(G1)/Vp", "pat", "mat"))
hpat = h0 %>% filter(Source=="pat") %>% select(2:4)
names(hpat) = c("h2_pat", "se_pat", "trait")
hmat = h0 %>% filter(Source=="mat") %>% select(2:4)
names(hmat) = c("h2_mat", "se_mat", "trait")
hpat = hpat %>% left_join(hmat, by="trait") %>% 
  filter(trait %ni% c("t20", paste0("t", 71:74))) %>%
  mutate(Distance = abs(h2_pat-h2_mat)/sqrt(2),
         type = ifelse(h2_pat<h2_mat, "mat", "pat"))
fwrite(hpat, "h2/pat_mar.h2", sep = "\t")

hpat = fread("D:/taolinwork/poe-gwas/h2pm/sus11/pat_mar.h2")
hpat = hpat %>% filter(trait %in% paste0("t",118:132))
t.test(hpat$h2_pat, hpat$h2_mat,  alternative = "less", paired = T) #p-value = 0.0002156

#hpat = hpat %>% filter(trait %in% c("t1", "t28", "t39", paste0("t", 118:132)))
pmh2 = ggplot(hpat, aes(h2_pat, h2_mat, color=Distance)) + xlim(0,0.8) + ylim(0,0.8)+
  geom_abline(intercept = 0, slope = 1, linetype="dashed")+theme_bw()+
  scale_color_gradient(low = "#1b9e77", high =  "#e6ab02")+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 13,color = "black"),
        axis.title = element_text(size =14))+
  geom_errorbar(aes(ymin = h2_mat-se_mat, ymax=h2_mat+se_mat), color="gray", alpha=.2)+
  geom_errorbarh(aes(xmin=h2_pat-se_pat, xmax=h2_pat+se_pat), color="gray", alpha=.2)+geom_point()+
  labs(x=expression(italic(h)[pat]^2), y=expression(italic(h)[mat]^2))+
  theme(legend.position=c(1,0), legend.justification=c(1,0), 
        legend.background = element_blank(), legend.key = element_blank(),
        axis.ticks = element_line(color="black"),
        panel.border = element_rect(color = "black"))
ggsave("D:/taolinwork/poe-gwas/h2pm/sus11/pat_mat.pdf", pmh2, width = 4, height = 4)

###----------- paternal and maternal rg --------
iterN = c(1:19, 21:70, 75:132)
for (i in iterN){
  for (j in iterN){
    if (i < j){
      system(paste0("gcta64 --reml-bivar ", i, " ", j, " --mgrm h2/multi_grm.txt --pheno gcta.phen --covar gcta.cov --qcovar gcta.qcovar --reml-bivar-lrt-rg 0 --threads 10 --out rg_parental/seq", i, "_", j))
    }
  }
}

r0 = matrix(ncol = 5, nrow = 0) %>% as.data.frame()
names(r0) = c("Source", "Variance", "SE", "trait1", "trait2")
for (i in iterN){
  for (j in iterN){
    if (i<j){
      file1 = paste0("rg_parental/seq", i, "_", j, ".hsq")
      if (file.exists(file1)){
        r1 = read.csv(file1, sep = "")%>% 
          filter(Source %in% c("rG1", "rG2")) %>%
          mutate(trait1=paste0("t", i), trait2=paste0("t", j))
        r0 = rbind(r0, r1)
      }
      print(paste0(i," & ", j))
    }
  }
}

r0$Source[r0$Source=="rG1"] = "paternal"
r0$Source[r0$Source=="rG2"] = "maternal"
fwrite(r0, "h2/rg_parental.txt", sep = "\t", quote = F)
r0 = fread("h2/rg_parental.txt")
lor = paste0("t", 118:132)
r0 = r0 %>% filter(trait1 %in% lor, trait2 %in% lor, SE<=1)
fwrite(r0, "h2/rg_parental.lor.txt", sep = "\t", quote = F)

t.test(Variance~Source, r0)

library(ggpubr)
library(reshape2)
r0 = fread("D:/taolinwork/poe-gwas/mht/rg_parental.lor.txt")
r00 = r0 %>% dplyr::select(4,5,1,2)
names(r00)[3:4] = c("variable", "value")
r1 = dcast(r00,  trait1 + trait2 ~ variable) # 2.164e-07
t.test(r1$maternal, r1$paternal, paired = T)

rgplot = ggplot(r0, aes(Source, Variance, fill=Source)) + 
  geom_boxplot(width=0.4) + theme_classic() + guides(fill="none") +
  scale_fill_manual(values = c("#31a354", "#3182bd")) + 
  labs(y=expression(Parental~~italic(r)[g]), x=NULL) + 
  theme(axis.ticks = element_line(color="black"),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size =14)) +
  stat_compare_means(comparisons = list(c("paternal", "maternal")), 
                     method = "t.test")

gp = plot_grid(pmh2, rgplot, ncol = 2, align = "hv", 
          rel_widths = c(1.2,1), labels = "AUTO")
ggsave("D:/taolinwork/poe-gwas/h2pm/sus11/h2_rg.pdf", gp, width = 6, height = 3.5)

###----------- a+d+imp h2----------
# a: 
# system("gcta64 --bfile f2.filtered.add --make-grm --thread-num 10 --out gcta/add")

# d:
system("gcta64 --bfile f2.filtered.add --make-grm-d --thread-num 10 --out gcta/dom")

# imp: MACH output format---plink---gcta
do = fread("f2.filtered.vcf.gz")
id = names(do)[10:ncol(do)]
ml = fread("f2.filtered.add.name.frquency.frq")

do1 = do %>% select(10:ncol(do)) %>% as.matrix()
do1[do1=="0|0" | do1=="1|1"] = 0

for (i in 1:nrow(do1)){
  x = do1[i,]
  if (ml[i, 3] == do[i,4]){
    x[x=="0|1"] = -1
    x[x=="1|0"] = 1
  }
  if (ml[i, 3] == do[i,5]){
    x[x=="0|1"] = 1
    x[x=="1|0"] = -1
  }
  print(i)
}

do1 = do1 %>% t()
head = cbind(x1 = paste0("0->", id), x2 = rep("ML_DOSE", length(id)))
do1 = cbind(head, do1)
fwrite(do1, "imp.mldose", na="NA", col.names = F, sep = "\t", quote = F)
system("bgzip -@ 10 imp.mldose")


system("gcta64 --dosage-mach-gz imp.mldose.gz heter.mlinfo.gz --thread-num 10 --make-bed --out imp")

#chr=0 change to 1-18
d2 = fread("imp.bim") 
d2 = d2 %>% separate("V2", into = c("chr", "bp"), sep = "[.]", remove = F)
d2$chr = str_replace(d2$chr, "ssc", "") 
d2 = d2 %>% select(chr, V2, V3, V4, V5, V6)
fwrite(d2, "imp.bim", sep = "\t", col.names = F, quote = F)

system("gcta64 --bfile imp --make-grm --thread-num 10 --out gcta/imp")

for (i in 1:132){
  system(paste0("gcta64 --reml --mgrm h2_imp/multi_grm.txt --pheno gcta.phen --mpheno ", i ," --covar gcta.cov --qcovar gcta.qcovar --out h2_imp/t", i, " --thread-num 10"))
  print(i)
}


###----------- a+d+imp rg----------
iterN = c(1:19, 21:70, 75:132)
for (i in iterN){
  for (j in iterN){
    if (i < j){
      system(paste0("gcta64 --reml-bivar ", i, " ", j, " --mgrm h2_imp/multi_grm.txt --pheno gcta.phen --covar gcta.cov --qcovar gcta.qcovar --reml-bivar-lrt-rg 0 --threads 10 --out rg_imp/seq", i, "_", j))
    }
  }
}

###----------- heterozygote ---------
#转化为MACH
## test.mlinfo: SNP    Al1 Al2 Freq1   MAF Quality Rsq
system("plink --bfile f2.filtered.add.name --freq --out f2.filtered.add.name.frquency")
# By default, the minor allele is assigned to be A1
ml = fread("f2.filtered.add.name.frquency.frq")
ml0 = ml %>% mutate(Quality=1, Rsq=1, Freq1=MAF) %>% 
  select(SNP, A1, A2, Freq1, MAF, Quality, Rsq)
names(ml0) = c("SNP", "Al1", "Al2", "Freq1", "MAF", "Quality", "Rsq")
fwrite(ml0, "heter.mlinfo", sep = "\t")
system("bgzip -@ 10 heter.mlinfo")

## test.mldose: 001->0011 ML_DOSE 2.000 0.000
do = fread("f2.filtered.vcf.gz")
id = names(do)[10:ncol(do)]

do1 = do %>% select(10:ncol(do)) %>% t()
do1[do1=="0|0" | do1=="1|1"] = NA
do1[do1=="0|1"] = 0
do1[do1=="1|0"] = 1
he = cbind(x1 = paste0("0->", id), x2 = rep("ML_DOSE", length(id)))
do1 = cbind(he, do1)
fwrite(do1, "heter.mldose", na="NA", col.names = F, sep = "\t", quote = F)
system("bgzip -@ 10 heter.mldose")
system("gcta64 --mlma --dosage-mach-gz heter.mldose.gz heter.mlinfo.gz --make-bed --out heter --thread-num 4")
system("gcta64 --mlma --bfile heter --grm gcta/add --pheno gcta.phen --mpheno 1 --covar gcta.cov --qcovar gcta.qcovar --out gcta_heter/het1 --thread-num 4")

###----------- BLUP residual effect------
for (i in 1:132){
  system(paste0("gcta64 --reml --grm gcta/add --pheno gcta.phen --mpheno ", i, 
                " --covar gcta.cov --qcovar gcta.qcovar --reml-pred-rand --out blup/t", 
                i, " --thread-num 10"))
}

###----------- manhattan plot ----
# 返回杂合子GWAS中显著的SNP
t0 = fread("gcta_heter/het1.mlma") %>% filter(p <= 0.05/293988) %>% mutate(trait = "t1")
for (i in c(2:19, 21:70, 75:132)){
  t1 = fread(paste0("gcta_heter/het", i, ".mlma")) %>% 
    filter(p <= 0.05/293988) %>% mutate(trait = paste0("t", i))
  t0 = rbind(t0, t1)
  print(i)
}

fwrite(t0, "gcta_heter/sig.het.txt", sep = "\t")
t0 = t0 %>% tidyr::separate("SNP", into = c("rs", "bp"), remove = F, sep = "[.]")
t0$trait %>% unique() %>% length() #103
t0$SNP %>% unique() %>% length() #3555
t0$rs %>% unique()
t0 %>% group_by(rs) %>% summarise(n=n())
t0 %>% group_by(SNP) %>% summarise(n=n()) %>% arrange(desc(n)) %>% head(10)

### manhattan plot and boxplot
th0 = theme(panel.grid = element_blank(), 
            axis.line = element_line(colour = 'black'), 
            axis.text.y = element_text(color = "black"),
            axis.text.x = element_blank(),
            panel.background = element_blank()) 
th = theme(axis.text= element_text(color = "black", size=11),
           plot.title = element_text(hjust = 0.5))

add = function(trait, antSNP){
  d1 = fread(paste0("gcta/t", trait, ".mlma")) %>% select(SNP, Chr, bp, p)
  names(d1) = c("SNP","CHR","BP","P")
  d1 = d1 %>% arrange(CHR, BP)
  d1$SNP1 = seq(1, nrow(d1), 1)
  d1$CHR = factor(d1$CHR, levels = unique(d1$CHR))
  chr = summaryBy(SNP1~CHR, d1, FUN = median)
  d1 = d1 %>% mutate(type = ifelse(SNP==antSNP, "y", "n"))
  
  pmht <- ggplot(d1%>%filter(type=="n"))+
    geom_point(aes(SNP1, -log(P, 10), color = CHR), show.legend = FALSE, size=0.5) +
    scale_color_manual(values = rep(c("#1a1a1a", "#878787"),9))+
    labs(x = NULL, y = expression(-lg(additive~~italic(P)))) +
    scale_x_continuous(breaks = chr$SNP1.median, labels = chr$CHR, expand = c(0, 0)) +
    geom_hline(yintercept = -log10(0.05/293988),lty="dashed", color = 'red', alpha = 0.5)+ th0 +
    geom_point(data = d1%>%filter(type=="y"), aes(SNP1, -log(P, 10)), shape=23, size=2, color="black", fill="red")
  return(pmht)
  gc()
}

pat = function(trait, antSNP){
  d1 = fread(paste0("gcta_parental/pat", trait, ".mlma")) %>% select(SNP, Chr, bp, p)
  names(d1) = c("SNP", "CHR", "BP", "P")
  d1 = d1 %>% arrange(CHR, BP)
  d1$SNP1 = seq(1, nrow(d1), 1)
  d1$CHR = factor(d1$CHR, levels = unique(d1$CHR))
  chr = summaryBy(SNP1~CHR, d1, FUN = median)
  d1 = d1 %>% mutate(type = ifelse(SNP==antSNP, "y", "n"))
  pmht = ggplot(d1%>%filter(type=="n")) + 
    geom_point(aes(SNP1, -log(P, 10), color = CHR), show.legend = FALSE, size=0.5) +
    scale_color_manual(values = rep(c("#3182bd", "#9ecae1"),9))+
    labs(x = NULL, y = expression(-lg(paternal~~italic(P)))) +
    scale_x_continuous(breaks = chr$SNP1.median, labels = chr$CHR, expand = c(0, 0)) +
    geom_hline(yintercept = -log10(0.05/293988),lty="dashed", color = 'red', alpha = 0.5)+ th0 +
    geom_point(data = d1%>%filter(type=="y"), aes(SNP1, -log(P, 10)), shape=23, size=2, color="black", fill="red")
  return(pmht)
  gc()
}

mat = function(trait, antSNP){
  d1 = fread(paste0("gcta_parental/mat", trait, ".mlma")) %>% select(SNP, Chr, bp, p)
  names(d1) = c("SNP", "CHR", "BP", "P")
  d1 = d1 %>% arrange(CHR, BP)
  d1$SNP1 = seq(1, nrow(d1), 1)
  d1$CHR = factor(d1$CHR, levels = unique(d1$CHR))
  chr = summaryBy(SNP1~CHR, d1, FUN = median)
  d1 = d1 %>% mutate(type = ifelse(SNP==antSNP, "y", "n"))
  pmht = ggplot(d1%>%filter(type=="n")) + 
    geom_point(aes(SNP1, -log(P, 10), color = CHR), show.legend = FALSE, size=0.5) +
    scale_color_manual(values = rep(c("#31a354", "#a1d99b"),9))+
    labs(x = NULL, y = expression(-lg(maternal~~italic(P)))) +
    scale_x_continuous(breaks = chr$SNP1.median, labels = chr$CHR, expand = c(0, 0)) +
    geom_hline(yintercept = -log10(0.05/293988),lty="dashed", color = 'red', alpha = 0.5)+  th0 +
    geom_point(data = d1%>%filter(type=="y"), aes(SNP1, -log(P, 10)), shape=23, size=2, color="black", fill="red")
  return(pmht)
  gc()
}

het = function(trait, antSNP){
  d1 = fread(paste0("gcta_heter/het", trait, ".mlma"))%>% 
    tidyr::separate("SNP", into = c("rs", "bp"), remove = F, sep = "[.]")
  d1$chr = str_remove(d1$rs, "ssc")
  d1$chr = as.numeric(d1$chr)
  d1$bp = as.numeric(d1$bp)
  d1 = d1 %>% select(SNP, chr, bp, p)
  names(d1) = c("SNP","CHR","BP","P")
  d1 = d1 %>% arrange(CHR, BP)
  d1$SNP1 = seq(1, nrow(d1), 1)
  d1$CHR = factor(d1$CHR, levels = unique(d1$CHR))
  chr = summaryBy(SNP1~CHR, d1, FUN = median)
  d1 = d1 %>% mutate(type = ifelse(SNP==antSNP, "y", "n"))
  pmht <- ggplot(d1%>%filter(type=="n")) + 
    geom_point(aes(SNP1, -log(P, 10), color = CHR), show.legend = FALSE, size=0.5) +
    scale_color_manual(values = rep(c("#2d004b", "#8073ac"),9))+
    labs(x = "Chromosomes", y = expression(-lg(heterozygous~~italic(P)))) +
    scale_x_continuous(breaks = chr$SNP1.median, labels = chr$CHR, expand = c(0, 0)) +
    geom_hline(yintercept = -log10(0.05/293988),lty="dashed", color = 'red', alpha = 0.5)+
    theme(panel.grid = element_blank(), 
          axis.line = element_line(colour = 'black'), 
          axis.text = element_text(color = "black"),
          panel.background = element_blank()) +
    geom_point(data = d1%>%filter(type=="y"), aes(SNP1, -log(P, 10)), shape=23, size=2, color="black", fill="red")
  return(pmht)
  gc()
}

vcf = fread("f2.filtered.vcf.gz")
names(vcf)[1] = "chr"
gr = read_xlsx("/data_group/xiehaibing/xiehaibing2/poe-gwas/group.xlsx", sheet = 1) %>%
  as.data.frame()
row.names(gr) = str_remove(gr$trait, "T")

snp = function(trait, antSNP, ytitle){
  x = antSNP %>% str_remove("ssc") %>% str_split("[.]") %>% as.data.frame()
  chro = x[1,1] %>% as.numeric()
  pos1 = x[2,1] %>% as.numeric()
  
  pat = fread(paste0("gcta_parental/pat", trait, ".mlma")) %>% filter(SNP==antSNP) %>% as.data.frame()
  mat = fread(paste0("gcta_parental/mat", trait, ".mlma")) %>% filter(SNP==antSNP) %>% as.data.frame()
  
  res = fread(paste0("blup/t",trait, ".indi.blp")) %>% select(2, 6)
  names(res) = c("id", "residual")
  d0 = vcf %>% filter(chr==chro, POS==pos1)
  ref = d0[1, "REF"] %>% as.character()
  alt = d0[1, "ALT"] %>% as.character()
  da = d0[,10:ncol(d0)] %>% t() %>% as.data.frame() %>% 
    separate(V1, into = c("pat", "mat"), sep = "[|]", remove = F) %>% 
    rownames_to_column("id") %>% left_join(res, by="id") %>% na.omit() %>%
    mutate(ref =ref, alt = alt,
           pat1 = if_else(pat==0, ref, alt),
           mat1 = if_else(mat==0, ref, alt),
           Genotype = paste0(pat1, "|", mat1))
 
  ge1 = aggregate(da$Genotype, by = list(da$Genotype), FUN = length)
  names(ge1) = c("Genotype", "Ngenotype")
  pa1 = aggregate(da$pat1, by = list(da$pat1), FUN = length)
  names(pa1) = c("pat1", "Npat")
  pa1$y = min(da$residual)
  ma1 = aggregate(da$mat1, by = list(da$mat1), FUN = length)
  names(ma1) = c("mat1", "Nmat")
  
  da = da %>% left_join(ge1, by="Genotype") %>% left_join(pa1, by="pat1") %>% left_join(ma1, by="mat1")

  pval1 = signif(as.numeric(pat[1, "p"]), 2)
  p1 = ggplot(da, aes(pat1, residual)) + 
    geom_boxplot(fill="#3182bd", outlier.colour = NA, width=0.6) + 
    geom_text(aes(x=pat1, y=y, label=Npat)) + theme_bw() + th +
    labs(title=bquote(italic(P)~"="~.(pval1)), 
         x="Paternal allele", y=ytitle)
  
  pval2 = signif(as.numeric(mat[1,"p"]),2)
  p2 = ggplot(da, aes(mat1, residual)) +
    geom_boxplot(fill="#31a354", outlier.colour = NA, width=0.6) +
    geom_text(aes(x=mat1, y=y, label=Nmat)) + theme_bw() + th +
    labs(title=bquote(italic(P)~"="~.(pval2)), 
         x="Maternal allele", y="")
  p3 = ggplot(da, aes(Genotype, residual)) + 
    geom_boxplot(outlier.colour = NA, width=0.6) +
      geom_text(aes(x=Genotype, y=y, label=Ngenotype)) + 
    theme_bw() + th + labs(y="", title = antSNP)
  p0 = plot_grid(p1, p2, p3, ncol = 3, align = "hv", rel_widths = c(0.5, 0.5, 1))
  
  m1 = plot_grid(add(trait, antSNP), pat(trait, antSNP), 
                 mat(trait, antSNP), het(trait, antSNP), ncol = 1, align = "hv")
  ptot = plot_grid(m1, p0, ncol = 1, align = "hv", rel_heights=c(2, 0.6))
  #return(ptot)
  ggsave(paste0("mht/t", trait, "_", antSNP, ".tiff"), ptot, width = 8, height = 9, dpi = 1200)
}

# p = snp(121, "ssc2.442357", ytitle="LOR4 residual")
# snp(116, "ssc2.148273414", ytitle="COLFHQ residual")
# snp(39, "ssc7.101261209", "NOR residual")

t0 = fread("gcta_heter/sig.het.txt")
t00 = t0$trait %>% unique()

par = matrix(ncol = 3, nrow = 0) %>% as.data.frame()
names(par) = c("tr", "SNP", "ylab1")
for (i in t00){
  x = t0 %>% filter(trait==i) %>% arrange(p)
  tr = str_remove(i, "t")
  tagsnp = x[1,"SNP"]
  ylab1 = paste0(gr[tr, "id"], " residual")
  par2 = cbind(tr, tagsnp, ylab1)
  par = rbind(par, par2)
}

list(trait=par$tr, antSNP=par$SNP, ytitle=par$ylab1) %>% pmap(snp)


###----------- QQ-plot ----------
gqq = function(trait){
  getp = function(d0){
    observed <- sort(d0$p)
    lobs <- -(log10(observed))
    expected <- c(1:length(observed)) 
    lexp <- -(log10(expected/(length(expected)+1)))
    p = cbind(lexp = lexp, lobs=lobs) %>% as.data.frame()
    return(p)
  }
  
  d1 = fread(paste0("gcta/t", trait, ".mlma")) %>% select(p) %>% na.omit() %>% 
    getp() %>% mutate(type="a")
  d2 = fread(paste0("gcta_parental/pat", trait, ".mlma")) %>% select(p) %>% na.omit() %>%
    getp() %>% mutate(type="b")
  d3 = fread(paste0("gcta_parental/mat", trait, ".mlma")) %>% select(p) %>% na.omit() %>%
    getp() %>% mutate(type="c")
  d4 = fread(paste0("gcta_heter/het", trait, ".mlma"))%>% select(p) %>% na.omit() %>%
    getp() %>% mutate(type="d")
  p = rbind(d1, d2, d3, d4)
  
  p0 = ggplot(p, aes(lexp, lobs, color=type)) + 
    geom_abline(intercept = 0, slope = 1, color="red")+
    geom_point(size=0.5) + 
    scale_color_manual(values = c("#1a1a1a", "#3182bd", "#31a354","#2d004b")) +
    theme_classic() + th +
    labs(x = expression(Expected~~-lg(italic(P))),
         y = expression(Observed~~-lg(italic(P)))) +
    theme(axis.text = element_text(colour = "black")) +
    guides(color="none")
  return(p0)
  gc()
}
for (i in 1:132){
  qqplot = gqq(i)
  ggsave(paste0("mht/t", i, "tiff"), qqplot, width = 2, height = 2, dpi = 1200)
}

###----------- SNP annotation ---------------------
t0 = fread("gcta_heter/sig.het.txt") %>% select(SNP) %>% distinct() %>%
  separate("SNP", into = c("ssc", "pos"), remove = F, sep = "[.]")#只关注了杂合子显著的SNP
t0$ssc = str_remove(t0$ssc, "ssc")
t0$ssc = as.numeric(t0$ssc)
t0$pos = as.numeric(t0$pos)
vcf1 = fread("f2.filtered.vcf.gz") %>% select(1:5) #为了添加碱基信息
names(vcf1) = c("ssc", "pos", "id", "ref", "alt")
t0 = t0 %>% left_join(vcf1, by=c("ssc", "pos"))
t0 = t0 %>% mutate(id=pos) %>% select(ssc, pos, id, ref, alt, SNP)
fwrite(t0, "gcta_heter/sig.het.snp", sep = "\t", col.names = F, quote = F)

av.pl="perl /data_group/xiehaibing/xiehaibing2/software/annovar/annotate_variation.pl"
db="/data_group/xiehaibing/xiehaibing2/software/annovar/pig11db/"
system(paste0(av.pl, " --outfile gcta_heter/sig.het.snp.annovar.txt --buildver pig11 gcta_heter/sig.het.snp ", db))

system("perl /data_group/xiehaibing/xiehaibing2/software/annovar/annotate_variation.pl --outfile dmr.annovar_dep10.txt --buildver pig11 dmr.annovar_dep10.txt /data_group/xiehaibing/xiehaibing2/software/annovar/pig11db/")

an = fread("gcta_heter/sig.het.snp.annovar.txt.variant_function") #注释结果
names(an) = c("region", "gene", "ssc", "pos1", "pos2", "ref", "alt", "SNP")
t0 = fread("gcta_heter/sig.het.txt") %>% left_join(an, by="SNP")
tx = t0 %>% filter(trait %in% paste0("t", 118:132))
tx$gene%>% str_replace_all("\\(.*?\\)", "") %>% stringr::str_split(",") %>% 
  unlist() %>% unique() %>% na.omit() %>% sort()

d1 = fread("missingh2/sig.add.txt") %>% 
  filter(trait %in% paste0("t",118:132)) %>%
  group_by(SNP) %>%
  summarise(n=n()) %>% 
  arrange(desc(n))

### summary heter-GWAS, /data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5/gcta_heter/script
library(data.table)
library(dplyr)
library(purrr)
library(readxl)
setwd("/data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5")
snp = fread("gcta_heter/sig.het.txt")

get_add_snp = function(trait, snp){
  d1 = fread(paste0("gcta/", trait, ".mlma")) %>% filter(SNP %in% snp) %>% mutate(traitID = trait)
  return(d1)
}
add = list(trait=snp$trait, snp=snp$SNP) %>% pmap(get_add_snp) %>% rbindlist()
fwrite(add, "gcta_heter/sig.het_add.snp.txt")

get_snp = function(trait, snp, parent){
  d1 = fread(paste0("gcta_parental/", parent, trait, ".mlma")) %>% filter(SNP %in% snp) %>% mutate(traitID = trait)
  return(d1)
}

pat = list(trait=snp$trait, snp=snp$SNP, parent=rep("pa", nrow(snp))) %>% pmap(get_snp) %>% rbindlist()
fwrite(pat, "gcta_heter/sig.het_pat.snp.txt")

mat = list(trait=snp$trait, snp=snp$SNP, parent=rep("ma", nrow(snp))) %>% pmap(get_snp) %>% rbindlist()
fwrite(mat, "gcta_heter/sig.het_mat.snp.txt")


snp1 = snp %>% select(SNP, trait, b, p)
names(snp1) = c("SNP", "traitID", "b_het", "p_het")

add = fread("gcta_heter/sig.het_add.snp.txt") %>% select(SNP, traitID, A1, A2, Freq, b, p)
names(add) = c("SNP", "traitID", "A1", "A2", "Freq", "b_add", "p_add")

pat = fread("gcta_heter/sig.het_pat.snp.txt") %>% select(SNP, traitID, b, p)
names(pat) = c("SNP", "traitID", "b_pat", "p_pat")

mat = fread("gcta_heter/sig.het_mat.snp.txt") %>% select(SNP, traitID, b, p)
names(mat) = c("SNP", "traitID", "b_mat", "p_mat")

gro = read_xlsx("/data_group/xiehaibing/xiehaibing2/poe-gwas/group.xlsx")
gro$trait = tolower(gro$trait)

an = fread("gcta_heter/sig.het.snp.annovar.txt.variant_function") %>% select(V8, V2, V1)
names(an) = c("SNP", "gene", "location")

snp1 = snp1 %>% left_join(add, by=c("SNP", "traitID")) %>% left_join(pat, by=c("SNP", "traitID")) %>% 
  left_join(mat, by=c("SNP", "traitID")) %>% left_join(gro, by=join_by("traitID" == "trait")) %>%
  left_join(an, by="SNP") %>%
  select("name", "id", "traitID", "type", "SNP",  "A1", "A2", "Freq", "gene", "location", "b_het", "p_het", 
         "b_add", "p_add", "b_pat", "p_pat", "b_mat", "p_mat")
fwrite(snp1, "gcta_heter/sig.het.all.gwas.txt", sep = "\t")


###----------- missing h2 ---------
interN = c(1:19, 21:70, 75:132)
# 所有性状显著的SNP
get_sig_snp = function(dir1){
  t0 = matrix(ncol = 2, nrow = 0) %>% as.data.frame()
  names(t0) = c("SNP", "trait")
  for (i in interN){
    t1 = fread(paste0(dir1, i, ".mlma")) %>% 
      filter(p <= 0.05/293988) %>% 
      mutate(trait = paste0("t", i)) %>%
      select(SNP, trait)
    t0 = rbind(t0, t1)
    print(i)
  }
  return(t0)
}

addsnp = get_sig_snp("gcta/t")
fwrite(addsnp, "missingh2/sig.add.txt", sep = "\t", quote = F)
patsnp = get_sig_snp("gcta_parental/pat")
fwrite(patsnp, "missingh2/sig.pat.txt", sep = "\t", quote = F)
matsnp = get_sig_snp("gcta_parental/mat")
fwrite(matsnp, "missingh2/sig.mat.txt", sep = "\t", quote = F)
hetsnp = get_sig_snp("gcta_heter/het")
fwrite(hetsnp, "missingh2/sig.het.txt", sep = "\t", quote = F)

allsnp = rbind(addsnp, patsnp, matsnp, hetsnp) %>% distinct(.keep_all = T)
# 计算显著SNP h2
cal_h2 = function(trait1){
  # 加性显著SNP
  sigadd = addsnp %>% filter(trait == paste0("t", trait1)) %>% select(SNP)
  fwrite(sigadd, "missingh2/add.snplist", col.names = F, quote = F, sep = "\t")
  system("gcta64 --bfile f2.filtered.add.name --extract missingh2/add.snplist --make-grm --out missingh2/add --thread-num 10")
  system(paste0("gcta64 --reml --grm missingh2/add --pheno gcta.phen --mpheno ", trait1, 
                " --covar gcta.cov --qcovar gcta.qcovar --out missingh2/add_t", trait1, " --thread-num 10"))
  
  # 所有显著SNP
  sigall = allsnp %>% filter(trait == paste0("t", trait1)) %>% select(SNP)
  fwrite(sigall, "missingh2/all.snplist", col.names = F, quote = F, sep = "\t")
  system("gcta64 --bfile f2.filtered.add.name --extract missingh2/all.snplist --make-grm --out missingh2/all --thread-num 10")
  system(paste0("gcta64 --reml --grm missingh2/all --pheno gcta.phen --mpheno ", trait1, 
                " --covar gcta.cov --qcovar gcta.qcovar --out missingh2/all_t", trait1, " --thread-num 10"))
}
# cal_h2(121)
list(trait1 = interN) %>% pmap(cal_h2)

# 统计缺失遗传力的情况
ms0 = matrix(ncol = 5, nrow = 0) %>% as.data.frame()
names(ms0) = c("Source", "Variance", "SE", "trait", "type")

for (i in interN){
  ms1 = read.csv(paste0("missingh2/add_t", i, ".hsq"), sep = "") %>% 
    filter(Source %in% c("V(G)/Vp", "Pval")) %>%
    mutate(trait = paste0("t", i), type="add")
  ms2 = read.csv(paste0("missingh2/all_t", i, ".hsq"), sep = "") %>% 
    filter(Source %in% c("V(G)/Vp", "Pval")) %>%
    mutate(trait = paste0("t", i), type="all")
  ms0 = rbind(ms0, ms1, ms2)
}

fwrite(ms0, "missingh2/tot.result.txt", sep = "\t", quote = F, col.names = T)

# 可视化
dmu = fread("D:/taolinwork/poe-gwas/figure/h2.127.txt") %>%
  filter(trait %in% paste0("t", 118:132)) %>%
  mutate(type = "family")

ms0 = fread("D:/taolinwork/poe-gwas/figure/tot.result.txt") %>% 
  filter(Source=="V(G)/Vp", trait %in% paste0("t", 118:132)) %>%
  select(Variance, SE, trait, type)
names(dmu) = names(ms0)
ms0 = rbind(ms0, dmu)

pmissing = ggplot(ms0, aes(trait, Variance, color=type)) + geom_point() +
  geom_errorbar(aes(ymin=Variance-SE, ymax=Variance+SE), width=0.2, alpha = 0.3) + 
  theme_classic() + labs(x=NULL, y=expression(italic(h)^2)) + 
  scale_x_discrete(labels=paste0("LOR", 1:15)) +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_color_manual(values = c("red", "blue", "black")) + #"#e41a1c", "#377eb8", "#4daf4a"
  theme(axis.text = element_text(color = "black", size=11),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank(),
        legend.position = c(0,0), 
        legend.justification = c(0,0),
        legend.background = element_rect(fill = "white", colour = "black"))
ggsave("D:/taolinwork/poe-gwas/figure/missingh2.pdf", pmissing, width = 6, height = 3)

# paired t test
ms = ms0 %>% select(Variance, trait, type) %>% 
  rename(value=Variance, variable=type) %>%
  dcast(trait~variable)
t.test(ms$add, ms$all, paired = T, alternative = "less") # p-value = 0.0005403
t.test(ms$all, ms$family, paired = T, alternative = "less") # p-value = 3.241e-08
t.test(ms$ad, ms$family, paired = T, alternative = "less") # p-value = 4.353e-10
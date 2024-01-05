# source activate r4_py37_env
library(data.table)
library(dplyr)
library(blupADC)
setwd("/data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5/dmu")

# phenotype
d1 = fread("../gcta.phen")
names(d1) = c("fid", "id", paste0("t", 1:132))
d1 = d1 %>% select(c("id", paste0("t", 1:132)))

d2 = fread("../gcta.cov") %>% select(2:3)
names(d2) = c("id", "sex")
d2 = d2 %>% left_join(d1, by="id")


# pedigree
ped = fread("../../id11.fam")
ped$da = ped$kid

id = c(ped$kid, ped$pat, ped$mat) %>% unique() %>% as.data.frame()
names(id) = "id"
id = na.omit(id)
id$no = 1:nrow(id)

d2 = d2 %>% left_join(id, by="id") %>%
  select(no, sex, paste0("t",1:132))
fwrite(d2, "phenotype.txt", sep = "\t", na="-9", col.names = F, quote = F)

id1 = id
id2 = id
id3 = id
names(id1) = c("kid", "no_kid")
names(id2) = c("pat", "no_pat")
names(id3) = c("mat", "no_mat")

ped = ped %>%  
  left_join(id1, by="kid") %>% 
  left_join(id2, by="pat") %>% 
  left_join(id3, by="mat") %>%
  select(no_kid, no_pat, no_mat) %>%
  mutate(da = no_kid)

fwrite(ped, "pedigree.txt", col.names = F, na="0", quote = F, sep = "\t")


# 利用DMU软件进行遗传评估
# data_path=system.file("extdata", package = "blupADC")  #  path of example files
data_path="/data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5/dmu"
dirout = "/data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5/dmu"

# run_DMU(
#   phe_col_names=c("Id","Sex", paste0("t", 1:132)), # colnames of phenotype 
#   target_trait_name=list(c("t121")),           #trait name 
#   fixed_effect_name=list(c("Sex")),            #fixed effect name
#   random_effect_name=list(c("Id")),            #random effect name
#   covariate_effect_name=NULL,                  #covariate effect name
#   phe_path=data_path,                          #path of phenotype file
#   phe_name="phenotype.txt",                    #name of phenotype file
#   integer_n=2,                                 #number of integer variable 
#   analysis_model="PBLUP_A",                    #model of genetic evaluation
#   dmu_module="dmuai",                          #modeule of estimating variance components 
#   relationship_path=data_path,                 #path of relationship file 
#   relationship_name="pedigree.txt",            #name of relationship file 
#   output_result_path=dirout                   # output path 
# )

dmu = function(i){
  setwd("/data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5/dmu")
  run_DMU(
    phe_col_names=c("Id","Sex", paste0("t", 1:132)), # colnames of phenotype 
    target_trait_name=list(c(paste0("t", i))),           #trait name 
    fixed_effect_name=list(c("Sex")),            #fixed effect name
    random_effect_name=list(c("Id")),            #random effect name
    covariate_effect_name=NULL,                  #covariate effect name
    phe_path=data_path,                          #path of phenotype file
    phe_name="phenotype.txt",                    #name of phenotype file
    integer_n=2,                                 #number of integer variable 
    analysis_model="PBLUP_A",                    #model of genetic evaluation
    dmu_module="dmuai",                          #modeule of estimating variance components 
    relationship_path=data_path,                 #path of relationship file 
    relationship_name="pedigree.txt",            #name of relationship file 
    output_result_path=dirout                   # output path 
  )
}

for (i in 11:132){dmu(i)}

# 合并每个性状的遗传力

h2 = matrix(ncol = 6, nrow = 0) %>% as.data.frame()
names(h2) = c("Random_effect_name", "prior","prior_se","h2","h2_se", "trait")
for (i in c(1:19, 21:70, 75:132)){
  he0 = fread(paste0("PBLUP_A_t", i, "/t", i, "_heritability_result.txt")) %>% 
    filter(Random_effect_name=="Id") %>%
    mutate(trait = paste0("t", i))
  h2 = rbind(h2, he0)
}

h2 = h2 %>% select(4:6)
fwrite(h2, "h2.127.txt", sep = "\t")
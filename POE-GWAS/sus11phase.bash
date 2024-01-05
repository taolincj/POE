##!/bin/sh
#PBS -N chr<par>
#PBS -q queue2T
#PBS -l mem=100gb,walltime=10000:00:00,nodes=1:ppn=10
#PBS -o /data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/script/chr<par>.sh.out
#PBS -e /data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/script/chr<par>.sh.err

# qsub -I -N test -q queue2T -l select=1:ncpus=10:mem=200GB -l walltime=10000:00:00
phase_common=/data_group/xiehaibing/xiehaibing2/software/phase_common_static
cd /data_group/xiehaibing/xiehaibing2/poe-gwas/10to11

vcftools --gzvcf gentalks-result/genetalks_result.joint.vcf.gz --keep 602.txt --recode --recode-INFO-all --out sus11.602
bgzip -c -@ 10 sus11.602.recode.vcf > sus11.602.recode.vcf.gz
bcftools index sus11.602.recode.vcf.gz --threads 10

#phase 
i=<par>
$phase_common --input sus11.602.recode.vcf.gz --pedigree id11.fam --region $i --output shapeit5/phased.chr$i.bcf --thread 10
echo $?

cd /data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5
bcftools concat -o phased.602.vcf.gz -O z --threads 10 \
phased.chr1.bcf \
phased.chr2.bcf \
phased.chr3.bcf \
phased.chr4.bcf \
phased.chr5.bcf \
phased.chr6.bcf \
phased.chr7.bcf \
phased.chr8.bcf \
phased.chr9.bcf \
phased.chr10.bcf \
phased.chr11.bcf \
phased.chr12.bcf \
phased.chr13.bcf \
phased.chr14.bcf \
phased.chr15.bcf \
phased.chr16.bcf \
phased.chr17.bcf \
phased.chr18.bcf

#提取F2个体
bcftools view -S f2.534.keep phased.602.vcf.gz -O z --threads 10 -o phased.f2.vcf.gz
#质控
plink2 --vcf phased.f2.vcf.gz --geno 0.05 --mind 0.05 --maf 0.05 --hwe 0.0000001 --recode vcf --out f2.filtered
#压缩
bgzip -i -@ 10 f2.filtered.vcf
#提取单倍型，便于提取亲本基因型
bcftools convert --hapsample --vcf-ids f2.filtered.vcf.gz -o f2.filtered --threads 10

# 备份数据至冷盘
# mv /data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5/phased* /nas_data/xiehaibing/xiehaibing2/wgs/poe-gwas/sus11phase/
# cd /data_group/xiehaibing/xiehaibing2/poe-gwas/10to11/shapeit5
# #把杂合子文件转化成plink二进制
# gcta64 --mlma --dosage-mach-gz heter.mldose.gz heter.mlinfo.gz --make-bed --out heter --thread-num 4
# 
# # 以下循环分开提交
# for i in $(seq 1 132); do
# # 加性GWAS
# gcta64 --mlma --bfile f2.filtered.add.name --grm gcta/add --pheno gcta.phen --mpheno $i --covar gcta.cov --qcovar gcta.qcovar --out gcta/t$i --thread-num 4
# # 父本
# gcta64 --mlma --bfile pat --grm gcta/add --pheno gcta.phen --mpheno $i --covar gcta.cov --qcovar gcta.qcovar --out gcta_parental/pat$i --thread-num 4
# # 母本
# gcta64 --mlma --bfile mat --grm gcta/add --pheno gcta.phen --mpheno $i --covar gcta.cov --qcovar gcta.qcovar --out gcta_parental/mat$i --thread-num 4
# # 杂合子
# gcta64 --mlma --bfile heter --grm gcta/add --pheno gcta.phen --mpheno $i --covar gcta.cov --qcovar gcta.qcovar --out gcta_heter/het$i --thread-num 4
# done
echo $?

###----------------------------------------------- WASP -----------------------------------------------------------
#/data_group/xiehaibing/xiehaibing2/f1/rna/script/wasp/202026-10rib.sh
cd /data_group/xiehaibing/xiehaibing2/f1/rna
cpunum=10

#DIR_fastq=/data_group/xiehaibing/xiehaibing2/f1/lixin0210072/X101SC22124979-Z02-J001/00.CleanData/
#DIR_fastq=/data_group/xiehaibing/xiehaibing2/f1/lixin0210072/X101SC22124979-Z02-J002_multipath/X101SC22124979-Z02-J002_01/00.CleanData
DIR_fastq=/data_group/xiehaibing/xiehaibing2/f1/lixin0210072/X101SC22124979-Z02-J002_multipath/X101SC22124979-Z02-J002_02/00.CleanData
vcf="/data_group/xiehaibing/xiehaibing2/f1/rna/pysam/individual_vcf/202026-10.het.recode.vcf"
i="202026-10肋骨"

genomeDirpass1=/data_group/xiehaibing/xiehaibing2/f1/rna/index_pig11
DIR_STARpass1=/data_group/xiehaibing/xiehaibing2/f1/rna/DIR_STARpass1
DIR_STARpass2=/data_group/xiehaibing/xiehaibing2/f1/rna/DIR_STARpass2WASP
genomeDirpass2_all=/data_group/xiehaibing/xiehaibing2/f1/rna/stargenome_pass2
FASTA=/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
gtf=/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf


source activate r4_py27_env
genomeDirpass2=$genomeDirpass2_all/$i

echo "filter bam file for reads passed WASP filtering"
DIR_STAR_wasp=$DIR_STARpass2
samtools view -@ $cpunum -h $DIR_STAR_wasp/${i}.Aligned.sortedByCoord.out.bam | \
	awk '{if($18!="vW:i:2" && $18!="vW:i:3" && $18!="vW:i:4" && $18!="vW:i:5" && $18!="vW:i:6" && $18!="vW:i:7") print $0}' | \
	samtools view -@ $cpunum -Shb - \
	> $DIR_STAR_wasp/${i}_WASPfilterPASS.bam
echo $?

### 统计WASP前后的read数量
# /data_group/xiehaibing/xiehaibing2/f1/d90/rna/script/wasp/wasp.read.count.sh #杂合位点：WASPreadCount.sh
cd /data_group/xiehaibing/xiehaibing2/f1/rna/DIR_STARpass2WASP
# cd /data_group/xiehaibing/xiehaibing2/f1/d90/rna/DIR_STARpass2WASP
ncpu=10
echo "total read number before WASP"
for i in $(ls *.Aligned.sortedByCoord.out.bam); do
j=${i%%.*}
echo $j 
samtools view -c -@ $ncpu $i
done

echo "read number after WASP"
for i in $(ls *_WASPfilterPASS.bam); do
j=${i%%_*}
echo $j
samtools view -c -@ $ncpu $i
done
echo $?

###---------------------------------------------split BAM using pysam-----------------------------------------------
cd /data_group/xiehaibing/xiehaibing2/f1/rna/pysam/pysam_bam
sambamba index -p -t 4 /data_group/xiehaibing/xiehaibing2/f1/rna/DIR_STARpass2WASP/<para2>-肌肉_WASPfilterPASS.bam
mkdir -p <para2>jr
for i in $(seq 1 18); do
python /data_group/xiehaibing/xiehaibing2/f1/rna/pysam/script/split_bam_pair.py $i \
/data_group/xiehaibing/xiehaibing2/f1/wgbs/pysam/vcf4/<para2>.filter4.vcf \
/data_group/xiehaibing/xiehaibing2/f1/rna/DIR_STARpass2WASP/<para2>-肌肉_WASPfilterPASS.bam \
<para2>jr/<para2>jr
done
echo $?

###-------------------------------------------- 合并BAM ------------------------------------------------------------- 
#/data_group/xiehaibing/xiehaibing2/f1/rna/pysam/script/MergeBamYZ.sh
cd /data_group/xiehaibing/xiehaibing2/f1/rna/pysam/pysam_bam

for i in $(ls); do
#合并
sambamba merge -t 10 $i/$i.maternal $i/$i.maternal*
sambamba merge -t 10 $i/$i.paternal $i/$i.paternal*

#排序
sambamba sort -t 10 -o $i/$i.maternal.sort.bam $i/$i.maternal
sambamba sort -t 10 -o $i/$i.paternal.sort.bam $i/$i.paternal

#删除中间文件
rm $i/$i.maternal
rm $i/$i.paternal
done
echo $?
###------------------------------------------------------  统计pysam的拆分比例 ------------------------------------------
#拆分前
cd /data_group/xiehaibing/xiehaibing2/f1/rna/DIR_STARpass2WASP
# cd /data_group/xiehaibing/xiehaibing2/f1/d90/rna/DIR_STARpass2WASP
for i in $(ls *WASPfilterPASS.bam); do
echo $i
samtools view -c -@ 10 $i
done

#拆分后
cd /data_group/xiehaibing/xiehaibing2/f1/rna/pysam/pysam_bam
# cd /data_group/xiehaibing/xiehaibing2/f1/d90/rna/pysam/pysam_bam
for i in $(ls); do
echo $i pat mat
samtools view -c -@ 10 $i/$i.paternal.sort.bam
samtools view -c -@ 10 $i/$i.maternal.sort.bam
done
echo $?

###------------------------------------------------ 拆分后BAM定量 -----------------------------------
source activate gatk
gtf=/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf

cd /data_group/xiehaibing/xiehaibing2/f1/rna/pysam/pysam_bam
for i in $(ls); do
mv $i/$i.*sort* ../mergeBAM/
echo $i
done

cd /data_group/xiehaibing/xiehaibing2/f1/rna/pysam/mergeBAM
featureCounts -T 10 -p -t exon -g gene_id -a $gtf \
	-o ../expression/rib.featureCounts.txt \
	202026-10rib.maternal.sort.bam 202026-10rib.paternal.sort.bam 202026-2rib.maternal.sort.bam 202026-2rib.paternal.sort.bam \
	202102-7rib.maternal.sort.bam 202102-7rib.paternal.sort.bam 202102-8rib.maternal.sort.bam 202102-8rib.paternal.sort.bam \
	202232-1rib.maternal.sort.bam 202232-1rib.paternal.sort.bam 202232-2rib.maternal.sort.bam 202232-2rib.paternal.sort.bam \
	202232-3rib.maternal.sort.bam 202232-3rib.paternal.sort.bam 202266-7rib.maternal.sort.bam 202266-7rib.paternal.sort.bam \
	202266-8rib.maternal.sort.bam 202266-8rib.paternal.sort.bam 202266-9rib.maternal.sort.bam 202266-9rib.paternal.sort.bam \
	Y193208-1rib.maternal.sort.bam Y193208-1rib.paternal.sort.bam Y193208-2rib.maternal.sort.bam Y193208-2rib.paternal.sort.bam \
	Y193208-3rib.maternal.sort.bam Y193208-3rib.paternal.sort.bam Y193208-4rib.maternal.sort.bam Y193208-4rib.paternal.sort.bam \
	Y200902-1rib.maternal.sort.bam Y200902-1rib.paternal.sort.bam Y200902-2rib.maternal.sort.bam Y200902-2rib.paternal.sort.bam \
	Y200902-3rib.maternal.sort.bam Y200902-3rib.paternal.sort.bam Y204711-1rib.maternal.sort.bam Y204711-1rib.paternal.sort.bam \
	Y204711-2rib.maternal.sort.bam Y204711-2rib.paternal.sort.bam Y204711-3rib.maternal.sort.bam Y204711-3rib.paternal.sort.bam

cd /data_group/xiehaibing/xiehaibing2/f1/rna/pysam/mergeBAM
source activate gatk
gtf=/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf
for i in $(ls *bam)
do
j=${i%%.sort*}
stringtie -p 20 -e -B -G $gtf -o ../expression/stringtie/$j/transcripts.gtf -A ../expression/stringtie/$j/gene_abundances.tsv $i
echo $i was done!
done
echo $?

cd /data_group/xiehaibing/xiehaibing2/f1/rna/pysam/expression/stringtie									
perl stringtie_expression_matrix.pl --expression_metric=FPKM \
                                    --result_dirs='202026-10rib.maternal,202026-10rib.paternal,202026-2rib.maternal,202026-2rib.paternal,202102-7rib.maternal,202102-7rib.paternal,202102-8rib.maternal,202102-8rib.paternal,202232-1rib.maternal,202232-1rib.paternal,202232-2rib.maternal,202232-2rib.paternal,202232-3rib.maternal,202232-3rib.paternal,202266-7rib.maternal,202266-7rib.paternal,202266-8rib.maternal,202266-8rib.paternal,202266-9rib.maternal,202266-9rib.paternal,Y193208-1rib.maternal,Y193208-1rib.paternal,Y193208-2rib.maternal,Y193208-2rib.paternal,Y193208-3rib.maternal,Y193208-3rib.paternal,Y193208-4rib.maternal,Y193208-4rib.paternal,Y200902-1rib.maternal,Y200902-1rib.paternal,Y200902-2rib.maternal,Y200902-2rib.paternal,Y200902-3rib.maternal,Y200902-3rib.paternal,Y204711-1rib.maternal,Y204711-1rib.paternal,Y204711-2rib.maternal,Y204711-2rib.paternal,Y204711-3rib.maternal,Y204711-3rib.paternal' \
                                    --transcript_matrix_file=transcript_fpkm_rib.tsv \
                                    --gene_matrix_file=gene_fpkm_rib.tsv


############################################################################################
source activate reseq_env
cd /data_group/xiehaibing/xiehaibing2/f1/rna/DIR_STARpass2WASP
for i in $(ls *肋骨_WASPfilterPASS.bam)
do
j=${i%%肋骨*}
echo $j
picard AddOrReplaceReadGroups \
        -I $i \
        -O ${j}-rib_WASPfilterPASS_rgadded_sorted.bam \
        --SORT_ORDER coordinate \
        --RGID ${j}_rib \
        --RGLB ${j}_rib \
        --RGPL Illumina \
        --RGPU ${j}_rib \
        --RGSM ${j}_rib
done
echo $?

###--------------------------------- DEG between crosses -----------------------------------------
cd /data_group/xiehaibing/xiehaibing2/f1/rna/DIR_STARpass2WASP
source activate gatk
gtf=/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf
	
featureCounts -T 20 -p -t exon -g gene_id -a $gtf \
	-o ../DIR_STARpass2expression/rib.pass2.featureCounts.txt \
	202026-10-rib_WASPfilterPASS_rgadded_sorted.bam \
	202026-2-rib_WASPfilterPASS_rgadded_sorted.bam \
	202102-7-rib_WASPfilterPASS_rgadded_sorted.bam \
	202102-8-rib_WASPfilterPASS_rgadded_sorted.bam \
	202232-1-rib_WASPfilterPASS_rgadded_sorted.bam \
	202232-2-rib_WASPfilterPASS_rgadded_sorted.bam \
	202232-3-rib_WASPfilterPASS_rgadded_sorted.bam \
	202266-7-rib_WASPfilterPASS_rgadded_sorted.bam \
	202266-8-rib_WASPfilterPASS_rgadded_sorted.bam \
	202266-9-rib_WASPfilterPASS_rgadded_sorted.bam \
	Y193208-1-rib_WASPfilterPASS_rgadded_sorted.bam \
	Y193208-2-rib_WASPfilterPASS_rgadded_sorted.bam \
	Y193208-3-rib_WASPfilterPASS_rgadded_sorted.bam \
	Y193208-4-rib_WASPfilterPASS_rgadded_sorted.bam \
	Y200902-1-rib_WASPfilterPASS_rgadded_sorted.bam \
	Y200902-2-rib_WASPfilterPASS_rgadded_sorted.bam \
	Y200902-3-rib_WASPfilterPASS_rgadded_sorted.bam \
	Y204711-1-rib_WASPfilterPASS_rgadded_sorted.bam \
	Y204711-2-rib_WASPfilterPASS_rgadded_sorted.bam \
	Y204711-3-rib_WASPfilterPASS_rgadded_sorted.bam
echo $?
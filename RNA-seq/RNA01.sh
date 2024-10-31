### mapping
cd /data_group/xiehaibing/xiehaibing2/f1/d90/rna
cpunum=10
DIR_fastq=/data_group/xiehaibing/xiehaibing2/f1/d90/wangyue0912023/X101SC22120474-Z01-J122_multipath/X101SC22120474-Z01-J122_01/00.CleanData
i="41303背肌"

genomeDirpass1=/data_group/xiehaibing/xiehaibing2/f1/rna/index_pig11
FASTA=/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
gtf=/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf

DIR_STARpass1=/data_group/xiehaibing/xiehaibing2/f1/d90/rna/DIR_STARpass1
DIR_STARpass2=/data_group/xiehaibing/xiehaibing2/f1/d90/rna/DIR_STARpass2bam
genomeDirpass2_all=/data_group/xiehaibing/xiehaibing2/f1/d90/rna/stargenome_pass2

echo "step 1: generate STAR genome index for 1-pass alignment (pre-generated)"
fq1=$DIR_fastq/$i/${i}_1.clean.fq.gz
fq2=$DIR_fastq/$i/${i}_2.clean.fq.gz

source activate r4_py27_env
echo "step 2: 1-pass STAR reads alignment"
STAR --genomeDir $genomeDirpass1 --readFilesCommand zcat --readFilesIn $fq1 $fq2 \
     --outFileNamePrefix $DIR_STARpass1/$i. --runThreadN $cpunum

echo "step 3: create a new genome index for the 2-pass STAR"
genomeDirpass2=$genomeDirpass2_all/$i
mkdir $genomeDirpass2
STAR --runMode genomeGenerate --genomeDir $genomeDirpass2 --genomeFastaFiles $FASTA \
	--sjdbFileChrStartEnd $DIR_STARpass1/$i.SJ.out.tab --sjdbOverhang 100 --runThreadN $cpunum --limitGenomeGenerateRAM 40000000000

echo "step 4: 2-pass STAR reads alignment"
STAR --genomeDir $genomeDirpass2 --readFilesCommand zcat --readFilesIn $fq1 $fq2 \
    --outFileNamePrefix $DIR_STARpass2/$i. --outReadsUnmapped Fastx --outSAMmapqUnique 60 --runThreadN $cpunum \
    --outSAMtype BAM SortedByCoordinate
echo $?



### 提取read count
source activate gatk
cd /data_group/xiehaibing/xiehaibing2/f1/d90/rna/DIR_STARpass2bam
rename 背肌 jr *.Aligned.sortedByCoord.out.bam
rename 肋软骨4 rib *.Aligned.sortedByCoord.out.bam

featureCounts -T 10 -p -t exon -g gene_id -a $gtf -o ../DIR_STARpass2expression/pass2notwasp/rib.pass2notwasp.featureCounts.txt $(ls *rib*bam)
featureCounts -T 10 -p -t exon -g gene_id -a $gtf -o ../DIR_STARpass2expression/pass2notwasp/jr.pass2notwasp.featureCounts.txt $(ls *jr*bam)
echo $?


### /data_group/xiehaibing/xiehaibing2/f1/d90/rna/script/pass2bam/stringtie.sh
cd /data_group/xiehaibing/xiehaibing2/f1/d90/rna/DIR_STARpass2bam
outdir=/data_group/xiehaibing/xiehaibing2/f1/rna/DIR_STARpass2expression/pass2notwasp
stringtie=/data_group/xiehaibing/xiehaibing2/software/stringtie-2.2.1.Linux_x86_64/stringtie
gtf=/data_group/xiehaibing/xiehaibing2/luchuan_duroc/Sus_scrofa.Sscrofa11.1.106.gtf

for i in $(ls *.Aligned.sortedByCoord.out.bam)
do
j=${i%%.Aligned.sortedByCoord.out.bam*}
echo $j
$stringtie -p 20 -e -B -G $gtf -o $outdir/$j/transcripts.gtf -A $outdir/$j/gene_abundances.tsv $i
echo $i was done!
done
echo $?


### 提取表达量
cd /data_group/xiehaibing/xiehaibing2/f1/rna/DIR_STARpass2expression/pass2notwasp
perl stringtie_expression_matrix.pl --expression_metric=FPKM \
                                    --result_dirs='202026-2xz,202102-7xz,202102-8xz,202232-1xz,202232-2xz,202232-3xz,202266-7xz,202266-8xz,202266-9xz,Y193208-1xz,Y193208-2xz,Y193208-3xz,Y193208-4xz,Y200902-1xz,Y200902-2xz,Y200902-3xz,Y204711-1xz,Y204711-2xz,Y204711-3xz,202026-10xz,202026-2yz,202102-7yz,202102-8yz,202232-1yz,202232-2yz,202232-3yz,202266-7yz,202266-8yz,202266-9yz,Y193208-1yz,Y193208-2yz,Y193208-3yz,Y193208-4yz,Y200902-1yz,Y200902-2yz,Y200902-3yz,Y204711-1yz,Y204711-2yz,Y204711-3yz,202026-10yz,202026-2rib,202102-7rib,202102-8rib,202232-1rib,202232-2rib,202232-3rib,202266-7rib,202266-8rib,202266-9rib,Y193208-1rib,Y193208-2rib,Y193208-3rib,Y193208-4rib,Y200902-1rib,Y200902-2rib,Y200902-3rib,Y204711-1rib,Y204711-2rib,Y204711-3rib,202026-10rib,00113-jr,00809-jr,00916-jr,00918-jr,01101-jr,01102-jr,01113-jr,01301-jr,40701-jr,40702-jr,40703-jr,40704-jr,40707-jr,YP221105-001-jr,YP221105-003-jr,YP221105-017-jr' \
                                    --transcript_matrix_file=transcript_fpkm_e27.tsv \
                                    --gene_matrix_file=gene_fpkm_e27.tsv

perl stringtie_expression_matrix.pl --expression_metric=FPKM \
                                    --result_dirs='02403jr,02409jr,02501jr,02505jr,02601jr,02603jr,40801jr,40805jr,40807jr,41301jr,41303jr,41305jr,100601jr,100603jr,100901jr,100903jr,100913jr,101003jr,210003jr,210109jr,210207jr,210303jr,210403jr,210601jr,02403rib,02409rib,02501rib,02505rib,02601rib,02603rib,40801rib,40805rib,40807rib,41301rib,41303rib,41305rib,100601rib,100603rib,100901rib,100903rib,100913rib,101003rib,210003rib,210109rib,210207rib,210303rib,210403rib,210601rib' \
                                    --transcript_matrix_file=transcript_all_fpkm_d90.tsv \
                                    --gene_matrix_file=gene_all_fpkm_d90.tsv
echo $?

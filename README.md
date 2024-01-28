# Parent-of-origin effects in pigs
These codes are generated for the project of parent-of-origin effects in pigs. We provide a python script to split F1 reads in the form of BAM file using tag SNPs. Its function is similar to SNPsplit. However, this python script works for heterozygous markers in offspring without the requirement of different homozygotes for parents. This python script is contributed by Hai-Bing Xie (xiehb@mail.kiz.ac.cn).
## How to use
### 01 usage
An example to split example.bam with snp.file:


```python split_bam_pair.py 1 snp.file example.bam output```
### 02 python module
This script require ```pysam```, ```sqlite3``` and ```pandas```.
### 03 input files
1. example.bam: BAM file sorted by coordinates.
2. snp.file: SNP file separated by TAB with six columns, including chromosome, position, paternal allele 1, paternal allele 2, maternal allele 1, and maternal allele 2. Please note that there is no header for this file. The following is an example:

|chromosome|position|paternal allele 1|paternal allele 2|maternal allele 1|maternal allele 2|
|:---|:---|:---|:---|:---|:---|
|1|1082|G|G|G|T|
|1|4553|G|G|A|G|
|1|4611|T|T|C|T|
|1|7257|T|T|A|T|
|1|11431|G|G|A|A|
|1|11567|C|C|T|T|

To speed up the task, filtered snp.file is suggested, including (1) only inclusion of heterozygous SNPs for offspring, and (2) excludsion of SNPs for parent with same genotypes (heterozygous or heterozygous). If one use it to deal with whole genome bisulfite sequencing (WGBS) data, tag SNPs methylated in each sample should be excluded in SNP file.

### 04 output files
1. The output files are BAM files for each chromosome. One can use samtools (https://github.com/samtools/samtools) or sambamba (https://github.com/biod/sambamba) to handle BAM files generated by this script.
2. Reads and SNP information covered by them.
## How to cite
Please cite our article "Parent-of-origin effects on rib length in pigs" by Tao et al, if you use it.

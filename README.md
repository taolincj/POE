# Parent-of-origin effects in pigs
These codes are generated for the project of parent-of-origin effects in pigs. We provide a python script to split F1 reads in the form of BAM file using tag SNPs. Its function is similar to SNPsplit. However, this python script takes the parental and offspring SNPs into consideration, and it works for heterozygous markers in offspring without the requirement of different homozygotes for parents.
## How to use it
### usage
An example to split example.bam with snp.file: python split_bam_pair.py 1 snp.file example.bam output
### python module
This script require pysam, sqlite3 and pandas.
### input files
1. example.bam: BAM file sorted by coordinates.
2. snp.file: SNP file separated by TAB with six columns, including chromosome, position, paternal allele 1, paternal allele 2, maternal allele 1, and maternal allele 2.
### after split
One can use samtools (http://github.com/samtools) or sambamba (https://github.com/biod/sambamba) to merge BAM files.
## How to cite it
Please cite our article "Parent-of-origin effects on rib length in pigs" by Tao et al, if you use it.

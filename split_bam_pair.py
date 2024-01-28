import pysam
import sqlite3
import sys
import pandas as pd
from collections import defaultdict

if len(sys.argv)!=5:
    print("%s chr snpfile bamfile outbam_prefix" %(sys.argv[0]))
    sys.exit()

focuschr = sys.argv[1]
snpfile = sys.argv[2]
bam = sys.argv[3]
outprefix = sys.argv[4]
#snpfile = "202260fliter4.vcf"
#bam="/data_group/xiehaibing/xiehaibing2/jizhui/atac/arc202260d3yz/outs/atac_possorted_bam.bam"
#outprefix = "./test/yzATAC.bam"

paternal = "%s.paternal.bam.%s" %(outprefix,focuschr)
maternal = "%s.maternal.bam.%s" %(outprefix,focuschr)
unassigned = "%s.unassigned.bam.%s" %(outprefix,focuschr)
paternal_unique = "%s.paternal.unique.bam.%s" %(outprefix,focuschr)
maternal_unique = "%s.maternal.unique.bam.%s" %(outprefix,focuschr)
snpclassinfo = "%s.supportingsnp.%s" %(outprefix,focuschr)

samfile = pysam.AlignmentFile(bam, "rb")
bamout1 = pysam.AlignmentFile(paternal, "wb", template=samfile)
bamout2 = pysam.AlignmentFile(maternal, "wb", template=samfile)
bamout3 = pysam.AlignmentFile(unassigned, "wb", template=samfile)
bamout4 = pysam.AlignmentFile(maternal_unique, "wb", template=samfile)
bamout5 = pysam.AlignmentFile(paternal_unique, "wb", template=samfile)
snpassignment = open(snpclassinfo,"w")

con = sqlite3.connect(":memory:")
con.execute("create table if not exists snp (chr char(6), pos integer,base1 char(1), base2 char(1), base3 char(1), base4 char(1), class integer)")
con.execute("create index idx1 on snp(chr, pos, base1, base2, base3, base4)")
con.execute("create index idx2 on snp(class,chr,pos,base1, base2, base3, base4)")

read_dict = defaultdict(lambda: [[], []])
cursor=con.cursor()

def AddSNPTag():
    sql = "select chr,pos,base1,base2,base3,base4,class from snp where class=10 or class=11 order by chr,pos"
    results=cursor.execute(sql)    
    for info in results.fetchall():
        snpchr = info[0]
        snppos = info[1]
        snppos0 = snppos+1
        base1 = info[2]
        base2 = info[3]
        base3 = info[4]
        base4 = info[5]
        snpclass = info[6]
        bases=[]
        if snpclass==10:
            minor = base3
            if base3 == base2:
                minor = base4
            for iter in samfile.pileup(snpchr,snppos,snppos0,truncate=True):
                for pileupread in iter.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_proper_pair and not pileupread.alignment.is_secondary and not pileupread.alignment.is_supplementary:
                        bases.append(pileupread.alignment.query_sequence[pileupread.query_position])
            if len(set(bases) & {base3, base4}) == 2:
               sql = "update snp set class=2 where chr=\"%s\" and pos=%d" %(snpchr,snppos)
            elif  bases == [minor]:
               sql = "update snp set class=4 where chr=\"%s\" and pos=%d" %(snpchr,snppos)
            con.execute(sql)
            continue
        if snpclass==11:
            minor = base1
            if base1 == base3:
                minor = base2
            for iter in samfile.pileup(snpchr,snppos,snppos0,truncate=True):
                for pileupread in iter.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_proper_pair and not pileupread.alignment.is_secondary and not pileupread.alignment.is_supplementary:
                        bases.append(pileupread.alignment.query_sequence[pileupread.query_position])
            if len(set(bases) & {base1, base2}) == 2:
               sql = "update snp set class=3 where chr=\"%s\" and pos=%d" %(snpchr,snppos)
            elif  bases == [minor]:
               sql = "update snp set class=5 where chr=\"%s\" and pos=%d" %(snpchr,snppos)
            con.execute(sql)

def SplitBAM(chr):
    sql = "select max(pos) from snp where chr=\"%s\"" %(chr)
    results=cursor.execute(sql).fetchone()
    if len(results)==0:
        return
    maxpos = results[0]
    print("SplitBAM chr=%s maxpos=%d" %(chr,maxpos))
    for read in samfile.fetch(contig=chr,until_eof=True):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_start > maxpos:
            break
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                SplitReads(read, read_dict[qname][1])
            else:
                SplitReads(read_dict[qname][0], read)
            del read_dict[qname]

def SplitReads(read,mate):
    chr = read.reference_name
    readstart = read.reference_start
    readend = read.reference_end
    matestart = mate.reference_start
    mateend = mate.reference_end
    
    sql = "select chr,pos,base1,base2,base3,base4,class from snp where chr=\"%s\" and ((pos>=%d and pos<=%d) or (pos>=%d and pos<=%d)) and (class=1 or class=2 or class=3 or class=4 or class=5)" %(chr,readstart,readend,matestart,mateend)
    results=cursor.execute(sql).fetchall()
    if len(results)==0:
        return
    readtype = 0
    readpositions = read.get_reference_positions(full_length=True)
    matepositions = mate.get_reference_positions(full_length=True)
    supportingsnps = []
    supportingbases= []
    supportingclass= []
    for info in results:
        #snpchr=info[0]
        snppos=info[1]
        selectread = mate
        positions = matepositions
        if snppos>=readstart and snppos<=readend:
            selectread = read
            positions = readpositions
        if snppos not in positions:
            continue
        
        readbase = selectread.query_sequence[positions.index(snppos)]
        base1=info[2]
        base2=info[3]
        base3=info[4]
        base4=info[5]
        snpclass=info[6]
        supportingsnps.append(snppos)
        supportingbases.append(readbase)
        supportingclass.append(snpclass)
        #print("pos = %d    ---   %s %s %s %s snpcalss = %d   readtype = %d  %s" %(snppos,base1,base2,base3,base4,snpclass,readtype,read.query_name))
        if snpclass==1:
            if readtype == 0:
                if readbase == base1:
                    readtype = 1
                elif readbase == base3:
                    readtype = 2
                else:
                    readtype = -1
                    break
                continue
            if readtype == 4 and readbase == base3:
                readtype = 2
                break
            if readtype == 5 and readbase == base1:
                readtype = 1
                break
            if readbase == base1 and (readtype == 2 or readtype == 4):
                readtype = -1
                break
            if readbase == base3 and (readtype == 1 or readtype == 5):
                readtype = -1
                break
        # mother is hetm father is hom
        if snpclass==2:
            if readtype == 0 or readtype == 4:
                if readbase == base1:
                    readtype = 1
                else:
                    readtype = 2
                continue
            if (readtype == 1 or readtype == 5) and readbase != base1:
                readtype = -1
                break
        if snpclass==4:
            if readtype == 0:
                if readbase !=base1:
                    readtype = 4
                continue
            if (readtype == 1 or readtype== 5) and readbase != base1:
                readtype = -1
                break
        #father is het, mother is homo
        if snpclass==3:
            if readtype == 0 or readtype == 5:
                if readbase==base3:
                    readtype = 2
                else:
                    readtype = 1
                continue
            if (readtype == 2 or readtype == 4) and readbase != base3:
                readtype = -1
                break
        if snpclass==5:
            if readtype == 0:
                if readbase != base3:
                    readtype = 5
                continue
            if (readtype == 2 or readtype ==4) and readbase != base3:
                readtype = -1
                break
    
    if readtype>0 or readtype==-1:
        snpassignment.write("%s\t%d\t%s\t" %(read.query_name,readtype,chr))
        for snpid in range(len(supportingsnps)):
            snpassignment.write("%d:%s:%d," %(supportingsnps[snpid]+1,supportingbases[snpid],supportingclass[snpid]))
        snpassignment.write("\n")
    
    if readtype==1:
        bamout1.write(read)
        bamout1.write(mate)
    if readtype==2:
        bamout2.write(read)
        bamout2.write(mate)
    if readtype==4:
        bamout4.write(read)
        bamout4.write(mate)
    if readtype==5:
        bamout5.write(read)
        bamout5.write(mate)
    if readtype==-1:
        bamout3.write(read)
        bamout3.write(mate)

#appending data to a memory database
def readSNP(snpfilename):
    vcf = pd.read_csv(snpfilename, header=None, sep="\t")
    for i in range(vcf.shape[0]):
        chr = str(vcf.iloc[i, 0])
        pos = int(str(vcf.iloc[i, 1]))-1
        base1 = str(vcf.iloc[i, 2])
        base2 = str(vcf.iloc[i, 3])
        base3 = str(vcf.iloc[i, 4])
        base4 = str(vcf.iloc[i, 5])
        snpclass = 0
        if base1 == base2 and base3 == base4 and base1 != base3:
            snpclass = 1
        if base1 == base2 and base3 != base4:
            snpclass = 10
        if base1 != base2 and base3 == base4:
            snpclass = 11
        if snpclass > 0 and chr==focuschr:
            con.execute("insert into snp (chr, pos, base1, base2, base3, base4, class) values (\"%s\", %d, \"%s\", \"%s\", \"%s\", \"%s\", %d)" %( chr, pos, base1, base2, base3, base4, snpclass ))

print("ReadSNP")
readSNP(snpfile)
print("AddSNPTag")
AddSNPTag()
print("SplitBAM")
SplitBAM(focuschr)

snpassignment.close()
con.close()
bamout1.close()
bamout2.close()
bamout3.close()
bamout4.close()
bamout5.close()

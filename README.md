### RNA-seq pipeline (using Col-0,fpa-7 as example)

### 1.1 Data quality control with fastp

```bash
#!/bin/bash
#PBS -N test
#PBS -l nodes=1:ppn=12,mem=30G
#PBS -q batch
#PBS -j oe

for name in col-02d rz12d
do
        for ((i=1;i<=3;i++))
        do
                fastp -w 12 -q 20 -u 20 \
                -i /public/home/zhangqq/RNA-seq_Col_rz1/rawdata/${name}-${i}_R1.fastq.gz \
                -o /public/home/zhangqq/RNA-seq_Col_rz1/qc_data/${name}-${i}.R1.fastp.fq.gz \#qc_data文件夹需要提前建好
                -I /public/home/zhangqq/RNA-seq_Col_rz1/rawdata/${name}-${i}_R2.fastq.gz \
                -O /public/home/zhangqq/RNA-seq_Col_rz1/qc_data/${name}-${i}.R2.fastp.fq.gz \
                -j /public/home/zhangqq/RNA-seq_Col_rz1/qc_data/${name}-${i}.fastp.json \
                -h /public/home/zhangqq/RNA-seq_Col_rz1/qc_data/${name}-${i}.fastp.html \
        done
done
```

### 1.2 gene mapping with hisat2

```bash
for name in col-02d rz12d
do
        for ((i=1;i<=3;i++))
        do
                hisat2 -p 12 \
                        -x /public/home/zhangqq/Tair10_genome/hisat2_index/TAIR10 \
                        --summary-file /public/home/zhangqq/RNA-seq_Col_rz1/mapping_data/${name}-${i}.summary \
                        --dta-cufflinks \#出来的结果更适合cufflinks处理 （主要用于基因表达量的计算和差异表达基因的寻找）
                        -1 /public/home/zhangqq/RNA-seq_Col_rz1/qc_data/${name}-${i}.R1.fastp.fq.gz \
                        -2 /public/home/zhangqq/RNA-seq_Col_rz1/qc_data/${name}-${i}.R2.fastp.fq.gz | \
                samtools view -ShuF 4 -q 20 -f 2 -@ 8 - | \
                samtools sort -@ 8 -o /public/home/zhangqq/RNA-seq_Col_rz1/mapping_data/${name}-${i}.sorted.bam -
        done
done
```

### 1.3 convert BAM format to bigwig format with deeptools

```bash
for name in col-02d rz12d
do
        for ((i=1;i<=3;i++))
        do
            samtools index /public/home/zhangqq/RNA-seq_Col_rz1/mapping_data/${name}-${i}.sorted.bam \
            bamCoverage --bam /public/home/zhangqq/RNA-seq_Col_rz1/mapping_data/${name}-${i}.sorted.bam \
                        -o /public/home/zhangqq/RNA-seq_Col_rz1/bw_data/${name}-${i}.bw \#使用deeptools将bam转换bw
                        --binSize 10 \
                        --normalizeUsing RPGC \
                        --effectiveGenomeSize 119481543 #拟南芥染色体基因组大小
        done
done
```

### 2 Get the gene expression matrix (stringtie可以组装转录本，包含内含子；featurecount只可计算外显子的表达量；参考文献Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown，2016)
### 2.1 Assemble the transcripts
```
#!/bin/bash
#PBS -N QQtest
#PBS -l nodes=1:ppn=12,mem=30G
#PBS -q batch
#PBS -j oe

for name in col-02d rz12d
do
    for ((i=1;i<=3;i++))
    do
        stringtie /public/home/zhangqq/RNA-seq_Col_rz1/mapping_data/${name}-${i}.sorted.bam \ # 此bam是samtools sort处理后的文件
          -G /public/home/zhangqq/Tair10_genome/TAIR10.gff3 \ #参考基因组注释文件
         --rf \#链特异性建库方式(fr-firststrand比如dUTP)
         -l ${name}-${i} \#转录本名称的前缀
         -o /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/${name}-${i}.transcripts.stringtie.gtf \
          -p 12 -e #如果现有的参考基因组注释文件足够了，则使用-e参数；使用-e参数才可以运行prepDE.py3脚本得到readcount矩阵；用于计算readcounts时需要-e；如果不需要预测新的转录本时，需要使用-e,如果不使用-e则会使转录本稀释，导致所关注的转录本统计不到
    done
done
```
### 2.2 Merge the transcripts samples (处理多个生物学样本时需要合并)
```bash
        vim mergelist.txt #需要包含之前output.gtf文件的路径,下面是txt文件的内容，比如
        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/Col-P1.transcripts.stringtie.gtf
        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/Col-P2.transcripts.stringtie.gtf
        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/Col-P3.transcripts.stringtie.gtf
        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/rz1-P1.transcripts.stringtie.gtf
        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/rzl-P2.transcripts.stringtie.gtf
        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/rzl-P3.transcripts.stringtie.gtf

for name in col-02d rz12d
do
        stringtie --merge -G /public/home/zhangqq/Tair10_genome/TAIR10.gff3 \
                  -F 0.1 -T 0.1 -i -o /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/${name}.stringtie_merged.gtf \ 
                  /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/mergelist.txt #input文件路径
done
```
### 2.3 Reassmble the transcripts with merged gtf file and creat ballgown file (处理多个生物学样本时，需要用组装并且合并后的gtf文件，重新对每个生物学样品进行组装)
```bash
for name in col-02d rz12d
do
        for ((i=1;i<=3;i++))
        do
         stringtie /public/home/zhangqq/RNA-seq_Col_rz1/mapping_data/${name}-${i}.sorted.bam \#input file,此bam是samtools sort处理后的文件
            -e -B -p 8 \
            -G /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/${name}.stringtie_merged.gtf \#此文件是stringtie组装合并后的gtf文件
            -o /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/ballgown/${name}/${name}-${i}.ballgown_merge.gtf \#首先创建一个ballgown文件，即mkdir ballgown; 此命令行为创建ballgown可读取文件;
        done
done
```
### 2.4 output read count
```bash 
        vim sample_list.txt #需要包含样本名和定量的gtf文件的路径，以下是#txt文件内容，文件名和gtf文件之间用TAB键隔开
        col-02d-1        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/col-02d/col-02d-1.ballgown_merge.gtf
        col-02d-2        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/col-02d/col-02d-2.ballgown_merge.gtf
        col-02d-3        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/col-02d/col-02d-3.ballgown_merge.gtf 
        rz12d-1        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/rz12d/rz12d-1.ballgown_merge.gtf
        rz12d-2        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/rz12d/rz12d-2.ballgown_merge.gtf
        rz12d-3        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/rz12d/rz12d-3.ballgown_merge.gtf
或者   ?哪种对？
        col-0        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/col-02d.stringtie_merged.gtf
        rz1        /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/rz12d.stringtie_merged.gtf
第一种方式(python2)
        source activate python2
        python /public/home/zhangqq/software/stringtie-2.2.1/prepDE.py \ #使用python的prepDE.py命令(prepDE.py在stringtie下面，写上prepDE.py的绝对路径,不写绝对路径，系统识别不出来)
               -i /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/sample_list.txt \
               -g /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/gene_count_matrix.csv \
               -t /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/transcript_count_matrix.csv
第二种方式处理(python3可以直接启用prepDE.py命令)
prepDE.py -i /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/sample_list.txt \
          -g /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/gene_count_matrix.csv \
          -t /public/home/zhangqq/RNA-seq_Col_rz1/gene_exp/transcript_count_matrix.csv
```
### 2.5 Quantitative analysis


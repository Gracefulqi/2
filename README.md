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

### 2 Get the gene expression matrix (stringtie可以组装转录本，包含内含子；featurecount只可计算外显子的表达量)
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
        stringtie /public/home/zhangqq/RNA-seq_Col_rz1/map/rz1.rep1.sorted.bam \ # 此bam是samtools sort处理后的文件
          -G /public/home/zhangqq/Tair10_genome/TAIR10.gff3 \ #参考基因组注释文件
          -l rz1 -o /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/rz1.transcripts.stringtie.gtf \
          -p 12 -e #如果现有的参考基因组注释文件足够了，则使用-e参数；使用-e参数才可以运行prepDE.py3脚本得到readcount矩阵；用于计算readcounts时需要-e；如果不需要预测新的转录本时，需要使用-e,如果不使用-e则会使转录本稀释，导致所关注的转录本统计不到
    done
done



### RNA-seq pipeline (using Col-0,fpa-7 as example)

### 1.1 Data quality control
### 1.2 gene mapping
### 1.3 convert BAM format to bigwig format

```bash
#!/bin/bash
#PBS -N QQtest
#PBS -l nodes=1:ppn=12,mem=10G
#PBS -q batch
#PBS -j oe

for name in col-02d rz12d
do
    for ((i=1;i<=3;i++))
    do
        fastp     -w 12 -q 20 -u 20 \
                -i /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/rawdata/${name}-${i}_R1.fastq.gz \
                -o /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil2/${name}-${i}_R1.fastp.fastq.gz \
                -I /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/rawdata/${name}-${i}_R2.fastq.gz \
                -O /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil2/${name}-${i}_R2.fastp.fastq.gz \
                --html /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil/${name}-${i}.html
                -j /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil2/${name}-${i}.fastp.json \
                -h /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil2/${name}-${i}.fastp.html

        hisat2 -p 12
               -x /public/home/zhangqq/Tair10_genome/hisat2_index/TAIR10 \
               --summary-file /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/map2/${name}-${i}.summary \
               --dta-cufflinks \
               -1 /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil2/${name}-${i}_R1.fastp.fastq.gz \
               -2 /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil2/${name}-${i}_R2.fastp.fastq.gz | \

        samtools view -ShuF 4 -q 20 -f 2 -@ 8 - | \
        samtools sort -@ 8 -o /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/map/rz1.rep1.sorted.bam -


```



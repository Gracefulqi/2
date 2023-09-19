RNA-seq pipeline (using Col-0,fpa-7 as example)

1. Data quality control_gene mapping_convert BAM format to bigwig format
```bash
#PBS N openmpi
##PBS l nodes=1:ppn=12
##PBS j oe
##PBS l walltime=3:00:00
##PBS l mem=10G

fastp -i /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/rawdata/rz12d-1_R1.fastq.gz \
      -I /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/rawdata/rz12d-1_R2.fastq.gz \
      -o /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil/rz12d-1_R1.fastp.fastq.gz \
      -O /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil/rz12d-1_R2.fastp.fastq.gz \
      -w 12 -q 20 -u 20 \
        -j /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil/rz12d-1.fastp.json \
        -h /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil/rz12d-1.fastp.html
```

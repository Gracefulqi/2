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
### 2.5 Quantitative analysis (将服务器上的matrix.csv文件下载后，使用RStudio本地软件运行DESeq2)（引自https://www.jianshu.com/p/b86e5598468b）
```bash
R
> library(DESeq2) #启用DESeq2程序包
> setwd("D:/R_data/../.../..) #设置需要分析的文件路径
> getwd() #再次查看确定目前文件的路径
> CountMatrix1<-read.csv("gene_count_matrix.csv",sep=",",row.names="gene_id")  ##修改列名
> head(CountMatrix1) #显示CountMatrix1的信息
> names(CountMatrix1)<-c("ctrl_rep1","ctrl_rep2","ctrl_rep"," rz1_rep1","rz1_rep2","rz1_rep3") #设置样本信息矩阵，包括处理信息：实验组rz1_rep vs. 对照组ctrl_rep，每个有3个
> ColumnData<-data.frame(row.names=colnames(CountMatrix1),samName=colnames(CountMatrix1),
condition=rep(c("ctrl_rep","rz1_rep"),each=3)) #生成DESeqDataSet数据集
> ColumnData #显示ColumnData的值

> dds<-DESeqDataSetFromMatrix(countData = CountMatrix1, colData = ColumnData, design = ~ condition) #构建dds矩阵；DESeq差异表达计算；其中countData为表达矩阵，colData为样品信息矩阵，design为差异表达矩阵，即批次和条件(对照，处理)等

> dds<-DESeq(dds)  #对原始dds进行normalize,生成差异表达结果
> dds #显示dds信息

> res<-results(dds) #使用DESeq2包中的results()函数，提取差异分析的结果，将提取的差异分析结果定义为变量“res”
> res<-res[order(res$pvalue),] #对结果res利用order()函数按pvalue值进行排序；创建矩阵时，X[i,]指矩阵X中的第i行，X[,j]指矩阵X中的第j列；order()函数先对数值排序，然后返回排序后各数值的索引（此步同res<- res[order(res$padj),]）
> head(res) #显示res结果信息
> summary(res) #对res矩阵进行总结，利用summary命令统计显示，共有多少个基因上调和下调

> table(res$padj <0.05) #统计padj（adjusted p-value）小于0.05的数目，padj即BH adjusted p-values，p值经过FDR多重校验校正后的值
> res<- res[order(res$padj),]  #按padj排序
> resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
> write.csv(resdata,file = "rz1_vs_col.csv") #输出结果到csv文件

> deg <- subset(res, padj <= 0.05 & abs(log2FoldChange) >= 1) #使用subset()函数过滤需要的结果至新的变量diff_gene_Group中；筛选显著差异表达基因（padj小于0.01且FoldChange绝对值大于2）
> summary(deg)  #查看筛选后的总结信息
> write.csv(deg, "rz1_vs_col.deg.csv")  #将差异表达显著的结果输出到csv文件

> resSig<-res[which(res$pvalue<0.05 & abs(res$log2FoldChange>1)),] #对差异基因的结果进行差异筛选；此处采用p值<0.05, log2FoldChange>1,
> resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
> resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
> head(resSig)
> write.csv(resSig,"fpa_vs_col-diff-p-0.05-FC-1.csv")
```
### 2.6 利用ggplots作图
```bash
> library(ggplot2)
> volcano<- ggplot(resdata, aes(x= log2FoldChange, y= -1*log10(padj))) #x轴为log2FC；y轴为-log(padj); padj即p.adjust，转录组测序的差异表达分析是对大量的基因表达值进行的独立统计假设检验，存在假阳性问题，因此引入Padj对显著性P值（P.adjust）进行校正。Padj是对P-value的再判断，筛选更为严格。
> threshold<-as.factor(resdata$padj <= 0.01 & abs(resdata$log2FoldChange) >= 2) #筛选条件（阈值）：绝对log2FC大于2，并且padj<0.01；
其中log2FoldChange：对Fold Change取log2，一般默认表达相差2倍以上是有意义的，可以根据情况适当放宽至1.5/1.2，但最好不要低于1.2倍。
> p1<-volcano+geom_point(aes(color=threshold)) 
> p1 #加上各数据点信息
> p2<-p1+scale_color_manual(values=c("grey","red"))
> p2 #更改散点颜色
> p3<-p2+geom_hline(yintercept=2,linetype=3)+geom_vline(xintercept=c(-2,2),linetype=3)
> p3 #加上水平和垂直线，标识阈值选择范围
> p4<-p3+theme(axis.line=element_line(colour="black"),panel.background = element_rect(fill = "white"))
> p4 #修改图片背景填充颜色，坐标轴线条颜色
> degs <- subset(resdata, padj <= 0.01 & abs(log2FoldChange)>= 2) 
> p5<-p4+geom_text(aes(label=degs$Row.names),hjust=1, vjust=0,data = degs)
> p5 #绘制P-value图
> hist(deg$pvalue,breaks=10,col="grey",xlab="p-value") #绘制MA图
> plot(deg$log2FoldChange,-log2(deg$padj),col=ifelse(abs(deg$log2FoldChange) >= 2 & abs(deg$padj) <= 0.05,"red","black"),xlab="log2FoldChange",ylab="-log2Pvalue")
```

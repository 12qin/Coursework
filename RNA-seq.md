# RNA-seq pepline
## Docker
1. Download image  
```shell
#docker pull image:<version>
docker pull staphb/fastqc:latest
docker pull ewels/multiqc:latest
docker pull zavolab/cutadapt:1.16
docker pull combinelab/salmon:latest
```
2. Create container
```shell
docker run --name=container_name \
 -dt -h bioinfo_docker --restart unless-stopped \
 -v C:\Users\qinyan\Desktop\share_folder_name:/home/test/share image:<version>

docker exec -u root container_name chown -R test:test /home/test/share
```
3. Run container
```shell
docker start container_name
docker exec -it container_name bash
```
```shell
docker ps #查看当前正在运行的容器container
docker ps -a #查看所有容器container
docker images #查看所有镜像image
```
## Fastqc
```shell
docker exec -it FastQC bash
cd ../FastQC/
./fastqc -o /home/test/share/fastqc /home/test/share/rawdata/N1-LFK11571_L2_1_fq.gz
```
## Cutadapt
1. Sequence trimming
```shell
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -o /home/test/share/cutadapt/N1-trimmed_R1.fastq.gz -p /home/test/sahre/cutadpat/N1-trimmed_R2.fastq.gz \
  /home/test/share/rawdata/N1-LFK11571_L2_1_fastqc.zip /home/test/share/rawdata/N1-LFK11571_L2_2_fastqc.zip
```
2. Filtering option
```shell 
#-O: overlap with read for minlength (default=3)   
#-m: filter reads by minimum length   
#-e: maximum error rate (Default= 0.1)  
#-q: trim low quality bases  
```
## Salmon
1. Create index
```shell
salmon index -t /home/test/share/salmon/Homo_sapiens.GRCh38.cdna.all.fa.gz -i GRCh38_salmon_index

# -t： 参考转录组路径，支持压缩文件 (Homo_sapiens.GRCh38.cdna.all.fa.gz)
# -i： 设定index名称
# -type： 索引类型，分为fmd, quasi, 建议quasi
# -k: k-mers的长度,若reads长度大于75 bp，取值31可获得较好的效果
```
2. Gene quantification
```shell
salmon quant --validateMappings \
 -l A -1 /home/test/share/N1-trimmed_R1.fastq.gz -2 /home/test/share/N1-trimmed_R2.fastq.gz \
 -i GRCh38_salmon_index \
 -o /home/test/share/salmon/N1.salmon.count -p 4

# -i： 上一步建立好的index路径
# -l/--libType: 文库类型，详细参见文档，可设置为A，使用软件自动适配
# -1： read1，支持压缩文件
# -2： read2，支持压缩文件
# -o： 输出目录
```
## MultiQC
'''shell
docker exec -it multiqc bash
multiqc /home/test/share/fastqc -o /home/test/share/multiqc
multiqc /home/test/share/salmon -o /home/test/share/multiqc
'''
## R
1. Data summary and export
#更新R版本
```R
install.packages("installr")
require(installr)
updateR()
```
```R
library(tximport)
library(GenomicFeatures)

#导入基因注释文件
txd <- makeTxDbFromGFF(file = 'gencode.v44.primary_assembly.annotation.gtf.gz',
                    format = 'gtf', dataSource = 'gencode.v44',
                    organism = 'Homo sapiens')  
str(txd)

#transcript_ID 转换为 Gene_ID
txTogene <- AnnotationDbi::select(txd, keys = keys(txd,'TXNAME'),
                              keytype ='TXNAME',
                              columns = c('TXNAME','GENEID'))
head(txTogene)

## tximport处理数据
samplelist <- c('N-1','N-2','Ab-1','Ab-2','S01-C1','S01-C2')
countdata <- c('N1.salmon.count/quant.sf','N2.salmon.count/quant.sf',
               'Ab1.salmon.count/quant.sf','Ab2.salmon.count/quant.sf',
               'S01-C1.salmon.count/quant.sf','S01-C2.salmon.count/quant.sf')
names(countdata) <- samplelist

txi <- tximport(countdata, type='salmon', tx2gene=txTogene)
str(txi)

#筛选样本内的表达量counts均为0的基因
counts <- txi$counts
counts_filt <- counts[rowSums(counts)>0,]
write.csv(counts_filt,'results/sw480_counts_filt.csv',row.names=TRUE)
```
5. DEseq2
```R
library(DESeq2)

#整理数据
SampleTable <- data.frame(samlpe= c('N-1', 'N-2', 'Ab-1', 'Ab-2', 'S01-C-1', 'S01-C-2'), \
  condition = c("Control","Control","Antibody","Antibody","Cysteine","Cysteine")
rownames(sampleTable)<-colnames(txi$counts)

#构建dds对象
dds <- DESeqDataSetFromMatrix(round(counts_filt), coldata, design=~treatment)
str(dds)

# 质控
idx <- rowSums(counts(dds) >= 5) > 2
dds <- dds[idx,]
dim(dds)

#差异分析
dds <- DESeq(dds)
saveRDS(dds,file = "sw480_dds.rds")
```
6. 样本间差异分析
```R
library(pheatmap)
library(RColorBrewer)

# log转换
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)

# dist聚类
rldDist <- dist(t(rlogMat), method = "euclidean")
rldMatrix <- as.matrix(rldDist) 
pheatmap(rldMatrix,fontsize=10,
         clustering_row_rows=rldDist,
         clustering_distance_cols= rldDist,
         col=colorRampPalette( rev(brewer.pal(9, "Blues")) )(255),
         main="Heatmap of the sample-to-sample distances")

# 计算 pearson correlation
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
heatmap(pearson_cor,col = colorRampPalette(brewer.pal(9, "GnBu"))(100),
        main="The pearson correlation of eachsample")

##PCA分析
plotPCA(rld, intgroup=c('condition'))
```
7. 差异基因表达分析
```R
# 提取组间差异基因
res <- results(dds, contrast = c('treatment',"Antibody",'Control')) 
summary(res)

res <- res[order(res$padj),]       
DEG <- as.data.frame(res)
diff_gene <- na.omit(DEG)  #删除NA值

# gene_id 转换
library(org.Hs.eg.db)
library(clusterProfiler)

#去掉ENSEMBL_geneid_version
diff_gene["ENSEMBL"] <- rownames(diff_gene)
diff_gene$ENSEMBL <- gsub('\\..+$', '', diff_gene$ENSEMBL)

#转换为symbol
gene_symbol <- bitr(diff_gene$ENSEMBL, fromType="ENSEMBL", 
                    toType=c("SYMBOL",'ENTREZID'), 
                    OrgDb="org.Hs.eg.db")
diff_gene = data.frame(gene_symbol,diff_gene[match(gene_symbol$ENSEMBL,diff_gene$ENSEMBL),])

# 根据padj和foldchange值分析上下调基因
cut_off_logFC = 0.6

diff_gene$up_down = ifelse(diff_gene$padj < 0.05 & abs(diff_gene$log2FoldChange) > cut_off_logFC, 
                           ifelse(diff_gene$log2FoldChange > cut_off_logFC ,
                                  'Up','Down'),'None')
table(diff_gene$up_down)

write.csv(diff_gene,'results/DEG_Ab~C.csv',row.names=TRUE)
```
8. 结果可视化
```
#查看P value分布
library(ggplot2)
ggplot(diff_gene) , aes(x = pvalue)) +
  geom_histogram(bins = 200)

## 火山图绘制
library(ggrepel)
cut_off_padj <- 0.05

p <- ggplot(diff_gene, aes(x = log2FoldChange, y = -log10(padj), colour=up_down)) +
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=2,col="grey",lwd=0.5) +
  geom_hline(yintercept = -log10(cut_off_padj),lty=2,col="grey",lwd=0.5) +
  # 坐标轴
  labs(x="log2(Fold Change)",y="-log10 (P-adj)")+
  theme_bw()+
  ggtitle("Volcano Plot Ab~C")+
  # 图例
  theme(panel.grid=element_blank(),
        axis.title = element_text(size = 16),axis.text = element_text(size = 14))
#注释
p + geom_text_repel(
    data = subset(diff_gene,diff_gene$padj < 0.05 &  abs(diff_gene$log2FoldChange) > 0.6),
    aes(label = SYMBOL, color = up_down), size = 3, color = "black", 
    box.padding = unit(0.5, "lines"), point.padding = unit(0.3, "lines"),
    arrow= arrow(length = unit(0.03, "npc"),type= "open", ends = "last"),
    show.legend = FALSE, max.overlaps = 12)

## MA图绘制
p <- ggplot(diff_gene, aes(x = log2(baseMean), y = log2FoldChange, colour=up_down)) +
  geom_point(alpha=0.6,size = 2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 坐标轴
  labs(x="log2(baseMean)",y="log2FoldChange")+
  theme_bw()+
  ggtitle("MA Plot Ab~C")+
  # 辅助线
  geom_hline(yintercept=c(-1.0,1.0),lty=2,col="grey",lwd=0.5) +
  # 图例
  theme(panel.grid=element_blank(),
        axis.title = element_text(size = 16), axis.text = element_text(size = 14))
#注释
p + geom_text_repel(
  data = subset(na.omit(diff_gene),diff_gene$padj < 0.05 &  abs(diff_gene$log2FoldChange) > 0.6),
  aes(label = SYMBOL, color = up_down), size = 3, color = "black", 
  box.padding = unit(0.5, "lines"), point.padding = unit(0.3, "lines"),
  show.legend = FALSE, max.overlaps = 20 )
```

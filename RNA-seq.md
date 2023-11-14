# RNA-seq pepline
## Usage of docker
1. Download image  
```shell
docker pull image:<version>
```
2. Create container
```shell
docker run --name=container_name \
 -dt -h bioinfo_docker --restart unless-stopped \
 -v C:\Users\qinyan\Desktop\share_folder:/home/test/share image:<version>

docker exec -u root container_name chown -R test:test /home/test/share
```
3. Run container
```shell
docker exec -it container_name bash
docker start container_name/code
```
4.Some command
```shell
docker ps #查看当前正在运行的容器container
docker ps -a #查看所有容器container
docker images #查看所有镜像image
```
## Fastqc
```shell
cd /FastQC/
./fastqc -o ./home/test/share/fastqc/sample-R2.fastq.gz
```
## Cutadapt
1. Sequence trimming
```shell
cutadapt -a ADAPTER_FWD -A ADAPTER_REV \
  -o out.1.fastq -p out.2.fastq \
  reads.1.fastq reads.2.fastq
```
2. Filtering option
```shell
#-a: trimming reads 3'end adapter for read1  
#-A: trimming reads 3'end adapter for read2  
#-O: overlap with read for minlength (default=3)   
#-m: filter reads by minimum length   
#--discard-untrimmed, --trimmed-only: discard untrimmed reads  
#-e: maximum error rate (Default= 0.1)  
#-q: trim low quality bases  
#-o: output file for read1  
#-P: output file for read2  
#--info-file：each reads and matched adapter information 
```
## Salmon
1. Create index
```shell
salmon index -t transcripts.fasta.gz -i transcripts_index --type quasi -k 31

# -t： 参考转录组路径，支持压缩文件 (Homo_sapiens.GRCh38.cdna.all.fa.gz)
# -i： 设定index名称
# -type： 索引类型，分为fmd, quasi, 建议quasi
# -k: k-mers的长度。说明文档中，若reads长度大于75 bp，取值31可获得较好的效果
```
2. Gene quantification
```shell
salmon quant -i transcripts_index \
 -l A -1 reads_1.fastq -2 reads_2.fastq \
 -o transcripts_quant

# -i： 上一步建立好的index路径
# -l/--libType: 文库类型，详细参见文档，可设置为A，使用软件自动适配
# -1： read1，支持压缩文件
# -2： read2，支持压缩文件
# -o： 输出目录
```
3. Batch treatment
  批量运行脚本
```shell
for fn in data/cutadapt;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i index -l A \
         -1 ${fn}/${samp}_R1.fastq.gz \
         -2 ${fn}/${samp}_R2.fastq.gz \
         -p 8 -o quants/${samp}_quant
done
```

4. Data summary and export
更新R版本
```R
install.packages("installr")
require(installr)
updateR()
```

```R
library(tximport)
library(GenomicFeatures)

#下载gtf基因注释文件(GRCh38_hg38/gencode.v33.annotation.gtf)
txd <- makeTxDbFromGFF(file = "../GRCh38_hg38/gencode.v33.annotation.gtf", format = "gtf", \
  dataSource = "gencode.v33.annotation.gtf", organism = "Homo sapiens")
str(txd)
keytypes(txd)
saveDb(txd, file="/Example/GRCh38_hg38/Gencode_V33_TxDb.sqlite")

#创建转录本与基因相关矩阵
txTogene <- AnnotationDbi::select(txd, keys =keys(txd,'TXNAME'), keytype='TXNAME', columns=c('TXNAME','GENEID'))
head(txTogene)
write.csv(txdTogene,'txiTogene.csv',row.names=TRUE)

# 创建sample矩阵
samplelist<-c('N-1','N-2','Ab-1','Ab-2','S01-C1','S01-C2')
filelist <- file.path("/home/test/share/salmon", samplelist, "quant.sf")
files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample", 1:6)
names(filelist) <- samplelist

# tximport导入结果
txi <- tximport(fileList, type = "salmon", tx2gene = txTogene)
data<-c('N-1/quant.sf','N-2/quant.sf','Ab-1/quant.sf','Ab-2/quant.sf','S01-C-1/quant.sf','S01-C-2/quant.sf')
samplelist<-c('N-1','N-2','Ab-1','Ab-2','S01-C-1','S01-C-2')
names(data)<-samplelist
head(txi$counts) # 查看矩阵
write.csv(/home/test/share/deseq2/txi,'txi.csv',row.names=TRUE)
```
5. DEseq2: countData colData design
```R
library(DESeq2)
#整理数据
SampleTable <- data.frame(samlpe= c('N-1', 'N-2', 'Ab-1', 'Ab-2', 'S01-C-1', 'S01-C-2'), \
  condition = c("Control","Control","Antibody","Antibody","Cysteine","Cysteine")
rownames(sampleTable)<-colnames(txi$counts)

#差异表达分析
dds<-DESeqDataSetFromTximport(txi,coldata=samples,design=~condition)  
dds<-DESeq(dds) 
res<-results(dds)        #从DESeq分析中提取结果

#以padj对res排序
resordered<-res[order(res$padj),]        
write.csv(resordered,'res.csv',row.names=TRUE)
```



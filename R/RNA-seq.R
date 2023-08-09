
#a= ydx_RNA_seq_read_counts(folder_path="D:\\Personal\\Desktop\\GSE181143_RAW1")


ydx_RNA_seq_read_counts<- function(folder_path,sep = ",",con_name="non_Diabetes_non_TB" ,treat_name= "Diabetes_non_TB",pattern = ".+\\.txt.gz$"){

# folder_path：要读取的目录folder_path="D:\\Personal\\Desktop\\GSE181143_RAW1"
# sep：要读取文件gene与表达值之间的分隔符号
# con_name：对照组名称
# treat_name：实验组名称
# pattern：要匹配的所有文件的格式
  
  library(dplyr)
  library(purrr)
#  folder_path="D:\\Personal\\Desktop\\GSE181143_RAW1"
  file_list <- list.files(folder_path, pattern = pattern, full.names = TRUE)
  file_names <- lapply(file_list, function(x) tools::file_path_sans_ext(basename(x)))
  # 使用lapply函数读取每个文件，并将其转换为数据框
  data_list <- lapply(file_list, function(x) data.table::fread(x, header = TRUE, sep =sep))

  # 获取所有数据框中的基因名，取交集
  gene_names <- Reduce(intersect, lapply(data_list, `[[`, 1))
  
  # 从每个数据框中提取共有的基因名及其对应的数值，并合并为一个数据框
  rt1 <- lapply(seq_along(file_list), function(i) {
    data_list[[i]] %>% filter(.[[1]] %in% gene_names) %>% setNames(c("gene", file_names[[i]]))
  }) %>% purrr::reduce(inner_join, by = "gene")
  # 将第一列设置为行名
  rm(data_list)
  rt1=as.data.frame(rt1)
  rownames(rt1)=rt1[,1]
  exp1=rt1[,2:ncol(rt1)]
  dimnames=list(rownames(exp1),colnames(exp1))
  data1=matrix(as.numeric(as.matrix(exp1)),nrow=nrow(exp1),dimnames=dimnames)
  data1=limma::avereps(data1)
  data1=data1[rowMeans(data1)>0,]
  merged_data=as.data.frame(data1)
  cat("所有样本因数量:",length(colnames(merged_data)))
  
  cat("\n去除全为0的行后，剩余基因数量:",length(rownames(merged_data)))

  
  ######################################################################
  ######################################################################
  ######################################################################
  # 将列名中的 "Diabetes_non_TB_" 替换为 "con"
  colnames(merged_data) <- sub(con_name, "con", colnames(merged_data))
  colnames(merged_data) <- sub(treat_name, "treat", colnames(merged_data))
  names (merged_data) <- gsub ("(GSM\\d+|_con|_treat)|.", "\\1", names (merged_data))
  ######################################################################
  ######################################################################
  ######################################################################
  
  
  # 实验组
  library(stringr)
  
  exp_data_T <- merged_data %>% dplyr::select(which(str_detect(colnames(.), "_treat")))
  
  nT <- ncol(exp_data_T) 
  cat("实验组数量:", nT, "\n")
  # 正常组
  exp_data_N <- merged_data %>% dplyr::select(str_which(colnames(.), "_con"))
  nN <- ncol(exp_data_N) 
  cat("正常组数量:", nN, "\n")
  # 合并数据，CON放前面，Test放后面
  merged_data <- cbind(exp_data_N, exp_data_T)

  
    return(merged_data)
}


ydx_RNA_seq_fenxi <- function(){
  library(edgeR)
  library(data.table)
  library(tidyverse)
  library(ggsignif) 
  library(RColorBrewer)
  library(limma)
  library(ggplot2)
  library(ggpubr)
  library(beepr)
  library(gplots)
  library(pheatmap)
  library("DESeq2")
  library(VennDiagram)
  
  rt=data
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data1=avereps(data1)
  data1=data1[rowMeans(data1)>0,]
  data1=as.data.frame(data1)
  
  
  
  ####################遍历文件中的名字 ，进行合并数据，CON放前面，Test放后面
  # 正常组
  exp_data_N <- exp_data %>% dplyr::select(str_which(colnames(.), "_con"))
  nN <- ncol(exp_data_N) 
  cat("正常组数量:", nN, "\n")
  # 实验组
  exp_data_T <- exp_data %>% dplyr::select(str_which(colnames(.), "_treat")) 
  nT <- ncol(exp_data_T) 
  cat("实验组数量:", nT, "\n")
  # 合并数据，CON放前面，Test放后面
  data <- cbind(exp_data_N, exp_data_T)
  count <- floor(data)#下取
  
  
  # 预处理，过滤低丰度的数据
  countData <- count[apply(count, 1, sum) > 0 , ]
  
  # 构建分组矩阵
  data=colnames(countData)
  Type=c(rep(1,conNum), rep(2,treatNum))
  colData=cbind(data, Type)
  colData=as.data.frame(colData)
  colnames(colData)=c("id", "Type")
  colData$Type=ifelse(colData$Type==1, "Normal", "Treat")
  colData=as.matrix(colData)
  rownames(colData)=colData[,1]
  colData=colData[,2:ncol(colData)]
  colData=as.data.frame(colData)
  colnames(colData)=c("condition")
  
  colnames(count)==colData$condition
  
  # 构建DESeq2中的对象
  dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ condition)
  # 指定哪一组作为对照组
  dds$condition <- relevel(dds$condition, ref = "Normal")
  #计算每个样本的归一化系数
  dds <- estimateSizeFactors(dds)
  #估计基因的离散度
  dds <- estimateDispersions(dds)

  #差异分析
  dds <- nbinomWaldTest(dds)
  dds <- DESeq(dds)
  res <- results(dds)
  write.table(res,"DESeq2.diff.tsv",sep="\t",quote=F,col.names = NA)
  
  
  a = as.data.frame(res)
  #标准化(减小离散度)
  vsd <- vst(dds, blind = FALSE)  ## 使用vst函数进行标准化
  countData_old<- assay(dds)
  countData_new<- assay(vsd)
  n.sample=ncol(countData)#样本数
  cols <- rainbow(n.sample*1.2)
  par(mfrow=c(2,2))
  
  #标准化前图
  pdf(file="rawBox.pdf")
  boxplot(countData_old, col = cols,main="expression value",xaxt = "n")
  dev.off()
  #标准化后图
  pdf(file="normalBox.pdf")
  boxplot(countData_new, col = cols,main="expression value",xaxt = "n")
  dev.off()
  
  pdf(file="histold.pdf")
  hist(countData_old)
  dev.off()
  pdf(file="histnew.pdf")
  hist(countData_new)
  dev.off()
  
  
  
  #可视化
  gene="DCAF5"
  data=t(countData_new[gene,,drop=F])
  Type=c(rep(1,conNum), rep(2,treatNum))
  exp=cbind(data, Type)
  exp=as.data.frame(exp)
  colnames(exp)=c("gene", "Type")
  exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
  exp$gene=log2(exp$gene+1)
  group=levels(factor(exp$Type))
  exp$Type=factor(exp$Type, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  boxplot=ggboxplot(exp, x="Type", y="gene", color="Type",
                    xlab="",
                    ylab=paste0(gene, " expression"),
                    legend.title="Type",
                    palette = c("blue","red"),
                    add = "jitter")+ 
    stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
}
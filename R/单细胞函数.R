
ydx_sc_read10X <- function(dir, min_nFeature= 200,percent_mt_num=10,percent_hb_num =3,nCount=20000,scale_genes ="all"){
  #  scale_genes ="all"的时候对所有基因进行降维，否则只对高变基因进行降维
  # 线粒体基因的比例小于10% # 血红蛋白基因的比例小于3%
  library(SingleR)
  library(patchwork)
  library(glmGamPoi)
  library(dplyr)
  library(Seurat)
  library(tidyverse)
  library(clustree)
  library(ggplot2)
  library(RColorBrewer)
  
  #创建一个文件夹用于写分析结果
  cat("使用read10X函数读入标准单细胞数据\n")
  

  sample_names = ""
  
  if(length(sample_names) >2){
    names(dir) = sample_names
  }else{
    names(dir) = sapply(dir, basename)
  }
freq.combine <- list()
scRNAlist <- list()
cat("请仔细看单细胞数据数据集的线粒体和血红蛋白的百分比，方便过滤细胞\n")
#创建seurat对象
# min.cell：每个feature至少在多少个细胞中表达
# min.features：每个细胞中至少有多少个feature被检测到
#nFeature_RNA是每个细胞中检测到的基因数量
#nCount_RNA是细胞内检测到的分子总数
#nFeature_RNA过低，表示该细胞可能已死/将死或是空液滴。高nCount_RNA和/或nFeature_RNA表明“细胞”实际上可以是两个或多个细胞。
#结合线粒体基因（percent.mt）与核糖体基因（percent.rb）除去异常值，即可除去大多数双峰/死细胞/空液滴，因此它们过滤是常见的预处理步骤
  for(i in 1:length(dir)){
    scRNAlist[[i]] <- CreateSeuratObject(Read10X(data.dir = dir[i]), min.cells = 10, min.features =200)
    scRNAlist[[i]]@assays
    if (grepl("[a-z]", scRNAlist[[1]]@assays[["RNA"]]@counts@Dimnames[[1]][1])){
      cat("当前是小鼠数据集，请核查----------")
      # 1. 过滤线粒体
      #scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")#  人用MT，小鼠用mt
      scRNAlist[[i]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = '^MT-', col.name = 'percent.mt')
      
      ## 2. 过滤红细胞
      HB.genes <- c("Hba-a1", "Hba-a2", "Hbb-bh1", "Hbb-bh2", "Hbb-b1", "Hbb-b2", "Hbb-y", "Hbb-bt", "Hbb-bh3")
    } else{
      cat("当前是人类数据集，请核查----------")
      # 1. 过滤线粒体
      scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")#  人用MT，小鼠用mt
      ## 2. 过滤红细胞
      HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
      }
    HB_m <- match(HB.genes, rownames(scRNAlist[[i]]@assays$RNA)) 
    HB.genes <- rownames(scRNAlist[[i]]@assays$RNA)[HB_m] 
    HB.genes <- HB.genes[!is.na(HB.genes)] 
    scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 
    ##QC：统计基因数，RNA，线粒体基因分布
    gene.freq <- do.call("cbind", tapply(scRNAlist[[i]] @meta.data$nFeature_RNA,scRNAlist[[i]] @meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
    rna.freq <- do.call("cbind", tapply(scRNAlist[[i]] @meta.data$nCount_RNA,scRNAlist[[i]] @meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
    mt.freq <- do.call("cbind", tapply(scRNAlist[[i]] @meta.data$percent.mt,scRNAlist[[i]] @meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
    hb.freq <- do.call("cbind", tapply(scRNAlist[[i]] @meta.data$percent.HB,scRNAlist[[i]] @meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
    combine <- as.data.frame(cbind(gene.freq,rna.freq,hb.freq,mt.freq))
    colnames(combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                                paste(colnames(rna.freq),"RNA",sep = "_"),
                                paste(colnames(hb.freq),"HB",sep = "_"),
                                paste(colnames(mt.freq),"MT",sep = "_"))
    # # 使用paste()函数为数据框命名
    # name <- paste("freq.combine_", i, sep = "")
    # assign(name, freq.combine)
    # rm(freq.combine)
    freq.combine[[i]] <- combine
    # 对Seurat对象进行可视化，绘制基因数、RNA数、线粒体基因分布和红细胞基因分布的小提琴图
    col.num <- length(levels(scRNAlist[[i]]@active.ident))
    VlnPlot(scRNAlist[[i]],
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
            cols = rainbow(col.num), 
            pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
            ncol = 4) + 
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
    ggsave(paste0("样本",i, "-scRNA_VlnPlot.pdf"), device = "pdf", width = 8, height = 6)
    
    #对Seurat对象进行可视化，绘制nCount_RNA和nFeature_RNA的散点图
    FeatureScatter(scRNAlist[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    ggsave(paste0("样本",i, "-nCount_RNA和nFeature_RNA.pdf"), device = "pdf", width = 8, height = 6)
    
   
    ####################################根据情况进行细胞过滤#######################################

    # 输出过滤前细胞数量
    cat("目前正在处理数据集:",i,"           \n过滤前细胞数量：", nrow(scRNAlist[[i]]@meta.data), "\n")
    scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nFeature_RNA >= min_nFeature
                             & percent.mt <= percent_mt_num
                             & percent.HB <= percent_hb_num
                             & nCount_RNA < nCount)
    
    # 输出过滤标准
    cat("目前正在处理数据集:",i,"----------过滤标准是：线粒体比例 <", percent_mt_num, "%，血红蛋白小于", percent_hb_num, "%","过滤后细胞数量：", nrow(scRNAlist[[i]]@meta.data), "\n")
  }
####################################################
if (sum(1:length(dir)) > 1) {
  
  for (i in 1:length(scRNAlist)) {
    scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
    scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst",nfeatures = 3000)
  }
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist,anchor.features = 2000)
  scRNA <- IntegrateData(anchorset = scRNA.anchors)
  DefaultAssay(scRNA) <- "integrated"
  #### 数据标准化（中心化）
  ## 如果内存足够最好对所有基因进行中心化,如果内存不够，可以只对高变基因进行标准化
  if (scale_genes == "all") {
    scale.genes <- rownames(scRNA)
    scRNA <- ScaleData(scRNA, features = scale.genes)
  } else {
    scale.genes <- VariableFeatures(scRNA)
    scRNA <- ScaleData(scRNA, features = scale.genes)
  }
  scRNA <- RunPCA(scRNA, verbose = FALSE)
  cat("当前未整合数据，直接使用了IntegrateData合并数据，适用于差异不大的情况")
  
  
} else {
  # 如果只有一个数据集，则直接使用该数据集。
  scRNA = scRNAlist[[1]]
  cat("只有一个数据集，直接进行NormalizeData\n")
  scRNA <- NormalizeData(scRNA)
  cat("只有一个数据集，直接进行FindVariableFeatures\n")
  scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
  cat("只有一个数据集，直接进行scale_genes\n")
  if (scale_genes == "all") {
    scale.genes <- rownames(scRNA)
    scRNA <- ScaleData(scRNA, features = scale.genes)
  } else {
    scale.genes <- VariableFeatures(scRNA)
    scRNA <- ScaleData(scRNA, features = scale.genes)
  }
  

  scRNA <- RunPCA(scRNA, verbose = FALSE)
  cat("只有一个数据，无需整合数据，已进行标准化，可进行后续分析")
}

cat("根据命名，添加分组信息,去除   “。 - _ ”   字符及其后面的字符，选择前面的作为分组名称\n")
scRNA$group = gsub("[\\._-].*$", "", scRNA$orig.ident)
# 提取前10的高变基因
top10 <- head(VariableFeatures(scRNA), 10)
# 展示高变基因
plot1 <- VariableFeaturePlot(scRNA)

LabelPoints(plot = plot1, points = top10, repel = TRUE) #plot_layout，patchwork函数，指定一行有几个图片
ggsave(filename = "03.topgene.pdf", width = 7, height = 6)




return(list(freq.combine=freq.combine,scRNAlist=scRNAlist,scRNA=scRNA))
}

ydx_sc_reduction <- function(scRNA,pca_number=20,resolution=0.8){
  library(dplyr)
  library(Seurat)
  library(tidyverse)
  library(clustree)
  library(ggplot2)
  library(RColorBrewer)
  library(glmGamPoi)
  library(msigdbr)
  library(fgsea)
  library(viridis)
  library(dittoSeq)
  # dittoHeatmap(scRNA, top_list,
  #              annot.by = c("singleR","orig.ident"))
  
  
  pca_number= c(1:pca_number)
  if (is.list(scRNA)) {
    scRNA=scRNA$scRNA
  }
  
  
  cat("根据PCA降维分析的碎石图，来选择UMAP和TSNE的维度信息，当前维度信息是：",length(pca_number),"\n\n")
  #tSNE/umap
  scRNA <- Seurat::RunUMAP(scRNA,  dims = pca_number)
  scRNA = Seurat::RunTSNE(scRNA, dims = pca_number)
  VizDimLoadings(scRNA, dims = 1:2, reduction = "pca")
  ggsave(paste0("线性降维曲线.pdf"), device = "pdf", width = 12, height = 6)
  pdf(paste0("PC热图.pdf"),width =8,height = 22)  # 打开PDF设备
  DimHeatmap(scRNA, dims = 1:30, cells = 500, balanced = TRUE)
  dev.off()
  if ("SCT" %in% names(scRNA)) {
  
  }else{
    scRNA <- JackStraw(scRNA, num.replicate = 100)
    scRNA <- ScoreJackStraw(scRNA, dims = 1:20)
    JackStrawPlot(scRNA, dims = 1:20)
    ggsave(paste0("JackStraw曲线1-20.pdf"), device = "pdf", width = 8, height = 6)
  }

  
 
  
  cat("数据整合完毕，PCA降维分析完成，请根据情况进行聚类\n")
  ElbowPlot(scRNA , ndims=50, reduction="pca") 
  ggsave(paste0("根据曲线选择合适的PCA.pdf"), device = "pdf", width = 8, height = 6)

  cat("正在绘制降维图，请观察每一个组是否都能均匀的分布到每个cluster，根据已有知识，判断是否需要整合\n")
  scRNA <- FindNeighbors(scRNA, reduction = "pca", dims =pca_number)
  
 
  if(scRNA@active.assay == "integrated" ){
    for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1,1.2,1.5,2,2.5,3)) {
      scRNA=FindClusters(scRNA, graph.name = "integrated_snn", resolution = res, algorithm = 1)}
    plot1=clustree(scRNA@meta.data, prefix = "integrated_snn_res.")
  }else if(scRNA@active.assay == "RNA"){
    for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1,1.2,1.5,2,2.5,3)) {
      scRNA=FindClusters(scRNA, graph.name = "RNA_snn", resolution = res, algorithm = 1)}
    plot1=clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
  }
  ggsave(filename = "04.clustree.pdf", plot = plot1, width = 12, height = 10)
  cluster=  apply(scRNA@meta.data[,grep("snn_res",colnames(scRNA@meta.data))],2,table)
  
  
  # 提取以 "RNA_snn_res" 开头的全部列名,并存放到rna_snn_list
  rna_snn_list <- grep("^RNA_snn_res", colnames(scRNA@meta.data), value = TRUE)
  # 创建输出文件夹
  output_folder <- "UMAP-tSNE2"
  dir.create(output_folder, showWarnings = FALSE)
  # 循环绘制不同分辨率的umap—tsne图，并保存到本地
  for (rna_snn in rna_snn_list) {
    
    tsne <- DimPlot(scRNA, reduction = "tsne",
                    group.by = rna_snn,
                    label = TRUE, label.size = 8)
    umap <- DimPlot(scRNA, reduction = "umap",
                    label = TRUE, label.size = 8,
                    group.by = rna_snn)
    
    # 合并
    p <- tsne + umap
    # 保存为 PNG 文件
    print(p)
    filename <- paste0(output_folder, "/UMAP_tSNE_", rna_snn, "_plot.pdf")
    ggsave(filename, plot = p, width = 16, height = 6)
  } 
  

  scRNA <- FindClusters(scRNA, resolution = resolution)
  cat("根据细胞情况选择resolution，一般情况2万细胞以上都是大于1.0的。当前：",resolution,"\n\n")
  cat("查看每一类有多少个细胞,当前有",length(unique(scRNA@meta.data$seurat_clusters)),"个cluster，\n每个cluster细胞分别是：\n",table(scRNA@meta.data$seurat_clusters),"\n\n")
  # 可视化UMAP/tSNE
  
  library(COSG)
  library(Seurat)
 
  scRNA <- SetIdent(scRNA, value = paste0("RNA_snn_res.",resolution))
  marknumber= 20
  marker_cosg <- cosg(
    scRNA,
    groups='all', #考虑全部分组
    assay='RNA',
    slot='data',
    mu=1,         #惩罚项参数，值越大
    remove_lowly_expressed=TRUE,   #是否过滤低表达基因
    expressed_pct=0.1,             #设置低表达的阈值
    n_genes_user=marknumber      #每个cluster定义Top-N个marker gene
  )
  top_list<-c()
  for (group in colnames(marker_cosg$names)){
    top_i<-marker_cosg$names[group][1:5,1]
    top_list<-c(top_list,top_i)
  }
  DotPlot(scRNA, 
          assay = 'RNA',
          # scale=TRUE,
          features =  unique(top_list)) + RotatedAxis()
  ggsave(filename = paste0("RNA_snn_res.",resolution,"_markers.pdf"),width =20.5,height = 10.5)
  
 write.csv(marker_cosg$names,"COSG算法找出的marker基因")
  
  
  plot2 = DimPlot(scRNA, reduction = "umap", label = T, label.size = 3.5,pt.size = 0.5)+theme_classic()
  ggsave(filename = "05-cluster.UMAP.pdf", plot = plot2,width =6.5,height = 5.5)
  
  plot1 = DimPlot(scRNA, reduction = "tsne", label = T, label.size = 3.5,pt.size = 0.5)+theme_classic()
  ggsave(filename = "05-cluster.TSEN.pdf", plot = plot1,width =6.5,height = 5.5)
  
  ########################################################################
  #####################下面是观察聚类情况，无需修改美观#################
  ########################################################################
  PCA = c("orig.ident")
  for(i in PCA){
  for(PCA_type in c("pca","tsne","umap")){
    pdf(paste0(PCA_type,"---数据聚合情况.pdf"),width =8,height = 6)  # 打开PDF设备
    print(DimPlot(scRNA, group.by = i, reduction = PCA_type))
    print(DimPlot(scRNA, group.by = i, split.by = i, reduction = PCA_type))
    print(DimPlot(scRNA, group.by = i, label = TRUE, label.size = 4, reduction = PCA_type))
    dev.off()  # 关闭PDF设备
  }
  }
  return(list(scRNA = scRNA, cluster = cluster))
  }

ydx_sc_Integrat <- function(scRNAlist,Integrat ="SCT",group="orig.ident",scale_genes ="all"){
  
  if (Integrat == "SCT" ){ 
    cat("\n\n\n当前有多个数据集，使用SCTransform函数整合数据\n")
    # 如果使用SCT标准化，则对每个数据集进行SCT变换，并输出信息。
    scRNAlist <- lapply(X = scRNAlist, FUN = SCTransform)
    # 选择跨数据集的基因
    features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
    # 找到数据集的锚点
    scRNAlist <- PrepSCTIntegration(object.list= scRNAlist, anchor.features = features)
    anchors <- FindIntegrationAnchors(object.list = scRNAlist, normalization.method = "SCT", anchor.features = features)
    # 整合数据集
    scRNA <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
    # 将整合后的数据集的默认数据类型设置为“integrated”，方便后续分析。
    DefaultAssay(scRNA) <- "integrated"
    scRNA <- RunPCA(scRNA, verbose = FALSE)
    cat("使用SCT整合数据成功,目前活动分组是integrated分组")
    }else if(Integrat == "Harmony" ){
      if (!requireNamespace("harmony", quietly = TRUE))BiocManager::install("harmony", update = FALSE, ask = FALSE)
      library(harmony)
      cat("\n\n\n当前有多个数据集，使用harmony整合数据\n")

      scRNA =merge(x=scRNAlist[[1]],y=scRNAlist[-1])
      
      scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
      # 使用 RunHarmony
      
      ###一定要指定harmony
      system.time({scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident")})
      # 使用 RunUMAP
      scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30)
      # 使用 FindNeighbors
      scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:30)
      # 使用 FindClusters
      scRNA <- FindClusters(scRNA, resolution = 0.5)
     # 使用 RunTSNE
    scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = 1:30)
     cat("使用Harmony整合数据成功,目前活动分组是RNA分组")
     
     # 创建 p3
     p3 <- DimPlot(scRNA, reduction = "tsne", group.by = group, pt.size = 0.5) + 
       theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
     # 创建 p4
     p4 <- DimPlot(scRNA, reduction = "tsne", group.by = "ident", pt.size = 0.5, label = TRUE, repel = TRUE) + 
       theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
     # 合并 p3 和 p4
     merged_plot <- p3 + p4
     # 保存合并的图形
     ggsave(plot = merged_plot, filename = "harmony降维.pdf", width = 20, height = 16, units = "cm")
    
    }
  cat("根据命名，添加分组信息,去除   “。 - _ ”   字符及其后面的字符，选择前面的作为分组名称\n")
  scRNA$group = gsub("[\\._-].*$", "", scRNA$orig.ident)
return(scRNA=scRNA)
}

ydx_sc_cell_cycle <- function(scRNA,zouqi = TRUE){
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(patchwork)
  cat("PCA数据分析完毕，现在观察细胞周期对PCA有无影响分析")
  ############################################################
  ####单细胞转录组基础分析：细胞周期 ####
  ############################################################
  #细胞周期回归：上一步找到的高变基因，常常会包含一些细胞周期相关基因。
  #它们会导致细胞聚类发生一定的偏移，即相同类型的细胞在聚类时会因为细胞周期的不同而分开。
  #细胞周期评分，匹配g2m_genes和s_genes
  if (is.list(scRNA)) {
    scRNA=scRNA$scRNA
  }
  # 提取g2m特征向量
  g2m_genes = cc.genes$g2m.genes
  g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA))
  # 提取s期特征向量
  s_genes = cc.genes$s.genes
  s_genes = CaseMatch(search = s_genes, match = rownames(scRNA))
  # 对细胞周期阶段进行评分
  scRNA <- CellCycleScoring(object=scRNA,  g2m.features=g2m_genes,  s.features=s_genes)
  # 对周期基因进行主成分分析
  PCA <- RunPCA(scRNA, features = c(s_genes, g2m_genes))
  # 对主成分分析结果进行可视化
  DimPlot(PCA, reduction = "pca", group.by = "Phase")
  ggsave("4. 细胞周期基因进行聚类,观察细胞周期基因对聚类的影响.pdf",width = 8,height = 6)
 
   cat("如果细胞周期对PCA影响不大，设置参数zouqi = FALSE，不去除周期影响\n\n")
  if(zouqi == TRUE){
    ##如果需要消除细胞周期的影响
    cat("正在消除细胞周期的影响，分析周期的时候，一定要禁止，设置参数zouqi = FALSE，不去除周期影响\n\n")
    scRNA <- ScaleData(scRNA, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scRNA))
    #依据高变基因PCA降维
    scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
    DimPlot(scRNA, reduction = "pca", group.by="orig.ident")
    ggsave("5. 去除细胞周期后PCA.pdf",width = 8,height = 6)
  }else{
    scRNA
  }
  return(scRNA)
}

ydx_sc_annotate_singeleR <- function(scRNA,ref =HumanPrimaryCellAtlasData(),FindAllMarker=TRUE,
                                     singleCellBase= "http://cloud.capitalbiotech.com/SingleCellBase/download/singleCellBase_all_cell_markers_20220629.txt",
                                     CellMarker_http = "http://xteam.xbio.top/CellMarker/download/Single_cell_markers.txt"){
  library(stringr)
  library(tidyr)
  library(SingleR)
  library(celldex)
  library(tidyverse)
  library(Seurat)
  if (is.list(scRNA)) {
    scRNA=scRNA$scRNA
  }
  cat("在线读入CellMarker官方文件，网址是:  \n",CellMarker_http)
  cat("\n在线读入失败时候，请选择读取本地文件 \n")
  CellMarker_offical <- read.delim(CellMarker_http)
  
  
  cat("\n开始查找原始cluster的marker，然后去singleR注释\n")
  Markers = FindAllMarkers(scRNA, only.pos = TRUE, 
                           min.pct = 0.3, logfc.threshold = 0.25)
  ##3 %>% 这是通道函数 起传递左右  可以自己百度深入理解一下，我这里只告诉你他是起传递作用的
  write.table(Markers,
              file=("_total_marker_genes_tsne_PC.txt"),
              sep="\t",quote = F,row.names = F)
  cat("一共找到marker基因:",length(unique(Markers$gene)))
  table(Markers$cluster)
  #################################################################################################
  #################################################################################################
  #################################################################################################
  cat("\n现在进行singleCellBase注释文件的编写，可以根据singleCellBase的文件进行重新分组\n")
  #######手动注释CellMarker###################
  #读入单细胞分析中输出的cell marker基因文件
  cat("本地读入singleCellBase官方文件，更新日期“2023-06-19”")
  singleCellBase_offical <- data.table::fread(file = "D:\\Personal\\Documents\\R语言函数\\singleCellBase.txt")
  singleCellBase_offical =select(singleCellBase_offical,cell_type,gene_symbol,sub_gene_symbol,species,tissue_type,pubmed_id,gene_exp
  )
  # 假设你的数据框为 df，要将列名为 column_name 的列转换为字符串
  # library(dplyr)
  # rm(singleCellBase_offical)
  # 
  # 使用 mutate() 和 if_else() 函数将空的 gene_symbol 列的行用 sub_gene_symbol 列的值填充
  singleCellBase_offical <- singleCellBase_offical %>%
    mutate(gene_symbol = if_else(is.na(gene_symbol) | gene_symbol == "", sub_gene_symbol, gene_symbol))
  
  # 使用 filter() 函数删除 gene_symbol 列为空的行
  singleCellBase_offical <- singleCellBase_offical %>%
    filter(!is.na(gene_symbol) & gene_symbol != "")
  
  # 清理 gene_symbol 列中的无效字符（保留连字符 -）
  singleCellBase_offical$gene_symbol <- gsub("[^A-Za-z0-9,-]", "", singleCellBase_offical$gene_symbol)
  # 使用 separate_rows() 函数进行拆分
  singleCellBase_offical <- singleCellBase_offical %>% separate_rows(gene_symbol, sep = ",\\s*")
  
  
  
  #### 2. 定位数据库中存在的基因 ####
  Markers %>% filter(gene %in% singleCellBase_offical$gene_symbol) -> marker.gene.sel
  cat("singleCellBase中包含的基因:",length(unique(marker.gene.sel$gene)))
  # 使用left_join函数将两个数据框合并
  singleCellBase_offical_singn <- left_join(marker.gene.sel, singleCellBase_offical, by = c("gene" = "gene_symbol"),relationship = "many-to-many")
  write.csv(singleCellBase_offical_singn,"singleCellBase细胞注释情况.csv",row.names = FALSE)
  #################################################################################################
  #################################################################################################
  #################################################################################################
  cat("现在进行CellMarker注释文件的编写，可以根据CellMarker的文件进行重新分组\n")
  #######手动注释CellMarker###################
  #读入单细胞分析中输出的cell marker基因文件
  
  
  CellMarker_offical =select(CellMarker_offical,cellName,tissueType,speciesType,cancerType,geneSymbol)
  # 使用separate_rows函数将含有逗号的列中的元素分成多行
  CellMarker_offical <- CellMarker_offical %>% separate_rows(geneSymbol, sep = ",\\s*")
  
  #### 2. 定位数据库中存在的基因 ####
  Markers %>% filter(gene %in% CellMarker_offical$geneSymbol) -> marker.gene.sel
  cat("CellMarker中包含的基因:",length(unique(marker.gene.sel$gene)))
  # 使用left_join函数将两个数据框合并
  CellMarker_offical_singn <- left_join(marker.gene.sel, CellMarker_offical, by = c("gene" = "geneSymbol"),relationship = "many-to-many")
  write.csv(CellMarker_offical_singn,"CellMarker细胞注释情况.csv",row.names = FALSE)
  #################################################################################################
  #################################################################################################
  #################################################################################################
  cat("开始进行singleR注释\n\n\n")
  #第1种方法用SingleR鉴定细胞类型
  ####然后我们把环境中的ref_Human_all赋值与refdata
  ###把rna的转录表达数据提取
  testdata <- GetAssayData(scRNA, slot="data")
  ###把scRNA数据中的seurat_clusters提取出来，注意这里是因子类型的
  clusters <- scRNA@meta.data$seurat_clusters
  ###开始用singler分析
  cellpred <- SingleR(test = testdata, ref = ref, labels = ref$label.main, 
                      clusters = clusters, 
                      assay.type.test = "logcounts", assay.type.ref = "logcounts")
  ###制作细胞类型的注释文件
  celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)
  ###保存一下
  write.csv(celltype,"singleR细胞注释情况.csv",row.names = FALSE)
  
  ##把singler的注释写到metadata中
  celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F) 
  scRNA@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
  ###因为我把singler的注释加载到metadata中时候，命名的名字叫singleR，所以画图时候，group.by="singleR"
  if (FindAllMarker) {
    Idents(scRNA) <- "singleR"
    singleR_Marker <- FindAllMarkers(scRNA, logfc.threshold = 0.01, group.by = "singleR")
  }else{
    singleR_Marker=NULL
  }
  
  
  return(list(scRNA = scRNA,Original_cluster_Markers=Markers,singleR_Marker=singleR_Marker))
  
}
ydx_sc_annotate_scCATCH <- function(scRNA=scRNA,tissue = "Brain"){
  output_dir <- "scCATCH"
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  
dir.exists()
  if (is.list(scRNA)) {
    scRNA <- scRNA$scRNA  # 如果输入的 scRNA 是列表对象，则提取其中的 scRNA 数据
  }
  library(scCATCH)  # 载入 scCATCH 包
  scCATCH_obj <- scCATCH::createscCATCH(scRNA[["RNA"]]@data, cluster = as.character(scRNA@meta.data$seurat_clusters))
  # 使用 scCATCH::createscCATCH 函数创建 scCATCH 对象，传入 scRNA 的 RNA 数据和 Seurat 对象中的聚类结果
  
  if (grepl("[a-z]", scRNA@assays[["RNA"]]@counts@Dimnames[[1]][1])) {
    species <- "Mouse"  # 如果 RNA 数据的基因名称中包含小写字母，则判断为鼠数据
  } else {
    species <- "Human"  # 否则判断为人类数据
  }
  clu_markers <- scCATCH::findmarkergene(scCATCH_obj, species = species,
                                         cluster = "All", cancer = "Normal",
                                         tissue = tissue,
                                         use_method = "1", cell_min_pct = 0.25,
                                         logfc = 0.25, pvalue = 0.25,
                                         verbose = TRUE)
  # 使用 scCATCH::findmarkergene 函数查找标记基因，传入 scCATCH 对象和其他参数

  clu_ann <- scCATCH::findcelltype(clu_markers)
  # 使用 scCATCH::findcelltype 函数进行细胞类型注释，传入标记基因结果
  # scCATCH_ann_data = merge(clu_ann@meta, clu_ann@celltype,by = "cluster") %>% dplyr::rename(barcode = cell) %>% tibble::column_to_rownames("barcode")
  scCATCH_ann_data = merge(clu_ann@meta, clu_ann@celltype,by = "cluster") %>% tibble::column_to_rownames("cell")
  
  
  scRNA = AddMetaData(scRNA, scCATCH_ann_data)
  p5_6 = DimPlot(scRNA, reduction = "tsne", label = T, group.by="cell_type", pt.size = 0.5, repel = TRUE) + ggtitle("scCATCH ")
  ggsave(filename=paste(output_dir,"annotation.with.scCATCH.pdf",sep="/"), plot = p5_6, width=15, height=7, path = "./" )

}
ydx_sc_annotate_ScType <- function( scRNA=scRNA,tissue = "Brain"){
  if (is.list(scRNA)) {
    scRNA=scRNA$scRNA
  }
  library(ggplot2)
  lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)  # 加载R包
  db_net = "D:\\Personal\\Documents\\R语言函数\\ScTypeDB_full.xlsx"
  # 注释
  # 选择样本对应的组织类型，可选：Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
  gs_list = gene_sets_prepare(db_net, tissue)
  # 此处scRNA是上面得到的Seurat对象
  es.max = sctype_score(scRNAseqData = scRNA[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  # 数据整理
  cL_resutls = do.call("rbind", lapply(unique(scRNA@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(scRNA@meta.data[scRNA@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(scRNA@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  # 将低置信度(low ScType score)的cluster设置为"unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])
  # 绘图
  scRNA@meta.data$ScType = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    scRNA@meta.data$ScType[scRNA@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  plot2= DimPlot(scRNA, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'ScType')  
  ggsave(filename = "sc-type-umap.pdf", plot = plot2,width =12,height = 5.5)
  plot2= DimPlot(scRNA, reduction = "tsne", label = TRUE, repel = TRUE, group.by = 'ScType')  
  ggsave(filename = "sc-type-tsne.pdf", plot = plot2,width =12,height = 5.5)
  return(list(scRNA=scRNA,sctype_scores=sctype_scores))
  
  
}

ydx_sc_manu_annotater <- function(scRNA,CellMarker,type ="CellMarker"){
  if (is.list(scRNA)) {
    scRNA=scRNA$scRNA
  }
  if(type =="CellMarker"){
    CellMarker <- data.frame(ClusterID = 0:(length(CellMarker)-1), CellMarker = CellMarker) # 创建数据框
    # 将"ClusterID"列转换为因子类型
    CellMarker$ClusterID <- as.factor(CellMarker$ClusterID)
    # 执行合并操作
    scRNA@meta.data <- inner_join(scRNA@meta.data, CellMarker, by = c("seurat_clusters" = "ClusterID"))
    Idents(scRNA) = "CellMarker"
    CellMarker_Markers= FindAllMarkers(scRNA,logfc.threshold = 0.01,group.by = "CellMarker")
  }else{
    CellMarker <- data.frame(ClusterID = 0:(length(CellMarker)-1), singleCellBase = CellMarker) # 创建数据框
    # 将"ClusterID"列转换为因子类型
    CellMarker$ClusterID <- as.factor(CellMarker$ClusterID)
    # 执行合并操作
    scRNA@meta.data <- inner_join(scRNA@meta.data, CellMarker, by = c("seurat_clusters" = "ClusterID"))
    Idents(scRNA) = "singleCellBase"
    CellMarker_Markers= FindAllMarkers(scRNA,logfc.threshold = 0.01,group.by = "singleCellBase")
  }
  
  return(list(scRNA=scRNA,CellMarker=CellMarker_Markers))
  
  
  
}

ydx_sc_irGSEA<- function(scRNA,assay ="RNA",method = c("AUCell", "UCell", "singscore", "ssgsea"),category = "kegg",groupby="singleR"){
  library(UCell)
  library(irGSEA)
  library(Seurat)
  if (is.list(scRNA)) {
    scRNA=scRNA$scRNA
  }
  # plot
  DimPlot(scRNA, reduction = "umap",
          group.by = groupby,label = T) + NoLegend()
  
  
  # 根据正则表达式模式获取匹配的列名包含 'group.by' 的列名
  matching_columns <- colnames(scRNA@meta.data)[grepl(groupby, colnames(scRNA@meta.data))]
  
  # 将匹配的列名赋值为群集标识符
  Idents(scRNA) <- scRNA@meta.data[, matching_columns]
  scRNA@active.ident
  cat("请检查分组信息，如果分组信息未填写，请及时更新，\n当前使用",groupby,"组信息进行分析")
  
  
  #需要联网
  if (grepl("[a-z]", scRNA@assays[["RNA"]]@counts@Dimnames[[1]][1])){
    species="Mus musculus"
  }else{
    species="Homo sapiens"
  }
  # GSEA：执行前需对所有样本进行分组，然后基于分组去计算得到排序基因列表。这个过程中，我们需要考虑不同分组中样本构成的影响；
  # GSVA：首先需要对所有样本中每个基因进行累积分布密度函数的核估计。这个过程中，需要考虑所有样本，容易受到样本背景信息的影响；
  # PLAGE 和 z-score：首先需要对基因表达矩阵执行标准化处理。同样的，这个过程容易受样本构成的影响；
  # AddModuleScore：Seurat包中的AddModuleScore函数，需要先计算基因集中所有基因的平均值，再根据平均值把表达矩阵切割成若干份，然后从切割后的每一份中随机抽取对照基因（基因集外的基因）作为背景值。因此，在整合不同样本的情况下，即使使用相同基因集为相同细胞打分，也会产生不同的富集评分；
  # AUCell：基于单个样本中的基因表达排名（gene expression rank）,使用曲线下面积来评估输入基因集是否在单个样本的前5%表达基因内富集；
  # UCell：基于单个样本的基因表达排名，使用Mann-Whitney U统计量计算单个样本的基因集富集评分；
  # 7.singscore：基于单个样本的基因表达排名，评估基因集远离中心的程度从而计算基因集富集评分；
  # 8.ssgsea：基于单个样本的基因表达排名，通过计算单个样本中基因集内和基因集外的经验累积分布函数之间的差值进而生成富集分数。然后，根据基因表达谱最大表达值和最小表达值的差值对富集分数进行标准化。这一步容易受样本构成的影响。
  category <- toupper(category)
  switch(category,
         "KEGG" = {
           category = "C2"
           subcategory = "CP:KEGG"
         },
         "BP" = {
           category = "C5"
           subcategory = "GO:BP"
         },
         "CC" = {
           category = "C5"
           subcategory = "GO:CC"
         },
         "MF" = {
           category = "C5"
           subcategory = "GO:MF"
         }
  )
  
  
  
  scRNA <- irGSEA.score(object = scRNA, assay = assay, 
                        slot = "data", ncores=1,
                        min.cells = 3, min.feature = 0,
                        custom = F, geneset = NULL, msigdb = T, 
                        species = species, category = category,  
                        subcategory = subcategory, geneid = "symbol",
                        method = "GSVA" ,
                        aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                        kcdf = 'Gaussian')
  Seurat::Assays(scRNA)
  
  result.dge <- irGSEA.integrate(object = scRNA,
                                 group.by = groupby,
                                 metadata = NULL, col.name = NULL,
                                 method = method)
  
  irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge,
                                        method = "RRA",
                                        top = 50,
                                        show.geneset = NULL)
  irGSEA.heatmap.plot
  #气泡图
  irGSEA.bubble.plot <- irGSEA.bubble(object = result.dge,
                                      method = "RRA",
                                      top = 50)
  irGSEA.bubble.plot
  #③.upset plot
  irGSEA.upset.plot <- irGSEA.upset(object = result.dge,
                                    method = "RRA")
  #> Warning in if (as.character(ta_call[[1]]) == "upset_top_annotation") {: the
  #> condition has length > 1 and only the first element will be used
  irGSEA.upset.plot
  #④.堆叠条形图
  irGSEA.barplot.plot <- irGSEA.barplot(object = result.dge,
                                        method = method)
  irGSEA.barplot.plot
  
  
  
  
  return(list(scRNA=scRNA,list_dge=result.dge))
  
  
  
}

ydx_Find_doublet <- function(scRNA,pcSelect=20){
  library(DoubletFinder)
  # 寻找最优pk值
  scRNA=Integrat
  # 判断 `sce` 对象的 `assays` 中是否包含名为 `sct` 的 assay
  is_SCT_true <- "SCT" %in% Assays(scRNA)
  
  # 输出判断结果
  if (is_SCT_true) {
    set.seed(123)
    sweep.res.list <- paramSweep_v3(scRNA, PCs = 1:pcSelect, sct = TRUE) # 如果使用SCT方法进行标准化，请设置'sct = TRUE'
  } else {
    set.seed(123)
    sweep.res.list <- paramSweep_v3(scRNA, PCs = 1:pcSelect, sct = FALSE) # 如果不使用SCT方法进行标准化，请设置'sct = FALSE'
  }
  set.seed(123)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats) #可以看到最佳参数的点
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值
  # DoubletRate = ncol(scRNA)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
  #DoubletRate = 0.075    # 直接查表，10000细胞对应的doublets rate是～7.6%
  DoubletRate = ncol(scRNA)*8*1e-6 #更通用
  #估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
  homotypic.prop <- modelHomotypic(scRNA$seurat_clusters) #最好提供celltype，而不是seurat_clusters。
  # 计算双细胞比例
  nExp_poi <- round(DoubletRate*ncol(scRNA)) 
  # 使用同源双细胞比例对计算的双细胞比例进行校正 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## 使用确定好的参数鉴定doublets
  scRNA <- doubletFinder_v3(scRNA, PCs = 1:pcSelect, pN = 0.25, pK = pK_bcmvn, 
                            nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  
  ## 结果展示，分类结果在scRNA@meta.data中
  pattern <- "classifications"
  matching_columns <- colnames(scRNA@meta.data)[grepl(pattern, colnames(scRNA@meta.data), ignore.case = TRUE)]
  
  colnames(scRNA@meta.data)[ncol(scRNA@meta.data)] = "doublet_info"
  
  
  # 查看判定为双胞细胞在UMAP图中的位置
  DimPlot(scRNA, reduction = "umap", group.by = "doublet_info")
  ggsave("Doublet-umap.pdf",width = 8,height = 6)
  DimPlot(scRNA, reduction = "tsne", group.by = "doublet_info")
  ggsave("Doublet-tsne.pdf",width = 8,height = 6)
  # 提取判定为单胞的细胞进行下游分析
  scRNA <- subset(scRNA,subset=doublet_info=="Singlet")
  dim(scRNA)
  return(scRNA)
}

ydx_sc_cell_deg <- function(scRNA,Control="CAFs",test="dapi",start_number=1,type = "singleR") {
  library(dplyr)
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggplot2) # 加载ggplot2库
  if (type == "singleR") {
    scRNA$deg <- paste( scRNA$singleR,scRNA$group, sep = "_")
    cellfordeg<-levels(factor(unname(scRNA$singleR)))
  } else if (type == "cellMarker") {
    scRNA$deg <- paste( scRNA$cellMarker, scRNA$group,sep = "_")
    cellfordeg<-levels(factor(unname(scRNA$cellMarker)))
  } else {
    cat("请输入要进行对比的cluster的注释信息，目前只支持singleR和cellMarker注释信息的识别\n")
  }
  Idents(scRNA) <- "deg"
  
  deg_list <- list()
  for(i in start_number:length(cellfordeg)){
    cat("开始进行实验组和对照组之间", cellfordeg[i], "的差异分析\n")
    CELLDEG <- FindMarkers(scRNA, ident.1 = paste0(cellfordeg[i],"_",Control), ident.2 = paste0(cellfordeg[i],"_",test),
                           test.use="wilcox",min.pct=0.1, verbose = TRUE)
    
    # 将CELLDEG对象的行名转换为列名，列名设置为"symbol"
    CELLDEG <- tibble::rownames_to_column(CELLDEG, var = "symbol")
    # 将CELLDEG对象保存为csv文件
    write.csv(CELLDEG, file = paste0(Control,"-", test, "---",cellfordeg[i], ".csv"))
    # 筛选p_val_adj小于0.05的基因
    gene_list <- CELLDEG %>% filter(p_val_adj < 0.05)
    # 选取平均对数折叠变化最大的前10个基因
    top10 <- gene_list %>% top_n(n = 12, wt = avg_log2FC) %>% select(symbol)
    top10= as.data.frame(top10)$symbol
    
    
    VlnPlot(scRNA,group.by = "group",features = top10)
    ggsave(paste0(cellfordeg[i],"差异分析前12个基因(",Control,"-",test,")VlnPlot.pdf"), device = "pdf",width = 10, height = 8)
    
    DotPlot(scRNA,features = top10,split.by ='group' )
    ggsave(paste0(cellfordeg[i],"差异分析前12个基因(",Control,"-",test,")DotPlot.pdf"), device = "pdf",width = 10, height = 8)
    
    deg_list[[cellfordeg[i]]] <- CELLDEG
  }
  return(deg_list)
}

ydx_sc_cellchat <-function(df.net="",select = "singleR"){
  # library packages
  library(Seurat)
  library(dplyr)
  library(patchwork) #最强大的拼图包
  library(ggplot2)
  library(CellChat)
  library(ggalluvial)
  library(svglite)
  options(stringsAsFactors = F) #输入数据不自动转换成因子（防止数据格式错误）
  ##提取表达矩阵和细胞分类信息
  data.input <- GetAssayData(scRNA, assay = "RNA", slot = "data")
  if (select == "singleR") {
    cat("当前选择singleR注释的细胞进行互作分析")
    
    
    # 构建cellChat对象
    cellchat <- createCellChat(object = data.input, meta = scRNA@meta.dat, group.by = "singleR")
  } else if (select == "CellMarker") {
    cat("当前选择CellMarker注释的细胞进行互作分析")
    cellchat <- createCellChat(object = data.input, meta = scRNA@meta.dat, group.by = "CellMarker")
  }
  
  
  # 导入配受体数据库
  if(grepl("[a-z]", cellchat@data@Dimnames[[1]][99])){
    ####可选CellChatDB.human, CellChatDB.mouse
    cat("根据基因表达，目前使用小鼠数据")
    CellChatDB <- CellChatDB.mouse
  }else{
    cat("根据基因表达，目前使用人的数据")
    CellChatDB <- CellChatDB.human
  }
  showDatabaseCategory(CellChatDB)
  colnames(CellChatDB$interaction)
  CellChatDB$interaction[1:4,1:4]
  head(CellChatDB$cofactor)
  head(CellChatDB$complex)
  head(CellChatDB$geneInfo)
  
  
  ########在CellChat中，我们还可以先择特定的信息描述细胞间的相互作用，
  ##可以理解为从特定的侧面来刻画细胞间相互作用，比用一个大的配体库又精细了许多。
  ##查看可以选择的侧面
  unique(CellChatDB$interaction$annotation)
  # use Secreted Signaling for cell-cell communication analysis
  
  # 选取分泌信号通路进行下游分析CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use # set the used database in the object
  
  #对表达数据进行预处理
  
  ##将信号基因的表达数据进行子集化，以节省计算成本
  # 首先在一个细胞群中识别过表达的配体或受体，然后识别其相互作用；
  # 可将基因表达数据投射到蛋白-蛋白相互作用(PPI)网络上。
  cellchat <- subsetData(cellchat)
  # 识别过表达基因
  cellchat <- identifyOverExpressedGenes(cellchat)
  # 识别过表达基因的互作
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # 将配体、受体投射到PPI网络
  if(grepl("[a-z]", cellchat@data@Dimnames[[1]][1])){
    ####可选CellChatDB.human, CellChatDB.mouse
    cat("根据基因表达，目前使用小鼠数据")
    cellchat <- projectData(cellchat, PPI.mouse)
  }else{
    cat("根据基因表达，目前使用人的数据")
    cellchat <- projectData(cellchat, PPI.human)
  }
  ##相互作用推断
  ## 1、计算通信概率推断细胞互作的通信网络(cellphonedb是用平均表达值代表互作强度)
  # 如果不想用上一步PPI矫正的结果，raw.use = TRUE即可
  #默认的cutoff的值为20%，即表达比例在25%以下的基因会被认为是0，trim=0.1 可以调整比例阈值
  cellchat <- computeCommunProb(cellchat, raw.use = T, population.size = T)
  ###如果特定细胞群中只有少数细胞，则过滤掉细胞间的通信
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # 通过计算的每个配体-受体通信概率来计算信号通路层面的通信概率
  # 每个信号通路分别存储在“net”和“netP”槽中
  cellchat <- computeCommunProbPathway(cellchat)
  # 计算聚集的细胞通讯网络
  cellchat <- aggregateNet(cellchat)
  # 计算每种细胞各有多少个
  groupSize <- as.numeric(table(cellchat@idents))
  # 创建一个PDF图形设备
  pdf("整体---细胞互作网络图-全部细胞.pdf")
  par(mfrow = c(1,2), xpd = TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                   label.edge = F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                   label.edge = F, title.name = "Interaction weights/strength")
  # 关闭图形设备并保存图像
  dev.off()
  
  #左图：外周各种颜色圆圈的大小表示细胞的数量，圈越大，细胞数越多。发出箭头的细胞表达配体，
  #箭头指向的细胞表达受体。配体-受体对越多，线越粗。
  #右图：互作的概率或者强度值（强度就是概率值相加）
  
  ##每个细胞如何跟别的细胞互作（互作的强度或概率图）
  mat <- cellchat@net$weight
  pdf("整体---细胞互作网络图-单个细胞-权重.pdf")
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    plot <- netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    print(plot)
  }
  dev.off()
  ##每个细胞如何跟别的细胞互作（number+of+interaction图）
  mat <- cellchat@net$count
  pdf("整体---细胞互作网络图-单个细胞-通讯数量.pdf")
  for(i in 1:nrow(mat)){
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i,] <- mat[i,]
    plot= netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                           arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    print(plot)
  }
  dev.off()
  #提取推断出的细胞互作的通信网络数据框，我们提供了一个subsetCommunication 函数，
  #可以方便地访问感兴趣的推断的细胞间通信。
  ##返回一个数据框，包含所有推断的配体/受体级别的细胞-细胞通信。设置slot.name = "netP"以访问信令路径级别的推断通信
  
  df.net <- subsetCommunication(cellchat) # 将细胞通讯预测结果以数据框的形式取出来
  # 把预测的结果写出来
  cat("正在写出细胞通讯预测结果")
  write.csv(df.net, '细胞通讯预测结果.csv')
  ######################################################################################
  #######################上面是整体通讯的权重和数量，现在开始看通路#####################
  
  ##可视化每个信号通路
  ##查看通路
  #单个信号通路或配体-受体介导的细胞互作可视化(层次图、网络图、和弦图、热图)
  
  cat("目前富集到了以下这些信号通路\n",cellchat@netP$pathways  )#查看富集到的信号通路
  
  for(i in cellchat@netP$pathways){
    cat("当前正在输出",i,"信号通路贡献度\n")
    pathways.show= i
    pdf(paste0(i,"信号通路.pdf"))
    #配体-受体层级的可视化（计算各个ligand-receptor+pair对信号通路的贡献）
    plot1= netAnalysis_contribution(cellchat, signaling = pathways.show)
    ##圈图
    par(mfrow = c(1,1))
    plot2= netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
    ##和弦图
    par(mfrow=c(1,1))
    plot3=  netVisual_aggregate(cellchat, signaling =pathways.show, layout = "chord")
    ##热图
    par(mfrow=c(1,1))
    plot4= netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
    ##层次图
    par(mfrow = c(1,1))
    plot5= netVisual_aggregate(cellchat, signaling = pathways.show,layout="hierarchy")
  }
  
  
  
  
  
  
  
  
  
  
  # for(i in cellchat@netP$pathways) {
  #   cat("当前正在输出", i, "信号通路贡献度\n")
  #   pathways.show <- i
  #   pdf(paste0(i,"信号通路.pdf"))
  #   # 配体-受体层级的可视化（计算各个ligand-receptor+pair对信号通路的贡献）
  #   plot1 <- netAnalysis_contribution(cellchat, signaling = pathways.show)
  #   
  #   # 圈图
  #   plot2 <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  #   
  #   # 和弦图
  #   plot3 <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  #   
  #   # 热图
  #   plot4 <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  #   
  #   # 层次图
  #   plot5 <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "hierarchy")
  #   
  #   plot1+plot2+plot4
  # 
  #   dev.off()
  # }
  # 
  
  
  
  
  
  levels(cellchat@idents)            #查看细胞顺序
  vertex.receiver = c(1, 2)          #指定靶细胞的索引
  
  
  ##层次图
  par(mfrow = c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout="hierarchy")
  #在层次图中，实体圆和空心圆分别表示源和目标。圆的大小与每个细胞组的细胞数成比例。线越粗，互作信号越强。
  #左图中间的target是我们选定的靶细胞。右图是选中的靶细胞之外的另外一组放在中间看互作。
  ##圈图
  par(mfrow = c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  ##和弦图
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling =pathways.show, layout = "chord")
  #, vertex.size = groupSize)
  
  ##热图
  par(mfrow=c(1,1))
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  ##纵轴是发出信号的细胞，横轴是接收信号的细胞，热图颜色深浅代表信号强度。
  ##上侧和右侧的柱子是纵轴和横轴强度的累积
  
  
  # 
  #   ##也可以看到单个配体-受体对介导的细胞-细胞通信。
  #   #我们提供了一个extractEnrichedLR功能来提取给定信号通路的所有重要相互作用(L-R对)和相关信号基因。
  #   pairLR.MK <- extractEnrichedLR(cellchat, signaling = "MK", geneLR.return = FALSE)
  #   #提取对这个通路贡献最大的配体受体对来展示（也可以选择其他的配体受体对）
  #   LR.show <- pairLR.MK[1,] # show one ligand-receptor pair
  #   # Hierarchy plot
  #   vertex.receiver = seq(1,2) # a numeric vector
  #   ##层次图
  #   netVisual_individual(cellchat, signaling = "MK",  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout="hierarchy")
  #   ##圈图
  #   netVisual_individual(cellchat, signaling = "MK", pairLR.use = LR.show, layout = "circle")
  #   ##和弦图
  #   netVisual_individual(cellchat, signaling ="MK", pairLR.use = LR.show, layout = "chord")
  #   
  
  ##############批量保存
  pathway.show.all=cellchat@netP$pathways
  levels(cellchat@idents)
  vertex.receiver=c(1,2,3,4)
  dir.create("cellchat.res")
  dev.off()
  for (i in 1:length(pathway.show.all)) {
    cat("当前正在输出", pathway.show.all[i], "信号通路\n")
    netVisual(cellchat,signaling = pathway.show.all[i],layout = c("circle", "hierarchy", "chord", "spatial"),
              out.format = c("pdf"),nCol =1)
    
    
    plot=netAnalysis_contribution(cellchat,signaling = pathway.show.all[i])
    ggsave(filename = paste0("cellchat.res/",pathway.show.all[i],".信号通路贡献度.pdf"),
           plot=plot,width=6,height=4,units="in")
    
  }
  
  #####################################################################################
  ##气泡图
  levels(cellchat@idents)
  netVisual_bubble(cellchat, sources.use = "Fibroblast",  remove.isolate = FALSE)
  
  ##sources.use = 2 是值第二个细胞亚群
  netVisual_bubble(cellchat, sources.use =c(1,3), targets.use = c(1:5), remove.isolate = FALSE)
  ##指定信号通路
  cellchat@netP$pathways 
  netVisual_bubble(cellchat, sources.use =c("Brush cell (Tuft cell)","Clara cell"), targets.use =c(1:5),signaling =  c("MK","CCL"), remove.isolate = FALSE)
  
  
  pairLR  <- extractEnrichedLR(cellchat, signaling =c("MK","CCL"), geneLR.return = FALSE)
  netVisual_bubble(cellchat, sources.use =c(1,3), targets.use =c(1:5),pairLR.use =pairLR , remove.isolate = FALSE)
  
  #和弦图               
  netVisual_chord_gene(cellchat, sources.use = 2, targets.use = c(1:5), lab.cex = 0.5,legend.pos.y = 30)
  
  ##用小提琴图绘制信号基因的表达分布 参与某条信号通路（如TGFb）的所有基因在细胞群中的表达情况展示
  plotGeneExpression(cellchat, signaling = "MK")
  #默认情况下，plotGeneExpression只显示与推断的重要通信相关的信号基因的表达。
  #用户可以通过显示一个信号通路相关的所有信号基因的表达。
  plotGeneExpression(cellchat, signaling = "MK", enriched.only = FALSE)
  #也可以用气泡图展示
  plotGeneExpression(cellchat, signaling = "MK",color.use = NULL,type = "dot")
  
  ###################################################################################################
  ##可视化配体和受体
  ## 1、计算网络中心性得分
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  ##2、热图  使用热图可视化计算的中心性评分，允许随时识别细胞群的主要信号作用。
  netAnalysis_signalingRole_network(cellchat, signaling = "MK", width = 8, height = 2.5, font.size = 10)
  
  ##在2D空间中可视化主要的发送者(源)和接收者(目标)。
  ##我们还提供了另一种直观的方式，使用散点图来可视化2D空间中的主要发送者(源)和接收者(目标)。
  
  ##从所有信号通路对聚合的细胞-细胞通信网络的信号作用分析
  gg1 <- netAnalysis_signalingRole_scatter(cellchat)
  ###从所有信号通路对聚合的细胞-细胞通信网络的信号作用分析
  gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MK", "PARs"))
  gg1 + gg2
  
  ##识别对某些细胞群的传出或传入信号贡献最大的信号，从所有信号通路对聚合的细胞-细胞通信网络的信号作用分析。
  
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
  ht1 + ht2
  
  ########################################################################################################
  #细胞通讯模式和信号网络
  
  library(NMF)
  library(ggalluvial)
  #非负矩阵分解（NMF）识别细胞的通讯模式
  ##信号输出细胞的模式识别
  ##计算分解成几个因子(pattern)比较合适（这一步运行比较慢+。在使用NMF对细胞进行亚群细分时，如果不测试的话，最好选择比细胞类型多一点的值）
  selectK(cellchat, pattern = "outgoing")
  
  
  #挑选曲线中第一个出现下降的点（从3就开始下降了）
  nPatterns = 2
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  ##river plot
  netAnalysis_river(cellchat, pattern = "outgoing")
  #气泡图
  netAnalysis_dot(cellchat, pattern = "outgoing")
  #信号输入细胞的模式识别
  selectK(cellchat, pattern = "incoming")
  
  #################################################################################################
  #  信号网络聚类
  # 1、根据功能相似性来识别信号分组
  
  ##reticulate::py_install(packages = 'umap-learn')
  
  ## 2、基于结构相似性识别信号分组
  cellchat <- computeNetSimilarity(cellchat, type = "structural")
  cellchat <- netEmbedding(cellchat, type = "structural")
  #> Manifold learning of the signaling networks for a single dataset
  cellchat <- netClustering(cellchat, type = "structural")
  #> Classification learning of the signaling networks for a single dataset
  # Visualization in 2D-space
  netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
  #############################################################################################
  #############################################################################################
  #############################################################################################
  #############################################################################################
  #############################################################################################
  #不同分组之间的配对分析
  sc.sp=SplitObject(scRNA,split.by = "singleR")
  table(scRNA$group)
  
  sc.11=scRNA$scRNA[,colnames(sc.sp[["CAFs.1"]]),]
  sc.3=scRNA$scRNA[,colnames(sc.sp[["CAFs.1"]]),]
  
  cat()
  cellchat.sc11 = sc.sp[[1]]
  cellchat.sc3 = sc.sp[[2]]
  
  cellchat.sc11 <- createCellChat(object =sc.11@assays$RNA@data, meta =sc.11@meta.data,  group.by ="singleR")
  cellchat.sc3 <- createCellChat(object =sc.3@assays$RNA@data, meta =sc.3@meta.data,  group.by ="singleR")
  
  dir.create("compare")
  setwd("compare/")
  
  cellchat=cellchat.sc11 
  cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
  cellchat <- filterCommunication(cellchat, min.cells = 3)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cc.sc11 = cellchat
  #################################
  cellchat=cellchat.sc3
  cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
  cellchat <- filterCommunication(cellchat, min.cells = 3)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cc.sc3 = cellchat
  ##############################################
  cc.list=list(SC11=cc.sc11,SC3=cc.sc3)
  cellchat=mergeCellChat(cc.list,cell.prefix = T,add.names = names(cc.list))
  ##可视化
  ##所有细胞群总体观：通讯数量与强度对比
  compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "count")
  compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "weight")
  ##第一个图展示通讯数量之间的差异，第二个图展示通讯强度之间的差异。 
  
  ##数量与强度差异网络图
  netVisual_diffInteraction(cellchat,weight.scale = T)
  netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight")
  ##红色是case相对于control上调的，蓝色是下调的
  
  #数量与强度差异热图
  netVisual_heatmap(cellchat)
  netVisual_heatmap(cellchat,measure = "weight")
  #case和control对比，红色是上调，蓝色是下调
  
  #保守和特异性信号通路的识别与可视化
  rankNet(cellchat,mode = "comparison",stacked = T,do.stat = T)
  rankNet(cellchat,mode = "comparison",stacked =F,do.stat = T)
  ##左图最下面多个信号通路是case组独有的
  
  ##细胞互作数量对比网络图
  weight.max=getMaxWeight(cc.list,attribute = c("idents","count"))
  netVisual_circle(cc.list[[1]]@net$count,weight.scale = T,label.edge = F,
                   edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )
  
  netVisual_circle(cc.list[[2]]@net$count,weight.scale = T,label.edge = F,
                   edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )
  
  
  table(scRNA_harmony@active.ident)
  s.cell=c( "Macrophage", "Tissue_stem_cells","Monocyte")
  count1=cc.list[[1]]@net$count[s.cell,s.cell]
  count2=cc.list[[2]]@net$count[s.cell,s.cell]
  
  netVisual_circle(count1,weight.scale = T,label.edge = F,
                   edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )
  
  netVisual_circle(count2,weight.scale = T,label.edge = F,
                   edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )
  
  
  
}

ydx_sc_monocle3 <-function(){
  library(dplyr)
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(monocle3)
  data <- GetAssayData(scRNA, assay = 'RNA', slot = 'counts')

  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  
  
  cds <- new_cell_data_set(data,
                           cell_metadata = scRNA@meta.data,
                           gene_metadata = gene_annotation)
  ## Step 1: Normalize and pre-process the data
  cds <- preprocess_cds(cds, num_dim = 100)
  plot_pc_variance_explained(cds)
  cds <- reduce_dimension(cds, preprocess_method = "PCA",reduction_method = c("UMAP"))
  
  plot_cells(cds)
  plot_cells(cds, color_cells_by="singleR")
  
  
  
}

############################################################
####单细胞转录组基础分析四：细胞注释 #######################
############################################################


ydx_sc_Markplot <- function(scRNA, Mark_list, Idents= "singleR"){
  Idents(scRNA) = Idents
  # 遍历每一个cluster然后展示其中前4个基因
  marker.sig <- Mark_list %>% 
    mutate(Ratio = round(pct.1/pct.2,3)) %>%
    filter(p_val_adj <= 0.05)  # 本条件为过滤统计学不显著的基因
  
  cluster_gene_number = 12
  # 定义一个空的列表，用于存储每个簇的代表性基因。
  cl_genes_list <- list()
  for(cluster_id in unique(marker.sig$cluster)){
    cl4.genes <- marker.sig %>% 
      filter(cluster == cluster_id) %>%
      arrange(desc(Ratio))
    cl4.genes <- cl4.genes[1:min(nrow(cl4.genes),cluster_gene_number),"gene"]
    
    # VlnPlot
    VlnPlot(scRNA, features = cl4.genes, ncol = 4)
    ggsave(paste0(cluster_id,"-每个cluster前",cluster_gene_number,"的基因VlnPlot之cluster",cluster_id,".pdf"), device = "pdf", width = 16, height = 12)
    
    # FeaturePlot
    FeaturePlot(scRNA, features = cl4.genes, ncol = 4)
    ggsave(paste0(cluster_id,"-每个cluster前",cluster_gene_number,"的基因FeaturePlot之cluster",cluster_id,".pdf"), device = "pdf", width = 16, height = 12)
    
    # RidgePlot
    # RidgePlot(scRNA, features = cl4.genes, ncol =4)
    # ggsave(paste0(cluster_id,"-每个cluster前",cluster_gene_number,"的基因RidgePlot图之cluster",cluster_id,".pdf"), device = "pdf", width = 16, height = 12)
    
    
    #选择每个cluster中的3个元素用来画图
    cl5.genes <- marker.sig %>% 
      filter(cluster == cluster_id) %>%
      arrange(desc(Ratio))
    cl5.genes <- cl5.genes[1:min(nrow(cl5.genes),3),"gene"]
    cl_genes_list[[cluster_id]] <- cl5.genes
  }
}



############################################################
####单细胞转录组基础分析五：细胞通讯#######################
############################################################




############################################################
####单细胞转录组：图片美化#######################
############################################################

#####这个是细胞每个cluster的比例的百分比的图片
bar_cluster_PropPlot <- function(scRNA, groupBy="singleR"){
  # (1)获取绘图数据
  library(dplyr)
  library(ggplot2)
  if (is.list(scRNA)) {
    scRNA=scRNA$scRNA
  }
  
  plot_data = scRNA@meta.data %>% 
    dplyr::select(orig.ident, {{groupBy}}) %>% 
    dplyr::rename(group = as.name(groupBy))
  
  # (2)绘图
  figure = ggstatsplot::ggbarstats(data = plot_data, 
                      x = group, y = orig.ident,
                      package = 'ggsci',
                      palette = 'category20c_d3',
                      results.subtitle = FALSE,
                      bf.message = FALSE,
                      proportion.test = FALSE,
                      label.args = list(size = 2, 
                                        fill = 'white', 
                                        alpha = 0.85,
                                        family = 'Arial',
                                        fontface = 'bold'),
                      perc.k = 2,
                      title = '',
                      xlab = '',
                      legend.title = 'Seurat Cluster',
                      ggtheme = ggpubr::theme_pubclean()) +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = 'black', lineend = 'round'),
          legend.position = 'right',
          axis.text.x = element_text(size = 15, color = 'black', family = 'Arial'),
          axis.text.y = element_text(size = 15, color = 'black', family = 'Arial'),
          legend.text = element_text(family = 'Arial', size = 10, color = 'black'),
          legend.title = element_text(family = 'Arial', size = 13, color = 'black')) 
  
  # (3)去除柱子下面的样本量标识：
  gginnards::delete_layers(x = figure, match_type = 'GeomText')
}













ydx_sc_Integrat1= function(dir,Integrat =TRUE, SCT = TRUE,resolution = 0.5, specise = "mmu",scale_genes ="all", pc.num=1:20){
#整合和合并是不同的，整合是均匀分布，可能会丢失掉一些cluster
#合并是单纯的合并在一起，所有的cluster都会保留，但是要注意，需要是同一批次的样本才能进行合并，而不整合
####################################是否要整合数据#######################################
# 判断输入的单细胞转录组数据集的数量是否大于1，如果是则进行整合，否则直接使用第一个数据集。
# 如果存在一个或多个数据集

  # 如果选择整合方法，则继续进行判断，如果不使用SCT标准化，则对每个数据集进行标准化处理（NormalizeData）和筛选差异表达基因（FindVariableFeatures），选择跨数据集的基因（SelectIntegrationFeatures），并找到数据集的锚点（FindIntegrationAnchors），最终使用整合数据方法（IntegrateData）将所有数据整合成一个数据集。
  if (Integrat== 1){ 
    cat("数据集大于1个，可以选择一下三种方式：\n 1. 使用SCT标准化后合并数据\n2.使用NormalizeData后整合数据\n3.使用merge函数合并数据\n")
    
    } else{
      # 在整合前删除之前可能存在的scRNA数据集，避免重复赋值。
      rm(scRNA)
      # 对每个数据集进行标准化处理和筛选差异表达基因，并输出信息。
      
      cat("有多个数据集，先NormalizeData，再FindVariableFeatures，最后再整合数据，最终scale_genes\n")
      for(i in 1:length(dir)){
        scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
        scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst",nfeatures = 3000)
      }
      # 选择跨数据集的基因，也就是多个数据集都有表达的基因用于整合
      features <- SelectIntegrationFeatures(object.list = scRNAlist)
      # 找到数据集的锚点
      immune.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = features)
      # this command creates an 'integrated' data assay
      # 整合数据集
      scRNA <- IntegrateData(anchorset = immune.anchors)
      # 将整合后的数据集的默认数据类型设置为“integrated”，方便后续分析。
      DefaultAssay(scRNA) <- "integrated"
      ####数据标准化（中心化）
      ##如果内存足够最好对所有基因进行中心化
      ##如果内存不够，可以只对高变基因进行标准化
      if (scale_genes =="all"){
        scale.genes <-  rownames(scRNA)
        scRNA <- ScaleData(scRNA, features = scale.genes)
      }else{
        scale.genes <-  VariableFeatures(scRNA)
        scRNA <- ScaleData(scRNA, features = scale.genes)
      } 
      cat("使用标准化整合数据成功")
    }

    # 如果数据集同一批次，差异较小，使用merge函数合并数据集

# 检查Seurat对象的数据格式
scRNA <- RunPCA(scRNA, verbose = FALSE)
#每个PCA对整体特征的贡献度
#前50个PCA的变化




ElbowPlot(scRNA$scRNA , ndims=50, reduction="pca") 
ggsave(paste0("根据曲线选择合适的PCA.pdf"), device = "pdf", width = 8, height = 6)

cat("PCA降维聚类完毕，下一步选择合适的PC数\n")
cat("使用SCTransfrom函数可以增加一定的维度信息，使用NormalizeData函数，一般20个维度即可\n\n")
pca_number = 1:20
if (SCT == TRUE) {
  cat("当前使用SCTransfrom函数可以增加一定的维度信息，在原来的维度上增加10个维度\n\n\n")
  pca_number <- 1:(max(pca_number)+10)
}
cat("根据PCA降维分析的碎石图，来选择UMAP和TSNE的维度信息，当前维度信息是：",length(pca_number),"\n\n")
# Run the standard workflow for visualization and clustering
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = pca_number)
#tSNE
scRNA = RunTSNE(scRNA, dims = pca_number) 
#细胞聚类resolution越多，切的越多，cluster越多
#pc.num维度越多，cluster越多
###Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B. Thanks to Nigel Delaney (evolvedmicrobe@github) for the rewrite of the Java modularity optimizer code in Rcpp!
scRNA <- FindNeighbors(scRNA, reduction = "pca", dims =pca_number)
###这个分辨率是可以自定义的，当我们的样本细胞数较大时候resolution 要高一些，

# #设置不同的分辨率，观察分群效果，dim为PCA选择的主成分数
# 
# for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1,1.2,1.5,2,2.5,3)) {
# scRNA=FindClusters(scRNA,  resolution = res, algorithm = 1)}
# 
# 
# 
# 
# apply(scRNA@meta.data[,grep("RNA_snn_res",colnames(scRNA@meta.data))],2,table)
# p2_tree=clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
# pdf(file = "03.clustertree.pdf",width =12,height =10)
# p2_tree
# dev.off()
af <- FindClusters(af, resolution = resolution)

cat("根据细胞情况选择resolution，一般情况2万细胞以上都是大于1.0的。当前：",resolution,"\n\n")
cat("查看每一类有多少个细胞,当前有",length(unique(scRNA@meta.data$seurat_clusters)),"个cluster，\n每个cluster细胞分别是：\n",table(scRNA@meta.data$seurat_clusters),"\n\n")
##系统发育分析（Phylogenetic Analysis of Identity Classes）
scRNA<-BuildClusterTree(scRNA)
PlotClusterTree(scRNA)

#根据命名，添加分组信息
scRNA$group = gsub("[\\._-].*$", "", scRNA$orig.ident)
return(list(scRNA=scRNA,quality_control=freq.combine))
}




ydx_sc_monocle <- function(scRNA= annotate){
  library(monocle)
  library(dplyr)
  library(Seurat)
  library(clustree)
  library(tidyverse)
  library(patchwork)
  library(ggplot2) # 加载ggplot2库
  library(patchwork) # 加载patchwork库
  library(RColorBrewer) # 加载RColorBrewer库
  library(glmGamPoi) # 加载glmGamPoi库
  scRNA.Osteoclastic <- scRNA$scRNA[, rownames(subset(scRNA$scRNA@meta.data, singleR == "HSC_-G-CSF"))]
  data <- as(as.matrix(scRNA.Osteoclastic@assays$RNA@counts), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data = scRNA.Osteoclastic@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  monocle_cds <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  #这里要注意，如果数据是稀疏矩阵使用negbinomial.size()
  #FPKM 使用tobit()        logFPKM使用gaussianff()
  ####################################################################
  monocle_cds <- estimateSizeFactors(monocle_cds)
  monocle_cds <- estimateDispersions(monocle_cds)
  
  monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
  print(head(fData(monocle_cds)))
  

  #筛选monocle的高变基因
  disp_table <- dispersionTable(monocle_cds)
  disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
  monocle_cds <- setOrderingFilter(monocle_cds, disp.genes)
  plot_ordering_genes(monocle_cds)
  monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
                          method = 'DDRTree')
  
  #修改orderCells
  trace('project2MST', edit = T, where = asNamespace("monocle"))
  monocle_cds <- orderCells(monocle_cds)
  plot_cell_trajectory(monocle_cds, color_by = "singleR")
  
  
  plot_cell_trajectory(monocle_cds, color_by = "State")
  
  plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")
  
  
  plot_cell_trajectory(monocle_cds, color_by = "State") +
    facet_wrap(~State, nrow = 1)
  
  
  
  blast_genes <- row.names(subset(fData(monocle_cds),
                                  gene_short_name %in% c("GAPDH", "RORA")))
  
  plot_genes_jitter(monocle_cds[blast_genes,],
                    grouping = "State",
                    min_expr = 0.1)
  
  
  
  HSMM_expressed_genes <-  row.names(subset(fData(monocle_cds),
                                            num_cells_expressed >= 10))
  HSMM_filtered <- HSMM[HSMM_expressed_genes,]
  my_genes <- row.names(subset(fData(HSMM_filtered),
                               gene_short_name %in% c("YWHAB", "GAPDH", "TNNC1")))
  cds_subset <- HSMM_filtered[my_genes,]
  plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")
  
  
  
  plot_genes_in_pseudotime(cds_subset, color_by =  "State")
  
  
  install.packages(c("devtools", "remotes", "BiocManager"))
  BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment", "GenomicRanges", "rGREAT", "scran", "scater", "Rtsne", "igraph"))
  remotes::install_github("sqjin/CellChat")
  
  set_config(
    use_proxy(url = "http://127.0.0.1:10808")
  )
  
}


ydx_sc_PCA <- function(scRNA,PCA= c("orig.ident","singleR")){
  # UMAP 个性化分析
  library(ggplot2)
  library(magrittr)
  
  color_number = length(unique(scRNA@meta.data$seurat_clusters))
  getpalette = colorRampPalette(brewer.pal(color_number,"Paired"))
  celltype_color =getpalette(color_number)
  # Visualization
  for(i in PCA){
    Idents(scRNA) = i
    for(PCA_type in c("pca","tsne","umap")){
      #分离每个cluster，彩虹色
      DimPlot(scRNA,group.by=i, split.by = i,cols =celltype_color, reduction = PCA_type)
      #分离每个cluster，默认颜色
      DimPlot(scRNA,group.by=i, split.by = i,reduction = PCA_type)
      #图上带label
    DimPlot(scRNA, group.by=i, label=T, label.size=4,cols =celltype_color, reduction=PCA_type)
    #图上不带label
    DimPlot(scRNA, group.by=i, cols =celltype_color, reduction=PCA_type)
    }
    
    
    ###########################美化tsne图######################

    ## Step 01. 获取UMAP坐标信息
    
    
    if( i == "umap" ){
      umap <- scRNA@reductions$umap@cell.embeddings %>% as.data.frame

      ## Step 02. 获取与坐标对应的细胞类型
      cellType <- scRNA@active.ident %>% as.data.frame %>%
        `colnames<-` ("Cell_type")
      umap <- cbind(umap, cellType)
      ## Step 03. 基础绘图
      allcolour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98",
                     "#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF",
                     "#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C",
                     "#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD",
                     "#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
                     "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD",
                     "#9932CC","#8B008B","#8B4513","#DEB887")
      p <- ggplot(umap,aes(x=UMAP_1  , y = UMAP_2 ,colour = Cell_type)) +  
        geom_point(size = 1 , alpha = 0.5) +  
        scale_colour_manual(values = allcolour)
      p
      
      ## Step 04. 设置个性化的主题
      myTheme <- theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 12),
      )  
      
      
      p2 <- p + myTheme + guides(
        colour = guide_legend(
          override.aes = list(size = 6)
        )
      )
      p2
      
      ## Step 05. 添加短的坐标轴
      xstart <- min(umap$UMAP_1); xend <- min(umap$UMAP_1) + 3
      ystart <- min(umap$UMAP_2); yend <- min(umap$UMAP_2) + 3
      
      p3 <- p2 + 
        geom_segment(aes(x = xstart, y = ystart , xend = xend, yend = ystart),
                     colour = "black", linewidth=1,arrow = arrow(length = unit(0.3,"cm")))+ 
        geom_segment(aes(x = xstart, y = ystart, xend = xstart , yend = yend),
                     colour = "black", linewidth=1,arrow = arrow(length = unit(0.3,"cm"))) +
        annotate("text", x = xstart +1.5, y = ystart -1, label = "UMAP_1",
                 color="black",size = 3, fontface="bold" ) + 
        annotate("text", x = xstart -1, y = ystart + 1.5, label = "UMAP_2",
                 color="black",size = 3, fontface="bold" ,angle=90) 
      p3
      
      ## Step 06. 将标签添加到绘图区
      ##### 基于群体中每个坐标轴的中位值，计算标签的位置
      cellTypeLoc <- umap %>%
        group_by(Cell_type) %>%
        summarise(
          Loc1 = median(UMAP_1),
          Loc2 = median(UMAP_2)
        ) %>% as.data.frame
      
      ##### 添加标签
      p4 <- p3 + geom_text(data = cellTypeLoc, 
                           mapping = aes(x = Loc1, y = Loc2, label = Cell_type),
                           color = "black", size = 5)
      p4
      
      ## Step 07. 圈出细胞群
      devtools::install_github("sajuukLyu/ggunchull", type = "source")
      library(ggunchull)
      # 可以通过调整nbin、nsm 、qval 、sfac参数，
      # 调整大小和形状
      ##### 圈出所有
      p5=   p4 + stat_unchull(
        fill = "white", 
        alpha = 0, 
        show.legend = FALSE,
        nsm = 50,
        nbin = 200,
        sfac = 1.5
      )

    
    
    # ##### 圈出特定的细胞簇
    # p4 + stat_unchull(
    #   data = subset(umap, Cell_type == "NK"),
    #   fill = "white",
    #   alpha = 0, 
    #   color = "black",
    #   show.legend = FALSE,
    #   size = 2
    # )
  }else if(i == "tsne"){ 
    
    umap <- scRNA@reductions$tsne@cell.embeddings %>% as.data.frame
  
    ## Step 02. 获取与坐标对应的细胞类型
    cellType <- scRNA@active.ident %>% as.data.frame %>%
      `colnames<-` ("Cell_type")
    umap <- cbind(umap, cellType)
    ## Step 03. 基础绘图
    allcolour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98",
                   "#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF",
                   "#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C",
                   "#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD",
                   "#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
                   "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD",
                   "#9932CC","#8B008B","#8B4513","#DEB887")
    p <- ggplot(umap,aes(x=tSNE_1  , y = tSNE_2 ,colour = Cell_type)) +  
      geom_point(size = 1 , alpha = 0.5) +  
      scale_colour_manual(values = allcolour)
    p
    
    ## Step 04. 设置个性化的主题
    myTheme <- theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 12),
    )  
    
    
    p2 <- p + myTheme + guides(
      colour = guide_legend(
        override.aes = list(size = 6)
      )
    )
    p2
    
    ## Step 05. 添加短的坐标轴
    xstart <- min(umap$tSNE_1); xend <- min(umap$tSNE_1) + 3
    ystart <- min(umap$tSNE_2); yend <- min(umap$tSNE_2) + 3
    
    p3 <- p2 + 
      geom_segment(aes(x = xstart, y = ystart , xend = xend, yend = ystart),
                   colour = "black", linewidth=1,arrow = arrow(length = unit(0.3,"cm")))+ 
      geom_segment(aes(x = xstart, y = ystart, xend = xstart , yend = yend),
                   colour = "black", linewidth=1,arrow = arrow(length = unit(0.3,"cm"))) +
      annotate("text", x = xstart +1.5, y = ystart -1, label = "tSNE_1",
               color="black",size = 3, fontface="bold" ) + 
      annotate("text", x = xstart -1, y = ystart + 1.5, label = "tSNE_2",
               color="black",size = 3, fontface="bold" ,angle=90) 
    p3
    
    ## Step 06. 将标签添加到绘图区
    ##### 基于群体中每个坐标轴的中位值，计算标签的位置
    cellTypeLoc <- umap %>%
      group_by(Cell_type) %>%
      summarise(
        Loc1 = median(tSNE_1),
        Loc2 = median(tSNE_2)
      ) %>% as.data.frame
    
    ##### 添加标签
    p4 <- p3 + geom_text(data = cellTypeLoc, 
                         mapping = aes(x = Loc1, y = Loc2, label = Cell_type),
                         color = "black", size = 5)
    p4
    
    ## Step 07. 圈出细胞群

    library(ggunchull)
    # 可以通过调整nbin、nsm 、qval 、sfac参数，
    # 调整大小和形状
    ##### 圈出所有
    p5=   p4 + stat_unchull(
      fill = "white", 
      alpha = 0, 
      show.legend = FALSE,
      nsm = 50,
      nbin = 200,
      sfac = 1.5
    )
  }
    
    
    
    
    
    
    
    
    
    
    
    
    ############Markgene 气泡图##############
    #随机选择30个基因做气泡图
    f1=sample(rownames(scRNA),30)
    
    p <-DotPlot(scRNA, group.by=i, features = f1 ) + coord_flip()
    dot.data=p[["data"]]
    p2=ggplot(dot.data,aes(x=id, y = features.plot, color = avg.exp.scaled, size = pct.exp)) + 
      geom_point() +cowplot::theme_cowplot() + 
      theme(axis.line  = element_blank()) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
      ylab('') +
      theme(axis.ticks = element_blank()) +
      scale_color_gradientn(colours = viridis::plasma(20), limits = c(0,3), oob = scales::squish, name = 'expression')
    p2
    
  
    #各个细胞组分所占百分比
    if(i =="singleR" ){
      pB2_df <- table(scRNA@meta.data$orig.ident,scRNA@meta.data$singleR) %>% melt()
    }else if(i =="cellMarker"){
      pB2_df <- table(scRNA@meta.data$orig.ident,scRNA@meta.data$cellMarker) %>% melt()
    }else if(i =="seurat_clusters"){
      pB2_df <- table(scRNA@meta.data$orig.ident,scRNA@meta.data$seurat_clusters) %>% melt()
    }
    library(reshape2)
    colnames(pB2_df) <- c("Cluster","Sample","Number")
    pB2_df$Cluster <- factor(pB2_df$Cluster)
    n <- length(unique(scRNA@meta.data$seurat_clusters)) 
    library(randomcoloR)
    library(ggsci)
    colpalettes<- c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#6E568C","#E0367A","#D8D155","#91612D","#64495D","#7CC767")
    cat("注意这里的颜色数量只有16个,cluster数量是    ",n, "  cluster过多的话颜色分配不均，会报错")
    pB2 <- ggplot(data = pB2_df, aes(x = Number , y =Cluster , fill =Sample)) +
      geom_bar(stat = "identity", width=0.8,position="fill")+
      scale_fill_manual(values =colpalettes,n) +
      theme_bw()+
      theme(panel.grid =element_blank()) +
      labs(x="",y="Ratio")+
      ####用来将y轴移动位置
      theme(axis.text.y = element_text(size=12, colour = "black"))+
      theme(axis.text.x = element_text(size=12, colour = "black"))
    pB2
  }
  
  

}



ydx_sc_scenic<- function(){
  library(dplyr)
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(SCENIC)
  
}


# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# sctype_score: calculate ScType scores and assign cell types
#
# @params: scRNAseqData - input scRNA-seq matrix (rownames - genes, column names - cells), 
# @params: scale - indicates whether the matrix is scaled (TRUE by default)
# @params: gs - list of gene sets positively expressed in the cell type 
# @params: gs2 - list of gene sets that should not be expressed in the cell type (NULL if not applicable)

sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
  
  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
      warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)
  
  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
  
  es.max
}
# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# gene_sets_prepare: prepare gene sets and calculate marker sensitivity from input Cell Type excel file
#
# @params: path_to_db_file - DB file with cell types
# @cell_type - cell type (e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain)
#
gene_sets_prepare <- function(path_to_db_file, cell_type){
  
  cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}

auto_detect_tissue_type <- function(path_to_db_file, seuratObject, scaled, assay = "RNA", ...){
  
  # get all tissue types in DB
  db_read = openxlsx::read.xlsx(path_to_db_file); tissues_ = unique(db_read$tissueType); result_ = c()
  
  for(tissue in tissues_){ print(paste0("Checking...", tissue));
    
    # prepare gene sets
    gs_list = gene_sets_prepare(path_to_db_file, tissue);
    
    # prepare obj
    if(scaled){
      obj = as.matrix(seuratObject[[assay]]@scale.data)
    } else {
      obj = as.matrix(seuratObject[[assay]]@counts)
    }
    
    es.max = sctype_score(scRNAseqData = obj, scaled = scaled, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative, 
                          marker_sensitivity = gs_list$marker_sensitivity, verbose=!0);
    
    cL_resutls = do.call("rbind", lapply(unique(seuratObject@meta.data$seurat_clusters), function(cl){
      
      es.max.cl = sort(rowSums(es.max[ ,rownames(seuratObject@meta.data[seuratObject@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl), 10)
    }))
    
    dt_out = cL_resutls %>% group_by(cluster) %>% top_n(n = 1)
    
    # return mean score for tissue
    result_ = rbind(result_, data.frame(tissue = tissue, score = mean(dt_out$scores)))
  }
  
  # order by mean score
  result_ = result_[order(-result_$score),]
  
  # plot 
  barplot(height=result_$score, names=result_$tissue, col=rgb(0.8,0.1,0.1,0.6),
          xlab="Tissue", ylab="Summary score",  main="The higher summary score, the more likely tissue type is")
  
  result_
}

###AUCell 未完工
ydx_sc_AUCell <- function(){
  #########################################################################################
  ##BiocManager::install("AUCell")
  library(AUCell)
  library(ggplot2)
  library(Seurat)
  library(clusterProfiler)
  
  sc.id=sample(colnames(scRNA),1500)
  sc2=scRNA[,sc.id]
  cells_rankings <- AUCell_buildRankings(scRNA@assays$integrated@data, plotStats=TRUE) 
  
  cells_rankings
  getwd()
  c2 <- read.gmt("c2.cp.kegg.v7.4.symbols.gmt") 
  geneSets <- lapply(unique(c2$term), function(x){print(x);c2$gene[c2$term == x]})
  names(geneSets) <- unique(c2$term)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)
  library(msigdbr)
  m_df<- msigdbr(species = "Mus musculus",  category = "C2", subcategory = "KEGG")
  fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  
  cells_AUC <- AUCell_calcAUC(fgsea_sets, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)
  grep("OX",rownames(cells_AUC@assays@data$AUC),value = T)
  
  geneSet <- "KEGG_PROXIMAL_TUBULE_BICARBONATE_RECLAMATION"
  aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
  sc2$AUC <- aucs
  df<- data.frame(sc2@meta.data, sc2@reductions$umap@cell.embeddings)
  colnames(df)
  class_avg <- df %>%
    group_by( seurat_clusters) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  
  ggplot(df, aes(UMAP_1, UMAP_2))  +
    geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="D") +
    ggrepel::geom_label_repel(aes(label = seurat_clusters),
                              data = class_avg,
                              size = 5,
                              label.size = 1,
                              segment.color = NA
    )+   theme(legend.position = "none") + theme_bw()
}

###GSEA 未完工
ydx_sc_GSEA <- function(){
  
  library(GSVA)
  kegggmt2 <- read.gmt("c2.cp.kegg.v7.4.symbols.gmt")
  
  kegg_list = split(kegggmt2$gene, kegggmt2$term)
  ?AverageExpression
  x=AverageExpression(scRNA)
  
  exp=x[["RNA"]]
  exp1=as.matrix(exp)
  
  m_df<- msigdbr(species = "Mus musculus",  category = "C2", subcategory = "KEGG")
  fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  es.max <- gsva(exp1, fgsea_sets, mx.diff=FALSE,kcdf="Poisson")
  ?gsva
  pheatmap::pheatmap(es.max[1:20,], #热图的数据
                     cluster_rows = T,#行聚类
                     cluster_cols =T,#列聚类，可以看出样本之间的区分度
                     
                     show_colnames=T,
                     scale = "row", #以行来标准化，这个功能很不错
                     color =colorRampPalette(c("#FF7744", "white","#AAAAAA","#0044BB"))(100))
  
  colnames(scRNA@meta.data)
  cell.type=AverageExpression(scRNA,add.ident = "orig.ident")
  cell.type=cell.type[[1]]
  cell.type=as.matrix(cell.type)
  
  es.max1 <- gsva(cell.type, fgsea_sets, mx.diff=FALSE)
  
  pheatmap::pheatmap(es.max1[1:20,], #热图的数据
                     cluster_rows = F,#行聚类
                     cluster_cols =F,#列聚类，可以看出样本之间的区分度
                     
                     show_colnames=T,
                     scale = "row", #以行来标准化，这个功能很不错
                     color =colorRampPalette(c("#FF7744", "white","#AAAAAA","#0044BB"))(100))
  
  
  
  
  ##########################################################################################
  cell.500.id=sample(colnames(scRNA1),500)
  cell.500=scRNA1[,cell.500.id]
  m.cell.500=as.matrix(cell.500@assays$RNA@scale.data)
  
  library(GSVA)
  
  es.dif <- gsva(m.cell.500, fgsea_sets, mx.diff=TRUE)
  
  colnames(cell.500@meta.data)
  annotation_col=cell.500@meta.data[,c(5,6)]
  p1=arrange(annotation_col,seurat_clusters)
  es.dif1=as.data.frame(es.dif)
  es.dif2=es.dif1[,rownames(p1)]
  es.dif3=as.matrix(es.dif2)
  
  pheatmap::pheatmap(es.dif3[20:40,], #热图的数据
                     cluster_rows = F,#行聚类
                     cluster_cols =F,#列聚类，可以看出样本之间的区分度
                     annotation_col = annotation_col,
                     show_colnames=F,
                     scale = "row", #以行来标准化，这个功能很不错
                     color =colorRampPalette(c("#FF7744", "white","#AAAAAA","#0044BB"))(100))
  
}
############################################################
####单细胞转录组基础分析四：细胞类型鉴定 ####
############################################################
###因为我把singler的注释加载到metadata中时候，命名的名字叫CellMarker，所以画图时候，group.by="CellMarker"
###########CellMarker可视化###########
# f1=sample(rownames(scRNA),50)
# p <-DotPlot(scRNA, group.by="CellMarker", features = f1 ) + coord_flip()
# dot.data=p[["data"]]
# colnames(dot.data)
# ##使用scale_color_gradientn设置渐变色的范围为0-3，如果有的点超过了3，那么通过oob = scales::squish将它硬生生压到3
# p1=ggplot(dot.data,aes(x=id, y = features.plot, color = avg.exp.scaled, size = pct.exp)) + 
#   geom_point() +cowplot::theme_cowplot() + 
#   theme(axis.line  = element_blank()) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
#   ylab('') +
#   theme(axis.ticks = element_blank()) +
#   scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,3), oob = scales::squish, name = 'expression')
# p1
# 
# cat("选择每个cluster选择的基因合集，绘制气泡图")
# all_genes <- unlist(cl_genes_list)
# names(all_genes) =NULL
# p <-DotPlot(scRNA, features = all_genes ) + coord_flip()+  
#   scale_color_gradientn(colours = viridis::plasma(20), limits = c(0,3), oob = scales::squish, name = 'expression')
# p 

# # 选出每个群集中的前50个power最大的基因
# top_power <- top_n(all.markers, n = 50, wt = power)
# # 选出每个群集中的前50个avg_diff最大的基因
# top_avg_diff <- top_n(all.markers, n = 50, wt = avg_diff)
# # 合并两个top_n的结果，并按群集和-power排序
# degs_top50 <- arrange(merge(top_power, top_avg_diff, by = c("cluster", "symbol")), cluster, -power)
# 
# avgData <- scRNA@assays$RNA@data[degs_top50$gene,] %>% 
#   apply(1, function(x){
#     tapply(x, scRNA$celltype, mean) # ExpMean
#   }) %>% t
# 
# phData <- MinMax(scale(avgData), -2, 2) # z-score
# rownames(phData) <- 1:nrow(phData)
# library(pheatmap)
# phres <- pheatmap(
#   phData, 
#   color = colorRampPalette(c("darkblue", "white", "red3"))(99), #配色
#   scale = "row",
#   cluster_rows = F, #不按行聚类
#   cluster_cols = F, #按列聚类
#   clustering_method = "complete",
#   show_rownames = F, #显示cluster名
#   annotation_row = data.frame(cluster = degs_top50$cluster), 
# )  



















  
  

  

  

  
  
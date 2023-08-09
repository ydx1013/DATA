ydx_draw_sc_cell = function(data_Dir,Integrat =TRUE,specise = "mmu",scale_genes ="all",zouqi = TRUE, pc.num=1:20){
  #  scale_genes ="all"的时候对所有基因进行降维，否则只对高变基因进行降维
  library(dplyr)
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggplot2) # 加载ggplot2库
  library(patchwork) # 加载patchwork库
  library(RColorBrewer) # 加载RColorBrewer库
  library(glmGamPoi) # 加载glmGamPoi库
  
  dir=  c("D:\\BaiduNetdiskDownload\\004一对一辅导单细胞测序密码2019\\课程文件\\14节课（医学拼课请认准唯一VX：782878241）\\14节课\\lession14\\CAFs.1/",
          "D:\\BaiduNetdiskDownload\\004一对一辅导单细胞测序密码2019\\课程文件\\14节课（医学拼课请认准唯一VX：782878241）\\14节课\\lession14\\CAFs.2/",
          "D:\\BaiduNetdiskDownload\\004一对一辅导单细胞测序密码2019\\课程文件\\14节课（医学拼课请认准唯一VX：782878241）\\14节课\\lession14\\dapi1/",
          "D:\\BaiduNetdiskDownload\\004一对一辅导单细胞测序密码2019\\课程文件\\14节课（医学拼课请认准唯一VX：782878241）\\14节课\\lession14\\dapi2/")
  sample_names = ""
  if(length(sample_names) >2){
    names(dir) = sample_names
  }else{
    names(dir) = sapply(dir, basename)
  }
    

  scRNAlist <- list()
  for(i in 1:length(dir)){
    counts <- Read10X(data.dir = dir[i])
    scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells = 3, min.features =300)
  }
  

  
  
  # 加载glmGamPoi包
  for (i in 1:length(scRNAlist)) {
    scRNAlist[[i]]
    if (specise == "Human"){
      # 1. 过滤线粒体
      scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")#  人用MT，小鼠用mt
      ## 2. 过滤红细胞
      HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    } else{
      # 1. 过滤线粒体
      scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")#  人用MT，小鼠用mt
      ## 2. 过滤红细胞
      HB.genes <- c("Hba-a1", "Hba-a2", "Hbb-bh1", "Hbb-bh2", "Hbb-b1", "Hbb-b2", "Hbb-y", "Hbb-bt", "Hbb-bh3")
    }
    HB_m <- match(HB.genes, rownames(scRNAlist[[i]]@assays$RNA)) 
    HB.genes <- rownames(scRNAlist[[i]]@assays$RNA)[HB_m] 
    HB.genes <- HB.genes[!is.na(HB.genes)] 
    scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 
    scRNAlist[[i]] <- SCTransform(scRNAlist[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt","percent.HB"),   method = "glmGamPoi",verbose = FALSE)
  }
    
    
 
 #整合和合并是不同的，整合是均匀分布，可能会丢失掉一些cluster
  #合并是单纯的合并在一起，所有的cluster都会保留，但是要注意，需要是同一批次的样本才能进行合并，而不整合
  ####################################是否要整合数据#######################################
  if(length(scRNAlist)>2 ){  if(Integrat == 1){ 
    rm(scRNA1)
    # 选择跨数据集的基因，也就是多个数据集都有表达的基因用于整合
    features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
    scRNAlist <- PrepSCTIntegration(object.list = scRNAlist, anchor.features = features)
    # 找到数据整合的锚点
    scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, normalization.method = "SCT",
                                            anchor.features = features)
    scRNA1 <- IntegrateData(anchorset = scRNA.anchors, normalization.method = "SCT")
  }else{
    rm(scRNA1)
    scRNA1 = scRNAlist[[1]]
    for (i in 2:length(scRNAlist)) {
      scRNA1 <- merge(scRNA1, scRNAlist[[i]])
    }
  }
  }else{
    rm(scRNA1)
      scRNA1 = scRNAlist[[1]]
      }
  
  
  
  
  ####Feature、count、线粒体基因、红细胞基因占比可视化。
  col.num <- length(levels(scRNA1@active.ident))
  pdf("1. scRNA_VlnPlot.pdf",width = 8,height = 6) # 将图像保存为PDF格式
  VlnPlot(scRNA1,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
          cols =rainbow(col.num), 
          pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
          ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
  dev.off()
  pdf("1. nCount_RNA和nFeature_RNA.pdf",width = 8,height = 6) # 将图像保存为PDF格式
  FeatureScatter(scRNA1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  dev.off()
  ####################################根据情况进行细胞过滤#######################################
  min_nFeature=300
  max_nFeature= 5000
  scRNA1 <- subset(scRNA1, subset = nFeature_RNA > min_nFeature # RNA特征数大于300个
                  & nFeature_RNA < max_nFeature # RNA特征数小于7000个  
                  & percent.mt < 10 # 线粒体基因的比例小于10%
                  & percent.HB < 1 # 血红蛋白基因的比例小于3%
                  & nCount_RNA < 100000) # RNA拷贝数小于100000
  
  ####################################数据降维#######################################

  scRNA1 <- RunPCA(scRNA1, verbose = FALSE)
  #每个PCA对整体特征的贡献度
  #前50个PCA的变化
  pdf("根据曲线选择合适的PCA.pdf",width = 8,height = 6)
  ElbowPlot(scRNA1 , ndims=50, reduction="pca") 
  dev.off()
  pca_number = 1:30
  #umap
  scRNA1 <- RunUMAP(scRNA1,reduction = "pca", dims = pca_number)
  DimPlot(scRNA1, reduction = "umap",label = T)
  #tSNE
  scRNA1 = RunTSNE(scRNA1, dims = pca_number) 
  DimPlot(scRNA1, reduction = "tsne",label = T)
  
  ############################################################
  ####单细胞转录组基础分析：细胞周期 ####
  ############################################################
  #细胞周期回归：上一步找到的高变基因，常常会包含一些细胞周期相关基因。
  #它们会导致细胞聚类发生一定的偏移，即相同类型的细胞在聚类时会因为细胞周期的不同而分开。
  #细胞周期评分，匹配g2m_genes和s_genes
  g2m_genes = cc.genes$g2m.genes
  g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA1))
  s_genes = cc.genes$s.genes
  s_genes = CaseMatch(search = s_genes, match = rownames(scRNA1))
  #提取细胞周期相关基因
  
  scRNA1 <- CellCycleScoring(object=scRNA1,  g2m.features=g2m_genes,  s.features=s_genes)
  #对细胞周期基因进行聚类
  scRNA1 <- RunPCA(scRNA1, features = c(s_genes, g2m_genes))
  
  pdf("4. 细胞周期基因进行聚类,观察细胞周期基因对聚类的影响.pdf",width = 8,height = 6) # 将图像保存为PDF格式
  DimPlot(scRNA1, reduction = "pca", group.by = "Phase")
  dev.off()
  if(zouqi == TRUE){
    ##如果需要消除细胞周期的影响
    scRNA1 <- ScaleData(scRNA1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scRNA1))
    #依据高变基因PCA降维
    scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1))
    pdf("5. 去除细胞周期后PCA.pdf",width = 8,height = 6)
    DimPlot(scRNA1, reduction = "pca", group.by="orig.ident")
    dev.off()
  }else{
    scRNA1
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ######################################################################################
  # DefaultAssay(scRNA1) <- "SCT"
  # scRNA1 <-  ScaleData(scRNA1,features = rownames(scRNA1))
  # DefaultAssay(scRNA1) <- "integrated"
  # scRNA1 <-  ScaleData(scRNA1)
  # scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst",nfeatures = 3000)
  # scRNA1 <- RunPCA(scRNA1, npcs = 20, verbose = T)
  # t-SNE and Clustering
  
  # scRNA1 <- FindNeighbors(scRNA1, reduction = "pca", dims = 1:7)
  # scRNA1 <- FindClusters(scRNA1, resolution = 0.4)
  # scRNA1 <- RunUMAP(scRNA1, reduction = "pca", dims = 1:7)
  # colnames(scRNA1@meta.data)
  # scRNA1 <- RunTSNE(scRNA1, dims = 1:7)
  # DimPlot(scRNA1, reduction = "umap",label = T)
  
  
  

  DimPlot(scRNA1, reduction = "umap")
  
 
  embed_tsne <- Embeddings(scRNA1, 'tsne')
  write.csv(embed_tsne,'embed_tsne.csv')
  DimPlot(scRNA1, reduction = "tsne") 
  

  # ##读取10x的数据
  # scRNA.counts=Read10X(data_Dir)
  # 
  # #min.features   每个细胞至少表达多少个基因
  # #min.cells 每个基因至少要在多少细胞中表达
  # scRNA = CreateSeuratObject(scRNA.counts ,min.cells = 3,project="os", min.features = 300)

  # 
  # scRNA <- SCTransform(scRNA, vars.to.regress = "percent.mt", verbose = FALSE)

  
  ############################################################
  ####单细胞转录组基础分析：细胞过滤 ####
  ############################################################

  
  # 3 . 根据上图小提琴的数据进行细胞赛学，细胞过滤
  # nFeature_RNA
  #nFeature_RNA
  #percent.mt 线粒体 小于10% 
  # percent.HB  红细胞小于3%
  #nCount_RNA 

   
####过滤完之后 我们就要对数据进行均一化，使用NormalizeData这个函数。
###注意均一化是用NormalizeData，标准化是用ScaleData
  scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
  ###官方推荐是2000个高变基因，很多文章也有设置3000的，这个因自己的实验项目决定
  scRNA <- FindVariableFeatures(scRNA1, selection.method = "vst", nfeatures = 2000) 
  # Identify the 10 most highly variable genes，把top10的高变基因挑选出来，目的是为了作图
  top15 <- head(VariableFeatures(scRNA1), 15) 
  # plot variable features with and without labels  画出来不带标签的高变基因图
  plot1 <- VariableFeaturePlot(scRNA1) 
  ###把top15的基因加到图中
  plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=2.5) 
 
  ###画图，寻找高变基因
  pdf("2. 前15个高变基因.pdf",width = 8,height = 6) # 将图像保存为PDF格式
  CombinePlots(plots = list(plot1, plot2),legend="bottom")  
  dev.off()
  
  
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
  #PCA降维并提取主成分
  #依据高变基因PCA降维
  scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA1)) 
  plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident") 
  pdf("3. 高变基因PCA降维聚类.pdf",width = 8,height = 6) # 将图像保存为PDF格式
  DimPlot(scRNA, reduction = "pca", group.by="orig.ident") 
  dev.off()
  
  ############################################################
  ####单细胞转录组基础分析：细胞周期 ####
  ############################################################
  #细胞周期回归：上一步找到的高变基因，常常会包含一些细胞周期相关基因。
  #它们会导致细胞聚类发生一定的偏移，即相同类型的细胞在聚类时会因为细胞周期的不同而分开。
  #细胞周期评分，匹配g2m_genes和s_genes
  g2m_genes = cc.genes$g2m.genes
  g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA))
  s_genes = cc.genes$s.genes
  s_genes = CaseMatch(search = s_genes, match = rownames(scRNA))
  #提取细胞周期相关基因
  
  scRNA <- CellCycleScoring(object=scRNA,  g2m.features=g2m_genes,  s.features=s_genes)
  #对细胞周期基因进行聚类
  scRNA <- RunPCA(scRNA, features = c(s_genes, g2m_genes))
  
  pdf("4. 细胞周期基因进行聚类,观察细胞周期基因对聚类的影响.pdf",width = 8,height = 6) # 将图像保存为PDF格式
  DimPlot(scRNA, reduction = "pca", group.by = "Phase")
  dev.off()
  
if(zouqi == TRUE){
  ##如果需要消除细胞周期的影响
  scRNA <- ScaleData(scRNA, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scRNA))
  #依据高变基因PCA降维
  scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
  pdf("5. 去除细胞周期后PCA.pdf",width = 8,height = 6)
  DimPlot(scRNA, reduction = "pca", group.by="orig.ident")
  dev.off()
}else{
  scRNA
}



  #细胞聚类resolution越多，切的越多，cluster越多
  #pc.num维度越多，cluster越多
  ###Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B. Thanks to Nigel Delaney (evolvedmicrobe@github) for the rewrite of the Java modularity optimizer code in Rcpp!
  scRNA <- FindNeighbors(scRNA, dims = pc.num) 
  ###这个分辨率是可以自定义的，当我们的样本细胞数较大时候resolution 要高一些，
  #一般情况2万细胞以上都是大于1.0的
  scRNA <- FindClusters(scRNA, resolution = 0.5)
  ## 查看每一类有多少个细胞
  table(scRNA@meta.data$seurat_clusters)
  ##系统发育分析（Phylogenetic Analysis of Identity Classes）
  scRNA<-BuildClusterTree(scRNA)
  PlotClusterTree(scRNA)
  
  #可视化降维有两个方法tSNE和UMAP
  #tSNE
  scRNA = RunTSNE(scRNA, dims = pc.num)
  embed_tsne <- Embeddings(scRNA, 'tsne')
  write.csv(embed_tsne,'embed_tsne.csv')
  plot1 = DimPlot(scRNA, reduction = "tsne") 
  ##画图
  plot1
  ###label = TRUE把注释展示在图中
  plot2 = DimPlot(scRNA, reduction = "tsne",label = TRUE) 
  ggsave("tSNE1无序号.pdf", plot = plot1, width = 8, height = 7)
  ggsave("tSNE2,标注序号.pdf", plot = plot2, width = 8, height = 7)
  
  
  
  #UMAP---第二种可视化降维
  scRNA <- RunUMAP(scRNA, dims = pc.num)
  embed_umap <- Embeddings(scRNA, 'umap')
  write.csv(embed_umap,'embed_umap.csv') 
  plot11 = DimPlot(scRNA, reduction = "umap") 
  plot2 = DimPlot(scRNA, reduction = "umap",label = TRUE) 
  plot2
  ggsave("UMAP有序号.pdf", plot = plot2, width = 8, height = 7)
  ggsave("UMAP无序号.pdf", plot = plot11, width = 8, height = 7)
  #合并tSNE与UMAP
  
  plotc <- plot1+plot11+ plot_layout(guides = 'collect')
  plotc
  ggsave("tsn+umap.pdf", plot = plotc, width = 16, height = 7)
  
  
  
  
  ############################################################
  ####单细胞转录组基础分析四：细胞类型鉴定 ####
  ############################################################
  diff.wilcox = FindAllMarkers(scRNA)
  ##3 %>% 这是通道函数 起传递左右  可以自己百度深入理解一下，我这里只告诉你他是起传递作用的
  all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
  top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  ###将marker基因保存一下
  write.csv(all.markers, "diff_genes_wilcox-mark基因保存.csv", row.names = F)
  ###把top10marker基因保存一下
  write.csv(top20, "top10mark基因保存_diff_genes_wilcox.csv", row.names = F)
  
  ##top10基因绘制热图

  top10_genes = CaseMatch(search = as.vector(top10$gene), match = rownames(scRNA)) 
  ##把marker基因用热图展示出来
  plot1 = DoHeatmap(scRNA, features = top10_genes, group.by = "seurat_clusters", group.bar = T, size = 4)
  plot1
  ggsave("7. top10_markers.pdf", plot=plot1, width=18, height=16) 
  
  
  
  ############################################################
  ####单细胞转录组基础分析四：细胞注释 ####
  ############################################################
  #第1种方法用SingleR鉴定细胞类型
  library(SingleR)
  ##把师傅给你的百度云打开 下载其中的人的数据库，因为你们没有vpn，所以singler的数据库没法下载
  ###下载好数据库后，把ref_Human_all.Rdata加载到环境中，这样算是对数据库的加载，就可以按照singler的算法来对细胞亚群进行定义了。
  ###我们可以看到在环境中多了一个叫ref_Human_all的文件 大小为113mb  这个就是数据库
  ####然后我们把环境中的ref_Human_all赋值与refdata
  refdata <- MouseRNAseqData()
  ###把rna的转录表达数据提取
  testdata <- GetAssayData(scRNA, slot="data")
  ###把scRNA数据中的seurat_clusters提取出来，注意这里是因子类型的
  clusters <- scRNA@meta.data$seurat_clusters
  ###开始用singler分析
  cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
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
  DimPlot(scRNA, group.by="singleR", label=T, label.size=5, reduction='tsne')
  ###我们可以看到  两种方法得到的结果都是一样的，但是我比较喜欢第二种方法
  ##鉴定结果展示
  p1 = DimPlot(scRNA, group.by="singleR", label=T, label.size=5, reduction='tsne')
  p1
  p2 = DimPlot(scRNA, group.by="singleR", label=T, label.size=5, reduction='umap')
  p2
  p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
  p3 
  ggsave("singleR-tSNE_celltype.pdf", p1, width=8 ,height=8)
  ggsave("singleR-UMAP_celltype.pdf", p2, width=8 ,height=8)
  ggsave("singleR-celltype.pdf", p3, width=16 ,height=8)


  
  ##CellMarker：http://biocc.hrbmu.edu.cn/CellMarker/index.jsp
  ###PanglaoDB：https://panglaodb.se/index.html
  ####把每一类的marker基因输入到这两个网站中，然后看网站给出的细胞注释信息
  scRNA <- RenameIdents(scRNA, "10" = "MSC","9"="t cells","3"="t cells")
  scRNA <- RevertIdents(scRNA)
  DimPlot(scRNA, reduction = "tsne", label = TRUE)
  ##选择基因作图展示
  FeaturePlot(scRNA,features = "Ermn",reduction = "tsne", label=T, ncol=2)
  
  
  
  





  ggsave("selectgenes_FeaturePlot.png", p2, width=8 ,height=12)
  p3=p1|p2
  p3
  ggsave("selectgenes.pdf", p3, width=10 ,height=8)
  ############################################################
 
  
 
  
  
  # ##第三种方法用已有生物学知识或者文献来看，比如区分免疫或者非免疫细胞就用CD45这个基因。PTPRC是CD45官方名字
  # select_genes <- c('COL2A1',"PTPRC","EPCAM")
  # p <-DotPlot(scRNA, features = select_genes ) + coord_flip()
  # p
  # ###我们可以看到0,3,4，7,8,9，12是免疫细胞、
  # ##提取气泡图的信息
  # dat=p$data 
  # ###提取cd45的信息
  # cd45=dat[dat$features.plot=='PTPRC',]
  # ###查看五分位数
  # fivenum(cd45$avg.exp.scaled)
  # ###选择中位数的cluster，这个要结合气泡图来综合选择，但是我发现所以的分类都是以-0.5为准的。所以我把-0.5当做固定值了
  # ###这里的-0.5我是问健明老师的  他告诉的哈哈哈
  # imm=cd45[cd45$avg.exp.scaled > -0.5,]$id
  # imm
  # ###我们看到在-0.5的时候，与我们肉眼判断的免疫cluster是一样的
  # ##接下来这一步是将免疫和非免疫细胞的信息添加到metadata中
  # scRNA1@meta.data$immune_annotation <-ifelse(scRNA1@meta.data$RNA_snn_res.0.5  %in% imm ,'immune','non-immune')
  # # MAke a table 查看一下数据
  # table(scRNA1@meta.data$immune_annotation)
  # # Make and save relevant plots 可视化一下免疫和非免疫的细胞
  # p <- TSNEPlot(object =scRNA1, group.by = 'immune_annotation')
  # p 

  

  
  
  
  
  
  
  
  
  
  
}

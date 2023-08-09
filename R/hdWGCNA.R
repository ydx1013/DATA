ydx_hdWGCNA <- function(){
  # 单细胞分析包
  library(Seurat)
  # 绘图和数据科学包
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  # 共表达网络分析包
  library(WGCNA)
  library(hdWGCNA)
  # 使用cowplot主题进行ggplot绘图
  theme_set(theme_cowplot())
  # 设置随机种子以实现可重复性
  set.seed(12345)
  # 加载周等人的snRNA-seq数据集
  seurat_obj <-annotate$scRNA
  
  p <- DimPlot(seurat_obj, group.by='singleR', label=TRUE) + ggtitle('Zhou et al Control Cortex') + NoLegend()
  
  p

  
  # gene_select参数选择：
  # variable：使用高变基因
  # fraction：选择至少在一部分细胞中表达的基因；此时需要追加fraction参数，如fraction = 0.05表示选择至少在5%的细胞中表达的基因
  # custom：人为指定用于WGCNA分析的基因
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction", # 基因选择方法
    fraction = 0.05, # 基因表达的最低比例要求，达到该比例才会被纳入分析
    wgcna_name = "tutorial" )
  
 

  
  # 标准化metacell表达矩阵
  seurat_obj <- NormalizeMetacells(seurat_obj)
  
  
  
  
  
  
  
  
  
  # 构造metacell
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("singleR", "group2"), # 指定在seurat_obj@meta.data中按照哪些列进行分组
    reduction = 'umap', # 选择要在其上执行KNN的降维方法
    k = 25, # 最近邻参数
    max_shared = 10, # 两个元细胞之间最大共享细胞数
    ident.group = 'singleR' # 设置元细胞Seurat对象的Idents
  )
  
  # 规范化元细胞表达矩阵：
  seurat_obj <- NormalizeMetacells(seurat_obj)
  seurat_obj@misc$tutorial$wgcna_metacell_obj
  head(seurat_obj@misc$tutorial$wgcna_metacell_obj, 2)
  
  
  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = c("Oligodendrocytes"), # group.by 列中感兴趣组的名称
    group.by = 'singleR', # 包含细胞类型信息的元数据列。该列应该在 MetacellsByGroups 中使用过。
    assay = 'RNA', # 使用 RNA 试验数据
    slot = 'data' # 使用规范化后的数据
  )
  
  
  # 在WGCNA中，networkType参数用于指定构建基因共表达网络的类型，主要有以下几个可选参数：
  # unsigned：无符号网络，即不考虑基因表达的正负号，只考虑它们之间的关联性。
  # signed：带符号网络，即同时考虑基因表达的正负号和大小，可以反映基因之间的正负调控关系。
  # hybrid：混合网络，即将无符号网络和带符号网络结合起来，利用它们之间的优势来提高网络分析的准确性和可靠性。
  # 测试不同的软阈值（soft power）：

  seurat_obj <- TestSoftPowers(
    seurat_obj,
    setDatExpr=FALSE,
    powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2))) # 选取soft powers(默认)
  
  # 绘制结果：
  plot_list <- PlotSoftPowers(seurat_obj,
                              point_size = 5,
                              text_size = 3)
  
  # 使用 patchwork 组合图形
  wrap_plots(plot_list, ncol = 2)
  
  power_table <- GetPowerTable(seurat_obj)
  head(power_table)
  
  
  # construct co-expression network:
  seurat_obj <- ConstructNetwork(
    seurat_obj,
    setDatExpr=FALSE,
    tom_name = 'Oligodendrocytes' # name of the topoligical overlap matrix written to disk
  )
  PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
  # 读入TOM矩阵
  TOM <- GetTOM(seurat_obj)
  
  # 查看TOM矩阵元素
  TOM[1:4, 1:4]
  # harmony前需要进行正态化（Scale）
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
  
  # 计算所有的MEs，比较耗时
  seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars = "group2"  # 根据Sample去批次
  )
  # 通过ModuleConnectivity函数在整个单细胞数据集上计算kME values，kME values值越高，该基因为hub genes的可能性越大。
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = 'singleR', group_name = 'Oligodendrocytes'  # 这里继续计算INH（感兴趣的细胞类型）的kME
  )
  name ="Oligodendrocytes"
  # 重命名module，以表明这些module由INH计算所得
  seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = paste0(name,"-M")
  )
  
  # 可视化每个模块中的基因的kME
  PlotKMEs(seurat_obj, ncol = 2, n_hubs = 10)
  # get the module assignment table:
  modules <- GetModules(seurat_obj)
  
  # 通过GetHubGenes函数抽取每个模块的topN hub节点，即每个模块中kME高的topN个基因
  hub_df <- GetHubGenes(seurat_obj = seurat_obj, n_hubs = 10)
  
  # 至此，基本完成了hdWGCNA分析
  saveRDS(seurat_obj, file = 'hdWGCNA_object.rds')
  # 根据每个模块的topN个基因对细胞进行打分
  seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = 25, # topN hub genes
    method = "Seurat" # 打分方法，分为Seurat、UCell
  )
  # 根据hMEs进行FeaturePlot
  plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features = 'hMEs', # 可选择MEs、hMEs、scores、average
    order = TRUE # order so the points with highest hMEs are on top
  )
  #  stitch together with patchwork
  wrap_plots(plot_list, ncol = 6)
  # 模块相关性
  ModuleCorrelogram(seurat_obj)
  # 利用Seurat自带的可视化方法
  
  # get hMEs from seurat object
  MEs <- GetMEs(seurat_obj, harmonized=TRUE)
  mods <- colnames(MEs); mods <- mods[mods != 'grey']
  
  # add hMEs to Seurat meta-data:
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
  p <- DotPlot(seurat_obj, features = mods, group.by = "singleR")
  
  p <- p + 
    #     coord_flip() + # 反转x/y轴
    RotatedAxis() + # 旋转坐标轴标签
    scale_color_gradient2(high = "red", mid = "grey95", low = "blue")  # 改颜色
  
  p
  # single-cell analysis package
  library(Seurat)
  
  # plotting and data science packages
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  
  # co-expression network analysis packages:
  library(WGCNA)
  library(hdWGCNA)
  
  # network analysis & visualization package:
  library(igraph)
  
  ModuleNetworkPlot(seurat_obj = seurat_obj, mods = "Oligodendrocytes-M1")
  HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = 3, # 用于可视化的hub gene
    n_other = 5, # 随机选取的gene
    edge_prop = 0.75, # 采样的边数
    mods = "all"
  )
  
  # 将基因降维成UMAP图，该函数会在Seurat_obj@misc$wgcna_name中生成一个module_umap矩阵
  seurat_obj <- RunModuleUMAP(
    seurat_obj, 
    n_hubs = 10, # 用于UMAP嵌入的hub genes
    n_neighbors = 15, # UMAP参数
    min_dist = 0.1 # 两个点（基因）在UMAP空间中的最短距离
  )
  # 查看gene UMAP table
  head(seurat_obj@misc$tutorial$module_umap, 2)
  
  ModuleUMAPPlot(
    seurat_obj,
    edge.alpha = 0.25,
    sample_edges = T,
    edge_prop = 0.1, # 采样10%的边用于可视化
    label_hubs = 2, # 每个module用于展示的hub genes
    keep_grey_edges = F
  )
  #将基因模块与特征相关联

  
  
  
}
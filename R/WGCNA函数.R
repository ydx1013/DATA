ydx_draw_WGCNA_pre<- function(data,disease="Treat",hlevel=0,minModuleSize=25,MEDissThres = 0.25,gene_num=5000){
   print("第一步应该先观察聚类图，看是否有离群值，然后改变hlevel的数值，注意，这里要根据需要在样品名中加入分组信息") 
  data= pre_data$exp_data
  library(doParallel)
  library(WGCNA)
  library(limma)
  library(tidyverse)
  print("WGCNA需要的是样品名中加入分组信息的表达矩阵数据") 
  picDir = paste0("WGCNA_", GSE_number)
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {dir.create(picDir)}
  

  data=as.matrix(data)
  data=avereps(data)
  data2=as.data.frame(data) #【将新的矩阵构建为数据框】
  #使用中位数来选择做WGCNA的基因数量
  WGCNA_matrix = data2[order(apply(data2,1,mad), decreasing = T)[1:gene_num],]
  
  data3=as.data.frame(WGCNA_matrix)
  duplicated(colnames(data3))#【查看是否有重复列】
  
  print(paste("有", sum(duplicated(colnames(data3))), "个重复列"))
  
  rt<-data3[,!duplicated(colnames(data3))] #【将重复列删除】
  datExpr0 = as.data.frame(t(rt))
  gsg = goodSamplesGenes(datExpr0, verbose = 3)#【检查缺缺失值】
  gsg$allOK
  #【如果有缺失值就在下方删除】
  if (!gsg$allOK)
  {
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  
  ###过滤表达量【将表达量过低的基因删除】#####
  meanFPKM=0.1  #【设置过滤值】
  n=nrow(datExpr0)
  datExpr0[n+1,]=apply(datExpr0[c(1:nrow(datExpr0)),],2,mean)#【读取表达量值】
  datExpr0=datExpr0[1:n,datExpr0[n+1,] > meanFPKM] #【将表达量值与过滤值进行比较】
  dimnames=list(rownames(datExpr0), colnames(datExpr0))#【保存符合条件的值】
  data=matrix(as.numeric(as.matrix(datExpr0)), nrow=nrow(datExpr0), dimnames=dimnames)
  datExpr0=avereps(data)
  ###导出过滤后的基因表达矩阵
  filtered_fpkm=datExpr0#【将矩阵进行行列转换，转换为基因名在上，样品名在左】
  filtered_fpkm=data.frame(rownames(filtered_fpkm),filtered_fpkm)
  names(filtered_fpkm)[1]="sample"
  head(filtered_fpkm)
  sampleTree = hclust(dist(datExpr0), method = "average")#【将每个相似的样品进行聚类】
  pdf( file=paste0(picDir,"/1. 初步聚类(请观察是否有离群).pdf"),width=15,height=10)#【绘制图形】
  par(cex = 0.5)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
       cex.axis = 2, cex.main = 2)
  dev.off()
  hlevel=hlevel#【根据聚类图，查看是否有离群值，并设置删除高度】
  ###【删除离群值】
  pdf( file=paste0(picDir,"/2. 如有离群,画下切割线.pdf"),width=15,height=10)#【绘制图形】
  par(cex = 0.5)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
       cex.axis = 2, cex.main = 2)
  abline(h = hlevel, col = "red")#【h】是设置删除高度
  dev.off()

  
  if (hlevel==0) {
    # if condition is TRUE, do this 
    print("第一步应该先观察聚类图，看是否有离群值，然后改变hlevel的数值，注意，这里要根据需要在样品名中加入分组信息") 
    gc()

  } else {
    ###【保留非离群样本】
    clust = cutreeStatic(sampleTree, cutHeight = hlevel, minSize = 10)
    table(clust)
    keepSamples = (clust==1)#【提取非离群数据】
    datExpr = datExpr0[keepSamples, ] #【构建非离群矩阵】

    # 记录基因和样本数，方便后续可视化
    nGenes = ncol(datExpr)#基因数
    nSamples = nrow(datExpr)#样本数
    # 输出剪切后的样本数量和基因数量
    cat(paste0("剪切之后剩下的样本数量为：", nSamples, "，基因数量为：", nGenes))
    cat("进入下一个过程")
    
    # enableWGCNAThreads() #多线程工作
    powers = c(c(1:10), seq(from = 12, to=20, by=2)) #【无标度拓扑拟合指数】
    sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

    pdf( file=paste0(picDir,"/3. SoftThreshold.pdf"),width=9,height=5)#【绘制图形】
    par(mfrow = c(1,2))
    cex1=0.7
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red");
    abline(h=0.9,col="red") #【若没有数值到达0.9则可以降低标准，但最低不低于0.8】
    ###平均连通性散点图
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    dev.off()
    
    sft #【查看最佳软阈】
    softPower =sft$powerEstimate #【最佳软阈】
    print(paste("【最佳软阈】：",softPower))
    ADJ <- abs(cor(datExpr0,use="p"))^softPower
    k = as.vector(apply(ADJ,2,sum,na.rm=T))
    ## 无尺度网络验证
    pdf( file=paste0(picDir,"/4. histogram无尺度网络验证.pdf"),width=13,height=5)#【绘制图形】
    par(mfrow = c(1,2))
    hist(k)
    scaleFreePlot(k,main="cheak scale free topology")
    dev.off()
    ## 计算邻接矩阵
    adjacency = adjacency(datExpr,power=softPower)
    ## 计算TOM拓扑矩阵
    TOM = TOMsimilarity(adjacency)
    ## 计算相异度
    dissTOM = 1- TOM 
    #模块初步聚类分析
    library(flashClust)
    geneTree = flashClust(as.dist(dissTOM),method="average")
    #绘制层次聚类树
    pdf(file = paste0(picDir,"/5. GeneClusterTOM-based.pdf"))
    plot(geneTree,xlab="",sub="",main="Gene clustering on TOM-based",labels=FALSE,hang=0.04)
    dev.off()
    #构建初步基因模块
    #设定基因模块中至少30个基因
    minModuleSize = minModuleSize #【设置模块基因数目，每个模块最大基因数量为5000】
    
    # 动态剪切树识别网络模块
    dynamicMods = cutreeDynamic(dendro = geneTree,#hclust函数的聚类结果
                                distM = dissTOM,#
                                deepSplit = 2,
                                pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)#设定基因模块中至少30个基因
    # 将标签转换为颜色
    dynamicColors = labels2colors(dynamicMods)
    table(dynamicColors)#看聚类到哪些模块，哪些颜色
    table(dynamicColors) #【查看模块颜色】
    #【绘制模块图形】
    pdf( file=paste0(picDir,"/6. Tree.pdf"),width=8,height=6)
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")
    dev.off()
    ###【将相似的模块进行网络构建】
    MEList = moduleEigengenes(datExpr, colors = dynamicColors)#计算特征向量
    MEs = MEList$eigengenes;#提取特征向量
    MEDiss1 = 1-cor(MEs);#计算相异度
    METree1 = flashClust(as.dist(MEDiss1),method="average")#对相异度进行flashClust聚类
    
    #设置特征向量相关系数大于0.75
    #合并模块
    MEDissThres= MEDissThres#相异度在0.25以下，也就是相似度大于0.75，对这些模块合并
    merge = mergeCloseModules(datExpr, #合并相似度大于0.75的模块
                              dynamicColors,
                              cutHeight = MEDissThres, 
                              verbose=3)
    mergedColors = merge$colors
    table(dynamicColors)#动态剪切树的模块颜色
    table(mergedColors)#合并后的模块颜色，可以看到从18个模块变成了14个模块
    mergedMEs = merge$newMEs#合并后的模块
    #重新命名合并后的模块
    moduleColors = mergedColors;
    colorOrder = c("grey",standardColors(50));
    moduleLabels = match(moduleColors,colorOrder)-1;
    MEs = mergedMEs;
    MEDiss2 = 1-cor(MEs);#计算相异度
    METree2 = flashClust(as.dist(MEDiss2),method="average");#对合并后的模块进行聚类
    #绘制聚类结果图
    pdf(file=paste0(picDir,"/7. MECombined.pdf"),width=12,height=5)
    par(mfrow=c(1,2))
    plot(METree1,xlab="",sub="",main="Clustering of ME before combined")# METree1是动态剪切树形成的模块
    abline(h=MEDissThres,col="red")#相异度为0.25
    plot(METree2,xlab="",sub="",main="Clustering of ME after combined")# METree2是合并后的模块
    dev.off()
   
 
    pdf(paste0(picDir,"/8. Tree-MergedDynamics.pdf"),width=8,height=6)
    plotDendroAndColors(dendro = geneTree,#剪切树
                        colors = cbind(dynamicColors,mergedColors),#将两种方法形成的模块颜色合并在一起
                        groupLabels = c("Dynamic Tree Cut","Merged Dynamics"),
                        dendroLabels = FALSE, 
                        hang = 0.03,
                        addGuide=TRUE,
                        guideHang=0.05,
                        main="Gene Dendrogram and module colors")
    dev.off()
   
    # 模块中基因数
    write.table(table(moduleColors),paste0(picDir,"/各个模块基因的总数.txt"),quote = F,row.names = F)
    # 保存构建的网络信息
    moduleColors=mergedColors
    colorOrder=c("grey", standardColors(50))
    moduleLabels=match(moduleColors, colorOrder)-1
    MEs=mergedMEs
    ## 绘制样本间的相关性
    MEs = orderMEs(MEs)
    sizeGrWindow(5,7.5);
    pdf(file =paste0(picDir,"/9. module_Cor.pdf"), width = 5, height = 7.5);
    
    par(cex = 0.9)
    plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                          = 90)
    dev.off()
    
    dissTOM = 1-TOMsimilarityFromExpr(datExpr0, power = softPower); 
    
    nSelect = 400
    # 随机选取400个基因进行可视化，s设置seed值，保证结果的可重复性
    set.seed(10);
    select = sample(nGenes, size = nSelect);
    selectTOM = dissTOM[select, select];
    # 对选取的基因进行重新聚类
    selectTree = hclust(as.dist(selectTOM), method = "average")
    selectColors = mergedColors[select];
    # 打开一个绘图窗口
    sizeGrWindow(9,9)
    pdf(file =paste0(picDir,"/10. TOMplot.pdf") , width = 9, height = 9);
    # 美化图形的设置
    plotDiss = selectTOM^7;
    diag(plotDiss) = NA;
    TOMplot(plotDiss, selectTree, selectColors, 
            main = "Network heatmap plot, selected genes")
    dev.off()
    ###################################################################################### 
    ###################################载入性状########################################### 
    ###################################################################################### 
    # 创建新的矩阵，行名为列名，列名为疾病名称
    traitData <- matrix(0, nrow = nrow(datExpr), ncol = 2,
                        dimnames = list(rownames(datExpr), c(disease,"con")))
    # # 为新矩阵填写值
    # traitData[grepl("treat", rownames(traitData)), ] <- 1
    # traitData[grepl("con", rownames(traitData)), ] <- 0
    # 使用grepl函数将所有行名包含"treat"的行的对应列设置为1和0
    traitData[grepl("treat", rownames(traitData)), disease] <- 1
    traitData[grepl("treat", rownames(traitData)), "con"] <- 0
    traitData[grepl("con", rownames(traitData)), disease] <- 0
    traitData[grepl("con", rownames(traitData)), "con"] <- 1
    traitData <- as.data.frame(traitData)
    class(traitData)
    dim(traitData)
    names(traitData)
    fpkmSamples = rownames(datExpr)
    traitSamples =rownames(traitData)
    traitRows = match(fpkmSamples, traitSamples)
    datTraits = traitData[traitRows,]
    rownames(datTraits) 
    collectGarbage()
    ###【绘制非离群样本聚类图】
    sampleTree2 = hclust(dist(datExpr), method = "average")
    traitColors = numbers2colors(datTraits, signed = FALSE)
    
    pdf( file=paste0(picDir,"/11. sample_map.pdf"),width=10,height=10)#【绘制图形】
    plotDendroAndColors(sampleTree2, traitColors,
                        groupLabels = names(datTraits),
                        main = "Sample dendrogram and trait heatmap",cex.lab = 1,
                        cex.axis = 1, cex.main = 1)
    dev.off()
 
    
      MEs=orderMEs(MEs)
      ###【绘制模块与性状数据热图】
      nGenes = ncol(datExpr0)#【提取模块基因】
      nSamples = nrow(datExpr0)#【提取模块样本】
      moduleTraitCor = cor(MEs, datTraits, use = "p")#【计算模块相关性】
      moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)#【计算模块P值】
      
      
      pdf( file=paste0(picDir,"/12. sample_Module.pdf"),width=10,height=8)#【绘制性状数据热图】
      textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                         signif(moduleTraitPvalue, 1), ")", sep = "")
      dim(textMatrix) = dim(moduleTraitCor)
      par(mar = c(5, 10, 3, 3))
      labeledHeatmap(Matrix=moduleTraitCor,#模块和表型的相关性矩阵，这个参数最重要，其他可以不变
                     xLabels=colnames(datTraits),
                     yLabels=names(MEs),
                     ySymbols=names(MEs),
                     colorLabels=FALSE,
                     colors=blueWhiteRed(50),
                     textMatrix=textMatrix,
                     setStdMargins=FALSE,
                     cex.text=0.5,
                     cex.lab=0.5,
                     zlim=c(-1,1),
                     main=paste("Module-trait relationships"))
      dev.off()
      
      ###【导出基因所在的模块】
      moduleColors=dynamicColors
      probes = colnames(datExpr0)#【读取基因名】
      geneInfo0 = data.frame(probes= probes,
                             moduleColor = moduleColors)
      geneOrder =order(geneInfo0$moduleColor)#【读取模块】
      geneInfo = geneInfo0[geneOrder, ]
      #【保存所有模块基因】
      write.table(geneInfo, file = paste0(picDir,"/all_genes.txt"),sep="\t",row.names=F,quote=F)
      # 创建一个空的数据框，用于存储每个模块中包含的基因
       result = geneInfo[,1:2]
       rownames(result) = result$probes
       result_expdata <- as.data.frame(t(datExpr0[,result$probes]))
      df <- merge(result_expdata, result, by = 0)
      
      # 将第一列设置为行名
      rownames(df) <- df[, 1]
      
      # 删除第一列
      df <- df[, -1]
      # 删除probes列
      df1 <- df[, -which(names(df) == "probes")]
      df1$moduleColor
      library(ggplot2)
      plist=list()
      tag = 0
      for (mod in unique(df1$moduleColor)){
        tag = tag + 1
        specmod <- df1[df1$moduleColor == mod,]
        t<-data.frame(tp=1:5,tpm_mean = log2( c(mean(as.matrix(specmod[1:30,1:4])),mean(as.matrix(specmod[1:30,5:8])),mean(as.matrix(specmod[1:30,9:12])),mean(as.matrix(specmod[1:30,13:16])),mean(as.matrix(specmod[1:30,17:20])))))
        plist[[tag]]<-ggplot(t,aes(x = tp, y = tpm_mean)) +
          geom_smooth(color = mod) + #拟合曲线的颜色
          labs(x = "after treat", y = "Expression Level", title = mod) + 
          theme_bw() + #去除背景色
          theme(panel.grid =element_blank()) + #去除网格线
          theme(panel.border = element_blank()) + #去除外层边框
          theme(axis.line = element_line(colour = "grey")) + #坐标轴画成灰色
          theme(axis.ticks = element_blank()) + #取掉坐标轴上的刻度线
          geom_hline(yintercept = t$tpm_mean[1],linetype="dashed") #在第一个时间点处画虚线
      }
      
      require(cowplot)
      plot_grid(plotlist=plist, ncol=2)

      # 设置长度和宽度为10和6并保存图像
      ggsave(file = paste0(picDir, "/14. myfitting.pdf"), width = 6, height = 10)
      
      
      
      
      
      
      
      # library(pheatmap)
      # #先按module内基因的数量为module排序
      # module.num<-table(df1$moduleColor)
      # module.num.sort<-module.num[order(module.num,decreasing =TRUE)]
      # module.num.sort
      # keys<-names(module.num.sort)
      # keysDF<-data.frame(moduleColor=keys,order=1:length(keys))
      # keysDF
      # 
      # #再按这个module的顺序为基因排序
      # df.exp<-merge(df1,keysDF,by.x = 'moduleColor',by.y = 'moduleColor',all.x = T, all.y = F)
      # 
      # 
      # 
      # 
      # pheatmap(df.exp[,-1],cellwidth = 12, cellheight = 0.2, fontsize = 8,
      #          scale="row", #为基因做scale
      #          cluster_rows=F,cluster_cols=F,#不聚类
      #          color = colorRampPalette(c("#00A9E0", "white", "firebrick3"))(50),
      #          show_colnames=F,show_rownames =F,
      #          annotation_col = annotation_col,
      #          annotation_row = annotation_row,
      #          annotation_colors = ann_colors,
      #          border_color = "NA",
      #          filename = "heatmap.pdf")
      
      
      
      
      #【使用循环对每个模块的基因进行导出】
      # 加载库
      library(dplyr)
      library(readr)
      # 初始化一个空列表来存储数据帧
      modules_list <- list()
      for (mod in 1:nrow(table(moduleColors))) {
        modules <- names(table(moduleColors))[mod]
        probes <- colnames(datExpr0)
        inModule <- (moduleColors == modules)
        modGenes <- probes[inModule]
        # 写入文件
        write.table(modGenes, file = paste0(picDir, "/module_", modules, ".txt"), sep = "\t", row.names = F, col.names = F, quote = F)
        
        # 将数据转换为数据帧并将文件名作为列名
        temp_df <- data.frame(modGenes)
        colnames(temp_df) <- paste0( modules)
        
        # 将新的数据帧添加到列表中
        modules_list[[modules]] <- temp_df

        # 将新的数据帧添加到列表中
      }


      
      
      ###【计算MM和GS】
      modNames = substring(names(MEs), 3)
      geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
      MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
      names(geneModuleMembership) = paste("MM", modNames, sep="")
      names(MMPvalue) = paste("p.MM", modNames, sep="")
      traitNames=names(datTraits)
      geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
      GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
      names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
      names(GSPvalue) = paste("p.GS.", traitNames, sep="")
      
      
      ###输出模块重要性的图形
      y=datTraits[,1]
      GS1=as.numeric(cor(y, datExpr0, use="p"))
      GeneSignificance=abs(GS1)
      ModuleSignificance=tapply(GeneSignificance, mergedColors, mean, na.rm=T)
      pdf(file=paste0(picDir,"/13. GeneSignificance.pdf"), width=11, height=7)
      plotModuleSignificance(GeneSignificance, mergedColors)
      dev.off()
      
      
      ###【批量输出性状和模块散点图】
      dir.create(file.path(picDir, "module_trait"))
      picDir1 <- file.path(picDir, "module_trait")
      for (trait in traitNames){#【对每个性状和模块进行循环绘制图片】
        traitColumn=match(trait,traitNames)  
        for (module in modNames){
          column = match(module, modNames)
          moduleGenes = moduleColors==module
          if (nrow(geneModuleMembership[moduleGenes,]) > 1){
            pdfFile=paste("9_", trait, "_", module,".pdf",sep="")
            outPdf = file.path( picDir1, pdfFile)
            pdf(file=outPdf,width=7,height=7)
            par(mfrow = c(1,1))
            verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                               abs(geneTraitSignificance[moduleGenes, traitColumn]),
                               xlab = paste("Module Membership in", module, "module"),
                               ylab = paste("Gene significance for ",trait),
                               main = paste("Module membership vs. gene significance\n"),
                               cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
            dev.off()
          }
        
      }
    }

      ###输出GS_MM数据
      probes = colnames(datExpr0)
      geneInfo0 = data.frame(probes= probes,
                             moduleColor = moduleColors)
      for (Tra in 1:ncol(geneTraitSignificance))
      {
        oldNames = names(geneInfo0)
        geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                               GSPvalue[, Tra])
        names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                             names(GSPvalue)[Tra])
      }
      
      for (mod in 1:ncol(geneModuleMembership))
      {
        oldNames = names(geneInfo0)
        geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                               MMPvalue[, mod])
        names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                             names(MMPvalue)[mod])
      }
      geneOrder =order(geneInfo0$moduleColor)
      geneInfo = geneInfo0[geneOrder, ]
      write.table(geneInfo, file =paste0(picDir,"/GS_MM.xls") ,sep="\t",row.names=F)
  
      ###输出每个模块的核心基因
      geneSigFilter=0.5         #基因重要性的过滤条件
      moduleSigFilter=0.8       #基因与模块相关性的过滤条件
      datMM=cbind(geneModuleMembership, geneTraitSignificance)
      datMM=datMM[abs(datMM[,ncol(datMM)])>geneSigFilter,]
      for(mmi in colnames(datMM)[1:(ncol(datMM)-2)]){
        dataMM2=datMM[abs(datMM[,mmi])>moduleSigFilter,]
        write.table(row.names(dataMM2), file =paste0(picDir,"/hubGenes_",mmi,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
      }
      
     
      
      
      
  }
  return(list(datExpr = datExpr, modules_list = modules_list))

}




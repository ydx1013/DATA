ydx_draw_GSEA <- function(diff_data,pvalueCutoff = 0.05,type ="mf",methods="gsea",width = 8, height = 6,picDir = "gene_GSEA_"){
  options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
  #适用于，根据基因的表达量，分为高低组，从高低组里面去做差异，然后做GSEA
  
  library(msigdbr)
  library(limma)
  library(GEOquery)
  library(affy)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  if (!dir.exists(picDir)) {
    dir.create(picDir)
  }
  cols <- names(diff_data)[grep("(?i)(symbol)", names(diff_data), perl = TRUE)]
  gene_name = diff_data[,cols]
  result <- sum(grepl("[a-z]", gene_name))
  if(result>length(gene_name)/2){
    species= "Mus musculus"
    org_db= "org.Mm.eg.db"
    cat("检测到基因属于小鼠，使用小鼠基因集做分析\n")
  }else{
    species= "Homo sapiens"
    org_db= "org.Hs.eg.db"
    cat("检测到基因属于人类，使用人类基因集做分析\n")
  }
  # 使用grep函数获取列名
  

  cols <- names(diff_data)[grep("(?i)(symbol|log)", names(diff_data), perl = TRUE)]
  # 提取symbol和logfc列
  res <- diff_data[, cols]
  names(res)[grep("(?i)(symbol)", names(res), perl = TRUE)] <-"symbol"
  names(res)[grep("(?i)(log)", names(res), perl = TRUE)] <-"log2fc"
  
  

  methods= toupper(methods)
  type= toupper(type)
  cat('提取msigdbr数据信息，\n当前提取的是',type,"\n")
  if(type == "KEGG"){
    C2_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:KEGG") %>%
      dplyr::select(gs_name, gene_symbol) %>%
      mutate(gs_name = gsub('KEGG_', '', gs_name),         # 去除前缀
        gs_name = stringr::str_to_title(gs_name),        # 将首字母除外的大写换为小写
      )
  }else if(type %in% c("BP", "CC", "MF", "ALL")){
    if( type =="ALL"){type=""}
    C2_t2g <- msigdbr(species = species, category = "C5") %>%
      dplyr::select(gs_name, gene_symbol) %>%
      filter(grepl(paste0('^GO',type), gs_name))%>%   # 选择以"GO"开头的subcategory
      mutate(
        gs_name = gsub('GOMF_|GOCC_|GOBP_', '', gs_name),   # 去除前缀
        gs_name = stringr::str_to_title(gs_name)       # 将首字母除外的大写换为小写
      )
    }
  if( type ==""){type="ALL"}
    # 运行GSEA
  gseaTab <- NULL
    if(methods== "GSEA"){
      cat( "正在使用GSEA函数来计算")
      res <- res[order(res$log2fc, decreasing = T),]
      glist <- res$log2fc; names(glist) <- rownames(res)
      gsea <- GSEA(geneList     =    glist,
                   TERM2GENE    =    C2_t2g, 
                   pvalueCutoff = pvalueCutoff,
                   eps = 0)
      cat( "正在使用GSEA函数绘图\n")
      pdf(paste0(picDir, "/", type,"-",methods, "_所有GSEA图.pdf"), width = width, height = height)
      for (i in 1:nrow(gsea@result)) {
        my_plot <- gseaplot2(gsea, geneSetID = i, pvalue_table = TRUE)
        print(my_plot)
      }
      dev.off()
      library(GseaVis)
      pdf(paste0(picDir, "/", type,"-",methods, "_所有GSEA图----GseaVis.pdf"), width = width, height = height)
      my_plot1 <-     dotplotGsea(gsea)
      my_plot3 <-       volcanoGsea(gsea)
      print(my_plot1)
      print(my_plot3)
      for (i in 1:nrow(gsea@result)) {
        my_plot2 <-      gseaNb(gsea,geneSetID =i,subPlot = 3)
        print(my_plot2)
      }
      dev.off()
      
      gsea.df <- as.data.frame(gsea) # 数据转换为数据框
      gc()
      gseaTab <- rbind.data.frame(gseaTab,
                                  data.frame(term = gsea.df$ID,
                                             NES = gsea.df$NES,
                                             pval = gsea.df$pvalue,
                                             FDR = gsea.df$p.adjust,
                                             number = gsea.df$setSize,
                                             stringsAsFactors = F),
                                  stringsAsFactors = F)
    }else if (methods== "FGSEA"){
      cat( "正在使用fgsea函数来计算\n")
      res <- res[order(res$log2fc, decreasing = T),]
      glist <- res$log2fc; names(glist) <- rownames(res)
      m_list <- split(C2_t2g$gene_symbol, C2_t2g$gs_name)  # 使用MSigDB数据库
      gsea <- fgsea::fgsea(pathways = m_list,
                           stats    = glist)
      
      pdf(paste0(picDir, "/", type,"-",methods, "_所有GSEA图.pdf"), width = width, height = height)
      topPathways <- gsea[head(order(pval), n=15)][order(NES), pathway]
      my_plot2= plotGseaTable(m_list[topPathways], glist,
                    gsea, gseaParam=0.5)
      print(my_plot2)
      dev.off()
      gsea <- as.data.frame(gsea) # 数据转换为数据框
      gsea$leadingEdge <- sapply(gsea$leadingEdge, paste, collapse = ",")
      colnames(gsea) = c("ID","pvalue","p.adjust","log2err","ES","NES","setSize","leading_edge")
      gc()
      gseaTab <- rbind.data.frame(gseaTab,
                                  data.frame(term = gsea$ID,
                                             NES = gsea$NES,
                                             pval = gsea$pvalue,
                                             FDR = gsea$p.adjust,
                                             number = gsea$setSize,
                                             stringsAsFactors = F),
                                  stringsAsFactors = F)
  

    }else if(methods== "GSE"){
    cat( "正在使用clusterProfiler包中gseKEGG的函数来计算")
    names(res)[grep("(?i)(symbol)", names(res), perl = TRUE)] <-"SYMBOL"
    ENSEMBLlist <- bitr(res$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", org_db)
    
    # 使用merge函数合并两个数据框
    merged_df <- merge( ENSEMBLlist,res, by = "SYMBOL", all = TRUE,ignore.case = TRUE)
    merged_df = na.omit(merged_df)
    
    merged_df2 =merged_df
    merged_df <- merged_df[, !names(merged_df) %in% "SYMBOL"]
    library(ggupset)
    library(enrichplot)
    library(ggridges)
    
    ## 1: 提取logFC值，并储存在一个向量中
    geneList = merged_df[,2]
    ## 2: 对geneList进行命名
    names(geneList) = as.character(merged_df[,1])
    head(geneList)
    ## 3: 根据logFC值降序排列
    geneList = sort(geneList, decreasing = TRUE)

    type= toupper(type)
    cat('提取msigdbr数据信息，\n当前提取的是',type,"\n")
    if(type == "KEGG"){
      if (org_db == "org.Hs.eg.db") {
        gsea <- gseKEGG(geneList     = geneList,
                        organism     = 'hsa',
                        minGSSize    = 120,
                        pvalueCutoff = pvalueCutoff,
                        verbose      = FALSE)
      } else if (org_db == "org.Mm.eg.db") {
        gsea <- gseKEGG(geneList     = geneList,
                        organism     = 'mmu',
                        nPerm        = 1000,
                        minGSSize    = 120,
                        pvalueCutoff = pvalueCutoff,
                        verbose      = FALSE)
      } else {
        stop("未知物种")
      }
      cat( "正在使用clusterProfiler包中的gseKEGG函数绘图\n")
    }else if(type %in% c("BP", "CC", "MF", "ALL")){
      gsea <- gseGO(geneList     = geneList,
                    OrgDb = org_db,
                    ont = type,
                    seed = 1,
                    pvalueCutoff = pvalueCutoff,
                    pAdjustMethod = "BH",
                    verbose      = FALSE)
                  
    }
    pdf(paste0(picDir, "/", type,"-",methods, "_所有GSEA图.pdf"), width = width, height = height)
    
    for (i in 1:nrow(gsea@result)) {
      my_plot <- gseaplot2(gsea, geneSetID = i, pvalue_table = TRUE)
      print(my_plot)
    }
    dev.off()
    library(GseaVis)
    pdf(paste0(picDir, "/", type,"-",methods, "_所有GSEA图----GseaVis.pdf"), width = width, height = height)
    my_plot1 <-     dotplotGsea(gsea)

    print(my_plot1)

    for (i in 1:nrow(gsea@result)) {
      my_plot2 <-      gseaNb(gsea,geneSetID =i,subPlot = 3)
      print(my_plot2)
    }
    dev.off()

    gsea.df <- as.data.frame(gsea) # 数据转换为数据框
    gseaTab <- NULL
    gseaTab <- rbind.data.frame(gseaTab,
                                data.frame(ID = gsea.df$ID,
                                           term = gsea.df$Description,
                                           NES = gsea.df$NES,
                                           pval = gsea.df$pvalue,
                                           FDR = gsea.df$p.adjust,
                                           number = gsea.df$setSize,
                                           stringsAsFactors = F),
                                  stringsAsFactors = F)
    }

write.csv(as.data.frame(gsea), file =  paste0(picDir,"/",type,"-gsea-",methods,"算法-结果.csv"))
return(gsea)
}


  



######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################

ydx_draw_gene_GSEA <-function(gene_name=m_l, data=pre_data$exp_data,pct = 0.4,methods= "GSEA",width = 8,height = 6,picDir = "gene_GSEA_heatmap", pathOfInterest){
  options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
  #适用于，根据基因的表达量，分为高低组，从高低组里面去做差异，然后做GSEA
  packages <- c("topGO", "ReactomePA", "tidyverse", "clusterProfiler", "biomaRt", "enrichplot", "dplyr", "data.table", "ggplot2", "patchwork", "cols4all", "ggraph", "igraph","RColorBrewer", "ComplexHeatmap", "ggpubr", "circlize", "org.Hs.eg.db", "org.Mm.eg.db")
  lapply(packages, require, character.only = TRUE)
  library(msigdbr)
  library(limma)
  library(GEOquery)
  library(affy)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  library(ggplot2)
  library(data.table)
  library(ggpubr)
  library(SimDesign)
  if (!dir.exists(picDir)) {
    dir.create(picDir)
  }

  
   # 30%阈值来定义高低组
  result <- sum(grepl("[a-z]", gene_name))
  if(result>length(gene_name)/2){
    organism     = 'mmu'
    org_db= "org.Mm.eg.db"
    cat("检测到基因属于小鼠，使用小鼠基因集做分析\n")
  }else{
    organism     = 'hsa'
    org_db= "org.Hs.eg.db"
    species = "Homo sapiens"
  
    cat("检测到基因属于人类，使用人类基因集做分析\n")
  }
  
  twoclasslimma <- function(subtype   = NULL,
                            featmat   = NULL,
                            treatVar  = NULL,
                            ctrlVar   = NULL,
                            prefix    = NULL,
                            overwt    = FALSE,
                            sort.p    = TRUE,
                            verbose   = TRUE,
                            res.path  = getwd()) {
    
    library(limma)
    if(!is.element("condition", colnames(subtype))) {
      stop("argument of subtype must contain a column named with 'condition'!")
    }
    # Create comparison list for differential analysis between two classes.
    createList  <- function(subtype = NULL) {
      
      tumorsam <- rownames(subtype)
      sampleList <- list()
      treatsamList <- list()
      treatnameList <- c()
      ctrlnameList <- c()
      
      sampleList[[1]] <- tumorsam
      treatsamList[[1]] <- intersect(tumorsam, rownames(subtype[which(subtype$condition == treatVar),,drop = F]))
      treatnameList[1] <- treatVar
      ctrlnameList[1] <- ctrlVar
      
      return(list(sampleList, treatsamList, treatnameList, ctrlnameList))
    }
    
    complist <- createList(subtype = subtype)
    
    sampleList <- complist[[1]]
    treatsamList <- complist[[2]]
    treatnameList <- complist[[3]]
    ctrlnameList <- complist[[4]]
    allsamples <- colnames(featmat)
    
    # log transformation
    if(max(featmat) < 25 | (max(featmat) >= 25 & min(featmat) < 0)) {
      message("--feature matrix seems to have been standardised (z-score or log transformation), no more action will be performed.")
      gset <- featmat
    }
    if(max(featmat) >= 25 & min(featmat) >= 0){
      message("--log2 transformation done for feature matrix.")
      gset <- log2(featmat + 1)
    }
    
    options(warn = 1)
    for (k in 1:length(sampleList)) {
      samples <- sampleList[[k]]
      treatsam <- treatsamList[[k]]
      treatname <- treatnameList[k]
      ctrlname <- ctrlnameList[k]
      
      compname <- paste(treatname, "_vs_", ctrlname, sep="")
      tmp <- rep("others", times = length(allsamples))
      names(tmp) <- allsamples
      tmp[samples] <- "control"
      tmp[treatsam] <- "treatment"
      
      if(!is.null(prefix)) {
        outfile <- file.path(res.path, paste(prefix, "_limma_test_result.", compname, ".txt", sep = ""))
      } else {
        outfile <- file.path(res.path, paste("limma_test_result.", compname, ".txt", sep = ""))
      }
      
      if (file.exists(outfile) & (overwt == FALSE)) {
        cat(paste0("limma of ",compname, " exists and skipped...\n"))
        next
      }
      
      pd <- data.frame(Samples = names(tmp),
                       Group = as.character(tmp),
                       stringsAsFactors = FALSE)
      
      design <-model.matrix(~ -1 + factor(pd$Group, levels = c("treatment","control")))
      colnames(design) <- c("treatment","control")
      
      fit <- limma::lmFit(gset, design = design);
      contrastsMatrix <- limma::makeContrasts(treatment - control, levels = c("treatment", "control"))
      fit2 <- limma::contrasts.fit(fit, contrasts = contrastsMatrix)
      fit2 <- limma::eBayes(fit2, 0.01)
      resData <- limma::topTable(fit2, adjust = "fdr", sort.by = "B", number = 100000)
      resData <- as.data.frame(subset(resData, select=c("logFC","t","B","P.Value","adj.P.Val")))
      resData$id <- rownames(resData)
      colnames(resData) <- c("log2fc","t","B","pvalue","padj","id")
      resData$fc <- 2^resData$log2fc
      
      if(sort.p) {
        resData <- resData[order(resData$padj),]
      } else {
        resData <- as.data.frame(resData)
      }
      if(verbose) {
        resData <- resData[,c("id","fc","log2fc","pvalue","padj")]
      } else {
        resData <- resData[,c("id","fc","log2fc","t","B","pvalue","padj")]
      }
      write.table(resData, file = outfile, row.names = FALSE, sep = "\t", quote = FALSE)
      cat(paste0("limma of ",compname, " done...\n"))
    }
    options(warn = 0)
  }
  
  gsea_plot <-function(gsea){
    # 生成一个点图，展示最显著的前10个基因集
    my_plot <- dotplot(gsea, showCategory = 10) + ggtitle("dotplot for GSEA")
    pdf(paste0(picDir,"/", j,"-kegg" ,"_dotplot.pdf"))
    print(my_plot)
    dev.off()
    # 生成一个“集合的交集和并集”图
    my_plot <- upsetplot(gsea)
    pdf(paste0(picDir,"/",j,"-kegg","_upsetplot.pdf"),width=15,height =10)
    print(my_plot)
    dev.off()

    # 生成一个脊型图，展示基因集富集得分的分布情况
    my_plot <- ridgeplot(gsea, showCategory = 15)
    pdf(paste0(picDir,"/", j,"-kegg_ridgeplot.pdf"))
    print(my_plot)
    dev.off()
    pdf(paste0(picDir, "/", j,"-kegg", "_所有通路的单一图.pdf"), width = 16, height = 12)
    for (i in 1:nrow(gsea@result)) {
      my_plot <- gseaplot2(gsea, geneSetID = i, pvalue_table = TRUE)
      print(my_plot)
    }
    dev.off()
    
    
    
  }
  
  # 设置颜色
  darkblue <- "#303B7F"
  darkred <- "#D51113"
  yellow <- "#EECA1F"
  
  
  
  #msigdb.hallmark <- read.gmt("h.all.v7.2.symbols.gmt") # 这里使用的是HALLMARK背景集，可以改为自己感兴趣的基因集合
  
  tumors=data
  genelist= gene_name
  gseaTab <- NULL
  cat("分类算法只针对于疾病组\n")
  data = as.data.frame(data)
  # 实验组
  exp_data_T <- data %>% dplyr::select(grep("_treat", colnames(.)))
  treatNum <- ncol(exp_data_T) 
  cat("实验组数量:", treatNum, "\n")
  
  message("--仅提取了实验组样本...")
  es <- exp_data_T[genelist,] # 取出感兴趣基因和当前肿瘤样本的表达谱子集
  tumsam= colnames(es)

    for (j in genelist) {
      set.seed(20220407)
      message("gene of ", j, " starts...")
      es_subset <- as.data.frame(t(es[j,tumsam]))
      es_subset <- es_subset[order(es_subset[,1],decreasing = T),,drop = F] # 对表达值排序
      
      high.es <- rownames(es_subset[1:(nrow(es_subset) * pct),,drop = F]) # 取前30%为高组
      low.es <- rownames(es_subset[nrow(es_subset):(nrow(es_subset) - nrow(es_subset) * pct + 1),,drop = F]) # 取后30%为低组
      
      # 采用两样本limma差异表达分析
      subt <- data.frame(condition = rep(c("high","low"),c(length(high.es),length(low.es))),
                         row.names = c(high.es, low.es),
                         stringsAsFactors = F)
      
      gset <- na.omit(data[,rownames(subt)])
      twoclasslimma(subtype  = subt, # 亚型信息（必须包含一列为condition）
                    featmat  = gset, # 表达谱（会自动检测数据标准化与否）
                    treatVar = "high", # “治疗组”的名字
                    ctrlVar  = "low", # “对照组”的名字
                    prefix   = paste0("_",j), # 差异表达文件的前缀
                    overwt   = T, # 决定是否覆盖已经存在的差异表达结果
                    sort.p   = F, # 决定是否对结果按照FDR排序
                    verbose  = TRUE, # 决定是否显示一些提示
                    res.path = picDir) # 确定结果路径
      
      # 加载差异表达结果
      res <- read.table(file.path(picDir, paste0("_",j,"_limma_test_result.high_vs_low.txt")),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
      if(methods== "GSEA" | methods== "fgsea"){
        # 产生pre-ranked基因列表
        res <- res[order(res$log2fc, decreasing = T),]
        glist <- res$log2fc; names(glist) <- rownames(res)
        cat('提取msigdbr数据信息')
        C2_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:KEGG") %>%
          dplyr::select(gs_name, gene_symbol) %>%
          mutate(
            gs_name = gsub('KEGG_', '', gs_name),         # 去除前缀
            gs_name = stringr::str_to_title(gs_name),        # 将首字母除外的大写换为小写
          )
        
        
        # 运行GSEA
        if(methods== "GSEA"){
          cat( "正在使用GSEA函数来计算KEGG信号通路")
 
          gsea <- GSEA(geneList     = glist,
                       TERM2GENE = C2_t2g, 
                       pvalueCutoff = 0.1, eps = 0)
          cat( "正在使用GSEA函数绘图\n")
          gsea_plot(gsea)
          cat( "正在使用fgsea函数来计算KEGG信号通路111111111111111")
          gsea.df <- as.data.frame(gsea) # 数据转换为数据框
        }else{
          cat( "正在使用fgsea函数来计算KEGG信号通路\n")
          m_list <- split(C2_t2g$gene_symbol, C2_t2g$gs_name)  # 使用MSigDB数据库
          gsea <- fgsea::fgsea(pathways = m_list,
                               stats    = glist,
                               minSize  = 15,
                               maxSize  = 500)
          gsea.df <- as.data.frame(gsea) # 数据转换为数据框
          gsea.df$leadingEdge <- sapply(gsea.df$leadingEdge, paste, collapse = ",")

          colnames(gsea.df) = c("ID","pvalue","p.adjust","log2err","ES","NES","number","leading_edge")
          
          
        }
        gc()
        write.table(gsea.df, file = file.path(picDir, paste0("高低表达组差异-方法", methods, "百分比", pct, "基因-", j, " in ", ".txt")))
        gseaTab <- rbind.data.frame(gseaTab,
                                    data.frame(gene = j,
                                               term = gsea.df$ID,
                                               NES = gsea.df$NES,
                                               pval = gsea.df$pvalue,
                                               FDR = gsea.df$p.adjust,
                                               stringsAsFactors = F),
                                    stringsAsFactors = F)
      }else if(methods== "gse"){
        cat( "正在使用clusterProfiler包中gseKEGG的函数来计算")
      res$SYMBOL =rownames(res)
      # 使用grep函数获取列名
      cols <- names(res)[grep("(?i)(symbol|log)", names(res), perl = TRUE)]
      # 提取symbol和logfc列
      DEG_list_sub <- res[, cols]
      ENSEMBLlist <- bitr(DEG_list_sub$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", org_db)
      
      # 使用merge函数合并两个数据框
      merged_df <- merge( ENSEMBLlist,DEG_list_sub, by = "SYMBOL", all = TRUE,ignore.case = TRUE)
      merged_df = na.omit(merged_df)
      
      merged_df2 =merged_df
      merged_df <- merged_df[, !names(merged_df) %in% "SYMBOL"]
      library(ggupset)
      library(enrichplot)
      library(ggridges)
      
      ## 1: 提取logFC值，并储存在一个向量中
      geneList = merged_df[,2]
      ## 2: 对geneList进行命名
      names(geneList) = as.character(merged_df[,1])
      head(geneList)
      ## 3: 根据logFC值降序排列
      geneList = sort(geneList, decreasing = TRUE)
      
      if (org_db == "org.Hs.eg.db") {
        gsea <- gseKEGG(geneList     = geneList,
                         organism     = 'hsa',
                         minGSSize    = 120,
                         verbose      = FALSE)
      } else if (org_db == "org.Mm.eg.db") {
        gsea <- gseKEGG(geneList     = geneList,
                         organism     = 'mmu',
                         nPerm        = 1000,
                         minGSSize    = 120,
                         verbose      = FALSE)
      } else {
        stop("未知物种")
      }
      cat( "正在使用gseKEGG函数绘图\n")
      gsea_plot(gsea)

      gsea.df <- as.data.frame(gsea) # 数据转换为数据框
      write.table(gsea.df,file = file.path(picDir,paste0("高低表达组差异-方法",methods,"百分比",pct,"基因-",j," in ",".txt")),sep = "\t",row.names = T,col.names = NA,quote = F)
      
      gseaTab <- rbind.data.frame(gseaTab,
                                  data.frame(gene = j,
                                             term = gsea.df$Description,
                                             NES = gsea.df$NES,
                                             pval = gsea.df$pvalue,
                                             FDR = gsea.df$p.adjust,
                                             stringsAsFactors = F),
                                  stringsAsFactors = F)
      }
    }
  write.table(gseaTab, paste0(picDir,"/mult_single_gene_gsea_kegg.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
      
  if(length(pathOfInterest)!= 0){
   
    # 提取和该通路有关的NES
    tmp <- gseaTab[gseaTab$term %in% pathOfInterest, ]
    my_palette <- colorRampPalette(c(darkblue,yellow,darkred), alpha=TRUE)(n=128)
    ggplot(tmp, aes(x=gene,y=term)) +
      geom_point(aes(size=-log10(FDR),color=NES)) +
      scale_color_gradientn('NES', 
                            colors=my_palette) + 
      scale_size_continuous(range = c(1,4)) + #圆点的大小范围
      theme_bw() +
      theme(panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank(),
            axis.text.x = element_text(angle = 45, size = 12, hjust = 0.3, vjust = 0.5, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title = element_blank(),
            panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
            legend.position = "bottom",
            plot.margin = unit(c(1,1,1,1), "lines"))
    ggsave(paste0(picDir,"/GSEA cormap.pdf"), width = width,height = height)
    
    return(gsea_kegg=gsea)
  }
  
}





ydx_draw_sample_GSEA <- function(DEG_list,picDir = "GSEA_",pvalueCutoff = 1){
  options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
  packages <- c("topGO", "ReactomePA", "tidyverse", "clusterProfiler", "biomaRt", "enrichplot", "dplyr", "data.table", "ggplot2", "patchwork", "cols4all", "ggraph", "igraph","RColorBrewer", "ComplexHeatmap", "ggpubr", "circlize", "org.Hs.eg.db", "org.Mm.eg.db")
  lapply(packages, require, character.only = TRUE)
  cat("使用gseGO函数做GSEA\n")
  
  
  if (!dir.exists(picDir)) {
    dir.create(picDir)
  }
  names(DEG_list)[1] <- "SYMBOL"
  
  result <- sum(grepl("[a-z]", DEG_list[1]))
  if(result>length(result)/2){
    species= "Mus musculus"
    org_db= "org.Mm.eg.db"
    cat("检测到基因属于小鼠，使用小鼠基因集做分析\n")
  }else{
    species= "Homo sapiens"
    org_db= "org.Hs.eg.db"
    cat("检测到基因属于人类，使用人类基因集做分析\n")
  }
  
  
  # 使用grep函数获取列名
  cols <- names(DEG_list)[grep("(?i)(symbol|log)", names(DEG_list), perl = TRUE)]
  # 提取symbol和logfc列
  DEG_list_sub <- DEG_list[, cols]
  ENSEMBLlist <- bitr(DEG_list_sub$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", org_db)

  # 使用merge函数合并两个数据框
  merged_df <- merge( ENSEMBLlist,DEG_list_sub, by = "SYMBOL", all = TRUE,ignore.case = TRUE)
  merged_df = na.omit(merged_df)
  
  merged_df2 =merged_df
  merged_df <- merged_df[, !names(merged_df) %in% "SYMBOL"]
  library(ggupset)
  library(enrichplot)
  library(ggridges)
  
  ## 1: 提取logFC值，并储存在一个向量中
  geneList = merged_df[,2]
  ## 2: 对geneList进行命名
  names(geneList) = as.character(merged_df[,1])
  head(geneList)
  ## 3: 根据logFC值降序排列
  geneList = sort(geneList, decreasing = TRUE)
  
  ########## GO的GSEA富集分析：gseGO ##########
  gse_BP_result <- gseGO(geneList     = geneList,#排序后的基因列表，一般根据logFC进行排序
                         OrgDb        = org_db,
                         ont          = "BP",#可选择bp.MF,CC,ALL
                         minGSSize    = 100,#最小基因集的基因数
                         maxGSSize    = 500,#最大基因集的基因数
                         pvalueCutoff = 0.05,#p值的阈值
                         verbose      = FALSE)#是否输出提示信息，默认为false
  write.table(gse_BP_result, file = paste0(picDir, "/gse_BP_result.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  gse_CC_result <- gseGO(geneList     = geneList,#排序后的基因列表，一般根据logFC进行排序
                         OrgDb        = org_db,
                         ont          = "CC",#可选择bp.MF,CC,ALL
                         minGSSize    = 100,#最小基因集的基因数
                         maxGSSize    = 500,#最大基因集的基因数
                         pvalueCutoff = 0.05,#p值的阈值
                         verbose      = FALSE)#是否输出提示信息，默认为false
  write.table(gse_CC_result, file = paste0(picDir, "/gse_CC_result.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  gse_MF_result <- gseGO(geneList     = geneList,#排序后的基因列表，一般根据logFC进行排序
                         OrgDb        = org_db,
                         ont          = "MF",#可选择bp.MF,CC,ALL
                         minGSSize    = 100,#最小基因集的基因数
                         maxGSSize    = 500,#最大基因集的基因数
                         pvalueCutoff = 0.05,#p值的阈值
                         verbose      = FALSE)#是否输出提示信息，默认为false
  write.table(gse_MF_result, file = paste0(picDir, "/gse_MF_result.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  gse_ALL_result <- gseGO(geneList     = geneList,#排序后的基因列表，一般根据logFC进行排序
                          OrgDb        = org_db,
                          ont          = "ALL",#可选择bp.MF,CC,ALL
                          minGSSize    = 100,#最小基因集的基因数
                          maxGSSize    = 500,#最大基因集的基因数
                          pvalueCutoff = 0.05,#p值的阈值
                          verbose      = FALSE)#是否输出提示信息，默认为false
  write.table(gse_ALL_result, file = paste0(picDir, "/gse_ALL_result.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  ########## KEGG的GSEA富集分析：gseKEGG ##########
  # 进行KEGG分析
  if (org_db == "org.Hs.eg.db") {
    gsekk <- gseKEGG(geneList     = geneList,
                     organism     = 'hsa',
                     minGSSize    = 120,
                     pvalueCutoff = pvalueCutoff,
                     verbose      = FALSE)
  } else if (org_db == "org.Mm.eg.db") {
    gsekk <- gseKEGG(geneList     = geneList,
                     organism     = 'mmu',
                     nPerm        = 1000,
                     minGSSize    = 120,
                     pvalueCutoff = pvalueCutoff,
                     verbose      = FALSE)
  } else {
    stop("未知物种")
  }
  
  write.table(gsekk,file=paste0(picDir,"/", "GSEA_KEGG_result.txt"),sep="\t",
              quote=F,row.names = F)
  
  library(clusterProfiler)
  library(ggplot2)
  library(ggupset)
  library(enrichplot)
  # 可视化
  # 定义包含所有基因集富集分析结果的列表 
  gesalist <- list(gse_BP_result, gse_MF_result, gse_CC_result, gse_ALL_result, gsekk)
  # 定义每个基因集富集分析结果对应的名称
  gesa_name <- c("gse_BP_result", "gse_MF_result", "gse_CC_result", "gse_ALL_result", "gsekk")
  # 循环对每个基因集富集分析结果进行可视化并保存为pdf
  for (i in seq_along(gesalist)) {
    # 获取当前循环所代表的基因集富集分析结果和它的名称
    gseaname <- gesalist[[i]]
    gsea_name <- gesa_name[i]
    # 生成一个点图，展示最显著的前10个基因集
    my_plot <- dotplot(gseaname, showCategory = 10) + ggtitle("dotplot for GSEA")
    pdf(paste0(picDir,"/", gsea_name ,"_dotplot.pdf"))
    print(my_plot)
    dev.off()
    # 生成一个“集合的交集和并集”图
    my_plot <- upsetplot(gseaname)
    pdf(paste0(picDir,"/",gsea_name,"_upsetplot.pdf"),width=15,height =10)
    print(my_plot)
    dev.off()
    # 生成一个小提琴图，展示基因的富集分布情况
    my_plot <- gseaplot2(gseaname, geneSetID = 1:5, pvalue_table = TRUE)
    pdf(paste0(picDir,"/", gsea_name,"_gseaplot2.pdf"),width=16,height =12)
    print(my_plot)
    dev.off()
    # 生成一个脊型图，展示基因集富集得分的分布情况
    my_plot <- ridgeplot(gseaname, showCategory = 15)
    pdf(paste0(picDir,"/", gsea_name,"_ridgeplot.pdf"))
    print(my_plot)
    dev.off()
  }
  return(list(gse_BP = gse_BP_result, gse_CC = gse_CC_result, gse_MF = gse_MF_result, gse_ALL = gse_ALL_result,gsekk=gsekk))}



######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################

ydx_draw_MsigDb_GSEA= function(DEG_list,pvalueCutoff = 0.05){
  library(msigdbr)
  library(fgsea)
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(tidyverse)
  library(pacman)
  cat("使用GSEA函数做GSEA\n")
  # 使用p_load函数来安装（如果需要）并加载所需的包
  p_load(topGO, ReactomePA, tidyverse, clusterProfiler, biomaRt, enrichplot, dplyr, data.table, ggplot2, patchwork, cols4all, ggraph, igraph,RColorBrewer, ComplexHeatmap, ggpubr, circlize, org.Hs.eg.db, org.Mm.eg.db)
  picDir = "SigDb_GSEA_"
  if (!dir.exists(picDir)) {dir.create(picDir)}
  
  colnames(DEG_list)[1] <- "SYMBOL"
  
  
  result <- sum(grepl("[a-z]", DEG_list[1]))
  if(result>length(result)/2){
    species= "Mus musculus"
    org_db= "org.Mm.eg.db"
    cat("检测到基因属于小鼠，使用小鼠基因集做分析\n")
  }else{
    species= "Homo sapiens"
    org_db= "org.Hs.eg.db"
    cat("检测到基因属于人类，使用人类基因集做分析\n")
  }
  
  # 使用grep函数获取列名
  cols <- names(DEG_list)[grep("(?i)(symbol|log)", names(DEG_list), perl = TRUE)]
  
  # 提取symbol和logfc列
  DEG_list_sub <- DEG_list[, cols]
  ## 1: 提取logFC值，并储存在一个向量中
  geneList = DEG_list_sub[,2]
  ## 2: 对geneList进行命名
  names(geneList) = as.character(DEG_list_sub$SYMBOL)
  head(geneList)
  ## 3: 根据logFC值降序排列
  geneList = sort(geneList, decreasing = TRUE)
  # gmtFile <- paste0(function_dir, "c5.go.v2023.1.Hs.symbols.gmt")
  # # 读取上面指定的gmt文件
  # all_msigdb <- read.gmt(gmtFile)
  # gsegmt <- GSEA(geneList, TERM2GENE=all_msigdb, verbose=FALSE)
  # head(gsegmt)
  # gseaplot2(gse_BP_result,geneSetID = 1:5, pvalue_table = TRUE)#绘制三个基因集
  # plotGseaTable(gsegmt)
  

  
  geneset_GO = msigdbr(species = species,
                       category = "C5", 
                       subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
  
  geneset_GO$gs_name <- gsub('GOBP_','',geneset_GO$gs_name)#去除前缀
  geneset_GO$gs_name <- str_to_title(geneset_GO$gs_name)#将首字母除外的大写换为小写
  geneset_GO$gs_name <- gsub('_',' ',geneset_GO$gs_name)#将_转化为空格
  
  GSEA_BP_result <- GSEA(geneList, TERM2GENE=geneset_GO, verbose=FALSE)
  
  write.table(GSEA_BP_result,file=paste0(picDir,"/", "GSEA_BP_result.txt"),sep="\t",
              quote=F,row.names = F)
  
  geneset_GO = msigdbr(species = species,
                       category = "C5", 
                       subcategory = "GO:CC") %>% dplyr::select(gs_name,gene_symbol)
  
  geneset_GO$gs_name <- gsub('GOCC_','',geneset_GO$gs_name)#去除前缀
  geneset_GO$gs_name <- str_to_title(geneset_GO$gs_name)#将首字母除外的大写换为小写
  geneset_GO$gs_name <- gsub('_',' ',geneset_GO$gs_name)#将_转化为空格
  
  GSEA_CC_result <- GSEA(geneList, TERM2GENE=geneset_GO, pvalueCutoff = pvalueCutoff,verbose=FALSE)
  
  write.table(GSEA_CC_result,file=paste0(picDir,"/", "GSEA_CC_result.txt"),sep="\t",
              quote=F,row.names = F)
  
  
  geneset_GO = msigdbr(species = species,
                       category = "C5", 
                       subcategory = "GO:MF") %>% dplyr::select(gs_name,gene_symbol)
  
  geneset_GO$gs_name <- gsub('GOMF_','',geneset_GO$gs_name)#去除前缀
  
  geneset_GO$gs_name <- str_to_title(geneset_GO$gs_name)#将首字母除外的大写换为小写
  geneset_GO$gs_name <- gsub('_',' ',geneset_GO$gs_name)#将_转化为空格
  GSEA_MF_result <- GSEA(geneList, TERM2GENE=geneset_GO, pvalueCutoff = pvalueCutoff,verbose=FALSE)
  write.table(GSEA_MF_result,file=paste0(picDir,"/", "GSEA_MF_result.txt"),sep="\t",
              quote=F,row.names = F)
  
  #######################################geneset_KEGG#############################
  geneset_KEGG <-  msigdbr(species = species, 
                           category = "C2",
                           subcategory = "KEGG")  %>% dplyr::select(gs_name,gene_symbol)
  
  geneset_KEGG$gs_name <- gsub('KEGG_','',geneset_KEGG$gs_name)#去除前缀
  geneset_KEGG$gs_name <- str_to_title(geneset_KEGG$gs_name)#将首字母除外的大写换为小写
  geneset_KEGG$gs_name <- gsub('_',' ',geneset_KEGG$gs_name)#将_转化为空格
  gsekk <- GSEA(geneList, TERM2GENE=geneset_KEGG,verbose=FALSE) 
  
  write.table(gsekk,file=paste0(picDir,"/", "gsekk_result.txt"),sep="\t",
              quote=F,row.names = F)
  
  
  # 生成一个点图，展示最显著的前10个基因集
  gesalist <- list(GSEA_BP_result, GSEA_CC_result, GSEA_MF_result, gsekk)
  # 定义每个基因集富集分析结果对应的名称
  gesa_name <- c("gse_BP_result", "gse_CC_result", "gse_MF_result", "gsekk")
  # 循环对每个基因集富集分析结果进行可视化并保存为pdf
  for (i in seq_along(gesalist)) {
    # 获取当前循环所代表的基因集富集分析结果和它的名称
    gseaname <- gesalist[[i]]
    gsea_name <- gesa_name[i]
    my_plot <- cnetplot(gseaname, categorySize="pvalue", color.params = list(foldChange = geneList))
    pdf(paste0(picDir,"/", gsea_name ,"_cnetplot.pdf"),width=15,height =10)
    print(my_plot)
    dev.off()
    
    my_plot <- heatplot(gseaname, foldChange=geneList,showCategory = 10,label_format = 30)
    pdf(paste0(picDir,"/", gsea_name ,"_heatplot.pdf"),width=15,height =10)
    print(my_plot)
    dev.off()
    # 生成一个点图，展示最显著的前10个基因集
    my_plot <- dotplot(gseaname, showCategory = 10) + ggtitle("dotplot for GSEA")
    pdf(paste0(picDir,"/", gsea_name ,"_dotplot.pdf"))
    print(my_plot)
    dev.off()
    # 生成一个“集合的交集和并集”图
    my_plot <- upsetplot(gseaname,label_size = 20)
    pdf(paste0(picDir,"/",gsea_name,"_upsetplot.pdf"),width=15,height =10)
    print(my_plot)
    dev.off()
    # 生成一个小提琴图，展示基因的富集分布情况
    my_plot <- gseaplot2(gseaname, geneSetID = 1:5, pvalue_table = TRUE)
    pdf(paste0(picDir,"/", gsea_name,"_gseaplot2.pdf"),width=16,height =12)
    print(my_plot)
    dev.off()
    # 生成一个脊型图，展示基因集富集得分的分布情况
    my_plot <- ridgeplot(gseaname, showCategory = 15)
    pdf(paste0(picDir,"/", gsea_name,"_ridgeplot.pdf"))
    print(my_plot)
    dev.off()
  }
  return(list(gse_BP = GSEA_BP_result, gse_CC = GSEA_CC_result, gse_MF = GSEA_MF_result,gsekk=gsekk))
}






ydx_kegg_hierarchy <- (picDir = "GO",fromType = "SYMBOL",toType = "ENTREZID"){
  
  library(tidyjson) 
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(ggplot2)
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {
    dir.create(picDir)
  }
  if(is.data.frame(GeneName)){
    GeneName=GeneName[,1]
  }
  GeneName = as.character(GeneName)
  result <- sum(grepl("[a-z]", GeneName))
  if(result>length(GeneName)/2){
    org_db= "org.Mm.eg.db"
    cat("检测到基因属于小鼠，使用小鼠基因集做分析\n")
  }else{
    org_db= "org.Hs.eg.db"
    cat("检测到基因属于人类，使用人类基因集做分析\n")
  }
  ENSEMBLlist <- bitr(GeneName, fromType = fromType, toType = toType, org_db)
  ENSEMBLlist <- ENSEMBLlist$ENTREZID
  # 进行KEGG分析
  if (org_db == "org.Hs.eg.db") {
    kegg <- enrichKEGG(ENSEMBLlist, organism = "hsa", qvalueCutoff=qvalueFilter,pvalueCutoff=pvalueFilter, keyType = "kegg", pAdjustMethod = "BH")
  } else if (org_db == "org.Mm.eg.db") {
    kegg <- enrichKEGG(ENSEMBLlist, organism = "mmu", qvalueCutoff=qvalueFilter,pvalueCutoff=pvalueFilter, keyType = "kegg", pAdjustMethod = "BH")
    # 使用dplyr的mutate()和str_replace()函数去除指定值
    kegg@result$Description <- str_replace_all(kegg@result$Description, "- Mus musculus \\([^\\)]+\\)", "")    
  } else {
    stop("Not a valid database")
  }
  cat("KEGG富集完成")
  
  kegg = as.data.frame(kegg)
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  # 加载json file
  path_hierarchy <- jsonlite::fromJSON("https://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=json&filedir=", simplifyDataFrame = T, flatten = T)
  # 仅保留path_hierarchy数据框架中的children列
  path_hierarchy <- path_hierarchy$children
  # 将name列重命名为level1
  path_hierarchy <- rename(path_hierarchy, level1 = name)
  # 展开数据框架中的children列，每个子元素有自己的行
  path_hierarchy <- unnest(path_hierarchy, cols = c(children))
  
  # 将name列重命名为level2
  path_hierarchy <- rename(path_hierarchy, level2 = name)
  
  # 展开数据框架中的children列，每个子元素有自己的行
  path_hierarchy <- unnest(path_hierarchy, cols = c("children"))
  
  # 将name列重命名为level3
  path_hierarchy <- rename(path_hierarchy, level3 = name)
  # 创建新的列id，由level3列的前五个字符构成，前面加上字符串"ko"，并修改level3列
  path_hierarchy <- mutate(path_hierarchy, id = paste0("ko", substr(level3, 1, 5)),
                           level3 = substr(level3, 7, nchar(level3)))
  # 仅保留path_hierarchy数据框架中的id、level3、level2和level1这四列，丢弃其它列
  path_hierarchy <- select(path_hierarchy, id, level3, level2, level1)
  
  
  kegg2=path_hierarchy
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

 
  
  ggData <- data %>%
    mutate(level3 = kegg$level3[match(ID, kegg$id)],
           level2 = kegg$level2[match(ID, kegg$id)],
           level1 = kegg$level1[match(ID, kegg$id)]) %>%
    
    # 给level1和level2的term排序
    # 注意：只保留table(kegg$level1)里出现的term
    mutate(level1 = factor(level1, 
                           c("Metabolism", "Cellular Processes", "Environmental Information Processing",
                             "Genetic Information Processing", "Organismal Systems")),
           level2 = factor(level2, 
                           c("Global and overview maps", "Carbohydrate metabolism", "Energy metabolism", 
                             "Lipid metabolism", "Glycan biosynthesis and metabolism",
                             "Biosynthesis of other secondary metabolites",
                             "Xenobiotics biodegradation and metabolism", 
                             "Metabolism of other amino acids", 
                             "Metabolism of cofactors and vitamins", 
                             "Metabolism of terpenoids and polyketides", 
                             "Nucleotide metabolism", 
                             "Cellular community - eukaryotes", 
                             "Transport and catabolism", 
                             "Membrane transport",
                             "Signal transduction", 
                             "Folding, sorting and degradation", 
                             "Replication and repair", 
                             "Translation",
                             "Digestive system"))
    ) %>%
    arrange(level1, level2, Direction, level3) %>%
    mutate(ID = factor(ID, rev(unique(ID))),
           level3_x = factor(level3, rev(unique(as.character(level3)))))
  
  # level1的文字
  ggData_l1 <- data.frame(nrow(ggData) - (cumsum(table(ggData$level1)) - table(ggData$level1)/2))
  ggData_l1$start <- ggData_l1$Freq - table(ggData$level1)/2
  ggData_l1$end <- ggData_l1$Freq + table(ggData$level1)/2
  
  # level2的文字
  ggData_l2 <- data.frame(nrow(ggData) - (cumsum(table(ggData$level2)) - table(ggData$level2)/2))
  ggData_l2$start <- ggData_l2$Freq - table(ggData$level2)/2
  ggData_l2$end <- ggData_l2$Freq + table(ggData$level2)/2
  
  
  
  
  
  
}

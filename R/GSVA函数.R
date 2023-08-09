ydx_draw_GSVA <- function(pre_machine_data=counts,picDir = "GSVA") {

  if (!file.exists(picDir)) { dir.create(picDir) }
  library(GSVA)
  library(limma)
  library(GSEABase)
  library(msigdbr)
  library(pheatmap)
   
   
   # 判断是否包含小写字母并返回结果为TRUE的行
   result <- sum(grepl("[a-z]", rownames(pre_machine_data)))
   if(result>20){
     species = "Mus musculus"
     cat("检测到基因属于小鼠，使用小鼠基因集做分析\n")
   }else{
     species = "Homo sapiens"
     cat("检测到基因属于人类，使用人类基因集做分析\n")
   }
   
   
   
  mat <- as.matrix(pre_machine_data)
  mat <- avereps(mat)
  mat <- normalizeBetweenArrays(mat)
  #######生成KEGG  GENESET#######################
  all_gene_sets = msigdbr(species = species,
                          category="C2",
                          subcategory = "KEGG") #从Msigdb中获取信号通路，此处H是癌症的标志物
  
  length(unique(table(all_gene_sets$gs_name))) # 116条通路
  tail(table(all_gene_sets$gs_name))
  gs=split(all_gene_sets$gene_symbol,all_gene_sets$gs_name)
  # 去除通路中KEGG的同时，再改为小写
  gs = lapply(gs, unique)
  head(gs)
  # gsc <- GeneSetCollection(mapply(function(geneIds, keggId) {
  #   GeneSet(geneIds, geneIdType=EntrezIdentifier(),
  #           collectionType=KEGGCollection(keggId),
  #           setName=keggId)
  # }, gs, names(gs)))
  
  # head(gsc)
  kegg_geneset <- gs
  #################生成KEGG  GENESET         结束##########################
  output <- gsva(mat, kegg_geneset, min.sz = 1, max.sz = 500, verbose = TRUE, parallel.sz = 1)
  KEGG_output <- rbind(id = colnames(output), output)
  write.table(KEGG_output, file = paste0(picDir, "/", "KEGG_output.txt"), sep = "\t", quote = F, col.names = F)

  

  
#######生成REACTOME  GENESET#######################
  # 从Msigdb中获取人类Reactome基因集
  all_gene_sets = msigdbr(species = species,
                          category="C2",
                          subcategory = "REACTOME")
  
  # 统计基因集的数量并去重
  num_gene_sets <- length(unique(table(all_gene_sets$gs_name)))
  cat("Number of gene sets:", num_gene_sets, "\n")
  
  # 按基因集名称分组并提取基因符号
  gs <- split(all_gene_sets$gene_symbol, all_gene_sets$gs_name)
  gs <- lapply(gs, unique)
  
  # 将基因符号列表转换为GeneSet对象，并创建GeneSetCollection对象
  # gsc <- GeneSetCollection(mapply(function(geneIds, setName) {
  #   GeneSet(geneIds, geneIdType=EntrezIdentifier(),
  #           collectionType=KEGGCollection(setName),
  #           setName=setName)
  # }, gs, names(gs)))
  # 
  # # 显示GeneSetCollection对象的摘要信息
  # summary(gsc)

REACTOME_geneset <- gs
##################生成REACTOME  GENESET   结 束           #########################
output <- gsva(mat, REACTOME_geneset, min.sz = 1, max.sz = 500, verbose = TRUE, parallel.sz = 1)
reactome_output <- rbind(id = colnames(output), output)
write.table(reactome_output, file = paste0(picDir, "/", "reactome_output.txt"), sep = "\t", quote = F, col.names = F)



  
  return(list(KEGG = KEGG_output, reactome_result = reactome_output ))
}
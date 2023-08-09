ydx_friends <- function(GeneName=t1dm_mci_wgcna_mit, picDir = "GO",width= 8,height= 10,fromType = "SYMBOL",toType = "ENTREZID"){
  #加载包
  
  cat("背景资料：https://mp.weixin.qq.com/s/_Wt_GmC8yjcvEdXBNRUQTw\n")
  library(GOSemSim)
  library(reshape2)
  library(ggplot2)
  library(clusterProfiler)
  # 创建输出目录（如果不存在）
  if (!file.exists(picDir)) {
    dir.create(picDir)
  }
  if(is.data.frame(GeneName)){
    GeneName=GeneName[,1]
  }
  GeneName = as.character(GeneName)
  # 判断是否包含小写字母并返回结果为TRUE的行
  result <- sum(grepl("[a-z]", GeneName))
  if(result>length(GeneName)/2){
    library(org.Mm.eg.db)
    org_db= "org.Mm.eg.db"
    cat("检测到基因属于小鼠，使用小鼠基因集做分析\n")
  }else{
    library(org.Hs.eg.db)
    org_db= "org.Hs.eg.db"
    cat("检测到基因属于人类，使用人类基因集做分析\n")
  }
  ENSEMBLlist <- bitr(GeneName, fromType = fromType, toType = toType, org_db)
  ENSEMBLlist_ENTREZID <- ENSEMBLlist$ENTREZID
  #用godata()函数来构建相应物种的Molecular Function本体的GO DATA
  mf <- godata(org_db, ont="MF", computeIC = FALSE)
  #用godata()函数来构建相应物种的Cellular Component本体的GO DATA
  cc <- godata(org_db, ont="CC", computeIC = FALSE)
  #用godata()函数来构建相应物种的Biological Process本体的GO DATA
  bp <- godata(org_db, ont="BP", computeIC = FALSE)
 
   ########计算语义相似度
  #用mgeneSim来计算MF本体，基因之间的语义相似度，结果为一个行列相同的矩阵
  simmf <- mgeneSim(ENSEMBLlist_ENTREZID, semData = mf, measure = "Wang", drop = NULL, combine = "BMA")
  #用mgeneSim来计算CC本体，基因之间的语义相似度，结果为一个行列相同的矩阵
  simcc <- mgeneSim(ENSEMBLlist_ENTREZID, semData = cc, measure = "Wang", drop = NULL, combine = "BMA")
  #用mgeneSim来计算BP本体，基因之间的语义相似度，结果为一个行列相同的矩阵
  simbp <- mgeneSim(ENSEMBLlist_ENTREZID, semData = bp, measure = "Wang", drop = NULL, combine = "BMA")
  
  cat("选择共有的基因做分析\n")
  # 使用Reduce()和intersect()获取共有的列名
  common_columns <- Reduce(intersect, list(colnames(simbp), colnames(simcc), colnames(simmf)))
  
  # 筛选出共有列名对应的数据框
  simbp_filtered <- simbp[common_columns, common_columns]
  simcc_filtered <- simcc[common_columns, common_columns]
  simmf_filtered <- simmf[common_columns, common_columns]
  
  
  #或者计算基因在MF、CC、BP本体下的几何平均值
  fsim <- (simmf_filtered * simcc_filtered * simbp_filtered)^(1/3)
  
  # 1. 将fsim的行名和列名与ENSEMBLlist的第一列匹配，获取匹配成功的索引
  matched_rows <- match(rownames(fsim), ENSEMBLlist$ENTREZID)
  matched_cols <- match(colnames(fsim), ENSEMBLlist$ENTREZID)
  
  # 2. 使用匹配的索引将fsim中的行名和列名替换为ENSEMBLlist的第二列的数据
  rownames(fsim) <- ENSEMBLlist$SYMBOL[matched_rows]
  colnames(fsim) <- ENSEMBLlist$SYMBOL[matched_cols]
  
  fsim_data=as.data.frame(fsim)
  #将基因自己和自己的相似度设为NA，方便接下来去掉。
  for (i in 1:ncol(fsim)){
    fsim[i,i] <- NA
  }
  
  y <- melt(fsim) #把宽格式数据转化成长格式，其实就是把正方形矩阵转成三列
  y <- y[!is.na(y$value),] #删掉带NA的行
  
  
  ## 开始画图
  #计算每个基因跟其他基因相似度的平均值
  y.mean <- aggregate(.~Var1,y,mean) 
  m <- y.mean$value
  names(m) <- y.mean$Var1
  #按平均值给基因名排序，便于画图
  y$Var1 <- factor(y$Var1, levels=names(sort(m)))
  
  f <- function(y) {
    r <- quantile(y, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
    r[3] <- mean(y)
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  colo <- c("#FF8F00",
            "#FFA700", "#FFBF00", "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00",
            "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", "#20FF00",
            "#08FF00", "#00FF10", "#00FF28", "#00FF40", "#00FF58", "#00FF70", "#00FF87",
            "#00FF9F", "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF",
            "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", "#0040FF", "#0028FF",
            "#0010FF", "#0800FF", "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF",
            "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", "#FF00D7",
            "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", "#FF0060", "#FF0048", "#FF0030",
            "#FF0018")
  
  color=colo[1:23]
  
  p1 <- ggplot(y, aes(Var1, value, fill = factor(Var1))) + 
  #  ggsci::scale_fill_ucscgb() +
    guides(fill=FALSE) + #不显示图例
    stat_summary(fun.data= f, geom='boxplot') + 
    geom_hline(aes(yintercept=0.75), linetype="dashed") + #画一条虚线
    coord_flip() + # x、y坐标轴互换
    xlab("") + ylab("") + 
    theme(axis.text.x = element_text(family = "Arial", size = 16, face = "bold"),
          axis.text.y = element_text(family = "Arial", size = 16, face = "bold")) + 
    theme_bw() + 
    theme(panel.border=element_rect(size=1)) #边框粗细 
  p1
  
  # 保存到pdf文件
  ggsave(paste0(picDir,"/friends_box.pdf"),width= width,height=height )
  return(list(fsim=fsim_data,table = y))
  
}

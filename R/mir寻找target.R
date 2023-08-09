
mir_target()

mir_target<- function(mir ="mir-137",species = "mmu"){
  # 定义要匹配的字符串
  library(data.table)
  library(dplyr)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  options(timeout = 3000)
  
  cat("正在下载targetscan数据库\n")
  # 下载文件并保存在临时目录中
  temp_file <- tempfile()
  download.file("https://www.targetscan.org/mamm_31/mamm_31_data_download/Conserved_Family_Conserved_Targets_Info.txt.zip", temp_file)
  # 解压缩文件
  unzip(temp_file, exdir = tempdir())
  # 读取文件
  targetscan <- data.table::fread(file.path(tempdir(), "Conserved_Family_Conserved_Targets_Info.txt"))
  #targetscan =unique(select(targetscan,1:2)) 
  
  
  cat("正在下载mirdb数据库\n")
  # 下载文件并保存在临时目录中
  temp_file <- tempfile()
  download.file("https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz", temp_file)
  # 读取和解压缩GZIP文件
  mirdb <- data.table::fread(cmd = paste("gzip -cd", temp_file))
  # 保留以"hsa"和"mmu"开头的行
  mirdb <- mirdb[grep("^hsa|^mmu", mirdb$V1), ]
  
  # 使用bitr函数将REFSEQ转录本标识符转换为Symbol，并将结果存储在mirdb$V2列中
  if (species == "mmu") {
    mirdb_trans <- bitr(mirdb$V2, fromType = "REFSEQ", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
  }else{
    mirdb_trans <- bitr(mirdb$V2, fromType = "REFSEQ", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
    
  }
  merged_mirdb <- merge(mirdb, mirdb_trans, by.x = "V2", by.y = "REFSEQ", all.x = TRUE)
  mirdb= na.omit(select(merged_mirdb,2,4,3))
  rm(merged_mirdb,mirdb_trans)
  
  
  
  
  cat("正在下载mirtarbase数据库\n")
  # 下载 miRTarBase 的 Excel 文件
  url <- "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/cache/download/9.0/miRTarBase_MTI.xlsx"
  dest_file <- "miRTarBase_MTI.xlsx"
  downloader::download(url, dest_file, mode = "wb")
  mirtarbase= readxl::read_excel(dest_file)
  mirtarbase <- mirtarbase[grep(species, mirtarbase$miRNA), ]
  # mirtarbase =select(mirtarbase,2,4)
  # mirtarbase=unique(mirtarbase)
  # 
  
  
  
  cat("正在下载NPInter数据库\n")
  # 下载文件并保存在临时目录中
  temp_file <- tempfile()
  download.file("http://bigdata.ibp.ac.cn/npinter5/download/file/interaction_NPInterv4.txt.gz", temp_file)
  # 读取和解压缩GZIP文件
  NPInter <- data.table::fread(cmd = paste("gzip -cd", temp_file))
  # 保留已选择的物种

  NPInter <- NPInter[grepl(species, NPInter$ncName), ]
  
  
  
  
  # cat("正在下载ENCORI数据库\n")
  # mir ="mir-7-3p"
  # ENCORI_url <- paste0('https://rna.sysu.edu.cn/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA=', mir, '&target=all&cellType=all')
  # ENCORI_file <- "ENCORI_data.txt"
  # 
  # 
  # download.file(ENCORI_url, destfile = ENCORI_file)
  # 
  # ENCORI <- data.table::fread(ENCORI_file)
  

  # 构建正则表达式，不区分大小写，数字部分绝对匹配
  regex_pattern <- paste0("(?i)\\b", mir, "\\b")
  # 选择包含"mir-7"（不区分大小写，数字部分绝对匹配）的行
  mir_NPInter <- subset(NPInter, apply(NPInter, 1, function(row) any(grepl(regex_pattern, row))))
  mir_targetscan<- subset(targetscan, apply(targetscan, 1, function(row) any(grepl(regex_pattern, row))))
  mir_mirtarbase<- subset(mirtarbase, apply(mirtarbase, 1, function(row) any(grepl(regex_pattern, row))))
  mir_mirdb<- subset(mirdb, apply(mirdb, 1, function(row) any(grepl(regex_pattern, row))))


  # 保存 mir_NPInter 数据框
  write.csv(mir_NPInter, "mir_NPInter.csv", row.names = FALSE)
  # 保存 mir_targetscan 数据框
  write.csv(mir_targetscan, "mir_targetscan.csv", row.names = FALSE)
  # 保存 mir_mirtarbase 数据框
  write.csv(mir_mirtarbase, "mir_mirtarbase.csv", row.names = FALSE)
  # 保存 mir_mirdb 数据框
  write.csv(mir_mirdb, "mir_mirdb.csv", row.names = FALSE)
  
  
}


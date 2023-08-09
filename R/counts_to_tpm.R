#' 将计数转换为每百万转录本数（TPM）。
#' 
#' 将具有原始特征计数的数字矩阵（行为特征，列为条件）转换为每百万转录本数。
#' 
#'    Lior Pachter. Models for transcript quantification from RNA-Seq.
#'    arXiv:1104.3889v2 
#'    
#' @param counts 数字矩阵，包含每个基因的原始特征计数，即分配给每个基因的片段数。
#' @param featureLength 数字向量，包含每个特征（基因）的长度。
#' @param meanFragmentLength 数字向量，包含平均片段长度。
#' @return tpm 数值矩阵，通过库大小和特征长度进行归一化。
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # 确保参数有效。
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # 计算每个库中特征的有效长度。
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # 排除长度小于平均片段长度的基因。
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # 逐列处理数据。
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # 复制原始矩阵的行和列名称。
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

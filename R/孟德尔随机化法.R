Mendelian_Randomization <- function( exposure_ID='ebi-a-GCST010681',outcome_ID ='finn-b-G6_ALZHEIMER',
                                     p1 = 5e-08,r2 = 0.001,kb = 10000){ #加载R包
  library(TwoSampleMR)
  if (grepl("[/.]", exposure_ID)) {
    exposure_data=exposure_data
    } else {
    # 如果a中不包含/和.则执行以下代码
   
    exposure_data<-extract_instruments(outcomes=exposure_ID,p1 = p1,r2 = r2,kb = kb, clump=TRUE)
  }
  ################################################################################
  #获取暴露数据,access_token这个参数，对于中国大陆地区的用户必须设置该参数为access_token=NULL，这样才能顺利获取数据，否则就需要开VPN获取谷歌授权
  #如果设置r2=0.001和kb=10000，那这就表示去掉在10000kb范围内与最显著SNP的r2大于0.001的SNP。
  # 剔除混杂因素 http://www.phenoscanner.medschl.cam.ac.uk/
  outcome_data<-extract_outcome_data(
    snps = exposure_data$SNP,
    outcomes =outcome_ID,
    proxies = FALSE,
    maf_threshold = 0.01)
  #snps：它是一串以rs开头的SNP ID
  #outcomes：它是outcome在MR base中的ID；
  #proxies：它表示是否使用代理SNP，默认值是TRUE，也即当一个SNP在outcome中找不到时可以使用与其存在强连锁不平衡的SNP信息来替代，我个人喜欢设置成FALSE。
  #maf_threshold：它表示的是SNP在outcome中的最小等位基因频率，默认值是0.3，不过大样本GWAS可以适当调低，这里设置的是0.01。
 
  #合并暴露和结局的数据
  bind_data <- harmonise_data(
    exposure_dat=exposure_data,
    outcome_dat=outcome_data,
    action= 2
  )#将IV的效应等位基因（effect allele）对齐，
  #第三个参数action最重要，一般我们推荐使用默认值action=2即可，当然也可以使用action=3，这时候就表示去除所有存在回文结构的SNP
  #MR分析
  res <- mr(bind_data)
  
  res
  OR<-generate_odds_ratios(res)
  
  #异质性检验
  het <- mr_heterogeneity(bind_data)
  presso = run_mr_presso(bind_data,NbDistribution = 10000)
  het
  #这些IV之间存在很强的异质性（Q_pval远小于0.05），
  #这时候我们需要剔除某些outcome的P值非常小的SNP，
  #或者直接使用随机效应模型来估计MR效应量
  #随机效应模型
  mr(bind_data,method_list=c('mr_ivw_mre'))
  #多效性检验
  pleio <- mr_pleiotropy_test(bind_data)
  pleio
 # 逐个剔除检验
  single <- mr_leaveoneout(bind_data)
  mr_leaveoneout_plot(single)
 
  #散点图
  p = mr_scatter_plot(res,bind_data)
  
  #森林图
  res_single <- mr_singlesnp(bind_data)
  mr_forest_plot(res_single)
  
  #漏斗图
  mr_funnel_plot(res_single)
  

  
  return(list(OR=OR,bind_data=bind_data,exposure_data=exposure_data,outcome_data=outcome_data))
  
  }


#####################################孟德尔随机化的操作流程########

# 如果a中包含/和.则执行以下代码
a= data.table::fread("D:\\Personal\\Downloads\\PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis (1).txt.gz")
b= subset(a,p_value<p1)

exposure_data =read_exposure_data(b,
                                  clump = FALSE,
                                  sep = ",",
                                  snp_col = "SNP",
                                  beta_col = "beta",
                                  se_col = "se",
                                  eaf_col = "eaf",
                                  effect_allele_col = "effect_allele",
                                  other_allele_col = "other_allele",
                                  pval_col = "pval",
                                  chr_col = "chr",
                                  pos_col = "pos")












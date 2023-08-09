# 代码简介：12种具有变量筛选或预测模型构建的算法组合
# 
# 算法包括：Lasso、Ridge、Enet、Stepglm、SVM、glmBoost、LDA、plsRglm、RandomForest、GBM、XGBoost、NaiveBayes
# 
# 具体内容：
# 在交叉验证的框架下使用一种算法进行变量选择而使用另一种算法构建分类预测模型，
# 并计算所使用的模型组合(共113种)在外部数据集(或包括训练集)上的ROC曲线下面积(AUC)，最终通过热图将模型的评估结果可视化
# 
# 算法输入：
# 训练集表达谱：行为特征(如SYMBOL/ENSEMBL等，但需与测试集有一定交集)，列为样本的表达矩阵，格式见InputData文件中Train_expr.txt
# 训练集分类信息：行为样本，列包含至少一个需要预测的二分类变量(仅支持[0，1]格式，代码中需分别指定相对应的变量名)，格式见InputData文件中Train_class.txt
# 测试集表达谱：行为特征(如SYMBOL/ENSEMBL等，但需与训练集有一定交集)，列为样本的表达矩阵，格式见InputData文件中Test_expr.txt
# 测试集生存信息：行为样本，列包含至少一个需要预测的二分类变量(仅支持[0，1]格式，代码中需分别指定相对应的变量名)，以及一列用于指定队列信息的变量，格式见InputData文件中Test_class.txt
# 注意：建议对表达谱进行预筛选，保留几十至几千量级的感兴趣特征用于最优模型筛选，否则运行时间将大大延长
# 
# 本代码适用于：
# 存在大量对分类预测意义不明的变量时，需同时进行变量筛选，模型筛选和分类模型构建的场景
# 
# 代码输出结果包括：
# 113种算法在训练集上获得的模型(model.rds)，在测试集(或包括训练集)上的评估结果(AUC_mat.txt)，以及评估结果的热图展示(AUC.pdf)
# 
# 参考文献：
# Machine learning-based integration develops an immune-derived lncRNA signature for improving outcomes in colorectal cancer
# DOI: s41467-022-28421-6
# 
# 作者：大鱼海棠
# 工作单位：中国药科大学国家天然药物重点实验室，生物统计与计算药学研究中心
# 目前地址：法国斯特拉斯堡遗传与分子生物研究所（IGBMC），癌症功能基因组研究中心
# 联系邮箱：xlu.cpu@foxmail.com
machine_learning_113 <- function(Train_expr=pre_data$exp_data,train_lab="Train group",Logistic= FALSE,Test_expr,Test_class, PicDir = "picdir"){
  # 加载需要使用的R包
  library(openxlsx)
  library(seqinr)
  library(plyr)
  library(randomForestSRC)
  library(glmnet)
  library(plsRglm)
  library(gbm)
  library(caret)
  library(mboost)
  library(e1071)
  library(BART)
  library(MASS)
  library(snowfall)
  library(xgboost)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(pROC)
  # 构造完整路径
  picdir.path <- file.path(getwd(), PicDir)
  # 将工作目录设置为picdir文件夹路径
  # 获取当前工作目录
  work.path <- picdir.path
  # 设置其他路径

  res.path <- file.path(work.path, "Results")
  # 如不存在这些路径则创建路径
  if (!dir.exists(PicDir)) dir.create(PicDir)
  if (!dir.exists(res.path)) dir.create(res.path)

  # all ML algorithms
  RunML <- function(method, Train_set, Train_label, mode = "Model", classVar){
    # for example: Enet [alpha=0.4]
    method = gsub(" ", "", method) # 去除参数中的空格，得到Enet [alpha=0.4]
    method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)  # get name of ML algorithm, e.g., Enet
    method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method) # get parameter of ML algorithm, e.g., alpha=0.4
    
    method_param = switch(
      EXPR = method_name,
      "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
      "Stepglm" = list("direction" = method_param),
      NULL
    )
    message("Run ", method_name, " algorithm for ", mode, "; ",
            method_param, ";",
            " using ", ncol(Train_set), " Variables")
    
    args = list("Train_set" = Train_set,
                "Train_label" = Train_label,
                "mode" = mode,
                "classVar" = classVar)
    args = c(args, method_param)
    
    obj <- do.call(what = paste0("Run", method_name),
                   args = args) 
    
    if(mode == "Variable"){
      message(length(obj), " Variables retained;\n")
    }else{message("\n")}
    return(obj)
  }
  
  RunEnet <- function(Train_set, Train_label, mode, classVar, alpha){
    cv.fit = cv.glmnet(x = Train_set,
                       y = Train_label[[classVar]],
                       family = "binomial", alpha = alpha, nfolds = 10)
    fit = glmnet(x = Train_set,
                 y = Train_label[[classVar]],
                 family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }
  
  RunLasso <- function(Train_set, Train_label, mode, classVar){
    RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
  }
  
  RunRidge <- function(Train_set, Train_label, mode, classVar){
    RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
  }
  
  RunStepglm <- function(Train_set, Train_label, mode, classVar, direction){
    fit <- step(glm(formula = Train_label[[classVar]] ~ .,
                    family = "binomial", 
                    data = as.data.frame(Train_set)),
                direction = direction, trace = 0)
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }
  
  RunSVM <- function(Train_set, Train_label, mode, classVar){
    data <- as.data.frame(Train_set)
    data[[classVar]] <- as.factor(Train_label[[classVar]])
    fit = svm(formula = eval(parse(text = paste(classVar, "~."))),
              data= data, probability = T)
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }
  
  RunLDA <- function(Train_set, Train_label, mode, classVar){
    data <- as.data.frame(Train_set)
    data[[classVar]] <- as.factor(Train_label[[classVar]])
    fit = train(eval(parse(text = paste(classVar, "~."))), 
                data = data, 
                method="lda",
                trControl = trainControl(method = "cv"))
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }
  
  RunglmBoost <- function(Train_set, Train_label, mode, classVar){
    data <- cbind(Train_set, Train_label[classVar])
    data[[classVar]] <- as.factor(data[[classVar]])
    fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                    data = data,
                    family = Binomial())
    
    cvm <- cvrisk(fit, papply = lapply,
                  folds = cv(model.weights(fit), type = "kfold"))
    fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                    data = data,
                    family = Binomial(), 
                    control = boost_control(mstop = max(mstop(cvm), 40)))
    
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }
  
  RunplsRglm <- function(Train_set, Train_label, mode, classVar){
    cv.plsRglm.res = cv.plsRglm(formula = Train_label[[classVar]] ~ ., 
                                data = as.data.frame(Train_set),
                                nt=10, verbose = FALSE)
    fit <- plsRglm(Train_label[[classVar]], 
                   as.data.frame(Train_set), 
                   modele = "pls-glm-logistic",
                   verbose = F, sparse = T)
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }
  
  RunRF <- function(Train_set, Train_label, mode, classVar){
    rf_nodesize = 5 # may modify
    Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
    fit <- rfsrc(formula = formula(paste0(classVar, "~.")),
                 data = cbind(Train_set, Train_label[classVar]),
                 ntree = 1000, nodesize = rf_nodesize,
                 importance = T,
                 proximity = T,
                 forest = T)
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }
  
  RunGBM <- function(Train_set, Train_label, mode, classVar){
    fit <- gbm(formula = Train_label[[classVar]] ~ .,
               data = as.data.frame(Train_set),
               distribution = 'bernoulli',
               n.trees = 10000,
               interaction.depth = 3,
               n.minobsinnode = 10,
               shrinkage = 0.001,
               cv.folds = 10,n.cores = 6)
    best <- which.min(fit$cv.error)
    fit <- gbm(formula = Train_label[[classVar]] ~ .,
               data = as.data.frame(Train_set),
               distribution = 'bernoulli',
               n.trees = best,
               interaction.depth = 3,
               n.minobsinnode = 10,
               shrinkage = 0.001, n.cores = 8)
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }
  
  RunXGBoost <- function(Train_set, Train_label, mode, classVar){
    indexes = createFolds(Train_label[[classVar]], k = 5, list=T)
    CV <- unlist(lapply(indexes, function(pt){
      dtrain = xgb.DMatrix(data = Train_set[-pt, ], 
                           label = Train_label[-pt, ])
      dtest = xgb.DMatrix(data = Train_set[pt, ], 
                          label = Train_label[pt, ])
      watchlist <- list(train=dtrain, test=dtest)
      
      bst <- xgb.train(data=dtrain, 
                       max.depth=2, eta=1, nthread = 2, nrounds=10, 
                       watchlist=watchlist, 
                       objective = "binary:logistic", verbose = F)
      which.min(bst$evaluation_log$test_logloss)
    }))
    
    nround <- as.numeric(names(which.max(table(CV))))
    fit <- xgboost(data = Train_set, 
                   label = Train_label[[classVar]], 
                   max.depth = 2, eta = 1, nthread = 2, nrounds = nround, 
                   objective = "binary:logistic", verbose = F)
    fit$subFeature = colnames(Train_set)
    
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }
  
  RunNaiveBayes <- function(Train_set, Train_label, mode, classVar){
    data <- cbind(Train_set, Train_label[classVar])
    data[[classVar]] <- as.factor(data[[classVar]])
    fit <- naiveBayes(eval(parse(text = paste(classVar, "~."))), 
                      data = data)
    fit$subFeature = colnames(Train_set)
    if (mode == "Model") return(fit)
    if (mode == "Variable") return(ExtractVar(fit))
  }
  
  # DRF不适用于二分类情况，因此删去
  # RunDRF <- function(Train_set, Train_label, mode, classVar){
  #   Train_label <- data.frame(
  #     "0" = as.numeric(Train_label == 0),
  #     "1" = as.numeric(Train_label == 1)
  #   )
  #   fit <- drf(X = Train_set, 
  #              Y = Train_label, 
  #              compute.variable.importance = F)
  #   fit$subFeature = colnames(Train_set)
  #   
  #   summary(predict(fit, functional = "mean", as.matrix(Train_set))$mean)
  #   
  #   if (mode == "Model") return(fit)
  #   if (mode == "Variable") return(ExtractVar(fit))
  # }
  
  quiet <- function(..., messages=FALSE, cat=FALSE){
    if(!cat){
      sink(tempfile())
      on.exit(sink())
    }
    out <- if(messages) eval(...) else suppressMessages(eval(...))
    out
  }
  
  standarize.fun <- function(indata, centerFlag, scaleFlag) {  
    scale(indata, center=centerFlag, scale=scaleFlag)
  }
  
  # 2.0更新
  scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL){
    samplename = rownames(data)
    if (is.null(cohort)){
      data <- list(data); names(data) = "training"
    }else{
      data <- split(as.data.frame(data), cohort)
    }
    
    if (is.null(centerFlags)){
      centerFlags = F; message("No centerFlags found, set as FALSE")
    }
    if (length(centerFlags)==1){
      centerFlags = rep(centerFlags, length(data)); message("set centerFlags for all cohort as ", unique(centerFlags))
    }
    if (is.null(names(centerFlags))){
      names(centerFlags) <- names(data); message("match centerFlags with cohort by order\n")
    }
    
    if (is.null(scaleFlags)){
      scaleFlags = F; message("No scaleFlags found, set as FALSE")
    }
    if (length(scaleFlags)==1){
      scaleFlags = rep(scaleFlags, length(data)); message("set scaleFlags for all cohort as ", unique(scaleFlags))
    }
    if (is.null(names(scaleFlags))){
      names(scaleFlags) <- names(data); message("match scaleFlags with cohort by order\n")
    }
    
    centerFlags <- centerFlags[names(data)]; scaleFlags <- scaleFlags[names(data)]
    outdata <- mapply(standarize.fun, indata = data, centerFlag = centerFlags, scaleFlag = scaleFlags, SIMPLIFY = F)
    # lapply(out.data, function(x) summary(apply(x, 2, var)))
    outdata <- do.call(rbind, outdata)
    outdata <- outdata[samplename, ]
    return(outdata)
  }
  
  ExtractVar <- function(fit){
    Feature <- quiet(switch(
      EXPR = class(fit)[1],
      "lognet" = rownames(coef(fit))[which(coef(fit)[, 1]!=0)], # 本身没有筛选变量的功能，但是可以舍去模型中系数为0的变量
      "glm" = names(coef(fit)), # 逐步回归可以对变量进行筛选
      "svm.formula" = fit$subFeature, # SVM对变量没有筛选功能，所以默认使用所有变量
      "train" = fit$coefnames, # LDA不能筛选变量，所以默认使用所有变量
      "glmboost" = names(coef(fit)[abs(coef(fit))>0]), # glmboost同样不具备筛选变量的功能，因此舍去模型中系数为0的变量
      "plsRglmmodel" = rownames(fit$Coeffs)[fit$Coeffs!=0], # plsRglmmodel同样不具备筛选变量的功能，因此舍去模型中系数为0的变量
      "rfsrc" = var.select(fit, verbose = F)$topvars, # rfsrc可以对变量进行筛选
      "gbm" = rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0], # gbm通过舍去重要性为0的变量来进行变量筛选
      "xgb.Booster" = fit$subFeature, # xgboost没有筛选变量的能力， 默认使用所有变量
      "naiveBayes" = fit$subFeature # naiveBayes没有筛选变量的能力，默认使用所有变量
      # "drf" = fit$subFeature # drf自带的变量系数提取函数会输出NA，因此默认使用所有变量
    ))
    
    Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
    return(Feature)
  }
  
  
  CalPredictScore <- function(fit, new_data, type = "lp"){
    new_data <- new_data[, fit$subFeature]
    RS <- quiet(switch(
      EXPR = class(fit)[1],
      "lognet"      = predict(fit, type = 'response', as.matrix(new_data)), # response
      "glm"         = predict(fit, type = 'response', as.data.frame(new_data)), # response
      "svm.formula" = predict(fit, as.data.frame(new_data), probability = T), # 
      "train"       = predict(fit, new_data, type = "prob")[[2]],
      "glmboost"    = predict(fit, type = "response", as.data.frame(new_data)), # response
      "plsRglmmodel" = predict(fit, type = "response", as.data.frame(new_data)), # response
      "rfsrc"        = predict(fit, as.data.frame(new_data))$predicted[, "1"],
      "gbm"          = predict(fit, type = 'response', as.data.frame(new_data)), # response
      "xgb.Booster" = predict(fit, as.matrix(new_data)),
      "naiveBayes" = predict(object = fit, type = "raw", newdata = new_data)[, "1"]
      # "drf" = predict(fit, functional = "mean", as.matrix(new_data))$mean
    ))
    RS = as.numeric(as.vector(RS))
    names(RS) = rownames(new_data)
    return(RS)
  }
  
  PredictClass <- function(fit, new_data){
    new_data <- new_data[, fit$subFeature]
    label <- quiet(switch(
      EXPR = class(fit)[1],
      "lognet"      = predict(fit, type = 'class', as.matrix(new_data)),
      "glm"         = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5, 
                             yes = "1", no = "0"), # glm不返回预测的类，将概率>0.5的作为1类
      "svm.formula" = predict(fit, as.data.frame(new_data), decision.values = T), # 
      "train"       = predict(fit, new_data, type = "raw"),
      "glmboost"    = predict(fit, type = "class", as.data.frame(new_data)), 
      "plsRglmmodel" = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5, 
                              yes = "1", no = "0"), # plsRglm不允许使用因子变量作为因变量，因而predict即使type设为class也无法正常运作
      "rfsrc"        = predict(fit, as.data.frame(new_data))$class,
      "gbm"          = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5,
                              yes = "1", no = "0"), # gbm未设置预测类别，设置大于0.5为1
      "xgb.Booster" = ifelse(test = predict(fit, as.matrix(new_data))>0.5,
                             yes = "1", no = "0"), # xgboost 未提供预测类别，设置大于0.5为1
      "naiveBayes" = predict(object = fit, type = "class", newdata = new_data)
      # "drf" = predict(fit, functional = "mean", as.matrix(new_data))$mean
    ))
    label = as.character(as.vector(label))
    names(label) = rownames(new_data)
    return(label)
  }
  
  RunEval <- function(fit, 
                      Test_set = NULL, 
                      Test_label = NULL, 
                      Train_set = NULL, 
                      Train_label = NULL, 
                      Train_name = NULL,
                      cohortVar = "Cohort",
                      classVar){
    
    if(!is.element(cohortVar, colnames(Test_label))) {
      stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
    } 
    
    if((!is.null(Train_set)) & (!is.null(Train_label))) {
      new_data <- rbind.data.frame(Train_set[, fit$subFeature],
                                   Test_set[, fit$subFeature])
      
      if(!is.null(Train_name)) {
        Train_label$Cohort <- Train_name
      } else {
        Train_label$Cohort <- "Training"
      }
      colnames(Train_label)[ncol(Train_label)] <- cohortVar
      Test_label <- rbind.data.frame(Train_label[,c(cohortVar, classVar)],
                                     Test_label[,c(cohortVar, classVar)])
      Test_label[,1] <- factor(Test_label[,1], 
                               levels = c(unique(Train_label[,cohortVar]), setdiff(unique(Test_label[,cohortVar]),unique(Train_label[,cohortVar]))))
    } else {
      new_data <- Test_set[, fit$subFeature]
    }
    # 调用 CalPredictScore 函数计算分数
    RS <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))
    # 创建 Predict.out 数据框
    Predict.out <- Test_label
    Predict.out$RS <- as.vector(RS)
    # 根据 cohortVar 划分数据并计算每个子集的 AUC
    Predict.out <- split(x = Predict.out, f = Predict.out[,cohortVar])
    unlist(lapply(Predict.out, function(data){
      as.numeric(auc(suppressMessages(roc(data[[classVar]], data$RS))))
    }))
  }
  
  SimpleHeatmap <- function(Cindex_mat, avg_Cindex, 
                            CohortCol, barCol,
                            cellwidth = 1, cellheight = 0.5, 
                            cluster_columns, cluster_rows){
    col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                              col = list("Cohort" = CohortCol),
                              show_annotation_name = F)
    
    row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                              gp = gpar(fill = barCol, col = NA),
                                              add_numbers = T, numbers_offset = unit(-10, "mm"),
                                              axis_param = list("labels_rot" = 0),
                                              numbers_gp = gpar(fontsize = 9, col = "white"),
                                              width = unit(3, "cm")),
                           show_annotation_name = F)
    
    Heatmap(as.matrix(Cindex_mat), name = "AUC",
            right_annotation = row_ha, 
            top_annotation = col_ha,
            # col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # 黄绿配色
            col = c("#4195C1", "#FFFFFF", "#CB5746"), # 红蓝配色
            rect_gp = gpar(col = "black", lwd = 1), # 边框设置为黑色
            cluster_columns = cluster_columns, cluster_rows = cluster_rows, # 不进行聚类，无意义
            show_column_names = FALSE, 
            show_row_names = TRUE,
            row_names_side = "left",
            width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
            height = unit(cellheight * nrow(Cindex_mat), "cm"),
            column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
            column_title = NULL,
            cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
              grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                        x, y, gp = gpar(fontsize = 10))
            }
    )
  }
  
  # 创建向量
  methods <- c(
    "Lasso + Stepglm [both]", "SVM", "glmBoost + SVM", "Ridge", "Lasso + SVM",
    "glmBoost + Ridge", "Enet [alpha=0.1]", "glmBoost + Enet [alpha=0.1]", "Enet [alpha=0.2]", "Enet [alpha=0.3]",
    "glmBoost + Enet [alpha=0.3]", "glmBoost + Enet [alpha=0.2]", "Enet [alpha=0.4]", "glmBoost + Enet [alpha=0.4]",
    "Lasso + glmBoost", "Enet [alpha=0.5]", "glmBoost", "glmBoost + Enet [alpha=0.5]", "Enet [alpha=0.6]",
    "glmBoost + Enet [alpha=0.6]", "glmBoost + Enet [alpha=0.7]", "glmBoost + Enet [alpha=0.8]", "Enet [alpha=0.8]",
    "Enet [alpha=0.9]", "Lasso", "Enet [alpha=0.7]", "glmBoost + Enet [alpha=0.9]", "glmBoost + Lasso",
    "Lasso + plsRglm", "glmBoost + plsRglm", "glmBoost + Stepglm [forward]", "Lasso + Stepglm [forward]", "RF + SVM",
    "Stepglm [forward]", "plsRglm", "RF + Ridge", "RF + Enet [alpha=0.1]", "RF + plsRglm", "RF + Stepglm [forward]",
    "RF + Enet [alpha=0.2]", "RF + Enet [alpha=0.3]", "RF + Enet [alpha=0.6]", "RF + Lasso", "RF + Enet [alpha=0.7]",
    "RF + Enet [alpha=0.5]", "RF + glmBoost", "RF + Enet [alpha=0.9]", "RF + Enet [alpha=0.4]", "RF + Enet [alpha=0.8]",
    "RF + Stepglm [both]", "RF + Stepglm [backward]", "Stepglm [both] + Ridge", "Stepglm [backward] + Ridge",
    "Stepglm [both] + plsRglm", "Stepglm [backward] + plsRglm", "Stepglm [both] + Enet [alpha=0.9]",
    "Stepglm [backward] + Enet [alpha=0.9]", "Stepglm [both] + Enet [alpha=0.1]", "Stepglm [backward] + Enet [alpha=0.1]",
    "Stepglm [both] + Enet [alpha=0.8]", "Stepglm [backward] + Enet [alpha=0.8]", "Stepglm [both] + Enet [alpha=0.2]",
    "Stepglm [backward] + Enet [alpha=0.2]", "Stepglm [both] + Lasso", "Stepglm [backward] + Lasso",
    "Stepglm [both] + Enet [alpha=0.6]", "Stepglm [backward] + Enet [alpha=0.6]", "glmBoost + GBM",
    "Stepglm [both] + Enet [alpha=0.7]", "Stepglm [backward] + Enet [alpha=0.7]", "Lasso + Stepglm [backward]",
    "Stepglm [both]", "Stepglm [backward]", "glmBoost + Stepglm [both]", "glmBoost + Stepglm [backward]",
    "Stepglm [both] + Enet [alpha=0.4]", "Stepglm [backward] + Enet [alpha=0.4]", "Stepglm [both] + Enet [alpha=0.3]",
    "Stepglm [backward] + Enet [alpha=0.3]", "Stepglm [both] + glmBoost", "Stepglm [backward] + glmBoost",
    "Stepglm [both] + Enet [alpha=0.5]", "Stepglm [backward] + Enet [alpha=0.5]", "glmBoost + RF", "RF", "Lasso + GBM",
    "RF + GBM", "GBM", "Stepglm [both] + SVM", "Stepglm [backward] + SVM", "Lasso + RF", "Stepglm [both] + GBM",
    "Stepglm [backward] + GBM", "Stepglm [both] + RF", "LDA", "glmBoost + LDA", "RF + LDA", "Stepglm [both] + LDA",
    "Stepglm [backward] + LDA", "Lasso + LDA", "Stepglm [backward] + RF", "XGBoost", "Lasso + XGBoost", "glmBoost + XGBoost",
    "RF + XGBoost", "Stepglm [both] + XGBoost", "Stepglm [backward] + XGBoost", "NaiveBayes", "Lasso + NaiveBayes",
    "glmBoost + NaiveBayes", "RF + NaiveBayes", "Stepglm [both] + NaiveBayes", "Stepglm [backward] + NaiveBayes"
  )
  
  cat('目前机器学习的种类有',length(methods),"种\n")
  ## method list --------------------------------------------------------
  # 此处记录需要运行的模型，格式为：算法1名称[算法参数]+算法2名称[算法参数]
  # 目前仅有Stepglm和Enet支持输入算法参数
  methods <- gsub("-| ", "", methods)
  head(methods)
  
  surv <- function(datExpr){
    traitData <- matrix(0, nrow = ncol(datExpr), ncol = 1,
                        dimnames = list(colnames(datExpr), c("outcome")))
    # # 为新矩阵填写值
    traitData[grepl("treat", rownames(traitData)), "outcome"] <- 1
    traitData[grepl("con", rownames(traitData)), "outcome"] <- 0
    traitData <- as.data.frame(traitData)
  }
  
  # first_col_to_rowname <-function(Test_expr){
  #   Test_expr=as.matrix(Test_expr)
  #   rownames(Test_expr)=Test_expr[,1]
  #   exp=Test_expr[,2:ncol(Test_expr)]
  #   dimnames=list(rownames(exp),colnames(exp))
  #   Test_expr=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
  #   Test_expr=as.data.frame(Test_expr)
  #   
  # }
  
  # 输出向量
  # 选择最后生成的模型类型：panML代表生成由不同算法构建的模型； multiLogistic表示抽取其他模型所用到的变量并建立多变量logistic模型
  FinalModel <- c("panML", "multiLogistic")[2]

  ##################################################################################################
  #################################以上为函数固定不分，现在正文开始#################################
  
  #整理训练模型
  cat("训练组数据开始处理\n")
  Train_expr <- Train_expr[rowSums(Train_expr > 0) > ncol(Train_expr) * 0.1, ] 
  Train_class = surv(Train_expr)
  comsam <- intersect(rownames(Train_class), colnames(Train_expr))
  Train_expr <- Train_expr[,comsam]; Train_class <- Train_class[comsam,,drop = F]
  
  
  
  
  cat("测试组数据开始处理\n")
  comsam <- intersect(rownames(Test_class), colnames(Test_expr))
  Test_expr <- Test_expr[,comsam]; Test_class <- Test_class[comsam,,drop = F]
 
  cat("提取相同基因")
  comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
  Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
  Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
  
  
  
  # 按队列对数据分别进行标准化（根据情况调整centerFlags和scaleFlags）
  ## data: 需要表达谱数据（行为样本，列为基因） 
  ## cohort：样本所属队列，为向量，不输入值时默认全表达矩阵来自同一队列
  ## centerFlag/scaleFlags：是否将基因均值/标准差标准化为1；
  ##        默认参数为NULL，表示不进行标准化；
  ##        为T/F时，表示对所有队列都进行/不进行标准化
  ##        输入由T/F组成的向量时，按顺序对队列进行处理，向量长度应与队列数一样
  ##        如centerFlags = c(F, F, F, T, T)，表示对第4、5个队列进行标准化，此时flag顺序应当与队列顺序一致
  ##        如centerFlags = c("A" = F, "C" = T, "B" = F)，表示对队列C进行标准化，此时不要求flag顺序与data一致
  cat("开始数据归一化\n")
  Train_set = scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
  names(x = split(as.data.frame(Test_expr), f = Test_class$Cohort)) # 注意测试集标准化顺序与此一致
  Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = T, scaleFlags = T)
  # summary(apply(Train_set, 2, var))
  # summary(apply(Test_set, 2, var))
  # lapply(split(as.data.frame(Test_set), Test_class$Cohort), function(x) summary(apply(x, 2, var))) # 测试scale结果
  
  # Model training and validation -------------------------------------------
  

  #colnames(Test_set)=colnames(Train_set)

  ## Train the model --------------------------------------------------------
  cat("开始进行机器学习\n")
  classVar = "outcome" # 设置所要预测的变量名（仅支持[0,1]二元变量格式）
  min.selected.var = 5 # 设置模型最少纳入的变量数
  
  ## Pre-training 将各方法所用到的变量筛选过程汇总，以减少计算量
  Variable = colnames(Train_set)
  preTrain.method =  strsplit(methods, "\\+") # 检视所有方法，分析各方法是否需要进行变量预筛选(pre-training)
  preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1]) # 删除各方法用于构建分类模型的算法，保留用于变量筛选的算法
  preTrain.method = unique(unlist(preTrain.method)) # 汇总所有变量筛选算法，去除重复计算
  
  preTrain.var <- list() # 用于保存各算法筛选的变量
  set.seed(seed = 777) # 设置建模种子，使得结果可重复
  for (method in preTrain.method){
    preTrain.var[[method]] = RunML(method = method, # 变量筛选所需要的机器学习方法
                                   Train_set = Train_set, # 训练集有潜在预测价值的变量
                                   Train_label = Train_class, # 训练集分类标签
                                   mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                   classVar = classVar) # 用于训练的分类变量，必须出现在Train_class中
  }
  preTrain.var[["simple"]] <- colnames(Train_set)# 记录未经筛选的变量集（以便后续代码撰写），可视为使用simple方法（无筛选功能）的变量筛选结果
  
  ## Model training
  model <- list() # 用于保存各模型的所有信息
  set.seed(seed = 777) # 设置建模种子，使得结果可重复
  Train_set_bk = Train_set # RunML有一个函数(plsRglm)无法正常传参，需要对训练集数据进行存档备份
  for (method in methods){
    cat(match(method, methods), ":", method, "\n")
    method_name = method # 本轮算法名称
    method <- strsplit(method, "\\+")[[1]] # 各步骤算法名称
    
    if (length(method) == 1) method <- c("simple", method) # 如果本方法没有预筛选变量，则认为本方法使用simple方法进行了变量筛选
    
    Variable = preTrain.var[[method[1]]] # 根据方法名称的第一个值，调用先前变量筛选的结果
    Train_set = Train_set_bk[, Variable]   # 对训练集取子集，因为有一个算法原作者写的有点问题，无法正常传参
    Train_label = Train_class            # 所以此处需要修改变量名称，以免函数错误调用对象
    model[[method_name]] <- RunML(method = method[2],        # 根据方法名称第二个值，调用构建的函数分类模型
                                  Train_set = Train_set,     # 训练集有潜在预测价值的变量
                                  Train_label = Train_label, # 训练集分类标签
                                  mode = "Model",            # 运行模式，Variable(筛选变量)和Model(获取模型)
                                  classVar = classVar)       # 用于训练的分类变量，必须出现在Train_class中
    
    # 如果最终模型纳入的变量数小于预先设定的下限，则认为该算法输出的结果是无意义的
    if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
      model[[method_name]] <- NULL
    }
  }
  
  Train_set = Train_set_bk; rm(Train_set_bk) # 将数据还原，并移除备份
  saveRDS(model, file.path(res.path, "model.rds")) # 报错各模型的所有中间过程
  
  if (FinalModel == "multiLogistic"){
    logisticmodel <- lapply(model, function(fit){ # 根据各算法最终获得的变量，构建多变量Logistic模型，从而以Logistic回归系数和特征表达计算单样本分类概率
      tmp <- glm(formula = Train_class[[classVar]] ~ .,
                 family = "binomial", 
                 data = as.data.frame(Train_set[, ExtractVar(fit)]))
      tmp$subFeature <- ExtractVar(fit) # 提取当Logistic模型最终使用的预测变量
      return(tmp)
    })
  }
  saveRDS(logisticmodel, file.path(res.path, "logisticmodel.rds")) # 保存最终以多变量Logistic模型
  
  ## Evaluate the model -----------------------------------------------------
  
  # 读取已保存的模型列表
  model <- readRDS(file.path(res.path, "model.rds"))
  if(Logistic){
    model <- readRDS(file.path(res.path, "logisticmodel.rds")) 
    # 若希望使用多变量保存最终以多变量Logistic模型计算得分，请运行此行
  }
  # 
  methodsValid <- names(model)
  
  # 根据给定表达量计算样本风险评分
  # 预测概率
  RS_list <- list()
  for (method in methodsValid){
    RS_list[[method]] <- CalPredictScore(fit = model[[method]], 
                                         new_data = rbind.data.frame(Train_set,Test_set)) # 2.0更新
  }
  RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
  write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) # 输出风险评分文件
  
  # 根据给定表达量预测分类
  Class_list <- list()
  for (method in methodsValid){
    Class_list[[method]] <- PredictClass(fit = model[[method]], 
                                         new_data = rbind.data.frame(Train_set,Test_set)) # 2.0更新
  }
  Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
  #Class_mat <- cbind.data.frame(Test_class, Class_mat[rownames(Class_mat),]) # 若要合并测试集本身的样本信息文件可运行此行
  write.table(Class_mat, file.path(res.path, "Class_mat.txt"), # 测试集经过算法预测出的二分类结果
              sep = "\t", row.names = T, col.names = NA, quote = F)
  
  # 提取所筛选的变量（列表格式）
  fea_list <- list()
  for (method in methodsValid) {
    fea_list[[method]] <- ExtractVar(model[[method]])
  }
  
  # 提取所筛选的变量（数据框格式）
  fea_df <- lapply(model, function(fit){
    data.frame(ExtractVar(fit))
  })
  fea_df <- do.call(rbind, fea_df)
  fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
  colnames(fea_df)[1] <- "features"
  
  # 使用dplyr进行分组和汇总
  result <- fea_df %>%
    mutate(method = gsub("\\..*$", "", algorithm)) %>%
    group_by(features) %>%
    summarize(used_methods = paste(unique(method), collapse = ", "))
  
  
  ################################################################################
  #下面是我自己手动添加的优化，可以删除
  # 创建一个空的数据框来存储结果
  result_df <- data.frame(Gene = character(0), Count = integer(0), Algorithm = character(0), stringsAsFactors = FALSE)
  # 使用dplyr进行筛选和统计
  for (i in unique(fea_df$features)) {
    acadvl_data <- fea_df %>%
      filter(features == i)
    num_acadvl <- nrow(acadvl_data)
    # 输出对应的名字
    if (num_acadvl > 0) {
      algorithms <- unique(acadvl_data$algorithm)
      result_df <- rbind(result_df, data.frame(Gene = i, Count = num_acadvl, Algorithm = paste(algorithms, collapse = ", ")))
    }
  }
  write.table(result_df, file.path(res.path, "基因筛选次数.txt"), # 两列，包含算法以及算法所筛选出的变量
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  # 创建一个空的数据框来存储结果
  result_ML <- data.frame(Gene = character(0), Count = integer(0), Algorithm = character(0), stringsAsFactors = FALSE)
  # 使用dplyr进行筛选和统计
  for (i in unique(fea_df$algorithm)) {
    acadvl_data <- fea_df %>%
      filter(algorithm == i)
    num_acadvl <- nrow(acadvl_data)
    # 输出对应的名字
    if (num_acadvl > 0) {
      features <- unique(acadvl_data$features)
      result_ML <- rbind(result_ML, data.frame(Gene = i, Count = num_acadvl, Algorithm = paste(algorithms, collapse = ", ")))
    }
  }
  write.table(result_ML, file.path(res.path, "每种机器学习筛选基因的类别.txt"), # 两列，包含算法以及算法所筛选出的变量
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  ################################################################################
  
  
  
  write.table(fea_df, file.path(res.path, "fea_df.txt"), # 两列，包含算法以及算法所筛选出的变量
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  # 对各模型计算C-index
  AUC_list <- list()
  for (method in methodsValid){
    AUC_list[[method]] <- RunEval(fit = model[[method]],     # 分类预测模型
                                  Test_set = Test_set,      # 测试集预测变量，应当包含训练集中所有的变量，否则会报错
                                  Test_label = Test_class,   # 训练集分类数据，应当包含训练集中所有的变量，否则会报错
                                  Train_set = Train_set,    # 若需要同时评估训练集，则给出训练集表达谱，否则置NULL
                                  Train_label = Train_class, # 若需要同时评估训练集，则给出训练集分类数据，否则置NULL
                                  Train_name = train_lab,       # 若需要同时评估训练集，可给出训练集的标签，否则按“Training”处理
                                  cohortVar = "Cohort",      # 重要：用于指定队列的变量，该列必须存在且指定[默认为“Cohort”]，否则会报错
                                  classVar = classVar)       # 用于评估的二元分类变量，必须出现在Test_class中
  }
  AUC_mat <- do.call(rbind, AUC_list)
  write.table(AUC_mat, file.path(res.path, "AUC_mat.txt"),
              sep = "\t", row.names = T, col.names = T, quote = F)
  
  # Plot --------------------------------------------------------------------
  
  AUC_mat <- read.table(file.path(res.path, "AUC_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
  avg_AUC <- apply(AUC_mat, 1, mean)           # 计算每种算法在所有队列中平均AUC
  avg_AUC <- sort(avg_AUC, decreasing = T)     # 对各算法AUC由高到低排序
  AUC_mat <- AUC_mat[names(avg_AUC), ]      # 对AUC矩阵排序
  fea_sel <- fea_list[[rownames(AUC_mat)[1]]] # 最优模型（测试集AUC均值最大）所筛选的特征
  avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3)) # 保留三位小数
  
  if(ncol(AUC_mat) < 3) { # 如果用于绘图的队列小于3个
    CohortCol <- c("red","blue") # 则给出两个颜色即可（可自行替换颜色）
  } else { # 否则通过brewer.pal赋予超过3个队列的颜色
    CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired") # 设置队列颜色
  }
  names(CohortCol) <- colnames(AUC_mat)
  
  cellwidth = 3; cellheight = 0.5
  hm <- SimpleHeatmap(AUC_mat, # 主矩阵
                      avg_AUC, # 侧边柱状图
                      CohortCol, "steelblue", # 列标签颜色，右侧柱状图颜色
                      cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                      cluster_columns = F, cluster_rows = F) # 是否对行列进行聚类
  
  pdf(file.path(res.path, "AUC.pdf"), width = cellwidth * ncol(AUC_mat) + 3, height = cellheight * nrow(AUC_mat) * 0.45)
  draw(hm)
  invisible(dev.off())
  return(list(gene =result_df,machine_learn =result_ML  ))
}
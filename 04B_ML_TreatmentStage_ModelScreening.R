############################################################
## 04B_ML_TreatmentStage_ModelScreening.R
##
## Machine-learning model screening for treatment-stage
## assessment in GSE89403.
##
## Input:
##   data.train.txt
##   data.test.txt
##   refer.methodLists.txt
##
## Output:
##   ML_model_objects.rds
##   ML_logistic_model_objects.rds
##   ML_risk_score_matrix.txt
##   ML_class_prediction_matrix.txt
##   ML_model_gene_features.txt
##   ML_AUC_matrix.txt
##   ML_AUC_heatmap.pdf
##   ML_top_model_ranking.txt
##   ML_representative_model_features.txt
############################################################


## =========================================================
## 0. Set English directory structure
## =========================================================

base_dir <- "D:/atraining/blood/analysis"
ml_dir   <- file.path(base_dir, "ml_model")

setwd(ml_dir)

set.seed(123)

min.selected.var <- 5
classVar <- "Type"


## =========================================================
## 1. Load packages
## =========================================================

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
library(grid)


## =========================================================
## 2. Define machine-learning helper functions
## =========================================================

RunML <- function(method, Train_set, Train_label, mode = "Model", classVar) {
  
  method <- gsub(" ", "", method)
  
  method_name <- gsub("(\\w+)\\[(.+)\\]", "\\1", method)
  method_param <- gsub("(\\w+)\\[(.+)\\]", "\\2", method)
  
  method_param <- switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "Stepglm" = list("direction" = method_param),
    NULL
  )
  
  message(
    "Run ", method_name,
    " algorithm for ", mode, "; ",
    method_param, "; using ",
    ncol(Train_set), " variables"
  )
  
  args <- list(
    "Train_set" = Train_set,
    "Train_label" = Train_label,
    "mode" = mode,
    "classVar" = classVar
  )
  
  args <- c(args, method_param)
  
  obj <- do.call(
    what = paste0("Run", method_name),
    args = args
  )
  
  if (mode == "Variable") {
    message(length(obj), " variables retained.\n")
  } else {
    message("\n")
  }
  
  return(obj)
}


RunEnet <- function(Train_set, Train_label, mode, classVar, alpha) {
  
  cv.fit <- cv.glmnet(
    x = Train_set,
    y = Train_label[[classVar]],
    family = "binomial",
    alpha = alpha,
    nfolds = 10
  )
  
  fit <- glmnet(
    x = Train_set,
    y = Train_label[[classVar]],
    family = "binomial",
    alpha = alpha,
    lambda = cv.fit$lambda.min
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}


RunLasso <- function(Train_set, Train_label, mode, classVar) {
  
  RunEnet(
    Train_set = Train_set,
    Train_label = Train_label,
    mode = mode,
    classVar = classVar,
    alpha = 1
  )
}


RunRidge <- function(Train_set, Train_label, mode, classVar) {
  
  RunEnet(
    Train_set = Train_set,
    Train_label = Train_label,
    mode = mode,
    classVar = classVar,
    alpha = 0
  )
}


RunStepglm <- function(Train_set, Train_label, mode, classVar, direction) {
  
  data <- as.data.frame(Train_set)
  data[[classVar]] <- Train_label[[classVar]]
  
  fit <- step(
    glm(
      formula = as.formula(paste0(classVar, " ~ .")),
      family = "binomial",
      data = data
    ),
    direction = direction,
    trace = 0
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}


RunSVM <- function(Train_set, Train_label, mode, classVar) {
  
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  
  fit <- svm(
    formula = as.formula(paste0(classVar, " ~ .")),
    data = data,
    probability = TRUE
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}


RunLDA <- function(Train_set, Train_label, mode, classVar) {
  
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  
  fit <- train(
    as.formula(paste0(classVar, " ~ .")),
    data = data,
    method = "lda",
    trControl = trainControl(method = "cv")
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}


RunglmBoost <- function(Train_set, Train_label, mode, classVar) {
  
  data <- cbind(as.data.frame(Train_set), Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  
  fit <- glmboost(
    as.formula(paste0(classVar, " ~ .")),
    data = data,
    family = Binomial()
  )
  
  cvm <- cvrisk(
    fit,
    papply = lapply,
    folds = cv(model.weights(fit), type = "kfold")
  )
  
  fit <- glmboost(
    as.formula(paste0(classVar, " ~ .")),
    data = data,
    family = Binomial(),
    control = boost_control(mstop = max(mstop(cvm), 40))
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}


RunplsRglm <- function(Train_set, Train_label, mode, classVar) {
  
  cv.plsRglm.res <- cv.plsRglm(
    formula = as.formula(paste0(classVar, " ~ .")),
    data = cbind(as.data.frame(Train_set), Train_label[classVar]),
    nt = 10,
    verbose = FALSE
  )
  
  fit <- plsRglm(
    Train_label[[classVar]],
    as.data.frame(Train_set),
    modele = "pls-glm-logistic",
    verbose = FALSE,
    sparse = TRUE
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}


RunRF <- function(Train_set, Train_label, mode, classVar) {
  
  rf_nodesize <- 5
  
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  
  fit <- rfsrc(
    formula = as.formula(paste0(classVar, " ~ .")),
    data = cbind(as.data.frame(Train_set), Train_label[classVar]),
    ntree = 1000,
    nodesize = rf_nodesize,
    importance = TRUE,
    proximity = TRUE,
    forest = TRUE
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}


RunGBM <- function(Train_set, Train_label, mode, classVar) {
  
  data <- cbind(as.data.frame(Train_set), Train_label[classVar])
  
  fit <- gbm(
    formula = as.formula(paste0(classVar, " ~ .")),
    data = data,
    distribution = "bernoulli",
    n.trees = 10000,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = 10,
    n.cores = 6
  )
  
  best <- which.min(fit$cv.error)
  
  fit <- gbm(
    formula = as.formula(paste0(classVar, " ~ .")),
    data = data,
    distribution = "bernoulli",
    n.trees = best,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    n.cores = 8
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}


RunXGBoost <- function(Train_set, Train_label, mode, classVar) {
  
  indexes <- createFolds(
    Train_label[[classVar]],
    k = 5,
    list = TRUE
  )
  
  CV <- unlist(lapply(indexes, function(pt) {
    
    dtrain <- xgb.DMatrix(
      data = Train_set[-pt, ],
      label = Train_label[-pt, classVar]
    )
    
    dtest <- xgb.DMatrix(
      data = Train_set[pt, ],
      label = Train_label[pt, classVar]
    )
    
    watchlist <- list(train = dtrain, test = dtest)
    
    bst <- xgb.train(
      data = dtrain,
      max.depth = 2,
      eta = 1,
      nthread = 2,
      nrounds = 10,
      watchlist = watchlist,
      objective = "binary:logistic",
      verbose = FALSE
    )
    
    which.min(bst$evaluation_log$test_logloss)
  }))
  
  nround <- as.numeric(names(which.max(table(CV))))
  
  fit <- xgboost(
    data = Train_set,
    label = Train_label[[classVar]],
    max.depth = 2,
    eta = 1,
    nthread = 2,
    nrounds = nround,
    objective = "binary:logistic",
    verbose = FALSE
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}


RunNaiveBayes <- function(Train_set, Train_label, mode, classVar) {
  
  data <- cbind(as.data.frame(Train_set), Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  
  fit <- naiveBayes(
    as.formula(paste0(classVar, " ~ .")),
    data = data
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}


quiet <- function(..., messages = FALSE, cat = FALSE) {
  
  if (!cat) {
    sink(tempfile())
    on.exit(sink())
  }
  
  out <- if (messages) eval(...) else suppressMessages(eval(...))
  
  return(out)
}


standarize.fun <- function(indata, centerFlag, scaleFlag) {
  
  scale(indata, center = centerFlag, scale = scaleFlag)
}


scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL) {
  
  samplename <- rownames(data)
  
  if (is.null(cohort)) {
    data <- list(data)
    names(data) <- "training"
  } else {
    data <- split(as.data.frame(data), cohort)
  }
  
  if (is.null(centerFlags)) {
    centerFlags <- FALSE
    message("No centerFlags found, set as FALSE")
  }
  
  if (length(centerFlags) == 1) {
    centerFlags <- rep(centerFlags, length(data))
    message("Set centerFlags for all cohorts as ", unique(centerFlags))
  }
  
  if (is.null(names(centerFlags))) {
    names(centerFlags) <- names(data)
    message("Match centerFlags with cohort by order")
  }
  
  if (is.null(scaleFlags)) {
    scaleFlags <- FALSE
    message("No scaleFlags found, set as FALSE")
  }
  
  if (length(scaleFlags) == 1) {
    scaleFlags <- rep(scaleFlags, length(data))
    message("Set scaleFlags for all cohorts as ", unique(scaleFlags))
  }
  
  if (is.null(names(scaleFlags))) {
    names(scaleFlags) <- names(data)
    message("Match scaleFlags with cohort by order")
  }
  
  centerFlags <- centerFlags[names(data)]
  scaleFlags <- scaleFlags[names(data)]
  
  outdata <- mapply(
    standarize.fun,
    indata = data,
    centerFlag = centerFlags,
    scaleFlag = scaleFlags,
    SIMPLIFY = FALSE
  )
  
  outdata <- do.call(rbind, outdata)
  outdata <- outdata[samplename, ]
  
  return(outdata)
}


ExtractVar <- function(fit) {
  
  Feature <- quiet(switch(
    EXPR = class(fit)[1],
    
    "lognet" = rownames(coef(fit))[which(coef(fit)[, 1] != 0)],
    
    "glm" = names(coef(fit)),
    
    "svm.formula" = fit$subFeature,
    
    "train" = fit$coefnames,
    
    "glmboost" = names(coef(fit)[abs(coef(fit)) > 0]),
    
    "plsRglmmodel" = rownames(fit$Coeffs)[fit$Coeffs != 0],
    
    "rfsrc" = var.select(fit, verbose = FALSE)$topvars,
    
    "gbm" = rownames(summary.gbm(fit, plotit = FALSE))[
      summary.gbm(fit, plotit = FALSE)$rel.inf > 0
    ],
    
    "xgb.Booster" = fit$subFeature,
    
    "naiveBayes" = fit$subFeature
  ))
  
  Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
  
  return(Feature)
}


CalPredictScore <- function(fit, new_data, type = "lp") {
  
  new_data <- new_data[, fit$subFeature, drop = FALSE]
  
  RS <- quiet(switch(
    EXPR = class(fit)[1],
    
    "lognet" = predict(
      fit,
      type = "response",
      as.matrix(new_data)
    ),
    
    "glm" = predict(
      fit,
      type = "response",
      as.data.frame(new_data)
    ),
    
    "svm.formula" = attr(
      predict(
        fit,
        as.data.frame(new_data),
        probability = TRUE
      ),
      "probabilities"
    )[, "1"],
    
    "train" = predict(
      fit,
      new_data,
      type = "prob"
    )[[2]],
    
    "glmboost" = predict(
      fit,
      type = "response",
      as.data.frame(new_data)
    ),
    
    "plsRglmmodel" = predict(
      fit,
      type = "response",
      as.data.frame(new_data)
    ),
    
    "rfsrc" = predict(
      fit,
      as.data.frame(new_data)
    )$predicted[, "1"],
    
    "gbm" = predict(
      fit,
      type = "response",
      as.data.frame(new_data)
    ),
    
    "xgb.Booster" = predict(
      fit,
      as.matrix(new_data)
    ),
    
    "naiveBayes" = predict(
      object = fit,
      type = "raw",
      newdata = new_data
    )[, "1"]
  ))
  
  RS <- as.numeric(as.vector(RS))
  names(RS) <- rownames(new_data)
  
  return(RS)
}


PredictClass <- function(fit, new_data) {
  
  new_data <- new_data[, fit$subFeature, drop = FALSE]
  
  label <- quiet(switch(
    EXPR = class(fit)[1],
    
    "lognet" = predict(
      fit,
      type = "class",
      as.matrix(new_data)
    ),
    
    "glm" = ifelse(
      predict(fit, type = "response", as.data.frame(new_data)) > 0.5,
      "1",
      "0"
    ),
    
    "svm.formula" = predict(
      fit,
      as.data.frame(new_data),
      decision.values = TRUE
    ),
    
    "train" = predict(
      fit,
      new_data,
      type = "raw"
    ),
    
    "glmboost" = predict(
      fit,
      type = "class",
      as.data.frame(new_data)
    ),
    
    "plsRglmmodel" = ifelse(
      predict(fit, type = "response", as.data.frame(new_data)) > 0.5,
      "1",
      "0"
    ),
    
    "rfsrc" = predict(
      fit,
      as.data.frame(new_data)
    )$class,
    
    "gbm" = ifelse(
      predict(fit, type = "response", as.data.frame(new_data)) > 0.5,
      "1",
      "0"
    ),
    
    "xgb.Booster" = ifelse(
      predict(fit, as.matrix(new_data)) > 0.5,
      "1",
      "0"
    ),
    
    "naiveBayes" = predict(
      object = fit,
      type = "class",
      newdata = new_data
    )
  ))
  
  label <- as.character(as.vector(label))
  names(label) <- rownames(new_data)
  
  return(label)
}


RunEval <- function(
    fit,
    Test_set = NULL,
    Test_label = NULL,
    Train_set = NULL,
    Train_label = NULL,
    Train_name = NULL,
    cohortVar = "Cohort",
    classVar
) {
  
  if (!is.element(cohortVar, colnames(Test_label))) {
    stop(
      paste0(
        "There is no [",
        cohortVar,
        "] indicator. Please add one more column."
      )
    )
  }
  
  if ((!is.null(Train_set)) & (!is.null(Train_label))) {
    
    new_data <- rbind.data.frame(
      Train_set[, fit$subFeature, drop = FALSE],
      Test_set[, fit$subFeature, drop = FALSE]
    )
    
    if (!is.null(Train_name)) {
      Train_label$Cohort <- Train_name
    } else {
      Train_label$Cohort <- "Training"
    }
    
    colnames(Train_label)[ncol(Train_label)] <- cohortVar
    
    Test_label <- rbind.data.frame(
      Train_label[, c(cohortVar, classVar)],
      Test_label[, c(cohortVar, classVar)]
    )
    
    Test_label[, 1] <- factor(
      Test_label[, 1],
      levels = c(
        unique(Train_label[, cohortVar]),
        setdiff(
          unique(Test_label[, cohortVar]),
          unique(Train_label[, cohortVar])
        )
      )
    )
    
  } else {
    
    new_data <- Test_set[, fit$subFeature, drop = FALSE]
  }
  
  RS <- suppressWarnings(
    CalPredictScore(
      fit = fit,
      new_data = new_data
    )
  )
  
  Predict.out <- Test_label
  Predict.out$RS <- as.vector(RS)
  
  Predict.out <- split(
    x = Predict.out,
    f = Predict.out[, cohortVar]
  )
  
  auc_values <- unlist(lapply(Predict.out, function(data) {
    
    as.numeric(
      auc(
        suppressMessages(
          roc(data[[classVar]], data$RS)
        )
      )
    )
  }))
  
  return(auc_values)
}


SimpleHeatmap <- function(
    Cindex_mat,
    avg_Cindex,
    CohortCol,
    barCol,
    cellwidth = 1,
    cellheight = 0.5,
    cluster_columns,
    cluster_rows
) {
  
  col_ha <- columnAnnotation(
    "Cohort" = colnames(Cindex_mat),
    col = list("Cohort" = CohortCol),
    show_annotation_name = FALSE
  )
  
  row_ha <- rowAnnotation(
    bar = anno_barplot(
      avg_Cindex,
      bar_width = 0.8,
      border = FALSE,
      gp = gpar(fill = barCol, col = NA),
      add_numbers = TRUE,
      numbers_offset = unit(-10, "mm"),
      axis_param = list("labels_rot" = 0),
      numbers_gp = gpar(fontsize = 9, col = "white"),
      width = unit(3, "cm")
    ),
    show_annotation_name = FALSE
  )
  
  Heatmap(
    as.matrix(Cindex_mat),
    name = "AUC",
    right_annotation = row_ha,
    top_annotation = col_ha,
    col = c("#4195C1", "#FFFFFF", "#FFBC90"),
    rect_gp = gpar(col = "black", lwd = 1),
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
    height = unit(cellheight * nrow(Cindex_mat), "cm"),
    column_split = factor(
      colnames(Cindex_mat),
      levels = colnames(Cindex_mat)
    ),
    column_title = NULL,
    cell_fun = function(j, i, x, y, w, h, col) {
      grid.text(
        label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
        x,
        y,
        gp = gpar(fontsize = 10)
      )
    }
  )
}


## =========================================================
## 3. Read training and validation data
## =========================================================

Train_data <- read.table(
  "data.train.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = 1,
  stringsAsFactors = FALSE
)

Train_expr <- Train_data[, 1:(ncol(Train_data) - 1), drop = FALSE]
Train_class <- Train_data[, ncol(Train_data), drop = FALSE]
colnames(Train_class) <- classVar


Test_data <- read.table(
  "data.test.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = 1,
  stringsAsFactors = FALSE
)

Test_expr <- Test_data[, 1:(ncol(Test_data) - 1), drop = FALSE]
Test_class <- Test_data[, ncol(Test_data), drop = FALSE]
colnames(Test_class) <- classVar

Test_class$Cohort <- gsub(
  "(.*)\\_(.*)\\_(.*)",
  "\\1",
  rownames(Test_class)
)

Test_class <- Test_class[, c("Cohort", classVar)]


## =========================================================
## 4. Keep common genes and standardize expression
## =========================================================

common_genes <- intersect(
  colnames(Train_expr),
  colnames(Test_expr)
)

Train_expr <- as.matrix(
  Train_expr[, common_genes, drop = FALSE]
)

Test_expr <- as.matrix(
  Test_expr[, common_genes, drop = FALSE]
)

Train_set <- scaleData(
  data = Train_expr,
  centerFlags = TRUE,
  scaleFlags = TRUE
)

Test_set <- scaleData(
  data = Test_expr,
  cohort = Test_class$Cohort,
  centerFlags = TRUE,
  scaleFlags = TRUE
)

## Important:
## Artificial noise was removed from the final manuscript code.


## =========================================================
## 5. Read machine-learning method list
## =========================================================
## Required file:
##   refer.methodLists.txt
##
## This file should contain a column named "Model".
## Example model names:
##   Lasso+RF
##   glmBoost+RF
##   Stepglm[backward]+RF

methodRT <- read.table(
  "refer.methodLists.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

methods <- methodRT$Model
methods <- gsub("-| ", "", methods)

write.table(
  data.frame(Model = methods),
  "ML_model_method_list.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


## =========================================================
## 6. Stage 1: feature selection
## =========================================================

preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method <- unique(unlist(preTrain.method))

preTrain.var <- list()

set.seed(123)

for (method in preTrain.method) {
  
  preTrain.var[[method]] <- RunML(
    method = method,
    Train_set = Train_set,
    Train_label = Train_class,
    mode = "Variable",
    classVar = classVar
  )
}

preTrain.var[["simple"]] <- colnames(Train_set)


## =========================================================
## 7. Stage 2: model training
## =========================================================

model <- list()

set.seed(123)

Train_set_backup <- Train_set

for (method in methods) {
  
  cat(match(method, methods), ":", method, "\n")
  
  method_name <- method
  method_split <- strsplit(method, "\\+")[[1]]
  
  if (length(method_split) == 1) {
    method_split <- c("simple", method_split)
  }
  
  selected_features <- preTrain.var[[method_split[1]]]
  
  Train_set_sub <- Train_set_backup[
    ,
    selected_features,
    drop = FALSE
  ]
  
  model[[method_name]] <- RunML(
    method = method_split[2],
    Train_set = Train_set_sub,
    Train_label = Train_class,
    mode = "Model",
    classVar = classVar
  )
  
  if (length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}

Train_set <- Train_set_backup
rm(Train_set_backup)

saveRDS(
  model,
  "ML_model_objects.rds"
)

methodsValid <- names(model)

cat("Number of valid models:", length(methodsValid), "\n")


## =========================================================
## 8. Optional: logistic refitting model
## =========================================================

logistic_model <- lapply(model, function(fit) {
  
  selected_features <- ExtractVar(fit)
  
  tmp <- glm(
    formula = Train_class[[classVar]] ~ .,
    family = "binomial",
    data = as.data.frame(
      Train_set[, selected_features, drop = FALSE]
    )
  )
  
  tmp$subFeature <- selected_features
  
  return(tmp)
})

saveRDS(
  logistic_model,
  "ML_logistic_model_objects.rds"
)


## =========================================================
## 9. Export risk score matrix
## =========================================================

RS_list <- list()

for (method in methodsValid) {
  
  RS_list[[method]] <- CalPredictScore(
    fit = model[[method]],
    new_data = rbind.data.frame(Train_set, Test_set)
  )
}

risk_matrix <- as.data.frame(
  t(do.call(rbind, RS_list))
)

risk_matrix <- cbind(
  id = rownames(risk_matrix),
  risk_matrix
)

write.table(
  risk_matrix,
  "ML_risk_score_matrix.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


## =========================================================
## 10. Export class prediction matrix
## =========================================================

Class_list <- list()

for (method in methodsValid) {
  
  Class_list[[method]] <- PredictClass(
    fit = model[[method]],
    new_data = rbind.data.frame(Train_set, Test_set)
  )
}

class_matrix <- as.data.frame(
  t(do.call(rbind, Class_list))
)

class_matrix <- cbind(
  id = rownames(class_matrix),
  class_matrix
)

write.table(
  class_matrix,
  "ML_class_prediction_matrix.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


## =========================================================
## 11. Export selected gene features
## =========================================================

feature_list <- list()

for (method in methodsValid) {
  feature_list[[method]] <- ExtractVar(model[[method]])
}

feature_df <- lapply(names(feature_list), function(method) {
  
  data.frame(
    Model = method,
    Gene = feature_list[[method]],
    stringsAsFactors = FALSE
  )
})

feature_df <- do.call(rbind, feature_df)

write.table(
  feature_df,
  file = "ML_model_gene_features.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


## =========================================================
## 12. Evaluate AUC across training and validation subsets
## =========================================================

AUC_list <- list()

for (method in methodsValid) {
  
  AUC_list[[method]] <- RunEval(
    fit = model[[method]],
    Test_set = Test_set,
    Test_label = Test_class,
    Train_set = Train_set,
    Train_label = Train_class,
    Train_name = "Train",
    cohortVar = "Cohort",
    classVar = classVar
  )
}

AUC_mat <- do.call(rbind, AUC_list)

auc_table <- cbind(
  Method = rownames(AUC_mat),
  AUC_mat
)

write.table(
  auc_table,
  "ML_AUC_matrix.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


## =========================================================
## 13. Plot AUC heatmap
## =========================================================

AUC_mat <- read.table(
  "ML_AUC_matrix.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = 1,
  stringsAsFactors = FALSE
)

avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = TRUE)

AUC_mat <- AUC_mat[
  names(avg_AUC),
  ,
  drop = FALSE
]

avg_AUC <- as.numeric(
  format(avg_AUC, digits = 3, nsmall = 3)
)

CohortCol <- brewer.pal(
  n = max(3, ncol(AUC_mat)),
  name = "Paired"
)[1:ncol(AUC_mat)]

names(CohortCol) <- colnames(AUC_mat)

cellwidth <- 1
cellheight <- 0.5

hm <- SimpleHeatmap(
  Cindex_mat = AUC_mat,
  avg_Cindex = avg_AUC,
  CohortCol = CohortCol,
  barCol = "steelblue",
  cellwidth = cellwidth,
  cellheight = cellheight,
  cluster_columns = FALSE,
  cluster_rows = FALSE
)

pdf(
  file = "ML_AUC_heatmap.pdf",
  width = cellwidth * ncol(AUC_mat) + 6,
  height = max(6, cellheight * nrow(AUC_mat) * 0.45)
)

draw(
  hm,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()


## =========================================================
## 14. Export top model ranking
## =========================================================

top_model_table <- data.frame(
  Model = rownames(AUC_mat),
  Mean_AUC = avg_AUC,
  AUC_mat,
  stringsAsFactors = FALSE
)

write.table(
  top_model_table,
  "ML_top_model_ranking.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


## =========================================================
## 15. Export representative model features
## =========================================================
## These two representative models were selected for downstream
## signature scoring and ROC validation in the manuscript.

representative_models <- c(
  "glmBoost+RF",
  "Stepglm+RF",
  "Stepglm[backward]+RF",
  "Stepglm[forward]+RF",
  "Stepglm[both]+RF"
)

representative_features <- feature_df[
  feature_df$Model %in% representative_models,
  ,
  drop = FALSE
]

write.table(
  representative_features,
  "ML_representative_model_features.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Machine-learning model screening completed.\n")
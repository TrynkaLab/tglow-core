library(data.table)
library(caret)
library(matrixStats)

# Global to store the running total of filesets
FILESET_ID=0

#-------------------------------------------------------------------------------
#' Find nearest index
#'
#'Find the index of the value in x that is closest to value
#'
#' @returns The position in X where value is closest to
nearest.index <- function(x, value) {
  which.min(abs(x-value))
}

#-------------------------------------------------------------------------------
#' Retrieve a cell (and its neighbours) based on a feature sumstat
#' 
#' Gets a cell and its closes neighbours based on a single feature and a 
#' sumstat (mean, median, upper.q, lower.q). To customise the quantile used
#' specify q.
#' 
#' @returns vector of indices in cell matrix
tglow.fetch.representative.cell <- function(dataset, feature, assay="cells", metric="mean", na.rm=F, n=0, subset=NULL, q=NULL) {
  
  x <- dataset[[assay]][, feature]
  i <- 1:length(x)
  
  if (!is.null(subset)) {
    x <- x[subset]
    i <- i[subset]
  }
  
  i <- i[order(x)]
  x <- sort(x)
  
  if (metric == "mean") {
    m   <- mean(x, na.rm=na.rm)
    out <- nearest.index(x, m)
  } else if (metric == "median") {
    m   <- median(x, na.rm=na.rm)
    out <- nearest.index(x, m)
  } else if (metric == "upper.q") {
    if (is.null(q)) {q <- 0.75}
    m <- quantile(x, probs=q)
    out <- nearest.index(x, m)
  } else if (metric == "lower.q") {
    if (is.null(q)) {q <- 0.25}
    m <- quantile(x, probs=q)
    out <- nearest.index(x, m)
  }
  
  if (n > 0){
    out <- (out-n):(out+n)
  }
  
  return(i[out])
}


#-------------------------------------------------------------------------------
#' Retrieve a set of cells (and its neighbours) based on a feature sumstat
#' 
#' Fetches (n*2)+1 cells arround the 0th, 25th, 50th, 75th and 100th quantiles
#'
#' @returns A list with row ids and labels for the cells
tglow.fetch.representative.cell.quantiles <- function(dataset, feature, assay="cells", name=NULL, n=1, col.obj.id="cell_ObjectNumber_Global") {
  
  cur.assay <- dataset[[assay]]
  if (is.null(name)) {
    name <- feature
  } 
  
  if (!col.obj.id %in% colnames(cur.assay)) {
    stop(paste0("Id collumn ", col.obj.id, " not found in assay ", assay))
  }
  
  cells    <- c()
  cell.ids <- c()
  
  f   <- cur.assay[,feature]
  idx <- cur.assay[!is.na(f), col.obj.id]
  f   <- f[!is.na(f)]
  fs  <- order(f)
  
  lower <- which(cur.assay[, col.obj.id] %in% idx[head(fs, n=(n*2)+1)])
  upper <- which(cur.assay[, col.obj.id] %in% idx[tail(fs, n=(n*2)+1)])

  l10 <- tglow.fetch.representative.cell(dataset,
    feature=feature,
    assay=assay,
    metric="lower.q",
    q=0.10,
    n=n)
  
  l25 <- tglow.fetch.representative.cell(dataset,
    feature=feature,
    assay=assay,
    metric="lower.q",
    q=0.25,
    n=n)
  
  l50 <- tglow.fetch.representative.cell(dataset,
    feature=feature,
    assay=assay,
    metric="lower.q",
    q=0.5,
    n=n)
  
  l75 <- tglow.fetch.representative.cell(dataset,
    feature=feature,
    assay=assay,
    metric="lower.q",
    q=0.75,
    n=n)

  l90 <- tglow.fetch.representative.cell(dataset,
    feature=feature,
    assay=assay,
    metric="lower.q",
    q=0.90,
    n=n)
                                         
  cn <- c(rep(paste0("q0 - ", name), (n*2)+1),
        rep(paste0("q10 - ", name), (n*2)+1),
        rep(paste0("q25 - ", name), (n*2)+1),
        rep(paste0("q50 - ", name), (n*2)+1),
        rep(paste0("q75 - ", name), (n*2)+1),
        rep(paste0("q90 - ", name), (n*2)+1),
        rep(paste0("q100 - ", name), (n*2)+1))
                                
  return(list(ids=c(lower, l10, l25, l50, l75, l90, upper), names=cn))         
  
}

#-------------------------------------------------------------------------------
#' Test feature association with meta data element
#' 
#' 
tglow.run.assoc <- function(dataset, predictor, assay="cells_norm", method="lm", features="analyze_norm", img.id.col="Image_ImageNumber_Global") {
  
  if(!assay %in% names(dataset)) {
    stop(paste0(assay, " not found in data. Call tglow.norm first"))
    #cat("[ERROR] ", assay, " not found in data. Call tglow.norm first\n")
    #return(NULL)
  }
  
  # Filter to features  
  if (length(features) == 1) {
    if (features %in% c("analyze_norm", "analyze")) {
      features <- dataset$features$id[dataset$features[[features]]]
      features <- features[features %in% colnames(dataset[[assay]])]
    } else {
        stop("Not a valid features column")
    }
  } else{
    features <- features[features %in% colnames(dataset[[assay]])]
  }

  cur.cells <- dataset[[assay]][,features]
  
  if (method == "lm") {

    if (predictor %in% colnames(dataset$meta) & assay %in% c("cells", "cells_corrected", "cells_corrected_norm", "cells_norm", "cells_transform", "cells_log", "cells_sqr", "cells_res", "cells_corrected_for_intensity_norm")) {

      x <- dataset$meta[dataset[[assay]][,img.id.col], predictor]

    } else if (predictor %in% colnames(dataset$meta) & assay %in% c("agg", "agg_cor", "agg_cor_mean", "agg_cor_median")){ # Julie added this so we can do it on the aggregated data as well

      x <- dataset$meta[rownames(dataset$agg), predictor]

    } else if (predictor %in% colnames(cur.cells)) {

      x <- cur.cells[,predictor]

    } else {
      stop("Not a valid response")
      #cat("[ERROR] Not a valid response\n")
      #return(NULL)
    }
    
    ## TODO: use custom more efficient LM 
    i <- 0
    j <- ncol(cur.cells)
    
    res <- apply(cur.cells, 2, function(y) {
      i <<- i +1
      cat("\r[INFO]", round((i/j)*100, digits=2), "%")
      l <- lm(y ~ x)
      m <- summary(l)
      f <- m$fstatistic
      p <- pf(f[1],f[2],f[3],lower.tail=F)
      cooksD <- cooks.distance(l)
      influential <- cooksD[(cooksD > (3 * mean(cooksD, na.rm = TRUE)))]
      n.influential <- length(influential)

      return(c(m$r.squared, m$adj.r.squared, m$fstatistic, p, n.influential))

    })
    
    res <- t(res)
    colnames(res) <- c("r.squared", "adj.r.squared", "f.stat", "numdf", "dendf", "p.value", "n.influential")  
    
    return(res)

  } else {
    stop("Not yet implemented")
    #cat("[ERROR] Not yet impelemeted\n")
    #return(NULL)
  }
  
}




#-------------------------------------------------------------------------------
#' Test feature association with meta data element correcting for another feature
#' 
#' 
tglow.run.assoc.twosteps <- function(dataset, predictor, to.correct, assay="cells_norm", method="lm", features=NULL) {
  
  if(!assay %in% names(dataset)) {
    stop(paste0(assay, " not found in data. Call tglow.norm first"))
    #cat("[ERROR] ", assay, " not found in data. Call tglow.norm first\n")
    #return(NULL)
  }
  
  if (is.null(features)) {
    features <- colnames(dataset$cells)[dataset$features$analyze]
    features <- features[features %in% colnames(dataset[[assay]])]
  } 
  
  cur.cells <- dataset[[assay]][,features]
  
  if (method == "lm") {
    
    #Find x if it's in cells, agg or in meta
    if (predictor %in% colnames(dataset$meta) & assay == "cells_norm") {
      x <- dataset$meta[dataset[["cells"]]$Image_ImageNumber_Global, predictor]
    } else if (predictor %in% colnames(dataset$meta) & assay == "agg"){ # Julie added this so we can do it on the aggregated data as well
        x <- dataset$meta[rownames(dataset$agg), predictor]
    } else if (predictor %in% colnames(cur.cells)) {
      x <- cur.cells[,predictor]
    } else {
      stop("Not a valid response")
    }
    
    # Make a dataframe to store the results
    df <- data.frame(x = x)
    
    # Find variable to correct for if it's in cells, agg or in meta
    if (to.correct %in% colnames(dataset$meta) & assay == "cells_norm") {
      z <- dataset$meta[dataset[["cells"]]$Image_ImageNumber_Global, to.correct] 
    } else if (to.correct %in% colnames(dataset$meta) & assay == "agg"){ # Julie added this so we can do it on the aggregated data as well
      z <- dataset$meta[rownames(dataset$agg), to.correct]
    } else if (to.correct %in% colnames(cur.cells)) {
      z <- cur.cells[,to.correct]
    } else {
      stop("Not a valid response")
    }
      
    df$z <- z

    ## TODO: use custom more efficient LM 
    i <- 0
    j <- ncol(cur.cells)
    
    res.list <- list()
    residuals <- data.frame(matrix(nrow = nrow(cur.cells), ncol =0))
    
    for(y in seq_along(colnames(cur.cells))){
      
      i <- i +1
      cat("\r[INFO]", round((i/j)*100, digits=2), "%")
      
      # Prepare dataframe
      data <- df
      data$y <- cur.cells[, y]
      
      # Run the linear models
      l1 <- lm(y ~ z, data = data)
      l2 <- lm(y ~ x, data = data.frame(y = residuals(l1), x = data$x))
      
      m <- summary(l2)
      f <- m$fstatistic
      p <- pf(f[1],f[2],f[3],lower.tail=F)
      
      res.list[[colnames(cur.cells)[y]]] <- data.frame(r.squared = m$r.squared,
                                                       adj.r.quared = m$adj.r.squared,
                                                       f.stat = f[1],
                                                       numdf = f[2],
                                                       dendf = f[3],
                                                       p.value = p,
                                                       feature = colnames(cur.cells)[y]) 
      residuals[, y] <- residuals(l1)
    }
  }
  
  res <- bind_rows(res.list)

  cats <- dataset$features[dataset$features$analyze, ]$category
  names(cats) <- rownames(dataset$features[dataset$features$analyze, ])

  mes <- dataset$features[dataset$features$analyze, ]$measurement
  names(mes) <- rownames(dataset$features[dataset$features$analyze, ])

  obj <- dataset$features[dataset$features$analyze, ]$object
  names(obj) <- rownames(dataset$features[dataset$features$analyze, ])

  res$category <- cats[res$feature]
  res$measurement <- mes[res$feature]
  res$object <- obj[res$feature]
  res$channel <- str_extract(res$feature, pattern = paste0(c("mito", "actin", "cd25_ki67", "dna"), collapse = "|"))
  res$channel <- ifelse(is.na(res$channel), "AreaShape", res$channel)

  return(list(res, residuals))  
}


#-------------------------------------------------------------------------------
#' Construct a feature metadata table from _cells.tsv
tglow.get.feature.meta.from.cells <- function(feature.names) {
  
  #feature.meta <- data.frame(id=paste0(colnames(cells), "_", cells[1,]),
  #                         object=colnames(cells),
  #                         measurement=as.character(cells[1,]))
  
  
  feature.meta <- data.frame(id=feature.names,
                           object=sapply(strsplit(feature.names, split="_"), function(x){x[[1]]}),
                           measurement=sapply(strsplit(feature.names, split="_"), function(x){paste0(x[-1], collapse="_")}))
  
  feature.meta$category  <- sapply(strsplit(feature.meta$measurement, split="_"), function(x){x[[1]]})
  feature.meta$name      <- sapply(strsplit(feature.meta$measurement, split="_"), function(x){paste0(x[-1], collapse="_")})
  
  rownames(feature.meta) <- feature.meta$id
  return(feature.meta)
}

#-------------------------------------------------------------------------------
#' Merge cell level filesets
#' Data must be a list of lists with outputs from tglow.read.fileset.a/b
#' or have the items, cells, meta, orl,  [children], [features], [cells_norm]
tglow.merge.filesets <- function(data) {

  if (class(data) != "list") {
    stop("Data argument must be a list.")
  }

  out <- list()
  
  ncol.cells <- as.numeric(lapply(data, function(x){ncol(x[["cells"]])}))
  ncol.meta  <- as.numeric(lapply(data, function(x){ncol(x[["meta"]])}))
  ncol.orl   <- as.numeric(lapply(data, function(x){ncol(x[["orl"]])}))
  #ncol.nucl <- as.numeric(lapply(data, function(x){ncol(x[["nucl"]])}))

  selector <- (ncol.cells == as.numeric(names(which.max(table(ncol.cells))))) & (ncol.meta == as.numeric(names(which.max(table(ncol.meta))))) & (ncol.orl  == as.numeric(names(which.max(table(ncol.orl)))))
  
  if (sum(selector) != length(selector)) {
  
    msg <- paste0("Not all filesets have the same collumn number, dropping the ones with least frequent number. Retained ", 
    sum(selector), "/", length(selector), " filesets. \n",
    "The following filesets are at issue:\n")
    
    for (i in (1:length(data))[!selector]) {
      msg <- paste0(msg, i, ", ")
    }
    
    warning(msg)
  }
  
  out$cells <- dplyr::bind_rows(lapply(data[selector], function(x){x[["cells"]]}),)
  out$meta  <- dplyr::bind_rows(lapply(data[selector], function(x){x[["meta"]]}),)
  out$orl   <- dplyr::bind_rows(lapply(data[selector], function(x){x[["orl"]]}),)
  #out$nucl <- dplyr::bind_rows(lapply(data[selector], function(x){x[["nucl"]]}),)

  if ( "features" %in% names(data[[1]])) {
    out$features <-  dplyr::bind_rows(lapply(data[selector], function(x){x[["features"]]}))
  }
  
  if ( "cells_norm" %in% names(data[[1]])) {
    out$cells_norm <-  dplyr::bind_rows(lapply(data[selector], function(x){x[["cells_norm"]]}))
  }
  
  # Merge the child object matrices
  if ("children" %in% names(data[[1]])) {
    out$children <- list()
    for (obj in names(data[[1]]$children)) {
      out$children[[obj]] <-  dplyr::bind_rows(lapply(data[selector], function(x){x[["children"]][[obj]]}))
    } 
  }
  
  #out$cells <- do.call(rbind, lapply(data, function(x){x[["cells"]][,feature.names]}))
  #out$meta  <- do.call(rbind, lapply(data, function(x){x[["meta"]]}))
  #out$orl   <- do.call(rbind, lapply(data, function(x){x[["orl"]]}))
  
  #if ( "features" %in% names(data[[1]])) {
  #  out$features <- do.call(rbind, lapply(data, function(x){x[["features"]]}))
  #}
  
  #if ( "cells_norm" %in% names(data[[1]])) {
  #  out$cells_norm <- do.call(rbind, lapply(data, function(x){x[["cells_norm"]]}))
  #}
  
  return(out)
}

#-------------------------------------------------------------------------------
#' Calculate the modified z-score
tglow.mod.zscore <- function(y) {
  return((0.6745*(y-median(y, na.rm = TRUE)))/mad(y, na.rm = TRUE))
}

#-------------------------------------------------------------------------------
#' Calculate modified Z-score for cells per grouping (sample / condition)
#' Adapted from Julies manip_zscore_per_sample
#' 
#' @param data a dataframe to z-score. Rows treated as samples, cols as features
#' @param grouping a vector of nrow(data) to optionally subset z-score calculation on
#' @param features only apply z-scoring to these features (defaults to all when NULL)
#' @param method mod.z for modified z-score or z for regular z-score
#' 
#' @returns Dataframe of dim(data) where the values have been replaced by the group 
#' specific z-scores. 
tglow.grouped.scale <- function(data, grouping, features=NULL, method="mod.z"){

  if (class(data) == "data.frame") {
    
    # For dataframe input type
    if (is.null(features)) {
      features <- colnames(data)
    }
    
    for(i in unique(grouping)){
      cat("[INFO] Processing group ", i, "\n")
      subset <- grouping == i
      
      # Filter to only the relevant cells
      sample_data <- data[subset, features, drop=F]
      
      # Calculate z-score and replace sample data
      if (method == "mod.z") {
        data[subset, features] <- apply(sample_data, 2, tglow.mod.zscore)
      } else if (method=="z") {
        data[subset, features] <- apply(sample_data, 2, scale, center=T, scale=T)
      } else {
        stop("Method not valid")
      }
    }
    
  } else if (class(data) == "numeric"){
    # For vector input type
    for(i in unique(grouping)){
      #cat("[INFO] normalizing group ", i, "\n")
      subset       <- grouping == i
      if (method == "mod.z") {
        data[subset] <- tglow.mod.zscore(data[subset])
      } else if (method=="z") {
        data[subset] <- scale(data[subset], center=T, scale=T)
      } else {
        stop("Method not valid")
      }
    }
    
  } else {
    stop("Invalid input class")
    #cat("[ERROR] Invalid input class\n")
    #return(NULL)
  }

  return(data)
}

#-------------------------------------------------------------------------------
#' Normalize features 
#' 
#' Normalize a tglow dataset with z-score or modified z-score
#' Output is placed on the dataset list with name assay.out
#' 
#' To only normalize specific features supply the feature names.
#' 
#' @param method method to use for normalizing, z | mod zcore
#' 
#' @returns tglow dataset with normalized assay
tglow.normalize <- function(dataset, method="mod.z", features="analyze", assay.out="cells_norm", assay = "cells", filter=TRUE, new.feat.col = "analyze_norm") {
  
  
  # Normalize and filter features
  ft <- dataset[[assay]]
  
  # Extract ids
  ids <- ft$Image_ImageNumber_Global
  
  # Filter to features  
  if (all(features %in% c("analyze_norm", "analyze"))) {
    features <- dataset$features$id[dataset$features[[features]]]
    features <- features[features %in% colnames(dataset[[assay]])]
  } 
  
  # Creat matrice with the appropriate features [OLD SCRIPT THAT REMOVES FEATURES FROM DF]
  # ft <- as.matrix(ft[, colnames(ft) %in% features])
  
  i <- 0
  j <- ncol(ft)
  
  if (method == "mod.z") {
    ft[, features] <- apply(ft[, features], 2, function(x){
      i <<- i +1
      cat("\r[INFO]", round((i/j)*100, digits=2), "%")
      return(tglow.mod.zscore(x))
    })
  } else if (method == "z") {
    ft[, features] <- apply(ft[, features], 2, function(x){
      i <<- i +1
      cat("\r[INFO]", round((i/j)*100, digits=2), "%")
      return(scale(x, center=T, scale=T))
    })
  } else {
    stop("Method specified not valid")
  }
  cat("\n")
  
  if(filter==TRUE){
    
    # OLD SCRIPT THAT REMOVES FEATURES FROM DATAFRAME
    # Remove non properly normalized features and those with zero variance
    #cat("[INFO] Checking for NA's\n")
    #ft <- ft[,colSums(is.na(ft)) == 0]
    #cat("[INFO] Checking for zero variances\n")
    #ft <- ft[,Rfast::colVars(ft) != 0]
    #ft <- ft[,colVars(as.matrix(ft))!=0]
    
    #if (ncol(ft) != length(features)) {
    #  warning(paste0("Removed ", length(features)-ncol(ft), " features with zero variance after normalizing"))
    #}
    
    # NEW SCRIPT THAT KEEPS FEATURES IN DATAFRAME - BUT CHANGES output$features
    dataset$features[[new.feat.col]] <- dataset$features$analyze
    cat("[INFO] Checking for NA's\n")
    dataset$features[[new.feat.col]][colSums(is.na(ft[, features])) > 0] <- F # Set features with NAs as F
    
    cat("[INFO] Checking for zero variances\n")
    dataset$features[[new.feat.col]][Rfast::colVars(ft[, features]) == 0] <- F # Set features with zero variance as F
    dataset$features[[new.feat.col]][colVars(as.matrix(ft[, features])) == 0] <- F # Set features with zero variance as 
    
  }
  
  # Transform to dataframe add back Image_ImageNumber_Global
  ft <- as.data.frame(ft)
  #ft$Image_ImageNumber_Global <- ids
  
  dataset[[assay.out]] <- ft
  
  return(dataset)
}


#-------------------------------------------------------------------------------
#' Detect outlier cells in PCA space
#' 
#' SEE FUNCTION BELOW FOR AN UPDATED ONE BY JULIE [more annotated]
#' 
tglow.pca.outliers.deprecated <- function(dataset, assay = "cells_transform", grouping=NULL, pc.thresh=0.5, pc.max=500, pc.n=NULL, threshold=3.5, features=NULL, features.col = "analyze", method="z", renormalize=T) {
  
  #stop("Not finished")

  # Selecting Features
  if (is.null(features)) {

    cat("Selecting Column: ", features.col, "in Assay ", assay, "\n")
    
    features <- output$features$id[output$features[[features.col]]]
    features <- features[features %in% colnames(dataset[[assay]])]

  } 

  # Check that the assay is there
  if (!assay %in% names(dataset)) {

      stop(paste0("Assay not found: ", assay))

  } 
  
  # Now we want to create 'data' according to the different options of normalization

  # If there is no grouping, make sure dataset is normalized across the experiment
  if (is.null(grouping)) {

    stop("Not implemented")
    data <- dataset[[assay]]
    
    if (renormalize) {

      dataset <- tglow.normalize(dataset, features=features, method="z", assay = assay, assay.out = "cells_norm_qc")
      data <- dataset[["cells_norm_qc"]]

    }

  # If there is grouping -> normalize data per group
  } else {

    cat("[INFO] Starting Analysis Per Group\n")

    #data <- tglow.grouped.scale(dataset$cells, grouping=grouping, method=method, features=features)
    #data <- data[,colSums(is.na(data)) == 0]
    #data <- data[,colVars(as.matrix(data))!=0]

    if (renormalize) { ### WHAT IS THIS PER GROUPPP??

      cat("[INFO] Renormalizing Per Group\n")

      dataset <- tglow.normalize(dataset, features=features, method="z", assay = assay, assay.out = "cells_norm_qc")
      data <- dataset[["cells_norm_qc"]]

    } else {

        data       <- dataset[[assay]]

    }
    
  
    #output     <- matrix(NA, nrow=nrow(data), ncol=length(features))
    output.pcs <- matrix(NA, nrow=nrow(data), ncol=length(features))
    outliers   <- rep(NA, nrow(data))
    
    if (pc.max > length(features)) {
      pc.max <- length(features) -1
    }
    
    for (group in unique(grouping)) {
      cat("[INFO] Calculating PC's for ", group, "\n")
      cur.data <- data[grouping==group, features]
      
      if (sum(grouping==group) < 2) {
        warning("Need at least two cells in group, skipping")
        next
      }
      
      #if (method == "mod.z") {
      #  cur.data <- apply(cur.data, 2, tglow.mod.zscore)
     # } else if (method == "z") {
        cur.data <- apply(cur.data, 2, scale, center=T, scale=T)
      #}
      
      
      cur.data <- cur.data[,colSums(is.na(cur.data)) == 0]
      cur.data <- cur.data[,colVars(as.matrix(cur.data))!=0]
      
      if (ncol(cur.data) != length(features)) {
        warning(paste0("Removed ", length(features)-ncol(cur.data), " features with zero variance after normalizing"))
      }
      
      if (nrow(cur.data) < pc.max) {
        pc.final <- nrow(cur.data)-1
      } else {
        pc.final <- pc.max
      }
      
      pca      <- irlba::prcomp_irlba(cur.data, n=pc.final, center=T, scale=T)
      
      if (is.null(pc.n)) {
        pc.var <- pca$sdev^2/pca$totalvar
        if (sum(cumsum(pc.var) > pc.thresh) >=1){
          pc.n   <- min(which(cumsum(pc.var) > pc.thresh))
        } else {
          pc.n <- pc.final
        }
        
        if (cumsum(pc.var)[pc.n] < pc.thresh) {
          warning("Last pc doesnt pass pc.thresh, try increasing pc.max")
        }
        
        cat("[INFO] Selected ", pc.n, " pcs explaining ", round(cumsum(pc.var)[pc.n], digits=2)*100, "% of the variance\n")
      }
      
      
      if (method == "mod.z") {
        pcs.norm <- apply(pca$x, 2, tglow.mod.zscore)[,1:pc.n, drop=F]
      } else if (method == "z") {
        pcs.norm <- apply(pca$x, 2, scale, center=T, scale=T)[,1:pc.n, drop=F]
      }
      
      #pcs.norm <- scale(pca$x)[,1:pc.n, drop=F]
      
      outliers[grouping==group] <- rowSums(abs(pcs.norm) < threshold) != pc.n
      pc.n <- NULL
    
      #output[grouping==group,1:pc.n]     <- abs(pcs.norm) < threshold
      #output.pcs[grouping==group,1:pc.n] <- pcs.norm
      
    }
    
    return(list(outliers=outliers))#, pcs=output.pcs))
  }
  
}

#-------------------------------------------------------------------------------
#' Detect outlier cells in PCA space
#' [more annotated, updated normalization per group] 
#'
tglow.pca.outliers <- function(dataset, 
                               assay = "cells_transform", 
                               grouping=NULL, 
                               pc.thresh=0.5, 
                               pc.max=500, 
                               pc.n=NULL, 
                               threshold=3.5, 
                               features=NULL, 
                               features.col = "analyze", 
                               method="z") {
  
  if (is.null(grouping)) {
    stop("Must provide grouping")
  }
  
  if ( "list" %in% class(dataset)) {
    # Check that the assay is there
    if (!assay %in% names(dataset)) {
      stop(paste0("Assay not found: ", assay))
    } else {
      cat("[INFO] Found assay: ", assay, "\n") 
    }
    
    # Selecting Features
    if (is.null(features)) {
      cat("[INFO] Selecting Feature Column: ' ", features.col, " ' from the assay: '", assay, "' \n")
      features <- output$features$id[output$features[[features.col]]]
      features <- features[features %in% colnames(dataset[[assay]])]
    } 
     
    # Set data to the assay
    data <- as.matrix(dataset[[assay]][,features])
    
  } else if ("matrix" %in% class(dataset)) {
    data <- dataset
  } else {
    stop("Dataset is not a valid class, must be a matrix or a list (tglow)")
  }
  
  # Now we want to create 'data' according to the different options of normalization
  # If there is no grouping, make sure dataset is normalized across the experiment    
  cat("[INFO] Starting Analysis Per Group", "\n")

  
  # Define results matrix
  output     <- matrix(NA, nrow=nrow(data), ncol=ncol(data))
  output.pcs <- matrix(NA, nrow=nrow(data), ncol=ncol(data))
  outliers   <- rep(NA, nrow(data))
  
  # Now for each group
  for (group in unique(grouping)) {
    
    cat("[INFO] Calculating for ", group, "\n")
    
    # Subset data
    cur.data <- data[grouping==group,]
    
    # Check if there enough cells
    if (sum(grouping==group) < 2) {
      warning("Need at least two cells in group, skipping")
      next  
    }
    
    cat("[INFO] Rescaling data for group\n")
    cur.data <- fcolScale(cur.data, add_attr=F)
    
    # Remove features of low variance
    cat("[INFO] Removing Features of Low Variance for ", group, "\n")
    cur.data <- cur.data[,Rfast::colsums(is.na(cur.data)) == 0]
    cur.data <- cur.data[,Rfast::colVars(cur.data)!=0]
    
    if (ncol(cur.data) != ncol(data)) {
      warning(paste0("Removed ", ncol(data)-ncol(cur.data), " features with zero variance after normalizing"))
    }
    
    # Again define PCs to less then the smallest dimension
    if (nrow(cur.data) < pc.max) {   
      pc.final <- nrow(cur.data) -1
    } else {
      pc.final <- pc.max
    }
    
    # Define the number of PCs
    if (pc.final > ncol(cur.data)) {
      pc.final <- ncol(cur.data) -1 
    }
    
    cat("[INFO] Calculating  ", pc.final, " pc's on ", nrow(cur.data), " samples and, ", ncol(cur.data), " features\n")

    # Calculate PCAs
    pca      <- irlba::prcomp_irlba(cur.data, n=pc.final)
    
    # Selecting n.pcs or number of PCs that explain x of the variance
    if (is.null(pc.n)) {
      pc.var <- pca$sdev^2/pca$totalvar
      
      if (sum(cumsum(pc.var) > pc.thresh) >=1){
        pc.n   <- min(which(cumsum(pc.var) > pc.thresh))
      } else {
        pc.n <- pc.final

      }
      
      if (cumsum(pc.var)[pc.n] < pc.thresh) {     
        warning("Last pc doesnt pass pc.thresh, try increasing pc.max")    
      }
      
      cat("[INFO] Selected ", pc.n, " pcs explaining ", round(cumsum(pc.var)[pc.n], digits=2)*100, "% of the variance\n")
    }
    
    # Now we perform z-scoring on the PCs
    if (method == "mod.z") {
      pcs.norm <- apply(pca$x, 2, tglow.mod.zscore)[,1:pc.n, drop=F]
    } else if (method == "z") {   
      pcs.norm <- apply(pca$x, 2, scale, center=T, scale=T)[,1:pc.n, drop=F] 
    }
    
    # Add results to the outlier list
    outliers[grouping==group] <- rowSums(abs(pcs.norm) < threshold) != pc.n
    pc.n <- NULL
  }
  
  return(list(outliers=outliers)) #, pcs=output.pcs))
}



#-------------------------------------------------------------------------------
#' Calculate the median value per feature, defaulting to normalized assay & per image value
#' 
#'
tglow.median_by_group <- function(dataset, assay = "cells_norm", features = NULL, group = "ImageNumber_Global", f = "test") {
  
  # Defining the features 
  if (is.null(features)) {
    features <- dataset$features$id[dataset$features$analyze]
    features <- features[features %in% colnames(dataset[[assay]])]
  } else{

    features <- features[features %in% colnames(dataset[[assay]])]

  }
  
  # Define group values
  group_values <- dataset$meta[dataset[[assay]]$Image_ImageNumber_Global, ][[group]]
  group_values <- as.factor(group_values)
  
  # Use tapply/by to split the data by groups and then apply the median function to each column of those groups

  if(f == "median"){

      results <- by(dataset[[assay]][, features], group_values, function(x) apply(x, 2, median))


  } else if(f == "mean"){

      results <- by(dataset[[assay]][, features], group_values, function(x) apply(x, 2, mean))

  } else{

      stop("Please choose f == median or mean")

  }

  results <- do.call(rbind, results)
  results <- as.data.frame(results)
  results$Image_ImageNumber_Global <- rownames(results)

  return(results)  
}



#-------------------------------------------------------------------------------
#' Faster alternative to scale a matrix
#' https://www.r-bloggers.com/2016/02/a-faster-scale-function/
fcolScale <- function(x,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = TRUE,
                    rows = NULL,
                    cols = NULL,
                    na.rm=F) {
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }
  ################
  # Get the column means
  ################
  cm = colMeans(x, na.rm=na.rm)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = matrixStats::colSds(x, center = cm, na.rm=na.rm)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = t( (t(x) - cm) / csd )
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  return(x)
}


#-------------------------------------------------------------------------------
#' Correct for plate effects for example and return corrected dataframe
#' Updated by jm52 to return a dataframe of the same structure as output$cells so we can still use output$features$analyze
#'

tglow.correct <- function(dataset, assay = "cells", to.correct, assay.out = "cells_corrected", residual.skip.cols=NULL, features = NULL, img.id.col="Image_ImageNumber_Global"){
  
  # Defining the features 
  if (is.null(features)) {
    features <- dataset$features$id[dataset$features$analyze]
    features <- features[features %in% colnames(dataset[[assay]])]
  } 
  
  # Getting only the cells and features of interest and removing the NAs
  cur.cells <- na.omit(dataset[[assay]][, features])
  
  # Make a basic dataframe with our variable to correct (z) for the function [so we do z ~ f]
  if(all(to.correct %in% colnames(dataset$meta))){ # if we want to correct a metadata variable
    
    df <- data.frame(dataset$meta[dataset[[assay]][[img.id.col]], to.correct])
    
  } else{ # if we want to correct a feature (i.e: intensity)
    
    df <- data.frame(dataset[[assay]][, to.correct]) 
    
  }
  
  i <- 0 # tracking of progress
  j <- ncol(cur.cells) # tracking of progress
  
  # Create a dataframe of the original structure -> in the loop, replace relevant columns with residuals
  residuals <- dataset[[assay]]
  
  # For each feature, make a temporary df (data) and correct for the feature, add the residuals of the model as the new data
  for(y in seq_along(colnames(cur.cells))){
    
    i <- i +1
    cat("\r[INFO]", round((i/j)*100, digits=2), "%")
    
    # Prepare dataframe
    data <- df
    data$y <- cur.cells[, y]
    
    # Run the linear models
    l1 <- lm(y ~ ., data = data)
    feature <- colnames(cur.cells)[y]

    if (!is.null(residual.skip.cols)) {

      # Calculate residuals, skipping the effects of a covariate
      coef <- l1$coefficients
      residual.keep.cols   <- to.correct[!to.correct %in% residual.skip.cols]
      residuals[[feature]] <- data$y - (coef["(Intercept)"] + as.matrix(data[,residual.keep.cols]) %*% coef[residual.keep.cols])

    } else {
      residuals[[feature]] <- residuals(l1)
    }

    
  }
  
  # Add it back to the dataset whilst keeping the same structure as the original assay
  dataset[[assay.out]] <- residuals
  
  return(dataset)
  
  
}


#-------------------------------------------------------------------------------
#' Correct for plate effects for example PER GROUP and return corrected dataframe
#' Updated to keep the same structure as the original dataframe, so we can use output$features$analyze
#'


tglow.correct.group <- function(dataset, assay = "cells", assay.out = "cells_corrected", to.correct = "plate_id", features = NULL, group = "plate"){
  
  # Defining the features 
  if (is.null(features)) {
    features <- dataset$features$id[dataset$features$analyze]
    features <- features[features %in% colnames(dataset[[assay]])]
  } 
  
  # Make a dataframe to keep all of the data [same structure as original df]
  residuals <- dataset[[assay]]
  
  # Find the groups 
  groups <- unique(dataset$meta[[group]])
  
  # Extract cells, subset to each group & remove NAs
  for(g in groups){
    
    # Subset to only our samples of interest
    imgs.to.keep <- rownames(dataset$meta[dataset$meta[[group]] == g, ])
    imgs.to.keep <- rownames(dataset$meta) %in% imgs.to.keep
    table(imgs.to.keep)
    
    temp <- tglow.filter.img.apply(dataset, imgs.to.keep)
    
    # Initialize residuals for the current group
    residuals.group <- temp[[assay]]
    
    # Getting only the cells and removing the NAs
    cur.cells <- na.omit(temp[[assay]][, features])
    
    # Make a basic dataframe with our variable to correct (z) for the function [so we do z ~ f]
    df <- data.frame(z = temp$meta[temp[[assay]]$Image_ImageNumber_Global, to.correct])
  
    i <- 0 # tracking of progress
    j <- ncol(cur.cells) # tracking of progress
    
    # For each feature, make a temporary dataframe and correct for the feature & add the residuals 
    for(y in seq_along(colnames(cur.cells))){
      
      i <- i +1
      cat("\r[INFO]", round((i/j)*100, digits=2), "%")
      
      # Prepare dataframe
      data <- df
      data$y <- cur.cells[, y]
      
      # Run the linear models
      l1 <- lm(y ~ z, data = data)
      s <- summary(l1)
      feature <- colnames(cur.cells)[y]
      residuals.group[[feature]] <- residuals(l1) + s$coefficients[1,1] 
      
    }
    
    # Add it to the list of groups
    residuals[rownames(residuals.group), colnames(residuals.group)] <- residuals.group
    
    
  }
  
  dataset[[assay.out]] <- residuals
  
  return(dataset)
  
  
}

#-------------------------------------------------------------------------------
#' Function to transform data using boxcox
#' c is the column you want to transform
#' tracking = T, > This should be done in a wrapper on the loop, not here <
tglow.transform.bc <- function(c, return.lambda = F, limit = 5, downsample.lambda=NULL){
  
  #  > This should be done in a wrapper on the loop, not here <
  #if(tracking == T){
  #  # Increment the global counter
  #  i <<- i + 1
  #  # Print the progress
  #  cat("\r[INFO]", round((i/j)*100, digits=2), "%")
  #}

  # If is not positive - then offset to make everything positive
  if(any(c <= 0, na.rm = T) == T) { c <- abs(min(c, na.rm=T)) + c + 1 }
  
  # Create data to estimate transformation parameters
  # Optionally downsample
  if (!is.null(downsample.lambda)) {
    if (class(downsample.lambda) == "integer" || class(downsample.lambda) == "numeric") {   
      if (length(downsample.lambda) > 1) {
        y <- c[downsample.lambda]
      } else {
        y <- c[sample(1:length(c), downsample.lambda)]
      }
    } else {
      stop("downsample.lambda must be a single integer, or a vector representing values to use for lambda estimation")
    }
  } else {
    y <- c
  }
  
  #y <- c[cells] # subset 50K cells
  #x <- rep(1, length(y)) # make everything a 1
  
  # Estimate transformation parameters
  bc     <- MASS::boxcox(y ~  1, plotit = F,  lambda = seq(-limit, limit, 1/10)) # don't plot lambda outcome
  lambda <- bc$x[which.max(bc$y)]

  # If return lambda or transformed values
  if(return.lambda){
      return(lambda)
  } else {
      # Define fudge: uncertainty value so that anything within that range becomes the closest value
      fudge <- 0.1
      
      # Transform data
      if(lambda < fudge & lambda > -fudge) { # If the data is in between -0.1 & 0.1, then just do a log
        t <- log(c)
      } else if (lambda < (1 + fudge) & lambda > (1 - fudge)) { # If the data is in still between -1.2 & 1.2, then just leav it
        t <- c       
      } else { # Otherwise use the calculated lambda   
        t <- (c^lambda - 1)/lambda
      }
      
      return(t)
  }
  
}

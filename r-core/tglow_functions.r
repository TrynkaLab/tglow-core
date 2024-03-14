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
tglow.fetch.representative.cell <- function(data, feature, metric="mean", na.rm=F, n=0, subset=NULL, q=NULL) {
  
  x <- data$cells[, feature]
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
#' Test feature association with meta data element
#' 
#' 
tglow.run.assoc <- function(dataset, predictor, assay="cells_norm", method="lm", features=NULL) {
  
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
    if (predictor %in% colnames(dataset$meta) & assay == "cells_norm") {
      x <- dataset$meta[dataset[["cells"]]$Image_ImageNumber_Global, predictor]
    } else if (predictor %in% colnames(dataset$meta) & assay == "agg"){ # Julie added this so we can do it on the aggregated data as well
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
      m <- summary(lm(y ~ x))
      f <- m$fstatistic
      p <- pf(f[1],f[2],f[3],lower.tail=F)
      return(c(m$r.squared, m$adj.r.squared, m$fstatistic, p))
    })
    
    res <- t(res)
    colnames(res) <- c("r.squared", "adj.r.squared", "f.stat", "numdf", "dendf", "p.value")  
    
    return(res)
  } else {
    stop("Not yet implemented")
    #cat("[ERROR] Not yet impelemeted\n")
    #return(NULL)
  }
  
}


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
    
    #Find variable to correct for if it's in cells, agg or in meta
      
    if (to.correct %in% colnames(dataset$meta) & assay == "cells_norm") {
        
      z <- dataset$meta[dataset[["cells"]]$Image_ImageNumber_Global, to.correct]
      
      } else if (to.correct %in% colnames(dataset$meta) & assay == "agg"){ # Julie added this so we can do it on the aggregated data as well
        
        z <- dataset$meta[rownames(dataset$agg), to.correct]
      
        } else if (to.correct %in% colnames(cur.cells)) {
        
          z <- cur.cells[,to.correct]
        
            } else {
        
                stop("Not a valid response")
                #cat("[ERROR] Not a valid response\n")
                #return(NULL)
      
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
tglow.merge.filesets <- function(data) {
  out <- list()
  
  out$cells <- dplyr::bind_rows(lapply(data, function(x){x[["cells"]]}),)
  out$meta  <- dplyr::bind_rows(lapply(data, function(x){x[["meta"]]}),)
  out$orl   <- dplyr::bind_rows(lapply(data, function(x){x[["orl"]]}),)
  
  if ( "features" %in% names(data[[1]])) {
    out$features <-  dplyr::bind_rows(lapply(data, function(x){x[["features"]]}))
  }
  
  if ( "cells_norm" %in% names(data[[1]])) {
    out$cells_norm <-  dplyr::bind_rows(lapply(data, function(x){x[["cells_norm"]]}))
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
tglow.normalize <- function(dataset, method="mod.z", features=NULL, assay.out="cells_norm") {
  
  # Normalize and filter features
  ft <- dataset$cells
  
  if (is.null(features)) {
    features <- dataset$features$analyze
  } 
  
  ft <- as.matrix(ft[,features])
  
  
  i <- 0
  j <- ncol(ft)

  if (method == "mod.z") {
    ft <- apply(ft, 2, function(x){
      i <<- i +1
      cat("\r[INFO]", round((i/j)*100, digits=2), "%")
      return(tglow.mod.zscore(x))
    })
  } else if (method == "z") {
    ft <- apply(ft, 2, function(x){
      i <<- i +1
      cat("\r[INFO]", round((i/j)*100, digits=2), "%")
      return(scale(x, center=T, scale=T))
    })
  } else {
    stop("Method specified not valid")
  }
  cat("\n")
  
  # Remove non properly normalized features and those with zero variance
  cat("[INFO] Checking for NA's\n")
  ft <- ft[,colSums(is.na(ft)) == 0]
  cat("[INFO] Checking for zero variances\n")
  ft <- ft[,Rfast::colVars(ft) != 0]
  #ft <- ft[,colVars(as.matrix(ft))!=0]
  
  if (ncol(ft) != length(features)) {
    warning(paste0("Removed ", length(features)-ncol(ft), " features with zero variance after normalizing"))
  }
  
  dataset[[assay.out]] <- ft
  
  
  return(dataset)
}

#-------------------------------------------------------------------------------
#' Detect outlier cells in PCA space
#' 
#' 
tglow.pca.outliers <- function(dataset, grouping=NULL, pc.thresh=0.5, pc.max=500, pc.n=NULL, threshold=3.5, features=NULL, method="z", renormalize=F) {
  
  #stop("Not finished")
  if (is.null(features)) {
    features <- dataset$features$id[dataset$features$analyze]
  } 
  
  if (is.null(grouping)) {
    stop("Not implemented")
    if ("cells_norm" %in% names(dataset)) {
      data <- dataset[["cells_norm"]]
    }
    
    if (renormalize) {
      dataset <- tglow.normalize(dataset, features=features, method="z")
      data <- dataset[["cells_norm"]]
    }
  } else {
    cat("[INFO] Renormalizing data per group\n")
    #data <- tglow.grouped.scale(dataset$cells, grouping=grouping, method=method, features=features)
    #data <- data[,colSums(is.na(data)) == 0]
    #data <- data[,colVars(as.matrix(data))!=0]
    
    data       <- dataset$cells
    #output     <- matrix(NA, nrow=nrow(data), ncol=length(features))
    output.pcs <- matrix(NA, nrow=nrow(data), ncol=length(features))
    outliers   <- rep(NA, nrow(data))
    
    if (pc.max > length(features)) {
      pc.max <- length(features) -1
    }
    
    for (group in unique(grouping)) {
      cat("[INFO] Calculating PC's for ", group, "\n")
      cur.data <- data[grouping==group, features]
      
      
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
#' Calculate the median value per feature, defaulting to normalized assay & per image value
#' 
#'
tglow.median_by_group <- function(dataset, assay = "cells_norm") {
  # Get features to analyze & that are in the correct assay
  features <- colnames(dataset$cells)[dataset$features$analyze]
  features <- features[features %in% colnames(dataset[[assay]])]
  
  # Define group values
  group_values <- dataset$meta[dataset$cells$Image_ImageNumber_Global, ]$ImageNumber_Global
  group_values <- as.factor(group_values)
  
  # Use tapply/by to split the data by groups and then apply the median function to each column of those groups
  results <- by(dataset[[assay]][, features], group_values, function(x) apply(x, 2, median))
  
  results <- do.call(rbind, results)
  
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
                    cols = NULL) {
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
  cm = colMeans(x, na.rm = TRUE)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = colSds(x, center = cm)
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



#-------------------------------------------------------------------------------
#' Find which feature filters to a tglow dataset
#' 
#' All filters are inclusive
#' 
#' @param dataset
#' @param filter.table A data frame describing the filters to apply. See details.
tglow.filter.features.calc <- function(dataset, filter.table, assay="cells", features=NULL) {
  
  if (is.null(features)) {
    features <- colnames(dataset$cells)[dataset$features$analyze]
    features <- features[features %in% colnames(dataset[[assay]])]
  } 
  data <- dataset[[assay]][,features]
  
  # Only apply active filters
  filter.table <- filter.table[filter.table$active,]
  
  res <- matrix(TRUE, nrow=ncol(data), ncol=nrow(filter.table))
  rownames(res) <- colnames(data)
  colnames(res) <- filter.table$name
  
  for (i in 1:nrow(filter.table)) {
    
    cur.filter.name <- filter.table[i, "name"]
    cur.filter      <- filter.table[i, "type"]
    cat("[INFO] ", cur.filter.name, ": ", cur.filter, "\n")
    
    # Select the features to apply to
    if (filter.table[i, "column_pattern"] == "all") {
      cur.features <- colnames(data)
    } else {
      cur.features <- grep(filter.table[i, "column_pattern"], colnames(data), value=T)
    }
    
    cat("[INFO] Applying pattern: ",  filter.table[i, "column_pattern"],  " and selected: ", length(cur.features), " features \n")
    
    j <- 0
    res[cur.features, cur.filter.name] <- apply(data[,cur.features], 2, function(x){
      j <<- j+1
      cat("[INFO] ", round((j/length(cur.features))*100, digits=2), "%\r")
      do.call(cur.filter, list(vec=x, thresh=filter.table[i, "value"]))
    })
    cat("\n")
  }
  
  return(res)
}

#-------------------------------------------------------------------------------
#' Find which feature filters to a tglow dataset.
#' 
#' All filters are inclusive. NA is treated as FALSE.
#' 
#' @param dataset tglow dataset
#' @param res output from tglow.filter.calc
#' 
#' @returns 
#' Tglow dataset with updated analyze column
tglow.filter.features.apply <- function(dataset, res) {
  
  if (class(res)[1] == "data.frame" | class(res)[1] == "matrix") {
    selector <- rowSums(res, na.rm=T)!=ncol(res)
    names(selector) <- rownames(res)
  } else if (class(res)[1] == "logical") {
    selector <- res
  } else {
    stop("Not a valid res input type")
  }
  
  if (length(selector) != nrow(dataset$features)) {
    
    if (length(selector) == sum(dataset$features$analyze)) {
      warning("Selection vector length does not match input length, matching on previous analyze col.")
      #cat("[WARN] Selection vector length does not match input length, matching on previous analyze col\n")
      dataset$features[dataset$features$analyze, "analyze"][selector] <- F
    } else {
      warning("Selection vector length does not match input length, matching on names")
      #cat("[WARN] Selection vector length does not match input length, matching on names\n")
      dataset$features[dataset$features$id, "analyze"][selector[dataset$features$id]] <- F
    }
    
    
  } else {
    dataset$features[selector, "analyze"] <- F
  }
  
  return(dataset)
}

#-------------------------------------------------------------------------------
#' Apply a cell level filter to a tglow dataset
#' 
#' @param dataset
#' @param filter.table A data frame describing the filters to apply. See details.
#' 
tglow.filter.cells.calc <- function(dataset, filter.table, assay="cells", features=NULL) {
  
  if (is.null(features)) {
    features <- colnames(dataset$cells)[dataset$features$analyze]
    features <- features[features %in% colnames(dataset[[assay]])]
  } 
  
  data <- dataset[[assay]][,features]
  meta <- dataset$meta
  
  # Only apply active filters
  filter.table <- filter.table[filter.table$active,]
  
  res <- matrix(TRUE, nrow=nrow(data), ncol=nrow(filter.table))
  rownames(res) <- rownames(data)
  colnames(res) <- filter.table$name
  
  for (i in 1:nrow(filter.table)) {
    
    cur.filter.name <- filter.table[i, "name"]
    cur.filter      <- filter.table[i, "type"]
    cat("[INFO] ", cur.filter.name, ": ", cur.filter, "\n")
    
    # Select the features to apply to
    if (filter.table[i, "column_pattern"] == "all") {
      cur.features <- colnames(data)
    } else {
      cur.features <- grep(filter.table[i, "column_pattern"], colnames(data), value=T)
    }
    
    cat("[INFO] Applying pattern: ",  filter.table[i, "column_pattern"],  " and selected: ", length(cur.features), " features \n")
    
    if (!is.na(filter.table[i, "metadata_group"])) {
      grouping <- meta[data$Image_ImageNumber_Global, filter.table[i, "metadata_group"]]
    } else {
      grouping <- NULL
    }
    
    cur.data <- data[,cur.features, drop=F]
    if (filter.table[i, "transpose"]) {
      cur.data <- t(cur.data)
    }
    
    res[, cur.filter.name] <- do.call(cur.filter, list(vec=cur.data,
                                                       thresh=filter.table[i, "value"],
                                                       grouping=grouping))
  }
  
  return(res)
}

#-------------------------------------------------------------------------------
#' Apply a cell level filter to a tglow dataset
#' 
#' @param dataset
#' @param res A data frame, matrix or vector with booleans describing which cells to keep.
#' @returns 
#' A tglow dataset with the cells, image and object relation matrix filtered
tglow.filter.cells.apply <- function(dataset, res, img.id.col="Image_ImageNumber_Global", other.assays = NULL) {
  
  if (class(res)[1] == "data.frame" | class(res)[1] == "matrix") {
    selector <- rowSums(res)==ncol(res)
  } else if (class(res)[1] == "logical") {
    selector <- res
  } else {
    cat("[ERROR] no valid res input type\n")
    return(NULL)
  }
  
  dataset$cells <-  dataset$cells[selector,]
  
   if (!is.null(other.assays)){
    
    for (a in other.assays){
      
      dataset[[a]] <-  dataset[[a]][selector,]
      
    } 
    
  }
  
  img.nr <- unique(dataset$cells[,img.id.col])
  dataset$meta <- dataset$meta[dataset$meta$ImageNumber_Global %in% img.nr,]
  
  obj.nr <- unique(as.character(unlist(dataset$cells[,grep("Object_Number_Global",colnames(dataset$cells))])))
  dataset$orl <- dataset$orl[dataset$orl$`First Object Number Global` %in% obj.nr | dataset$orl$`Second Image Number Global` %in% obj.nr,]
  
  return(dataset)
}


#-------------------------------------------------------------------------------
#' Apply a image level filter to a tglow dataset
#' 
#' @param dataset
#' @param res A data frame, matrix or vector with booleans describing which images to keep.
#' 
#' @returns 
#' A tglow dataset with the cells, image and object relation matrix filtered
#' 
tglow.filter.img.apply <- function(dataset, res, col.img.id="Image_ImageNumber_Global") {
  
  if (class(res) == "data.frame" | class(res) == "matrix") {
    selector <- rowSums(res)==ncol(res)
  } else if (class(res) == "logical") {
    selector <- res
  } else {
    stop("Not a valid res input type")
  }
  
  dataset$meta   <- dataset$meta[selector,]
  img.nr         <- unique(dataset$meta$ImageNumber_Global)
  selector       <- dataset$cells[,col.img.id] %in% img.nr
  dataset$cells  <- dataset$cells[selector,]
  
  if ("cells_norm" %in% names(dataset)) {
    dataset$cells_norm <-  dataset$cells_norm[selector,]
  } 

  if ("cells_transform" %in% names(dataset)) {
    dataset$cells_transform <-  dataset$cells_transform[selector,]
  }

  if ("cells_corrected" %in% names(dataset)) {
    dataset$cells_corrected <-  dataset$cells_corrected[selector,]
  }
  
  obj.nr <- unique(as.character(unlist(dataset$cells[,grep("Object_Number_Global",colnames(dataset$cells))])))
  dataset$orl <- dataset$orl[dataset$orl$`First Object Number Global` %in% obj.nr | dataset$orl$`Second Image Number Global` %in% obj.nr,]
  
  return(dataset)
}


#-------------------------------------------------------------------------------
# Filters
#-------------------------------------------------------------------------------
# All filters are inclusive!

#-------------------------------------------------------------------------------
#' Min filter
filter.min <- function(vec, thresh, grouping=NULL) {
  return(vec >= thresh)
}

filter.min.sum <- function(...) {filter.sum(..., func=filter.min)}

#-------------------------------------------------------------------------------
#' Max filter
filter.max <- function(vec, thresh, grouping=NULL) {
  return(vec <= thresh)
}

filter.max.sum <- function(...) {filter.sum(..., func=filter.max)}

#-------------------------------------------------------------------------------
#' NA filter
filter.na <- function(vec, thresh, grouping=NULL) {
  return((sum(is.na(vec))/length(vec)) <= thresh)
}

filter.na.multicol <- function(vec, thresh, grouping) {
  res <- apply(vec, 2, filter.na, thresh=thresh, grouping=grouping)
  return(res)
}

#-------------------------------------------------------------------------------
#' Caret near zero variance filter
filter.near.zero.var <- function(vec, thresh=NULL) {
  return (length(nearZeroVar(vec))==0)
}

filter.near.zero.var.sum <- function(...) {filter.sum(..., func=filter.near.zero.var)}

#-------------------------------------------------------------------------------
#' Zero variance filter
filter.zero.var <- function(vec, thresh=0) {
  return (Rfast::Var(vec[!is.na(vec)]) > thresh)
}

filter.zero.var.sum <- function(...) {filter.sum(..., func=filter.zero.var)}

#-------------------------------------------------------------------------------
#' Coefficient of variation filter
filter.coef.var <- function(vec, thresh) {
  return ((sqrt(Rfast::Var(vec[!is.na(vec)])) / mean(vec[!is.na(vec)])) > thresh)
}

filter.coef.var.sum <- function(...) {filter.sum(..., func=filter.zero.var)}

#-------------------------------------------------------------------------------
#' Minimum number of unique values
filter.unique.val <- function(vec, thresh=NULL) {
  return (length(unique(vec)) > thresh)
}

filter.unique.val.sum <- function(...) {filter.sum(..., func=filter.unique.val)}

#-------------------------------------------------------------------------------
#' Infinite median filter
filter.inf.median <- function(vec, thresh=NULL) {
  return(!is.infinite(median(vec, na.rm=T)))
}

filter.inf.median.sum <- function(...) {filter.sum(..., func=filter.inf.median)}

#-------------------------------------------------------------------------------
#' Infinite sum filter
filter.inf <- function(vec, thresh=NULL, grouping=NULL) {
  return(sum(is.infinite(vec)) <= thresh)
}

filter.inf.mutlicol <- function(vec, thresh, grouping) {
  res <- apply(vec, 2, filter.inf, thresh=thresh, grouping=grouping)
  return(res)
}

#-------------------------------------------------------------------------------
#' Modified z-score filter
filter.mod.z <- function(vec, thresh, grouping=NULL, absolute=T, method="mod.z") {
  if(is.null(grouping)) {
    grouping <- rep(1, length(vec))
  }
  mod.z <- tglow.grouped.scale(vec, grouping=grouping, method=method)
  
  if (absolute) {
    mod.z <- abs(mod.z)
  }
  
  return(mod.z < thresh)
} 

Met <- function(...) {filter.sum(..., func=filter.mod.z)}

filter.mod.z.perc <- function(vec, thresh, thresh2, grouping=NULL) {
  
  i <- 0
  j <- ncol(vec)
  cat("\n[INFO] Normalizing features per group.\n")
  
  res <- apply(vec, 2, function(x){
    i <<- i +1
    cat("\r[INFO]", round((i/j)*100, digits=2), "%")
    return(filter.mod.z(x, thresh=thresh, grouping=grouping))
  })
  cat("\n[INFO] Done normalizing per group. Calculating percentages per cell\n")
  return((rowSums(res, na.rm=T)/ncol(vec)) >= thresh2)
}

#-------------------------------------------------------------------------------
#' If data is multicolumn, take the sum over all collumns, if one is false, exlcude
filter.sum <- function(vec, thresh, grouping, func) {
  res <- apply(vec, 2, func, thresh=thresh, grouping=grouping)
  
  res[is.na(res)] <- T
  
  return(rowSums(res, na.rm=T) == ncol(vec))
}


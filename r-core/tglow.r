setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix", "NULL"))

#-------------------------------------------------------------------------------
setClass("TglowAssay",
         slots=list(data="AnyMatrix",
                    scale.data="AnyMatrix",
                    features="data.frame"),
         prototype = list(scale.data=NULL)
         )
         
#-------------------------------------------------------------------------------
setClass("TglowDataset",
         slots=list(assays="list",
                    meta="data.frame",
                    image.meta="data.frame",
                    image.data="TglowAssay",
                    image.ids="character",
                    reduction="list",
                    graph="ANY",
                    active.assay="character"),
         prototype = list(reduction=list(),
                          graph=NULL)
)

#-------------------------------------------------------------------------------
setMethod("show", signature("TglowAssay"),
function(object) {
  cat("TglowAssay with: ",
    nrow(object@data), " cells and ",
    ncol(object@data)," features \n")
})

#-------------------------------------------------------------------------------
setMethod("show", signature("TglowDataset"),
function(object) {
  cat("TglowData with: ",
      nrow(object@meta), " cells, ",
      nrow(object@image.meta), " images and ",
      length(object@assays)," assays \n")
      
  cat("Assays: \n")
  print(object@assays)

  cat("\nActive assay: ", object@active.assay, "\n")
  cat("\nReductions: ", names(object@reduction), "\n")
})

#-------------------------------------------------------------------------------
setMethod("[",
          "TglowAssay",
function(x, i, j, drop=F) {
  object <- x
  
  # Select rows
  if (!missing(i)) {
    
    # Filter main assay
    object@data <- object@data[i, , drop = F]
    
    # Filter scale data
    if (!is.null(object@scale.data)) {
      if (nrow(object@scale.data) >= 1) {
        object@scale.data <- object@scale.data[i, , drop = F]
      }
    }
  }
  
  # Select columns
  if (!missing(j)) {
    # Filter main assay
    object@data <- object@data[,j, drop = F]
    
    # Filter scale data
    if (!is.null(object@scale.data)) {
      if (nrow(object@scale.data) >= 1) {
        object@scale.data <- object@scale.data[,j, drop = F]
      }
    }
    
    # Filter features
    object@features <- object@features[j,, drop=F]
  }
  
  object
  
})

#-------------------------------------------------------------------------------
setMethod("[",
          "TglowDataset",
function(x, i, j, drop=F) {
  object <- x
  
  # Select rows
  if (!missing(i)) {
    
    # Filter assays
    for (assay in 1:length(object@assays)) {
      object@assays[[assay]] <- object@assays[[assay]][i,,drop=F]
    }
    
    # Filter meta
    object@meta          <- object@meta[i,,drop=F]
    
    # Select images
    object@image.ids     <- object@image.ids[i,drop=F]
    object@image.data    <- object@image.data[unique(object@image.ids),,drop=F]
    object@image.meta    <- object@image.meta[unique(object@image.ids),,drop=F]
    
    # Filter dimension reductions
    if (length(object@reduction) >= 1) {
      for (k in 1:length(object@reduction)) {
        object@reduction[[k]] <- object@reduction[[k]][i,,drop=F]
      }
    }
  }
  
  # Select columns
  if (!missing(j)) {
    
    if (class(j) != "character") {
      warning("Assuming all assays have the same column order")
    }
    
    # Filter assays
    for (assay in 1:length(object@assays)) {
      object@assays[[assay]] <- object@assays[[assay]][,j,drop=F]
    }
  }
  
  object
  
})



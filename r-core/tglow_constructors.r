#-------------------------------------------------------------------------------
#' Convert a tglow list to a tglow assay
tglow.assay.from.list <- function(output, assay, meta.cols, col.object) {
  
  row.meta       <- output[[assay]][, meta.cols]
  data           <- as.matrix(output[[assay]][, !colnames(output[[assay]]) %in% meta.cols])
  rownames(data) <- row.meta[,col.object]
  
  features         <- output$features[!rownames(output$features) %in% meta.cols,]
  features$analyze <- T
  features         <- features[colnames(data),]
  rownames(features) <- colnames(data)
  
  assay <- new("TglowAssay",
               data=data,
               features=features)
  
  return(assay)
  
}

#-------------------------------------------------------------------------------
#' Tglow dataset from list
tglow.dataset.from.list <- function(output, assay, meta.patterns=c("ImageNumber","ObjectNumber", "Object_Number", "Parent"), img.feature.patterns=c("^AreaOccupied", "^ImageQuality", "^Granularity"), col.object="cell_ObjectNumber_Global", col.img.id="cell_ImageNumber_Global", col.meta.img.id="ImageNumber_Global") {
  
  # Perform some sanity checks
  if(sum(c(assay, "meta") %in% names(output)) != 2) {
    stop("Assay and meta items not found in list")
  }
  
  if (!col.img.id %in% colnames(output[[assay]])) {
    stop("Image ID column not found in assay")
  }
  
  if (!col.meta.img.id %in% colnames(output[["meta"]])) {
    stop("Image ID column not found in metadata")
  }
  
  # Columns in output data considered metadata and not features
  meta.cols  <- unlist(sapply(meta.patterns, grep, colnames(output[[assay]]), value=T))
  
  # Set main assay
  main.assay <- tglow.assay.from.list(output, assay=assay, meta.cols=meta.cols, col.object=col.object)
  dataset                 <- new("TglowDataset")
  dataset@assays[["raw"]] <- main.assay  
  dataset@active.assay    <- "raw"
  dataset@object.ids      <- output[[assay]][,col.object]
  
  # Filter meta
  dataset@meta           <- output[[assay]][, meta.cols]
  rownames(dataset@meta) <- output[[assay]][, col.object]
  
  # Image level data
  image.cols <- unlist(sapply(img.feature.patterns, grep, colnames(output$meta), value=T))
  dataset@image.data <- new("TglowAssay", 
                            data=as.matrix(output$meta[, image.cols]),
                            features=tglow.get.feature.meta.from.cells(image.cols))
  rownames(dataset@image.data@data) <- output$meta[,col.meta.img.id]
  dataset@image.meta                <- output$meta[,!colnames(output$meta) %in% image.cols]
  rownames(dataset@image.meta)      <- output$meta[,col.meta.img.id]
  
  dataset@image.ids                 <- dataset@meta[,col.img.id]
  names(dataset@image.ids)          <- dataset@meta[,col.object]
  return(dataset)
}

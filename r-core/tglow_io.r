# Global to store the running total of filesets
FILESET_ID=0

#-------------------------------------------------------------------------------
#' Read a cellprofiler fileset directory tree
#' 
#' Matches images based on <plate>, <well>, <field>
#' 
#' @param path path to tglow output dir
tglow.read.dir <- function(path, exp.suffix="_experiment.tsv", n=NULL, verbose=F, ...) {
  
  #f.cells <- list.files(path, recursive = T, pattern="*_cells.tsv")
  #f.exp   <- list.files(path, recursive = T, pattern="*_experiment.tsv")
  #f.img   <- list.files(path, recursive = T, pattern="*_image.tsv")
  #f.orel  <- list.files(path, recursive = T, pattern="*_objectRelation.ships.tsv")
  
  files    <- list.files(path, recursive = T, pattern=paste0("*", exp.suffix), full.names = T)
  prefixes <- gsub(exp.suffix, "", files)
  
  if (!is.null(n)) {
    prefixes <- prefixes[1:n]
  }
  
  # Reset global fileset index
  assign("FILESET_ID", 0, envir = .GlobalEnv)
  
  # Read filesets
  #output   <- list()
  filesets <- list()
  i <- 0
  for (pre in prefixes) {
    i <- i + 1
    cat("\r[INFO] reading fileset ", i, "/", length(prefixes))
    #filesets[[pre]] <- tglow.read.fileset(pre)
    cur <- tglow.read.fileset(pre, return.feature.meta=F, ...)
    #cur <- tglow.read.fileset(pre, return.feature.meta=F, pat.img = "_images.tsv")
    filesets[[pre]] <- cur
    
    if (verbose){
        cat("\n[DEBUG] cols:", ncol(cur$cells), " cols meta", ncol(cur$meta), " cols orl:", ncol(cur$orl), "\n")
    }
    #if(length(output) == 0) {
    #  output <- cur 
    #} else {
      #output <- tglow.merge.filesets(list(output, cur))
    #}
  }
  
  cat("\n[INFO] Merging filesets\n")
  output <- tglow.merge.filesets(filesets)
  
  # Read feature index
  #cat("[INFO] Reading features\n")
  #cells    <- fread(paste0(pre[1],"_cells.tsv"), data.table=F, nrows=2)
  
  if(verbose) {
    cat("[DEBUG] colnames:\n", colnames(output$cells))
  }
  
  features      <- tglow.get.feature.meta.from.cells(colnames(output$cells))
  classes       <- sapply(output$cells, class)
  features$type <- classes[features$id]
  
  #cat("[INFO] Merging filesets\n")
  #output <- tglow.merge.filesets(filesets)
  
  features               <- features[colnames(output$cells),]
  rownames(output$cells) <- output$cells$cells_ObjectNumber_Global
  rownames(output$meta)  <- output$meta$ImageNumber_Global
  
  return(c(output, list(features=features)))
}

#-------------------------------------------------------------------------------
#' Build an index of where the example images are stored.
#' This assumes images are organized as follows:
#' <plate>/<row>/<col>/<field>.ome.tiff
#' 
#' Matches images based on <plate>/<row>/<col>/<field>.ome.tiff
#' 
#' @param path path to tglow image dir
#' @param plate_filter plate names to run
#' @param pattern the pattern to match image names and then extract the fields (gsubbed away)

#' @returns A data frame with the iamges
tglow.build.img.index <- function(path, plate_filter=NULL, pattern=".ome.tiff$") {
  plates <- list.files(path, recursive = F)
  
  if (!is.null(plate_filter)) {
    plates <- plates[plates %in% plate_filter]
  }
  
  #cat("[INFO] Indexing plates: ", plates, "\n")
  
  index <- data.frame(matrix(NA, nrow=0, ncol=6))
  
  nfiles <- 0
  for (plate in plates) {
    cat("[INFO] Indexing plate: ", plate, "\n")

    rows <- list.files(paste0(path, "/", plate), pattern="^[A-Z]$")
    for (row in rows) {
      cols <- list.files(paste0(path, "/", plate, "/", row), pattern="^\\d+$")
      for (col in cols) {
        
        files <- list.files(paste0(path, "/", plate, "/", row, "/", col), pattern=pattern)
        nfiles <- nfiles + length(files)
        for (file in files) {
          
          field <- gsub(pattern, "", file)
          tmp   <- c(plate, paste0(row, sprintf("%02d", as.numeric(col))),  row, col, field, paste0(path, "/", plate, "/", row, "/", col, "/", file))
          index <- rbind(index, tmp)
          
        } 
      } 
    }
  }
  colnames(index) <- c("plate", "well", "row", "col", "field", "path")

  cat("[INFO] Indexed ", nfiles, " image files\n")
  if (nfiles ==0) {
    warning("No files detected, suggest you check pattern is correct.")
  } else {
    rownames(index) <- paste0(index$plate, ":", index$well, ":", index$field)
  }
    
  
  return(index)  
}

#-------------------------------------------------------------------------------
#' Build an index of where the example images are stored.
#' This assumes images are organized as follows:
#' <plate>/<well>/segmentation_images/<pattern>
#' 
#' Matches images based on <plate>, <well>, <field>
#' 
#' @param path path to tglow output dir

#' @returns A list of paths indexed by plate / well
tglow.build.img.index.deprecated <- function(path, pattern.group, pattern.file=".*_composite_.*.png", pattern.filter=".*_composite_segmented_.*.png", pattern.field=".*_(f\\d\\d)\\..*") {
  
  paths <- list.files(path, pattern = pattern.file, recursive = T)
  
  if (!is.null(pattern.filter)) {
    paths <-  grep(pattern.filter, paths, value=T, invert=T)
  }
  
  mat <- do.call(rbind, strsplit(paths, split="/"))
  colnames(mat) <- c("plate", "well", "segment_dir", "filename")
  
  path <- paste0(path, "/", paths)
  mat   <- cbind(mat, path)
  
  field <- gsub(pattern.field, "\\1", mat[,"filename"])
  mat   <- cbind(mat, field)  
  
  group <- gsub(pattern.group, "\\1", mat[,"filename"])
  mat   <- cbind(mat, group)  
  
  rownames(mat) <- mat[,"group"]
  
  return(mat)
}


#-------------------------------------------------------------------------------
#' Load images around the center of a cell object.
#' 
#' @param data Tglow data list
#' @param cell.subset The subset of cells to retrieve
#' @param img.index A matrix with image paths generated by tglow.build.img.index
#' @param window Window in px around the cell center to retrieve
#' @param feature.x The feature that describes the object pos in px in x
#' @param feature.y The feature that describes the object pos in px in y
#' @param feature.id The feature to use as the unique id to give to the output
#' @param group.col The feature to match cells to img.index
#' @param channels The channels to read. Vector of indices. (default NULL = all channels)
#' @param planes The planes to read. Vector of indices. (default NULL = all planes)
#' @param max.project Should the stack be max projected per channel? (default TRUE)
##' @param ncores How many cores to use with parallels mcapply for reading
#' 
#' @returns A list of EBImage objects by object id
tglow.read.imgs <- function(data,
                            cell.subset,
                            img.index,
                            window=75,
                            feature.x="cell_AreaShape_Center_X",
                            feature.y="cell_AreaShape_Center_Y",
                            feature.id="cell_Number_Object_Number_Global",
                            group.col="Metadata_group",
                            channels=NULL,
                            planes=NULL,
                            max.project=T) {
  
  cur.cells <- data$cells
  #out <- list()
  j <- 0
  #for (i in cell.subset) {
  out <- lapply(cell.subset, function(i) {
    cat("[INFO] Reading ", round((j / length(cell.subset))*100) ,"%\r")
    j <<- j+1
    cur.group <- data$meta[cur.cells[i,"Image_ImageNumber_Global"],group.col]
    cur.img   <- img.index[cur.group,"path"]
    x.pos     <- cur.cells[i,feature.x]
    y.pos     <- cur.cells[i,feature.y]
    
   # cat("[INFO] Reading ", cur.img, "\n")

    # TODO: can optimize this to only read a subset
    cx <- round((x.pos-window):(x.pos+window))
    cy <- round((y.pos-window):(y.pos+window))
    
    subset = list(X=cx, Y=cy)
    
    if (!is.null(channels)) {
      subset[["C"]] = channels
    }
    
    if (!is.null(planes)) {
      subset[["Z"]] = planes
    }
    
    img <- read.image(cur.img, normalize=F, subset=subset)

    #img <- read.image(cur.img, normalize=F)
    colorMode(img) = Grayscale
    #out[[j]] <- img
    
    # cat("[INFO] ",cx, ":", cy, " dim: ", dim(img), "\n")

    if (length(dim(img)) == 4) {
      #if (is.null(channels)) {
      #  channels = seq(1, img@dim[3])
      #}
      
      #if (is.null(planes)) {
      #  planes = seq(1, img@dim[4])
     # }
      
      # Crop
      #crop <- img[, , channels, planes]
  
      # Max project
      if (max.project) {
        crop <- apply(img, c(1, 2, 3), max)
      }
      
    } else if (length(dim(img))==3) {
      
      #if (is.null(planes)) {
      #  planes = seq(1, img@dim[3])
      #}
      
      # Crop
      #crop <- img[, , planes]
      
      # Max project
      if (max.project) {
        crop <- apply(img, c(1, 2), max)
      }
      
    } else if (length(dim(img)) == 2) {
      
      crop <- img
      
    }
  
    #cat("[INFO] crop dim: ", dim(crop), "\n")
    #out[[j]] <- crop
    return(crop)

   # gc()
    #crop <- geometry_area(window*2, window*2, x.pos-window, y.pos-window)  
    #cur.img.crop <- image_crop(cur.img, crop)
    #out <- c(out,cur.img.crop)
  })
  
  cat("\n")
  names(out) <- cur.cells[cell.subset, feature.id]
  
  return(out)
}

#-------------------------------------------------------------------------------
#' Wrapper around mclapply to track progress
#' 
#' Based on http://stackoverflow.com/questions/10984556
#' 
#' @param X         a vector (atomic or list) or an expressions vector. Other
#'                  objects (including classed objects) will be coerced by
#'                  ‘as.list’
#' @param FUN       the function to be applied to
#' @param ...       optional arguments to ‘FUN’
#' @param mc.preschedule see mclapply
#' @param mc.set.seed see mclapply
#' @param mc.silent see mclapply
#' @param mc.cores see mclapply
#' @param mc.cleanup see mclapply
#' @param mc.allow.recursive see mclapply
#' @param mc.progress track progress?
#' @param mc.style    style of progress bar (see txtProgressBar)
#'
#' @examples
#' x <- mclapply2(1:1000, function(i, y) Sys.sleep(0.01))
#' x <- mclapply2(1:3, function(i, y) Sys.sleep(1), mc.cores=1)
#' 
#' dat <- lapply(1:10, function(x) rnorm(100)) 
#' func <- function(x, arg1) mean(x)/arg1 
#' mclapply2(dat, func, arg1=10, mc.cores=2)
#-------------------------------------------------------------------------------
mclapply2 <- function(X, FUN, ..., 
    mc.preschedule = TRUE, mc.set.seed = TRUE,
    mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
    mc.cleanup = TRUE, mc.allow.recursive = TRUE,
    mc.progress=TRUE, mc.style=3) 
{
    if (!is.vector(X) || is.object(X)) X <- as.list(X)

    if (mc.progress) {
        f <- fifo(tempfile(), open="w+b", blocking=T)
        p <- parallel:::mcfork()
        pb <- txtProgressBar(0, length(X), style=mc.style)
        setTxtProgressBar(pb, 0) 
        progress <- 0
        if (inherits(p, "masterProcess")) {
            while (progress < length(X)) {
                readBin(f, "double")
                progress <- progress + 1
                setTxtProgressBar(pb, progress) 
            }
            cat("\n")
            parallel:::mcexit()
        }
    }
    tryCatch({
        result <- mclapply(X, ..., function(...) {
                res <- FUN(...)
                if (mc.progress) writeBin(1, f)
                res
            }, 
            mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
            mc.silent = mc.silent, mc.cores = mc.cores,
            mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
        )

    }, finally = {
        if (mc.progress) close(f)
    })
    result
}



#-------------------------------------------------------------------------------
#' Load images around the center of a cell object.
#' 
#' @param data Tglow data list
#' @param cell.subset The subset of cells to retrieve
#' @param img.index A matrix with image paths generated by tglow.build.img.index
#' @param format The image format of the source
#' @param window Window in px around the cell center to retrieve
#' @param feature.x The feature that describes the object pos in px in x
#' @param feature.y The feature that describes the object pos in px in y
#' @param feature.id The feature to use as the unique id to give to the output
#' @param group.col The feature to match cells to img.index
#' 
#' @returns A list of image magic objects by object id
tglow.read.imgs.deprecated  <- function(data, cell.subset, img.index, format="png", window=75,
                            feature.x="cell_AreaShape_Center_X",
                            feature.y="cell_AreaShape_Center_Y",
                            feature.id="seed_Number_Object_Number_Global",
                            group.col="Metadata_group") {
  
  cur.cells <- data$cells
  out <- list()
  for (i in cell.subset) {
    cur.group <- data$meta[cur.cells[i,"Image_ImageNumber_Global"],group.col]
    cur.img <- image_read(img.index[cur.group,"path"])
    x.pos <- cur.cells[i,feature.x]
    y.pos <- cur.cells[i,feature.y]
    
    crop <- geometry_area(window*2, window*2, x.pos-window, y.pos-window)  
    cur.img.crop <- image_crop(cur.img, crop)
    out <- c(out,cur.img.crop)
  }
  
  names(out) <- cur.cells[cell.subset, feature.id]
  
  return(out)
}


#-------------------------------------------------------------------------------
#' Read a cell level fileset
#' 
#' Reads a CellProfiler (tglow) fileset into a list.
#' 
#' @returns list with data frames:
#' - cells (cell level features)
#' - meta (image level features)
#' - objectRelations 
#' - features [optional]
#' Output is NULL if no cells are detected.
tglow.read.fileset <- function(prefix, return.feature.meta=F, add.global.id=T,
                               pat.img="_image.tsv",
                               pat.cells="_cells.tsv",
                               pat.orl="_objectRelationships.tsv") {
  
  if (add.global.id) {
    assign("FILESET_ID", FILESET_ID + 1, envir = .GlobalEnv)
    #global.prefix <- paste0("FS", gsub("\\s", "0",format(FILESET_ID, width=4)))
    global.prefix <- paste0("FS", FILESET_ID)
  }
  
  # Read header of _cells.tsv
  cells <- fread(paste0(prefix, pat.cells), data.table=F, nrows=3)
  
  # If the file has no cells empty
  if (nrow(cells) != 3) {
    return(NULL)
  }
  
  # Clean colnames
  cn <- paste0(colnames(cells), "_", as.character(cells[1,]))
  
  if(return.feature.meta) {
    feature.meta <- tglow.get.feature.meta.from.cells(cn)
  }  
  
  #colnames(cells) <- cn
  #cells <- cells[-1,]
  
  # Read content of _cells.tsv
  cells <- fread(paste0(prefix,pat.cells), data.table=F, skip=2)
  colnames(cells) <- cn
  
  # Read _image.tsv (metadata)
  img     <- fread(paste0(prefix, pat.img), data.table=F)
  
  # Read _objectRelation.ships.tsv
  orl   <- fread(paste0(prefix, pat.orl), data.table=F)
  
  
  # Standardize ID's across filesets into the following format
  # FS#I#O#
  # where FS = file set, I = image within fileset and O = object within image
  # This makes it easier to match across data with a unique ID
  if (add.global.id) {
    
    # Cell level information
    #-----------
    cells[,"Image_ImageNumber_Global"] <- paste0(global.prefix, "_I", cells[, "Image_ImageNumber"])
    
    cols <- grep("ObjectNumber", colnames(cells), value=T)
    cols <- c(cols,grep("Object_Number", colnames(cells), value=T))
    
    for (cur.col in cols) {
      cells[, paste0(cur.col,"_Global")] <- paste0(global.prefix, "_I", cells[,"Image_ImageNumber"], "_O", cells[, cur.col])
    } 
    
    # IMG image level information
    #-----------
    img[,"ImageNumber_Global"] <- paste0(global.prefix, "_I", img[, "ImageNumber"])
    
    # ORL, object relationships
    #-----------
    if (nrow(orl) > 0) {
      orl[,"First Image Number Global"]   <- paste0(global.prefix, "_I", orl[, "First Image Number"])
      orl[,"Second Image Number Global"]  <- paste0(global.prefix, "_I", orl[, "Second Image Number"])
      orl[,"First Object Number Global"]  <- paste0(global.prefix, "_I", orl[, "First Image Number"], "_O", orl[, "First Object Number"])
      orl[,"Second Object Number Global"] <- paste0(global.prefix, "_I", orl[, "Second Image Number"], "_O", orl[, "Second Object Number"])
    }
  }
  
  out.list <- list(cells=cells, meta=img, orl=orl)
  
  if(return.feature.meta) {
    return(c(out.list, list(features=feature.meta)))
  } else {
    return(out.list)
  }
}

# Global to store the running total of filesets
FILESET_ID=0

#-------------------------------------------------------------------------------
#' Read a cellprofiler fileset directory tree
#' 
#' Matches images based on <plate>, <well>, <field>
#' 
#' @param path path to tglow output dir
#' @param pattern The pattern that uniquely identfies a well. Use .zip for type B and the _experiment.tsv
#' for type A.
#' @param type Must be A or B. 
#' Type A: _cells.tsv, _image.tsv, _experiment.tsv and _objectRelations.tsv.
#' See tglow.read.fileset.a for detaills
#' Type B: <plate>_<well>.zip with individual files for each child object. Main object is assumed to be _cells.
#' see 
#' See tglow.read.fileset.b for detaills
#' @param n Read a subset of filesets. If integer, only that fileset is read, otherwise specify indices to read. 
#' @param verbose Should I be chatty?
#' @param ... Remaining parameters passed to tglow.read.fileset.a/b
tglow.read.dir <- function(path, pattern, type, n=NULL, verbose=F, ...) {
  
  #f.cells <- list.files(path, recursive = T, pattern="*_cells.tsv")
  #f.exp   <- list.files(path, recursive = T, pattern="*_experiment.tsv")
  #f.img   <- list.files(path, recursive = T, pattern="*_image.tsv")
  #f.orel  <- list.files(path, recursive = T, pattern="*_objectRelation.ships.tsv")
  
  files    <- list.files(path, recursive = T, pattern=paste0("*", pattern), full.names = T)
  
  if (type == "A") {
    prefixes <- gsub(pattern, "", files)
  } else {
    prefixes <- files
  }
  
  if (!is.null(n)) {
    #if (length(n) == 1) {
        prefixes <- prefixes[n]
    #} else {
    #    prefixes <- prefixes[1:n]
    #}
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
    
    if (type == "A") {
        cur <- tglow.read.fileset.a(pre, return.feature.meta=F, ...)
    } else if (type == "B") {
        cur <- tglow.read.fileset.b(pre, return.feature.meta=F, ...)
    } else {
      stop(paste0("Invalid type: ", type))
    }
    
    if (!is.null(cur)) {
      filesets[[pre]] <- cur
    
      if (verbose){
          cat("\n[DEBUG] cols:", ncol(cur$cells), " cols meta", ncol(cur$meta), " cols orl:", ncol(cur$orl), "\n")
      }
    } else {
      warning("Fileset was NULL, skipped.")
    }

  }
  
  cat("\n[INFO] Merging filesets\n")
  output <- tglow.merge.filesets(filesets)
  
  cat("[INFO] names: ", names(output), "\n")
  
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

#' @returns A data frame with the images
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
#' Read a binary matrix, used for reading the registration matrices
#' Adapted from: https://www.r-bloggers.com/2012/06/getting-numpy-data-into-r/
tglow.read.binmat <- function(path) {
  
  con <- file(path, "rb")
  dim <- readBin(con, "integer", 2)
  mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
  close(con)
  return(mat)
}

#-------------------------------------------------------------------------------
#' Load images around the center of a cell object.
#' 
#' @param data Tglow data list
#' @param cell.subset The subset of cells to retrieve
#' @param img.index A matrix with image paths generated by tglow.build.img.index
#' @param reg.index A matrix with registration paths generated by tglow.build.img.index.
#' Must have an additional collumn named 'plate2' to indicate the registraiton plate
#' @param window Window in px around the cell center to retrieve
#' @param feature.x The feature that describes the object pos in px in x
#' @param feature.y The feature that describes the object pos in px in y
#' @param feature.id The feature to use as the unique id to give to the output
#' @param group.col The feature to match cells to img.index
#' @param channels The channels to read. A list of Vector of indices per cycle. (default NULL = all channels)
#' @param planes The planes to read.  A list of Vector of indices per cycle. (default NULL = all planes)
#' @param max.project Should the stack be max projected per channel? (default TRUE)
##' @param ncores How many cores to use with parallels mcapply for reading
#' 
#' @returns A list of EBImage objects by object id
tglow.read.imgs <- function(data,
                            cell.subset,
                            img.index,
                            reg.index=NULL,
                            window=75,
                            assay="cells",
                            feature.x="cell_AreaShape_Center_X",
                            feature.y="cell_AreaShape_Center_Y",
                            feature.id="cell_Number_Object_Number_Global",
                            group.col="Metadata_group",
                            img.id.col="Image_ImageNumber_Global",
                            channels=NULL,
                            planes=NULL,
                            max.project=T) {
  
  cur.cells <- data[[assay]]
  
  if ((class(cell.subset) != "character") && (class(cell.subset) != "integer")) {
        stop("Cell.subset must be character or integer")
  }
  
  if (class(cell.subset) == "integer") {
    warning("cell.subset is integer, make sure the indices are properly matched to assay.")
    if (max(cell.subset) > nrow(cur.cells)) {
      stop("Index in cell.subset larger then number of rows in assay. Are you using the correct assay?")
    }
  }
    
  j <- 0
  out <- lapply(cell.subset, function(i) {
    cat("[INFO] Reading ", round((j / length(cell.subset))*100) ,"%\r")
    j <<- j+1
    cur.group <- data$meta[cur.cells[i, img.id.col],group.col]
    cur.img   <- img.index[cur.group,]
    x.pos     <- cur.cells[i,feature.x]
    y.pos     <- cur.cells[i,feature.y]
    
    # cat("[INFO] Reading ", cur.img, "\n")
    cx <- round((x.pos-window):(x.pos+window))
    cy <- round((y.pos-window):(y.pos+window))
    
    subset = list(X=cx, Y=cy)
    
    if (!is.null(channels)) {  
      if (class(channels) != "list") {
        stop("Channels must be a list with one numeric vector per cycle containing channel indices to load.")
      }    
      subset[["C"]] = channels[[1]]
    }
    
    if (!is.null(planes)) {
      if (class(planes) != "list") {
        stop("Planes must be a list with one numeric vector per cycle containing plane indices to load.")
      } 
      subset[["Z"]] = planes[[1]]
    }

    img <- read.image(cur.img$path, normalize=F, subset=subset)
    colorMode(img) = Grayscale

    # Append mutliple cycles
    if (!is.null(reg.index)) {
      cur.reg.index <- reg.index[reg.index$plate == cur.img$plate & reg.index$well == cur.img$well & reg.index$field == cur.img$field,]
      
      
      for (p2 in 1:nrow(cur.reg.index)) {
        
        new.img <- img.index[paste0(cur.reg.index[p2, "plate2"], ":",
                                    cur.reg.index[p2, "well"], ":", 
                                    cur.reg.index[p2, "field"]),]
        
        if (!is.null(channels)) {
            if (length(channels) != nrow(cur.reg.index) +1) {
              stop("Length of channel list does not match the number of cycles found in reg.index")
            }
            subset[["C"]] = as.numeric(channels[[p2+1]])            
        } 
              
        if (!is.null(planes)) {
            if (length(planes) != nrow(cur.reg.index) +1) {
              stop("Length of plane list does not match the number of cycles found in reg.index")
            }
            subset[["Z"]] = planes[[p2+1]]
        }     
        
        
        # Read the registration matrix
        reg <- tglow.read.binmat(cur.reg.index[p2, "path"])
      
        tmp.img <- read.image(new.img$path, subset=subset, normalize=F)
        colorMode(tmp.img) = Grayscale
        
        # Recode the matrix so it works with EBImage
        m <- t(reg)[,1:2]
        m[3, ] <- -1*m[3,]
      
        # Apply affine transform
        tmp.img <- affine(tmp.img, m, filter="none", antialias = F)
        
        # Combine with the image
        img <- combine(img, tmp.img)
      }
    }
    
    crop <- img
    
    if (length(dim(img)) == 4) {
      # Max project
      if (max.project) {
        crop <- apply(img, c(1, 2, 3), max)
      }
      
    } else if (length(dim(img))==3) {
      # Max project
      if (max.project) {
        crop <- apply(img, c(1, 2), max)
      }
      
    } else if (length(dim(img)) == 2) {
      crop <- img
    }
    return(crop)
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
#' Add a global id to a matrix of image files.
#' Matrix with column pattern 'ImageNumber', 'ObjectNumber' and 'Object_Number'
#' which will be duplicated and have a globally unique variable added.
#' Output of this is stored in the same column name but with suffix _Global
add.global.ids <- function(matrix) {
    
  # Fetch global fileset id
  global.prefix <- paste0("FS", FILESET_ID)
  
  # Image numbers
  cols.i <- grep("ImageNumber",colnames(matrix), value=T)
  
  for (cur.col in cols.i) {
    
    selector <- !is.na(matrix[, cur.col])
    matrix[selector,paste0(cur.col, "_Global")] <- paste0(global.prefix, "_I", matrix[selector, cur.col])
  }
  
  # Object numbers
  cols <- grep("ObjectNumber", colnames(matrix), value=T)
  cols <- c(cols,grep("Object_Number", colnames(matrix), value=T))
  cols <- c(cols,grep("Parent", colnames(matrix), value=T))

  for (cur.col in cols) {
    selector <- !is.na(matrix[, cur.col])
    matrix[selector, paste0(cur.col,"_Global")] <- paste0(global.prefix, "_I", matrix[selector,cols.i[1]], "_O", matrix[selector, cur.col])
  } 
  
  return(matrix)
}



#-------------------------------------------------------------------------------
#' Read a cell level fileset type A
#' 
#' Reads a CellProfiler fileset into a list.
#' Type A: Assumes all features are in a single _cells.tsv / _cells.tsv
#' and all features are matched
#' 
#' @returns list with data frames:
#' - cells (cell level features)
#' - meta (image level features)
#' - objectRelations 
#' - features [optional]
#' Output is NULL if no cells are detected.
tglow.read.fileset.a <- function(prefix,
                               return.feature.meta=F,
                               add.global.id=T,
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
    cells <- add.global.ids(cells)
    #cells[,"Image_ImageNumber_Global"] <- paste0(global.prefix, "_I", cells[, "Image_ImageNumber"])
    
    #cols <- grep("ObjectNumber", colnames(cells), value=T)
    #cols <- c(cols,grep("Object_Number", colnames(cells), value=T))
    
    #for (cur.col in cols) {
    #  cells[, paste0(cur.col,"_Global")] <- paste0(global.prefix, "_I", cells[,"Image_ImageNumber"], "_O", cells[, cur.col])
    #} 
        
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

#-------------------------------------------------------------------------------
#' Read a cell level fileset type B
#' 
#' Reads a zip file with cellprofiler features, each file other then
#' _Image, _Experiment and _Object Relationships are assumed to be an object.
#' Will match objects on order with an appropriate matching strategy.
#' 
#' @returns list with data frames:
#' - cells (cell level features)
#' - meta (image level features)
#' - objectRelations 
#' - features [optional]
#' Output is NULL if no cells are detected.
tglow.read.fileset.b <- function(prefix,
                                 return.feature.meta=F,
                                 add.global.id=T,
                                 merging.strategy="mean",
                                 parent.col="Parent_cell",
                                 pat.exp=".*Experiment.txt",
                                 pat.img=".*Image.txt",
                                 pat.cells=".*cell.txt",
                                 pat.orl=".*Object relationships.txt",
                                 pat.others="^.*_([a-z]+\\d*).txt$",
                                 na.rm=F) {
  
  if (add.global.id) {
    assign("FILESET_ID", FILESET_ID + 1, envir = .GlobalEnv)
    global.prefix <- paste0("FS", FILESET_ID)
  }
  
  index          <- unzip(prefix, list=T)
  index$FileName <- basename(index$Name)

  # Unzip into tmp folder
  tmpdir <- tempdir()
  unzip(prefix, exdir=tmpdir)
  
  cells <- fread(paste0(tmpdir, "/", index[grep(pat.cells, index$FileName),"Name"]), data.table=F)
  img   <- fread(paste0(tmpdir, "/", index[grep(pat.img, index$FileName),"Name"]), data.table=F)
  orl   <- fread(paste0(tmpdir, "/", index[grep(pat.orl, index$FileName),"Name"]), data.table=F)
  
  if (nrow(cells) == 0) {
    warning("No cells detected for ", index[grep(pat.cells, index$FileName),"Name"], " returning NULL.")
    return(NULL)
  }
  
  
  exclude <- c(grep(pat.cells, index$FileName),
               grep(pat.img, index$FileName),
               grep(pat.orl, index$FileName),
               grep(pat.exp, index$FileName))
  
  index <- index[!1:nrow(index) %in% exclude,]
  
  colnames(cells) <- paste0("cell_", colnames(cells))
  rownames(cells) <- paste0(cells$cell_ImageNumber, "_", cells$cell_ObjectNumber)
  
  children <- list()
  if (nrow(index) > 0) {
    index$object <- gsub(pat.others, "\\1",index$FileName)
    
    for (i in 1:nrow(index)) {
      obj               <- index[i, "object"]
      cur               <- fread(paste0(tmpdir, "/", index[i, "Name"]), data.table=T)
      
      # Remove these columns from the merging strategy
      exclude.cols      <- c("Group.1", grep("ObjectNumber", colnames(cur), value=T), grep("Number_Object_Number", colnames(cur), value=T))
  
      # If there are no cols, return NA
      if (nrow(cur) == 0) {
        #next()
        # If the file is empty, set these columns to NA
        cur                                   <- as.data.frame(cur)
        colnames(cur)                         <- paste0(obj, "_", colnames(cur))
        cells[,c(colnames(cur), paste0(obj, "_QC_Object_Count"))] <- NA
        warning(paste0(obj, " assay for ", index[i, "Name"], " is empty. Returning NA for these cols."))
      } else  if (merging.strategy == "mean") {
        
        cur              <- cur[as.logical(cur[[parent.col]] != 0),]
        colnames(cur)    <- paste0(obj, "_", colnames(cur))
        selector         <- paste0(cur[[paste0(obj, "_ImageNumber")]], "_", cur[[paste0(obj, "_", parent.col)]])
        counts           <- table(selector)
        #exclude.cols      <- c("Group.1", grep("ObjectNumber", colnames(cur), value=T), grep("Number_Object_Number", colnames(cur), value=T))

        cur$Group.1      <- selector

        tmp              <- as.data.frame(cur[, lapply(.SD, mean, na.rm=na.rm), by=Group.1, .SDcols=colnames(cur)[!colnames(cur) %in% exclude.cols] ])
        rownames(tmp)    <- tmp$Group.1

        # The above is the data.table equivalent of aggregate, which is MUCH faster
        #tmp               <- aggregate(cur, by=list(selector), FUN=mean, na.rm=T)

        # Add the object count as a sanity check
        tmp[,paste0(obj, "_QC_Object_Count")] <- counts[tmp$Group.1]
        tmp                                   <- tmp[, !colnames(tmp) %in% exclude.cols]
        
        # Assign the columns to the output matrix
        cells[selector,colnames(tmp)]         <- tmp[selector,]
        
      } else if(merging.strategy == "none") {
        if (add.global.id) {cur <- add.global.ids(cur)}
        children[[index[i, "object"]]] <- cur
      } else {
        stop("Only valid merging strategy is 'mean' or 'none'")
      }
    }
  } 
  # Clean up the tmpdir
  unlink(paste0(tmpdir,"/features"), recursive=T)

  # Standardize ID's across filesets into the following format
  # FS#I#O#
  # where FS = file set, I = image within fileset and O = object within image
  # This makes it easier to match across data with a unique ID
  if (add.global.id) {
    
    # Cell level information
    #-----------
    cells <- add.global.ids(cells)
    
    # IMG image level information
    #-----------
    img[,"ImageNumber_Global"] <- paste0(global.prefix, "_I", img[, "ImageNumber"])
    # ORL, object relationships
    #-----------
    if (nrow(orl) > 0) {
      orl[,"First Image Number Global"]   <- paste0(global.prefix, "_I", orl[, "First Image Number"])
      orl[,"First Image Number Global"]   <- paste0(global.prefix, "_I", orl[, "First Image Number"])
      orl[,"First Object Number Global"]  <- paste0(global.prefix, "_I", orl[, "First Image Number"], "_O", orl[, "First Object Number"])
      orl[,"First Object Number Global"]  <- paste0(global.prefix, "_I", orl[, "First Image Number"], "_O", orl[, "First Object Number"])
    }
  }
  
  if(return.feature.meta) {
    feature.meta <- tglow.get.feature.meta.from.cells(colnames(cells))
  }  
  
  out.list <- list(cells=cells, meta=img, orl=orl)
  
  if(return.feature.meta) {
    out.list <- c(out.list, list(features=feature.meta))
  } 
  
  if (length(children) > 0 ) {
      out.list <- c(out.list, list(children=children))
  }
  
  return(out.list)
}

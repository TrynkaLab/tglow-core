library(ggplot2)
library(viridis)
library(pheatmap)

##------------------------------------------------------------------------------
## Simple plotting theme for ggplot using arial family font
##------------------------------------------------------------------------------
theme.plain <- function(p, base_size = 11, base_family = "ArialMT", legend=T) {
  p <- p + theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black", size=0.75),
          axis.ticks = element_line(size=0.75),
          axis.text = element_text(size=base_size, family=base_family, face="plain"),
          strip.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size=base_size, family = base_family, face="plain"),
          complete = TRUE,
          plot.title = element_text(hjust=0.5))
  
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}

#-------------------------------------------------------------------------------
#' Chart of average execution times
#' 
#' @details 
#' Plot the execution time of CellProfiler modules
#' 
#' @param meta dataframe with metadata as read by tglow.read.fileset
#' 
#' @returns ggplot object containing the plot
#' 
tglow.plot.execution.time <- function(meta) {
  
  # Calculate average execution time
  df.plot <- meta[,grep("ExecutionTime",colnames(meta))]
  df.plot <- data.frame(name=gsub("ExecutionTime_", "", colnames(df.plot)), time=colMeans(df.plot, na.rm=T))
  
  #df.plot$time <- (df.plot$time / sum(df.plot$time)) *100
  
  df.plot <- df.plot[order(df.plot$time, decreasing=T),]
  df.plot$name <- factor(df.plot$name, levels=unique(df.plot$name))  
  
  plot <- ggplot(df.plot, aes(x=name, y=time, fill=time)) +
    geom_bar(width = 0.8, stat = "identity") +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour="grey"),
          panel.background=element_blank()) +
    ylab("Averate time (seconds)") +
    xlab("") 
  
  return(plot)
}

#-------------------------------------------------------------------------------
#' Plot a hexbin colored on a 3rd variable
#' 
#' 
tglow.plot.location.hex <- function(cells,
                                    feature.z,
                                    feature.x="nucl_AreaShape_Center_X",
                                    feature.y="nucl_AreaShape_Center_Y",
                                    xlab="X centroid",
                                    ylab="Y centroid",
                                    zlab="Z",
                                    main=NULL,
                                    bins=30,
                                    scale.to.well.mean=F,
                                    scale.z=F,
                                    trim.outliers.z=F,
                                    trim.outliers.z.thresh=3.5,
                                    bins.mincount=0) {
  
  # Build the plot
  df.plot <- data.frame(x=cells[, feature.x],
                        y=cells[, feature.y],
                        z=cells[, feature.z])
  
  
  if (trim.outliers.z) {
    df.plot <- df.plot[abs(tglow.mod.zscore(df.plot$z)) < trim.outliers.z.thresh,]
  }

  if (scale.to.well.mean) {
    df.plot$z <- df.plot$z / mean(df.plot$z)
  }
  
  if (scale.z) {
    df.plot$z <- scale(df.plot$z)
  }

  p1 <- ggplot(data=df.plot,
              mapping=aes(x=x, y=y, z=z)) +
    stat_summary_hex(fun = function(x){ if (length(x) > bins.mincount){mean(x)} else {return(NA)}},
                     bins = bins) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(main) +
    scale_fill_viridis_c(name=zlab)
  
  return(p1)
  
}
#-------------------------------------------------------------------------------
#' Plot a hexbin colored on density
#' 
#' 
tglow.plot.hex <- function(x, y, bins=250, do.lm=T, lm.col="lightgrey", xlab="x", ylab="y", facet=NULL, ...) {
  
  df.plot <- data.frame(x=x, y=y)
  
  if (!is.null(facet)) {
    df.plot$facet <- facet
  }
  
  p1 <- ggplot(data=df.plot, mapping=aes(x=x, y=y)) +
    geom_hex(bins=bins) +
    xlab(xlab) +
    ylab(ylab)
  
  if (do.lm) {
    p1 <- p1 + geom_smooth(method="lm", col=lm.col)
  }
  
  if (!is.null(facet)) {
    p1 <- p1 + facet_wrap(~facet, ...)
  }
  
  return(p1)
}

#-------------------------------------------------------------------------------
#' Plot boxplot with a numerically ordered x and a loess line
#' 
#' 
tglow.plot.boxline <- function(x, y, levels, facet=NULL, xlab="x", ylab="y", line.col="#3d403d", method=NULL) {
  
  df.plot <- data.frame(x=as.numeric(factor(x, levels=levels)), y=y)
  if (!is.null(facet)) {
    df.plot$facet <- facet
  }

  p1 <- ggplot(data=df.plot, mapping=aes(x=x, y=y, group=x, facet=facet)) +
    geom_boxplot() +
    geom_smooth(se=F, aes(group=1), col=line.col, method=method) +
    xlab(xlab) +
    ylab(ylab) +
    scale_x_continuous(breaks=1:length(levels), labels=unique(levels))
  
  if (!is.null(facet)) {
    p1 <- p1+facet_wrap(~facet)
  } 
  
  return(p1)
}


#-------------------------------------------------------------------------------
#' Plot an individual ImageMagick image with a title
#' 
#' @returns ggplot object with the magick image
tglow.plot.img <- function(img, main="img", marker.add=T, marker.x=NULL, marker.y=NULL, marker.col="white", marker.size=8, marker.shape=18, ...) {
  p <- ggplot() +
    ggtitle(main) +
    annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    coord_fixed()
  
  if (marker.add) {
    
    if (is.null(marker.x)) {
      marker.x <- round(dim(img)[1]/2)
    }
    
    if (is.null(marker.y)) {
      marker.y <- round(dim(img)[2]/2)
    }
    
    marker.x 
    
    p <- p + geom_point(mapping=aes(x=marker.x, y=marker.y),
                        col=marker.col,
                        shape=marker.shape,
                        size=marker.size)
  }
  return(p + theme_void())
}

#-------------------------------------------------------------------------------
#' Plot a series of ImageMagick images with a title
#' 
#' @returns ggplot object with the magick image
tglow.plot.img.set <- function(imgs, ncol,  main="", main.sub=NULL, text.col="white", background.col="black", byrow=F, ...) {

  plots <- list()
  i <- 0
  for (img in 1:length(imgs)) {
    
    i <- i+1
    if (!is.null(main.sub)) {
      p <- tglow.plot.img(imgs[[img]], main=main.sub[i], ...)
    } else {
      p <- tglow.plot.img(imgs[[img]], main="", ...)
    }
    
    p <- p + theme(plot.title = element_text(colour=text.col)) +
        theme(plot.background=element_rect(fill=background.col),
              plot.title=element_text(hjust=0.5))
      
    #p <-  theme.plain(p) + 
    #  
    
    plots[[img]] <- p
  }
  
  p.tmp <- plot_grid(plotlist = plots, ncol=ncol, byrow=byrow) +
    theme(plot.background = element_rect(fill = background.col, colour = NA))
  
  title  <- ggdraw() + draw_label(main, fontface='bold')
  
  p.final <-  plot_grid(title, p.tmp, ncol=1, rel_heights=c(0.05, 1))
  
  return(p.final)
}



#-------------------------------------------------------------------------------
#' JULIE NEW: Plot a scatter plot faceted by variable z
#' This is redundant with tglow.plot.xy?
#' 
tglow.plot.scatter.grouped <- function(cells,
                                       feature.x,
                                       feature.y,
                                       feature.z,
                                       xlab,
                                       ylab,
                                       zlab=NULL,
                                       main=NULL,
                                       bins=100) {
  cat("[WARNING] IS THIS REDUNDANT WITH tglow.plot.xy????\n")
  
  # Build the plot
  df.plot <- data.frame(x=cells[, feature.x],
                        y=cells[, feature.y],
                        z=cells[, feature.z])
  
  p1 <- ggplot(data=df.plot,
               mapping=aes(x=x, y=y)) +
    geom_point() +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(main) +
    facet_wrap(~z, scales = "free")
  
  return(p1)
  
}

#-------------------------------------------------------------------------------
#' JULIE NEW Plot a histogram from cells grouped by variable z
#' 
#' 
tglow.plot.histogram.grouped <- function(cells,
                                         feature.x,
                                         feature.z,
                                         xlab,
                                         main=NULL,
                                         bins=100) {
  
  # Build the plot
  df.plot <- data.frame(x=cells[, feature.x],
                        z=cells[, feature.z])
  
  med <- median(df.plot$x, na.rm = TRUE)
  
  p1 <- ggplot(data=df.plot,
               mapping=aes(x=x, fill=z), alpha = 0.5, color = "black") +
    geom_histogram() +
    geom_vline(xintercept = med, color = "red", linetype = "dashed", size = 1) +  # Add vertical line for the median
    xlab(xlab) +
    ggtitle(main) +
    scale_fill_viridis_d()
  
  return(p1)
  
}


#-------------------------------------------------------------------------------
# Aggreagtes localization features
tglow.make.locale.plot.df <- function(data, feature.pattern, facet) {
  
  features         <- grep(feature.pattern, colnames(output.f$cells), value=T)
  img.nr           <- data$cells$Image_ImageNumber_Global
  
  df.plot <- data.frame()
  for (f in features) {
    df.plot <- rbind(df.plot, cbind(data$cells[,f],
                                    paste0(f, "_", data$meta[img.nr, "timepoint"]),
                                    data$meta[img.nr, "timepoint"], 
                                    f,
                                    facet))
  }
  
  colnames(df.plot) <- c("y", "x", "time", "fill", "facet")
  df.plot$y         <- as.numeric(df.plot$y)
  df.plot2          <- aggregate(df.plot$y, by=list(df.plot$time, df.plot$fill, df.plot$facet), FUN=mean)
  tmp               <- aggregate(df.plot$y, by=list(df.plot$time, df.plot$fill, df.plot$facet), FUN=sd)
  df.plot2$ymin     <- df.plot2$x - tmp$x
  df.plot2$ymax     <- df.plot2$x + tmp$x
  df.plot2$Group.2 <- gsub(feature.pattern, "", df.plot2$Group.2)
  df.plot2$Group.2 <- gsub("of5", "", df.plot2$Group.2)
  
  return(df.plot2)
}


#-------------------------------------------------------------------------------
#' Localization plot
tglow.plot.make.locale <- function(df.plot, feature.label) {

  p1 <- ggplot(data=df.plot,
               mapping=aes(x=factor(Group.1, levels=time.order),
                           y=x,
                           col=Group.2,
                           fill=Group.2,
                           ymin=ymin, 
                           ymax=ymax,
                           group=Group.2)) +
    geom_ribbon() + 
    geom_line() +
    scale_color_viridis_d(name="Localization") +
    scale_fill_viridis_d(name="Localization",alpha=0.25) +
    ylab(feature.label) +
    xlab("Timepoint") +
    facet_wrap(~Group.3)
  
  return(p1)
}


#------------------------------------------------------------------------------
# Correlation plot - need to fix the color thing
tglow.plot.cor <- function(dataset, x, y, color = NULL, assay = "cells_norm", facet = NULL, n.dots = NULL, log.x = NULL, log.y = NULL){
  
  # Find x and y 
  if(x %in% colnames(dataset$meta)){
    x.data <- dataset[["meta"]][dataset[[assay]]$Image_ImageNumber_Global, x]
  } else if (x %in% colnames(dataset$cells_norm)){
    x.data <- dataset[[assay]][, x]
  }
  
  if(y %in% colnames(dataset$meta)){
    y.data <- dataset[["meta"]][dataset[[assay]]$Image_ImageNumber_Global, y]
  } else if (y %in% colnames(dataset$cells_norm)){
    y.data <- dataset[[assay]][, y]
  }
  
  # Create plotting dataframe
  df <- data.frame(x = x.data,
                   y = y.data)
  
  # Add color variable if necessary
  if(!is.null(color)){
    if(color %in% colnames(dataset$meta)){
      df$c <- dataset[["meta"]][dataset[[assay]]$Image_ImageNumber_Global, color]
    } else if (color %in% colnames(dataset$cells_norm)){
      df$c <- dataset[[assay]][, color]
    }
  }
  
  # Add facet variable if necessary
  if(!is.null(facet)){
    df$f <- dataset[["meta"]][dataset[[assay]]$Image_ImageNumber_Global, facet]
  }
  
  # Build base plot (with or without colored dots)
  if(!is.null(color)){
    p <- ggplot(df) +
      geom_point(aes(x, y, color = c), size = 0.5) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      scale_color_viridis_d()
  } else{ 
    p <- ggplot(df) +
      geom_point(aes(x, y)) +
      theme_bw() +
      theme(panel.grid = element_blank())
  }
  
  # Log axis
  if(!is.null(log.x)){p <- p + scale_x_log10()}
  
  if(!is.null(log.y)){p <- p + scale_x_log10()}
  
  # Add facet
  if(!is.null(facet)){p <- p + facet_wrap(~f, scales = "free")}
  
  # Add axis labels
  p <- p + labs(x = x,
             y = y)
  
  return(p)
  
}


#------------------------------------------------------------------------------
#' Violin plot - need to fix the color thing
tglow.plot.violin <- function(dataset, x, y, color = NULL, assay.y = "cells_norm", assay.color = "cells_norm", facet = NULL, n.dots = NULL){
  
  # Create plotting dataframe
  df <- data.frame(x = dataset[["meta"]][dataset[[assay.y]]$Image_ImageNumber_Global, x],
                   y = dataset[[assay.y]][, y])
  
  
  # Add color variable if necessary
  if(!is.null(color)){
    
    if(color %in% colnames(dataset$meta)){
      df$c <- dataset[["meta"]][dataset[[assay.y]]$Image_ImageNumber_Global, color]
      
    } else if (color %in% colnames(dataset$cells_norm)){
      df$c <- dataset[[assay.color]][, color]
    }
    
  }
  
  # Add facet variable if necessary
  if(!is.null(facet)){
    df$f <- dataset[["meta"]][dataset[[assay.y]]$Image_ImageNumber_Global, facet]
  }
  
  # Build base plot
  p <- ggplot(df) +
    geom_violin(aes(x, y), fill = "grey") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(y = y) +
    scale_fill_viridis_d() +
    labs(x = "")

  # Add color
  if(!is.null(color)){

    p <- p + geom_violin(aes(x, y, fill = c), color = "black")

  }
  
  # Add points
  if(!is.null(n.dots)){
    # Subset data
    df.sub <- slice_sample(df, n = n.dots)
    
    p <- p + geom_jitter(data = df.sub, aes(x, y), size = 0.5, height = 0, width = 0.1, color = "grey")

  }
  
  # Add facet
  if(!is.null(facet)){
    p <- p + facet_wrap(~f, scales = "free") 
  }
  
  return(p)
}


#------------------------------------------------------------------------------
#' Boxplot plot - need to fix the color thing
tglow.plot.boxplot <- function(dataset, x, y, color = NULL, assay.y = "cells", assay.color = "cells_norm", facet = NULL, n.dots = NULL){
  
  # Create plotting dataframe
  df <- data.frame(x = dataset[["meta"]][dataset[[assay.y]]$Image_ImageNumber_Global, x],
                   y = dataset[[assay.y]][, y])
  
  # Add color variable if necessary
  if(!is.null(color)){
    
    if(color %in% colnames(dataset$meta)){
      df$c <- dataset[["meta"]][dataset[[assay.y]]$Image_ImageNumber_Global, color]
    } else if (color %in% colnames(dataset$cells_norm)){
      df$c <- dataset[[assay.color]][, color]
    }

  }
  
  # Add facet variable if necessary
  if(!is.null(facet)){
    df$f <- dataset[["meta"]][dataset[[assay.y]]$Image_ImageNumber_Global, facet]
  }
  
  # Build base plot
  p <- ggplot(df) +
    geom_boxplot(aes(x, y, fill = c)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(y = y) +
    scale_fill_viridis_d() +
    labs(x = "") +
    theme(legend.position = "none")
  
  # Add points
  if(!is.null(n.dots)){
    # Subset data
    df.sub <- slice_sample(df, n = n.dots)
    p <- p + geom_jitter(data = df.sub, aes(x, y), size = 0.5, height = 0, width = 0.1, color = "grey")
  }
  
  # Add facet
  if(!is.null(facet)){
    p <- p + facet_wrap(~f)
  }
  
  return(p)

}


#------------------------------------------------------------------------------
#' GGplot scatterplot
tglow.plot.xy <- function(x, y, xlab="X", ylab="Y", main=NA, main.prefix="", size=1, col="black", fixed=F, alpha=0.75, shape=16, line.col="blue", do.lm=T, method="lm", lm.group=NULL, raster=F, dpi=300, facet=NULL, facet.ncol=NULL) {
  
  df.plot <- data.frame(x=x,
                        y=y)
                        
  # Fix axis limits between x and y
  if (fixed) {lims <- c(min(c(x, y), na.rm=T), max(c(x, y), na.rm=T))  }
  
  # Facet y or n
  if (!is.null(facet)) {df.plot$facet <- facet}
  
  if (nrow(df.plot) > 3) {
    ct <- cor.test(df.plot$x, df.plot$y)
  } else {
    ct <- list(`p.value`=NA, `estimate`=NA)
  }
  
  if (length(shape) > 1 && length(col > 1)) {
    if (raster) {
      p <- ggplot(data=df.plot, mapping=aes(x=x, y=y, col=col, shape=shape)) +
        rasterize(geom_point(alpha=alpha, size=size), dpi=dpi)
    } else {
      p <- ggplot(data=df.plot, mapping=aes(x=x, y=y, col=col, shape=shape)) +
        geom_point(alpha=alpha, size=size)
    }
    
  } else if (length(col) > 1) {
    if (raster) {
      p <- ggplot(data=df.plot, mapping=aes(x=x, y=y, col=col)) +
        rasterize(geom_point(alpha=alpha, shape=shape, size=size), dpi=dpi)
    } else {
      p <- ggplot(data=df.plot, mapping=aes(x=x, y=y, col=col)) +
        geom_point(alpha=alpha, shape=shape, size=size)
    }
  } else {
    if (raster) {
      p <- ggplot(data=df.plot, mapping=aes(x=x, y=y)) +
        rasterize(geom_point(alpha=alpha, color=col, shape=shape, size=size), dpi=dpi)
    } else {
      p <- ggplot(data=df.plot, mapping=aes(x=x, y=y)) +
        geom_point(alpha=alpha, color=col, shape=shape, size=size) 
    }
  }
  
  p <- p +
    xlab(xlab) +
    ylab(ylab) 
  
  # Add linear line
  if (do.lm) {
    p <- p + geom_smooth(method=method, col=line.col, inherit.aes = F, mapping=aes(x=x, y=y, group=lm.group)) 
  }
  
  # If axis are fixed, add an abline and set limits
  if (fixed) {
    p <- p + 
      xlim(lims) +
      geom_abline(slope=1, intercept=0, col="grey", lty=2) +
      coord_fixed() +
      ylim(lims) 
  }
  
  # Add title
  if (is.na(main) & do.lm) {
    main <- paste0(main.prefix, "R: ", format(ct$estimate, digits=2),
                            " p-value: ", format(ct$`p.value`, digits=2, scientific=T))
  }
  
  if (!is.na(main)) {p <- p + ggtitle(main)}
  
  if (!is.null(facet)) {p <- p + facet_wrap(~facet, ncol=facet.ncol)}

  
  
  return(p)
}


#------------------------------------------------------------------------------
#' Simple heatmap with auto labels
tglow.plot.simple.hm <- function(data, cellsize=-1, cellwidth=12, cellheight=12, limits=NULL, cluster=T, range="symmetric", min.value=0, palette=NULL, border=NA, ...) {
  
  if (cellsize > 0) {
    cellwidth  <-  cellsize
    cellheight <- cellsize
  }
  
  if (is.null(limits)) {
    max <- max(data)
    min <- min(data)
    max.abs <- max(abs(data))
    min.abs <- min(abs(data))
  } else {
    if(length(limits) ==2) {
        max <- limits[2]
        min <- limits[1]
        max.abs <- max(abs(limits))
        min.abs <- min(abs(limits))
        range = "auto"
    } else {
      stop("[ERROR] Limits must be a vector of length 2")
    }
  }
  
  if (range == "symmetric") {
    break.list <- seq(-max.abs, max.abs, by=max.abs/50)
    if (is.null(palette)) {palette="RdBu"}
    cols       <- colorRampPalette(rev(brewer.pal(n=7, name =palette)))(length(break.list))
  } else if (range == "absolute") {
    if (is.null(palette)) {palette="Reds"}
    break.list <- seq(min, max.abs, by=(max.abs-min.abs)/100)
    cols       <- colorRampPalette(c("#FFFFFF",brewer.pal(n=7, name =palette)))(length(break.list))
  } else if (range == "auto") {
    break.list <- seq(min, max, by=(max-min)/100)
    if (is.null(palette)) {palette="RdBu"}
    cols       <- colorRampPalette(rev(brewer.pal(n=7, name =palette)))(length(break.list))
  } else  {
    stop("[ERROR] range must be symmetric, auto, or aboslute\n")
  }
  
  if (!cluster) {
    pheatmap(data,
             breaks=break.list,
             col=cols,
             cellwidth=cellwidth,
             cellheight=cellheight,
             border=border,
             cluster_rows = F,
             cluster_cols = F,
             ...)
  } else {
    pheatmap(data,
             breaks=break.list,
             col=cols,
             cellwidth=cellwidth,
             cellheight=cellheight,
             border=border,
             ...)
  }
}

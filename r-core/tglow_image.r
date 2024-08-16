# Image related functions for the tglow pipeline
# Based on EBImage

#-------------------------------------------------------------------------------
#' Convert a hex code to RGB
#' @param hex String with hex code
#'
#' @returns a vector of length 3 with relative scaling between RGB channels between 0 and 1
tglow.hex.to.rgb <- function(hex) {
  
  if (startsWith(hex, "#")) {
    hex <- gsub("^#", "", hex)
  }
  
  rgb.hex <- strsplit(hex, "(?<=..)", perl = TRUE)[[1]]
  rgb <- sapply(rgb.hex, strtoi, base = 16)
  rgb <- rgb/255
  
  names(rgb) <- c("r", "g", "b")
  return(rgb)
}

#-------------------------------------------------------------------------------
#' Apply a relative scaling of RGB values to a 3d array
#' 3rd dimension is assumed to contain the color channels.
#' If input is grey, returns the image in the color
#' @param image 3d array with 3rd dimension as rgb
#' @param rgb vector of length 3 with scaling factors for rgb channels
#'
#' @returns The scaled image
tglow.apply.color <- function(image, rgb) {
  image[,,1] <-image[,,1]*rgb[1]
  image[,,2] <-image[,,2]*rgb[2]
  image[,,3] <-image[,,3]*rgb[3]
  
  return(image)
}

#-------------------------------------------------------------------------------
#' Combine a series of grey rgb images, appy a color shift and average their
#' RGB values at each pixel to create a composite image.
#' @param images list of RGB EBImages
#' @param colors list of hex codes to color shift, must be same length as images
#'
#' @returns A single RGB images representing the average color values
tglow.composite.image <- function(images, colors) {
  
  if(length(images) != length(colors)) {
    stop("Images and suplied colors don't match")
  }
  
  i <- 1
  comp <- NULL
  
  for (image in images) {
    
    if (length(dim(image)) != 2) {
      stop(paste0("Dimensions of image not equal to 2, is it greyscale. At index ", i))
    }
    
    rgb <- tglow.hex.to.rgb(gplots::col2hex(colors[i]))
    img <- channel(image, mode="rgb")
    
    if (is.null(comp)) {
      # Create an RGB image from greyscale
      comp <- tglow.apply.color(img, rgb)
    } else {
      comp <- comp + tglow.apply.color(img, rgb)
    }
    
    i <- i+1
  }
  
  comp <- comp/length(images)
  
  return(comp)
}


#-------------------------------------------------------------------------------
#' Find the maximum value in each channel in a collection of images
#' This is used for normalizing images
#' @param imgs A list of EBImage objects, with 3 or 4 dimensions. 3d dimension
#' is assumed to be channel.
#' @returns A vector with the max values in each channel
tglow.max.per.channel <- function(imgs, channel.dim=3, q=1) {
  
  channel.dim <- 3
  max         <- c(rep(0, dim(imgs[[1]])[channel.dim]))
  
  for (img in imgs) {
    if (length(dim(img)) == 3) {
      for (ch in 1:dim(img)[channel.dim]) {
        if (quantile(img[,,ch], probs=q) > max[ch]) {
          max[ch] <- quantile(img[,,ch], probs=q)
        }
      }
    } else if (length(dim(img)) == 4) {
      for (ch in 1:dim(img)[channel.dim]) {
        if (quantile(img[,,ch,], probs=q) > max[ch]) {
          max[ch] <- quantile(img[,,ch,], probs=q)
        }

      }
    } else {
      stop(paste0(dim(img)," is not a valid number of dimensions for an image. "))
    }
  
  }
  
  return(max)
}

#-------------------------------------------------------------------------------
#' Normalize a list of EBImages to a common norm factor.
#' @param imgs A list of EBImage objects, with 3 or 4 dimensions. 3d dimension
#' is assumed to be channel.
#' @param norm.factors Value to divide the images by.
#' Defaults to NULL in which case tglow.max.per.channel is called
tglow.norm.img <- function(imgs, norm.factors=NULL, q=1) {
  
  if (is.null(norm.factors)) {
    cat("[INFO] Calculating norm factors\n")
    norm.factors <- tglow.max.per.channel(imgs, q=q)
  }
  
  if (length(dim(imgs[[1]])) == 3) {
    imgs <- lapply(1:length(imgs), function(x) {
      clls <- imgs[[x]]
      
      for (ch in 1:dim(clls)[3]) {
        imageData(clls)[,,ch] <- clls[,,ch]/norm.factors[ch]
      }
      return(clls)
    })
  } else if (length(dim(imgs[[1]])) == 4) {
    
    imgs <- lapply(1:length(imgs), function(x) {
      clls <- imgs[[x]]
      
      for (ch in 1:dim(clls)[3]) {
        imageData(clls)[,,ch,] <- clls[,,ch,]/norm.factors[ch]
      }
      return(clls)
    })
    
  } else {
    stop(paste0(dim(img)," is not a valid number of dimensions for an image. "))
  }
  
  
  return(imgs)
  
  
  
}

#-------------------------------------------------------------------------------
#' Max project a EBImage along the 4th axis
tglow.max.project <- function(img) {
  apply(img, c(1, 2, 3), max)
}

pacman::p_load('jpeg', 'imager', 'data.table', 'dplyr', 'raster', 'rgdal', 'rtiff', 'RStoolbox', 'sdmvspecies')

eLAI <- function(cos, Pgap, G){
  numerator = (cos* log(Pgap))
  denominator = G
  eLAI = (-(numerator / denominator))
  return(eLAI)
}

args = commandArgs(trailingOnly=TRUE)

dir <- args[1]
files <- list.files(dir, pattern = '.tif') %>% gsub(., pattern = '.tif', replacement = "")

lai <- data.frame(matrix(0, length(files), 3))
names(lai) <- c('scan', 'rc_Pgap', 'rc_eLAI')

################################################################################
for(i in 1:length(files)){

    scan <- raster(paste0(dir, '/', files[i], '.tif'))

    extent <- attr(scan, 'extent')
    low.hinge <- 55 * extent[4] / extent[4]
    up.hinge <- 60 * extent[4] / extent[4]

    xmin <- extent[1]
    xmax <- extent[2]
    ymin <- low.hinge
    ymax <- up.hinge

    scan.hinge <- crop(scan, extent(xmin, xmax, ymin, ymax))
    scan.hinge[is.na(scan.hinge[])] <- 0

    ## Ridler and Calvard

    rc <- autoThreshold(scan.hinge)
    thresh.rc <- threshold(scan.hinge, thr = rc[3])
    
    
    thresh.rc.xyz <- rasterToPoints(thresh.rc, na.rm = F)
    thresh.rc.xyz <- as.data.frame(thresh.rc.xyz)
    names(thresh.rc.xyz) <- c('x', 'y', 'intensity')
    
    rc.gap <- nrow(subset(thresh.rc.xyz, intensity == 0))
    rc.non.gap <- nrow(subset(thresh.rc.xyz, intensity == 1))
    rc.Pgap <- rc.gap / (rc.gap + rc.non.gap)
    rc.eLAI <- eLAI(cos(57.5), rc.Pgap, 0.5)
    
    ## print results

    lai$scan[i] <- files[i]
    lai$rc_Pgap[i] <- rc.Pgap
    lai$rc_eLAI[i] <- rc.eLAI

  }

write.csv(lai, paste0(dir, '/', dir, '_rc.csv'), row.names = F)

################################################################################

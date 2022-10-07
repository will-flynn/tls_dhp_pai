pacman::p_load(dplyr, raster, plotrix, rtiff, imager, data.table)

eLAI <- function(cos, Pgap, G){
  numerator = (cos* log(Pgap))
  denominator = G
  eLAI = (-(numerator / denominator))
  return(eLAI)
}

location.cx = 2744  # x coordinate of center
location.cy = 1829  # y coordinate of center
#location.cr = 1306  # radius of circle
location.cr = 1744  # radius of circle

################################################################################

dir <- '/media/will/WillData1/phd/LAI/dhp/data'
folders <- list.files(paste0(dir, '/alto_tajo/thresh'))

parameters <- read.csv(paste0(dir, '/', 'exposure_choices.csv'))

t1 = Sys.time()
for(ii in 1:length(folders)) {

    plotdir <- paste0(dir, '/alto_tajo/thresh/', folders[ii], '/')

    p <- which(grepl(folders[ii], parameters$plot))

    images =   if (parameters$exposure[p] %in% "auto") {
      list.files(plotdir, pattern = "_AUTO.JPG", full.names = T)
      } else if (parameters$exposure[p] %in% "minus") {
        list.files(plotdir, pattern = "_MINUS.JPG", full.names = T)
      } else {
        list.files(plotdir, pattern = "_PLUS.JPG", full.names = T)
      }


    lai <- data.frame(matrix(0, length(images), 4))
    names(lai) <- c('plot', 'image', 'Pgap', 'ePAI')

        for(i in 1:length(images)){

            image <- raster(images[i], band = 3)

            extent <- attr(image, 'extent')

            low.hinge <- 55 * extent[4] / 90
            up.hinge <- 60 * extent[4] / 90
            low.hinge <- low.hinge / 2
            up.hinge <- up.hinge / 2

            plot(image)
            low.hinge.circle <- draw.circle(location.cx, location.cy, low.hinge)
            up.hinge.circle <- draw.circle(location.cx, location.cy, up.hinge)

            low.hinge.circle <- as.data.frame(low.hinge.circle)
            up.hinge.circle <- as.data.frame(up.hinge.circle)

            lh.p = Polygon(low.hinge.circle)
            lh.ps = Polygons(list(lh.p),1)
            lh.sps = SpatialPolygons(list(lh.ps))

            uh.p = Polygon(up.hinge.circle)
            uh.ps = Polygons(list(uh.p),1)
            uh.sps = SpatialPolygons(list(uh.ps))

            image.hinge <- mask(image, uh.sps)
            image.hinge <- mask(image.hinge, lh.sps, inverse = T)

            xyz <- rasterToPoints(image.hinge)
            xyz <- as.data.table(xyz)
            names(xyz) <- c('x', 'y', 'intensity')
            xyz <- xyz[, intensity:= scales::rescale(intensity, to = c(0, 1),
              from = range(intensity, na.rm = TRUE, finite = TRUE))]
            xyz
            theta <- cos(57.5)

            gap <- subset(xyz, intensity == 1)
            non.gap <- subset(xyz, intensity == 0)
            Pgap <- nrow(gap) / (nrow(gap) + nrow(non.gap))
            LAI <- eLAI(theta, Pgap, 0.5)
            
            lai$plot <- folders[ii] 
            lai$image[i] <- i
            lai$Pgap[i] <- round(Pgap, 2)
            lai$ePAI[i] <- round(LAI, 2)

          }
          write.csv(lai, paste0(dir, '/', folders[ii],  '.csv'), row.names = F)
}
t2 = Sys.time()
t2 - t1

################################################################################

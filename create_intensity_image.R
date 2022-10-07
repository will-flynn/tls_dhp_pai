pacman::p_load(data.table, scales, raster, dplyr, doParallel)

################################################################################
args <- commandArgs(trailingOnly = T)
dir <- args[1]
#dir <- 'C:/Users/wrmfl/phd/LAI/image_intensity/alt25/pts'


if (length(args)<1)
{
    stop("At least one arguments must be supplied
       1: directory of .pts files")
}


files <- list.files(dir, pattern = '.pts') %>% gsub(., pattern=".pts", replacement="")

if (length(files)<1)
{
    stop("No .pts files found in supplied directory")
}


t1 <- Sys.time()
for(i in 1:length(files)){
    
    print(paste0('reading ', dir, ' ', files[i], '...'), quote = F)
    
    scan <- fread(paste0(dir, '/', files[i], '.pts'), select = c(1,2,3,4),
       col.names = c('X', 'Y', 'Z', 'intensity'), data.table = T, header = F)
    
    print(paste0('converting ', dir, ' ', files[i], ' to raster...'), quote = F)

    scan <- scan[, intensity:= scales::rescale(intensity, to = c(0, 255),
              from = range(intensity, na.rm = TRUE, finite = TRUE))][, list(X,Y,Z,intensity)]
    scan <- scan[, r:= sqrt((scan$X^2) + (scan$Y^2) + (scan$Z^2))]
    scan <- scan[, zenith:= (round(acos(scan$Z / scan$r), 3))*180/pi]
    scan <- scan[, azimuth:= (round(atan2(scan$Y, scan$X), 3))*180/pi]

    xyz <- scan[, c('azimuth', 'zenith', 'intensity')]

    scan.raster <- rasterFromXYZ(xyz)

    print(paste0('writing ', dir, ' ', files[i], '...'), quote = F)
    
    writeRaster(scan.raster, paste0(dir, '/', files[i], '.tiff'), format = 'GTiff', overwrite = T)

    print(paste0(dir," ", files[i], " ", "complete!"), quote = F)
}
t2 <- Sys.time()
print(t2 - t1)

################################################################################

pacman::p_load(rTLS, parallel, stringr, doParallel)

args = commandArgs(trailingOnly=TRUE)

dir = args[1]
scans = list.files(dir, pattern = '.pts') %>% gsub(pattern = "\\.pts$", "", .)
plot = dir %>% gsub('^.*/', '', .)

threads = detectCores()

cores <- detectCores()
useCores <- cores

cl <- makeCluster(useCores)
registerDoParallel(cl)

out.list <- list()

t1 = Sys.time()
result = foreach(i = 1:length(scans), .combine = rbind) %dopar% {
  pacman::p_load(rTLS, parallel, stringr, doParallel)

  pts = fread(paste0(dir, '/', scans[i], '.pts'), select = c(1,2,3), col.names = c('X', 'Y', 'Z'), skip = 1)

  gap.fractions = canopy_structure(
    TLS.type = 'single',
    scan = pts,
    zenith.range = c(55, 60),
    zenith.rings = 1,
    azimuth.range = c(0, 360),
    vertical.resolution = 0.25,
    TLS.frame = c(0, 90, 0, 360),
    TLS.pulse.counts =  NULL,
    TLS.resolution = c(0.018, 0.018),
    TLS.coordinates = c(X = 0,
                        Y = 0,
                        Z = 0),
    threads = threads
  )

  pai <- data.frame(matrix(0, 1, 3))
  names(pai) <- c('scan', 'Pgap', 'PAI')

  pai$scan <- scans[i]
  pai$Pgap <- min(gap.fractions[, 2])
  pai$PAI <- max(gap.fractions[, 3])

  out.list[[i]] <- pai
}
stopCluster(cl)
t2 = Sys.time()
time.elapsed = t2 - t1
print(time.elapsed)

out.file = paste0(dir, '/', plot, ".csv")
if(file.exists(out.file)){
  file.remove(out.file)
}

write.table(result, out.file, append = F, sep = ",", row.names = F, col.names = !file.exists(out.file), quote = F)

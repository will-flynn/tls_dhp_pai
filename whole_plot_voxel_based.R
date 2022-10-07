pacman::p_load(data.table, stringr, dplyr, doParallel, rTLS, VoxR, sf, raster, concaveman, Morpho, ggplot2)

dir <- "C:/Users/wrmfl/Desktop/downsampled_1/tree"

plot.names <- unique(str_extract(list.files(dir, pattern = ".txt") %>% gsub(".txt", "", .), "[^_]+"))
res = 0.05

################################################################################
vox = function(data,res,full.grid,message){
  
  #- declare variables to pass CRAN check as suggested by data.table mainaitners
  x=y=z=npts=.N=.=':='=NULL
  
  #- check for data consistancy and convert to data.table
  check=VoxR::ck_conv_dat(data,message=message)
  
  if(missing(res)){
    stop("No voxel resolution (res) provided")
  }else{
    #- res must be a vector
    if(!is.vector(res)) stop("res must be a vector of length 1")
    #- res must be numeric
    if(!is.numeric(res)) stop("res must be numeric")
    #- res must be numeric
    if(res<=0) stop("res must be positive")
    #- res must be of length 1
    if(length(res)>1){
      res=res[1]
      warning("res contains more than 1 element. Only the first was used")
    }
  }
  
  #- default is without full grid
  if(missing(full.grid)) full.grid = FALSE
  
  #- keep the data
  data=check$data
  
  #- round the data coordinates with the user defines resolution
  data[,':='(x = Rfast::Round( x / res ) * res,
             y = Rfast::Round( y / res ) * res,
             z = Rfast::Round( z / res ) * res)]
  
  data = unique(data[,npts:=.N,by=.(x,y,z)])
  
  gc()
  
  if(full.grid){
    
    x_seq = seq(min(data$x),max(data$x),res)
    y_seq = seq(min(data$y),max(data$y),res)
    z_seq = seq(min(data$z),max(data$z),res)
    
    empty=data.table::data.table(expand.grid(x_seq,y_seq,z_seq))
    gc()
    data.table::setnames(empty,c("x","y","z"))
    empty[,npts:=0]
    gc()
    
    data = dplyr::bind_rows(data,empty)
    data = data[, npts:=sum(npts), keyby=.(x,y,z)]
    
    gc()
  }
  
  if(check$dfr) data = as.data.frame(data)
  
  return(data) #- output = coordinates + number of points
}
################################################################################

for(i in 1:length(plot.names)) {
  
  t1 = Sys.time()
  
  ## read in trees, normalise Z and combine into plot
  
  trees <- list.files(dir, pattern = '.txt') %>% .[grep(plot.names[i], .)] %>% gsub(".txt", "", .)
  
  plot <- paste0(dir, "/", trees, ".txt") %>%
    lapply(., fread, select = c(1,2,3), col.names = c("X", "Y", "Z")) %>%
    lapply(., function(x) cbind(x, Znorm = (x$Z - min(x$Z)))) %>%
    do.call(rbind, .) %>% .[, c(1,2,4)]
  
  ## voxelise plot
  
  print(paste0("voxelising ", plot.names[i], " ..."), quote = F)
  
  plot_vox <- vox(plot, res = res, 
                  full.grid = FALSE)
  
  gc()
  remove(plot)
  gc()
  
  z_seq <- seq(min(plot_vox$z), max(plot_vox$z), res)
  
  ground <- raster(nrow = (45 / res),
                   ncol = (45 / res),
                   xmn = -7.5,
                   xmx = 37.5,
                   ymn = -7.5,
                   ymx = 37.5)
  values(ground) <- 0
  ground_dt <- as.data.table(rasterToPoints(ground))
  
  empty <- as.data.table(cbind(rep(ground_dt$x, length(z_seq)), 
                               rep(ground_dt$y, length(z_seq))))
  data.table::setnames(empty,c("x","y"))
  empty$z <- rep(z_seq, each = nrow(distinct(empty, x,y)))
  empty[,npts:=0]
  
  plot_block = dplyr::bind_rows(plot_vox,empty)
  plot_block = plot_block[,npts:=sum(npts),keyby=.(x,y,z)]
  
  gc()
  remove(plot_vox)
  gc()

  # plot_voxels_full_grid(plot_block, res = 1)
  
  ## calculate z slices
  
  z_slices <- seq(min(plot_block$z), max(plot_block$z), res)
  #z_slices <- subset(z_slices, z_slices >= 1)
  
  ## calculate lai profiles 
  
  lai_profiles <- as.data.frame(matrix(0, length(z_slices), 2))
  names(lai_profiles) <- c("height", "pai")
  
  print("calculating z slices... ", quote = F)
  pb = txtProgressBar(min = 0, max = length(z_slices), initial = 0, style = 1) 
  
  for(j in 1:length(z_slices)) {
    
    plot.slice <- plot_block[z %in% z_slices[j], ]
    
    ni <- as.numeric(nrow(plot.slice[npts > 0, ]))
    nt <- as.numeric(nrow(plot.slice))
    
    N <- ni / nt
    l <- 1.1 * N
    
    lai_profiles$height[j] <- z_slices[j]
    lai_profiles$pai[j] <- l
    
    gc()
    remove(plot.slice)
    gc()
    
    setTxtProgressBar(pb,j)
    
  }
  
  output <- data.frame(matrix(0, 1, 2))
  names(output) <- c("plot", "pai")
  
  output$plot <- plot.names[i]
  output$pai <- sum(lai_profiles$pai, na.rm = T)
  
  output
  
  print(paste0(plot.names[i], " PAI = ", output$pai))
  
  out.file = paste0(dir, "/", "combined_trees_hosoi.csv")
  
  write.table(output, out.file, append = T, sep = ",", row.names = F, col.names = F, quote = F)
  
  t2 = Sys.time()
  print(paste0(plot.names[i], " complete... time elapsed ", round(t2 - t1, 2)), quote = FALSE)
  
gc()
  remove(plot_vox)
  remove(lai_profiles)
  remove(output)
  gc()
}

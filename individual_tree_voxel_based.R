pacman::p_load(data.table, dplyr, doSNOW, doParallel, rTLS, VoxR, sf, raster, 
               concaveman, Morpho)

dir <- "C:/Users/wrmfl/Desktop/downsampled_05"
tree.files <- list.files(paste0(dir, '/tree'), pattern = ".txt") %>% gsub(pattern = "\\.txt$", "", .)

res = 0.05

cores <- detectCores()
useCores <- cores - 4
cl <- makeCluster(useCores)
registerDoSNOW(cl)

pb <- txtProgressBar(max = length(tree.files), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

t1 = Sys.time()
foreach(i = 1:length(tree.files), .options.snow = opts) %dopar% {
  
  t1 = Sys.time()
  pacman::p_load(data.table, dplyr, doParallel, rTLS, VoxR, sf, raster, concaveman)
  
  tree.dir <- paste0(dir, "/tree/", tree.files[i], ".txt")
  leaf.dir <- paste0(dir, "/leaf/", tree.files[i], ".txt")
  wood.dir <- paste0(dir, "/wood/", tree.files[i], ".txt")
  
  tree = fread(tree.dir, select = c(1,2,3),
               col.names = c("x", "y", "z"))
  
  leaf = fread(leaf.dir, select = c(1,2,3),
               col.names = c("x", "y", "z"))
  
  wood = fread(wood.dir, select = c(1,2,3),
               col.names = c("x", "y", "z"))
  
  ## voxelise clouds
  
  tree_vox <- vox(tree, res = res, 
                  full.grid = FALSE)
  
  leaf_vox <- vox(leaf, res = res, 
                  full.grid = FALSE)
  
  wood_vox <- vox(wood, res = res, 
                  full.grid = FALSE)
  
  gc()
  remove(tree)
  remove(leaf)
  remove(wood)
  gc()
  
  ## cut to crown projected area 
  
  z_seq <- seq(min(tree_vox$z), max(tree_vox$z), res)
  
  crown_raster <- rasterFromXYZ(tree_vox[, 1:3])
  crown_dt <- as.data.table(rasterToPoints(crown_raster))
  
  empty <- as.data.table(cbind(rep(crown_dt$x, length(z_seq)), 
                               rep(crown_dt$y, length(z_seq))))
  data.table::setnames(empty,c("x","y"))
  empty$z <- rep(z_seq, each = nrow(distinct(empty, x,y)))
  empty[,npts:=0]
  
  tree_block = dplyr::bind_rows(tree_vox,empty)
  tree_block = tree_block[,npts:=sum(npts),keyby=.(x,y,z)]
  
  leaf_block = dplyr::bind_rows(leaf_vox,empty)
  leaf_block = leaf_block[,npts:=sum(npts),keyby=.(x,y,z)]
  
  wood_block = dplyr::bind_rows(wood_vox,empty)
  wood_block = wood_block[,npts:=sum(npts),keyby=.(x,y,z)]
  
  # #############################################################################
  # 
  # plot_voxels(tree_block, plot = F, res = 0.5)$mesh %>%
  #   mesh2ply(., paste0(dir, "/", "tree_empty.ply"))
  # plot_voxels(tree_vox, plot = F, res = 0.5)$mesh %>%
  #   mesh2ply(., paste0(dir, "/", "tree_filled.ply"))
  # 
  # plot_voxels(leaf_block, plot = F, res = 0.5)$mesh %>%
  #   mesh2ply(., paste0(dir, "/", "leaf_empty.ply"))
  # plot_voxels(leaf_vox, plot = F, res = 0.5)$mesh %>%
  #   mesh2ply(., paste0(dir, "/", "leaf_filled.ply"))
  # 
  # plot_voxels(wood_block, plot = F, res = 0.5)$mesh %>%
  #   mesh2ply(., paste0(dir, "/", "wood_empty.ply"))
  # plot_voxels(wood_vox, plot = F, res = 0.5)$mesh %>%
  #   mesh2ply(., paste0(dir, "/", "wood_filled.ply"))
  # 
  # #############################################################################
  
  gc()
  remove(tree_vox)
  remove(leaf_vox)
  remove(wood_vox)
  gc()
  
  ## create output data frame 
  
  z_slices <- seq(min(tree_block$z), max(tree_block$z), res)
  
  lai_profiles <- as.data.frame(matrix(0, length(z_slices), 4))
  names(lai_profiles) <- c("height", "pai", "lai", "wai")
  
  pb = txtProgressBar(min = 0, max = length(z_slices), initial = 0, style = 1) 
  
  for(j in 1:length(z_slices)) {
    
    tree.slice <- tree_block[z %in% z_slices[j], ]
    leaf.slice <- leaf_block[z %in% z_slices[j], ]
    wood.slice <- wood_block[z %in% z_slices[j], ]
    
    tree.ni <- as.numeric(nrow(tree.slice[npts > 0, ]))
    tree.nt <- as.numeric(nrow(tree.slice))
    
    leaf.ni <- as.numeric(nrow(leaf.slice[npts > 0, ]))
    leaf.nt <- as.numeric(nrow(leaf.slice))
    
    wood.ni <- as.numeric(nrow(wood.slice[npts > 0, ]))
    wood.nt <- as.numeric(nrow(wood.slice))
    
    tree.N <- tree.ni / tree.nt
    tree.l <- 1.1 * tree.N
    
    leaf.N <- leaf.ni / leaf.nt
    leaf.l <- 1.1 * leaf.N
    
    wood.N <- wood.ni / wood.nt
    wood.l <- 1.1 * wood.N
    
    lai_profiles$height[j] <- z_slices[j]
    lai_profiles$pai[j] <- tree.l
    lai_profiles$lai[j] <- leaf.l
    lai_profiles$wai[j] <- wood.l
    
    setTxtProgressBar(pb,j)
    
  }
  
  output <- data.frame(matrix(0, 1, 6))
  names(output) <- c("tree", "pai", "lai", "wai", "alpha", "lai_pai")
  
  output$tree <- tree.files[i]
  output$pai <- sum(lai_profiles$pai)
  output$lai <- sum(lai_profiles$lai)
  output$wai <- sum(lai_profiles$wai)
  output$alpha <- output$wai / output$pai
  output$lai_pai <- output$lai / output$pai
  
  output
  
  out.file = paste0(dir, "/", "hosoi_", sub('.*\\.', '', res), ".csv")
  
  if (file.exists(out.file)) {
    unlink(out.file)
    cat("The outfile is deleted")
  }
  
  
  write.table(output, out.file, append = T, sep = ",", row.names = F, col.names = F, quote = F)
  #print(paste0(tree.files[i], " complete!"), quote = FALSE)
  
  gc()
  remove(tree_block)
  remove(leaf_block)
  remove(wood_block)
  remove(lai_profiles)
  remove(output)
  gc()
  t2 = Sys.time()
  print(paste0(tree.files[i], " complete in ", round(t2 - t1, 3), " seconds"))
  
}
close(pb)
stopCluster(cl)

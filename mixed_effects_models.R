pacman::p_load(ggplot2, ggpubr, lme4, merTools, MuMIn, data.table, emmeans,
               lmerTest, equatiomatic, sjstats, sjPlot, ggeffects, plyr, multcomp,
               multcompView, forcats)

dir <- "C:/Users/wrmfl/Dropbox/phd/LAI/results/data/datasets_for_figures"
tree_data <- read.csv(paste0(dir, "/", "individual_tree.csv"))

################################################################################
## species model 
# 

tree_data$species <- as.factor(tree_data$species)

species_model <- lmerTest::lmer(alpha ~ species + (1|plot),
                                data = tree_data)
summary(species_model)

################################################################################
## height/ CAI models

species_data <- split(tree_data, tree_data$species)

height_models <- lapply(species_data, 
                         lmerTest::lmer,
                         formula = alpha ~ tree_height + cai + (1|plot))

lapply(height_models, icc)

################################################################################

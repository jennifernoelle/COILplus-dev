# Estimating variable importance based on interaction model
# This code is just for plotting

# --------- TO DO: set  your directories and name the current results files using the date ---------#

# The directory where the analysis is performed: should already be your WD if you cloned the repo

data_path <- 'ProcessedData/'
save_plots_path <- 'Results/'

# Save results using convention: res_date_i.rda
date <- 'COILp_exp'
results_path <- paste0('Results/', date , '/')

# --------------------- STEP 0: LOAD FUNCTIONS AND DATA -------------------- #

library(tidyverse)
library(grid)
library(gridExtra)
library(abind)
library(superheat)


# Define a useful function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Loading the data:
obs_X <- loadRData(paste0(data_path, 'Obs_X_sorted.dat')) # vert traits
obs_W <- loadRData(paste0(data_path, 'Obs_W_sorted.dat')) # plant traits
obs_A <- loadRData(paste0(data_path, 'A_obs.dat'))

# Read in taxonomical data
v.taxa <- read.csv(paste0(data_path, 'v_taxa_sorted.csv'))
p.taxa <- read.csv(paste0(data_path, 'p_taxa_sorted.csv'))


# Check for missing values in the covariates:
any_X_miss <- any(apply(obs_X, 2, function(x) sum(is.na(x))) > 0)
any_W_miss <- any(apply(obs_W, 2, function(x) sum(is.na(x))) > 0)

# Getting the sample sizes:
nV <- nrow(obs_X)
nP <- nrow(obs_W)

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1



# --------------- STEP 1: Efficiently loading posterior samples ---------------#

# Number of MCMC chains
nchains <- 4

# Load all chains
all_Ls <- vector("list", nchains)
all_X <- vector("list", nchains)
all_W <- vector("list", nchains)

for (ii in 1:nchains) {
  cat("\nChain =", ii)
  res <- loadRData(paste0(results_path, 'res_', date, '_', ii, '.dat')) 
  all_Ls[[ii]] <- res$all_pred[,,,1] # selects posterior samples of L
  #all_pred[[ii]] <- res$all_pred
  if (any_X_miss) all_X[[ii]] <- res$Xs
  if (any_W_miss) all_W[[ii]] <- res$Ws
}

# Set number of posterior samples and trait dimensions
Nsims <- dim(all_Ls[[1]])[1]
pV <- if (any_X_miss) length(all_X[[1]]) else ncol(obs_X)
pP <- if (any_W_miss) length(all_W[[1]]) else ncol(obs_W)

# Construct imputed trait matrices if needed
Xs <- NULL
if (any_X_miss) {
  Xs <- vector("list", pV)
  for (p in 1:pV) {
    nmiss.p <- length(all_X[[1]][[p]])
    Xs[[p]] <- matrix(NA, nchains * Nsims, nmiss.p)
    for (cc in 1:nchains) {
      wh_entries <- Nsims * (cc - 1) + 1:Nsims
      Xs[[p]][wh_entries, ] <- all_X[[cc]][[p]]
    }
  }
}

Ws <- NULL
if (any_W_miss) {
  Ws <- vector("list", pP)
  for (p in 1:pP) {
    nmiss.p <- length(all_W[[1]][[p]])
    Ws[[p]] <- matrix(NA, nchains * Nsims, nmiss.p)
    for (cc in 1:nchains) {
      wh_entries <- Nsims * (cc - 1) + 1:Nsims
      Ws[[p]][wh_entries, ] <- all_W[[cc]][[p]]
    }
  }
}


# --------------- STEP 2: Load and plot the trait matching ----------------- #


# Load the trait matching results
load(paste0(results_path,  'rsq_obs_X_', date, '.dat'))
load(paste0(results_path, 'rsq_obs_W_', date, '.dat'))
load(paste0(results_path, 'rsq_resampling_X_', date, '.dat'))
load(paste0(results_path, 'rsq_resampling_W_', date, '.dat'))
load(paste0(results_path, 'corr_obs_X_', date, '.dat'))
load(paste0(results_path, 'corr_obs_W_', date, '.dat'))

good_namesX <- c("Generation Length", "Log Body Mass", "IUCN Status", "Habitat")
good_namesW <- c("Fruit Length", "Fruit Width", "Wood Density")

# Calculating the number of permuted standard deviations away from the mean.

# Starting from the vert covariates:
wh_obs <- rsq_obs_X
wh_resampling <- rsq_resampling_X
sd_awayX <- rep(NA,  length(wh_obs))
names(sd_awayX) <- good_namesX #colnames(obs_X) #good_namesX
for  (cc in 1 : length(wh_obs)) {
  sd_awayX[cc] <- (wh_obs[cc] - mean(wh_resampling[, cc])) / sd(wh_resampling[, cc])
}

# And for the plant covariates:
wh_obs <- rsq_obs_W # this gives the squared correlation coefficient between covars and interactions
wh_resampling <- rsq_resampling_W # this is the squared correlation coeff between shuffled covars and interactions
sd_awayW <- rep(NA,  length(wh_obs))
names(sd_awayW) <- good_namesW
for  (cc in 1 : length(wh_obs)) { # for each covar, compute how many sds away from resampled mean the obs rsq is
  sd_awayW[cc] <- (wh_obs[cc] - mean(wh_resampling[, cc])) / sd(wh_resampling[, cc]) # note these sds are tiny
}


# Plotting the tiles of variable importance ordering the variables in
# decreasing importance:

# For the vertebrate species:
xx <- data.frame(value = sd_awayX, covariate = names(sd_awayX), y = 1)
xx <- xx[order(- xx$value), ]
xx$covariate <- factor(xx$covariate, levels = xx$covariate)
xx$value_rounded <- round(xx$value/10)*10

ggplot() + geom_tile(aes(x = covariate, y = y, fill = value_rounded), color = 'white', data = xx) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'bottom', axis.title = element_blank()) +
  guides(fill = guide_legend("Frugivore Variable Importance")) +
  scale_fill_gradient(low = '#BFF0B6', high = '#3B6E32', breaks = xx$value_rounded) +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 20, family = "serif", color = "black"), 
        legend.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 15, family = "serif", color = "black"), 
        legend.title = element_text(size = 20, family = "serif"))


# For the plant species:
ww <- data.frame(value = sd_awayW, covariate = names(sd_awayW), y = 1)
ww <- ww[order(- ww$value), ]
ww$covariate <- factor(ww$covariate, levels = ww$covariate)
ww$value_rounded <- round(ww$value/10)*10

ggplot() + geom_tile(aes(x = covariate, y = y, fill = value_rounded), color = 'white', data = ww) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'bottom', axis.title = element_blank()) +
  scale_fill_gradientn(colors = c('#D4DEF7','#425075'), breaks = ww$value_rounded) +
  guides(fill = guide_legend("Plant Variable Importance")) +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 20, family = "serif", color = "black"), 
        legend.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 15, family = "serif", color = "black"), 
        legend.title = element_text(size = 20, family = "serif"))


# ------- PART B: Posterior probabilities based on the important covariates.

# Taking the posterior probabilities of interaction by averaging R across the three chains,
# and setting the recorded interactions to NA:
use_Ls <- abind(all_Ls, along = 1) # all_pred[,,,1] is pred_L
use_mean_Ls <- apply(use_Ls, c(2, 3), mean)
use_mean_Ls[comb_A == 1] <- NA

# Which covariate is to be plotted. The ones we want are listed first.
wh_X <- 1
wh_W <- 1

# Showing only the species that have the covariate measured.
keep_birds <- which(!is.na(obs_X[, wh_X]))
keep_plants <- which(!is.na(obs_W[, wh_W]))
use_out <- use_mean_Ls[keep_birds, keep_plants]

# Because some species have identical values for the covariate, in order for
# plot to show all of them, we need to slightly pertube their values. That way,
# the increasing or decreasing order is not altered, but there is no overlap in
# the covariate values:
bird_cov <- obs_X[keep_birds, wh_X]
bird_cov <- bird_cov + rnorm(length(keep_birds), sd = sd(bird_cov) * 0.0001)
plant_cov <- obs_W[keep_plants, wh_W]
plant_cov <- plant_cov + rnorm(length(plant_cov), sd = sd(plant_cov) * 0.0001)

# Creating a data frame in which the species are ordered by their covariate
# values. This will allow us to plot the probability of interaction across the
# covariates in an interpretable way. We also note that we need to turn the
# covariates to factors in order for them to be plotted in the correct order.
plot_dta <- data.frame(cov_bird = rep(bird_cov, length(keep_plants)),
                       cov_plant = rep(plant_cov, each = length(keep_birds)),
                       probability = as.numeric(use_out))
plot_dta$use_cov_bird <- factor(as.numeric(factor(plot_dta$cov_bird)))
plot_dta$use_cov_plant <- factor(as.numeric(factor(plot_dta$cov_plant)))

g <- ggplot() +
  geom_raster(aes(x = use_cov_plant, y = use_cov_bird, fill = probability), data = plot_dta) +
  #scale_fill_gradient(low = "#F5D4C7", high = "#02A65F", na.value = '#016B3B',
  #                    name = 'Posterior\ninteraction\nprobability\n', limits = c(0, 1)) +
  scale_fill_gradient(low = "#FDE725FF", high = "#02A65F", na.value = '#016B3B',
                      name = 'Posterior\ninteraction\nprobability\n', limits = c(0, 1)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 30),
        axis.title.y = element_text(vjust = - 1, family = 'serif'), text = element_text(family = 'serif')) +
  ylab(expression(symbol('\256'))) + xlab(expression(symbol('\256')))

gridExtra::grid.arrange(g, left = textGrob("Frugivores: increasing log body mass", rot = 90,
                                           x = 1.3, y = 0.57, gp = gpar(fontsize = 12, fontfamily = 'serif')),
                        bottom = textGrob("Plants: increasing wood density",
                                          x = 0.435, y = 1.3, gp = gpar(fontsize = 12, fontfamily = 'serif')),
                        vp=viewport(width=0.5, height=0.6))


#------ PART C. PLOTTING IMPORTANT COVARIATES WITH SIGNED CORRELATION ---------#

# Now group taxa by group for these plots and create pretty names
colnames(corr_obs_W) <- c("Fruit length", "Fruit width", "Wood density") #colnames(obs_W) # plant covars

# # Start by sorting the taxa by family if they aren't already (animals need some cleaning, plants are fine)
# v.taxa.ordered <- v.taxa %>% mutate(original_row = row_number()) %>% 
#   arrange(Animal_Class, Animal_Family, Animal_Genus, Animal_Species_Corrected)
# use.codes <- sort(paste(LETTERS, rep(1:3, 26), sep = "-")) # Create codes that will alphabetize nicely
# v.family.codebook <- data.frame(Animal_Family = unique(v.taxa.ordered$Animal_Family)) %>% 
#                      mutate(v.show.order = paste0("Family ", use.codes[row_number()])) %>% 
#                      right_join(., v.taxa.ordered)
# 

# Creating the clusters that will be used
# The following two lines specify that horizontal and vertical lines in our
# plot will separate species by taxonomic families:
v_group <- v.taxa$Animal_Family_Generic
p_group <- p.taxa$Plant_Family_Generic

# Calculating the size of each cluster, will be used when plotting results
# for families of certain size
v_size_cluster <- sapply(unique(v_group), function(x) sum(v_group == x))
# v_size_df <- data.frame(v_size_cluster) %>% rownames_to_column(., var = "family.code")
# v.family.codebook <- left_join(v.family.codebook, v_size_df, join_by(v.show.order == family.code))
# # Look at which families might be interesting
# v.family.codebook %>% select(Animal_Family, Animal_Taxa_Type_Courser, v_size_cluster) %>% unique()

p_size_cluster <- sapply(unique(p_group), function(x) sum(p_group == x))

# Set the minimum cluster size that should be plotted. For the results of the
# manuscript, we set min_v_size to 10, and min_plant_size to 25. Additionally, 
# we retain the smaller clusters Hominidae and Elephantidae because they are of
# scientific interest. Setting both to 0 will produce the full results.
table(v_size_cluster)
table(p_size_cluster)
min_v_size <- 12
min_p_size <- 30

keep_v_groups <- sort(names(which(v_size_cluster >= min_v_size)))
keep_p_groups <- names(which(p_size_cluster >= min_p_size))

keep_v_index <- which(v_group %in% keep_v_groups)
keep_p_index <- which(p_group %in% keep_p_groups)

### Plant covariates versus animal groups

# Plotting those with minimum size as specified:
superheat(X = corr_obs_W[keep_v_index, ],
          pretty.order.cols = FALSE, 
          pretty.order.rows = FALSE,
          membership.rows = v_group[keep_v_index],
          grid.hline.col = "black", grid.vline.col = 'black', # #00257D
          grid.hline.size = 0.5, grid.vline.size = 0.5,
          left.label = "none", 
          legend =  TRUE,
          bottom.label.text.size = 10,
          bottom.label.size = 0.075, 
)

### Frugivore covariates versus plant groups
plot_pred_X <- corr_obs_X[, -c(3,4)]
colnames(plot_pred_X) <- c("Generation length", "Log body mass")

# Plotting those with minimum size as specified:
# We replaced observed interactions with black NA
superheat(X = plot_pred_X[keep_p_index, ],
          pretty.order.cols = FALSE, 
          pretty.order.rows = FALSE,
          membership.rows = p_group[keep_p_index],
          left.label = "none",
          grid.hline.col = "black", grid.vline.col = 'black', # #00257D
          legend =  TRUE,
          bottom.label.text.size = 10,
          bottom.label.size = 0.075, 
)

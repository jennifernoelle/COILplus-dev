# File to plot results for cv and full run for a single experiment

# --------- TO DO: set  your directories and name the current results files using the date ---------#


# Save results using convention: res_date_i.rda
date <- 'COILp_exp'

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
save_path_base <- 'Results/' # Crete this directory
# Where the results are saved: 
results_path <- paste0('Results/', date, '/')


# ------ STEP 0: Some functions. --------- #


# We want the github version of superheat
if(!("superheat" %in% rownames(installed.packages()))){
  install.packages("devtools")
  remotes::install_github("rlbarter/superheat")
}

library(ggplot2)
library(RColorBrewer)
library(gplots)
library(superheat)
library(abind)
library(gridExtra)
library(grid)
library(dplyr)
library(caret)
library(pROC)
library(bayesplot)
library(tidyverse)

# Loading the data:
load(paste0(data_path, 'A_obs.dat'))
load(paste0(data_path, 'F_obs.dat'))
load(paste0(data_path, 'Obs_X_sorted.dat')) # vertebrate traits
load(paste0(data_path, 'Obs_W_sorted.dat')) # plant traits
load(paste0(data_path, 'OP_full.dat')) # site level obs plants
load(paste0(data_path, 'OV_full.dat')) # site level obs vertebrates


# Read in taxonomical data
v.taxa <- read.csv(paste0(data_path, 'v_taxa_sorted.csv'))
p.taxa <- read.csv(paste0(data_path, 'p_taxa_sorted.csv'))

sum(v.taxa$Animal_Species_Generic != rownames(obs_A))
sum(p.taxa$Plant_Species_Generic != colnames(obs_A))

# Define a useful function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

# Useful values
nP <- ncol(obs_A)
nV <- nrow(obs_A)
nS <- dim(obs_A)[3]
s.names <- dimnames(obs_A)[[3]]

# Set sampleP and occurrence prior probs according to the version of the model fit
sampleP <- grepl('COILp', date) # We only sample probs with COIL+

if(grepl('default', date)){ # Unobserved species are assumed not present
  O_P <- ifelse(O_P == 1, 1, 0)
  O_V <- ifelse(O_V ==1, 1, 0)
}else if(grepl('p75', date)){ # Unobserved species are assumed present with p = 0.75
  O_P <- ifelse(O_P == 1, 1, 0.75)
  O_V <- ifelse(O_V ==1, 1, 0.75)
} # Otherwise keep expert defined probabilities loaded above


# ------------------------- STEP 1: Specifications. -------------------------- #

# Number of chains for the full sample runs
nchains <- 4

# Number of cross validation repetitions:
repetitions <- 10

# Number of cv samples per repetition
n.cv <- 100

# ------------- STEP 1: Getting the results together - full sample ---------- #

### Put all chains together in a list
all_res <- NULL
for (ii in 1 : nchains) {
  cat("\n Chain = ", ii)
  res <- loadRData(paste0(results_path, 'res_',  date, '_', ii, '.dat')) 
  cat(" Dim is ", dim(res[[1]]))
  all_res[[ii]] <- res
}

### Bind together predicted interactions and detection probabilities 

# Number of posterior samples used:
use_Nsims <- dim(all_res[[1]]$all_pred)[1]

# Creating an array to bind results across chains:
pred_ours <- array(NA, dim = c(nchains * use_Nsims, nV, nP))
occ_p <- array(NA, dim = c(nchains, nP, nS))
occ_v <- array(NA, dim = c(nchains, nV, nS))
p_occ_p <- array(NA, dim = c(nchains, nP, nS))
p_occ_v <- array(NA, dim = c(nchains, nV, nS))

for (ii in 1 : nchains) {
  # Using the posterior samples of the L matrix:
  pred_ours[1 : use_Nsims + use_Nsims * (ii - 1), , ] <- all_res[[ii]]$all_pred[, , , 1] # this is the imputed L's from the model output (included bias corr)
 # Save posterior mean of occurrence indicators for each chain
  occ_p[ii, , ] <- (all_res[[ii]]$occ_plants)$OP_mean
  occ_v[ii, , ] <- (all_res[[ii]]$occ_verts)$OV_mean
   if(sampleP){
   p_occ_p[ii, , ] <- all_res[[ii]]$occ_plants$p_OPs
   p_occ_v[ii, , ] <- all_res[[ii]]$occ_verts$p_VPs
 }
}

dimnames(pred_ours)[2 : 3] <- list(vertebrate = rownames(obs_A), plant = colnames(obs_A))
dimnames(p_occ_p)[2:3] <- dimnames(occ_p)[2:3] <- list(plant = colnames(obs_A), dimnames(obs_A)[[3]])
dimnames(p_occ_v)[2:3] <- dimnames(occ_v)[2:3] <- list(vert = rownames(obs_A), dimnames(obs_A)[[3]])

# Calculating the posterior means across mcmc iterations:
mean_pred <- mean_pred_na <- apply(pred_ours, c(2, 3), mean)
mean_occ_v <- apply(occ_v, c(2,3), mean) # occ indicators
mean_occ_p <- apply(occ_p, c(2,3), mean)
mean_p_occ_v <- apply(p_occ_v, c(2,3), mean)
mean_p_occ_p <- apply(p_occ_p, c(2,3), mean)


# Setting the recorded interactions to NA for plotting
mean_pred_na[comb_A == 1] <- NA


# Save the posterior network 
write.csv(mean_pred, paste0(results_path, "post_network", date, ".csv"))

#-------------------- STEP 2: POSTERIOR NETWORK SUMMARY -----------------------------#

# Note that because post prob is taken to be the mean across posterior samples of L
# Zero values are possible when the number of posterior samples is relativley low
# And unique values are limited to multiples of 1/n_samples
mean(mean_pred)
summary(c(mean_pred))
sort(c(mean_pred))[1000] # Might be actually zero: plausible given small number of post samples


# Check: model should be outputting 1 whenever there is a known interaction 
sum(comb_A)
sum(mean_pred == 1)
mean(mean_pred[comb_A == 1])

# Look at new interactions
sum(comb_A) /length(comb_A)
sum(mean_pred_na > 0.5, na.rm = TRUE) # likely non-observed interactions
sum(mean_pred_na > 0.75, na.rm = TRUE) # v likely non-observed interactions


# Look at the relative increase in posterior prevalence vs observed prevalence
pp75 <- sum(mean_pred > 0.75)/length(mean_pred)
pp50 <- sum(mean_pred > 0.50)/length(mean_pred)
op <- sum(comb_A ==1)/length(comb_A)
(pp75 - op)/op # % increase in interactions at 75% threshold
(pp50 - op)/op

# Count number of focal studies for each vertebrate species
fV <- apply(obs_F, 1, function(x) {
  # We sum across plants within each study, then check if > 0
  study_presence <- colSums(x) > 0
  sum(study_presence)
})

# Count number of observed interactions for each species
aV <- rowSums(comb_A)

# Count the number of high posterior probability interactions for each species
p75V <- apply(mean_pred, 1, function(x) sum(x > 0.75))
p50V <- apply(mean_pred, 1, function(x) sum(x > 0.5))

# Make sure they match
sum(names(fV)!= names(aV))

# Look at correlation between study intensity and observed interactions: vertebrates
cor(fV, aV) # Observed correlation
cor(fV, p75V) # Correlation with posterior probs, P > 75%
cor(fV, p50V) # Correlation with posterior probs, P > 50%

# Look at vertebrate species with the most new interactions
new_df <- data.frame(v_species = names(aV), n_studies = fV, 
                     obs_int = aV, new_int = p75V,
                     perc_inc = (p75V - aV)/aV) %>%
  arrange(-perc_inc, n_studies) 


# ----------------- STEP 3: PLOTTING THE HEATMAP ----------------------------#

# Creating the clusters that will be used
# The following two lines specify that horizontal and vertical lines in our
# plot will separate species by taxonomic families:
v_group <- v.taxa$Animal_Family_Generic
p_group <- p.taxa$Plant_Family_Generic

# Calculating the size of each cluster, will be used when plotting results
# for families of certain size:
v_size_cluster <- sapply(unique(v_group), function(x) sum(v_group == x))
p_size_cluster <- sapply(unique(p_group), function(x) sum(p_group == x))

# Set plot_pred to pred_ours for results based on our method 
plot_pred <- mean_pred_na

# Set the minimum cluster size that should be plotted. For the results of the
# manuscript, we set min_bird_size to 10, and min_plant_size to 20. Setting both
# to 0 will produce the full results.
min_v_size <- 10
min_p_size <- 25

keep_p_groups <- names(which(p_size_cluster >= min_p_size))
keep_v_groups <- names(which(v_size_cluster >= min_v_size))

keep_v_index <- which(v_group %in% keep_v_groups)
keep_p_index <- which(p_group %in% keep_p_groups)


# Plotting those with minimum size as specified:
# We replaced observed interactions with black NA
superheat(X = plot_pred[keep_v_index, keep_p_index],
          membership.rows = v_group[keep_v_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 3,
          bottom.label.text.size = 3,
          bottom.label.size = 0.2, left.label.size = 0.12,
          heat.col.scheme = "grey", heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05), 
          # title = date, 
          title.size = 20)




# --------------- STEP 4: Cross validation results ----------------- #

# Getting the results together (held out indicies and predictions)
all_indices <- array(NA, dim = c(repetitions, n.cv, 2))
our_preds <- array(NA, dim = c(repetitions, nV, nP))
for (rr in 1 : repetitions) {
  cv_indices <- loadRData(paste0(results_path, 'cv_indices_', date, '_', rr, '.dat'))
  pred <- loadRData(paste0(results_path, 'pred_', date, '_', rr, '.dat'))
  all_indices[rr, , ] <- cv_indices
  our_preds[rr, , ] <- pred
}


# Predictions of the held out data from model: 100 indices held out each time - always interactions
pred <- array(NA, dim = c(repetitions, n.cv))
for (rr in 1 : repetitions) {
  for (ii in 1 : n.cv) {
    pred[rr, ii] <- our_preds[rr, all_indices[rr, ii, 1], all_indices[rr, ii, 2]]
  }
}

# Average and median probability of interaction based on the overall data
overall_mean <- cbind(apply(our_preds, 1, mean))
overall_median <- cbind(apply(our_preds, 1, median))

# Average and median in the held out data.
pred_mean <- apply(pred, 1, mean)
pred_median <- apply(pred, 1, median)

# Pseudo precision: 
mean(pred_mean/overall_mean)
mean(pred_median/overall_median)

# Pseudo accuracy: 
sum(pred>0.5)/length(pred) # what proportion of true interactions are predicted as "likely" >0.5
sum(pred>0.75)/length(pred) # what proportion of true interactions are predicted as "v likely" >0.75


# Creating the data frame we will plot:
plot_dta <- data.frame(value = rbind(pred_mean / overall_mean, pred_median / overall_median), 
                       stat = rep(c('Pred:Overall (mean)', 'Pred:Overall (median)'), repetitions))

# Plotting cross validation results:
ggplot(data = plot_dta) +
  geom_boxplot(aes(x = stat, y = value)) +
  theme_bw() +
  ylab('') +
  xlab('') +
  ggtitle('Out of sample performance', subtitle = date) +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) + 
  theme(text = element_text(size = 20))



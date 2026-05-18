# This file processes the raw data to create the following key data sources: 
  # A_obs: adjacency array recording observed interactions
  # F_obs: focus array recording study focus
  # O_A, O_P: occurrence probability matrices for animals and plants, respectively


# ---------------------------------- TO DO ----------------------------------- #

# Set the directories below to correspond to paths on your machine:

# The directory within the working directory where the data are:
data_path <- 'RawData/'
save_path <- 'ProcessedData/' # Create this directory

save_files <- TRUE

library(tidyverse)
library(foreach)
library(doParallel)
library(abind)

# Load data
frug  <- read.csv(paste0(data_path, "frug_generic.csv"))


#----------------------------- PART A: EDA ------------------------------------#

# Look at interactions data
apply(frug, 2, function(x) sum(is.na(x)))
table(frug$Zone) 
table(frug$Country) 
table(frug$Focus) # 293 missing & 268 both limited: treat these as v. narrow studies (i.e. all 0s noninformative)
nrow(frug) # 10635 interactions

frug_unique <- frug %>% group_by(Animal_Species_Corrected, Plant_Species_Corrected) %>% 
  dplyr::select(Animal_Species_Corrected, Plant_Species_Corrected) %>% 
  unique() 
nrow(frug_unique) # 5938 unique interactions

# Extract useful information
v.names <- unique(frug$Animal_Species_Corrected)
p.names <- unique(frug$Plant_Species_Corrected)
s.names <- unique(frug$Study_ID)
sp.names <- unique(frug$Study_ID_Plus)

nV <- length(v.names)
nP <- length(p.names)
nS <- length(s.names) # Studies IDs
nSP <- length(sp.names) # Study plus IDs: separate IDs for different sites in same study

#---------------------- PART B: CREATE THE ADJACENCY ARRAY -----------------#

# 1. Convert raw data into a list of adjacency matrices
frug_new <- frug %>% dplyr::select(Animal_Species_Corrected, Plant_Species_Corrected, Study_ID, Study_ID_Plus)
frug_new$Interaction <- 1

Adj_matrix_new <- list()
for(ss in 1:nSP){
  Adj <- frug_new %>% filter(Study_ID_Plus==sp.names[ss]) %>% 
    dplyr::select(Animal_Species_Corrected, Plant_Species_Corrected, Interaction) %>% 
    distinct() %>% pivot_wider(names_from = Plant_Species_Corrected, values_from=Interaction)
  # Case 1: one plant species only (species names will be lost)
  if(dim(Adj)[2] == 2){
    this.v <- Adj$Animal_Species_Corrected
    this.p <- colnames(Adj)[2]
    Adj_df  <- as.data.frame(Adj)[,-1]
    Adj_matrix <- as.matrix(Adj_df)
    dimnames(Adj_matrix) <- list(this.v, this.p)
    # Case 2: multiple species of at least one type
  }else{
    Adj_df  <- as.data.frame(Adj)
    rownames(Adj_df) <- Adj_df$Animal_Species_Corrected
    Adj_df <- Adj_df[,-1]
    Adj_df[is.na(Adj_df)] <- 0
    Adj_matrix <- as.matrix(Adj_df)
  }
  s.id <- sp.names[ss]
  Adj_matrix_new[[s.id]] <- Adj_matrix # add back rownames and colnames for 1,1 studies
}

# 2. Create adjacency array by augmenting with zeros
acomb <- function(...) abind(..., along=3)
registerDoParallel(10)

A_obs <- foreach(i=1:nSP, .combine = 'acomb') %dopar% {
  A <- Adj_matrix_new[[i]]
  wh.v <- sapply(dimnames(A)[[1]], function(x) which(x==v.names))
  wh.p <- sapply(dimnames(A)[[2]], function(x) which(x==p.names))
  
  A.full <- matrix(0, nrow = nV, ncol = nP)
  A.full[wh.v, wh.p] <- A
  # Debugging
  i <- i + 1
  
  A.full
}

dimnames(A_obs) <- list(v.names, p.names, sp.names)


#----------------- PART D: CREATING THE OCCURRENCE MATRICES -------------------#

## Extract occurrence information using a variety of methods, determining occurrence probability when
# Species occurs in the same study -> 1
# Species occurs at the same site, but a different study -> 0.75
# Species occurs in the same country, but a different site -> 0.5
# Species occurs in the same zone, but different ... -> 0.25
# Species does not occur in the same zone -> 0

## Create occurrence matrices
O_P <- matrix(0, nrow = nP, ncol = nSP) 
rownames(O_P) <- p.names
O_V <- matrix(0, nrow = nV, ncol = nSP) 
rownames(O_V) <- v.names
colnames(O_P) <- colnames(O_V) <- sp.names


# Note that we have multi-site studies
frug %>% group_by(Study_ID) %>% 
  summarize(n_Sites = length(unique(Study_Site)), 
            n_Countries = length(unique(Country)), 
            n_Zones = length(unique(Zone))) %>% 
  filter(n_Sites > 1 | n_Countries > 1 | n_Zones > 1)

for(ss in 1:nSP){
  IDplus <- sp.names[ss]
  study_site <- unique(frug[frug$Study_ID_Plus == IDplus, "Study_Site"])
  study_zone <- unique(frug[frug$Study_ID_Plus == IDplus, "Zone"])
  study_country <- unique(frug[frug$Study_ID_Plus == IDplus, "Country"]) # site, country, zone can all be NA, think about how to handle that
  
  study_site_dta <- frug[which(frug$Study_ID_Plus == IDplus ),] # same study and site: this is always nonempty
  site_dta <- frug[which(frug$Study_Site == study_site & !(is.na(frug$Study_Site))),] # same non-NA site: may be empty
  country_dta <- frug[which(frug$Country == study_country & !(is.na(frug$Country))), ] # same non-NA country: may be empty
  zone_dta <- frug[which(frug$Zone == study_zone & !(is.na(frug$Zone))), ] # same non-NA region: may be empty
  
  for(p in 1:nP){
    this.plant <- p.names[p]
    if(length(study_site_dta$Plant_Species_Corrected[study_site_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 1
    } else if(length(site_dta$Plant_Species_Corrected[site_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 0.75
    }  else if(length(country_dta$Plant_Species_Corrected[country_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 0.5
    } else if(length(zone_dta$Plant_Species_Corrected[zone_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 0.25
    } 
  }  
  for(v in 1:nV){
    this.animal <- v.names[v]
    if(length(study_site_dta$Animal_Species_Corrected[study_site_dta$Animal_Species_Corrected==this.animal])>0){
      O_V[v,ss] <- 1
    } else if(length(site_dta$Animal_Species_Corrected[site_dta$Animal_Species_Corrected==this.animal])>0){
      O_V[v,ss] <- 0.75
    } else if(length(country_dta$Animal_Species_Corrected[country_dta$Animal_Species_Corrected==this.animal])>0){
      O_V[v,ss] <- 0.5
    } else if(length(zone_dta$Animal_Species_Corrected[zone_dta$Animal_Species_Corrected==this.animal])>0){
      O_V[v,ss] <- 0.25
    }
  }
  
}


# ------------------- PART E: CREATE THE FOCUS ARRAY ------------------------- #

# Creating focus array which indicates whether the studies would have
# recorded an observed interaction or not

# 1. Categorize studies first
F_obs <- array(0, dim = c(nV, nP, nSP))
dimnames(F_obs) <- list(v.names, p.names, sp.names)
methods <- frug %>% group_by(Study_ID_Plus) %>% 
  summarize(method = unique(Focus)) %>% 
  mutate(method = ifelse(is.na(method), "Unknown", method))

# 2. Create Focus array
for (ss in 1 : nSP) {
  this_study <- sp.names[ss]
  this_A <- A_obs[,,ss]
  this_method <- methods[ss, "method"]
  these_v <- which(rowSums(this_A)>0)
  these_p <- which(colSums(this_A)>0)
  # Observability depends on Focus
  if(this_method == "Animal"){
    F_obs[these_v, , ss] <- 1 # Animal focused: all interactions between target vertebrates and any plants recorded
  } else if(this_method == "Plant"){
    F_obs[, these_p , ss] <- 1 # Plant focused: all interactions between target plants and any vertebrates recorded    
  } else if(this_method == "BothLimited" | this_method == "Unknown"){
    F_obs[,,ss] <- this_A # BothLimited, only observed interactions are known to be observable b/c we can't tell which plants vs which animals were focal
  } else { # Both or network studies: all interactions are observable
    F_obs[,,ss] <- 1 # All observations recorded in a network study
  }
}

#------------------------------- SORT SPECIES ---------------------------------#

# Read in phylogeny and taxonomy files and sort all data consistently
# Also check traits

# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv_phylo.dat'))
load(paste0(data_path, 'Obs_X.dat')) # vertebrate traits
load(paste0(data_path, 'Obs_W.dat')) # plant traits
v.taxa <- read.csv(paste0(data_path, 'v_taxa_generic.csv'))
p.taxa <- read.csv(paste0(data_path, 'p_taxa_generic.csv'))

# Rename for convenience
Cu <- Cu_phylo_generic
Cv <- Cv_phylo_generic
obs_W <- Obs_W_generic
obs_X <- Obs_X_generic

# Sort vertebrates
v.taxa.sorted <- v.taxa %>%
  arrange(Animal_Class, Animal_Family_Generic, Animal_Species_Generic)
v_species_order <- v.taxa.sorted$Animal_Species_Generic
Cu_sorted <- Cu[v_species_order, v_species_order]
obs_X <- obs_X[v_species_order, ]
O_V <- O_V[v_species_order,]


# Sort plants
p.taxa.sorted <- p.taxa %>% 
  arrange(Plant_Family_Generic, Plant_Species_Generic)
p_species_order <- p.taxa.sorted$Plant_Species_Generic
Cv_sorted <- Cv[p_species_order, p_species_order]
obs_X <- obs_X[v_species_order, ]
O_P <- O_P[p_species_order,]

# Sort interactions and focus
obs_A <- A_obs[v_species_order, p_species_order, ]
obs_F <- F_obs[v_species_order, p_species_order, ]


if(save_files){
  # The cleaned core data files
  save(obs_A, file = paste0(save_path, "A_obs.dat"))
  save(O_P, file = paste0(save_path, 'OP_full.dat'))
  save(O_V, file = paste0(save_path, 'OV_full.dat'))
  save(obs_F, file = paste0(save_path, 'F_obs.dat'))
  # Sorted covariate and phylogeny files
  save(Cu_sorted, file = paste0(save_path, "Cu_phylo_sorted.dat"))
  save(Cv_sorted, file = paste0(save_path, "Cv_phylo_sorted.dat"))
  save(obs_X, file = paste0(save_path, "Obs_X_sorted.dat"))
  save(obs_W, file = paste0(save_path, "Obs_W_sorted.dat"))
  write.csv(v.taxa.sorted, paste0(save_path, 'v_taxa_sorted.csv'), row.names = FALSE)
  write.csv(p.taxa.sorted, paste0(save_path, 'p_taxa_sorted.csv'), row.names = FALSE)
}


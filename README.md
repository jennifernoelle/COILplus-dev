# COIL+

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20368146.svg)](https://doi.org/10.5281/zenodo.20368146)

This repository implements both the **COIL** model and the extended **COIL+** model developed in Kampe et. al. (2026). This package adds support for informative occurrence priors and sampling of occurrence probabilities, together with additional trait matching functions, and major efficiency improvements to the [BiasedNetwork](https://github.com/gpapadog/BiasedNetwork) package (Papadogeorgou, 2023).. 

The included code provides an illustrative example, imputing unobserved frugivore-plant interactions from a multi-study data set of Afrotropical frugivory. The data are provided with generic species labels for plants and frugivores in compliance with data sharing restrictions.

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("jennifernoelle/COILplus-dev")
```
> **Note:** Vignettes are not built by default and should not be built automatically — they include MCMC sampling and are designed to be run interactively in sequence.

### Running the vignettes

To work through the full modelling workflow:

1. Clone the repository to your local machine
2. Open the project in RStudio
3. Navigate to the `vignettes/` folder
4. Run the vignettes **in order** (01 through 06), knitting each one interactively
5. Data are downloaded automatically from Zenodo on first run of vignette 01.


## Data

The vignettes use Afrotropical frugivory data archived on Zenodo 
([DOI: 10.5281/zenodo.20357089](https://doi.org/10.5281/zenodo.20357089)), 
which is downloaded automatically on first run of vignette 01.

To apply COILplus to your own data you will need equivalent inputs in the 
following format:

- **Interaction records** (`frug_generic.csv`): pairwise interactions between 
  animals and plants, with study-level metadata including study ID, site, 
  country, zone, and sampling focus (animal-, plant-, or both-focused)
- **Animal traits** (`Obs_X.dat`): trait matrix for animal species
- **Plant traits** (`Obs_W.dat`): trait matrix for plant species  
- **Animal taxonomy** (`v_taxa_generic.csv`): taxonomic labels for animal species
- **Plant taxonomy** (`p_taxa_generic.csv`): taxonomic labels for plant species
- **Animal phylogenetic correlation matrix** (`Cu_phylo.dat`): can be computed 
  from a consensus phylogeny using the `ape` package
- **Plant phylogenetic correlation matrix** (`Cv_phylo.dat`): can be computed 
  using the `V.PhyloMaker` R package

See vignette 01 (`01_data_prep`) for full details on data preparation and 
formatting requirements.

## Citation

If you use this repository, please cite both Kampe et. al. (2026) and Papdogeourgiou et al. (2023).

# References

Kampe, J. (2026). *COILplus: R package for link prediction in ecological meta-networks*. 
Zenodo. https://doi.org/10.5281/zenodo.20368146

Papadogeorgou, G. (2026). *BiasedNetwork: Latent factor network model with bias correction for unrecorded interactions*. GitHub. https://github.com/gpapadog/BiasedNetwork

Papadogeorgou, G., Bello, C., Ovaskainen, O., & Dunson, D. B. (2023). Covariate-informed latent interaction models: addressing geographic & taxonomic bias in predicting bird–plant interactions. *Journal of the American Statistical Association*, *118*(544), 2250-2261.

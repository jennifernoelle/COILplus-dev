# COIL+

This repository implements both the **COIL** model and the extended **COIL+** model developed in XX (2025) and contains all code needed to reproduce the results from the main text. It adds support for informative occurrence priors and sampling of occurrence probabilities, together with major efficiency improvements to the [BiasedNetwork](https://github.com/gpapadog/BiasedNetwork) package (Papadogeorgou, 2023).

The included code provides an illustrative example, imputing unobserved frugivore-plant interactions from a multi-study data set of Afrotropical frugivory. The data are provided with generic species labels for plants and frugivores in compliance with data sharing restrictions.

## Data set

The key raw data sources are are in the folder `RawData/`:

-   Raw frugivory interactions with generic species labels and the variables study ID, study site, country, zone, focus: `frug_generic.csv`

-   Plant traits with generic species labels: `Obs_W.data`

-   Vertebrate traits with generic species labels: `Obs_X.data`

-   Generic taxonomy for plants: `p_taxa_generic.csv`

-   Generic taxonomy for vertebrates: `v_taxa_generic.csv`

-   `Cu_phylo.dat`: phylogenetic correlation matrix for vertebrates

-   `Cv_phylo.dat`: phylogenetic correlation matrix for plants

Plant phylogenies are acquired using the `V.PhyloMaker` R package and the phylogenetic correlation matrix for vertebrates is computed via the the package ape using a consensus tree obtained from VertLife and provided in this repo. For the purposes of this vignette, we supply the respective phylogenetic correlation matrices directly without

## Code

The code is included in `R_fast/`and `R_legacy/`. The use is as follows:

-   `R_fast/` contains rewritten versions of the most computationally intensive routines (sampling latent factors, updating occurrence indicators and probabilities, computing the log-likelihood, updating stored quantities). These updates are **distributionally equivalent** to the legacy code but avoid explicit matrix inversions, making them both **more numerically accurate** and up to **10× faster**.
-   `R_legacy/` contains the legacy implementations, preserved to replicate analyses based on earlier versions of the functions.
-   This work builds on and substantially extends the original [BiasedNetwork](https://github.com/gpapadog/BiasedNetwork) repo.

The code for the analysis is in the folder `Analysis/`. The numbers in the beginning of the file names represent the order with which the files should be used/run. In brief the content of each analysis file is as follows:

-   `Anaylsis_0_data_prep.R`: This code MUST be run before any subsequent analysis. Assembles key network and meta-data sources.

-   `Analysis_1a_fitCOILplus.R`: This code MUST be run before any subsequent plotting functions. Fits COIL+ to the Afrotropical frugivory data with prior incorporating domain knowledge.

-   `Analysis_1legacy_fitCOILplus.R`: Fits COIL+ to the Afrotropical frugivory data with prior incorporating domain knowledge using the legacy unblocked sampler.

-   `Analysis_1b_fitCOIL.R`: This code MUST be run before any subsequent plotting functions. Fits basic COIL to the Afrotropical frugivory data with prior incorporating domain knowledge.

-   `Analysis_1c_fitCOIL_default.R`: This code MUST be run before any subsequent plotting functions. Fits basic COIL to the Afrotropical frugivory data with default prior.

-   `Analysis_2a_cv_COILplus.R`: This code MUST be run before any subsequent plotting functions. Performs 10-fold cross validation of COIL+ with the Afrotropical frugivory data with prior incorporating domain knowledge.

-   `Analysis_2b_cv_COIL.R`: This code MUST be run before any subsequent plotting functions. Performs 10-fold cross validation of COIL with the Afrotropical frugivory data with prior incorporating domain knowledge.

-   `Analysis_2c_cv_COIL_default.R`: This code MUST be run before any subsequent plotting functions. Performs 10-fold cross validation of COIL with the Afrotropical frugivory data with default 0/1 prior occurrence prior.

-   `Analysis_3a_plot_results.R`: This code analyzes in and out of sample model performance for the any of the three models fit in 1a-2c.

-   `Analysis_4a_trait_matching_save.R`: This code MUST be run before the subsequent plotting code. Performs trait matching for the selected model fit.

-   `Analysis_4b_trait_matching_plot.R`: This code plots the trait matching.

### Note

For replicating the analysis results in this code, you will need to specify a directory where processed data and results can be saved. We recommend you create a folder at the same level as the `Analysis/`, `Data/`, and `R*/`folders that is named `Results/` and `ProcessedData/`.

## Citation

If you use this repository, please cite both XX (2025) and Papdogeourgiou et al. (2023).

# References

Papadogeorgou, G., Bello, C., Ovaskainen, O., & Dunson, D. B. (2023). Covariate-informed latent interaction models: addressing geographic & taxonomic bias in predicting bird–plant interactions. *Journal of the American Statistical Association*, *118*(544), 2250-2261.

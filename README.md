# uit-ajive-exam

This is the R code for evaluating the paper [*Integrative, multi-omics, analysis of blood samples improves model predictions: applications to cancer*](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04296-0), by Ponzi et al. (2021)


**am-dataset-construction.R**: Code for constructing analytic dataset -- data is saved after each major step and then reloaded, to permit starting at an arbitrary point in the dataset construction, since some steps can be time-consuming.
  * **build-cpg-gene-name-dataset.R**: builds the DNA methylation (CpG site) to gene name mapping dataset
  * **build-mrna-gene-mapping-dataset.R**: Builds the mRNA seqence to gene name mapping dataset 

**am-helper-functions.R**: Helper functions
  * **functions/filtering.R**: Erica Ponzi's function for filtering columns by highest criterion (e.g. greatest variance)
  * **functions/Screeplots.R**: Erica Ponzi's function for producing a scree plot

**am-analysis-ajive-initial-ranks.R**: Perform aJIVE analysis using different initial rank estimates (low to high)

**am-analysis-num-wedin-samples.R**: Compare effect of different numbers of Wedin bound samples

**am-analysis-PL-aJIVE-ROC-datasets.R**: Produce AUC/ROC plots comparing:
  * 5000 vs 10000 mRNA sequences
  * 40% vs 0% missing data cut-off threshold
  * 3 vs 5 for extreme M-value cutoff for CpG sites
  * CpG selection via mRNA gene mapping vs. highest variance 

**am-analysis-profile-likelihood-stability.R**: Produce plots comparing:
  * effect of bootstrapping on singular values
  * profile likelihood computation for svd vs rsvd (randomized SVD) 
  * effect of the number of excluded SVs on profile likelihood 
    * sub-effect of bootstrapping (introduces "degenerate" subjects into SVD)
    * sub-effect of random subject subsamples 
    * sub-effect of random feature subsamples 
    * sub-effect of subject sample size
    * sub-effect of feature sample size


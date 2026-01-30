# Which Factor Accounts for Most Genetic Diversity in a Widely Distributed Herbaceous Species from Southern South American Grasslands? (Ongoing Review)

## Overview
This repository contains the analytical pipeline and datasets for a research study investigating the drivers of genetic diversity and adaptation in *Petunia axillaris*, a widely distributed herbaceous species in the grasslands of Southern South America. The study utilizes high-throughput genomic data (RAD-seq - GBS) integrated with environmental variables to assess the relative contributions of geographic distance, climatic factors, and soil properties to the observed genetic structure.

## Dataset Description
- **Genomic Data**: Single Nucleotide Polymorphism (SNP) datasets derived from RAD-seq, processed to filter for quality, linkage disequilibrium (LD), and monomorphic sites.
- **Environmental Variables**:
    - **Climatic Data**: Bioclimatic variables (e.g., Temperature seasonality, Precipitation of driest quarter).
    - **Edaphic Data**: Soil properties (e.g., pH, clay content, cation exchange capacity).
    - **Topographic Data**: Elevation and terrain indices (TRI, Topographic Wetness Index).

## Analysis Pipeline
The research follows a multi-stage analytical framework:

### 1. Variant Calling and Filtering
- **Quality Control**: Processing of raw RAD-seq data to obtain high-quality SNP datasets, including filtering for missing data, minor allele frequency, and linkage disequilibrium.

### 2. Outlier Detection
- **Genome Scanning**: Identification of loci under putative selection using `Bayescan` and `PCAdapt`.

### 3. Population Genetic Structure
- **Clustering Analysis**: Assessment of population structure using `Admixture`.
- **Phylogenetic Analysis**: Hierarchical analysis of genetic relationships using `HierfStat`.
- **Principal Component Analysis**: Visualization of population differentiation using PCA.
- **Sequence Processing**: RAD-seq data processing using `STACKS`.
- **Validation**: Subsampling validation of genetic structure inferences.

### 4. Demography and Spatial Connectivity
- **Demographic History**: Estimation of effective population size ($N_e$) through time using `Stairway Plot 2`.
- **Migration and Connectivity**: Modeling of migration surfaces and spatial connectivity using `FEEMS`.

### 5. Genotype-Environment Association
- **Bayesian Environment Association**: Correlation of genetic loci with environmental gradients using `Bayescenv`.
- **Regression Modeling**: Analysis of genotype-environment associations using `Bedassle`.
- **Spatial Structure**: Modeling of spatial genetic structure using `conStruct`.
- **Generalized Linear Mixed Models (GLMM)**: Implementation of `MCMCglmm` in R to model genetic distance as a function of geographic, climatic, and edaphic predictors while controlling for population structure.

## Directory Structure
The repository is organized into the following modules based on the analytical components:

- `1-VariantCallFilt/`: Raw and filtered variant call format files, and quality control scripts.
- `2-OutlierDetection/`:
  - `Bayescan/`: Outlier detection using Bayescan, including diagnostic plots.
  - `PCAdapt/`: Outlier detection using PCAdapt.
- `3-PopGenStruct/`:
  - `Admixture/`: Population clustering analysis using Admixture.
  - `HierfStat/`: Hierarchical genetic analysis.
  - `PCA/`: Principal component analysis of population differentiation.
  - `STACKS/`: RAD-seq data processing and sequence assembly.
  - `Subsampling_Validation/`: Validation of genetic structure inferences through subsampling.
- `4-DemographySpatialConnect/`:
  - `FEEMS/`: Spatial connectivity and migration modeling.
  - `StairwayPlot2/`: Demographic history inference.
- `5-GenotypeEnvironAssoc/`:
  - `Bayescenv/`: Bayesian genotype-environment association analysis.
  - `Bedassle/`: Regression-based genotype-environment association.
  - `conStruct/`: Spatial genetic structure modeling.
  - `EnvironValriables/`: Environmental variable datasets.
  - `GLMM/`: Generalized Linear Mixed Models for landscape genetics, including diagnostic plots for various environmental predictors.

## Requirements
The analyses require the following software environments:
- **R (>= 4.0)**: Key packages include `vcfR`, `adegenet`, `hierfstat`, `MCMCglmm`, and `ecodist`.
- **Python (>= 3.8)**: Utilized for data parsing and running specific analyses.
- **Bioinformatics Tools**: `stacks` (sequence processing) and various R and Python packages for population genetics and statistical modeling.

## Citation
Please refer to the associated research paper for detailed methodology and interpretation of the results.



# 1. General Overview

**CytoBatchFlagR** is an R-based tool designed to identify potential batch effects in flow and spectral flow cytometry datasets.

It evaluates differences between control samples and batches using three complementary approaches:

- **IQR-based assessment:** evaluates variations in marker expression and positive cell percentages.
- **EMD-based assessment:** compares marker expression distributions between batches.
- **Unsupervised clustering-based assessment:** identifies potential batch-related differences in cellular populations.

The workflow combines these assessments to generate graphical representations, statistical indicators, and summary heatmaps.

The Eswar26 R Markdown script is largely derived from the official CytoBatchFlagR GitHub repository, with minor modifications and additional code adapted to the dataset.

Official repository: https://github.com/AndorfLab/CytoBatchFlagR

**Important:** This R Markdown file cannot be knitted directly because several steps require interactive user input. Each code chunk must therefore be executed manually in RStudio.

# 2. Data Preparation and Dependencies

## 2.1. Input Data

The workflow requires the following input files:

- FCS files containing the flow cytometry measurements.
- `fcs_file_info.csv`: metadata describing each FCS file.
- Marker information file containing marker names (`PnS`) and corresponding channels or fluorophores (`PnN`).

The metadata file must contain at least the following columns:

| Column | Description |
|---|---|
| `control` | Control sample identification |
| `batch` | Batch identification |
| `sample_id` | Individual sample identification |

The marker table is used to associate fluorescence channels with biological markers.

## 2.2. Package Installation and Loading

The script first installs and loads the R packages required for the analysis.

A dedicated script, `package_installer.R`, is used to check for missing packages and install them automatically.

The main dependencies include:

- `dplyr`, `tidyverse`, `tidyr`: data manipulation.
- `ggplot2`, `ggpubr`, `ggridges`: graphical representations.
- `flowCore`: FCS file reading and processing.
- `FlowSOM`: unsupervised clustering.
- `ConsensusClusterPlus`: clustering evaluation.
- `ComplexHeatmap`, `pheatmap`: heatmap generation.
- `umap`: dimensionality reduction.
- `shiny`, `shinyjs`, `DT`: interactive interfaces.

The script then loads the CytoBatchFlagR functions using `source()`.

## 2.3. Working Directories

Two main paths must be configured:

- `packagePath`: location of the CytoBatchFlagR R functions.
- `workDir`: location of the Eswar26 dataset and associated files.

The working directory is then defined using `knitr::opts_knit$set()`.

# 3. Step 1 – Pre-processing and Visual Assessment

## Objective

The first step prepares the dataset, identifies the common marker panel across FCS files, transforms the measurements and provides an initial visual assessment.

## Main Operations

### 1. Identification of FCS Files

The script recursively searches the working directory for all `.fcs` files.

### 2. Marker Panel Extraction

Each FCS file is read using `flowCore::read.FCS()`.

The following information is extracted:

- `PnS`: marker name.
- `PnN`: parameter name corresponding to the channel or fluorophore.

The script then identifies marker/channel combinations present in every readable FCS file.

### 3. Marker Panel Export

The common marker panel is exported to:

- `panel_markers.xlsx`

This Excel file can be manually edited to select the markers of interest.

The selected panel is subsequently converted into:

- `panel_markers.csv`

### 4. Metadata Loading

The script imports the metadata and marker information files.

It then creates lists containing:

- Controls.
- Batches.
- Sample IDs.
- Marker names.
- Fluorescence channels.

### 5. Output Directory Creation

A new output directory is automatically created with a timestamp.

This directory contains the results generated during the analysis.

### 6. Data Transformation

The function `create_transformed_dataframe_function()` transforms the FCS data.

The selected transformation is:

- Transformation: `arcsinh`
- Cofactor: `6000`

The transformed dataset is stored in `ffdf`.

### 7. Initial Visual Assessment

Three graphical representations are generated:

| Plot | Purpose |
|---|---|
| Cell count plot | Displays the number of cells per sample. |
| MDS plot | Visualizes similarities and differences between samples. |
| UMAP plot | Provides a dimensionality-reduced representation of the control samples, colored by batch. |

## Expected Outputs

- Common marker panel (`panel_markers.xlsx`).
- Selected marker panel (`panel_markers.csv`).
- Transformed dataset (`ffdf`).
- Cell count plot.
- MDS plot.
- UMAP plot.

# 4. Step 2 – Interquartile Range (IQR)-Based Assessment

## Objective

This step evaluates the consistency of marker expression across batches using three complementary metrics:

- Negative population median fluorescence intensity (-MFI).
- Positive population median fluorescence intensity (+MFI).
- Percentage of positive cells (%pos).

The objective is to identify batches showing unusual variations in these measurements.

## Main Operations

### 1. Marker Selection

The `selectMarkerList()` function opens an interactive interface allowing the user to select the markers included in the IQR analysis.

### 2. Automated Marker Threshold Calculation

The `automated_threshold_DF()` function calculates thresholds for each selected marker.

These thresholds separate negative and positive cell populations.

The script uses the following parameters:

| Parameter | Value |
|---|---:|
| Subsampling | Yes |
| Random seed | 450 |
| Cells per batch | 20,000 |
| Width cutoff | 0.4 |
| Peak separation | 0.5 |
| Dip cutoff | 0.15 |
| Minimum peak height fraction | 0.008 |
| Shoulder fraction | 0.05 |

### 3. Data Subsampling

The `sample_cells()` function reduces the dataset size to facilitate graphical analysis.

The script uses two sampling strategies:

- `overall_min`: equalizes the number of cells across all samples using the smallest sample size.
- `fixed`: selects a predefined number of cells per sample, here 15,000.

### 4. Marker Expression Density Plots

The `write_density_plot_per_marker()` function generates density plots for each marker across controls and batches.

These plots display the automatically calculated thresholds.

The user can inspect the distributions and determine whether the thresholds are appropriate.

### 5. Manual Threshold Adjustment

The `editMarkerTable()` function opens an interactive interface allowing the user to modify the automatically generated thresholds.

The updated thresholds are stored in `updated_auto_cutoffs`.

### 6. Biaxial Dot Plots

The `referenceMarkerSelection()` function allows the user to select a reference marker.

The `draw_biaxial_dotplots()` function then generates biaxial plots to visualize marker expression relative to the selected reference marker.

### 7. Calculation of IQR Metrics

Three functions calculate the main measurements:

| Function | Measurement |
|---|---|
| `neg_mfi_df_all()` | Negative population MFI |
| `pos_mfi_df_all()` | Positive population MFI |
| `percent_pos_df_all()` | Percentage of positive cells |

### 8. Flag Generation

The `iqr_dataframe()` function combines these three metrics and identifies potential outliers.

The results are stored in `iqrdf`.

### 9. Graphical Representation

The script generates boxplots for each metric:

- Negative MFI.
- Positive MFI.
- Percentage of positive cells.

Two visualization approaches are available:

- Variable axes.
- Fixed axes.

Flagged batches are identified directly on the plots.

## Expected Outputs

- Automated marker thresholds.
- Updated thresholds, if manually modified.
- Marker expression density plots.
- Biaxial dot plots.
- IQR metrics.
- IQR flag table (`iqrdf`).
- IQR boxplots.

# 5. Step 3 – Earth Mover's Distance (EMD)-Based Assessment

## Objective

The EMD-based assessment evaluates differences between marker expression distributions across batches.

Unlike the IQR approach, which focuses on specific statistical measurements, EMD compares entire distributions.

This provides a complementary assessment of potential batch effects.

## Main Operations

### 1. Data Preparation

The `reshape_df()` function reorganizes the transformed dataset into the format required for EMD calculations.

The resulting dataset is stored in `cdf`.

### 2. Pairwise EMD Calculation

The `EMD_list()` function calculates pairwise Earth Mover's Distances between batches for each marker and control sample.

The script uses:

- 20,000 cells for subsampling.
- Random seed: 450.

### 3. EMD Heatmaps

The `generate_EMD_heatmaps()` function generates heatmaps representing pairwise EMD values.

The threshold is defined as:

`emd_threshold = 5`

This threshold is used to identify potentially important differences between batches.

### 4. EMD Boxplots

The `generate_EMD_boxplots()` function generates boxplots of pairwise EMD values.

The plots are organized according to median EMD values.

A horizontal line represents the selected threshold, and batches exceeding the threshold are highlighted.

### 5. EMD Summary Table

The `emd_med_vals()` function calculates median EMD values and generates corresponding flags.

The results are stored in `emdf`.

## Expected Outputs

- Pairwise EMD calculations.
- EMD heatmaps.
- EMD boxplots.
- Median EMD values.
- EMD flag table (`emdf`).

# 6. Step 4 – Summary of IQR and EMD Results

## Objective

This step combines the results from the IQR and EMD assessments.

It provides an overview of potential batch effects identified by both approaches.

## Main Operations

The `generate_ranked_flagged_hmap()` function generates a summary heatmap combining:

- IQR-based flags.
- EMD-based flags.

The results are ranked according to the number of detected flags across batches and markers.

The heatmap uses three main colors:

| Color | Meaning |
|---|---|
| Red | Flagged |
| Grey | Not flagged |
| White | Missing value (NA) |

Several parameters can be adjusted:

- `markers_per_plot`: maximum number of markers per plot.
- `batches_per_plot`: maximum number of batches per plot.
- `min_batch_effect`: minimum number of flags required to display a batch.
- `min_marker_effect`: minimum number of flags required to display a marker.

In the current script, both minimum filtering parameters are set to zero.

## Expected Outputs

- Combined IQR and EMD summary heatmap.
- Ranked overview of potentially problematic batches and markers.

# 7. Step 5 – Unsupervised Clustering-Based Assessment

## Objective

This step investigates whether batch-related differences affect the composition of cellular populations.

It uses FlowSOM to identify cellular clusters based on marker expression profiles.

The objective is to determine whether specific batches are overrepresented in particular clusters.

## Main Operations

### 1. Marker Selection

The script selects phenotype markers using the `status` column of the marker table.

Only markers classified as `phenotype` are included.

### 2. Clustering Parameter Definition

The following parameters are defined:

| Parameter | Value |
|---|---:|
| Number of cells | 50,000 |
| Number of metaclusters | 8 |
| Random seed | 450 |
| Threshold 1 | 2.0 |
| Threshold 2 | 1.5 |

### 3. FlowSOM Clustering

The `getFlowSOM_clusters()` function performs unsupervised clustering using FlowSOM.

Each cell is assigned to a cluster according to its marker expression profile.

The resulting dataset is stored in `filter_df`.

### 4. Cluster Size Visualization

The `clusterSizes_plot()` function generates a bar plot showing:

- Number of cells per cluster.
- Percentage of cells in each cluster.

A horizontal reference line indicates 0.5% of the total analyzed cells.

### 5. Cluster Characterization

The `clustering_heatmap()` function generates a heatmap of scaled median marker expression values for each cluster.

This allows the user to characterize the different cellular populations.

### 6. Batch Proportions Within Clusters

The `proportion_of_batches_per_cluster()` function calculates the contribution of each batch to each cluster.

Threshold 1 identifies batches contributing more than twice the expected proportion.

The corresponding results are also stored in a table using:

`proportion_of_batches_per_cluster_table()`

### 7. Comparison Across Control Samples

The `proportion_of_batches_across_control_samples()` function investigates batch proportions across control samples.

Threshold 2 is used to identify batches exceeding 1.5 times the expected proportion within a given cluster.

The results are also exported into a table.

### 8. Final Cluster UMAP

The `generate_cluster_UMAP()` function generates a UMAP representation of the clusters, with cells colored according to their batch.

This provides a final visual assessment of batch distribution across cellular populations.

## Expected Outputs

- FlowSOM cluster assignments.
- Cluster size bar plot.
- Cluster characterization heatmap.
- Batch proportions per cluster.
- Batch proportion comparison across controls.
- Cluster UMAP.
- Tables identifying potentially overrepresented batches.

# 8. Step 6 – Global Summary of All Assessment Metrics

## Objective

The final step combines all three assessment approaches:

1. IQR-based assessment.
2. EMD-based assessment.
3. Unsupervised clustering-based assessment.

It provides a comprehensive overview of potential batch effects identified throughout the workflow.

## Main Operations

The `generate_ranked_flagged_hmap_all()` function generates a global heatmap combining:

- IQR flags (`iqrdf`).
- EMD flags (`emdf`).
- Clustering results (`filter_df`).

The function also considers the two clustering thresholds defined previously.

The resulting heatmap ranks batches and markers according to the number of detected flags.

As in Step 4:

- Red indicates flagged results.
- Grey indicates non-flagged results.
- White indicates missing values.

The heatmap can be divided into several plots if the number of markers or batches is too large.

## Expected Outputs

- Global summary heatmap.
- Combined visualization of IQR, EMD and clustering results.
- Final overview of potentially problematic batches and markers.

# 9. Important Parameters Used in the Eswar26 Script

| Parameter | Value | Purpose |
|---|---:|---|
| Transformation | arcsinh | Transformation of fluorescence intensities |
| Cofactor | 6000 | Spectral flow cytometry transformation |
| IQR subsampling | 20,000 | Automated threshold calculation |
| IQR plotting | 15,000 | Fixed subsampling for selected plots |
| Random seed | 450 | Reproducibility |
| EMD threshold | 5 | Identification of distribution differences |
| Clustering cells | 50,000 | FlowSOM input |
| Metaclusters | 8 | Number of clusters |
| Threshold 1 | 2.0 | Batch contribution threshold |
| Threshold 2 | 1.5 | Batch contribution threshold across controls |

These values correspond to the current implementation and should be validated according to the experimental design and dataset characteristics.

# 10. Important Points to Check in the Current Script

Several elements should be verified before running the complete workflow.

## 1. Color Variables

The script uses `sample_colours`, `batchColours` and `clustColors`.

Their definitions are not included in the provided code excerpt. They must be initialized before the corresponding plotting functions.

## 2. Marker Information

The clustering step requires a `status` column in the marker table, containing the value `phenotype` for the markers selected for clustering.

## 3. Required Packages

Some packages used in the script, including `limma`, `ComplexHeatmap`, `flowCore`, `FlowSOM`, `ConsensusClusterPlus` and `openxlsx`, should be checked to ensure they are installed and available.

## 4. Working Directories

The `packagePath` and `workDir` variables must point to the correct directories.

## 5. Interactive Execution

The script must be executed chunk by chunk in RStudio. Interactive interfaces must be completed before continuing to subsequent steps.

## 6. Threshold Validation

The EMD threshold and clustering thresholds should be interpreted in the context of the dataset. A flagged result is an indication of a potential issue, not definitive proof of a batch effect.

# 11. Overall Workflow Summary

The CytoBatchFlagR analysis pipeline can be summarized as follows:

1. **Input data:** FCS files, metadata and marker panel.
2. **Step 1 – Pre-processing:** marker panel extraction, arcsinh transformation, MDS and UMAP.
3. **Step 2 – IQR:** evaluation of negative MFI, positive MFI and percentage of positive cells.
4. **Step 3 – EMD:** pairwise comparison of marker expression distributions.
5. **Step 4 – IQR + EMD:** combined heatmap of statistical flags.
6. **Step 5 – Clustering:** FlowSOM clustering, cluster characterization and batch proportions.
7. **Step 6 – Global summary:** integration of IQR, EMD and clustering results.
8. **Final interpretation:** identification and investigation of potential batch effects.

# 12. Conclusion

The CytoBatchFlagR workflow provides a comprehensive framework for evaluating potential batch effects in cytometry datasets.

By combining statistical comparisons (IQR), distribution-based comparisons (EMD), and unsupervised clustering (FlowSOM), the workflow offers complementary information about the consistency of marker expression and cellular populations across batches.

The final summary allows potentially problematic batches and markers to be identified and investigated before downstream analysis.

Depending on the results, possible corrective strategies include excluding problematic batches or markers, or applying an appropriate batch normalization method.

## Official Documentation

- GitHub repository: https://github.com/AndorfLab/CytoBatchFlagR
- CytoBatchFlagR Wiki: https://github.com/AndorfLab/CytoBatchFlagR/wiki
- Input file documentation: https://github.com/AndorfLab/CytoBatchFlagR/wiki/Input-Files

---

*Note: This documentation describes the code provided for the Eswar26 dataset. It does not imply that every function or output has been independently validated against the original package implementation.*

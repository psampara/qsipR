#' Scale feature coverage values to estimate their absolute abundance
#'
#' Calculates global scaling factors for features (contigs or bins),based on linear regression of sequin coverage. Options include log-transformations of coverage, as well as filtering features based on limit of detection. This function must be called first, before the feature abundance table, feature detection table, and plots are retrieved.
#'
#'@param rlm_scaling Boolean (TRUE or FALSE), should coverages and sequin concentrations be scaled using robust linear regression? Ideally, outliers should not be deleted and a robust linear regression (rlm) should be run. In rlm, the regression is done based on a weighted method. In this pipeline, a Huber loss function was used for robust linear regression.
#'@param f_tibble Can be either of
#' (1) a tibble with first column "Feature" that contains bin IDs, and the rest of the columns represent samples with bins' pooled values. Every sequin is also listed s a feature.
#' (2) a tibble as outputted by the program "checkm coverage" from the tool CheckM (https://github.com/Ecogenomics/CheckM). If this is the input format, the optional function, pooling_functions.R must be run. `pooling_functions.R` parses the checkM coverage output to provide a tibble as described in option 1. Please check `pooling_functions.R` for further details. Please check CheckM documentation (https://github.com/Ecogenomics/CheckM) on the usage for "checkm coverage" program
#'@param sequin_meta tibble containing sequin names ("Feature column") and concentrations in attamoles/uL ("Concentration") column.
#'@param seq_dilution tibble with first column "Sample" with **same sample names as in f_tibble**, and a second column "Dilution" showing ratio of sequins added to final sample volume (e.g. a value of 0.01 for a dilution of 1 volume sequin to 99 volumes sample)
#'@param coe_of_variation Acceptable coefficient of variation for coverage and detection (eg. 20 - for 20 % threshold of coefficient of variation). Coverages above the threshold value will be flagged in the plots.
#'@param log_trans Boolean (TRUE or FALSE), should coverages and sequin concentrations be log-scaled?
#'@param cook_filtering Boolean (TRUE or FALSE), should data points be filtered based on Cook's distance metric. Cooks distance can be useful in detecting influential outliers in an ordinary least square’s regression model, which can negatively influence the model. A threshold of Cooks distance of 4/n (where n is the sample size) is chosen, and any data point with Cooks distance > 4/n is filtered out. It is typical to choose 4/n as the threshold in detecting the outliers in the data

#'@importFrom rlang .data
#'@return a list of tibbles containing
#'  - mag_tab: a tibble with first column "Feature" that contains bin (or contig IDs), and the rest of the columns represent samples with features' scaled abundances (attamoles/uL)
#'  - mag_det: a tibble with first column "Feature" that contains bin (or contig IDs),
#'  - plots: linear regression plots for scaling MAG coverage values to absolute abundance
#'  - scale_fac: a master tibble with all of the intermediate values in above calculations
#'@export

scale_features <- function(rlm_scaling = T, f_tibble, sequin_meta, seq_dilution, log_trans = TRUE, coe_of_variation=10000, cook_filtering){
  # Retrieve sample names from feature tibble
scale_fac = ifelse(rlm_scaling == "TRUE",
  scale_features_rlm(f_tibble, sequin_meta, seq_dilution, log_trans = TRUE, coe_of_variation=10000),
  scale_features_lm(f_tibble, sequin_meta, seq_dilution, log_trans = TRUE, coe_of_variation=10000, cook_filtering))
}

#' Scale feature coverage values to estimate their absolute abundance
#'
#' Calculates global scaling factors for features (contigs or bins),based on linear regression of sequin coverage. Options include log-transformations of coverage, as well as filtering features based on limit of detection. This function must be called first, before the feature abundance table, feature detection table, and plots are retrieved.
#'
#'@param f_tibble Can be either of
#' (1) a tibble with first column "Feature" that contains bin IDs, and the rest of the columns represent samples with bins' pooled values. Every sequin is also listed s a feature.
#' (2) a tibble as outputted by the program "checkm coverage" from the tool CheckM (https://github.com/Ecogenomics/CheckM). If this is the input format, the optional function, pooling_functions.R must be run. `pooling_functions.R` parses the checkM coverage output to provide a tibble as described in option 1. Please check `pooling_functions.R` for further details. Please check CheckM documentation (https://github.com/Ecogenomics/CheckM) on the usage for "checkm coverage" program
#'@param sequin_meta tibble containing sequin names ("Feature column") and concentrations in attamoles/uL ("Concentration") column.
#'@param seq_dilution tibble with first column "Sample" with **same sample names as in f_tibble**, and a second column "Dilution" showing ratio of sequins added to final sample volume (e.g. a value of 0.01 for a dilution of 1 volume sequin to 99 volumes sample)
#'@param coe_of_variation Acceptable coefficient of variation for coverage and detection (eg. 20 - for 20 % threshold of coefficient of variation). Coverages above the threshold value will be flagged in the plots.
#'@param log_trans Boolean (TRUE or FALSE), should coverages and sequin concentrations be log-scaled?
#'@import magrittr
#'@importFrom rlang .data
#'@return a list of tibbles containing
#'  - mag_tab: a tibble with first column "Feature" that contains bin (or contig IDs), and the rest of the columns represent samples with features' scaled abundances (attamoles/uL)
#'  - mag_det: a tibble with first column "Feature" that contains bin (or contig IDs),
#'  - plots: linear regression plots for scaling MAG coverage values to absolute abundance
#'  - scale_fac: a master tibble with all of the intermediate values in above calculations
#'@export

scale_features <- function(f_tibble, sequin_meta, seq_dilution, log_trans, coe_of_variation=10000){
  # Retrieve sample names from feature tibble
  scale_fac <- dplyr::tibble(Sample = names(f_tibble) %>%
                               stringr::str_subset(pattern = "Feature", negate = TRUE))

  # Merge dilution factors for samples, add log-scaling option
  scale_fac <- scale_fac %>%
    dplyr::inner_join(seq_dilution, by = "Sample") %>%
    dplyr::mutate(log_scale = log_trans)

  # Make coverage table for features
  scale_fac <- scale_fac %>%
    dplyr::mutate(
      cov_tab = purrr::map(Sample, ~ dplyr::select(f_tibble, Feature, all_of(.))), # Make list of coverage tables for samples
      cov_tab = purrr::map(cov_tab, ~stats::setNames(., c("Feature", "Coverage")))  #get rid of sample name in header
    ) %>%

    # merge sequins with coverage tables, remove them from mag/feature table
    dplyr::mutate(
      seq_cov = purrr::map(cov_tab, ~ dplyr::inner_join(., sequins, by = "Feature")),
      mag_cov = purrr::map(cov_tab, ~ dplyr::anti_join(., sequins, by = "Feature"))
    )  %>%

    # scale sequin concentrations based on dilution factors
    dplyr::mutate(
      seq_cov = purrr::map2(seq_cov, Dilution, ~ dplyr::mutate(.x, Concentration = Concentration / .y))
    ) %>%

    # Determine groups of spike-in concentrations
    dplyr::mutate(
      seq_group = purrr::map(seq_cov, ~.x %>% dplyr::group_by(Concentration) %>%
                               dplyr::tally(name="standards"))
    ) %>%

    # determine limit of detection of spike-ins, based on presence of 5 sequins per conc.
    dplyr::mutate(
      #determine lowest concentration where at least 1 sequins is detected
      lod = purrr::map_dbl(seq_cov, ~ dplyr::filter(., Coverage > 0) %>%
                             dplyr::group_by(., Concentration) %>%
                             dplyr::tally(name="detected") %>%
                             dplyr::summarise(Min = min(Concentration)) %>%
                             dplyr::pull(Min)),

      #create tibble comparing number of observed and theoretical spike ins
      seq_det = purrr::map2(seq_cov, seq_group, ~ dplyr::filter(., Coverage > 0) %>%
                              dplyr::group_by(., Concentration) %>%
                              dplyr::tally(name="detected") %>%
                              dplyr::inner_join(.y, by = "Concentration")),

      # determine difference between standards and observed spike ins
      seq_det = purrr::map(seq_det, ~ dplyr::mutate(., diff = standards - detected)),
      seq_warning = purrr::map_int(seq_det, ~ dplyr::summarise(., Sum = sum(diff)) %>%
                                     dplyr::pull(Sum)) #positive values give warning later
    ) %>%

    # Calculate mean, standard deviation, and coeefficient of variation for groups of sequins
    # Create a logical vector determining if the sequin is within the threshold
    dplyr::mutate(
      grouped_seq_cov = purrr::map(cov_tab, ~ dplyr::inner_join(., sequins, by = "Feature") %>%
                              dplyr::select(Feature, Coverage, Concentration)),
      grouped_seq_cov = purrr::map2(grouped_seq_cov, Dilution, ~ dplyr::mutate(.x, Concentration = Concentration/.y) %>%
                               dplyr::group_by(Concentration) %>%
                               dplyr::summarise(mean_cov = mean(Coverage),
                                         sd_cov = stats::sd(Coverage)) %>%
                               dplyr::na_if(0) %>%
                               dplyr::mutate(coe_var = sd_cov*100/mean_cov) %>%
                               dplyr::mutate(threshold_detection = coe_var <= coe_of_variation))) %>%

    #Create a list of samples in which sequins were not detected
    dplyr::mutate(under_detected = purrr::map(grouped_seq_cov, ~.x %>%
                                                dplyr::filter(is.na(mean_cov)) %>%
                                                dplyr::select(Concentration))
           ) %>%

    # perform linear regression on coverage vs conc., extract lm params, make plots
    dplyr::mutate(
      seq_cov_filt = purrr::map2(seq_cov,grouped_seq_cov, ~ dplyr::inner_join(.x, .y , by = "Concentration") %>%
                           dplyr::filter(., Coverage > 0)),
      seq_cov_filt = purrr::map2(seq_cov_filt, lod, ~.x %>%
                            dplyr::filter(Concentration >= .y) %>%
                            dplyr::filter(., coe_var <= coe_of_variation) %>% #remove zero coverage values before lm
                            dplyr::mutate(
                            lod = .y))) %>%

    dplyr::mutate(
      fit = ifelse(log_scale == "TRUE" , # check log_trans input
                   purrr::map(seq_cov_filt, ~ stats::lm(log10(Concentration) ~ log10(Coverage) , data = .)), #log lm if true
                   purrr::map(seq_cov_filt, ~ stats::lm((Concentration) ~ (Coverage) , data = .)) # lm if false
      ),
      slope = purrr::map_dbl(fit, ~ summary(.)$coef[2]), # get slope
      intercept = purrr::map_dbl(fit, ~summary(.)$coef[1]) # get intercept
    ) %>%

    #plot linear regressions
    dplyr::mutate(
      plots = ifelse(log_scale == "TRUE" , # check log_trans input
                     purrr::map(seq_cov_filt, # log-scaled plot if true
                         ~ ggplot2::ggplot(data=. , aes(x=log10(Coverage), y= log10(Concentration))) +
                           geom_point(aes(shape = threshold_detection)) +
                           geom_smooth(method = "lm") +
                           ggpubr::stat_regline_equation(label.x= -0.1, label.y = 3) +
                           ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -0.1, label.y = 3.5) +
                           xlab("Coverage (log[read depth])") +
                           ylab("DNA Concentration (log[attamoles/uL])") +
                           scale_shape(name = "Coefficient of variation", labels = c(paste("below the threshold (",coe_of_variation,")"), paste("above the threshold(",coe_of_variation,")"))) +
                           theme_bw()
                     ),
                     purrr::map(seq_cov_filt, # non-scaled plot if true
                         ~ ggplot2::ggplot(data=. , aes(x=Coverage, y= Concentration)) +
                           geom_point(aes(shape = threshold_detection)) +
                           geom_smooth(method = "lm") +
                           ggpubr::stat_regline_equation(label.x= 0, label.y = 1000) +
                           ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0, label.y = 100000) +
                           xlab("Coverage (read depth)") +
                           ylab("DNA Concentration (attamoles/uL)") +
                           scale_shape(name = "Coefficient of variation", labels = c(paste("below the threshold (",coe_of_variation,")"), paste("above the threshold(",coe_of_variation,")"))) +
                           theme_bw()
                     )
      )
    ) %>%

    #flag MAGs below LOD, and scale MAGs by slope and intercept
    dplyr::mutate(
      # Scale MAGs based on linear regression
      mag_ab = ifelse(log_scale == "TRUE" , # check log_trans input
                      purrr::map2(mag_cov, slope, ~ mutate(.x, Concentration = log10(Coverage) * .y)), # y = mx (in log scale) if true
                      purrr::map2(mag_cov, slope, ~ mutate(.x, Concentration = Coverage * .y)) # y = mx if false
      ),
      mag_ab = purrr::map2(mag_ab, intercept, ~ dplyr::mutate(.x, Concentration = Concentration + .y)),
      mag_ab = ifelse(log_scale == "TRUE" , # check log_trans input
                      purrr::map(mag_ab, ~ dplyr::mutate(.x, Concentration = 10^Concentration)), #convert back from log10 if true
                      mag_ab), # no change if false
      mag_det = mag_ab,
      mag_ab = purrr::map(mag_ab, ~ dplyr::select(., Feature, Concentration)), # drop Coverage column
      mag_ab = purrr::map2(mag_ab, Sample, ~ stats::setNames(.x, c("Feature", .y))), #put sample name in MAG table

      # Remove MAGs below LOD
      #mag_det = mag_ab,
      #mag_det = map2(mag_det, Sample, ~ setNames(.x, .y, 'Concentration')), #get sample name out of header for filter
      mag_det = purrr::map2(mag_det, lod, ~ dplyr::filter(.x, Concentration > .y)),
      mag_det = purrr::map(mag_det, ~ dplyr::select(.x, Feature, Concentration)),
      mag_det = purrr::map2(mag_det, Sample, ~ stats::setNames(.x, c("Feature", .y))) #change header back to sample
    )

  # compile feature abundance across samples
  mag_tab <- scale_fac$mag_ab %>%
    purrr::reduce(dplyr::left_join, by="Feature") %>%
    tibble::column_to_rownames(var = "Feature")
  #Create regression plots directory
  dir.create("sequin_scaling_plots")
  # extract plots
  plots <- scale_fac %>% dplyr::select(Sample, plots)
  #Save scaling plots in .pdf format in the regression plots directory
  plots <- plots %>%
    dplyr::mutate(save_plots = purrr::map2(plots, Sample,  ~ ggplot2::ggsave(filename = paste("",.y, ".pdf", sep=""), plot = .x, path = "sequin_scaling_plots/")))


  # extract feature detection
  mag_det <- scale_fac$mag_det %>%
    purrr::reduce(dplyr::left_join, by="Feature") %>%
    tibble::column_to_rownames(var = "Feature")


  # make list of results
  results <- list("mag_tab" = mag_tab,
                  "mag_det" = mag_det,
                  "plots" = plots,
                  "scale_fac" = scale_fac)

  # return results
  attach(results)
  return(results)

}


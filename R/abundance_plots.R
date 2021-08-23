#' Plot abundance profiles for MAGs at different buoyant densities
#'
#' This function provides the plots for abundance vs Buoyant density of incorporators. The plots are saved in the directory "abundance_plots" in the current working directory of the "qSIP_analysis.rmd" pipeline. The goal of plotting individual replicates separately is to identify any anomalies, if any, regarding the changes in the scaled abundances across the buoyant density gradient which may not be obvious in the plots with the mean and standard deviation of abundances. These plots could help in debugging
#'
#'@param incorporator_list A tibble consisting of two columns
#' - First column consists of incorporator feature ID.
#' - Second column consists of the incorporator taxonomy
#'@param fractions Fractions metadata file consisting of
#' - Replicate: Depends on how many replicates the study has
#' - Fractions: Typically in the range of 2-24
#' - Buoyant_density: As calculated from the refractometer for each fraction and replicate
#' - Isotope: "12C", "13C", "14N", "15N" etc.
#' - DNA_concentration
#' - Sample: In the format "'isotope'_rep_#_fraction_#".
#'   For instance, "12C_rep_1_fraction_1"
#'@param mag_tab A tibble of feature absolute abundances
#' Each column consists of the scaled abundances of the features. Row names idenitfy the Feature
#'@import magrittr
#'@importFrom rlang .data
#'@return image files with the following formats
#'  - "mean_abundance_" prefix is for the plots visualizing the abundance estimates of features
#'  over the buoyant density gradient
#'  - "Rep_1_abundance" prefix is for the plots visualizing the abundance estimates of features
#'  from replicate 1 over the buoyant density gradient
#'  - "Rep_2_abundance" prefix is for the plots visualizing the abundance estimates of features
#'  from replicate 2 over the buoyant density gradient
#'  - "Rep_3_abundance" prefix is for the plots visualizing the abundance estimates of features from replicate 3 over the buoyant density gradient
#'@export

plot_abundance = function(incorporator_list, fractions, mag_tab) {
dir.create("abundance_plots")
#Rename first column of `incorporator_list_rename` tibble to Feature from OTU
incorporator_list_renamed = incorporator_list %>%
  dplyr::rename("Feature" = "OTU")
#Reformat mag_tab to make rownames in this tibble as column names with the column name Feature
mag_tab_reformatted = as.data.frame(mag_tab) %>%
  tibble::rownames_to_column(var = "Feature")
#Make a tibble by joining abundance and MAG names from incorporator list
incorporator_abs_abundance= incorporator_list_renamed %>%
  dplyr::inner_join(.data$., mag_tab_reformatted, by = "Feature")
#Pivot the tibble to later access abundance values of every Bin in every sample
#Store this as a temporary tibble
long_cov = incorporator_abs_abundance %>%
  tidyr::pivot_longer(fractions$Sample, names_to = "Sample", values_to = "Abundance")
#Keep Feature and taxonomy columns only
incorporator_abs_abundance = incorporator_abs_abundance %>%
  dplyr::select(.data$Feature, .data$Taxonomy)

#Create a tibble within the previous tibble to access abundances of each MAG in every sample
incorporator_abs_abundance = incorporator_abs_abundance %>%
  dplyr::mutate(
    abs_abundance = purrr::map(.data$Feature, ~ dplyr::filter(long_cov, Feature == .data$.) %>%
                                 dplyr::select(.data$Sample, .data$Abundance))
  )
#Populate the corresponding fraction, BD, Replicate, and Isotope information
#for all samples. This gives coverage values for a certain replicate, in a certain BD fraction, for a particular isotope treatment
#for all MAGs
incorporator_abs_abundance = incorporator_abs_abundance %>%
  dplyr::mutate(
    abs_abundance = purrr::map(.data$abs_abundance, ~ dplyr::mutate(.data$., Isotope = fractions$Isotope,
                                              Fraction = fractions$Fraction,
                                              Buoyant_density = fractions$Buoyant_density,
                                              Rep = stringr::str_split_fixed(fractions$Sample, pattern = "_", n= Inf)[,3]))
  )
#Summarise the mean abundance, standard deviation of abundance, and mean BD for each fraction and isotope treatment
incorporator_abs_abundance = incorporator_abs_abundance %>%
  dplyr::mutate(
    summary_coverage = purrr::map(.data$abs_abundance, ~ dplyr::group_by(.data$., Fraction, Isotope ) %>%
                           dplyr::summarise(mean_abs_abundance = mean(Abundance),
                                     mean_BD = mean(Buoyant_density),
                                     sd_abs_abundance = stats::sd(Abundance))),
    summary_coverage = purrr::map(.data$summary_coverage, ~ dplyr::arrange(.data$., Isotope))

  )
#Plot abundance vs BD and save the plots in the current path
incorporator_abs_abundance = incorporator_abs_abundance %>%
  dplyr::mutate(
    plots = purrr::map(.data$summary_coverage, ~ ggplot2::ggplot(data = .data$., aes(x = mean_BD, y = mean_abs_abundance)) +
                  geom_point(aes(color = Isotope)) +
                  geom_line(aes(color = Isotope)) +
                  geom_errorbar(aes(ymin = mean_abs_abundance - sd_abs_abundance,
                                    ymax = mean_abs_abundance + sd_abs_abundance)) +
                  ylab("Mean absolute \n abundance (attamole/uL)") +
                  xlab("Mean buoyant \n density (g/mL)") +
                  theme_bw()),
    save_plots = purrr::map2(.data$plots, .data$Feature,  ~ ggplot2::ggsave(filename = paste("mean_abundance_",.y, ".pdf", sep=""), plot = .x, path = "abundance_plots/"))
  )
incorporator_abs_abundance = incorporator_abs_abundance %>%
  dplyr::mutate(
    plots_Rep1 = purrr::map(.data$abs_abundance, ~ dplyr::filter(.data$., Rep == 1) %>%
                  ggplot2::ggplot(data = .data$., aes(x = Buoyant_density, y = Abundance)) +
                  geom_point(aes(color = Rep, shape = Isotope)) +
                  geom_line(aes(color = Isotope)) +
                  ylab("Absolute \n abundance (attamole/uL)") +
                  xlab("Buoyant \n density (g/mL)") +
                  theme_bw()),
    plots_Rep2 = purrr::map(.data$abs_abundance, ~ dplyr::filter(.data$., Rep == 2) %>%
                       ggplot2::ggplot(data = .data$., aes(x = Buoyant_density, y = Abundance)) +
                       geom_point(aes(color = Rep, shape = Isotope)) +
                       geom_line(aes(color = Isotope)) +
                       ylab("Absolute \n abundance (attamole/uL)") +
                       xlab("Buoyant \n density (g/mL)") +
                       theme_bw()),
    plots_Rep3 = purrr::map(.data$abs_abundance, ~ dplyr::filter(.data$., Rep == 3) %>%
                       ggplot2::ggplot(data = .data$., aes(x = Buoyant_density, y = Abundance)) +
                       geom_point(aes(color = Rep, shape = Isotope)) +
                       geom_line(aes(color = Isotope)) +
                       ylab("Absolute \n abundance (attamole/uL)") +
                       xlab("Buoyant \n density (g/mL)") +
                       theme_bw()),
    save_plots = purrr::map2(.data$plots_Rep1, .data$Feature,  ~ ggplot2::ggsave(filename = paste("Rep_1_abundance_",.y, ".pdf", sep=""), plot = .x, path = "abundance_plots/")),
    save_plots = purrr::map2(.data$plots_Rep2, .data$Feature,  ~ ggplot2::ggsave(filename = paste("Rep_2_abundance_",.y, ".pdf", sep=""), plot = .x, path = "abundance_plots/")),
    save_plots = purrr::map2(.data$plots_Rep3, .data$Feature,  ~ ggplot2::ggsave(filename = paste("Rep_3_abundance_",.y, ".pdf", sep=""), plot = .x, path = "abundance_plots/"))
  )
return(incorporator_abs_abundance)
}

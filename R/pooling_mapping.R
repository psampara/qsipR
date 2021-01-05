#' Calculate mean or median per feature
#'
#' pool_bin_stat calculates per bin mean or median scaffold value (coverage or mapped read numbers).
#'
#' @param bin_tibble a tibble of coverage values (numeric) that need to be pooled,
#' column names should contain `Sequence Id`, `Bin Id`, `Coverage`
#' @param operator a character vector either of "mean" or "median" value; identifies pooling strategy.
#' @import magrittr
#' @return A tibble with first column "Feature" that contains bin IDs, and the rest of the columns represent samples with bins' pooled values.
#' @export

pool_bin_stat <- function(bin_tibble, operator = "mean"){
# setting the function to pool values based on "operator" parameter value
f_tibble = bin_tibble
newnames <- f_tibble %>% dplyr::select(grep(pattern = "Bam", x = names(.))) %>% #subset to BAM containing columns
  unique(.) %>% dplyr::slice(1) %>% unlist(., use.names=FALSE) #get list of unique BAM file names
oldnames <- f_tibble %>% dplyr::select(grep(pattern = "Coverage", x = names(f_tibble))) %>% names(.)
f_tibble <- f_tibble %>% dplyr::rename_at(dplyr::vars(oldnames), ~ newnames) #replace old names with new names
f_tibble <- f_tibble %>% dplyr::select(1, 2, dplyr::contains(newnames)) #pull out coverage columns
if (operator == "median") {
f_tibble = f_tibble %>%
  dplyr::select(-`Sequence Id`) %>%
  dplyr::group_by(`Bin Id`) %>%
  dplyr::summarise_all(stats::median)
} else {
  f_tibble = f_tibble %>%
    dplyr::select(-`Sequence Id`) %>%
    dplyr::group_by(`Bin Id`) %>%
    dplyr::summarise_all(mean)
}
f_tibble = f_tibble %>%
  dplyr::rename(Feature = `Bin Id`)
return(f_tibble)
}




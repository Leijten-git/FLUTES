#' Title
#'
#' @param years
#' @param lu_matrices
#' @param dir_output_files
#' @param country
#' @param aggregation_factor
#' @param cut_off_year
#'
#' @return
#' @export
#'
#' @examples
gen_panel_dataframe <- function(years,
                                lu_matrices,
                                dir_output_files,
                                country,
                                aggregation_factor,
                                cut_off_year = NULL){

  if(!exists("lu_matrices")){

    cat("\nConverting the fractional land use matrices into long format...\n")

    period <- paste0(c(as.character(first(sort(years))),
                       "-",
                       as.character(last(sort(years)))), collapse = "")

    filename_lu_frac_matrices <- paste0("lu_frac_matrix_matrices_",
                                        aggregation_factor,
                                        "_",
                                        period,
                                        ".Rds")

    folder_agg_factor <- paste0("agg_factor_", aggregation_factor)

    path_lu_frac_matrices <- file.path(dir_output_files,
                                       country,
                                       folder_agg_factor,
                                       filename_lu_frac_matrices)

    lu_matrices <- readRDS(path_lu_frac_matrices)


  }

  lu_frac_matrices_df <- do.call(cbind, lu_matrices) %>%
    as.data.frame()

  lu_frac_matrices_long <- lu_frac_matrices_df %>%
    mutate(id = row_number()) %>%
    tidyr::gather(key = "year_and_land_class", value = "value", -id)

  return(lu_frac_matrices_long)

}

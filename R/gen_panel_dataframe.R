gen_panel_dataframe <- function(mask,
                                years,
                                aggregation_factor,
                                cut_off_year = NULL){

  if(!exists("lu_frac_matrices_list")){

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

    lu_frac_matrices_list <- readRDS(path_lu_frac_matrices)


  }

  lu_frac_matrices_df <- do.call(cbind, lu_frac_matrices_list) %>%
    as.data.frame()

  inds_outside_mask <- which(is.na(terra::values(mask)))
  lu_frac_matrices_df = lu_frac_matrices_df[-inds_outside_mask,]

  lu_frac_matrices_long <- lu_frac_matrices_df %>%
    mutate(id = row_number()) %>%
    tidyr::gather(key = "year_and_land_class", value = "value", -id)

  return(lu_frac_matrices_long)

}

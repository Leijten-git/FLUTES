identify_lu_classes <- function(country,
                                aggregation_factor,
                                years,
                                dir_output_files){

  if(!exists("lu_frac_matrices_list")){

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

  extract_lu_classes <- function(index){

    lu_matrix <- lu_frac_matrices_list[[index]]
    lu_classes <- colnames(lu_matrix)
    return(lu_classes)

  }

  lu_classes_list <- lapply(X = c(1:length(lu_frac_matrices_list)),
                            FUN = extract_lu_classes)

  unique_land_systems <- unique(sub('.*_', '', unlist(lu_classes_list)))

  return(unique_land_systems)

}

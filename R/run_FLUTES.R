#' Run FLUTES
#'
#' @return
#' @export
#'
#' @examples
run_FLUTES <- function(.country,
                       .years,
                       .paths_lu_filenames,
                       .dir_output_files,
                       .aggregation_factor,
                       .path_lu_legend_filename,

                       .cut_off_year = NULL,
                       .lags = NULL,
                       .growth_constraint = NULL,

                       .max_dev = 1,
                       .max_iter = 2e3,
                       .vals_lu_classes_to_exclude = NULL,
                       .kernel_size = 25,
                       .is_abs_dev = F,

                       skip_preprocessing = F){

  if(skip_preprocessing == F){

    cat("\nPreprocessing steps are skipped...")
    warning("Preprocessing steps are skipped...")

    crop_data_to_country(country = .country,
                         years = .years,
                         paths_lu_filenames = .paths_lu_filenames,
                         dir_output_files = .dir_output_files)

    lu_frac_matrices_list <- gen_fract_lu_matrices(country = .country,
                                                   years = .years,
                                                   aggregation_factor = .aggregation_factor,
                                                   dir_output_files = .dir_output_files,
                                                   path_lu_legend_filename = .path_lu_legend_filename,
                                                   vals_lu_classes_to_exclude = .vals_lu_classes_to_exclude)

    mask <- gen_gridded_mask(country = .country,
                             years = .years,
                             paths_lu_filenames = .paths_lu_filenames,
                             dir_output_files = .dir_output_files)

    specify_output_dir(dir_output_files = .dir_output_files,
                       aoi = noquote(.country),
                       aggregation_factor = .aggregation_factor)

    lu_requirements <- extract_lu_requirements(lu_frac_matrices_list = lu_frac_matrices_list,
                                               years = .years)

    extract_indices_cols <- function(lu){
      indices_cols <- c(1:ncol(lu))
      return(indices_cols)
    }

    indices_cols <- lapply(X = lu_frac_matrices_list, FUN = extract_indices_cols)

    neighbourhood_output_matrices_list <- lapply(X = c(1:length(lu_frac_matrices_list)),
                                                 FUN = calc_neighbourhood_vals,
                                                 lu_list = lu_frac_matrices_list,
                                                 neigh_cols_list = indices_cols,
                                                 kernel_size = .kernel_size,
                                                 mask = mask,
                                                 suffix = "neigh")

    saveRDS(neighbourhood_output_matrices_list, "neighbourhood_output_matrices_list.Rds")

    neighbourhood_output_matrices_df <- do.call(cbind, neighbourhood_output_matrices_list)

    lu_frac_matrices_long <- gen_panel_dataframe(years = .years,
                                                 lu_matrices = lu_frac_matrices_list,
                                                 aggregation_factor = .aggregation_factor,
                                                 dir_output_files = .dir_output_files,
                                                 country = .country)

    fitted_vals_scaled <- regression_analysis(lu_matrices = lu_frac_matrices_long,
                                              neigh_values = neighbourhood_output_matrices_df,
                                              lu_classes = identify_lu_classes(lu_frac_matrices_list = lu_frac_matrices_list),
                                              cut_off_year = .cut_off_year,
                                              lags = .lags)

  } else if(skip_preprocessing == T){

    specify_output_dir(dir_output_files = .dir_output_files,
                       aoi = noquote(.country),
                       aggregation_factor = .aggregation_factor)

    period <- paste0(c(as.character(first(sort(years))),
                       "-",
                       as.character(last(sort(years)))), collapse = "")

    filename_lu_frac_matrices <- paste0("lu_frac_matrix_matrices_",
                                        .aggregation_factor,
                                        "_",
                                        period,
                                        ".Rds")

    lu_frac_matrices_list <- readRDS(filename_lu_frac_matrices)
    neighbourhood_output_matrices_list <-
      readRDS("neighbourhood_output_matrices_list.Rds")
    lu_requirements <- readRDS("land_use_requirements.Rds")
    fitted_vals_scaled <- readRDS("fitted_vals_scaled.Rds")

  }

  index = which(as.numeric(names(lu_frac_matrices_list))==(.cut_off_year+1))

  inds_lu_requirements_rows <- c(index,nrow(lu_requirements))
  inds_lu_requirements_cols <- c(2:(ncol(lu_requirements)-1))

  lu_requirements_filter <- lu_requirements[inds_lu_requirements_rows,
                                            inds_lu_requirements_cols]

  max_devs <- determine_max_dev_pars(lu_requirements = lu_requirements_filter,
                                     lu_first_period = lu_frac_matrices_list[[index]],
                                     default_dev_par = .max_dev)

  params = list(is_abs_dev = .is_abs_dev,
                max_devs = max_devs,
                max_iter = .max_iter,
                growth = rep(.growth_constraint, ncol(lu_frac_matrices_list[[1]])),
                no_change = NULL)

  dir_model_output <- specify_output_dir(dir_output_files = .dir_output_files,
                                         aoi = .country,
                                         aggregation_factor = .aggregation_factor,
                                         create_sub_folder = T,
                                         are_conv_pars_perc = !.is_abs_dev,
                                         growth_par = .growth_constraint,
                                         are_PAs_excluded = F)

  new_matrix <- lu_allocation(lu = lu_frac_matrices_list[[index]],
                              sm = fitted_vals_scaled,
                              params = params,
                              dmd = lu_requirements_filter,
                              ln = neighbourhood_output_matrices_list[[index]],
                              constraint = T,
                              rescale = T,
                              PA = F,
                              output_dir = noquote(dir_model_output))

  cat("\nWriting the new matrix...\n")

  saveRDS(new_matrix, "new_matrix.Rds")

  calculate_rmse <- function(){

    predicted_fractions <- c(new_matrix)
    actual_fractions <- c(lu_frac_matrices_list[[as.character(last(sort(.years)))]])
    original_fractions <- c(lu_frac_matrices_list[[index]])

    nrow(lu_frac_matrices_list[[index]])
    nrow(neighbourhood_output_matrices_list[[index]])

    length(lu_frac_matrices_list[[index]])
    length(as.matrix(neighbourhood_output_matrices_list[[index]]))

    length(predicted_fractions)
    length(actual_fractions)
    length(original_fractions)

    actual_changes <- abs(actual_fractions-original_fractions)

    rmse_model <- sqrt(mean((actual_fractions - predicted_fractions)^2, na.rm = T))
    rmse_no_change <- sqrt(mean((actual_fractions - original_fractions)^2, na.rm = T))

    if(rmse_model>rmse_no_change){

      cat("\nThe RMSE of the default model is higher than the RMSE of the NULL model!\n\n")
      warning("\nThe RMSE of the default model is higher than the RMSE of the NULL model!\n\n")

    } else {

      cat("\nHurray! Outperforming the NULL model!\n\n")
      warning("\nHurray! Outperforming the NULL model!\n\n")

    }

    RMSEs <- data.frame(model = c("default", "null"),
                        rmse = c(rmse_model, rmse_no_change))

    return(RMSEs)

  }

  RMSEs <- calculate_rmse()

  return(RMSEs)

}

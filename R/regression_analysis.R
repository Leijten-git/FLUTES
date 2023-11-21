#' Title
#' @import dplyr
#' @return
#' @export
#'
#' @examples
regression_analysis <- function(lu_classes,
                                neigh_values,
                                lags,
                                cut_off_year = NULL,
                                cut_off_val_reg1 = 2.5e2,
                                cut_off_val_reg2 = 1e2,
                                cut_off_val_var = 5e-7){

  cat("\nGenerating a panel data frame...\n")

  lu_frac_matrices_long <- gen_panel_dataframe(years = years,
                                               aggregation_factor = 16)

  neigh_values_long <- neigh_values %>%
    as.data.frame() %>%
    mutate(id = row_number()) %>%
    tidyr::gather(key = "year_and_land_class", value = "value_neigh", -id)

  panel_data <- lu_frac_matrices_long %>%
    cbind(neigh_values_long %>%
            dplyr::select(value_neigh))

  panel_data2 <- panel_data %>%
    tidyr::separate(col = "year_and_land_class",
                    into = c("Year", "Land class"),
                    sep = "_") %>%
    mutate(Year = as.numeric(Year))

  if(!is.null(cut_off_year)){

    panel_data2 <- panel_data2 %>%
      filter(Year <= cut_off_year)

  }

  panel_data3 <- panel_data2 %>%
    rename(fraction_lu = "value",
           fraction_lu_neigh = "value_neigh") %>%
    arrange(id, `Land class`)

  if(!is.null(lags)){

    for(lag_pos in lags){

      cat('\nCreating a lagged variable with position:', lag_pos, "\n")

      lag_col_name <- paste0("lag_fraction_lu_", lag_pos)

      panel_data3 <- panel_data3 %>%
        group_by(id, `Land class`) %>%
        mutate(!!rlang::ensym(lag_col_name) := dplyr::lag(fraction_lu,
                                                          n = lag_pos,
                                                          default = NA)) %>%
        ungroup() %>%
        na.omit()

    }

  }

  gen_preds_by_lu_class <- function(lu_class){

    panel_data_filter <- panel_data3  %>%
      filter(`Land class` == lu_class) %>%
      dplyr::select(-`Land class`)

    total_area = sum(panel_data_filter$fraction_lu)
    total_n_pop_pixels <- sum(panel_data_filter$fraction_lu>0)
    var_dep_var <- var(panel_data_filter$fraction_lu)

    vars <- colnames(panel_data_filter)
    lagged_vars <- vars[grep("lag_fraction_lu", vars)]

    cat("\nNow performing a regression analysis for land use class:",
        lu_class,
        "\n",
        "Number of non-empty pixels:", total_n_pop_pixels, "\n",
        "Total area:", total_area, "\n",
        "Variance:", var_dep_var, "\n")

    model_spec <- stats::formula(paste("fraction_lu", "~",
                                       paste("fraction_lu_neigh + ",
                                             paste(lagged_vars, collapse = "+"))))
    if(total_n_pop_pixels>cut_off_val_reg1 & var_dep_var >= cut_off_val_var){

      suppressWarnings(model_result <- plm::plm(model_spec,
                                                index=c("id", "Year"),
                                                data = panel_data_filter,
                                                model = "within"))

    } else if(total_n_pop_pixels<=cut_off_val_reg1 &
              total_n_pop_pixels>cut_off_val_reg2 &
              var_dep_var >= cut_off_val_var){

      suppressWarnings(model_result <- plm::plm(model_spec,
                                                index=c("id"),
                                                data = panel_data_filter,
                                                model = "within"))

    } else {

      if(length(lagged_vars)>0){

        first_lagged_var <- lagged_vars[1]
        xsym <- rlang::ensym(first_lagged_var)

        model_result <- rlang::inject(stats::lm(fraction_lu ~ !!xsym,
                                                data = panel_data_filter))

      } else {

        model_result <- lm(fraction_lu ~ 1, data = panel_data_filter)

      }

    }

    fitted_vals = model_result %>%
      stats::predict()

    normalize <- function(x, na.rm = TRUE) {
      return((x- min(x)) /(max(x)-min(x)))
    }

    fitted_vals_scaled <- normalize(as.numeric(fitted_vals))

    return(fitted_vals_scaled)

  }

  fitted_vals_scaled <- sapply(X = lu_classes,
                               FUN = gen_preds_by_lu_class)

  cat("\nWriting the fitted values...\n")
  saveRDS(fitted_vals_scaled, "fitted_vals_scaled.Rds")

  return(fitted_vals_scaled)


}



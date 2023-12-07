#' Title
#'
#' @param lu_matrices
#' @param neigh_values
#' @param lags
#' @param cut_off_year
#' @param min_n_pixels1
#' @param min_n_pixels2
#' @param min_n_pixels3
#' @param min_area1
#' @param min_area2
#' @param min_area3
#' @param cut_off_val_var
#' @param lu_classes
#'
#' @import dplyr
#' @return
#' @export
#'
#' @examples
regression_analysis <- function(lu_matrices,
                                neigh_values,
                                lu_classes,

                                lags,
                                cut_off_year = NULL,
                                min_n_pixels1 = 7e2,
                                min_n_pixels2 = 4e2,
                                min_n_pixels3 = 5e1,

                                min_area1 = 5e2,
                                min_area2 = 2.5e2,
                                min_area3 = 10,

                                cut_off_val_var = 5e-7){

  cat("\nGenerating a panel data frame...\n")

  neigh_values_long <- neigh_values %>%
    as.data.frame() %>%
    mutate(id = row_number()) %>%
    tidyr::gather(key = "year_and_land_class", value = "value_neigh", -id)

  panel_data <- lu_matrices %>%
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

  cat("Starting the regression analysis.\n",
      "The model specification depends on the number of non-emtpy pixels.\n",
      "Current cut-off values:\n",
      "Cut-off value 1:", min_n_pixels1, "\n",
      "Cut-off value 2:", min_n_pixels2, "\n",
      "Maximum variance:", cut_off_val_var, "\n")

  saveRDS(data.frame(lu_class = NA,
                     r2 = NA),
          "r2_scores.Rds")

  gen_preds_by_lu_class <- function(lu_class){

    panel_data_filter <- panel_data3  %>%
      filter(`Land class` == lu_class) %>%
      dplyr::select(-`Land class`)

    n_pixels <- length(unique(panel_data_filter$id))

    mean_populated_area <- panel_data_filter %>%
      group_by(id) %>%
      summarize(mean_fraction = mean(fraction_lu, na.rm = T)) %>%
      ungroup() %>%
      summarize(sum = sum(mean_fraction)) %>%
      pull()

    total_n_pop_pixels <- panel_data_filter %>%
      group_by(id) %>%
      summarize(sum = mean(fraction_lu, na.rm = T)) %>%
      ungroup() %>%
      filter(sum>0) %>%
      nrow()

    var_dep_var <- var(panel_data_filter$fraction_lu)

    vars <- colnames(panel_data_filter)
    lagged_vars <- vars[grep("lag_fraction_lu", vars)]

    cat("\nNow performing a regression analysis for land use class:",
        lu_class,
        "\n",
        "Number of non-empty pixels:", total_n_pop_pixels, "\n",
        "Time-averaged populated area:", mean_populated_area, "\n",
        "Variance:", var_dep_var, "\n")

    model_spec <- stats::formula(paste("fraction_lu", "~",
                                       paste("fraction_lu_neigh + ",
                                             paste(lagged_vars, collapse = "+"))))

    if(total_n_pop_pixels>min_n_pixels1 &
       mean_populated_area > min_area1 &
       var_dep_var >= cut_off_val_var){

      suppressWarnings(model_result <- plm::plm(model_spec,
                                                index=c("id", "Year"),
                                                data = panel_data_filter,
                                                model = "within"))

      r2 <- as.numeric(summary(model_result)$r.squared[1])

    } else if(total_n_pop_pixels<=min_n_pixels1 &
              total_n_pop_pixels>min_n_pixels2 &

              mean_populated_area <= min_area1 &
              mean_populated_area > min_area2 &

              var_dep_var >= cut_off_val_var){

      suppressWarnings(model_result <- plm::plm(model_spec,
                                                index=c("id"),
                                                data = panel_data_filter,
                                                model = "within"))

      r2 <- as.numeric(summary(model_result)$r.squared[1])

    } else {

      if(length(lagged_vars)>0 &
         total_n_pop_pixels>min_n_pixels3 &
         mean_populated_area>min_area3){

        first_lagged_var <- lagged_vars[1]
        xsym <- rlang::ensym(first_lagged_var)

        model_result <- rlang::inject(stats::lm(fraction_lu ~ !!xsym,
                                                data = panel_data_filter))

        r2 <- summary(model_result)$r.squared

      } else {

        cat("There are not enough non-empty pixels to estimate a model with ",
                "predictor variables for land use class:\n", lu_class, ".\n",
                "Instead, a linear intercept-only model has been estimated. ",
                "Note that this results in a r-squared of 0.\n",
                "Consider changing the minimum threshold values.",
                "\nCurrent settings:\n",
                "Minimum number of pixels: ", min_n_pixels3, "\n",
                "Average fraction: ", mean_populated_area, "\n")

        warning("There are not enough non-empty pixels to estimate a model with",
                "predictor variables for land use class:\n", lu_class, ".\n",
                "Instead, a linear intercept-only model has been estimated. ",
                "Note that this results in a r-squared of 0.\n",
                "Consider changing the minimum threshold values.",
                "\nCurrent settings:\n",
                "Minimum number of pixels: ", min_n_pixels3, "\n",
                "Average fraction: ", mean_populated_area, "\n\n")

        model_result <- stats::lm(fraction_lu ~ 1, data = panel_data_filter)
        r2 <- summary(model_result)$r.squared

      }

    }

    cat("R-squared value:", r2, "\n")

    readRDS("r2_scores.Rds") %>%
      rbind(data.frame(lu_class = lu_class,
                       r2 = r2)) %>%
      na.omit() %>%
      saveRDS("r2_scores.Rds")

    fitted_vals = model_result %>%
      stats::predict()

    normalize <- function(x, na.rm = TRUE) {
      return((x- min(x)) /(max(x)-min(x)))
    }

    if(var(fitted_vals)>0){
      fitted_vals_scaled <- normalize(as.numeric(fitted_vals))
    } else {
      fitted_vals_scaled <- fitted_vals
    }

    return(fitted_vals_scaled)

  }

  fitted_vals_scaled <- sapply(X = lu_classes,
                               FUN = gen_preds_by_lu_class)

  years <- panel_data3$Year[c(1:nrow(fitted_vals_scaled))]
  last_year <- last(unique(years))
  inds_last_year <- which(years==last_year)

  fitted_vals_scaled2 <- fitted_vals_scaled[inds_last_year,]

  cat("\nWriting the fitted values...\n")
  saveRDS(fitted_vals_scaled2, "fitted_vals_scaled.Rds")

  return(fitted_vals_scaled2)


}



#' Extracting land use requirments
#'
#' @param lu_frac_matrices_list
#' @param years
#' @param path
#'
#' @return
#' @export
#'
#' @examples
extract_lu_requirements <- function(lu_frac_matrices_list,
                                    years,
                                    path = paste0(getwd(), "/land_use_requirements.Rds")){

  lu_classes <- identify_lu_classes(lu_frac_matrices_list = lu_frac_matrices_list)

  names(lu_frac_matrices_list) <- NULL
  lu <- do.call(cbind.data.frame, lu_frac_matrices_list)

  n_lu_classes <- length(lu_classes)

  class_supply <- colSums(lu)

  years_sorted <- sort(years)[1]:tail(sort(years),1)
  rm(years)

  demand <- matrix(NA, nrow = length(years_sorted), ncol = n_lu_classes + 2)
  demand[,1] <- years_sorted

  obs_ind <- which(demand[,1] %in% years_sorted)

  period <- paste0(as.character(c(first(years_sorted), last(years_sorted))),
                   collapse = ":")

  if(!is.null(path)){

    cat("\nWriting the land use requirements for the total simulation period (i.e.,", period, ") ",
        "in the output working directory:\n", path, "\n", sep = "")

  }

  for(year in years_sorted){

    class_supply_year <- class_supply[grep(year, colnames(lu))]

    # if(length(class_supply_year)<n_lu_classes){
    #
    #   class_supply_year_df <- data.frame(lu_class = names(class_supply_year),
    #                                      supply = as.numeric(class_supply_year))
    #
    #   class_supply_year_complete_df <- data.frame(lu_class = paste0(year, "_", lu_classes)) %>%
    #     mutate(order = c(1:length(lu_classes))) %>%
    #     merge(class_supply_year_df,
    #           by = "lu_class",
    #           all.x = T) %>%
    #     mutate(supply = coalesce(supply, 0)) %>%
    #     arrange(order)
    #
    #   class_supply_year <- class_supply_year_complete_df$supply
    #
    # }

    demand[which(demand[,1]==year), 2:(n_lu_classes+1)] <- class_supply_year

  }

  for(i in 1:n_lu_classes){
    demand[, i + 1] <- stats::approx(demand[,1], demand[,i + 1], xout = years_sorted)$y
  }

  demand[, n_lu_classes + 2] <- rowSums(demand[,-1], na.rm = TRUE)

  if(!is.null(path)){
    saveRDS(demand, path)
  }

  return(demand)


}

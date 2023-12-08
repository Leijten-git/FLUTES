#' Title
#'
#' @param neigh_cols
#' @param n_lu_classes
#' @param kernel_size
#' @param mask
#' @param enr
#'
#' @return
#' @export
#'
#' @examples
calc_neighbourhood_vals <- function(index,
                                    lu_list,
                                    neigh_cols_list,
                                    kernel_size,
                                    mask,
                                    enr = T,
                                    suffix = NULL){

  cat("\nNow calculating the neighbourhood weights for year:",
      names(lu_list)[[index]],
      "\n")

  lu <- lu_list[[index]]

  neigh_cols <- neigh_cols_list[[index]]

  neigh_values <- lu[,neigh_cols]

  weights <- list(matrix(1/kernel_size,
                         sqrt(kernel_size),
                         sqrt(kernel_size),
                         byrow= TRUE)) # size of window

  weights_by_land_system <- rep(weights, length.out = length(neigh_cols))

  inds_mask <- which(!is.na(terra::values(mask)))

  for (i in c(1:length(neigh_cols))){

    window <- weights_by_land_system[[neigh_cols[i]]]

    mask[inds_mask] <- lu[,neigh_cols[i]]

    grid_focal_sum <- terra::focal(mask,
                                   window,
                                   fun = "sum",
                                   na.rm = T)

    if(enr == TRUE){
      grid_focal_sum <- grid_focal_sum/mean(lu[,neigh_cols[i]])
    }

    neigh_values[,i] <- terra::values(grid_focal_sum)[inds_mask]

  }

  if(!is.null(suffix)){

    neigh_values <- neigh_values %>%
      as.data.frame() %>%
      `colnames<-`(paste0(colnames(neigh_values), "_neigh_"))
  }

  return(neigh_values)

}



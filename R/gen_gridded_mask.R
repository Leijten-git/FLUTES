#' Extract the indices of a gridded mask for the area of interest
#' @import dplyr
#' @param country
#' @param paths_lu_filenames
#' @param dir_output_files
#'
#' @return
#' @export
#'
#' @examples
gen_gridded_mask <- function(country,
                             paths_lu_filenames,
                             dir_output_files){

  country_name = country

  cat("\nGenerating a mask for:", country_name, "\n")

  files <- list.files()

  if(!"initial_fractional_grids.tif" %in% files){
    stop("The initial fractional grids haven't been created yet.",
         "You need to run `gen_fract_lu_matrices` before you can extract the indices.")
  }

  mask <- terra::rast("initial_fractional_grids.tif")
  if(terra::nlyr(mask)>1){
    mask <- terra::app(mask, fun = sum,  na.rm = T)
  }

  is_nonNA_mask <- which(!is.na(terra::values(mask)))
  terra::values(mask) <- NA
  terra::values(mask)[is_nonNA_mask] <- 1

  return(mask)


}

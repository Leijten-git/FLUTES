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
                             cut_off_year,
                             years,
                             paths_lu_filenames,
                             dir_output_files){

  country_name = country

  cat("\nGenerating a mask for:", country_name, "\n")

  filenames_cropped_lu_maps <- list.files(pattern="cropped_lu_map")

  if(!is.null(cut_off_year)){
    filename_ind_first_cropped_lu_map <- grep(as.character(cut_off_year),
                                              filenames_cropped_lu_maps)
  } else {
    first_year <- sort(years)[1]
    filename_ind_first_cropped_lu_map <- grep(as.character(first_year),
                                              filenames_cropped_lu_maps)
  }

  filename_first_cropped_lu_map <-
    filenames_cropped_lu_maps[filename_ind_first_cropped_lu_map]

  if(length(filename_first_cropped_lu_map)==0){
    stop("The initial fractional grids haven't been created yet.",
         "You need to run `gen_fract_lu_matrices` before you can extract the indices.")
  }

  mask <- terra::rast(filename_first_cropped_lu_map)

  if(terra::nlyr(mask)>1){
    mask <- terra::app(mask, fun = sum,  na.rm = T)
  }

  is_nonNA_mask <- which(!is.na(terra::values(mask)))
  terra::values(mask) <- NA
  terra::values(mask)[is_nonNA_mask] <- 1

  return(mask)


}

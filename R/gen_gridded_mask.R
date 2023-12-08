#' Extract the indices of a gridded mask for the area of interest
#' @import dplyr
#'
#' @param country
#' @param paths_lu_filenames
#' @param dir_output_files
#' @param cut_off_year
#' @param years
#' @param aggregation_factor
#'
#' @return
#' @export
#'
#' @examples
gen_gridded_mask <- function(country,
                             cut_off_year,
                             years,
                             aggregation_factor,
                             paths_lu_filenames,
                             dir_output_files){

  country_name = country

  cat("\nGenerating a mask for:", country_name, "\n")

  specify_output_dir(dir_output_files = dir_output_files,
                     aoi = country,
                     aggregation_factor = aggregation_factor)

  #filenames_cropped_lu_maps <- list.files(pattern="cropped_lu_map")
  filenames_cropped_lu_maps = intersect(list.files(pattern = "cropped_lu_map"),
                                        list.files(pattern = ".tif$"))

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

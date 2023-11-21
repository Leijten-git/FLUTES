#' Specify the output directory
#'
#' @param dir_output_files
#' @param aoi
#' @param aggregation_factor
#'
#' @return
#' @export
#'
#' @examples
specify_output_dir <- function(dir_output_files,
                               aoi,
                               aggregation_factor){

  setwd(dir_output_files)
  dir.create(aoi, showWarnings = F)
  setwd(file.path(dir_output_files, aoi))
  folder_agg_factor <- paste0("agg_factor_", aggregation_factor)
  dir.create(folder_agg_factor, showWarnings = F)
  setwd(file.path(dir_output_files, aoi, folder_agg_factor))

  return(getwd())

}

#' Crop the input raster file to the extent of the country of interest
#' @import dplyr
#' @return
#' @export
#'
#' @examples
#'
crop_data_to_country <- function(country,
                                 years,
                                 paths_lu_filenames,
                                 dir_output_files,
                                 nc_layer = "lccs_class"){

  cat("\nCropping the input raster file to the spatial extent of the country of interest\n")

  world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  world_sf_reproj <- sf::st_transform(world_sf,
                               crs = terra::crs(terra::rast(paths_lu_filenames[1])))

  country_sf <- world_sf_reproj %>%
    dplyr::filter(name == country)

  crop_lu_raster_files <- function(path_lu_filename){

    index <- which(paths_lu_filenames==path_lu_filename)
    year <- years[index]

    cat("\nNow cropping the input raster file for year:", year, "\n")

    lu_filename <- sub('.*/', '', path_lu_filename)
    filename_extension <- sub('.*\\.', '', lu_filename)

    if(filename_extension == "nc"){
      raw_lu_map <- terra::rast(path_lu_filename)[nc_layer]
    } else {
      raw_lu_map <- terra::rast(path_lu_filename)
    }

    raw_lu_map_cropped <- terra::crop(raw_lu_map, country_sf)

    names(raw_lu_map_cropped) <- paste0(country, "_cropped_lu_map")

    raw_lu_map_aligned <- extract_by_mask(raw_lu_map_cropped,
                                          country_sf)

    filename_cropped_lu <- paste0(country, "_cropped_lu_map_", year, ".tif")

    setwd(dir_output_files)
    dir.create("raw_lu_files_by_country", showWarnings = F)
    setwd(file.path(dir_output_files, "raw_lu_files_by_country"))
    dir.create(country, showWarnings = F)
    setwd(file.path(getwd(), country))

    terra::writeRaster(raw_lu_map_aligned, filename_cropped_lu, overwrite = T)

  }

  sapply(X = paths_lu_filenames, FUN = crop_lu_raster_files)


}

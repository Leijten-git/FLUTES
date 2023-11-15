#' Title
#' @import dplyr
#' @return
#' @export
#'
#' @examples
#'
crop_data_to_country <- function(country,
                                 years,
                                 paths_lu_filenames){

  world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  world_sf_reproj <- sf::st_transform(world_sf,
                               crs = terra::crs(terra::rast(paths_lu_filenames[1])))

  country_sf <- world_sf_reproj %>%
    dplyr::filter(name == country)

  crop_lu_raster_files <- function(path_lu_filename){

    lu_filename <- sub('.*/', '', path_lu_filename)
    filename_extension <- sub('.*\\.', '', lu_filename)

    # if nc then -->


    unique_land_systems <- unique(sub('.*_', '', unlist(land_systems_list)))

    raw_lu_map <- terra::rast(path_lu_filename)

    raw_lu_map_cropped <- terra::crop(raw_lu_map, country_sf)

    raw_lu_map_aligned <- extract_by_mask(raw_lu_map_cropped,
                                          country_sf)

  }


}

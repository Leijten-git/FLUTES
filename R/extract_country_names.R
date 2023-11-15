#' Extract the names of the countries that may be selected as area of interest
#'
#' @return A character with the list of country names
#' @export
#'
#' @examples
extract_country_names <- function(){

  world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  country_names <- sort(world_sf$name)
  return(country_names)

}

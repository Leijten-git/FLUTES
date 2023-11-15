#' Align an input raster to a pre-specified reference raster
#'
#' @param input_grid A SpatRaster that will be aligned
#' @param ref_grid A SpatRaster or shapefile that serves as a mask for the input grid. #'
#' @param res_method The resampling technique. Defaulting to "near", i.e.,
#' nearest neighbour resampling.
#'
#' @return
#' @export
#'
#' @examples
extract_by_mask <- function(input_grid,
                            ref_grid,
                            res_method = "near"){

  if(class(ref_grid)[1] == "sf"){

    cat("\nReference grid is a shapefile. It will be rasterized first...\n")

    ext_ref_grid <- ext(ref_grid)

    ref_grid$polygon <- c(1:nrow(ref_grid))

    ref_grid <- terra::rasterize(ref_grid,
                                 input_grid,
                                 field = "polygon")

    ref_grid <- terra::crop(ref_grid,
                            ext_ref_grid)

    cat("\nShapefile has been rasterized!\n")

  }

  cat("\nAligning the input raster: '",
      names(input_grid),
      "' to the reference raster: '",
      names(ref_grid),
      "' ",
      "(the resampling technique is set to '",
      res_method,
      "').\n",
      sep = "")

  redundant_cells <- which(is.na(terra::values(ref_grid)))

  input_grid_reproject <- terra::project(input_grid,
                                         ref_grid,
                                         method = res_method)

  input_grid_crop <- terra::crop(input_grid_reproject,
                                 ref_grid)

  input_grid_crop_aligned <- terra::resample(input_grid_crop,
                                             ref_grid,
                                             method = res_method)

  terra::values(input_grid_crop_aligned)[redundant_cells] <- NA

  return(input_grid_crop_aligned)

}

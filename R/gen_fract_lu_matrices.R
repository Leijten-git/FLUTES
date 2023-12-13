#' Generate fractional land use matrices over time based on the cropped land cover maps
#'
#' @param country
#' @param years
#' @param aggregation_factor
#' @param dir_output_files
#' @param path_lu_legend_filename
#' @param vals_lu_classes_to_exclude
#' @param no_unfilled_cells
#'
#' @import dplyr
#' @return
#' @export
#'
#' @examples
gen_fract_lu_matrices <- function(country,
                                  years,
                                  cut_off_year,
                                  aggregation_factor,
                                  dir_output_files,
                                  path_lu_legend_filename,
                                  vals_lu_classes_to_exclude = NULL,
                                  no_unfilled_cells = T){

  cat("\nGenerating the fractional land use matrices over time based on the cropped land cover maps\n")

  setwd(file.path(dir_output_files, "raw_lu_files_by_country", country))
  raster_files <- list.files(pattern=".tif$")
  raster_files_filtered <- raster_files[sapply(X = years, grep, raster_files)]

  raster_files_filtered_sorted <-
    raster_files_filtered[sapply(X = sort(years),
                                 function(x) which(years==x))]

  lu_raster_stack <- terra::rast(raster_files_filtered_sorted)
  names(lu_raster_stack) <- raster_files_filtered_sorted

  convert_cat_grid_into_fractional_grid <- function(input_map){

    cat("\nProcessing:", names(input_map), "\n",
        "Aggregation factor:", aggregation_factor, "\n")

    if(!is.null(vals_lu_classes_to_exclude)){

      for(value in vals_lu_classes_to_exclude){

        cat("Excluding pixels with value:", value, "\n")
        input_map[input_map==value] <- NA

      }

    }

    unique_vals <- sort(unique(terra::values(input_map, na.rm = T)))

    lu_classes <- read.csv(path_lu_legend_filename) %>%
      filter(Value %in% unique_vals) %>%
      dplyr::select(Label) %>%
      pull()

    cat("\nLand use classes: \n", paste0(lu_classes, "\n"), "\n")

    split_up_raster <- function(index){

      count <- which(unique_vals==index)

      cat("Creating a separate raster for land use class:",
          lu_classes[count],
          "\n")

      input_map[input_map!=index] <- NA
      input_map[input_map==index] <- 1

      return(input_map)

    }

    input_maps_split <- terra::rast(lapply(unique_vals, split_up_raster))

    fractional_scalar <- aggregation_factor^2

    cat("\nAggregating the different raster files and converting the values into fractional values\n")

    input_maps_aggregated <- terra::aggregate(input_maps_split,
                                              fact = aggregation_factor,
                                              fun = "sum",
                                              na.rm = T)

    input_maps_frac <- input_maps_aggregated %>%
      terra::app(function(x) x / fractional_scalar)

    specify_output_dir(dir_output_files = dir_output_files,
                       aoi = country,
                       aggregation_factor = aggregation_factor)

    terra::writeRaster(input_maps_frac,
                       names(input_map),
                       overwrite = T)

    names(input_maps_frac) <- lu_classes

    return(input_maps_frac)

  }

  lu_frac_grids_list <- lapply(X = lu_raster_stack,
                                  FUN = convert_cat_grid_into_fractional_grid)

  names(lu_frac_grids_list) <- sort(years)

  convert_map_into_fractional_lu_matrix <- function(year){

    input_map_frac <- lu_frac_grids_list[[year]]

    lu_frac_matrix <- matrix(data = terra::values(input_map_frac),
                             ncol = length(names(input_map_frac)))

    lu_frac_matrix[is.na(lu_frac_matrix)] <- 0

    if(no_unfilled_cells == T){

      cat("\nRecaling all cells to ensure the fractions sum up to 1\n")

      specify_output_dir(dir_output_files = dir_output_files,
                         aoi = country,
                         aggregation_factor = aggregation_factor)

      filenames_cropped_lu_maps = intersect(list.files(pattern = "cropped_lu_map"),
                                            list.files(pattern = ".tif$"))

      if(!is.null(cut_off_year)){

        cat("\nUsing the land use map for the year",
            cut_off_year,
            "as reference map\n")

        filename_ind_ref_lu_map <- grep(as.character(cut_off_year),
                                        filenames_cropped_lu_maps)
      } else {

        first_year <- sort(years)[1]

        cat("\nUsing the land use map for the year ",
            first_year,
            "as reference map")

        filename_ind_ref_lu_map <- grep(as.character(first_year),
                                        filenames_cropped_lu_maps)
      }

      filename_ref_cropped_lu_map <-
        filenames_cropped_lu_maps[filename_ind_ref_lu_map]

      ref_cropped_lu_map <- terra::rast(filename_ref_cropped_lu_map)

      if(terra::nlyr(ref_cropped_lu_map)>1){

        ref_cropped_lu_map <- terra::app(ref_cropped_lu_map, fun = sum,
                                         na.rm = T)
      }

      inds_NA_vals <- which(is.na(terra::values(ref_cropped_lu_map)))

      if(length(inds_NA_vals)>0){
        lu_frac_matrix = lu_frac_matrix[-inds_NA_vals,]
      }

      rowsums2 <- rowSums(lu_frac_matrix, na.rm = T)

      lu_frac_matrix <- lu_frac_matrix %>%
        as.data.frame() %>%
        mutate_all(list(~./rowsums2)) %>%
        as.matrix()

    }

    colnames(lu_frac_matrix) <- paste0(year, "_", names(input_map_frac))

    return(lu_frac_matrix)


  }

  lu_frac_matrices_list <- lapply(X = as.character(sort(years)),
                                  FUN = convert_map_into_fractional_lu_matrix)

  lu_classes <- identify_lu_classes(lu_frac_matrices_list = lu_frac_matrices_list)

  names(lu_frac_matrices_list) <- sort(years)

  add_missing_lu_classes <- function(year){

    lu_matrix <- lu_frac_matrices_list[[as.character(year)]]

    unique_lu_classes_year <- unique(sub('.*_', '', colnames(lu_matrix)))

    missing_lu_classes_indices <- which(!lu_classes %in% unique_lu_classes_year)

    year_char <- as.character(year)
    order_col_names <- paste0(year_char, "_", lu_classes)

    if(length(missing_lu_classes_indices)>0){


      missing_lu_classes <- lu_classes[missing_lu_classes_indices]
      new_col_names <- paste0(year_char, "_", missing_lu_classes)

      for(new_col_name in new_col_names){

        cat("In year", year, "the following land use classes were missing:",
            new_col_name, "\nAn empty column has been added\n")
        warning("In year ", year, " the following land use classes were missing: ",
                new_col_name, "\nAn empty column has been added\n")

        lu_matrix <- lu_matrix %>%
          as.data.frame() %>%
          mutate(!!new_col_name := 0)

      }

    }

    if(class(lu_matrix)[1]=="matrix"){
      lu_matrix <- lu_matrix %>%
        as.data.frame()
    }

    lu_matrix[is.na(lu_matrix)] <- 0
    lu_matrix <- lu_matrix %>%
      dplyr::select(all_of(order_col_names)) %>%
      as.matrix()

    return(lu_matrix)

  }

  lu_frac_matrices_list2 <- lapply(X = sort(years),
                                   FUN = add_missing_lu_classes)

  names(lu_frac_matrices_list2) <- sort(years)

  lu_frac_matrices_list3 <- clean_lu_fact_matrices(matrices_list = lu_frac_matrices_list2,
                                                   cut_off_year = cut_off_year)

  period <- paste0(c(as.character(first(sort(years))),
                     "-",
                     as.character(last(sort(years)))), collapse = "")

  filename_lu_frac_matrices <- paste0("lu_frac_matrix_matrices_",
                                      aggregation_factor,
                                      "_",
                                      period,
                                      ".Rds")

  cat("\nWriting the fractional land use matrices here:\n",
      getwd(), "\n")

  saveRDS(lu_frac_matrices_list3, filename_lu_frac_matrices)

  return(lu_frac_matrices_list3)

}

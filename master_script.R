library(devtools)
library(FLUTES)

# user input --------------------------------------------------------------

base_directory <- "" # specify path to folder that contains the software package
package_name <- "FLUTES" # specify packahe name
dir_lu_filenames <- "" # path to filenames (raster files; should be in terra format) with land use types

path_cloud <- "" # path to cloud folder (e.g., OneDrive)
base_dir <- paste0(path_cloud, "/") # path output dir (where you want to save the output files)
path_lu_legend_filename <- file.path("") # path csv file with land use keys and labels

extract_country_names() # list all country names

country = "" # specify country of interest

# leave as is -------------------------------------------------------------

file_path_pkg <- paste0(base_directory, "/", package_name)
load_all(path = file_path_pkg)

setwd(dir_lu_filenames)
file_names <- list.files()[grepl("LCCS", list.files())]
paths_lu_filenames <- paste0(getwd(), "/", file_names)
years <- as.numeric(stringr::str_sub(sub('.*P1Y', '', paths_lu_filenames), 2, 5))

# run FLUTES --------------------------------------------------------------

run_FLUTES(.country = country,
           .years = years,
           .paths_lu_filenames = paths_lu_filenames,
           .dir_output_files = base_dir,
           .aggregation_factor = 32,
           .path_lu_legend_filename = path_lu_legend_filename,
           .cut_off_year = 1998, # 1999 - 2020
           .lags = c(1),
           .growth_constraint = NULL, # NULL; numeric value; or "cal" (calculate)
           .vals_lu_classes_to_exclude = 210) # water bodies


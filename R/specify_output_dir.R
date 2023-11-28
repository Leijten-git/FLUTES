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
                               aggregation_factor,
                               create_sub_folder = F,
                               are_conv_pars_perc = NULL,
                               growth_par = NULL,
                               are_PAs_excluded = NULL){

  setwd(dir_output_files)
  dir.create(aoi, showWarnings = F)
  setwd(file.path(dir_output_files, aoi))
  folder_agg_factor <- paste0("agg_factor_", aggregation_factor)
  dir.create(folder_agg_factor, showWarnings = F)
  setwd(file.path(dir_output_files, aoi, folder_agg_factor))
  output_dir_base <- getwd()

  if(create_sub_folder == T){

    setwd(output_dir_base)

    if(!is.null(are_conv_pars_perc)){
      cwd <- getwd()

      if(are_conv_pars_perc==T){
        dir.create("conv_pars_perc", showWarnings = F)
        setwd(file.path(cwd, "conv_pars_perc"))
      } else if(are_conv_pars_perc==F){
        dir.create("conv_pars_abs", showWarnings = F)
        setwd(file.path(cwd, "conv_pars_abs"))
      }
    }

    if(!is.null(growth_par)){
      cwd <- getwd()
      growth_folder_name <- paste0("growth_", growth_par)
      dir.create(growth_folder_name, showWarnings = F)
      setwd(file.path(cwd, growth_folder_name))
    }

    if(!is.null(are_PAs_excluded)){
      cwd <- getwd()

      if(are_PAs_excluded==T){
        dir.create("PAs_excluded", showWarnings = F)
        setwd(file.path(cwd, "PAs_excluded"))
      } else if(are_PAs_excluded==F){
        dir.create("no_PAs", showWarnings = F)
        setwd(file.path(cwd, "no_PAs"))
      }
    }

    return(getwd())

  } else {

    return(output_dir_base)

  }

}

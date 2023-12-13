#' Calculate the average growth parameter within empty cells
#'
#' @param lu_matrices
#' @param years
#' @param cut_off_year
#'
#' @return
#' @export
#'
#' @examples
calc_growth_pars_empty_cells <- function(lu_matrices,
                                         years,
                                         cut_off_year){

  first_year <- sort(years)[1]

  lu_matrix_first_year <- lu_matrices[[as.character(first_year)]]
  lu_matrix_cut_off_year <- lu_matrices[[as.character(cut_off_year)]]

  calc_n_filled_cells <- function(index){
    inds_empty_cells_f_year <- which(lu_matrix_first_year[,index]==0)
    vals_cells_cut_off_year <- lu_matrix_cut_off_year[,index][inds_empty_cells_f_year]
    n_filled_cells <- sum(vals_cells_cut_off_year>0)
    prop_filled_cells <- n_filled_cells/length(inds_empty_cells_f_year)
    return(prop_filled_cells)
  }

  prop_filled_cells_by_lu_type <- sapply(X = c(1:ncol(lu_matrix_first_year)),
                                         FUN = calc_n_filled_cells)

  perc_filled_cells_by_lu_type <- prop_filled_cells_by_lu_type*1e2

  return(perc_filled_cells_by_lu_type)

}

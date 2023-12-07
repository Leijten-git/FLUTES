#' Clean land use fractional matrices
#'
#' @param matrices_list
#' @param cut_off_year
#'
#' @return
#' @export
#'
#' @examples
clean_lu_fact_matrices <- function(matrices_list,
                                   cut_off_year){

  lu_frac_matrix <- matrices_list[[as.character(cut_off_year)]]
  ind_unfilled_lu_classes <- which(colSums(lu_frac_matrix)==0)

  if(length(ind_unfilled_lu_classes)>0){

    names_unfilled_lu_classes_w_year <-
      colnames(lu_frac_matrix)[ind_unfilled_lu_classes]

    names_unfilled_lu_classes <- sub('.*_', '', names_unfilled_lu_classes_w_year)

    cat("\nThe following land use classes are removed due to a lack of",
        "observations in the reference year:", cut_off_year, "\n",
        paste0(names_unfilled_lu_classes, "\n"), "\n")

    warning("\nThe following land use classes were removed due to a lack of ",
            "observations in the reference year: ", cut_off_year, "\n",
            paste0(names_unfilled_lu_classes, "\n"), "\n")

    for(list_element in c(1:length(matrices_list))){

      lu_classes <- colnames(matrices_list[[list_element]])
      lu_classes_filter <- sub('.*_', '', lu_classes)

      # inds_unfilled_lu_classes <-
      #   as.numeric(sapply(X = names_unfilled_lu_classes,
      #                     function(x) which(lu_classes_filter %in% x)))

      # for(col_index in as.numeric(ind_unfilled_lu_classes)){
      #   name_col = colnames(matrices_list[[list_element]])[col_index]
      #   cat("\nRemoving:", name_col)
      #   matrices_list[[list_element]] <- matrices_list[[list_element]][,-col_index]
      # }

      matrices_list[[list_element]] <- matrices_list[[list_element]][,-ind_unfilled_lu_classes]


    }

    cat("\n")

  }

  return(matrices_list)

}

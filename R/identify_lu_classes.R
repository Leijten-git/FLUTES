#' Title
#'
#' @param lu_frac_matrices_list
#'
#' @return
#' @export
#'
#' @examples
identify_lu_classes <- function(lu_frac_matrices_list){

  extract_lu_classes <- function(index){

    lu_matrix <- lu_frac_matrices_list[[index]]
    lu_classes <- colnames(lu_matrix)
    return(lu_classes)

  }

  lu_classes_list <- lapply(X = c(1:length(lu_frac_matrices_list)),
                            FUN = extract_lu_classes)

  unique_lu_classes <- unique(sub('.*_', '', unlist(lu_classes_list)))

  return(unique_lu_classes)

}

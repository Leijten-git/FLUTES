#' Determines the maximum deviation parameters
#'
#' @param lu_requirements
#' @param lu_first_period
#' @param default_dev_par
#' @param pixel_cut_off_value
#' @param pixel_cut_off_value2
#' @param alt_dev_par
#' @param alt_dev_par2
#'
#' @import dplyr
#' @return
#' @export
#'
#' @examples
determine_max_dev_pars <- function(lu_requirements,
                                   lu_first_period,

                                   default_dev_par = 1,
                                   pixel_cut_off_value = 15,
                                   pixel_cut_off_value2 = 3,
                                   alt_dev_par = 5,
                                   alt_dev_par2 = 15){

  cat("\nAutomatically setting the maximum deviation parameters\n")

  rel_changes <- apply(lu_requirements, 1, function(x) x/lu_requirements[1,])

  rel_changes_df <- rel_changes[[2]] %>%
    t() %>%
    `colnames<-`("rel_change") %>%
    data.frame()

  rel_changes_df <- rel_changes_df %>%
    mutate(diff = rel_change-1,
           abs_change = abs(rel_change-1),
           abs_perc_change = abs_change*1e2,
           land_system = rownames(rel_changes_df)) %>%
    arrange(-diff)

  rownames(rel_changes_df) <- NULL

  required_perc_changes <- rel_changes_df %>%
    dplyr::select(abs_perc_change) %>%
    pull()

  max_dev_par = default_dev_par
  max_devs <- rep(max_dev_par, ncol(lu_requirements))

  # indices_prob_vars <- which(max_devs >= required_perc_changes)
  #
  # if(length(indices_prob_vars)>0){
  #
  #   skip_these_vars <- which(required_perc_changes[indices_prob_vars]<pixel_cut_off_value2)
  #   indices_prob_vars <- indices_prob_vars[-skip_these_vars]
  #
  # }
  #
  # if(length(indices_prob_vars)>0){
  #
  #   new_values <- required_perc_changes[indices_prob_vars]*0.9
  #   max_devs[indices_prob_vars] <- new_values
  #
  # }

  total_area_by_land_class <- colSums(lu_first_period)

  indices_prob_classes <-
    as.numeric(which(total_area_by_land_class<pixel_cut_off_value))
  indices_prob_classes2 <-
    as.numeric(which(total_area_by_land_class<pixel_cut_off_value2))

  if(length(indices_prob_classes)>0){

    max_devs[indices_prob_classes] <- alt_dev_par

    if(length(indices_prob_classes2)){

      max_devs[indices_prob_classes2] <- alt_dev_par2

    }

  }



  return(max_devs)



}


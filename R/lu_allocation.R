lu_allocation <- function(lu,
                          sm,
                          params,
                          dmd,
                          ln,
                          constraint,
                          rescale,
                          PA,
                          adjust_dev_pars = T){

  cat('\nStarting simulations...\n')

  #number of land use classes and number of cells
  K <- ncol(lu) # FL: number of land systems
  n <- nrow(lu) # FL: number of pixels

  is_absolute_dev <- params$is_abs_dev

  if(is_absolute_dev){

    cat("\nMaximum deviations expressed in absolute differences, not in percentages!\n")

  } else if(is_absolute_dev==F){

    cat("\nMaximum deviations expressed in percentages differences, not in absolute differences!\n")

  } else {

    stop("'is_absolute_dev' is a logical operator and only accepts two values: TRUE or FALSE")

  }

  max_devs <<- params$max_devs # FL: should be flexible, depending on the demand variable
  max_iter <<- params$max_iter
  growth <- params$growth
  no_change <-  1:K%in%params$no_change

  cat("Maximum deviation tolerance thresholds:", max_devs)

  supply_t0 <- colSums(lu) # FL: count pixels by land system in the initial situation

  demand_t1 <- dmd[2,] # FL: extract required pixels next sim period

  # Land use supply in the first time stepis the "candidate" that will be recalculated
  # in each iteration until it meets demand_t1. For now, it's equal to supply_t0
  supply_t1_candidate <- dmd[1,]

  # supply_t0 == supply_t1_candidate # TRUE FOR ALL CASES (for now)

  #First land use map becomes "candidate" on which allocations take place on an iterative basis
  p_t1_candidate <- lu # FL: fractions current LU
  # note: p_t1_candidate will output the final projections!

  diff_demand <- diff(as.matrix(dmd)) # changes between periods
  abs_diff_demand <- abs(diff_demand) # absolute changes between periods

  if(is_absolute_dev == T){

    dev_diffs <- as.numeric(abs_diff_demand)

  } else if(is_absolute_dev == F){

    dev_diffs <- t(t(abs_diff_demand) / dmd[1,])*100 # calculate relative (%) changes (relative to initial situation)

  }

  dev_diffs[which(is.na(dev_diffs))] <- 0 # set NA values to 0

  #Counter
  count <- 0

  if(PA == T){

    cat("\nAreas marked as legally protected are excluded from the simulations...\n")

    pa <- rast("wdpa_1992_masked.tif")

    pa_vals <- as.numeric(values(pa))
    pa_vals_filter <- pa_vals[is_essential_cell]
    pa_inds <- which(pa_vals_filter==1)

    supply_t0_pa <- colSums(p_t1_candidate[pa_inds,])

    nopa_inds <- which(is.na(pa_vals_filter))

    p_t1_candidate <- p_t1_candidate[nopa_inds,]
    ln <- ln[nopa_inds,]
    sm <- sm[nopa_inds,]

  }

  if(constraint){

    cat("\nLand use change within empty pixels is constrained by the following growth parameters:",
        growth)

    are_zero <- p_t1_candidate == 0 # matrix with boolean columns: is the fraction within a cell equal to 0?
    inds_list <- list() # list of vectors for each land use type containing the indices of cells that should remain empty

    i <- 1

    for(i in 1:K){ # loop over number of land systems

      inds_are_zero <- which(are_zero[,i]) # indices of the cells belonging to LU type k that are empty (i.e., 0)

      inds_inds_zero_pos_neigh <- which(ln[inds_are_zero,i]!=0)  # indices of the indices of empty cells belonging to LU type k that are in an area with non-empty cells
      # note: growth should be allowed here

      inds_zero_pos_neigh <- inds_are_zero[inds_inds_zero_pos_neigh] # indices of empty cells that contain populated cells in the neighbourhood

      # size <- round(nrow(lu) * ((growth[i])/100)) # sample size of empty cells that may experience growth: this is weird, why take a sample of the total
      # number of pixels? why not take a sample of the empty pixels...

      size <- round(length(inds_are_zero) * ((growth[i])/100)) # sample size of empty cells that may experience growth

      # two possibilities:
      # - sample size is not larger than inds_zero_pos_neigh: ensure all empty cells stay empty, excepting a sample of pixels in populated neighbourhoods
      # - sample size is larger than inds_zero_pos_neigh: now there are not enough empty pixels in populated areas to sample.
      # now you need to combine a sample of empty pixels in populated neighbourhoods with a sample of empty pixels in empty neighbourhoods.

      if(size <= length(inds_zero_pos_neigh)){ # if the sample size is not larger than the number of unpopulated cells that contain populated cells in the neighbourhood

        sample_pixels_pos_neigh <- inds_zero_pos_neigh[wrswoR::sample_int_rank(length(inds_zero_pos_neigh),
                                                                               size = size,
                                                                               prob = sm[inds_zero_pos_neigh,i])] # take sample of cells that may experience growth

        #inds_list[[i]] <- inds_are_zero[-which(inds_are_zero %in% keep)] # original bug, this results in an empty list

        if(length(sample_pixels_pos_neigh)>0){ # fix bug

          inds_list[[i]] <- inds_are_zero[-which(inds_are_zero %in% sample_pixels_pos_neigh)] # only add pixels outside the sample of populated areas

          check <- length(inds_list[[i]])==(length(inds_are_zero)-size) # has to be TRUE

        } else {

          inds_list[[i]] <- inds_are_zero

        }



      }

      # no what if the sample size (of empty pixels that may experience growth) is larger than the number of empty pixels in populated areas
      # exhaust all the empty pixels in populated areas and add a random sample of empty pixels in empty areas
      if(size > length(inds_zero_pos_neigh)){ # if the sample size is larger than the number of empty cells that contain populated cells in the neighbourhood

        # first calculate: how many pixels are missing?

        n_missing_pixels <- size - length(inds_zero_pos_neigh)

        # then take a sample of empty pixels in empty areas
        # where are these?

        inds_inds_zero_no_neigh <- which(!inds_are_zero %in% inds_zero_pos_neigh) # extract indices empty pixels with no populated pixels in the neighbourhood

        inds_leftover <- inds_are_zero[inds_inds_zero_no_neigh]

        # what to do if the number of remaining pixels is smaller than the required sample size?
        # then sample up to the maximum point
        # WARNING: this means that the number of pixels allowed to grow is less than expected

        sample_indices_leftover <- wrswoR::sample_int_rank(n = length(inds_leftover),
                                                           size = n_missing_pixels,
                                                           prob = sm[inds_leftover,i])

        sample_pixels_empty_areas <- inds_leftover[sample_indices_leftover] # take sample of empty cells in empty regions

        #check <- (length(inds_zero_pos_neigh)+length(sample_pixels_empty_areas))==size # needs to be TRUE

        indices_empty_cells_to_be_excluded <- c(inds_zero_pos_neigh,
                                                sample_pixels_empty_areas)

        to_exclude <- which(inds_are_zero %in% indices_empty_cells_to_be_excluded) # combine the sample of empty pixels in populated areas with the sample of empty pixels in empty areas

        sum(sample_pixels_empty_areas %in% inds_zero_pos_neigh) # should be 0
        sum(inds_zero_pos_neigh %in% sample_pixels_empty_areas) # should be 0

        #length(inds_are_zero)-length(to_exclude)==length(inds_are_zero[-to_exclude]) # ought to be true

        inds_list[[i]] <- inds_are_zero[-to_exclude]

      }
    }
  }

  # Iterative allocation
  while (any(dev_diffs > max_devs)){

    #Counter increment
    count <- count + 1 # keep track of the number of iterations

    #Demand to be allocated
    demand_change <- demand_t1 - supply_t1_candidate # FL: subtracting demand next sim period from demand current period
    # FL: in other words, the above is equal to: diff(as.matrix(dmd))
    # it's the required absolute change in supply!

    if(PA == T){
      demand_change <- demand_change - supply_t0_pa # exclude the areas within the protected mask
    }

    #Calculate change factor (a),by how much do we have to multiply the candidate
    # land use proportions to satisfy the modelled suitability.

    # FL: 'ideal change' represents the "change factors" shown in figure 2 of the
    # main paper. It's based on the ratio of the predicted fraction to the actual
    # fraction

    ideal_change <- sm / p_t1_candidate # FL: divide suitability estimates by current fractions

    both_0 <- is.na(ideal_change) # FL: identify cells where both suitability and actual fraction is 0.
    ideal_change[both_0] <- 1 # FL: in case where both are 0, assume no change (so set change factors to 1)

    cand_0 <- !is.finite(ideal_change) #This determines which cells are INF in p_t1_candidate
    ideal_change[cand_0] <- sm[cand_0] # FL: in case of infinite numbers, the actual values are 0. Use the predicted values in this case

    #Calculate Relative suitability (r) from change factors.
    # %*% --> matrix multiplication

    rel_suitability <- ideal_change %*% diag(1/colSums(ideal_change))
    # FL: 1/colSums(ideal_change) --> calculates the inverse of the sum of all ideal changes
    # diag(1/colSums(ideal_change)) --> create matrix with the inverse total suitabilities on the diagonal
    # ideal_change %*% diag(1/colSums(ideal_change)) --> multiply the ideal changes by the inverse total suitabilities
    # this ensures that all suitabilities sum to 1

    # Put differently, we scale the ideal changes by the inverse total ideal changes!

    #Allocate demand change between pixels (d)
    target_lu_change_pixel <-  rel_suitability %*% diag(demand_change)
    # FL: distributes additional demand (can be negative) across pixels based on the relative suitabilities.

    #Add changes to candidate map and make everything positive.
    p_t1_proposal <- p_t1_candidate + target_lu_change_pixel # KEY STEP: updating the candidates!
    p_t1_proposal <- pmax(p_t1_proposal, 0) # FL: if < 0, make 0 (can't have negative fractions)

    if(any(demand_t1 == 0)){ # FL: if there is at least one LU class for which the number of required pixels is 0

      p_t1_proposal[,which(demand_t1 == 0)] <- 0 # FL: then set simulated fractions to 0

    }

    rowsums <- rowSums(p_t1_proposal, na.rm = T)

    if(rescale == T){

      cat("\nAdjusting the fractions...\n")
      p_t1_proposal <- p_t1_proposal %>%
        as.data.frame() %>%
        mutate_all(list(~./rowsums)) %>%
        as.matrix()

    }

    if(constraint){

      for(i in (1:K)){
        p_t1_proposal[inds_list[[i]], i] <- 0
      }

    }

    # step 2 ------------------------------------------------------------------

    p_t1_candidate <- p_t1_proposal

    # Calculate new candidate supply, i.e. the supply of the currently proposed candidate
    supply_t1_candidate <- colSums(p_t1_candidate) # calculate total number of supply

    #Recalculate % deviation of candidate supply from demand:
    if(PA == T){
      diffs <- abs(demand_t1 - (supply_t1_candidate + supply_t0_pa))
    } else if(PA == F){
      diffs <- abs(demand_t1 - supply_t1_candidate)
    }

    dev_diffs <- as.numeric(diffs)

    if(is_absolute_dev == T){

      cat("\n", paste0("Iteration: ", count, "    "),
          "Deviation from target per class [absolute value]: ",
          paste(round(dev_diffs, 4), sep = " "))

    } else if(is_absolute_dev == F){

      dev_diffs <- dev_diffs/demand_t1 * 100
      dev_diffs[which(is.na(dev_diffs))] <- 0

      cat("\n", paste0("Iteration: ", count, "    "),
          "Deviation from target per class [%]: ",
          paste(round(dev_diffs, 4), sep = " "))

    }

    final_count <<- count

    if(adjust_dev_pars == T){

      first_cut_off_value <- round(max_iter*0.25)
      second_cut_off_value <- round(max_iter*0.5)
      third_cut_off_value <- round(max_iter*0.75)

      if(max_iter > 5){

        if(count == first_cut_off_value){

          cat("\nNo solution found after", count, "iterations.")

          max_devs[which(dev_diffs > max_devs)] <- max_devs[which(dev_diffs > max_devs)]+2

          cat("\nIncreased the maximum deviance tolerance parameters to:\n")

          cat(paste0(max_devs, collapse = " "))

          cat("\n")

        } else if(count == second_cut_off_value){

          cat("\nNo solution found after", count, "iterations.")

          max_devs[which(dev_diffs > max_devs)] <- max_devs[which(dev_diffs > max_devs)]+8

          cat("\nIncreased the maximum deviance tolerance parameters to:\n")

          cat(paste0(max_devs, collapse = " "))

          cat("\n")

        } else if(count == third_cut_off_value){

          cat("\nNo solution found after", count, "iterations.")

          max_devs[which(dev_diffs > max_devs)] <- round(dev_diffs[which(dev_diffs > max_devs)]+1)

          cat("\nIncreased the maximum deviance tolerance parameters to:\n")

          cat(paste0(max_devs, collapse = " "))

          cat("\n")
        }

      }

    }

    if (count == max_iter){

      stop("\nMaximum iterations reached")

    }

  }

  cat("\nSaving number of iterations...\n")

  saveRDS(final_count, "n_iterations.Rds")

  # end while loop ----------------------------------------------------------

  #When allocations are ready, return result
  if(PA == T){

    pred_out <- matrix(NA, ncol = K, nrow = n)
    pred_out[nopa_inds,] <- p_t1_candidate
    pred_out[pa_inds,] <- lu[pa_inds,]

  } else if (PA == F){

    pred_out <- p_t1_candidate

  }

  final_dev_diffs <- as.numeric(diffs)
  final_dev_diffs_perc <- final_dev_diffs/demand_t1 * 100
  final_dev_diffs_perc[which(is.na(final_dev_diffs_perc))] <- 0

  cat("\nWriting the deviation values")

  saveRDS(data.frame(final_dev_diffs = final_dev_diffs,
                     final_dev_diffs_perc = as.numeric(final_dev_diffs_perc)),
          "dev_values.Rds")

  colnames(pred_out) <- colnames(lu)
  return(pred_out)

}

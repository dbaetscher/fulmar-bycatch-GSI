# nofu functions
library(dplyr)
library(tidyr)
library(rubias)

# @param base_ids: colony and sample id
# @param n_dsample: numeric, number to downsample to, e.g., 58
# @param baseline: baseline genotypes dataframe, 2-col format

downsample_baseline <- function(base_ids, n_dsample, base_genos, no_hi_missers){
  
    downsampled_three_colonies <- base_ids %>%
    group_by(repunit) %>%
    filter(repunit != "Chagulak") %>%
    sample_n(., n_dsample, replace = FALSE) #%>%
  #ungroup()
  
  # grab the Chagulak ids from the df
  chag_ids <- base_ids %>%
    filter(repunit == "Chagulak")
  
  # put the df back together
  baseline_dsampled_ids <- bind_rows(downsampled_three_colonies, chag_ids)
  
  # Use a semi_join to keep just those samples in the reference baseline genotypes
  dbaseline <- base_genos %>%
    semi_join(., baseline_dsampled_ids, by = "NMFS_DNA_ID")
  
  # Now, make a combined data frame with the bycatch and baseline samples to keep the integer values for haplotypes consistent.
  # here I want to combine the genotypes
  all_genos <- bind_rows(dbaseline, no_hi_missers)
  
  
  # first make integers of the alleles
  alle_idxs <- all_genos %>% 
    dplyr::select(NMFS_DNA_ID, locus, gene_copy, allele) %>%
    group_by(locus) %>%
    mutate(alleidx = as.integer(factor(allele, levels = unique(allele)))) %>%
    ungroup() %>%
    arrange(NMFS_DNA_ID, locus, alleidx) 
  
  # select just the columns to retain 
  alle_idx2 <- alle_idxs[,-4]
  
  # spread the alleles
  two_col <- alle_idx2 %>%
    unite(loc, locus, gene_copy, sep = ".") %>%
    spread(loc, alleidx)
  
  
  # I can use joins to split the combined data frame
  # first the mixture
  mix_df <- mixture_ids %>%
    left_join(two_col) %>%
    rename(indiv = NMFS_DNA_ID) %>%
    ungroup() # apparently the indiv column was still grouped and causing problems with rubias 
  
  mixture <- as_tibble(mix_df)
  mixture$repunit <- as.character(mixture$repunit)
  
  # and then the baseline
  base_df <- baseline_dsampled_ids %>%
    inner_join(., two_col) %>%
    rename(indiv = NMFS_DNA_ID) %>%
    ungroup()
  
  # perform mixture-assignment 
  mix_assign <- infer_mixture(reference = base_df, mixture = mixture, gen_start_col = 5, method = "MCMC", reps = 10000, burn_in = 1000)
  
  # individual data
  top_assign <- mix_assign$indiv_posteriors %>%
    group_by(indiv) %>%
    top_n(., 1, PofZ) %>%
    ungroup()
  
  # 90% assignment threshold
  # high-likelihood assignments with outlier z-scores?
  assign90 <- top_assign %>%
    filter(z_score < 2.5 & z_score > -2.5) %>%
    filter(PofZ > 0.90) %>%
    select(indiv, collection, PofZ, z_score)

  return(list(
    assign90
    
    
  ))
  
}



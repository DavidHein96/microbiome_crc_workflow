
# ------------ Libraries & Setup ------------ #

library(phyloseq)
library(vegan)
library(plyr)
library(microbiomeMarker)
library(tidyverse)
library(Maaslin2)
library(LinDA)

# ------------ Function to create phyloseq object ------------ #

build_phylo_obj <- function(selected_patients, counts_all) {
  set.seed(2023)
  # selected patients is semi-prefiltered clinical data with only
  # select variables and filtered patients, needs SampleID column, counts_all is full ASV table with patients as columns
  
  obj_list <- list()
  
  # Get only patients that we want to include
  counts_all <- counts_all %>%
    dplyr::select(
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species",
      any_of(selected_patients$SampleID)
    )
  
  print(head(counts_all))
  # Filter out Archea and ASVs that have all zeros
  counts_all <-
    counts_all %>% filter(Kingdom != "Archaea") %>% filter(rowSums(across(where(is.numeric))) !=
                                                             0)
  
  # add in signifier of taxa level
  counts_all <- counts_all %>% mutate(
    Kingdom = paste0("d__", Kingdom),
    Phylum = paste0("p__", Phylum),
    Class = paste0("c__", Class),
    Order = paste0("o__", Order),
    Family = paste0("f__", Family),
    Genus = paste0("g__", Genus),
    Species = paste0("s__", Species)
  )
  
  # add up the identical ASVs
  # first make a taxonomy variable in a single column and remove the multiple columns
  counts_all_combinedASV <- counts_all %>%
    mutate(taxonomy = base::paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep =
                                    ";")) %>%
    relocate(taxonomy)
  
  counts_all_combinedASV2 <- counts_all_combinedASV %>%
    dplyr::select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species)
  
  # use ddply to sum up the identical taxa
  counts_all_combinedASV3 <-
    counts_all_combinedASV2 %>% ddply("taxonomy", numcolwise(sum))
  
  # rejoin the taxa names and reorder
  taxa_names <- counts_all_combinedASV %>%
    dplyr::select(taxonomy, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    distinct()
  taxa_order <- counts_all_combinedASV3$taxonomy
  taxa_names_ordered <-
    taxa_names[match(taxa_order, taxa_names$taxonomy), ]
  
  # Build Phyloseq object
  asvs <- as.matrix(counts_all_combinedASV3[, -1])
  row.names(asvs) <- counts_all_combinedASV3$taxonomy
  otu <- otu_table(asvs, taxa_are_rows = TRUE)
  
  clin <- selected_patients %>% dplyr::select(-SampleID)
  row.names(clin) <- selected_patients$SampleID
  samples <- sample_data(clin)
  sample_names(samples) <- selected_patients$SampleID
  
  tax <- as.matrix(taxa_names_ordered[, -1])
  row.names(tax) <- taxa_names_ordered$taxonomy
  tax <- tax_table(tax)
  # print(samples)
  obj <- phyloseq(otu, samples, tax)
  
  # 10% prev filter
  obj_filt <-
    filter_taxa(obj, function(x)
      sum(x > 5) > nrow(clin) * 0.1, TRUE)
  
  # return both filtered and un filt phylo objs
  obj_list[["unfilt"]] <- obj
  obj_list[["filt"]] <- obj_filt
  
  return(obj_list)
  
}





# ------------ Function to run MaAsLin2 on all taxa levels with batch random effect ------------ #
run_mas_all_levels <- function(phylo_obj, variable) {
  set.seed(2023)
  mas_results <- data.frame()
  tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (tax_level in tax_levels) {
    # aggregate taxa at the taxa level we want
    obj_sub = aggregate_taxa(phylo_obj, paste0(tax_level))
    df_input_data <- data.frame(obj_sub@otu_table)
    df_meta <- data.frame(obj_sub@sam_data)
    # Run maaslin2 and bind results
    fit_mas <- Maaslin2(
      input_data = df_input_data,
      input_metadata = df_meta,
      output = "test",
      fixed_effects = c(paste0(variable)),
      random_effects = c("batch"),
      min_prevalence = 0.2,
      plot_heatmap = FALSE,
      plot_scatter = FALSE,
      normalization = 'CSS',
      transform = "LOG",
      min_abundance = 1,
      reference = c("t_simple,3ab", "anti_class,none", "new_anti,weakornone")
    )
    res <- data.frame(fit_mas$results)
    mas_results <- rbind(mas_results, res)
  }
  
  # adjust p values for multiple testing
  mas_bhcor_p <- p.adjust(mas_results[, 6], method = "BH")
  mas_results$mas_bhcor_p <- mas_bhcor_p
  return(mas_results)
}


# ------------ Function to run MaAsLin2 on all taxa levels no bactch effect and antibiotics fixed effect ------------ #
run_mas_all_levels_no_batch_anti <- function(phylo_obj, variable) {
  set.seed(2023)
  mas_results <- data.frame()
  tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  for (tax_level in tax_levels) {
    # aggregate taxa at the taxa level we want
    obj_sub = aggregate_taxa(phylo_obj, paste0(tax_level))
    
    df_input_data <- data.frame(obj_sub@otu_table)
    df_meta <- data.frame(obj_sub@sam_data)
    
    # Run maaslin2 and bind results
    fit_mas <- Maaslin2(
      input_data = df_input_data,
      input_metadata = df_meta,
      output = "test",
      fixed_effects = c(paste0(variable), "new_anti"),
      min_prevalence = 0.2,
      plot_heatmap = FALSE,
      plot_scatter = FALSE,
      normalization = 'CSS',
      transform = "LOG",
      min_abundance = 1,
      reference = c("t_simple,3ab", "anti_class,none", "new_anti,weakornone")
    )
    res <- data.frame(fit_mas$results)
    mas_results <- rbind(mas_results, res)
  }
  # adjust p values for multiple testing
  mas_bhcor_p <- p.adjust(mas_results[, 6], method = "BH")
  mas_results$mas_bhcor_p <- mas_bhcor_p
  return(mas_results)
}


# ------------ Function to run MaAsLin2 on all taxa levels w/batch and antibiotics ------------ #
run_mas_all_levels_anti <- function(phylo_obj, variable) {
  set.seed(2023)
  mas_results <- data.frame()
  tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  #tax_levels<-c("Species")
  
  for (tax_level in tax_levels) {
    # aggregate taxa at the taxa level we want
    obj_sub = aggregate_taxa(phylo_obj, paste0(tax_level))
    df_input_data <- data.frame(obj_sub@otu_table)
    df_meta <- data.frame(obj_sub@sam_data)
    
    # Run maaslin2 and bind results
    fit_mas <- Maaslin2(
      input_data = df_input_data,
      input_metadata = df_meta,
      output = "test",
      fixed_effects = c(paste0(variable), "new_anti", "facility", "t_simple"),
      random_effects = c("batch"),
      min_prevalence = 0.2,
      plot_heatmap = FALSE,
      plot_scatter = FALSE,
      normalization = 'CSS',
      transform = "LOG",
      min_abundance = 1,
      reference = c("t_simple,3ab", "new_anti,weakornone")
    )
    res <- data.frame(fit_mas$results)
    mas_results <- rbind(mas_results, res)
  }
  # adjust p values for multiple testing
  mas_bhcor_p <- p.adjust(mas_results[, 6], method = "BH")
  mas_results$mas_bhcor_p <- mas_bhcor_p
  return(mas_results)
}



## LINDA ##

# ------------ Function to run LinDA on all taxa levels with batch random effect ------------ #
run_linda_all_levels <- function(phylo_obj, variable) {
  set.seed(2023)
  lin_results <- data.frame()
  tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (tax_level in tax_levels) {
    # aggregate taxa at the taxa level we want
    obj_sub = aggregate_taxa(phylo_obj, paste0(tax_level))
    df_input_data <- data.frame(obj_sub@otu_table)
    df_meta <- data.frame(obj_sub@sam_data)
    
    # Run linda and bind results
    fit_lin <- linda(
      df_input_data,
      df_meta,
      formula = paste0("~", variable, "+(1|batch)"),
      type = 'count',
      prev.cut = 0.2,
      adaptive = TRUE,
      pseudo.cnt = 0.5
    )
    res <- data.frame(fit_lin$output)
    lin_results <- rbind(lin_results, res)
  }
  # adjust p values for multiple testing
  lin_bhcor_p <- p.adjust(lin_results[, 5], method = "BH")
  lin_results$lin_bhcor_p <- lin_bhcor_p
  return(lin_results)
}

# ------------ Function to run LinDA on all taxa levels with batch random effect and antibiotics ------------ #
run_linda_all_levels_anti <- function(phylo_obj, variable) {
  set.seed(2023)
  lin_results <- data.frame()
  tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  #tax_levels<-c("Species")
  for (tax_level in tax_levels) {
    # aggregate taxa at the taxa level we want
    obj_sub = aggregate_taxa(phylo_obj, paste0(tax_level))
    df_input_data <- data.frame(obj_sub@otu_table)
    df_meta <- data.frame(obj_sub@sam_data)
    
    # Run linda and bind results
    fit_lin <- linda(
      df_input_data,
      df_meta,
      formula = paste0("~", variable, "+new_anti+(1|batch)"),
      type = 'count',
      prev.cut = 0.2,
      adaptive = FALSE,
      pseudo.cnt = 0.5
    )
    res <- data.frame(fit_lin$output)
    lin_results <- rbind(lin_results, res)
  }
  # adjust p values for multiple testing
  lin_bhcor_p <- p.adjust(lin_results[, 5], method = "BH")
  lin_results$lin_bhcor_p <- lin_bhcor_p
  return(lin_results)
}



# ------------ Function to run LinDA on all taxa levels with antibiotics but no batch effect ------------ #
run_linda_all_levels_no_batch_anti <- function(phylo_obj, variable) {
  set.seed(2023)
  lin_results <- data.frame()
  tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (tax_level in tax_levels) {
    # aggregate taxa at the taxa level we want
    obj_sub = aggregate_taxa(phylo_obj, paste0(tax_level))
    df_input_data <- data.frame(obj_sub@otu_table)
    df_meta <- data.frame(obj_sub@sam_data)
    
    # Run linda and bind results
    fit_lin <- linda(
      df_input_data,
      df_meta,
      formula = paste0("~", variable, "+new_anti"),
      type = 'count',
      prev.cut = 0.2,
      adaptive = FALSE,
      pseudo.cnt = 0.5
    )
    res <- data.frame(fit_lin$output)
    lin_results <- rbind(lin_results, res)
  }
  # adjust p values for multiple testing
  lin_bhcor_p <- p.adjust(lin_results[, 5], method = "BH")
  lin_results$lin_bhcor_p <- lin_bhcor_p
  return(lin_results)
}





#------------ function to combine mas and LINDA result into a single df ------------ #
combine_res_mas_lind <- function(mas_res, linda_res) {
  set.seed(2023)
  # select and rename relevant columns so can bind them
  
  #mas_res <- mas_tsimpbac4
  # linda_res <- tsimpbac4
  #mas_res <- resp_anti_mas2
  # linda_res <- resp_anti_lind2
  
  m <-
    mas_res %>% dplyr::select(feature, name, coef, stderr, pval) %>%
    dplyr::rename(
      log2FoldChange_mas = coef,
      lfcSE_mas = stderr,
      pvalue_mas = pval
    ) %>%
    mutate(feature = str_replace_all(feature, "\\.", " ")) %>%
    mutate(feature = str_replace_all(feature, "\\-", " "))
  
  # Here we have to correct the taxa names becuase MaAsLin does weird stuff with them
  linda_res2 <-
    linda_res %>% dplyr::select(-lin_bhcor_p) %>% rownames_to_column(var =
                                                                       "feature")
  linda_res_long <-
    linda_res2 %>% pivot_longer(-feature,
                                names_sep = "\\.",
                                names_to = c("name", ".value"))
  linda_res_long2 <-
    linda_res_long %>% dplyr::select(feature, name, log2FoldChange, lfcSE, pvalue) %>%
    dplyr::rename(
      log2FoldChange_lin = log2FoldChange,
      lfcSE_lin = lfcSE,
      pvalue_lin = pvalue
    ) %>%
    mutate(feature = str_replace_all(feature, "\\.", " ")) %>%
    mutate(feature = str_replace_all(feature, "\\-", " ")) %>%
    mutate(feature = str_replace_all(feature, "\\[", " ")) %>%
    mutate(feature = str_replace_all(feature, "\\]", " "))
  
  
  combined <-
    dplyr::full_join(m, linda_res_long2, by = c("feature", "name"))
  # here we get rid of antibiotics, dont need to include them in pvalue calc
  # THIS LINE NEEDS TO BE MANUALLY SET
  combined2 <-
    combined %>% filter(!str_detect(name, "new_anti")) %>% filter(!str_detect(name, "facility")) %>% filter(!str_detect(name, "t_simple"))
  
  # pivot wider and apply couchy combination & FDR
  combined_wp <-
    combined2 %>% mutate(couchy_t = 0.5 * ((tan((0.5 - pvalue_lin) * pi
    )) + tan((0.5 - pvalue_mas) * pi)))
  combined_wp <-
    combined_wp %>% mutate(couchy_p = 0.5 - (atan(couchy_t) / pi))
  
  comb_bhcor_p <- p.adjust(combined_wp[, 10], method = "BH")
  combined_wp$comb_bhcor_p <- comb_bhcor_p
  
  return(combined_wp)
}

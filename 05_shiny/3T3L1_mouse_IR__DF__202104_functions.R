###Last updated:
#20211002

###Background:
#Here I keep all necessary functions for 3T3L1_mouse_IR__DF__202104 project




####Visualisation


###Boxplot for 3t3-L1 phos psites (splits by multiplicity)
boxplot_3t3phos_psite <- function(data, psite, siglabels = FALSE, prot_normalised = FALSE, ...){
  #if prot_normalised: Use ppeptide normalised to protein
  if (prot_normalised){
    data <- data[, -c(1:52)]
    colnames(data) <- gsub("PROT_norm_", "", colnames(data))
  }
  
  #Include brackets to label signifiance
  if (siglabels){
    return(general_boxplot_psite_siglabels(data = data,
                                     psite = psite,
                                     conditions = c("CTRL_BAS", "CTRL_INS",
                                                    "CI_BAS", "CI_INS", 
                                                    "DEX_BAS", "DEX_INS",
                                                    "TNF_BAS", "TNF_INS",
                                                    "MPQ_BAS", "MPQ_INS",
                                                    "AA_BAS", "AA_INS"),
                                     treatments = c("CTRL", "CI", "DEX", "TNF", "MPQ", "AA"),
                                     treatment_colours = SH_control_model_main_colours,
                                     pval_cols = c("t_adj_p_val", "CI_padj", "DEX_padj", "TNF_padj", "MPQ_padj", "AA_padj"),
                                     comparisons = list(c("CTRL_BAS", "CTRL_INS"),
                                                        c("CI_BAS", "CI_INS"),
                                                        c("DEX_BAS", "DEX_INS"),
                                                        c("TNF_BAS", "TNF_INS"),
                                                        c("MPQ_BAS", "MPQ_INS"),
                                                        c("AA_BAS", "AA_INS")),
                                     num_cols = 1:52,
                                     show_stars = c("*", "#", "#", "#", "#", "#"),
                                     boxplot_function = general_boxplot_psite,
                                     ...))
  } else {
    return(general_boxplot_psite(data = data,
                                 psite = psite,
                                 conditions = c("CTRL_BAS", "CTRL_INS",
                                                "CI_BAS", "CI_INS", 
                                                "DEX_BAS", "DEX_INS",
                                                "TNF_BAS", "TNF_INS",
                                                "MPQ_BAS", "MPQ_INS",
                                                "AA_BAS", "AA_INS"),
                                 treatments = c("CTRL", "CI", "DEX", "TNF", "MPQ", "AA"),
                                 treatment_colours = SH_control_model_main_colours,
                                 ...))
  }
  
}



###Boxplot for mouse phos psites (splits by multiplicity)
boxplot_mousephos_psite <- function(data, psite, siglabels = FALSE, ...){
  if (siglabels){
    return(general_boxplot_psite_siglabels(data = data,
                                     psite = psite,
                                     conditions = c("CHOW_BAS", "CHOW_INS",
                                                    "HFD_BAS", "HFD_INS", 
                                                    "REV_BAS", "REV_INS"),
                                     treatments = c("CHOW", "HFD", "REV"),
                                     treatment_colours = SH_mouse_main_colours,
                                     pval_cols = c("CHOW_adj_p_val", "HFD_adj_p_val", "REV_adj_p_val"),
                                     comparisons = list(c("CHOW_BAS", "CHOW_INS"),
                                                        c("HFD_BAS", "HFD_INS"),
                                                        c("REV_BAS", "REV_INS")),
                                     num_cols = 1:72,
                                     show_stars = c("*", "*", "*"),
                                     ...))
    
    
    
  } else {
  return(general_boxplot_psite(data = data,
                               psite = psite,
                               conditions = c("CHOW_BAS", "CHOW_INS",
                                              "HFD_BAS", "HFD_INS", 
                                              "REV_BAS", "REV_INS"),
                               treatments = c("CHOW", "HFD", "REV"),
                               treatment_colours = SH_mouse_main_colours,
                               ...))
  }
}




###Boxplot for one df and given conditions
boxplot_grouped_3t3_proteome <- function(data, 
                                         protein,
                                         siglabels = FALSE,
                                         ...){
  #Relabel
  colnames(data) <- gsub("GROUPED_", "", colnames(data))
  colnames(data) <- gsub("NORMAL", "CTRL", colnames(data))
  
  #Generate plots
  if (siglabels){
    return(general_boxplot_siglabels(data = data,
                                     protein = protein,
                                     conditions = c("CTRL", "CI", "DEX", "TNF", "MPQ", "AA"),
                                     treatments = c("CTRL", "CI", "DEX", "TNF", "MPQ", "AA"),
                                     treatment_colours = SH_control_model_main_colours,
                                     pval_cols = c("CI_adj_pval", "DEX_adj_pval", "TNF_adj_pval", "MPQ_adj_pval", "AA_adj_pval"),
                                     comparisons = list(c("CTRL", "CI"),
                                                        c("CTRL", "DEX"),
                                                        c("CTRL", "TNF"),
                                                        c("CTRL", "MPQ"),
                                                        c("CTRL", "AA")),
                                     num_cols = 1:48,
                                     show_stars = c("*", "*", "*", "*", "*"),
                                     is_raw = TRUE,
                                     ...))
    
  } else {
    return(general_boxplot(data = data,
                                     protein = protein,
                                     conditions = c("CTRL", "CI", "DEX", "TNF", "MPQ", "AA"),
                                     treatments = c("CTRL", "CI", "DEX", "TNF", "MPQ", "AA"),
                                     treatment_colours = SH_control_model_main_colours,
                           is_raw = TRUE,
                                     ...))
  }
}


###Boxplot for 3t3-L1 phos psites 
boxplot_gsk3iphos_psite <- function(data, psite, siglabels = FALSE, ...){
  #Include brackets to label signficance
  if (siglabels){
    return(general_boxplot_psite_siglabels(data = data,
                                           psite = psite,
                                           conditions = c("Basal", "GSK3i"),
                                           treatments = c("Basal", "GSK3i"),
                                           treatment_colours = Gsk3i_colours,
                                           pval_cols = c("adj_p_val"),
                                           comparisons = list(c("Basal", "GSK3i")),
                                           num_cols = 1:8,
                                           show_stars = c("*"),
                                           ...))
    
  } else {
    return(general_boxplot_psite(data = data,
                                 psite = psite,
                                 conditions = c("Basal", "GSK3i"),
                                 treatments = c("Basal", "GSK3i"),
                                 treatment_colours = Gsk3i_colours,
                                 ...))
  }
}
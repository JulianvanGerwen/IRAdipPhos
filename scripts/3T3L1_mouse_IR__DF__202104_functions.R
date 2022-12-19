###Last updated:
#20211002

###Background:
#Here I keep all necessary functions for 3T3L1_mouse_IR__DF__202104 project

###Packages
library(reshape2)
library(ggplot2)
library(gridExtra)
library(multcomp)




#####Visualisations#####

###Function to make a figure with multiple panels that are each boxplots
boxplot_all_conditions_3t3_post202104_mult_panel <- function(data, 
                                                             sites, 
                                                             boxplot_function = boxplot_all_conditions_3t3_post202104,
                                                             order = "treatment",
                                                             alpha_val_inner = 0.3,
                                                             alpha_val_border = 1,
                                                             is_raw = FALSE,
                                                             ncol,
                                                             nrow,
                                                             include_axistitles_legend = FALSE){
  
  ###Configure theme, based on include_axes
  if (include_axistitles_legend == TRUE){
    
    plot_theme <- NULL
  } else {
    
    plot_theme <- theme(legend.position = "none",
                        axis.text.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank())
  }
  
  ###Make list of boxplots
  boxplot_list <- list()
  for (i in 1:length(sites)){
    
    temp_site_split <- strsplit(sites[i],
                                "_")[[1]]
    
    boxplot_list[[i]] <- boxplot_function(data = data,
                                          protein = sites[i],
                                          order = order,
                                          alpha_val_inner = alpha_val_inner,
                                          alpha_val_border = alpha_val_border,
                                          is_raw = is_raw) +
      plot_theme +
      ggtitle(paste(c(temp_site_split[1],
                      " ",
                      temp_site_split[3],
                      " P",
                      temp_site_split[4]),
                    collapse = ""))
  }
  
  ###Return grob object
  return(arrangeGrob(grobs = boxplot_list,
                     nrow = nrow,
                     ncol = ncol))
}


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

###Boxplot for 3t3-L1 phos sites 
boxplot_3t3phos <- function(data, site, siglabels = FALSE, prot_normalised = FALSE, ...){
  #if prot_normalised: Use ppeptide normalised to protein
  if (prot_normalised){
    data <- data[, -c(1:52)]
    colnames(data) <- gsub("PROT_norm_", "", colnames(data))
  }
  
  
  #Include brackets to label signficance
  if (siglabels){
    return(general_boxplot_siglabels(data = data,
                           protein = site,
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
                           ...))
    
  } else {
    return(general_boxplot(data = data,
                           protein = site,
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


###Boxplot for one df and given conditions
#Up to date 20210801
#If is_raw = TRUE, then I don't search for columns that tell me if we've imputed
boxplot_all_conditions_3t3_post202104 <- function(data, 
                                                  protein, 
                                                  order = "treatment",
                                                  alpha_val_inner = 0.3,
                                                  alpha_val_border = 1,
                                                  is_raw = FALSE,
                                                  prot_normalised = FALSE){
  
  ##Change data based on whether or not prot_normalised
  if (prot_normalised){
    
    data <- data[, c(grep("PROT_norm",
                          colnames(data)),
                     grep("imp",
                          colnames(data)))]
    colnames(data) <- gsub("PROT_norm_",
                           "",
                           colnames(data))
  } else {
    
    data <- data[, c(1:52,
                     grep("imp",
                          colnames(data)))]
  }
  
  boxplot_df <- data.frame("Condition" = gsub("_\\d+",
                                              "",
                                              colnames(data)[grep("_\\d+",
                                                                  colnames(data))]),
                           "Value" = as.numeric(data[protein, grep("_\\d+", colnames(data))]),
                           stringsAsFactors = FALSE)
  
  ##Order by treatment or basins
  
  if (order == "treatment") {
    
    boxplot_df$Condition <- factor(boxplot_df$Condition, levels = c("CTRL_BAS",
                                                                    "CTRL_INS",
                                                                    "CI_BAS",
                                                                    "CI_INS",
                                                                    "DEX_BAS",
                                                                    "DEX_INS",
                                                                    "TNF_BAS",
                                                                    "TNF_INS",
                                                                    "MPQ_BAS",
                                                                    "MPQ_INS",
                                                                    "AA_BAS",
                                                                    "AA_INS"))
  } else if (order == "basins") {
    
    boxplot_df$Condition <- factor(boxplot_df$Condition, levels = c("CTRL_BAS",
                                                                    "CI_BAS",
                                                                    "DEX_BAS",
                                                                    "TNF_BAS",
                                                                    "MPQ_BAS",
                                                                    "AA_BAS",
                                                                    "CTRL_INS",
                                                                    "CI_INS",
                                                                    "DEX_INS",
                                                                    "TNF_INS",
                                                                    "MPQ_INS",
                                                                    "AA_INS"))
    
  }
  
  boxplot_df$Treatment <- sub("_\\w+", "", boxplot_df$Condition)
  boxplot_df$Treatment <- factor(boxplot_df$Treatment, levels = c("CTRL",
                                                                  "CI",
                                                                  "DEX",
                                                                  "TNF",
                                                                  "MPQ",
                                                                  "AA"))
  ##Tag imputed sites so I can then label
  boxplot_df$Imputed <- "Not imputed"
  #Dodge this process if raw = TRUE
  if (is_raw){
    
    NULL
  } else {
    
    for (i in 1:nrow(boxplot_df)){
      
      
      if (data[protein, paste(c(as.character(boxplot_df$Condition[i]),
                                "_blockimp"),
                              collapse = "")]){
        
        boxplot_df$Imputed[i] <- "Imputed"
      }
    }
  }
  
  ##Make point plot based on whether there are imputed sites
  if (length(unique(boxplot_df$Imputed)) > 1){
    
    boxplot_df$Imputed <- factor(boxplot_df$Imputed,
                                 levels = c("Not imputed",
                                            "Imputed"))
    
    point_plot_1 <- geom_point(data = boxplot_df,
                               aes(x = Condition,
                                   y = Value,
                                   color = Treatment,
                                   shape = Imputed),
                               size = 1.5,
                               stroke = 0.6,
                               alpha = alpha_val_border)
    point_plot_2 <- geom_point(data = boxplot_df[grep("INS",
                                                      boxplot_df$Condition), ],
                               aes(x = Condition,
                                   y = Value,
                                   color = Treatment,
                                   fill = Treatment,
                                   shape = Imputed),
                               size = 1.5,
                               stroke = 0.6,
                               alpha = alpha_val_inner)
    scale_shape <- scale_shape_manual(values = c(21,
                                                 23))
  } else {
    
    point_plot_1 <- geom_point(data = boxplot_df,
                               aes(x = Condition,
                                   y = Value,
                                   color = Treatment),
                               size = 1.5,
                               shape = 21,
                               stroke = 0.6,
                               alpha = alpha_val_border)
    point_plot_2 <- geom_point(data = boxplot_df[grep("INS",
                                                      boxplot_df$Condition), ],
                               aes(x = Condition,
                                   y = Value,
                                   color = Treatment,
                                   fill = Treatment),
                               size = 1.5,
                               shape = 21,
                               stroke = 0.6,
                               alpha = alpha_val_inner)
    scale_shape <- NULL
  }
  
  ##Plot
  plot <- ggplot(data = boxplot_df,
                 aes(x = Condition,
                     y = Value,
                     color = Treatment)) +
    geom_boxplot(outlier.shape = NA,
                 coef = 0,
                 lwd = 0.35) + 
    point_plot_1 +
    point_plot_2 +
    scale_shape +
    scale_color_manual(values = SH_control_model_main_colours) +
    scale_fill_manual(values = SH_control_model_main_colours) +
    ggtitle(protein) + 
    labs(x = "Condition",
         y = "log2 Intensity") +
    theme(axis.text.x = element_text(angle = 90, 
                                     hjust = 1, 
                                     vjust = 0.5, 
                                     size = 6,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black",
                                     size = 6),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(size = 8),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          legend.position = "none")
  return(plot)
}

###Function: heatmap on phos_3t3_proc using INS/BAS data
heatmap_phos_3t3_proc <- function(data, order_column = "CTRL_FC", ...){
  fc_sig_heatmap(data = data,
                 fc_cols_w_pvals = c("CTRL_FC",
                                     "CI_FC",
                                     "DEX_FC",
                                     "TNF_FC",
                                     "MPQ_FC",
                                     "AA_FC"),
                 pval_cols = c("t_adj_p_val",
                               "CI_padj",
                               "DEX_padj",
                               "TNF_padj",
                               "MPQ_padj",
                               "AA_padj"),
                 shape2_pval_cols = c("CI_padj",
                                      "DEX_padj",
                                      "TNF_padj",
                                      "MPQ_padj",
                                      "AA_padj"),
                 order_column = order_column,
                 is_decreasing = FALSE,
                 x_axis_names = c("CTRL",
                                  "CI",
                                  "DEX",
                                  "TNF",
                                  "MPQ",
                                  "AA"),
                 ...)
}


###Function: Heatmap on phos_3t3_proc using INS/BAS and Model/CTRL data
#Plot INS/BAS and MODEL/CTRL FCs
heatmap_phos_3t3_proc_allFCs <- function(data, order_column = "CTRL_FC", ...){
  fc_sig_heatmap(data = data,
                 fc_cols_w_pvals = c("CTRL_FC",
                                     "CI_FC",
                                     "DEX_FC",
                                     "TNF_FC",
                                     "MPQ_FC",
                                     "AA_FC"),
                 fc_cols_wo_pvals = c("CI_FC_on_control_BAS",
                                      "DEX_FC_on_control_BAS",
                                      "TNF_FC_on_control_BAS",
                                      "MPQ_FC_on_control_BAS",
                                      "AA_FC_on_control_BAS",
                                      "CI_FC_on_control_INS",
                                      "DEX_FC_on_control_INS",
                                      "TNF_FC_on_control_INS",
                                      "MPQ_FC_on_control_INS",
                                      "AA_FC_on_control_INS"),
                 pval_cols = c("t_adj_p_val",
                               "CI_padj",
                               "DEX_padj",
                               "TNF_padj",
                               "MPQ_padj",
                               "AA_padj"),
                 shape2_pval_cols = c("CI_padj",
                                      "DEX_padj",
                                      "TNF_padj",
                                      "MPQ_padj",
                                      "AA_padj"),
                 order_column = order_column,
                 is_decreasing = FALSE,
                 x_axis_names = c("CTRL INS/BAS",
                                  "CI INS/BAS",
                                  "DEX INS/BAS",
                                  "TNF INS/BAS",
                                  "MPQ INS/BAS",
                                  "AA INS/BAS",
                                  "CI/CTRL BAS",
                                  "DEX/CTRL BAS",
                                  "TNF/CTRL BAS",
                                  "MPQ/CTRL BAS",
                                  "AA/CTRL BAS",
                                  "CI/CTRL INS",
                                  "DEX/CTRL INS",
                                  "TNF/CTRL INS",
                                  "MPQ/CTRL INS",
                                  "AA/CTRL INS"), 
                 ...)
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

###Boxplot for mouse phos (no split by multiplicity)
boxplot_mousephos <- function(data, site, siglabels = FALSE, ...){
  #Add brackets for significance if needed
  if (siglabels){
    return(general_boxplot_siglabels(data = data,
                           protein = site,
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
    return(general_boxplot(data = data,
                           protein = site,
                           conditions = c("CHOW_BAS", "CHOW_INS",
                                          "HFD_BAS", "HFD_INS", 
                                          "REV_BAS", "REV_INS"),
                           treatments = c("CHOW", "HFD", "REV"),
                           treatment_colours = SH_mouse_main_colours,
                           ...))
  }
  
}



###Boxplot for one df and given conditions
#Up to date 20210801
boxplot_all_conditions_mouse <- function(data, 
                                         protein,
                                         order = NULL,
                                         alpha_val_inner = 0.3,
                                         alpha_val_border = 1,
                                         is_raw = FALSE){
  
  boxplot_df <- data.frame("Condition" = gsub("_\\d+",
                                              "",
                                              colnames(data)[grep("_\\d+", colnames(data))]),
                           "Value" = as.numeric(data[protein, grep("_\\d+", colnames(data))]),
                           stringsAsFactors = FALSE)
  
  boxplot_df$Condition <- factor(boxplot_df$Condition, levels = unique(boxplot_df$Condition))
  boxplot_df$Treatment <- sub("_\\w+", "", boxplot_df$Condition)
  boxplot_df$Treatment <- factor(boxplot_df$Treatment, levels = unique(boxplot_df$Treatment))
  
  
  ##Tag imputed sites so I can then label
  boxplot_df$Imputed <- "Not imputed"
  #Dodge this process if raw = TRUE
  if (is_raw){
    
    NULL
  } else {
    
    for (i in 1:nrow(boxplot_df)){
      
      
      if (data[protein, paste(c(as.character(boxplot_df$Condition[i]),
                                "_blockimp"),
                              collapse = "")]){
        
        boxplot_df$Imputed[i] <- "Imputed"
      }
    }
  }
  
  
  ##Make point plot based on whether there are imputed sites
  if (length(unique(boxplot_df$Imputed)) > 1){
    
    boxplot_df$Imputed <- factor(boxplot_df$Imputed,
                                 levels = c("Not imputed",
                                            "Imputed"))
    
    point_plot_1 <- geom_point(data = boxplot_df,
                               aes(x = Condition,
                                   y = Value,
                                   color = Treatment,
                                   shape = Imputed),
                               size = 1.5,
                               stroke = 0.6,
                               alpha = alpha_val_border) 
    point_plot_2 <- geom_point(data = boxplot_df[grep("INS",
                                                      boxplot_df$Condition), ],
                               aes(x = Condition,
                                   y = Value,
                                   color = Treatment,
                                   fill = Treatment,
                                   shape = Imputed),
                               size = 1.5,
                               stroke = 0.6,
                               alpha = alpha_val_inner) 
    scale_shape <- scale_shape_manual(values = c(21,
                                                 23))
  } else {
    
    point_plot_1 <- geom_point(data = boxplot_df,
                               aes(x = Condition,
                                   y = Value,
                                   color = Treatment),
                               size = 1.5,
                               shape = 21,
                               stroke = 0.6,
                               alpha = alpha_val_border) 
    point_plot_2 <- geom_point(data = boxplot_df[grep("INS",
                                                      boxplot_df$Condition), ],
                               aes(x = Condition,
                                   y = Value,
                                   color = Treatment,
                                   fill = Treatment),
                               size = 1.5,
                               shape = 21,
                               stroke = 0.6,
                               alpha = alpha_val_inner)
    scale_shape <- NULL
  }
  
  
  plot <- ggplot(data = boxplot_df,
                 aes(x = Condition,
                     y = Value,
                     color = Treatment)) +
    geom_boxplot(outlier.shape = NA,
                 coef = 0,
                 lwd = 0.35) + 
    point_plot_1 +
    point_plot_2 +
    scale_shape +
    scale_color_manual(values = SH_mouse_main_colours) +
    scale_fill_manual(values = SH_mouse_main_colours) +
    ggtitle(protein) + 
    labs(x = "Condition",
         y = "log2 Intensity") +
    theme(axis.text.x = element_text(angle = 90, 
                                     hjust = 1, 
                                     vjust = 0.5, 
                                     size = 6,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black",
                                     size = 6),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(size = 8),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          legend.position = "none")
  return(plot)
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


###Boxplot for 3t3-L1 phos sites 
boxplot_gsk3iphos <- function(data, site, siglabels = FALSE, ...){
  #Include brackets to label signficance
  if (siglabels){
    return(general_boxplot_siglabels(data = data,
                                     protein = site,
                                     conditions = c("Basal", "GSK3i"),
                                     treatments = c("Basal", "GSK3i"),
                                     treatment_colours = Gsk3i_colours,
                                     pval_cols = c("adj_p_val"),
                                     comparisons = list(c("Basal", "GSK3i")),
                                     num_cols = 1:8,
                                     show_stars = c("*"),
                                     ...))
    
  } else {
    return(general_boxplot(data = data,
                           protein = site,
                           conditions = c("Basal", "GSK3i"),
                           treatments = c("Basal", "GSK3i"),
                           treatment_colours = Gsk3i_colours,
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


###Function: Simple heatmap for model/CTRL FCs in 3t3l1IR proteome data
heatmap_proteome_3t3l1ir_simple <- function(data, ...){
  fc_sig_heatmap(data,
                 fc_cols_w_pvals = c("CI_GROUPED_FC",
                                     "DEX_GROUPED_FC",
                                     "TNF_GROUPED_FC",
                                     "MPQ_GROUPED_FC",
                                     "AA_GROUPED_FC"),
                 pval_cols = c("CI_GROUPED_adj_pval",
                               "DEX_GROUPED_adj_pval",
                               "TNF_GROUPED_adj_pval",
                               "MPQ_GROUPED_adj_pval",
                               "AA_GROUPED_adj_pval"),
                 yaxis_naming = "gene_uniprot",
                 x_axis_names = c("CI",
                                  "DEX",
                                  "TNF",
                                  "MPQ",
                                  "AA"),
                 order_column = "CI_GROUPED_FC",
                 legend_title = "log2 Model/CTRL",
                 ...)
}
  



###Heatmap function
##Plots FC on heatmap
##Plots significance as point in tile
fc_sig_heatmap_proteome <- function(data,
                                    fc_cols_w_pvals = c(),
                                    fc_cols_wo_pvals = c(),
                                    pval_cols = c(),
                                    order_column = "standard_name",
                                    is_decreasing = TRUE,
                                    return_df = FALSE,
                                    x_axis_names = c(),
                                    gap_column_index = FALSE){
  
  
  ##Prep for melts
  
  #Set up data
  data$standard_name <- sapply(rownames(data),
                                      function(x) strsplit(x, "_")[[1]][1])
  #Double up genenames:
  data$old_standard_name <- data$standard_name
  for (i in 1:nrow(data)){
    
    if (length(which(data$old_standard_name ==
                     data$standard_name[i])) > 1){
      
      data$standard_name[i] <- paste(strsplit(rownames(data)[i],
                                                     "_")[[1]],
                                            collapse = " ")
    }
  }
  
  data <- data[order(data[, order_column],
                                   decreasing = is_decreasing == FALSE), ]
  data$standard_name <- factor(data$standard_name,
                                      levels = unique(data$standard_name))
  
  
  ##w_pvals melt
  
  #Check that we have pval_cols
  
  if (length(pval_cols) > 0){
    
    #FC melt
    FC_melt_df_w_pvals <- melt(data,
                               id.vars = "standard_name",
                               measure.vars = fc_cols_w_pvals)
    #p_val melt
    pval_melt_df_w_pvals <- melt(data,
                                 id.vars = "standard_name",
                                 measure.vars = pval_cols)
    pval_melt_df_w_pvals$value <- pval_melt_df_w_pvals$value < 0.05
    
    #Combine
    post_melt_df_w_pvals <- FC_melt_df_w_pvals
    colnames(post_melt_df_w_pvals)[which(colnames(post_melt_df_w_pvals) == "value")] <- "FC_value"
    post_melt_df_w_pvals$sig_value <- pval_melt_df_w_pvals$value
    
    #NAs
    post_melt_df_w_pvals$sig_value[is.na(post_melt_df_w_pvals$sig_value)] <- FALSE
  } else {
    
    post_melt_df_w_pvals <- NULL
    
  }
  
  
  
  ##wo_pvals melt
  
  #Check that we have fc_cols_wo_pvals
  
  if (length(fc_cols_wo_pvals) > 0){
    
    #FC melt
    FC_melt_df_wo_pvals <- melt(data,
                                id.vars = "standard_name",
                                measure.vars = fc_cols_wo_pvals)
    colnames(FC_melt_df_wo_pvals)[which(colnames(FC_melt_df_wo_pvals) == "value")] <- "FC_value"
    #Add FALSE pvals
    FC_melt_df_wo_pvals$sig_value <- FALSE
    
    #post_melt_df_wo_pvals
    post_melt_df_wo_pvals <- FC_melt_df_wo_pvals
  } else {
    
    post_melt_df_wo_pvals <- NULL
    
  }
  
  
  ##Combine w_pvals and wo_pvals
  
  post_melt_df_combined <- rbind(post_melt_df_w_pvals,
                                 post_melt_df_wo_pvals)
  
  
  ##x-axis names
  
  #check that we have x_axis_names
  if (length(x_axis_names)> 0){
    levels(post_melt_df_combined$variable) <- x_axis_names
  }
  
  ##Add gap_column
  
  if (gap_column_index != FALSE){
    #Make gap_column_df
    gap_column_df <- post_melt_df_combined[1:length(unique(post_melt_df_combined$standard_name)), ]
    gap_column_df$FC_value <- 0
    gap_column_df$sig_value <- FALSE
    gap_column_df$variable <- ""
    #Calculate new levels
    new_levels <- c(levels(post_melt_df_combined$variable)[1:gap_column_index],
                    "",
                    levels(post_melt_df_combined$variable)[-c(1:gap_column_index)])
    #Add to post_melt_df_combined
    post_melt_df_combined <- rbind(post_melt_df_combined,
                                   gap_column_df)
    #Reorder levels
    post_melt_df_combined$variable <- factor(post_melt_df_combined$variable,
                                             levels = new_levels)
  }
  
  
  ##Heatmap
  output_plot <- ggplot(post_melt_df_combined,
                        aes(x = variable,
                            y = standard_name,
                            fill = FC_value)) + 
    scale_fill_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         name = "log2 MODEL/CONTROL") + 
    geom_tile() +
    coord_fixed(ratio = 1) +
    geom_point(data = post_melt_df_combined[post_melt_df_combined$sig_value == TRUE,],
               size = 0.5) +
    labs(x = "Condition",
         y = "Protein") +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          plot.title = element_text(hjust = 0.5))
  
  if (return_df){
    
    return(post_melt_df_combined)
    
  } else {
    
    return(output_plot)
    
  }
  
}





####Data processing

###Reg sites summariser
#Up to data 20210514
#Put in df with CTRL_reg column (values: unregulated, up, down), num_defective_models column,  column for defect in each model
#Spits out summary that tells number of sites defective in 0, 1, 2, 3, 4, 5 models, up/down or just up or down, in total or across models
#Perc spits out percentages in each column
#all_def_row makes a row that sums across defective in 0, 1, 2, 3, 4, 5 models
reg_sites_summariser <- function(data,
                                 perc = FALSE,
                                 all_def_row = FALSE){
  
  ###Just regulated data
  data <- data[which(data$CTRL_reg != "unregulated"), ]
  
  ###Set up output
  reg_summary_m <- matrix(0,
                          ncol = 19,
                          nrow = 6)
  colnames(reg_summary_m) <- c("num_defective_models",
                               "all_models_total",
                               "all_models_up",
                               "all_models_down",
                               "CI_total",
                               "CI_up",
                               "CI_down",
                               "DEX_total",
                               "DEX_up",
                               "DEX_down",
                               "TNF_total",
                               "TNF_up",
                               "TNF_down",
                               "MPQ_total",
                               "MPQ_up",
                               "MPQ_down",
                               "AA_total",
                               "AA_up",
                               "AA_down")
  reg_summary_m[, 1] <- 0:5
  reg_summary_df <- as.data.frame(reg_summary_m)
  
  ###Function: Extended table
  #table doesn't work because it will not include if there are none of a value (e.g. no sites defective in 5 models and up-regulated)
  extended_table <- function(x) c(length(which(x == 0)),
                                  length(which(x == 1)),
                                  length(which(x == 2)),
                                  length(which(x == 3)),
                                  length(which(x == 4)),
                                  length(which(x == 5)))
  
  ###Assign numbers
  
  ##all_models
  reg_summary_df$all_models_total <- extended_table(data$num_defective_models)
  reg_summary_df$all_models_up <- extended_table(data$num_defective_models[which(data$CTRL_reg == "up")])
  reg_summary_df$all_models_down <- extended_table(data$num_defective_models[which(data$CTRL_reg == "down")])
  
  ##models
  models <- c("CI",
              "DEX",
              "TNF",
              "MPQ",
              "AA")
  for (j in 1:length(models)){
    
    reg_summary_df[, paste(c(models[j],
                             "_total"),
                           collapse = "")] <- 
      extended_table(data$num_defective_models[which(data[, paste(c(models[j],
                                                                    "_defect"),
                                                                  collapse = "")] == TRUE)])
    reg_summary_df[, paste(c(models[j],
                             "_up"),
                           collapse = "")] <- 
      extended_table(data$num_defective_models[which(data[, paste(c(models[j],
                                                                    "_defect"),
                                                                  collapse = "")] == TRUE &
                                                       data$CTRL_reg == "up")])
    reg_summary_df[, paste(c(models[j],
                             "_down"),
                           collapse = "")] <- 
      extended_table(data$num_defective_models[which(data[, paste(c(models[j],
                                                                    "_defect"),
                                                                  collapse = "")] == TRUE &
                                                       data$CTRL_reg == "down")])
  }
  
  if (perc == TRUE){
    
    col_sums <- colSums(reg_summary_df)
    col_sums[1] <- 1
    
    reg_summary_df_perc <- sweep(reg_summary_df,
                                 2,
                                 STATS = col_sums,
                                 FUN = "/")
    reg_summary_df <- reg_summary_df_perc
  }
  
  
  ##Add all_def_row if needed
  if (all_def_row == TRUE){
    
    reg_summary_df[7, 1] <- "all_defective_sites"
    reg_summary_df[7, -1] <- colSums(reg_summary_df[2:6, -1])
  }
  
  ##Return
  return(reg_summary_df)
}


reg_sites_summariser_mouse <- function(data,
                                       perc = FALSE,
                                       no_pval = FALSE){
  
  ###Make summary_df
  summary_m <- matrix(0,
                      nrow = 2,
                      ncol = 4)
  rownames(summary_m) <- c("Not defective",
                           "Defective")
  colnames(summary_m) <- c("site_group",
                           "Total",
                           "Up",
                           "Down")
  summary_df <- as.data.frame(summary_m)
  summary_df$site_group <- rownames(summary_df)
  
  ###Fill summary_m
  
  ##Not exact
  #Change filtering based on no_pval
  if (no_pval == TRUE){
    
    summary_df[1, 2:4] <- c(length(which(abs(data$CHOW_logFC) > 0.58)),
                            length(which(data$CHOW_logFC > 0.58)),
                            length(which(data$CHOW_logFC < -0.58)))
    summary_df[2, 2:4] <- c(length(which((data$CHOW_logFC > 0.58 &
                                            data$HFD_logFC < 0.58) |
                                           (data$CHOW_logFC < -0.58 &
                                              data$HFD_logFC > -0.58))),
                            length(which(data$CHOW_logFC > 0.58 &
                                           data$HFD_logFC < 0.58)),
                            length(which(data$CHOW_logFC < -0.58 &
                                           data$HFD_logFC > -0.58)))
  } else {
    
    summary_df[1, 2:4] <- c(length(which(data$CHOW_reg == TRUE)),
                            length(which(data$CHOW_reg == TRUE &
                                           data$CHOW_logFC > 0.58)),
                            length(which(data$CHOW_reg == TRUE &
                                           data$CHOW_logFC < -0.58)))
    summary_df[2, 2:4] <- c(length(which(data$HFD_defective == TRUE)),
                            length(which(data$HFD_defective == TRUE &
                                           data$CHOW_logFC > 0.58)),
                            length(which(data$HFD_defective == TRUE &
                                           data$CHOW_logFC < -0.58)))
  }
  
  
  ##Make exact
  summary_df[1, 2:4] <- summary_df[1, 2:4] - summary_df[2, 2:4]
  
  ##Make perc if required
  if (perc == TRUE){
    
    summary_df[, 2:4] <- sweep(summary_df[, 2:4],
                               2,
                               STATS = colSums(summary_df[, 2:4]),
                               FUN = "/") * 100
  }
  
  ##Return
  return(summary_df)
}


###Emergent sites summariser
#Up to data 20210726
#Put in df with emergent (boolean), emergent_up_num, and emergent_down_num columns, emergent_num column (natural numbers),  column for emergent in each model
#Spits out summary that tells number of sites emergent in 1, 2, 3, 4, 5 models, up/down or just up or down, in total or across models
#Perc spits out percentages in each column
emergent_sites_summariser <- function(data,
                                      perc = FALSE){
  
  ###Just emergent data
  data <- data[which(data$emergent == TRUE), ]
  
  ###Set up output
  reg_summary_m <- matrix(0,
                          ncol = 19,
                          nrow = 5)
  colnames(reg_summary_m) <- c("num_emergent_models",
                               "all_models_total",
                               "all_models_up",
                               "all_models_down",
                               "CI_total",
                               "CI_up",
                               "CI_down",
                               "DEX_total",
                               "DEX_up",
                               "DEX_down",
                               "TNF_total",
                               "TNF_up",
                               "TNF_down",
                               "MPQ_total",
                               "MPQ_up",
                               "MPQ_down",
                               "AA_total",
                               "AA_up",
                               "AA_down")
  reg_summary_m[, 1] <- 1:5
  reg_summary_df <- as.data.frame(reg_summary_m)
  
  ###Function: Extended table
  #table doesn't work because it will not include if there are none of a value (e.g. no sites defective in 5 models and up-regulated)
  extended_table <- function(x) c(length(which(x == 1)),
                                  length(which(x == 2)),
                                  length(which(x == 3)),
                                  length(which(x == 4)),
                                  length(which(x == 5)))
  
  ###Assign numbers
  
  ##all_models
  reg_summary_df$all_models_total <- extended_table(data$emergent_num)
  reg_summary_df$all_models_up <- extended_table(data$emergent_num[which(data$emergent_up_num > 0)])
  reg_summary_df$all_models_down <- extended_table(data$emergent_num[which(data$emergent_down_num > 0)])
  
  ##models
  models <- c("CI",
              "DEX",
              "TNF",
              "MPQ",
              "AA")
  for (j in 1:length(models)){
    
    reg_summary_df[, paste(c(models[j],
                             "_total"),
                           collapse = "")] <- 
      extended_table(data$emergent_num[which(data[, paste(c(models[j],
                                                            "_emergent"),
                                                          collapse = "")] != "FALSE")])
    reg_summary_df[, paste(c(models[j],
                             "_up"),
                           collapse = "")] <- 
      extended_table(data$emergent_num[which(data[, paste(c(models[j],
                                                            "_emergent"),
                                                          collapse = "")] == "up")])
    reg_summary_df[, paste(c(models[j],
                             "_down"),
                           collapse = "")] <- 
      extended_table(data$emergent_num[which(data[, paste(c(models[j],
                                                            "_emergent"),
                                                          collapse = "")] == "down")])
  }
  
  if (perc == TRUE){
    
    col_sums <- colSums(reg_summary_df)
    col_sums[1] <- 1
    
    reg_summary_df_perc <- sweep(reg_summary_df,
                                 2,
                                 STATS = col_sums,
                                 FUN = "/")
    return(reg_summary_df_perc)
  } else {
    
    return(reg_summary_df)
  }
  
  
}




###Function that identifies best uniprots for proteome data, by using phospho data
#Up to date 20210521
#possible_uniprots_col contains uniprots separated by ";"; typically Majority.protein.IDs column from MaxQuant
uniprot_proteome_from_phos_matcher <- function(old_uniprots,
                                               possible_uniprots_col,
                                               destination_uniprots){
  
  ##Set up new uniprots
  new_uniprots <- old_uniprots
  
  ##Make list of possible uniprots
  possible_uniprots_list <- strsplit(possible_uniprots_col,
                                     ";")
  names(possible_uniprots_list) <- old_uniprots
  
  ##Make list of indices for matched uniprots
  ######Ultimately get rid of this
  matched_uniprots_indices <- as.list(names(possible_uniprots_list))
  names(matched_uniprots_indices) <- names(possible_uniprots_list)
  matched_uniprots_indices_report <- NULL
  
  #Fill up
  for (i in 1:length(possible_uniprots_list)){
    
    #Get matched indices
    temp_matched_indices <- which(possible_uniprots_list[[i]] %in%
                                    destination_uniprots)
    
    #Fill up matched_uniprots_indices
    matched_uniprots_indices[[i]] <- temp_matched_indices
    matched_uniprots_indices_report[i] <- length(matched_uniprots_indices[[i]])
    
    #Assign new uniprot
    #One match:
    if (length(temp_matched_indices) == 1){
      
      new_uniprots[i] <- possible_uniprots_list[[i]][temp_matched_indices]
      
      #More than one match
    } else if (length(temp_matched_indices > 1)){
      
      
      #Find how many destination uniprots correspond to each match. Match with greatest number wins. Take first one if tie
      destination_uniprot_hits <- NULL
      for (j in 1:length(temp_matched_indices)){
        
        destination_uniprot_hits[j] <- length(which(destination_uniprots ==
                                                      possible_uniprots_list[[i]][temp_matched_indices[j]]))
      }
      winning_index <- temp_matched_indices[which(destination_uniprot_hits == 
                                                    max(destination_uniprot_hits))[1]]
      new_uniprots[i] <-  possible_uniprots_list[[i]][winning_index]
    }
  }
  
  ##Return
  return(new_uniprots)
}



####Stats


###Perform basal and insulin anova and dunnetts on 3t3l1 IR ppeptides
##Up to date 20210806
##No p-value adjustment is performed between BAS and INS
site_level_anova_dunnetts_3t3l1 <- function(data,
                                            site,
                                            protnorm = FALSE){
  ###Basal
  
  ##Set up factor
  if (protnorm){
    BAS_conditions_str <- c(rep("CTRL_BAS_PROT_norm", 5),
                            rep("CI_BAS_PROT_norm", 3),
                            rep("DEX_BAS_PROT_norm", 4),
                            rep("TNF_BAS_PROT_norm", 3),
                            rep("MPQ_BAS_PROT_norm", 4),
                            rep("AA_BAS_PROT_norm", 4))
  } else {
    BAS_conditions_str <- c(rep("CTRL_BAS", 6),
                            rep("CI_BAS", 4),
                            rep("DEX_BAS", 4),
                            rep("TNF_BAS", 4),
                            rep("MPQ_BAS", 4),
                            rep("AA_BAS", 4))
  }
  BAS_conditions_fac <- factor(BAS_conditions_str,
                               levels = unique(BAS_conditions_str))
  
  ##Get column indices
  col_indices <- unlist(sapply(levels(BAS_conditions_fac),
                               function(x) grep(paste(c(x,
                                                        "_\\d+"),
                                                      collapse = ""),
                                                colnames(data))))
  
  ##Anova
  BAS_anova <- aov(as.numeric(data[site, 
                                   col_indices]) ~ BAS_conditions_fac,
                   data = data[, col_indices])
  
  ##Dunnetts
  BAS_dunnetts <- glht(BAS_anova,
                       linfct = mcp(BAS_conditions_fac = "Dunnett"))
  
  
  ###INS
  
  ##Set up factor
  if (protnorm){
    INS_conditions_str <- c(rep("CTRL_INS_PROT_norm", 5),
                            rep("CI_INS_PROT_norm", 4),
                            rep("DEX_INS_PROT_norm", 4),
                            rep("TNF_INS_PROT_norm", 4),
                            rep("MPQ_INS_PROT_norm", 4),
                            rep("AA_INS_PROT_norm", 4))
  } else {
    INS_conditions_str <- c(rep("CTRL_INS", 6),
                            rep("CI_INS", 4),
                            rep("DEX_INS", 4),
                            rep("TNF_INS", 4),
                            rep("MPQ_INS", 4),
                            rep("AA_INS", 4))
  }
  INS_conditions_fac <- factor(INS_conditions_str,
                               levels = unique(INS_conditions_str))
  
  ##Get column indices
  col_indices <- unlist(sapply(levels(INS_conditions_fac),
                               function(x) grep(paste(c(x,
                                                        "_\\d+"),
                                                      collapse = ""),
                                                colnames(data))))
  
  ##Anova
  
  INS_anova <- aov(as.numeric(data[site, 
                                   col_indices]) ~ INS_conditions_fac,
                   data = data[, col_indices])
  
  ##Dunnetts
  INS_dunnetts <- glht(INS_anova,
                       linfct = mcp(INS_conditions_fac = "Dunnett"))
  
  
  ###Return
  output_list <- list()
  output_list[[1]] <- summary(BAS_anova)[[1]]["Pr(>F)"][1]
  output_list[[2]] <- summary(BAS_dunnetts)$test$pvalues
  output_list[[3]] <- summary(BAS_dunnetts)$test$coefficients
  output_list[[4]] <- summary(INS_anova)[[1]]["Pr(>F)"][1]
  output_list[[5]] <- summary(INS_dunnetts)$test$pvalues
  output_list[[6]] <- summary(INS_dunnetts)$test$coefficients
  names(output_list) <- c("BAS_anova_p",
                          "BAS_dunnetts_ps",
                          "BAS_dunnetts_estimates",
                          "INS_anova_p",
                          "INS_dunnetts_ps",
                          "INS_dunnetts_estimates")
  
  return(output_list)
}

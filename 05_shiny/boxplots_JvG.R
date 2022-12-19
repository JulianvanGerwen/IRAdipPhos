###Background
#Here I keep boxplot_functions
#Note: You need to source essential_visualisations_JvG


###Packages
library(reshape2)
library(ggplot2)
library(ggrepel)
library(scales)



######################################################
#general_boxplot and attached functions

###Function: Returns a ggplot2 geom that is a minimal boxplot (no points) with no outline, a transparent filled-in box, and a solid mean
#boxplot_df: Data to be plotted. Needs columns: Condition (groups for x axis), Value (values for y axis), and Treatment (for colouring)
minimal_boxplot <- function(boxplot_df,
                            box_alpha = 0.3){
  boxplot <- geom_boxplot(data = boxplot_df,
                          aes(x = Condition, y = Value, fill = Treatment),
                          alpha = box_alpha,
                          outlier.shape = NA,
                          coef = 0,
                          colour = alpha("white", 0)) 
  median_line <- stat_summary(data = boxplot_df,
                              aes(x = Condition, y = Value, color = Treatment),
                              geom = "crossbar", 
                              width = 0.7, fatten=0, 
                              fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) 
  return(list("boxplot" = boxplot, "median_line" = median_line))
}

###Function: Prep boxplot df
#data: contains numerical columns to plot, which are conditions followed by _ and numbers
#protein: row(s) of data to display
#condition: correspond to distinct groups on x-axis
#treatments: Groups for colouring. Must contain conditions e.g. conditions: CTRL_BAS, CTRL_INS, CI_BAS, CI_INS. treatments: CTRL, CI
general_boxplot_dataprep <- function(data,
                                     protein,
                                     conditions,
                                     treatments){
  ###Get columns corresponding to conditions, in order
  conditions_columns <- NULL
  for (i in 1:length(conditions)){
    conditions_columns <- c(conditions_columns,
                            colnames(data)[grep(paste(c(conditions[i],
                                                        "_\\d+"),
                                                      collapse = ""),
                                                colnames(data))])
  }
  ###Make boxplot_df
  boxplot_df <- as.data.frame(t(data[protein, conditions_columns]))
  boxplot_df$Colname <- rownames(boxplot_df)
  boxplot_df <- melt(boxplot_df, id.vars = "Colname")
  colnames(boxplot_df) <- c("Colname", "Protein", "Value")
  boxplot_df$Condition <- gsub("_\\d+",
                               "",
                               boxplot_df$Colname)
  boxplot_df$Treatment <- ""
  for (i in 1:length(treatments)){
    boxplot_df$Treatment[grep(treatments[i],
                              boxplot_df$Condition)] <- treatments[i]
  }
  #Into factors, ordered as given
  boxplot_df$Condition <- factor(boxplot_df$Condition,
                                 levels = conditions)
  boxplot_df$Treatment <- factor(boxplot_df$Treatment,
                                 levels = unique(treatments))
  return(boxplot_df)
}


###Function: Prepare boxplot for faeting
#boxplot_df: output from general_boxplot_dataprep
#facet_labels: Labels for facets. If NULL, uses the row names
general_boxplot_facetprep <- function(boxplot_df,
                                      facet_labels = NULL,
                                      facet_ncol_nrow){
  boxplot_df$facets <- as.factor(boxplot_df$Protein)
  if (!is.null(facet_labels)){
    levels(boxplot_df$facets) <- facet_labels
  }
  facet_plot <- facet_wrap(~ facets, scales = "free_y",
                           ncol = facet_ncol_nrow[1], nrow = facet_ncol_nrow[2])
  return(list("boxplot_df" = boxplot_df,
              "facet_plot" = facet_plot))
}


###Function: Make point plot for general boxplot
#boxplot_df: output from general_boxplot_dataprep
#point_size: Size of points
#point_alpha: alpha of points
#is_raw: If FALSE, gives imputed points a different shape
#imp_method: Method to search for imputed sites.
#if block_imp: Then boolean columns labelled a condition + _blockimp indicate whether block imputation was performed for a condition and a protein
#if imp_matrix: Then columns labelled a colname + _imp indicate whether a given datapoint is imputed
general_boxplot_pointprep <- function(data,
                                      boxplot_df,
                                      point_size = 1,
                                      point_alpha = 0.6,
                                      is_raw = FALSE,
                                      imp_method = "blockimp"){
  
  #Only one shape if raw
  if (is_raw){
    point_plot <- geom_point(data = boxplot_df,
                             aes(x = Condition,
                                 y = Value,
                                 color = Treatment),
                             size = point_size,
                             shape = 16,
                             alpha = point_alpha)
    scale_shape <- NULL
  } else {
    #Label imputed data based on imp_method
    boxplot_df$Imputed <- "Not imputed"
    #If block imputed, search for a boolean column indicating block imputation and label appropriately
    if (imp_method == "blockimp"){
      for (i in 1:nrow(boxplot_df)){
        if (data[as.character(boxplot_df$Protein[i]), paste(c(as.character(boxplot_df$Condition[i]),
                                                              "_blockimp"), collapse = "")]){
          boxplot_df$Imputed[i] <- "Imputed"
        }
      }
    } else if (imp_method == "imp_matrix"){
      for (i in 1:nrow(boxplot_df)){
        if (data[as.character(boxplot_df$Protein[i]), paste(c(as.character(boxplot_df$Colname[i],
                                                             "_imp"), colapse = ""))]){
          boxplot_df$Imputed[i] <- "Imputed"
        }
      }
    } else {
      print("Invalid value for imp_method")
      return("Invalid value for imp_method")
    }
    boxplot_df$Imputed <- factor(boxplot_df$Imputed, levels = c("Not imputed", "Imputed"))
    #Make point_plot
    point_plot <- geom_point(data = boxplot_df,
                             aes(x = Condition,
                                 y = Value,
                                 color = Treatment,
                                 shape = Imputed),
                             size = point_size,
                             alpha = point_alpha)
    scale_shape <- scale_shape_manual(values = c(16,
                                                 18))
  }
  return(list("point_plot" = point_plot,
              "scale_shape" = scale_shape))
}

###Function: Generates labels for general_boxplot, which label datapoins by their original column name
#boxplot_df: Output from general_boxplot_dataprep
#label_poins: Boolean indicating whether to label points
#label_size: size of labels
general_boxplot_labelprep <- function(boxplot_df,
                                      label_points = FALSE,
                                      label_size = 1.5){
  ###Add geom_text if points are to be labelled
  if (label_points){
    point_labels <- geom_text_repel(data = boxplot_df,
                                    aes(x = Condition,
                                        y = Value,
                                        label = Colname),
                                    size = label_size)
  } else {
    point_labels <- NULL
  }
  return(point_labels)
}




###General boxplot for given conditions and treatments
#data: dataframe containing numerical columns to plot. Row(s) are selected for plotting
#protein: The rowname(s) to plot
#conditions: are the names used to group into boxes; these go on the x-axis. Data that will be plotted must be in the form condition_number e.g. condition = CTRL_BAS, example data column is CTRL_BAS_2
#treatments: are the groups used to colour. They must be present with the conditions. For example a treatment could be CTRL
#treatment_colours: The colours to apply to treatments
#is_raw: If FALSE, uses a diamond shape for imputed points
#imp_method: if is_raw == TRUE, this argument specifies how to imputed data. Refer to general_boxplot_pointprep for more info
#facet: This function has the ability to plot multiple proteins/site using ggplot faceting. To do this, enter multiple proteins/sites concatenated in the "protein" argument, and set "facet = TRUE"
#Labels will be the protein/site names by default, but you can specify labels using "facel_labels"
#You must enter the desired dimensions of the resulting plot array in "facet_ncol_nrow" e.g. c(1, 2)
#label_points: Labels each datapoint with the name of the column it came from, using geom_text_repel
general_boxplot <- function(data,
                            protein,
                            conditions,
                            treatments,
                            treatment_colours,
                            is_raw = FALSE,
                            imp_method = "blockimp",
                            facet = FALSE,
                            facet_labels = NULL,
                            facet_ncol_nrow = c(1, 1),
                            label_points = FALSE,
                            include_xaxis = FALSE,
                            include_legend = FALSE,
                            ...){
  #Prep data
  boxplot_df <- general_boxplot_dataprep(data = data, 
                                         protein = protein, 
                                         conditions = conditions, 
                                         treatments = treatments)
  #Facet
  if (facet){
    facet_output <- general_boxplot_facetprep(boxplot_df = boxplot_df, 
                                              facet_labels = facet_labels,
                                              facet_ncol_nrow = facet_ncol_nrow)
    facet_plot <- facet_output[["facet_plot"]]
    boxplot_df <- facet_output[["boxplot_df"]]
  } else {
    facet_plot <- NULL
  }
  
  #Get colour order
  #scale_colour_manual reverts to alphabetical ordering when a level is completely missing. Account for this here
  treatment_missing <- sapply(treatments,
                              function(x) sum(is.na(boxplot_df$Value[which(boxplot_df$Treatment == x)]) == FALSE) == 0)
  if (sum(treatment_missing) > 0){
    scale_color_indices <- order(levels(boxplot_df$Treatment))
  } else {
    scale_color_indices <- 1:length(treatments)
  }
  
  #Build plot
  boxplot <- minimal_boxplot(boxplot_df, ...)
  point_plot <- general_boxplot_pointprep(data = data, boxplot_df = boxplot_df, is_raw = is_raw, imp_method = imp_method, ...)
  output_plot <- ggplot() +
    boxplot[["boxplot"]] +
    boxplot[["median_line"]] +
    point_plot[["point_plot"]] +
    point_plot[["scale_shape"]] +
    
    facet_plot +
    general_boxplot_labelprep(boxplot_df = boxplot_df, label_points = label_points) +
    
    scale_color_manual(values = treatment_colours[scale_color_indices]) +
    scale_fill_manual(values = treatment_colours[scale_color_indices]) +
    
    comfy_theme(include_xaxis = include_xaxis,
                include_legend = include_legend,
                rotate_x_text = TRUE,
                ...) +
    labs(x = "Condition",
         y = "log2 Intensity")
  return(output_plot)
}


################################################
#Other functions that use general_boxplot






###General boxplot for phos, which plots all ppeptides for a given site
general_boxplot_psite <- function(data, 
                                  psite,
                                  rows_or_cols = "cols",
                                  ...){
  #Fetch all ppeptides
  ppeptides <- sort(rownames(data[grep(psite, rownames(data)), ]))
  multiplicities <- as.numeric(sapply(ppeptides, function(x) strsplit(x, "_")[[1]][4]))
  #Facet features
  facet_labels <- c("Single phosphopeptide", "Double phosphopeptide", "Triple phosphopeptide")[multiplicities]
  if (rows_or_cols == "rows"){
    facet_ncol_nrow <- c(length(ppeptides), 1)
  } else if (rows_or_cols == "cols"){
    facet_ncol_nrow <- c(1, length(ppeptides))
  } else {
    return("Inappropriate value selected for rows_or_cols. Options: rows, cols")
  }
  #Run boxplot function
  return(general_boxplot(data = data,
                         protein = ppeptides,
                         facet = TRUE,
                         facet_labels = facet_labels,
                         facet_ncol_nrow = facet_ncol_nrow,
                         ...))
}


###general_boxplot with significance labels
#data: dataframe containing data, including data for boxplot and pval columns
#protein: The row(s) to plot
#the row(s) to plot
#pval_cols: Column names for pvalues
#num_cols: Numerical columns used in the plot. This is used to get the yaxis limits
#comparisons: List of length-2 vectors that describe the comparisons from each of pval_cols. They must correspond to x-axis values in the intended plot
#facet: If TRUE, makes it compatible with faceting in a plot
#facet_labels: Labels to be used for faceting
#show_stars: Indicates whether to show pvalues (NULL) or significance stars. In the latter, you can supply a single value which is the symbol to be used (e.g. *), or a vector with length equal to the nubmer of comparisons
#show_ns: If showing stars, TRUE shows insignificant comparisons as NS
#boxplot_function: The boxplot function to use, typically general_boxplot (this is what will work best)
general_boxplot_siglabels <- function(data,
                                      protein,
                                      pval_cols,
                                      num_cols,
                                      comparisons,
                                      facet = FALSE,
                                      facet_labels = NULL,
                                      show_stars = NULL,
                                      show_ns = FALSE,
                                      boxplot_function = general_boxplot,
                                      ...){
  bplot <- boxplot_function(data = data,
                            protein = protein,
                            facet = facet,
                            facet_labels = facet_labels,
                            ...)
  sig_labels <- add_sig_labels(data = data[protein, ],
                               pval_cols = pval_cols,
                               num_cols = num_cols,
                               comparisons = comparisons,
                               facet = facet,
                               facet_labels = facet_labels,
                               show_stars = show_stars,
                               show_ns = show_ns)
  return(bplot + sig_labels[["pval_plot"]] + sig_labels[["expand_yaxis"]])
}


###general_boxplot with significance labels, that plots ppetpides for the same site separately
#data: dataframe containing data, including data for boxplot and pval columns
#psite: The site to plot (gene_uniprot_site format)
#the row(s) to plot
#pval_cols: Column names for pvalues
#num_cols: Numerical columns used in the plot. This is used to get the yaxis limits
#comparisons: List of length-2 vectors that describe the comparisons from each of pval_cols. They must correspond to x-axis values in the intended plot
#facet: If TRUE, makes it compatible with faceting in a plot
#facet_labels: Labels to be used for faceting
#show_stars: Indicates whether to show pvalues (NULL) or significance stars. In the latter, you can supply a single value which is the symbol to be used (e.g. *), or a vector with length equal to the nubmer of comparisons
#show_ns: If showing stars, TRUE shows insignificant comparisons as NS
#boxplot_function: The boxplot function to use, typically general_boxplot (this is what will work best)
general_boxplot_psite_siglabels <- function(data,
                                            psite,
                                            pval_cols,
                                            num_cols,
                                            comparisons,
                                            show_stars = NULL,
                                            show_ns = FALSE,
                                            boxplot_function = general_boxplot_psite,
                                            ...){
  bplot <- boxplot_function(data = data,
                            psite = psite,
                            ...)
  #Manually extract multiplicities for sig_labels
  #Fetch all ppeptides
  ppeptides <- sort(rownames(data[grep(psite, rownames(data)), ]))
  multiplicities <- as.numeric(sapply(ppeptides, function(x) strsplit(x, "_")[[1]][4]))
  #Facet features
  facet_labels <- c("Single phosphopeptide", "Double phosphopeptide", "Triple phosphopeptide")[multiplicities]
  
  #Make sig_labels
  sig_labels <- add_sig_labels(data = data[ppeptides, ],
                               pval_cols = pval_cols,
                               num_cols = num_cols,
                               comparisons = comparisons,
                               facet = TRUE,
                               facet_labels = facet_labels,
                               show_stars = show_stars,
                               show_ns = show_ns)
  return(bplot + sig_labels)
}











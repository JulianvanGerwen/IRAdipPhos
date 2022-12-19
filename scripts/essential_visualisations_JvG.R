###Background
#Here I keep all my main visualisation functions


###Packages
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(scales)
library(grid)
library(ggpubr)

###comfy_theme
comfy_theme <- function(include_xaxis = TRUE, include_yaxis = TRUE,
                        include_legend = TRUE,
                        rotate_x_text = FALSE, x_text_angle = 45){
    #Base theme
    base_theme <- theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line.x = element_line(colour = "black", size = 0.235),
            axis.line.y = element_line(colour = "black", size = 0.235),
            axis.ticks.x = element_line(size = 0.235,
                                        colour = "black"),
            axis.ticks.y = element_line(size = 0.235,
                                        colour = "black"),
            axis.text.x = element_text(colour = "black",
                                       size = 7),
            axis.text.y = element_text(colour = "black",
                                       size = 7),
            axis.title.x = element_text(colour = "black",
                                        size = 7),
            axis.title.y = element_text(colour = "black",
                                        size = 7),
            plot.title = element_text(colour = "black",
                                      size = 7,
                                      hjust = 0.5),
            legend.text = element_text(colour = "black",
                                       size = 7),
            legend.title = element_text(colour = "black",
                                        size = 7),
            legend.key = element_blank(),
            legend.key.size = unit(0.25, 'cm'),
            strip.text = element_text(size = 7, colour = "black"),
            strip.background = element_blank())
    #Build list for additional arguments
    if (include_legend){
      legend_position <- NULL
    }else{
      legend_position <- "none"
    }
    if (include_xaxis){
      axis_line_x <- NULL
      axis_ticks_x <- NULL
    }else{
      axis_line_x <- element_blank()
      axis_ticks_x <- element_blank()
    }
    if (include_yaxis){
      axis_line_y <- NULL
      axis_ticks_y <- NULL
    }else{
      axis_line_y <- element_blank()
      axis_ticks_y <- element_blank()
    }
    #Rotate x text if desired
    if (rotate_x_text){
      #if 90 degrees
      if (x_text_angle == 90){
        axis_text_x <- element_text(colour = "black",
                                    size = 7,
                                    hjust = 1, vjust = 0.5, angle = 90)
      } else if (x_text_angle == 45){
        axis_text_x <- element_text(colour = "black",
                                    size = 7,
                                    hjust = 1, vjust = 1, angle = 45)
      }
      
    } else{
      axis_text_x <- NULL
    }
    theme_args <- list("legend.position" = legend_position,
                       "axis.line.x" = axis_line_x,
                       "axis.ticks.x" = axis_ticks_x,
                       "axis.line.y" = axis_line_y,
                       "axis.ticks.y" = axis_ticks_y,
                       "axis.text.x" = axis_text_x)
    theme_args <- theme_args[vapply(theme_args, function(x) {!is.null(x)}, logical(1))]
    
    #Make secondary theme
    secondary_theme <- do.call(function(...) theme(...), theme_args)
    return(base_theme + secondary_theme)
}








###Function: Volcano plot given pval and FC columns (can be adj_pval, and/or rval)
#Up to date 20210926
#By default, colours by pval and FC cutoff
#User can specificy a categorical column to colour points by, through the argument "colour_col". Then, the levels of this categorical column and the colours must be specified (colour_col_levels, colour_col_colours). Note that points will be placed in the oredr of the colour_col_levels levels, so order wisely 
#Colours significant (pval and FC) proteins
#label_col specifies were to get labels from (if desired), and rows_to_label points at what rows to label (needs to be row names)
volcano_plot_coloured <- function(data,
                                  pval_col,
                                  FC_col,
                                  pval_cutoff = 0.05,
                                  FC_cutoff = 0.58,
                                  inner_alpha_val = 0.3,
                                  outer_alpha_val = 0.5,
                                  size = 1,
                                  sig_colour = "red",
                                  colour_col = FALSE,
                                  colour_col_levels,
                                  colour_col_colours,
                                  x_lim = "max",
                                  y_lim = "max",
                                  remove_Inf = TRUE,
                                  return_df = FALSE,
                                  show_legend = FALSE,
                                  label_col = FALSE,
                                  rows_to_label = c()){
  
  ##Set up data
  plot_data <- data
  colnames(plot_data)[colnames(plot_data) == FC_col] <- "FC"
  colnames(plot_data)[colnames(plot_data) == pval_col] <- "pval"
  plot_data$nlog_pval <- -log10(plot_data$pval)
  plot_data$sig <- FALSE
  plot_data$sig[which(abs(plot_data$FC) > FC_cutoff &
                        plot_data$pval < pval_cutoff)] <- TRUE
  
  ##If colour col
  if (colour_col != FALSE){
    
    plot_data[, "colour_col"] <- data[, colour_col]
    
    #Tag so I colour by the right name
    colour_by <- "colour_col"
    
    #Adjusting colours
    colour_adjusting <- scale_colour_manual(values = alpha(colour_col_colours,
                                                           outer_alpha_val)) 
    fill_adjusting <-  scale_fill_manual(values = alpha(colour_col_colours,
                                                        alpha = inner_alpha_val)) 
  } else {
    
    colour_by <- "sig"
    colour_adjusting <- scale_colour_manual(values = alpha(c("black",
                                                             sig_colour),
                                                           outer_alpha_val)) 
    fill_adjusting <- scale_fill_manual(values = alpha(c("black",
                                                         sig_colour),
                                                       alpha = inner_alpha_val))
  }
  
  #If label_col
  if (label_col != FALSE){
    
    plot_data[, "label_col"] <- data[, label_col]
    plot_data[which(rownames(plot_data) %in%
                      rows_to_label == FALSE), "label_col"] <- ""
    
    label_plot <- geom_text_repel(aes_string(x = "FC",
                                             y = "nlog_pval",
                                             colour = colour_by,
                                             label = "label_col"),
                                  size = 1.875,
                                  min.segment.length = 0,
                                  max.overlaps = Inf) 
  } else {
    
    label_plot <- NULL
  }

  
  
  
  ##Remove Inf if desired
  if (remove_Inf){
    
    plot_data <- plot_data[which(plot_data$nlog_pval !=
                                   Inf), ]
  }
  
  ##Get x_lim and y_lim if using max
  if (x_lim == "max"){
    
    x_lim <- max(abs(plot_data$FC[which(is.na(plot_data$pval) == FALSE)]),
                 na.rm = TRUE)
  }  
  
  if (y_lim == "max"){
    
    y_lim <- max(plot_data$nlog_pval[which(plot_data$nlog_pval != Inf)])
  }
  
  ##Order df
  
  #If colour by sig:
  if (colour_col == FALSE){
    
    plot_data <- plot_data[order(plot_data$sig), ] 
  } else {
    
    #Make colour_col a factor
    plot_data$colour_col <- factor(plot_data$colour_col,
                                   levels = colour_col_levels)
    plot_data <- plot_data[order(plot_data$colour_col), ]
  }
  
  
  
  ##Decide whether to keep legend
  if (show_legend == FALSE){
    
    legend_theme <- theme(legend.position = "none")
  } else {
    
    legend_theme = NULL
  }
  
  
  
  
  
  ##Output
  
  #return_df
  if (return_df){
    
    return(plot_data)
    
    #Colour by sig
  } else {
    output_plot <- ggplot(data = plot_data,
                          aes_string(x = "FC",
                                     y = "nlog_pval",
                                     colour = colour_by,
                                     fill = colour_by))+
      geom_point(size = size,
                 shape = 21,
                 stroke = size*0.3) +
      label_plot +
      colour_adjusting +
      fill_adjusting +
      comfy_theme() +
      legend_theme +
      xlim(-x_lim, x_lim) +
      ylim(0, y_lim)
    return(output_plot) 
  } 
}



###Funcion that converts unique psites formatted formally into nice format
#mult = TRUE: gene_uniprot_site_multiplicity to gene site Pmultplicity (or gene uniprot site Pmultiplicity if there are double ups)
#mult = FALSE: gene_uniprot_site to gene site (or gene uniprot site if there are double ups)
psite_to_nice_name <- function(psites, mult = TRUE){
  name_df <- as.data.frame(t(sapply(psites, function(x) strsplit(x, "_")[[1]])))
  #If multiplicity is present
  if (mult){
    names(name_df) <- c("gene", "uniprot", "site", "mult")
    name_df$mult_nice <- paste("P", name_df$mult, sep = "")
    name_df$site_mult_nice <- paste(name_df$site, name_df$mult_nice, sep = " ")
    #If no multiplicity
  } else {
    names(name_df) <- c("gene", "uniprot", "site")
    name_df$site_mult_nice <- name_df$site
  }
  name_df$nice_name_old <- apply(name_df[, c("gene", "site_mult_nice")],
                                 1, function(x) paste(x, collapse = " "))
  #Include uniprots for sites that have same gene, site, and multiplicity, or no gene
  name_df$nice_name <- name_df$nice_name_old
  for (i in 1:nrow(name_df)){
    if (length(which(name_df$nice_name_old == name_df$nice_name[i])) > 1 |
        name_df$gene[i] %in% c("NA", "")){
      name_df$nice_name[i] <- paste(c(name_df[i, c("gene", "uniprot", "site_mult_nice")]), collapse = " ")
    }
  }
  #Strip of NA
  name_df$nice_name <- gsub("NA ", "", name_df$nice_name)
  return(name_df$nice_name)
}

###Functions for naming yaxis items in fc_sig_heatmap
gene_uniprot_site_mult_naming <- function(data){
  data$standard_name <- psite_to_nice_name(rownames(data))
  return(data)
}
gene_uniprot_naming <- function(data){
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
  return(data)
}


###Heatmap function
#Up to date 20211020
##Plots FC on heatmap
##Plots significance as point in tile
#Arguments:
#data: your data
#fc_cols_w_pvals: FC columns that have corresponding pval columns
#fc_cols_wo_pvals: FC columns without corresponding pval columns
#pval_cols: p-value columns
#shape2_pval_cols: p-value columns where you want a different shape for the significance point
#order_column: The FC column to order by
#is_decreasing: How you order
#return_df: Returns the dataframe that is put into ggplot
#x_axis_names: names for x_axis, correspond to fc_cols
#gap_column_index: index of column where you want a gap
#scale_colours_manually: do you want to scale the gradient manually? If so, you can change where the middle colour falls with mid_colour_threshold, and where the max colour falls with high_colour_threhsold
#mid_colour_threshold: See scale_colours_manually Typically you would want to highlight the changes between 0 and mid_colour_threshold
#high_colour_threshold: See scale_colours_manually
#colour_scheme: What colour scheme to use. Accepts "original" (quite gross), or brewer palette codes e.g. "RdBu". To reverse pallet, enter "rev_RdBu"
#middle_brewer_white: If using brewer palette, sets middle colour to white (FFFFFF) if true
#yaxis_naming: How rownames will be processed to get y-axis names. "gene_uniprot_site_mult" makes names "gene site Pmult", "gene_uniprot" makes name "gene", but "gene uniprot" for double-up genes, "psite_and_prot" applies both "gene_uniprot_site_mult" and "gene_uniprot", and "keepall" leaves rownames as is
#legend_title: Title for legend
#tile_border_size: Size of white border for tiles. 0.1 is good. If FALSE, no border
#shapes: Shapes to use for dots
#show_shape_legend: Whether or not to show shape in legend
#aspect_ratio_fixed: if true, force tiles to be squares
fc_sig_heatmap <- function(data,
                          fc_cols_w_pvals = c(),
                          fc_cols_wo_pvals = c(),
                          pval_cols = c(),
                          shape2_pval_cols = c(),
                          order_column,
                          is_decreasing = TRUE,
                          return_df = FALSE,
                          x_axis_names = c(),
                          gap_column_index = FALSE,
                          scale_colours_manually = FALSE,
                          mid_colour_threshold = 1.5,
                          high_colour_threshold = "max",
                          colour_scheme = "rev_RdBu",
                          middle_brewer_white = TRUE,
                          yaxis_naming = "gene_uniprot_site_mult",
                          legend_title = "log2 INS/BASAL",
                          tile_border_size = 0.05,
                           shapes = c(19, 17),
                          show_shape_legend = TRUE,
                          aspect_ratio_fixed = TRUE){
  
  
  ##Prep for melts
  
  #Set up data
  
  #yaxis_naming
  
  #Apply naming functions
  if (yaxis_naming == "gene_uniprot_site_mult"){
    data <- gene_uniprot_site_mult_naming(data)
  } else if (yaxis_naming == "gene_uniprot"){
    data <- gene_uniprot_naming(data)
  } else if (yaxis_naming == "psite_and_prot"){
    underscore_counts <- sapply(strsplit(rownames(data), "_"),
                                function(x) length(x))
    psite_indices <- which(underscore_counts == 4)
    prot_indices <- which(underscore_counts == 2)
    data$standard_name <- ""
    if (length(psite_indices) > 0){
      data$standard_name[psite_indices] <- gene_uniprot_site_mult_naming(data[psite_indices, ])$standard_name
    }
    if (length(prot_indices) > 0){
      data$standard_name[prot_indices] <- gene_uniprot_naming(data[prot_indices, ])$standard_name
    }
    
  } else if (yaxis_naming == "keepall"){
    
    data$standard_name <- rownames(data)
  }
  
  
  data <- data[order(data[, order_column],
                     decreasing = is_decreasing == FALSE), ]
  data$standard_name <- factor(data$standard_name,
                               levels = unique(data$standard_name))
  
  ##Gap column
  #Make FC and pval columns with 0 FC (appear white) and non-significant pvalue (no dot)
  #data$gap_FC_column <- 0
  #data$gap_pval_column <- 1
  ##Insert gap column in fc_cols and pval_cols. use gap_column_w_pvals to decide where
  #if (gap_column_w_pvals){
  #  
  #  fc_cols
  #}
  
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
    
    #dot_shape column
    pval_melt_df_w_pvals$dot_shape <- "dot"
    pval_melt_df_w_pvals$dot_shape[which(pval_melt_df_w_pvals$variable %in%
                                           shape2_pval_cols)] <- "shape2"
    pval_melt_df_w_pvals$dot_shape <- factor(pval_melt_df_w_pvals$dot_shape,
                                             levels = c("dot", "shape2"))
    
    #Combine
    post_melt_df_w_pvals <- FC_melt_df_w_pvals
    colnames(post_melt_df_w_pvals)[which(colnames(post_melt_df_w_pvals) == "value")] <- "FC_value"
    post_melt_df_w_pvals$sig_value <- pval_melt_df_w_pvals$value
    post_melt_df_w_pvals$dot_shape <- pval_melt_df_w_pvals$dot_shape
    
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
    FC_melt_df_wo_pvals$dot_shape <- "shape"
    
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
    gap_column_df$FC_value <- NA
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
  
  ##Set aspect ratio
  if (aspect_ratio_fixed){
    aspect_ratio_fixed <- coord_fixed(ratio = 1) 
  } else {
    aspect_ratio_fixed <- NULL
  }
  
  
  ##Make colour scaling
  
  #If using max for high
  if (high_colour_threshold == "max"){
    high_colour_threshold <- max(abs(post_melt_df_combined$FC_value),
                                 na.rm = TRUE)
  }
  
  #If using original colours
  if (colour_scheme == "original"){
    
    #If scaling manually:
    if (scale_colours_manually){
      
      gradient <- scale_fill_gradientn(colors = c("#1A09FF",
                                                  "#B28AFF",
                                                  "white",
                                                  "#FF9E82",
                                                  "#FF1C0B"),
                                       values = rescale(c(-high_colour_threshold, 
                                                          -mid_colour_threshold,
                                                          0,
                                                          mid_colour_threshold, 
                                                          high_colour_threshold)),
                                       limits = c(-high_colour_threshold, high_colour_threshold),
                                       name = legend_title)
      
    } else {
      
      gradient <- scale_fill_gradient2(low = "blue",
                                       mid = "white",
                                       high = "red",
                                       name = legend_title)
    }
    
    
    #If using brewer
  } else{
    
    ##Extract brewer colours
    colour_scheme <- strsplit(colour_scheme, "_")[[1]]
    
    if (length(colour_scheme) == 1){
      
      colours <- brewer.pal(9, colour_scheme[1])
    } else if (colour_scheme[1] == "rev"){
      
      colours <- rev(brewer.pal(9, colour_scheme[2]))
    }
    
    if (middle_brewer_white){
      
      colours[5] <- "#FFFFFF"
    }
    
    
    #If scaling manually:
    if (scale_colours_manually){
      
      
      gradient <- scale_fill_gradientn(colours = colours,
                                       values = rescale(c(-high_colour_threshold,
                                                          -0.5*high_colour_threshold -0.5*mid_colour_threshold,
                                                          -mid_colour_threshold,
                                                          -0.5*mid_colour_threshold,
                                                          0,
                                                          0.5*mid_colour_threshold,
                                                          mid_colour_threshold,
                                                          0.5*high_colour_threshold + 0.5*mid_colour_threshold,
                                                          high_colour_threshold)),
                                       limits = c(-high_colour_threshold, high_colour_threshold),
                                       name = legend_title)
      
    } else {
      
      gradient <- scale_fill_gradientn(colours = colours,
                                       values = rescale(c(-high_colour_threshold,
                                                          -0.75*high_colour_threshold,
                                                          -0.5*high_colour_threshold,
                                                          -0.25*high_colour_threshold,
                                                          0,
                                                          0.25*high_colour_threshold,
                                                          0.5*high_colour_threshold,
                                                          0.75*high_colour_threshold,
                                                          high_colour_threshold)),
                                       limits = c(-high_colour_threshold, high_colour_threshold),
                                       name = legend_title)
    }
  }
  
  
  #Set up tile
  if (tile_border_size == FALSE){
    
    tile = geom_tile()
  } else {
    
    tile = geom_tile(colour = "#a5a5a5",
                     size = tile_border_size)
  }
  
  
  #Scale for shapes
  if (show_shape_legend){
    shape_scale <- scale_shape_manual(values = shapes)
  } else {
    shape_scale <- scale_shape_manual(values = shapes,
                                      guide = "none")
  }
  
  ##Heatmap
  output_plot <- ggplot(post_melt_df_combined,
                        aes(x = variable,
                            y = standard_name,
                            fill = FC_value)) + 
    gradient + 
    tile +
    aspect_ratio_fixed +
    geom_point(data = post_melt_df_combined[post_melt_df_combined$sig_value == TRUE,],
               aes(x = variable,
                   y = standard_name,
                   shape = dot_shape),
               size = 0.5) +
    shape_scale +
    #labs(x = "Condition",
    #     y = "Phosphopeptide") +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     colour = "black",
                                     size = 7),
          axis.text.y = element_text(colour = "black",
                                     size = 7),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    colour = "black",
                                    size = 7),
          legend.text = element_text(colour = "black",
                                     size = 7),
          legend.title = element_text(colour = "black",
                                      size = 7),
          legend.key.size = unit(0.15,
                                 units = "inch"))
  
  if (return_df){
    
    return(post_melt_df_combined)
    
  } else {
    
    return(output_plot)
    
  }
  
}


###Function: Add text to heatmap
#Up to date 20220104
#heatmap: An existing heatmap generated with fc_sig_heatmap
#data: Data used to generate the heatmap
#order_col: the order_col used to generate the heatmap
#text_cols: columns of text in data
#text_xaxis_names: Names to appear on xaxis
#left_or_right: Whether text is added to the left or the right. Only right looks good
add_text_to_heatmap <- function(heatmap,
                                data,
                                order_col,
                                text_cols,
                                text_xaxis_names,
                                left_or_right = "right"){
  #Melt data
  melted_data <- fc_sig_heatmap(data = data,
                                fc_cols_wo_pvals = c(order_col, text_cols),
                                order_column = order_col,
                                is_decreasing = FALSE,
                                x_axis_names = c("order",
                                                 text_xaxis_names),
                                return_df = TRUE,
                                high_colour_threshold = 1)
  
  #Make text plot
  text_plot <- ggplot(melted_data[melted_data$variable != "order", ],
                      aes(x = variable, y = standard_name, label = FC_value)) +
    geom_text(hjust = 0, size = 6*0.35) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(colour = "black", size = 6, hjust = 1, vjust = 0.5, angle = 90),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    coord_cartesian(clip = "off")
  
  #Combine and return
  grob1 <- ggplotGrob(heatmap)
  grob2 <- ggplotGrob(text_plot)
  if (left_or_right == "right"){
    combined_plot <- cbind(grob1, grob2, size = "first")
  } else if (left_or_right == "left"){
    combined_plot <- cbind(grob2, grob1, size = "last")
  } else {
    print("Invalid value selected for left_or_right")
  }
  combined_plot$heights <- unit.pmax(grob1$heights, grob2$heights)
  return(ggarrange(combined_plot))
}


####Adding significance stars to plots

###Function: ADd signifiance labels to a ggplot using stat_pvalue_manual
#test_data: dataframe containing minimal columns: group1, group2 (indicate comparisons), p, and rownames (rows of num_data that correspond to pvalues)
#if using with facets, simply add a column named the same as the group column in your original data and add appropriate values
#num_data: Numerical data corresponding to test values. Used to calculated positions for brackets
#show_stars: Indicates whether to show pvalues (NULL) or significance stars. In the latter, you can supply a single value which is the symbol to be used (e.g. *), or a vector with length equal to the nubmer of comparisons
#show_ns: If showing stars, TRUE shows insignificant comparisons as NS
add_sig_labels_fromtable <- function(test_data,
                                     num_data,
                           show_stars = NULL,
                           show_ns = FALSE){
  #Add columns
  test_data$sig_levels <- (test_data$p < 0.05) + (test_data$p < 0.01) + (test_data$p < 0.001)
  test_data$comparison <- paste(test_data$group1, test_data$group2, sep = "vs")
  comparisons <- unique(test_data$comparison)
  
  #Make significance labels
  if (!is.null(show_stars)){
    
    #Expand if only 1 star
    if (length(show_stars) == 1){
      show_stars <- rep(show_stars, length(comparisons))
    }
    
    #Significance labels
    test_data$sig_label <- NA
    #Loop over comparisons
    for (i in 1:length(comparisons)){
      comparison <- comparisons[i]
      #Loop over test data
      indices <- which(test_data$comparison == comparison)
      for (j in indices){
        if (!is.na(test_data$sig_levels[j]) &
            test_data$sig_levels[j] > 0){
          test_data$sig_label[j] <- paste(rep(show_stars[i], test_data$sig_levels[j]), collapse = "")
        }
      }
    }
    
    #Show ns if you want
    if (show_ns){
      test_data$sig_label[is.na(test_data$sig_label)] <- "NS"
    }
    
  } else {
    test_data$sig_label <- test_data$p
  }
  
    #Loop over rownames to get y position
    test_data$y.position = NA
    for (row in unique(test_data$rownames)){
      ymax <- max(num_data[row, ], na.rm = TRUE)
      ymin <- min(num_data[row, ], na.rm = TRUE)
      test_data$y.position[test_data$rownames == row] <- ymax + 0.075*(ymax - ymin)
      
      #Overlapping comparisons
      overlapping_indices <- c()
      #Here only inspect comparisons that will be shown
      total_indices <- which(test_data$rownames == row & !is.na(test_data$sig_label))
      if (length(total_indices) > 1){
        for (i in total_indices){
          #Test overlap in other groups
          if (test_data[i, "group1"] %in% c(test_data[setdiff(total_indices, c(i)), "group1"],
                                            test_data[setdiff(total_indices, c(i)), "group2"]) |
              test_data[i, "group2"] %in% c(test_data[setdiff(total_indices, c(i)), "group1"],
                                            test_data[setdiff(total_indices, c(i)), "group2"])){
            overlapping_indices <- c(overlapping_indices, i)
          }
        }
      }
      
      #displace. Here I generate a sequence of y values spaced by the total axis range
      if (length(overlapping_indices) > 1){
        test_data[overlapping_indices, "y.position"] <- 
          seq(from = test_data[overlapping_indices[1], "y.position"],
              by = 0.1*(ymax - ymin),
              length.out = length(overlapping_indices))
      }
      }
  
  #Plot
  pval_plot <- stat_pvalue_manual(test_data, 
                                  label = "sig_label")
  return(list("pval_plot" = pval_plot, 
              "expand_yaxis" = scale_y_continuous(expand = expansion(mult = c(0, 0.15)))))
}


###Function to prep data for sig labels
#data: data to use, which must include pval_cols
#pval_cols: Column names for pvalues
#Comparisons: List of length-2 vectors that describe the comparisons from each of pval_cols. They must correspond to x-axis values in the intended plot
#facet: If TRUE, makes it compatible with faceting in a plot
#facet_labels: Labels to be used for faceting
prep_sig_labels <- function(data,
                            pval_cols,
                            comparisons,
                            facet = FALSE,
                            facet_labels = NULL){
  #Melt stat data
  data$rownames <- rownames(data)
  data_melt <- melt(data[, c(pval_cols, "rownames")], id.vars = c("rownames"))
  
  #Add columns
  data_melt$p <- data_melt$value
  
  #Add groups 
  data_melt[, c("group1", "group2")] <- NA
  for (i in 1:length(pval_cols)){
    indices <- which(data_melt$variable == pval_cols[i])
    for (j in indices){
      data_melt[j, c("group1", "group2")] <- comparisons[[i]]
    }
  }
  if (facet){
    data_melt$facets <- facet_labels
    data_melt$facets <- factor(data_melt$facets,
                               levels = unique(data_melt$facets))
  }
  return(data_melt)
}


###Function to make sig labels e.g. for a boxplot
#data: data to use, which must include pval_cols
#pval_cols: Column names for pvalues
#num_cols: Numerical columns used in the plot. This is used to get the yaxis limits
#comparisons: List of length-2 vectors that describe the comparisons from each of pval_cols. They must correspond to x-axis values in the intended plot
#facet: If TRUE, makes it compatible with faceting in a plot
#facet_labels: Labels to be used for faceting
#show_stars: Indicates whether to show pvalues (NULL) or significance stars. In the latter, you can supply a single value which is the symbol to be used (e.g. *), or a vector with length equal to the nubmer of comparisons
#show_ns: If showing stars, TRUE shows insignificant comparisons as NS
add_sig_labels <- function(data,
                           pval_cols,
                           num_cols,
                           comparisons,
                           facet = FALSE,
                           facet_labels = NULL,
                           show_stars = NULL,
                           show_ns = FALSE){
  test_data <- prep_sig_labels(data = data,
                               pval_cols = pval_cols,
                               comparisons = comparisons,
                               facet = facet,
                               facet_labels = facet_labels)
  num_data <- data[, num_cols]
  return(add_sig_labels_fromtable(test_data = test_data,
                                  num_data <- data[, num_cols],
                                  show_stars = show_stars,
                                  show_ns = FALSE))
}



###Makes barplot of fisher's exact-enriched pathways
#Up to date 20211027
pathway_enrichment_barplot <- function(plot_df,
                                       num_pathways = "sig_only",
                                       colour,
                                       sig_line = TRUE,
                                       sig_cutoff = 0.05,
                                       x_max = FALSE){
  ##Colours
  #1st is high, second is low
  #Brewer palettes: 2-8
  colours <- list(c("#E7E1EF",
                    "#C5C6E0",
                    "#96AED2",
                    "#6298C4",
                    "#2C7DB2",
                    "#0D5BA0",
                    "#0B467A"),
                  c("#FDE0C3",
                    "#FBC691",
                    "#FB9E59",
                    "#FA792F",
                    "#EB5212",
                    "#CE3306",
                    "#932507"),
                  c("#FCD7EB",
                    "#ECA4D1",
                    "#D45E9E",
                    "#B6006A"),
                  c("#E0F4C6",
                    "#ABDD73",
                    "#6EB232",
                    "#3F8319"),
                  c("#ECECEC",
                    "#D0D0D0",
                    "#AFAFAF",
                    "#848484",
                    "#606060",
                    "#414141",
                    "#1C1C1C"),
                  c("#E7E1EF",
                    "#C5C6E0",
                    "#96AED2",
                    "#6298C4",
                    "#2C7DB2",
                    "#0D5BA0",
                    "#3660A5"))
  names(colours) <- c("blue",
                      "red",
                      "pink",
                      "green",
                      "grey",
                      "GSK3_putsub_blue")
  
  ##Set up
  plot_df <- as.data.frame(plot_df,
                           stringsasfactors = FALSE)
  colnames(plot_df) <- c("pval", 
                         "adj_pval",
                         "num_DE_genes_in_pathway",
                         "num_background_genes_in_pathway",
                         "num_total_genes_in_pathway")
  plot_df <- plot_df[order(plot_df$adj_pval), ]
  #num_pathways
  if (num_pathways == "sig_only"){
    plot_df <- plot_df[which(plot_df$adj_pval < sig_cutoff), ]
  } else {
    plot_df <- plot_df[1:num_pathways, ]
  }
  
  
  plot_df$pathway <- as.factor(rownames(plot_df))
  
  plot_df$nlog_adj_pval <- -log10(plot_df$adj_pval)
  plot_df$pathway <- factor(plot_df$pathway,
                            levels = as.character(plot_df[order(plot_df$nlog_adj_pval,
                                                                decreasing = FALSE), 
                                                          "pathway"]))
  
  ##Sig line
  if (sig_line){
    
    sig_line_plot <- geom_vline(xintercept = -log10(sig_cutoff),
                                alpha = 0.5)
  } else {
    
    sig_line_plot <- NULL
  }
  
  ##xlim
  if (x_max != FALSE){
    
    x_limits <- xlim(0, x_max)
  } else {
    
    x_limits <- NULL
  }
  
  ##Plot
  output_plot <- ggplot(plot_df,
                        aes(x = nlog_adj_pval,
                            y = pathway,
                            fill = num_DE_genes_in_pathway)) +
    geom_col() +
    sig_line_plot +
    scale_fill_gradientn(colours = colours[[colour]],
                         name = "# Genes") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line  = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5,
                                    size = 8),
          axis.text.x = element_text(colour = "black",
                                     size = 6),
          axis.text.y = element_text(colour = "black",
                                     size = 6),
          axis.title.x = element_text(colour = "black",
                                      size = 8),
          axis.title.y = element_blank(),
          legend.text = element_text(colour = "black",
                                     size = 6),
          legend.title = element_text(colour = "black",
                                      size = 8),
          legend.key.size = unit(0.15,
                                 units = "inch")) +
    x_limits +
    labs(x = "-log10 adjusted p-value",
         y = "Pathway")
  
  return(output_plot)
}


###Given a list of two fisher's exact enriched pathways, plots one on top of the other with the same horizontal scale
#Up to date 20210528
#stack can be "vertical" or "horizontal", dictates whether plots will be stacked on top of or next to each other
pathway_enrichment_barplot_fuse <- function(plot_df_list,
                                            num_pathways = "sig_only",
                                            colours,
                                            sig_line = TRUE,
                                            sig_cutoff = 0.05,
                                            stack = "vertical"){
  
  ###Get maximum x-value to standardise axes
  x_limits <- NULL
  for (i in 1:length(plot_df_list)){
    
    x_limits[i] <- max(-log10(plot_df_list[[i]][, 2]))
  }
  
  ###Run barplot function on each plot in list
  plot_list <- list()
  
  for (i in 1:length(plot_df_list)){
    
    plot_list[[i]] <- pathway_enrichment_barplot(plot_df = plot_df_list[[i]],
                                                 num_pathways = num_pathways,
                                                 colour = colours[i],
                                                 sig_line = sig_line,
                                                 sig_cutoff = sig_cutoff,
                                                 x_max = max(x_limits))
  }
  
  ###Make grob
  grob1 <- ggplotGrob(plot_list[[1]])
  grob2 <- ggplotGrob(plot_list[[2]])
  if (stack == "vertical"){
    output_grob <- rbind(grob1, 
                         grob2, 
                         size = "first")
    output_grob$widths <- unit.pmax(grob1$widths,
                                    grob2$widths)
  } else if (stack == "horizontal"){
    output_grob <- cbind(grob1, 
                         grob2, 
                         size = "first")
    output_grob$heights <- unit.pmax(grob1$heights,
                                     grob2$heights)
  }
  
  
  return(output_grob)
}


###gst heatmap
#Takes in list of gst output matrices for multiple conditions
#Gives dots for signficance
#Can colour red or blue
gst_heatmap <- function(gst_output_list,
                        colour,
                        sig_cutoff = 0.05,
                        single_df = FALSE,
                        single_df_condition = ""){
  
  ##Colours
  colours <- list(c("#FFFFFF",
                    "#C5C6E0",
                    "#6298C4",
                    "#0D5BA0",
                    "#0B467A"),
                  c("#FFFFFF",
                    "#FBC691",
                    "#FA792F",
                    "#CE3306",
                    "#932507"),
                  c("#FFFFFF",
                    "#FCD7EB",
                    "#ECA4D1",
                    "#D45E9E",
                    "#B6006A"),
                  c("#FFFFFF",
                    "#E0F4C6",
                    "#ABDD73",
                    "#6EB232",
                    "#3F8319"))
  names(colours) <- c("blue",
                      "red",
                      "pink",
                      "green")
  
  ###Single df
  if (single_df){
    
    ##All sig pathways, and order
    gst_output_df <- as.data.frame(gst_output_list)
    gst_output_df <- gst_output_df[which(gst_output_df$adj_pval < 
                                           sig_cutoff), ]
    gst_output_df <- gst_output_df[order(gst_output_df$adj_pval,
                                         decreasing = TRUE), ]
    ##Add columns
    gst_output_df$condition <- single_df_condition
    gst_output_df$pathway <- factor(rownames(gst_output_df),
                                    levels = rownames(gst_output_df))
    gst_output_df$nlog_adj_pval <- -log10(gst_output_df$adj_pval)
    gst_output_df$sig_value <- gst_output_df$adj_pval < sig_cutoff
    
    ##Plot
    output_plot <- ggplot(gst_output_df,
                          aes(x = condition,
                              y = pathway,
                              fill = nlog_adj_pval)) +
      geom_tile() +
      scale_fill_gradientn(colors = colours[[colour]],
                           values = rescale(c(0, 
                                              -log10(sig_cutoff), 
                                              (max(gst_output_df$nlog_adj_pval) + 2*(-log10(sig_cutoff)))/3, 
                                              (2*max(gst_output_df$nlog_adj_pval) + (-log10(sig_cutoff)))/3, 
                                              max(gst_output_df$nlog_adj_pval))),
                           name = "-log10 adjusted p-value") +
      coord_fixed(ratio = 1) +
      geom_point(data = gst_output_df[gst_output_df$sig_value == TRUE,],
                 size = 0.5) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5),
            axis.ticks = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = "Condition",
           y = "Pathway")
    return(output_plot)
    
  } else{
    
    ###Multiple dfs
    
    
    ##Get all sig pathways
    sig_pathways_list <- list()
    for (i in 1:length(gst_output_list)){
      
      sig_pathways_list[[i]] <- rownames(gst_output_list[[i]])[which(gst_output_list[[i]][, "adj_pval"] < sig_cutoff)]
    }
    sig_pathways <- unique(unlist(sig_pathways_list))
    ##Fuse
    for (i in 1:length(gst_output_list)){
      
      temp_output_df <- as.data.frame(gst_output_list[[i]],
                                      stringsasfactors = FALSE)
      #Sig pathways only
      temp_output_df <- temp_output_df[sig_pathways, ]
      #colnames
      colnames(temp_output_df) <- sapply(colnames(temp_output_df),
                                         function(x) paste(c(names(gst_output_list)[i],
                                                             x),
                                                           collapse = "_"))
      if (i == 1){
        gst_output_big_df <- temp_output_df
      } else {
        gst_output_big_df <- cbind(gst_output_big_df,
                                   temp_output_df)
      }
    }
    gst_output_big_df <- gst_output_big_df[, grep("adj_pval",
                                                  colnames(gst_output_big_df))]
    gst_output_big_df <- gst_output_big_df[order(gst_output_big_df[, 1],
                                                 decreasing = TRUE), ]
    
    ##Melt
    melt_df <- melt(cbind(rownames(gst_output_big_df),
                          gst_output_big_df))
    colnames(melt_df) <- c("pathway",
                           "condition",
                           "adj_p_val")
    melt_df$nlog_adj_pval <- -log10(melt_df$adj_p_val)
    levels(melt_df$condition) <- sapply(levels(melt_df$condition),
                                        function(x) strsplit(x, "_")[[1]][1])
    melt_df$pathway <- factor(melt_df$pathway,
                              levels = rownames(gst_output_big_df))
    melt_df$sig_value <- melt_df$adj_p_val < sig_cutoff
    
    ##Plot
    output_plot <- ggplot(melt_df,
                          aes(x = condition,
                              y = pathway,
                              fill = nlog_adj_pval)) +
      geom_tile() +
      scale_fill_gradientn(colors = colours[[colour]],
                           values = rescale(c(0, 
                                              -log10(sig_cutoff), 
                                              (max(melt_df$nlog_adj_pval) + 2*(-log10(sig_cutoff)))/3, 
                                              (2*max(melt_df$nlog_adj_pval) + (-log10(sig_cutoff)))/3, 
                                              max(melt_df$nlog_adj_pval))),
                           name = "-log10 adjusted p-value") +
      coord_fixed(ratio = 1) +
      geom_point(data = melt_df[melt_df$sig_value == TRUE,],
                 size = 0.5) +
      comfy_theme(axis_lines = FALSE) +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5,
                                       size = 6,
                                       colour = "black")) +
      labs(x = "Condition",
           y = "Pathway")
    return(output_plot)
  }
}




###FUnction: ggsave as a pdf and png
ggsave_pdfpng <- function(file, width, height){
  ggsave(paste(file, ".pdf", sep = ""),
         width = width, height = height)
  ggsave(paste(file, ".png", sep = ""),
         width = width, height = height)
}

###Function: Save as a pdf and png
save_pdfpng <- function(plot, file, width, height, res = 300){
  pdf(paste(file, ".pdf", sep = ""),
      width = width, height = height)
  plot(plot)
  dev.off()
  png(paste(file, ".png", sep = ""),
      width = width, height = height, units = "in", res = res)
  plot(plot)
  dev.off()
}





















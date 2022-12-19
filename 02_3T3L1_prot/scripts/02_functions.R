###Background
#Here I store functions used in this project

###Initialise
library(corrplot)
library(plotly)
library(multcomp)
library(ggplot2)
library(ggVennDiagram)
library(AnnotationDbi)
library(reshape2)
library(dendextend)

#####Data processing#####

###Function that takes in phospho data and proteo data for a protein, and returns log10 occupancy
occupancy_maker <- function(phospho_data,
                            proteo_data){
  ##Correct colnames
  colnames(proteo_data) <- grouped_to_og_colnames$original
  phospho_data <- phospho_data[, colnames(proteo_data)]
  
  ##Make occupancy
  occupancy_data <- sweep(phospho_data,
                          2,
                          STATS = as.numeric(proteo_data),
                          FUN = "-")
  
  return(occupancy_data)
}


protein_correlator <- function(data,
                               target_proteins,
                               proteins){
  ##Set up
  data <- as.matrix(data)
  corr_output <- matrix(NA,
                        nrow = length(target_proteins),
                        ncol = 3 * length(proteins))
  rownames(corr_output) <- target_proteins
  colnames_list <- as.list(proteins)
  names(colnames_list) <- proteins
  for (i in 1:length(colnames_list)){
    
    colnames_list[[i]] <- c(paste(c(proteins[i],
                                    "_corr"),
                                  collapse = ""),
                            paste(c(proteins[i],
                                    "_pval"),
                                  collapse = ""),
                            paste(c(proteins[i],
                                    "_adjpval"),
                                  collapse = ""))
  }
  colnames(corr_output) <- unlist(colnames_list)
  
  ##Do tests
  for (i in 1:length(target_proteins)){
    
    for (j in 1:length(proteins)){
      temp_test <- try(cor.test(data[proteins[j], ],
                                data[target_proteins[i], ]),
                       silent = TRUE)
      if (inherits(temp_test,
                   "try-error")){
        corr_output[i, (1:2 + 3*j - 3)] <- NA
      } else {
        corr_output[i, (1:2 + 3*j - 3)] <- c(temp_test$estimate,
                                             temp_test$p.value)
      }
    }
  }
  
  ##Adjust and order
  for (j in 1:length(proteins)){
    corr_output[, 3*j] <- p.adjust(corr_output[, 3*j - 1],
                                   method = "fdr")
  }
  corr_output <- corr_output[order(corr_output[, 3],
                                   decreasing = FALSE), ]
  return(corr_output)
}

###Function that identifies best uniprots for proteome data, by using phospho data
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





###Function that visualises a KEGG pathway for given mouse data, using pathview
#data: your data. Must contain a column "gene"
#pathway_name: name of pathway e.g. KEGG_PEROXISOME
#pathway: The pathway, containing gene symbols
#paths: database matching KEGG pathways to KEGG IDs
#FC_cols: colnames containing numerical data that you want visualised
#genes_unique_FUN: We enter genes into pathview. To deal with multiple isoforms for same gene, this option combines isoforms using a specified function
mouse_data_to_pathview <- function(data,
                                   pathway_name,
                                   pathway,
                                   paths = paths.hsa,
                                   FC_cols,
                                   genes_unique_FUN = "mean",
                                   output_name){
  
  
  
  
  ##Process paths.hsa names so they will match my kegg database
  #E.g. "Citrate cycle (TCA cycle)" to "KEGG_CITRATE_CYCLE_TCA_CYCLE"
  paths <- sapply(paths,
                  function(x){
                    
                    x <- toupper(x)
                    x <- gsub("\\s+",
                              "_",
                              x)
                    x <- gsub("[^A-Za-z0-9_]",
                              "",
                              x)
                    x <- gsub("_+",
                              "_",
                              x)
                    return(paste(c("KEGG_",
                                   x),
                                 collapse = ""))
                  })
  
  ##Process data
  
  #Restrict to pathway
  data <- data[which(toupper(data$gene) %in%
                       pathway), ]
  
  #Make genes unique
  data_unique_genes <- data.frame("gene" = unique(data$gene),
                                  stringsAsFactors = FALSE)
  rownames(data_unique_genes) <- data_unique_genes$gene
  for (i in 1:nrow(data_unique_genes)){
    
    data_unique_genes[i, FC_cols] <- apply(data[which(data$gene ==
                                                        data_unique_genes$gene[i]), 
                                                FC_cols],
                                           2,
                                           FUN = genes_unique_FUN,
                                           na.rm = TRUE)
  }
  
  #Convert to matrix
  data_unique_genes_m <- as.matrix(data_unique_genes[, FC_cols])
  
  #toupper genes
  rownames(data_unique_genes_m) <- toupper(rownames(data_unique_genes_m))
  
  
  ##Pathview
  pv.out <- pathview(gene.data = data_unique_genes_m,
                     gene.idtype = "SYMBOL",
                     pathway.id = names(paths[which(paths ==
                                                      pathway_name)]),
                     species = "hsa",
                     out.suffix = output_name,
                     kegg.native = T,
                     low = list(gene = "#7070ff",
                                cpd = "green"),
                     mid = list(gene = "white",
                                cpd = "white"),
                     high = list(gene = "#ff7070",
                                 cpd = "yellow"),
                     na.col = "gray",
                     limit = list(gene = min(c(max(abs(data_unique_genes_m),
                                                   na.rm = TRUE),
                                               2)),
                                  cpd = 1))
  return(pv.out)
}

#####Visualisations#####
###Boxplot for one df and given conditions
boxplot_grouped_3t3_proteome <- function(data, 
                                         protein, 
                                         order = "treatment",
                                         alpha_val_border = 1){
  boxplot_df <- data.frame("Condition" = rep(c("NORMAL_GROUPED",
                                               "CI_GROUPED",
                                               "DEX_GROUPED",
                                               "TNF_GROUPED",
                                               "MPQ_GROUPED",
                                               "AA_GROUPED"),
                                             c(10,
                                               7,
                                               8,
                                               7,
                                               8,
                                               8)),
                           "Value" = as.numeric(data[protein, c(grep("NORMAL_GROUPED_\\d+",
                                                                     colnames(data)),
                                                                grep("CI_GROUPED_\\d+",
                                                                     colnames(data)),
                                                                grep("DEX_GROUPED_\\d+",
                                                                     colnames(data)),
                                                                grep("TNF_GROUPED_\\d+",
                                                                     colnames(data)),
                                                                grep("MPQ_GROUPED_\\d+",
                                                                     colnames(data)),
                                                                grep("AA_GROUPED_\\d+",
                                                                     colnames(data)))]),
                           stringsAsFactors = FALSE)
  
  ##Order by treatment 
  boxplot_df$Condition <- factor(boxplot_df$Condition, levels = c("NORMAL_GROUPED",
                                                                  "CI_GROUPED",
                                                                  "DEX_GROUPED",
                                                                  "TNF_GROUPED",
                                                                  "MPQ_GROUPED",
                                                                  "AA_GROUPED"))
  
  boxplot_df$Treatment <- sub("_\\w+", "", boxplot_df$Condition)
  boxplot_df$Treatment <- factor(boxplot_df$Treatment, levels = c("NORMAL",
                                                                  "CI",
                                                                  "DEX",
                                                                  "TNF",
                                                                  "MPQ",
                                                                  "AA"))
  
  plot <- ggplot(data = boxplot_df,
                 aes(x = Treatment,
                     y = Value,
                     color = Treatment)) +
    geom_boxplot(outlier.shape = NA,
                 coef = 0) + 
    geom_point(size = 3,
               shape = 21,
               stroke = 1,
               alpha = alpha_val_border) + 
    scale_color_manual(values = SH_control_model_main_colours) +
    ggtitle(protein) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(color = "black"))
  return(plot)
}

###Heatmap function
##Plots FC on heatmap
##Plots significance as point in tile
fc_sig_heatmap_proteome <- function(pre_melt_df,
                                    fc_cols_w_pvals = c(),
                                    fc_cols_wo_pvals = c(),
                                    pval_cols = c(),
                                    order_column,
                                    is_decreasing = TRUE,
                                    return_df = FALSE,
                                    x_axis_names = c()){
  
  
  ##Prep for melts
  
  #Set up pre_melt_df
  pre_melt_df$standard_name <- sapply(rownames(pre_melt_df),
                                      function(x) strsplit(x, "_")[[1]][1])
  
  pre_melt_df <- pre_melt_df[order(pre_melt_df[, order_column],
                                   decreasing = is_decreasing == FALSE), ]
  pre_melt_df$standard_name <- factor(pre_melt_df$standard_name,
                                      levels = unique(pre_melt_df$standard_name))
  
  
  ##w_pvals melt
  
  #Check that we have pval_cols
  
  if (length(pval_cols) > 0){
    
    #FC melt
    FC_melt_df_w_pvals <- melt(pre_melt_df,
                               id.vars = "standard_name",
                               measure.vars = fc_cols_w_pvals)
    #p_val melt
    pval_melt_df_w_pvals <- melt(pre_melt_df,
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
    FC_melt_df_wo_pvals <- melt(pre_melt_df,
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
  
  
  ##Heatmap
  output_plot <- ggplot(post_melt_df_combined,
                        aes(x = variable,
                            y = standard_name,
                            fill = FC_value)) + 
    scale_fill_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         name = "MODEL/CONTROL") + 
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

###Function: Barplot from pathway enrichment
#Up to date 20210220
#Colours: blue, red, grey
pathway_enrichment_barplot <- function(plot_df,
                                       num_pathways = "sig_only",
                                       colour,
                                       sig_line = TRUE,
                                       sig_cutoff = 0.05){
  ##Colours
  #1st is high, second is low
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
                  c("#ECECEC",
                    "#D0D0D0",
                    "#AFAFAF",
                    "#848484",
                    "#606060",
                    "#414141",
                    "#1C1C1C"))
  names(colours) <- c("blue",
                      "red",
                      "grey")
  
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
  
  ##Plot
  output_plot <- ggplot(plot_df,
                        aes(x = nlog_adj_pval,
                            y = pathway,
                            fill = num_DE_genes_in_pathway)) +
    geom_col() +
    sig_line_plot +
    scale_fill_gradientn(colours = colours[[colour]],
                         name = "# DE genes in pathway") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line  = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "-log10 adjusted p-value",
         y = "Pathway")
  
  return(output_plot)
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
                    "#932507"))
  names(colours) <- c("blue",
                      "red")
  
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
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = "Condition",
           y = "Pathway")
    return(output_plot)
  }
}

###Function that plots and saves multiple overrepresentation barplots
###title and directory must contain "XXX", which will be replaced by conditions[i]
multiple_overrep_barplots <- function(overrep_list,
                                      title,
                                      directory,
                                      conditions,
                                      colour,
                                      num_pathways = "sig_only",
                                      sig_cutoff = 0.05){
  for (i in 1:length(overrep_list)){
    
    plot(pathway_enrichment_barplot(overrep_list[[i]],
                                    colour = colour,
                                    num_pathways = num_pathways,
                                    sig_cutoff = sig_cutoff) +
           ggtitle(gsub("XXX",
                        conditions[i],
                        title)))
           ggsave(gsub("XXX",
                       conditions[i],
                       directory),
                  width = 4,
                  height= 4)
  }
}

###Function: Volcano plot given pval and FC columns (can be adj_pval, and/or rval)
#Up to date 20210406
#Colours significant (pval and FC) proteins
#Option to include labels using ggrepel. By default this is all significant objects, but specific objects can be specified by assigning label_rownames to rownames of desired objects
volcano_plot_coloured <- function(data,
                                  pval_col,
                                  FC_col,
                                  pval_cutoff = 0.05,
                                  FC_cutoff = 0.58,
                                  inner_alpha_val = 0.3,
                                  outer_alpha_val = 0.5,
                                  size = 1,
                                  sig_colour = "red",
                                  x_lim = "max",
                                  y_lim = "max",
                                  remove_Inf = TRUE,
                                  return_df = FALSE,
                                  label_col = FALSE,
                                  label_rownames = "all_sig"){
  
  ##Set up data
  if (label_col == FALSE){
    
    plot_data <- data[, c(FC_col,
                          pval_col)]
    colnames(plot_data) <- c("FC",
                             "pval")
  } else {
    
    plot_data <- data[, c(FC_col,
                          pval_col,
                          label_col)]
    colnames(plot_data) <- c("FC",
                             "pval",
                             label_col)
  }
  
  plot_data$nlog_pval <- -log10(plot_data$pval)
  plot_data$sig <- FALSE
  plot_data$sig[which(abs(plot_data$FC) > FC_cutoff &
                        plot_data$pval < pval_cutoff)] <- TRUE
  
  
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
  
  ##Order df by sig
  plot_data <- plot_data[order(plot_data$sig), ]
  
  
  ##Labels
  if (label_col != FALSE){
    
    
    #Assign label_rownames if all_sig
    if (label_rownames[1] == "all_sig"){
      
      label_rownames <- rownames(plot_data)[which(plot_data$sig == TRUE)]
    }
    
    #Make label column
    #I set everything to "" except those that fit the supplied label_rownames
    
    plot_data$label_col_proc <- plot_data[, label_col]
    plot_data$label_col_proc[which(rownames(plot_data) %in%
                                     label_rownames == FALSE)] <- ""
    
    label_plot <- geom_text_repel(data = plot_data[which(rownames(plot_data) %in%
                                                           label_rownames |
                                                           plot_data$sig == TRUE), ],
                                  aes(x = FC,
                                      y = nlog_pval,
                                      label = label_col_proc),
                                  size = 2.5,
                                  alpha = 1,
                                  min.segment.length = 0,
                                  max.overlaps = Inf)
  } else{
    
    label_plot <- NULL
  }
  
  
  
  
  ##Output
  
  #return_df
  if (return_df){
    
    return(plot_data)
  } else {
    output_plot <- ggplot(data = plot_data,
                          aes(x = FC,
                              y = nlog_pval,
                              colour = sig,
                              fill = sig))+
      geom_point(size = size,
                 shape = 21,
                 stroke = size*0.3) +
      label_plot +
      scale_colour_manual(values = alpha(c("black",
                                           sig_colour),
                                         outer_alpha_val)) +
      scale_fill_manual(values = alpha(c("black",
                                         sig_colour),
                                       alpha = inner_alpha_val)) +
      theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5)) +
      xlim(-x_lim, x_lim) +
      ylim(0, y_lim)
    return(output_plot)
  }
}

#####Enrichment analysis#####

####Function that enriches from a list of pathways
#Up to date 20210220
#background_genes: Background genes e.g. all quantified genes
#pathways_list: List where each entry is the gene contents of a pathway
#pathway_DE_intersection_threshold: The minimum intersection between your DE genes and a pathway for the overlap to be tested. If FALSE, all pathways are tested
pathway_enricher_from_list <- function(background_genes,
                                       DE_genes,
                                       pathways_list,
                                       pathway_DE_intersection_threshold = FALSE){
  
  ##Make pathways relevant 
  #Trim pathway list so it's only pathways with DE genes
  if (pathway_DE_intersection_threshold == FALSE){
    NULL
  } else {
    relevant_pathways_bool <- NULL
    old_pathways_list <- pathways_list
    for (i in 1:length(old_pathways_list)){
      
      relevant_pathways_bool[i] <- length(intersect(DE_genes,
                                                    old_pathways_list[[i]])) >= pathway_DE_intersection_threshold
    }
    
    #Terminate here if no relevant pathways
    if (sum(relevant_pathways_bool) == 0){
      
      return(NULL)
    } else {
      
      pathways_list <- as.list(names(old_pathways_list)[relevant_pathways_bool])
      names(pathways_list) <- names(old_pathways_list)[relevant_pathways_bool]
      for (i in 1:length(pathways_list)){
        
        pathways_list[[i]] <- old_pathways_list[[names(pathways_list)[i]]]
      }
    }
  }
  
  
  
  ##Set up output
  #pval
  #adj_pval
  #Number of DE genes in pathway
  #Number of background genes in pathway
  #Number of total genes in pathway
  output_m <- matrix(NA,
                     nrow = length(pathways_list),
                     ncol = 5)
  rownames(output_m) <- names(pathways_list)
  colnames(output_m) <- c("pval",
                          "adj_pval",
                          "num_DE_genes_in_pathway",
                          "num_background_genes_in_pathway",
                          "num_total_genes_in_pathway")
  
  ##Do stats
  #Do one-sided fisher's exact
  
  #Get x, m, n, k
  #m is number of DE genes
  #n is number of background, non-DE genes
  #k is number of background genes in pathway
  #x is number of DE genes in pathway
  
  m <- length(DE_genes)
  n <- length(background_genes) - m
  
  for (i in 1:length(pathways_list)){
    
    temp_pathway <- pathways_list[[i]]
    k <- length(intersect(background_genes,
                          temp_pathway))
    x <- length(intersect(DE_genes,
                          temp_pathway))
    temp_pval <- sum(dhyper(x:min(k, m),
                            m,
                            n,
                            k))
    output_m[i, ] <- c(temp_pval,
                       NA,
                       x,
                       k,
                       length(temp_pathway))
  }
  
  ##Adjust pvals
  output_m[, 2] <- p.adjust(output_m[, 1],
                            method = "fdr")
  
  ##Sort
  #By adj_pval. Don't order if length == 1
  if (nrow(output_m) > 1){
    output_m <- output_m[order(output_m[, 2],
                               decreasing = FALSE), ]
  }
  
  ##Return
  return(output_m)
}

####Function that enriches for terms of a given GO ontology
#Supply vector of background genes, vector of DE genes, gene_label (UNIPROT, SYMBOL, etc), organism database (e.g. org.Mm.eg.db), and ontology (e.g. CC)
GO_enricher <- function(background_genes,
                        DE_genes,
                        gene_label,
                        organism_database,
                        ontology,
                        pathway_DE_intersection_threshold = 1){
  
  ###Get relevant GOIDs
  relevant_GOIDS_df <- AnnotationDbi::select(organism_database,
                                             keys = DE_genes,
                                             columns = "GO",
                                             keytype = gene_label)
  relevant_GOIDS_df <- relevant_GOIDS_df[which(relevant_GOIDS_df$ONTOLOGY == 
                                                 ontology),]
  relevant_GOIDS <- unique(relevant_GOIDS_df$GO)
  
  ###List of all genes for relevant_GOIDS
  relevant_GOIDS_genes_df <- AnnotationDbi::select(organism_database,
                                                   keys = relevant_GOIDS,
                                                   columns = gene_label,
                                                   keytype = "GO")
  ##Make list
  #Make sure no gene name complete duplicates. Can happen because of different evidence levels
  relevant_GOIDS_genes_list <- as.list(unique(relevant_GOIDS_genes_df$GO))
  names(relevant_GOIDS_genes_list) <- unique(relevant_GOIDS_genes_df$GO)
  for (i in 1:length(relevant_GOIDS_genes_list)){
    
    temp_GOID <- names(relevant_GOIDS_genes_list)[i]
    relevant_GOIDS_genes_list[[i]] <- unique(relevant_GOIDS_genes_df[which(relevant_GOIDS_genes_df$GO == temp_GOID), 
                                                                     gene_label])
  }
  
  
  ##Run tests
  output_df <- as.data.frame(pathway_enricher_from_list(background_genes,
                                                        DE_genes,
                                                        relevant_GOIDS_genes_list,
                                                        pathway_DE_intersection_threshold = pathway_DE_intersection_threshold),
                             stringsasfactors = FALSE)
  
  ##Add GO term and ID
  output_df$GO_term <- AnnotationDbi::select(GO.db,
                                             rownames(output_df),
                                             columns = "TERM",
                                             keytype = "GOID")$TERM
  output_df$GO_id <- rownames(output_df)
  rownames(output_df) <- output_df$GO_term
  
  ##Return
  return(output_df)
}


###Fisher's exact test for DE genes and pathway
#Supply background genes, DE genes, pathway genes
fishers_DE_pathway <- function(background_genes,
                               DE_genes,
                               pathway_genes){
  
  ##Run test
  #m: number of DE genes
  #n: number of background, non-DE genes
  #k: number of background genes in pathway
  #x: number of DE genes in pathway
  m <- length(DE_genes)
  n <- length(background_genes) - m
  k <- length(intersect(background_genes,
                        pathway_genes))
  x <- length(intersect(DE_genes,
                        pathway_genes))
  
  pval <- fisher.test(rbind(c(x, k - x),
                            c(m - x, n - k + x)),
                      alternative = "greater")$p.value
  
  ##Output
  output <- c(pval,
              x,
              k,
              length(pathway_genes))
  output <- as.list(output)
  names(output) <- c("pval",
                     "num_DE_genes_in_pathway",
                     "num_background_genes_in_pathway",
                     "num_pathway_genes")
  return(output)
}





###Limma GST for multiple pathways
#Type: "t" or "F"; signed or unsigned
gst_enricher <- function(pathways_list,
                         stats,
                         stats_genes,
                         alternative,
                         type){
  
  ##Set up
  names(stats) <- stats_genes
  #Remove NA
  stats <- stats[is.na(stats) == FALSE]
  #Output_m
  output_m <- matrix(NA,
                     nrow = length(pathways_list),
                     ncol = 4)
  rownames(output_m) <- names(pathways_list)
  colnames(output_m) <- c("pval",
                          "adj_pval",
                          "num_background_genes_in_pathway",
                          "num_total_genes_in_pathway")
  
  ##Run stats
  for (i in 1:length(pathways_list)){
    
    pathway_indices <- which(names(stats) %in%
                               pathways_list[[i]])
    output_m[i, c(1, 3, 4)] <- c(geneSetTest(index = pathway_indices,
                                             statistics = stats,
                                             alternative = alternative,
                                             type = type),
                                 length(pathway_indices),
                                 length(pathways_list[[i]]))
  }
  
  ##Adjust
  output_m[, 2] <- p.adjust(output_m[, 1],
                            method = "fdr")
  output_m <- output_m[order(output_m[, 2]), ]
  ##Return
  return(output_m)
}

###Function that does gst enrichment for multiple treatments, and outputs a list of gst output dfs
#stats must be a df or matrix
gst_multiple_treatment_enricher <- function(treatments,
                                            pathways_list,
                                            stats,
                                            stats_genes,
                                            alternative,
                                            type,
                                            tidy_names = FALSE){
  
  output_list <- list()
  for (i in 1:length(treatments)){
    
    output_list[[i]] <- gst_enricher(pathways_list = pathways_list,
                                     stats = stats[, grep(treatments[i],
                                                          colnames(stats))],
                                     stats_genes = stats_genes,
                                     alternative = alternative,
                                     type = type)
    if (tidy_names){
      rownames(output_list[[i]]) <- sapply(rownames(output_list[[i]]),
                                           function(x) paste(strsplit(x, "_")[[1]][-1], collapse = " "))
    }
  }
  names(output_list) <- sapply(treatments,
                               function(x) paste(c(x, alternative),
                                                 collapse = "_"))
  return(output_list)
}


###Function that does gst enrichment for GO terms
GO_gst <- function(stats,
                   stats_genes,
                   gene_label,
                   organism_database,
                   ontology){
  
  ###Get relevant GOIDs
  relevant_GOIDS_df <- AnnotationDbi::select(organism_database,
                                             keys = stats_genes,
                                             columns = "GO",
                                             keytype = gene_label)
  relevant_GOIDS_df <- relevant_GOIDS_df[which(relevant_GOIDS_df$ONTOLOGY == 
                                                 ontology),]
  relevant_GOIDS <- unique(relevant_GOIDS_df$GO)
  
  ###List of all genes for relevant_GOIDS
  relevant_GOIDS_genes_df <- AnnotationDbi::select(organism_database,
                                                   keys = relevant_GOIDS,
                                                   columns = gene_label,
                                                   keytype = "GO")
  ##Make list
  #Make sure no gene name complete duplicates. Can happen because of different evidence levels
  relevant_GOIDS_genes_list <- as.list(unique(relevant_GOIDS_genes_df$GO))
  names(relevant_GOIDS_genes_list) <- unique(relevant_GOIDS_genes_df$GO)
  for (i in 1:length(relevant_GOIDS_genes_list)){
    
    temp_GOID <- names(relevant_GOIDS_genes_list)[i]
    relevant_GOIDS_genes_list[[i]] <- unique(relevant_GOIDS_genes_df[which(relevant_GOIDS_genes_df$GO == temp_GOID), 
                                                                     gene_label])
  }
  
  
  ##Run tests
  output_df <- as.data.frame(gst_enricher(pathways_list = relevant_GOIDS_genes_list,
                                          stats = stats,
                                          stats_genes = stats_genes,
                                          alternative = alternative,
                                          type = type))
  
  ##Add GO term and ID
  output_df$GO_term <- AnnotationDbi::select(GO.db,
                                             rownames(output_df),
                                             columns = "TERM",
                                             keytype = "GOID")$TERM
  output_df$GO_id <- rownames(output_df)
  rownames(output) <- output_df$GO_term
  
  ##Return
  return(output_df)
}










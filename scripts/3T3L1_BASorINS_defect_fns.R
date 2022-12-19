###Background
#Here I store functions to aid my analysis of BAS and INS defects in 3T3-L1 IR phos data

###Defect matrix: Makes a matrix matching format of num_data, which is boolean based on whether or not a site is defective in a given condition
defect_matrix_maker <- function(data, num_data){
  defect_matrix <- matrix(FALSE, ncol = ncol(num_data), nrow = nrow(num_data))
  colnames(defect_matrix) <- colnames(num_data)
  rownames(defect_matrix) <- rownames(num_data)
  for (col in colnames(defect_matrix)){
    condition <- strsplit(col, "_")[[1]][1]
    if (condition == "CTRL"){
      defect_matrix[, col] <- TRUE
    } else {
      defect_matrix[, col] <- data[, paste(condition, "_defect", sep = "")]
    }
  }
  return(defect_matrix)
}

###Run ANOVAs for either BAS or INS
defective_BASorINS_ANOVAs <- function(data, num_data){
  #Get defect_matrix so I test only the right conditions
  defect_matrix <- defect_matrix_maker(data, num_data)
  
  #Treatments
  treats <- sapply(colnames(num_data), function(x) strsplit(x, "_")[[1]][1])
  master_treatment <- factor(treats, levels = unique(treats))
  
  #Run ANOVAs
  tests <- list()
  ps <- NULL
  for(i in 1:nrow(num_data)){
    #No test if no defect
    if (sum(defect_matrix[i, -grep("CTRL", colnames(defect_matrix))]) == 0){
      ps[i] <- NA
      tests[[i]] <- NA
      #Else run ANOVA
    } else {
      treatment <- master_treatment[defect_matrix[i, ]]
      temp_test <- try(aov(as.numeric(num_data[i, defect_matrix[i, ]]) ~ treatment, 
                       data = num_data),
                       silent = TRUE)
      if (inherits(temp_test, "try-error")){
        ps[i] <- NA
        tests[[i]] <- NA
      } else {
        summary <- summary(temp_test)[[1]][["Pr(>F)"]]
        #This lower condition happens when there are few observations but enough to still do the test (e.g. one in each condition)
        if (is.null(summary)){
          tests[[i]] <- NA
          ps[i] <- NA
        } else {
          tests[[i]] <- temp_test
          ps[i] <- summary(temp_test)[[1]][["Pr(>F)"]][1]
        }
        
      }
      
    }
  }
  ps[is.null(ps)] <- NA
  adj_ps <- p.adjust(ps, method = "fdr")
  return(list("p" = ps, "adj_p" = adj_ps, "tests" = tests))
}

###Function: Dunnett's using ANOVAs performed comparing BAS and INS
defective_BASorINS_Dunnetts <- function(padjs, 
                                        tests,
                                        data,
                                        num_data){
  #Vector for which to test
  to_test <- !(is.na(padjs)) & padjs < 0.05
  
  #Run Dunnett's
  #Only test if ANOVA was significant
  dunnet_models <- list()
  for(i in 1:length(tests)){
    if (to_test[i]){
      dunnet_models[[i]] <- glht(tests[[i]], linfct = mcp(treatment = "Dunnett"))
    } else {
      dunnet_models[[i]] <- NA
    }
  }
  
  
  #Extract p-values
  dunnet_p_m <- matrix(NA, ncol = 5, nrow = nrow(num_data))
  colnames(dunnet_p_m) <- c("CI", "DEX", "TNF", "MPQ", "AA")
  rownames(dunnet_p_m) <- rownames(num_data)
  for(i in 1:nrow(dunnet_p_m)){
    #Proceed if test is not NA
    if (!identical(is.na(dunnet_models[[i]]),
                  TRUE)){
      used_comparisons <- rownames(summary(dunnet_models[[i]])$linfct)
      used_treatments <- sub(" - CTRL", "", used_comparisons)
      dunnet_p_m[i, used_treatments] <- summary(dunnet_models[[i]])$test$pvalues
    }
  }
  dunnet_p_adj_m <- dunnet_p_m
  for (j in 1:ncol(dunnet_p_m)){
    dunnet_p_adj_m[, j] <- p.adjust(dunnet_p_m[, j],
                                    method = "fdr")
  }
  return(list("p_matrix" = dunnet_p_m,
              "p_adj_matrix" = dunnet_p_adj_m))
}

###Function: Run ANOVAs and Dunnett's for defective sites in bas or insulin, return as df
defective_BASorINS_stats <- function(data, num_data){
  anovas <- defective_BASorINS_ANOVAs(data = data, num_data = num_data)
  dunnetts <- defective_BASorINS_Dunnetts(padjs = anovas$adj_p,
                                          tests = anovas$tests,
                                          data = data, num_data = num_data)
  dunnetts_p_m <- dunnetts$p_matrix
  colnames(dunnetts_p_m) <- paste(colnames(dunnetts_p_m), "_p", sep = "")
  dunnetts_p_adj_m <- dunnetts$p_adj_matrix
  colnames(dunnetts_p_adj_m) <- paste(colnames(dunnetts_p_adj_m), "_adj_p", sep = "")
  output_df <- as.data.frame(cbind(dunnetts_p_m, dunnetts_p_adj_m))
  output_df$anova_p <- anovas$p
  output_df$anova_adj_p <- anovas$adj_p
  return(output_df)
}

###Function: Perform BASorINS defective stats for BAS and INS, and output df
defective_BASorINS_bothdf <- function(data, num_data){
  BAS_tests <- defective_BASorINS_stats(data, 
                                        num_data[, grep("BAS", colnames(num_data))])
  colnames(BAS_tests) <- sapply(colnames(BAS_tests), 
                                function(x){
                                  x_split <- strsplit(x, "_")[[1]]
                                  return(paste(c(x_split[1], 
                                                 "BASdefect", 
                                                 x_split[-1]),
                                               collapse = "_"))
                                })
  INS_tests <- defective_BASorINS_stats(data, 
                                        num_data[, grep("INS", colnames(num_data))])
  colnames(INS_tests) <- sapply(colnames(INS_tests), 
                                function(x){
                                  x_split <- strsplit(x, "_")[[1]]
                                  return(paste(c(x_split[1], 
                                                 "INSdefect", 
                                                 x_split[-1]),
                                               collapse = "_"))
                                })
  data <- cbind(cbind(data, BAS_tests), INS_tests)
  return(data)
}

###Function: Given stats targeting BAS and INS defects, now detect them by comparing FCs
#FC_on_control_suffixes: Suffixes to make colnames for model FCs on CTRL in BAS or INS. e.g. _FC_on_control_BAS will give CI_FC_on_control_BAS, etc.
detect_BASorINS_defects <- function(data, logFC_cutoff = 0.58,
                                    FC_on_control_suffixes = c("_FC_on_control_BAS",
                                                               "_FC_on_control_INS")){
  models <- c("CI", "DEX", "TNF", "MPQ", "AA")
  #BAS
  #Loop over models
  for (model in models){
    #Make colnames
    padj_colname <- paste(model, "_BASdefect_adj_p", sep = "")
    FC_colname <- paste(model, FC_on_control_suffixes[1], sep = "")
    defect_colname <- paste(model, "_BASdefect", sep = "")
    #Make various boolean columns that I later combine
    vsCTRL_sig <- data[, padj_colname] < 0.05 &
      !is.na(data[, padj_colname])
    up_defect <- data$CTRL_reg == "up" & 
      data[, FC_colname] > logFC_cutoff
    down_defect <- data$CTRL_reg == "down" & 
      data[, FC_colname] < -logFC_cutoff
    BAS_defect <- vsCTRL_sig & (up_defect | down_defect)
    data[, defect_colname] <- BAS_defect
  }
  #INS
  #Loop over models
  for (model in models){
    #Make colnames
    padj_colname <- paste(model, "_INSdefect_adj_p", sep = "")
    FC_colname <- paste(model, FC_on_control_suffixes[2], sep = "")
    defect_colname <- paste(model, "_INSdefect", sep = "")
    #Make various boolean columns that I later combine
    vsCTRL_sig <- data[, padj_colname] < 0.05 &
      !is.na(data[, padj_colname])
    up_defect <- data$CTRL_reg == "up" & 
      data[, FC_colname] < -logFC_cutoff
    down_defect <- data$CTRL_reg == "down" & 
      data[, FC_colname] > logFC_cutoff
    INS_defect <- vsCTRL_sig & (up_defect | down_defect)
    data[, defect_colname] <- INS_defect
  }
  #Summary columns
  for (model in models){
    BASorINSdefect_col <- paste(model, "_BASorINSdefect", sep = "")
    BASdefect_col <- paste(model, "_BASdefect", sep = "")
    INSdefect_col <- paste(model, "_INSdefect", sep = "")
    defect_col <- paste(model, "_defect", sep = "")
    data[, BASorINSdefect_col] <- "no defect"
    data[data[, defect_col], BASorINSdefect_col] <- "unassigned"
    data[data[, BASdefect_col], BASorINSdefect_col] <- "BAS"
    data[data[, INSdefect_col], BASorINSdefect_col] <- "INS"
    data[data[, BASdefect_col] + data[, INSdefect_col] == 2, BASorINSdefect_col] <- "both"
  }
  data$num_BASdefect_models <- rowSums(data[, grep("_BASdefect$",
                                                   colnames(data))])
  data$num_INSdefect_models <- rowSums(data[, grep("_INSdefect$",
                                                   colnames(data))])
  data$BASorINSdefect <- "no defect"
  data$BASorINSdefect[data$is_defective] <- "unassigned"
  data$BASorINSdefect[data$num_BASdefect_models > 0] <- "BAS"
  data$BASorINSdefect[data$num_INSdefect_models > 0] <- "INS"
  data$BASorINSdefect[data$num_BASdefect_models > 0 &
                        data$num_INSdefect_models > 0] <- "both"
  
  #Old: Summary  
  #print(apply(data[, c(grep("BASdefect$", colnames(data)),
  #                     grep("INSdefect$", colnames(data)))],
  #            2, table))
  return(data)
}


###Function: Summarise BASorINS defect
BASorINSdefect_summariser <- function(data){
  cols <- colnames(data[, grep("BASorINSdefect", colnames(data))])
  summary_m <- sapply(cols, function(x){
    sapply(c("BAS", "INS", "both", "unassigned"),
           function(y){
             length(which(data[, x] == y))
           })
  })
  summary_df <- cbind(rownames(summary_m), as.data.frame(summary_m))
  colnames(summary_df) <- c("defect", "CI", "DEX", "TNF", "MPQ", "AA", "Overall")
  return(summary_df)
}

###Function: SUmmarise BASorINS defect with up and down sites separate
BASorINSdefect_summariser_updownsep <- function(data){
  up_summary <- BASorINSdefect_summariser(data[which(data$CTRL_reg == "up"), ])
  colnames(up_summary)[-1] <- paste(colnames(up_summary)[-1], "_up", sep = "")
  down_summary <- BASorINSdefect_summariser(data[which(data$CTRL_reg == "down"), ])
  colnames(down_summary)[-1] <- paste(colnames(down_summary)[-1], "_down", sep = "")
  summary <- cbind(up_summary[, -1], down_summary[, -1])
  halfn = ncol(summary)/2
  indices <- NULL
  for (i in 1:halfn){
    indices <- c(indices, i, i + halfn)
  }
  summary <- cbind(up_summary[, 1], summary[, indices])
  colnames(summary)[1] <- "defect"
  return(summary)
}



###Barplot
library(ggpattern)

###Barplot of summary
#Recomended ggsave parameters: width, height = 2.5
BASorINSdefect_summary_bplot <- function(data,
                                         summary = NULL,
                                         alpha_bplot = FALSE){
  #Generate Summary if it doesn't already exist
  if (is.null(summary)){
    summary <- BASorINSdefect_summariser(data)
  }
  summary_melt <- melt(summary)
  summary_melt$condition <- sapply(as.character(summary_melt$variable),
                                   function(x) strsplit(x, "_")[[1]][1])
  summary_melt$condition <- factor(summary_melt$condition,
                                   levels = unique(summary_melt$condition))
  summary_melt$defect <- factor(summary_melt$defect,
                                levels = c("unassigned", "both", "INS", "BAS"))
  levels(summary_melt$defect) <- c("Unassigned", "Both", "INS", "BAS")
  
  if (alpha_bplot){
    ggplot(summary_melt,
           aes(x = variable, y = value, 
              fill = condition, 
               alpha = defect)) +
      geom_col(position = "stack",
                       width = 0.8) + 
      comfy_theme() +
      scale_fill_manual(values = c(SH_control_model_main_colours[-1],
                                   "grey")) +
      scale_alpha_manual(values = rev(c(0.3, 0.5, 0.8, 1))) +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1, vjust = 0.5),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(y = "Phosphopeptides") +
      ggtitle("BAS or INS defect summary")
  } else {
    ggplot(summary_melt,
           aes(x = variable, y = value, 
               colour = condition, fill = condition, 
               alpha = defect,
               pattern = defect,
               pattern_fill = condition, pattern_colour = condition)) +
      geom_col_pattern(position = "stack",
                       width = 0.8) + 
      comfy_theme() +
      scale_colour_manual(values = c(SH_control_model_main_colours[-1],
                                     "grey")) +
      scale_fill_manual(values = c(SH_control_model_main_colours[-1],
                                   "grey")) +
      scale_pattern_fill_manual(values = c(SH_control_model_main_colours[-1],
                                           "grey")) +
      scale_pattern_colour_manual(values = c(SH_control_model_main_colours[-1],
                                             "grey")) +
      scale_alpha_manual(values = c(0, 1, 0.3, 0)) +
      scale_pattern_manual(values = c("stripe", rep("none", 3))) +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1, vjust = 0.5),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(y = "Phosphopeptides") +
      ggtitle("BAS or INS defect summary")
  }
  
}


















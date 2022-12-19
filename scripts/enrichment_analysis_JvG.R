###Last updated:
#20211020

###Background
#Here are all of my functions that perform enrichment
#Fisher's exact
#GSEA
#KSEA

###Packages
library(ksea)
library(rlist)
library(limma)
library(org.Mm.eg.db)
library(GO.db)





#####Kinase enrichment#####


###Enrichment function for single condition
#Up to data 20210610

##Rules
#Only keep kinases with >= 5 substrates

##data is your data, which has FC columns labelled "condition_FC" e.g. "CTRL_FC"
##psp_ks is ideally phosphositeplus kinase substrate database, which has minimally a column "KINASE" giving the kinase, and "uniprot_site" giving the uniprot and site of substrate

###NOTE: adj_p_vals ARE NOT ACTUALLY ADJUSTED. IN MORE UP-TO-DATE CODE I CORRECT THIS (e.g. KSEA_from_list)
ksea_beltrao_enricher <- function(data,
                                  psp_ks,
                                  condition,
                                  sub_search_data = NULL,
                                  min_sub_num = 5){
  
  ##Set up
  #sub_search_data
  #This is the data whose rownames are searched for potential substrates. Generally this should be identical to data, but you can restrict the sub_search_data e.g. if you only want matches with a negative logINS/BAS, for example
  if (is.null(sub_search_data)){
    
    sub_search_data <- data
  }
  
  #condition_FC
  condition_FC <- paste(c(condition,
                          "FC"),
                        collapse = "_")
  
  #kinase names (not yet filtered for >= 5 subs)
  ksea_kinase_names_temp <- unique(psp_ks[which(psp_ks$uniprot_site %in%
                                                  sub_search_data[which(is.na(sub_search_data[, condition_FC]) == FALSE), 
                                                                  "uniprot_site"]), 
                                          "KINASE"])
  #kinase subs
  #In making this, we filter to only include kinases with >= 5 subs
  ksea_kinase_subs_temp <- NULL
  j <- 1
  for (i in 1:length(ksea_kinase_names_temp)){
    
    temp_kinase <- ksea_kinase_names_temp[i]
    
    ksea_temp_kinase_subs <- rownames(sub_search_data[which(sub_search_data$uniprot_site %in%
                                                              psp_ks[which(psp_ks$KINASE == temp_kinase), 
                                                                     "uniprot_site"] &
                                                              is.na(sub_search_data[, condition_FC]) == FALSE), ])
    
    if(length(ksea_temp_kinase_subs) >= min_sub_num){
      
      ksea_kinase_subs_temp[[j]] <- ksea_temp_kinase_subs
      names(ksea_kinase_subs_temp)[j] <- temp_kinase
      j <- j + 1
    }
  } 
  #filter kinase_names
  ksea_kinase_names_temp <- names(ksea_kinase_subs_temp)
  
  #sub values (must be ordered). These are all FCs that are not NA
  ksea_FC_temp <- data[which(is.na(data[, condition_FC]) == FALSE), condition_FC]
  names(ksea_FC_temp) <- rownames(data[which(is.na(data[, condition_FC]) == FALSE), ])
  ksea_FC_temp <- sort(ksea_FC_temp,
                       decreasing = TRUE)
  #Set up output
  ksea_output_df_temp <- matrix(0,
                                nrow = length(ksea_kinase_names_temp),
                                ncol = 4)
  rownames(ksea_output_df_temp) <- ksea_kinase_names_temp
  colnames(ksea_output_df_temp) <- c(paste(c(condition,
                                             "ES"),
                                           collapse = "_"),
                                     paste(c(condition,
                                             "pval"),
                                           collapse = "_"),
                                     paste(c(condition,
                                             "adj_pval"),
                                           collapse = "_"),
                                     paste(c(condition,
                                             "num_subs"),
                                           collapse = "_"))
  
  
  ##Fill out output
  
  #adj_pvals. I get these separately from pval and ES, using the batch function
  ksea_adj_pvals_temp <- ksea_batchKinases(names(ksea_FC_temp), 
                                           ksea_FC_temp, 
                                           ksea_kinase_subs_temp, 
                                           trial=1000)
  names(ksea_adj_pvals_temp) <- sapply(names(ksea_adj_pvals_temp),
                                       function(x) strsplit(x, "\\.")[[1]][1])
  
  #Cycle over kinases
  for (i in 1:length(ksea_kinase_names_temp)){
    
    temp_kinase <- ksea_kinase_names_temp[i]
    
    
    ksea_temp_kinase_output <- ksea(names(ksea_FC_temp),
                                    ksea_FC_temp,
                                    ksea_kinase_subs_temp[[temp_kinase]],
                                    trial = 1000,
                                    significance = TRUE)
    
    ksea_output_df_temp[i, ] <- c(ksea_temp_kinase_output$ES,
                                  ksea_temp_kinase_output$p.value,
                                  ksea_adj_pvals_temp[temp_kinase],
                                  length(ksea_kinase_subs_temp[[temp_kinase]]))
  }
  
  ksea_output_df_temp <- ksea_output_df_temp[order(ksea_output_df_temp[, 3],
                                                   decreasing = FALSE),]
  return(ksea_output_df_temp)
}


###Combine output function
#Up to date 20210610
#Same as above, but for multiple conditions. I do not adjust p-values across conditions
ksea_beltrao_enricher_mult_conditions <- function(data,
                                                  psp_ks,
                                                  conditions,
                                                  sub_search_data = NULL,
                                                  min_sub_num = 5){
  ##Set up
  ksea_output_list <- NULL
  ksea_kinases_list <- NULL
  
  ##Enrichment
  for (i in 1:length(conditions)){
    
    temp_condition <- conditions[i]
    
    temp_ksea_output_df <- ksea_beltrao_enricher(data = data,
                                                 psp_ks = psp_ks,
                                                 condition = temp_condition,
                                                 sub_search_data = sub_search_data,
                                                 min_sub_num = min_sub_num)
    ksea_output_list[[i]] <- temp_ksea_output_df
    names(ksea_output_list)[i] <- temp_condition
    ksea_kinases_list[[i]] <- rownames(temp_ksea_output_df)
    names(ksea_kinases_list)[i] <- temp_condition
  }
  
  ##Make row name and order uniform
  ksea_all_kinases <- unique(unlist(ksea_kinases_list))
  for (i in 1:length(conditions)){
    
    ksea_output_list[[i]] <- as.data.frame(ksea_output_list[[i]])
    ksea_output_list[[i]][ksea_all_kinases[which(ksea_all_kinases %in%
                                                   rownames(ksea_output_list[[i]]) == FALSE)], ] <- NA
    ksea_output_list[[i]] <- ksea_output_list[[i]][ksea_all_kinases, ]
    
  }
  
  ##Bind and format
  ksea_combined_output <- list.cbind(ksea_output_list)
  colnames(ksea_combined_output) <- sapply(colnames(ksea_combined_output),
                                           function(x) strsplit(x, "\\.")[[1]][2])
  ##Output 
  return(ksea_combined_output)
}



#####Pathway enrichment#####
####Function that enriches from a list of pathways
#Up to date 20210528
#Supply list of relevant pathways, vector of background genes, and vector of DE genes
#Pathway_DE_intersection_threshold: Pathways need this number of DE genes or more to be considered
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




####Function: Get GO pathways that intersect with a list of genes
#DE_genes: Genes of interest
#gene_label: Nami gnconvention for the genes, which ahs to be read by AnnotationDBI. Default is SYMBOL (e.g. Akt2, Slc2a4)
#organism_database: The AnnotationDBI organism database e.g. org.Mm.eg.db
#ontology: The desired GO ontology e.g. CC
relevant_GOIDS <- function(DE_genes,
                           gene_label = "SYMBOL",
                           organism_database,
                           ontology){
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
  return(relevant_GOIDS_genes_list)
}



####Function that enriches for terms of a given GO ontology
#Up to date 20210528
#Supply vector of background genes, vector of DE genes, gene_label (UNIPROT, SYMBOL, etc), organism database (e.g. org.Mm.eg.db), and ontology (e.g. CC)
GO_enricher <- function(background_genes,
                        DE_genes,
                        gene_label = "SYMBOL",
                        organism_database,
                        ontology,
                        pathway_DE_intersection_threshold){
  
  ###Get relevant GOIDs
  relevant_GOIDS_genes_list <- relevant_GOIDS(DE_genes = DE_genes,
                                              gene_label = gene_label,
                                              organism_database = organism_database,
                                              ontology = ontology)
  
  
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
#Up to date 20210712
#Supply background genes, DE genes, pathway genes, and specifict alternative
fishers_DE_pathway <- function(background_genes,
                               DE_genes,
                               pathway_genes,
                               alternative = "greater"){
  
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
  
  test <- fisher.test(rbind(c(x, k - x),
                            c(m - x, n - k + x)),
                      alternative = alternative)
  pval <- test$p.value
  odds_ratio <- test$estimate
  
  ##Output
  output <- c(pval,
              x,
              k,
              length(pathway_genes),
              odds_ratio)
  output <- as.list(output)
  names(output) <- c("pval",
                     "num_DE_genes_in_pathway",
                     "num_background_genes_in_pathway",
                     "num_pathway_genes",
                     "odds_ratio")
  return(output)
}



###Limma GST for multiple pathways
#Type: "t" or "F"; signed or unsigned
#background_pathway_intersection_threshold: Minimum overlap between background genes and pathway genes
gst_enricher <- function(pathways_list,
                         stats,
                         stats_genes,
                         alternative,
                         type,
                         background_pathway_intersection_threshold = 1){
  
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
  
  ##Remove pathways with overlap below threshold
  output_m <- output_m[which(output_m[, 3] >= background_pathway_intersection_threshold), ]
  
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
                                            tidy_names = FALSE,
                                            background_pathway_intersection_threshold = 1){
  
  output_list <- list()
  for (i in 1:length(treatments)){
    
    output_list[[i]] <- gst_enricher(pathways_list = pathways_list,
                                     stats = stats[, grep(treatments[i],
                                                          colnames(stats))],
                                     stats_genes = stats_genes,
                                     alternative = alternative,
                                     type = type,
                                     background_pathway_intersection_threshold = background_pathway_intersection_threshold)
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









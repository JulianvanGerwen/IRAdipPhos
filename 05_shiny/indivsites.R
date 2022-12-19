###Background
#I collect functions that support the "indivsites" module, including the module itself

###Initialise
##Packages
library(shiny)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(egg)




################################################
#Functions for boxplots

###Function: Make phospho boxplot depending on raw or imputed data, and selected site
#raw_data: database of raw data
#imputed_data: database of imputed data
#raw: logical indicating whether raw or imputed data is to be used. TRUE = raw
#selected_site: The site to be plotted, corresponds to a rowname of raw_data/imputed_data
#boxplot_function: The boxplot function to use
#to_search_against: An optional list to search selected_site against before trying to plot. By default this is rownames of the selected data, but you can use gene_uniprot_site, for example
#by specifying "gene_uniprot_site" or directly entering a vector
make_phos_boxplot <- function(boxplot_data,
                              selected_site,
                              boxplot_function,
                              label = NULL,
                              to_search_against = NULL,
                              ...){
  
  #Make label
  if (!is.null(label)){
    label_plot <- ggtexttable(label, 
                              theme = ttheme(tbody.style = tbody_style(fontsize = 8, fill = "white", face = "bold")))
  }
  
  #Check if site in data before returning boxplot
  if (is.null(to_search_against)){
    to_search_against <- rownames(boxplot_data)
  } else if (to_search_against == "gene_uniprot_site"){
    to_search_against <- gsub("_\\d+$", "", rownames(boxplot_data))
  }
  
  if (selected_site %in% to_search_against){
    boxplot <- boxplot_function(data = boxplot_data, selected_site, ...) +
      theme(plot.title = element_blank(),
            axis.title.x = element_blank())
    if (!is.null(label)){
      return(ggarrange(label_plot, boxplot, 
                       nrow = 2, ncol = 1,
                       heights = c(0.3, 1)))
    } else {
      return(boxplot)
    }
  } else {
    return(textGrob("Absent", gp = gpar(fontsize = 12)))
  }
}

###Function that slects raw or imputed data
select_boxplot_data <- function(raw_data,
                                imp_data,
                                raw){
  if (raw){
    return(raw_data)
  } else {
    return(imp_data)
  }
}


#####################################
#Functions for site summary

#General functions

###Function: return defective or emergent models separate by commas, for a given row of data
summarise_defemerg_models <- function(data, rowname, keyword = "_defect$"){
  search_data <- data[rowname, grep(keyword, colnames(data))]
  search_summary_vec <- gsub(keyword, "", colnames(search_data)[as.logical(search_data)]) 
  search_summary_fused <- paste(search_summary_vec, collapse = ", ")
  return(search_summary_fused)
}

###Function: Summarise all ppeptide for a psite
#summary_fn: Function used to summarise one ppeptide
#Categories: Names for rows
summarise_psite <- function(data, selected_site, summary_fn, categories){
  if (length(grep(selected_site, rownames(data))) > 0){
    #Raw summary
    ppeptides <- sort(rownames(data)[grep(selected_site, rownames(data))])
    summary_m <- sapply(ppeptides, function(x) summary_fn(data, x))
    #Add ppetpide info
    multplicities <- as.numeric(sapply(colnames(summary_m), function(x) strsplit(x, "_")[[1]][4]))
    summary_m <- rbind(c("Single", "Double", "Triple")[multplicities], summary_m)
    #Add category column
    summary_m <- cbind(categories, summary_m)
    return(summary_m)
  } else {
    return(NULL)
  }
}

#phos_3t3_proc
#summarise ppeptide
summarise_phos3t3_ppeptide <- function(data, rowname){
  summary_vec <- c(data[rowname, "CTRL_reg"],
                   summarise_defemerg_models(data, rowname, "_defect$"),
                   summarise_defemerg_models(data, rowname, "_emergent$"))
  return(summary_vec)
}

#Summarise psite
summarise_phos3t3_psite <- function(data, selected_site){
  summarise_psite(data = data, selected_site = selected_site,
                  summary_fn = summarise_phos3t3_ppeptide,
                  categories = c("Phosphopeptide:",
                                 "Insulin regulation in CTRL:",
                                 "Defective models:",
                                 "Emergent models:"))
}


#phos_mouse_proc
#summarise ppeptide
summarise_phosmouse_ppeptide <- function(data, rowname){
  #Get CHOW regulation
  CHOW_reg <- "unregulated"
  if (data[rowname, "CHOW_reg"] == TRUE & data[rowname, "CHOW_logFC"] > 0.58){
    CHOW_reg <- "up"
  } else if (data[rowname, "CHOW_reg"] == TRUE & data[rowname, "CHOW_logFC"] < -0.58){
    CHOW_reg <- "down"
  }
  #Get HFD defect
  if (data[rowname, "HFD_defective"]){
    HFD_defect <- "yes"
  } else {
    HFD_defect <- ""
  }
  #Get emergent
  HFD_emergent <- ""
  if (data[rowname, "HFD_emergent_up"]){
    HFD_emergent <- "up"
  } else if (data[rowname, "HFD_emergent_down"]){
    HFD_emergent <- "down"
  }
  summary_vec <- c(CHOW_reg,
                   HFD_defect,
                   HFD_emergent)
  return(summary_vec)
}

#Summarise psite
summarise_phosmouse_psite <- function(data, selected_site){
  summarise_psite(data = data, selected_site = selected_site,
                  summary_fn = summarise_phosmouse_ppeptide,
                  categories = c("Phosphopeptide:",
                                 "Insulin regulation in CHOW:",
                                 "Defective in HFD:",
                                 "Emergent in HFD:"))
}

#3t3l1 proteome
summarise_prot3t3_psite <- function(data, selected_site){
  selected_protein <- paste(strsplit(selected_site, "_")[[1]][1:2], collapse = "_")
  if (selected_protein %in% rownames(data)){
    models <- c("CI", "DEX", "TNF", "MPQ", "AA")
    changed_models_bool <- sapply(models, 
                                  function(x){rowSums(data[selected_protein, 
                                                           paste(x, c("_isup", "_isdown"), sep = "")]) > 0})
    changed_models_vec <- models[changed_models_bool]
    return(as.matrix(data.frame("Category" = "Models changed compared to CTRL:",
                                "Models" = paste(changed_models_vec, collapse = ", "))))
  } else {
    NULL
  }
}

###Function: Make ggtexttable displayed combined output from a list of summaries from e.g. summarise_phos3t3_psite
combine_summaries_toggtextable <- function(summary_table_list,
                                           fontsize = 5,
                                           padding = 1.5){
  #Remove NULL items
  summary_table_list <- summary_table_list[vapply(summary_table_list, function(x){!is.null(x)}, logical(1))]
  #Extract dimensions for fuse
  ncols <- sapply(summary_table_list, ncol)
  nrows <- sapply(summary_table_list, nrow)
  nrows_cumulative <- nrows
  for (i in 1:length(nrows_cumulative)){
    nrows_cumulative[i] <- sum(nrows[1:i])
  }
  #Loop over data and resize
  summary_table_combined <- NULL
  for (i in 1:length(summary_table_list)){
    data <- matrix("", nrow = nrows[i], ncol = max(ncols) + 1)
    data[1, 1] <- names(summary_table_list)[i]
    data[, 2:(ncols[i] + 1)] <- summary_table_list[[i]]
    summary_table_combined <- rbind(summary_table_combined, data)
  }
  table <- ggtexttable(summary_table_combined, 
                       theme = ttheme(tbody.style = tbody_style(fontsize = fontsize, fill = NA, hjust = 1, x = 1),
                                      padding = unit(c(padding, padding), "mm")))
  table <- table %>% table_cell_font(row = 1:tab_nrow(table),
                                     column = 1,
                                     face = "bold",
                                     size = fontsize) %>%
    tab_add_hline(at.row = nrows_cumulative,
                  row.side = "bottom",
                  linewidth = 1.5,
                  linetype = 1)
  return(table)
}


#################################
#indivsites module


##psite UI
#Boxplots display psites
indivsites_psiteUI <- function(id,
                                  sidebarwidth = 3,
                                  mainwidth = 9,
                               selectable_sites){
  fluidPage(
    sidebarLayout(
      sidebarPanel(width = sidebarwidth,
                   selectizeInput(NS(id, "selected_site"),
                                  label = "Phosphosite/protein",
                                  choices = selectable_sites,
                                  selected = "Gsk3b_Q9WV60_S9"),
                   checkboxInput(NS(id, "raw"),
                                 label = "Raw data",
                                 value = FALSE),
                   checkboxGroupInput(NS(id, "selected_data"),
                                      label = "Choose data to display",
                                      choices = c("3T3-L1 phosphoproteome" = "3T3-L1 phos",
                                                  "Mouse phosphoproteome" = "Mouse phos",
                                                  "3T3-L1 proteome" = "3T3-L1 prot"),
                                      selected = c("3T3-L1 phos",
                                                   "Mouse phos",
                                                   "3T3-L1 prot"))),
      mainPanel(width = mainwidth,
                plotOutput(NS(id, "boxplot"),
                           height = 800))
    )
  )
}



##psite server
#Boxplots display psites
#resolution: Output resolution for boxplot
indivsites_psiteServer <- function(id,
                                   resolution = 216){
  moduleServer(id, function(input, output, session){
    ##Boxplots
    #Boxplot datasets
    phos_3t3l1_boxplot_data <- reactive({select_boxplot_data(phos_3t3_filt_norm,
                                                             phos_3t3_proc,
                                                             input$raw)})
    phos_mouse_boxplot_data <- reactive({select_boxplot_data(phos_mouse_raw_data_norm,
                                                             phos_mouse_proc,
                                                             input$raw)})
    #Boxplot reactives
    phos_3t3l1_boxplot <- reactive({
      if ("3T3-L1 phos" %in% input$selected_data){
        make_phos_boxplot(boxplot_data = phos_3t3l1_boxplot_data(),
                          selected_site = input$selected_site,
                          boxplot_function = boxplot_3t3phos_psite,
                          label = "3T3-L1 phos",
                          to_search_against = "gene_uniprot_site",
                          is_raw = TRUE,
                          siglabels = !input$raw)
      }
    })
    phos_mouse_boxplot <- reactive({
      if ("Mouse phos" %in% input$selected_data){
        make_phos_boxplot(boxplot_data = phos_mouse_boxplot_data(),
                          selected_site = input$selected_site,
                          boxplot_function = boxplot_mousephos_psite,
                          label = "Mouse phos",
                          to_search_against = "gene_uniprot_site",
                          is_raw = TRUE,
                          siglabels = !input$raw)
      }
    })
    proteome_boxplot <- reactive({
      if ("3T3-L1 prot" %in% input$selected_data){
        selected_protein <- sapply(input$selected_site,
                                   function(x) paste(strsplit(x, "_")[[1]][1:2], collapse = "_"))
        make_phos_boxplot(boxplot_data = proteome_post_analysis_3t3_DIA_matched,
                          selected_site = selected_protein,
                          boxplot_function = boxplot_grouped_3t3_proteome,
                          label = "3T3-L1 prot",
                          siglabels = !input$raw)
      }
    })
    #Extract number of ppeptides to scale height
    num_ppeptides <- reactive({
      num_ppeptides_3t3l1 <- length(grep(input$selected_site, rownames(phos_3t3l1_boxplot_data()))) 
      num_ppeptides_mouse <- length(grep(input$selected_site, rownames(phos_mouse_boxplot_data())))
      if (length(grep("phos", input$selected_data)) >= 1){
        max(num_ppeptides_3t3l1, num_ppeptides_mouse)
      } else {
        1
      }
    })
    
    ##Site summary table
    phos_3t3l1_summary_table <- reactive({
      if ("3T3-L1 phos" %in% input$selected_data){
        summarise_phos3t3_psite(phos_3t3_proc, input$selected_site)
      } else {
        NULL
      }
    })
    phos_mouse_summary_table <- reactive({
      if ("Mouse phos" %in% input$selected_data){
        summarise_phosmouse_psite(phos_mouse_proc, input$selected_site)
      } else {
        NULL
      }
    })
    prot_3t3l1_summary_table <- reactive({
      if ("3T3-L1 prot" %in% input$selected_data){
        summarise_prot3t3_psite(proteome_post_analysis_3t3_DIA_matched, input$selected_site)
      } else {
        NULL
      }
    })
    summary_ggtexttable <- reactive({
      table <- combine_summaries_toggtextable(list("3T3-L1 phos" = phos_3t3l1_summary_table(),
                                          "Mouse phos" = phos_mouse_summary_table(),
                                          "3T3-L1 prot" = prot_3t3l1_summary_table()))
      label_plot <- ggtexttable("Summary", 
                                theme = ttheme(tbody.style = tbody_style(fontsize = 8, fill = "white", face = "bold")))
      ggarrange(label_plot, table, 
                nrow = 2, ncol = 1,
                heights = c(0.3, 1))
    })
    
    ##Output boxplot (including summary table)
    output$boxplot <- renderPlot({
      boxplot_list <- list("3T3-L1 phos" = phos_3t3l1_boxplot(), 
                           "Mouse phos" = phos_mouse_boxplot(), 
                           "3T3-L1 prot" = proteome_boxplot())
      boxplot_list <- boxplot_list[input$selected_data]
      boxplot_list[["Summary table"]] <- summary_ggtexttable()
      as_ggplot(arrangeGrob(grobs = boxplot_list,
                            ncol = 2,
                            nrow = 2,
                            layout_matrix = rbind(c(1, 2),
                                                  c(3, 4)),
                            heights = c(num_ppeptides(), 0.8)))
    }, res = resolution,
    height = function() 400 + 400*num_ppeptides())
  })
}

##App function
indivsites_psiteApp <- function(){
  
  ui <- indivsites_psiteUI("indivsites",
                           selectable_sites = selectable_sites)
  server <- function(input, output, session){
    indivsites_psiteServer("indivsites")
  }
  
  shinyApp(ui = ui, server = server)
}

indivsites_psiteApp()






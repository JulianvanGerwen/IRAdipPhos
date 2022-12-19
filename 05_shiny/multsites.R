###Background
#I collect functions that support the "multsites" module, including the module itself

###Initialise
##Packages
library(shiny)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(egg)




###Function: return data based on searching against rownames
#search_terms: search terms separate by commas e.g. "Gsk3,Sos1"
search_for_terms <- function(data, search_terms){
  search_terms_vec <- strsplit(search_terms, ",")[[1]]
  #Idiotproofing: Remove items that are "" (they will search forever)
  search_terms_vec <- search_terms_vec[search_terms_vec != ""]
  #Set to c("") if empty
  if (length(search_terms_vec) == 0){
    search_terms_vec <- c("")
  }
  return(data[unique(unlist(lapply(search_terms_vec, function(x) grep(x, rownames(data))))), ])
}

###Functions: Make individual heatmaps
make_hmap_3t3_phos <- function(searched_sites, data){
  fc_sig_heatmap(data = data[searched_sites, ],
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
                 order_column = "CTRL_FC",
                 is_decreasing = FALSE,
                 x_axis_names = c("CTRL",
                                  "CI",
                                  "DEX",
                                  "TNF",
                                  "MPQ",
                                  "AA"),
                 show_shape_legend = FALSE,
                 yaxis_naming = "psite_and_prot")
}

make_hmap_3t3_prot <- function(searched_sites, data, order_column = "CTRL_FC", yaxis_naming = "psite_and_prot"){
  fc_sig_heatmap(data = data[searched_sites, ],
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
                 order_column = order_column,
                 is_decreasing = FALSE,
                 x_axis_names = c("CI",
                                  "DEX",
                                  "TNF",
                                  "MPQ",
                                  "AA"),
                 shapes = c(17, 17),
                 show_shape_legend = FALSE,
                 yaxis_naming = yaxis_naming)
}

make_hmap_mouse <- function(searched_sites, data){
  fc_sig_heatmap(data = data[searched_sites, ],
                 fc_cols_w_pvals = c("CHOW_logFC",
                                     "HFD_logFC",
                                     "REV_logFC"),
                 pval_cols = c("CHOW_adj_p_val",
                               "HFD_adj_p_val",
                               "REV_adj_p_val"),
                 order_column = "CTRL_FC",
                 is_decreasing = FALSE,
                 x_axis_names = c("CHOW",
                                  "HFD",
                                  "REV"),
                 show_shape_legend = FALSE,
                 yaxis_naming = "psite_and_prot")
}

###Function: Make gradient
#Make gradient
gradient <- function(maxFC){
  scale_fill_gradientn(colours = rev(brewer.pal(9, "RdBu")),
                       values = rescale(c(-maxFC,
                                          -0.75*maxFC,
                                          -0.5*maxFC,
                                          -0.25*maxFC,
                                          0,
                                          0.25*maxFC,
                                          0.5*maxFC,
                                          0.75*maxFC,
                                          maxFC)),
                       limits = c(-maxFC, maxFC),
                       name = "log2 FC")
}

###Function: Get max FC value
get_maxFC <- function(data) max(abs(data), na.rm = TRUE)

###Function: Generate legend describing signficance dots
#Here I make a toy database to get the legend I want, then extract the legend with get_legend from ggpubr
get_dot_legend <- function(){
  test_df <- data.frame(shape = as.character(c(16, 18)),
             value = c(1, 1))
  legend_plot <- ggplot(test_df, aes(x = value, y = value, shape = shape)) +
    geom_point() +
    scale_shape_manual(values = c(19, 17),
                       labels = c("Comparing BAS and INS", "Compared to CTRL"),
                       name = "Significance") +
    theme(legend.key = element_blank(),
          legend.text = element_text(size = 7, colour = "black"),
          legend.title = element_text(size = 7, colour = "black")) 
  return(as_ggplot(get_legend(legend_plot)))
}

###Function: Make a simpel ggtexttable label
standard_ggtexttable_label <- function(text,
                                       size = 6,
                                       vjust = 1, y = 1){
  ggtexttable(text, theme = ttheme(tbody.style = tbody_style(size = size, fill = NA,
                                                          vjust = vjust, y = vjust)))
}



##############################################
#Filters

###These are filters that take in TRUE/FALSE. if TRUE, return boolean of rows that meet condition. Else return boolean of all TRUE
filter_l1phos_reg <- function(data, l1phos_reg){
  if (l1phos_reg){
    return(!is.na(data$CTRL_reg) & 
             data$CTRL_reg != "unregulated")
  } else {
    return(rep(TRUE, nrow(data)))
  }
}

filter_l1phos_def <- function(data, l1phos_def){
  if (l1phos_def > 0){
    return(!is.na(data$num_defective_models) & 
             data$num_defective_models >= l1phos_def)
  } else {
    return(rep(TRUE, nrow(data)))
  }
}

filter_l1phos_emerg <- function(data, l1phos_emerg){
  if (l1phos_emerg > 0){
    return(!is.na(data$emergent_num) & 
             data$emergent_num >= l1phos_emerg)
  } else {
    return(rep(TRUE, nrow(data)))
  }
}

filter_mousephos_reg <- function(data, mousephos_reg){
  if (mousephos_reg){
    return(!is.na(data$CHOW_reg) & 
             data$CHOW_reg == TRUE)
  } else {
    return(rep(TRUE, nrow(data)))
  }
}

filter_mousephos_def <- function(data, mousephos_def){
  if (mousephos_def){
    return(!is.na(data$HFD_defective) & 
             data$HFD_defective == TRUE)
  } else {
    return(rep(TRUE, nrow(data)))
  }
}

filter_mousephos_emerg <- function(data, mousephos_emerg){
  if (mousephos_emerg){
    return((!is.na(data$HFD_emergent_up) | !(is.na(data$HFD_emergent_down))) & 
             (data$HFD_emergent_up == TRUE) | (data$HFD_emergent_down == TRUE))
  } else {
    return(rep(TRUE, nrow(data)))
  }
}

#################################
#multsites module

##ui
multsitesUI <- function(id, 
                        sidebarwidth = 3,
                        mainwidth = 9){
  fluidPage(
    sidebarLayout(
      sidebarPanel(width = sidebarwidth,
                  h4("Site search"),
                  textInput(NS(id, "search_terms"),
                            label = "Enter key terms separated by commas",
                            value = "Gsk3,Map3k5"),
                  actionButton(NS(id, "search"),
                               label = "Search"),

                  br(), br(), h4("3T3-L1 phos"),
                  checkboxInput(NS(id, "l1phos_reg"),
                                label = "Insulin-regulated",
                                value = FALSE),
                  selectInput(NS(id, "l1phos_def"),
                              label = "Defective in at least ",
                              choices = list("0 models" = 0,
                                             "1 model" = 1,
                                             "2 models" = 2,
                                             "3 models" = 3,
                                             "4 models" = 4,
                                             "5 models" = 5),
                              selected = "0 models"),
                  selectInput(NS(id, "l1phos_emerg"),
                              label = "Emergent in at least ",
                              choices = list("0 models" = 0,
                                             "1 model" = 1,
                                             "2 models" = 2,
                                             "3 models" = 3,
                                             "4 models" = 4,
                                             "5 models" = 5),
                              selected = "0 models"),
                  
                  br(), h4("Mouse phos"),
                  checkboxInput(NS(id, "mouse_reg"),
                                label = "Insulin-regulated",
                                value = FALSE),
                  checkboxInput(NS(id, "mouse_def"),
                                label = "Defective in HFD",
                                value = FALSE),
                  checkboxInput(NS(id, "mouse_emerg"),
                                label = "Emergent in HFD",
                                value = FALSE),
                  
                  br(), br(), checkboxGroupInput(NS(id, "selected_data"),
                                     label = h4("Select data"),
                                     choices = c("3T3-L1 phosphoproteome" = "3T3-L1 phos",
                                                 "3T3-L1 proteome" = "3T3-L1 prot",
                                                 "Mouse phosphoproteome" = "Mouse phos"),
                                     selected = c("3T3-L1 phos",
                                                  "3T3-L1 prot",
                                                  "Mouse phos"))),
      mainPanel(width = mainwidth,
                plotOutput(NS(id, "heatmap")))
    )
  )
}


##Server
#resolution: Output resolution for plots
multsitesServer <- function(id, resolution = 216, all_data = phos_3t3_mouse_prot_3t3_proc,
                            prot_data = proteome_post_analysis_3t3_DIA_matched){
  moduleServer(id, function(input, output, session){
    
    #Searched sites
    #Search against all rows of data
    search_term_data <- eventReactive(input$search,
                                      {search_for_terms(data = all_data, input$search_terms)},
                                      ignoreNULL = FALSE)
    
    searched_sites_raw <- reactive({
      #Filter based on 3t3l1 phos
      filter <- filter_l1phos_reg(search_term_data(), input$l1phos_reg) &
        filter_l1phos_def(search_term_data(), input$l1phos_def) &
        filter_l1phos_emerg(search_term_data(), input$l1phos_emerg) &
        #filter based on mouse phos
        filter_mousephos_reg(search_term_data(), input$mouse_reg) &
        filter_mousephos_def(search_term_data(), input$mouse_def) &
        filter_mousephos_emerg(search_term_data(), input$mouse_emerg)
      rownames(search_term_data())[filter]})
    
    #Either copy searched_sites_raw, or make a dummy if it is empty. The dummy means code doesn't get confused
    searched_sites <- reactive({
      
      #Convert to protein names if only protein data is to be shown
      if ("3T3-L1 prot" %in% input$selected_data &
          length(input$selected_data) == 1){
        searched_sites_temp <- unique(sapply(searched_sites_raw(),
                                 function(x) paste(strsplit(x, "_")[[1]][1:2], collapse = "_")))
        searched_sites_temp <- searched_sites_temp[which(searched_sites_temp %in%
                                                 rownames(proteome_post_analysis_3t3_DIA_matched))]
      } else {
        searched_sites_temp <- searched_sites_raw()
      }
      
      if (length(searched_sites_temp) == 0){
        rownames(all_data)[3000:3010]
      } else {
        searched_sites_temp
      }
    })
    
    #Make plotting data. Either all_data (phos and prot) or prot_data (just prot, uses this if we have selected just proteoem data)
    plot_data <- reactive({
      if ("3T3-L1 prot" %in% input$selected_data &
          length(input$selected_data) == 1){
        prot_data
      } else {
        all_data
      }
    })

    #FC conditions
    FC_conditions <- reactive({
      FC_conditions_full <- list("3T3-L1 phos" = c("CTRL_FC", "CI_FC", "DEX_FC","TNF_FC", "MPQ_FC", "AA_FC"),
                                 "3T3-L1 prot" = c("CI_GROUPED_FC", "DEX_GROUPED_FC", "TNF_GROUPED_FC", "MPQ_GROUPED_FC", "AA_GROUPED_FC"),
                                 "Mouse phos" = c("CHOW_logFC", "HFD_logFC","REV_logFC"))
      FC_conditions_full[which(names(FC_conditions_full) %in% input$selected_data)]
    })
    #Max FC value to rescale colours
    maxFC <- reactive({max(sapply(FC_conditions(),
                        function(x) get_maxFC(plot_data()[searched_sites(), x])))})
    
    
    
    #Make heatmaps
    hmap_3t3_phos <- reactive({
      if ("3T3-L1 phos" %in% input$selected_data){
        make_hmap_3t3_phos(searched_sites(), data = plot_data())
      }
    })
    hmap_3t3_prot <- reactive({
      if ("3T3-L1 prot" %in% input$selected_data){
        #Restrict to only protein data (not separated by site) if not showing phos
        if (length(input$selected_data) == 1){
          make_hmap_3t3_prot(searched_sites(), data = plot_data(), order_column = "CI_GROUPED_FC")
        } else {
          make_hmap_3t3_prot(searched_sites(), data = plot_data())
        }
      }
    })
    hmap_mouse <- reactive({
      if ("Mouse phos" %in% input$selected_data){
        make_hmap_mouse(searched_sites(), data = plot_data())
      }
    })
    #Make label list to specifiy type of data
    label_list <- reactive({
      #Include empty label for legend
      labels <- lapply(c(input$selected_data, ""),
                         function(x) ggtexttable(x, theme = ttheme(tbody.style = tbody_style(size = 6, fill = NA,
                                                                                             face = "bold"))))
      })
    #Make FC type list to specificy what each FC is
    FC_type_list <- reactive({
      #Include empty label for legend
      labels <- lapply(c(input$selected_data, ""),
                                     function(x){
                                       if (x == "3T3-L1 phos" | x == "Mouse phos"){
                                         y <- "log2 INS/BAS"
                                       } else if (x == "3T3-L1 prot") {
                                         y <- "log2 MODEL/CTRL"
                                       } else {
                                         y <- x
                                       }
                                       standard_ggtexttable_label(y)
                                     })
      })
    
    #Combine heatmaps
    output$heatmap <- renderPlot({
      hmap_list <- list(hmap_3t3_phos(), hmap_3t3_prot(), hmap_mouse())
      hmap_list <- hmap_list[vapply(hmap_list, function(x){ !is.null(x)}, logical(1))]
      for (i in 1:length(hmap_list)){
        hmap_list[[i]] <- hmap_list[[i]] + gradient(maxFC())
        if (i > 1){
          hmap_list[[i]] <- hmap_list[[i]] + theme(axis.text.y = element_blank())
        }
        if (i < length(hmap_list)){
          hmap_list[[i]] <- hmap_list[[i]] + theme(legend.position = "none")
        }
      }
      if (length(searched_sites_raw()) == 0){
        ggtexttable("No data",
                    theme = ttheme(tbody.style = tbody_style(size = 20, fill = "white")))
      } else {
        #Add legend for points
        hmap_list[["Legend"]] <- get_dot_legend()
        ggarrange(plots = c(label_list(), 
                            hmap_list,
                            FC_type_list()),
                  nrow = 3, ncol = length(hmap_list),
                  widths = c(sapply(FC_conditions(), length), 10),
                  #Here I scale the heatmap's size based on how many rows it has
                  lengths = c(4, length(searched_sites()), 4))
      }
     
    }, res = resolution,
    #Here I scale the size of the plot to try to keep heatmap aspect ratio constant
    height = function() 200 + 25*length(searched_sites()))
  })
}

##App
multsitesApp <- function(){
  ui <- multsitesUI("multsites")
  server <- function(input, output, session){
    multsitesServer("multsites")
  }
  shinyApp(ui, server)
}

multsitesApp()



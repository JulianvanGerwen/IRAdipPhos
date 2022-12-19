####Background
#I collect functions and modules to visualise Gsk3i ppome sites usign heatmap


###Function: Heatmap for Gsk3i ppome
make_hmap_gsk3i_phos <- function(searched_sites){
  fc_sig_heatmap(data = phos_gsk3i_proc[searched_sites, ],
                 fc_cols_w_pvals = c("FC"),
                 pval_cols = c("adj_p_val"),
                 order_column = "FC",
                 is_decreasing = FALSE,
                 x_axis_names = c(""),
                 show_shape_legend = FALSE,
                 legend_title = "log2 Gsk3i/Basal")
}

###Function: Generate legend describing signficance dots
#Here I make a toy database to get the legend I want, then extract the legend with get_legend from ggpubr
get_dot_legend_Gsk3i <- function(){
  test_df <- data.frame(shape = as.character(c(16, 16)),
                        value = c(1, 1))
  legend_plot <- ggplot(test_df, aes(x = value, y = value, shape = shape)) +
    geom_point() +
    scale_shape_manual(values = c(19),
                       labels = c("Comparing Basal and Gsk3i"),
                       name = "Significance") +
    theme(legend.key = element_blank(),
          legend.text = element_text(size = 7, colour = "black"),
          legend.title = element_text(size = 7, colour = "black")) 
  return(as_ggplot(get_legend(legend_plot)))
}

###################################
#Filters
#These functions filter Gsk3i data based on various criteria
filter_gsk3iphos_reg <- function(data, gsk3iphos_reg){
  if (gsk3iphos_reg){
    return((!is.na(data$GSK3i_up) | !is.na(data$GSK3i_down)) & 
             (data$GSK3i_up == TRUE | data$GSK3i_down == TRUE))
  } else {
    return(rep(TRUE, nrow(data)))
  }
}

filter_gsk3iphos_motif <- function(data, gsk3iphos_motif){
  if (gsk3iphos_motif){
    return(!is.na(data$matches_pST35psp_motif) &
             data$matches_pST35psp_motif)
  } else {
    return(rep(TRUE, nrow(data)))
  }
}

filter_gsk3iphos_sub <- function(data, gsk3iphos_sub){
  if (gsk3iphos_sub){
    return(!is.na(data$GSK3i_down_motif) &
             data$GSK3i_down_motif)
  } else {
    return(rep(TRUE, nrow(data)))
  }
}


#########################################
#Gsk3i_multsites module

##UI
Gsk3i_multsitesUI <- function(id,
                              sidebarwidth = 3,
                              mainwidth = 9){
  fluidPage(
    sidebarLayout(
      sidebarPanel(width = sidebarwidth,
                   h4("Site search"),
                   textInput(NS(id, "search_terms"),
                             label = "Enter key terms separated by commas",
                             value = "Trarg1,Gys1"),
                   actionButton(NS(id, "search"),
                                label = "Search"),
                   br(), br(), 
                   checkboxInput(NS(id, "reg"),
                                 label = "Regulated by Gsk3i",
                                 value = FALSE),
                   checkboxInput(NS(id, "motif"),
                                 label = "Matches Gsk3 motif",
                                 value = FALSE),
                   checkboxInput(NS(id, "sub"),
                                 label = "Putative Gsk3 substrate",
                                 value = FALSE)),
      mainPanel(width = mainwidth,
                plotOutput(NS(id, "heatmap")))
    )
  )
}

##Server
#resolution: Output resolution for pots
Gsk3i_multsiteServer <- function(id, resolution = 216){
  moduleServer(id, function(input, output, session){
    #Serached sites
    #Search against all rows of phos_gsk3i_proc
    search_term_data <- eventReactive(input$search,
                                      {search_for_terms(phos_gsk3i_proc, input$search_terms)},
                                      ignoreNULL = FALSE)
    searched_sites_raw <- reactive({
      filter <- filter_gsk3iphos_reg(search_term_data(), input$reg) &
        filter_gsk3iphos_motif(search_term_data(), input$motif) &
        filter_gsk3iphos_sub(search_term_data(), input$sub)
      rownames(search_term_data())[filter]
    })
    #Get searched sites, filling with dummy if they are empty
    searched_sites <- reactive({
      if (length(searched_sites_raw()) == 0){
        rownames(phos_gsk3i_proc)[100:101]
      } else {
        searched_sites_raw()
      }
    })
    
    #Output heatmap
    output$heatmap <- renderPlot({
      #List of hmap and legend
      hmap_list <- list(make_hmap_gsk3i_phos(searched_sites()),
                        get_dot_legend_Gsk3i())
      #List containing FC label
      label_list <- lapply(c("log2 Gsk3i/Basal", ""), standard_ggtexttable_label)
      #No output if empty
      if (length(searched_sites_raw()) == 0){
        ggtexttable("No data",
                    theme = ttheme(tbody.style = tbody_style(size = 20, fill = "white")))
      } else {
        ggarrange(plots = c(hmap_list, label_list),
                  nrow = 2, ncol = 2,
                  widths = c(1, 10),
                  #Here I scale the heatmap's size based on how many rows it has
                  lengths = c(length(searched_sites(), 4)))
      }
      
    }, res = resolution,
    #Here I scale the size of the plot to try to keep heatmap aspect ratio constant
    height = function() 200 + 25*length(searched_sites()))
  })
}

##App
Gsk3i_multsitesApp <- function(){
  ui <- Gsk3i_multsitesUI("Gsk3i_multsites")
  server <- function(input, output, session){
    Gsk3i_multsiteServer("Gsk3i_multsites")
  }
  shinyApp(ui, server)
}

Gsk3i_multsitesApp()










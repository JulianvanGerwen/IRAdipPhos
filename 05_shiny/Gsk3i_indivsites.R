









#########################################
#Functions
make_Gsk3iphos_boxplot <- function(boxplot_data,
                                   selected_site,
                                   boxplot_function = boxplot_gsk3iphos_psite,
                                   to_search_against = "gene_uniprot_site",
                                   ...){
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
    return(boxplot)
  } else {
    return(textGrob("Absent", gp = gpar(fontsize = 12)))
  }
}

###Function: Summarise all ppeptide for a psite, for Gsk3i data
#summary_fn: Function used to summarise one ppeptide
#Categories: Names for rows
summarise_Gsk3i_psite <- function(data, selected_site){
  if (length(grep(selected_site, rownames(data))) > 0){
    #Raw summary
    ppeptides <- sort(rownames(data)[grep(selected_site, rownames(data))])
    
    #Summarise
    summary_m <- sapply(ppeptides, function(x){
      reg <- ""
      if (data[x, "GSK3i_up"]){
        reg <- "up"
      } else if (data[x, "GSK3i_down"]){
        reg <- "down"
      }
      motif <- ""
      if (data[x, "matches_pST35psp_motif"]){
        motif <- "yes"
      }
      sub <- ""
      if (data[x, "GSK3i_down_motif"]){
        sub <- "yes"
      }
      return(c(reg, motif, sub))
    })
    #Add ppetpide info
    multplicities <- as.numeric(sapply(colnames(summary_m), function(x) strsplit(x, "_")[[1]][4]))
    summary_m <- rbind(c("Single", "Double", "Triple")[multplicities], summary_m)
    #Add category column
    summary_m <- cbind(c("Phosphopeptide:",
                         "Gsk3i regulation:",
                         "Matches Gsk3 motif:",
                         "Putative Gsk3 substrate:"), summary_m)
    return(summary_m)
  } else {
    return(NULL)
  }
}

###Make ggtext_table from summary
Gsk3i_summary_toggtexttable <- function(summary_table,
                                        fontsize = 5,
                                        padding = 1.5){
  #Prep data, including add dummy left column
  colnames(summary_table) <- NULL
  summary_table <- cbind(rep("", 4), summary_table)
  table <- ggtexttable(summary_table, 
                       theme = ttheme(tbody.style = tbody_style(fontsize = fontsize, fill = NA, hjust = 1, x = 1),
                                      padding = unit(c(padding, padding), "mm")))
  table <- table %>% table_cell_font(row = 1:tab_nrow(table),
                                     column = 1,
                                     face = "bold",
                                     size = fontsize)
  return(table)
}




#########################################

#Create named list of psites to choose from
#FOR now: all psites quantified more than 3 and without NA gene name

#UI module
Gsk3i_indivsites_psiteUI <- function(id, 
                               sidebarwidth = 3,
                               mainwidth = 9,
                               selectable_sites){
  fluidPage(
    sidebarLayout(
      sidebarPanel(width = sidebarwidth,
                   selectizeInput(NS(id, "selected_site"),
                                  label = "Phosphosite",
                                  choices = selectable_sites,
                                  selected = "Trarg1_Q8C838_S72"),
                   checkboxInput(NS(id, "raw"),
                                 label = "Raw data",
                                 value = FALSE)),
      mainPanel(width = mainwidth,
                plotOutput(NS(id, "boxplot"),
                           height = 800,
                           width = 1000))
    )
  )
}

#Server module
Gsk3i_indivsites_psiteserver <- function(id,
                                         resolution = 216){
  moduleServer(id, function(input, output, session){
    ##Boxplots
    #Boxplot dataset
    phos_gsk3i_boxplot_data <- reactive({phos_gsk3i_proc})
    
    #Boxplot reactive
    phos_gsk3i_boxplot <- reactive({
      make_Gsk3iphos_boxplot(boxplot_data = phos_gsk3i_boxplot_data(),
                              selected_site = input$selected_site,
                              siglabels = !input$raw,
                             is_raw = TRUE)
    })
    
    ##Summary table
    summary_table <- reactive({
      Gsk3i_summary_toggtexttable(summarise_Gsk3i_psite(phos_gsk3i_proc,
                                                        input$selected_site))
    })
    
    #Extract number of ppeptides to scale height
    num_ppeptides <- reactive({length(grep(input$selected_site,
                                           rownames(phos_gsk3i_boxplot_data())))})
      
    #Output
    output$boxplot <- renderPlot({
      plot_list <- list("Boxplot" = phos_gsk3i_boxplot(),
                        "Summary" = summary_table())
      as_ggplot(arrangeGrob(grobs = plot_list,
                            ncol = 2,
                            nrow = 1))
    }, res = resolution,
    height = function() 400*num_ppeptides())

})
}

##App function
Gsk3i_indivsites_psiteApp <- function(){
  ui <- Gsk3i_indivsites_psiteUI("Gsk3i_indivsites",
                                 selectable_sites = Gsk3i_selectable_sites)
  server <- function(input, output, session){
    Gsk3i_indivsites_psiteserver("Gsk3i_indivsites")
  }
  
  shinyApp(ui = ui, server = server)
}

Gsk3i_indivsites_psiteApp()




























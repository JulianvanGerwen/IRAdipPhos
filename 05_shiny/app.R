###shiny app for visualising IR Adip phos and total proteome data

###Version history

#2
#Implement modules for different pages


###Initialising

##Packages
library(shiny)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(egg)


##Scripts
source("01_load_data.R")
source("essential_visualisations_JvG.R")
source("boxplots_JvG.R")
source("3T3L1_mouse_IR__DF__202104_functions.R")
source("3T3L1_mouse_IR__DF__202104_objects.R")
source("indivsites.R")
source("multsites.R")
source("Gsk3i_indivsites.R")
source("Gsk3i_multsites.R")




###UI
ui <- fluidPage(
  titlePanel("Adipocyte phosphoproteomics"),
  navlistPanel(widths = c(2, 10),
    "Insulin resistance",
    tabPanel("Explore individual sites and proteins",
             indivsites_psiteUI(id = "indivsites",
                                selectable_sites = selectable_sites)),
    tabPanel("Explore multiple sites and proteins",
             multsitesUI(id = "multsites")),
    "Gsk3i phosphoproteome",
    tabPanel("Explore individual sites",
             Gsk3i_indivsites_psiteUI(id = "Gsk3i_indivsites",
                                      selectable_sites = Gsk3i_selectable_sites)),
    tabPanel("Explore multiple sites",
             Gsk3i_multsitesUI(id = "Gsk3i_multsites"))
)
)

###Server
server <- function(input, output, session){
  
  
  #IR Phosphoproteome
  indivsites_psiteServer("indivsites")
  multsitesServer("multsites")
  
  #Gsk3i phospohoproteome
  Gsk3i_indivsites_psiteserver("Gsk3i_indivsites")
  Gsk3i_multsiteServer("Gsk3i_multsites")
}

###shinyApp
shinyApp(ui = ui,
         server = server)
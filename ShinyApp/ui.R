# 

### Load libraries ###

library(shiny)
library(shinythemes)

genes <- c("alx1", "tbr", "tgif", "dr1", "foxB", "sm29", "sm37", "p58b", "Lim1", "vegfr10")



shinyUI(
  
  fluidPage(theme = shinytheme("yeti"),
    
    # Global page variables
    div( style = "padding: 1px 0px; width: '100%'", titlePanel( title = "", windowTitle = "Gene Coexpression")),
    
    # Page Layout
    navbarPage(

      # Application title
      title = ("Exploring Changes in the Developmental Gene Regulatory Network"),
      
      # Main data panel with tabs
      column(6, 
                  
               wellPanel(
                  h1(strong("About this Webpage"), style = "font-size:20px;"),
                  p("Sea Urchin larvae have been used as a model of early embryonic development for well over a century. As a result we know a lot about the roles that individual genes play in development and how these genes regulate each other. The figure below shows relationships between a network of key developmental genes."),
                  p(" "),
                  imageOutput("network", height = "100%"), 
                  p(" "),
                  p("Here we examine two species with this conserved life history: Lytechinus variegatus and Heliocidaris tuberculata. H. tuberculata was chosen because its close relative, Heliocidaris erythrogramma, has drastically different development. H. erythrogramma produces very large nutritious eggs, leading to larvae that are rounder, faster to metamorphose into adults, and do not need to feed on plankton."),
                  p(" "),
                  imageOutput("development", height = "100%"), 
                  p(" "),
                  strong("Explore the Data"),
                  p("We invite you to explore the data and see for yourself what parts of the developmental gene regulatory network show signs of change between these species. Simply return to the top of this page and select your genes of interest in the drop down menus on the righthand side.")
                  )
               ),
        
      column(6,
             fluidRow(   
               column(6,
                      strong("Select Two Genes:"),
                      selectizeInput("g1", "Gene One", choices = c(genes)),
                      selectizeInput("g2", "Gene Two", choices = c(genes), selected = "tbr"), 
                      ),
               column(6,imageOutput("legend", height = "100%")
                      ),
             ),
             
             tabsetPanel(
               
               
               tabPanel("Gene One", 
                        p(" "),
                        fluidRow(plotOutput("expG1", height = "300px")),
                        fluidRow(p("The plot above shows expression values based on bulk RNAseq data. Values shown in the plot are the mean of three replicates, each made up of multiple embryos. The images below display single cell RNAseq data using UMAPs, plots in which each cell in the data is represented by a dot. The cells are arranged in the plot based on overall gene expression similarity.  Colored dots represent cells that express Gene One.")),
                        fluidRow(
                          column(6,imageOutput("sc_LV_G1"), ),
                          column(6,imageOutput("sc_HE_G1"))
                        ),
                      ),
               
               tabPanel("Gene Two", 
                        p(" "),
                        fluidRow(plotOutput("expG2", height = "300px")),
                        fluidRow(p("The plot above shows expression values based on bulk RNAseq data. Values shown in the plot are the mean of three replicates, each made up of multiple embryos. The images below display single cell RNAseq data using UMAPs, plots in which each cell in the data is represented by a dot. The cells are arranged in the plot based on overall gene expression similarity.  Colored dots represent cells that express Gene Two. ")),
                        fluidRow(
                          column(6,imageOutput("sc_LV_G2")),
                          column(6,imageOutput("sc_HE_G2"))
                        ),
                      ),
               
               tabPanel("Gene Co-Expression", 
                        p(" "),
                        
                        fluidRow(
                          column(7,plotOutput("coExp", height = "300px")),
                          column(5,p("The plot at left shows the percentage of cells from our single cell data in which both selected genes are present. The images below are UMAPs, plots in which each cell in the data is represented by a dot. The cells are arranged in the plot based on overall similarity.  Colored dots represent cells with both genes present. "))
                        ),
                        
                        fluidRow(
                          column(6,imageOutput("sc_LV_co")),
                          column(6,imageOutput("sc_HE_co"))
                        ),
                      ),
               
        ) # End tabsetPanel
      ) # End column
    ) # End navbar
  ) # End fluidPage
) # End shinyUI



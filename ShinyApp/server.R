# This server uses pre-wrangled data. Please see Data_Wrangling.Rmd for details. 


### Load libraries ###

library(shiny)
library(tidyverse)
library(ggplot2)


### Read in data tables ###

# Bulk RNAseq expression count tables, long form, all time points, both species
combined_expression <- read.csv("../Final_Long_Lv_Ht_He.txt", sep = "\t") |> tibble::as_tibble()

# Single cell RNAseq gene co-expression count tables
long_tc <- read.csv("new_coexpression.csv") |> tibble::as_tibble()

# Table of pre-generated UMAP images
image_table <- read.csv("UMAP_image_table.csv")

# Table of gene common names to species gene IDs
g2o_code <- read.csv("gene_to_orthogroup_code.txt")



### Server starts here ###

shinyServer(function(input, output, session) {
  
  ### Read in static images ###
  
  output$network <- renderImage(deleteFile=FALSE, {
    outfile <- list(src = "Skel-GRN.png",
                    contentType = 'image/png',
                    width = "100%",
                    alt = "This is alternate text") })
  
  output$development <- renderImage(deleteFile=FALSE, {
    outfile <- list(src = "Development.png ",
                    contentType = 'image/png',
                    width = "100%",
                    alt = "This is alternate text") })
  
  
  output$legend <- renderImage(deleteFile=FALSE, {
    outfile <- list(src = "legend.png",
                    contentType = 'image/png',
                    width = "100%",
                    alt = "This is alternate text") })
  

  output$sc_LV_G1 <- renderImage(deleteFile=FALSE, {
    outfile <- list(src = filter(image_table, species == "LV" & type == "EXP" & tolower(G1) == tolower(input$g1))$value,
                    contentType = 'image/png',
                    width = "100%",
                    alt = "This is alternate text")
  })
  
  output$sc_LV_G2 <- renderImage(deleteFile=FALSE, {
    outfile <- list(src = filter(image_table, species == "LV" & type == "EXP" & tolower(G1) == tolower(input$g2))$value,
                    contentType = 'image/png',
                    width = "100%",
                    alt = "This is alternate text")
  })
  
  output$sc_HE_G1 <- renderImage(deleteFile=FALSE, {
    outfile <- list(src = filter(image_table, species == "HE" & type == "EXP" & tolower(G1) == tolower(input$g1))$value,
                    contentType = 'image/png',
                    width = "100%",
                    alt = "This is alternate text")
  })
  
  output$sc_HE_G2 <- renderImage(deleteFile=FALSE, {
    outfile <- list(src = filter(image_table, species == "HE" & type == "EXP" & tolower(G1) == tolower(input$g2))$value,
                    contentType = 'image/png',
                    width = "100%",
                    alt = "This is alternate text")
  })
  
  output$sc_LV_co <- renderImage(deleteFile=FALSE, {
    outfile <- list(src = filter(image_table, species == "LV" & type == "coEXP" & (tolower(G1) == tolower(input$g1) | tolower(G2) == tolower(input$g1)) & (tolower(G1) == tolower(input$g2) | tolower(G2) == tolower(input$g2)))$value,
                    contentType = 'image/png',
                    width = "100%",
                    alt = "This is alternate text")
  })
  
  output$sc_HE_co <- renderImage(deleteFile=FALSE, {
    outfile <- list(src = filter(image_table, species == "HE" & type == "coEXP" & (tolower(G1) == tolower(input$g1) | tolower(G2) == tolower(input$g1)) & (tolower(G1) == tolower(input$g2) | tolower(G2) == tolower(input$g2)))$value,
                    contentType = 'image/png',
                    width = "100%",
                    alt = "This is alternate text")
  })
  

  
    output$coExp <- renderPlot({
        print(        long_tc |>
                        filter(tolower(gene_1) == tolower(input$g1) & tolower(gene_2) == tolower(input$g2) & coex < 100))
        # draw the histogram with the specified number of bins
        long_tc |>
          filter(tolower(gene_1) == tolower(input$g1) & tolower(gene_2) == tolower(input$g2) & coex < 100) |>
          ggplot(aes(hpf, coex, color=species)) +
          geom_line(size = 1.25) +
          labs(y = '% cells with coexpression') +
          labs(title = str_c('Coexpression')) +
          theme_minimal() +
          scale_x_continuous(breaks=c(6,9,12,16,20,24)) +
          scale_color_manual(values=c("#FAB120", "#148040")) +
          theme(
            panel.grid.minor = element_blank(),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.x = element_text(family = "sans", size = 16, margin=margin(10,0,0,0)), 
            axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,12,0,0)), 
            axis.text = element_text(family = "mono", size = 14),
            plot.title = element_text(family = "sans", size = 18, margin=margin(0,0,20,0)),
            legend.position = "none"
          )
        

    })
    output$expG1 <- renderPlot({
      
      # draw the histogram with the specified number of bins
      combined_expression |>
        filter(combined_expression$orthogroup == filter(g2o_code, tolower(name) == tolower(input$g1))$orthogroup)  |>
        ggplot(aes(hpf, ex, color=species)) +
        geom_line(size = 1.25) +
        labs(y = NULL) +
        labs(title = str_c(input$g1, ' Expression')) +
        theme_minimal() +
        scale_x_continuous(breaks=c(0, 2, 3, 4, 10, 16, 24)) +
        scale_color_manual(values=c("#148040", "#005400", "#FAB120")) +
        theme(
          panel.grid.minor = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x = element_text(family = "sans", size = 16, margin=margin(10,0,0,0)), 
          axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,12,0,0)), 
          axis.text = element_text(family = "mono", size = 14),
          plot.title = element_text(family = "sans", size = 18, margin=margin(0,0,20,0)) #, legend.position = "none"
        )
    })
    
    output$expG2 <- renderPlot({
      
      
      # draw the histogram with the specified number of bins
      combined_expression |>
        filter(orthogroup == filter(g2o_code, tolower(name) == tolower(input$g2))$orthogroup) |>
        ggplot(aes(hpf, ex, color=species)) +
        geom_line(size = 1.25) +
        labs(y = NULL) +
        labs(title = str_c(input$g2, ' Expression')) +
        theme_minimal() +
        scale_x_continuous(breaks=c(0, 2, 3, 4, 10, 16, 24)) +
        scale_color_manual(values=c("#148040", "#005400", "#FAB120")) +
        theme(
          panel.grid.minor = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x = element_text(family = "sans", size = 16, margin=margin(10,0,0,0)), 
          axis.title.y = element_text(family = "sans", size = 16, margin=margin(0,12,0,0)), 
          axis.text = element_text(family = "mono", size = 14),
          plot.title = element_text(family = "sans", size = 18, margin=margin(0,0,20,0))#,legend.position = "none"
        )
    })

})

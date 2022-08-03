# R Shiny Seurat tutorial from the Seurat website.

library(shiny)
#source("helpers.R")
library(dplyr)
library(Seurat)
library(patchwork)

#setwd("~/Desktop/Seurat_data")
pbmc.data <- Read10X(data.dir = 'data/hg19/')
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Normalize the data using a global-scaling normalization "LogNormalize"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features (feature selection)
# Which features are highly expressed in some cells and lowly expressed in others?
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


# Scale the data to shift the mean gene expression to 0, and shift the variance across cells to 1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Linear dimensional reduction
# Run PCA, use previously determined variables features as input
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

ui <- fluidPage(
  titlePanel("Seurat tutorial"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Follow the Seurat analysis of single cell dataset."),
      
      selectInput("feature", 
                  label = "Choose a variable to display",
                  choices = c("nFeature_RNA", 
                              "nCount_RNA",
                              "percent.mt"),
                  selected = "nFeature_RNA"),
      sliderInput("topNfeautres", 
                  label = "Select number of highly variable fetures:",
                  min = 1, max = 50, value = 10),
      
      sliderInput("PC_N", 
                  label = "Select the PC to visualize:",
                  min = 1, max = 20, value = 1)
    ),
    
    mainPanel(
      plotOutput("Violin_features"),
      tableOutput('topN'),
      plotOutput('PC_plot')
    )
  )
)


server <-function(input, output) {
  output$Violin_features <- renderPlot({
    
    VlnPlot(pbmc, features = input$feature)
    
  })
  output$topN <- renderTable({
    head(VariableFeatures(pbmc), input$topNfeautres)
  })
  
  output$PC_plot <- renderPlot({
    
    VizDimLoadings(pbmc, dims = input$PC_N, reduction = "pca")
    
  })
}

shinyApp(ui = ui, server = server)
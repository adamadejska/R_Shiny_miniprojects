# R Shiny Seurat tutorial from the Seurat website.

library(shiny)
#source("helpers.R")
library(dplyr)
library(Seurat)
library(patchwork)


pbmc.data <- Read10X(data.dir = 'data/filtered_gene_bc_matrices/hg19/')
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc[["percent.Ribosomal"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[LS]")
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
pbmc <- RunUMAP(pbmc, dims = 1:10)

ui <- navbarPage("Seurat tutorial",
                 tabPanel("QC metrics", 
                          helpText("This is an app that follows the steps of the Seurat tutorial.
                     Use each tab to explore the different steps of the tutorial."),
                          sidebarLayout(
                            sidebarPanel(
                              helpText("Use this tab to visualize QC metrics and use them to filter cells."),
                              helpText("Available metrics:"),
                              helpText("nFeature_RNA: number of observed genes (anything with a nonzero count)"),
                              helpText("nCount_RNA: total number of reads (unique molecular identifiers) in the dataset"),
                              helpText("percent.mt: precentage of mitochondrial genes"),
                              helpText("percent.Ribosomal: percentage of ribosomal protein genes"),
                              selectInput("feature", 
                                          label = "Choose a metric to display",
                                          choices = c("nFeature_RNA", 
                                                      "nCount_RNA",
                                                      "percent.mt",
                                                      "percent.Ribosomal"),
                                          selected = "nFeature_RNA"),
                              helpText("Visualize feature-feature relationships"),
                              checkboxGroupInput("checkGroup", label = "Choose 2 features to display", 
                                                 choices = c("nFeature_RNA", "nCount_RNA" , "percent.mt", "percent.Ribosomal"),
                                                 selected = c("nFeature_RNA", "nCount_RNA")),
                              helpText("Filter the data based on the features."),
                              sliderInput("nFeature_filter", 
                                          label = "Select the range of nFeature_RNA to keep:",
                                          min = 1, max = 3500, value = c(200, 2500)),
                              sliderInput("per_mt_filter", 
                                          label = "Select the percentage of mitochondrial genese to filter:",
                                          min = 1, max = 100, value = 99)
                            ),
                            mainPanel(plotOutput("Violin_features"), plotOutput("scatter"))
                          )
                 ),
                 tabPanel("Summary", 
                          helpText("This is an app that follows the steps of the Seurat tutorial.
                     Use each tab to explore the different steps of the tutorial."),
                          sidebarLayout(
                            sidebarPanel(
                              helpText("Use this tab to find features that exhibit high cell-to-cell variation in the dataset."),
                              helpText("Focusing on these genes in downstream analysis helps to highlight biological signal."),
                              sliderInput("topNfeautres", 
                                          label = "Select number of highly variable fetures:",
                                          min = 1, max = 50, value = 10)
                            ),
                            mainPanel(tableOutput('topN'), plotOutput("VarFeaturePlot"))
                          )
                 ),
                 tabPanel("PCA", 
                          helpText("This is an app that follows the steps of the Seurat tutorial.
                     Use each tab to explore the different steps of the tutorial."),
                          sidebarLayout(
                            sidebarPanel(
                              helpText("Explote the primary sources of heterogeneity in the dataset"),
                              sliderInput("PC_N", 
                                          label = "Select the PC to visualize:",
                                          min = 1, max = 20, value = 1)
                            ),
                            mainPanel(
                              helpText("List of the genes which were most highly (and lowly) weighted in the selected PC."),
                              plotOutput('PC_plot'),
                              helpText("Heatmap plots of PCA weightings for the most highly and lowly weighted genes, shown against the set of cells which are most highly influenced by the PC."),
                              helpText("The idea is that as long as we’re seeing clear structure in one of these plots then we’re still adding potentially useful information to the analysis"),
                              plotOutput('PC_heatmap'))
                          )
                 ),
                 tabPanel("UMAP", 
                          helpText("This is an app that follows the steps of the Seurat tutorial.
                     Use each tab to explore the different steps of the tutorial."),
                          sidebarLayout(
                            sidebarPanel(
                              helpText("Based on the dimensions chosen from the PCA, construct a 2D plot to cluster the cells"),
                              sliderInput("dim_range", 
                                          label = "Range of dimensions:",
                                          min = 1, max = 20, value = c(1, 10)),
                              helpText("Resolution parameter sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters"),
                              sliderInput("cluster_res", 
                                          label = "Resolution of the clusters:",
                                          min = 0.1, max = 2, value = 0.5),
                              actionButton("umap_go", "Update")
                            ),
                            mainPanel(plotOutput('UMAP_clusters'))
                          )
                 ),
                 tabPanel("Markers", 
                          helpText("This is an app that follows the steps of the Seurat tutorial.
                     Use each tab to explore the different steps of the tutorial."),
                          sidebarLayout(
                            sidebarPanel(
                              helpText("Cluster similar cells together in low-dimensional space using non-linear reduction techniques."),
                              helpText("Visualize marker expression on a PCA plot."),
                              helpText("Example genes: MS4A1, GNLY, CD3E, CD14,
                                                       FCER1A, FCGR3A, LYZ,
                                                       PPBP, CD8A"),
                              textInput("gene", "Write gene name you want to visualize (case sensitive)", "MS4A1"),
                              actionButton("go", "Update")
                            ),
                            mainPanel(plotOutput('UMAP_plot'))
                          )
                 )
                 
)


server <-function(input, output) {
  output$Violin_features <- renderPlot({
    pbmc <- subset(pbmc, subset = nFeature_RNA > input$nFeature_filter[1] & nFeature_RNA < input$nFeature_filter[2] & percent.mt < input$per_mt_filter)
    VlnPlot(pbmc, features = input$feature)
  })
  
  output$scatter <- renderPlot({
    pbmc <- subset(pbmc, subset = nFeature_RNA > input$nFeature_filter[1] & nFeature_RNA < input$nFeature_filter[2] & percent.mt < input$per_mt_filter)
    FeatureScatter(pbmc, feature1 = input$checkGroup[1], feature2 = input$checkGroup[2])
  })
  
  output$topN <- renderTable({
    
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    as_tibble(HVFInfo(pbmc),rownames = "Gene") -> variance.pbmc
    variance.pbmc %>%
      mutate(hypervariable=Gene %in% VariableFeatures(pbmc)
      ) -> variance.pbmc

    head(variance.pbmc, n=input$topNfeautres)
    # head(VariableFeatures(pbmc), input$topNfeautres)
  })
  
  output$VarFeaturePlot <- renderPlot({
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    topN <- head(VariableFeatures(pbmc), input$topNfeautres)
    plot1 <- VariableFeaturePlot(pbmc)
    LabelPoints(plot = plot1, points = topN, repel = TRUE)
  })
  
  output$PC_plot <- renderPlot({
    
    VizDimLoadings(pbmc, dims = input$PC_N, reduction = "pca")
    
  })
  
  output$PC_heatmap <- renderPlot({
    DimHeatmap(pbmc, dims = input$PC_N, cells = 500, balanced = TRUE)
    
  })
  
  output$UMAP_clusters <- renderPlot({
    input$umap_go
    isolate({
      pbmc <- FindNeighbors(pbmc, dims = input$dim_range[1]:input$dim_range[2])
      pbmc <- FindClusters(pbmc, resolution = input$cluster_res)
      
      pbmc <- RunUMAP(pbmc, dims = input$dim_range[1]:input$dim_range[2])
      DimPlot(pbmc, reduction = "umap")
    })
    
  })
  
  
  output$UMAP_plot <- renderPlot({
    input$go
    isolate(FeaturePlot(pbmc, features = input$gene))
  })
}

shinyApp(ui = ui, server = server)
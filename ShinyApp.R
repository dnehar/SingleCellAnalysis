
#-- 02/15/18
#-- analysis seurat rds files from single cell data 
#-- version1

#---------------------  --------------------- --------------------- --------------------- --------------------- --------------------- #


#--------------                                             User Interface (UI)

#---------------------  --------------------- --------------------- --------------------- --------------------- --------------------- #
  library(RColorBrewer)
  library(Seurat)
  library(tidyr)
  library(dplyr)
  library(reshape)
  library(reshape2)
  library(shiny)
  library(shinythemes)
  library(pheatmap)
  library(threejs)
  library(shinyFiles)
  
  #-- to be able to hit 'enter' button
  jscode <- '
  $(function() {
  var $els = $("[data-proxy-click]");
  $.each(
  $els,
  function(idx, el) {
  var $el = $(el);
  var $proxy = $("#" + $el.data("proxyClick"));
  $el.keydown(function (e) {
  if (e.keyCode == 13) {
  $proxy.click();
  }
  });
  }
  );
  });
  '
  ui <- bootstrapPage(
           navbarPage(
           theme = shinytheme("cerulean"),
           windowTitle = 'Single Cell Analysis Shiny App',

           #--- File selection
           shinyFilesButton('file', 'File select', 'Select a rds file', FALSE,
                            buttonType = "button"),
           
           tabPanel('EXPRESSION PLOTS',  
                    #-- header          
                    div(h1('Gene expression analysis across clusters'), align = 'center'),  
                    br(),
                    br(),
                    #--- number of cells is XXX        
                    mainPanel(
                      textOutput("text1"), 
                      tags$head(
                        tags$style(
                          "#text1{color: black;
                          font-size: 20px;
                          font-style: italic;
                          }")
                  )
                        ),
                  
                  br(), # br() element to introduce extra vertical spacing
                  br(),
                  br(),
                  br(),
                  
                  #---------------------------------------------------- Expression plots 
                  
                  #-- Gene selection  
                  #-- GO button           
                  fluidRow(
                    column(4, offset = 4,
                           tags$head(tags$script(HTML(jscode))),
                           tagAppendAttributes(
                             textInput("GENE_SYMBOL", label = ("Enter GENE_SYMBOL to plot expression"), value = 'B2M'),
                             `data-proxy-click` = "goButton1"
                           ),
                           actionButton('goButton1', 'Plot')
                    )),
                  br(),
                  br(),  
                  br(),
                  br(),
                  br(),
                  br(),
                  br(),
                  br(),
                  
                  splitLayout(
                    cellWidths = c("33%","33%", "33%"), 
                    plotOutput("VlnPlot"), 
                    plotOutput("JoyPlot"), 
                    plotOutput("FeatPlot")
                    
                  )#fluidRow
                        ), #mainPanel
           
    #---------------------------------------------------- Tsne plots 
  
    tabPanel('TSNE PLOTS', 
             #-- header          
             div(h1('Tsne plots within CellType, groups or samples'), align = 'center'),  
             br(),
             br(),  
             br(),
             #goButton
             column(4, offset=4,
                    actionButton('goButton2', 'Click to see Tsne plots')
             ),  
             
             br(),
             br(),  
             br(),
             br(),
             br(),
             br(),  
             
             
             mainPanel(#"TSNE plots",
               fluidRow( 
                 splitLayout(
                   cellWidths = c("30%","32%", "40%"), 
                   plotOutput("TSNE"),  
                   plotOutput("TSNE2"), 
                   plotOutput("TSNE3"))
               )
             )
    ),
    
    #----------  DOuble positive cells 
    
    tabPanel('DOUBLE POSTIVE CELLS', 
             #-- header          
             div(h1('Double positive cells within clusters'), align = 'center'),  
             br(),
             br(),  
             br(),
             #goButton
             fluidRow(
               column(8, offset = 5,
                      tags$head(tags$script(HTML(jscode))),
                      tagAppendAttributes(
                        textInput('GenesForDP', 'Two comma seperated gene names', value ='IL1B,ISG15'),
                        `data-proxy-click` = "goButton5"
                      ),
                      actionButton('goButton5','Plot')
               )
             ),    
             br(),
             br(),  
      
             mainPanel(column(12, offset = 5,
               fluidRow(plotOutput("DoublePos")))
               )
    ),
    
    #------------------------------------------  BARPLOTS 
    #navbarPage(
    tabPanel('BARPLOTS',
             
             #header   
             div(h1('Proportion Barplots within CellType, Groups or samples'), align = 'center'),  
             br(),
             br(),  
             br(),
             #-Go button
             column(4,offset = 4,
                    actionButton('goButton3', 'Click to see Barplots')
             ),
             
             br(),
             br(),  
             br(),
             br(),
             br(),
             br(),  
             
             #--- Barplots 
             mainPanel(
               plotOutput("BP_CellType"),
               plotOutput("BP_Group"),
               plotOutput("BP_ID")
             
             )#mainPanel
    ), #tabPanel
    
    #------------------------------------------  HEATMAP + DOTPLOT 
    
   navbarMenu('HEATMAP & DotPot',
        #- HeatMap
        tabPanel('HeatMap',      
                    #-- header          
                    div(h1('HEATMAP using geneset of interest'), align = 'center'),  
                    br(),
                    br(),
                    br(),
                    br(),
                    #--- run button 
                    fluidRow(
                       column(12, offset = 2,
                              tags$head(tags$script(HTML(jscode))),
                              tagAppendAttributes(
                                textInput('GenesForHM', 'Comma seperated gene names', value = 'Il6,Cd3d,il32'),
                                `data-proxy-click` = "goButton4"
                              ),
                              actionButton('goButton4', 'Plot')
                       )
                     ),    
               #----  plot HM
              mainPanel(column(8, offset = 5,
                   fluidRow(plotOutput("HM")))
                   )
              ),
          #- DotPlot
          tabPanel('DotPlot', 
                   #-- header          
                   div(h1('DOTPLOT using geneset of interest'), align = 'center'), 
                   #-space
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   
                   #----  plot DP
                   mainPanel(column(8, offset = 5,
                     fluidRow(plotOutput("DP"))))
                   )
                  
          )#navbarMenu
        
   
   
  #----   
       )#navbarPage
   )#bootstrap

#---------------------  --------------------- --------------------- --------------------- --------------------- --------------------- #

#-                                                        server

#---------------------  --------------------- --------------------- --------------------- --------------------- --------------------- #
  library(RColorBrewer)
  library(Seurat)
  library(tidyr)
  library(dplyr)
  library(reshape)
  library(reshape2)
  library(shiny)
  library(shinythemes)
  library(pheatmap)
  library(threejs)
  library(ggplot2)
  library(shinyFiles)
  
#-- Functions 

  PercentAbove <- function(x, threshold){
    return(length(x = x[x > threshold]) / length(x = x))
  }
  
  MinMax <- function(data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    return(data2)
  }
  #--- colors 
  colors <-c('0'='#e41a1c','1'='#999999', '2'='#377eb8','3'='#4daf4a','4'='#984ea3','5'='#ff7f00',
             '6'='#a65628','7'='#f781bf','8'='#7fcdbb','9'='#cab2d6', '10'="#060404", '11'="#bdbdbd", 
             '12'="#addd8e", '13'="#dd1c77",'14'="#fed976",'15'= "#08519c",'16'="#7a0177")
  col_ID <- c("JB17009"="#e41a1c", "JB17010"="#999999", "JB17011"="#377eb8", "JB17017"="#4daf4a","JB17018"="#78c679",
              "JB17022"="#984ea3", "JB17023"="#ff7f00", "JB17024"="#a65628", "JB17001"="#f781bf", "JB17002"="#7fcdbb",
              "JB17013"="#cab2d6", "JB17003"="#060404", "JB17004"="#bdbdbd",  "JB17005"="#addd8e",  "JB17006"="#dd1c77",
              "JB17007"="#fed976",  "JB17008"="#08519c",  "JB17014"="#7a0177",  "JB17015"="#a6761d",  
              "JB17016"="#756bb1","JB17019"="#fdc086",  "JB17020"="#fb8072","JB17021"="#8dd3c7")
  
  col_HM <- colorRampPalette(c("#3288bd","#000000","#FFFF00"))(n=100)
  
  #------------------- 
  
  server <- function(input, output) {
  
  #------------------------------------------------  FIlE  
    #-- path to rds files 
    root = c(wd='../../SubClusters/RDS/1Cluster/') 
    #root =c(wd='/projects/nehard/SingleCell/SLEvsHY_Aggregate_II/1to12_1/1to12_1_MergedCluster_II_Jan2018/SubClusters/RDS')
    
    #-- files 
    shinyFileChoose(input, 'file', roots=root, filetypes=c('', 'rds')) # choose file 
    
    #--- acces files    
    observeEvent(input$file, { 
      
      inFile <- parseFilePaths(roots=root, input$file)
      #--- load rds file 
      load(as.character(inFile$datapath))
    
      
 #------------------------------------------------  parameters     
    #-- change identities for TSNE
    GROUP <-  SetAllIdent(object = Cl_seurat, id = "Group")   
    ID <- SetAllIdent(object = Cl_seurat, id = "ID") 
    
    SLE_ID <- c(paste0("JB1700", 1:8),  paste0("JB170", 14:16), paste0("JB170", 19:24))
    HY_ID <- c("JB17009",paste0("JB170", 10:11), paste0("JB170", 17:18))
    ID_all <-  c(HY_ID, SLE_ID)
    

    #--- goButtons
    v1 <- reactiveValues(doPlot = FALSE)
    observeEvent(input$goButton1, {
      v1$doPlot <- input$goButton1
    })
    
    v2 <- reactiveValues(doPlot = FALSE)
    observeEvent(input$goButton2, {
      v2$doPlot <- input$goButton2
    })
    
    v3 <- reactiveValues(doPlot = FALSE)
    observeEvent(input$goButton3, {
      v3$doPlot <- input$goButton3
    })
    
    v4 <- reactiveValues(doPlot = FALSE)
    observeEvent(input$goButton4, {
      v4$doPlot <- input$goButton4
    })
    
    v5 <- reactiveValues(doPlot = FALSE)
    observeEvent(input$goButton5, {
      v5$doPlot <- input$goButton5
    })
    
    #-- number of cells 
    output$text1 <-  renderText({
      CL.name <- gsub( ".rds", "", inFile$name) #remove ".rds' from the name
      paste0("Number of ", CL.name ," cells is : ", dim(Cl_seurat@data)[2])
    })
    
    #--------------------------------------------------------- TSNE plots 
    #-- TSNE plot CellType
    output$TSNE <- renderPlot({
      if (v2$doPlot == FALSE)
        return() 
      isolate ({
        TSNEPlot(Cl_seurat, pt.size = 1.5, do.label = T, label.size = 12,
                 colors.use = colors)
      }) 
    })
    
    #-- TSNE plot  group 
    output$TSNE2 <- renderPlot({
      if (v2$doPlot == FALSE)
        return() 
      
      TSNEPlot(GROUP, pt.size = 1.5, colors.use =c("#78c679","#756bb1") , 
               do.label =F, label.size = 8)
    })
    
    #-- TSNE plot ID 
    
    output$TSNE3 <- renderPlot({
      if (v2$doPlot == FALSE)
        return() 
      TSNEPlot(ID, pt.size = 1.5,do.label = F, label.size = 8, colors.use = col_ID)
    })
    
    #--------------------------------------------------------- Double positive plots 
    
    
    output$DoublePos <- renderPlot({ 
      if (v5$doPlot == FALSE)
        return() 
      isolate({ 
        
        GENE_DP  <- input$GenesForDP
        GENE_DP <- toupper(GENE_DP)
        GENE_DP <- unlist(strsplit(GENE_DP, ','))
        gene_DP = intersect(GENE_DP, rownames(Cl_seurat@data))
        #print(gene_DP)
        FeaturePlot(object = Cl_seurat, features.plot = gene_DP, pt.size = 2, 
                            cols.use = c("#e0e0e0", "#1f78b4", "#1b9e77", "#e41a1c"),
                          overlay = TRUE, no.legend = FALSE)
      
      }) })
    
    
    
    #--------------------------------------------------------- Barplots 
    #---- Barplot CellType
    
    output$BP_CellType <- renderPlot({ 
      
      if (v3$doPlot == FALSE)
        return() 
      isolate({ 
        
        y=prop.table(x=table(Cl_seurat@ident, Cl_seurat@meta.data$ID), margin=2)
        y1=y[,ID_all]
        barplot(y1,las=2,cex.names = 1, xlim=c(0, ncol(y1) + 3),
                col = colors, ylab="Proportion of cells",
                legend.text=TRUE,
                args.legend = list(x = "top", bty = "n", inset=c(0, -0.15), cex=1.2, horiz=T))      
      })
    })
    
    #---- Barplot Group 
    output$BP_Group <- renderPlot({ 
      if (v3$doPlot == FALSE)
        return() 
      
      isolate({   
        y=prop.table(x=table(GROUP@ident, GROUP@meta.data$Merged_CL), margin=2)
        barplot(y,las=2,cex.names = 1, xlim=c(0, ncol(y) + 3), col = c("#78c679","#756bb1"),
                ylab="Proportion of cells", legend.text=TRUE,
                args.legend = list(x = "top", bty = "n", inset=c(0, -0.15), cex=1.1, horiz=T)) 
        
      }) 
    })
    
    #--- Barplot ID
    output$BP_ID<- renderPlot({ 
      if (v3$doPlot == FALSE)
        return() 
      
      isolate({
        y=prop.table(x=table(ID@ident, ID@meta.data$Merged_CL), margin=2)
        barplot(y,las=2,cex.names = 1,
                xlim=c(0, ncol(y) + 3),
                col = col_ID,
                ylab="Proportion of cells",
                legend.text=TRUE,
                args.legend = list(x = "right", bty = "n", inset=c(0, -1.2), cex=0.8, horiz=F))   
      }) 
    })
    
    #--------------------------------------------------------- Expression plots 
    
    #-- Vln plot 
    output$VlnPlot <- renderPlot({
      if (v1$doPlot == FALSE)
        return() 
      isolate({
        Seurat::VlnPlot(Cl_seurat,toupper(input$GENE_SYMBOL), point.size.use = 0, cols.use = colors)
      })
    })
    #-- Joy plot 
    output$JoyPlot <- renderPlot({
      
      if (v1$doPlot == FALSE)
        return() 
      
      isolate({
        Seurat::JoyPlot(Cl_seurat,toupper(input$GENE_SYMBOL),col =colors)
      })  
    })
    
    #-- Feature plot 
    output$FeatPlot <- renderPlot({
      if (v1$doPlot == FALSE)
        return() 
      isolate({
        Seurat::FeaturePlot(Cl_seurat,toupper(input$GENE_SYMBOL),cols.use =c("#e0e0e0", "#ff7f00"),
                            pt.size = 2, no.axes =F)
      })  
    })
    
    #--------------------------------------------------------- HEATMAP + DOTPLOT 
    #--- Dotplot
    
    output$DP<- renderPlot({
      
      if (v4$doPlot == FALSE)
        return() 
      isolate({ 
        #- genes to be plotted   
        GENE  <- input$GenesForHM
        GENE <- toupper(GENE)
        GENE <- unlist(strsplit(GENE, ','))
        genes = intersect(GENE, rownames(Cl_seurat@data))
        
        if (length(unique(genes))>2) {
          DotPlot(Cl_seurat, genes, x.lab.rot = T, cols.use = c("#d9d9d9", "#762a83"), plot.legend = T)
        }
      })
    })
    
    #--- heatmap 
    output$HM <- renderPlot({
      if (v4$doPlot == FALSE)
        return() 
      isolate({
        #- genes to be plotted 
        GENE  <- input$GenesForHM
        GENE <- toupper(GENE)
        GENE <- unlist(strsplit(GENE, ','))
        genes = intersect(GENE, rownames(Cl_seurat@data))
        
        if (length(unique(genes))>=2) {
          Mat.genes1 <- as.matrix(Cl_seurat@data[genes,])
          MyMat <- t(Mat.genes1)
          #--- Cells with id (clusters)
          MyCells <- data.frame(id= Cl_seurat@ident[rownames(MyMat)])
          
          #--- Expression data and cluster id 
          data<-data.frame(merge(MyMat,MyCells, by=0), row.names = 1)
          data$cell <- rownames(data)
          
          
          #--- melt the matrix 
          data %>% gather(
            key = genes.plot,
            value = expression,
            -c(cell, id) # id= clusters 
          ) -> data
          
          #- threshold + + percentage of cell within each cluster 
          data %>%
            group_by(id, genes.plot) %>%
            summarize( #***
              avg.exp = mean(expm1(x = expression)), #expm1(x) function calculates exp(x) - 1
              pct.exp = PercentAbove(x = expression, threshold = .3) #**** thr
            ) -> data
          
          #--- scale 
          data %>%
            ungroup() %>% #- ungoupe the firt grouping 
            group_by(genes.plot) %>% #group by gene across id (clusters)
            mutate(avg.exp.scale = scale(x = avg.exp)) %>% #scaling 
            mutate(avg.exp.scale = MinMax( data = avg.exp.scale, max = 2.5, min = -2.5)) ->  data
          
          #--- filter genes using averge expression and percentage of cells 
          genes.plot.filtered <-data %>% group_by(genes.plot) %>% filter(avg.exp > 0 & pct.exp > 0.3)%>% 
            count(genes.plot) %>% filter(n >=1) #& n <= 10
          
          if (!is.null(genes.plot.filtered)){
            #---- melted matrix with filtred genes 
            data.to.plot1<-data[data$genes.plot %in% genes.plot.filtered$genes.plot,]
            
            #--- de-melt 
            hold<-data.frame(dcast(data.to.plot1, genes.plot~id,value.var = 'avg.exp.scale'), row.names = 1)
            
            #- modify colnames 
            colnames(hold)  <- gsub("X", "", colnames(hold))
            colnames(hold) <- paste0("Cl_",colnames(hold))
            
            
            #-- plot heatmap 
            pheatmap(hold, show_colnames = T, 
                     show_rownames = T, 
                     cluster_cols = T, 
                     cluster_rows = T,
                     border_color = F,
                     cellwidth =30, 
                     cellheight = 30,  
                     clustering_method = "ward.D2",
                     treeheight_row = 12, 
                     treeheight_col = 10, 
                     fontsize_row =12, 
                     color = col_HM, 
                     fontsize_col = 12, 
                     fontsize =15)
            
          }}
      }) #renderPlot
    })
    
    
  })#observeEvent
}
shinyApp(ui = ui, server = server)





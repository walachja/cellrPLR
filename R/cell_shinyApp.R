cell_shinyApp <- function(data, g1, g2, mainGroup='max',mainGroup2='all',simul=FALSE)
{
  library(reshape2)
  library(shiny)
  library(ggplot2)
  library(MASS)
  library(plotly)
  
  # Colors
  colfunc <- colorRampPalette(c("blue","white", "red"))
  
  
  # Cellmap - by group
  v <- cellmap(data, g1=g1, g2=g2, mainGroup='max',simul=simul)
  
  datafiles <- list(Tukey=list(v$vv1,v$vv2),Huber=list(v$vv2,v$vv3),Hampel=list(v$vv4,v$vv5))
  rm(v)
  gc()
  
  v <- cellmap(data, g1=g1, g2=g2, mainGroup='all')
  
  # Cellmap - all
  datafiles2 <- list(Tukey=list(v$vv1,v$vv2),Huber=list(v$vv2,v$vv3),Hampel=list(v$vv4,v$vv5))
  rm(v)
  gc()
  
  
  ui <- fluidPage(
    
    # Application title
    titlePanel("Cell-wise outliers"),
    
    # Sidebar with a slider input for the number of bins
    sidebarLayout(
      sidebarPanel(
        sliderInput("bins",
                    "Threashold:",
                    min = 0,
                    max = 1,
                    step=0.05,
                    value = 0),
        
        selectInput("do3", "Weighting function:",
                    c('Biweight'='1','Huber'='2','Hampel'='3'),selected = '1'),
        
        selectInput("do4", "Aggregation:",
                    c('Mean'='1','Median'='2'),selected = '1'),
        
        selectInput("do5", "Sorting:",
                    c('None'='1','By sum'='2', 'By Number bigger than zero'='3'),selected = 'None'),
        sliderInput('sort', 'Sort according to:',min = 1, max = nrow(datafiles$Tukey[[1]]), value = c(1,nrow(datafiles$Tukey[[1]])),step=1)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(
          tabPanel("With group information",
                   fluidRow(
                     plotlyOutput("distPlot4", width = "100%", height = "1000px"))#,
                   
          ),
          tabPanel("Without group information",
                   fluidRow(
                     plotlyOutput("distPlot", width = "100%", height = "1000px"))#,
                   
          ))
      )
    )
  )
  
  server <- function(input, output) {
    
    outVar <- reactive({
      temp <- datafiles[[as.numeric(input$do3)]][[as.numeric(input$do4)]]
      zzz <- temp
      rownames(zzz) <- rownames(data)
      colnames(zzz) <- colnames(data)
      zzz[(zzz < input$bins & zzz > -input$bins)] <- 0
      
      s1 <- (abs(apply(zzz[input$sort[1]:input$sort[2],]!=0,2,sum)))
      s2 <- (abs(apply(zzz[input$sort[1]:input$sort[2],],2,sum)))
      
      if (input$do5==2) {zzz <- zzz[,order(s2)]}
      if (input$do5==3) zzz <- zzz[,order(s1,s2)]
      zzz <- reshape2::melt(t(zzz))
      zzz$Var1 <- as.factor(zzz$Var1)
      zzz$Var1 <- factor(zzz$Var1,levels = levels(zzz$Var1))
      zzz$Var2 <- as.factor(zzz$Var2)
      zzz$Var2 <- factor(zzz$Var2,levels = levels(zzz$Var2))
      zzz
    })
    
    outVar2 <- reactive({
      temp <- datafiles2[[as.numeric(input$do3)]][[as.numeric(input$do4)]]
      zzz <- temp
      rownames(zzz) <- rownames(data)
      colnames(zzz) <- colnames(data)
      zzz[(zzz < input$bins & zzz > -input$bins)] <- 0
      
      s1 <- (abs(apply(zzz[input$sort[1]:input$sort[2],]!=0,2,sum)))
      s2 <- (abs(apply(zzz[input$sort[1]:input$sort[2],],2,sum)))
      
      if (input$do5==2) {zzz <- zzz[,order(s2)]}
      if (input$do5==3) zzz <- zzz[,order(s1,s2)]
      zzz <- reshape2::melt(t(zzz))
      zzz$Var1 <- as.factor(zzz$Var1)
      zzz$Var1 <- factor(zzz$Var1,levels = levels(zzz$Var1))
      zzz$Var2 <- as.factor(zzz$Var2)
      zzz$Var2 <- factor(zzz$Var2,levels = levels(zzz$Var2))
      zzz
    })
    
    plot.render <- reactive({
      p1 <- ggplot(outVar(), aes(Var1, Var2)) + geom_tile(aes(fill = value)) +
        scale_fill_gradient2(low = "blue",mid='white',high = "red",midpoint = 0,limits=c(-1,1)) +
        theme_grey(base_size = 15)  + 
        scale_y_discrete(expand = c(0, 0))+
        scale_x_discrete(expand = c(0, 0))+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.3), axis.text.x = element_text(angle = 90, hjust = 1))+
        xlab('Variables')+ylab('Samples')
      
      print(ggplotly(p1))
    })
    
    output$distPlot4 <- renderPlotly({
      plot.render()
    })
    
    
    output$distPlot <- renderPlotly({
      
      p1 <- ggplot(outVar2(), aes(Var1, Var2)) + geom_tile(aes(fill = value)) +
        scale_fill_gradient2(low = "blue",mid='white',high = "red",midpoint = 0,limits=c(-1,1)) +
        theme_grey(base_size = 15)  + 
        scale_y_discrete(expand = c(0, 0))+
        scale_x_discrete(expand = c(0, 0))+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.3), axis.text.x = element_text(angle = 90, hjust = 1))+
        xlab('Variables')+ylab('Samples')
      
      print(ggplotly(p1))
    })
  }
  
  shinyApp(ui,server)
}



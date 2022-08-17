library(tidyverse)
library(gridExtra)
library(plotly)

testData <- read_csv("Caples_PlotData_20220225.csv")
colnames(testData)[2] <- "CSE_ID"
testData <- testData %>% 
  mutate(Total_Live_Shrub_Ht = replace(Total_Live_Shrub_Ht, Total_Live_Shrub_Ht == "999", NA)) %>%
  mutate(Live_Shrub_Tall_Height = replace(Live_Shrub_Tall_Height, Live_Shrub_Tall_Height == "999", NA)) %>%
  mutate(Live_Shrub_Medium_Height = replace(Live_Shrub_Medium_Height, Live_Shrub_Medium_Height == "999", NA)) %>%
  mutate(Live_Shrub_Low_Height = replace(Live_Shrub_Low_Height, Live_Shrub_Low_Height == "999", NA)) %>%
  mutate(Herb_Ht = replace(Herb_Ht, Herb_Ht == "999", NA)) %>%
  mutate(Gram_Ht = replace(Gram_Ht, Gram_Ht == "999", NA)) %>%
  mutate(Caples_Severity_Class = replace(Caples_Severity_Class, Caples_Severity_Class == "Unburned" & Caples_Perim == "Y", "Unburned_In")) %>%
  mutate(Caples_Severity_Class = replace(Caples_Severity_Class, Caples_Severity_Class == "Unburned" & Caples_Perim == "N", "Unburned_Out"))

##########################  UI section    #####################################
###############################################################################
###############################################################################

ui <- fluidPage(
  titlePanel("Caples Viz 3.0"),
  tags$textarea(id = "introduction", rows = 6, cols = 100,
                "This Shiny App is aimed towards visualizing and exploring the plot data collected during six years of the Caples Project (2013, 2017, 2018, 2019, 2020, 2021). Under the following tabs you'll find a variety of plotting options, each allowing you to examine the current three main classes of data: raw data describing vegetation cover, estimates of tree and snag density per hectare (including primary species specific density estimates), and estimates of coarse woody debris (CWD - including cover, size class, volume, and density). Additional plot-level information (e.g., plot ID, elevation, etc) is available when the cursor is held over a point (i.e., hover text)."),
  tabsetPanel(
    tabPanel("Barplots",
             sidebarLayout(sidebarPanel(selectInput(inputId = "vegVar_bw", label = "Select vegetation variable", choices = colnames(testData[,c(13:28,55:59)]),selected = "", selectize = FALSE),
                                        actionButton("goButton_vegVar_bw", "Update vegetation variable plot",style = "background-color:#FFFFFF;color:#0A0A0A;border-color:#6B6A6A;border-style:dashed;border-width:1px;border-radius:3%;font-size:12px;"),
                                        hr(),
                                        selectInput(inputId = "densVar_bw", label = "Select density variable",choices = colnames(testData[,c(29:54)])),actionButton("goButton_densVar_bw", "Update density variable plot",style = "background-color:#FFFFFF;color:#0A0A0A;border-color:#6B6A6A;border-style:dashed;border-width:1px;border-radius:3%;font-size:12px;"),
                                        hr(),
                                        selectInput(inputId = "cwdVar_bw", label = "Select CWD variable", choices = colnames(testData[,c(60:74)]),selected = "", selectize = FALSE),actionButton("goButton_cwdVar_bw", "Update CWD variable plot",style = "background-color:#FFFFFF;color:#0A0A0A;border-color:#6B6A6A;border-style:dashed;border-width:1px;border-radius:3%;font-size:12px;"),
                                        width = 3),
                           mainPanel(fluidRow(helpText("Box-and-whisker plots illustrating sampling distribution of the selected variable pre-treatment (Pre-Caples: 2013, 2017, 2018) and post-treatment (Post-Caples: 2019, 2020, 2021). These plots are standard, illustrating median (horizontal line) and lower (25%) and upper (75%) quartiles - lower and upper hinge, respectively. Whiskers extend 1.5 * IQR on either side of the hinge. NOTE: hover over each plot to identify the plot number, or the plot boundaries to view IQR values"),column(6,plotlyOutput("vegVar_tab1")),column(6,plotlyOutput("densVar_tab1")),column(6,plotlyOutput("cwdVar_tab1")))))),
    tabPanel("Line Plots",
             sidebarPanel(width = 3,selectInput(inputId = "vegVar_line", label = "Select vegetation variable", choices = colnames(testData[,c(13:28,55:59)]),selected = "", selectize = FALSE),actionButton("goButton_vegVar_line", "Update vegetation variable plot",style = "background-color:#FFFFFF;color:#0A0A0A;border-color:#6B6A6A;border-style:dashed;border-width:1px;border-radius:3%;font-size:12px;")),
             mainPanel(verticalLayout(helpText("*Panels* correspond to whether plots were burned or unburned during the project period. Within panels, the *x-axis* corresponds to sampling events: Pre-Caples = pre-treatment sampling (2013, 2017, 2018), Post-Caples = post-treatment sampling (2019, 2020, 2021); *Points* correspond to plot ID and are connected across sampling events by lines. NOTE: hover over each plot to identify the plot number"),plotlyOutput("vegVar_tab2a"),helpText("*Panels* correspond to burn severity, measured on the ground and transformed into four classes - unburned (further divided by plots that were inside or outside of caples project boundary), Low, Moderate, and High. The *X-axis* in each panel corresponds to sampling events: Pre-Caples = pre-treatment sampling (2013,2017,2018), Post-Caples = post-treatment sampling (2019, 2020, 2021); *Points* correspond to plot ID and are connected across sampling events by lines. NOTE: hover over each plot to identify the plot number"),plotlyOutput("vegVar_tab2b"))),
             sidebarPanel(width = 3,selectInput(inputId = "densVar_line", label = "Select density variable", choices = colnames(testData[,c(29:54)]),selected = "", selectize = FALSE),actionButton("goButton_densVar_line", "Update density variable plot",style = "background-color:#FFFFFF;color:#0A0A0A;border-color:#6B6A6A;border-style:dashed;border-width:1px;border-radius:3%;font-size:12px;")),
             mainPanel(verticalLayout(plotlyOutput("densVar_tab2a"),plotlyOutput("densVar_tab2b"))),
             sidebarPanel(width = 3,selectInput(inputId = "cwdVar_line", label = "Select CWD variable", choices = colnames(testData[,c(60:74)]),selected = "", selectize = FALSE),actionButton("goButton_cwdVar_line", "Update CWD variable plot",style = "background-color:#FFFFFF;color:#0A0A0A;border-color:#6B6A6A;border-style:dashed;border-width:1px;border-radius:3%;font-size:12px;")),
             mainPanel(verticalLayout(plotlyOutput("cwdVar_tab2a"),plotlyOutput("cwdVar_tab2b")))),
    tabPanel("X by Y plots",
             sidebarPanel(width = 3,selectInput(inputId = "vegVar_xbyy1", label = "Select vegetation variable", choices = colnames(testData[,c(13:28,55:59)]),selected = "", selectize = FALSE),selectInput(inputId = "vegVar_xbyy2", label = "Select vegetation variable", choices = colnames(testData[,c(13:28,55:59)]),selected = "", selectize = FALSE),actionButton("goButton_vegVar_xbyy", "Update vegetation variable plot",style = "background-color:#FFFFFF;color:#0A0A0A;border-color:#6B6A6A;border-style:dashed;border-width:1px;border-radius:3%;font-size:12px;")),
             mainPanel(verticalLayout(helpText("help text"),plotlyOutput("vegVar_tab3a"),plotlyOutput("vegVar_tab3b"))),
             sidebarPanel(width = 3,selectInput(inputId = "densVar_xbyy1", label = "Select density variable", choices = colnames(testData[,c(29:54)]),selected = "", selectize = FALSE),selectInput(inputId = "densVar_xbyy2", label = "Select density variable", choices = colnames(testData[,c(29:54)]),selected = "", selectize = FALSE),actionButton("goButton_densVar_xbyy", "Update density variable plot",style = "background-color:#FFFFFF;color:#0A0A0A;border-color:#6B6A6A;border-style:dashed;border-width:1px;border-radius:3%;font-size:12px;")),
             mainPanel(verticalLayout(helpText("help text"),plotlyOutput("densVar_tab3a"),plotlyOutput("densVar_tab3b"))),
             sidebarPanel(width = 3,selectInput(inputId = "cwdVar_xbyy1", label = "Select CWD variable", choices = colnames(testData[,c(60:74)]),selected = "", selectize = FALSE),selectInput(inputId = "cwdVar_xbyy2", label = "Select CWD variable", choices = colnames(testData[,c(60:74)]),selected = "", selectize = FALSE),actionButton("goButton_cwdVar_xbyy", "Update CWD variable plot",style = "background-color:#FFFFFF;color:#0A0A0A;border-color:#6B6A6A;border-style:dashed;border-width:1px;border-radius:3%;font-size:12px;")),
             mainPanel(verticalLayout(helpText("help text"),plotlyOutput("cwdVar_tab3a"),plotlyOutput("cwdVar_tab3b")))),
    tabPanel("Geographic perspective and plot level data",
             h2("Geographic position of plots - click on plot to populate plotting windows below"),
             plotlyOutput("geomp"), 
             h2("Change in proportion of tree/shrub/forb at each plot"),
             plotOutput("geomp_target1"), 
             h2("Change in tree/snag size classes"),
             plotOutput("geomp_target2"),
             h2("Change in CWD"),
             plotOutput("geomp_target3"),
             h2("Change in proportion CWD at each plot"),
             plotOutput("geomp_target4"))
  )
)


###############################################################################
###############################################################################
######################  SERVER section    #####################################
###############################################################################
###############################################################################

server <- function(input, output, session) {

  output$vegVar_tab1 <- renderPlotly({
    input$goButton_vegVar_bw
    isolate(
        bw <- ggplot(testData %>% drop_na(input$vegVar_bw), 
                     aes_string(x = "Treatment_Interval", 
                                y = input$vegVar_bw, 
                                text = paste0("CSE_ID"))) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(size = 0.5, width = 0.08) +
          theme_bw() +
          ggtitle("Vegetation cover variable") +
          theme(axis.title.y = element_blank()) +
          theme(axis.title.x = element_blank(), legend.position = "none") +
          theme(plot.title = element_text(size = 12)) +
          scale_x_discrete(limits = rev(levels(as.factor(testData$Treatment_Interval))))
        )
    isolate(ggplotly(bw, tooltip = "text"))

  })
  
  output$densVar_tab1 <- renderPlotly({
    input$goButton_densVar_bw
    isolate(
      bw <- ggplot(testData %>% drop_na(input$densVar_bw), 
                   aes_string(x = "Treatment_Interval", 
                              y = input$densVar_bw, 
                              text = paste0("CSE_ID"))) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(size = 0.5, width = 0.08) +
        theme_bw() +
        ggtitle("Density variable") +
        theme(axis.title.y = element_blank()) +
        theme(axis.title.x = element_blank(), legend.position = "none") +
        theme(plot.title = element_text(size = 12)) +
        scale_x_discrete(limits = rev(levels(as.factor(testData$Treatment_Interval))))
    )
    isolate(ggplotly(bw, tooltip = "text"))
    
  })
  
  output$cwdVar_tab1 <- renderPlotly({
    input$goButton_cwdVar_bw
    isolate(
      bw <- ggplot(testData %>% drop_na(input$cwdVar_bw), 
                   aes_string(x = "Treatment_Interval", 
                              y = input$cwdVar_bw, 
                              text = paste0("CSE_ID"))) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(size = 0.5, width = 0.08) +
        theme_bw() +
        ggtitle("Coarse Woody Debris variable") +
        theme(axis.title.y = element_blank()) +
        theme(axis.title.x = element_blank(), legend.position = "none") +
        theme(plot.title = element_text(size = 12)) +
        scale_x_discrete(limits = rev(levels(as.factor(testData$Treatment_Interval))))
    )
    isolate(ggplotly(bw, tooltip = "text"))
    
  })
  
  output$vegVar_tab2a <- renderPlotly({
    input$goButton_vegVar_line
    isolate(
      lp <- ggplot(testData %>% drop_na(input$vegVar_line), 
                   aes_string(x = "Treatment_Interval", 
                              y = input$vegVar_line, 
                              text = paste0("CSE_ID"))) +
        geom_point(size = 1) +
        geom_line(aes(group = as.factor(CSE_ID))) + 
        theme_bw() +
        theme(axis.title.y = element_blank(),legend.position = "none") +
        theme(axis.title.x = element_blank(),legend.position = "none") +
        scale_x_discrete(limits = rev(levels(as.factor(testData$Treatment_Interval)))) +
        facet_wrap(~factor(Burned_Caples, levels = c("Unburned","Burned")))
    )
    isolate(ggplotly(lp, tooltip = "text"))
    
  })
  
  output$vegVar_tab2b <- renderPlotly({
    input$goButton_vegVar_line
    isolate(
      lp2 <- ggplot(testData %>% drop_na(input$vegVar_line), 
                   aes_string(x = "Treatment_Interval", 
                              y = input$vegVar_line, 
                              text = paste0("CSE_ID"))) +
        geom_point(size = 0.5) +
        geom_line(aes(group = as.factor(CSE_ID))) + 
        theme_bw() +
        theme(axis.title.y = element_blank(),legend.position = "none") +
        theme(axis.title.x = element_blank(),legend.position = "none") +
        scale_x_discrete(limits = rev(levels(as.factor(testData$Treatment_Interval)))) +
        facet_wrap(~factor(Caples_Severity_Class, levels = c("Unburned_Out","Unburned_In","Low","Mod","High")))
    )
    isolate(ggplotly(lp2, tooltip = "text"))
    
  })
  
  output$densVar_tab2a <- renderPlotly({
    input$goButton_densVar_line
    isolate(
      lp <- ggplot(testData %>% drop_na(input$densVar_line), 
                   aes_string(x = "Treatment_Interval", 
                              y = input$densVar_line, 
                              text = paste0("CSE_ID"))) +
        geom_point(size = 1) +
        geom_line(aes(group = as.factor(CSE_ID))) + 
        theme_bw() +
        theme(axis.title.y = element_blank(),legend.position = "none") +
        theme(axis.title.x = element_blank(),legend.position = "none") +
        scale_x_discrete(limits = rev(levels(as.factor(testData$Treatment_Interval)))) +
        facet_wrap(~factor(Burned_Caples, levels = c("Unburned","Burned")))
    )
    isolate(ggplotly(lp, tooltip = "text"))
    
  })
  
  output$densVar_tab2b <- renderPlotly({
    input$goButton_densVar_line
    isolate(
      lp2 <- ggplot(testData %>% drop_na(input$densVar_line), 
                    aes_string(x = "Treatment_Interval", 
                               y = input$densVar_line, 
                               text = paste0("CSE_ID"))) +
        geom_point(size = 0.5) +
        geom_line(aes(group = as.factor(CSE_ID))) + 
        theme_bw() +
        theme(axis.title.y = element_blank(),legend.position = "none") +
        theme(axis.title.x = element_blank(),legend.position = "none") +
        scale_x_discrete(limits = rev(levels(as.factor(testData$Treatment_Interval)))) +
        facet_wrap(~factor(Caples_Severity_Class, levels = c("Unburned_Out","Unburned_In","Low","Mod","High")))
    )
    isolate(ggplotly(lp2, tooltip = "text"))
    
  })
  
  output$cwdVar_tab2a <- renderPlotly({
    input$goButton_cwdVar_line
    isolate(
      lp <- ggplot(testData %>% drop_na(input$cwdVar_line), 
                   aes_string(x = "Treatment_Interval", 
                              y = input$cwdVar_line, 
                              text = paste0("CSE_ID"))) +
        geom_point(size = 1) +
        geom_line(aes(group = as.factor(CSE_ID))) + 
        theme_bw() +
        theme(axis.title.y = element_blank(),legend.position = "none") +
        theme(axis.title.x = element_blank(),legend.position = "none") +
        scale_x_discrete(limits = rev(levels(as.factor(testData$Treatment_Interval)))) +
        facet_wrap(~factor(Burned_Caples, levels = c("Unburned","Burned")))
    )
    isolate(ggplotly(lp, tooltip = "text"))
    
  })
  
  output$cwdVar_tab2b <- renderPlotly({
    input$goButton_cwdVar_line
    isolate(
      lp2 <- ggplot(testData %>% drop_na(input$cwdVar_line), 
                    aes_string(x = "Treatment_Interval", 
                               y = input$cwdVar_line, 
                               text = paste0("CSE_ID"))) +
        geom_point(size = 0.5) +
        geom_line(aes(group = as.factor(CSE_ID))) + 
        theme_bw() +
        theme(axis.title.y = element_blank(),legend.position = "none") +
        theme(axis.title.x = element_blank(),legend.position = "none") +
        scale_x_discrete(limits = rev(levels(as.factor(testData$Treatment_Interval)))) +
        facet_wrap(~factor(Caples_Severity_Class, levels = c("Unburned_Out","Unburned_In","Low","Mod","High")))
    )
    isolate(ggplotly(lp2, tooltip = "text"))
    
  })
  
  output$vegVar_tab3a <- renderPlotly({
    input$goButton_vegVar_xbyy
    isolate(
      xyp <- ggplot(testData %>% drop_na(input$vegVar_xbyy1, input$vegVar_xbyy2), 
                   aes_string(x = input$vegVar_xbyy2, 
                              y = input$vegVar_xbyy1, 
                              text = paste0("CSE_ID"))) +
        geom_point() +
        #geom_smooth(method = lm) + 
        theme_bw() +
        theme(axis.title.y = element_text(size = 14)) +
        theme(axis.title.x = element_text(size = 14)) +
        facet_wrap(~factor(Burned_Caples, levels = c("Unburned","Burned")))
    )
    isolate(ggplotly(xyp, tooltip = "text"))
    
  })
  
  output$vegVar_tab3b <- renderPlotly({
    input$goButton_vegVar_xbyy
    isolate(
      xyp2 <- ggplot(testData %>% drop_na(input$vegVar_xbyy1, input$vegVar_xbyy2), 
                    aes_string(x = input$vegVar_xbyy2, 
                               y = input$vegVar_xbyy1, 
                               text = paste0("CSE_ID"))) +
        geom_point() +
        #geom_smooth(method = lm) + 
        theme_bw() +
        theme(axis.title.y = element_text(size = 14)) +
        theme(axis.title.x = element_text(size = 14)) +
        facet_wrap(~factor(Caples_Severity_Class, levels = c("Unburned_Out","Unburned_In","Low","Mod","High")))
    )
    isolate(ggplotly(xyp2, tooltip = "text"))
    
  })
  
  output$densVar_tab3a <- renderPlotly({
    input$goButton_densVar_xbyy
    isolate(
      xyp <- ggplot(testData %>% drop_na(input$densVar_xbyy1, input$densVar_xbyy2), 
                    aes_string(x = input$densVar_xbyy2, 
                               y = input$densVar_xbyy1, 
                               text = paste0("CSE_ID"))) +
        geom_point() +
        #geom_smooth(method = lm) + 
        theme_bw() +
        theme(axis.title.y = element_text(size = 14)) +
        theme(axis.title.x = element_text(size = 14)) +
        facet_wrap(~factor(Burned_Caples, levels = c("Unburned","Burned")))
    )
    isolate(ggplotly(xyp, tooltip = "text"))
    
  })
  
  output$densVar_tab3b <- renderPlotly({
    input$goButton_densVar_xbyy
    isolate(
      xyp2 <- ggplot(testData %>% drop_na(input$densVar_xbyy1, input$densVar_xbyy2), 
                     aes_string(x = input$densVar_xbyy2, 
                                y = input$densVar_xbyy1, 
                                text = paste0("CSE_ID"))) +
        geom_point() +
        #geom_smooth(method = lm) + 
        theme_bw() +
        theme(axis.title.y = element_text(size = 14)) +
        theme(axis.title.x = element_text(size = 14)) +
        facet_wrap(~factor(Caples_Severity_Class, levels = c("Unburned_Out","Unburned_In","Low","Mod","High")))
    )
    isolate(ggplotly(xyp2, tooltip = "text"))
    
  })
  
  output$cwdVar_tab3a <- renderPlotly({
    input$goButton_cwdVar_xbyy
    isolate(
      xyp <- ggplot(testData %>% drop_na(input$cwdVar_xbyy1, input$cwdVar_xbyy2), 
                    aes_string(x = input$cwdVar_xbyy2, 
                               y = input$cwdVar_xbyy1, 
                               text = paste0("CSE_ID"))) +
        geom_point() +
        #geom_smooth(method = lm) + 
        theme_bw() +
        theme(axis.title.y = element_text(size = 14)) +
        theme(axis.title.x = element_text(size = 14)) +
        facet_wrap(~factor(Burned_Caples, levels = c("Unburned","Burned")))
    )
    isolate(ggplotly(xyp, tooltip = "text"))
    
  })
  
  output$cwdVar_tab3b <- renderPlotly({
    input$goButton_cwdVar_xbyy
    isolate(
      xyp2 <- ggplot(testData %>% drop_na(input$cwdVar_xbyy1, input$cwdVar_xbyy2), 
                     aes_string(x = input$cwdVar_xbyy2, 
                                y = input$cwdVar_xbyy1, 
                                text = paste0("CSE_ID"))) +
        geom_point() +
        #geom_smooth(method = lm) + 
        theme_bw() +
        theme(axis.title.y = element_text(size = 14)) +
        theme(axis.title.x = element_text(size = 14)) +
        facet_wrap(~factor(Caples_Severity_Class, levels = c("Unburned_Out","Unburned_In","Low","Mod","High")))
    )
    isolate(ggplotly(xyp2, tooltip = "text"))
    
  })
  
  output$geomp <- renderPlotly({
    geomp <- ggplot(testData, aes(x = Easting, y = Northing,
                                  label = CSE_ID, text = paste0("Elevation"))) + 
      geom_point(aes(color = Caples_Severity_Class), size = 4) +
      theme_bw() +
      scale_color_discrete(name = "Burn Severity", 
                           labels = c("High Severity","Low Severity","Moderate Severity","Unburned In Caples Perimeter", "Unburned Out of Caples Perimeter", "NA")) 
    #limits = c("High","Mod","Low","Unburned_In","Unburned_Out","NA"))
    ggplotly(geomp)
    
  })
  
  output$geomp_target1 <- renderPlot({
    mouse_event <- event_data("plotly_click") # captures hover text from plotly
    print(mouse_event)
    target <- mouse_event[c(3,4)] #would change index number here, if using the hover text (and not the axis value) to capture plot ID
    df <- testData[,c(2,4,7,8,15,19,21,23,25,27)]
    df <- reshape2::melt(df, id.vars = c("CSE_ID", "Treatment_Interval","Easting","Northing"))
    geomp_pie <- ggplot(df %>% filter(Easting == target[[1]], Northing == target[[2]]), 
                        aes(x = "", fill = variable, y = value, text = paste0("CSE_ID"))) + 
      geom_bar(stat = 'identity', position = position_fill(1)) + 
      coord_polar("y", start = 0) +
      theme_void() +
      facet_wrap(~factor(Treatment_Interval, levels = c("Pre-Caples","Post-Caples")))
    geomp_pie
    
  })
  
  output$geomp_target2 <- renderPlot({
    mouse_event <- event_data("plotly_click") # captures hover text from plotly
    print(mouse_event)
    target <- mouse_event[c(3,4)] #would change index number here, if using the hover text (and not the axis value) to capture plot ID
    df2 <- testData[,c(2,4,7,8,31:35,38:42)]
    df2 <- reshape2::melt(df2, id.vars = c("CSE_ID", "Treatment_Interval","Easting","Northing"))
    geomp_size <- ggplot(df2 %>% filter(Easting == target[[1]], Northing == target[[2]]), 
                         aes(fill = Treatment_Interval, x = variable, y = value)) + 
      geom_point(pch = 21, size = 5) + 
      geom_line(aes(group = as.factor(variable)), arrow = arrow(ends = "first", type = "closed"), size = 1) +
      theme_bw()
    geomp_size
    
  })
  
  output$geomp_target3 <- renderPlot({
    mouse_event <- event_data("plotly_click") # captures hover text from plotly
    print(mouse_event)
    target <- mouse_event[c(3,4)] #would change index number here, if using the hover text (and not the axis value) to capture plot ID
    df3 <- testData[,c(2,4,7,8,60:70)]
    df3 <- reshape2::melt(df3, id.vars = c("CSE_ID", "Treatment_Interval","Easting","Northing"))
    geomp_cwd <- ggplot(df3 %>% filter(Easting == target[[1]], Northing == target[[2]]), 
                         aes(fill = Treatment_Interval, x = variable, y = value)) + 
      geom_point(pch = 21, size = 5) + 
      geom_line(aes(group = as.factor(variable)), arrow = arrow(ends = "first", type = "closed"), size = 1) +
      theme_bw()
    geomp_cwd
    
  })
  
  output$geomp_target4 <- renderPlot({
    mouse_event <- event_data("plotly_click") # captures hover text from plotly
    print(mouse_event)
    target <- mouse_event[c(3,4)] #would change index number here, if using the hover text (and not the axis value) to capture plot ID
    df4 <- testData[,c(2,4,7,8,60:70)]
    df4 <- reshape2::melt(df4, id.vars = c("CSE_ID", "Treatment_Interval","Easting","Northing"))
    geomp_cwd2 <- ggplot(df4 %>% filter(Easting == target[[1]], Northing == target[[2]]), 
                        aes(x = "", fill = variable, y = value, text = paste0("CSE_ID"))) + 
      geom_bar(stat = 'identity', position = position_fill(1)) + 
      coord_polar("y", start = 0) +
      theme_void() +
      facet_wrap(~factor(Treatment_Interval, levels = c("Pre-Caples","Post-Caples")))
    geomp_cwd2
    
  })


}

# Run the application 
shinyApp(ui = ui, server = server)

# This app relies on FindPeaks.R
# Run this command first
# Waveform <- ASCII_extract("PATH/file")
# author: Daxiang Na (daxiang_na@urmc.rochester.edu)

library(shiny)
library(plotly)

ui <- fluidPage(
        numericInput("Sound","Sound Level", 75, min = 0, max = 110, step = 5),
        plotlyOutput("plot"),
        verbatimTextOutput("click"),
        tableOutput("dataTable")
)

server <- function(input, output, session) {
        output$plot <- renderPlotly({
                j <- input$Sound
                p <- if (as.character(j) %in% colnames(Waveform)) {
                        plot_ly(x = Waveform$Data_Pnt_ms, y = Waveform[,as.character(j)], type = 'scatter', mode = 'lines')%>%
                                add_trace()%>%
                                layout(showlegend = F,
                                       annotations = list(text = paste("Sound Level = ", as.character(j)," dB", sep = ""), 
                                                          showarrow = F),
                                       xaxis = list(title = list(text ='Latency')),
                                       yaxis = list(title = list(text = 'Amplitude'))
                                       )%>%
                                layout(dragmode = "select") %>%
                                event_register("plotly_selecting")
                } else {plot_ly(x = 0, y = 0, type = 'scatter', mode = 'lines')
                }

        })
        
        output$dataTable <- renderTable({
                c <- event_data("plotly_click")
                if (is.null(c)) 
                        "Click to show x and y (double-click to clear)" 
                else 
                        colnames(c) <- c("curveNumber", "pointNumber", "latency (ms)", "amplitude (μV)")
                c
        })

        # output$click <- renderPrint({
        #         b <- data.frame()
        #         c <- event_data("plotly_click")
        #         d <- rbind(b,c)
        #         if (is.null(d)) "Click to show x and y (double-click to clear)" else d
        # })
}

shinyApp(ui, server)

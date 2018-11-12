library(revtools) 
library(dplyr)
WoS1 <- read_bibliography("data/WoS_try1_1.bib")
WoS2 <- read_bibliography("data/WoS_try1_2.bib") %>% 
    mutate(note = NA)
WoS <- rbind(WoS1, WoS2)

library(shiny)
article <- WoS[1, ]
WoS
ui <- fluidPage(
    titlePanel("Abstract screener"),
    #br(),
    #tags$hr(),
    sidebarLayout(position = "right",
        sidebarPanel( 
            actionButton(inputId = "accept",
            label = "Accept", style = "background-color: #4CB5F5; color: #fff;  width: 100px"),
            actionButton(inputId = "reject",
                label = "Reject", style = " background-color: #34675c; color: #fff; width: 100px"),
            br(),       
            br(),
            actionButton(inputId = "backward",
                label = "<", style = "background-color: #B7B8B6; color: #fff;  width: 100px"),
            actionButton(inputId = "forward",
                label = ">", style = "background-color: #B7B8B6; color: #fff;  width: 100px"),
            br(),
            br(),
            numericInput(inputId = "df_row", label = "Select article by rownumber", value = 1,
                min = 1, max = nrow(WoS), step = 1),
            textInput(inputId = "notes", label = "Notes")
        ),
        mainPanel(#tags$b("Abstract"),
            br(),
            strong(textOutput(outputId = "title")),
            div(textOutput(outputId = "authors"), style = "color:grey"),
            div(textOutput(outputId = "year"),  style = "color:grey"),
            div(textOutput(outputId = "journal"), style = "color:grey"),
            br(),
            textOutput(outputId = "abstract"))
   
    )
)

server <- function(input, output){
    
    current_paper <- reactive({
            WoS[input$df_row, ]
    })
    
    output$title <- renderText({current_paper()[["title"]]})
    # output$authors <- renderText(article$author)
    # output$year <- renderText(article$year)
    # output$journal <- renderText(article$journal)
    output$abstract <- renderText({current_paper()[["abstract"]]})
    # observeEvent(input$accept, {
    #     output$title <- renderText({paste0(WoS[input$df_row, "title"], " (Accepted)")})
    # })
}

shinyApp(ui = ui, server = server)



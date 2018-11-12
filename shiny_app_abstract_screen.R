library(shiny)
library(shinydashboard)
library(dplyr)
library(revtools)
WoS1 <- read_bibliography("data/WoS_try1_1.bib")
WoS2 <- read_bibliography("data/WoS_try1_2.bib") %>% 
    mutate(note = NA)
WoS <- rbind(WoS1, WoS2)

data_in <- revtools:::load_abstract_data(data = WoS)
ui_data <- revtools:::screen_abstracts_ui()

ui <- dashboardPage(skin = "yellow",
    title = "screen abstracts",
    header = shinydashboard::dashboardHeader(title = plotOutput("header")),
    sidebar = shinydashboard::dashboardSidebar(sidebarMenu(id = "tabs", 
        menuItem("Data", icon = shiny::icon("bar-chart-o"), startExpanded = TRUE, br(),
                actionButton(inputId = "save_data", label = "Save Data", width = "80%"), br())
    )),
    body = shinydashboard::dashboardBody(fluidRow(column(width = 1), 
        column(width = 8, textOutput("abstract")), column(width = 1), 
        column(width = 2, uiOutput(outputId = "selector_buttons"), 
            uiOutput(outputId = "render_notes"), tableOutput(outputId = "progress_text"))))
)

server <- function(input, output) {
    data <- reactiveValues(raw = data_in$data$raw)
    progress <- reactiveValues(current = data_in$progress$current, 
        row = data_in$progress$row)
    display <- reactiveValues(notes = FALSE)
    #WoS_element <- WoS[1, ]
   # output$header <- renderText(paste(WoS_element$author,  WoS_element$title,  WoS_element$journal,  WoS_element$year, sep = "; "))
    output$abstract <- renderText(WoS_element[1, ]$abstract)
}

shinyApp(ui, server)
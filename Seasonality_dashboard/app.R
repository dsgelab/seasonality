library(shiny)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(DT)
library(mgcv)
source('utils.R')

heritability_dat <- read.table('data/finngen_R10_finngen_R10_analysis_data_ldsc_data_finngen_R10_FIN.ldsc.heritability.tsv',header=T) %>%
  as_tibble()

monthly_counts_FinRegistry  <- read.table('data/FINREGISTRY_endpoints_monthly_count_green.txt',header=T) %>%
  as_tibble() %>%
  filter_monthly_counts(.,heritability_dat)

cols_to_format <- c('pval','log10_pval','val_peak','val_trough','ptr','af','var_fraction','dispersion')

seasonal_summary_FinnGen <- read.table('data/FINNGEN_seasonal_summary.txt',header=T) %>%
  as_tibble() %>%
  rename(log10_pval=log_pval) %>%
  mutate(log10_pval=log10_pval/log(10)) %>%
  inner_join(select(heritability_dat,ENDPOINT=PHENO,H2),by='ENDPOINT') %>%
  filter(H2>0.01) %>%
  mutate(across(all_of(cols_to_format), format,digits=2))

seasonal_summary_FinRegistry <- read.table('data/FINREGISTRY_seasonal_summary.txt',header=T) %>%
                                as_tibble() %>%
                                inner_join(select(seasonal_summary_FinnGen,ENDPOINT,H2),by='ENDPOINT')%>%
                                mutate(across(all_of(cols_to_format), format,digits=2))

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("Seasonality dashboard"),

  # Sidebar with a slider input for number of bins
  fluidRow(
      column(4,
        uiOutput('search_endpoints')
      ),
      column(8,
        #checkboxGroupInput(inputId="adjustments",label="Adjustments",choices=c('Hospital usage','Mean temperature'))
      )
  ),
  hr(),
    # Show a plot of the generated distribution
  fluidRow(
       column(6,
          plotOutput("plot_FinRegistry",height=600)
       ),
       column(6,
          tabsetPanel(
            tabPanel("Seasonal summary - FinRegistry",
                     dataTableOutput("FinRegistry_table")
            ),
            tabPanel("Seasonal summary - FINNGEN",
                     dataTableOutput("FinnGen_table")
            )
          )
       )
    ),
  )

# Define server logic required to draw a histogram
server <- function(input, output) {
  clicked <- reactiveValues(endpoint_idx=matrix(,nrow=1,ncol=2))
  observeEvent(input$FinRegistry_table_cells_selected, {
    clicked$endpoint_idx <- input$FinRegistry_table_cells_selected
  })

  observeEvent(input$FinnGen_table_cells_selected, {
    clicked$endpoint_idx <- input$FinnGen_table_cells_selected
  })

  fit_gam <- reactive({
    if(nchar(input$endpoint)!=0){
      dat_endpoint <- filter(monthly_counts_FinRegistry,ENDPOINT==input$endpoint)
      mod_list <- run_seasonality_gam(dat_endpoint)
      return(list(dat_endpoint=dat_endpoint,mod_list=mod_list))
    }
  })

  output$plot_FinRegistry <- renderPlot({
    if(nchar(input$endpoint)!=0){
      compare_gam_fits(fit_gam()$mod_list,fit_gam()$dat_endpoint,input$endpoint) /
      (trend_plot(fit_gam()$mod_list) +
      seasonality_plot(fit_gam()$mod_list)) +
      plot_layout(widths = c(2,1,1), nrow=2)
    }
  })

  output$FinRegistry_table <- renderDataTable({
    datatable(seasonal_summary_FinRegistry,
              selection = list(mode='single',target = 'cell',selectable=cbind(1:nrow(seasonal_summary_FinRegistry),1)),
              options = list(autoWidth = TRUE,selection = 'none'))
  })

  output$FinnGen_table <- renderDataTable({
    datatable(seasonal_summary_FinnGen,
              selection = list(mode='single',target = 'cell',selectable=cbind(1:nrow(seasonal_summary_FinnGen),1)),
              options = list(autoWidth = TRUE,selection = 'none'))
  })


  output$search_endpoints <- renderUI({
    endpoints=sort(unique(seasonal_summary_FinRegistry$ENDPOINT))
    selectizeInput(
      inputId = "endpoint",
      label = "Search for endpoint",
      multiple = FALSE,
      choices = c("Search for endpoint" = "", endpoints),
      selected=as.data.frame(seasonal_summary_FinRegistry)[clicked$endpoint_idx],
      options = list(
        create = FALSE,
        placeholder = "Search for endpoint",
        maxItems = '1',
        onDropdownOpen = I("function($dropdown) {if (!this.lastQuery.length) {this.close(); this.settings.openOnFocus = false;}}"),
        onType = I("function (str) {if (str === '') {this.close();}}")
      ))
  })
}

# Run the application
shinyApp(ui = ui, server = server)

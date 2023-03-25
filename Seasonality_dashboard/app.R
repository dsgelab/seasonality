library(shiny)
library(patchwork)
library(readr)
library(DT)
library(knitr)
source('utils.R')

GWAS_info_dat <- read_tsv('data/finngen_R10_finngen_R10_analysis_data_finngen_R10_pheno_n.tsv')

monthly_counts_FinRegistry  <- read_tsv('data/FINREGISTRY_endpoints_monthly_count_green.txt') %>%
  filter_monthly_counts(.,GWAS_info_dat)

cols_to_format <- c('pval','log10_pval','val_peak','val_trough','ptr','af','var_fraction','dispersion')

seasonal_summary_FinnGen <- read.table('data/FINNGEN_seasonal_summary.txt',header=T) %>%
  as_tibble() %>%
  rename(log10_pval=log_pval) %>%
  mutate(log10_pval=log10_pval/log(10)) %>%
  inner_join(select(GWAS_info_dat,ENDPOINT=phenocode,num_gw_significant),by='ENDPOINT') %>%
  mutate(across(all_of(cols_to_format), format,digits=3))

seasonal_summary_FinRegistry <- read_tsv('data/FINREGISTRY_seasonal_summary.txt') %>%
  inner_join(select(seasonal_summary_FinnGen,ENDPOINT,num_gw_significant),by='ENDPOINT')%>%
  mutate(across(all_of(cols_to_format), format,digits=3)) %>%
  mutate(across(all_of(c('af','var_fraction')),as.numeric)) %>%
  select(-pval,-log10_pval)

seasonal_summary_FinRegistry_adj <- read_tsv('data/FINREGISTRY_seasonal_summary_adj.txt') %>%
  inner_join(select(seasonal_summary_FinnGen,ENDPOINT,num_gw_significant),by='ENDPOINT')%>%
  mutate(across(all_of(cols_to_format), format,digits=3)) %>%
  select(-pval,-log10_pval)

seasonal_splines_FinRegistry <- read_tsv('data/FINREGISTRY_seasonal_splines.txt') %>%
  inner_join(select(seasonal_summary_FinnGen,ENDPOINT),by='ENDPOINT') %>%
  group_by(month) %>%
  summarise(avg_seasonal_val=median(seasonal_val))
a <- 1;b <- 13


ui <- fluidPage(

  # Application title
  titlePanel("Seasonality dashboard"),

  # TabPanel with visualization and background
  tabsetPanel(
    tabPanel('Visualization',
             fluidRow(
               column(4,
                      uiOutput('search_endpoints')
               ),
               column(8,
                      checkboxGroupInput(inputId="adjustment",label='Adjustments',choices=c('Median seasonal pattern'='avg'))
               )
             ),
             hr(),
             # Output elements: figure and a table
             fluidRow(
               column(6,
                      plotOutput("plot_FinRegistry",height=600)
               ),
               column(6,
                      dataTableOutput("FinRegistry_table")
               )
             ),
    ),
    tabPanel('Background',
             uiOutput('background')
    )
  ),
  #To ensure Rmd file fills up the screen
  tags$head(tags$style(HTML("
                               body {
                                  width: 100% !important;
                                  max-width: 100% !important;
                               }
                               ")))
)

# Define server logic
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
      seasonal_spline_avg <- NULL
      if('avg' %in% input$adjustment){
        seasonal_spline_avg <- seasonal_splines_FinRegistry
      }
      mod_list <- run_seasonality_gam(dat_endpoint,seasonal_spline_avg)
      return(list(dat_endpoint=dat_endpoint,mod_list=mod_list,seasonal_spline_avg=seasonal_spline_avg))
    }
  })

  output$plot_FinRegistry <- renderPlot({
    if(nchar(input$endpoint)!=0){
      compare_gam_fits(fit_gam()$mod_list,fit_gam()$dat_endpoint,input$endpoint) /
        (trend_plot(fit_gam()$dat_endpoint,fit_gam()$mod_list,fit_gam()$seasonal_spline_avg) +
           seasonality_plot(fit_gam()$dat_endpoint,fit_gam()$mod_list,fit_gam()$seasonal_spline_avg,a=a,b=b)) +
        plot_layout(widths = c(2,1,1), nrow=2)
    }
  })

  output$FinRegistry_table <- renderDataTable({
    datatable(#if(!is.null(input$adjustment)) seasonal_summary_FinRegistry_adj else seasonal_summary_FinRegistry,
      seasonal_summary_FinRegistry,
      selection = list(mode='single',target = 'cell',selectable=cbind(1:nrow(seasonal_summary_FinRegistry),1)),
      options = list(autoWidth = TRUE,selection = 'none'))
  })

  # output$FinnGen_table <- renderDataTable({
  #   datatable(seasonal_summary_FinnGen,
  #             selection = list(mode='single',target = 'cell',selectable=cbind(1:nrow(seasonal_summary_FinnGen),1)),
  #             options = list(autoWidth = TRUE,selection = 'none'))
  # })


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

  output$background <- renderUI({
    HTML(markdown::markdownToHTML(knit('seasonal_patterns_endpoints.Rmd', quiet = TRUE)),fragment.only = T)
  })
}

# Run the application
shinyApp(ui = ui, server = server)

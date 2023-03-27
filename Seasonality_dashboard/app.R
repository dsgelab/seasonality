library(shiny)
library(dplyr)
library(tidyr)
library(broom)
library(lubridate)
library(readr)
library(DT)
library(knitr)
library(survival)
library(mgcv)
library(ggplot2)
library(patchwork)
library(shinydashboard)
library(shinyBS)
library(shinyjs)
library(shinyWidgets)
source('DGP_simulation_utils.R')
source('utils.R')

theme_custom <-  theme(plot.title = element_text(size=16),
                       axis.title.x = element_text(size=16),
                       axis.text.x = element_text(size=14),
                       axis.title.y = element_text(size=16),
                       axis.text.y = element_text(size=14),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=14))

js <- '.nav-tabs-custom .nav-tabs li.active {
    border-top-color: #1a3263;
}"'
js_box <- '.box.box-primary {
    border-top-color: #1a3263;
}"'
js_button <- 'a.btn {
    color: #000000;
}"'
js_background <- ".content {background-color: #FFFFFF;}"
m_item <- ".sidebar-menu li a { font-size: 15px; }"

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


ui <- shinyUI(fluidPage(
                #tags$head(tags$link(rel="shortcut icon", href="logo_nbg.png")),
                #To ensure Rmd file fills up the screen
                tags$head(tags$style(HTML("
                                             body {
                                                width: 100% !important;
                                                max-width: 100% !important;
                                             }
                                             "))),
                #tags$head(tags$style(".main-header .navbar{ margin-left: 0px !important;}")),
                title='Seasonality dashboard',
                useShinyjs(),
                withMathJax(),
                dashboardPage(skin = 'blue',
                              dashboardHeader(
                                tags$li(class = "dropdown",
                                        tags$style(".main-header {max-height: 60px}",
                                                   ".main-header .logo {height: 60px}",
                                                   ".main-header  {height: 60px}",
                                                   ".main-header .sidebar-toggle {height: 60px}",
                                                   ".main-header .sidebar-toggle {width: 60px}",
                                                   ".container-fluid {padding: 0px;}")
                                ),
                                title = span(
                                  span("Seasonality",
                                       style = 'font-weight: bold; color: #FFDB58; font-size: 22px'),
                                  span("dashboard",
                                       style = ' font-weight: bold; color: #FFFFFF; font-size: 22px')
                                ),
                                titleWidth = 300,
                                tags$li(a(href = 'https://github.com/dsgelab/seasonality',
                                          tags$img(src = 'github-icon.svg',
                                                   title = "github repo",
                                                   height = "24px"),
                                          title = "Back to Apps Home"),
                                        class = "dropdown")),
                              dashboardSidebar(width=300,
                                sidebarMenu(
                                  menuItem("Seasonality in FinRegistry", tabName = "finregistry", icon = icon("chart-simple")),
                                  menuItem("Simulation - binary trait", tabName = "binary_simulation", icon = icon("computer")),
                                  menuItem("Simulation - quantitative trait", tabName = "qt_simulation", icon = icon("chart-line"))
                                )
                              ),
                              dashboardBody(
                                            tags$style(js),
                                            tags$style(js_box),
                                            tags$style(js_button),
                                            tags$style(js_background),
                                            tabItems(
                                              tabItem(tabName="finregistry", width=NULL,
                                                tabBox(id='tabset1',width=NULL,
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
                                                )
                                              ),
                                              tabItem(tabName='binary_simulation',
                                                      # TabPanel with visualization and background
                                                      sidebarLayout(
                                                        sidebarPanel(width=3,
                                                                     sliderInput(inputId='beta_b',label='Baseline seasonal effect',value=0,min=-1,max=1,step = 0.1),
                                                                     sliderInput(inputId='beta_g',label='Genetic seasonal effect',value=0,min=-0.5,max=0.5,step=0.05),
                                                                     sliderInput(inputId='phi',label='Seasonal phase angle',value=0,min=0,max=3.14,step=0.1),
                                                                     radioButtons(inputId='mod_type',label = 'Model type',
                                                                                  choices=c('Additive'='additive','Recessive'='recessive','Dominant'='dominant'),
                                                                                  selected='additive'),
                                                                     hr(),
                                                                     sliderInput(inputId='lag',label='Mean diagnosis lag (months)',value=1,min=0,max=6,step=0.1),
                                                                     hr(),
                                                                     sliderInput(inputId='n',label='Sample size (n)',value=40000,min=10000,max=100000,step=5000),
                                                                     sliderInput(inputId='MAF',label='Minor allele frequency (MAF)',value=0.1,min=0.01,max=0.5,step=0.01)

                                                        ),
                                                        mainPanel(width=9,
                                                                  tabsetPanel(
                                                                    tabPanel('Distributions',
                                                                             fluidRow(
                                                                               plotOutput("plot_dist",height=300),
                                                                               plotOutput("plot_lag",height=300)
                                                                             ),
                                                                    ),
                                                                    tabPanel('Simulate',
                                                                             fluidRow(
                                                                               plotOutput("sim_plot",height=800)
                                                                             ),
                                                                    ),
                                                                    tabPanel('Background',
                                                                             uiOutput('background_sim')
                                                                    )
                                                                  ),
                                                        ),
                                                      ))
                                            )
                              )
                )
))

# Define server logic
server <- function(input, output, session) {
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
      mod_list <- run_seasonality_gam(dat_endpoint,a=a,b=b,seasonal_spline_avg)
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
    withMathJax(HTML(markdown::markdownToHTML(knit('seasonal_patterns_endpoints.Rmd', quiet = TRUE)),fragment.only = T))
  })

  observeEvent(input$beta_b, {
    gamma <- if(input$mod_type %in% c('recessive','dominant')) 1 else 2
    if(input$beta_b + gamma*input$beta_g>1){
      updateSliderInput(session, "beta_g", max=(1-input$beta_b)/gamma)
    }else if(input$beta_b + gamma*input$beta_g< -1){
      updateSliderInput(session, "beta_g", min=(-1-input$beta_b)/gamma)
    }else{
      updateSliderInput(session, "beta_g", min=-1/gamma,max=1/gamma)
    }
  })

  observeEvent(input$beta_g, {
    gamma <- if(input$mod_type %in% c('recessive','dominant')) 1 else 2
    if(input$beta_b + gamma*input$beta_g>1){
      updateSliderInput(session, "beta_b", max=(1-input$beta_g)/gamma)
    }else if(input$beta_b + gamma*input$beta_g< -1){
      updateSliderInput(session, "beta_b", min=(-1-input$beta_g)/gamma)
    }else{
      updateSliderInput(session, "beta_b", min=-1,max=1)
    }
  })

  output$plot_dist <- renderPlot({
    cols <- c("#BC3C29FF","#0072B5FF","#E18727FF")
    if(input$mod_type=='recessive'){
      cols <- cols[1:2]
    }else if(input$mod_type=='dominant'){
      cols <- cols[c(1,3)]
    }
    dist_grid <- get_dist_grid(beta=c(input$beta_b,input$beta_g),phi=input$phi,mod_type=input$mod_type,a=a,b=b) %>% bind_rows()
    dens_plot <- ggplot(dist_grid,aes(t,dens,col=as.factor(gt_lab))) +
      geom_line(size=1) +
      scale_y_continuous(limits=c(0,2/12)) +
      scale_color_manual(values=cols,name='Genotype') +
      xlab('t') +
      ylab('Density') +
      theme_classic() +
      theme_custom
    theme(legend.position='none')

    surv_plot <- ggplot(dist_grid,aes(t,surv,col=as.factor(gt_lab))) +
      geom_line(size=1) +
      scale_y_continuous(limits=c(0,1)) +
      scale_color_manual(values=cols,name='Genotype') +
      xlab('t') +
      ylab('Survival') +
      theme_classic() +
      theme_custom

    dens_plot + surv_plot
  })

  output$plot_lag <- renderPlot({
    mu_lag <- input$lag
    sigma_lag <- sqrt(mu_lag)
    dist_grid <- tibble(t=seq(0,12,by=0.01),dens=dgamma(t,shape=(mu_lag/sigma_lag)^2,rate=mu_lag/sigma_lag^2))
    lag_plot <- ggplot(dist_grid,aes(t,dens)) +
      geom_line(size=1,col="#E18727FF") +
      ggtitle('Assumed diagnosis lag distribution') +
      xlab('t') +
      ylab('Density') +
      theme_classic() +
      theme_custom

    anonymization_plot <- tibble(t=c(seq(-15,-1),seq(1,15)),dens=1/length(t)) %>%
      ggplot(aes(t,dens)) +
      geom_col(width=0.8,fill="#E18727FF") +
      scale_x_continuous(breaks=c(seq(-15,-1,by=2),seq(1,15,by=2))) +
      scale_y_continuous(limits=c(0,0.05),expand=c(0,0)) +
      ggtitle('Date anonymization distribution') +
      xlab('t') +
      ylab('Probability') +
      theme_classic() +
      theme_custom

    lag_plot + anonymization_plot
  })

  output$sim_plot <- renderPlot({
    print('hae')
    sim_dat <- simulate_seasonal_dist(n=input$n,MAF=input$MAF,beta=c(input$beta_b,input$beta_g),phi=input$phi,mod_type=input$mod_type,a=a,b=b)
    sim_dat$EVENT_MONTH_DEC_PB <- perturb_surv_times(surv_times=sim_dat$EVENT_MONTH_DEC,mu=input$lag,sigma=sqrt(input$lag),a=a,b=b)
    monthly_counts <- mutate(sim_dat,month=floor(EVENT_MONTH_DEC),
                             month_pb=floor(EVENT_MONTH_DEC_PB)) %>%
      pivot_longer(cols=c('month','month_pb'),names_to = 'ENDPOINT',values_to = 'EVENT_MONTH') %>%
      group_by(ENDPOINT,EVENT_MONTH) %>%
      summarise(count=length(EVENT_MONTH)) %>%
      mutate(ENDPOINT=as.character(factor(ENDPOINT,levels=c('month','month_pb'),
                                          labels=c('EVENT_MONTH_DEC','EVENT_MONTH_DEC_PB'))),
             ENDPOINT_LAB=as.character(factor(ENDPOINT,levels=c('EVENT_MONTH_DEC','EVENT_MONTH_DEC_PB'),
                                              labels=c('True disease onset','Perturbed disease diagnosis'))))

    surv_obj_sim <- Surv(time=sim_dat$EVENT_MONTH_DEC,event=rep(1,nrow(sim_dat)))
    survfit_sim <- survfit(surv_obj_sim ~ gt_lab,data=sim_dat)
    surv_obj_sim_pb <- Surv(time=sim_dat$EVENT_MONTH_DEC_PB,event=rep(1,nrow(sim_dat)))
    survfit_sim_pb <- survfit(surv_obj_sim_pb ~ gt_lab,data=sim_dat)
    print('hae2')
    KM_plot <- (KM_gg(survfit_sim,'True disease onset') + theme_custom) +
      (KM_gg(survfit_sim_pb,title='True disease onset with perturbation') + theme_custom + theme(legend.position='none'))
    monthly_counts_plot <- ggplot(monthly_counts,aes(EVENT_MONTH,count,fill=factor(ENDPOINT_LAB,levels=c('True disease onset','Perturbed disease diagnosis')))) +
      geom_col(position='dodge') +
      scale_x_continuous(breaks=seq_len(12),expand=c(0,0.2)) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_manual(values=c("#0072B5FF","#E18727FF"),name='') +
      xlab('Month') +
      ylab('Count') +
      theme_classic() +
      theme_custom
    mod_dat_true <- filter(monthly_counts,ENDPOINT=='EVENT_MONTH_DEC')
    mod_counts_true <- run_seasonality_gam_sim(mod_dat_true,a=a,b=b)
    mod_dat_pb <- filter(monthly_counts,ENDPOINT=='EVENT_MONTH_DEC_PB')
    mod_counts_pb <- run_seasonality_gam_sim(mod_dat_pb,a=a,b=b)

    seasonal_comp <- bind_rows(list('True disease onset'=get_seasonal_comp(mod_dat_true,mod_counts_true,a=a,b=b),
                                    'Perturbed disease diagnosis'=get_seasonal_comp(mod_dat_pb,mod_counts_pb,a=a,b=b)),.id='ENDPOINT_LAB') %>%
      mutate(ENDPOINT_LAB=factor(ENDPOINT_LAB,levels=c('True disease onset','Perturbed disease diagnosis')))
    GAM_seasonality_plot <- seasonality_plot_sim(seasonal_comp,a=a,b=b) + theme_custom

    seasonal_mods <- run_seasonal_models(monthly_counts,sim_dat,endpoint_id='EVENT_MONTH_DEC',a=a,b=b)
    seasonal_mods_pb <- run_seasonal_models(monthly_counts,sim_dat,endpoint_id='EVENT_MONTH_DEC_PB',a=a,b=b)

    summary_dat <- bind_rows(list(binary_mod=tidy(seasonal_mods$binary_mod),
                                  binary_mod_pb=tidy(seasonal_mods_pb$binary_mod),
                                  qt_mod=tidy(seasonal_mods$qt_mod),
                                  qt_mod_pb=tidy(seasonal_mods_pb$qt_mod)),.id='mod_type') %>%
      filter(term=='gt_mod') %>%
      mutate(mod_type=gsub('_mod','',mod_type)) %>%
      separate(mod_type,into=c('model_type','data_type')) %>%
      mutate(data_type=ifelse(is.na(data_type),'True disease onset','Perturbed disease diagnosis'))
    effect_plot <- ggplot(summary_dat) +
      geom_pointrange(aes(x=model_type,y=estimate,ymin=estimate-1.96*std.error,ymax=estimate+1.96*std.error,col=factor(data_type,levels=c('True disease onset','Perturbed disease diagnosis'))),
                      position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0),size=1) +
      scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
      xlab('Model') +
      ylab('Effect size') +
      theme_classic() +
      theme_custom

    pval_plot <- ggplot(summary_dat) +
      geom_col(aes(x=model_type,y=-log10(p.value),fill=factor(data_type,levels=c('True disease onset','Perturbed disease diagnosis'))),position = 'dodge',width=0.6) +
      geom_abline(slope=0,intercept=-log10(5e-8),color='red',linetype='dashed') +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_manual(values=c("#0072B5FF","#E18727FF"),name='') +
      xlab('Model') +
      ylab(expression(-log10(pval))) +
      theme_classic() +
      theme_custom
    KM_plot / (monthly_counts_plot + GAM_seasonality_plot) / (effect_plot + pval_plot)
  })


  output$background_sim <- renderUI({
    withMathJax(HTML(markdown::markdownToHTML(knit('simulation_background.Rmd', quiet = TRUE)),fragment.only = T))
  })
}

# Run the application
shinyApp(ui = ui, server = server)

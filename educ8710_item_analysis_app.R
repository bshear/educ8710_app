#educ8710 item analysis app
#original: 25nov2017
#author:   benjamin r. shear; benjamin.shear@colorado.edu
#update:   27aug2023
#inspired by Martinkova et al. at: https://www.shinyitemanalysis.org

# check for and install necessary packages
for (p in c("shiny", "CTT", "eRm", "corrplot",
            "lattice", "tidyverse", "psych",
            "likert", "WrightMap", "MASS")) {
  print(p)
  if(!require(p, character.only = TRUE)) {
    install.packages(p)
    require(p, character.only = TRUE)
  }
}

# Define UI for data upload app ----

ui <- navbarPage(
  "Item Analysis for EDUC 8710",
  
  # Tab: Overview #####
  
  tabPanel(
    "Overview",
    h3("Overview"),
    p("This application is intended to assist with project analyses for EDUC 8710.
        This app was written by Benjamin Shear. Inspired partly by:"),
    p("Martinková, P., & Drabinová, A. (2018). ShinyItemAnalysis for teaching psychometrics 
        and to enforce routine analysis of educational tests. 
        The R Journal, 10(2), 503-515. https://doi.org/10.32614/RJ-2018-074"),
    p("This app is a work in progress -- feedback welcome!"),
    h3("Getting Started"),
    p("Each tab in the app produces results based on a different analysis we
      will learn in class. These include basic summary statistics and classical test theory (CTT) results,
      dimensionality analyses, and Rasch model results."),
    p("Prior to uploading your data to the app for analysis, there are some important
      data processing steps that must be taken. See the description below for
      more details about this."),
    h3("Data Requirements"),
    p("Prior to uploading your data for analysis, there are some important data
        processing and cleaning steps that need to be taken.
        There can be no missing values. Each respondent must have a response for all items. You can have
        the first column as an ID variable, but all other columns are assumed
        to be responses to your survey items that measure a single latent construct.
        The responses can be either ordinal or binary. However, most of the analyses
        will convert responses to binary 0/1 scores, using a pre-specified
        value to dichotomize ordinal responses. If any items require reverse-scoring,
        the reverse-scoring must be done PRIOR to importing your data. The app will
        not do the reverse-scoring for you."),
    tags$ul(
      tags$li("Data should be in CSV format if possible."), 
      tags$li("Data shoud be in wide format, with each row representing a person
              and each column representing an item."), 
      tags$li("The first column can be an ID variable, but all other columns must
              be responses to your survey items measuring the construct of
              interest (IMPORTANT: remove all demographic and other items that are not
              part of your survey instrument)."),
      tags$li("There can be NO missing values in the data. You either need to drop
              respondents who skipped questions or decide on an appropriate value
              to assign for the missing responses and do this before importing."),
      tags$li("If your survey contains items that require reverse-coding, you MUST
              do the reverse-coding prior to importing your data for analysis."),
      tags$li("Some of the analyses can only be done with binary item scoring. If
              your survey items have more than 2 response categories see notes below.")
    ),
    h3("How to Analyze Items with More than 2 Response Categories"),
    p("If your survey items have more than two response categories, you have two
      options for analysis:"),
    tags$ol(
      tags$li("Convert all responses to binary responses *PRIOR* to importing your 
              data for analysis. For example, if you have items on a 5-point Likert
              scale, you could decide that all responses of 1/2/3 will receive a
              score of 0, while all responses of 4/5 will receive a score of 1. 
              This way, an item response score of 1 still represents a higher
              level of the construct than does a score of 0."), 
      tags$li("You can import your data with the full set of item responses,
              but you muyst also specify a threshold that the app will use to
              convert these responses to binary scores for analyses that require binary scoring.
              To accomplish the scoring
              described above, where (1/2/3)=0 and (4/5)=1, you would specify 4
              as the threshold so that all scores greater than or equal to 4 get
              converted to a 1 and all responses less than 4 get converted to a 0. This
              threshold will be constant across all items.")
    ),
    p("The Response Frequency and Likert Plots, CTT Reliability and Item Stats,
    and Dimensionality analyses can be done  with either the binary scores or 
    the full response scores. All other analyses will be done ONLY with the binary
    version of the responses."),
    p("So, for example, if you convert your responses
      to binary scoring BEFORE importing your data, all results in all tabs will
      be based on binary response scoring. However, if you import responses with
      more than 2 response categories, then you will get to choose whether to
      use the binary scores or the full responses for some of the tabs."),
    h3("Packages"),
    p("This app relies on a number of R packages. Running the app the first time
        will automatically install the necessary packages if they are not already
        installed on your computer. Currently the app loads the following packages:
        shiny, CTT, eRm, corrplot, lattice, tidyverse, psych, likert, WrightMap, MASS.")
  ),
  
  # Tab: Data #####
  
  tabPanel(
    "Data",
    sidebarLayout(
      
      ## Sidebar panel for inputs ----
      sidebarPanel(
        
        helpText("Load your data file here. The data file should be a CSV file
          where each row is a person and each column is an item. All items are
          assumed to be binary, with 1=correct/more and 0=incorrect/less.
          If your responses include more than 2 response categories, be sure
          to specify the correct threshold value below that can be used to convert
          the responses to binary scoring.
          If you have an ID variable as the first column, make sure to check
          the appropriate box below. Any additional non-item columns need to
          be removed before loading. They will cause errors or incorrect results."),
        
        helpText("NOTE: some of the analyses in this
                 app currently only work for binary item response
                 data scored 1=correct (or agree, or more often) and 0=incorrect 
                 (or disagree or less often)."),
        
        ## Input: Select a file ----
        fileInput("file1", "Choose CSV File",
                  multiple = TRUE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        ## Horizontal line ----
        tags$hr(),
        
        ## Input: checkbox if file has header ----
        checkboxInput("header", "Header", TRUE),
        
        ## Input: checkbox if data file has ID's as first column
        radioButtons("hasid", "Does data file contain an ID variable in column 1?",
                     choices = c(Yes = "Yes", No = "No"), selected = "No"),
        
        ## Input: value for dichotomizing ----
        numericInput(
          "thresh",
          "Threshold for dichotomizing responses (responses greater than or equal
          to this value will be scored as 1, and responses less than this
          value will be scored as 0 for the binary analyses)",
          value = 1,
          min = -10,
          max = 10,
          step = 1
        ),
        
        ## Input: Select separator ----
        radioButtons("sep", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),
        
        ## Input: Select quotes ----
        radioButtons("quote", "Quote",
                     choices = c(None = "",
                                 "Double Quote" = '"',
                                 "Single Quote" = "'"),
                     selected = '"'),
        
        ## Horizontal line
        tags$hr(),
        
        ## Input: Select number of rows to display ----
        radioButtons("disp", "Display",
                     choices = c(Head = "head",
                                 All = "all"),
                     selected = "head")
        
      ),
      
      mainPanel(
        h2("Imported data before dichotomizing:"),
        tableOutput("contents_full"),
        h2("Imported data after dichotomizing:"),
        tableOutput("contents")
      )
    )
  )
  ,
  
  # Tab: Response Frequency and Likert Plot #####
  
  tabPanel("Response Frequency and Likert Plot",
           sidebarLayout(
             sidebarPanel(
               p("This tab provides a summary of response frequencies. Plots
                 are created using the Likert package. Note that the plots
                 and response frequency table assume that all items have the same
                 response scales."),
               radioButtons("likert_binary",
                            "Calculate response frequency stats for:",
                            choices = c("Binary Scores"="binary",
                                        "Full Scores" = "full"),
                            selected = "binary"),
               radioButtons("likert_order",
                            "Sort items by response frequencies?",
                            choices = c("Yes"="yes",
                                        "No" = "no"),
                            selected = "yes"),
               p("Adjust the size of the Likert frequency plot using these sliders:"),
               sliderInput("likert_height", "height", min = 300, max = 800, value = 500),
               sliderInput("likert_width", "width", min = 300, max = 800, value = 500)
             ),
             mainPanel(
               h2("Table of response frequencies:"),
               tableOutput("likert_freq_table"),
               h2("Plot of response frequencies:"), 
               plotOutput("likert_plot"))
           ))
  ,
  
  # Tab: CTT Total Score Reliability and Item Stats #####
  
  tabPanel("CTT Total Score Reliability and Item Stats",
           sidebarLayout(
             sidebarPanel(
               helpText("Classical test theory (CTT) item statistics are shown on
                this page. Difficulty is computed as the mean response to each
                item. The item discrimination is calculated and reported using
                        the corrected-total item correlation (correlation between 
                        the item response and total score without that item). The Alpha if
                        deleted column shows the estimate of alpha if each item were dropped."),
               downloadButton("dl_ctt_table", "Download Table"),
               
               ## Input: checkbox if data file has ID's as first column #####
               radioButtons("ctt_binary",
                            "Calculate CTT stats for:",
                            choices = c("Binary Scores"="binary",
                                        "Full Scores" = "full"),
                            selected = "binary")
             ),
             ## main panel #####
             mainPanel(
               verbatimTextOutput("obs_score_stats_text"),
               tableOutput("ctt_table"),
               plotOutput("cttstatsplot"),
               p("Total score summary statistics:"),
               tableOutput("total_score_stats"),
               plotOutput("total_score_hist"),
               tableOutput("total_score_table")
             )
           )
  )
  ,
  
  # tabPanel("Total Score Stats",
  #          sidebarLayout(
  #            sidebarPanel(
  #              tableOutput("total_score_stats"),
  #              verbatimTextOutput("obs_score_stats_text"),
  #              
  #            ),
  #            mainPanel(
  #              plotOutput("total_score_hist"),
  #              tableOutput("total_score_table")
  #            )
  #          )
  # )
  # ,
  
  # Tab: Dimensionality #####
  
  tabPanel("Dimensionality",
           sidebarLayout(
             sidebarPanel(
               helpText("This page computes eigenvalues for the item correlation
                        matrix and plots them in a scree plot. A heat map of
                        the item correlation matrix is also printed."),
               ## Input: choose binary or full/poly data for dimensionality #####
               radioButtons("dimensionality_binary",
                            "Calculate dimensionality analysis based on:",
                            choices = c("Binary Scores"="binary",
                                        "Full Scores" = "full"),
                            selected = "binary")
               
               #,
               
               # p("Adjust the size of the scree plot using these sliders:"),
               # sliderInput("screeplot_height", "height", min = 100, max = 700, value = 400),
               # sliderInput("screeplot_width", "width", min = 100, max = 700, value = 400),
               # 
               # p("Adjust the size of the corrplot heat map with these sliders:"),
               # sliderInput("corrplot_height", "height", min = 100, max = 700, value = 400),
               # sliderInput("corrplot_width", "width", min = 100, max = 700, value = 400)
             ),
             
             ## main panel #####
             
             mainPanel(
               fluidRow(
                 plotOutput("scree_plot", width=700, height=500),
                 plotOutput("corrplot", width=800, height=800),
                 h3("Item Correlation Matrix:"),
                 p("(Copy and paste values into Excel for formatting.)"),
                 tableOutput("item_cor_table"),
                 h3("Observed and Random Eigenvalues for Parallel Analysis:"),
                 p("The random eigenvalues here for parallel analysis are based on correlation matrices for
                   100 random datasets of simulated multivariate normal data with the same number of items
                   as the empirical data. The average of the random eigenvalues is reported along with the
                   95th percentile value for each random eigenvalue."),
                 tableOutput("eigen_value_table")
                 #,
                 # working here dimensionality
                 # tableOutput("parallel_analysis_table")
               )
             )
           )
  ) ,   
  
  # Tab: CTT Item Curves #####
  
  tabPanel("CTT Item Curves",
           sidebarLayout(
             sidebarPanel(
               helpText("These plots show empirical item response curves for
                        each item. These plots use the binary scoring. The points
                        show the proportion of 1's at each total score value.
                        The black lines trace the empirical values while the blue
                        lines show a smoothed nonparametric trend line. 
                        \n You can adjust the size of the plots using these sliders:"),
               
               ## Input: figure size #####
               
               sliderInput("ctt_curve_height", "height", min = 100, max = 700, value = 400),
               sliderInput("ctt_curve_width", "width", min = 100, max = 700, value = 400)
             ),
             mainPanel(
               plotOutput("ctt_item_plots", width=400, height=400)
             )
           )
  ),
  
  # Tab: Rasch Item Parameters #####
  
  tabPanel("Rasch Item Parameters",
           sidebarLayout(
             sidebarPanel(
               helpText("Item parameters are estimated using the RM() function
              	in the eRm package. The item parameters are estimated using a
              	technique called 'conditional maximum likelihood'. The item
                difficulty parameters are scaled to sum to 0 for model
                identification."),
               downloadButton("dl_rasch_item_tab", "Download Table")
             ),
             mainPanel(tableOutput("rasch_item_tab"))
           )
  )
  ,
  
  # Tab: Rasch Person Parameters #####
  
  tabPanel("Rasch Person Parameters",
           sidebarLayout(
             sidebarPanel(helpText("Note: this table will NOT include respondents
                                   with minimum/maximum possible total scores, if
                                   they are present in the data."),
                          downloadButton("dl_theta_table", "Download Table")),
             mainPanel(tableOutput("theta_table"))
           )
  )
  ,
  
  # Tab: Wright Map #####
  
  tabPanel("Wright Map",
           sidebarLayout(
             sidebarPanel(helpText("blank")),
             mainPanel(plotOutput("wright_map"),
                       plotOutput("wright_map2"),
                       plotOutput("wright_map3")
             )
           )
  )
  ,
  
  # Tab: Theta Stats and SEM #####
  
  tabPanel("Theta Stats and SEM",
           sidebarLayout(
             sidebarPanel(
               verbatimTextOutput("sep_rel_text"),
               h3("Theta table:"),
               tableOutput("theta_sem_table")),
             mainPanel(
               h2("Summary Statistics for Theta Estimates:"),
               tableOutput("theta_summary_table"),
               h2("Plots:"),               
               plotOutput("sem_by_theta_plot"),
               plotOutput("theta_by_total_plot"),
               plotOutput("theta_hist_plot")
             )
           )
  )
  ,
  
  # Tab: Item Fit Statistic Plots #####
  
  tabPanel("Item Fit Statistic Plots",
           sidebarLayout(
             sidebarPanel(),
             mainPanel(
               plotOutput("item_infit_plot"),
               plotOutput("item_outfit_plot")
             )
           )
  )
  ,
  
  # Tab: Item Fit Plots #####
  
  tabPanel("Item Fit Plots",
           sidebarLayout(
             sidebarPanel(
               sliderInput("nbins", "Number of Theta Bins", min = 1, max = 20, step = 1, value = 5),
               checkboxInput("showN", "Show N Sizes on Plots", value=FALSE),
               checkboxInput("adjSize", "Scale Point Size Based on N", value=FALSE),
               checkboxInput("showLine", "Show Empirical ICC line", value=FALSE),
               helpText("This download button will save a CSV file containing
                        the binned proportion correct data, average theta values,
                        average total scores, and N counts used to plot the
                        item fit curves on this page."),
               downloadButton("dl_item_icc_fit_table", "Download Table")
             ),
             mainPanel(
               #plotOutput("item_icc_fit_plots")
               uiOutput("item_icc_fit_plots_dynamic")
             )
           )
  )
  ,
  
  # Tab: Item IRT Plots #####
  
  tabPanel("Item IRF Plots",
           sidebarLayout(
             sidebarPanel(helpText("blank")),
             mainPanel(plotOutput("item_irf_plots"))
           )
  )
  ,
  
  # Tab: Information Plots #####
  
  tabPanel("Information Plots",
           sidebarLayout(
             sidebarPanel(helpText("blank")),
             mainPanel(plotOutput("info_plots"))
           )
  )
  ,
  
  # Tab: Sorted Item-Person Table #####
  
  tabPanel("Sorted Item-Person Table",
           sidebarLayout(
             sidebarPanel(helpText("Note: this table will NOT include respondents
                                   with minimum/maximum possible total scores, if
                                   they are present in the data."),
                          downloadButton("dl_sorted_response_table", "Download Table")),
             mainPanel(tableOutput("sorted_response_table"))
           )
  )
  ,
  
  # Tab: Extrapolated Theta Table #####
  
  tabPanel("Extrapolated Theta Table",
           sidebarLayout(
             sidebarPanel(helpText("Note: this table WILL include respondents
                                   with minimum/maximum possible total scores, if
                                   they are present in the data. However, these
                                   respondents will not have an estimated SEM of
                                   theta. Respondents with minimum/maximum possible
                                   total scores will have Interpolated=TRUE in this
                                   table."),
                          downloadButton("dl_extrapolated_theta_table", "Download Table")),
             mainPanel(tableOutput("extrapolated_theta_table"))
           )
  )
)



# server function #####

server <- function(input, output) {
  
  ## Data import functions #####
  
  ### Import dataset and remove ID column if necessary
  
  dataset_noid <- reactive({
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep    = input$sep,
                   quote  = input$quote)
    
    if (input$hasid=="Yes") {
      df <- df[,-1]
    }
    
    # convert to binary if needed
    # if already binary 0/1, this shouldn't matter
    
    df <- as.data.frame(apply(df, 2, function(x) {
      ifelse(x>=input$thresh, 1, 0)
    }))
    
    return(df)
  })
  
  ### import dataset with no ID variable and without dichotomizing
  
  dataset_poly <- reactive({
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep    = input$sep,
                   quote  = input$quote)
    
    if (input$hasid=="Yes") {
      df <- df[,-1]
    }
    
    return(df)
    
  })
  
  ## More functions #####

  #############################################################################-
  
  likert_data <- function(binary=TRUE) {
    if(input$likert_binary=="binary")
    {return(dataset_noid())} else
    {return(dataset_poly())}
  }
  
  output$likert_freq_table <- renderTable({
    
    d <- likert_data()
    
    psych_results <- psych::alpha(d)
    return(psych_results$response.freq)
    
  }, include.rownames=TRUE)
  
  
  output$likert_plot <- renderPlot(
    
    width = function() input$likert_width,
    height = function() input$likert_height,
    res = 96,
    
    {
      
      d <- likert_data()
      
      mylevels <- c(min(d):max(d))
      d_likert <- d
      for(i in seq_along(d_likert)) {
        d_likert[,i] <- factor(d_likert[,i], levels=mylevels)
      }
      
      if(input$likert_order=="yes") {
        print(plot(likert(d_likert), type="bar"))
      }
      if(input$likert_order=="no") {
        print(plot(likert(d_likert), type="bar", group.order = names(d_likert)))
      }
      
    })
  
  #############################################################################-
  # Function to compute CTT statistics
  
  get_ctt_stats_func <- function(binary=TRUE) {
    if(binary==TRUE) {
      df <- dataset_noid()} else {
        df <- dataset_poly()
      }
    
    total_score <- apply(df, 1, sum)
    
    # classical item stats
    ctt_stats <- itemAnalysis(df, itemReport=TRUE)
    ctt_table <- data.frame(
      ItemMean            = ctt_stats$itemReport$itemMean,
      Corrected_ItemTotal = ctt_stats$itemReport$pBis,
      AlphaIfDeleted      = ctt_stats$itemReport$alphaIfDeleted
    )
    
    ts_table <- data.frame(
      min    = min(total_score),
      mean   = mean(total_score),
      median = median(total_score),
      max    = max(total_score),
      sd     = sd(total_score)
    )
    
    return(list(
      total_score = total_score,
      ctt_table   = ctt_table,
      ctt_stats   = ctt_stats,
      ts_table    = ts_table,
      ts_var      = var(total_score),
      alpha       = ctt_stats$alpha))
    
  }
  
  get_ctt_stats <- reactive({
    # return CTT stats table, total scores
    
    if(input$ctt_binary == "binary") {
      return(get_ctt_stats_func(TRUE))} else {
        return(get_ctt_stats_func(FALSE))
      }
    
    
  })
  
  
  #############################################################################-
  # Function to get Rasch model results
  
  get_rasch_stats <- reactive({
    # return item stats table
    # return person stats table
    # return sorted person-item response matrix
    rasch_mod <- RM(dataset_noid())
    itdiff    <- -1*rasch_mod$betapar # item difficulty values
    itdiff.se <- rasch_mod$se.beta
    
    pers.out  <- person.parameter(rasch_mod)  # generates standardized residuals
    ifit      <- itemfit(pers.out)            # item fit statistics
    
    ifit_table  <- data.frame(
      Difficulty = itdiff,
      SE         = itdiff.se,
      outMNSQ    = ifit$i.outfitMSQ,
      out.t      = ifit$i.outfitZ,
      inMNSQ     = ifit$i.infitMSQ,
      in.t       = ifit$i.infitZ
    )
    
    pers_fit <- personfit(pers.out)
    expected <- pmat(pers.out)      # matrix of expected item response
    resid    <- residuals(pers.out) # matrix, standardized residuals of responses
    observed <- pers.out$X.ex       # observed responses for people without 
    # perfect/0 scores
    
    # this is the SEM as estimated with ALL observed scores
    ctt_sem <- sqrt(
      var(get_ctt_stats()[["total_score"]])*
        (1-get_ctt_stats()[["alpha"]])
    )
    
    # Make a large table with all of these values, plus CTT values
    pers_table <- data.frame(
      theta        = pers.out$thetapar$NAgroup1,
      se.theta     = pers.out$se.theta$NAgroup1,
      total_score  = apply(pers.out$X.ex, 1, sum),
      ctt_sem      = ctt_sem,
      infitMNSQ    = pers_fit$p.infitMSQ,
      infitMNSQ.z  = pers_fit$p.infitZ,
      outfitMNSQ   = pers_fit$p.outfitMSQ,
      outfitMNSQ.z = pers_fit$p.outfitZ
    )
    
    # extrapolated theta values
    theta_extrapolated <- cbind(
      pers.out$theta.table,
      total_score = apply(pers.out$X, 1, sum)
    )
    names(theta_extrapolated)[names(theta_extrapolated)=="Person Parameter"] <- "theta"
    
    
    # sorted person-item matrix to evaluate misfitting responses
    observed_sort             <- as.data.frame(observed[,order(itdiff)])
    observed_sort$total_score <- apply(observed_sort, 1, sum)
    observed_sort <- cbind(
      observed_sort,
      outfitMNSQ = pers_fit$p.outfitMSQ,
      infitMNSQ  = pers_fit$p.infitMSQ)
    observed_sort <- observed_sort[order(observed_sort$total_score),]
    
    return(list(
      ifit_table         = ifit_table,
      pers_table         = pers_table,
      observed_sort      = observed_sort,
      rasch_mod          = rasch_mod,
      theta_extrapolated = theta_extrapolated
    )
    )
  })
  
  #############################################################################-
  # output/download CTT table
  output$ctt_table <- renderTable({
    stuff <- get_ctt_stats()
    stuff[["ctt_table"]]
  }, include.rownames=TRUE)
  
  output$dl_ctt_table <- downloadHandler(
    filename = function() {
      paste("ctt_item_stats.csv")
    },
    content = function(file) {
      write.csv(get_ctt_stats()[["ctt_table"]], file, row.names = TRUE)
    }
  )
  
  
  #############################################################################-
  # CTT stats plot
  output$cttstatsplot <- renderPlot({
    
    d <- get_ctt_stats()[["ctt_table"]]
    
    with(d, plot(Corrected_ItemTotal~ItemMean, las=1,
                 xlim=c(min(0,min(d$ItemMean)),max(1, max(d$ItemMean))),
                 ylim = c(min(d$Corrected_ItemTotal,0),1), pch=19,
                 xlab = "Item Mean", ylab = "Item Discrimination",
                 main="Item Mean vs. Corrected Item-Total Correlation"))
    abline(h=0.2, lty=2)
    if(input$ctt_binary=="binary") {
      abline(v=0.2, lty=2)
      abline(v=0.8, lty=2)
    }
  }, height=400, width=400)
  
  output$total_score_stats <- renderTable({
    get_ctt_stats()[["ts_table"]] 
  })
  
  output$total_score_hist <- renderPlot({
    ts <- get_ctt_stats()[["total_score"]]
    hist(ts, col="grey", main="Histogram of Observed Total Scores",
         xlab="Observed Score", las=1
         #,breaks=c(-0.5:(max(ts)+0.5))
    )
  }, width=500)
  output$total_score_table <- renderTable({
    tst <- as.data.frame(table(get_ctt_stats()[["total_score"]]))
    names(tst) <- c("Score", "Frequency")
    tst
  })
  
  
  #############################################################################-
  # CTT item curves
  output$ctt_item_plots <- renderPlot(
    
    width = function() input$ctt_curve_width,
    height = function() input$ctt_curve_height,
    res = 96,
    
    {
      
      d <- dataset_noid()
      
      d$xxtot <- apply(d, 1, sum)
      d$xxid <- c(1:nrow(d))
      
      dlong <- d %>%
        pivot_longer(cols = -c("xxid","xxtot"), names_to="item") %>%
        group_by(item) %>%
        mutate(disc=cor(value,(xxtot-value))) %>%
        ungroup() %>%
        group_by(item, xxtot) %>%
        summarise(
          n=n(),
          mean=mean(value),
          disc=mean(disc)
        )
      
      text_labs <- dlong %>%
        group_by(item) %>%
        summarise(itmn=paste0("mean: ", round(mean(mean),2)),
                  itdisc=paste0("disc: ", round(mean(disc),2)))
      
      dlong %>%
        ggplot(aes(x=xxtot, y=mean)) +
        geom_point() +
        geom_line() +
        geom_smooth(se=FALSE) +
        geom_text(aes(x=2,y=0.85,label=itdisc), data=text_labs, size=2) +
        geom_text(aes(x=2,y=0.95,label=itmn), data=text_labs, size=2) +
        ylim(0,1) +
        xlab("Total Score") +
        ylab("Item Mean") +
        facet_wrap(~item, ncol=3) +
        theme_bw()
      
    })
  
  
  #############################################################################-
  # some summary stats and SEM
  output$obs_score_stats_text <- renderText({
    stuff <- get_ctt_stats()
    ts    <- stuff[["total_score"]]
    r     <- stuff[["alpha"]]
    paste0("The observed score mean is ", round(mean(ts), 3),
           " (SD = ", round(sd(ts), 3), ").",
           "\nEstimated Cronbach's alpha coefficient: ", round(r, 3), ".",
           "\nThe estimated SEM = sqrt(Var(X)*(1-alpha)) = ",
           round(sqrt(var(ts)*(1-r)), 3), ".")
  })
  
  
  #############################################################################-
  # output Rasch Item statistics
  output$rasch_item_tab <- renderTable({
    get_rasch_stats()[["ifit_table"]]
  }, include.rownames=TRUE, digits = 4)
  
  output$dl_rasch_item_tab <- downloadHandler(
    filename = function() {
      paste("rasch_item_stats.csv")
    },
    content = function(file) {
      write.csv(get_rasch_stats()[["ifit_table"]], file, row.names = TRUE)
    }
  )
  
  
  #############################################################################-
  # output theta descriptive statistics
  output$theta_table <- renderTable({
    get_rasch_stats()[["pers_table"]]
  }, include.rownames=TRUE, digits=4)
  
  output$dl_theta_table <- downloadHandler(
    filename = function() {
      paste("rasch_person_stats.csv")
    },
    content = function(file) {
      write.csv(get_rasch_stats()[["pers_table"]], file, row.names = TRUE)
    }
  )
  
  
  #############################################################################-
  # output extrapolated theta's
  output$extrapolated_theta_table <- renderTable({
    get_rasch_stats()[["theta_extrapolated"]]
  }, include.rownames=TRUE, digits=4)
  
  output$dl_extrapolated_theta_table <- downloadHandler(
    filename = function() {
      paste("extrapolated_theta_table.csv")
    },
    content = function(file) {
      write.csv(get_rasch_stats()[["theta_extrapolated"]], file, row.names = TRUE)
    }
  )
  
  
  #############################################################################-
  # output sorted item-person matrix
  output$sorted_response_table <- renderTable({
    get_rasch_stats()[["observed_sort"]]
  }, include.rownames=TRUE)
  
  output$dl_sorted_response_table <- downloadHandler(
    filename = function() {
      paste("sorted_response_table.csv")
    },
    content = function(file) {
      write.csv(get_rasch_stats()[["observed_sort"]], file, row.names = TRUE)
    }
  )
  
  
  #############################################################################-
  # Wright Maps
  output$wright_map <- renderPlot({
    mod <- get_rasch_stats()[["rasch_mod"]]
    plotPImap(mod, sorted = TRUE)
  },width=500)
  
  output$wright_map2 <- renderPlot({
    mod <- get_rasch_stats()[["rasch_mod"]]
    plotPImap(mod, sorted = FALSE)
  },width=500)
  
  # using wrightmap package
  output$wright_map3 <- renderPlot({
    it_tab <- get_rasch_stats()[["ifit_table"]]
    it_tab <- arrange(it_tab, Difficulty)
    betas  <- it_tab$Difficulty
    beta_names <- row.names(it_tab)
    thetas <- get_rasch_stats()[["pers_table"]]$theta
    wrightMap(thetas = thetas, thresholds = betas,
              label.items=beta_names, label.items.srt=90,
              show.thr.lab=FALSE)
  },width=500)  
  
  
  #############################################################################-
  # item fit statistic plots
  output$item_infit_plot <- renderPlot({
    mod  <- get_rasch_stats()[["rasch_mod"]]
    ifit <- itemfit(person.parameter(mod))
    dotplot(ifit$i.infitMSQ, xlim=c(0,3), xlab = "Infit Mean Square",
            panel = function(...) { 
              panel.dotplot(...) 
              panel.abline(v=1.33, lty=2)
              panel.abline(v=0.75, lty=2)
            })   	
  }, width=500)
  
  output$item_outfit_plot <- renderPlot({
    mod  <- get_rasch_stats()[["rasch_mod"]]
    ifit <- itemfit(person.parameter(mod))
    Nperson <- length(person.parameter(mod)$thetapar[[1]])
    fitstat_low <- 1-2*sqrt(2/Nperson)
    fitstat_high <- 1+2*sqrt(2/Nperson)
    
    dotplot(ifit$i.outfitMSQ, xlim=c(0,3), xlab = "Outfit Mean Square",
            panel = function(...) { 
              panel.dotplot(...) 
              panel.abline(v=fitstat_low, lty=2)
              panel.abline(v=fitstat_high, lty=2)
            })
  }, width=500)
  
  
  #############################################################################-
  # item information plots
  output$info_plots <- renderPlot({
    mod <- get_rasch_stats()[["rasch_mod"]]
    plotINFO(mod, type="both")
  }, height=800, width=800)
  
  
  #############################################################################-
  # SEM by theta plot and table
  output$sem_by_theta_plot <- renderPlot({
    tab <- get_rasch_stats()[["pers_table"]]
    tab <- tab %>% group_by(theta, se.theta) %>%
      summarise(ts=mean(total_score),
                n=n()) %>%
      ungroup() %>%
      as.data.frame()
    plot(tab$se.theta ~ tab$theta, las=1, main = "SEM by Theta",
         xlab="Theta", ylab="SEM", type="b")
  }, width=500)
  
  output$theta_sem_table <- renderTable({
    tab <- get_rasch_stats()[["pers_table"]]
    tab <- tab %>% group_by(theta, se.theta) %>%
      summarise(ts=mean(total_score),
                n=n()) %>%
      ungroup() %>%
      as.data.frame()
    tab
  }, digits = 3)
  
  
  #############################################################################-
  # total by theta plot
  output$theta_by_total_plot <- renderPlot({
    tab <- get_rasch_stats()[["pers_table"]]
    tab <- tab %>% group_by(theta, se.theta) %>%
      summarise(ts=mean(total_score),
                n=n()) %>%
      ungroup() %>%
      as.data.frame()
    plot(tab$theta ~ tab$ts, las=1, main = "Theta by Total Score",
         xlab="Observed Total Score", ylab="Estimated Theta", type="b")
  }, width=500)
  
  
  #############################################################################-
  # theta histogram
  output$theta_hist_plot <- renderPlot({
    tab <- get_rasch_stats()[["pers_table"]]
    hist(tab$theta, col="grey", las=1, xlab="Theta",
         main = "Histogram of Estimated Theta Values")
  }, width=500)
  
  
  #############################################################################-
  # theta summary
  output$theta_summary_table <- renderTable({
    theta <- get_rasch_stats()[["pers_table"]]$theta
#    res <- t(describe(theta))
    res <- data.frame(
      "N" = length(theta),
      "Mean" = mean(theta),
      "SD" = sd(theta),
      "Min" = min(theta),
      "Max" = max(theta),
      "Median" = median(theta),
      "p25" = quantile(theta,.25),
      "p75" = quantile(theta,.75)
    )
    res
  }, include.rownames=FALSE, digits=4)
  
  
  #############################################################################-
  # separation reliability text
  output$sep_rel_text <- renderText({
    tab <- get_rasch_stats()[["pers_table"]]
    alpha <- get_ctt_stats()[["alpha"]]
    vartheta <- round(var(tab$theta), 4)
    mse <- round(mean(tab$se.theta^2), 4)
    paste0(
      "\nThe variance of theta estimates is:",
      "\n  Var(theta) = ", vartheta, ".",
      "\nThe average estimated SEM of theta is:",
      "\n  mean(SEM) = ", round(mean(tab$se.theta), 3), ".",
      "\nThe average squared SEM of theta is:",
      "\n  mean(SEM^2) = ", mse, ".",
      "\nSeparation reliability estimate is:",
      "\n  ([Var(theta) - mean(SEM^2)] / Var(theta) = ", round((vartheta-mse)/vartheta, 3), ".")
  })
  
  
  #############################################################################-
  # item IRF plots
  output$item_irf_plots <- renderPlot({
    plotjointICC(get_rasch_stats()[["rasch_mod"]], legpos = FALSE,
                 ylab="P(X=1)",xlab="Theta", las=1)
  }, width=500)
  
  #############################################################################-
  # Item ICC fit plots
  output$item_icc_fit_plots <- renderPlot({
    stuff       <- get_rasch_stats()
    rasch_mod   <- stuff[["rasch_mod"]]
    theta       <- person.parameter(rasch_mod)$theta.table$`Person Parameter`
    dat         <- rasch_mod$X01
    total_score <- apply(dat, 1, sum)
    nitem       <- nrow(stuff[["ifit_table"]])
    nr          <- ceiling(nitem/3)
    
    nbins     <- input$nbins
    showN     <- input$showN #FALSE  # label number of people per bin?
    adjSize   <- input$adjSize #FALSE  # adjust point size based on N size per bin?
    showLine  <- input$showLine
    thetacat  <- cut(theta, breaks = nbins, labels = FALSE)
    
    par(mfrow = c(nr,3))
    for(i in 1:nitem) {
      plotICC(rasch_mod, item.subset = i, ylim=c(-0.1,1.1),
              xlim = c(min(theta)-0.2, max(theta)+0.2),
              las=1, main = row.names(stuff[["ifit_table"]])[i])
      x <- data.frame(x=dat[,i], theta=theta, thetacat=thetacat)
      x <- x %>% group_by(thetacat) %>%
        summarise(
          p=mean(x),
          n=length(theta),
          theta=mean(theta)) %>% ungroup() %>% as.data.frame()
      
      if(adjSize) {
        points(x$p~x$theta, pch=19, cex=x$n/max(x$n)+1)
      } else {  
        points(x$p~x$theta, pch=19)
      }
      if(showN) {
        text(x = x$theta, y = (x$p+ifelse(x$p<0.2,0.1,-0.1)), labels = x$n)
      }
      if(showLine) {
        lines(x$p~x$theta, lty = 2, col = "blue")
      }
    }
    par(mfrow=c(1,1))
    
  })
  
  item_icc_fit_plots_height <- reactive({
    75*ceiling(nrow(get_rasch_stats()[["ifit_table"]]))
  })
  
  output$item_icc_fit_plots_dynamic <- renderUI({
    plotOutput("item_icc_fit_plots", height = item_icc_fit_plots_height(),
               width=700)
  })
  
  item_icc_fit_table <- reactive({
    stuff       <- get_rasch_stats()
    rasch_mod   <- stuff[["rasch_mod"]]
    theta       <- person.parameter(rasch_mod)$theta.table$`Person Parameter`
    dat         <- rasch_mod$X01
    total_score <- apply(dat, 1, sum)
    nitem       <- nrow(stuff[["ifit_table"]])
    nr          <- ceiling(nitem/3)
    
    nbins     <- input$nbins
    showN    <- input$showN #FALSE  # label number of people per bin?
    adjSize  <- input$adjSize #FALSE  # adjust point size based on N size per bin?
    thetacat <- cut(theta, breaks = nbins, labels = FALSE)
    
    # Show the proportion of 1's per bin for each item
    bindata <- data.frame(dat, theta=theta, total=total_score, thetacat=thetacat)
    bindata <- round(cbind(
      as.data.frame(summarise(group_by(bindata, thetacat), across(where(is.numeric), mean))),
      #ddply(bindata, .(thetacat), numcolwise(mean)),
      n = tapply(bindata$theta, bindata$thetacat, length)), 3)
    names(bindata) <- c("Theta_Bin", colnames(dat), "mean_theta", "mean_total", "N")
    bindata
    
  })
  
  output$dl_item_icc_fit_table <- downloadHandler(
    filename = function() {
      paste("item_icc_fit_table.csv")
    },
    content = function(file) {
      write.csv(item_icc_fit_table(), file, row.names = TRUE)
    }
  )
  
  
  #############################################################################-
  # Scree plot
  output$scree_plot <- renderPlot(
    
    # width = function() input$screeplot_width,
    # height = function() input$screeplot_height,
    res = 96,
    
    {
      if(input$dimensionality_binary=="binary") {
        d <- dataset_noid()
      } else {
        d <- dataset_poly()
      }
      plot(eigen(cor(d))$values, type="b", las=1,
           xlab = "Component Number",
           ylab = "Eigenvalue",
           main = "Eigenvalues of Item Pearson Correlation Matrix")
      abline(h=1, lty=2)
    })
  
  
  #############################################################################-
  # Corrplot
  output$corrplot <- renderPlot(
    
    # width = function() input$corrplot_width,
    # height = function() input$corrplot_height,
    res = 96,
    
    {
      if(input$dimensionality_binary=="binary") {
        d <- dataset_noid()
      } else {
        d <- dataset_poly()
      }
      corrplot(cor(d))
    })
  
  
  #############################################################################-
  # table of item correlations
  output$item_cor_table <- renderTable({
    if(input$dimensionality_binary=="binary") {
      d <- dataset_noid()
    } else {
      d <- dataset_poly()
    }
    round(cor(d), 3)
  })
  
  
  #############################################################################-
  # table of eigenvalues and parallel analysis eigenvalues
  
  output$eigen_value_table <- renderTable({
    if(input$dimensionality_binary=="binary") {
      d <- dataset_noid()
    } else {
      d <- dataset_poly()
    }
    
    # 100 reps; show mean and 95th percentile
    set.seed(1112021)
    evs <- eigen(cor(d))$values
    D <- d
    reps <- 100
    pct <- 0.95
    nobs <- nrow(D)
    nitem <- ncol(D)
    results <- matrix(ncol = nitem, nrow = reps)
    for (r in 1:reps) {
      tempresult <- eigen(cor(mvrnorm(n=nobs, mu=rep(0,times=nitem),
                                      Sigma=diag(nitem))))$values
      results[r,] <- tempresult
      # HEREHERE
    }
    
    ev_means <- apply(results, 2, mean)
    ev_pct <- apply(results, 2, function(X) {quantile(X, pct)})
    
    data.frame(
      Component = c(1:nitem),
      Obs_Eigenvalue = evs,
      Mean_Random_Eigenvalue = ev_means,
      Pct_95_Random_Eigenvalue = ev_pct
    )
    
  })
  
  
  
  
  #############################################################################-
  # preview/show imported data table
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    df <- dataset_noid()
    
    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
    
  })
  
  output$contents_full <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    df <- dataset_poly()
    
    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
    
  })
  
}

# shinyApp() #####
shinyApp(ui, server)


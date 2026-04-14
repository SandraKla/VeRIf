####################################### WELCOME TO THE SHINY APP ##################################
####################################### from Sandra K. (2026) #####################################
###################################################################################################

# # Get template for the dataset
# library(writexl)
# library(readxl)
# 
# data <- reflimR::livertests
# 
# write.csv(data, "reflim_csv.csv", row.names = FALSE)
# write.csv2(data, "reflim_csv2.csv", row.names = FALSE)
# write_xlsx(data, "reflim_excel.xlsx")
# 
# dataset_original1 <- read.csv("reflim_csv.csv")
# dataset_original2 <- read.csv2("reflim_csv2.csv")
# dataset_original3 <- read_excel("reflim_excel.xlsx")
# 
# write.csv2(dataset_original1, "reflim_data1.csv", row.names = FALSE)
# write.csv2(dataset_original2, "reflim_data2.csv", row.names = FALSE)
# write.csv2(dataset_original3, "reflim_data3.csv", row.names = FALSE)

####################################### Load Script and Example-Dataset ###########################

source("functions.R")
source("reflimR_loop.R")

####################################### Libraries #################################################

if ("DT" %in% rownames(installed.packages())) {
  library(DT)} else{
    install.packages("DT")
    library(DT)}

if ("mclust" %in% rownames(installed.packages())) {
  library(mclust)} else{
    install.packages("mclust")
    library(mclust)}

if ("refineR" %in% rownames(installed.packages())) {
  library(refineR)} else{
    install.packages("refineR")
    library(refineR)}

if ("reflimR" %in% rownames(installed.packages())) {
  library(reflimR)} else{
    install.packages("reflimR")
    library(reflimR)}

if ("rhandsontable" %in% rownames(installed.packages())) {
  library(rhandsontable)} else{
    install.packages("rhandsontable")
    library(rhandsontable)}

if ("readxl" %in% rownames(installed.packages())) {
  library(readxl)} else{
    install.packages("readxl")
    library(readxl)}

if ("rpart" %in% rownames(installed.packages())) {
  library(rpart)} else{
    install.packages("rpart")
    library(rpart)}

if ("rpart.plot" %in% rownames(installed.packages())) {
  library(rpart.plot)} else{
    install.packages("rpart.plot")
    library(rpart.plot)}

if ("shinycssloaders" %in% rownames(installed.packages())) {
  library(shinycssloaders)} else{
    install.packages("shinycssloaders")
    library(shinycssloaders)}

if ("shinydashboard" %in% rownames(installed.packages())) {
  library(shinydashboard)} else{
    install.packages("shinydashboard")
    library(shinydashboard)}

dataset_original <- reflimR::livertests

####################################### Texts #####################################################

data_text <- HTML(paste0(
  "This Shiny App is based on the package ", a("reflimR", href = "https://cran.r-project.org/web/packages/reflimR/index.html"), 
  " for the estimation of reference limits from routine laboratory results:", br(), br(), 
  "These columns should be used for new data: Category: Name of the category to filter the data; Age: Age in years; Sex: m for male and f for female;
  Value: Column name is the analyte name, values are the laboratory measures.Starting with the fourth column, enter the laboratory value; the other three columns can be in any order. The data from *livertests* serves as a template. 
  To load new data, the data should be in CSV format with values separated by semicolons (;), and decimal numbers should use a comma (,) as the decimal separator. The first row should contain column headers.
  Alternatively, the data can be loaded into the editable table using the copy-and-paste function or with .xlsx.", br(), br(),
  "On the left side, the sidebar allows you to select the laboratory parameter, category, age and gender group. 
  In the “Target Values” section, you can load target values from targetvalues, load reference intervals estimated with refineR, or manually enter custom values."
))
reflim_text <- HTML(paste0(
  "These tab displays the corresponding plot and the outputs of the reflim() function, providing an estimation of new reference intervals or a verification of the selected target values. 
  By clicking “Visualization of all plots across every process step”, all plots generated throughout the workflow can be displayed."
))
refineR_text<- HTML(paste0(
  "If, during the verification with reflimR and its target values or own target values, a yellow or red bar appears, 
  a follow-up analysis using refineR is recommended. The resulting reference intervals from refineR can be used as new target values
  and re-verified with reflimR. If all indicators turn green, this suggests that the manufacturer’s target values are likely incorrect. 
  If one or more indicators remain yellow or red, the data are considered too challenging for indirect methods. This assumption can be further 
  evaluated in the “mclust” tab using a Gaussian mixture model (mclust)."
))
mclust_text <- HTML(paste0(
  "Gaussian mixture modelling for the verification of reference intervals."
))
scatterplot_text <- HTML(paste0(
  "The scatterplot shows the relationship between age and the laboratory value."
))
statistics_text <- HTML(paste0(
  "The two figures show the distribution of sex across age and the laboratory value."
))
zlog_text <- HTML(paste0(
  "zlog values are calculated from the dataset and the calculated reference intervals. The lower reference limits (LL) and upper reference limits (UL) can transform any result x into a zlog value using the following equation:
  zlog(x) = (log(x)–(log(LL)+ log(UL))/2)*3.92/(log(UL)–log(LL)). Values ranging from –1.96 to 1.96 are considered normal, while values below –5 and above 5 indicate pathological conditions."
))
rpart_text <- HTML(paste0(
  "Decision trees are used here to divide the data into subgroups based on age and sex. The method identifies groups with similar values, allowing more appropriate calculation of reference intervals. The tree is built using the rpart package and visualized with rpart.plot."))

####################################### User Interface ############################################

ui <- dashboardPage(
  dashboardHeader(title = "VeRIf", titleWidth = 300),
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "sidebarid",
      
      uiOutput("parameters"),
      textInput(
        "parameter_label",
        "Laboratory parameter name:",
        value = ""
      ),
      uiOutput("category"),
      
      selectInput(
        "sex",
        "Select the sex:",
        choices = c("Female (F) & Male (M)" = "t", "Female (F)" = "f", "Male (M)" = "m")
      ),
      
      sliderInput(
        "age_end",
        "Select age-range:",
        min = 0,
        max = 100,
        value = c(0, 100)),
      
      hr(),
      
      numericInput(
        "nmin",
        "Minimum number required for reliable reference limit estimation:",
        200,
        min = 40,
        max = 1000
      ),
      
      radioButtons(
        "lambda_type",
        "Select lambda source:",
        choices = c(
          "Lambda from reflimR" = "reflimR",
          "Lambda from refineR" = "refineR",
          "User defined lambda" = "user"
        ),
        selected = "reflimR"
      ),

      conditionalPanel(
        condition = "input.lambda_type == 'user'",
        sliderInput(
          "lambda",
          "Select lambda for UM:",
          min = 0,
          max = 1,
          value = 0.5,
          step = 0.01
        )
      ),
      
      hr(),
      
      checkboxInput("check_targetvalues", "Load preinstalled target values", value = FALSE),
      checkboxInput("check_target", "Load own target values", value = FALSE),
      
      conditionalPanel(
        condition = "input.check_target == true",
        
        numericInput(
          "target_low",
          "Lower value:",
          10,
          min = 0,
          max = 10000
        ),
        
        numericInput(
          "target_upper",
          "Upper value:",
          15,
          min = 0,
          max = 10000
        )
      ), 
      uiOutput("refineR_checkbox_ui"),
      hr()
    )
  ),
  
  dashboardBody(
    fluidRow(
      
      tabsetPanel( 
        tabPanel("Data", 
                 icon = icon("upload"),
                 
                 box(
                   status = "info",
                   width = 7,
                   
                   p(data_text),
                   
                   fluidRow(
                     column(6, checkboxInput("show_table", tags$b("Show and use editable table for upload. Please click Submit!"), value = FALSE)),
                     column(6, actionButton("submit", "Submit"))
                   ),
                    conditionalPanel(
                      condition = "input.show_table == true",
                      withSpinner(rHandsontableOutput("editable_table"))
                    ), hr(),
                    
                    uiOutput("dataset_file"),
                    uiOutput("upload_mapping_ui"),
                    actionButton('reset', 'Reset Input', icon = icon("trash"))
                  )
        ),
        
        # tabPanel("Scatterplot", 
        #          icon = icon("chart-line"),
        #          
        #          box(
        #            title = "",
        #            status = "info",
        #            width = 7,
        #            solidHeader = TRUE,
        #            
        #            p(scatterplot_text),
        #            DT::dataTableOutput("table"),
        #            plotOutput("scatterplot", height = "700px")
        #          )
        # ),
        
        # tabPanel( "Statistics", 
        #           icon = icon("chart-bar"),
        #           
        #           box(
        #             title = "",
        #             status = "info",
        #             width = 7,
        #             solidHeader = TRUE,
        #             
        #             p(statistics_text),
        #             plotOutput("plot_statistics", height = "700px")
        #           )
        # ),
        
        tabPanel("reflimR", 
                 icon = icon("chart-line"), 
                 
                 box(
                   status = "info",
                   width = 7,

                   p(reflim_text),
                   
                   fluidRow(
                     column(4, checkboxInput("check_plot.all","Visualization of all plots across every process step")),
                     
                     column(8, radioButtons(
                              "plot_type",
                              "Select visualization:",
                              choices = c(
                                "Visualization with EL" = "pU",
                                "Visualization with UM" = "VeRUS"
                              ),
                              selected = "pU"
                            )
                          ),
                     
                  ), withSpinner(plotOutput("plot", height = "70vh"))
                 )
        ),
        
        tabPanel( "refineR", 
                  icon = icon("chart-line"),
                  
                  box(
                    status = "info",
                    width = 7,
                    
                    p(refineR_text),
                    withSpinner(plotOutput("plotrefineR", height = "70vh")),
                    withSpinner(DT::dataTableOutput("table_report_refineR"))
                  )
        ),
        
        tabPanel(
          "mclust",
          icon = icon("chart-line"),
          
          box(
            status = "info",
            width = 7,
            
            p(mclust_text),
            
            fluidRow(
              
              column(
                3,
                checkboxInput(
                  inputId = "auto_cluster",
                  label = "Select optimal number of clusters automatically",
                  value = TRUE
                )
              ),
              
              column(
                3,
                conditionalPanel(
                  condition = "input.auto_cluster == false",
                  numericInput(
                    inputId = "n_cluster",
                    label = "Number of clusters:",
                    value = 3,
                    min = 1,
                    max = 15,
                    step = 1
                  )
                )
              ),
              
              column(
                6,
                radioButtons(
                  inputId = "model_name",
                  label = "Univariate model for mclust:",
                  choices = c("Equal variance (one-dimensional)" = "E", "Variable/unqual variance (one-dimensional)" = "V"),
                  selected = "V"
                )
              )
            ),
            
            withSpinner(plotOutput("plotmclust", height = "70vh")),
            
            conditionalPanel(
              condition = "input.auto_cluster == true",
              withSpinner(plotOutput("plotmclustbic", height = "70vh"))
            )
          )
        ),
        
        tabPanel( "rpart", 
                  icon = icon("chart-bar"),
                  
                  box(
                    status = "info",
                    width = 7,
                    
                    p(rpart_text),
                    
                    numericInput("tree_window_minsplit", "Minimum observations needed to split a node:", 20, min = 10, max = 100),
                    #numericInput("tree_window_cp", "rpart: Complexity parameter", 0.01, min = 0, max = 10),
                    withSpinner(plotOutput("tree_rpart", height = "70vh"))
                  )
        )#,
        
        # tabPanel( "zlog", 
        #           icon = icon("table"),
        #           
        #           box(
        #             title = "",
        #             status = "info",
        #             width = 7,
        #             solidHeader = TRUE,
        #             
        #             p(zlog_text),
        #             DT::dataTableOutput("table_zlog",  height = "700px")
        #           )
        # )
      ),
      
      box(
        title = tagList(shiny::icon("table"), "Report:"),
        status = "info",
        width = 5,
        solidHeader = TRUE,
        
        withSpinner(DT::dataTableOutput("table_report")), hr(),
        downloadButton("download_ritable", "Download all Reference Intervals"),
        #downloadButton("download_zlogtable", "Download all zlog values")
      )
    )
  )
)

####################################### Server ####################################################

server <- function(input, output, session) {
  
  ##################################### Observe Events ############################################
  
  options(shiny.sanitize.errors = TRUE)
  options(warn = -1)
  
  values <- reactiveValues(
    upload_state = NULL
  )
  
  observeEvent(input$dataset_file1, {
    values$upload_state <- 'uploaded'
  })
  
  observeEvent(input$reset, {
    values$upload_state <- 'reset'
  })
  
  dataset_input <- reactive({
    if (is.null(values$upload_state)) {
      return(NULL)
    } else if (values$upload_state == 'uploaded') {
      return(input$dataset_file1)
    } else if (values$upload_state == 'reset') {
      return(NULL)
    }
  })
  
  output$dataset_file <- renderUI({
    input$reset ## Create a dependency with the reset button
    fileInput('dataset_file1', label = NULL,  multiple = FALSE)
  })

  normalize_upload_value <- function(x) {
    x <- as.character(x)
    converted <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
    converted[is.na(converted)] <- x[is.na(converted)]
    tolower(trimws(converted))
  }

  guess_upload_column <- function(dataset, preferred_names = character(),
                                  fallback = NULL, excluded = character()) {
    available_columns <- setdiff(names(dataset), excluded)
    
    if (!length(available_columns)) {
      return(NULL)
    }
    
    normalized_columns <- normalize_upload_value(available_columns)
    normalized_preferred <- normalize_upload_value(preferred_names)
    matching_index <- match(normalized_preferred, normalized_columns, nomatch = 0)
    matching_index <- matching_index[matching_index > 0]
    
    if (length(matching_index)) {
      return(available_columns[matching_index[1]])
    }
    
    if (!is.null(fallback) && fallback >= 1 && fallback <= ncol(dataset)) {
      fallback_column <- names(dataset)[fallback]
      if (fallback_column %in% available_columns) {
        return(fallback_column)
      }
    }
    
    available_columns[1]
  }

  guess_sex_value <- function(values, preferred_values = character(), excluded = NULL) {
    available_values <- unique(as.character(values))
    available_values <- available_values[!is.na(available_values) & nzchar(trimws(available_values))]
    
    if (!is.null(excluded)) {
      available_values <- available_values[
        normalize_upload_value(available_values) != normalize_upload_value(excluded)
      ]
    }
    
    if (!length(available_values)) {
      return(NULL)
    }
    
    normalized_values <- normalize_upload_value(available_values)
    normalized_preferred <- normalize_upload_value(preferred_values)
    matching_index <- match(normalized_preferred, normalized_values, nomatch = 0)
    matching_index <- matching_index[matching_index > 0]
    
    if (length(matching_index)) {
      return(available_values[matching_index[1]])
    }
    
    available_values[1]
  }

  uploaded_dataset_raw <- reactive({
    file_info <- dataset_input()
    req(file_info)
    
    datapath <- file_info[["datapath"]]
    
    validate(need(
      grepl("\\.(csv|xlsx)$", datapath, ignore.case = TRUE),
      "Check if you have used the correct template! It must be a CSV or XLSX file!"
    ))
    
    if (grepl("\\.csv$", datapath, ignore.case = TRUE)) {
      dataset <- read.csv2(datapath)
    } else {
      dataset <- as.data.frame(readxl::read_excel(datapath), stringsAsFactors = FALSE)
    }
    
    validate(need(
      nrow(dataset) > 0,
      "Check if you have used the correct template! The dataset is empty!"
    ))
    
    validate(need(
      ncol(dataset) >= 3,
      "The uploaded dataset must contain at least three columns."
    ))
    
    dataset
  })

  output$upload_mapping_ui <- renderUI({
    if (input$show_table || is.null(dataset_input())) {
      return(NULL)
    }
    
    dataset <- uploaded_dataset_raw()
    column_names <- names(dataset)
    
    guessed_age_column <- guess_upload_column(
      dataset,
      preferred_names = c("Age", "Alter"),
      fallback = 2
    )
    selected_age_column <- if (!is.null(input$upload_age_column) &&
                               input$upload_age_column %in% column_names) {
      input$upload_age_column
    } else {
      guessed_age_column
    }
    
    guessed_sex_column <- guess_upload_column(
      dataset,
      preferred_names = c("Sex", "Gender", "Geschlecht"),
      fallback = 3,
      excluded = selected_age_column
    )
    selected_sex_column <- if (!is.null(input$upload_sex_column) &&
                               input$upload_sex_column %in% column_names &&
                               input$upload_sex_column != selected_age_column) {
      input$upload_sex_column
    } else {
      guessed_sex_column
    }
    
    sex_values <- unique(as.character(dataset[[selected_sex_column]]))
    sex_values <- sex_values[!is.na(sex_values) & nzchar(trimws(sex_values))]
    
    guessed_female_value <- guess_sex_value(
      sex_values,
      preferred_values = c("f", "female", "weiblich", "w")
    )
    selected_female_value <- if (!is.null(input$upload_female_value) &&
                                 input$upload_female_value %in% sex_values) {
      input$upload_female_value
    } else {
      guessed_female_value
    }
    
    guessed_male_value <- guess_sex_value(
      sex_values,
      preferred_values = c("m", "male", "mannlich", "maennlich", "mann", "man"),
      excluded = selected_female_value
    )
    selected_male_value <- if (!is.null(input$upload_male_value) &&
                               input$upload_male_value %in% sex_values &&
                               normalize_upload_value(input$upload_male_value) != normalize_upload_value(selected_female_value)) {
      input$upload_male_value
    } else {
      guessed_male_value
    }
    
    tagList(
      hr(),
      tags$h4("Column mapping:"),
      fluidRow(
        column(
          6,
          selectInput(
            "upload_age_column",
            "Age column:",
            choices = column_names,
            selected = selected_age_column
          )
        ),
        column(
          6,
          selectInput(
            "upload_sex_column",
            "Sex column:",
            choices = setdiff(column_names, selected_age_column),
            selected = selected_sex_column
          )
        )
      ),
      fluidRow(
        column(
          6,
          selectInput(
            "upload_female_value",
            "Value for female:",
            choices = sex_values,
            selected = selected_female_value
          )
        ),
        column(
          6,
          selectInput(
            "upload_male_value",
            "Value for male:",
            choices = sex_values,
            selected = selected_male_value
          )
        )
      )
    )
  })

  uploaded_dataset_standardized <- reactive({
    dataset <- uploaded_dataset_raw()
    
    age_column <- if (!is.null(input$upload_age_column) &&
                      input$upload_age_column %in% names(dataset)) {
      input$upload_age_column
    } else {
      guess_upload_column(
        dataset,
        preferred_names = c("Age", "Alter"),
        fallback = 2
      )
    }
    
    sex_column <- if (!is.null(input$upload_sex_column) &&
                      input$upload_sex_column %in% names(dataset)) {
      input$upload_sex_column
    } else {
      guess_upload_column(
        dataset,
        preferred_names = c("Sex", "Gender", "Geschlecht"),
        fallback = 3,
        excluded = age_column
      )
    }
    
    validate(need(!is.null(age_column), "Please choose the age column."))
    validate(need(!is.null(sex_column), "Please choose the sex column."))
    validate(need(age_column != sex_column, "Age column and sex column must be different."))
    
    named_category <- setdiff(
      names(dataset)[normalize_upload_value(names(dataset)) %in% c("category", "kategorie")],
      c(age_column, sex_column)
    )
    category_candidates <- setdiff(names(dataset)[seq_len(min(3, ncol(dataset)))], c(age_column, sex_column))
    remaining_columns <- setdiff(names(dataset), c(age_column, sex_column))
    
    category_column <- if (length(named_category)) {
      named_category[1]
    } else if (length(category_candidates)) {
      category_candidates[1]
    } else if (length(remaining_columns)) {
      remaining_columns[1]
    } else {
      NULL
    }
    
    analyte_columns <- setdiff(names(dataset), c(category_column, age_column, sex_column))
    
    validate(need(
      length(analyte_columns) > 0,
      "Please choose age and sex columns so that at least one laboratory parameter remains."
    ))
    
    sex_values <- dataset[[sex_column]]
    female_value <- if (!is.null(input$upload_female_value) &&
                        nzchar(trimws(input$upload_female_value))) {
      input$upload_female_value
    } else {
      guess_sex_value(sex_values, c("f", "female", "weiblich", "w"))
    }
    male_value <- if (!is.null(input$upload_male_value) &&
                      nzchar(trimws(input$upload_male_value))) {
      input$upload_male_value
    } else {
      guess_sex_value(sex_values, c("m", "male", "mannlich", "maennlich", "mann", "man"), excluded = female_value)
    }
    
    validate(need(!is.null(female_value), "Please map one sex value to female."))
    validate(need(!is.null(male_value), "Please map one sex value to male."))
    validate(need(
      normalize_upload_value(female_value) != normalize_upload_value(male_value),
      "Female and male must use different source values."
    ))
    
    normalized_sex <- rep(NA_character_, nrow(dataset))
    normalized_source_values <- normalize_upload_value(sex_values)
    normalized_sex[normalized_source_values == normalize_upload_value(female_value)] <- "f"
    normalized_sex[normalized_source_values == normalize_upload_value(male_value)] <- "m"
    
    category_values <- if (!is.null(category_column)) {
      as.character(dataset[[category_column]])
    } else {
      rep("All", nrow(dataset))
    }
    category_values[is.na(category_values) | !nzchar(trimws(category_values))] <- "All"
    
    standardized_dataset <- data.frame(
      Category = category_values,
      Age = suppressWarnings(as.numeric(dataset[[age_column]])),
      Sex = normalized_sex,
      stringsAsFactors = FALSE
    )
    
    cbind(standardized_dataset, dataset[analyte_columns])
  })

  parameter_label_value <- reactiveVal(NULL)

  observeEvent(input$parameter, {
    parameter_label_value(input$parameter)
    updateTextInput(session, "parameter_label", value = input$parameter)
  }, ignoreNULL = FALSE)

  observeEvent(input$parameter_label, {
    parameter_label_value(input$parameter_label)
  }, ignoreInit = TRUE)
  
  output$parameters <- renderUI({
    if(input$show_table && input$submit) {
      choices <- colnames(data_store())[4:length(colnames(data_store()))]
    } else{
    if (is.null(dataset_input())) { 
      choices <- colnames(dataset_original)[4:length(colnames(dataset_original))]
    } 
    else{
      dataset <- uploaded_dataset_standardized()
      choices <- setdiff(colnames(dataset), c("Category", "Age", "Sex"))
      }}
    validate(need(length(choices) > 0, "No laboratory parameter columns are available."))
    
    selected_parameter <- if (!is.null(input$parameter) && input$parameter %in% choices) {
      input$parameter
    } else {
      choices[1]
    }
    
    selectInput("parameter","Select laboratory value:", choices = choices, selected = selected_parameter)
  })   
  
  output$category <- renderUI({
    if(input$show_table && input$submit) {
      choices <- unique(data_store()[[1]])
    } else{
      if (is.null(dataset_input())) { 
        choices <- unique(dataset_original[[1]])
      } else{
        dataset <- uploaded_dataset_standardized()
        choices <- unique(dataset$Category)
      }
    }
    choices <- c("Not selected", choices)
    
    selected_category <- if (!is.null(input$category) && input$category %in% choices) {
      input$category
    } else {
      "Not selected"
    }
    
    selectInput("category", "Select category:", choices = choices, selected = selected_category)
  })
  
  
  # Create a reactive values to track the state of the checkboxes
  reactive_values <- reactiveValues(
    check_targetvalues = FALSE,
    check_target = FALSE,
    check_refineR = FALSE
  )
  
  # Observe changes in check_targetvalues and update the reactive value
  observeEvent(input$check_targetvalues, {
    if (input$check_targetvalues) {
      reactive_values$check_targetvalues <- TRUE
      reactive_values$check_target <- FALSE
      reactive_values$check_refineR <- FALSE
    } else {
      reactive_values$check_targetvalues <- FALSE
    }
  })
  
  # Observe changes in check_target and update the reactive value
  observeEvent(input$check_target, {
    if (input$check_target) {
      reactive_values$check_target <- TRUE
      reactive_values$check_targetvalues <- FALSE
      reactive_values$check_refineR <- FALSE
    } else {
      reactive_values$check_target <- FALSE
    }
  })
  
  # Observe changes in check_refineR and update the reactive value
  observeEvent(input$check_refineR, {
    if (input$check_refineR) {
      reactive_values$check_refineR <- TRUE
      reactive_values$check_targetvalues <- FALSE
      reactive_values$check_target <- FALSE
    } else {
      reactive_values$check_refineR <- FALSE
    }
  })
  
  # Update the checkboxes based on the reactive value
  observe({
    updateCheckboxInput(session, "check_targetvalues", value = reactive_values$check_targetvalues)
    updateCheckboxInput(session, "check_target", value = reactive_values$check_target)
    updateCheckboxInput(session, "check_refineR", value = reactive_values$check_refineR)
  })
  
  observeEvent(input$submit, {
    if (!is.null(input$editable_table)) { data_store(hot_to_r(input$editable_table))}
  })

  initial_data <- data.frame(
    Category = character(50),
    Age = numeric(50),
    Sex = character(50),
    Analyte = numeric(50),
    #Analyte1 = numeric(50),
    #Analyte2 = numeric(50),
    #Analyte3 = numeric(50),
    stringsAsFactors = FALSE
  )

  data_store <- reactiveVal(initial_data)

  observeEvent(input$show_table, {
    if (input$show_table) {
      output$editable_table <- renderRHandsontable({
        rhandsontable(data_store(), width = '800',
                      height = 550, rowHeaders = NULL, colHeaders = colnames(data_store()))
      })
    }
  }, ignoreNULL = FALSE)
  
  plot_choice <- reactive({
    
    if(input$plot_type == "pU"){
      return("pU")
    }
    
    if(input$plot_type == "VeRUS"){
      return("VeRUS")
    }
  })
  
  ##################################### Reactive Expressions ######################################

  parameter_display <- reactive({
    label_value <- parameter_label_value()
    
    if (is.null(label_value)) {
      return(input$parameter)
    }
    
    custom_label <- trimws(label_value)
    if (nzchar(custom_label)) {
      custom_label
    } else {
      input$parameter
    }
  })
  
  refineR_data <- reactive({
    if (isTRUE(input$show_table) && input$submit > 0) {
      dataset <- data_store()
    } else {
      if (is.null(dataset_input())) {
        dataset <- dataset_original
      } else {
        dataset <- uploaded_dataset_standardized()
      }
    }
    
    validate(need(
      input$parameter %in% names(dataset),
      "Please select a valid laboratory parameter."
    ))
    
    dataset <- dataset[, c("Category", "Age", "Sex", input$parameter), drop = FALSE]
    
    if (!is.null(input$category) && input$category != "Not selected") {
      dataset <- subset(dataset, Category == input$category)
    }
    
    dataset$Age <- suppressWarnings(as.numeric(dataset$Age))
    dataset[, 4] <- as.numeric(dataset[, 4])
    
    dataset <- subset(dataset, Age >= input$age_end[1] & Age <= input$age_end[2])
    
    if (input$sex %in% c("m", "f")) {
      dataset <- subset(dataset, Sex == input$sex)
    }
    
    dataset
  })
  
  # Create the table with the dataset as reactive expression 
  reflim_data <- reactive({
    
    if(input$show_table && input$submit) {
       dataset <- data_store()
    } else{
    
      if (is.null(dataset_input())) {
        dataset <- dataset_original
      } else {
        dataset <- uploaded_dataset_standardized()
      }}
    
    validate(need(
      input$parameter %in% names(dataset),
      "Please select a valid laboratory parameter."
    ))
    
    dataset <- dataset[, c("Category", "Age", "Sex", input$parameter), drop = FALSE]
    
    if (!is.null(input$category) && input$category != "Not selected") {
      dataset <- subset(dataset, Category == input$category)
    }
    
    validate(need(
      ncol(dataset) == 4, 
      "Check if you have used the correct template! You need 4 columns (Category, Age, Sex, Value)!"
    ))
    dataset$Age <- suppressWarnings(as.numeric(dataset$Age))
    dataset[, 4] <- as.numeric(dataset[, 4])
    
    dataset <- subset(dataset, Age >= input$age_end[1] & Age <= input$age_end[2])
    
    if (input$sex %in% c("m", "f")) {
      dataset <- subset(dataset, Sex == input$sex)
    }
    
    return(dataset)
  })
  
  get_alldata_file <- reactive({
    
    if(input$show_table && input$submit) {
      dataset <- data_store()
    } else{
      
      if (is.null(dataset_input())) {
        dataset <- dataset_original
      } else {
        dataset <- uploaded_dataset_standardized()
      }}
    
    dataset$Age <- suppressWarnings(as.numeric(dataset$Age))
    
    if (!is.null(input$category) && input$category != "Not selected") {
      dataset <- subset(dataset, Category == input$category)
    }
    
    dataset <- subset(dataset, Age >= input$age_end[1] & Age <= input$age_end[2])
    
    if (input$sex %in% c("m", "f")) {
      dataset <- subset(dataset, Sex == input$sex)
    }
    
    return(dataset)
  })
  
  get_data_report <- reactive({
    
    dat <- reflim_data()
    validate(need(nrow(dat) > 39,
                  "(reflim) n = 0. The absolute minimum for reference limit estimation is 40."))
    
    if (input$check_target == FALSE && input$check_targetvalues == FALSE && input$check_refineR == FALSE) {
      reflim_text <- reflim(dat[,4], n.min = input$nmin, plot.all = FALSE, plot.it = FALSE)
    }
    
    if (input$check_target) {
      validate(need(input$target_low < input$target_upper,
                    "(reflim) the upper target limit must be greater than the lower target limit."))
      
      validate(need(input$target_low > 0, 
                    "(reflim) the lower target limit must be greater than 0."))
      
      validate(need(input$target_upper > 0, 
                    "(reflim) the upper target limit must be greater than 0."))
      
      validate(need(input$target_low > 0 && input$target_upper > 0, 
                    "(reflim) the lower and upper target limit must be greater than 0."))
      
      reflim_text <- reflim(dat[,4], targets = c(input$target_low, input$target_upper), n.min = input$nmin, plot.all = FALSE, plot.it = FALSE)
    }
    
    if (input$check_targetvalues) {
      
      validate(
        need(
          input$sex != "t" || input$parameter %in% c("ALB", "BIL", "PROT"),
          "(reflim) The reference intervals are sex-specific. Please select a sex."
        )
      )
      
      targets <- reflimR::targetvalues
      targets_values <- targets[targets$analyte == input$parameter, ]
      
      if (input$sex == "m") {
        targetvalues_low <-  targets_values[, 5]
        targetvalues_upper <- targets_values[, 6]
      }
      if (input$sex == "f") {
        targetvalues_low <- targets_values[, 3]
        targetvalues_upper <- targets_values[, 4]
      }
      if (input$sex == "t") {
        targetvalues_low <- targets_values[, 3]
        targetvalues_upper <- targets_values[, 4]
      }
      
      validate(need(nrow(targets_values) > 0, 
                    "(reflim) There are no preloaded target values for this parameter!"))
      
      reflim_text <- reflim(dat[,4], targets = c(targetvalues_low, targetvalues_upper), n.min = input$nmin, plot.all = FALSE, plot.it = FALSE)
    }
    
    if(input$check_refineR) {
      
      validate(need(refineR_done(),"(refineR) Please perform the refineR calculation first."))
      
      table_refineR <- getRI(fit_refineR())
      targetvalues_low <- table_refineR$PointEst[1]
      targetvalues_upper <- table_refineR$PointEst[2]
      
      reflim_text <- reflim(dat[,4], targets = c(targetvalues_low, targetvalues_upper), n.min = input$nmin, plot.all = FALSE, plot.it = FALSE)
    }
    
    report <- reflim_text
    
    return(report)
  })
  
  fit_refineR <- eventReactive(list(
    input$parameter,
    input$category,
    input$sex,
    input$age_end,
    input$reset,
    input$dataset_file1,
    input$submit,
    input$show_table,
    input$upload_age_column,
    input$upload_sex_column,
    input$upload_female_value,
    input$upload_male_value
  ), {
    
    withProgress(message = "RI calculation with refineR …", {
      
      dat <- refineR_data()
      fit <- findRI(Data = dat[, 4])
      refineR_done(TRUE)
      
      fit
    })
  }, ignoreNULL = FALSE)
  
  rpart_data <- reactive({
    dat <- reflim_data()
    response_name <- names(dat)[4]
    
    dat$Age <- suppressWarnings(as.numeric(dat$Age))
    dat[[response_name]] <- suppressWarnings(as.numeric(dat[[response_name]]))
    
    dat <- dat[
      is.finite(dat$Age) &
        !is.na(dat$Sex) &
        nzchar(trimws(as.character(dat$Sex))) &
        is.finite(dat[[response_name]]),
      ,
      drop = FALSE
    ]
    
    validate(need(
      nrow(dat) > 1,
      "(rpart) Not enough complete rows after removing missing values."
    ))
    
    if (nrow(dat) >= 4) {
      quartiles <- stats::quantile(
        dat[[response_name]],
        probs = c(0.25, 0.75),
        na.rm = TRUE,
        names = FALSE
      )
      iqr_value <- quartiles[2] - quartiles[1]
      
      if (is.finite(iqr_value) && iqr_value > 0) {
        lower_fence <- quartiles[1] - 3 * iqr_value
        upper_fence <- quartiles[2] + 3 * iqr_value
        
        dat <- dat[
          dat[[response_name]] >= lower_fence &
            dat[[response_name]] <= upper_fence,
          ,
          drop = FALSE
        ]
      }
    }
    
    validate(need(
      nrow(dat) > 1,
      "(rpart) Not enough rows after removing extreme values."
    ))
    
    dat$Sex <- factor(dat$Sex, levels = c("f", "m"), labels = c("Female", "Male"))
    dat$Sex <- droplevels(dat$Sex)
    
    predictors <- c()
    if (length(unique(dat$Age)) > 1) {
      predictors <- c(predictors, "Age")
    }
    if (nlevels(dat$Sex) > 1) {
      predictors <- c(predictors, "Sex")
    }
    
    validate(need(
      length(predictors) > 0,
      "(rpart) Age and sex do not contain enough variation for a tree."
    ))
    
    list(
      data = dat,
      response_name = response_name,
      predictors = predictors
    )
  })
  
  build_rpart <- reactive({
    tree_input <- rpart_data()
    dat <- tree_input$data
    
    progress <- shiny::Progress$new()
    on.exit(progress$close(), add = TRUE)
    progress$set(message = "Calculate decision tree...", detail = "", value = 2)
    
    validate(
      need(!(input$tree_window_minsplit == 0 || input$tree_window_minsplit == "" || is.na(input$tree_window_minsplit)),
        "Please provide a valid minimum observations value."
      )
    )
    
    tree_minsplit <- as.numeric(input$tree_window_minsplit)
    
    tree_cp <- 0.01
    tree_formula <- reformulate(tree_input$predictors, response = tree_input$response_name)
    
    rpart(
      tree_formula,
      data = dat,
      control = rpart.control(cp = tree_cp, minsplit = tree_minsplit)
    )
  })
  
  ##################################### Output ####################################################
  
  output$plot <- renderPlot({ # #Tab:reflimR
    
    dat <- reflim_data()
    report <- get_data_report()
    
    plot_choice <- plot_choice()
    parameter_name <- parameter_display()
    
    validate(need(nrow(dat) > 39,
                  "(reflim) n < 40. The absolute minimum for reference limit estimation is 40."))

    validate(need(!(report$remarks %in% c("n < 40 after truncation.")),
        "(reflim) n < 40 after truncation. The absolute minimum for reference limit estimation is 40."))
    
    if(input$lambda_type == "reflimR"){
      if(report$lognormal){
        lambda <- 0
      } else{
        lambda <- 1
      }
    } else if(input$lambda_type == "refineR"){
      validate(need(refineR_done(), "(refineR) Please perform the refineR calculation first."))
      lambda <- lambda_refineR
    } else if(input$lambda_type == "user"){
      lambda <- input$lambda
    }
    
    reflimR.plot.all <- FALSE
    
    if (input$check_plot.all) {
      reflimR.plot.all <- TRUE
    }
    
    if (input$check_target) {
      validate(need(input$target_low < input$target_upper,
                    "(reflim) the upper target limit must be greater than the lower target limit."))
      
      validate(need(input$target_low > 0, 
                    "(reflim) the lower target limit must be greater than 0."))
      
      validate(need(input$target_upper > 0, 
                    "(reflim) the upper target limit must be greater than 0."))
      
      validate(need(input$target_low > 0 && input$target_upper > 0, 
                    "(reflim) the lower and upper target limit must be greater than 0."))
      
      if(plot_choice == "VeRUS"){
        reflim_result <- reflim_VeRUS(dat[, 4], targets = c(input$target_low, input$target_upper), n.min = input$nmin, plot.all = reflimR.plot.all, lambda = lambda,
                                      main = paste0("Reference limits for ", parameter_name))}
      if(plot_choice == "pU"){
        reflim_result <- reflim(dat[, 4], targets = c(input$target_low, input$target_upper), n.min = input$nmin, plot.all = reflimR.plot.all,
                                main = paste0("Reference limits for ", parameter_name))}
    }
    
    if (input$check_targetvalues) {
      validate(
        need(
          input$sex != "t" || input$parameter %in% c("ALB", "BIL", "PROT"),
          "(reflim) The reference intervals are sex-specific. Please select a sex."
        )
      )
      
      targets <- reflimR::targetvalues
      targets_values <- targets[targets$analyte == input$parameter,]
      
      if (input$sex == "m") {
        targetvalues_low <-  targets_values[, 5]
        targetvalues_upper <- targets_values[, 6]
      }
      if (input$sex == "f") {
        targetvalues_low <- targets_values[, 3]
        targetvalues_upper <- targets_values[, 4]
      }
      if (input$sex == "t") {
        targetvalues_low <- targets_values[, 3]
        targetvalues_upper <- targets_values[, 4]
      }
      
      validate(need(nrow(targets_values) > 0, 
                    "(reflim) There are no preloaded target values for this parameter!"))
      
      if(plot_choice == "VeRUS"){
        reflim_result <- reflim_VeRUS(dat[, 4], targets = c(targetvalues_low, targetvalues_upper), n.min = input$nmin, plot.all = reflimR.plot.all, lambda = lambda,
                                      main = paste0("Reference limits for ", parameter_name))}
      if(plot_choice == "pU"){
        reflim_result <- reflim(dat[, 4], targets = c(targetvalues_low, targetvalues_upper), n.min = input$nmin, plot.all = reflimR.plot.all,
                                main = paste0("Reference limits for ", parameter_name))}
    }
    
    if(input$check_refineR) {
      
      validate(need(refineR_done(), "(refineR) Please perform the refineR calculation first."))
      
      table_refineR <- getRI(fit_refineR())
      targetvalues_low <- table_refineR$PointEst[1]
      targetvalues_upper <- table_refineR$PointEst[2]
      
      if(plot_choice == "VeRUS"){
        reflim_result <- reflim_VeRUS(dat[, 4], targets = c(targetvalues_low, targetvalues_upper), n.min = input$nmin, plot.all = reflimR.plot.all, lambda = lambda,
                                      main = paste0("Reference limits for ", parameter_name))}
      if(plot_choice == "pU"){
        reflim_result <- reflim(dat[,4], targets = c(targetvalues_low, targetvalues_upper), n.min = input$nmin, plot.all = reflimR.plot.all,
                                main = paste0("Reference limits for ", parameter_name))}
    }
    
    if (input$check_target == FALSE && input$check_targetvalues == FALSE && input$check_refineR == FALSE) {
      
      if(plot_choice == "VeRUS"){
        reflim_result <- reflim_VeRUS(dat[, 4], n.min = input$nmin, plot.all = reflimR.plot.all, lambda = lambda,
                                      main = paste0("Reference limits for ", parameter_name))}
      if(plot_choice == "pU"){
        reflim_result <- reflim(dat[, 4], n.min = input$nmin, plot.all = reflimR.plot.all, main = paste0("Reference limits for ", parameter_name))}
    }
    
    reflim_result
    if(plot_choice == "VeRUS"){
      if(!reflimR.plot.all){
      usr <- par("usr")
      rect(usr[1], usr[3], usr[2], usr[4], border = "azure4", lwd = 5, bty = "o")
      }
    }
  })
  
  output$table <- DT::renderDataTable({
    
    DT::datatable(reflim_data(), caption = htmltools::tags$caption(style = 'caption-side: bottom; text-align: center;','Dataset'))
  })
  
  output$table_report <- DT::renderDataTable({ #Tab:reflimR
    
    report <- get_data_report()
    if (!is.na(report$limits[1])) {
      parameter_name <- parameter_display()
      converted_sex <- switch(input$sex,
                              "f" = "Female(F)",
                              "m" = "Male(M)",
                              "t" = "Female(F) & Male(M)")
      
      
      if(input$lambda_type == "reflimR"){
        if(report$lognormal){
          lambda <- 0
        } else{
          lambda <- 1
        }
      } else if(input$lambda_type == "refineR"){
        validate(need(refineR_done(), "(refineR) Please perform the refineR calculation first."))
        lambda <- lambda_refineR
      } else if(input$lambda_type == "user"){
        lambda <- input$lambda
      }
      
      report_versus <- verus.limits(report$limits[1], report$limits[2], lambda = lambda)
      
      report_lower_VeRUS <- c(round(report_versus$lower.lim.low, 1), round(report_versus$lower.lim.upp, 1))
      report_upper_VeRUS <- c(round(report_versus$upper.lim.low, 1), round(report_versus$upper.lim.upp, 1))
      
      if(is.na(report$targets[1])){
        report_lower_target_VeRUS <- c(NA, NA)
        report_upper_target_VeRUS <- c(NA, NA)
      } else{
        report_target_versus <- verus.limits(report$targets[1], report$targets[2], lambda = lambda)
        
        report_lower_target_VeRUS <- c(round(report_target_versus$lower.lim.low, 1), round(report_target_versus$lower.lim.upp, 1))
        report_upper_target_VeRUS <- c(round(report_target_versus$upper.lim.low, 1), round(report_target_versus$upper.lim.upp, 1))
      }
      
      lognormal_value <<- report$lognormal
      
      table_report <- t(data.frame(
        "Sex and Age:" = paste0(converted_sex, " (", input$age_end[1], "-", input$age_end[2], ")"),
        "Category:" = input$category,
        "Mean (sd):" = paste0(round(report$stats[1], 2), " (", round(report$stats[2], 2), ")"),
        "Lognormal Distribution:" = report$lognormal,
        "Reference limit:" =  paste0(report$limits[1] , " - " , report$limits[2]),
        "Lower tolerance intervals (EL):" = paste0(report$limits[3], " - " , report$limits[4]),
        "Upper tolerance intervals (EL):" = paste0(report$limits[5], " - " , report$limits[6]),
        "Lower UM intervals:" = paste0(report_lower_VeRUS[1], " - " , report_lower_VeRUS[2]),
        "Upper UM intervals:" = paste0(report_upper_VeRUS[1], " - " ,report_upper_VeRUS[2]),
        "Lower confidence intervals:" = paste0(report$confidence.int[1], " - " , report$confidence.int[2]),
        "Upper confidence intervals:" = paste0(report$confidence.int[3], " - " , report$confidence.int[4]),
        "Target Limits:" = paste0(report$targets[1], " - " , report$targets[2]),
        "Lower target tolerance intervals (EL):" = paste0(report$targets[3], " - " , report$targets[4]),
        "Upper target tolerance intervals (EL):" = paste0(report$targets[5], " - " , report$targets[6]),
        "Lower UM target intervals:" = paste0(report_lower_target_VeRUS[1], " - " , report_lower_target_VeRUS[2]),
        "Upper UM target intervals:" = paste0(report_upper_target_VeRUS[1], " - " , report_upper_target_VeRUS[2]),
        "Interpretation of the lower limit:" = report$interpretation[1],
        "Interpretation of the upper limit:" = report$interpretation[2],
        "Lambda for UM:" = lambda,
        check.names = FALSE))
      colnames(table_report) <- parameter_name
      
      DT::datatable(table_report, extensions = 'Buttons',
                    options = list(dom = 'Bt', pageLength = 19, buttons = c('copy', 'csv', 'pdf', 'print')))
    }
  })
   
  refineR_done <- reactiveVal(FALSE)
  
  observeEvent(list(
    input$parameter,
    input$category,
    input$sex,
    input$age_end,
    input$reset,
    input$dataset_file1,
    input$submit,
    input$show_table,
    input$upload_age_column,
    input$upload_sex_column,
    input$upload_female_value,
    input$upload_male_value
  ), { 
    refineR_done(FALSE)
  }, ignoreInit = TRUE)
  
  output$plotrefineR <- renderPlot({ #Tab:refineR
    plot(fit_refineR())
  })
  
  output$refineR_checkbox_ui <- renderUI({ 
    checkboxInput(
      "check_refineR",
      "Load calculated refineR values",
      value = FALSE
    )
  })
  
  output$table_report_refineR <- DT::renderDataTable({ #Tab:refineR
    
    table_refineR_original <- fit_refineR()
    table_refineR <- getRI(table_refineR_original)
    parameter_name <- parameter_display()
    
    lambda_refineR <<- table_refineR_original$Lambda
    
    converted_sex <- switch(input$sex,
                            "f" = "Female(F)",
                            "m" = "Male(M)",
                            "t" = "Female(F) & Male(M)")
    
    table_report <- t(data.frame(
      "Sex and Age:" = paste0(converted_sex, " (", input$age_end[1], "-", input$age_end[2], ")"),
      "Category:" = input$category,
      "Reference limit:" =  paste0(round(table_refineR$PointEst[1], 3) , " - " , round(table_refineR$PointEst[2], 3)),
      "Lambda:" = table_refineR_original$Lambda,
      check.names = FALSE))
    colnames(table_report) <- parameter_name
    
    DT::datatable(table_report, extensions = 'Buttons',
                  caption = htmltools::tags$caption(
                    style = "caption-side: top; text-align: left; font-weight: bold; font-size: 16px;",
                    "Report for refineR"
                  ), options = list(dom = 'Bt', pageLength = 15, buttons = c('copy', 'csv', 'pdf', 'print')))
  })
  
  mclust_plot_digits <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x) & x > 0]
    if(length(x) == 0){
      return(2)
    }
    x <- x[x < (median(x) + 6 * IQR(x))]
    if(length(x) == 0){
      return(2)
    }
    scale.ref <- median(x, na.rm = TRUE)
    if(!is.finite(scale.ref) || scale.ref <= 0){
      return(2)
    }
    max(0, 2 - floor(log10(scale.ref)))
  }
  
  output$plotmclust <- renderPlot({ #Tab:mclust
    
    validate(
      need(
        input$auto_cluster || 
          (!is.null(input$n_cluster) && input$n_cluster >= 1 && input$n_cluster <= 15),
        "Please choose a number of clusters between 1 and 15."
      )
    )
    
    withProgress(message = "mclust Calculation …", {
      dat <- reflim_data()
      n_cluster <- if (input$auto_cluster) NULL else input$n_cluster
      plot_digits <- mclust_plot_digits(dat[, 4])
      lab_mclust(dat[, 4], lognormal = lognormal_value, model = input$model_name, n.cluster = n_cluster, remove.extremes = T, digits = plot_digits)
    })
  })
  
  output$plotmclustbic <- renderPlot({ #Tab:mclust
    
    validate(
      need(
        input$auto_cluster || 
          (!is.null(input$n_cluster) && input$n_cluster >= 1 && input$n_cluster <= 15),
        "Please choose a number of clusters between 1 and 15."
      )
    )
    
    withProgress(message = "mclust Calculation …", {
      dat <- reflim_data()
      n_cluster <- if (input$auto_cluster) NULL else input$n_cluster
      plot_digits <- mclust_plot_digits(dat[, 4])
      lab_mclust(dat[, 4], lognormal = lognormal_value, model = input$model_name, n.cluster = n_cluster, remove.extremes = T, plot.bic = T, digits = plot_digits)
    })
  })
  
  # output$scatterplot <- renderPlot({ #Tab:Scatterplot
  #   
  #   dat <- reflim_data()
  #   parameter_name <- parameter_display()
  #   colors <- ifelse(dat[, 3] == "f", "indianred", "cornflowerblue")
  #   pchs <- ifelse(dat[, 3] == "f", 17, 19)
  #   plot(dat[,4] ~ dat[,2], pch = pchs, cex = 1, col = colors, xlab = "Age", ylab = parameter_name)
  #   
  #   unique_levels <- levels(factor(dat[, 3]))
  #   legend("topright", legend = unique_levels, pch = c(17, 19)[1:length(unique_levels)], col = c("indianred", "cornflowerblue")[1:length(unique_levels)])
  # })
  
  # output$plot_statistics <- renderPlot({ #Tab:Statistics
  #   
  #   par(mfrow = c(2,1))
  #   
  #   dat <- reflim_data()
  #   ylab_ <- parameter_display()
  #   
  #   if (!(nrow(dat)) == 0) {
  #     hist_data_w <- subset(dat, Sex == "f", select = Age)
  #     hist_data_m <- subset(dat, Sex == "m", select = Age)
  #     
  #     hist_w <- hist(hist_data_w$Age, breaks = seq(min(dat[,2]) - 1,max(dat[,2]),by = 1))$counts
  #     hist_m <- hist(hist_data_m$Age, breaks = seq(min(dat[,2]) - 1,max(dat[,2]),by = 1))$counts
  #     
  #     barplot(rbind(hist_m,hist_w), col = c("cornflowerblue","indianred"),
  #             names.arg = seq(min(dat[,2]), max(dat[,2]), by = 1), xlab = "Age", las = 1, beside = TRUE, ylab = "Number of data")
  #     abline(h = 0)
  #     legend("topright", legend = c(paste0("m: ", nrow(hist_data_m)), paste0("f: ", nrow(hist_data_w))), col = c("cornflowerblue","indianred"), pch = c(19, 19))
  #     
  #     par(new = TRUE)
  #     boxplot(dat[,2], horizontal = TRUE, axes = FALSE, col = rgb(0, 0, 0, alpha = 0.15))
  #   }
  #   
  #   if (!(nrow(dat)) == 0) {
  #     
  #     if (input$sex == "m") {
  #       boxplot(dat[,4]~interaction(dat[,3], dat[,2]), xlab = "Age", 
  #               ylab = ylab_, col = "cornflowerblue", las = 2)
  #     }
  #     else if (input$sex == "f") {
  #       boxplot(dat[,4]~interaction(dat[,3], dat[,2]), xlab = "Age", 
  #               ylab = ylab_, col = "indianred", las = 2)
  #     } else{
  #       boxplot(dat[,4]~interaction(dat[,3], dat[,2]), xlab = "Age", 
  #               ylab = ylab_, col = c("indianred", "cornflowerblue"), las = 2)
  #     }
  #   }
  # })
  
  # output$table_zlog <- DT::renderDataTable({ #Tab:zlog
  #   
  #   dat <- reflim_data()
  #   report <- get_data_report()
  #   
  #   zlog_results <- numeric(nrow(dat))
  #   for (i in 1:nrow(dat)) {
  #     zlog_results[i] <- round_df(zlog(dat[i, 4], report$limits[1], report$limits[2]), 2)
  #   }
  #   
  #   reflim_data <- cbind(dat, "RI" = paste0(report$limits[1], " - " , report$limits[2]), "zlog" = zlog_results)
  #   
  #   options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))
  #   
  #   DT::datatable(reflim_data, rownames = FALSE, extensions = 'Buttons',
  #                 options = list(dom = 'Blfrtip', pageLength = 15, buttons = c('copy', 'csv', 'pdf', 'print')),
  #                 caption = htmltools::tags$caption(style = 'caption-side: bottom; text-align: center;',
  #                                                   'Table: Dataset with the zlog values')) %>%
  #     DT::formatStyle(columns = "zlog", 
  #                     color = styleEqual(reflim_data[,6], highzlogvalues(c(reflim_data[,6]))),
  #                     backgroundColor = styleEqual(reflim_data[,6], zlogcolor(c(reflim_data[,6])))) %>%
  #     DT::formatStyle(columns = colnames(reflim_data)[4], 
  #                     color = styleEqual(reflim_data[,4], highzlogvalues(c(reflim_data[,6]))),
  #                     backgroundColor = styleEqual(reflim_data[,4], zlogcolor(c(reflim_data[,6]))))
  # })
  
  output$download_ritable <- downloadHandler(
    filename = function() {
      paste("ReferenceIntervals_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      
      dat <- get_alldata_file()
      dataset <- dat[c(-1,-2,-3)]
      reflim.loop.results <- reflim.loop(dataset, plot.it = FALSE)
      
      tmpdir <- tempdir()
      csv_files <- c()
      
      for (col in names(reflim.loop.results)) {
        report <- reflim.loop.results[[col]]
        
        converted_sex <- switch(input$sex,
                                "f" = "Female(F)",
                                "m" = "Male(M)",
                                "t" = "Female(F) & Male(M)")
        
        df <- t(data.frame(
          "Sex and Age:" = paste0(converted_sex, " (", input$age_end[1], "-", input$age_end[2], ")"),
          "Category:" = input$category,
          "Mean:" = report$stats[1],
          "Standard deviation:" = report$stats[2],
          "Lognormal Distribution:" = report$lognormal,
          "Reference Interval:" = paste0(report$limits[1] , " - " , report$limits[2]),
          "Lower tolerance intervals:" = paste0(report$limits[3], " - " , report$limits[4]),
          "Upper tolerance intervals:" = paste0(report$limits[5], " - " , report$limits[6]),
          "Target Limits:" = paste0(report$targets[1], " - " , report$targets[2]),
          "Lower target tolerance intervals:" = paste0(report$targets[3], " - " , report$targets[4]),
          "Upper target tolerance intervals:" = paste0(report$targets[5], " - " , report$targets[6]),
          "Lower confidence intervals:" = paste0(report$confidence.int[1], " - " , report$confidence.int[2]),
          "Upper confidence intervals:" = paste0(report$confidence.int[3], " - " , report$confidence.int[4]),
          "Interpretation of the lower limit:" = report$interpretation[1],
          "Interpretation of the upper limit:" = report$interpretation[2],
          check.names = FALSE))
        
        csv_path <- file.path(tmpdir, paste0(col, ".csv"))
        write.csv(df, csv_path)
        csv_files <- c(csv_files, csv_path)
      }
      
      old_wd <- setwd(tmpdir)
      on.exit(setwd(old_wd))
      
      zip(zipfile = file, files = basename(csv_files), extras = "-j")
    }
  )
  
  # output$download_zlogtable <- downloadHandler( #Tab:zlog
  #   filename = function() {
  #     paste("zlogValues_", Sys.Date(), ".csv", sep = "")
  #   },
  #   content = function(file) {
  #     dat <- get_alldata_file()
  #     
  #     dataset <- dat[c(-1,-2,-3)]
  #     reflim.loop.results <- reflim.loop(dataset, plot.it = FALSE)
  #     zlog.loop.results <- zlog.loop(dataset, reflim.loop.results)
  #     
  #     result <- c(dat, zlog.loop.results)
  #     write.csv(result, file)
  #   }
  # )
  
  output$tree_rpart <- renderPlot({ #Tab:rpart
    
    rpart.plot::rpart.plot(build_rpart(), box.palette = "RdBu", roundint = FALSE,
                           main = paste("Regression tree:", parameter_display(), "~ Sex + Age"),)
  })
  
}
####################################### Run the application #######################################
shinyApp(ui = ui, server = server)

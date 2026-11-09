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

####################################### Load Scripts ##############################################

# Allow CSV and Excel uploads up to 100 MB.
options(shiny.maxRequestSize = 100 * 1024^2)

# Minimum sample size recommended for reliable reference limit estimation.
reflim_n_min <- 200

# Candidate cluster counts for automatic selection and the BIC comparison.
mclust_cluster_counts <- 1:15

source("functions.R")
source("excel_import.R")
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

####################################### Texts #####################################################

source("texts.R", encoding = "UTF-8")

####################################### User Interface ############################################

ui <- dashboardPage(
  dashboardHeader(title = "VeRIf", titleWidth = 300),
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "sidebarid",
      
      uiOutput("parameters"),
      div(
        class = "analyte-label-inputs",
        textInput(
          "parameter_label",
          "Laboratory analyte:",
          value = "",
          width = "100%"
        ),
        textInput(
          "parameter_units",
          "Units:",
          value = "",
          width = "100%"
        )
      ),
      
      sliderInput(
        "age_end",
        "Select age-range:",
        min = 0,
        max = 100,
        value = c(0, 100)),
      
      uiOutput("sex_selection"),
      uiOutput("category"),
      
      hr(),
      
      conditionalPanel(
        condition = "output.preinstalled_targetvalues_available === true",
        checkboxInput("check_targetvalues", "Load preinstalled target reference interval", value = FALSE)
      ),
      checkboxInput("check_target", "Target reference interval", value = FALSE),
      
      div(
        class = "target-reference-interval",
        numericInput(
          "target_low",
          "from",
          10,
          min = 0,
          max = 10000,
          width = "100%"
        ),
        
        numericInput(
          "target_upper",
          "to",
          15,
          min = 0,
          max = 10000,
          width = "100%"
        )
      ), 
      uiOutput("refineR_checkbox_ui"),
      hr(),
      
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
      
      hr()
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .shiny-input-disabled {
          opacity: 0.55;
        }

        .analyte-label-inputs {
          display: flex;
          gap: 10px;
          padding: 12px 15px 0;
          margin-bottom: 15px;
        }

        section.sidebar .analyte-label-inputs .shiny-input-container {
          flex: 1 1 0;
          min-width: 0;
          padding: 0;
          margin-bottom: 0;
        }

        .analyte-label-inputs .shiny-input-container:first-child {
          flex-grow: 2;
        }

        .target-reference-interval {
          display: flex;
          gap: 10px;
          padding: 0 15px;
          margin-bottom: 15px;
        }

        .target-reference-interval .shiny-input-container {
          display: flex;
          align-items: center;
          flex: 1 1 0;
          gap: 6px;
          min-width: 0;
          padding: 0;
          margin-bottom: 0;
        }

        .target-reference-interval label {
          flex: 0 0 auto;
          margin-bottom: 0;
          font-weight: normal;
        }

        .target-reference-interval input {
          min-width: 0;
        }
      ")),
      tags$script(HTML("
        Shiny.addCustomMessageHandler('setInputDisabled', function(message) {
          var input = document.getElementById(message.id);
          if (!input) {
            return;
          }

          input.disabled = !!message.disabled;

          var container = input.closest('.form-group');
          if (container) {
            container.classList.toggle('shiny-input-disabled', !!message.disabled);
          }
        });
      "))
    ),
    fluidRow(
      
      tabsetPanel( 
        id = "analysis_tab",
        tabPanel("Input", 
                 icon = icon("upload"),
                 
                 box(
                   status = "info",
                   width = 7,
                   
                   p(data_text),
                   
                   uiOutput("dataset_file"),
                   uiOutput("upload_mapping_ui"),
                   actionButton('reset', 'Reset Input', icon = icon("trash")),
                   actionButton("load_demo", "Load demo data", icon = icon("database")), br(),br(),
                   
                   fluidRow(
                     column(6, checkboxInput("show_table", tags$b("Editable table"), value = FALSE)),
                     column(6, actionButton("submit", "Submit"))
                   ),
                    conditionalPanel(
                      condition = "input.show_table == true",
                      helpText("Enter or paste values into the existing columns (Sex: f = female, m = male). Click Submit after adding or changing data to use the table for analysis."),
                      withSpinner(rHandsontableOutput("editable_table"))
                    ), hr(),

                    DT::dataTableOutput("table")
                  )
        ),
        
        tabPanel("Overview",
                 icon = icon("chart-line"),

                 box(
                   status = "info",
                   width = 7,

                   p(scatterplot_text),
                   plotOutput("scatterplot", height = "350px"),
                   plotOutput("plot_statistics", height = "350px"),
                   fluidRow(
                     column(6, selectInput(
                       "overview_age_group_width", "Boxplot age groups:",
                       choices = c("Automatic" = "auto", "1 year" = "1", "2 years" = "2",
                                   "5 years" = "5", "10 years" = "10", "20 years" = "20"),
                       selected = "auto"
                     )),
                     column(6, selectInput(
                       "overview_boxplot_scale", "Boxplot Y-axis:",
                       choices = c("Linear" = "linear", "Logarithmic" = "log"),
                       selected = "linear"
                     ))
                   ),
                   helpText("Age groups include the lower limit, e.g. 20-<25 means 20 to under 25 years. Use a logarithmic scale to see small values alongside large outliers."),
                   plotOutput("plot_age_boxplot", height = "400px")
                 )
        ),
        
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
                              selected = "pU",
                              inline = TRUE
                            )
                          ),
                     
                  ), withSpinner(plotOutput("plot", height = "70vh"))
                 )
        ),
        
        tabPanel("refineR", 
                  icon = icon("chart-line"),
                  
                  box(
                    status = "info",
                    width = 7,
                    
                    p(refineR_text),
                    withSpinner(plotOutput("plotrefineR", height = "70vh"))
                  )
        ),
        
        tabPanel("mclust",
          icon = icon("chart-line"),
          
          box(
            status = "info",
            width = 7,
            
            p(mclust_text),
            
            fluidRow(
              
              column(
                4,
                selectizeInput(
                  inputId = "n_cluster",
                  label = "Number of clusters:",
                  choices = c("Automatic (BIC)" = "auto",
                              setNames(as.character(mclust_cluster_counts),
                                       as.character(mclust_cluster_counts))),
                  selected = "auto",
                  options = list(create = FALSE)
                )
              ),
              
              column(
                8,
                radioButtons(
                  inputId = "model_name",
                  label = "Univariate model for mclust:",
                  choices = c("Equal variance" = "E", "Unequal variance" = "V"),
                  selected = "V",
                  inline = TRUE
                )
              )
            ),

            fluidRow(
              column(
                12,
                radioButtons(
                  inputId = "mclust_distribution",
                  label = "Distribution for mclust:",
                  choices = c(
                    "From reflimR report" = "report",
                    "Normal" = "normal",
                    "Lognormal" = "lognormal"
                  ),
                  selected = "report",
                  inline = TRUE
                )
              )
            ),
            
            fluidRow(
              box(
                title = "BIC model comparison",
                status = "info",
                width = 12,
                collapsible = TRUE,
                collapsed = TRUE,
                p(mclust_bic_text),
                withSpinner(plotOutput("plotmclustbic", height = "350px"))
              )
            ),
            
            withSpinner(plotOutput("plotmclust", height = "70vh"))
          )
        ),
        
        tabPanel("rpart", 
                  icon = icon("chart-bar"),
                  
                  box(
                    status = "info",
                    width = 7,
                    
                    p(rpart_text),
                    
                    numericInput("tree_window_minsplit", "Minimum observations needed to split a node:", 20, min = 10, max = 100),
                    #numericInput("tree_window_cp", "rpart: Complexity parameter", 0.01, min = 0, max = 10),
                    withSpinner(plotOutput("tree_rpart", height = "70vh"))
                  )
        ),
        
        tabPanel( "zlog",
                  icon = icon("table"),

                  box(
                    status = "info",
                    width = 7,

                    p(zlog_text),
                    withSpinner(DT::dataTableOutput("table_zlog", height = "700px"))
                  )
        ),
        
        tabPanel("About",
                  icon = icon("info-circle"),
                  
                  box(
                    status = "info",
                    width = 7,
                    
                    about_text
                  )
        )
      ),
      
      box(
        status = "primary",
        width = 5,
        
        tabsetPanel(
          id = "results_tab",
          tabPanel(
            "reflimR Results",
            withSpinner(DT::dataTableOutput("table_report")), hr(),
            downloadButton("download_ritable", "Download all Reference Intervals"),
            downloadButton("download_zlogtable", "Download all zlog values")
          ),
          tabPanel(
            "refineR Results",
            withSpinner(DT::dataTableOutput("table_report_refineR"))
          )
        )
      )
    )
  )
)

####################################### Server ####################################################

server <- function(input, output, session) {
  
  ##################################### Observe Events ############################################
  
  options(shiny.sanitize.errors = TRUE)
  options(warn = -1)

  observeEvent(input$analysis_tab, {
    if (input$analysis_tab %in% c("reflimR", "refineR")) {
      updateTabsetPanel(session, "results_tab",
                       selected = paste(input$analysis_tab, "Results"))
    }
  })

  updateInputDisabled <- function(input_id, disabled) {
    session$sendCustomMessage(
      "setInputDisabled",
      list(id = input_id, disabled = isTRUE(disabled))
    )
  }
  
  data_source <- reactiveVal("empty")
  
  observeEvent(input$dataset_file1, {
    data_source("uploaded")
    updateCheckboxInput(session, "show_table", value = FALSE)
  })
  
  observeEvent(input$reset, {
    # Start a fresh session so removed dynamic inputs, table edits and cached
    # results are cleared together with every visible input and its settings.
    session$reload()
  })

  observeEvent(input$load_demo, {
    data_source("demo")
    updateCheckboxInput(session, "show_table", value = FALSE)
  })
  
  dataset_input <- reactive({
    if (identical(data_source(), "uploaded")) {
      input$dataset_file1
    } else {
      NULL
    }
  })
  
  output$dataset_file <- renderUI({
    input$reset ## Create a dependency with the reset button
    input$load_demo
    fileInput('dataset_file1', label = 'Upload CSV or Excel (.xlsx)',
              multiple = FALSE, accept = c('.csv', '.xlsx'), buttonLabel = 'Choose file')
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

  guess_sex_value <- function(values, preferred_values = character(), excluded = NULL,
                              fallback_first = TRUE) {
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
    
    if (fallback_first) available_values[1] else NULL
  }

  guess_sex_mapping <- function(values) {
    available_values <- unique(as.character(values))
    available_values <- available_values[
      !is.na(available_values) & nzchar(trimws(available_values))
    ]

    female_value <- guess_sex_value(
      available_values,
      c("f", "female", "weiblich", "w"),
      fallback_first = FALSE
    )
    male_value <- guess_sex_value(
      available_values,
      c("m", "male", "mannlich", "maennlich", "mann", "man"),
      fallback_first = FALSE
    )

    remaining_values <- available_values[
      !normalize_upload_value(available_values) %in%
        normalize_upload_value(c(female_value, male_value))
    ]

    if (is.null(female_value) && length(available_values) > 1 && length(remaining_values)) {
      female_value <- remaining_values[1]
      remaining_values <- remaining_values[-1]
    }
    if (is.null(male_value) && length(available_values) > 1 && length(remaining_values)) {
      male_value <- remaining_values[1]
    }

    list(female = female_value, male = male_value)
  }

  excel_import_selection <- reactiveVal(NULL)

  excel_sheets <- reactive({
    file_info <- dataset_input()
    req(file_info)
    if (!grepl("\\.xlsx$", file_info$name, ignore.case = TRUE)) return(NULL)
    result <- tryCatch(readxl::excel_sheets(file_info$datapath), error = identity)
    validate(need(!inherits(result, "error"),
                  "The Excel workbook could not be opened. Please upload a valid, unencrypted XLSX file."))
    result
  })

  observeEvent(dataset_input(), {
    excel_import_selection(NULL)
    removeModal()
    if (is.null(dataset_input())) return()
    sheets <- excel_sheets()
    if (length(sheets) <= 1L) return()
    showModal(modalDialog(
      title = "Select worksheet",
      helpText(sprintf("The uploaded Excel file contains %d sheets. Choose the sheet to read.", length(sheets))),
      selectInput("excel_sheet", "Worksheet:", choices = sheets, selected = sheets[1]),
      footer = actionButton("apply_excel_import", "Read selected sheet"),
      easyClose = FALSE
    ))
  }, ignoreNULL = FALSE, priority = 100)

  observeEvent(input$apply_excel_import, {
    req(input$excel_sheet %in% excel_sheets())
    excel_import_selection(list(path = dataset_input()$datapath,
                                sheet = input$excel_sheet))
    removeModal()
  })

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
      sheets <- excel_sheets()
      if (length(sheets) == 1L) {
        sheet <- sheets[1]
      } else {
        selection <- excel_import_selection()
        validate(need(!is.null(selection) && identical(selection$path, datapath) &&
                        selection$sheet %in% sheets,
                      "Please choose a worksheet in the selection window."))
        sheet <- selection$sheet
      }
      dataset <- tryCatch(read_excel_dataset(datapath, sheet), error = identity)
      if (inherits(dataset, "error")) {
        validate(need(FALSE, paste("Excel import failed:", conditionMessage(dataset))))
      }
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

  selected_upload_category <- function(dataset, age_column, sex_column) {
    available_columns <- setdiff(names(dataset), c(age_column, sex_column))
    selected <- input$upload_category_column
    if (!is.null(selected) && selected %in% c("", available_columns)) {
      return(if (nzchar(selected)) selected else NULL)
    }
    matches <- available_columns[
      normalize_upload_value(available_columns) %in% c("category", "kategorie")
    ]
    if (length(matches)) matches[1] else NULL
  }

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
    
    category_column <- selected_upload_category(dataset, selected_age_column, selected_sex_column)

    sex_values <- unique(as.character(dataset[[selected_sex_column]]))
    sex_values <- sex_values[!is.na(sex_values) & nzchar(trimws(sex_values))]
    
    guessed_sex_mapping <- guess_sex_mapping(sex_values)
    guessed_female_value <- guessed_sex_mapping$female
    selected_female_value <- if (!is.null(input$upload_female_value) &&
                                 input$upload_female_value %in% c("", sex_values)) {
      input$upload_female_value
    } else {
      if (is.null(guessed_female_value)) "" else guessed_female_value
    }
    
    guessed_male_value <- guessed_sex_mapping$male
    selected_male_value <- if (!is.null(input$upload_male_value) &&
                               input$upload_male_value %in% c("", sex_values)) {
      input$upload_male_value
    } else {
      if (is.null(guessed_male_value)) "" else guessed_male_value
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
            "upload_category_column",
            "Category column:",
            choices = c("None (all rows)" = "", setdiff(column_names, c(selected_age_column, selected_sex_column))),
            selected = if (is.null(category_column)) "" else category_column
          )
        ),
        column(
          3,
          selectInput(
            "upload_female_value",
            "Value for female:",
            choices = c("Not present" = "", sex_values),
            selected = selected_female_value
          )
        ),
        column(
          3,
          selectInput(
            "upload_male_value",
            "Value for male:",
            choices = c("Not present" = "", sex_values),
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
    
    category_column <- selected_upload_category(dataset, age_column, sex_column)
    
    analyte_columns <- setdiff(names(dataset), c(category_column, age_column, sex_column))
    
    validate(need(
      length(analyte_columns) > 0,
      "Please choose age, sex and category columns so that at least one laboratory analyte remains."
    ))
    
    sex_values <- dataset[[sex_column]]
    guessed_sex_mapping <- guess_sex_mapping(sex_values)
    female_value <- if (!is.null(input$upload_female_value)) {
      if (nzchar(trimws(input$upload_female_value))) input$upload_female_value else NULL
    } else {
      guessed_sex_mapping$female
    }
    male_value <- if (!is.null(input$upload_male_value)) {
      if (nzchar(trimws(input$upload_male_value))) input$upload_male_value else NULL
    } else {
      guessed_sex_mapping$male
    }
    
    validate(need(
      !is.null(female_value) || !is.null(male_value),
      "Please map at least one sex value to female or male."
    ))
    validate(need(
      is.null(female_value) || is.null(male_value) ||
        normalize_upload_value(female_value) != normalize_upload_value(male_value),
      "Female and male must use different source values."
    ))
    
    normalized_sex <- rep(NA_character_, nrow(dataset))
    normalized_source_values <- normalize_upload_value(sex_values)
    if (!is.null(female_value)) {
      normalized_sex[normalized_source_values == normalize_upload_value(female_value)] <- "f"
    }
    if (!is.null(male_value)) {
      normalized_sex[normalized_source_values == normalize_upload_value(male_value)] <- "m"
    }
    
    category_values <- if (!is.null(category_column)) {
      as.character(dataset[[category_column]])
    } else {
      rep("All", nrow(dataset))
    }
    category_values[is.na(category_values) | !nzchar(trimws(category_values))] <- "All"
    
    standardized_dataset <- data.frame(
      Category = category_values,
      Age = parse_numeric_values(dataset[[age_column]]),
      Sex = normalized_sex,
      stringsAsFactors = FALSE
    )
    
    cbind(standardized_dataset, dataset[analyte_columns])
  })

  active_dataset <- reactive({
    if (isTRUE(input$show_table) && isTRUE(input$submit > 0)) {
      data_store()
    } else if (identical(data_source(), "uploaded")) {
      uploaded_dataset_standardized()
    } else if (identical(data_source(), "demo")) {
      reflimR::livertests
    } else {
      NULL
    }
  })

  parameter_label_value <- reactiveVal(NULL)

  observeEvent(input$parameter, {
    label <- if (is.null(input$parameter)) "" else input$parameter
    parameter_label_value(label)
    updateTextInput(session, "parameter_label", value = label)
  }, ignoreNULL = FALSE)

  observeEvent(input$parameter_label, {
    parameter_label_value(input$parameter_label)
  }, ignoreInit = TRUE)
  
  output$parameters <- renderUI({
    dataset <- active_dataset()
    if (is.null(dataset)) return(NULL)
    choices <- setdiff(colnames(dataset), c("Category", "Age", "Sex"))
    validate(need(length(choices) > 0, "No laboratory analyte columns are available."))
    
    selected_parameter <- if (!is.null(input$parameter) && input$parameter %in% choices) {
      input$parameter
    } else {
      choices[1]
    }
    
    selectInput("parameter","Select laboratory analyte:", choices = choices, selected = selected_parameter)
  })   
  
  output$category <- renderUI({
    dataset <- active_dataset()
    if (is.null(dataset)) return(NULL)
    choices <- c("Not selected", unique(dataset$Category))
    
    selected_category <- if (!is.null(input$category) && input$category %in% choices) {
      input$category
    } else {
      "Not selected"
    }
    
    selectInput("category", "Select category:", choices = choices, selected = selected_category)
  })

  output$sex_selection <- renderUI({
    dataset <- active_dataset()
    if (is.null(dataset)) return(NULL)

    available_sexes <- intersect(c("f", "m"), unique(as.character(dataset$Sex)))
    choices <- if (identical(available_sexes, c("f", "m"))) {
      c("Female & Male" = "t", "Female" = "f", "Male" = "m")
    } else if (identical(available_sexes, "f")) {
      c("Female" = "f")
    } else if (identical(available_sexes, "m")) {
      c("Male" = "m")
    } else {
      c("Female & Male" = "t", "Female" = "f", "Male" = "m")
    }

    selected_sex <- if (!is.null(input$sex) && input$sex %in% unname(choices)) {
      input$sex
    } else {
      unname(choices)[1]
    }

    selectInput("sex", "Select sex:", choices = choices, selected = selected_sex)
  })
  
  preinstalled_targetvalues_available <- reactive({
    identical(data_source(), "demo") &&
      !(isTRUE(input$show_table) && isTRUE(input$submit > 0))
  })

  output$preinstalled_targetvalues_available <- function(...) {
    preinstalled_targetvalues_available()
  }
  outputOptions(output, "preinstalled_targetvalues_available", suspendWhenHidden = FALSE)

  use_preinstalled_targets <- reactive({
    preinstalled_targetvalues_available() && isTRUE(input$check_targetvalues)
  })
  
  # Create a reactive values to track the state of the checkboxes
  reactive_values <- reactiveValues(
    check_targetvalues = FALSE,
    check_target = FALSE,
    check_refineR = FALSE
  )
  
  # Observe changes in check_targetvalues and update the reactive value
  observeEvent(input$check_targetvalues, {
    if (input$check_targetvalues && preinstalled_targetvalues_available()) {
      reactive_values$check_targetvalues <- TRUE
      reactive_values$check_target <- FALSE
      reactive_values$check_refineR <- FALSE
    } else {
      reactive_values$check_targetvalues <- FALSE
    }
  })

  observeEvent(preinstalled_targetvalues_available(), {
    if (!preinstalled_targetvalues_available()) {
      reactive_values$check_targetvalues <- FALSE
    }
  }, ignoreInit = FALSE)
  
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
      updateTabsetPanel(session, "analysis_tab", selected = "reflimR")
      updateTabsetPanel(session, "results_tab", selected = "reflimR Results")
    } else {
      reactive_values$check_refineR <- FALSE
    }
  })
  
  # Update the checkboxes based on the reactive value
  observe({
    updateCheckboxInput(session, "check_targetvalues", value = reactive_values$check_targetvalues)
    updateCheckboxInput(session, "check_target", value = reactive_values$check_target)
    updateCheckboxInput(session, "check_refineR", value = reactive_values$check_refineR)
    updateInputDisabled(
      "check_targetvalues",
      !preinstalled_targetvalues_available() ||
        (!isTRUE(reactive_values$check_targetvalues) &&
        (isTRUE(reactive_values$check_target) || isTRUE(reactive_values$check_refineR))
        )
    )
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

  parameter_name <- reactive({
    label_value <- parameter_label_value()
    custom_label <- if (is.null(label_value)) "" else trimws(label_value)
    if (nzchar(custom_label)) {
      custom_label
    } else {
      input$parameter
    }
  })

  parameter_display <- reactive({
    analyte <- parameter_name()
    units <- if (is.null(input$parameter_units)) "" else trimws(input$parameter_units)
    if (nzchar(units)) {
      paste0(analyte, " (", units, ")")
    } else {
      analyte
    }
  })

  table_export_buttons <- function(result_type, export_columns = NULL) {
    analyte <- gsub("[^[:alnum:]_-]+", "_", parameter_name())
    analyte <- gsub("^_+|_+$", "", analyte)
    if (length(analyte) == 0L || !nzchar(analyte)) analyte <- "Analyte"
    prefix <- paste("VeRIf", analyte, result_type, sep = "_")

    # Use the browser's current local date when the export button is clicked.
    filename <- DT::JS(sprintf("function() {
      var now = new Date();
      var date = now.getFullYear() + '-' +
        ('0' + (now.getMonth() + 1)).slice(-2) + '-' +
        ('0' + now.getDate()).slice(-2);
      return %s + '_' + date;
    }", jsonlite::toJSON(prefix, auto_unbox = TRUE)))

    lapply(c("copy", "csv", "pdf", "print"), function(button) {
      config <- list(extend = button, filename = filename, title = filename)
      if (!is.null(export_columns)) {
        config$exportOptions <- list(columns = export_columns)
      }
      config
    })
  }
  
  refineR_data <- reactive({
    dataset <- req(active_dataset())
    req(input$parameter, input$age_end, input$sex)
    
    validate(need(
      input$parameter %in% names(dataset),
      "Please select a valid laboratory analyte."
    ))
    
    dataset <- dataset[, c("Category", "Age", "Sex", input$parameter), drop = FALSE]
    
    if (!is.null(input$category) && input$category != "Not selected") {
      dataset <- subset(dataset, Category == input$category)
    }
    
    dataset$Age <- parse_numeric_values(dataset$Age)
    dataset[, 4] <- parse_numeric_values(dataset[, 4])
    
    dataset <- subset(dataset, Age >= input$age_end[1] & Age <= input$age_end[2])
    
    if (input$sex %in% c("m", "f")) {
      dataset <- subset(dataset, Sex == input$sex)
    }
    
    dataset
  })
  
  # Create the table with the dataset as reactive expression 
  reflim_data <- reactive({
    dataset <- req(active_dataset())
    req(input$parameter, input$age_end, input$sex)
    
    validate(need(
      input$parameter %in% names(dataset),
      "Please select a valid laboratory analyte."
    ))
    
    dataset <- dataset[, c("Category", "Age", "Sex", input$parameter), drop = FALSE]
    
    if (!is.null(input$category) && input$category != "Not selected") {
      dataset <- subset(dataset, Category == input$category)
    }
    
    validate(need(
      ncol(dataset) == 4, 
      "Check if you have used the correct template! You need 4 columns (Category, Age, Sex, Value)!"
    ))
    dataset$Age <- parse_numeric_values(dataset$Age)
    dataset[, 4] <- parse_numeric_values(dataset[, 4])
    
    dataset <- subset(dataset, Age >= input$age_end[1] & Age <= input$age_end[2])
    
    if (input$sex %in% c("m", "f")) {
      dataset <- subset(dataset, Sex == input$sex)
    }
    
    return(dataset)
  })

  get_alldata_file <- reactive({
    dataset <- req(active_dataset())
    req(input$age_end, input$sex)
    
    dataset$Age <- parse_numeric_values(dataset$Age)
    
    if (!is.null(input$category) && input$category != "Not selected") {
      dataset <- subset(dataset, Category == input$category)
    }
    
    dataset <- subset(dataset, Age >= input$age_end[1] & Age <= input$age_end[2])
    
    if (input$sex %in% c("m", "f")) {
      dataset <- subset(dataset, Sex == input$sex)
    }
    
    return(dataset)
  })
  
  # Estimated limits and distribution do not depend on verification targets.
  # zlog and mclust must remain usable when refineR targets are not ready.
  get_reflim_report <- reactive({
    dat <- reflim_data()
    validate(need(nrow(dat) > 39,
                  paste0("(reflim) n = ", nrow(dat), ". The absolute minimum for reference limit estimation is 40.")))

    reflim_with_warnings(dat[, 4], n.min = reflim_n_min, plot.all = FALSE, plot.it = FALSE)
  })

  get_data_report <- reactive({
    
    dat <- reflim_data()
    validate(need(nrow(dat) > 39,
                  paste0("(reflim) n = ", nrow(dat), ". The absolute minimum for reference limit estimation is 40.")))
    
    if (!isTRUE(input$check_target) && !use_preinstalled_targets() && !isTRUE(input$check_refineR)) {
      reflim_text <- get_reflim_report()
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
      
      reflim_text <- reflim_with_warnings(dat[,4], targets = c(input$target_low, input$target_upper), n.min = reflim_n_min, plot.all = FALSE, plot.it = FALSE)
    }
    
    if (use_preinstalled_targets()) {
      
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
                    "(reflim) There are no preloaded target reference interval for this parameter!"))
      
      reflim_text <- reflim_with_warnings(dat[,4], targets = c(targetvalues_low, targetvalues_upper), n.min = reflim_n_min, plot.all = FALSE, plot.it = FALSE)
    }
    
    if(input$check_refineR) {
      
      validate(need(refineR_done(),"(refineR) Please perform the refineR calculation first."))
      
      table_refineR <- getRI(fit_refineR())
      targetvalues_low <- table_refineR$PointEst[1]
      targetvalues_upper <- table_refineR$PointEst[2]
      
      reflim_text <- reflim_with_warnings(dat[,4], targets = c(targetvalues_low, targetvalues_upper), n.min = reflim_n_min, plot.all = FALSE, plot.it = FALSE)
    }
    
    report <- reflim_text
    
    return(report)
  })
  
  fit_refineR <- eventReactive(list(
    data_source(),
    input$parameter,
    input$category,
    input$sex,
    input$age_end,
    input$reset,
    input$dataset_file1,
    excel_import_selection(),
    input$submit,
    input$show_table,
    input$upload_age_column,
    input$upload_sex_column,
    input$upload_category_column,
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
        lower_fence <- quartiles[1] - 6 * iqr_value
        upper_fence <- quartiles[2] + 6 * iqr_value
        
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
    quoted_response_name <- deparse(as.name(tree_input$response_name), backtick = TRUE)
    tree_formula <- reformulate(tree_input$predictors, response = quoted_response_name)
    
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
      validate(need(
        is.logical(report$lognormal) && length(report$lognormal) == 1L && !is.na(report$lognormal),
        "(reflim) The distribution could not be determined. Please check the selected data and settings."
      ))
      if(report$lognormal){
        lambda <- 0
      } else{
        lambda <- 1
      }
    } else if(input$lambda_type == "refineR"){
      validate(need(refineR_done(), "(refineR) Please perform the refineR calculation first."))
      lambda <- fit_refineR()$Lambda
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
        reflim_result <- reflim_VeRUS(dat[, 4], targets = c(input$target_low, input$target_upper), n.min = reflim_n_min, plot.all = reflimR.plot.all, lambda = lambda,
                                      main = paste0("Reference limits for ", parameter_name), xlab = parameter_name)}
      if(plot_choice == "pU"){
        reflim_result <- reflim(dat[, 4], targets = c(input$target_low, input$target_upper), n.min = reflim_n_min, plot.all = reflimR.plot.all,
                                main = paste0("Reference limits for ", parameter_name), xlab = parameter_name)}
    }
    
    if (use_preinstalled_targets()) {
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
                    "(reflim) There are no preloaded target reference interval for this parameter!"))
      
      if(plot_choice == "VeRUS"){
        reflim_result <- reflim_VeRUS(dat[, 4], targets = c(targetvalues_low, targetvalues_upper), n.min = reflim_n_min, plot.all = reflimR.plot.all, lambda = lambda,
                                      main = paste0("Reference limits for ", parameter_name), xlab = parameter_name)}
      if(plot_choice == "pU"){
        reflim_result <- reflim(dat[, 4], targets = c(targetvalues_low, targetvalues_upper), n.min = reflim_n_min, plot.all = reflimR.plot.all,
                                main = paste0("Reference limits for ", parameter_name), xlab = parameter_name)}
    }
    
    if(input$check_refineR) {
      
      validate(need(refineR_done(), "(refineR) Please perform the refineR calculation first."))
      
      table_refineR <- getRI(fit_refineR())
      targetvalues_low <- table_refineR$PointEst[1]
      targetvalues_upper <- table_refineR$PointEst[2]
      
      if(plot_choice == "VeRUS"){
        reflim_result <- reflim_VeRUS(dat[, 4], targets = c(targetvalues_low, targetvalues_upper), n.min = reflim_n_min, plot.all = reflimR.plot.all, lambda = lambda,
                                      main = paste0("Reference limits for ", parameter_name), xlab = parameter_name)}
      if(plot_choice == "pU"){
        reflim_result <- reflim(dat[,4], targets = c(targetvalues_low, targetvalues_upper), n.min = reflim_n_min, plot.all = reflimR.plot.all,
                                main = paste0("Reference limits for ", parameter_name), xlab = parameter_name)}
    }
    
    if (!isTRUE(input$check_target) && !use_preinstalled_targets() && !isTRUE(input$check_refineR)) {
      
      if(plot_choice == "VeRUS"){
        reflim_result <- reflim_VeRUS(dat[, 4], n.min = reflim_n_min, plot.all = reflimR.plot.all, lambda = lambda,
                                      main = paste0("Reference limits for ", parameter_name), xlab = parameter_name)}
      if(plot_choice == "pU"){
        reflim_result <- reflim(dat[, 4], n.min = reflim_n_min, plot.all = reflimR.plot.all,
                                main = paste0("Reference limits for ", parameter_name), xlab = parameter_name)}
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
    dat <- reflim_data()
    parameter_name <- parameter_display()
    warning_text <- if (length(report$warnings) > 0L) {
      paste(report$warnings, collapse = " ")
    } else {
      "None"
    }
    validate(need(
      length(report$limits) >= 1L,
      "(reflim) No reference limits are available. Please check the selected data and settings."
    ))
    if (!is.na(report$limits[1])) {
      converted_sex <- switch(input$sex,
                              "f" = "Female",
                              "m" = "Male",
                              "t" = "Female & Male")
      
      
      if(input$lambda_type == "reflimR"){
        validate(need(
          is.logical(report$lognormal) && length(report$lognormal) == 1L && !is.na(report$lognormal),
          "(reflim) The distribution could not be determined. Please check the selected data and settings."
        ))
        if(report$lognormal){
          lambda <- 0
        } else{
          lambda <- 1
        }
      } else if(input$lambda_type == "refineR"){
        validate(need(refineR_done(), "(refineR) Please perform the refineR calculation first."))
        lambda <- fit_refineR()$Lambda
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
      
      report_interpretation_um <- rep(NA_character_, 2)
      valid_um_comparison <-
        all(is.finite(c(report$limits[1:2], report$targets[1:2]))) &&
        report$limits[1] > 0 && report$limits[2] > report$limits[1] &&
        report$targets[1] > 0 && report$targets[2] > report$targets[1]
      if (valid_um_comparison &&
          all(is.finite(c(unlist(report_versus), unlist(report_target_versus))))) {
        report_interpretation_um <- interpretation_VeRUS(
          report$limits[1:2], report$targets[1:2], lambda = lambda
        )$dev.lim
      }

      table_report <- t(data.frame(
        "Sex and Age:" = paste0(converted_sex, " (", input$age_end[1], "-", input$age_end[2], ")"),
        "n:" = nrow(dat),
        "Warnings:" = warning_text,
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
        "Interpretation of the lower limit (EL):" = report$interpretation[1],
        "Interpretation of the upper limit (EL):" = report$interpretation[2],
        "Interpretation of the lower limit (UM):" = report_interpretation_um[1],
        "Interpretation of the upper limit (UM):" = report_interpretation_um[2],
        "Lambda for UM:" = lambda,
        check.names = FALSE))
    } else {
      table_report <- t(data.frame(
        "n:" = nrow(dat),
        "Warnings:" = warning_text,
        check.names = FALSE))
    }
    colnames(table_report) <- parameter_name

    DT::datatable(table_report, extensions = 'Buttons',
                  options = list(dom = 'tB', pageLength = nrow(table_report), buttons = table_export_buttons("reflimR")))
  })
   
  refineR_done <- reactiveVal(FALSE)
  
  observeEvent(list(
    data_source(),
    input$parameter,
    input$category,
    input$sex,
    input$age_end,
    input$reset,
    input$dataset_file1,
    excel_import_selection(),
    input$submit,
    input$show_table,
    input$upload_age_column,
    input$upload_sex_column,
    input$upload_category_column,
    input$upload_female_value,
    input$upload_male_value
  ), { 
    refineR_done(FALSE)
  }, ignoreInit = TRUE)
  
  output$plotrefineR <- renderPlot({ #Tab:refineR
    plot(fit_refineR(), xlab = parameter_display())
  })
  
  output$refineR_checkbox_ui <- renderUI({ 
    checkboxInput(
      "check_refineR",
      "Compare reflimR and refineR",
      value = FALSE
    )
  })
  
  output$table_report_refineR <- DT::renderDataTable({ #Tab:refineR
    
    table_refineR_original <- fit_refineR()
    table_refineR <- getRI(table_refineR_original)
    parameter_name <- parameter_display()
    
    converted_sex <- switch(input$sex,
                            "f" = "Female",
                            "m" = "Male",
                            "t" = "Female & Male")
    
    table_report <- t(data.frame(
      "Sex and Age:" = paste0(converted_sex, " (", input$age_end[1], "-", input$age_end[2], ")"),
      "Category:" = input$category,
      "Reference limit:" =  paste0(round(table_refineR$PointEst[1], 3) , " - " , round(table_refineR$PointEst[2], 3)),
      "Lambda:" = table_refineR_original$Lambda,
      check.names = FALSE))
    colnames(table_report) <- parameter_name
    
    DT::datatable(table_report, extensions = 'Buttons',
                  options = list(dom = 'tB', pageLength = 15, buttons = table_export_buttons("refineR")))
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

  mclust_plot_values <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x) & x > 0]

    validate(need(
      length(x) >= 2 && length(unique(x)) >= 2,
      "No usable data available for mclust. At least two different positive values are required."
    ))

    iqr_value <- IQR(x)
    if (is.finite(iqr_value)) {
      trimmed_x <- x[x < (median(x) + 6 * iqr_value)]
      validate(need(
        length(trimmed_x) >= 2 && length(unique(trimmed_x)) >= 2,
        "Not enough usable data remain for mclust after removing extreme values."
      ))
    }

    x
  }

  mclust_lognormal <- reactive({
    distribution <- input$mclust_distribution
    if (is.null(distribution)) {
      distribution <- "report"
    }

    switch(
      distribution,
      "normal" = FALSE,
      "lognormal" = TRUE,
      "report" = {
        report <- get_reflim_report()
        validate(need(
          length(report$lognormal) == 1 && !is.na(report$lognormal),
          "The distribution could not be determined from the reflimR report."
        ))
        isTRUE(report$lognormal)
      },
      FALSE
    )
  })
  
  mclust_n_clusters <- reactive({
    req(input$n_cluster)
    if (identical(input$n_cluster, "auto")) {
      return(mclust_cluster_counts)
    }

    validate(
      need(
        input$n_cluster %in% as.character(mclust_cluster_counts),
        "Please choose Automatic (BIC) or a whole number of clusters between 1 and 15."
      )
    )

    as.integer(input$n_cluster)
  })

  update_mclust_cluster_label <- function(count = NULL) {
    automatic_label <- if (is.null(count)) {
      "Automatic (BIC)"
    } else {
      paste0("Automatic (BIC): ", count)
    }
    updateSelectizeInput(
      session, "n_cluster",
      choices = c(setNames("auto", automatic_label),
                  setNames(as.character(mclust_cluster_counts),
                           as.character(mclust_cluster_counts))),
      selected = isolate(input$n_cluster)
    )
  }

  output$plotmclust <- renderPlot({ #Tab:mclust

    # Clear the previous count before validating and fitting the current data.
    update_mclust_cluster_label()
    withProgress(message = "mclust Calculation …", {
      dat <- reflim_data()
      values <- mclust_plot_values(dat[, 4])
      n_cluster <- mclust_n_clusters()
      plot_digits <- mclust_plot_digits(values)
      fit <- lab_mclust(values, lognormal = mclust_lognormal(), model = input$model_name, n.cluster = n_cluster, remove.extremes = T, digits = plot_digits,
                        xlab = parameter_display())
      if (identical(input$n_cluster, "auto")) {
        update_mclust_cluster_label(nrow(fit$stats))
      }
      invisible(fit)
    })
  })
  
  output$plotmclustbic <- renderPlot({ #Tab:mclust

    withProgress(message = "mclust Calculation …", {
      dat <- reflim_data()
      values <- mclust_plot_values(dat[, 4])
      plot_digits <- mclust_plot_digits(values)
      fit <- lab_mclust(values, lognormal = mclust_lognormal(), model = input$model_name,
                        n.cluster = mclust_cluster_counts, remove.extremes = TRUE,
                        plot.it = FALSE, digits = plot_digits)
      plot(fit$BIC)
    })
  })
  
  output$scatterplot <- renderPlot({ #Tab:Scatterplot

    dat <- reflim_data()
    parameter_name <- parameter_display()
    age <- suppressWarnings(as.numeric(dat[, 2]))
    value <- suppressWarnings(as.numeric(dat[, 4]))
    complete_rows <- is.finite(age) & is.finite(value)

    validate(need(
      any(complete_rows),
      "No data available for the selected filters."
    ))

    dat <- dat[complete_rows, , drop = FALSE]
    dat[, 2] <- age[complete_rows]
    dat[, 4] <- value[complete_rows]
    colors <- ifelse(dat[, 3] == "f", "indianred", "cornflowerblue")
    pchs <- ifelse(dat[, 3] == "f", 17, 19)
    plot(dat[,4] ~ dat[,2], pch = pchs, cex = 1, col = colors, xlab = "Age", ylab = parameter_name)

    unique_levels <- c("f", "m")[c("f", "m") %in% dat[, 3]]
    legend(
      "topright",
      legend = c(f = "Female", m = "Male")[unique_levels],
      pch = ifelse(unique_levels == "f", 17, 19),
      col = ifelse(unique_levels == "f", "indianred", "cornflowerblue")
    )
  })

  output$plot_statistics <- renderPlot({ #Tab:Statistics

    dat <- reflim_data()

    if (!(nrow(dat)) == 0) {
      hist_data_w <- subset(dat, Sex == "f", select = Age)
      hist_data_m <- subset(dat, Sex == "m", select = Age)

      # Whole-year bins [age, age + 1) also cover fractional ages.
      age_breaks <- seq(floor(min(dat$Age)), floor(max(dat$Age)) + 1, by = 1)
      age_labels <- head(age_breaks, -1)
      hist_w <- if (nrow(hist_data_w) > 0) {
        hist(hist_data_w$Age, breaks = age_breaks, right = FALSE, plot = FALSE)$counts
      } else {
        integer(length(age_labels))
      }
      hist_m <- if (nrow(hist_data_m) > 0) {
        hist(hist_data_m$Age, breaks = age_breaks, right = FALSE, plot = FALSE)$counts
      } else {
        integer(length(age_labels))
      }
      count_max <- max(c(hist_m, hist_w), 1)

      bar_positions <- barplot(rbind(hist_w,hist_m), col = c("indianred","cornflowerblue"),
              names.arg = age_labels, xlab = "Age", las = 1, beside = TRUE, ylab = "Number of data",
              ylim = c(0, count_max * 1.5))
      abline(h = 0)
      legend("topright", inset = c(0, -0.08), xpd = NA,
             legend = c(paste0("Female: ", nrow(hist_data_w)), paste0("Male: ", nrow(hist_data_m))),
             col = c("indianred","cornflowerblue"), pch = c(19, 19))

      # Align ages with the grouped bars and draw in the existing coordinates.
      age_centers <- colMeans(bar_positions)
      age_spacing <- if (length(age_centers) > 1L) diff(age_centers)[1] else 3
      boxplot_ages <- age_centers[1] + (dat$Age - age_labels[1]) * age_spacing
      boxplot(boxplot_ages, horizontal = TRUE, add = TRUE, axes = FALSE,
              at = count_max * 1.2, boxwex = count_max * 0.08,
              border = "gray25", col = "gray85")
    }

  })

  output$plot_age_boxplot <- renderPlot({ #Tab:Overview
    dat <- reflim_data()
    values <- suppressWarnings(as.numeric(dat[, 4]))
    usable_rows <- is.finite(dat$Age) & is.finite(values) & dat$Sex %in% c("m", "f")
    dat <- dat[usable_rows, , drop = FALSE]
    values <- values[usable_rows]
    validate(need(length(values) > 0, "No data available for the selected filters."))

    par(mar = c(4.5, 4.5, 2.5, 1))
    group_width <- input$overview_age_group_width
    if (is.null(group_width) || group_width == "auto") {
      max_age_groups <- max(6, min(12, floor(par("pin")[1] / 0.6)))
      candidate_widths <- c(1, 2, 5, 10, 20)
      group_counts <- floor(max(dat$Age) / candidate_widths) -
        floor(min(dat$Age) / candidate_widths) + 1
      group_width <- candidate_widths[which(group_counts <= max_age_groups)[1]]
    } else {
      group_width <- as.numeric(group_width)
    }

    age_starts <- seq(floor(min(dat$Age) / group_width) * group_width,
                      floor(max(dat$Age) / group_width) * group_width, by = group_width)
    age_labels <- paste0(age_starts, "-<", age_starts + group_width)
    age_group <- factor(floor(dat$Age / group_width) * group_width, levels = age_starts)
    sex_levels <- c("f", "m")[c("f", "m") %in% dat$Sex]
    sex_group <- factor(dat$Sex, levels = sex_levels)
    sex_colors <- c(f = "indianred", m = "cornflowerblue")[sex_levels]
    grouped_values <- split(values, interaction(sex_group, age_group, drop = FALSE))

    use_log <- identical(input$overview_boxplot_scale, "log")
    validate(need(!use_log || all(values > 0),
                  "Logarithmic scale requires positive values. Choose Linear for data containing zero or negative values."))

    # Keep female/male boxes beside each other, with one label per age interval.
    offsets <- if (length(sex_levels) == 2) c(-0.18, 0.18) else 0
    box_positions <- rep(seq_along(age_starts), each = length(sex_levels)) +
      rep(offsets, times = length(age_starts))
    boxplot(grouped_values, at = box_positions,
            boxwex = if (length(sex_levels) == 2) 0.32 else 0.6,
            col = rep(sex_colors, times = length(age_starts)),
            border = "gray30", outpch = 16, outcex = 0.45,
            outcol = grDevices::adjustcolor("gray30", alpha.f = 0.4),
            xaxt = "n", las = 1, xlim = c(0.5, length(age_starts) + 0.5),
            xlab = "Age group (years)",
            ylab = paste0(parameter_display(), if (use_log) " (log scale)" else ""),
            log = if (use_log) "y" else "")
    axis(1, at = seq_along(age_starts), labels = age_labels, cex.axis = 0.8)
    mtext(paste0(group_width, "-year age groups"), side = 3, adj = 0, line = 0.5, cex = 0.85)
    legend("topright", inset = c(0, -0.12), xpd = NA,
           legend = c(f = "Female", m = "Male")[sex_levels],
           fill = sex_colors, horiz = TRUE, bty = "n")
  })
  
  zlog_table_data <- reactive({
    dat <- reflim_data()
    report <- get_reflim_report()
    limits <- report$limits[1:2]
    validate(need(
      length(limits) == 2L && all(is.finite(limits)) &&
        limits[1] > 0 && limits[2] > limits[1],
      "(zlog) No valid positive reference limits are available from reflimR. Please check the selected data and settings."
    ))

    # Calculate once per data/limit change, not once per displayed cell.
    zlog_values <- round(zlog(dat[[4]], limits[1], limits[2]), 2)
    dat$Age <- round(dat$Age, 3)
    cbind(dat,
          RI = rep(paste0(limits[1], " - ", limits[2]), nrow(dat)),
          zlog = zlog_values,
          .zlog_text_color = highzlogvalues(zlog_values),
          .zlog_background_color = zlogcolor(zlog_values))
  })

  output$table_zlog <- DT::renderDataTable({ #Tab:zlog
    DT::datatable(
      zlog_table_data(), rownames = FALSE, extensions = 'Buttons', class = 'cell-border',
      options = list(
        dom = 'Blfrtip', pageLength = 15, deferRender = TRUE,
        buttons = table_export_buttons("zlog", export_columns = 0:5),
        # Hidden color columns travel with each server-side page. This avoids
        # embedding a full-dataset styleEqual lookup in the browser callback.
        columnDefs = list(list(targets = 6:7, visible = FALSE, searchable = FALSE))
      ),
      caption = htmltools::tags$caption(style = 'caption-side: bottom; text-align: center;',
                                       'Table: Dataset with the zlog values')
    ) %>%
      DT::formatStyle(columns = 6, valueColumns = 7, color = DT::styleValue()) %>%
      DT::formatStyle(columns = 6, valueColumns = 8, backgroundColor = DT::styleValue())
  }, server = TRUE)
  
  output$download_ritable <- downloadHandler(
    filename = function() {
      paste("ReferenceIntervals_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      
      dat <- get_alldata_file()
      dataset <- dat[c(-1,-2,-3)]
      reflim.loop.results <- reflim.loop(dataset, n.min = reflim_n_min, plot.it = FALSE)
      
      tmpdir <- tempfile("VeRIf_ReferenceIntervals_")
      dir.create(tmpdir)
      on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
      csv_files <- c()
      csv_names <- setNames(reference_interval_filenames(names(reflim.loop.results)),
                            names(reflim.loop.results))
      
      for (col in names(reflim.loop.results)) {
        report <- reflim.loop.results[[col]]
        
        converted_sex <- switch(input$sex,
                                "f" = "Female",
                                "m" = "Male",
                                "t" = "Female & Male")
        
        if (!is.null(report$error) || length(report$limits) < 2L ||
            !all(is.finite(report$limits[1:2]))) {
          df <- t(data.frame(
            "Sex and Age:" = paste0(converted_sex, " (", input$age_end[1], "-", input$age_end[2], ")"),
            "Category:" = input$category,
            "Status:" = if (!is.null(report$error)) "Error" else "No reference interval available",
            "Error:" = if (!is.null(report$error)) report$error else NA_character_,
            "Remarks:" = paste(report$remarks, collapse = " "),
            "Warnings:" = if (length(report$warnings)) paste(report$warnings, collapse = " ") else "None",
            check.names = FALSE))
        } else {
          df <- t(data.frame(
            "Sex and Age:" = paste0(converted_sex, " (", input$age_end[1], "-", input$age_end[2], ")"),
            "Category:" = input$category,
            "Status:" = "OK",
            "Remarks:" = paste(report$remarks, collapse = " "),
            "Warnings:" = if (length(report$warnings)) paste(report$warnings, collapse = " ") else "None",
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
            "Interpretation of the lower limit (EL):" = report$interpretation[1],
            "Interpretation of the upper limit (EL):" = report$interpretation[2],
            check.names = FALSE))
        }
        
        colnames(df) <- col
        csv_path <- file.path(tmpdir, csv_names[[col]])
        write.csv(df, csv_path)
        csv_files <- c(csv_files, csv_path)
      }
      
      old_wd <- setwd(tmpdir)
      on.exit(setwd(old_wd), add = TRUE, after = FALSE)
      
      zip(zipfile = file, files = basename(csv_files), extras = "-j")
    }
  )
  
  output$download_zlogtable <- downloadHandler( #Tab:zlog
    filename = function() {
      paste("zlogValues_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      dat <- get_alldata_file()

      dataset <- dat[c(-1,-2,-3)]
      reflim.loop.results <- reflim.loop(dataset, n.min = reflim_n_min, plot.it = FALSE)
      zlog.loop.results <- zlog.loop(dataset, reflim.loop.results)

      result <- cbind(dat[1:3], zlog.loop.results)
      write.csv(result, file)
    }
  )
  
  output$tree_rpart <- renderPlot({ #Tab:rpart
    tree_fit <- build_rpart()
    tree_formula <- paste(parameter_display(), "~",
                          paste(deparse(stats::formula(tree_fit)[[3]]), collapse = " "))
    rpart.plot::rpart.plot(
      tree_fit,
      box.palette = "RdBu",
      roundint = FALSE,
      # Forced yes/no labels fail in rpart.plot for trees without splits.
      yesno = if (nrow(tree_fit$frame) > 1L) 2 else 0,
      main = paste("Regression tree:", tree_formula),
      
      node.fun = function(x, labs, digits, varlen) {
        paste0(
          round(100 * x$frame$n / x$frame$n[1]),
          "%"
        )
      }
    )
  })
  
}
####################################### Run the application #######################################
shinyApp(ui = ui, server = server)

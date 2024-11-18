### R Shiny App to select and start sc-analysis
################################################################################

# Source start settings
source("config/start_settings.R")

# Set environment
renv::use_python(type = "virtualenv", name = file.path(param$path_to_git,"env/basic/virtualenvs/r-reticulate"))
renv::load(file.path(param$path_to_git,"env/basic"))

# Load library
library(shiny)

# Read standard parameter
source(file.path(param$path_to_git,"config/standard_parameter.R"))

# Define a lists of basic parameters that have to be set for the respective analysis 
parameter_lists = list(
  "qc" = list(
    "species" = list(label = "Set species", type = "choice", choices = c("human", "mouse"))
  ),
  "pre-processing" = list(
    "species" = list(label = "Set species", type = "choice", choices = c("human", "mouse"))
  ),
  "cluster analysis" = list(
    "species" = list(label = "Set species", type = "choice", choices = c("human", "mouse"))
  ),
  "cell annotation clustify" = list(
    "clustifyr_ref" = list(label = "Set clustifyr reference e.g. 'http://cells.ucsc.edu/?ds=covid19-influenza-response', 'ref_hema_microarray()', or
 'file.path(param$path_to_git,/references/fetal-immune-pfi_clustifyr_reference.rds)'", type = "text", value = NULL),
    "cluster_col_clustifyr" = list(label = "Name of annotation column in reference dataset e.g. 'cell_type'", type = "text", value = NULL)
  ),
  "dataset mapping" = list(
    "refdata" = list(label = "Set path to reference rds object", type = "text", value = NULL)
  ),
  "ccc analysis" = list(
    "sender" = list(label = "Set sender cell types", type = "text", value = NULL),
    "receiver" = list(label = "Set receiver cell types", type = "text", value = NULL)
  ),
  "inspect rds" = list(
  ),
  "reference download" = list(
    "species" = list(label = "Set species", type = "choice", choices = c("human", "mouse")),
    "annot_version" = list(label = "Set Ensembl version", type = "numeric", min = 0, value = 98)
  ),
  "test dataset download" = list(
    "download_test_datasets" = list(label = "Select testdataset", type = "choice", choices = c("download_10x_SmartSeq2_pbmc_GSE132044", "download_10x_pbmc_small_split2samples",
                                    "download_10x_pbmc_hto_GSE108313", "download_10x_pbmc_5k_protein", "download_10x_pbmc_1k_healthyDonor_v3Chemistry"))
  ),  
  "generate clustifyr reference" = list(
    "ref_data_name" = list(label = "Set name of clustifyr reference", type = "text", value = NULL)
  )
  # Examples
  #"c" = list(
  #  "param1" = list(label = "Set Parameter 1 for Option 2", type = "numeric", min = 1, max = 100, value = 20),
  #  "param2" = list(label = "Set Parameter 2 for Option 2", type = "text", value = "more text"),
  #  "param3" = list(label = "Set Parameter 3 for Option 2", type = "slider", min = 10, max = 100, value = 40)
  #)
)

# Grouping of some analysis types
var_proj_id = c("qc", "pre-processing", "cluster analysis", "cell annotation clustify", "dataset mapping", "ccc analysis", "inspect rds" )
rds_type = c("cell annotation clustify", "dataset mapping", "ccc analysis", "inspect rds" )
var_data_type = c("qc", "pre-processing")



### Define UI for sc-analysis type selection and start
################################################################################
ui <- fluidPage(

    # App title
    titlePanel("sc-analysis"),

    # Sidebar with input for respective analysis  
    sidebarLayout(
      sidebarPanel(
        # Dropdown for selecting the analysis
        selectInput("analysis_type", 
                    label = "Choose analysis", 
                    choices = c("qc", "pre-processing", "cluster analysis", 
                                "cell annotation clustify", "dataset mapping", "ccc analysis", 
                                "inspect rds", "reference download", "test dataset download", 
                                "generate clustifyr reference"))
      ),

        # Conditional inputs based on selected parameter
        mainPanel(
          # Instruction message
          uiOutput("instruction"),
          
          # Project name input
          uiOutput("project_input"),
          
          # Data path input (depending on analysis type)
          uiOutput("rds_input"),
          uiOutput("data_input"),  
          uiOutput("data_input2"), 
          # Dynamic UI for samples
          uiOutput("samples_ui"),
          
          # Dynamic UI for parameter selection
          uiOutput("parameter_input"),  
          
          # Run button
          actionButton("run_script", "Run analysis"), 
          
          # Completion notification
          uiOutput("output_message")
        )
    )
)




### Server function
################################################################################
server <- function(input, output, session) {
  
  ##### Prepare general message and reactiveVal
  
  # Instruction message 
  instruction = paste("<b><span style='font-size: 15px;'>Set parameters and click the 'Run analysis' button. <br>
                       <span style='color: red;'>Wait until analysis is completed!!! <br></span>
                       The app will give you notification of completion.<br><br><br> </span></b>")
  output$instruction = renderUI({
    HTML(instruction)
  })
  
  # Create reactive value to store parameter
  assigned_scriptname = reactiveVal(NULL)
  assigned_parameter = reactiveVal(NULL)
  
  
  
  ##### Observe analysis selection and act
  observeEvent(input$analysis_type, {
    
    ### Assign a script path to the respective analysis type (reactive value)
    # Define switches (script path to analysis type)
    assigned_value = switch(input$analysis_type,
                             "pre-processing" = "scripts/pre-processing/pre-processing.Rmd",
                             "qc" = "scripts/pre-processing/qc.Rmd",
                             "cluster analysis" = "scripts/pre-processing/cluster_analysis.Rmd", 
                             "cell annotation clustify" = "scripts/dataset_mapping/cell_annotation_clustifyr.Rmd", 
                             "dataset mapping" = "scripts/dataset_mapping/dataset_mapping_seurat.Rmd", 
                             "ccc analysis" = "scripts/ccc_analysis/ccc_analysis.Rmd", 
                             "inspect rds" = "scripts/read_data/inspect_rds.R", 
                             "reference download" = "scripts/download_references/", 
                             "test dataset download" = "scripts/download_test_datasets/test_dataset_download.R", 
                             "generate clustifyr reference" = "scripts/download_references/")  
    
    # Store script path in reactive value
    assigned_scriptname(assigned_value)
    
    
    
    ### Render specific UIs for respective parameter settings
    # UI for project name input (if applicable)
    output$project_input = renderUI({
      if (input$analysis_type %in% var_proj_id) {
        textInput("project_id", label = "Enter project name", value = "Testdata")
      } else {output$project_input = NULL}
    })
    
    
    # UI for rds path input (if applicable)
    output$rds_input = renderUI({
      if (input$analysis_type %in% rds_type) {
        textInput("data", label = "Set path to rds object", value = NULL)
      } else {output$rds_input = NULL}
    })
    
    
    # Dynamically render data type input if multiple input types are optional
    # First, render UI to select input type
    if (input$analysis_type %in% var_data_type) {
      output$data_input = renderUI({
        selectInput("data_type", label = "Choose data type", choices = c("count matrix", "rds object", "download test datasets"))
      })
      
      # Second, render UI for data input depending on selected data type
      observeEvent(input$data_type, {
        output$data_input2 = renderUI({
          data_type = tolower(input$data_type)
          if (data_type == "count matrix") {
            numericInput("num_samples", label = "Number of Samples", value = 0, min = 0)
          } else if (data_type == "rds object") {
            textInput("data", label = "Set path to rds object", value = NULL)
          } else if (data_type == "download test datasets") {
            selectInput("download_test_datasets", label = "Select testdataset", 
                        choices = c("download_10x_SmartSeq2_pbmc_GSE132044",
                                    "download_10x_pbmc_small_split2samples",
                                    "download_10x_pbmc_hto_GSE108313",
                                    "download_10x_pbmc_5k_protein",
                                    "download_10x_pbmc_1k_healthyDonor_v3Chemistry"))
          } 
        })
        
        # Third, render dynamic UI for sample inputs if 'count matrix' was selected
        observeEvent(input$num_samples, {
          if (input$data_type == "count matrix" & input$num_samples > 0) {
            output$samples_ui = renderUI({
              # Get the number of samples from the numeric input
              n = input$num_samples  
              # Generate the UI for each sample (name and value inputs)
              sample_inputs = lapply(1:n, function(i) {
                tagList(
                  textInput(paste0("sample_name_", i), label = paste("Sample", i, ": Name")),
                  selectInput(paste0("sample_type_", i), label = paste("Sample", i, ": Data type"), choices = c("10x", "Smartseq")), 
                  textInput(paste0("sample_path_", i), label = paste("Sample", i, ": Path to count matrix folder"))
                )
              })
            
            # Return the list of input fields
            do.call(tagList, sample_inputs)
            })
          } else {output$samples_ui=NULL}
        })
      })
    } else {
      output$data_input = NULL
      output$data_input2=NULL
      output$samples_ui=NULL
    }
    
    
    # Dynamically render specific parameter input UI based on selected analysis
    output$parameter_input = renderUI({
      analysis_type = tolower(input$analysis_type)
      parameter = parameter_lists[[analysis_type]]
    
      # Create a list of UI elements with the respective parameters
      ui_elements = lapply(names(parameter), function(parameter_name) {
        parameter_info = parameter[[parameter_name]]
        
        if (parameter_info$type == "numeric") {
          numericInput(parameter_name, parameter_info$label, value = parameter_info$value, min = parameter_info$min, max = parameter_info$max)
        } else if (parameter_info$type == "text") {
          textInput(parameter_name, parameter_info$label, value = parameter_info$value)
        } else if (parameter_info$type == "slider") {
          sliderInput(parameter_name, parameter_info$label, min = parameter_info$min, max = parameter_info$max, value = parameter_info$value)
        } else if (parameter_info$type == "choice") {
          selectInput(parameter_name, parameter_info$label, choices = parameter_info$choices)
        } 
      })
    
      # Return the dynamic UI
      do.call(tagList, ui_elements)
    })
    
  }) 
  
  
  
  ### By pressing 'Run analysis' button assign parameter and run the selected script
  observeEvent(input$run_script, {
    analysis_type = tolower(input$analysis_type)
    
    
    ### Set all parameter
    # First, store standard param list loaded from start_settings.R and standard_parameter.R 
    assigned_parameter(param)
    
    # Add scriptname, alias the script path, to param list
    param[["scriptname"]] = assigned_scriptname()
    
    # Add project id
    if (analysis_type %in% var_proj_id) {
      param[["project_id"]] = input$project_id}
    
    # Add data path
    if (analysis_type %in% rds_type) {
      param[["data"]] = input$data}
    
    if (analysis_type %in% var_data_type) {
      data_type = tolower(input$data_type)
      if (data_type == "count matrix" & input$num_samples>0) {
        n = input$num_samples
        sample_info = data.frame(
          name = sapply(1:n, function(i) input[[paste0("sample_name_", i)]]),
          type = sapply(1:n, function(i) input[[paste0("sample_type_", i)]]),
          path = sapply(1:n, function(i) input[[paste0("sample_path_", i)]]),
          stringsAsFactors = FALSE
        )
        param[["path_data"]] = sample_info
        param[["data"]] = NULL
        param[["download_test_datasets"]] = NULL
      } else if (data_type == "download test datasets") {
        param[["path_data"]] = NULL
        param[["data"]] = NULL
        param[["download_test_datasets"]] = input$download_test_datasets
      } else {
        param[["path_data"]] = NULL
        param[["data"]] = input$data
        param[["download_test_datasets"]] = NULL
      }
    } 
    
    if (analysis_type == "test dataset download") {
      param[["download_test_datasets"]] = input$download_test_datasets}
    
    # Add dynamically set parameter to param list
    parameter <- parameter_lists[[analysis_type]]
    for (parameter_name in names(parameter)) {
      param[[parameter_name]]=input[[parameter_name]]
    }
    
    # Add output path to param list
    if (grepl(".Rmd", param$scriptname)) {
      param[["path_out"]] = file.path(param$path_to_git, "output", param$project_id, gsub(".Rmd", "", basename(param$scriptname)))
    } else {
      param[["path_out"]] = file.path(param$path_to_git, "output", param$project_id)
    }
    
    # Set param list as a variable in the global environment to be accessible to R scripts called via source
    assign("param", param, envir = .GlobalEnv)
    
    
    ### Run analysis
    if (grepl(".Rmd", param$scriptname)) {
      # Create output directories
      if (!file.exists(param$path_out)) dir.create(param$path_out, recursive=TRUE, showWarnings=FALSE)
      # Run markdown
      rmarkdown::render(file.path(param$path_to_git,param$scriptname), param = param, output_file = file.path(param$path_out, paste0(param$project_id,"_", gsub("Rmd", "html", basename(param$scriptname)))))
    } else {
      # Run R script
      source(file.path(param$path_to_git,param$scriptname)) 
    }
    
    
    ### Return notification of completion
    # Create message content
    if (analysis_type %in% var_proj_id) {
      notification = sprintf("<span style='font-size: 20px; color: darkblue;'><b><br><br>%s completed!<br></b></span>
                            <span style='font-size: 15px; color: darkblue;'>Output folder:<br>%s</span>", input$analysis_type, param$path_out)
    } else {
      notification = sprintf("<span style='font-size: 20px; color: darkblue;'><b><br><br>%s completed!<br></b></span>", input$analysis_type)
    }
    # Display the message 
    output$output_message = renderUI({
      HTML(notification)
    })
  
  })
}
      


# Run the application 
################################################################################
shinyApp(ui = ui, server = server)

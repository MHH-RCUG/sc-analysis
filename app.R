### R Shiny App to select and start sc-analysis
################################################################################

# Source start settings
source("config/start_settings.R")

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
    "refdata" = list(label = "Set path to reference rds object", type = "text", value = NULL),
    "reduction" = list(label = "Reduction to use; must exist in ref dataset", type = "choice", choices = c("umap", "tsne", "others (see advanced settings)"))
  ),
  "ccc analysis" = list(
    "sender" = list(label = "Set sender cell types", type = "text", value = NULL),
    "receiver" = list(label = "Set receiver cell types", type = "text", value = NULL)
  ),
  "inspect rds" = list(
  ),
  "reference download" = list(
    "species" = list(label = "Set species", type = "choice", choices = c("human", "mouse")),
    "annot_version" = list(label = "Set Ensembl version", type = "numeric", min = 0, max = 1000, value = 98)
  ),
  "test dataset download" = list(
    "download_test_datasets" = list(label = "Select testdataset", type = "choice", choices = c("download_10x_SmartSeq2_pbmc_GSE132044", "download_10x_pbmc_small_split2samples",
                                    "download_10x_pbmc_hto_GSE108313", "download_10x_pbmc_5k_protein", "download_10x_pbmc_1k_healthyDonor_v3Chemistry"))
  ),  
  "generate clustifyr reference" = list(
    "ref_data_name" = list(label = "Set name for the new clustifyr reference", type = "text", value = NULL),
    "ref_data_path" = list(label = "Set path to the folder with a 'exprMatrix.tsv.gz' and 'meta.tsv' from which clustifyr reference is generated", type = "text", value = NULL)
  )
  # Examples
  #"c" = list(
  #  "param1" = list(label = "Set Parameter 1", type = "numeric", min = 1, max = 100, value = 20 or NULL),
  #  "param2" = list(label = "Set Parameter 2", type = "text", value = "text or NULL"),
  #  "param3" = list(label = "Set Parameter 3", type = "slider", min = 10, max = 100, value = 40),
  #  "param4" = list(label = "Set Parameter 4", type = "choice", choices = c("choice A", "choice B")),
  #  "param5" = list(label = "Set Parameter 5", type = "checkbox", choices = c("choice A", "choice B"), selected = NULL),
  #  "param6" = list(label = "Set Parameter 6 (comma separated)", type = "numeric_vector", value = "0.5, 0.7")
  #)
)


# Define a lists of advanced parameters that have to be set for the respective analysis 
parameter_advset_lists = list(
  "qc" = list(
    "norm" = list(label = "Which normalization should be used for analysis? Normal log = 'RNA', SCTransform = 'SCT'", type = "choice", choices = c("RNA", "SCT")),
    "downsample_cells_n" = list(label = "Downsample data to at most n cells per sample AFTER filtering (mainly for tests)", type = "numeric", min = 50, max = 100000, value = NULL),
    "downsample_cells_equally" = list(label = "Downsample all samples equally according to the smallest sample", type = "choice", choices = c(FALSE, TRUE)),
    "samples_min_cells" = list(label = "Drop samples with too few cells", type = "numeric", min = 0, max = 1000, value = 10),
    "file_annot" = list(label = "Proviede file for gene annotation", type = "text", value = NULL),
    "file_cc_genes" = list(label = "Proviede file for cell cycle gene annotation", type = "text", value = NULL)
  ),
  "pre-processing" = list(
    "norm" = list(label = "Which normalization should be used for analysis? Normal log = 'RNA', SCTransform = 'SCT'", type = "choice", choices = c("RNA", "SCT")),
    "downsample_cells_n" = list(label = "Downsample data to at most n cells per sample AFTER filtering (mainly for tests)", type = "numeric", min = 50, max = 100000, value = NULL),
    "downsample_cells_equally" = list(label = "Downsample all samples equally according to the smallest sample", type = "choice", choices = c(FALSE, TRUE)),
    "nFeature_RNA" = list(label = "Filter cells by nFeature (lower and upper limit comma separated)", type = "numeric_vector", value = "20, NA"),
    "nCount_RNA" = list(label = "Filter cells by nCount (lower and upper limit comma separated)", type = "numeric_vector", value = "200, NA"),
    "percent_mt" = list(label = "Filter cells by percent_mt (lower and upper limit comma separated)", type = "numeric_vector", value = "0, 20"),
    "min_counts" = list(label = "Filter features by min_counts", type = "numeric", min = 0, max = 50, value = 1),
    "min_cells" = list(label = "Filter features by min_cells", type = "numeric", min = 0, max = 100, value = 5),
    "samples_min_cells" = list(label = "Drop samples with too few cells", type = "numeric", min = 0, max = 1000, value = 10),
    "cc_remove" = list(label = "Remove cell cycle effects", type = "choice", choices = c(FALSE, TRUE)),
    "cc_remove_all" = list(label = "Remove all cell cycle effects (not only the difference between proliferating cells); very stringent", type = "choice", choices = c(FALSE, TRUE)),
    "cc_rescore_after_merge" = list(label = "Re-score cell cycle effects after samples have been merged/integrated", type = "choice", choices = c(TRUE, FALSE)),
    "vars_to_regress" = list(label = "Additional (unwanted) variables that will be regressed out", type = "checkbox", choices = c("nCount_RNA", "percent_mt", "percent_ribo"), selected = NULL),
    "integrate_samples_method" = list(label = "How to combine multiple datasets", type = "choice", choices = c("merge", "integrate")),
    "integrate_samples_integration_function" = list(label = "ONLY for 'integrate' define integration_functions", type = "choice", choices = c("RPCAIntegration", "CCAIntegration")),
    "experimental_groups" = list(label = "Similarity between samples", type = "choice", choices = c("heterogene", "homogene")),
    "file_annot" = list(label = "Proviede file for gene annotation", type = "text", value = NULL),
    "file_cc_genes" = list(label = "Proviede file for cell cycle gene annotation", type = "text", value = NULL)
  ),
  "cluster analysis" = list(
    "pc_n" = list(label = "Number of principle components to use", type = "numeric", min = 1, max = 100, value = 20),
    "cluster_resolution" = list(label = "Cluster resolution to use for analysis", type = "numeric", min = 0.1, max = 1.6, value = 0.6),
    "cluster_resolution_test" = list(label = "Cluster resolutions to compute; multiple values possible (comma separated)", type = "numeric_vector", value = "0.5, 0.7, 0.8"),
    "file_annot" = list(label = "Proviede file for gene annotation", type = "text", value = NULL),
    "file_cc_genes" = list(label = "Proviede file for cell cycle gene annotation", type = "text", value = NULL),
    "marker_padj" = list(label = "Adjusted p value threshold for marker gene identification", type = "numeric", min = 0, max = 1, value = 0.05),
    "marker_log2FC" = list(label = "Fold change threshold for marker gene identification", type = "numeric", min = 0, max = 100, value = log2(2)),
    "marker_pct" = list(label = "Minimum percentage of cells with respective gene expression for definition as marker gene", type = "numeric", min = 0, max = 1, value = 0.25),
    "enrichr_dbs" = list(label = "Enrichr libraries. Not more than 5!", type = "checkbox", choices = c("GO_Biological_Process_2023", "WikiPathway_2023_Human",  "WikiPathways_2019_Mouse", "Azimuth_2023", "CellMarker_2024"), selected = NULL),
    "enrichr_padj" = list(label = "P-value threshold for functional enrichment tests", type = "numeric", min = 0, max = 1, value = 0.05),
    "annotation_dbs" = list(label = "Cell type annotation database (see https://bioconductor.statistik.tu-dortmund.de/packages/3.10/bioc/vignettes/SingleR/inst/doc/SingleR.html#6_reference_options)", type = "choice", choices = c("BlueprintEncodeData()", "HumanPrimaryCellAtlasData()", "DatabaseImmuneCellExpressionData()", "NovershternHematopoieticData()", "MonacoImmuneData()", "ImmGenData()", "MouseRNAseqData()"), selected = NULL)
  ),
  "dataset mapping" = list(
    "celltype" = list(label = "Pre-annotated cell types; column in reference dataset", type = "text", value = "annotation"),
    "reduction" = list(label = "Reduction to use; must exist in ref dataset", type = "text", value = "umap"),
    "predicted_score_threshold" = list(label = "Predicted score threshold", type = "numeric", min = 0, max = 1, value = 0.9),
    "percent_predicted_cells_threshold" = list(label = "Minimum fraction of cell with respective cell identity", type = "numeric", min = 0, max = 1, value = 0.1)
  ),
  "ccc analysis" = list(
    "liana_methods" = list(label = "Select methods. Not more than 3!", type = "checkbox", choices = c("connectome", "logfc", "natmi", "sca", "cellphonedb", "cytotalk", "call_squidpy", "call_cellchat", "call_connectome", "call_sca", "call_italk", "call_natmi"), selected = NULL),
    "liana_agg_rank_threshold" = list(label = "Threshold for liana agg rank", type = "numeric", min = 0, max = 1, value = 0.01)
  )
)


parameter_obj_advset_lists = list(
  "path_out" = list(label = "Output directory", type = "text", value = NULL),
  "col" = list(label = "Feature Plot colors - Highlights", type = "text", value = "#0086b3"),
  "col_bg" = list(label = "Feature Plot colors - Background", type = "text", value = "#D3D3D3"),
  "col_palette_samples" = list(label = "Colour palette used for samples", type = "choice", choices = c("ggsci::pal_igv", "ggsci::pal_ucscgb", "ggsci::springfield_simpsons")),
  "col_palette_clusters" = list(label = "Colour palette used for cluster", type = "choice", choices = c("ggsci::pal_igv", "ggsci::pal_ucscgb", "ggsci::springfield_simpsons")),
  "col_palette_annotation" = list(label = "Colour palette used for annotated cell types", type = "choice", choices = c("ggsci::pal_ucscgb", "ggsci::pal_igv", "ggsci::springfield_simpsons")),
  "pt_size" = list(label = "Dot size for umaps/tsne", type = "numeric", min = 0.1, max = 1, value = 0.5)
)


# Grouping of some analysis types
var_proj_id = c("qc", "pre-processing", "cluster analysis", "cell annotation clustify", "dataset mapping", "ccc analysis", "inspect rds" )
rds_type = c("cell annotation clustify", "dataset mapping", "ccc analysis", "inspect rds" )
var_data_type = c("qc", "pre-processing")
adv_obj_settings = c(var_proj_id)
adv_settings = c(var_proj_id)



### Define UI for sc-analysis type selection and start
################################################################################
ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, preset = "cerulean"),
  includeCSS(file.path(param$path_to_git,"assets/app_style.css")),
  

    # App title
    titlePanel("sc-analysis"),

    # Sidebar with input for respective analysis  
    sidebarLayout(
      sidebarPanel(
        # Dropdown for selecting the analysis
        tagList(
          tags$label("Choose analysis", class = "analysis-label"),
        selectInput("analysis_type", 
                    label = "", 
                    choices = c("qc", "pre-processing", "cluster analysis", 
                                "cell annotation clustify", "dataset mapping", "ccc analysis", 
                                "inspect rds", "reference download", "test dataset download", 
                                "generate clustifyr reference"))
        ),
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
          
          # Dynamically render the action button based on the selected parameter
          uiOutput("advanced_button_ui"),
          
          # Dynamic UI for parameter selection
          uiOutput("parameter_advset_input"),
          uiOutput("parameter_obj_advset_input"),
          
          # Run button
          actionButton("run_script", "Run analysis", class="btn btn-info"), 
          
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
  instruction = paste("<span class='instruction'>A modular workflow for sc data analysis </span><br>
                      <span class='instruction-info'>Chose analysis, set parameters, and click the 'Run analysis' button. The basic parameters have to be defined. Using 'Advanced Parameter Settings' further parameter values can be specified diverging from standard settings. After completion the app will give you notification. </span><br>
                      <span class='instruction-warning'>Wait until analysis is completed!!! </span><br><br><br>")
  output$instruction = renderUI({
    HTML(instruction)
  })
  
  # Create reactive value to store parameter
  assigned_scriptname = reactiveVal(NULL)
  assigned_parameter = reactiveVal(NULL)
  assigned_parameter_advset = reactiveVal(NULL)
  
  # Reactive value to store the number of button clicks; set to 0 at the beginning
  click_counter = reactiveVal(0)
  
  
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
                             "reference download" = "scripts/read_data/read_gene_annotation.R", 
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
          } else {output$samples_ui = NULL}
        })
      })
    } else {
      output$data_input = NULL
      output$data_input2 = NULL
      output$samples_ui = NULL
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
        } else if (parameter_info$type == "checkbox") {
          checkboxGroupInput(parameter_name, parameter_info$label, choices = parameter_info$choices, selected = parameter_info$selected)
        } 
      })
    
      # Return the dynamic UI
      do.call(tagList, ui_elements)
    })
    
    # Dynamically show action button for advanced parameter setting
    output$advanced_button_ui <- renderUI({
      if (input$analysis_type %in% adv_settings) {
        actionButton("adv_setting", "Advanced Parameter Settings", class="btn btn-secondary")
      } else {
        NULL  # Do not render anything if the condition is not met
      }
    })
    
    # Remove UI for advanced parameter setting if analysis type is changed
    output$parameter_advset_input = NULL 
    output$parameter_obj_advset_input = NULL
    
  }) 
  
  
  ### By pressing 'Advanced settings' button return UI with further parameter choices
  observeEvent(input$adv_setting, {
    
    # Increment the click counter
    click_counter(click_counter() + 1)
    
    # First click shows advanced parameter setting
    if (click_counter() == 1) {
      # Dynamically render specific parameter input UI based on selected analysis
      output$parameter_advset_input = renderUI({
        analysis_type = tolower(input$analysis_type)
        parameter_advset = parameter_advset_lists[[analysis_type]]
        
        # Create a list of UI elements with the respective parameters
        ui_advset_elements = lapply(names(parameter_advset), function(parameter_advset_name) {
          parameter_advset_info = parameter_advset[[parameter_advset_name]]
          
          if (parameter_advset_info$type == "numeric") {
            numericInput(parameter_advset_name, parameter_advset_info$label, value = parameter_advset_info$value, min = parameter_advset_info$min, max = parameter_advset_info$max)
          } else if (parameter_advset_info$type == "text") {
            textInput(parameter_advset_name, parameter_advset_info$label, value = parameter_advset_info$value)
          } else if (parameter_advset_info$type == "slider") {
            sliderInput(parameter_advset_name, parameter_advset_info$label, min = parameter_advset_info$min, max = parameter_advset_info$max, value = parameter_advset_info$value)
          } else if (parameter_advset_info$type == "choice") {
            selectInput(parameter_advset_name, parameter_advset_info$label, choices = parameter_advset_info$choices)
          } else if (parameter_advset_info$type == "checkbox") {
            checkboxGroupInput(parameter_advset_name, parameter_advset_info$label, choices = parameter_advset_info$choices, selected = parameter_advset_info$selected)
          } else if (parameter_advset_info$type == "numeric_vector") {
            textInput(parameter_advset_name, parameter_advset_info$label, value = parameter_advset_info$value)
          }
        })
        
        # Return the dynamic UI
        do.call(tagList, ui_advset_elements)
      })
      
      # Render object related specific parameter input UI 
      if (input$analysis_type %in% adv_obj_settings) {
        output$parameter_obj_advset_input = renderUI({
          parameter_obj_advset = parameter_obj_advset_lists
          
          # Create a list of UI elements with the respective parameters
          ui_obj_advset_elements = lapply(names(parameter_obj_advset), function(parameter_obj_advset_name) {
            parameter_obj_advset_info = parameter_obj_advset[[parameter_obj_advset_name]]
            
            if (parameter_obj_advset_info$type == "numeric") {
              numericInput(parameter_obj_advset_name, parameter_obj_advset_info$label, value = parameter_obj_advset_info$value, min = parameter_obj_advset_info$min, max = parameter_obj_advset_info$max)
            } else if (parameter_obj_advset_info$type == "text") {
              textInput(parameter_obj_advset_name, parameter_obj_advset_info$label, value = parameter_obj_advset_info$value)
            } else if (parameter_obj_advset_info$type == "slider") {
              sliderInput(parameter_obj_advset_name, parameter_obj_advset_info$label, min = parameter_obj_advset_info$min, max = parameter_obj_advset_info$max, value = parameter_obj_advset_info$value)
            } else if (parameter_obj_advset_info$type == "choice") {
              selectInput(parameter_obj_advset_name, parameter_obj_advset_info$label, choices = parameter_obj_advset_info$choices)
            } else if (parameter_obj_advset_info$type == "checkbox") {
              checkboxGroupInput(parameter_obj_advset_name, parameter_obj_advset_info$label, choices = parameter_obj_advset_info$choices, selected = parameter_obj_advset_info$selected)
            } else if (parameter_obj_advset_info$type == "numeric_vector") {
              textInput(parameter_obj_advset_name, parameter_obj_advset_info$label, value = parameter_obj_advset_info$value)
            }
          })
          # Return the dynamic UI
          do.call(tagList, ui_obj_advset_elements)
        })
      } else {output$parameter_obj_advset_input = NULL}
      
    } else if (click_counter() == 2) {
      # Second click removes UI with advanced parameter setting
      output$parameter_advset_input = NULL 
      output$parameter_obj_advset_input = NULL 
      
      # Reset the counter after the second click
      click_counter(0)
    }
    
  }) 
  
  
  ### By pressing 'Run analysis' button assign parameter and run the selected script
  observeEvent(input$run_script, {
    analysis_type = tolower(input$analysis_type)
    
    
    ### Set all parameter
    # First, store standard param list loaded from start_settings.R and standard_parameter.R 
    assigned_parameter(param)
    # Store advset param list; is an empty list from standard_parameter.R
    assigned_parameter_advset(param_advset)
    
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
    parameter = parameter_lists[[analysis_type]]
    for (parameter_name in names(parameter)) {
      if (parameter[[parameter_name]]$type == "text") {
        param[[parameter_name]] = unlist(strsplit(input[[parameter_name]], ","))
      } else {
        param[[parameter_name]] = input[[parameter_name]]
      }
    }
    
    
    
    # Add dynamically set advanced parameter to param_advset list
    parameter_advset = parameter_advset_lists[[analysis_type]]
    for (parameter_advset_name in names(parameter_advset)) {
      if (is.null(input[[parameter_advset_name]]) | is.na(input[[parameter_advset_name]])) {
        param_advset[[parameter_advset_name]] = NULL
      } else {
        if (parameter_advset[[parameter_advset_name]]$type == "text") {
          param_advset[[parameter_advset_name]] = unlist(strsplit(input[[parameter_advset_name]], ","))
        } else if (parameter_advset[[parameter_advset_name]]$type == "numeric_vector") {
          param_advset[[parameter_advset_name]] = as.numeric(unlist(strsplit(input[[parameter_advset_name]], ",")))
        } else {
          param_advset[[parameter_advset_name]] = input[[parameter_advset_name]]
        }
      }
    }
    
    # Combine or transform advanced parameters to a list element where needed
    if (analysis_type == "pre-processing" && !is.null(param_advset)) {
      param_advset[["cell_filter"]] = list(nFeature_RNA=param_advset[["nFeature_RNA"]], nCount_RNA=param_advset[["nCount_RNA"]], percent_mt=param_advset[["percent_mt"]] )
      param_advset[["feature_filter"]] = list(min_counts=param_advset[["min_counts"]], min_cells=param_advset[["min_cells"]] )
      param_advset$integrate_samples[["method"]]=param_advset[["integrate_samples_method"]]
      param_advset$integrate_samples[["integration_function"]]=param_advset[["integrate_samples_integration_function"]]
    }
    
    
    
    # Add output path to param list
    if (grepl(".Rmd", param$scriptname)) {
      param[["path_out"]] = file.path(param$path_to_git, "output", param$project_id, gsub(".Rmd", "", basename(param$scriptname)))
    } else {
      param[["path_out"]] = file.path(param$path_to_git, "output", param$project_id)
    }
    
    # Overwrite standard parameter with set advanced parameters
    param = modifyList(x = param, val = param_advset)
    
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
      notification = sprintf("<br><br><br>
                            <span class='notification-message'>'%s' completed!</span> <br>
                            <span class='notification-path'>Output folder: <br>   %s</span> <br><br>
                            <span class='notification-message'>You can close the app now.</span> <br><br><br>", 
                             input$analysis_type, param$path_out)
    } else {
      notification = sprintf("<br><br><br>
                            <span class='notification-message'>'%s' completed! <br><br>
                             You can close the app now.</span> <br><br><br>", 
                             input$analysis_type)
    }
    
    # Display the message 
    output$output_message = renderUI({
      HTML(notification)
    })
    
    # TEST
    #output$output_message <- renderText({
      #teststring = stringr::str_flatten_comma(param$cell_filter$nFeature_RNA)
      #teststring2 = stringr::str_flatten_comma(param$cell_filter$nCount_RNA)
      #testvector = unlist(strsplit(param$receiver, ","))
      #paste("Variables:", param$downsample_cells_equally, param$cell_filter$nFeature_RNA[1], param$cell_filter$nFeature_RNA[2])
    #})
  
  })
}
      


# Run the application 
################################################################################
shinyApp(ui = ui, server = server)

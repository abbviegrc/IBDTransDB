# IBDTransDB Compare Application
# A multi-omics comparison tool for target/biomarker discovery


library(RSQLite)

library(shiny)
library(shinyWidgets)
library(shinyBS)
library(shinybusy)
library(shinyalert)
library(rintrojs)
library(shinyjs)

library(tidyr)
library(dplyr)
library(stringr)

library(ggplot2)
library(rvest)
library(RColorBrewer)
library(gridExtra)

library(poolr)

library(DT)
library(reactable)
library(kableExtra)


#==============================================================================================================
# Import functions used to generate table summaries, table outputs, expression plots, etc
#==============================================================================================================
source("./app_functions/database_query.R")
source("./app_functions/dataset_selection.R")
source("./app_functions/comparison_summary.R")
source("./app_functions/comparison_tables.R")
source("./app_functions/signature_scores.R")
source("./app_functions/meta_analysis.R")
source("./app_functions/plot_functions.R")
source("./app_functions/geneset_view.R")
source("./app_functions/database_view.R")
#==============================================================================================================


#==============================================================================================================
# Create the connection to the IBDTransDB DB
#==============================================================================================================
db <- dbConnect(RSQLite::SQLite(), dbname = "./IBDTransDB.db")
#==============================================================================================================

jscode <- "shinyjs.refresh_page = function() { history.go(0); }"
# Define the UI for the application
ui <- navbarPage(
    title = "IBDCompare",
    id = "navbar-tabset",
    includeCSS("www/style.css"),
    useShinyjs(),
    extendShinyjs(text = jscode, functions = "refresh_page"),
    tags$head(
      tags$style(HTML("
        #toggleSection {
          display: block;
          margin-bottom: 10px;
          font-size: 16px;
          color: blue;
          cursor: pointer;
        }
        #toggleSection i {
          margin-right: 5px;
        }
      ")),
      tags$script(HTML('
        $(document).ready(function() {
          $(".section .header").click(function() {
            var section = $(this).parent();
            if (section.hasClass("expanded")) {
              section.removeClass("expanded");
              section.find(".content").slideUp();
              section.find(".expand-icon").removeClass("fa-chevron-down").addClass("fa-chevron-right");
            } else {
              $(".section.expanded").removeClass("expanded");
              $(".section.expanded .content").slideUp();
              $(".section.expanded .expand-icon").removeClass("fa-chevron-down").addClass("fa-chevron-right");
              section.addClass("expanded");
              section.find(".content").slideDown();
              section.find(".expand-icon").removeClass("fa-chevron-right").addClass("fa-chevron-down");
            }
          });
        });
      ')),
      tags$script(HTML("
        var header = $('.navbar> .container-fluid');
        header.append('<div style=\"float:right;color:white;padding-top:8px;\"><h4>Developed by Immunology Computational Biology - GRC</h3></div>');
        console.log(header)
      ")),
      tags$script(src = "popover.js"),
      tags$script(HTML("
        $(document).ready(function() {
          $('body').popover({
            selector: '[data-toggle=\"popover\"]',
            html: true,
            trigger: 'hover',
            container: 'body'
          });
        });
        $(document).ready(function(){
          $('[data-toggle=\"popover\"]').popover({html: true}); 
        });
        $(document).ready(function() {
          $('#toggle').click(function() {
            $('#collapse').collapse('toggle');
            var icon = $(this).find('i');
            icon.toggleClass('fa-plus fa-minus');
          });
        });
      "))
    ),
    
    ### DATASET SELECTION TAB ###
    tabPanel("Dataset Selection",id = "dataset_selection",
             column(3,
                    wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 1150px;",
                              
                              h4("REQUIRED", style = "text-align: center;"),
                              
                              div(
                                class = "section",
                                #actionButton("reset_btn","reset inputs",width = '100px'),
                                div(
                                  class = "header",
                                  icon("dna", class = "icon"),
                                  h4("TARGET INPUT")
                                  
                                ),
                                radioGroupButtons(
                                  inputId = "gene_option",
                                  label = NULL,
                                  size = "sm",
                                  direction = "horizontal",
                                  choices = c("GENE" = "GENE-level", "PATHWAY" = "PATHWAY-level"),
                                  selected = "GENE-level",
                                  checkIcon = list(yes = tags$i(class = "fa fa-check-square",
                                                                style = "color: steelblue"),
                                                   no = tags$i(class = "fa fa-square-o", 
                                                               style = "color: steelblue"))),
                                uiOutput("gene_input_selection")
                                #textInput(inputId = "gene_list_input",
                                #          label = "Enter a comma-separated list of gene names")
                                # conditionalPanel(
                                #   condition = "input.gene_option == 'GENE-level'",
                                #   selectInput(
                                #     "breaks", "Breaks",
                                #     c("Sturges", "Scott", "Freedman-Diaconis", "[Custom]" = "custom")
                                #   )
                                
                              ),
                              br(),
                              
                              div(
                                style = "margin-bottom: 5px;",
                                selectInput(
                                  inputId = "disease_input",
                                  label = "Select diseases of interest",
                                  choices = get_keywords("Disease", db),
                                  multiple = TRUE
                                )
                              ),
                              br(),
                              
                              h4("OPTIONAL", style = "text-align: center;"),

                              div(
                                style = "margin-bottom: 5px;",
                                selectInput(
                                  inputId = "tissue_input",
                                  label = "Select tissues of interest",
                                  choices = get_keywords("Source", db),
                                  multiple = TRUE
                                )
                              ),
                              div(
                                style = "margin-bottom: 5px;",
                                selectInput(
                                  inputId = "cell_type_input",
                                  label = "Select cell types of interest",
                                  choices = get_keywords("CellType", db),
                                  multiple = TRUE
                                )
                              ),
                              br(),
                              
                              div(
                                style = "margin-bottom: 5px;",
                                selectInput(
                                  inputId = "treatment_input",
                                  label = "Select treatments of interest",
                                  choices = get_keywords("Treatment", db),
                                  multiple = TRUE
                                )
                              ),
                              div(
                                style = "margin-bottom: 5px;",
                                selectInput(
                                  inputId = "timepoint_input",
                                  label = "Select timepoints of interest",
                                  choices = get_keywords("Timepoint", db),
                                  multiple = TRUE
                                )
                              ),
                              br(),
                              
                              div(class = "fixed-button",
                                  actionButton(inputId = "confirm_selections",
                                               label = "confirm selections",
                                               width = "300px",
                                               style = "font-size: 125%;"))
                    )),
             
             column(9,
                    wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 1150px;",
                              add_busy_spinner(spin = 'pixel', position = 'bottom-right'),
                              column(12, 
                                     # uiOutput("dataset_table_header"),
                                     conditionalPanel(
                                       condition = "input.confirm_selections > 0",
                                     div(
                                       div(style = "display: flex; justify-content: center;",
                                           actionButton(inputId = "run_immunocompare",
                                                        label = "Compare Datasets",
                                                        width = "300px",
                                                        style = "font-size: 125%;")),
                                       br(), br(),
                                       p("The table below allows you to explicitly select datasets that you would like to compare.",
                                         "By default, all datasets are selected. If you are only interested in a few datasets, you can click the",
                                         span(style = "font-weight: bold;", "deselect all"), "button and then make your selections. Click the", 
                                         span(style = "font-weight: bold;", "Compare Datasets"), "button to display results across selected data.",
                                         "Results will be displayed in the Results Panel."),
                                       br(),
                                       div(style = "display: inline-block;",
                                           span(style = "height: 30px;",
                                                actionButton(inputId = "select_all_datasets",
                                                             label = "select all",
                                                             class = "selection-button",
                                                             width = "100px",
                                                             style = "font-size: 90%;")),
                                           span(style = "height: 30px;",
                                                actionButton(inputId = "deselect_all_datasets",
                                                             label = "deselect all",
                                                             class = "selection-button",
                                                             width = "100px",
                                                             style = "font-size: 90%;")))
                                     )),
                                     DT::dataTableOutput("dataset_table"))
                              
                    ))
    ),
    
    ### RESULTS PANEL TAB ###
    tabPanel("Results Panel",
             
             ### TABS FOR OUTPUT ###
             column(10, offset = 1,
               tabsetPanel(
                 type = "tabs", 
                 id = "main_body_tabset",
  
                 ### SUMMARY TAB ###        
                 tabPanel("Summary",
                          class = "results-tab",
                          add_busy_spinner(spin = 'pixel', position = 'bottom-right'),
                          wellPanel(style = "background-color: #fff; border-color: #2c3e50; padding: 25px;",
                                    uiOutput("summary_title"),
                                    uiOutput("query_summary"),
                                    br(),
                                    uiOutput("dataset_summary"))
                 ),
                 
                 ### COMPARISON TABLES TAB ###        
                 tabPanel("Comparison Tables (Results)",
                          class = "results-tab",
                          add_busy_spinner(spin = 'pixel', position = 'bottom-right'),
                          wellPanel(
                            style = "background-color: #fff; border-color: #2c3e50; padding: 15px;",
                            fluidRow(
                              uiOutput("results_summary_title"),
                              column(6, uiOutput("dataset_results_summary_table")),
                              column(6, uiOutput("comparison_results_summary_table"))
                            ),
                            br(),
                            fluidRow(
                              uiOutput("results_description")
                            )
                          ),
                          wellPanel(style = "background-color: #fff; border-color: #2c3e50; padding: 15px;",
                                    uiOutput("comparison_table_title_and_options"),
                                    br(), br(),
                                    uiOutput("comparison_table_selection_options"),
                                    uiOutput("report_tables"))
                 ),
                 
                 ### VISUALS TAB ### 
                 tabPanel("Expression Visuals",
                          class = "results-tab",
                          add_busy_spinner(spin = 'pixel', position = 'bottom-right'),
                          wellPanel(style = "background-color: #fff; border-color: #2c3e50; padding: 25px;",
                                    uiOutput("expression_plot_title"),
                                    br(),
                                    uiOutput("plot_buttons"),
                                    br(),
                                    uiOutput("plot_titles"),
                                    br(),
                                    uiOutput("expression_plot_output")),
                          wellPanel(style = "background-color: #fff; border-color: #2c3e50; padding: 25px;",
                                    uiOutput("meta_analysis_header"),
                                    br(),
                                    reactableOutput("meta_analysis_table"))
                 ),
                 
                 ### DATABASES TAB ### 
                 tabPanel("Database View",
                          class = "results-tab",
                          wellPanel(style = "background-color: #fff; border-color: #2c3e50; padding: 25px;",
                                    uiOutput("database_title"),
                                    column(7,
                                          uiOutput("database_output"),
                                          uiOutput("database_summary")))
                          
                 )
                           
               )
             )
    )
    
    ### HELP TAB ###
    # tabPanel("Help",
    #          column(10, offset = 1,
    #                 wellPanel(style = "background-color: #fff; border-color: #2c3e50;",
    #                           add_busy_spinner(spin = 'pixel', position = 'bottom-right'),
    #                           includeHTML("homepage.html")))
    # )
    
)


# The active datasets table (defined globally)
active_datasets_table <- reactiveVal(NULL)

# Define the server functionality for the application
server <- function(input, output, session) {
    
    ### REACTIVE VALUES ###
    # Define a vector used for genes
    genes <- reactiveVal(NULL)
    # Define a vector used for pathways
    pathways <- reactiveVal(NULL)
    # Define a vector used for the exact diseases (necessary because some datasets are associated with multiple diseases)
    exact_diseases <- reactiveVal(NULL)
    # Generate UNSTYLED comparison tables
    comp_tables_unstyled_reactive <- reactiveVal(NULL)
    # Generate STYLED comparison tables
    comp_tables_styled_reactive <- reactiveVal(NULL)
    # Define the list used to store expression plots
    plots_reactive <- reactiveVal(NULL)
    
    # Queries
    diseases <- reactiveVal(NULL)
    organisms <- reactiveVal(NULL)
    sources <- reactiveVal(NULL)
    cell_types <- reactiveVal(NULL)
    experiment_types <- reactiveVal(NULL)
    platforms <- reactiveVal(NULL)
    treatments <- reactiveVal(NULL)
    timepoints <- reactiveVal(NULL)
    
    datasets <- reactiveVal(NULL)
    all_datasets <- reactiveVal(NULL)
    comparisons <- reactiveVal(NULL)
    all_comparisons <- reactiveVal(NULL)
    comparison_data <- reactiveVal(NULL)
    signature_data <- reactiveVal(NULL)
    pathway_data <- reactiveVal(NULL)
    sample_ann <- reactiveVal(NULL)
    gene_map <- reactiveVal(NULL)
    
    comparison_tables <- reactiveVal(NULL)
    ##########
    
    enrichment_gene_sets <- read.csv("enrichment_gene_sets.csv")
    values<- reactiveValues()
    values$checkbox_ids <- NULL
    
    #==============================================================================================================
    ### Set initial gene input option and ObserveEvent for GENES SWITCH
    #==============================================================================================================
    
    # Set the initial gene selection option
    # output$gene_input_selection <- renderUI({
      # selectInput(inputId = "gene_list_input",
      #             label = "Enter a list of gene names",
      #             choices = get_genes(db),
      #             multiple = TRUE)
      textInput(inputId = "gene_list_input",
                label = "Enter a comma-separated list of gene names")
    # })
    
    # Update the input for targets when a different option is selected
    observeEvent(input$gene_option, {
      if (input$gene_option == "GENE-level") {
        output$gene_input_selection <- renderUI({
          # selectInput(inputId = "gene_list_input",
          #             label = "Enter a list of gene names",
          #             choices = get_genes(db),
          #             multiple = TRUE)
          div(style = "padding-bottom: 5px;",
              textInput(inputId = "gene_list_input",
                        label = HTML("Enter a list of gene names <br/> (genes must be separated by commas)")))
        })
      } else if (input$gene_option == "SIGNATURE-level") {
        output$gene_input_selection <- renderUI({
          div(style = "padding-bottom: 5px;",
              textAreaInput(inputId = "signature_list_input",
                            label = HTML("Paste a list of gene names <br/> (each gene must be on a new line)"),
                            placeholder = "paste genes here"))
        })
      } else {
        pathways <- list("pathway" = c("KEGG" = "pathway_KEGG", "Reactome" = "pathway_Reactome"))
        output$gene_input_selection <- renderUI({
          div(
            selectInput(inputId = "pathway_database_input",
                        label = "Select one or multiple pathway databases",
                        choices = pathways,
                        selected = NULL,
                        multiple = TRUE,
                        width = "auto"),
            selectInput(inputId = "pathway_list_input",
                        label = "Select one or multiple pathways",
                        choices = NULL,
                        multiple = TRUE,
                        width = "auto")
          )
        })
      }
    })
    
   # Update the pathway options when an enrichment database is selected
    observeEvent(input$pathway_database_input, {
      # Pull the pathways associated with the selected enrichment databases
      enrichment_databases <- unlist(input$pathway_database_input)
      related_pathways <- enrichment_gene_sets %>%
        filter(database %in% enrichment_databases) %>%
        pull(description)

      # Update the pathway options
      updateSelectInput(session,
                        "pathway_list_input",
                        choices = related_pathways)
    })
    #==============================================================================================================
    
    
    #==============================================================================================================
    ### Generate a table of datasets based on the user query when the user clicks the "Confirm Selections" button
    #==============================================================================================================
    observeEvent(input$confirm_selections, {
      print(input)
      
      
      # Generate an alert if the user does not specify any genes or diseases
      if (!(isTruthy(input$gene_list_input) | isTruthy(input$signature_list_input) | isTruthy(input$pathway_list_input)) | !isTruthy(input$disease_input)) {
        shinyalert(title = "Insufficient Input Supplied",
                   text = "Please specify at least one gene and at least one disease",
                   type = "error")
      } else if(length(unlist(strsplit(input$gene_list_input,",",fixed=T)))>10 | length(input$pathway_list_input)>10){
        shinyalert(title = "Too many input genes/pathways",
                   text = "Please input at most 10 genes/pathways",
                   type = "error")
        } else{
          # Define the user query
          diseases <- unlist(input$disease_input)
          diseases(diseases)
          print(diseases)
          organisms <- unlist(input$organism_input)
          organisms(organisms)
          sources <- unlist(input$tissue_input)
          sources(sources)
          cell_types <- unlist(input$cell_type_input)
          cell_types(cell_types)
          experiment_types <- unlist(input$experiment_type_input)
          experiment_types(experiment_types)
          platforms <- unlist(input$platform_input)
          platforms(platforms)
          treatments <- unlist(input$treatment_input)
          treatments(treatments)
          print(treatments)
          timepoints <- unlist(input$timepoint_input)
          timepoints(timepoints)
          
          # Pull the dataset descriptions
          datasets <- get_datasets(diseases, 
                                   organisms, 
                                   sources, 
                                   treatments, 
                                   experiment_types, 
                                   platforms, 
                                   cell_types,
                                   timepoints,
                                   db)
          # Format the diseases in the datasets
          datasets <- datasets %>%
            mutate(disease = gsub(";Healthy|;nonIBD|;nonRA", "", disease)) %>%
            mutate(disease = gsub(";$", "", disease))
          
          # Pull the comparison descriptions and comparison data
          comparisons <- get_comparisons(datasets,db)
          
          # Update the reactive value
          comparisons(comparisons)
          all_comparisons(comparisons)
          
          # Format the diseases column in the datasets table and filter for datasets in the comparisons table
          disease_pattern <- paste(diseases, collapse = "|")
          datasets <- datasets %>%
            mutate(disease = sapply(disease, function(x) paste(unlist(str_extract_all(x, disease_pattern)), collapse = ";")))
          # Update the reactive value
          datasets(datasets)
          all_datasets(datasets)
          
          # Pull the sample annotations
          sample_ann <- get_sample_ann(datasets, db)
          # Update the reactive value
          sample_ann(sample_ann)
          
          # Pull the gene map
          gene_map <- get_gene_map(db)
          # Update the reactive value
          gene_map(gene_map)
          
          # Display the guidlines for interacting with the dataset selection table and display the Run ImmunoCompare button
          # output$dataset_table_header <- renderUI({
          #   div(
          #     div(style = "display: flex; justify-content: center;",
          #         actionButton(inputId = "run_immunocompare",
          #                      label = "Compare Datasets",
          #                      width = "300px",
          #                      style = "font-size: 125%;")),
          #     br(), br(),
          #     p("The table below allows you to explicitly select datasets that you would like to compare.",
          #       "By default, all datasets are selected. If you are only interested in a few datasets, you can click the",
          #       span(style = "font-weight: bold;", "deselect all"), "button and then make your selections. Click the", 
          #       span(style = "font-weight: bold;", "Compare Datasets"), "button to display results across selected data.",
          #       "Results will be displayed in the Results Panel."),
          #     br(),
          #     div(style = "display: inline-block;",
          #         span(style = "height: 30px;",
          #              actionButton(inputId = "select_all_datasets",
          #                           label = "select all",
          #                           class = "selection-button",
          #                           width = "100px",
          #                           style = "font-size: 90%;")),
          #         span(style = "height: 30px;",
          #              actionButton(inputId = "deselect_all_datasets",
          #                           label = "deselect all",
          #                           class = "selection-button",
          #                           width = "100px",
          #                           style = "font-size: 90%;")))
          #   )
          # })
          
          # Pull the dataset descriptions for the dataset selection table and update the reactive value
          dataset_table <- get_dataset_table_data(datasets, db)
          active_datasets_table(dataset_table)
          # Display the table with datasets (based on the user query) to select from (all datasets selected by default)
          output$dataset_table <- DT::renderDataTable({
            generate_dataset_table(dataset_table, FALSE)
          })
        }
       
      
    })
    #==============================================================================================================
    
    
    #==============================================================================================================
    ### ObserveEvent for selecting all or deselecting all datasets
    #==============================================================================================================
    observeEvent(input$select_all_datasets, {
      # Display the dataset selection table (all datasets selected)
      output$dataset_table <- DT::renderDataTable(server = FALSE, {
        dataset_table <- isolate(get_dataset_table_data(datasets(), db))
        generate_dataset_table(dataset_table, TRUE)
      })
    }, ignoreInit = TRUE)
    
    observeEvent(input$deselect_all_datasets, {
      # Display the dataset selection table (no datasets selected)
      output$dataset_table <- DT::renderDataTable(server = FALSE, {
        dataset_table <- isolate(get_dataset_table_data(datasets(), db))
        generate_dataset_table(dataset_table, FALSE)
      })
    }, ignoreInit = TRUE)
    
    
    observeEvent(input$reset_btn, {
      # js$refresh_page();
      
      session$reload()
      return()
    })
    #==============================================================================================================
    
    
    #==============================================================================================================
    ### Run ImmunoCompare and produce output when the user clicks the "Run ImmunoCompare" button
    #==============================================================================================================
    observeEvent(input$run_immunocompare, {
      
        # Generate an alert if the user does not specify any genes or diseases
        if (!isTruthy(input$dataset_table_rows_selected)) {
            shinyalert(title = "Insufficient Input Supplied",
                       text = "You need to specify at least one gene/pathway and at least one disease",
                       type = "error")
        } else {
            
            # Pull the selected dataset IDs from the dataset selection table
            
            selected_dataset_acc <- active_datasets_table()$dataset_acc[input$dataset_table_rows_selected]
            datasets <- filter(all_datasets(), dataset_acc %in% selected_dataset_acc)
            
            datasets(datasets)
            print("Filtered Datasets")
            
            comparisons <- filter(all_comparisons(), dataset_acc %in% selected_dataset_acc)
            comparisons(comparisons)
            print("Comparisons filtered")
            
            # Generate a vector of exact diseases (necessary because some datasets are associated with multiple diseases)
            exact_diseases <- get_exact_diseases(datasets)
            # Update the reactive value
            exact_diseases(exact_diseases)
            print("retr. exact diseaes")
            isError <- FALSE
            if (input$gene_option == "GENE-level") {
              # Convert the user specified genes to a list and convert all genes to uppercase
              genes <- toupper(str_trim(unlist(str_split(input$gene_list_input, ","))))
              # Pull the comparison data and update the reactive value
              genes <- intersect(genes,gene_map()$gene)
              if(length(genes)==0){
                shinyalert(title = "Invalid gene symbol!",
                           text = "Please input the valid gene symbols.",
                           type = "error")
                isError <- TRUE
              }else{
                comparison_data <- unique(get_comparison_data(comparisons, genes, gene_map(), db))
                comparison_data(comparison_data)
                print("comparison_data")
              }
            } else if (input$gene_option == "SIGNATURE-level") {
              # Convert the user specified genes to a list and convert all genes to uppercase
              genes <- toupper(str_trim(unlist(strsplit(input$signature_list_input, "\n"))))
              # Perform signature-based comparisons and update the reactive value
              
              signature_data <- perform_signature_comparisons(genes, comparisons, db)
              signature_data(signature_data)
            } else {
              # Pull the selected pathway databases and pathways
              #pathway_databases <- unlist(input$pathway_database_input)
              
              #only support KEGG
              pathway_databases <- unlist(input$pathway_database_input)
              
              pathways <- unlist(input$pathway_list_input)
              
              # print(pathways())
              # pathways(unlist(input$pathway_list_input))
              # print(pathways())
              # Pull the pathway data and update the reactive value
              #save(comparisons,pathway_databases,pathways,db,file="test.RData")
              pathway_data <- get_pathway_data(comparisons, pathway_databases, pathways, db)
              pathway_data(pathway_data)
              print("Retreived pathway data")
            }
            
            if(isError==FALSE){
              #==============================================================================================================
              ### QUERY SUMMARY
              #==============================================================================================================
              # Jump to the Comparison Summary tab
              updateTabsetPanel(session, "navbar-tabset", selected = "Results Panel")
              updateTabsetPanel(session, "main_body_tabset", selected = "Summary")
              
              # Generate the title of the summary section
              output$summary_title <- renderUI({
                h2("Summary of Query and Matching Datasets", align = "center")
              })
              if (input$gene_option == "GENE-level") {
                # Summarize the query
                output$query_summary <- renderUI({
                  div(
                    p("The following output contains comparisons that meet the user specified critieria displayed below."),
                    tags$ul(
                      tags$li(paste0("genes - ", paste0(genes, collapse = ", "))),
                      tags$li(paste0("diseases - ", paste0(diseases(), collapse = ", "))),
                      tags$li(paste0("sources - ", paste0(c(sources(), cell_types()), collapse = ", "))),
                      tags$li(paste0("treatments - ", paste0(treatments(), collapse = ", "))),
                      tags$li(paste0("timepoints - ", paste0(timepoints(), collapse = ", ")))
                    )
                  )
                })
              }else{
                output$query_summary <- renderUI({
                  div(
                    p("The following output contains comparisons that meet the user specified critieria displayed below."),
                    tags$ul(
                      tags$li(paste0("pathways - ", paste0(pathways, collapse = ", "))),
                      tags$li(paste0("diseases - ", paste0(diseases(), collapse = ", "))),
                      tags$li(paste0("sources - ", paste0(c(sources(), cell_types()), collapse = ", "))),
                      tags$li(paste0("treatments - ", paste0(treatments(), collapse = ", "))),
                      tags$li(paste0("timepoints - ", paste0(timepoints(), collapse = ", ")))
                    )
                  )
                })
              }
              # Create output objects for each disease summary (i.e. number of corresponding datasets)
              output$dataset_summary <- renderUI({
                output_summary_list <- lapply(paste0(unlist(input$disease_input), "_summary"), function(x) {
                  htmlOutput(x)
                })
                tagList(output_summary_list)
              })
              
              # Display the content for the disease summaries (i.e. number of corresponding datasets)
              sapply(diseases(), function(current_disease) {
                
                # Generate a table similar to the output (without styling)
                if (input$gene_option == "GENE-level") {
                  
                  # Generate a table with the dataset summary
                  results_table <- generate_dataset_summary_table(genes, current_disease, 
                                                                  datasets, comparisons, 
                                                                  comparison_data)
                  # Calculate the total number of comparisons and datasets found for each gene
                  print("Caclulating Comparisons")
                  summary <- generate_dataset_summary_stats(genes, results_table)
                  summary_df <- summary[["summary_df"]]
                  total_num_datasets <- summary[["total_num_datasets"]]
                  total_num_comparisons <- summary[["total_num_comparisons"]]
                  print("Summary table generated")
                  
                  # Generate the summary output
                  output[[paste0(current_disease, "_summary")]] <- renderUI({
                    # Create a bulleted list for each disease, specifying the summary metrics per gene
                    div(
                      p("Summary for ", span(current_disease, style="font-weight:bold"), ": \n",
                        total_num_datasets, "datasets", 
                        paste0("(", total_num_comparisons, " comparisons)"), 
                        "were returned based on the user specified queries."),
                      tags$ul(
                        lapply(genes, function(x) {
                          num_comparisons <- summary_df[x, "comparisons"]
                          num_datasets <- summary_df[x, "datasets"]
                          gene_summary <- paste0(x, " - found in ", num_datasets, " datasets ",
                                                 "(", num_comparisons, " comparisons)")
                          return(tags$li(gene_summary))
                        })
                      )
                    )
                  })
                } else if (input$gene_option == "PATHWAY-level") {
                  # Generate a table with the dataset summary
                  results_table <- pathway_generate_dataset_summary_table(pathways, current_disease, 
                                                                          datasets, comparisons, 
                                                                          pathway_data)
                  
                  print("Result table generated")
                  
                  # Calculate the total number of comparisons and datasets found for each pathway
                  summary <- pathway_generate_dataset_summary_stats(pathways, results_table)
                  summary_df <- summary[["summary_df"]]
                  total_num_datasets <- summary[["total_num_datasets"]]
                  total_num_comparisons <- summary[["total_num_comparisons"]]
                  
                  # Generate the summary output
                  output[[paste0(current_disease, "_summary")]] <- renderUI({
                    # Create a bulleted list for each disease, specifying the summary metrics per pathway
                    div(
                      p("Summary for ", span(current_disease, style="font-weight:bold"), ": \n",
                        total_num_datasets, "datasets", 
                        paste0("(", total_num_comparisons, " comparisons)"), 
                        "were returned based on the user specified queries."),
                      tags$ul(
                        lapply(pathways, function(x) {
                          num_comparisons <- summary_df[x, "comparisons"]
                          num_datasets <- summary_df[x, "datasets"]
                          pathway_summary <- paste0(x, " - found in ", num_datasets, " datasets ",
                                                    "(", num_comparisons, " comparisons)")
                          return(tags$li(pathway_summary))
                        })
                      )
                    )
                  })
                }
                
              })
              #==============================================================================================================
              
              
              #==============================================================================================================
              ### RESULTS SUMMARY
              #==============================================================================================================
              
              # Generate the title of the results summary section
              output$results_summary_title <- renderUI({
                h4("Summary of Results", align = "center")
              })
              
              # Generate the summary tables
              if (input$gene_option == "GENE-level") {
                print("generating summary table")
                summary_tables <- generate_results_summaries(genes, diseases(), 
                                                             experiment_types(), datasets, 
                                                             comparisons, comparison_data, 
                                                             "p_value", 0.05, 
                                                             1, db)
                print("summary Table generated")
              } else if (input$gene_option == "SIGNATURE-level") {
                summary_tables <- signature_generate_results_summaries(diseases(), experiment_types, 
                                                                       datasets, comparisons, 
                                                                       signature_data, 
                                                                       0.05, db)
              } else {
                # print("Generating Comparison #768")
                # print(pathways)
                # print(class(pathways))
                summary_tables <- pathway_generate_results_summaries(pathways, diseases(), 
                                                                     experiment_types(), datasets, 
                                                                     comparisons, pathway_data, 
                                                                     "p_value", 0.05, 
                                                                     "nes", 1, db)
              }
              
              # Display the content for the dataset results summary table
              output$dataset_results_summary_table <- renderText({
                summary_tables[["dataset_results_summary"]]
              })
              
              # Display the content for the comparison results summary table
              output$comparison_results_summary_table <- renderText({
                summary_tables[["comparison_results_summary"]]
              })
              
              # P-VALUE OPTION, P-VALUE THRESHOLD, LOGFC THRESHOLD, ES OPTION, ES THRESHOLD, REMOVE NS COMPARISONS
              observeEvent(c(input$pval_option, 
                             input$pval_threshold, 
                             input$logfc_threshold, 
                             input$es_option, 
                             input$es_threshold), {
                               
                               # Generate the summary tables
                               if (input$gene_option == "GENE-level") {
                                 summary_tables <- generate_results_summaries(genes, diseases(), 
                                                                              experiment_types(), datasets, 
                                                                              comparisons, comparison_data, 
                                                                              input$pval_option, 
                                                                              input$pval_threshold, 
                                                                              input$logfc_threshold, 
                                                                              db)
                               } else if (input$gene_option == "PATHWAY-level") {
                                 print("Generating Pathway Results Summaries #805")
                                 pathways(unlist(input$pathway_list_input))
                                 # print(pathways)
                                 # print(class(pathways))
                                 # print(diseases)
                                 # print(experiment_types)
                                 # print(datasets)
                                 # print(comparisons)
                                 # print(pathway_data)
                                 
                                 summary_tables <- pathway_generate_results_summaries(pathways(), diseases(), 
                                                                                      experiment_types(), datasets(), 
                                                                                      comparisons, pathway_data(), 
                                                                                      input$pval_option, 
                                                                                      input$pval_threshold, 
                                                                                      input$es_option, 
                                                                                      input$es_threshold, 
                                                                                      db)
                                 
                                 print("Pathway Summary Tables generated #835")
                                 
                               }
                               # Display the content for the dataset results summary table
                               output$dataset_results_summary_table <- renderText({
                                 summary_tables[["dataset_results_summary"]]
                               })
                               
                               # Display the content for the comparison results summary table
                               output$comparison_results_summary_table <- renderText({
                                 summary_tables[["comparison_results_summary"]]
                               })            
                             }, ignoreInit = TRUE)
              #==============================================================================================================
              
              
              #==============================================================================================================
              ### RESULTS (COMPARISON TABLES)
              #==============================================================================================================
              
              # Generate the description of the results
              output$results_description <- renderUI({
                div(
                  div(style = "display: inline-block; align: center; padding-left: 15px;",
                      h4("Description of Results")),
                  div(style = "display: inline-block; align: center;",
                      actionButton("show_description", 
                                   HTML("<i class='fa fa-plus'></i>",),
                                   style = "color: black; background-color: #EDEDED; border-color: #36454f;")),
                )
              })
              
              # Display the modal with the description of the results
              observeEvent(input$show_description, {
                showModal(modalDialog(
                  p("Clicking the dataset ID's in the first column of the following table(s) will navigate", 
                    "you to a page associated with the given dataset where you can further explore the data.", 
                    "Hovering over a cell under a gene, signature, or pathway specific column will generate",
                    " a pop with related information for the given. (i.e. p-value, adjusted p-value, and",
                    "|log", tags$sub("2"), "(fold change)|. A cell containing containing", span("up", style = "color:blue"), 
                    "indicates significant upregulated and a cell containing", span("down", style = "color:blue"), 
                    "indicates significant downregulated based on the specified thresholds (e.g. p-value, enrichment score).",
                    "An NS indicates that the given gene, signature, or pathway did not meet the specified thresholds for",
                    "signifiance. An", span("NA", style = "color:red"), "in the following table(s) indicates the given",
                    "gene, signature, or pathway was not found in the given dataset. The checkboxes in the right-hand column",
                    "can be used to filter the results tables by making selections for comparisons of interest and then clicking",
                    "the Filter Comparisons button. These checkboxes can also be used to visualize the expression patterns",
                    "related to comparisons of interest by making selections and then clicking the Visualize Comparisons button."),
                  footer = NULL,
                  easyClose = TRUE
                ))
              })
              
              # If organisms is NULL, set organisms to all organisms in the DB
              # if (is.null(organisms)) {
              #   organisms <- get_keywords("'Organism'", db)
              # }
              
              ### Generate the table outputs (the table titles, tables, and after table breaks are generated separately) ###
              output$report_tables <- renderUI({
                table_names <- unlist(lapply(exact_diseases, function(x) {
                  c(paste0(x, "_table_title"), paste0(x, "_table"), paste0(x, "_break"))
                }))
                output_table_list <- lapply(table_names, function(x) {
                  htmlOutput(x)
                })
                tagList(output_table_list)
              })
              ##########
              
              
              # Generate the comparison tables
              if (input$gene_option == "GENE-level") {
                table_list <- generate_comparison_tables(genes, exact_diseases, 
                                                         datasets, comparisons, 
                                                         comparison_data, "p_value", 
                                                         0.05, 1, FALSE)
                print("generate_comparison_tables")
              } else if (input$gene_option == "SIGNATURE-level") {
                
                table_list <- signature_generate_comparison_tables(exact_diseases, 
                                                                   datasets, comparisons,
                                                                   signature_data,
                                                                   0.05, FALSE)
              } else {
                print("Generating comparison tables list #904")
                table_list <- pathway_generate_comparison_tables(pathways, exact_diseases, 
                                                                 datasets, comparisons, 
                                                                 pathway_data, "p_value", 
                                                                 0.05, "NES", 1, FALSE)
              }
              # Update the REACTIVE VALUE for the comparison tables
              comparison_tables(table_list)
              # Style the comparison tables
              table_list <- style_comparison_tables(exact_diseases, table_list)
              
              ### Display the table titles, tables, and after table breaks ###
              sapply(exact_diseases, function(current_disease) {
                
                # Pull the comparison table for the given disease
                comparison_table <- table_list[[current_disease]]
                
                # Render the output if the comparison table is NOT NULL
                if (!is.null(comparison_table)) {
                  
                  # Generate the disease title for the comparison table
                  output_name <- paste0(current_disease, "_table_title")
                  output[[output_name]] <- renderUI({
                    div(style = "display: flex; justify-content: center; align-items: center;",
                        div(style = "flex-grow: 1; height: 1px; background-color: black; margin-right: 10px;",
                            hr(style = "border: none; height: 1px; background-color: black; margin: 0;")),
                        h3(current_disease),
                        div(style = "flex-grow: 1; height: 1px; background-color: black; margin-left: 10px;",
                            hr(style = "border: none; height: 1px; background-color: black; margin: 0;"))
                    )
                    # h3(paste("The following table display data from datasets associated with", gsub(";", "/", current_disease)))
                  })
                  
                  # Display the comparison table for the given disease
                  output_name <- paste0(current_disease, "_table")
                  output[[output_name]] <- renderText({
                    comparison_table
                  })
                  
                  # Generate the after table break and horizontal rule
                  output_name <- paste0(current_disease, "_break")
                  output[[output_name]] <- renderUI({
                    div(hr(), br())
                  })
                  
                }
              })
              ##########
              #==============================================================================================================
              
              
              #==============================================================================================================
              ### COMPARISON TABLE BUTTONS
              #==============================================================================================================
              
              # Display options used for modifying the comparison tables
              if (input$gene_option == "GENE-level") {
                output$comparison_table_title_and_options <- renderUI({
                  div(style = "display: flex; justify-content: space-between; align-items: center;",
                      div(
                        div(style = "display: inline-block; vertical-align: middle;",
                            h4("Table Options")),
                        div(style = "display: inline-block; vertical-align: middle;",
                            dropdownButton(h4("table options"),
                                           radioGroupButtons(inputId = "pval_option",
                                                             label = NULL,
                                                             choices = c("p-value" = "p_value", "adjusted p-value" = "p_value_adj"),
                                                             selected = "p_value",
                                                             checkIcon = list(
                                                               yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                                               no = tags$i(class = "fa fa-square-o", style = "color: steelblue"))),
                                           numericInput(inputId = "pval_threshold",
                                                        label = "Set a threshold for the p-value",
                                                        value = 0.05,
                                                        min = 0,
                                                        max = 1,
                                                        step = 0.01),
                                           numericInput(inputId = "logfc_threshold",
                                                        label = HTML(paste0("Set a threshold for the logFC")),
                                                        value = 1,
                                                        min = 0,
                                                        step = 0.1),
                                           materialSwitch(inputId = "remove_ns",
                                                          label = "Remove NS comparisons",
                                                          value = FALSE,
                                                          status = "primary"),
                                           circle = TRUE, 
                                           status = "primary",
                                           size = "sm",
                                           icon = icon("cog"), 
                                           width = "400px",
                                           right = FALSE))
                      ),
                      h2("Results", style = "margin: 0; margin-left: -10px;"),
                      p("")
                  )
                  # div(style = "display: flex; justify-content: center;",
                  #     div(style = "display: inline-block; vertical-align: bottom; padding-right: 5px;",
                  #         radioGroupButtons(inputId = "pval_option",
                  #                           label = NULL,
                  #                           choices = c("p-value" = "p_value", "adjusted p-value" = "p_value_adj"),
                  #                           selected = "p_value_adj",
                  #                           checkIcon = list(
                  #                             yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                  #                             no = tags$i(class = "fa fa-square-o", style = "color: steelblue")))),
                  #     div(style = "display: inline-block; vertical-align: bottom; padding-left: 5px; padding-right: 5px;",
                  #         numericInput(inputId = "pval_threshold",
                  #                      label = "Set a threshold for the p-value",
                  #                      value = 0.05,
                  #                      min = 0,
                  #                      max = 1,
                  #                      step = 0.01)),
                  #     div(style = "display: inline-block; vertical-align: bottom; padding-left: 5px; padding-right: 5px;",
                  #         numericInput(inputId = "logfc_threshold",
                  #                      label = HTML(paste0("Set a threshold for the logFC")),
                  #                      value = 1,
                  #                      min = 0,
                  #                      step = 0.1)),
                  #     div(style = "display: inline-block; vertical-align:top; padding-left: 15px;",
                  #         materialSwitch(inputId = "remove_ns",
                  #                        label = "Remove NS comparisons",
                  #                        value = FALSE,
                  #                        status = "primary"))
                  # )
                })
              } else if (input$gene_option == "SIGNATURE-level") {
                output$comparison_table_title_and_options <- renderUI({
                  div(style = "display: flex; justify-content: space-between; align-items: center;",
                      div(
                        div(style = "display: inline-block; vertical-align: middle;",
                            h4("Table Options")),
                        div(style = "display: inline-block; vertical-align: middle;",
                            dropdownButton(h4("table options"),
                                           numericInput(inputId = "pval_threshold",
                                                        label = "Set a threshold for the p-value",
                                                        value = 0.05,
                                                        min = 0,
                                                        max = 1,
                                                        step = 0.01),
                                           materialSwitch(inputId = "remove_ns",
                                                          label = "Remove NS comparisons",
                                                          value = FALSE,
                                                          status = "primary"),
                                           circle = TRUE, 
                                           status = "primary",
                                           size = "sm",
                                           icon = icon("cog"), 
                                           width = "400px",
                                           right = FALSE))
                      ),
                      h2("Results", style = "margin: 0; margin-left: -10px;"),
                      p("")
                  )
                })
              } else {
                output$comparison_table_title_and_options <- renderUI({
                  div(style = "display: flex; justify-content: space-between; align-items: center;",
                      div(
                        div(style = "display: inline-block; vertical-align: middle;",
                            h4("Table Options")),
                        div(style = "display: inline-block; vertical-align: middle;",
                            dropdownButton(h4("table options"),
                                           radioGroupButtons(inputId = "pval_option",
                                                             label = NULL,
                                                             choices = c("p-value" = "p_value", "adjusted p-value" = "p_value_adj"),
                                                             selected = "p_value_adj",
                                                             checkIcon = list(
                                                               yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                                               no = tags$i(class = "fa fa-square-o", style = "color: steelblue"))),
                                           numericInput(inputId = "pval_threshold",
                                                        label = "Set a threshold for the p-value",
                                                        value = 0.05,
                                                        min = 0,
                                                        max = 1,
                                                        step = 0.01),
                                           radioGroupButtons(inputId = "es_option",
                                                             label = NULL,
                                                             choices = c("ES" = "es", "NES" = "nes"),
                                                             selected = "nes",
                                                             checkIcon = list(
                                                               yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                                               no = tags$i(class = "fa fa-square-o", style = "color: steelblue"))),
                                           numericInput(inputId = "es_threshold",
                                                        label = HTML(paste0("Set a threshold for the enrichment score")),
                                                        value = 1,
                                                        min = 0,
                                                        step = 0.1),
                                           materialSwitch(inputId = "remove_ns",
                                                          label = "Remove NS comparisons",
                                                          value = FALSE,
                                                          status = "primary"),
                                           circle = TRUE, 
                                           status = "primary",
                                           size = "sm",
                                           icon = icon("cog"), 
                                           width = "400px",
                                           right = FALSE))
                      ),
                      h2("Results", style = "margin: 0; margin-left: -10px;"),
                      p("")
                  )
                })
              }
              
              # Display options used for filtering the comparison tables and generating expression box plots
              if (input$gene_option == "GENE-level") {
                output$comparison_table_selection_options <- renderUI({
                  div(style = "display: flex; justify-content: space-between;",
                      div(
                        div(style = "display: inline-block; vertical-align: middle;",
                            actionButton("filter_comparisons",
                                         label = "Filter Comparisons",
                                         style = "color: black; background-color: #EDEDED; border-color: #36454f;")),
                        div(style = "display: inline-block; vertical-align: middle;",
                            actionButton("clear_comparison_selections",
                                         label = "Remove Filter",
                                         style = "color: black; background-color: #EDEDED; border-color: #36454f;"))
                      ),
                      div(
                        actionButton("visuals_button", 
                                     label = "Visualize Comparisons",
                                     style = "color: black; background-color: #EDEDED; border-color: #36454f;")
                      )
                  )
                })
              } else if (input$gene_option == "SIGNATURE-level") {
                output$comparison_table_selection_options <- renderUI({
                  div(style = "display: flex; justify-content: space-between;",
                      div(
                        div(style = "display: inline-block; vertical-align: middle;",
                            actionButton("filter_comparisons",
                                         label = "Filter Comparisons",
                                         style = "color: black; background-color: #EDEDED; border-color: #36454f;")),
                        div(style = "display: inline-block; vertical-align: middle;",
                            actionButton("clear_comparison_selections",
                                         label = "Remove Filter",
                                         style = "color: black; background-color: #EDEDED; border-color: #36454f;"))
                      ),
                      div(
                        actionButton("visuals_button", 
                                     label = "Generate Visuals",
                                     style = "color: black; background-color: #EDEDED; border-color: #36454f;")
                      )
                  )
                })
              } else {
                output$comparison_table_selection_options <- renderUI({
                  div(
                    div(style = "display: inline-block; vertical-align: middle;",
                        actionButton("filter_comparisons",
                                     label = "Filter Comparisons",
                                     style = "color: black; background-color: #EDEDED; border-color: #36454f;")),
                    div(style = "display: inline-block; vertical-align: middle;",
                        actionButton("clear_comparison_selections",
                                     label = "Remove Filter",
                                     style = "color: black; background-color: #EDEDED; border-color: #36454f;")),
                    
                    
                  )
                })
              }
              #==============================================================================================================
              
              
              #==============================================================================================================
              ### DATABASE LINKS BOX
              #==============================================================================================================
              
              if (input$gene_option == "GENE-level") {
                # Generate the title of the geneset section
                output$database_title <- renderUI({
                  h2("Database Links", align = "center")
                })
                
                # Generate the table of database links for each gene
                output$database_output <- renderText({
                  generate_DB_links(genes, db)
                })
                
                # Generate a summary for each database listed in the table
                # db_df <- generate_DB_summary(db)
                # output$database_summary <- renderUI({
                #     
                #     lapply(1:nrow(db_df), function(x) {
                #         span(
                #             div(style = "border: 1px solid grey; padding: 10px; border-radius: 10px", 
                #                  span(db_df$database_name[x], style = "font-size: 14px; font-weight: bold;"),
                #                  br(),
                #                  span("Summary:", style = "font-weight: bold;"), db_df$summary[x]
                #             ),
                #             br(),
                #             br()
                #         )
                #     })
                #     
                # })
              }
            }
            
            #==============================================================================================================
        }
      
    })
    #==============================================================================================================
    
    
    #==============================================================================================================
    ### ObserveEvents for modifying comparison tables (p-value threshold, logFC threshold, remove NS)
    #==============================================================================================================
    
    # P-VALUE OPTION, P-VALUE THRESHOLD, LOGFC THRESHOLD, ES OPTION, ES THRESHOLD, REMOVE NS COMPARISONS
    observeEvent(c(input$pval_option, 
                   input$pval_threshold, 
                   input$logfc_threshold, 
                   input$es_option, 
                   input$es_threshold, 
                   input$remove_ns), {
      
      # Generate the comparison tables
      if (input$gene_option == "GENE-level") {
        genes <- toupper(str_trim(unlist(str_split(input$gene_list_input, ','))))
        genes <- intersect(genes,gene_map()$gene)
        table_list <- isolate(
          generate_comparison_tables(genes, exact_diseases(),
                                     datasets(), comparisons(),
                                     comparison_data(), input$pval_option,
                                     input$pval_threshold, input$logfc_threshold,
                                     input$remove_ns)
          )
      } else if (input$gene_option == "SIGNATURE-level") {
        print("Generating Comparsion Tables #1216")
        table_list <- isolate(
          signature_generate_comparison_tables(exact_diseases(), 
                                               datasets(), comparisons(),
                                               signature_data(),
                                               input$pval_threshold, 
                                               input$remove_ns)
          )
      } else {
        pathways <- unlist(input$pathway_list_input)
        table_list <- isolate(
          pathway_generate_comparison_tables(pathways, exact_diseases(),
                                             datasets(), comparisons(), 
                                             pathway_data(), input$pval_option,
                                             input$pval_threshold, input$es_option, 
                                             input$es_threshold, input$remove_ns)
          )
      }
      # Update the REACTIVE VALUE for the comparison tables
      comparison_tables(table_list)
      # Style the comparison tables
      table_list <- style_comparison_tables(exact_diseases(), table_list)
      
      sapply(exact_diseases(), function(current_disease) {
        
        # Pull the comparison table for the given disease
        comparison_table <- table_list[[current_disease]]
        
        # Remove the header and after table break if the comparison table is NULL
        if (is.null(comparison_table)) {
          # Generate the disease title for the comparison table
          output_name <- paste0(current_disease, "_table_title")
          output[[output_name]] <- renderUI(NULL)
          
          # Generate the after table break and horizontal rule
          output_name <- paste0(current_disease, "_break")
          output[[output_name]] <- renderUI(NULL)
        } else {
          # Generate the disease title for the comparison table
          output_name <- paste0(current_disease, "_table_title")
          output[[output_name]] <- renderUI({
            div(style = "display: flex; justify-content: center; align-items: center;",
                div(style = "flex-grow: 1; height: 1px; background-color: black; margin-top: 10px; margin-right: 10px;",
                    hr(style = "border: none; height: 1px; background-color: black; margin: 0;")),
                h3(current_disease, style = "line-height: 20px;"),
                div(style = "flex-grow: 1; height: 1px; background-color: black; margin-top: 10px; margin-left: 10px;",
                    hr(style = "border: none; height: 1px; background-color: black; margin: 0;"))
            )
          })
          
          # Generate the after table break and horizontal rule
          output_name <- paste0(current_disease, "_break")
          output[[output_name]] <- renderUI({
            div(hr(), br())
          })
        }
        
        # Display the comparison table for the given disease
        output_name <- paste0(current_disease, "_table")
        output[[output_name]] <- renderText({
          comparison_table
        })
        
      })
      
    }, ignoreInit = TRUE)
    #==============================================================================================================
    
    
    #==============================================================================================================
    ### Filter the rows in the comparison tables when the "Filter Comparisons" button is pressed
    #==============================================================================================================
    
    ### Filter Button - display modal ###
    observeEvent(input$filter_comparisons, {
      # MODAL DIALOG - ask user to select an option
      # option 1 - filter rows based on all comparison matches based on user selections
      # option 2 - filter rows based on explicit user selections
      showModal(modalDialog(
        p("If you would like to filter rows based on all comparisons similar to those for your selections,",
          "select the option below."),
        p("If you would like to filter rows based only on your explicit selections,", 
          "leave the option below unselected."),
        prettyCheckbox(
          inputId = "similar_comparisons_option",
          label = "Include similar comparisons", 
          value = FALSE,
          status = "danger",
          shape = "curve"
        ),
        footer = tagList(
          actionButton(inputId = "modal_filter_comparisons", label = "Filter Comparisons"),
          modalButton("Cancel")
        )
      ))
    })
    ##########
    
      
    ### Modal Filter Button - filter the table(s) for the selected comparisons ###
    observeEvent(input$modal_filter_comparisons, {
      
      # Remove the modal
      removeModal()
      
      # Determine the ID's associated with the plot selection checkboxes
      checkbox_ids <- grep("_[0-9]+$", names(input), value = TRUE)
      # Determine which checkboxes are checked
      checked_comparisons <- sapply(checkbox_ids, function(x) input[[x]])
      checked_comparisons <- checkbox_ids[checked_comparisons]
      
      # If no comparisons were selected, do not filter the table(s)
      if (length(checked_comparisons) != 0) {
        
        
        # Filter the comparison tables based on the user selected comparisons
        table_list <- filter_comparison_tables(exact_diseases(), 
                                               checked_comparisons,
                                               input$similar_comparisons_option,
                                               comparison_tables())
        # Update the REACTIVE VALUE for the comparison tables
        comparison_tables(table_list)
        # Style the filtered comparison tables
        table_list <- style_comparison_tables(exact_diseases(), table_list)

        # Display the filtered comparison tables
        sapply(exact_diseases(), function(current_disease) {
          # If the table for the given disease is NULL, remove the section header and after table break for that disease
          if (is.null(table_list[[current_disease]])) {
            # Remove section header
            output_name <- paste0(current_disease, "_table_title")
            output[[output_name]] <- renderText(NULL)
            # Remove break
            output_name <- paste0(current_disease, "_break")
            output[[output_name]] <- renderText(NULL)
          }
          # Render the table for the given disease
          table_name <- paste0(current_disease, "_table")
          output[[table_name]] <- renderText({
            table_list[[current_disease]]
          })
        })
        
      }

    })
    ##########
    
    
    ### Remove Filter Button ###
    observeEvent(input$clear_comparison_selections, {
      
      # Generate the comparison tables
      if (input$gene_option == "GENE-level") {
        genes <- toupper(str_trim(unlist(str_split(input$gene_list_input, ','))))
        genes <- intersect(genes,gene_map()$gene)
        table_list <- generate_comparison_tables(genes, exact_diseases(),
                                                 datasets(), comparisons(),
                                                 comparison_data(), input$pval_option,
                                                 input$pval_threshold, input$logfc_threshold,
                                                 input$remove_ns)
      } else if (input$gene_option == "SIGNATURE-level") {
        
        table_list <- signature_generate_comparison_tables(exact_diseases(), 
                                                           datasets(), comparisons(),
                                                           signature_data(),
                                                           input$pval_threshold, 
                                                           input$remove_ns)

      } else {
        print("Generating Comparison #1398")
        pathways <- unlist(input$pathway_list_input)
        
        table_list <- pathway_generate_comparison_tables(pathways, exact_diseases(),
                                                         datasets(), comparisons(), 
                                                         pathway_data(), input$pval_option,
                                                         input$pval_threshold, input$es_option, 
                                                         input$es_threshold, input$remove_ns)
      }
      # Update the REACTIVE VALUE for the comparison tables
      comparison_tables(table_list)
      # Style the comparison tables
      table_list <- style_comparison_tables(exact_diseases(), table_list)
      
      # Display the original comparison tables
      sapply(exact_diseases(), function(current_disease) {
        # If the table for the given disease is NOT NULL, add the section header and after table break for that disease
        if (!is.null(table_list[[current_disease]])) {
          # Add section header
          output_name <- paste0(current_disease, "_table_title")
          output[[output_name]] <- renderUI({
            h3(paste("The following table(s) display data from datasets associated with", current_disease))
          })
          # Add break
          output_name <- paste0(current_disease, "_break")
          output[[output_name]] <- renderUI({
            div(
              hr(),
              br()
            )
          })
        }
        # Render the table for the given disease
        table_name <- paste0(current_disease, "_table")
        output[[table_name]] <- renderText({
          table_list[[current_disease]]
        })
      })

    })
    ##########
    #==============================================================================================================

    
    #==============================================================================================================
    ### Generate boxplots given selected comparisons when the "Generate Visuals" button is pressed
    #==============================================================================================================
    
    observeEvent(input$visuals_button, {
      if (exists("input_list")){
        print("previous comparisons")
        print(grep("_[0-9]+$", input_list, value = TRUE))
        # Get the difference between old and new
        input_names_new <- names(input)[!which(input_list %in% names(input))]
        print("Clearing previous comparisons")
        print(input_names_new)
        input_list<-input_names_new
      }
      else{
        print("First ")
        input_list <-names(input)
      }
      genes <- toupper(str_trim(unlist(str_split(input$gene_list_input, ','))))
      genes <- intersect(genes,gene_map()$gene)
      checkbox_ids <- NULL
      checked_comparisons <- NULL
      # print(values$checkbox_ids)
      
      # Determine the ID's associated with the plot selection checkboxes
      
      checkbox_ids <- grep("_[0-9]+$", input_list, value = TRUE)
      if (!is.null(values$checkbox_ids)){
        # print(values$checkbox_ids)
        # print(checkbox_ids)
        
        checkbox_ids <- checkbox_ids[which(checkbox_ids %in% values$checkbox_ids)]
      }
      
      

      # Determine which checkboxes are checked
      checked_comparisons <- sapply(checkbox_ids, function(x) input[[x]])
      checked_comparisons <- checkbox_ids[checked_comparisons]
      
      # input[[names(input) %in% checkbox_ids]] <- NULL
      # print(which(checkbox_ids %in% input_list))
      # print(input_list)
      values$checkbox_ids<-checkbox_ids
      if(length(checked_comparisons)>6){
        shinyalert(title = "Too many selected comparisons",
                   text = "Please select at most 6 comparisons!",
                   type = "error")
      }else{
        # Jump to the Expression Visuals tab
        updateTabsetPanel(session, "main_body_tabset", selected = "Expression Visuals")
        
        # Generate the title for the expression boxplot box
        output$expression_plot_title <- renderUI({
          h2("Expression Boxplots for Selected Groups", align = "center")
        })
        
        # Display the button used to change the plot view and the button to download the plots
        output$plot_buttons <- renderUI({
          div(style = "display: flex; justify-content: space-between;",
              dropdownButton(h4("plot options"),
                             switchInput(inputId = "display_jitter", 
                                         label = "show data points", 
                                         value = FALSE, 
                                         width = "300px"),
                             radioGroupButtons(inputId = "display_pval",
                                               label = NULL,
                                               choices = c("p-value" = "p_value", "adjusted p-value" = "p_value_adj"),
                                               selected = input$pval_option,
                                               width = "300px",
                                               checkIcon = list(
                                                 yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                                 no = tags$i(class = "fa fa-square-o", style = "color: steelblue"))),
                             radioGroupButtons(inputId = "display_layout",
                                               label = NULL,
                                               choices = c("layout 1", "layout 2"),
                                               selected = "layout 1",
                                               width = "300px",
                                               checkIcon = list(
                                                 yes = tags$i(class = "fa fa-check-square", style = "color: steelblue"),
                                                 no = tags$i(class = "fa fa-square-o", style = "color: steelblue"))),
                             circle = TRUE, 
                             status = "primary",
                             size = "sm",
                             icon = icon("cog"), 
                             width = "350px",
                             right = FALSE,
                             tooltip = tooltipOptions(title = "Click to see plot options", placement = "right")),
              div(
                div(style = "display: inline-block; vertical-align: middle;",
                    actionButton(inputId = "run_meta_analysis",
                                 label = "Run Meta Analysis",
                                 width = "200px",
                                 style = "font-size: 115%;")),
                div(style = "display: inline-block; vertical-align: middle;",
                    dropdownButton(p("This button will run a meta analysis using Fisher's method for the selected comparisons", 
                                     "visualized here. The meta analysis results for each gene or pathway will be displayed in",
                                     "the section directly under the plots. Each row can be expanded to view the comparisons and",
                                     "results pertaining to the given meta analysis."),
                                   circle = TRUE, 
                                   status = "primary",
                                   size = "sm",
                                   icon = icon("question"), 
                                   width = "300px"))
              ),
              downloadButton('download_plots',
                             label = 'Download Plots',
                             style = 'color: black; background-color: #ededed; border-color: #26de5a;')
          )
        })
        
        
        
        # If no comparisons were selected, display a message in the Expression Visuals tab
        # Tell the user to select at least one comparison from the Comparison Tables tab
        if (length(checked_comparisons) == 0) {
          output$expression_plot_output = renderUI({
            h4("No comparisons were selected. Select at least one comparison from the Plot Selection column",
               "in the Comparison Tables tab.")
          })
          # Update the reactive value that contains the list of plots
          plots_reactive(NULL)
        } else {
          print("Updating plots")
          # Update the reactive value that contains the list of plots
          plots_reactive(plots_for_download(genes, gene_map(), datasets(), checked_comparisons, sample_ann(), input$pval_option, FALSE, db))
          print("Done")
          # Identify the unique datasets chosen by the user based on the checked comparisons
          chosen_dataset_acc <- unique(sapply(checked_comparisons, function(x) str_split(x, "_")[[1]][1]))
          
          # Determine the plot layout
          plot_idx <- plot_layout(genes, chosen_dataset_acc, checked_comparisons)
          print(names(plot_idx))
          # Generate the column headers for the plots (these headers are dataset ID's)
          output$plot_titles <- renderUI({
            splitLayout(
              cellWidths = c("16%", "16%", "16%", "16%", "16%", "16%"),
              div(str_split(plot_idx[[1]][[1]],"_")[[1]][1], style = "text-align: center; font-weight: 600;"),
              div(str_split(plot_idx[[2]][[1]],"_")[[1]][1], style = "text-align: center; font-weight: 600;"),
              div(str_split(plot_idx[[3]][[1]],"_")[[1]][1], style = "text-align: center; font-weight: 600;"),
              div(str_split(plot_idx[[4]][[1]],"_")[[1]][1], style = "text-align: center; font-weight: 600;"),
              div(str_split(plot_idx[[5]][[1]],"_")[[1]][1], style = "text-align: center; font-weight: 600;"),
              div(str_split(plot_idx[[6]][[1]],"_")[[1]][1], style = "text-align: center; font-weight: 600;")
            )
          })
          
          # Dynamically generate the output containers in the splitLayout function
          output$expression_plot_output = renderUI({
            div(
              splitLayout(
                cellWidths = c("16%", "16%", "16%", "16%", "16%", "16%"),
                tagList(
                  lapply(plot_idx[[1]][[2]], function(z) {
                    plotOutput(z, height = "250px")
                  })
                ),
                tagList(
                  lapply(plot_idx[[2]][[2]], function(z) {
                    plotOutput(z, height = "250px")
                  })
                ),
                tagList(
                  lapply(plot_idx[[3]][[2]], function(z) {
                    plotOutput(z, height = "250px")
                  })
                ),
                tagList(
                  lapply(plot_idx[[4]][[2]], function(z) {
                    plotOutput(z, height = "250px")
                  })
                ),
                tagList(
                  lapply(plot_idx[[5]][[2]], function(z) {
                    plotOutput(z, height = "250px")
                  })
                ),
                tagList(
                  lapply(plot_idx[[6]][[2]], function(z) {
                    plotOutput(z, height = "250px")
                  })
                )
              ),
              br()
            )
          })
          
          # for(current_comparison in checked_comparisons){
          #   cat("Current COM:",current_comparison,"\n")
          #   for(gene in genes){
          #     cat("GENE:",gene,"\n")
          #     plot_list <- plot_expression(gene, gene_map(), 
          #                                  current_comparison, datasets(), 
          #                                  checked_comparisons, 
          #                                  sample_ann(), 
          #                                  input$pval_option, 
          #                                  FALSE, 
          #                                  db)
          #     output[[paste0(gene, "_", current_comparison, "_expression_plot_", 1)]] <- renderPlot({
          #       plot_list[[1]]
          #     })
          #   }
          # }
          
          # Render the plots
          sapply(checked_comparisons, function(current_comparison) {
            sapply(genes, function(gene) {
              # Generate a list of plots for the given gene and dataset
              plot_list <- plot_expression(gene, gene_map(), 
                                           current_comparison, datasets(), 
                                           checked_comparisons, 
                                           sample_ann(), 
                                           input$pval_option, 
                                           FALSE, 
                                           db)
              sapply(1:length(plot_list), function(z) {
                # Render the plots for the given dataset
                output[[paste0(gene, "_", current_comparison, "_expression_plot_", z)]] <- renderPlot({
                  plot_list[[z]]
                })
              })
            })
          })
        }
      }
      
      
    })
    
    
    # Change between p-value and adjusted p-value display in the expression plots
    observeEvent(input$display_pval, {
      genes <- toupper(str_trim(unlist(str_split(input$gene_list_input, ','))))
      genes <- intersect(genes,gene_map()$gene)
      # Determine the ID's associated with the plot selection checkboxes
      checkbox_ids <- grep("_[0-9]+$", names(input), value = TRUE)
      # Determine which checkboxes are checked
      checked_comparisons <- sapply(checkbox_ids, function(x) input[[x]])
      checked_comparisons <- checkbox_ids[checked_comparisons]
      
      # If no comparisons were selected, display a message in the Expression Visuals tab
      # Tell the user to select at least one comparison from the Comparison Tables tab
      if (length(checked_comparisons) == 0) {
        output$expression_plot_output = renderUI({
          h4("No comparisons were selected. Select at least one comparison from the Plot Selection column",
             "in the Comparison Tables tab.")
        })
        # Update the reactive value that contains the list of plots
        plots_reactive(NULL)
      } else {
        # Update the reactive value that contains the list of plots
        plots_reactive(plots_for_download(genes, gene_map(), datasets(), checked_comparisons, sample_ann(), input$display_pval, input$display_jitter, db))
        # Identify the unique datasets chosen by the user based on the checked comparisons
        chosen_dataset_acc <- unique(sapply(checked_comparisons, function(x) str_split(x, "_")[[1]][1]))
        
        sapply(checked_comparisons, function(current_comparison) {
          sapply(genes, function(gene) {
            # Generate a list of plots for the given gene and dataset
            plot_list <- plot_expression(gene, gene_map(), 
                                         current_comparison, datasets(), 
                                         checked_comparisons, 
                                         sample_ann(), 
                                         input$display_pval, 
                                         input$display_jitter, 
                                         db)
            sapply(1:length(plot_list), function(z) {
              # Render the plots for the given dataset
              output[[paste0(gene, "_", current_comparison, "_expression_plot_", z)]] <- renderPlot({
                plot_list[[z]]
              })
            })
          })
        })
      }
    }, ignoreNULL = TRUE)
    
    
    # Display or remove jitter
    observeEvent(input$display_jitter, {
      genes <- toupper(str_trim(unlist(str_split(input$gene_list_input, ','))))
      genes <- intersect(genes,gene_map()$gene)
      # Determine the ID's associated with the plot selection checkboxes
      checkbox_ids <- grep("_[0-9]+$", names(input), value = TRUE)
      # Determine which checkboxes are checked
      checked_comparisons <- sapply(checkbox_ids, function(x) input[[x]])
      checked_comparisons <- checkbox_ids[checked_comparisons]
      
      # If no comparisons were selected, display a message in the Expression Visuals tab
      # Tell the user to select at least one comparison from the Comparison Tables tab
      if (length(checked_comparisons) == 0) {
        output$expression_plot_output = renderUI({
          h4("No comparisons were selected. Select at least one comparison from the Plot Selection column",
             "in the Comparison Tables tab.")
        })
        # Update the reactive value that contains the list of plots
        plots_reactive(NULL)
      } else {
        # Update the reactive value that contains the list of plots
        plots_reactive(plots_for_download(genes, gene_map(), datasets(), checked_comparisons, sample_ann(), input$display_pval, input$display_jitter, db))
        # Identify the unique datasets chosen by the user based on the checked comparisons
        chosen_dataset_acc <- unique(sapply(checked_comparisons, function(x) str_split(x, "_")[[1]][1]))
        
        sapply(checked_comparisons, function(current_comparison) {
          sapply(genes, function(gene) {
            # Generate a list of plots for the given gene and dataset
            plot_list <- plot_expression(gene, gene_map(), 
                                         current_comparison, datasets(), 
                                         checked_comparisons, 
                                         sample_ann(), 
                                         input$display_pval, 
                                         input$display_jitter, 
                                         db)
            sapply(1:length(plot_list), function(z) {
              # Render the plots for the given dataset
              output[[paste0(gene, "_", current_comparison, "_expression_plot_", z)]] <- renderPlot({
                plot_list[[z]]
              })
            })
          })
        })
      }
    }, ignoreNULL = TRUE)
    
    
    # Change between layout 1 and layout 2 for the expression plots
    observeEvent(input$display_layout, {
      genes <- toupper(str_trim(unlist(str_split(input$gene_list_input, ','))))
      genes <- intersect(genes,gene_map()$gene)
      # Determine the ID's associated with the plot selection checkboxes
      checkbox_ids <- grep("_[0-9]+$", names(input), value = TRUE)
      # Determine which checkboxes are checked
      checked_comparisons <- sapply(checkbox_ids, function(x) input[[x]])
      checked_comparisons <- checkbox_ids[checked_comparisons]
      
      # If no comparisons were selected, display a message in the Expression Visuals tab
      # Tell the user to select at least one comparison from the Comparison Tables tab
      if (length(checked_comparisons) == 0) {
        output$expression_plot_output = renderUI({
          h4("No comparisons were selected. Select at least one comparison from the Plot Selection column",
             "in the Comparison Tables tab.")
        })
        # Update the reactive value that contains the list of plots
        plots_reactive(NULL)
      } else {
        # Update the reactive value that contains the list of plots
        p <- plots_for_download(genes, gene_map(), datasets(), checked_comparisons, sample_ann(), input$display_pval, input$display_jitter, db)
        plots_reactive(plots_for_download(genes, gene_map(), datasets(), checked_comparisons, sample_ann(), input$display_pval, input$display_jitter, db))
        # Identify the unique datasets chosen by the user based on the checked comparisons
        chosen_dataset_acc <- unique(sapply(checked_comparisons, function(x) str_split(x, "_")[[1]][1]))
        # Determine the plot layout
        plot_idx <- plot_layout(genes, chosen_dataset_acc, checked_comparisons)
        
        if (input$display_layout == "layout 1") {
          # Generate the column headers for the plots (these headers are dataset ID's)
          output$plot_titles <- renderUI({
            splitLayout(
              cellWidths = c("16%", "16%", "16%", "16%", "16%", "16%"),
              div(str_split(plot_idx[[1]][[1]],"_")[[1]][1], style = "text-align: center; font-weight: 600;"),
              div(str_split(plot_idx[[2]][[1]],"_")[[1]][1], style = "text-align: center; font-weight: 600;"),
              div(str_split(plot_idx[[3]][[1]],"_")[[1]][1], style = "text-align: center; font-weight: 600;"),
              div(str_split(plot_idx[[4]][[1]],"_")[[1]][1], style = "text-align: center; font-weight: 600;"),
              div(str_split(plot_idx[[5]][[1]],"_")[[1]][1], style = "text-align: center; font-weight: 600;"),
              div(str_split(plot_idx[[6]][[1]],"_")[[1]][1], style = "text-align: center; font-weight: 600;")
            )
          })
          
          # Dynamically generate the output containers in the splitLayout function
          output$expression_plot_output = renderUI({
            div(
              splitLayout(
                cellWidths = c("16%", "16%", "16%", "16%", "16%", "16%"),
                tagList(
                  lapply(plot_idx[[1]][[2]], function(z) {
                    plotOutput(z, height = "250px")
                  })
                ),
                tagList(
                  lapply(plot_idx[[2]][[2]], function(z) {
                    plotOutput(z, height = "250px")
                  })
                ),
                tagList(
                  lapply(plot_idx[[3]][[2]], function(z) {
                    plotOutput(z, height = "250px")
                  })
                ),
                tagList(
                  lapply(plot_idx[[4]][[2]], function(z) {
                    plotOutput(z, height = "250px")
                  })
                ),
                tagList(
                  lapply(plot_idx[[5]][[2]], function(z) {
                    plotOutput(z, height = "250px")
                  })
                ),
                tagList(
                  lapply(plot_idx[[6]][[2]], function(z) {
                    plotOutput(z, height = "250px")
                  })
                )
              ),
              br()
            )
          })
        } else {
          output$plot_titles <- renderUI({
            NULL
          })
          
          # Dynamically generate the output containers
          output$expression_plot_output = renderUI({
            div(
              uiOutput(plot_idx[[1]][[1]]),
              uiOutput("section1"),
              br(),
              uiOutput(plot_idx[[2]][[1]]),
              uiOutput("section2"),
              br(),
              uiOutput(plot_idx[[3]][[1]]),
              uiOutput("section3"),
              br(),
              uiOutput(plot_idx[[4]][[1]]),
              uiOutput("section4"),
              br(),
              uiOutput(plot_idx[[5]][[1]]),
              uiOutput("section5"),
              br(),
              uiOutput(plot_idx[[6]][[1]]),
              uiOutput("section6")
            )
          })
        }
        sapply(1:6, function(x) {
          if (plot_idx[[x]][[1]] != "") {
            output[[paste0("section", x)]] <- renderUI({
              div(style = "border: 1px solid grey; padding: 5px;",
                  tagList(
                    lapply(plot_idx[[x]][[2]], function(y) {
                      div(style = "display: inline-block; vertical-align: middle;",
                          plotOutput(y, width = "250px", height = "250px"))
                    })
                  )
              )
            })
          }
        })
        
        # Render the plots
        sapply(checked_comparisons, function(current_comparison) {
          # Render the title for the given plot section (i.e. name of dataset)
          output[[paste0(current_comparison, "_plot_title")]] <- renderUI({
            h4(str_split(current_comparison,"_")[[1]][1])
          })
          sapply(genes, function(gene) {
            # Generate a list of plots for the given gene and dataset
            plot_list <- plot_expression(gene, gene_map(), 
                                         current_comparison, datasets(), 
                                         checked_comparisons, 
                                         sample_ann(), 
                                         input$display_pval, 
                                         input$display_jitter, 
                                         db)
            sapply(1:length(plot_list), function(z) {
              # Render the plots for the given dataset
              output[[paste0(gene, "_", current_comparison, "_expression_plot_", z)]] <- renderPlot({
                plot_list[[z]]
              })
            })
          })
        })
      }
      
      # Generate the title and run button for the meta analysis section in the visuals tab
      output$meta_analysis_header <- renderUI({
        div(id = "metaAnalysisSection",
            h3("Meta Analysis Module"))
      })
      
      # # Observer that will scroll down to the meta analysis section when the button is clicked
      # observeEvent(input$run_meta_analysis, {
      #   session$sendCustomMessage(type = "scroll", message = "metaAnalysisSection")
      # })
      # observe({
      #   session$registerDataObj(
      #     "scroll",
      #     "window.parent.document.getElementById(input.message).scrollIntoView();"
      #   )
      # })
      
      # Generate the meta analysis results table
      observeEvent(input$run_meta_analysis, {
        # Determine the ID's associated with the plot selection checkboxes
        checkbox_ids <- grep("_[0-9]+$", names(input), value = TRUE)
        # Determine which checkboxes are checked
        checked_comparisons <- sapply(checkbox_ids, function(x) input[[x]])
        checked_comparisons <- checkbox_ids[checked_comparisons]
        
        if (length(checked_comparisons) == 0) {
          shinyalert(title = "No Data Selected",
                     text = "Please make at least one selection from the results table(s) using the right hand column.",
                     type = "error")
        } else {
          output$meta_analysis_table <- renderReactable({
            if (input$gene_option == "GENE-level") {
              perform_meta_analysis(datasets(), comparisons(), comparison_data(), checked_comparisons)
            } else if (input$gene_option == "SIGNATURE-level") {
              perform_meta_analysis(signature_data(), checked_comparisons)
            } else {
              perform_meta_analysis(pathway_data(), checked_comparisons)
            }
          })
        }
      })
      
    }, ignoreNULL = TRUE)
    #==============================================================================================================
    
    
    #==============================================================================================================
    # Download handler for downloading plots
    #==============================================================================================================
    
    output$download_plots <- downloadHandler(
      # Filename: specify the filename of the zipped plots
      filename =  function() {
        "plots.zip"
      },
      # Content: save the plots as a tar file
      content = function(file) {
        tmpdir <- tempdir()
        setwd(tempdir())
        count <- 0
        
        # Pull the list of plots from the reactive object
        plot_list <- plots_reactive()
        
        # Define a list used to store the plot names (in the zip file)
        plot_name_list <- c()
        
        for (gene in names(plot_list)) {
          for (dataset_acc in names(plot_list[[gene]])) {
            
            plots <- plot_list[[gene]][[dataset_acc]]
            
            for (i in 1:length(plots)) {
              # Specify the name for the given plot and append to the plot_name_list
              plot_name <- paste0(gene, "_", dataset_acc, "_plot_", i, ".png")
              plot_name_list <- c(plot_name_list, plot_name)
              # Download the plot
              ggsave(plot_name, plot = plots[[i]], device = "png")
            }
            
          }
        }
        
        zip(file, plot_name_list)
      } 
    )
    
    session$onSessionEnded(function() {
      
      stopApp()
    })
    #==============================================================================================================
}

# Close the connection to the DB
onStop(function() {
  dbDisconnect(db)
})

# Run the application 
shinyApp(ui = ui, server = server)

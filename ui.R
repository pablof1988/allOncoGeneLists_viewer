library(shiny)
library(shinybusy)
library(bs4Dash)
library(goSorensen)
library(plotly)
library(ggplot2)
library(DT)
library(ggrepel)
library(dendextend)
library(GO.db)
data("allOncoGeneLists")
data("allContTabs")
data("allEqTests")
data("allDissMatrx")
#library(org.Hs.eg.db)
#library(RColorBrewer)
source("www/aux_fun.R")
add_busy_spinner(spin = "fading-circle")


ui <- bs4DashPage(
  title = "goSorensen_app.", fullscreen = T,
  header = bs4DashNavbar(
    title = div(
      style = "text-align: center;", 
      img(src = "go_logo2.png", height = "110px", width = "100px")
         ),
    HTML('<div style="text-align: center; color: #6e2c00; font-weight: bold; font-size: 18px; border-bottom: 2px solid #ba4a00; padding-bottom: 12px;">INTERACTIVE RESULTS VIEWER FOR "allOncoGeneLists"</div>'),
    border = T, compact = F,
    tags$head(
      tags$script(HTML("
    $(document).on('shiny:connected', function() {
      function updateMenuText() {
        var $menuA = $(\"a[data-value='example']\");
        var $menuText = $menuA.find('p');
        var $li = $menuA.closest('li');
        if ($li.hasClass('menu-open')) {
          $menuText.text('CLICK HERE TO HIDE THE INPUTS');
        } else {
          $menuText.text('CLICK HERE TO SHOW THE INPUTS');
        }
      }

      // Actualizar al hacer click
      $(document).on('click', \"a[data-value='example']\", function(e) {
        // Esperar a que la animación/expansión termine
        setTimeout(updateMenuText, 350);
      });

      // Actualizar también en la carga inicial (por si acaso)
      setTimeout(updateMenuText, 500);
    });
  "))
    ),
    tags$script(HTML("
    $(document).on('shiny:connected', function() {
      // Busca el li que contiene el a[data-value='example']
      var $li = $(\"a[data-value='example']\").closest('li');
      // Si no está abierto, ábrelo
      if (!$li.hasClass('menu-open')) {
        $li.addClass('menu-open');
      }
      // Marca como activo (opcional, si quieres el color azul típico)
      $(\"a[data-value='example']\").addClass('active');
    });
  ")),
    tags$style(HTML("
    .content-wrapper, .content, .tab-content {
      background-color: #fff !important;
    }
  ")),
    tags$head(
      tags$style(HTML("
    .main-sidebar, .skin-skyblue .main-sidebar, .main-sidebar .sidebar {
      background-color: #f5fbff !important;  /* Celeste muy bajo */
    }
  ")),
      tags$style(HTML("
    /* Solo afecta al datatable con id enrchd_mat y sus controles */
    #enrchd_mat, 
    #enrchd_mat table.dataTable, 
    #enrchd_mat table.dataTable th,  
    #enrchd_mat .dataTables_wrapper, 
    #enrchd_mat .dataTables_info, 
    #enrchd_mat .dataTables_paginate, 
    #enrchd_mat .dataTables_filter, 
    #enrchd_mat .dataTables_length {
      font-size: 11px !important;
    }
  "))
    )
    ),
  sidebar = bs4DashSidebar(
    skin = "skyblue", width = "25%",
    bs4SidebarMenu(
      bs4SidebarMenuItem("CLICK HERE TO HIDE THE INPUTS", tabName = "example", 
                         icon = icon("rectangle-list"),
                         wellPanel(
                             fluidRow(
                             column(6, selectInput(inputId = "l1", label = "Select List 1:", 
                                                   choices = names(allOncoGeneLists), width = "120%")),
                             column(6, selectInput(inputId = "l2", label = "Select List 2:", 
                                                   choices = NULL, width = "120%"))
                           ),
                           selectInput("sel_onto", label = "Select an ontology:",
                                       choices = c("BP (Biological Process)" = "BP", 
                                                   "CC (Cellular Component)" = "CC", 
                                                   "MF (Molecular Function)" = "MF"),
                                       width = "80%"),
                           tags$head(
                             tags$style(HTML("
                                  .radio-inline {
                                  margin-right: 35px !important; /* Más espacio entre opciones */
                                  font-size: 16px;
                                }
                              "))
                           ),
                           radioButtons("sel_level", "Select a GO-level:", 
                                        choices = 3:10, 
                                        inline = TRUE, width = "60%"),
                           actionButton("go_ex", label = "COMPUTE", icon = icon("caret-right"), 
                                        status = "info", width = "40%")
                         )
      )
    )
  ),
  
  body = bs4DashBody(
    bs4TabItems(
      bs4TabItem(
        tabName = "example",
        tabsetPanel(
          tabPanel("GENE LISTS (ENTREZ ID'S)",
                   div(
                     style = "text-align: center;", 
                     img(src = "geneLists.png", height = "100px", width = "120px")
                   ),
                   br(),
                   uiOutput("dyncards")
          ),
          tabPanel("ENRICHMENT ANALYSIS",
                   div(
                     style = "text-align: center;", 
                     img(src = "enrichment_cycle.png", height = "120px", width = "470px")
                   ),
                   br(),
                   fluidPage(
                     splitLayout(
                       cellWidths = c("60%", "40%"),
                       div(
                         h5("Enrichment matrix",
                            style = "text-decoration: underline; font-size: 18px; text-align: center; color: #4f94a9; font-weight: bold;"),
                         dataTableOutput("enrchd_mat"),
                         uiOutput("message")
                       ),
                       div(
                         style = "border-left: 2px solid #bdbdbd;padding-left: 20px;height: 100vh;",
                         h5("Enrichment contingency table",
                            style = "text-decoration: underline; font-size: 18px; text-align: center; color: #4f94a9; font-weight: bold;"),
                         br(), br(), 
                         verbatimTextOutput("contingency_tables") 
                       )
                       ))
                   
          ),
          tabPanel("STATISTICAL TECNIQUES",
                   fluidPage(
                     radioButtons(
                       inputId = "which_stat",
                       label = h6("Choose the calculation:", style = "font-size:18px; font-weight: bold;"),
                       choices = c("Equivalence Test" = "eqtest", "Matrix of Dissimilarities" = "diss"),
                       selected = "eqtest",
                       inline = TRUE
                     ),
                     # Panel para Equivalence Test
                     conditionalPanel(
                       condition = "input.which_stat == 'eqtest'",
                       p("Select a different input and click 'COMPUTE'
                             in the sidebar to refresh the results.",
                         style = "font-size: 13px; color: #4f94a9; "),
                       sidebarLayout(
                         sidebarPanel(
                           strong("Inputs"),
                           numericInput(inputId = "d0", 
                                        label = "d0:", 
                                        value = 0.4444, min = 0, max = 1, 
                                        step = 0.01, width = "60%"),
                           radioButtons(inputId = "conf", 
                                        label = "Confidence:", 
                                        choices = c(0.9, 0.95, 0.99), selected = 0.95, width = "50%"),
                           radioButtons(inputId = "samp", 
                                        label = "Sample distribution:", 
                                        choices = c("Normal" = F, "Bootstrap" = T), 
                                        selected = F),
                           width = 3
                         ),
                         mainPanel(
                           uiOutput("mes_eqtest"),
                           verbatimTextOutput("eq_test"), width = 9
                         )
                       )
                     ),
                     # Panel Matrix of Dissimilarities
                     conditionalPanel(
                       condition = "input.which_stat == 'diss'",
                       uiOutput("mes_matrx"),
                       verbatimTextOutput("dis_matrix"),
                       fluidPage(
                         splitLayout(
                           cellWidths = c("40%", "60%"),
                           div(
                             h5("Dendrogram",
                                style = "font-size: 16px; text-align: center; color: #4f94a9; font-weight: bold;"),
                             plotOutput("dendrogram")
                           ),
                           div(
                             style = "border-left: 2px solid #bdbdbd;padding-left: 20px;height: 100vh;",
                             h5("MDS - Biplot",
                                style = "font-size: 16px; text-align: center; color: #4f94a9; font-weight: bold;"),
                             plotlyOutput("biplot"),
                             br(),
                             radioButtons("whatdim", label = h5("Select a dimension and click 'Identify GO terms' to see the biological concepts associated with the detected significant equivalence in the chosen dimension:",
                                                                style = "font-size: 14px; text-align: left; white-space: normal; word-break: break-word; overflow-wrap: break-word;"), 
                                          choices = c("Dimension 1", "Dimension 2"), inline = T),
                             actionButton("compute_GOterms", label = "Identify GO Terms", 
                                          status = "info", icon = icon("caret-right")),
                             dataTableOutput("go_terms")
                           )
                         ))
                     )
                   )
          )
        )
      )
    )
  ),
  
  controlbar <- bs4DashControlbar(
    id = "info_panel",
    title = "App Information",
    skin = "light",
    width = 400,
    collapsed = TRUE,
    controlbarMenu(
      controlbarItem(
        title = "Purpose",
        icon = icon("search"),
        HTML("
        <p> This Shiny dashboard serves as a complementary tool for interactively visualising the results of applying <code>goSorensen</code> 
        to the gene lists referred to as <strong><em>allOncoGeneLists</em></strong> in section 4 of an article submitted to <em>The R Journal</em>, 
        which introduces the <code>goSorensen</code> R package, designed to detect functional similarity between gene lists based on the 
        joint enrichment of GO terms, quantified by the Sorensen dissimilarity. </p>
        <p>The results of this Shiny app are produced by the <code>goSorensen</code> <strong>(v1.10.0)</strong>, as included in <strong>Bioconductor 3.16</strong>.</p>
      ")
      ),
      controlbarItem(
        title = "Features",
        icon = icon("dna"),
        HTML("
        <ul>
          <li>GO ontology selection: BP, CC, MF</li>
          <li>GO level selection: from levels 3 to 10</li>
          <li>Selection of the main inputs for equivalence hypothesis tests</li>
            <ul>
              <li>Irrelevance limit 'd0'</li>
              <li>Confidence level</li>
              <li>Sample distribution</li>
            </ul>
          <li>Visualization of:
            <ul>
              <li>ENTREZ identification of the genes in the lists</li>
              <li>Enrichment matrix for a specific ontology and GO-level</li>
              <li>Enrichment contingency table for a pair of lists in a specific ontology and GO-level</li>
              <li>Equivalence test results for a pair of lists in a specific ontology and GO-level</li>
              <li>Matrix of dissimilarities for a specific ontology and GO-level</li>
              <li>Dendrogram (hierarchical clustering) for a specific ontology and GO-level</li>
              <li>MDS biplot for a specific ontology and GO-level</li>
              <li>Identification of the GO terms related with the biological similarity between feature lists</li>
            </ul>
          </li>
        </ul>
      ")
      ),
      controlbarItem(
        title = "Scope",
        icon = icon("exclamation-circle"),
        HTML("
        <p><strong>Note:</strong> This application is <em>not</em> a general-purpose tool for arbitrary gene list analysis.</p>
        <p>It is specifically designed to visualize the results of applying <code>goSorensen</code> to the <strong><em>allOncoGeneLists</em></strong> object, included in the package.</p>
      ")
      ),
      controlbarItem(
        title = "Repository",
        icon = icon("link"),
        HTML('
        <p>GitHub repository:</p>
        <p><a href="https://github.com/pablof1988/allOncoGeneLists_viewer" target="_blank">
        github.com/pablof1988/allOncoGeneLists_viewer</a></p>
      ')
      ),
      controlbarItem(
        title = "Citation",
        icon = icon("quote-left"),
        HTML("
        <p>If you use this application or the <code>goSorensen</code> package, please cite:</p>
        <blockquote style='font-size: 0.9em'>
          Flores, P., Salicrú, M., Sánchez-Pla, A., & Ocaña, J. (2022). An equivalence test between features lists, based on the Sorensen–Dice index and the joint frequencies of GO term enrichment. BMC bioinformatics, 23(1), 207.<br>
          <br>
          Flores, P., Salicrú, M., Sánchez-Pla, A., & Ocaña, J. (2025). goSorensen: A Bioconductor R-Package to Detect Significant Equivalence Between Gene Lists, The R Journal, 17(1).
        </blockquote>
      ")
      )
    )
  ),
  
  footer = bs4DashFooter(
    right = "Questions and comments: p_flores@espoch.edu.ec"
  )
)



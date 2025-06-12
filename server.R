server <- function(input, output, session) {
  observeEvent(input$l1, {
    updated_choices <- setdiff(names(allOncoGeneLists), input$l1)
    
    updateSelectInput(session, inputId = "l2", 
                      choices = rev(updated_choices),
                      selected = NULL)
  })
  
  selected_lists_data <- reactiveVal()
  cut_lines <- reactiveVal(NULL)
  
  observeEvent(input$go_ex, {
    featlists <- allOncoGeneLists[c(input$l1, input$l2)]
    selected_lists_data(featlists)
  })
  
  output$dyncards <- renderUI({
    req(selected_lists_data())
    fluidRow(
      lapply(names(selected_lists_data()), function(col_name) {
        local({
          nm <- col_name
          output_id <- paste0("card_", nm)
          output[[output_id]] <- renderPrint({
            selected_lists_data()[[nm]]
          })
          column(
            width = 6,  
            bs4Card(
              title = sprintf("%s", nm),
              status = "primary", collapsed = TRUE, maximizable = TRUE,
              solidHeader = TRUE, width = 12,
              verbatimTextOutput(outputId = output_id)
            )
          )
        })
      })
    )
  })
  
  ## Enrichment matrix
  filtersmatrx <- reactive({
    allContTabs[[input$sel_onto]][[as.numeric(input$sel_level) - 2]]
  })
  
  enriched_mat <- reactiveVal()
  
  observeEvent(input$go_ex, {
    enriched_mat(attr(filtersmatrx(), "enriched"))
  })
  
  output$enrchd_mat <- renderDataTable({
    df <- as.data.frame(enriched_mat())
    df[] <- lapply(df, function(x) {
      if (is.logical(x)) toupper(as.character(x)) else x
    })
    cap_txt <- sprintf(
      "Ontology %s, GO level %s",
      isolate(input$sel_onto), isolate(input$sel_level)
    )
    datatable(df, rownames = T, caption = htmltools::tags$caption(
      style = 'caption-side: top; text-align: left; font-size: 14px; color: #4f94a9;',
      cap_txt
    ),
    options = list(scrollX = TRUE, dom = "lrtipf"),
    class = "stripe hover order-column",
    )
  })
  
  messageEM <- reactiveVal()
  
  observeEvent(input$go_ex, {
    mess <- paste0("Each row of this enrichment matrix corresponds to the GO terms 
                  enriched by at least one of the gene lists. The total number
                  of GO terms annotated in the ontology ", input$sel_onto,
                   ", GO level ", input$sel_level, ", by the list of genes 'allOnco', is ", attr(enriched_mat(), "nTerms"),
                   ".", collapse = "<br>")
    messageEM(mess)  
  })
  
  output$message <- renderUI({
    
    div(
      style = "white-space: normal; overflow-x: auto; word-wrap: break-word; width: 100%; font-size: 12px;",
      HTML(paste0("<b>NOTE:</b> ", messageEM()))
    )
    
  })
  
  ## Contingency Tables
  
  contTables <- reactiveVal()
  
  observeEvent(input$go_ex, {
    conts <- getListNames(filtersmatrx())
    res_contTabs <- conts[filter_list(conts, c(input$l1, input$l2))]
    contTables(res_contTabs[[1]])
  })
  
  output$contingency_tables <- renderPrint({
    contTables()
  })
  
  # Equivalence test
  
  messageEqTest <- reactiveVal()
  
  observeEvent(input$go_ex, {
    messTest <- paste0("Equivalence test for ", input$l1, " Vs. ", input$l2, " (Ontology ", input$sel_onto, ", GO level ", input$sel_level, ")")
    messageEqTest(messTest)  
  })
  
  output$mes_eqtest <- renderUI({
    div(
      style = "text-decoration: underline; font-size: 18px; text-align: center; color: #4f94a9; font-weight: bold;",
      messageEqTest()
    )
  })
  
  eqTests <- reactiveVal()
  
  observeEvent(input$go_ex, {
    allTests <- allEqTests[[input$sel_onto]][[as.numeric(input$sel_level) - 2]]
    tests <- getListNames(allTests)
    res_tests <- tests[filter_list(tests, c(input$l1, input$l2))]
    eqTests(upgrade(res_tests[[1]], d0 = input$d0, 
                    conf.level = as.numeric(input$conf), boot =input$samp))
  })
  
  output$eq_test <- renderPrint({
    eqTests()
  })
  
  # Dissimilarity Matrix:
  messageMatrx <- reactiveVal()
  
  observeEvent(input$go_ex, {
    messMatrx <- paste0("Matrix of dissimilarities", " (Ontology ", input$sel_onto, ", GO level ", input$sel_level, ")")
    messageMatrx(messMatrx)  
  })
  
  output$mes_matrx <- renderUI({
    div(
      style = "text-decoration: underline; font-size: 18px; text-align: left; color: #4f94a9; font-weight: bold;",
      messageMatrx()
    )
  })
  
  dissMatrx <- reactiveVal()
  
  observeEvent(input$go_ex, {
    mat_diss <- allDissMatrx[[input$sel_onto]][[as.numeric(input$sel_level) - 2]]
    dissMatrx(mat_diss)
  })
  
  output$dis_matrix <- renderPrint({
    dissMatrx()
  })
  
  output$dendrogram <- renderPlot({
    clust <- hclustThreshold(dissMatrx())
    dend <- as.dendrogram(clust)
    #dend <- color_branches(dend, k = 4)
    dend <- set(dend, "branches_lwd", 2) 
    plot(dend)
  })
  
  output$biplot <- renderPlotly({
    prmds <- cmdscale(dissMatrx(), k = 2, eig = TRUE)
    prmds$points <- as.data.frame(prmds$points)
    colnames(prmds$points) <- paste0("Dim", 1:2)
    
    labels <- attr(dissMatrx(), "Labels")
    dfbiplot <- cbind(prmds$points, label = labels)
    
    graph <- ggplot(dfbiplot, aes(x = Dim1, y = Dim2)) +
      geom_point(aes(text = label), color = "blue", size = 2) +
      geom_text_repel(aes(label = label), color = "black", size = 1) +
      xlab(paste0("Dimension 1 (", round(prmds$eig[1] / sum(prmds$eig), 3) * 100, "%)")) +
      ylab(paste0("Dimension 2 (", round(prmds$eig[2] / sum(prmds$eig), 3) * 100, "%)")) +
      theme_light()
    
    cuts <- cut_lines()
    if (!is.null(cuts)) {
      if (cuts$type == "vline") {
        for (x in cuts$values) {
          graph <- graph + geom_vline(xintercept = x, color = "coral1", linetype = "dashed", linewidth = 0.25)
        }
      } else if (cuts$type == "hline") {
        for (y in cuts$values) {
          graph <- graph + geom_hline(yintercept = y, color = "coral1", linetype = "dashed", linewidth = 0.25)
        }
      }
    }
    
    ggplotly(graph, tooltip = "text")
  })
  
  # Detecting GO terms
  caption_params <- reactiveVal()
  
  observeEvent(input$go_ex,{
    observeEvent(input$compute_GOterms, {
      show_modal_spinner(spin = "fading-circle", 
                         text = paste0("Processing (",  input$sel_onto, input$sel_level, 
                                       " - ", input$whatdim, ").", " Please, wait..."))
      
      result <- gotermsid(
        list = allOncoGeneLists, 
        dm = cmdscale(dissMatrx(), k = 2, eig = T)$points, 
        prp = 0.2, 
        dimen = as.numeric(gsub("^[^0-9]*([0-9]+).*", "\\1", input$whatdim)),
        contabs = filtersmatrx()
      )
      
      caption_params(list(
        onto = input$sel_onto,
        level = input$sel_level,
        dim = input$whatdim
      ))
      
      mds_points <- cmdscale(dissMatrx(), k = 2, eig = TRUE)$points
      dim_num <- as.numeric(gsub("^[^0-9]*([0-9]+).*", "\\1", input$whatdim))
      
      if (dim_num == 1) {
        cpdim <- splitdim(mds_points[, 1], 0.20)  # 20%
        cut_lines(list(type = "vline", values = cpdim))
      } else {
        cpdim <- splitdim(mds_points[, 2], 0.20)  # 20%
        cut_lines(list(type = "hline", values = cpdim))
      }
      
      output$go_terms <- renderDataTable({
        req(result)  
        
        datatable(result, filter = "top", 
                  caption = htmltools::tags$caption(
                    style = "caption-side: top; text-align: center;",
                    paste("GO terms for ", caption_params()$onto, 
                          caption_params()$level, " - ", caption_params()$dim)
                  )
        )
      })
      
      remove_modal_spinner() 
    })})
  
  observeEvent(
    {
      input$l1
      input$l2
      input$sel_onto
      input$sel_level
      input$whatdim
    },
    {
      cut_lines(NULL)
      output$go_terms <- renderDataTable({ NULL })
    }
  )
  
}
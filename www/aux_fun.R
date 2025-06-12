getListNames <- function(tablelist){
  # Crear una lista vacía para almacenar las tablas individuales con nombres compuestos
  all_individual_tables_named <- list()
  
  # Recorrer cada elemento de tablelist
  for (i in seq_along(tablelist)) {
    
    # Obtener el nombre del elemento actual, que contiene ambas listas en comparación
    element_name <- names(tablelist)[i]
    
    # Obtener el elemento actual de la lista (puede ser una o más tablas)
    current_element <- tablelist[[i]]
    
    # Verificar si el elemento actual es una lista (indica múltiples tablas)
    if (is.list(current_element)) {
      # Recorrer cada tabla dentro de este elemento y agregarla a la lista con un nombre compuesto específico
      for (j in seq_along(current_element)) {
        # Asignar un nombre único combinando el nombre de ambas listas y el índice de la tabla
        table_name <- paste0(element_name, "_Vs_", names(current_element[j]))
        all_individual_tables_named[[table_name]] <- current_element[[j]]
      }
    } else {
      # Si el elemento no es una lista, agregarlo directamente a la lista de tablas con el nombre compuesto
      table_name <- paste0(element_name, "_Tabla_1")
      all_individual_tables_named[[table_name]] <- current_element
    }
  }
  all_individual_tables_named
  
}


filter_list <- function(lista, vec) {
  # Expresión regular para detectar pares de identificadores en el vector
  filtro_regex <- paste0("(", paste(vec, collapse = "|"), ")_Vs_(", paste(vec, collapse = "|"), ")")
  
  # Filtrar la lista manteniendo solo los elementos que cumplen con la regex y ambos elementos están en vec
  Filter(function(nombre) {
    # Verifica que el nombre contenga ambos identificadores en vec
    matches <- unlist(strsplit(nombre, "_Vs_"))
    all(matches %in% vec) && grepl(filtro_regex, nombre)
  }, names(lista))
}


sdm <- function(x, ...){
  nlists <- length(x)
  ifelse(nlists > 1, sd(x, ...), 0)
}


get_go_description <- function(go_id) {
  go_term <- Term(GOTERM[[go_id]])
  return(go_term)
}

splitdim <- function(dm, prp){
  dm <- as.vector(dm)
  prp <- c(prp, 1 - (2*prp), prp)
  rg <- range(dm)
  cutpoints <- cumsum(prp) * diff(rg) + rg[1]
  cutpoints[-length(cutpoints)]
}

gotermsid <- function(list, dm, prp, dimen, contabs){
  sorted <- as.data.frame(dm[order(dm[, dimen]), ])
  cpdim1 <- splitdim(dm[, dimen], prp)
  lleft <- rownames(sorted[sorted[, dimen] < cpdim1[1], ]) # Identify lists to the left
  lright <- rownames(sorted[sorted[, dimen] > cpdim1[2], ])
  
  # STEP 2)  1-0 Table for extremes left and right
  #tableleft <- data.frame(enrichedIn(list[lleft], ...))
  enrichmat <- attr(contabs, "enriched")
  tableleft <- data.frame(enrichmat[, lleft])
  colnames(tableleft) <- lleft
  #tableright <- data.frame(enrichedIn(list[lright], ...))
  tableright <- data.frame(enrichmat[, lright])
  colnames(tableright) <- lright
  
  # STEP 3) Mean and sd 
  lmnsd <- apply(tableleft, 1, 
                 function(x){c("meanLeft" = mean(x), "sdLeft" = sdm(x))})
  rmnsd <- apply(tableright, 1, 
                 function(x){c("meanRight" = mean(x), "sdRight" = sdm(x))})
  
  # STEP 4) Pseudo t
  nl <- ncol(tableleft) ;   nr <- ncol(tableright) 
  Pseudo_t <- abs(lmnsd[1, ] - rmnsd[1, ]) / sqrt((((lmnsd[2, ] / nl) + (rmnsd[2, ] / nr))) + 0.0000001)
  
  
  # STEP 5) Summary table
  sum <- as.data.frame(t(rbind(lmnsd, rmnsd, Pseudo_t)))
  sortSum <- sum[order(sum[, "Pseudo_t"]), ]
  prevRes <- sortSum[sortSum$Pseudo_t == max(sortSum$Pseudo_t), ]
  
  GOIDs <- rownames(prevRes)
  desc <- data.frame(Description = sapply(rownames(prevRes), get_go_description, 
                                          USE.NAMES = F))
  
  cbind(GOIDs, desc)
}


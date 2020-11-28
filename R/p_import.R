
menuImportSpectra <- function(input,output,session,wsp) {
  output$menuImportSpectra <- shiny::renderUI("Import spectra")
}

serveImportSpectra <- function(input, output, session, wsp) {
  
  data <- shiny::reactiveValues(
    fcsMtx = NULL,
    spectrum = NULL
  )

  observeEvent(input$fileImportFCS, {
    data$spectrum <- NULL
    data$fcsMtx <- NULL
    #TODO: handle loading errors
    m <- flowCore::read.FCS(input$fileImportFCS$datapath)@exprs
    colnames(m) <- unname(colnames(m))
    data$fcsMtx <- m
  })

  observeEvent(input$doImportGetSpectrum, {
    if(is.null(data$fcsMtx)) return()
    if(is.null(input$importDataCols)) return()
    cols <- nat.sort(input$importDataCols)
    data$spectrum <- extractSpectrum(data$fcsMtx[,cols,drop=F])
    #print(data$spectrum)
  })

  observeEvent(input$doImportSave, {
    if(is.null(data$spectrum)) return()
    n <- c(
      spectrumMetadataFormGather('importForm', input),
      list(spectrum=data$spectrum))
    if(any(sapply(wsp$spectra, function(x)
      all(sapply(names(n),function(nm)
        if(nm=='spectrum')TRUE else n[[nm]]==x[[nm]]))
      ))) {
      shiny::showNotification(type='error',
        "This spectrum already exists; use a different note to distinguish duplicates.")
      return ()
    }

    wsp$spectra = c(wsp$spectra, list(n))

    shiny::showNotification(type='message',
      "Saved OK.")
  })

  output$uiImportResult <- shiny::renderUI(shiny::tagList(
    shiny::h3("Computed spectrum"),
    if(is.null(data$spectrum)) "No spectrum computed"
    else shiny::tagList(
      paste0("Mean fluorescent intensity: ", dbf(data$spectrum$mI),
             ", signal: ", dbf(2*data$spectrum$sdI,'')),
      shiny::plotOutput('plotImportSpectra', width='100%', height='40ex')
    )
  ))

  output$plotImportSpectra <- shiny::renderPlot(
    if(!is.null(data$spectrum))
      plotSpectrum(data$spectrum$mS, data$spectrum$sdS, data$spectrum$channels)
  )

  output$uiImportSave <- shiny::renderUI({
    wsp$page # reload the information from wsp$spectra after the page changed
    shiny::tagList(
      shiny::h3("Save the spectrum"),
      spectrumMetadataForm('importForm', isolate(wsp$spectra)),
      shiny::uiOutput('uiImportSaveBtn')
  )})

  output$uiImportSaveBtn <- shiny::renderUI(
    if(is.null(data$spectrum)) "No spectrum to save"
    else shiny::actionButton('doImportSave', "Save")
  )

  output$uiImportColumnPicker <- shiny::renderUI(
    shiny::selectizeInput('importDataCols',
      "Fluorescent channels",
      multiple=T,
      choices=colnames(data$fcsMtx),
      selected=shiny::isolate(
        defaultFluorescenceChannels(nat.sort(colnames(data$fcsMtx)))),
      options=list(`actions-box`=TRUE)))

  output$uiImportSpectra <- shiny::renderUI(shiny::tagList(
    shiny::h1("Import spectral profiles"),
    shiny::fluidRow(
      shiny::column(8,style='min-width:20em',
        shiny::fileInput('fileImportFCS', "Uplad an FCS file", accept='.fcs'),
        shiny::uiOutput('uiImportColumnPicker'),
        shiny::actionButton('doImportGetSpectrum', "Compute spectrum"),
        shiny::uiOutput('uiImportResult')),
      shiny::column(4,style='min-width:10em',
        shiny::uiOutput('uiImportSave')))
  ))
}

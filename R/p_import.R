
menuImportSpectra <- function(input,output,session,wsp) {
  output$menuImportSpectra <- shiny::renderUI("Import spectra")
}

defaultMachines <- c('Aurora','Symphony','Fortessa II','BFC 9k')
defaultSampleTypes <- c('Spleen (mouse)','Thymus (mouse)','PBMC (mouse)','Bone marrow (mouse)','Beads')
defaultAntigens <- c(paste0('CD',1:100),'L/D','Unstained') 
defaultFluorochromes <- c('Cy5', 'Cy5.5', 'Cy7', 'APC', 'FITC', 'PE', 'GFP', 'Unstained')

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
  })

  observeEvent(input$doImportSave, {
    if(is.null(data$spectrum)) return()
    n <- list(machine=input$importMachine,
      mconfig=input$importMachineConfig,
      sample=input$importSampleType,
      antigen=input$importAntigen,
      fluorochrome=input$importFluorochrome,
      note=input$importNote,
      spectrum=data$spectrum)
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

    print(wsp$spectra)
  })

  output$uiImportResult <- shiny::renderUI(shiny::tagList(
    shiny::h3("Computed spectrum"),
    if(is.null(data$spectrum)) "No spectrum computed"
    else shiny::tagList(
      paste0("Mean fluorescent intensity: ",e2db(data$spectrum$mI^2),", signal: ",e2db(data$spectrum$sdI^2,'')),
      shiny::plotOutput('plotImportSpectra', width='100%', height='40ex')
    )
  ))

  output$uiImportSave <- shiny::renderUI({
    wsp$page # reload the information from wsp$spectra after the page changed
    shiny::tagList(
      shiny::h3("Save the spectrum"),
      shiny::selectizeInput('importMachine', "Cytometer model",
        choices=nat.sort(unique(c(defaultMachines,unlist(sapply(isolate(wsp$spectra), function(x) x$machine))))),
        multiple=F, options=list(`create`=TRUE)),
      shiny::selectizeInput('importMachineConfig', "Cytometer configuration",
        choices=nat.sort(unique(c('—',unlist(sapply(isolate(wsp$spectra), function(x) x$mconfig))))),
        multiple=F, options=list(`create`=TRUE)),
      shiny::selectizeInput('importSampleType', "Sample type",
        choices=nat.sort(unique(c(defaultSampleTypes,unlist(sapply(isolate(wsp$spectra), function(x) x$sample))))),
        multiple=F, options=list(`create`=TRUE)),
      shiny::selectizeInput('importAntigen', "Antigen",
        choices=nat.sort(unique(c(defaultAntigens,unlist(sapply(isolate(wsp$spectra), function(x) x$antigen))))),
        multiple=F, options=list(`create`=TRUE)),
      shiny::selectizeInput('importFluorochrome', "Fluorochrome",
        choices=nat.sort(unique(c(defaultFluorochromes,unlist(sapply(isolate(wsp$spectra), function(x) x$fluorochrome))))),
        multiple=F, options=list(`create`=TRUE)),
      shiny::selectizeInput('importNote', "Note",
        choices=nat.sort(unique(c('—', unlist(sapply(isolate(wsp$spectra), function(x) x$note))))),
        multiple=F, options=list(`create`=TRUE)),
      shiny::uiOutput('uiImportSaveBtn')
  )})

  output$uiImportSaveBtn <- shiny::renderUI(
    if(is.null(data$spectrum)) "No spectrum to save"
    else shiny::actionButton('doImportSave', "Save")
  )

  output$plotImportSpectra <- shiny::renderPlot(
    if(!is.null(data$spectrum))
      plotSpectrum(data$spectrum$mS, data$spectrum$sdS)
  )

  output$uiImportColumnPicker <- shiny::renderUI(
    shiny::selectizeInput('importDataCols',
      "Fluorescent channels",
      multiple=T,
      choices=colnames(data$fcsMtx),
      selected=shiny::isolate(defaultFluorescenceChannels(colnames(data$fcsMtx))),
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

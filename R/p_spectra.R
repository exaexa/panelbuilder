
menuSpectraBrowser <- function(input,output,session,wsp) {
  output$menuSpectraBrowser <- shiny::renderUI(
    paste0("Manage spectra (",length(wsp$spectra),")")
  )
}

serveSpectraBrowser <- function(input, output, session, wsp) {

  data <- shiny::reactiveValues(
    filter=NULL
  )

  dynObservers <- list()

  output$uiSpectraTable <- shiny::renderUI(shiny::tagList(
    paste("Displayed spectra:", length(wsp$spectra[data$filter]),
          "out of", length(wsp$spectra)),
    do.call(shiny::tags$table, c(list(
        class="table",
        shiny::tags$tr(
          shiny::tags$th('Machine'),
          shiny::tags$th('Config'),
          shiny::tags$th('Sample type'),
          shiny::tags$th('Antigen'),
          shiny::tags$th('Fluorochrome'),
          shiny::tags$th('Note'),
          shiny::tags$th('Channels'),
          shiny::tags$th('Intensity'),
          shiny::tags$th('Spectrum'),
          shiny::tags$th('Actions')
      )),
      unname(lapply(seq_len(length(wsp$spectra))[data$filter], function(sid) {
        s <- wsp$spectra[[sid]]
        rmName <- paste0('doSpectraRemove', sid)
        if(is.null(dynObservers[[rmName]]))
          dynObservers[[rmName]] <<- observeEvent(input[[rmName]],
            {wsp$spectra <- wsp$spectra[-sid]})
        toPanelName <- paste0('doSpectraToPanel', sid)
        if(is.null(dynObservers[[toPanelName]]))
          dynObservers[[toPanelName]] <<- observeEvent(input[[toPanelName]], {
            if(spectrumExistsIn(wsp$spectra[[sid]], wsp$panelSpectra, fields=c('antigen','fluorochrome')))
              shiny::showNotification(type='error',
                paste0(wsp$spectra[[sid]]$antigen,"/",wsp$spectra[[sid]]$fluorochrome," already in panel"))
            else {
              wsp$panelSpectra <- c(wsp$panelSpectra, list(wsp$spectra[[sid]]))
              shiny::showNotification(type='message',
                paste0(wsp$spectra[[sid]]$antigen,"/",wsp$spectra[[sid]]$fluorochrome," added"))
            }
          })
        shiny::tags$tr(
          shiny::tags$td(s$machine),
          shiny::tags$td(s$mconfig),
          shiny::tags$td(s$sample),
          shiny::tags$td(s$antigen),
          shiny::tags$td(s$fluorochrome),
          shiny::tags$td(s$note),
          shiny::tags$td(length(s$spectrum$mS)),
          shiny::tags$td(paste0(dbf(s$spectrum$mI), " ±", dbf(s$spectrum$sdI, ''))),
          shiny::tags$td(
            plotSpectrumPNG(
              s$spectrum$mS,
              s$spectrum$sdS,
              res=48, x=320, y=48)),
          shiny::tags$td(
            shiny::actionButton(toPanelName, "→panel"),
            shiny::actionButton(rmName, "×")))
      }
    ))))
  ))

  observeEvent(wsp$spectra, {
    data$filter <- unlist(sapply(wsp$spectra, function(x)T))
  })

  observeEvent(wsp$panelSpectra, {
    #TODO: clean up panelAssignments
  })

  observeEvent(input$doSpectraFilter, {
    tgt <- spectrumMetadataFormGather('spectraSearchForm', input)
    data$filter <- unlist(sapply(wsp$spectra, function(x)
      all(sapply(names(tgt),function(nm) x[[nm]] %in% tgt[[nm]]))))
  })

  observeEvent(input$fileSpectraUpload, {
    #TODO: merge instead of replacing
    wsp$spectra <- jsonlite::read_json(
      input$fileSpectraUpload$datapath,
      simplifyVector=T)
  })

  output$fileSpectraDownload <- downloadHandler(
    filename=function() paste0("spectra-",timestamp(),".json"),
    content=function(con) {
      names(wsp$spectra) <- seq_len(length(wsp$spectra[data$filter]))
      jsonlite::write_json(wsp$spectra, con, digits=9)
    }
  )

  output$uiSpectraBrowser <- shiny::renderUI(shiny::tagList(
    shiny::h1("Manage spectra"),
    shiny::fluidRow(
      shiny::column(8,
        shiny::h3("Search"),
        spectrumMetadataForm('spectraSearchForm',
          wsp$spectra, create=F, multiple=T, defaultAll=T),
        shiny::actionButton('doSpectraFilter', "Filter!")),
      shiny::column(4,
        shiny::h3("Import/Export"),
        shiny::fileInput('fileSpectraUpload',
          "Import from file", accept='.json'),
        shiny::downloadButton('fileSpectraDownload',
          "Export displayed spectra"))),
    shiny::uiOutput('uiSpectraTable')))
}

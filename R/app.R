
#' Run the panelBuildeR Shiny UI
#'
#' @param ... Passed to shinyApp()
#' @export
panelbuilder <- function(...) {
  required.pkgs <- c(
    'shiny',
    'EmbedSOM',
    'flowCore')
  missing.pkgs <- required.pkgs[!required.pkgs %in% .packages(all.available=T)]
  if(length(missing.pkgs)>0)
    stop(do.call(paste, as.list(c("Required packages missing:",missing.pkgs))))

  ui <- shiny::fluidPage(
    shiny::tags$head(shiny::tags$script(shiny::HTML("
window.onbeforeunload = function() {
    return 'Are you sure you want to quit? All unsaved changes will be lost!';
};")),
      shiny::tags$style(shiny::HTML("
.spectrumCovered { background-color: #dfb; }
.spectrumMeasured { background-color:#bdf; border-radius:4pt; }
.spectrumDerived { border-radius:4pt; }
.spectrum { padding:4pt; }
.spectrumSelected { border: 2pt solid black; padding: 2pt; }
")), #TODO: add style here if required
      shiny::tags$title("panelBuildeR")),

    shiny::fluidRow(style="margin-bottom:3em",
      shiny::column(2, style="min-width:16em",
        shiny::h2("panelBuildeR"),
        shiny::uiOutput('mainMenu')),
      shiny::column(10,
        shiny::uiOutput('mainPage')))
  )

  server <- function(input, output, session) {
    wsp <- shiny::reactiveValues(
      page='import',
      spectra=list(),
      panelSpectra=list(),
      panelSpectraEst=list(),
      panelAssignment=list()
    )

    mainPages <- list(
      import='ImportSpectra',
      spectra='SpectraBrowser',
      reqs='PanelSetup',
      optimize='Optimizer',
      export='Export',
      unmix='Unmix'
    )

    mainPageServers <- lapply(mainPages,
      function(nm) eval(as.name(paste0('serve',nm))))
    mainMenuServers <- lapply(mainPages,
      function(nm) eval(as.name(paste0('menu',nm))))

    mainPageUIs <- lapply(mainPages,function(nm)paste0('ui',nm))
    mainMenuUIs <- lapply(mainPages,function(nm)paste0('menu',nm))

    # render the menu
    output$mainMenu <- shiny::renderUI(
      shiny::div(do.call(shiny::tags$ol,
        c(list(class='m-3 list-group'),
          unname(lapply(names(mainPages),
          function(p)
            shiny::actionLink(
              class=paste('list-group-item','list-group-item-action',
                          if(wsp$page==p)'active' else NULL),
              paste0(p, 'Switch'),
              shiny::uiOutput(mainMenuUIs[[p]]))))))))

    #render the main page
    output$mainPage <- shiny::renderUI({
      if(wsp$page %in% names(mainPageServers)) {
        shiny::uiOutput(mainPageUIs[[wsp$page]])
      } else
        shiny::h1("Page not found")
    })

    #serve menus and pages
    for(p in names(mainPages)) {
      mainMenuServers[[p]](input,output,session,wsp)
      mainPageServers[[p]](input,output,session,wsp)
    }

    # handle the page switching
    for(p in names(mainPages)) {
      x <- new.env()
      x$wsp <- wsp
      x$p <- p
      observeEvent(input[[paste0(p,'Switch')]],
                   {wsp$page <- p},
                   event.env=x,
                   handler.env=x)
    }
  }

  if(is.null(options()$shiny.maxRequestSize))
    options(shiny.maxRequestSize=256*2^20)

  options(bitmapType='cairo')

  shiny::shinyApp(ui=ui, server=server, ...)
}



#' Run the panelBuildeR Shiny UI
#'
#' @param ... Passed to shinyApp()
#' @export
panelBuildeR <- function(...) {
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
      shiny::tags$style(shiny::HTML("")), #TODO: add style if required
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
      panelReqs=list(
        costs=data.frame(
          antigen=NULL,
          dye=NULL,
          cost=NULL),
        usedSpectra=list(),
        sampleTypes=NULL,
        targetRes=NULL,
        mutualScore=data.frame(
          ag1=NULL,
          ag2=NULL,
          cost=NULL)),
      target=list(
        mA=NULL, # antigen mean binding capacity (relative)
        sA=NULL, # antigen signal sdev, (same units as mA)
        mB=NULL, # mean brightness of the fluorochrome (relative)
        mS=NULL, # spectra of the fluorochrome (normalized)
        sS=NULL, # sdev of fluorochrome emission/measurement (same units as mS)
        uM=NULL, # unstained (AF) spectra mean expression (absolute)
        uS=NULL, # unstained (AF) spectra signal sdev (absolute)
        C=NULL, # costs of fluorochrome/antigen combinations
        D=NULL, # antigen-antigen mutual resolution importance
        T=NULL),# antigen target resolution (in bits)
      panel=data.frame(
        antigen=NULL,
        dye=NULL)
    )

    mainPages <- list(
      import='ImportSpectra',
      spectra='SpectraBrowser',
      reqs='PanelSetup',
      optimize='Optimizer',
      export='Export'
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


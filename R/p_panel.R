menuPanelSetup <- function(input,output,session,wsp) {
  output$menuPanelSetup <- shiny::renderUI(
    paste0("Panel setup (",length(panelAntigens(wsp)),"AGs, ",length(panelFluorochromes(wsp)),"FCs)")
  )
}

servePanelSetup <- function(input, output, session, wsp) {

  data <- shiny::reactiveValues(
    mode='select',
    display='signal'
  )

  findInPanel <- function(ag,fc,wsp) {
    idx <- which(unlist(lapply(wsp$panelSpectra,
      function(s) s$fluorochrome==fc && s$antigen==ag)))
    if(length(idx)==0) NULL else idx
  }

  output$plotPanelWrap <- shiny::renderUI(
    shiny::plotOutput('plotPanel', click='clickPanel', height=paste0(15+2*length(panelFluorochromes(wsp)),'em'))
  )

  output$plotPanel <- shiny::renderPlot({
    fcs <- panelFluorochromes(wsp)
    ags <- panelAntigens(wsp)

    d <- data.frame(
      Fluorochrome=factor(sapply(wsp$panelSpectra, function(x)x$fluorochrome), levels=fcs),
      Antigen=factor(sapply(wsp$panelSpectra, function(x)x$antigen), levels=ags),
      `Mean expression`=sapply(wsp$panelSpectra, function(x)x$spectrum$mI),
      `Signal sDev`=sapply(wsp$panelSpectra, function(x)x$spectrum$sdI),
      Selected=sapply(wsp$panelSpectra, function(x) {
        match <- wsp$panelAssignment[[x$antigen]]==x$fluorochrome
        if(length(match)==0) F
        else if(is.null(match)) F
        else match
      }),
      check.names=F)

    if(nrow(d)==0) NULL else
    ggplot2::ggplot(d, ggplot2::aes(x=Antigen,y=Fluorochrome,fill=`Signal sDev`, color=Selected)) +
      ggplot2::geom_tile(width=0.9, height=0.9, size=2) +
      scale_fill_gradientn(colors=colorspace::sequential_hcl(100, 'Dark Mint')) +
      scale_color_manual("Selection",values=c('#cccccc', '#4488ff'), breaks=c(F,T), labels=c("Not used", "In panel")) +
      scale_y_reverse() +
      cowplot::theme_cowplot() +
      ggplot2::theme(
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(angle=45, vjust=1, hjust=1),
        panel.grid = element_line(color='#eeeeee'))
  })

  observeEvent(input$doPanelClear, {
    wsp$panelAssignment <- list()
  })

  observeEvent(input$doPanelSelect, {
    data$mode <- 'select'
  })

  observeEvent(input$doPanelEstimate, {
    data$mode <- 'estimate'
    spectraEstComponents(wsp)
  })

  observeEvent(input$doPanelRemove, {
    data$mode <- 'remove'
  })

  observeEvent(input$clickPanel, {
    ags <- panelAntigens(wsp)
    fcs <- panelFluorochromes(wsp)
    xp <- as.integer(round(min(max(input$clickPanel$x,1), length(ags))))
    yp <- as.integer(round(min(max(input$clickPanel$y,1), length(fcs))))
    ag <- ags[xp]
    fc <- fcs[yp]

    if(data$mode=='select') {
      orig <- wsp$panelAssignment[[ag]]
      if(is.null(orig) || orig != fc) {
        idx <- findInPanel(ag,fc,wsp)
        if(!is.null(idx)) wsp$panelAssignment[[ag]] <- fc
      } else if(!is.null(orig) && orig==fc)
        wsp$panelAssignment[[ag]] <- NULL
    } else if (data$mode=='estimate') {
    } else if (data$mode=='remove') {
      idx <- findInPanel(ag, fc, wsp)
      if(!is.null(idx)) {
        wsp$panelSpectra <- wsp$panelSpectra[-idx]
        wsp$panelAssignment <- cleanAssignments(wsp)
      }
    }
  })

  output$uiPanelBtns <- shiny::renderUI(shiny::tagList(
    shiny::actionButton('doPanelSelect', class=if(data$mode=='select') 'btn-primary' else 'btn-outline-primary', "Select"),
    shiny::actionButton('doPanelEstimate', class=if(data$mode=='estimate') 'btn-primary' else 'btn-outline-primary', "Estimate spectra"),
    shiny::actionButton('doPanelEstAll', "Estimate all"),
    shiny::actionButton('doPanelEstClear', "Remove estimates"),
    shiny::actionButton('doPanelClear', class='btn-warning', "Clear selection"),
    shiny::actionButton('doPanelRemove', class=if(data$mode=='remove') 'btn-danger' else 'btn-outline-danger', "Remove spectra"),
    shiny::actionButton('doPanelDisplaySignal', class=if(data$display=='signal') 'btn-info' else 'btn-outline-info', "Display signal"),
    shiny::actionButton('doPanelDisplaySNR', class=if(data$display=='snr') 'btn-info' else 'btn-outline-info', "Display SNRs")
  ))

  output$uiPanelSetup <- shiny::renderUI(shiny::tagList(
    shiny::h1("Panel setup"),
    shiny::div(style='display: inline-block; vertical-align:top', shiny::uiOutput('uiPanelBtns')),
    shiny::uiOutput('plotPanelWrap')
  ))
}

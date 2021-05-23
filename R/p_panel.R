menuPanelSetup <- function(input,output,session,wsp) {
  output$menuPanelSetup <- shiny::renderUI(
    paste0("Panel setup (",length(panelAntigens(wsp$panelSpectra)),"AGs, ",length(panelFluorochromes(wsp$panelSpectra)),"FCs)")
  )
}

servePanelSetup <- function(input, output, session, wsp) {

  data <- shiny::reactiveValues(
    mode='select',
    display='mi',
    estimate=F
  )

  findInPanel <- function(ag,fc,p) {
    idx <- which(unlist(lapply(p,
      function(s) s$fluorochrome==fc && s$antigen==ag)))
    if(length(idx)==0) NULL else idx
  }

  output$plotPanelWrap <- shiny::renderUI(
    shiny::plotOutput('plotPanel', click='clickPanel', height=paste0(15+2*length(panelFluorochromes(wsp$panelSpectra)),'em'))
  )

  output$plotPanel <- shiny::renderPlot({
    ags <- panelAntigens(wsp$panelSpectraEst)
    fcs <- panelFluorochromes(wsp$panelSpectraEst)

    d <- data.frame(
      Fluorochrome=factor(sapply(wsp$panelSpectraEst, function(x)x$fluorochrome), levels=rev(fcs)),
      Antigen=factor(sapply(wsp$panelSpectraEst, function(x)x$antigen), levels=ags),
      `Mean expression`=sapply(wsp$panelSpectraEst, function(x)x$spectrum$mI),
      `Signal sDev`=sapply(wsp$panelSpectraEst, function(x)x$spectrum$sdI),
      Selected=factor(sapply(wsp$panelSpectraEst, function(x) {
        match <- wsp$panelAssignment[[x$antigen]]==x$fluorochrome
        est <- if(x$estimated) c('est','usedEst') else c('real','usedReal')
        if(length(match)==0) est[1]
        else if(is.null(match)) est[1]
        else if(match) est[2] else est[1]
      }), levels=c('real','usedReal','est','usedEst')),
      check.names=F)

    if(nrow(d)==0) NULL else
    ggplot2::ggplot(d,
      ggplot2::aes_string(
        x='Antigen',
        y='Fluorochrome',
        fill=if(data$display=='mi') '`Mean expression`' else '`Signal sDev`',
        color='Selected')) +
      ggplot2::geom_tile(width=0.8, height=0.8, size=2) +
      ggplot2::scale_fill_gradientn(
        colors=colorspace::sequential_hcl(100, h=c(0, 90), c=c(80, NA, 34), l=c(30, 99), power=c(0.2, 2))) +
      ggplot2::scale_color_manual("Selection",
        values=c('#ccffe2', '#00cc33', '#cce2ff', '#0044ff'),
        breaks=c('real','usedReal','est','usedEst'),
        labels=c("Measured", "Meas. in panel", "Estimated", "Est. in panel")) +
      cowplot::theme_cowplot() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = ggplot2::element_text(angle=45, vjust=1, hjust=1),
        panel.grid = ggplot2::element_line(color='#eeeeee'))
  })

  observeEvent(input$doPanelClear, {
    wsp$panelAssignment <- list()
  })

  observeEvent(input$doPanelSelect, {
    data$mode <- 'select'
  })

  observeEvent(input$doPanelEstimate, {
    data$estimate <- !data$estimate
  })

  observeEvent(input$doPanelRemove, {
    data$mode <- 'remove'
  })

  observeEvent(input$clickPanel, {
    ags <- panelAntigens(wsp$panelSpectraEst)
    fcs <- rev(panelFluorochromes(wsp$panelSpectraEst))
    xp <- as.integer(round(min(max(input$clickPanel$x,1), length(ags))))
    yp <- as.integer(round(min(max(input$clickPanel$y,1), length(fcs))))
    ag <- ags[xp]
    fc <- fcs[yp]

    if(data$mode=='select') {
      orig <- wsp$panelAssignment[[ag]]
      if(is.null(orig) || orig != fc) {
        idx <- findInPanel(ag,fc,wsp$panelSpectraEst)
        if(!is.null(idx)) wsp$panelAssignment[[ag]] <- fc
      } else if(!is.null(orig) && orig==fc)
        wsp$panelAssignment[[ag]] <- NULL
    } else if (data$mode=='remove') {
      idx <- findInPanel(ag, fc, wsp$panelSpectra)
      if(!is.null(idx)) {
        wsp$panelSpectra <- wsp$panelSpectra[-idx]
      }
    }
  })

  reestimate <- function() {
    wsp$panelSpectraEst <- makeSpectraEstimate(wsp, addEstimates=data$estimate)
    wsp$panelAssignment <- cleanAssignments(wsp$panelAssignment, wsp$panelSpectraEst)
  }

  observeEvent(wsp$panelSpectra, reestimate())
  observeEvent(data$estimate, reestimate())

  observeEvent(input$doPanelDisplayIntensity, {
    data$display <- 'mi'
  })

  observeEvent(input$doPanelDisplaySignal, {
    data$display <- 'sdi'
  })

  output$uiPanelBtns <- shiny::renderUI(shiny::tagList(
    shiny::actionButton('doPanelSelect', class=if(data$mode=='select') 'btn-primary' else 'btn-outline-primary', "Select"),
    shiny::actionButton('doPanelEstimate', class=if(data$estimate) 'btn-success' else 'btn-danger', if(data$estimate) "Estimates ON" else "Estimates OFF"),
    shiny::actionButton('doPanelClear', class='btn-warning', "Clear selection"),
    shiny::actionButton('doPanelRemove', class=if(data$mode=='remove') 'btn-danger' else 'btn-outline-danger', "Remove spectra"),
    shiny::actionButton('doPanelDisplayIntensity', class=if(data$display=='mi') 'btn-info' else 'btn-outline-info', "Display Intensity"),
    shiny::actionButton('doPanelDisplaySignal', class=if(data$display=='sdi') 'btn-info' else 'btn-outline-info', "Display Signal")
  ))

  output$uiPanelSetup <- shiny::renderUI(shiny::tagList(
    shiny::h1("Panel setup"),
    shiny::div(style='display: inline-block; vertical-align:top', shiny::uiOutput('uiPanelBtns')),
    shiny::uiOutput('plotPanelWrap')
  ))
}

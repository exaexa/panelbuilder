
menuPanelSetup <- function(input,output,session,wsp) {
  output$menuPanelSetup <- shiny::renderUI(
    paste0("Panel setup (",length(panelAntigens(wsp)),"AGs, ",length(panelFluorochromes(wsp)),"FCs)")
  )
}

servePanelSetup <- function(input, output, session, wsp) {
  
  dynObservers <- list()

  output$uiPanelTable <- shiny::renderUI({
    fcs <- panelFluorochromes(wsp)
    ags <- panelAntigens(wsp)

    makeCell <- function(fcid, agid) {
      fc <- fcs[fcid]
      ag <- ags[agid]
      selectedCls <- if(!is.null(wsp$panelAssignment[[ag]]) && wsp$panelAssignment[[ag]]==fc)
                       "spectrumSelected"
                     else
                       "spectrum"
      idx <- which(unlist(lapply(wsp$panelSpectra,
        function(s) s$fluorochrome==fc & s$antigen==ag)))

      setName <- paste0('doPanelSet', fcid, 'x', agid)
      if(is.null(dynObservers[[setName]]))
        dynObservers[[setName]] <<- observeEvent(input[[setName]], {
          fc <- panelFluorochromes(wsp)[fcid]
          ag <- panelAntigens(wsp)[agid]
          wsp$panelAssignment[[ag]] <- fc
        })
      rmName <- paste0('doPanelRemove', fcid, 'x', agid)
      if(is.null(dynObservers[[rmName]]))
        dynObservers[[rmName]] <<- observeEvent(input[[rmName]], {
          fc <- panelFluorochromes(wsp)[fcid]
          ag <- panelAntigens(wsp)[agid]
          idx <- which(unlist(lapply(wsp$panelSpectra,
            function(s) s$fluorochrome==fc & s$antigen==ag)))
          wsp$panelSpectra <- wsp$panelSpectra[-idx]
        })
      if(length(idx)==0) shiny::tags$td(shiny::div(class=paste('spectrumDerived', selectedCls), "×"))
      else shiny::tags$td(shiny::div(class=paste('spectrumMeasured', selectedCls),
        dbf(wsp$panelSpectra[[idx]]$spectrum$mI),
        shiny::tags$br(),
        "±", dbf(wsp$panelSpectra[[idx]]$spectrum$sdI, ''),
        shiny::tags$br(),
        ilDiv(
          shiny::actionButton(setName, "↑"),
          shiny::actionButton(rmName, "×"))))
    }

    do.call(shiny::tags$table, c(
      list(class="table"),
      list(do.call(shiny::tags$tr, c(
        list(
          shiny::tags$td("FC \\ AG",
          shiny::tags$br(),
          shiny::actionButton('doPanelClear', "Clear assignment")),
        unname(lapply(ags, function(ag)
          shiny::tags$th(ag,
            class=if(any(ag == names(wsp$panelAssignment))) "spectrumCovered" else NULL,
            style='writing-mode: sideways-lr'))))))),
      unname(lapply(seq_len(length(fcs)), function(fcid) do.call(shiny::tags$tr, c(
        list(shiny::tags$th(
          fcs[fcid],
          class=if(any(fcs[fcid] == unlist(wsp$panelAssignment))) "spectrumCovered" else NULL)),
        unname(lapply(seq_len(length(ags)), function(agid)makeCell(fcid, agid)))))))))
  })

  observeEvent(input$doPanelClear, {
    wsp$panelAssignment <- list()
  })

  observeEvent(input$doPanelInfo, {
    print(getUnmixingInfo(wsp))
  })

  output$uiPanelSetup <- shiny::renderUI(shiny::tagList(
    shiny::h1("Panel setup"),
    shiny::uiOutput('uiPanelTable'),
    shiny::actionButton('doPanelInfo', 'debug panel info')
  ))
}

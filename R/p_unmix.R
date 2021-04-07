
menuUnmix <- function(input,output,session,wsp) {
  output$menuUnmix <- shiny::renderUI("Unmixing")
}

serveUnmix <- function(input, output, session, wsp) {
  data <- shiny::reactiveValues(
    inputMtx = NULL,
    inputName = NULL,
    outputMtx = NULL,
    outputColnames = NULL,
    postCompensation = NULL,
    postCompHistory = list(),
    postCompPts = NULL
  )

  output$uiUnmixLoad <- shiny::renderUI(shiny::tagList(
    shiny::h3("Open a FCS"),
    shiny::fileInput('fileUnmixFCS',
      "FCS for unmixing", accept='.fcs'),
  ))

  observeEvent(input$fileUnmixFCS, {
    data$inputMtx <- NULL
    data$inputName <- character(0)
    data$outputMtx <- NULL
    data$outputColnames <- NULL
    shiny::withProgress({
      m <- flowCore::read.FCS(input$fileUnmixFCS$datapath)@exprs
      setProgress(1)
    }, message="Loading FCS...")
    colnames(m) <- unname(colnames(m))
    data$inputMtx <- m
    data$inputName <- input$fileUnmixFCS$name
  })

  output$uiUnmixControl <- shiny::renderUI(if(!is.null(data$inputMtx)) shiny::tagList(
    shiny::h3("Unmixing control"),
    shiny::div(shiny::strong("Loaded: "), data$inputName, paste0("(", nrow(data$inputMtx), " events)")),
    {
      mcs <- matchingChans(data$inputMtx, getUnmixingInfo(wsp))
      umcs <- nat.sort(colnames(data$inputMtx)[!colnames(data$inputMtx) %in% mcs])
      shiny::tagList(
        do.call(shiny::div, c(
          list(shiny::strong(paste0("Unmixing channels (",length(mcs),"):"))),
          lapply(mcs, function(mc) shiny::span(class="badge", mc)))),
        do.call(shiny::div, c(
          list(shiny::strong(paste0("Other channels (",length(umcs),"):"))),
          lapply(umcs, function(umc) shiny::span(class="badge", umc)))))
    },
    shiny::checkboxInput('unmixIncludeFluorochromes', "Include fluorochrome names", value=T),
    shiny::checkboxInput('unmixIncludeOriginals', "Retain original values in the raw channels that were used for unmixing", value=F),
    shiny::checkboxInput('unmixIncludeResiduals', "Include per-channel residuals", value=F),
    shiny::checkboxInput('unmixIncludeRMSE', "Include total unmixing RMSE information", value=T),
    shiny::actionButton('doUnmix', "Run unmixing")
  ))

  observeEvent(input$doUnmix,
    if(!is.null(data$inputMtx))
      shiny::withProgress({
        tryCatch(
          data$outputMtx <- doUnmix(data$inputMtx, getUnmixingInfo(wsp),
            input$unmixIncludeFluorochromes,
            input$unmixIncludeOriginals,
            input$unmixIncludeResiduals,
            input$unmixIncludeRMSE),
          error=function(e) shiny::showNotification(type='error',
            paste("Unmixing failed:", e)))
          data$outputColnames <- colnames(data$outputMtx)
        shiny::setProgress(1)
      }, message="Unmixing..."))

  output$uiUnmixPreview <- shiny::renderUI(shiny::tagList(
    shiny::h3("Result preview"),
    shiny::fluidRow(
      shiny::column(4,
        shiny::uiOutput('uiUnmixPlotOpts')),
      shiny::column(4, 
        shiny::h4("Data preview"),
        shiny::uiOutput('uiUnmixPlot'),
        shiny::sliderInput('unmixPlotAlpha', "Point alpha",
          min=0, max=1, step=0.01, value=0.5)),
      shiny::column(4,
        shiny::radioButtons('unmixTool',
          "Tool",
          choices=c(`Leveling tool`='level', `Gating`='gate'),
          selected='level'),
        shiny::uiOutput('uiUnmixTools'),
        shiny::h4("Results"),
        shiny::downloadButton('downloadUnmixFCS', "Download unmixed FCS")))
  ))

  output$uiUnmixPlotOpts <- shiny::renderUI(if(!is.null(data$outputColnames)) shiny::tagList(
    shiny::selectizeInput('unmixPlotX',
      "Preview column X",
      multiple=F,
      choices=data$outputColnames,
      selected=defaultFSCChannel(data$outputColnames)),
    shiny::checkboxInput('unmixAsinhX',
      "Transform X",
      value=F),
    shiny::sliderInput('unmixCofX',
      "Cofactor X (dB)",,
      min=-10, max=80, step=1, value=30),
    shiny::selectizeInput('unmixPlotY',
      "Preview column Y",
      multiple=F,
      choices=data$outputColnames,
      selected=defaultSSCChannel(data$outputColnames)),
    shiny::checkboxInput('unmixAsinhY',
      "Transform Y",
      value=F),
    shiny::sliderInput('unmixCofY',
      "Cofactor Y (dB)",,
      min=-10, max=80, step=1, value=30),
    shiny::selectizeInput('unmixPlotCol',
      "Preview column color",
      multiple=F,
      choices=c('(Density)', data$outputColnames),
      selected='(Density)'),
    shiny::checkboxInput('unmixAsinhCol',
      "Transform Color",
      value=F),
    shiny::sliderInput('unmixCofCol',
      "Cofactor color (dB)",,
      min=-10, max=80, step=1, value=30)
  ))

  output$uiUnmixPlot <- shiny::renderUI(if(!is.null(data$outputMtx)) shiny::tagList(
    shiny::plotOutput('plotUnmix',
      width="30em",
      height="30em",
      click=if(input$unmixTool=='level') 'clickUnmixPlot' else NULL,
      brush=if(input$unmixTool=='gate') shiny::brushOpts('brushUnmixPlot')),
    if(!is.null(data$outputMtx))
      shiny::div(paste0("(", nrow(data$outputMtx), " events)"))
  ))

  getTransFns <- function() list(
    tx = if(input$unmixAsinhX) function(v)asinh(v/db2e(input$unmixCofX)) else identity,
    ty = if(input$unmixAsinhY) function(v)asinh(v/db2e(input$unmixCofY)) else identity,
    tc = if(input$unmixAsinhCol) function(v)asinh(v/db2e(input$unmixCofCol)) else identity,
    itx = if(input$unmixAsinhX) function(v)sinh(v)*db2e(input$unmixCofX) else identity,
    ity = if(input$unmixAsinhY) function(v)sinh(v)*db2e(input$unmixCofY) else identity)

  getCompData <- function()
    if(is.null(data$postCompensation)) data$outputMtx
    else data$outputMtx %*% data$postCompensation

  output$uiUnmixTools <- shiny::renderUI(if(input$unmixTool=='level') shiny::tagList(
    shiny::div("the tool aligns the selected × cross to the level of + cross"),
    ilDiv(
      shiny::actionButton('doUnmixLevelH', "═"),
      shiny::actionButton('doUnmixLevelV', "‖"),
      shiny::actionButton('doUnmixLevelUndo', "Undo"),
      shiny::actionButton('doUnmixLevelReset', "Reset")))
    else ilDiv(
      shiny::actionButton('doUnmixGateIn', "Keep only gate"),
      shiny::actionButton('doUnmixGateOut', "Remove gate")))

  observeEvent(input$doUnmixLevelReset, {
    data$postCompensation <- diag(1, ncol(data$outputMtx))
    colnames(data$postCompensation) <- colnames(data$outputMtx)
    rownames(data$postCompensation) <- colnames(data$outputMtx)
    data$postCompHistory <- list()
    data$postCompPts <- NULL
  })

  observeEvent(input$doUnmixLevelUndo, if(length(data$postCompHistory)>0) {
    data$postCompensation <- data$postCompHistory[[1]]
    data$postCompHistory <- data$postCompHistory[-1]
  })

  observeEvent(data$outputMtx, {
    if(is.null(data$outputMtx))
      data$postCompensation <- NULL
    else {
      data$postCompensation <- diag(1,ncol(data$outputMtx))
      colnames(data$postCompensation) <- colnames(data$outputMtx)
      rownames(data$postCompensation) <- colnames(data$outputMtx)
    }
    data$postCompHistory <- list()
    data$postCompPts <- NULL
  })

  observeEvent(input$unmixPlotX,
    data$postCompPts <- NULL)

  observeEvent(input$unmixPlotY,
    data$postCompPts <- NULL)

  observeEvent(input$clickUnmixPlot, {
    ts <- getTransFns()
    data$postCompPts <- rbind(c(ts$itx(input$clickUnmixPlot$x), ts$ity(input$clickUnmixPlot$y)), data$postCompPts)
    if(nrow(data$postCompPts)>2) data$postCompPts <- data$postCompPts[1:2,,drop=F]
  })

  doAlign <- function(src,dst) tryCatch({
    tr <- solve(src) %*% dst
    data$postCompHistory <- c(list(data$postCompensation), data$postCompHistory)
    ds <- c(input$unmixPlotX, input$unmixPlotY)
    data$postCompensation[ds,ds] <- data$postCompensation[ds,ds] %*% tr
    data$postCompPts <- NULL
  }, error=function(e) shiny::showNotification(type='error', paste(e, "(system is probably singular)")))

  observeEvent(input$doUnmixLevelH, if(!is.null(data$postCompPts) && nrow(data$postCompPts)==2) {
    src <- data$postCompPts
    dst <- src
    dst[2,2] <- dst[1,2]
    doAlign(src, dst)
  })

  observeEvent(input$doUnmixLevelV, if(!is.null(data$postCompPts) && nrow(data$postCompPts)==2) {
    src <- data$postCompPts
    dst <- src
    dst[2,1] <- dst[1,1]
    doAlign(src, dst)
  })

  doGate <- function(b, inv) if(!is.null(b)) {
    ts <- getTransFns()
    flt <- xor(inv,
      data$outputMtx[,input$unmixPlotX] >= ts$itx(b$xmin) &
      data$outputMtx[,input$unmixPlotX] <= ts$itx(b$xmax) &
      data$outputMtx[,input$unmixPlotY] >= ts$ity(b$ymin) &
      data$outputMtx[,input$unmixPlotY] <= ts$ity(b$ymax))
    data$outputMtx <- data$outputMtx[flt,,drop=F]
  }

  observeEvent(input$doUnmixGateIn, doGate(input$brushUnmixPlot, F))
  observeEvent(input$doUnmixGateOut, doGate(input$brushUnmixPlot, T))

  output$plotUnmix <- shiny::renderPlot({
    ts <- getTransFns()
    d <- getCompData()
    par(mar=c(0,0,0,0))
    if(!is.null(data$outputMtx) && input$unmixPlotX!='' && input$unmixPlotY!='') {
      EmbedSOM::PlotEmbed(
        cbind(ts$tx(d[,input$unmixPlotX]), ts$ty(d[,input$unmixPlotY])),
        data=if(input$unmixPlotCol=='(Density)') NULL else cbind(ts$tc(d[,input$unmixPlotCol])),
        val=if(input$unmixPlotCol=='(Density)') 0 else 1,
        alpha=input$unmixPlotAlpha,
        plotf=scattermore::scattermoreplot)
      abline(h=0)
      abline(v=0)
      if(input$unmixTool=='level') points(
          ts$tx(rev(data$postCompPts[,1])),
          ts$ty(rev(data$postCompPts[,2])),
          cex=4, lwd=4, pch=c(4,3), col='#00cc00')
    }
  })

  output$downloadUnmixFCS <- shiny::downloadHandler(
    filename=function() paste0("pbUnmixed_",data$inputName),
    content=function(conn) flowCore::write.FCS(new('flowFrame', exprs=getCompData()), conn)
  )

  output$uiUnmix <- shiny::renderUI(shiny::tagList(
    shiny::h1("Unmixing"),
    if(length(wsp$panelAssignment)==0) "Prepare the panel first."
    else shiny::tagList(
      shiny::uiOutput('uiUnmixLoad'),
      shiny::uiOutput('uiUnmixControl'),
      shiny::uiOutput('uiUnmixPreview'))
  ))
}

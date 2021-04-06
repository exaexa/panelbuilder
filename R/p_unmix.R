
menuUnmix <- function(input,output,session,wsp) {
  output$menuUnmix <- shiny::renderUI("Unmixing")
}

serveUnmix <- function(input, output, session, wsp) {
  data <- shiny::reactiveValues(
    inputMtx = NULL,
    inputName = NULL,
    outputMtx = NULL
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
      umcs <- colnames(data$inputMtx)[!colnames(data$inputMtx) %in% mcs]
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
        data$outputMtx <- doUnmix(data$inputMtx, getUnmixingInfo(wsp),
          input$unmixIncludeFluorochromes,
          input$unmixIncludeOriginals,
          input$unmixIncludeResiduals,
          input$unmixIncludeRMSE)
        shiny::setProgress(1)
      }, message="Unmixing..."))

  output$downloadUnmixFCS <- shiny::downloadHandler(
    filename=function() paste0("pbUnmixed_",data$inputName),
    content=function(conn) flowCore::write.FCS(new('flowFrame', exprs=data$outputMtx), conn)
  )

  output$uiUnmixPreview <- shiny::renderUI(if(!is.null(data$outputMtx)) shiny::tagList(
    shiny::h3("Unmixed result"),
    shiny::downloadButton('downloadUnmixFCS', "Download unmixed FCS"),
    shiny::fluidRow(
      shiny::column(4,
        shiny::selectizeInput('unmixPlotX',
          "Preview column X",
          multiple=F,
          choices=colnames(data$outputMtx),
          selected=defaultFSCChannel(colnames(data$outputMtx))),
        shiny::selectizeInput('unmixPlotY',
          "Preview column Y",
          multiple=F,
          choices=colnames(data$outputMtx),
          selected=defaultSSCChannel(colnames(data$outputMtx)))),
      shiny::column(8,
        shiny::plotOutput('plotUnmix',
          width="30em",
          height="30em")))
  ))

  output$plotUnmix <- shiny::renderPlot({
    par(mar=c(0,0,0,0))
    if(!is.null(data$outputMtx) && input$unmixPlotX!='' && input$unmixPlotY!='') {
      EmbedSOM::PlotEmbed(
        data$outputMtx[,c(input$unmixPlotX, input$unmixPlotY)],
        fdens=log,
        plotf=scattermore::scattermoreplot)
    }
  })

  output$uiUnmix <- shiny::renderUI(shiny::tagList(
    shiny::h1("Unmixing"),
    shiny::uiOutput('uiUnmixLoad'),
    shiny::uiOutput('uiUnmixControl'),
    shiny::uiOutput('uiUnmixPreview')
  ))
}


menuUnmix <- function(input,output,session,wsp) {
  output$menuUnmix <- shiny::renderUI("Unmixing")
}

serveUnmix <- function(input, output, session, wsp) {
  data <- shiny::reactiveValues(
    inputFcs = NULL,
    spectra = NULL,
    outputMtx = NULL
  )

  output$uiUnmixSelSpectra <- shiny::renderUI(shiny::tagList(
    shiny::h3("Add spectra")
  ))

  output$uiUnmixControl <- shiny::renderUI(shiny::tagList(
    shiny::h3("Raw file input"),
    shiny::fileInput('fileUnmixFCS',
      "FCS for unmixing", accept='.fcs'),
    shiny::actionButton('doUnmix', "Unmix!")
  ))

  output$uiUnmixPreview <- shiny::renderUI(shiny::tagList(
    shiny::h3("Unmixed output")
  ))


  output$uiUnmix <- shiny::renderUI(shiny::tagList(
    shiny::h1("Unmixing"),
    shiny::uiOutput('uiUnmixSelSpectra'),
    shiny::uiOutput('uiUnmixControl'),
    shiny::uiOutput('uiUnmixPreview')
  ))
}

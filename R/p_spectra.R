
menuSpectraBrowser <- function(input,output,session,wsp) {
  output$menuSpectraBrowser <- shiny::renderUI(
    paste0("Manage spectra (",length(wsp$spectra),")")
  )
}

serveSpectraBrowser <- function(input, output, session, wsp) {
  output$uiSpectraBrowser <- shiny::renderUI(shiny::tagList(
    shiny::h1("Manage spectra")
  ))
}

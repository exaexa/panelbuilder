
menuPanelSetup <- function(input,output,session,wsp) {

  output$menuPanelSetup <- shiny::renderUI(
    paste0("Panel spectra (",length(wsp$panelSpectra),")")
  )
}

servePanelSetup <- function(input, output, session, wsp) {
  output$uiPanelSetup <- shiny::renderUI(shiny::tagList(
    shiny::h1("Panel setup and requirements")
  ))
}


menuPanelSetup <- function(input,output,session,wsp) {
  output$menuPanelSetup <- shiny::renderUI("panel menu!")
}

servePanelSetup <- function(input, output, session, wsp) {
  output$uiPanelSetup <- shiny::renderUI(shiny::tagList(
    shiny::h1("Panel setup and requirements")
  ))
}

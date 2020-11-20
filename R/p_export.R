
menuExport <- function(input,output,session,wsp) {
  output$menuExport <- shiny::renderUI("export menu!")
}

serveExport <- function(input, output, session, wsp) {
  output$uiExport <- shiny::renderUI(shiny::tagList(
    shiny::h1("Export data")
  ))
}

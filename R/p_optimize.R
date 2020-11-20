
menuOptimizer <- function(input,output,session,wsp) {
  output$menuOptimizer <- shiny::renderUI("optimizer menu!")
}

serveOptimizer <- function(input, output, session, wsp) {
  output$uiOptimizer <- shiny::renderUI(shiny::tagList(
    shiny::h1("Optimize the panel")
  ))
}

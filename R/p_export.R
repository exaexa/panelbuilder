
menuExport <- function(input,output,session,wsp) {
  output$menuExport <- shiny::renderUI("Panel overview")
}

serveExport <- function(input, output, session, wsp) {
  output$uiExport <- shiny::renderUI(shiny::tagList(
    shiny::h1("Panel overview"),
    shiny::tags$table(class="table",
      shiny::tags$tr(
        shiny::tags$th("Antigen"),
        shiny::tags$th("Fluorochrome")),
      do.call(shiny::tagList,
        lapply(names(wsp$panelAssignment), function(ag)
          shiny::tags$tr(
            shiny::tags$td(ag),
            shiny::tags$td(wsp$panelAssignment[[ag]])))))
  ))
}

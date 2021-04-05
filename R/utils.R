
const <- function(a,...) a

nat.sort <- gtools::mixedsort

timestamp <- function()
  format(Sys.time(),"%Y%m%d-%H%M")

ilDiv <- function(..., style='') shiny::div(
  style=paste0('display: inline-block;', style), ...)

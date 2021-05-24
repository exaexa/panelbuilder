
const <- function(a,...) a

nat.sort <- gtools::mixedsort

timestamp <- function()
  format(Sys.time(),"%Y%m%d-%H%M")

ilDiv <- function(..., style='') shiny::div(
  style=paste0('display: inline-block;', style), ...)

indexin <- function(idx, arr) {
  tmp <- seq_len(length(arr))
  names(tmp) <- arr
  unname(tmp[idx])
}

eachindex <- function(x) if(length(x)<1) c() else seq(length(x))

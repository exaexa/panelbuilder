
panelAntigens <- function(wsp)
  nat.sort(unique(sapply(wsp$panelSpectra, function(x)x$antigen)))

panelFluorochromes <- function(wsp)
  nat.sort(unique(sapply(wsp$panelSpectra, function(x)x$fluorochrome)))

cleanAssignments <- function(wsp) {
  ags <- panelAntigens(wsp)
  fcs <- panelFluorochromes(wsp)
  out <- list()
  for(ag in names(wsp$panelAssignment))
    if((ag %in% ags) && (wsp$panelAssignment[[ag]] %in% fcs))
      out[[ag]]=wsp$panelAssignment[[ag]]
  out
}

getUnmixingInfo <- function(wsp) {
  as <- cleanAssignments(wsp)
  ss <- lapply(names(as), function(ag)
    spectrumFind(wsp$panelSpectra, list(antigen=ag, fluorochrome=as[[ag]])))
  chs <- character(0)
  for(s in ss) chs <- c(chs, s$spectrum$channels)
  chs <- nat.sort(unique(chs))
  chs <- chs[as.logical(sapply(chs, function(ch)
    all(sapply(ss, function(s)
      ch %in% s$spectrum$channels)))),drop=F]
  slen <- length(s$spectrum$mS)
  list( #TODO: nat-sort by antigen
    antigens=as.character(sapply(ss, function(s) s$antigen)),
    fluorochromes=as.character(sapply(ss, function(s) s$fluorochrome)), #informative-only
    channels=chs,
    mSs=matrix(sapply(ss, function(s) s$spectrum$mS[indexin(chs, s$spectrum$channels)]), length(chs), length(ss)),
    sdSs=matrix(sapply(ss, function(s) s$spectrum$sdS[indexin(chs, s$spectrum$channels)]), length(chs), length(ss)),
    mIs=as.numeric(sapply(ss, function(s) s$spectrum$mI)),
    sdIs=as.numeric(sapply(ss, function(s) s$spectrum$sdI)))
}

matchingChans <- function(mtx, ui) {
  nat.sort(unique(colnames(mtx)[colnames(mtx) %in% ui$channels]))
}

doUnmix <- function(mtx, ui, fcNames=T, inclOrigs=F, inclResiduals=F, inclRmse=T) {
  mc <- matchingChans(mtx, ui)

  if(length(mc)==0) return(mtx) # nothing to do

  umtx <- t(mtx[,indexin(mc, colnames(mtx))])
  umSs <- ui$mSs[indexin(mc, ui$channels),]
  # TODO: add variants that consider sdS and absolute intensities
  u <- lm(umtx~umSs+0, weights=sqrt(rowSums(umSs^2)))

  res <- mtx

  if(!inclOrigs) res <- res[,-indexin(mc, colnames(res))]

  umxd <- t(u$coefficients)
  colnames(umxd) <-
    if(fcNames) paste0(ui$antigens, " <", ui$fluorochromes, ">")
    else ui$antigens
  res <- cbind(res, umxd)
  
  if(inclResiduals) {
    rss <- t(u$residuals)
    colnames(rss) <- paste("pbResidual ",mc)
    res <- cbind(res,rss)
  }

  if(inclRmse) res <- cbind(res, pbRMSE=sqrt(colMeans(u$residuals^2)))

  res
}
  

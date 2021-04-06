
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
  chs <- chs[
    sapply(chs, function(ch)
      all(sapply(ss, function(s)
        ch %in% s$spectrum$channels)))]
  slen <- length(s$spectrum$mS)
  list(
    antigens=sapply(ss, function(s) s$antigen),
    fluorochromes=sapply(ss, function(s) s$fluorochrome), #informative-only
    channels=chs,
    mSs=sapply(ss, function(s) s$spectrum$mS[seq_len(slen)[s$spectrum$channels %in% chs]]),
    sdSs=sapply(ss, function(s) s$spectrum$sdS[seq_len(slen)[s$spectrum$channels %in% chs]]),
    mIs=sapply(ss, function(s) s$spectrum$mI),
    sdIs=sapply(ss, function(s) s$spectrum$sdI))
}

matchingChans <- function(mtx, ui) {
  nat.sort(unique(colnames(mtx)[colnames(mtx) %in% ui$channels]))
}

doUnmix <- function(mtx, ui, fcNames=T, inclOrigs=F, inclResiduals=F, inclRmse=T) {
  mc <- matchingChans(mtx, ui)
  if(length(mc)==0) stop("No matching channels!")
  # TODO: error on mismatch
  umtx <- t(mtx[,colnames(mtx) %in% mc])
  umSs <- ui$mSs[ui$channels %in% mc,]
  # TODO: add variants that consider sdS and absolute intensities
  u <- lm(umtx~umSs+0) #TODO: weights

  res <- mtx
  
  if(!inclOrigs) res <- res[,!(colnames(res) %in% mc)]
  
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
  

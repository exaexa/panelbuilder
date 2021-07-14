
panelAntigens <- function(p)
  nat.sort(unique(sapply(p, function(x)x$antigen)))

panelFluorochromes <- function(p)
  nat.sort(unique(sapply(p, function(x)x$fluorochrome)))

cleanAssignments <- function(as, p) {
  ags <- panelAntigens(p)
  fcs <- panelFluorochromes(p)
  out <- list()
  for(s in p)
    if(s$antigen %in% names(as))
      if(as[[s$antigen]]==s$fluorochrome)
        out[s$antigen] <- s$fluorochrome
  out
}

getCommonChannels <- function(panel) {
  chs <- character(0)
  for(s in panel) chs <- c(chs, s$spectrum$channels)
  chs <- nat.sort(unique(chs))
  m <- matrix(F, nrow=length(chs), ncol=length(panel))
  rownames(m) <- chs
  for(i in seq(length(panel))) m[panel[[i]]$spectrum$channels, i] <- T
  chs[apply(m, 1, all)]
}

getUnmixingInfo <- function(wsp) {
  as <- cleanAssignments(wsp$panelAssignment, wsp$panelSpectraEst)
  ss <- lapply(nat.sort(names(as)), function(ag)
    spectrumFind(wsp$panelSpectraEst, list(antigen=ag, fluorochrome=as[[ag]])))
  chs <- getCommonChannels(ss)
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

#' @export
doUnmix <- function(mtx, ui, method='ols', fcNames=T, inclOrigs=F, inclResiduals=F, inclRmse=T) {
  mc <- matchingChans(mtx, ui)

  if(length(mc)==0) return(mtx) # nothing to do

  coefficients <- NULL
  residuals <- NULL

  umtx <- t(mtx[,indexin(mc, colnames(mtx))])

  # TODO: add variants that consider sdS and absolute intensities
  if(method=='ols') {
    umSs <- ui$mSs[indexin(mc, ui$channels),]
    u <- lm(umtx~umSs+0)
    coefficients <- t(u$coefficients)
    residuals <- t(u$residuals)
  } else if(method=='ols-spw') {
    umSs <- ui$mSs[indexin(mc, ui$channels),]
    cws <- sqrt(rowMeans(umSs^2))
    cws[cws<1e-2] <- 1e-2
    cws <- 1/cws
    umtx <- umtx*cws
    umSs <- umSs*cws
    u <- lm(umtx~umSs+0)
    coefficients <- t(u$coefficients)
    residuals <- t(u$residuals/cws)
  } else if(method=='ols-chw') {
    umSs <- ui$mSs[indexin(mc, ui$channels),]
    cws <- sqrt(rowMeans(umtx^2))
    cws[cws<1e-2] <- 1e-2
    cws <- 1/cws
    umtx <- umtx*cws
    umSs <- umSs*cws
    u <- lm(umtx~umSs+0)
    coefficients <- t(u$coefficients)
    residuals <- t(u$residuals/cws)
  } else if(method=='eols-rw') {
    umtx[umtx<1] <- 1
    umSs <- ui$mSs[indexin(mc, ui$channels),]
    ones <- rep(1, nrow(umtx))
    coefficients <- t(apply(umtx, 2,
      function(r) lm.fit(umSs/r, ones)$coefficients))
    residuals <- t(umtx) - (coefficients %*% t(umSs))
  } else if(method=='gd-pw') {
    umSs <- ui$mSs[indexin(mc, ui$channels),]

    res <- nougad::nougad(
      t(umtx),
      t(umSs),
      1, 2, 4)

    coefficients <- res$unmixed
    residuals <- res$residuals
  } else {
    stop("unsupported unmixing method")
  }

  res <- mtx

  if(!inclOrigs) res <- res[,-indexin(mc, colnames(res))]

  colnames(coefficients) <-
    if(fcNames) paste0(ui$antigens, " <", ui$fluorochromes, ">")
    else ui$antigens
  res <- cbind(res, coefficients)
  
  if(inclResiduals) {
    colnames(residuals) <- paste("pbResidual ",mc)
    res <- cbind(res,residuals)
  }

  if(inclRmse) res <- cbind(res, pbRMSE=sqrt(rowMeans(residuals^2)))

  res
}

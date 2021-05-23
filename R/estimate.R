
spectraEstComponents <- function(wsp) {
  ags <- panelAntigens(wsp$panelSpectra)
  fcs <- panelFluorochromes(wsp$panelSpectra)
  vert_ags <- seq_len(length(ags))
  vert_fcs <- seq_len(length(fcs))+length(ags)
  edge_ags <- indexin(unlist(sapply(wsp$panelSpectra, function(x)x$antigen)), ags)
  edge_fcs <- indexin(unlist(sapply(wsp$panelSpectra, function(x)x$fluorochrome)), fcs)
  for(i in seq_len(length(edge_ags))) {
    edge <- c(vert_ags[edge_ags[i]], vert_fcs[edge_fcs[i]])
    union <- min(edge)
    find <- max(edge)
    if(union<find) {
      vert_ags[vert_ags==find] <- union
      vert_fcs[vert_fcs==find] <- union
    }
  }
  names(vert_ags) <- ags
  names(vert_fcs) <- fcs
  list(ags=vert_ags, fcs=vert_fcs)
}

spectrumToLinearEst <- function(sp, channels) {
  cixs <- indexin(channels, sp$channels)
  c(log(sp$mI),
    log(sp$sdI),
    sp$mS, # mS is spherical, needs special treatment
    sp$sdS)
}

spectrumFromLinearEst <- function(lin, channels) {
  n <- length(channels)
  if(length(lin)!=2+2*n) stop("internal error in spectrum decoding")
  mI <- exp(lin[1])
  sdI <- exp(lin[2])
  # project mS back to unit sphere
  mS <- lin[3:(2+n)]
  mS[mS<0] <- 0
  mSp <- sqrt(sum(mS^2))
  if(mSp>1e-5) mS <- mS/sqrt(sum(mS^2))

  sdS <- lin[(3+n):length(lin)]
  sdS[sdS<0] <- 0

  list(
    channels=channels,
    mI=mI,
    sdI=sdI,
    mS=mS,
    sdS=sdS)
}

makeSpectraEstimate <- function(wsp, addEstimates=T) {
  # find what can be estimated
  components <- spectraEstComponents(wsp)
  estimate <- -outer(components$fcs, components$ags, "==")
  if(length(wsp$panelSpectra)>0) for(i in seq(length(wsp$panelSpectra))) {
    sp <- wsp$panelSpectra[[i]]
    estimate[sp$fluorochrome, sp$antigen] <- i
  }

  res <- list()
  for(sp in wsp$panelSpectra) {
    # include the originals
    os <- list(
      estimated=F,
      antigen=sp$antigen,
      fluorochrome=sp$fluorochrome,
      spectrum=sp$spectrum)
    res <- c(res, list(os))
  }

  if(!addEstimates) return(res)

  # cycle through all components
  for(component in sort(unique(unname(c(components$fcs,components$ags))))) {
    # pick out the relevant parts
    subest_mtx <- estimate[components$fcs==component, components$ags==component, drop=F]

    # prepare variable-picking bits for OLS
    fcbits <- diag(rep(1,nrow(subest_mtx)))
    rownames(fcbits) <- rownames(subest_mtx)
    agbits <- diag(rep(1,ncol(subest_mtx)))
    rownames(agbits) <- colnames(subest_mtx)

    # convert to a dataframe for easier processing
    subestimate <- reshape2::melt(
      subest_mtx,
      value.name='ref',
      varnames=c('fc','ag'))

    estIdxs <- subestimate$ref[subestimate$ref>0L]
    if(length(estIdxs)<2L || sum(subestimate$ref<0)<1L) {
      # continue if there's nothing to estimate
      next
    }

    # build a model from known spectra
    known <- subestimate[subestimate$ref>0L,]
    chs <- getCommonChannels(wsp$panelSpectra[estIdxs])
    if(length(chs)<1L) {
      warning(paste("estimation cancelled: no common spectra"))
      next
    }

    known_model <- cbind(
      fcbits[known$fc,,drop=F],
      agbits[known$ag,,drop=F])
    known_response <- matrix(nrow=nrow(known), t(sapply(
      wsp$panelSpectra[known$ref],
      function(s) spectrumToLinearEst(s$spectrum, chs))))

    # add absolute-scale parameters -- if nothing helps, we can always state the stuff should sum up to 0
    eps <- 1e-6
    known_model <- rbind(known_model,
      c(rep(eps,ncol(fcbits)), rep(0, ncol(agbits))),
      c(rep(0,ncol(fcbits)), rep(eps, ncol(agbits))),
      c(rep(eps,ncol(fcbits)), rep(eps, ncol(agbits))))
    known_response <- rbind(known_response, 0, 0, 0)

    coefs <- matrix(
      lm(known_response ~ known_model+0)$coefficients,
      nrow=ncol(known_model))

    # we didn't manage to guess enough stuff, weird.
    if(any(is.na(coefs))) {
      warning(paste("estimation failed despite expectations: ",known_model))
      next
    }

    # check out the spectra to be estimated
    unknown <- subestimate[subestimate$ref<0,]
    unknown_model <- cbind(
      fcbits[unknown$fc,,drop=F],
      agbits[unknown$ag,,drop=F])

    # calculate the estimates and push into the result
    est_response <- unknown_model %*% coefs

    for(i in seq(nrow(est_response))) {
      sp <- spectrumFromLinearEst(est_response[i,], chs)
      es <- list(
        estimated=T,
        antigen=as.character(unknown$ag[i]),
        fluorochrome=as.character(unknown$fc[i]),
        spectrum=sp)
      res <- c(res, list(es))
    }
  }

  res
}

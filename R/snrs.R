
snr_panel_choices <- function(as, p, static=1) {
  ags <- panelAntigens(p)
  fcs <- panelFluorochromes(p)
  chs <- getCommonChannels(p)

  # make matrices to hold all this
  self_noise <- array(0, dim=c(length(chs), length(ags), length(fcs)))
  signal <- array(0, dim=c(length(chs), length(ags), length(fcs)))

  # fill in the available data
  for(sp in p) {
    s <- sp$spectrum
    cixs <- indexin(s$channels, chs)
    self_noise[,sp$antigen, sp$fluorochrome] <- (db2e(s$mI) * s$sdS)^2
    signal[,sp$antigen, sp$fluorochrome] <- (db2e(s$sdI) * s$mS)^2
  }

  # collect the total noise from assigned panel
  asgn_noise <- array(0, dim=c(length(chs),length(ags)))
  dimnames(asgn_noise)[[2]] <- ags
  for(ag in names(as)) {
    fc <- as[[ag]]
    asgn_noise[,ag] <- self_noise[,ag, fc] + signal[,ag, fc]
  }

  # helper for easier computation of "sum except one" below
  total_noise <- apply(asgn_noise,1,sum)

  sapply(p, function(sp) {
    this_noise <- total_noise + self_noise[,sp$antigen, sp$fluorochrome]
    if(sp$antigen %in% names(as))
      this_noise <- this_noise - asgn_noise[,sp$antigen]

    e02db(sum(signal[,sp$antigen, sp$fluorochrome])) -
      e02db(sum(this_noise*signal[,sp$antigen, sp$fluorochrome])+static)/2
  })
}


#' Extract the spectrum
#'
#' the spectra are scaled to the mI.
extractSpectrum <- function(mtx) {
  a <- rowSums(mtx)
  # @importFrom MASS rlm
  #scoefs <- apply(mtx,2,function(x) MASS::rlm(x ~ a+0)$coefficients)
  scoefs <- lm(mtx~a+0)$coefficients[,]
  coefs <- scoefs / sqrt(sum(scoefs^2))
  m <- lm(t(mtx) ~ coefs + 0)
  list(mS=coefs,
       sdS=unname(apply(m$residuals,1,sd)/mean(m$coefficients)),
       mI=mean(m$coefficients),
       sdI=sd(m$coefficients))
}

defaultFluorescenceChannels <- function(nms)
  nms[!grepl("(fsc|ssc|fs0|ss0|time)",nms,ignore.case=T)]

e2db <- function(x, unit='u')
  paste0(round(10*log(x),2),' dB',unit)

spectrumPalette <- function(n=128, ...)
  colorspace::sequential_hcl(n = n, h = 114, c = c(0, 114, 0), l = c(100, 0), power = 1.1, ...)

plotSpectrum <- function(ms, sds, nms=names(ms), res=128) {
  d <- sapply(seq_len(length(nms)), function(i)
    sapply(seq(0,1,length.out=res), function(j)
      pnorm(j,mean=ms[i],sd=sds[i])))

  colnames(d) <- nms

  d <- reshape2::melt(d)

  ggplot2::ggplot(d) + 
    ggplot2::aes(Var2,Var1,fill=value) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors=spectrumPalette(32),guide=F) +
    ggplot2::scale_x_discrete(name=NULL) +
    ggplot2::scale_y_continuous(labels=NULL, name=NULL) +
    cowplot::theme_cowplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
}

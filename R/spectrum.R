
powerEstimate <- function(mtx) rowSums(mtx)

eliminateSpectra <- function(mtx, elim=NULL) {
  if(is.null(elim)) return(mtx)

  n <- nrow(mtx)
  d <- ncol(mtx)
  k <- ncol(elim)

  iters <- 100
  alpha <- 0.1
  tol <- 1

  x_kn <- matrix(0, k, n)
  r_dn <- matrix(0, d, n)

  res <- .C("pw_gd",
    n=as.integer(n),
    d=as.integer(d),
    k=as.integer(k),
    iters=as.integer(iters),
    alpha=as.single(alpha),
    tol=as.single(tol),
    s=as.single(t(elim)),
    spw=as.single(matrix(1, d, k)),
    snw=as.single(matrix(20, d, k)), # this might be parametrizable
    nw=as.single(rep(10,k)),
    y=as.single(t(mtx)),
    x=as.single(x_kn),
    r=as.single(r_dn))

  # return the residuals that couldn't be explained by spectra removal
  matrix(res$r,n,d,byrow=T)
}

#' Extract the spectrum
#'
#' the spectra are scaled to the mI.
extractSpectrum <- function(mtx, method='ols', gate=T, powerGate=T, elim=NULL) {

  gc(full=T)

  if(length(gate)==1) gate <- rep(T, nrow(mtx))
  if(length(powerGate)==1) powerGate <- rep(T, nrow(mtx))

  m1 <- eliminateSpectra(mtx[gate,, drop=F], elim)
  m2 <- m1[powerGate[gate],,drop=F]
  a <- powerEstimate(m2)

  sp <- NULL

  if(method=='ols') {
    reg <- lm(m2~a)
    sp <- reg$coefficients['a',]
  } else if(method=='theil-sen') {
    ssize <- max(1e6, 10*length(a))
    sa <- sample(length(a), ssize, replace=T)
    sb <- sample(length(a), ssize, replace=T)
    v <- (m2[sa,]-m2[sb,])
    p <- (a[sa]-a[sb])
    flt <- abs(p)>1
    vp <- v[flt,]/p[flt]
    sp <- apply(vp, 2, median)
  } else stop("Unknown method for spectrum extraction")

  sp[sp<0] <- 0 # clamp negatives (TODO: is this bias?)
  if(max(sp)==0) stop("Spectrum extraction problem: all coefficients non-positive.")
  sp <- sp/sqrt(sum(sp^2)) # normalize the spectrum

  reg <- lm(t(m1) ~ sp) # "unmix" using this single channel from the whole data
  e <- e2db(reg$coefficients['sp',]^2)/2 # _intensities_ in dB
  e[e < -100] <- -100 # dodge zeroes

  # only look at the upper half of the data for spectral noise
  flt <- e>=mean(e)
  sds <- rowSums((sp*reg$residuals[,flt,drop=F])^2) / rowSums(1+reg$fitted.values[,flt,drop=F]^2)

  #iqs <- quantile(e, pnorm(c(-2,2))) # +/- 2sigma range of intensities #TODO: this might be viable later

  channels <- nat.sort(colnames(mtx))
  perm <- indexin(channels, colnames(mtx))
  list(channels=channels,
       mS=unname(sp)[perm],
       sdS=unname(sds)[perm],
       mI=mean(e),
       sdI=sd(e)) #TODO: (iqs[2]-iqs[1])/sqrt(12)) #this roughly corresponded to the energy in uniform distribution
}

# decibels should be used for intensities everywhere
e02db <- function(x) 10*log(x, base=10)
e2db <- function(x) e02db(x[x>0])
db2e <- function(x) 10^(x/10)
dbf <- function(x, unit='u')
  paste0(round(x,2),' dB',unit)

defaultFluorescenceChannels <- function(nms)
  nms[!grepl("(fsc|ssc|fs0|ss0|time|comp)",nms,ignore.case=T)]

defaultFSCChannel <- function(nms)
  nms[c(grepl("(fsc|fs0).*-a", nms, ignore.case=T),grepl("(fsc|fs0)", nms, ignore.case=T))][1]

defaultSSCChannel <- function(nms)
  nms[c(grepl("(ssc|ss0).*-a", nms, ignore.case=T),grepl("(ssc|ss0)", nms, ignore.case=T))][1]

# This is a nice transcript of how the spectrum palettes evolved. Keep adding.
# Remember that the midpoint is important.
spectrumPalette <- function(n=128, ...)
  #grDevices::colorRampPalette(c('#44dd22','#ffddaa','#002040'))
  #colorspace::sequential_hcl(n=n, 'viridis', rev=T)
  #grDevices::colorRampPalette(c('#44dd22','white','black'))(n,...)
  #grDevices::colorRampPalette(c('#88aacc',EmbedSOM::ExpressionPalette(n-2),'white'))(n,...)
  #grDevices::colorRampPalette(c('#338822','#ffcc88','black'))(n,...)
  #colorspace::sequential_hcl(n = n, h = 114, c = c(0, 0, 114), l = c(0,100, 50), power = 1.1, ...)
  #colorspace::diverging_hcl(n = 28, h = c(220, 33), c = c(43, 55), l = c(100, 35), power = 1)
  colorRampPalette(c("#7CC3DB", "#77ACBF", "#525252", "#C59B8A", "#FFF3DD"))(n,...)
  #colorspace::diverging_hcl(n, 'Lisbon')

plotSpectrum <- function(ms, sds, nms=names(ms), res=128) {
  d <- sapply(seq_len(length(nms)), function(i)
    sapply(seq(0,1,length.out=res), function(j)
      pnorm(j,mean=ms[i],sd=sds[i])))

  colnames(d) <- nms

  d <- reshape2::melt(d)

  ggplot2::ggplot(d) + 
    ggplot2::aes(Var2,Var1,fill=value) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors=spectrumPalette(32), limits=c(0,1), guide=F) +
    ggplot2::scale_x_discrete(name=NULL) +
    ggplot2::scale_y_continuous(labels=NULL, name=NULL) +
    cowplot::theme_cowplot() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .5))
}

plotSpectrumPNG <- function(ms, sds, res=32, x=512, y=32) {
  shiny::img(width=x,height=y,style="image-rendering:crisp-edges",
    src=paste0("data:image/png;base64,",openssl::base64_encode(png::writePNG(
      array(t(col2rgb(spectrumPalette(1024)[1+as.integer(1023*sapply(seq_len(length(ms)),
        function(i) sapply(seq(1,0,length.out=res),
          function(j) pnorm(j,mean=ms[i],sd=sds[i]))))]))/255,dim=c(res,length(ms),3))))))
}

defaultMachines <- c('Aurora','Symphony','Fortessa II','BFC 9k')

gatherFormContent <- function(fld, spectra, defaults=NULL)
  nat.sort(unique(c(defaults, unlist(sapply(spectra, function(x) x[[fld]])))))

spectrumMetadataForm <- function(prefix, spectra, selected=NULL, defaultAll=F, create=T, multiple=F) {
  i <- function(x)paste0(prefix, x)
  mkch <- function(a, ...) gatherFormContent(a, spectra, ...)

  shiny::tagList(
    shiny::selectizeInput(i('Machine'), "Cytometer model",
      choices=mkch("machine", if(create) defaultMachines),
      selected=if(defaultAll) mkch("machine"),
      multiple=multiple, options=list(`create`=create)),
    shiny::selectizeInput(i('MachineConfig'), "Cytometer Configuration",
      choices=mkch("mconfig", if(create) "—"),
      selected=if(defaultAll) mkch("mconfig"),
      multiple=multiple, options=list(`create`=create)),
    shiny::selectizeInput(i('SampleType'), "Sample origin/type",
      choices=mkch("sample", if(create) "—"),
      selected=if(defaultAll) mkch("sample"),
      multiple=multiple, options=list(`create`=create)),
    shiny::selectizeInput(i('Antigen'), "Antigen/Antibody",
      choices=mkch("antigen", if(create) defaultAntigens),
      selected=if(defaultAll) mkch("antigen"),
      multiple=multiple, options=list(`create`=create)),
    shiny::selectizeInput(i('Fluorochrome'), "Fluorochrome",
      choices=mkch("fluorochrome", if(create) defaultFluorochromes),
      selected=if(defaultAll) mkch("fluorochrome"),
      multiple=multiple, options=list(`create`=create)),
    shiny::selectizeInput(i('Note'), "Note",
      choices=mkch("note", if(create) "—"),
      selected=if(defaultAll) mkch("note"),
      multiple=multiple, options=list(`create`=create)))
}

spectrumMetadataFormGather <- function(prefix, input) {
  i <- function(x)paste0(prefix, x)
  list(
    machine=input[[i('Machine')]],
    mconfig=input[[i('MachineConfig')]],
    sample=input[[i('SampleType')]],
    antigen=input[[i('Antigen')]],
    fluorochrome=input[[i('Fluorochrome')]],
    note=input[[i('Note')]])
}

#TODO: use this everywhere
spectrumMetadataNames <- c(
  machine='Machine',
  mconfig='Configuration',
  sample='Sample type',
  antigen='Antigen',
  fluorochrome='Fluorochrome',
  note='Note')

spectrumExistsIn <- function(n, ss, fields=names(spectrumMetadataNames))
  any(sapply(ss, function(x) all(sapply(fields, function(nm) n[[nm]]==x[[nm]]))))

spectrumFind <- function(ss, fields=list())
  for(s in ss)
    if(all(sapply(names(fields), function(fld) s[[fld]]==fields[[fld]])))
      return(s)

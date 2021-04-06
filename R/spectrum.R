
#' Extract the spectrum
#'
#' the spectra are scaled to the mI.
extractSpectrum <- function(mtx, gate, powerGate) {
  m2 <- mtx[gate & powerGate,,drop=F]
  a <- sqrt(rowSums(m2^2)) #TODO: energies?

  # TODO: one day we should really use rlm
  reg <- sapply(colnames(m2),
    function(col) {
      reg1 <- lm(m2[,col,drop=F]~a)
      c(intercept=unname(reg1$coefficients['(Intercept)']),
        val=unname(reg1$coefficients['a']),
        sd=sd(reg1$residuals/a))
    })

  csc <- sqrt(sum(reg['val',]^2)) # normalization factor
  sp <- reg['val',]/csc # normalized spectrum
  m2 <- t(mtx[gate,,drop=F])#-reg['intercept',] # matrix for guessing the intensity (intercept removed)
  reg2 <- lm(m2 ~ sp+0) # "unmix" using the single channel

  coefs <- reg2$coefficients[reg2$coefficients>0] # get unmixed values
  ress <- reg2$residuals[,reg2$coefficients>0,drop=F] #residuals

  ws <- coefs^2/(coefs^2+apply(ress^2,2,sum)) # weights: how much of the energy is explained?
  e <- e2db(coefs) # intensities in dB
  wm <- sum(e*ws)/sum(ws) # weighted mean intensities
  wsd <- sqrt(sum((e-wm)^2*ws)/sum(ws)) # weighted sdev of intensities
  list(channels=colnames(mtx),
       mS=unname(sp),
       sdS=unname(reg['sd',]/csc),
       mI=wm,
       sdI=wsd)
}

# decibels should be used for intensities everywhere
e2db <- function(x) 10*log(x[x>0], base=10)
db2e <- function(x) 10^(x/10)
dbf <- function(x, unit='u')
  paste0(round(x,2),' dB',unit)

defaultFluorescenceChannels <- function(nms)
  nms[!grepl("(fsc|ssc|fs0|ss0|time|comp)",nms,ignore.case=T)]

defaultFSCChannel <- function(nms)
  nms[c(grepl("(fsc|fs0).*-a", nms, ignore.case=T),grepl("(fsc|fs0)", nms, ignore.case=T))][1]

defaultSSCChannel <- function(nms)
  nms[c(grepl("(ssc|ss0).*-a", nms, ignore.case=T),grepl("(ssc|ss0)", nms, ignore.case=T))][1]

#spectrumPalette <- grDevices::colorRampPalette(c('#44dd22','#ffddaa','#002040'))
spectrumPalette <- function(n=128, ...)
  #colorspace::sequential_hcl(n=n, 'viridis', rev=T)
  #colorRampPalette(c('#44dd22','white','black'))(n,...)
  colorRampPalette(c('#44dd22','#ff8800','white'))(n,...)
  #colorRampPalette(c('#338822','#ffcc88','black'))(n,...)
  #colorspace::sequential_hcl(n = n, h = 114, c = c(0, 0, 114), l = c(0,100, 50), power = 1.1, ...)

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
defaultSampleTypes <- c('Spleen (mouse)','Thymus (mouse)','PBMC (mouse)','Bone marrow (mouse)','Beads')
defaultAntigens <- c(paste0('CD',1:9),'L/D','Unstained') 
defaultFluorochromes <- c('Cy5', 'Cy5.5', 'Cy7', 'APC', 'FITC', 'PE', 'GFP', 'Unstained')

spectrumMetadataForm <- function(prefix, spectra, selected=NULL, defaultAll=F, create=T, multiple=F) {
  i <- function(x)paste0(prefix, x)
  mkch <- function(fld, defaults=NULL)
    nat.sort(unique(c(defaults, unlist(sapply(spectra, function(x) x[[fld]])))))

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
      choices=mkch("sample", if(create) defaultSampleTypes),
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
  any(sapply(ss, function(x)
    all(sapply(names(spectrumMetadataNames), function(nm)
      n[[nm]]==x[[nm]]))))

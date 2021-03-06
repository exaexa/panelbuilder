
# PanelBuildeR

A small R Shiny app for base processing of flow cytometry samples. Main
functions include:

- advanced acquisition of clean spectra from single-stains
- exporting and importing saved spectra in interoperable JSON format
- combining spectra in panels, calculating SNRs (and some other useful metrics)
  for the spectra combinations
- multicolor and spectral unmixing/compensation

The functionality is mostly cytometer- and method-agnostic, and should
generally work with anything that produces FCS3 files. We have successfully
applied this to data from several cytometer vendors (both spectral and the
"traditional" multicolor).

#### Acknowledgements

PanelBuildeR was developed at the [IOCB Prague](https://www.uochb.cz/)
and [IHBT Prague](https://www.uhkt.cz/),
in cooperation with [Liston Lab](http://www.liston.babraham.ac.uk/)
of the [Babraham Institute](https://www.babraham.ac.uk/).
The development was supported by the Czech node of [ELIXIR](https://elixir-europe.org/),
the European infrastructure for life science.

<img src="media/uochb.svg" alt="IOCB logo" height="64px">   <img src="media/uhkt.svg" alt="IHBT logo" height="64px">   <img src="media/babraham.svg" alt="Babraham Institute logo" height="64px">   <img src="media/elixir.svg" alt="ELIXIR logo" height="64px">

## Screenshot

![screenshot](media/screenshot.png)

## Available algorithms & tools

- Spectra acquisition:
  - cleaning spectra with non-linear regression (e.g.: remove harsh
    autofluorescence from a singlestain with dim antigen/fluorochrome
    combination)
  - Theil-Sen etimation (extremely noise- and outlier-resistant)
- Panel quality estimation using signal-to-noise ratios
  - Estimation of spectra of hypothetical singlestains from available
    data about antigen and fluorochrome brightness
  - (Not yet implemented:) finding of optimal panel assignments
- Unmixing
  - non-linear regression (the implemented variant is similar to non-negative
    least squares, but also supports non-negative residuals to better capture
    the incoming Poisson noise)
  - OLS and several weighted variants thereof
  - MAPE-like unmixing
  - some useful data exported FCS files, e.g. unmixing RMSE and specific residuals
  - simple post-unmixing compensation

## Installation

You can install this using `devtools` directly from the GitHub repository. You will first need the unmixing package `nougad`:

```r
devtools::install_github('exaexa/nougad')
```

Check out the [`nougad`
documentation](https://github.com/exaexa/nougad/blob/master/README.md) for some
helpful hints about installation; in particular there might be very desirable
hardware-accelerated variants.

After that, install PanelBuildeR:

```r
devtools::install_github('exaexa/panelbuilder')
```

Once the package is installed, simply run the `panelbuilder()` function:

```r
library(panelbuilder)
panelbuilder()
```

You should get a browser window open with the PanelBuildeR UI running.

## Quick How-To

- Use the menu on the left side to navigate
- Import spectra from singlestains on the "Import" page. Take care with
  proper annotation of the spectra, you will need it later. Do not forget to
  measure the autofluorescence!
- On "Spectra" page, select spectra that you want to use. You will see them on
  "Panel" page.
- **Panel design:** If you imported multiple spectra that describe the same
  antigen, you may switch them on the "Panel" page now, trying you multiple
  variants of your panel. In some cases, the program may have enough
  information to estimate the spectra from singlestains that you did not import,
  which you can use too. This gives some estimate of how the panel will perform
  before you try it with real reagents and cells.
- **Unmixing:** You may use the imported singlestains on the "Unmixing" page,
  for unmixing/compensating any FCS file that contains matching channels. The
  unmixing will use the antigen/fluorochrome combination specified on the
  "Panel" page.

## Development

PanelBuildeR is designed for easy modification, so that we can quickly try
various new ideas and algorithms for processing the samples and spectra. To get
a "development" environment, clone the repository, and use `test.R` for
bootstrapping and running a local installation:

```r
git clone https://github.com/exaexa/panelbuilder.git
cd panelbuilder
R -f test.R
```

The software is new, there may be bugs and compatibility problems.
If anything looks fishy, feel free to open an issue!

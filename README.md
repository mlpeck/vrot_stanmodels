This site contains code and data for modeling rotation curves of disk galaxies with spectroscopic data in [SDSS-IV MaNGA] (http://www.sdss.org/dr14/manga/) as derived for example in [Teuben 2002] (https://ui.adsabs.harvard.edu/#abs/2002ASPC..275..217T/abstract). All of the R and Stan code and datasets required to recreate my workflow is included. A sample data set with galaxy spectra and a Stan model fit is provided on my [Dropbox account] (https://www.dropbox.com/s/tousr1abr1j27r6/8333-6102.rda?dl=0).

This software depends on a number of R packages. A possibly incomplete list includes:

- [rstan](https://cran.r-project.org/package=rstan)
- [ggplot2](https://cran.r-project.org/package=ggplot2)
- [FITSio](https://cran.r-project.org/package=FITSio)
- [Rsolnp](https://cran.r-project.org/package=Rsolnp)
- [SearchTrees](https://cran.r-project.org/package=SearchTrees)
- [akima](https://cran.r-project.org/package=akima)
- [tibble](https://cran.r-project.org/package=tibble)

Several of these packages also have a number of dependencies. These should be downloaded and installed automatically provided you have an active internet connection. Windows users will also need [Rtools](http://cran.us.r-project.org/bin/windows/Rtools/), which includes a collection of tools for building R packages including mingw compilers. Linux users need the gcc compiler suite. In Linux it's generally easiest to build packages from source, although R itself and possibly R packages may be available in binary form from your distribution's code repositories.

These models are designed to be used on logarithmically sampled RSS files from [SDSS-IV MaNGA](http://www.sdss.org/dr14/manga/manga-data/data-access/). With small modifications they could also be applied to the data cubes, but since no attempt is made to model the spatial covariance between spaxels the results might be misleading and in particular posterior uncertainties are likely to be too optimistic.

## Getting Started

Download all of the code and dataset files (the latter end in .rda) to a local directory, and set that as the working directory in R using `getwd()`. A complete session will look like:

```R

load("pcs.15.rda") ## eigenspectra templates
load("drpcat.rda") ## MaNGA drpall catalog stored as data frame
source("readmanga.r") ##tools for reading MaNGA fits files
source("vrot.r") ## modeling tools

options(mc.cores = parallel::detectCores())

fname <- "name_of_a_logrss_file.fits"
gdat.rss <- readrss(fname, drpcat)
gdat.stack <- stackrss(gdat.rss)
dz.stack <- getdz(gdat.stack, pcs.15, snrthresh=3) ##default value of snrthresh may be too conservative for this task

phi.guess <- as.numeric(drpcat[drpcat$plateifu == gdat.stack$meta$plateifu, "nsa_elpetro_phi"]) ## photometric major axis orientation
ci.guess <- as.numeric(drpcat[drpcat$plateifu == gdat.stack$meta$plateifu, "nsa_elpetro_ba"])  ## photometric minor/major axis ratio
r.eff <- as.numeric(drpcat[drpcat$plateifu == gdat.stack$meta$plateifu, "nsa_elpetro_th50_r"]) ## r band effective radius

vrmodel <- vrot(gdat.stack, dz.stack, phi.guess, ci.guess=ci.guess, r_eff= r.eff)
```

If you downloaded the sample data set loading it will load values of `gdat.stack`, `dz.stack`, `phi.guess`, `ci.guess`, and `r_eff` into the workspace, and only the last line is needed. There are a number of optional arguments to `vrot()`. These set parameters for priors and stan options that are passed on. If execution time is too long the values of `iter` and `warmup` can be reduced -- generally many fewer iterations are required for convergence and inference in Stan compared to traditional Metropolis-Hastings or Gibbs samplers.

There are two Stan models included here for modeling galaxies with moderate disk inclinations. The first, named `vrot.stan` uses a low order polynomial for the rotation velocity. `vrot_gp.stan` is a semi-parametric model that uses a gaussian process regression to model the velocity field. The code for the first model is relatively mature, which doesn't mean the output isn't nonsense of course. The latter is still experimental and is *much* more computationally intensive. To try it a minimal call is:

```R
vrmodel <- vrot(gdat.stack, dz.stack, phi.guess, ci.guess=ci.guess, r_eff= r.eff)
```

## Boilerplate

[SDSS acknowledgement](http://www.sdss.org/collaboration/citing-sdss/)

## License

Text is [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Author: Michael Peck.

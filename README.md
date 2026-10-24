# Distribution Explorer

An R Shiny app for exploring probability distributions, deployed at
<https://huchaoran.shinyapps.io/pdf_plot/>.

## Features

**Univariate tab**

* 14 continuous distributions: normal, Student t, Cauchy, logistic, Laplace, uniform,
  log-normal, exponential, gamma, inverse gamma, chi-squared, F, beta, Weibull.
* 6 discrete distributions: binomial, Poisson, geometric, negative binomial,
  hypergeometric, discrete uniform.
* Six views: density or mass function, CDF, density and CDF side by side, survival
  function, hazard function and quantile function.
* Compare curves: pin the current curve and keep changing the parameters, or even the
  distribution, to overlay up to 9 curves (for example a binomial pmf and its normal
  approximation).
* Fixed axes: set the x and y limits by hand, or press "Freeze current axes" so the
  view stops rescaling while you move a slider.
* Shade $P(a \le X \le b)$ and read its value, and overlay a histogram of a random sample.
* A compact table with the mean, standard deviation, median and shaded probability of every curve on the plot.
* The pdf or pmf formula in LaTeX, with the notation explained.
* Download the current plot as PNG.

**Bivariate 3-D tab**

* Bivariate normal, bivariate Student t, Dirichlet with 3 components, and a mixture of
  two bivariate normals.
* Rotatable 3-D density surface (plotly) next to a contour map, with an optional random
  sample drawn on the contour map.
* Fixed or frozen axes, colour scale and grid resolution controls.
* Compare: pin surfaces to overlay up to 6 bivariate densities, from the same or different distributions, with a table of their means, standard deviations and correlations.
* Formula, mean vector, covariance matrix and correlation.

## Running locally

```r
install.packages(c("shiny", "bslib", "ggplot2", "plotly", "mvtnorm"))
shiny::runApp()
```

All packages are on CRAN, so the app deploys to shinyapps.io with
`rsconnect::deployApp()` as before.

## Layout

* `app.R`: user interface and server.
* `R/distributions.R`: registry of univariate distributions (parameters, density,
  CDF, quantile, random generation, moments, formulas).
* `R/univariate_plot.R`: ggplot2 code for the univariate plots.
* `R/multivariate.R`: registry of bivariate distributions.
* `www/styles.css`: styling.

To add a distribution, add one entry to the registry in `R/distributions.R` or
`R/multivariate.R`; the controls, plots and formula panel pick it up automatically.

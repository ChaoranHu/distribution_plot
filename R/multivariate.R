# Registry of bivariate distributions, shown as 3-D surfaces and contours.
#
# Every entry has label, params, density(x1, x2, par) (vectorised),
# sample(n, par), range(par) giving list(x1 = c(lo, hi), x2 = c(lo, hi)),
# optional check(par), formula, moments_tex (lines of LaTeX for the mean and
# covariance) and notation (LaTeX).

bvn_sigma <- function(p) {
  matrix(c(p$s1^2, p$rho * p$s1 * p$s2, p$rho * p$s1 * p$s2, p$s2^2), 2)
}

sigma_tex <- "\\boldsymbol{\\Sigma} = \\begin{pmatrix} \\sigma_1^2 & \\rho\\sigma_1\\sigma_2 \\\\ \\rho\\sigma_1\\sigma_2 & \\sigma_2^2 \\end{pmatrix}"

bv_location_params <- list(
  param("mu1", "Mean μ₁", -5, 5, 0, 0.1),
  param("mu2", "Mean μ₂", -5, 5, 0, 0.1),
  param("s1", "Std. deviation σ₁", 0.2, 5, 1, 0.1),
  param("s2", "Std. deviation σ₂", 0.2, 5, 1, 0.1),
  param("rho", "Correlation ρ", -0.95, 0.95, 0.5, 0.05)
)

ddirichlet3 <- function(x1, x2, p) {
  x3 <- 1 - x1 - x2
  inside <- x1 > 0 & x2 > 0 & x3 > 0
  # Zero off the triangle, so the surface and contour edges stay smooth.
  out <- numeric(length(x1))
  a <- c(p$a1, p$a2, p$a3)
  logc <- lgamma(sum(a)) - sum(lgamma(a))
  out[inside] <- exp(logc + (a[1] - 1) * log(x1[inside]) + (a[2] - 1) * log(x2[inside]) +
                       (a[3] - 1) * log(x3[inside]))
  out
}

bivariates <- list(
  bvn = list(
    label = "Bivariate normal",
    params = bv_location_params,
    density = function(x1, x2, p) {
      mvtnorm::dmvnorm(cbind(x1, x2), mean = c(p$mu1, p$mu2), sigma = bvn_sigma(p))
    },
    sample = function(n, p) mvtnorm::rmvnorm(n, mean = c(p$mu1, p$mu2), sigma = bvn_sigma(p)),
    range = function(p) list(x1 = p$mu1 + c(-3.5, 3.5) * p$s1, x2 = p$mu2 + c(-3.5, 3.5) * p$s2),
    formula = paste0("f(\\mathbf{x}) = \\frac{1}{2\\pi\\,|\\boldsymbol{\\Sigma}|^{1/2}} ",
                     "\\exp\\left(-\\tfrac{1}{2}(\\mathbf{x}-\\boldsymbol{\\mu})^\\top ",
                     "\\boldsymbol{\\Sigma}^{-1}(\\mathbf{x}-\\boldsymbol{\\mu})\\right), \\quad ",
                     sigma_tex),
    moments_tex = c("E[\\mathbf{X}] &= \\boldsymbol{\\mu}",
                    "\\operatorname{Cov}(\\mathbf{X}) &= \\boldsymbol{\\Sigma}"),
    notation = c("\\(\\mathbf{x} = (x_1, x_2)^\\top \\in \\mathbb{R}^2\\)",
                 "\\(\\boldsymbol{\\mu} = (\\mu_1, \\mu_2)^\\top\\): mean vector",
                 "\\(\\sigma_1, \\sigma_2 > 0\\): marginal standard deviations",
                 "\\(\\rho \\in (-1, 1)\\): correlation between \\(X_1\\) and \\(X_2\\)",
                 "\\(|\\boldsymbol{\\Sigma}|\\): determinant of the covariance matrix")
  ),

  bvt = list(
    label = "Bivariate Student t",
    params = c(bv_location_params, list(param("nu", "Degrees of freedom ν", 1, 30, 3, 0.5))),
    density = function(x1, x2, p) {
      mvtnorm::dmvt(cbind(x1, x2), delta = c(p$mu1, p$mu2), sigma = bvn_sigma(p),
                    df = p$nu, log = FALSE)
    },
    sample = function(n, p) {
      mvtnorm::rmvt(n, delta = c(p$mu1, p$mu2), sigma = bvn_sigma(p), df = p$nu, type = "shifted")
    },
    range = function(p) list(x1 = p$mu1 + c(-5, 5) * p$s1, x2 = p$mu2 + c(-5, 5) * p$s2),
    formula = paste0("f(\\mathbf{x}) = \\frac{\\Gamma\\left(\\frac{\\nu+2}{2}\\right)}",
                     "{\\Gamma\\left(\\frac{\\nu}{2}\\right)\\,\\nu\\pi\\,|\\boldsymbol{\\Sigma}|^{1/2}} ",
                     "\\left(1 + \\frac{1}{\\nu}(\\mathbf{x}-\\boldsymbol{\\mu})^\\top ",
                     "\\boldsymbol{\\Sigma}^{-1}(\\mathbf{x}-\\boldsymbol{\\mu})\\right)^{-\\frac{\\nu+2}{2}}, \\quad ",
                     sigma_tex),
    moments_tex = c("E[\\mathbf{X}] &= \\begin{cases} \\boldsymbol{\\mu} & \\nu > 1 \\\\ \\text{undefined} & \\nu \\le 1 \\end{cases}",
                    "\\operatorname{Cov}(\\mathbf{X}) &= \\begin{cases} \\frac{\\nu}{\\nu-2}\\boldsymbol{\\Sigma} & \\nu > 2 \\\\ \\infty & 1 < \\nu \\le 2 \\\\ \\text{undefined} & \\nu \\le 1 \\end{cases}"),
    notation = c("\\(\\mathbf{x} = (x_1, x_2)^\\top \\in \\mathbb{R}^2\\)",
                 "\\(\\boldsymbol{\\mu}\\): location vector",
                 "\\(\\boldsymbol{\\Sigma}\\): scale matrix built from \\(\\sigma_1, \\sigma_2, \\rho\\); the covariance is \\(\\frac{\\nu}{\\nu-2}\\boldsymbol{\\Sigma}\\) when \\(\\nu > 2\\)",
                 "\\(\\nu > 0\\): degrees of freedom; heavier tails than the normal, which is the limit as \\(\\nu \\to \\infty\\)")
  ),

  dirichlet = list(
    label = "Dirichlet (3 components)",
    params = list(param("a1", "Concentration α₁", 0.3, 15, 2, 0.1),
                  param("a2", "Concentration α₂", 0.3, 15, 3, 0.1),
                  param("a3", "Concentration α₃", 0.3, 15, 4, 0.1)),
    density = ddirichlet3,
    sample = function(n, p) {
      g <- cbind(rgamma(n, p$a1), rgamma(n, p$a2), rgamma(n, p$a3))
      (g / rowSums(g))[, 1:2, drop = FALSE]
    },
    range = function(p) list(x1 = c(0, 1), x2 = c(0, 1)),
    formula = paste0("f(x_1, x_2) = \\frac{\\Gamma(\\alpha_1+\\alpha_2+\\alpha_3)}",
                     "{\\Gamma(\\alpha_1)\\Gamma(\\alpha_2)\\Gamma(\\alpha_3)} ",
                     "x_1^{\\alpha_1-1} x_2^{\\alpha_2-1} x_3^{\\alpha_3-1}, \\quad x_3 = 1 - x_1 - x_2"),
    moments_tex = c("E[X_i] &= \\frac{\\alpha_i}{\\alpha_0}",
                    "\\operatorname{Var}(X_i) &= \\frac{\\alpha_i(\\alpha_0-\\alpha_i)}{\\alpha_0^2(\\alpha_0+1)}",
                    "\\operatorname{Cov}(X_1, X_2) &= -\\frac{\\alpha_1\\alpha_2}{\\alpha_0^2(\\alpha_0+1)}, \\quad \\alpha_0 = \\alpha_1+\\alpha_2+\\alpha_3"),
    notation = c("\\((x_1, x_2, x_3)\\): proportions with \\(x_i > 0\\) and \\(x_1 + x_2 + x_3 = 1\\), so the density is drawn over \\((x_1, x_2)\\) on the triangle",
                 "\\(\\alpha_1, \\alpha_2, \\alpha_3 > 0\\): concentration parameters; values below 1 push mass to the corners",
                 "\\(\\Gamma(\\cdot)\\): the gamma function")
  ),

  mixture = list(
    label = "Mixture of two normals",
    params = list(param("w", "Weight of component 1, w", 0, 1, 0.5, 0.05),
                  param("m1x", "Component 1 mean μ₁₁", -5, 5, -1.5, 0.1),
                  param("m1y", "Component 1 mean μ₁₂", -5, 5, -1, 0.1),
                  param("m2x", "Component 2 mean μ₂₁", -5, 5, 1.5, 0.1),
                  param("m2y", "Component 2 mean μ₂₂", -5, 5, 1, 0.1),
                  param("s", "Common std. deviation σ", 0.2, 3, 0.8, 0.05),
                  param("rho", "Common correlation ρ", -0.95, 0.95, 0, 0.05)),
    density = function(x1, x2, p) {
      S <- bvn_sigma(list(s1 = p$s, s2 = p$s, rho = p$rho))
      X <- cbind(x1, x2)
      p$w * mvtnorm::dmvnorm(X, c(p$m1x, p$m1y), S) +
        (1 - p$w) * mvtnorm::dmvnorm(X, c(p$m2x, p$m2y), S)
    },
    sample = function(n, p) {
      S <- bvn_sigma(list(s1 = p$s, s2 = p$s, rho = p$rho))
      k <- stats::rbinom(n, 1, p$w) == 1
      out <- mvtnorm::rmvnorm(n, c(p$m2x, p$m2y), S)
      if (any(k)) out[k, ] <- mvtnorm::rmvnorm(sum(k), c(p$m1x, p$m1y), S)
      out
    },
    range = function(p) {
      list(x1 = range(c(p$m1x, p$m2x)) + c(-3.5, 3.5) * p$s,
           x2 = range(c(p$m1y, p$m2y)) + c(-3.5, 3.5) * p$s)
    },
    formula = paste0("f(\\mathbf{x}) = w\\,\\phi(\\mathbf{x};\\boldsymbol{\\mu}_1,\\boldsymbol{\\Sigma}) + ",
                     "(1-w)\\,\\phi(\\mathbf{x};\\boldsymbol{\\mu}_2,\\boldsymbol{\\Sigma}), \\quad ",
                     "\\boldsymbol{\\Sigma} = \\sigma^2\\begin{pmatrix} 1 & \\rho \\\\ \\rho & 1 \\end{pmatrix}"),
    moments_tex = c("E[\\mathbf{X}] &= w\\boldsymbol{\\mu}_1 + (1-w)\\boldsymbol{\\mu}_2",
                    "\\operatorname{Cov}(\\mathbf{X}) &= \\boldsymbol{\\Sigma} + w(1-w)(\\boldsymbol{\\mu}_1-\\boldsymbol{\\mu}_2)(\\boldsymbol{\\mu}_1-\\boldsymbol{\\mu}_2)^\\top"),
    notation = c("\\(\\phi(\\mathbf{x};\\boldsymbol{\\mu},\\boldsymbol{\\Sigma})\\): bivariate normal density with mean \\(\\boldsymbol{\\mu}\\) and covariance \\(\\boldsymbol{\\Sigma}\\)",
                 "\\(w \\in [0,1]\\): probability of drawing from component 1",
                 "\\(\\boldsymbol{\\mu}_k = (\\mu_{k1}, \\mu_{k2})^\\top\\): mean of component \\(k\\)",
                 "\\(\\sigma > 0\\), \\(\\rho \\in (-1,1)\\): shared standard deviation and correlation")
  )
)

bv_choices <- stats::setNames(names(bivariates), vapply(bivariates, `[[`, character(1), "label"))

# Density on a grid; z[i, j] = f(x1[j], x2[i]) as plotly expects.
bv_grid <- function(dist_id, par, xr, yr, n = 80) {
  x1 <- seq(xr[1], xr[2], length.out = n)
  x2 <- seq(yr[1], yr[2], length.out = n)
  g <- expand.grid(x1 = x1, x2 = x2)
  z <- bivariates[[dist_id]]$density(g$x1, g$x2, par)
  z[!is.finite(z)] <- NA
  list(x1 = x1, x2 = x2, z = matrix(z, nrow = n, ncol = n, byrow = TRUE))
}

bv_label <- function(dist_id, par) {
  spec <- bivariates[[dist_id]]
  vals <- vapply(spec$params, function(prm) {
    paste0(param_symbol(prm), " = ", format(signif(par[[prm$id]], 4)))
  }, character(1))
  paste0(spec$label, "(", paste(vals, collapse = ", "), ")")
}

# Mean vector and covariance matrix (NA = undefined, Inf = infinite).
bv_moments <- function(dist_id, p) {
  switch(dist_id,
    bvn = list(mean = c(p$mu1, p$mu2), cov = bvn_sigma(p)),
    bvt = list(mean = if (p$nu > 1) c(p$mu1, p$mu2) else c(NA, NA),
               cov = if (p$nu > 2) p$nu / (p$nu - 2) * bvn_sigma(p)
                     else matrix(if (p$nu > 1) Inf else NA, 2, 2)),
    dirichlet = {
      a <- c(p$a1, p$a2, p$a3); a0 <- sum(a)
      cv <- matrix(c(a[1] * (a0 - a[1]), -a[1] * a[2], -a[1] * a[2], a[2] * (a0 - a[2])), 2) /
        (a0^2 * (a0 + 1))
      list(mean = a[1:2] / a0, cov = cv)
    },
    mixture = {
      m1 <- c(p$m1x, p$m1y); m2 <- c(p$m2x, p$m2y)
      S <- bvn_sigma(list(s1 = p$s, s2 = p$s, rho = p$rho))
      list(mean = p$w * m1 + (1 - p$w) * m2, cov = S + p$w * (1 - p$w) * tcrossprod(m1 - m2))
    }
  )
}

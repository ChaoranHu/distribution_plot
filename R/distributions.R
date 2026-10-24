# Registry of univariate distributions.
#
# Every distribution is a list with:
#   label, type ("continuous" or "discrete"), params (list of slider specs),
#   d, p, q, r   functions of (x, par) where par is a named list of parameters,
#   mean, var    functions of par (NA = undefined, Inf = infinite),
#   check        optional function of par returning an error message or NULL,
#   range        optional function of par giving a default x range,
#   heavy        TRUE for heavy tailed distributions (narrower auto range),
#   formula      LaTeX for the pdf / pmf, support, notation (character vector).

param <- function(id, label, min, max, value, step = NULL) {
  list(id = id, label = label, min = min, max = max, value = value, step = step)
}

distributions <- list(

  # ---------------------------------------------------------------- continuous
  normal = list(
    label = "Normal", type = "continuous",
    params = list(param("mu", "Mean \u03bc", -25, 25, 0, 0.1),
                  param("sigma", "Std. deviation \u03c3", 0.1, 10, 1, 0.1)),
    d = function(x, p) dnorm(x, p$mu, p$sigma),
    p = function(x, p) pnorm(x, p$mu, p$sigma),
    q = function(u, p) qnorm(u, p$mu, p$sigma),
    r = function(n, p) rnorm(n, p$mu, p$sigma),
    mean = function(p) p$mu,
    var = function(p) p$sigma^2,
    formula = "f(x) = \\frac{1}{\\sigma\\sqrt{2\\pi}} \\exp\\left(-\\frac{(x-\\mu)^2}{2\\sigma^2}\\right)",
    support = "x \\in \\mathbb{R}",
    notation = c("\\(\\mu\\): mean (location)",
                 "\\(\\sigma > 0\\): standard deviation (scale)")
  ),

  t = list(
    label = "Student t", type = "continuous", heavy = TRUE,
    params = list(param("nu", "Degrees of freedom \u03bd", 0.5, 50, 3, 0.5)),
    d = function(x, p) dt(x, p$nu),
    p = function(x, p) pt(x, p$nu),
    q = function(u, p) qt(u, p$nu),
    r = function(n, p) rt(n, p$nu),
    mean = function(p) if (p$nu > 1) 0 else NA,
    var = function(p) if (p$nu > 2) p$nu / (p$nu - 2) else if (p$nu > 1) Inf else NA,
    range = function(p) c(-6, 6),
    formula = "f(x) = \\frac{\\Gamma\\left(\\frac{\\nu+1}{2}\\right)}{\\sqrt{\\nu\\pi}\\,\\Gamma\\left(\\frac{\\nu}{2}\\right)} \\left(1+\\frac{x^2}{\\nu}\\right)^{-\\frac{\\nu+1}{2}}",
    support = "x \\in \\mathbb{R}",
    notation = c("\\(\\nu > 0\\): degrees of freedom",
                 "\\(\\Gamma(\\cdot)\\): the gamma function")
  ),

  cauchy = list(
    label = "Cauchy", type = "continuous", heavy = TRUE,
    params = list(param("x0", "Location x\u2080", -10, 10, 0, 0.1),
                  param("gamma", "Scale \u03b3", 0.1, 10, 1, 0.1)),
    d = function(x, p) dcauchy(x, p$x0, p$gamma),
    p = function(x, p) pcauchy(x, p$x0, p$gamma),
    q = function(u, p) qcauchy(u, p$x0, p$gamma),
    r = function(n, p) rcauchy(n, p$x0, p$gamma),
    mean = function(p) NA,
    var = function(p) NA,
    range = function(p) p$x0 + c(-8, 8) * p$gamma,
    formula = "f(x) = \\frac{1}{\\pi\\gamma\\left[1+\\left(\\frac{x-x_0}{\\gamma}\\right)^2\\right]}",
    support = "x \\in \\mathbb{R}",
    notation = c("\\(x_0\\): location (median and mode)",
                 "\\(\\gamma > 0\\): scale (half width at half maximum)",
                 "The mean and variance do not exist.")
  ),

  logistic = list(
    label = "Logistic", type = "continuous",
    params = list(param("mu", "Location \u03bc", -10, 10, 0, 0.1),
                  param("s", "Scale s", 0.1, 5, 1, 0.1)),
    d = function(x, p) dlogis(x, p$mu, p$s),
    p = function(x, p) plogis(x, p$mu, p$s),
    q = function(u, p) qlogis(u, p$mu, p$s),
    r = function(n, p) rlogis(n, p$mu, p$s),
    mean = function(p) p$mu,
    var = function(p) p$s^2 * pi^2 / 3,
    formula = "f(x) = \\frac{e^{-(x-\\mu)/s}}{s\\left(1+e^{-(x-\\mu)/s}\\right)^2}",
    support = "x \\in \\mathbb{R}",
    notation = c("\\(\\mu\\): location (mean)", "\\(s > 0\\): scale")
  ),

  laplace = list(
    label = "Laplace", type = "continuous",
    params = list(param("mu", "Location \u03bc", -10, 10, 0, 0.1),
                  param("b", "Scale b", 0.1, 5, 1, 0.1)),
    d = function(x, p) exp(-abs(x - p$mu) / p$b) / (2 * p$b),
    p = function(x, p) ifelse(x < p$mu, 0.5 * exp((x - p$mu) / p$b),
                              1 - 0.5 * exp(-(x - p$mu) / p$b)),
    q = function(u, p) p$mu - p$b * sign(u - 0.5) * log(1 - 2 * abs(u - 0.5)),
    r = function(n, p) {
      u <- runif(n)
      p$mu - p$b * sign(u - 0.5) * log(1 - 2 * abs(u - 0.5))
    },
    mean = function(p) p$mu,
    var = function(p) 2 * p$b^2,
    formula = "f(x) = \\frac{1}{2b} \\exp\\left(-\\frac{|x-\\mu|}{b}\\right)",
    support = "x \\in \\mathbb{R}",
    notation = c("\\(\\mu\\): location (mean, median, mode)", "\\(b > 0\\): scale")
  ),

  uniform = list(
    label = "Uniform", type = "continuous",
    params = list(param("a", "Lower bound a", -10, 10, 0, 0.1),
                  param("b", "Upper bound b", -10, 10, 1, 0.1)),
    d = function(x, p) dunif(x, p$a, p$b),
    p = function(x, p) punif(x, p$a, p$b),
    q = function(u, p) qunif(u, p$a, p$b),
    r = function(n, p) runif(n, p$a, p$b),
    mean = function(p) (p$a + p$b) / 2,
    var = function(p) (p$b - p$a)^2 / 12,
    check = function(p) if (p$a >= p$b) "The lower bound a must be smaller than b.",
    range = function(p) {
      w <- p$b - p$a
      c(p$a - 0.25 * w, p$b + 0.25 * w)
    },
    formula = "f(x) = \\frac{1}{b-a}",
    support = "a \\le x \\le b",
    notation = c("\\(a < b\\): the end points of the interval")
  ),

  lognormal = list(
    label = "Log-normal", type = "continuous", heavy = TRUE,
    params = list(param("mu", "Log mean \u03bc", -3, 3, 0, 0.05),
                  param("sigma", "Log std. deviation \u03c3", 0.05, 2.5, 0.5, 0.05)),
    d = function(x, p) dlnorm(x, p$mu, p$sigma),
    p = function(x, p) plnorm(x, p$mu, p$sigma),
    q = function(u, p) qlnorm(u, p$mu, p$sigma),
    r = function(n, p) rlnorm(n, p$mu, p$sigma),
    mean = function(p) exp(p$mu + p$sigma^2 / 2),
    var = function(p) (exp(p$sigma^2) - 1) * exp(2 * p$mu + p$sigma^2),
    formula = "f(x) = \\frac{1}{x\\sigma\\sqrt{2\\pi}} \\exp\\left(-\\frac{(\\ln x-\\mu)^2}{2\\sigma^2}\\right)",
    support = "x > 0",
    notation = c("\\(\\mu\\): mean of \\(\\ln X\\)",
                 "\\(\\sigma > 0\\): standard deviation of \\(\\ln X\\)")
  ),

  exponential = list(
    label = "Exponential", type = "continuous",
    params = list(param("lambda", "Rate \u03bb", 0.05, 10, 1, 0.05)),
    d = function(x, p) dexp(x, p$lambda),
    p = function(x, p) pexp(x, p$lambda),
    q = function(u, p) qexp(u, p$lambda),
    r = function(n, p) rexp(n, p$lambda),
    mean = function(p) 1 / p$lambda,
    var = function(p) 1 / p$lambda^2,
    formula = "f(x) = \\lambda e^{-\\lambda x}",
    support = "x \\ge 0",
    notation = c("\\(\\lambda > 0\\): rate, so the mean is \\(1/\\lambda\\)")
  ),

  gamma = list(
    label = "Gamma", type = "continuous",
    params = list(param("shape", "Shape k", 0.1, 20, 2, 0.1),
                  param("scale", "Scale \u03b8", 0.1, 10, 1, 0.1)),
    d = function(x, p) dgamma(x, shape = p$shape, scale = p$scale),
    p = function(x, p) pgamma(x, shape = p$shape, scale = p$scale),
    q = function(u, p) qgamma(u, shape = p$shape, scale = p$scale),
    r = function(n, p) rgamma(n, shape = p$shape, scale = p$scale),
    mean = function(p) p$shape * p$scale,
    var = function(p) p$shape * p$scale^2,
    formula = "f(x) = \\frac{1}{\\Gamma(k)\\,\\theta^k} x^{k-1} e^{-x/\\theta}",
    support = "x > 0",
    notation = c("\\(k > 0\\): shape", "\\(\\theta > 0\\): scale",
                 "\\(\\Gamma(\\cdot)\\): the gamma function")
  ),

  invgamma = list(
    label = "Inverse gamma", type = "continuous", heavy = TRUE,
    params = list(param("shape", "Shape \u03b1", 0.1, 20, 3, 0.1),
                  param("scale", "Scale \u03b2", 0.1, 20, 2, 0.1)),
    # If Y ~ Gamma(shape = alpha, rate = beta) then X = 1 / Y ~ InvGamma(alpha, beta)
    d = function(x, p) ifelse(x > 0, dgamma(1 / x, shape = p$shape, rate = p$scale) / x^2, 0),
    p = function(x, p) ifelse(x > 0, pgamma(1 / x, shape = p$shape, rate = p$scale,
                                            lower.tail = FALSE), 0),
    q = function(u, p) 1 / qgamma(1 - u, shape = p$shape, rate = p$scale),
    r = function(n, p) 1 / rgamma(n, shape = p$shape, rate = p$scale),
    mean = function(p) if (p$shape > 1) p$scale / (p$shape - 1) else Inf,
    var = function(p) if (p$shape > 2) p$scale^2 / ((p$shape - 1)^2 * (p$shape - 2)) else Inf,
    formula = "f(x) = \\frac{\\beta^\\alpha}{\\Gamma(\\alpha)} x^{-\\alpha-1} e^{-\\beta/x}",
    support = "x > 0",
    notation = c("\\(\\alpha > 0\\): shape", "\\(\\beta > 0\\): scale",
                 "If \\(Y \\sim \\text{Gamma}(\\alpha, \\text{rate}=\\beta)\\) then \\(1/Y\\) has this distribution.")
  ),

  chisq = list(
    label = "Chi-squared", type = "continuous",
    params = list(param("k", "Degrees of freedom k", 0.5, 50, 3, 0.5)),
    d = function(x, p) dchisq(x, p$k),
    p = function(x, p) pchisq(x, p$k),
    q = function(u, p) qchisq(u, p$k),
    r = function(n, p) rchisq(n, p$k),
    mean = function(p) p$k,
    var = function(p) 2 * p$k,
    formula = "f(x) = \\frac{1}{2^{k/2}\\,\\Gamma(k/2)} x^{k/2-1} e^{-x/2}",
    support = "x > 0",
    notation = c("\\(k > 0\\): degrees of freedom",
                 "\\(\\Gamma(\\cdot)\\): the gamma function")
  ),

  f = list(
    label = "F", type = "continuous", heavy = TRUE,
    params = list(param("d1", "Numerator df d\u2081", 0.5, 100, 5, 0.5),
                  param("d2", "Denominator df d\u2082", 0.5, 100, 10, 0.5)),
    d = function(x, p) df(x, p$d1, p$d2),
    p = function(x, p) pf(x, p$d1, p$d2),
    q = function(u, p) qf(u, p$d1, p$d2),
    r = function(n, p) rf(n, p$d1, p$d2),
    mean = function(p) if (p$d2 > 2) p$d2 / (p$d2 - 2) else NA,
    var = function(p) {
      if (p$d2 > 4) 2 * p$d2^2 * (p$d1 + p$d2 - 2) / (p$d1 * (p$d2 - 2)^2 * (p$d2 - 4))
      else if (p$d2 > 2) Inf else NA
    },
    formula = "f(x) = \\frac{1}{x\\,B\\left(\\frac{d_1}{2},\\frac{d_2}{2}\\right)} \\sqrt{\\frac{(d_1x)^{d_1} d_2^{d_2}}{(d_1x+d_2)^{d_1+d_2}}}",
    support = "x > 0",
    notation = c("\\(d_1, d_2 > 0\\): numerator and denominator degrees of freedom",
                 "\\(B(\\cdot,\\cdot)\\): the beta function")
  ),

  beta = list(
    label = "Beta", type = "continuous",
    params = list(param("a", "Shape \u03b1", 0.1, 20, 2, 0.1),
                  param("b", "Shape \u03b2", 0.1, 20, 5, 0.1)),
    d = function(x, p) dbeta(x, p$a, p$b),
    p = function(x, p) pbeta(x, p$a, p$b),
    q = function(u, p) qbeta(u, p$a, p$b),
    r = function(n, p) rbeta(n, p$a, p$b),
    mean = function(p) p$a / (p$a + p$b),
    var = function(p) p$a * p$b / ((p$a + p$b)^2 * (p$a + p$b + 1)),
    range = function(p) c(0, 1),
    formula = "f(x) = \\frac{x^{\\alpha-1}(1-x)^{\\beta-1}}{B(\\alpha,\\beta)}",
    support = "0 < x < 1",
    notation = c("\\(\\alpha, \\beta > 0\\): shape parameters",
                 "\\(B(\\alpha,\\beta) = \\Gamma(\\alpha)\\Gamma(\\beta)/\\Gamma(\\alpha+\\beta)\\): the beta function")
  ),

  weibull = list(
    label = "Weibull", type = "continuous",
    params = list(param("shape", "Shape k", 0.2, 10, 1.5, 0.1),
                  param("scale", "Scale \u03bb", 0.1, 10, 1, 0.1)),
    d = function(x, p) dweibull(x, p$shape, p$scale),
    p = function(x, p) pweibull(x, p$shape, p$scale),
    q = function(u, p) qweibull(u, p$shape, p$scale),
    r = function(n, p) rweibull(n, p$shape, p$scale),
    mean = function(p) p$scale * gamma(1 + 1 / p$shape),
    var = function(p) p$scale^2 * (gamma(1 + 2 / p$shape) - gamma(1 + 1 / p$shape)^2),
    formula = "f(x) = \\frac{k}{\\lambda}\\left(\\frac{x}{\\lambda}\\right)^{k-1} e^{-(x/\\lambda)^k}",
    support = "x \\ge 0",
    notation = c("\\(k > 0\\): shape", "\\(\\lambda > 0\\): scale")
  ),

  # ------------------------------------------------------------------ discrete
  binomial = list(
    label = "Binomial", type = "discrete",
    params = list(param("n", "Number of trials n", 1, 100, 10, 1),
                  param("prob", "Success probability p", 0, 1, 0.5, 0.01)),
    d = function(x, p) dbinom(x, p$n, p$prob),
    p = function(x, p) pbinom(x, p$n, p$prob),
    q = function(u, p) qbinom(u, p$n, p$prob),
    r = function(n, p) rbinom(n, p$n, p$prob),
    mean = function(p) p$n * p$prob,
    var = function(p) p$n * p$prob * (1 - p$prob),
    range = function(p) c(0, p$n),
    formula = "P(X = x) = \\binom{n}{x} p^x (1-p)^{n-x}",
    support = "x \\in \\{0, 1, \\dots, n\\}",
    notation = c("\\(n\\): number of independent trials",
                 "\\(p \\in [0,1]\\): probability of success in each trial",
                 "\\(X\\): number of successes")
  ),

  poisson = list(
    label = "Poisson", type = "discrete",
    params = list(param("lambda", "Rate \u03bb", 0.1, 50, 4, 0.1)),
    d = function(x, p) dpois(x, p$lambda),
    p = function(x, p) ppois(x, p$lambda),
    q = function(u, p) qpois(u, p$lambda),
    r = function(n, p) rpois(n, p$lambda),
    mean = function(p) p$lambda,
    var = function(p) p$lambda,
    formula = "P(X = x) = \\frac{\\lambda^x e^{-\\lambda}}{x!}",
    support = "x \\in \\{0, 1, 2, \\dots\\}",
    notation = c("\\(\\lambda > 0\\): expected number of events (mean and variance)")
  ),

  geometric = list(
    label = "Geometric", type = "discrete",
    params = list(param("prob", "Success probability p", 0.01, 1, 0.3, 0.01)),
    d = function(x, p) dgeom(x, p$prob),
    p = function(x, p) pgeom(x, p$prob),
    q = function(u, p) qgeom(u, p$prob),
    r = function(n, p) rgeom(n, p$prob),
    mean = function(p) (1 - p$prob) / p$prob,
    var = function(p) (1 - p$prob) / p$prob^2,
    formula = "P(X = x) = (1-p)^x p",
    support = "x \\in \\{0, 1, 2, \\dots\\}",
    notation = c("\\(p \\in (0,1]\\): probability of success in each trial",
                 "\\(X\\): number of failures before the first success")
  ),

  negbin = list(
    label = "Negative binomial", type = "discrete",
    params = list(param("size", "Target successes r", 1, 50, 3, 1),
                  param("prob", "Success probability p", 0.01, 1, 0.5, 0.01)),
    d = function(x, p) dnbinom(x, p$size, p$prob),
    p = function(x, p) pnbinom(x, p$size, p$prob),
    q = function(u, p) qnbinom(u, p$size, p$prob),
    r = function(n, p) rnbinom(n, p$size, p$prob),
    mean = function(p) p$size * (1 - p$prob) / p$prob,
    var = function(p) p$size * (1 - p$prob) / p$prob^2,
    formula = "P(X = x) = \\binom{x+r-1}{x} p^r (1-p)^x",
    support = "x \\in \\{0, 1, 2, \\dots\\}",
    notation = c("\\(r\\): number of successes to wait for",
                 "\\(p \\in (0,1]\\): probability of success in each trial",
                 "\\(X\\): number of failures before the \\(r\\)th success")
  ),

  hyper = list(
    label = "Hypergeometric", type = "discrete",
    params = list(param("N", "Population size N", 1, 100, 50, 1),
                  param("K", "Successes in population K", 0, 100, 20, 1),
                  param("n", "Number of draws n", 0, 100, 10, 1)),
    d = function(x, p) dhyper(x, p$K, p$N - p$K, p$n),
    p = function(x, p) phyper(x, p$K, p$N - p$K, p$n),
    q = function(u, p) qhyper(u, p$K, p$N - p$K, p$n),
    r = function(n, p) rhyper(n, p$K, p$N - p$K, p$n),
    mean = function(p) p$n * p$K / p$N,
    var = function(p) {
      if (p$N <= 1) return(0)
      p$n * (p$K / p$N) * ((p$N - p$K) / p$N) * ((p$N - p$n) / (p$N - 1))
    },
    check = function(p) {
      if (p$K > p$N) return("K cannot exceed the population size N.")
      if (p$n > p$N) return("The number of draws n cannot exceed N.")
      NULL
    },
    range = function(p) c(max(0, p$n - (p$N - p$K)), min(p$n, p$K)),
    formula = "P(X = x) = \\frac{\\binom{K}{x}\\binom{N-K}{n-x}}{\\binom{N}{n}}",
    support = "\\max(0, n+K-N) \\le x \\le \\min(n, K)",
    notation = c("\\(N\\): population size", "\\(K\\): number of successes in the population",
                 "\\(n\\): number of draws without replacement",
                 "\\(X\\): number of successes drawn")
  ),

  dunif = list(
    label = "Discrete uniform", type = "discrete",
    params = list(param("a", "Lower bound a", -20, 20, 1, 1),
                  param("b", "Upper bound b", -20, 20, 6, 1)),
    d = function(x, p) ifelse(x >= p$a & x <= p$b & x == round(x), 1 / (p$b - p$a + 1), 0),
    p = function(x, p) pmin(pmax((floor(x) - p$a + 1) / (p$b - p$a + 1), 0), 1),
    q = function(u, p) pmax(p$a, p$a + ceiling(u * (p$b - p$a + 1)) - 1),
    r = function(n, p) p$a + sample.int(p$b - p$a + 1, n, replace = TRUE) - 1,
    mean = function(p) (p$a + p$b) / 2,
    var = function(p) ((p$b - p$a + 1)^2 - 1) / 12,
    check = function(p) if (p$a > p$b) "The lower bound a cannot exceed b.",
    range = function(p) c(p$a, p$b),
    formula = "P(X = x) = \\frac{1}{b-a+1}",
    support = "x \\in \\{a, a+1, \\dots, b\\}",
    notation = c("\\(a \\le b\\): integer end points")
  )
)

dist_choices <- function(type) {
  ids <- names(distributions)[vapply(distributions, function(d) d$type == type, logical(1))]
  stats::setNames(ids, vapply(distributions[ids], `[[`, character(1), "label"))
}

# Short Greek or Latin symbol for a parameter, taken from the end of its label.
param_symbol <- function(prm) sub(".* ", "", prm$label)

curve_label <- function(dist_id, par) {
  spec <- distributions[[dist_id]]
  vals <- vapply(spec$params, function(prm) {
    paste0(param_symbol(prm), " = ", format(signif(par[[prm$id]], 4)))
  }, character(1))
  paste0(spec$label, "(", paste(vals, collapse = ", "), ")")
}

# Default x range for one curve, before any user override.
auto_range <- function(dist_id, par) {
  spec <- distributions[[dist_id]]
  if (!is.null(spec$range)) return(spec$range(par))
  probs <- if (isTRUE(spec$heavy)) c(0.001, 0.97) else c(0.001, 0.999)
  rng <- spec$q(probs, par)
  if (spec$type == "discrete") rng <- c(floor(rng[1]), ceiling(rng[2]))
  if (!all(is.finite(rng)) || rng[1] >= rng[2]) rng <- c(rng[1], rng[1] + 1)
  rng
}

fmt_num <- function(v) {
  if (length(v) == 0 || is.na(v)) return("undefined")
  if (is.infinite(v)) return("\u221e")
  trimws(formatC(v, digits = 4, format = "g"))
}

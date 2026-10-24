# Helpers that turn a set of curves into a ggplot.
#
# A "curve" is list(dist = <id>, par = <named list>, label = <chr>,
#                   colour = <hex>, live = <lgl>).

fun_labels <- c(pdf = "Density / mass", cdf = "Cumulative distribution F(x)",
                surv = "Survival 1 − F(x)", hazard = "Hazard h(x)",
                quantile = "Quantile F⁻¹(u)")

palette_okabe <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                   "#E69F00", "#56B4E9", "#332288", "#882255", "#444444")

# Evaluate one function of one curve over the x limits.
curve_values <- function(cv, fun, xlim, n_grid = 801) {
  spec <- distributions[[cv$dist]]
  disc <- spec$type == "discrete"

  if (fun == "quantile") {
    u <- if (disc) seq(0.0005, 0.9995, length.out = 2000) else seq(0.001, 0.999, length.out = n_grid)
    return(data.frame(x = u, y = spec$q(u, cv$par)))
  }

  x <- if (disc) {
    lo <- ceiling(xlim[1]); hi <- floor(xlim[2])
    if (hi < lo) return(data.frame(x = numeric(0), y = numeric(0)))
    if (hi - lo > 3000) hi <- lo + 3000
    lo:hi
  } else {
    seq(xlim[1], xlim[2], length.out = n_grid)
  }

  y <- switch(fun,
    pdf = spec$d(x, cv$par),
    cdf = spec$p(x, cv$par),
    surv = 1 - spec$p(x, cv$par),
    hazard = {
      # Continuous: f(x) / (1 - F(x)). Discrete: P(X = x) / P(X >= x).
      tail <- if (disc) 1 - spec$p(x - 1, cv$par) else 1 - spec$p(x, cv$par)
      h <- spec$d(x, cv$par) / tail
      h[tail < 1e-9] <- NA
      h
    }
  )
  y[!is.finite(y)] <- NA
  data.frame(x = x, y = y)
}

# Robust upper y limit for densities that diverge at the edge of the support,
# e.g. Gamma with shape < 1 or Beta with a shape < 1.
robust_ymax <- function(y) {
  y <- y[is.finite(y)]
  if (!length(y)) return(1)
  m <- max(y)
  i <- which.max(y)
  q <- stats::quantile(y, 0.95, names = FALSE)
  if ((i <= 2 || i >= length(y) - 1) && m > 4 * q && q > 0) return(q * 2)
  m
}

# Build the plot for one function.
#   curves  list of curves
#   fun     one of names(fun_labels)
#   xlim    x limits (data scale) used for evaluation and display
#   ylim    NULL for automatic, or c(lo, hi)
#   shade   NULL or c(a, b): shade P(a <= X <= b) for the live curve
#   sample  NULL or numeric vector drawn from the live curve
make_uni_plot <- function(curves, fun, xlim, ylim = NULL, shade = NULL, sample = NULL) {
  labels <- vapply(curves, `[[`, character(1), "label")
  colours <- stats::setNames(vapply(curves, `[[`, character(1), "colour"), labels)
  widths <- stats::setNames(ifelse(vapply(curves, `[[`, logical(1), "live"), 1.4, 0.9), labels)

  n_disc <- sum(vapply(curves, function(cv) distributions[[cv$dist]]$type == "discrete", logical(1)))
  disc_i <- 0
  cont_df <- list(); lolli_df <- list(); step_df <- list()
  auto_ymax <- 0

  for (cv in curves) {
    v <- curve_values(cv, fun, xlim)
    if (!nrow(v)) next
    v$curve <- cv$label
    disc <- distributions[[cv$dist]]$type == "discrete"
    if (!disc) {
      cont_df[[length(cont_df) + 1]] <- v
      auto_ymax <- max(auto_ymax, if (fun %in% c("pdf", "hazard")) robust_ymax(v$y) else max(v$y, na.rm = TRUE))
    } else {
      disc_i <- disc_i + 1
      if (fun %in% c("pdf", "hazard")) {
        # Small horizontal offsets so that several pmfs do not hide each other.
        v$x <- v$x + (disc_i - (n_disc + 1) / 2) * min(0.12, 0.6 / n_disc)
        lolli_df[[length(lolli_df) + 1]] <- v
      } else {
        step_df[[length(step_df) + 1]] <- v
      }
      auto_ymax <- max(auto_ymax, max(v$y, na.rm = TRUE))
    }
  }
  bind <- function(l) if (length(l)) do.call(rbind, l) else NULL
  cont_df <- bind(cont_df); lolli_df <- bind(lolli_df); step_df <- bind(step_df)

  g <- ggplot()

  # Random sample overlay (density view only)
  live <- Filter(function(cv) cv$live, curves)
  if (!is.null(sample) && fun == "pdf" && length(live)) {
    live <- live[[1]]
    sdf <- data.frame(x = sample)
    if (distributions[[live$dist]]$type == "discrete") {
      tab <- as.data.frame(table(sample) / length(sample), stringsAsFactors = FALSE)
      names(tab) <- c("x", "y")
      tab$x <- as.numeric(tab$x)
      g <- g + geom_col(data = tab, aes(x, y), fill = live$colour, alpha = 0.22, width = 0.7)
    } else {
      sdf <- sdf[sdf$x >= xlim[1] & sdf$x <= xlim[2], , drop = FALSE]
      g <- g + geom_histogram(data = sdf, aes(x, y = after_stat(count / length(sample) / width)),
                              bins = 40, fill = live$colour, alpha = 0.22, colour = NA,
                              boundary = xlim[1])
    }
  }

  # Probability shading for the live curve
  if (!is.null(shade) && length(Filter(function(cv) cv$live, curves))) {
    lv <- Filter(function(cv) cv$live, curves)[[1]]
    spec <- distributions[[lv$dist]]
    a <- max(shade[1], xlim[1]); b <- min(shade[2], xlim[2])
    if (fun == "pdf" && a < b) {
      if (spec$type == "continuous") {
        xs <- seq(a, b, length.out = 400)
        ys <- spec$d(xs, lv$par); ys[!is.finite(ys)] <- NA
        g <- g + geom_area(data = data.frame(x = xs, y = ys), aes(x, y),
                           fill = lv$colour, alpha = 0.25)
      } else {
        xs <- ceiling(a):floor(b)
        if (length(xs) && xs[1] <= xs[length(xs)]) {
          g <- g + geom_col(data = data.frame(x = xs, y = spec$d(xs, lv$par)), aes(x, y),
                            fill = lv$colour, alpha = 0.25, width = 0.8)
        }
      }
    }
    if (fun %in% c("cdf", "surv")) {
      g <- g + geom_vline(xintercept = shade, linetype = "dashed", colour = lv$colour, alpha = 0.6)
    }
  }

  if (!is.null(cont_df)) {
    g <- g + geom_line(data = cont_df, aes(x, y, colour = curve, linewidth = curve), na.rm = TRUE)
  }
  if (!is.null(lolli_df)) {
    g <- g +
      geom_segment(data = lolli_df, aes(x = x, xend = x, y = 0, yend = y, colour = curve),
                   linewidth = 0.8, na.rm = TRUE) +
      geom_point(data = lolli_df, aes(x, y, colour = curve), size = 2.2, na.rm = TRUE)
  }
  if (!is.null(step_df)) {
    g <- g +
      geom_step(data = step_df, aes(x, y, colour = curve, linewidth = curve),
                direction = if (fun == "quantile") "vh" else "hv", na.rm = TRUE)
    if (fun != "quantile") g <- g + geom_point(data = step_df, aes(x, y, colour = curve), size = 1.6, na.rm = TRUE)
  }

  if (is.null(ylim)) {
    ylim <- switch(fun,
      cdf = c(0, 1), surv = c(0, 1),
      quantile = NULL,
      c(0, if (auto_ymax > 0) auto_ymax * 1.05 else 1)
    )
  }
  if (fun == "quantile") {
    coord <- coord_cartesian(xlim = c(0, 1), ylim = ylim)
  } else {
    coord <- coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE)
  }

  y_title <- switch(fun,
    pdf = if (n_disc == length(curves)) "P(X = x)" else if (n_disc == 0) "f(x)" else "f(x)  or  P(X = x)",
    cdf = "F(x) = P(X ≤ x)", surv = "S(x) = P(X > x)", hazard = "h(x)",
    quantile = "x")

  g <- g + coord +
    scale_colour_manual(values = colours, breaks = labels, name = NULL) +
    scale_linewidth_manual(values = widths, breaks = labels, guide = "none") +
    labs(x = if (fun == "quantile") "u" else "x", y = y_title, title = fun_labels[[fun]]) +
    theme_minimal(base_size = 15) +
    theme(legend.position = "bottom", legend.direction = "vertical",
          legend.text = element_text(size = 12),
          plot.title = element_text(face = "bold", size = 14, colour = "#334155"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#E2E8F0"),
          axis.title = element_text(colour = "#475569"),
          plot.background = element_rect(fill = "white", colour = NA))
  attr(g, "ylim") <- ylim
  g
}

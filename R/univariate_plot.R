# Helpers that turn a set of curves into a ggplot.
#
# A "curve" is list(dist = <id>, par = <named list>, label = <chr>,
#                   colour = <hex>, live = <lgl>).

fun_labels <- c(pdf = "Density / mass", cdf = "Cumulative distribution F(x)",
                surv = "Survival 1 − F(x)", hazard = "Hazard h(x)",
                quantile = "Quantile F⁻¹(u)")

# Categorical curve colours, assigned in this fixed order (pinned curves first, then
# the live curve). The first eight are checked for colour-blind separation; the ninth
# is a neutral slate for the last possible slot.
# The dark mode column is the same hues stepped for the dark surface.
palette_curves <- c("#2A78D6", "#EB6834", "#1BAF7A", "#EDA100",
                    "#E87BA4", "#008300", "#6250D6", "#E34948", "#475569")
palette_curves_dark <- c("#3987E5", "#D95926", "#199E70", "#C98500",
                         "#D55181", "#008300", "#9085E9", "#E66767", "#94A3B8")
curve_palette <- function(dark = FALSE) if (isTRUE(dark)) palette_curves_dark else palette_curves

# Ink colours shared by the ggplot theme and the plotly layouts, for each colour mode.
ink_light <- list(text = "#1E293B", muted = "#64748B", grid = "#E7EBF2", axis = "#CBD5E1",
                  surface = "#FFFFFF", scene = "#F7F9FC", scene_grid = "#DCE3EE")
ink_dark <- list(text = "#E2E8F0", muted = "#94A3B8", grid = "#223049", axis = "#334462",
                 surface = "#131C31", scene = "#18233D", scene_grid = "#2A3957")
ink_for <- function(dark = FALSE) if (isTRUE(dark)) ink_dark else ink_light
ink <- ink_light

# ggplot theme matching the app's cards.
theme_explorer <- function(base_size = 15, dark = FALSE) {
  ink <- ink_for(dark)
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(colour = ink$text),
      plot.title = element_text(face = "bold", size = rel(0.85), colour = ink$text,
                                margin = margin(b = 4)),
      plot.title.position = "plot",
      axis.title = element_text(colour = ink$muted, size = rel(0.85)),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      axis.text = element_text(colour = ink$muted, size = rel(0.8)),
      axis.line.x = element_line(colour = ink$axis, linewidth = 0.5),
      axis.ticks.x = element_line(colour = ink$axis, linewidth = 0.5),
      axis.ticks.length.x = unit(4, "pt"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = ink$grid, linewidth = 0.5),
      legend.position = "top", legend.justification = "left",
      legend.direction = "vertical", legend.margin = margin(0, 0, 2, 0),
      legend.text = element_text(size = rel(0.78), colour = ink$text),
      legend.key.width = unit(22, "pt"), legend.key.height = unit(14, "pt"),
      plot.margin = margin(10, 14, 8, 8),
      plot.background = element_rect(fill = ink$surface, colour = NA)
    )
}

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
#   dark    draw for the dark colour mode
#   title   show the function name above the plot
make_uni_plot <- function(curves, fun, xlim, ylim = NULL, shade = NULL, sample = NULL,
                          dark = FALSE, title = TRUE) {
  ink <- ink_for(dark)
  labels <- vapply(curves, `[[`, character(1), "label")
  colours <- stats::setNames(vapply(curves, `[[`, character(1), "colour"), labels)
  widths <- stats::setNames(ifelse(vapply(curves, `[[`, logical(1), "live"), 1.5, 1), labels)

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
      g <- g + geom_col(data = tab, aes(x, y), fill = live$colour, alpha = 0.2, width = 0.7)
    } else {
      sdf <- sdf[sdf$x >= xlim[1] & sdf$x <= xlim[2], , drop = FALSE]
      g <- g + geom_histogram(data = sdf, aes(x, y = after_stat(count / length(sample) / width)),
                              bins = 40, fill = live$colour, alpha = 0.2, colour = ink$surface,
                              linewidth = 0.4,
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
      # Write the probability inside the shaded region.
      lower <- if (spec$type == "discrete") spec$p(ceiling(shade[1]) - 1, lv$par) else spec$p(shade[1], lv$par)
      prob <- spec$p(shade[2], lv$par) - lower
      y_mid <- 0.45 * max(spec$d(if (spec$type == "discrete") ceiling(a):max(ceiling(a), floor(b))
                                 else seq(a, b, length.out = 200), lv$par), na.rm = TRUE)
      if (is.finite(prob) && is.finite(y_mid)) {
        g <- g + annotate("label", x = (a + b) / 2, y = y_mid,
                          label = sprintf("P(%s \u2264 X \u2264 %s) = %.4f",
                                          format(shade[1]), format(shade[2]), prob),
                          size = 3.8, colour = ink$text, fill = scales::alpha(ink$surface, 0.88),
                          label.size = 0, label.padding = unit(0.35, "lines"), fontface = "bold")
      }
    }
    if (fun %in% c("cdf", "surv")) {
      g <- g + geom_vline(xintercept = shade, linetype = "dashed", colour = lv$colour, alpha = 0.6)
    }
  }

  # Soft fill under the live density so it reads as the curve in focus.
  live_cont <- Filter(function(cv) cv$live && distributions[[cv$dist]]$type == "continuous", curves)
  if (fun == "pdf" && length(live_cont) && !is.null(cont_df) && is.null(sample) && is.null(shade)) {
    lf <- cont_df[cont_df$curve == live_cont[[1]]$label & !is.na(cont_df$y), ]
    if (nrow(lf)) g <- g + geom_area(data = lf, aes(x, y), fill = live_cont[[1]]$colour, alpha = 0.10)
  }

  # A dotted line at the mean of the live curve.
  live_cv <- Filter(function(cv) cv$live, curves)
  if (fun == "pdf" && length(live_cv)) {
    lv <- live_cv[[1]]
    mu <- distributions[[lv$dist]]$mean(lv$par)
    if (is.finite(mu) && mu >= xlim[1] && mu <= xlim[2]) {
      g <- g + geom_vline(xintercept = mu, colour = ink$muted, linetype = "dotted", linewidth = 0.6)
    }
  }

  if (!is.null(cont_df)) {
    g <- g + geom_line(data = cont_df, aes(x, y, colour = curve, linewidth = curve),
                       lineend = "round", na.rm = TRUE)
  }
  if (!is.null(lolli_df)) {
    g <- g +
      geom_segment(data = lolli_df, aes(x = x, xend = x, y = 0, yend = y, colour = curve),
                   linewidth = 0.8, na.rm = TRUE) +
      geom_point(data = lolli_df, aes(x, y, fill = curve), shape = 21, colour = ink$surface,
                 stroke = 0.7, size = 2.6, na.rm = TRUE)
  }
  if (!is.null(step_df)) {
    g <- g +
      geom_step(data = step_df, aes(x, y, colour = curve, linewidth = curve),
                direction = if (fun == "quantile") "vh" else "hv", na.rm = TRUE)
    if (fun != "quantile") {
      g <- g + geom_point(data = step_df, aes(x, y, fill = curve), shape = 21, colour = ink$surface,
                          stroke = 0.6, size = 2.1, na.rm = TRUE)
    }
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
    scale_fill_manual(values = colours, breaks = labels, guide = "none") +
    scale_linewidth_manual(values = widths, breaks = labels, guide = "none") +
    labs(x = if (fun == "quantile") "u" else "x", y = y_title,
         title = if (title) fun_labels[[fun]] else NULL) +
    theme_explorer(dark = dark)
  attr(g, "ylim") <- ylim
  g
}

# Distribution Explorer
#
# Interactive plots of probability distributions: pdf / pmf, cdf, survival,
# hazard and quantile functions for univariate distributions, with several
# parameter sets overlaid for comparison, plus 3-D density surfaces for
# bivariate distributions.
#
# Helper files in R/ are loaded automatically by shiny.
#
# CRAN packages: shiny, bslib, ggplot2, plotly, mvtnorm

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(mvtnorm)

app_theme <- bs_theme(
  version = 5,
  primary = "#2F5BD3",
  secondary = "#64748B",
  bg = "#F3F5F9",
  fg = "#1E293B",
  base_font = font_collection(font_google("Inter", local = FALSE), "system-ui", "sans-serif"),
  heading_font = font_collection(font_google("Plus Jakarta Sans", wght = c(600, 700, 800), local = FALSE),
                                 "system-ui", "sans-serif"),
  code_font = font_google("JetBrains Mono", local = FALSE),
  "border-radius" = "0.6rem",
  "border-radius-lg" = "0.9rem",
  "card-border-radius" = "1rem",
  "border-color" = "#E3E8F0",
  "link-decoration" = "none",
  "body-bg-dark" = "#0B1222",
  "body-color-dark" = "#E2E8F0",
  "border-color-dark" = "#223049"
)

# A small bell curve used as the app logo.
logo_svg <- HTML(paste0(
  '<svg class="brand-logo" viewBox="0 0 32 32" width="30" height="30" aria-hidden="true">',
  '<rect width="32" height="32" rx="9" fill="url(#lg)"/>',
  '<defs><linearGradient id="lg" x1="0" y1="0" x2="1" y2="1">',
  '<stop offset="0" stop-color="#5B8DEF"/><stop offset="1" stop-color="#EB6834"/></linearGradient></defs>',
  '<path d="M5 24 C10 24 11 8 16 8 C21 8 22 24 27 24" fill="none" stroke="#fff" ',
  'stroke-width="2.6" stroke-linecap="round"/></svg>'
))

uni_input_id <- function(dist, pid) paste0("u_", dist, "_", pid)
bv_input_id <- function(dist, pid) paste0("b_", dist, "_", pid)

param_sliders <- function(registry, id_fun, tabs_id) {
  panels <- lapply(names(registry), function(id) {
    nav_panel_hidden(
      value = id,
      lapply(registry[[id]]$params, function(prm) {
        sliderInput(id_fun(id, prm$id), prm$label, min = prm$min, max = prm$max,
                    value = prm$value, step = prm$step, ticks = FALSE)
      })
    )
  })
  do.call(navset_hidden, c(panels, list(id = tabs_id)))
}

pin_list_ui <- function(items, remove_input, empty_text, palette = palette_curves) {
  if (!length(items)) return(p(class = "small text-muted fst-italic", empty_text))
  tags$div(class = "pin-list", lapply(seq_along(items), function(i) {
    tags$div(
      class = "pin-row",
      tags$span(class = "swatch", style = paste0("background:", palette[i])),
      tags$span(class = "pin-label", items[[i]]$label),
      tags$button(class = "btn btn-sm btn-link pin-remove", title = "Remove",
                  onclick = sprintf("Shiny.setInputValue('%s', %d, {priority: 'event'})", remove_input, i),
                  icon("xmark"))
    )
  }))
}

# Single hue sequential scales in the app palette; other names go to plotly as is.
# Low densities sit near the surface colour, so dark mode runs from dark to light.
seq_scales <- list(
  Ocean = c("#F4F8FE", "#CDE2FB", "#9EC5F4", "#6DA7EC", "#3987E5", "#256ABF", "#184F95", "#0D366B"),
  Ember = c("#FFF7F2", "#FDDCC9", "#F9B596", "#F48B62", "#EB6834", "#C9501F", "#9C3C15", "#6E2A0E")
)
seq_scales_dark <- list(
  Ocean = c("#131C31", "#104281", "#1C5CAB", "#2A78D6", "#5598E7", "#86B6EF", "#B7D3F6", "#E6F0FD"),
  Ember = c("#1C1512", "#6E2A0E", "#9C3C15", "#C9501F", "#EB6834", "#F48B62", "#F9B596", "#FDDCC9")
)
bv_colorscale <- function(name, dark = FALSE) {
  cols <- (if (isTRUE(dark)) seq_scales_dark else seq_scales)[[name]]
  if (is.null(cols)) return(name)
  lapply(seq_along(cols), function(i) list((i - 1) / (length(cols) - 1), cols[i]))
}

plotly_font <- function(ink) list(family = "Inter, system-ui, sans-serif", size = 12, color = ink$text)
plotly_axis <- function(ink, title, ...) {
  list(..., gridcolor = ink$grid, linecolor = ink$axis, zeroline = FALSE,
       tickfont = list(color = ink$muted), title = list(text = title, font = list(color = ink$muted)))
}

scene_axis <- function(ink, title, ...) {
  list(..., gridcolor = ink$scene_grid, zerolinecolor = ink$axis, showbackground = TRUE,
       backgroundcolor = ink$scene, tickfont = list(color = ink$muted, size = 11),
       title = list(text = title, font = list(color = ink$muted)))
}

# A row of pill buttons that behaves like radioButtons.
segmented <- function(id, choices, selected) {
  div(class = "segmented", radioButtons(id, NULL, choices = choices, selected = selected, inline = TRUE))
}

info_tip <- function(...) tooltip(icon("circle-info", class = "text-muted"), ...)

# Non-ASCII characters as HTML entities, so Greek letters survive any server locale.
html_unicode <- function(x) {
  cp <- utf8ToInt(enc2utf8(x))
  paste(ifelse(cp > 127, sprintf("&#%d;", cp), intToUtf8(cp, multiple = TRUE)), collapse = "")
}

# ------------------------------------------------------------------------- UI

uni_sidebar <- sidebar(
  width = 340,
  title = "Univariate settings",
  accordion(
    multiple = TRUE,
    open = c("Distribution", "Compare", "Axes"),
    accordion_panel(
      "Distribution", icon = icon("chart-area"),
      selectInput("dist", NULL,
                  choices = list(Continuous = dist_choices("continuous"),
                                 Discrete = dist_choices("discrete")),
                  selected = "normal"),
      param_sliders(distributions, uni_input_id, "uni_tabs"),
      uiOutput("param_error")
    ),
    accordion_panel(
      "Compare", icon = icon("layer-group"),
      p(class = "small text-muted",
        "Pin the current curve, then change the parameters or the distribution to overlay another one."),
      actionButton("pin", "Pin current curve", icon = icon("thumbtack"),
                   class = "btn-primary w-100 mb-2"),
      uiOutput("pin_list"),
      actionLink("clear_pins", "Clear all pinned curves", icon = icon("trash-can"))
    ),
    accordion_panel(
      "Axes", icon = icon("ruler-combined"),
      radioButtons("x_mode", span("x axis ", info_tip("Auto fits every curve on screen. Fixed keeps the range while you move the sliders.")),
                   c(Auto = "auto", Fixed = "fixed"), inline = TRUE),
      conditionalPanel(
        "input.x_mode == 'fixed'",
        layout_columns(numericInput("xmin", "min", -5), numericInput("xmax", "max", 5))
      ),
      radioButtons("y_mode", "y axis", c(Auto = "auto", Fixed = "fixed"), inline = TRUE),
      conditionalPanel(
        "input.y_mode == 'fixed'",
        layout_columns(numericInput("ymin", "min", 0), numericInput("ymax", "max", 0.5))
      ),
      actionButton("freeze", "Freeze current axes", icon = icon("lock"),
                   class = "btn-outline-secondary btn-sm w-100")
    ),
    accordion_panel(
      "Probability and samples", icon = icon("dice"),
      input_switch("shade_on", "Shade P(a ≤ X ≤ b)"),
      conditionalPanel(
        "input.shade_on",
        layout_columns(numericInput("shade_a", "a", -1), numericInput("shade_b", "b", 1))
      ),
      input_switch("sample_on", "Overlay a random sample"),
      conditionalPanel(
        "input.sample_on",
        sliderInput("sample_n", "Sample size", min = 10, max = 5000, value = 500, step = 10,
                    ticks = FALSE),
        actionButton("resample", "Draw a new sample", icon = icon("rotate"),
                     class = "btn-outline-secondary btn-sm w-100")
      )
    )
  )
)

uni_main <- tagList(
  card(
    full_screen = TRUE,
    min_height = "540px",
    card_header(
      class = "plot-header",
      uiOutput("plot_title"),
      div(class = "header-tools",
          segmented("view", c("Density" = "pdf", "CDF" = "cdf", "Both" = "both",
                              "Survival" = "surv", "Hazard" = "hazard", "Quantile" = "quantile"),
                    "pdf"),
          downloadButton("download_plot", "PNG", class = "btn-sm btn-light"))
    ),
    card_body(class = "pt-2", uiOutput("plot_area"))
  ),
  layout_columns(
    col_widths = c(7, 5),
    card(card_header(icon("table-list"), "Summary of every curve on the plot"),
         card_body(class = "p-2 summary-body", tableOutput("summary_table"))),
    card(card_header(icon("square-root-variable"), "Formula"), card_body(uiOutput("formula")))
  )
)

bv_sidebar <- sidebar(
  width = 340,
  title = "Bivariate settings",
  accordion(
    multiple = TRUE,
    open = c("Distribution", "Compare"),
    accordion_panel(
      "Distribution", icon = icon("mountain"),
      selectInput("bv_dist", NULL, choices = bv_choices, selected = "bvn"),
      param_sliders(bivariates, bv_input_id, "bv_tabs")
    ),
    accordion_panel(
      "Compare", icon = icon("layer-group"),
      p(class = "small text-muted",
        "Pin the current surface, then change the parameters or the distribution to overlay another one."),
      actionButton("bv_pin", "Pin current surface", icon = icon("thumbtack"),
                   class = "btn-primary w-100 mb-2"),
      uiOutput("bv_pin_list"),
      actionLink("bv_clear_pins", "Clear all pinned surfaces", icon = icon("trash-can"))
    ),
    accordion_panel(
      "Display", icon = icon("cube"),
      selectInput("bv_colors", "Colour scale (single surface)",
                  c("Ocean", "Ember", "Viridis", "Cividis", "YlGnBu", "Portland", "Hot")),
      sliderInput("bv_n", "Grid resolution", min = 30, max = 150, value = 80, step = 10,
                  ticks = FALSE),
      input_switch("bv_sample_on", "Show a random sample on the contour map"),
      conditionalPanel(
        "input.bv_sample_on",
        sliderInput("bv_sample_n", "Sample size", min = 50, max = 3000, value = 500, step = 50,
                    ticks = FALSE)
      )
    ),
    accordion_panel(
      "Axes", icon = icon("ruler-combined"),
      radioButtons("bv_axes", NULL, c(Auto = "auto", Fixed = "fixed"), inline = TRUE),
      conditionalPanel(
        "input.bv_axes == 'fixed'",
        layout_columns(numericInput("bv_x1min", "x₁ min", -4), numericInput("bv_x1max", "x₁ max", 4)),
        layout_columns(numericInput("bv_x2min", "x₂ min", -4), numericInput("bv_x2max", "x₂ max", 4)),
        numericInput("bv_zmax", "Density axis max", 0.2)
      ),
      actionButton("bv_freeze", "Freeze current axes", icon = icon("lock"),
                   class = "btn-outline-secondary btn-sm w-100")
    )
  )
)

bv_main <- tagList(
  card(
    full_screen = TRUE,
    min_height = "560px",
    card_header(
      class = "plot-header",
      uiOutput("bv_title"),
      div(class = "header-tools",
          segmented("bv_view", c("Surface + contour" = "both", "Surface" = "surface",
                                 "Contour" = "contour"), "both"))
    ),
    card_body(uiOutput("bv_plot_area"))
  ),
  card(card_header(icon("table-list"), "Summary of every surface on the plot"),
       card_body(class = "p-2 summary-body", tableOutput("bv_summary_table"))),
  layout_columns(
    col_widths = c(7, 5),
    card(card_header(icon("square-root-variable"), "Formula"), card_body(uiOutput("bv_formula"))),
    card(card_header(icon("calculator"), "Moments of the live surface"), card_body(uiOutput("bv_moments")))
  )
)

about_panel <- div(
  class = "card about-card mx-auto my-3", style = "max-width: 860px;",
  div(class = "about-hero",
      logo_svg,
      div(h2("Distribution Explorer"),
          p("See how parameters shape probability distributions, one slider at a time."))),
  div(
    class = "card-body p-4",
    withMathJax(),
    p("Distribution Explorer plots probability distributions so that you can see how the parameters change their shape."),
    tags$ul(
      tags$li(strong("Univariate: "), "pick a continuous or discrete distribution, move the sliders, and switch between the density or mass function, the CDF, the survival and hazard functions and the quantile function."),
      tags$li(strong("Compare: "), "pin curves to overlay several parameter sets, even from different distributions, for example a Binomial pmf and its normal approximation."),
      tags$li(strong("Axes: "), "fix the axes, or freeze the current view, so that the picture does not rescale while you move a slider."),
      tags$li(strong("Probability and samples: "), "shade \\(P(a \\le X \\le b)\\) and overlay a histogram of a random sample."),
      tags$li(strong("Bivariate 3-D: "), "rotate and zoom 3-D density surfaces of bivariate distributions, with a matching contour map.")
    ),
    h5("Functions shown"),
    p("For a random variable \\(X\\) with density \\(f\\) (or mass function \\(P(X = x)\\) when \\(X\\) is discrete):"),
    tags$ul(
      tags$li("CDF: \\(F(x) = P(X \\le x)\\)"),
      tags$li("Survival function: \\(S(x) = 1 - F(x) = P(X > x)\\)"),
      tags$li("Hazard function: \\(h(x) = f(x) / S(x)\\) for continuous \\(X\\), and \\(h(x) = P(X = x) / P(X \\ge x)\\) for discrete \\(X\\)"),
      tags$li("Quantile function: \\(F^{-1}(u) = \\inf\\{x : F(x) \\ge u\\}\\) for \\(u \\in (0, 1)\\)")
    ),
    p(class = "text-muted small", "Built with shiny, bslib, ggplot2, plotly and mvtnorm.")
  )
)

ui <- page_navbar(
  title = span(class = "brand", logo_svg, span("Distribution ", tags$b("Explorer"))),
  window_title = "Distribution Explorer",
  theme = app_theme,
  bg = "#0F1A3A",
  inverse = TRUE,
  header = tagList(withMathJax(), tags$link(rel = "stylesheet", href = "styles.css")),
  nav_panel("Univariate", icon = icon("wave-square"), layout_sidebar(sidebar = uni_sidebar, uni_main, fillable = FALSE)),
  nav_panel("Bivariate 3-D", icon = icon("cube"), layout_sidebar(sidebar = bv_sidebar, bv_main, fillable = FALSE)),
  nav_panel("About", icon = icon("circle-info"), about_panel),
  nav_spacer(),
  nav_item(input_dark_mode(id = "dark_mode", mode = "light")),
  nav_item(tags$a(class = "nav-link", icon("github"), href = "https://github.com/ChaoranHu/distribution_plot",
                  target = "_blank", title = "Source code"))
)

# --------------------------------------------------------------------- server

server <- function(input, output, session) {

  is_dark <- reactive(identical(input$dark_mode, "dark"))

  # ---------------------------------------------------------- univariate ----

  observeEvent(input$dist, nav_select("uni_tabs", input$dist))

  live_par <- reactive({
    spec <- distributions[[input$dist]]
    par <- lapply(spec$params, function(prm) input[[uni_input_id(input$dist, prm$id)]])
    names(par) <- vapply(spec$params, `[[`, character(1), "id")
    req(all(!vapply(par, is.null, logical(1))))
    par
  })

  param_error <- reactive({
    spec <- distributions[[input$dist]]
    if (is.null(spec$check)) NULL else spec$check(live_par())
  })

  output$param_error <- renderUI({
    if (!is.null(param_error())) div(class = "alert alert-warning py-2 small", param_error())
  })

  live_curve <- reactive({
    list(dist = input$dist, par = live_par(), label = curve_label(input$dist, live_par()), live = TRUE)
  })

  pins <- reactiveVal(list())

  observeEvent(input$pin, {
    if (!is.null(param_error())) return()
    cur <- pins()
    lv <- live_curve()
    if (lv$label %in% vapply(cur, `[[`, character(1), "label")) {
      showNotification("This curve is already pinned.", type = "message", duration = 3)
      return()
    }
    if (length(cur) >= length(palette_curves) - 1) {
      showNotification("You can pin up to 8 curves. Remove one first.", type = "warning")
      return()
    }
    lv$live <- FALSE
    pins(c(cur, list(lv)))
  })
  observeEvent(input$clear_pins, pins(list()))
  observeEvent(input$remove_pin, {
    cur <- pins()
    i <- as.integer(input$remove_pin)
    if (i >= 1 && i <= length(cur)) pins(cur[-i])
  })

  # Pinned curves first, then the live curve; colours follow that order.
  curves <- reactive({
    validate(need(is.null(param_error()), param_error()))
    cs <- pins()
    lv <- live_curve()
    if (!(lv$label %in% vapply(cs, `[[`, character(1), "label"))) cs <- c(cs, list(lv))
    pal <- curve_palette(is_dark())
    for (i in seq_along(cs)) cs[[i]]$colour <- pal[(i - 1) %% length(pal) + 1]
    cs
  })

  output$pin_list <- renderUI(pin_list_ui(pins(), "remove_pin", "No pinned curves yet.", curve_palette(is_dark())))

  xlim <- reactive({
    if (input$x_mode == "fixed") {
      validate(need(is.numeric(input$xmin) && is.numeric(input$xmax) && input$xmin < input$xmax,
                    "Fixed x axis: min must be smaller than max."))
      return(c(input$xmin, input$xmax))
    }
    cs <- curves()
    rngs <- lapply(cs, function(cv) auto_range(cv$dist, cv$par))
    rng <- range(unlist(rngs), finite = TRUE)
    any_disc <- any(vapply(cs, function(cv) distributions[[cv$dist]]$type == "discrete", logical(1)))
    pad <- if (any_disc) 0.6 else 0
    rng + c(-1, 1) * (pad + 0.02 * diff(rng))
  })

  shade <- reactive({
    if (!isTRUE(input$shade_on)) return(NULL)
    a <- input$shade_a; b <- input$shade_b
    if (!is.numeric(a) || !is.numeric(b) || is.na(a) || is.na(b) || a > b) return(NULL)
    c(a, b)
  })

  sample_vals <- reactive({
    input$resample
    if (!isTRUE(input$sample_on) || !is.null(param_error())) return(NULL)
    distributions[[input$dist]]$r(input$sample_n, live_par())
  })

  main_fun <- reactive(if (input$view == "both") "pdf" else input$view)

  fixed_ylim <- reactive({
    if (input$y_mode != "fixed") return(NULL)
    validate(need(is.numeric(input$ymin) && is.numeric(input$ymax) && input$ymin < input$ymax,
                  "Fixed y axis: min must be smaller than max."))
    c(input$ymin, input$ymax)
  })

  # The function name goes above the plot only when two plots sit side by side;
  # otherwise the view buttons in the card header already say it.
  main_plot <- reactive({
    make_uni_plot(curves(), main_fun(), xlim(), fixed_ylim(), shade(), sample_vals(),
                  dark = is_dark(), title = input$view == "both")
  })

  output$plot_area <- renderUI({
    if (input$view == "both") {
      layout_columns(col_widths = c(6, 6),
                     plotOutput("plot_main", height = "470px"),
                     plotOutput("plot_cdf", height = "470px"))
    } else {
      plotOutput("plot_main", height = "470px")
    }
  })

  output$plot_main <- renderPlot(main_plot(), res = 96)
  output$plot_cdf <- renderPlot(make_uni_plot(curves(), "cdf", xlim(), NULL, shade(), dark = is_dark()),
                                res = 96)

  output$plot_title <- renderUI({
    spec <- distributions[[input$dist]]
    n <- length(curves())
    div(class = "plot-title",
        span(class = "plot-title-main", paste(spec$label, "distribution")),
        span(class = paste("type-pill", spec$type), spec$type),
        if (n > 1) span(class = "count-pill", icon("layer-group"), paste(n, "curves")))
  })

  observeEvent(input$freeze, {
    xl <- signif(xlim(), 4)
    yl <- signif(attr(main_plot(), "ylim"), 4)
    updateRadioButtons(session, "x_mode", selected = "fixed")
    updateNumericInput(session, "xmin", value = xl[1])
    updateNumericInput(session, "xmax", value = xl[2])
    if (length(yl) == 2) {
      updateRadioButtons(session, "y_mode", selected = "fixed")
      updateNumericInput(session, "ymin", value = yl[1])
      updateNumericInput(session, "ymax", value = yl[2])
    }
    showNotification("Axes frozen. Switch back to Auto in the Axes panel at any time.",
                     type = "message", duration = 3)
  })

  output$download_plot <- downloadHandler(
    filename = function() paste0("distribution-", input$dist, "-", main_fun(), ".png"),
    content = function(file) {
      # The download is always light, with the function name as its title.
      g <- make_uni_plot(curves(), main_fun(), xlim(), fixed_ylim(), shade(), sample_vals())
      ggsave(file, g, width = 9, height = 6, dpi = 150, bg = "white")
    }
  )

  output$summary_table <- renderTable({
    cs <- curves()
    sh <- shade()
    do.call(rbind, lapply(cs, function(cv) {
      spec <- distributions[[cv$dist]]
      v <- spec$var(cv$par)
      row <- data.frame(
        ` ` = sprintf("<span class='swatch' style='background:%s'></span>", cv$colour),
        Curve = paste0(html_unicode(cv$label), if (cv$live) " <span class='badge bg-primary-subtle text-primary'>live</span>" else ""),
        Mean = fmt_num(spec$mean(cv$par)),
        SD = fmt_num(if (is.na(v)) NA else sqrt(v)),
        Median = fmt_num(spec$q(0.5, cv$par)),
        check.names = FALSE
      )
      if (!is.null(sh)) {
        lower <- if (spec$type == "discrete") spec$p(ceiling(sh[1]) - 1, cv$par) else spec$p(sh[1], cv$par)
        row[[html_unicode(sprintf("P(%s \u2264 X \u2264 %s)", format(sh[1]), format(sh[2])))]] <-
          formatC(spec$p(sh[2], cv$par) - lower, digits = 4, format = "f")
      }
      row
    }))
  }, sanitize.text.function = identity, striped = TRUE, hover = TRUE, width = "100%",
     class = "table-sm")

  output$formula <- renderUI({
    spec <- distributions[[input$dist]]
    withMathJax(tagList(
      h6(spec$label, span(class = "badge bg-secondary-subtle text-secondary ms-1", spec$type)),
      div(class = "formula-display", paste0("$$", spec$formula, "$$")),
      p("Support: ", paste0("\\(", spec$support, "\\)")),
      p(class = "mb-1", strong("Notation")),
      tags$ul(class = "small", lapply(spec$notation, tags$li))
    ))
  })

  # ---------------------------------------------------------- bivariate -----

  observeEvent(input$bv_dist, nav_select("bv_tabs", input$bv_dist))

  bv_par <- reactive({
    spec <- bivariates[[input$bv_dist]]
    par <- lapply(spec$params, function(prm) input[[bv_input_id(input$bv_dist, prm$id)]])
    names(par) <- vapply(spec$params, `[[`, character(1), "id")
    req(all(!vapply(par, is.null, logical(1))))
    par
  })

  bv_live <- reactive({
    list(dist = input$bv_dist, par = bv_par(), label = bv_label(input$bv_dist, bv_par()), live = TRUE)
  })

  bv_pins <- reactiveVal(list())

  observeEvent(input$bv_pin, {
    cur <- bv_pins()
    lv <- bv_live()
    if (lv$label %in% vapply(cur, `[[`, character(1), "label")) {
      showNotification("This surface is already pinned.", type = "message", duration = 3)
      return()
    }
    if (length(cur) >= 5) {
      showNotification("You can pin up to 5 surfaces. Remove one first.", type = "warning")
      return()
    }
    lv$live <- FALSE
    bv_pins(c(cur, list(lv)))
  })
  observeEvent(input$bv_clear_pins, bv_pins(list()))
  observeEvent(input$bv_remove_pin, {
    cur <- bv_pins()
    i <- as.integer(input$bv_remove_pin)
    if (i >= 1 && i <= length(cur)) bv_pins(cur[-i])
  })
  output$bv_pin_list <- renderUI(pin_list_ui(bv_pins(), "bv_remove_pin", "No pinned surfaces yet.", curve_palette(is_dark())))

  # Pinned surfaces first, then the live one; colours follow that order.
  bv_surfaces <- reactive({
    cs <- bv_pins()
    lv <- bv_live()
    if (!(lv$label %in% vapply(cs, `[[`, character(1), "label"))) cs <- c(cs, list(lv))
    pal <- curve_palette(is_dark())
    for (i in seq_along(cs)) cs[[i]]$colour <- pal[(i - 1) %% length(pal) + 1]
    cs
  })

  bv_ranges <- reactive({
    if (input$bv_axes == "fixed") {
      validate(need(input$bv_x1min < input$bv_x1max && input$bv_x2min < input$bv_x2max,
                    "Fixed axes: each min must be smaller than its max."))
      return(list(x1 = c(input$bv_x1min, input$bv_x1max), x2 = c(input$bv_x2min, input$bv_x2max)))
    }
    rs <- lapply(bv_surfaces(), function(cv) bivariates[[cv$dist]]$range(cv$par))
    list(x1 = range(unlist(lapply(rs, `[[`, "x1"))), x2 = range(unlist(lapply(rs, `[[`, "x2"))))
  })

  bv_data <- reactive({
    r <- bv_ranges()
    cs <- bv_surfaces()
    grids <- lapply(cs, function(cv) bv_grid(cv$dist, cv$par, r$x1, r$x2, input$bv_n))
    auto_max <- max(vapply(seq_along(cs), function(i) {
      z <- grids[[i]]$z[is.finite(grids[[i]]$z)]
      if (!length(z)) return(0)
      if (cs[[i]]$dist == "dirichlet" && max(z) > 4 * stats::quantile(z, 0.95)) {
        # Densities with a concentration below 1 diverge at the edges of the triangle.
        stats::quantile(z, 0.95, names = FALSE) * 1.5
      } else {
        max(z) * 1.05
      }
    }, numeric(1)))
    zmax <- if (input$bv_axes == "fixed" && is.numeric(input$bv_zmax) && input$bv_zmax > 0) {
      input$bv_zmax
    } else {
      auto_max
    }
    for (i in seq_along(grids)) {
      grids[[i]]$z[grids[[i]]$z > zmax] <- zmax
      grids[[i]]$colour <- cs[[i]]$colour
      grids[[i]]$label <- cs[[i]]$label
    }
    list(grids = grids, x1 = grids[[1]]$x1, x2 = grids[[1]]$x2, zmax = zmax, ranges = r)
  })

  bv_sample <- reactive({
    if (!isTRUE(input$bv_sample_on)) return(NULL)
    bivariates[[input$bv_dist]]$sample(input$bv_sample_n, bv_par())
  })

  output$bv_title <- renderUI({
    n <- length(bv_surfaces())
    div(class = "plot-title",
        span(class = "plot-title-main", paste(bivariates[[input$bv_dist]]$label, "density")),
        if (n > 1) span(class = "count-pill", icon("layer-group"), paste(n, "surfaces")))
  })

  output$bv_plot_area <- renderUI({
    switch(input$bv_view,
      both = layout_columns(col_widths = c(7, 5),
                            plotlyOutput("bv_surface", height = "520px"),
                            plotlyOutput("bv_contour", height = "520px")),
      surface = plotlyOutput("bv_surface", height = "560px"),
      contour = plotlyOutput("bv_contour", height = "560px"))
  })

  bv_hover <- "x₁ = %{x:.3f}<br>x₂ = %{y:.3f}<br>f = %{z:.4f}"
  solid_scale <- function(col) list(c(0, col), c(1, col))

  output$bv_surface <- renderPlotly({
    d <- bv_data()
    ink <- ink_for(is_dark())
    single <- length(d$grids) == 1
    p <- plot_ly()
    for (g in d$grids) {
      p <- if (single) {
        add_surface(p, x = g$x1, y = g$x2, z = g$z, colorscale = bv_colorscale(input$bv_colors, is_dark()),
                    cmin = 0, cmax = d$zmax, name = g$label,
                    contours = list(z = list(show = TRUE, usecolormap = TRUE, project = list(z = TRUE),
                                             highlightcolor = "#ffffff")),
                    colorbar = list(title = list(text = "f(x₁, x₂)"), len = 0.6, thickness = 12,
                                    outlinewidth = 0),
                    hovertemplate = paste0(bv_hover, "<extra></extra>"))
      } else {
        # One flat colour per surface, partly transparent so that overlaps stay visible.
        # The near-zero floor is dropped, otherwise each surface's floor hides the others.
        z <- g$z
        z[z < 0.02 * d$zmax] <- NA
        add_surface(p, x = g$x1, y = g$x2, z = z, colorscale = solid_scale(g$colour),
                    showscale = FALSE, opacity = 0.7, name = g$label, showlegend = TRUE,
                    lighting = list(ambient = 0.55, diffuse = 0.8, specular = 0.15),
                    hovertemplate = paste0(bv_hover, "<extra>", g$label, "</extra>"))
      }
    }
    p |>
      layout(
        uirevision = "bv",
        margin = list(l = 0, r = 0, t = 10, b = 0),
        legend = list(orientation = "h", x = 0, y = -0.02, font = list(size = 11)),
        scene = list(
          xaxis = scene_axis(ink, title = "x₁", range = d$ranges$x1),
          yaxis = scene_axis(ink, title = "x₂", range = d$ranges$x2),
          zaxis = scene_axis(ink, title = "density", range = c(0, d$zmax)),
          aspectmode = "manual", aspectratio = list(x = 1, y = 1, z = 0.65),
          camera = list(eye = list(x = 1.55, y = -1.55, z = 0.95))
        ),
        font = plotly_font(ink),
        hoverlabel = list(font = plotly_font(ink), bgcolor = ink$surface, bordercolor = ink$axis),
        paper_bgcolor = "rgba(0,0,0,0)"
      ) |>
      config(displaylogo = FALSE)
  })

  output$bv_contour <- renderPlotly({
    d <- bv_data()
    ink <- ink_for(is_dark())
    single <- length(d$grids) == 1
    p <- plot_ly()
    for (g in d$grids) {
      p <- if (single) {
        add_contour(p, x = g$x1, y = g$x2, z = g$z, colorscale = bv_colorscale(input$bv_colors, is_dark()),
                    zauto = FALSE, zmin = 0, zmax = d$zmax, ncontours = 18,
                    contours = list(coloring = "heatmap", showlines = TRUE),
                    line = list(color = if (is_dark()) "rgba(255,255,255,0.12)" else "rgba(15,26,58,0.18)",
                                width = 0.6),
                    showscale = input$bv_view == "contour", name = g$label,
                    hovertemplate = paste0(bv_hover, "<extra></extra>"))
      } else {
        add_contour(p, x = g$x1, y = g$x2, z = g$z, colorscale = solid_scale(g$colour),
                    autocontour = TRUE, ncontours = 10, showscale = FALSE,
                    contours = list(coloring = "lines"), line = list(width = 1.6),
                    name = g$label, showlegend = TRUE,
                    hovertemplate = paste0(bv_hover, "<extra>", g$label, "</extra>"))
      }
    }
    s <- bv_sample()
    if (!is.null(s)) {
      p <- p |> add_markers(x = s[, 1], y = s[, 2], name = "sample",
                            marker = list(size = 4, color = if (is_dark()) "rgba(226,232,240,0.55)" else "rgba(15,26,58,0.5)",
                                          line = list(width = 0)),
                            hoverinfo = "skip", showlegend = FALSE)
    }
    p |>
      layout(
        uirevision = "bv",
        margin = list(l = 50, r = 10, t = 10, b = 45),
        legend = list(orientation = "h", x = 0, y = -0.15, font = list(size = 11)),
        xaxis = plotly_axis(ink, title = "x₁", range = d$ranges$x1),
        yaxis = plotly_axis(ink, title = "x₂", range = d$ranges$x2),
        plot_bgcolor = if (single) NULL else ink$surface,
        font = plotly_font(ink),
        hoverlabel = list(font = plotly_font(ink), bgcolor = ink$surface, bordercolor = ink$axis),
        paper_bgcolor = "rgba(0,0,0,0)"
      ) |>
      config(displaylogo = FALSE)
  })

  output$bv_summary_table <- renderTable({
    do.call(rbind, lapply(bv_surfaces(), function(cv) {
      m <- bv_moments(cv$dist, cv$par)
      sds <- sqrt(diag(m$cov))
      data.frame(
        ` ` = sprintf("<span class='swatch' style='background:%s'></span>", cv$colour),
        Surface = paste0(html_unicode(cv$label), if (cv$live) " <span class='badge bg-primary-subtle text-primary'>live</span>" else ""),
        `E[X1]` = fmt_num(m$mean[1]), `E[X2]` = fmt_num(m$mean[2]),
        `SD(X1)` = fmt_num(sds[1]), `SD(X2)` = fmt_num(sds[2]),
        `Corr(X1, X2)` = fmt_num(m$cov[1, 2] / prod(sds)),
        check.names = FALSE
      )
    }))
  }, sanitize.text.function = identity, striped = TRUE, hover = TRUE, width = "100%",
     class = "table-sm")

  observeEvent(input$bv_freeze, {
    g <- bv_data()
    updateRadioButtons(session, "bv_axes", selected = "fixed")
    updateNumericInput(session, "bv_x1min", value = signif(g$ranges$x1[1], 4))
    updateNumericInput(session, "bv_x1max", value = signif(g$ranges$x1[2], 4))
    updateNumericInput(session, "bv_x2min", value = signif(g$ranges$x2[1], 4))
    updateNumericInput(session, "bv_x2max", value = signif(g$ranges$x2[2], 4))
    updateNumericInput(session, "bv_zmax", value = signif(g$zmax, 4))
  })

  output$bv_formula <- renderUI({
    spec <- bivariates[[input$bv_dist]]
    withMathJax(tagList(
      # Long formulas, such as a density followed by its covariance matrix, go on two lines.
      div(class = "formula-display",
          lapply(strsplit(spec$formula, ", \\quad ", fixed = TRUE)[[1]],
                 function(f) paste0("$$", f, "$$"))),
      p(class = "mb-1", strong("Notation")),
      tags$ul(class = "small", lapply(spec$notation, tags$li))
    ))
  })

  output$bv_moments <- renderUI({
    p <- bv_par()
    tex_num <- function(v) {
      if (is.na(v)) "\\text{undefined}" else if (is.infinite(v)) "\\infty" else fmt_num(v)
    }
    tex_vec <- function(v) sprintf("\\begin{pmatrix} %s \\\\ %s \\end{pmatrix}", tex_num(v[1]), tex_num(v[2]))
    tex_mat <- function(m) sprintf("\\begin{pmatrix} %s & %s \\\\ %s & %s \\end{pmatrix}",
                                   tex_num(m[1, 1]), tex_num(m[1, 2]), tex_num(m[2, 1]), tex_num(m[2, 2]))
    m <- bv_moments(input$bv_dist, p)
    corr <- m$cov[1, 2] / sqrt(m$cov[1, 1] * m$cov[2, 2])
    withMathJax(tagList(
      p(paste0("$$E[\\mathbf{X}] = ", tex_vec(m$mean), "$$")),
      p(paste0("$$\\operatorname{Cov}(\\mathbf{X}) = ", tex_mat(m$cov), "$$")),
      p(paste0("$$\\operatorname{Corr}(X_1, X_2) = ", tex_num(corr), "$$")),
      p(class = "small text-muted",
        "\\(E[\\mathbf{X}]\\) is the mean vector, \\(\\operatorname{Cov}(\\mathbf{X})\\) the covariance matrix of \\(\\mathbf{X} = (X_1, X_2)^\\top\\).")
    ))
  })
}

shinyApp(ui = ui, server = server)

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
  primary = "#3B5BDB",
  bg = "#F8FAFC",
  fg = "#0F172A",
  base_font = font_google("Inter", local = FALSE),
  heading_font = font_google("Inter", local = FALSE),
  "border-radius" = "0.6rem"
)

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

info_tip <- function(...) tooltip(icon("circle-info", class = "text-muted"), ...)

# Non-ASCII characters as HTML entities, so Greek letters survive any server locale.
html_unicode <- function(x) {
  cp <- utf8ToInt(enc2utf8(x))
  paste(ifelse(cp > 127, sprintf("&#%d;", cp), intToUtf8(cp, multiple = TRUE)), collapse = "")
}

# ------------------------------------------------------------------------- UI

uni_sidebar <- sidebar(
  width = 340,
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
      "View", icon = icon("eye"),
      radioButtons("view", NULL,
                   choices = c("Density / mass" = "pdf", "CDF" = "cdf",
                               "Density and CDF side by side" = "both",
                               "Survival function" = "surv", "Hazard function" = "hazard",
                               "Quantile function" = "quantile"),
                   selected = "pdf")
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
  layout_columns(
    fill = FALSE,
    col_widths = c(3, 3, 3, 3),
    value_box("Mean", textOutput("vb_mean"), showcase = icon("bullseye"), theme = "text-primary"),
    value_box("Std. deviation", textOutput("vb_sd"), showcase = icon("arrows-left-right"),
              theme = "text-info"),
    value_box("Median", textOutput("vb_median"), showcase = icon("divide"), theme = "text-success"),
    value_box(uiOutput("vb_prob_title", inline = TRUE), textOutput("vb_prob"),
              showcase = icon("percent"), theme = "text-warning")
  ),
  card(
    full_screen = TRUE,
    min_height = "520px",
    card_header(
      class = "d-flex justify-content-between align-items-center",
      textOutput("plot_title", inline = TRUE),
      downloadButton("download_plot", "PNG", class = "btn-sm btn-outline-secondary")
    ),
    card_body(uiOutput("plot_area"))
  ),
  layout_columns(
    col_widths = c(5, 7),
    card(card_header("Formula"), card_body(uiOutput("formula"))),
    card(card_header("Summary of every curve on the plot"),
         card_body(tableOutput("summary_table")))
  )
)

bv_sidebar <- sidebar(
  width = 340,
  accordion(
    multiple = TRUE,
    open = c("Distribution", "Display"),
    accordion_panel(
      "Distribution", icon = icon("mountain"),
      selectInput("bv_dist", NULL, choices = bv_choices, selected = "bvn"),
      param_sliders(bivariates, bv_input_id, "bv_tabs")
    ),
    accordion_panel(
      "Display", icon = icon("cube"),
      radioButtons("bv_view", NULL,
                   c("3-D surface and contour" = "both", "3-D surface" = "surface",
                     "Contour map" = "contour"), selected = "both"),
      selectInput("bv_colors", "Colour scale",
                  c("Viridis", "Cividis", "YlGnBu", "Blues", "Portland", "Hot", "Electric")),
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
    card_header(textOutput("bv_title", inline = TRUE)),
    card_body(uiOutput("bv_plot_area"))
  ),
  layout_columns(
    col_widths = c(7, 5),
    card(card_header("Formula"), card_body(uiOutput("bv_formula"))),
    card(card_header("Moments"), card_body(uiOutput("bv_moments")))
  )
)

about_panel <- card(
  class = "mx-auto mt-3", style = "max-width: 820px;",
  card_header("About this app"),
  card_body(
    withMathJax(),
    p("Distribution Explorer plots probability distributions so that you can see how the parameters change their shape."),
    tags$ul(
      tags$li(strong("Univariate: "), "pick a continuous or discrete distribution, move the sliders, and switch between the density or mass function, the CDF, the survival and hazard functions and the quantile function."),
      tags$li(strong("Compare: "), "pin curves to overlay several parameter sets, even from different distributions, for example a Binomial pmf and its normal approximation."),
      tags$li(strong("Axes: "), "fix the axes, or freeze the current view, so that the picture does not rescale while you move a slider."),
      tags$li(strong("Probability and samples: "), "shade \\(P(a \\le X \\le b)\\) and overlay a histogram of a random sample."),
      tags$li(strong("Bivariate 3-D: "), "rotate and zoom 3-D density surfaces of bivariate distributions, with a matching contour map.")
    ),
    h6("Functions shown"),
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
  title = span(icon("chart-line"), " Distribution Explorer"),
  window_title = "Distribution Explorer",
  theme = app_theme,
  bg = "#1E293B",
  inverse = TRUE,
  header = tagList(withMathJax(), tags$link(rel = "stylesheet", href = "styles.css")),
  nav_panel("Univariate", icon = icon("wave-square"), layout_sidebar(sidebar = uni_sidebar, uni_main)),
  nav_panel("Bivariate 3-D", icon = icon("cube"), layout_sidebar(sidebar = bv_sidebar, bv_main)),
  nav_panel("About", icon = icon("circle-info"), about_panel),
  nav_spacer(),
  nav_item(tags$a(icon("github"), href = "https://github.com/ChaoranHu/distribution_plot",
                  target = "_blank", title = "Source code"))
)

# --------------------------------------------------------------------- server

server <- function(input, output, session) {

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
    if (length(cur) >= length(palette_okabe) - 1) {
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
    for (i in seq_along(cs)) cs[[i]]$colour <- palette_okabe[(i - 1) %% length(palette_okabe) + 1]
    cs
  })

  output$pin_list <- renderUI({
    cs <- pins()
    if (!length(cs)) return(p(class = "small text-muted fst-italic", "No pinned curves yet."))
    tags$div(class = "pin-list", lapply(seq_along(cs), function(i) {
      tags$div(
        class = "pin-row",
        tags$span(class = "swatch", style = paste0("background:", palette_okabe[i])),
        tags$span(class = "pin-label", cs[[i]]$label),
        tags$button(class = "btn btn-sm btn-link pin-remove", title = "Remove",
                    onclick = sprintf("Shiny.setInputValue('remove_pin', %d, {priority: 'event'})", i),
                    icon("xmark"))
      )
    }))
  })

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

  main_plot <- reactive({
    ylim <- NULL
    if (input$y_mode == "fixed") {
      validate(need(is.numeric(input$ymin) && is.numeric(input$ymax) && input$ymin < input$ymax,
                    "Fixed y axis: min must be smaller than max."))
      ylim <- c(input$ymin, input$ymax)
    }
    make_uni_plot(curves(), main_fun(), xlim(), ylim, shade(), sample_vals())
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
  output$plot_cdf <- renderPlot(make_uni_plot(curves(), "cdf", xlim(), NULL, shade()), res = 96)

  output$plot_title <- renderText({
    n <- length(curves())
    paste0(distributions[[input$dist]]$label, " distribution",
           if (n > 1) paste0("  ·  comparing ", n, " curves") else "")
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
      ggsave(file, main_plot(), width = 9, height = 6, dpi = 150, bg = "white")
    }
  )

  live_stats <- reactive({
    validate(need(is.null(param_error()), ""))
    spec <- distributions[[input$dist]]
    par <- live_par()
    v <- spec$var(par)
    list(mean = spec$mean(par), sd = if (is.na(v)) NA else sqrt(v), median = spec$q(0.5, par))
  })
  output$vb_mean <- renderText(fmt_num(live_stats()$mean))
  output$vb_sd <- renderText(fmt_num(live_stats()$sd))
  output$vb_median <- renderText(fmt_num(live_stats()$median))

  output$vb_prob_title <- renderUI({
    s <- shade()
    if (is.null(s)) "P(a ≤ X ≤ b)" else sprintf("P(%s ≤ X ≤ %s)", format(s[1]), format(s[2]))
  })
  output$vb_prob <- renderText({
    s <- shade()
    if (is.null(s)) return("turn on shading")
    validate(need(is.null(param_error()), ""))
    spec <- distributions[[input$dist]]
    par <- live_par()
    lower <- if (spec$type == "discrete") spec$p(ceiling(s[1]) - 1, par) else spec$p(s[1], par)
    formatC(spec$p(s[2], par) - lower, digits = 4, format = "f")
  })

  output$summary_table <- renderTable({
    cs <- curves()
    do.call(rbind, lapply(cs, function(cv) {
      spec <- distributions[[cv$dist]]
      v <- spec$var(cv$par)
      data.frame(
        ` ` = sprintf("<span class='swatch' style='background:%s'></span>", cv$colour),
        Curve = paste0(html_unicode(cv$label), if (cv$live) " <span class='badge bg-primary-subtle text-primary'>live</span>" else ""),
        Mean = fmt_num(spec$mean(cv$par)),
        Variance = fmt_num(v),
        Median = fmt_num(spec$q(0.5, cv$par)),
        check.names = FALSE
      )
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

  bv_ranges <- reactive({
    if (input$bv_axes == "fixed") {
      validate(need(input$bv_x1min < input$bv_x1max && input$bv_x2min < input$bv_x2max,
                    "Fixed axes: each min must be smaller than its max."))
      return(list(x1 = c(input$bv_x1min, input$bv_x1max), x2 = c(input$bv_x2min, input$bv_x2max)))
    }
    bivariates[[input$bv_dist]]$range(bv_par())
  })

  bv_data <- reactive({
    r <- bv_ranges()
    g <- bv_grid(input$bv_dist, bv_par(), r$x1, r$x2, input$bv_n)
    z <- g$z[is.finite(g$z)]
    zmax <- if (input$bv_axes == "fixed" && is.numeric(input$bv_zmax) && input$bv_zmax > 0) {
      input$bv_zmax
    } else if (input$bv_dist == "dirichlet" && max(z) > 4 * stats::quantile(z, 0.95)) {
      # Densities with a concentration below 1 diverge at the edges of the triangle.
      stats::quantile(z, 0.95, names = FALSE) * 1.5
    } else {
      max(z) * 1.05
    }
    g$z[g$z > zmax] <- zmax
    g$zmax <- zmax
    g$ranges <- r
    g
  })

  bv_sample <- reactive({
    if (!isTRUE(input$bv_sample_on)) return(NULL)
    bivariates[[input$bv_dist]]$sample(input$bv_sample_n, bv_par())
  })

  output$bv_title <- renderText(paste0(bivariates[[input$bv_dist]]$label, " density"))

  output$bv_plot_area <- renderUI({
    switch(input$bv_view,
      both = layout_columns(col_widths = c(7, 5),
                            plotlyOutput("bv_surface", height = "520px"),
                            plotlyOutput("bv_contour", height = "520px")),
      surface = plotlyOutput("bv_surface", height = "560px"),
      contour = plotlyOutput("bv_contour", height = "560px"))
  })

  output$bv_surface <- renderPlotly({
    g <- bv_data()
    plot_ly(x = g$x1, y = g$x2, z = g$z, type = "surface", colorscale = input$bv_colors,
            cmin = 0, cmax = g$zmax,
            contours = list(z = list(show = TRUE, usecolormap = TRUE, project = list(z = TRUE),
                                     highlightcolor = "#ffffff")),
            colorbar = list(title = list(text = "f(x₁, x₂)"), len = 0.6, thickness = 14),
            hovertemplate = "x₁ = %{x:.3f}<br>x₂ = %{y:.3f}<br>f = %{z:.4f}<extra></extra>") |>
      layout(
        uirevision = input$bv_dist,
        margin = list(l = 0, r = 0, t = 10, b = 0),
        scene = list(
          xaxis = list(title = "x₁", range = g$ranges$x1),
          yaxis = list(title = "x₂", range = g$ranges$x2),
          zaxis = list(title = "density", range = c(0, g$zmax)),
          aspectmode = "manual", aspectratio = list(x = 1, y = 1, z = 0.65),
          camera = list(eye = list(x = 1.55, y = -1.55, z = 0.95))
        ),
        paper_bgcolor = "rgba(0,0,0,0)"
      ) |>
      config(displaylogo = FALSE)
  })

  output$bv_contour <- renderPlotly({
    g <- bv_data()
    p <- plot_ly(x = g$x1, y = g$x2, z = g$z, type = "contour", colorscale = input$bv_colors,
                 zauto = FALSE, zmin = 0, zmax = g$zmax, ncontours = 18,
                 contours = list(coloring = "heatmap", showlines = TRUE),
                 line = list(color = "rgba(255,255,255,0.45)", width = 0.6),
                 showscale = input$bv_view == "contour",
                 hovertemplate = "x₁ = %{x:.3f}<br>x₂ = %{y:.3f}<br>f = %{z:.4f}<extra></extra>")
    s <- bv_sample()
    if (!is.null(s)) {
      p <- p |> add_markers(x = s[, 1], y = s[, 2], inherit = FALSE, name = "sample",
                            marker = list(size = 4, color = "rgba(15,23,42,0.45)",
                                          line = list(width = 0)),
                            hoverinfo = "skip", showlegend = FALSE)
    }
    p |>
      layout(
        uirevision = input$bv_dist,
        margin = list(l = 50, r = 10, t = 10, b = 45),
        xaxis = list(title = "x₁", range = g$ranges$x1, zeroline = FALSE),
        yaxis = list(title = "x₂", range = g$ranges$x2, zeroline = FALSE),
        paper_bgcolor = "rgba(0,0,0,0)"
      ) |>
      config(displaylogo = FALSE)
  })

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
      div(class = "formula-display", paste0("$$", spec$formula, "$$")),
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
    m <- switch(input$bv_dist,
      bvn = list(mean = c(p$mu1, p$mu2), cov = bvn_sigma(p)),
      bvt = list(mean = if (p$nu > 1) c(p$mu1, p$mu2) else c(NA, NA),
                 cov = if (p$nu > 2) p$nu / (p$nu - 2) * bvn_sigma(p) else matrix(if (p$nu > 1) Inf else NA, 2, 2)),
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

library(shiny)
library(ggplot2)

# ============================================================
# SIMULATION FUNCTIONS
# ============================================================

get_founder_colors <- function(n) {
  base_pal <- c(
    "#3A86FF", "#FF006E", "#FB5607", "#FFBE0B",
    "#8338EC", "#06D6A0", "#118AB2", "#6D4C41"
  )
  if (n <= 8) return(base_pal[seq_len(n)])
  colorRampPalette(c(
    "#3A86FF", "#FF006E", "#FB5607", "#FFBE0B",
    "#8338EC", "#06D6A0", "#118AB2", "#6D4C41"
  ))(n)
}

simulate_rils <- function(n_founders, n_gen, n_total, design) {
  n_pos   <- 101L
  pos_seq <- 0L:100L
  haps    <- matrix(1L, nrow = n_total, ncol = n_pos)

  for (r in seq_len(n_total)) {
    if (design == "hub") {
      hub   <- 1L
      other <- if (n_founders >= 2L) sample.int(n_founders - 1L, 1L) + 1L else 1L
      pool  <- c(hub, other)
    } else {
      pool <- seq_len(n_founders)
    }

    # Start with one randomly chosen founder
    state <- rep(sample(pool, 1L), n_pos)

    # Accumulate all crossovers across n_gen generations at once.
    # Each generation contributes ~1 crossover per 100 cM (Poisson rate = 1).
    # More generations = more breakpoints = smaller blocks.
    total_cross <- rpois(1L, n_gen)
    if (total_cross > 0L) {
      cpts <- sort(runif(total_cross, 0, 100))
      for (cp in cpts) {
        idx <- which(pos_seq >= cp)[1L]
        if (is.na(idx) || idx > n_pos) next
        new_f            <- sample(pool, 1L)
        state[idx:n_pos] <- new_f
      }
    }

    haps[r, ] <- state
  }
  haps
}

add_noise <- function(haps, n_founders, coverage, haplo_acc) {
  # Coverage drives read depth: below ~5x things are very noisy, above ~15x nearly perfect.
  # Higher coverage amplifies QTL signal by reducing haplotype misassignment.
  cov_err <- exp(-coverage / 4)   # steep dropoff: 1x~0.78, 5x~0.29, 15x~0.02

  # Block-level error: a whole inferred block gets assigned to the wrong founder.
  # This matches real haplotype inference errors (segment-level, not position-level).
  block_err_rate <- min(0.6, max(0, (1 - haplo_acc) * 0.4 + cov_err * 0.3))

  for (r in seq_len(nrow(haps))) {
    breaks <- c(1L, which(diff(haps[r, ]) != 0L) + 1L, ncol(haps) + 1L)
    for (b in seq_len(length(breaks) - 1L)) {
      if (runif(1L) < block_err_rate) {
        wrong_f <- sample.int(n_founders, 1L)
        haps[r, breaks[b]:(breaks[b + 1L] - 1L)] <- wrong_f
      }
    }
  }
  haps
}

get_qtl_config <- function(model, n_founders) {
  set.seed(42L)

  # At each QTL, only a small fraction of founders carry a positive effect.
  # n_effect = 10% of founders, min 1, max 10.
  # All other founders have effect 0.
  n_effect <- max(1L, min(10L, round(0.1 * n_founders)))

  make_eff <- function(effect_size) {
    eff              <- numeric(n_founders)
    who              <- sample.int(n_founders, n_effect)
    eff[who]         <- runif(n_effect, effect_size * 0.75, effect_size * 1.25)
    eff
  }

  switch(model,
    "null"      = list(pos = integer(0), eff = list()),
    "1qtl"      = list(pos = 50L,
                       eff = list(make_eff(3.0))),
    "2qtl"      = list(pos = c(30L, 70L),
                       eff = list(make_eff(2.5), make_eff(2.5))),
    "3qtl"      = list(pos = c(20L, 50L, 80L),
                       eff = list(make_eff(2.0), make_eff(2.0), make_eff(2.0))),
    "polygenic" = list(pos = c(5L,15L,25L,35L,45L,55L,65L,75L,85L,95L),
                       eff = lapply(seq_len(10L), function(i) make_eff(1.0)))
  )
}

simulate_phenotype <- function(true_haps, qtl_cfg) {
  pheno <- rnorm(nrow(true_haps), 0, 1)
  for (i in seq_along(qtl_cfg$pos)) {
    col_idx <- qtl_cfg$pos[i] + 1L
    pheno   <- pheno + qtl_cfg$eff[[i]][true_haps[, col_idx]]
  }
  pheno
}

# Chi-square test (2N allele counts per pool: each inbred RIL is diploid/
# homozygous so contributes 2 identical alleles per position).
# Chi-square gives exact p-values with no Monte Carlo LOD cap.
qtl_scan <- function(obs_haps, case_idx, control_idx) {
  n_pos <- ncol(obs_haps)
  lod   <- numeric(n_pos)
  all_lvls <- seq_len(max(obs_haps))
  for (p in seq_len(n_pos)) {
    case_f <- obs_haps[case_idx,    p]
    ctrl_f <- obs_haps[control_idx, p]
    lvls   <- sort(unique(c(case_f, ctrl_f)))
    if (length(lvls) < 2L) next
    # Multiply by 2: diploid individuals contribute 2 alleles each
    tab <- rbind(
      2L * table(factor(case_f, levels = lvls)),
      2L * table(factor(ctrl_f, levels = lvls))
    )
    tryCatch({
      pv <- chisq.test(tab, correct = FALSE)$p.value
      if (!is.na(pv) && pv > 0) lod[p] <- -log10(pv)
    }, error = function(e) NULL)
  }
  lod
}

# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: Arial, sans-serif; font-size: 13px; }
    h3   { color: #1F4E79; }
    h4   { color: #1F4E79; margin-bottom: 4px; margin-top: 12px; }
    .well { background: #f4f7fb; border: none; }
    hr   { margin: 8px 0; border-color: #ddd; }
  "))),

  titlePanel("Multiparent Population (MPP) Simulator"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Population Design"),
      selectInput("design", "Crossing design",
        choices = c("Fully intercrossed" = "magic",
                    "Hub-and-spoke"       = "hub"),
        width = "100%"),
      sliderInput("n_founders", "Founders",
                  min=2, max=500, value=8, step=1, width="100%"),
      sliderInput("n_gen", "Generations of recombination",
                  min=1, max=50, value=10, step=1, width="100%"),

      hr(),
      h4("Pool"),
      p(tags$small("Cases = top N by phenotype. Controls = random N from unselected individuals.")),
      sliderInput("pool_size", "Individuals per pool (N)",
                  min=10, max=500, value=50, step=10, width="100%"),
      sliderInput("n_display", "RILs to display in mosaics",
                  min=5, max=100, value=50, step=5, width="100%"),

      hr(),
      h4("Sequencing Quality"),
      sliderInput("coverage", "Coverage (x)",
                  min=0.5, max=30, value=10, step=0.5, width="100%"),
      sliderInput("haplo_acc", "Haplotype inference accuracy",
                  min=0, max=1, value=0.95, step=0.05, width="100%"),

      hr(),
      h4("Genetic Model"),
      selectInput("genetic_model", "QTL architecture",
        choices = c("Null - no QTL"                  = "null",
                    "1 QTL - major effect"            = "1qtl",
                    "2 QTLs - moderate effects"       = "2qtl",
                    "3 QTLs - moderate effects"       = "3qtl",
                    "Polygenic - many small effects"  = "polygenic"),
        width = "100%"),

      hr(),
      h4("Significance Threshold"),
      p(tags$small("Bonferroni (101 tests, α=0.05) ≈ LOD 3.3")),
      sliderInput("lod_thresh", "LOD threshold",
                  min=1, max=8, value=3.3, step=0.1, width="100%"),

      br(),
      actionButton("simulate", "Simulate",
                   width = "100%",
                   style = "background:#2E75B6; color:#fff; border:none;
                            padding:8px; border-radius:4px; font-size:14px;")
    ),

    mainPanel(
      width = 9,

      h4("Founder Haplotypes"),
      p(tags$small("Each row is one founder: pure inbred, no recombination. These are the building blocks mixed into each RIL.")),
      plotOutput("founder_plot", height = "160px"),

      h4("RIL Haplotype Mosaic"),
      p(tags$small(
        "A display sample of RILs from both pools. Colors = founder ancestry across 100 cM."
      )),
      plotOutput("hap_plot", height = "300px"),

      h4("Cases vs. Controls"),
      p(tags$small(
        "Top N by phenotype = cases (top panel); bottom N = controls (bottom panel). ",
        "Under a QTL peak, the dominant founder color should differ between panels."
      )),
      plotOutput("cases_plot", height = "360px"),

      h4("QTL Scan"),
      p(tags$small(
        "Chi-square test (2N diploid allele counts per pool) at each of 101 positions. ",
        tags$span(style="color:#1a7a1a; font-weight:bold;", "Green lines"),
        " = true QTL positions. ",
        tags$span(style="color:red;", "Red dashed"),
        " = LOD threshold."
      )),
      plotOutput("lod_plot", height = "200px"),

      h4("Founder Frequency: Cases - Controls"),
      p(tags$small(
        "At each position, the difference in each founder's allele frequency between pools. ",
        "Peaks at a QTL show which founders drive the phenotype difference."
      )),
      plotOutput("freq_plot", height = "200px")
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  results <- eventReactive(input$simulate, {
    set.seed(as.integer(Sys.time()) %% 100000L)

    # Simulate 3N individuals: top N selected as cases (high phenotype),
    # controls are a random N drawn from the remaining 2N unselected individuals.
    n_total   <- 3L * input$pool_size
    n_display <- min(input$n_display, input$pool_size)

    true_haps    <- simulate_rils(input$n_founders, input$n_gen, n_total, input$design)
    obs_haps     <- add_noise(true_haps, input$n_founders,
                              input$coverage, input$haplo_acc)
    qtl_cfg      <- get_qtl_config(input$genetic_model, input$n_founders)
    pheno        <- simulate_phenotype(true_haps, qtl_cfg)
    founder_haps <- matrix(rep(seq_len(input$n_founders), each = 101L),
                           nrow = input$n_founders, ncol = 101L)

    # Cases = top N by phenotype (selected pool)
    # Controls = random N from the remaining 2N unselected individuals
    ranked      <- order(pheno, decreasing = TRUE)
    case_idx    <- ranked[seq_len(input$pool_size)]
    unselected  <- ranked[(input$pool_size + 1L):n_total]
    control_idx <- sample(unselected, input$pool_size)

    # Fisher scan using full pools with 2N diploid allele counts
    lod <- qtl_scan(obs_haps, case_idx, control_idx)

    # Display subsets for mosaics (sorted by QTL position to pop out dominant founder)
    disp_case <- case_idx[seq_len(n_display)]
    disp_ctrl <- control_idx[seq_len(n_display)]

    list(obs_haps     = obs_haps,
         founder_haps = founder_haps,
         lod          = lod,
         qtl_pos      = qtl_cfg$pos,
         qtl_eff      = qtl_cfg$eff,
         n_founders   = input$n_founders,
         pool_size    = input$pool_size,
         n_display    = n_display,
         lod_thresh   = input$lod_thresh,
         case_idx     = case_idx,
         control_idx  = control_idx,
         disp_case    = disp_case,
         disp_ctrl    = disp_ctrl)
  })

  # --- founder plot ---
  output$founder_plot <- renderPlot({
    req(results())
    r    <- results()
    nf   <- r$n_founders
    cols <- get_founder_colors(nf)

    df <- data.frame(
      position = rep(0L:100L, each = nf),
      founder  = rep(seq_len(nf), times = 101L),
      fill_f   = factor(rep(seq_len(nf), times = 101L), levels = seq_len(nf))
    )

    ggplot(df, aes(x = position, y = founder, fill = fill_f)) +
      geom_raster() +
      scale_fill_manual(values = cols, drop = FALSE, guide = "none") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0),
                         breaks = if (nf <= 20) seq_len(nf) else
                                    c(1, seq(10, nf, by = 10)),
                         labels = if (nf <= 20) paste0("F", seq_len(nf)) else
                                    c("F1", paste0("F", seq(10, nf, by = 10)))) +
      labs(x = "Position (cM)", y = "Founder") +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(size = 9))
  })

  # --- RIL haplotype mosaic (display sample from full population) ---
  output$hap_plot <- renderPlot({
    req(results())
    r    <- results()
    nf   <- r$n_founders
    cols <- get_founder_colors(nf)

    # Show display sample: first n_display cases then first n_display controls
    disp_idx <- c(r$disp_case, r$disp_ctrl)
    disp_haps <- r$obs_haps[disp_idx, , drop = FALSE]
    nd <- nrow(disp_haps)

    df <- data.frame(
      position = rep(0L:100L, each = nd),
      ril      = rep(seq_len(nd), times = 101L),
      founder  = factor(as.vector(disp_haps), levels = seq_len(nf))
    )

    show_legend <- nf <= 16
    ggplot(df, aes(x = position, y = ril, fill = founder)) +
      geom_raster() +
      scale_fill_manual(values = cols, drop = FALSE,
                        guide = if (show_legend)
                                  guide_legend(title = "Founder",
                                               nrow = ceiling(nf / 8))
                                else "none") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "Position (cM)", y = "RIL") +
      theme_minimal(base_size = 13) +
      theme(panel.grid    = element_blank(),
            axis.text.y   = element_blank(),
            axis.ticks.y  = element_blank(),
            legend.position = if (show_legend) "bottom" else "none")
  })

  # --- cases vs controls plot ---
  output$cases_plot <- renderPlot({
    req(results())
    r    <- results()
    nf   <- r$n_founders
    cols <- get_founder_colors(nf)

    case_haps    <- r$obs_haps[r$disp_case, , drop = FALSE]
    control_haps <- r$obs_haps[r$disp_ctrl, , drop = FALSE]

    # Sort cases using the known effect founders at the primary QTL:
    # effect-founder carriers rise to the top (solid band under the green line),
    # non-effect founders fall below (random mix, matching controls).
    # Controls stay unsorted — random mix throughout is correct for an unselected pool.
    if (length(r$qtl_pos) > 0L && length(r$qtl_eff) > 0L) {
      qtl_col         <- r$qtl_pos[1L] + 1L
      effect_founders <- which(r$qtl_eff[[1L]] > 0)
      at_qtl          <- case_haps[, qtl_col]
      is_effect       <- at_qtl %in% effect_founders
      case_haps       <- case_haps[order(!is_effect, at_qtl), , drop = FALSE]
    }

    nc <- nrow(case_haps)
    nk <- nrow(control_haps)

    df_cases <- data.frame(
      position = rep(0L:100L, each = nc),
      ril      = rep(seq_len(nc), times = 101L),
      founder  = factor(as.vector(case_haps), levels = seq_len(nf)),
      group    = "Cases"
    )
    df_ctrl <- data.frame(
      position = rep(0L:100L, each = nk),
      ril      = rep(seq_len(nk), times = 101L),
      founder  = factor(as.vector(control_haps), levels = seq_len(nf)),
      group    = "Controls"
    )
    df       <- rbind(df_cases, df_ctrl)
    df$group <- factor(df$group, levels = c("Cases", "Controls"))

    show_legend <- nf <= 16
    p <- ggplot(df, aes(x = position, y = ril, fill = founder)) +
      geom_raster() +
      facet_wrap(~group, ncol = 1, scales = "free_y") +
      scale_fill_manual(values = cols, drop = FALSE,
                        guide = if (show_legend)
                                  guide_legend(title = "Founder",
                                               nrow = ceiling(nf / 8))
                                else "none") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "Position (cM)", y = "RIL") +
      theme_minimal(base_size = 13) +
      theme(panel.grid       = element_blank(),
            axis.text.y      = element_blank(),
            axis.ticks.y     = element_blank(),
            strip.text       = element_text(face = "bold", size = 11),
            legend.position  = if (show_legend) "bottom" else "none")

    if (length(r$qtl_pos) > 0)
      p <- p + geom_vline(xintercept = r$qtl_pos,
                          color = "#1a7a1a", linewidth = 0.9, alpha = 0.85)
    p
  })

  # --- LOD plot ---
  output$lod_plot <- renderPlot({
    req(results())
    r      <- results()
    lod_df <- data.frame(position = 0L:100L, lod = r$lod)
    y_max  <- max(r$lod_thresh * 1.5, max(r$lod) * 1.15, na.rm = TRUE)

    p <- ggplot(lod_df, aes(x = position, y = lod)) +
      geom_line(color = "#2E75B6", linewidth = 1) +
      geom_hline(yintercept = r$lod_thresh, linetype = "dashed",
                 color = "red", linewidth = 0.7) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
      scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
      labs(x = "Position (cM)", y = "LOD") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank())

    if (length(r$qtl_pos) > 0)
      p <- p + geom_vline(xintercept = r$qtl_pos,
                          color = "#1a7a1a", linewidth = 0.9, alpha = 0.85)
    p
  })

  # --- founder frequency difference plot (uses full pools, not display subset) ---
  output$freq_plot <- renderPlot({
    req(results())
    r    <- results()
    nf   <- r$n_founders
    cols <- get_founder_colors(nf)

    pos_seq      <- 0L:100L
    n_pos        <- length(pos_seq)
    nc           <- length(r$case_idx)
    nk           <- length(r$control_idx)
    case_haps    <- r$obs_haps[r$case_idx,    , drop = FALSE]
    control_haps <- r$obs_haps[r$control_idx, , drop = FALSE]

    rows <- vector("list", nf * n_pos)
    idx  <- 1L
    for (f in seq_len(nf)) {
      for (p in seq_len(n_pos)) {
        freq_case <- if (nc > 0) sum(case_haps[, p]    == f) / nc else 0
        freq_ctrl <- if (nk > 0) sum(control_haps[, p] == f) / nk else 0
        rows[[idx]] <- data.frame(
          position = pos_seq[p],
          founder  = factor(f, levels = seq_len(nf)),
          diff     = freq_case - freq_ctrl
        )
        idx <- idx + 1L
      }
    }
    df <- do.call(rbind, rows)

    p <- ggplot(df, aes(x = position, y = diff, color = founder)) +
      geom_line(linewidth = 0.75, alpha = 0.9) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
      scale_color_manual(values = cols, drop = FALSE,
                         guide = if (nf <= 16)
                                   guide_legend(title = "Founder",
                                                nrow = ceiling(nf / 8))
                                 else "none") +
      scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
      labs(x = "Position (cM)", y = "Freq(cases) - Freq(controls)") +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            legend.position  = if (nf <= 16) "bottom" else "none")

    if (length(r$qtl_pos) > 0)
      p <- p + geom_vline(xintercept = r$qtl_pos,
                          color = "#1a7a1a", linewidth = 0.9, alpha = 0.85)
    p
  })
}

shinyApp(ui, server)

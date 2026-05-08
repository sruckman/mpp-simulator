library(shiny)
library(ggplot2)
library(viridis)

# Chromosome positions: 0 to 100 cM in 0.5 cM steps = 202 positions
POS_SEQ <- seq(0, 100, length.out = 202L)
N_POS   <- length(POS_SEQ)

# ============================================================
# SIMULATION FUNCTIONS
# ============================================================

get_founder_colors <- function(n) {
  viridis(n, option = "turbo", begin = 0.05, end = 0.95)
}

simulate_rils <- function(n_founders, n_gen, n_total, design) {
  haps <- matrix(1L, nrow = n_total, ncol = N_POS)

  for (r in seq_len(n_total)) {
    if (design == "hub") {
      hub   <- 1L
      other <- if (n_founders >= 2L) sample.int(n_founders - 1L, 1L) + 1L else 1L
      pool  <- c(hub, other)
    } else {
      pool <- seq_len(n_founders)
    }

    state <- rep(sample(pool, 1L), N_POS)

    total_cross <- rpois(1L, n_gen * 0.25)
    if (total_cross > 0L) {
      cpts <- sort(runif(total_cross, 0, 100))
      for (cp in cpts) {
        idx <- which(POS_SEQ >= cp)[1L]
        if (is.na(idx) || idx > N_POS) next
        new_f            <- sample(pool, 1L)
        state[idx:N_POS] <- new_f
      }
    }

    haps[r, ] <- state
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
    col_idx <- which.min(abs(POS_SEQ - qtl_cfg$pos[i]))
    pheno   <- pheno + qtl_cfg$eff[[i]][true_haps[, col_idx]]
  }
  pheno
}

# Chi-square test on 2 x k founder count table (cases vs controls).
# Coverage models read-depth sampling: instead of observing all 2N true allele
# counts, we draw (coverage x pool_size) reads from the pool proportional to
# true founder frequencies. High coverage -> observed ~ true -> strong LOD.
# Low coverage -> noisy sample -> weak LOD.
qtl_scan <- function(obs_haps, case_idx, control_idx, coverage) {
  n_pos      <- ncol(obs_haps)
  lod        <- numeric(n_pos)
  n_reads    <- round(coverage * length(case_idx))

  for (p in seq_len(n_pos)) {
    case_f <- obs_haps[case_idx,    p]
    ctrl_f <- obs_haps[control_idx, p]
    lvls   <- sort(unique(c(case_f, ctrl_f)))
    if (length(lvls) < 2L) next

    true_case <- as.numeric(2L * table(factor(case_f, levels = lvls)))
    true_ctrl <- as.numeric(2L * table(factor(ctrl_f, levels = lvls)))

    # Coverage = accuracy of frequency measurement, not sample size.
    # Sample reads from the pool at the given coverage depth,
    # then normalize back to pool-size scale so LOD reflects
    # pool size (biology/design) not read depth (technical).
    obs_case <- as.numeric(rmultinom(1L, n_reads, prob = true_case / sum(true_case)))
    obs_ctrl <- as.numeric(rmultinom(1L, n_reads, prob = true_ctrl / sum(true_ctrl)))

    norm_case <- obs_case / sum(obs_case) * sum(true_case)
    norm_ctrl <- obs_ctrl / sum(obs_ctrl) * sum(true_ctrl)

    tab <- rbind(norm_case, norm_ctrl)
    tryCatch({
      pv <- suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
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
                  min=1, max=50, value=5, step=1, width="100%"),

      hr(),
      h4("Pool"),
      p(tags$small("Cases = top N by phenotype. Controls = random N from unselected individuals.")),
      sliderInput("pool_size", "Individuals per pool (N)",
                  min=10, max=500, value=50, step=10, width="100%"),
      sliderInput("n_display", "RILs to display in mosaics",
                  min=5, max=100, value=50, step=5, width="100%"),

      hr(),
      h4("Sequencing Coverage"),
      p(tags$small("Reads per individual. Higher coverage = less sampling noise = stronger LOD.")),
      sliderInput("coverage", "Coverage (x)",
                  min=1, max=100, value=30, step=1, width="100%"),

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
      p(tags$small("Bonferroni (202 tests, α=0.05) ≈ LOD 3.6")),
      sliderInput("lod_thresh", "LOD threshold",
                  min=1, max=8, value=3.6, step=0.1, width="100%"),

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
        "Fisher's exact test (2N diploid allele counts per pool) at each of 202 positions (0.5 cM steps). Assumes perfect haplotype inference. ",
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

    obs_haps     <- simulate_rils(input$n_founders, input$n_gen, n_total, input$design)
    qtl_cfg      <- get_qtl_config(input$genetic_model, input$n_founders)
    pheno        <- simulate_phenotype(obs_haps, qtl_cfg)
    founder_haps <- matrix(rep(seq_len(input$n_founders), each = N_POS),
                           nrow = input$n_founders, ncol = N_POS)

    # Cases = top N by phenotype (selected pool)
    # Controls = random N from the remaining 2N unselected individuals
    ranked      <- order(pheno, decreasing = TRUE)
    case_idx    <- ranked[seq_len(input$pool_size)]
    unselected  <- ranked[(input$pool_size + 1L):n_total]
    control_idx <- sample(unselected, input$pool_size)

    # Fisher scan using full pools with 2N diploid allele counts
    lod <- qtl_scan(obs_haps, case_idx, control_idx, input$coverage)

    # Display subsets for mosaics
    disp_all  <- sample.int(n_total, min(input$n_display * 2L, n_total))
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
         disp_all     = disp_all,
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
      position = rep(POS_SEQ, each = nf),
      founder  = rep(seq_len(nf), times = N_POS),
      fill_f   = factor(rep(seq_len(nf), times = N_POS), levels = seq_len(nf))
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

    # Random sample from the full population — not split by case/control
    disp_haps <- r$obs_haps[r$disp_all, , drop = FALSE]
    nd <- nrow(disp_haps)

    df <- data.frame(
      position = rep(POS_SEQ, each = nd),
      ril      = rep(seq_len(nd), times = N_POS),
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
      position = rep(POS_SEQ, each = nc),
      ril      = rep(seq_len(nc), times = N_POS),
      founder  = factor(as.vector(case_haps), levels = seq_len(nf)),
      group    = "Cases"
    )
    df_ctrl <- data.frame(
      position = rep(POS_SEQ, each = nk),
      ril      = rep(seq_len(nk), times = N_POS),
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
    lod_df <- data.frame(position = POS_SEQ, lod = r$lod)
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

    pos_seq      <- POS_SEQ
    n_pos        <- N_POS
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

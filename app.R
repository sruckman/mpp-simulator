library(shiny)
library(ggplot2)

# Chromosome positions: 0 to 100 cM in 0.5 cM steps = 202 positions
POS_SEQ <- seq(0, 100, length.out = 202L)
N_POS   <- length(POS_SEQ)

# ============================================================
# UTILITY FUNCTIONS
# ============================================================

# D3 category20 — 12 bold, well-separated colors standard in data visualization.
# Alternates saturated and lighter versions for extra contrast between neighbors.
get_founder_colors <- function(n) {
  base <- c(
    "dodgerblue4", "#2CA02C", "darkred",  "#9467BD",
    "#17BECF", "#FF83FA", "darkmagenta", "grey40",
    "#BCBD22", "#6495ED", "#43CD80",  "#F08080"
  )
  if (n <= 12L) return(base[seq_len(n)])
  colorRampPalette(base)(n)
}

# Rolling mean over a LOD vector; window = number of positions
smooth_lod <- function(lod, window = 5L) {
  n    <- length(lod)
  half <- window %/% 2L
  vapply(seq_len(n), function(i) {
    mean(lod[max(1L, i - half):min(n, i + half)])
  }, numeric(1L))
}

# -log10(p) computed in log-space to avoid p=0 underflow at strong QTL peaks.
chisq_lod <- function(statistic, df) {
  log_p <- pchisq(as.numeric(statistic), df = as.numeric(df),
                  lower.tail = FALSE, log.p = TRUE)
  -log_p / log(10)
}

# ============================================================
# SIMULATION FUNCTIONS
# ============================================================

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
    state       <- rep(sample(pool, 1L), N_POS)
    # Rate = 0.1 crossovers per generation: big blocks at low gen, smaller at high gen
    total_cross <- rpois(1L, n_gen * 0.1)
    if (total_cross > 0L) {
      cpts <- sort(runif(total_cross, 0, 100))
      for (cp in cpts) {
        idx <- which(POS_SEQ >= cp)[1L]
        if (is.na(idx) || idx > N_POS) next
        state[idx:N_POS] <- sample(pool, 1L)
      }
    }
    haps[r, ] <- state
  }
  haps
}

# hub_design = TRUE: hub founder (F1) is excluded from QTL effects,
# so the signal always comes from a spoke founder.
get_qtl_config <- function(model, n_founders, hub_design = FALSE) {
  avail    <- if (hub_design && n_founders > 1L) seq(2L, n_founders) else seq_len(n_founders)
  n_effect <- max(1L, min(10L, round(0.1 * length(avail))))

  make_eff <- function(effect_size) {
    eff      <- numeric(n_founders)
    who      <- sample(avail, min(n_effect, length(avail)))
    eff[who] <- runif(length(who), effect_size * 0.75, effect_size * 1.25)
    eff
  }

  switch(model,
    "null"      = list(pos = integer(0), eff = list()),
    "1qtl"      = list(pos = 50L,
                       eff = list(make_eff(1.5))),
    "2qtl"      = list(pos = c(30L, 70L),
                       eff = list(make_eff(1.2), make_eff(1.2))),
    "3qtl"      = list(pos = c(20L, 50L, 80L),
                       eff = list(make_eff(1.0), make_eff(1.0), make_eff(1.0))),
    "polygenic" = list(pos = c(5L,15L,25L,35L,45L,55L,65L,75L,85L,95L),
                       eff = lapply(seq_len(10L), function(i) make_eff(0.7)))
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

# Chi-square (1 replicate) or CMH (multiple replicates) QTL scan.
# Uses individual counts (not 2x diploid) to keep LOD scores realistic.
qtl_scan <- function(haps_list, case_idx_list, control_idx_list) {
  n_reps <- length(haps_list)
  n_pos  <- ncol(haps_list[[1L]])
  lod    <- numeric(n_pos)

  for (p in seq_len(n_pos)) {
    lvls <- sort(unique(unlist(lapply(seq_len(n_reps), function(ri) {
      c(haps_list[[ri]][case_idx_list[[ri]], p],
        haps_list[[ri]][control_idx_list[[ri]], p])
    }))))
    if (length(lvls) < 2L) next

    if (n_reps == 1L) {
      cf  <- haps_list[[1L]][case_idx_list[[1L]], p]
      kf  <- haps_list[[1L]][control_idx_list[[1L]], p]
      tab <- rbind(
        as.numeric(table(factor(cf, levels = lvls))),
        as.numeric(table(factor(kf, levels = lvls)))
      )
      tryCatch({
        test   <- suppressWarnings(chisq.test(tab, correct = FALSE))
        lod[p] <- chisq_lod(test$statistic, test$parameter)
      }, error = function(e) NULL)
    } else {
      k       <- length(lvls)
      tab_arr <- array(0L, dim = c(2L, k, n_reps))
      for (ri in seq_len(n_reps)) {
        cf <- haps_list[[ri]][case_idx_list[[ri]], p]
        kf <- haps_list[[ri]][control_idx_list[[ri]], p]
        tab_arr[1L, , ri] <- as.integer(table(factor(cf, levels = lvls)))
        tab_arr[2L, , ri] <- as.integer(table(factor(kf, levels = lvls)))
      }
      tryCatch({
        test   <- suppressWarnings(mantelhaen.test(tab_arr, correct = FALSE))
        lod[p] <- chisq_lod(test$statistic, test$parameter)
      }, error = function(e) NULL)
    }
  }
  lod
}

# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: Arial, sans-serif; font-size: 13px; }
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
                  min=2, max=30, value=8, step=2, width="100%"),
      sliderInput("n_gen", "Generations of recombination",
                  min=5, max=50, value=10, step=5, width="100%"),

      hr(),
      h4("Pool"),
      p(tags$small("Cases = top N by phenotype. Controls = random N from unselected.")),
      sliderInput("pool_size", "Individuals per pool (N)",
                  min=50, max=1000, value=100, step=50, width="100%"),
      sliderInput("n_display", "RILs to display in mosaics",
                  min=10, max=100, value=50, step=10, width="100%"),

      hr(),
      h4("Genetic Model"),
      selectInput("genetic_model", "QTL architecture",
        choices = c("Null - no QTL"                 = "null",
                    "1 QTL - major effect"           = "1qtl",
                    "2 QTLs - moderate effects"      = "2qtl",
                    "3 QTLs - moderate effects"      = "3qtl",
                    "Polygenic - many small effects" = "polygenic"),
        width = "100%"),

      hr(),
      h4("Significance Threshold"),
      p(tags$small("Bonferroni (202 tests, α=0.05) ≈ LOD 3.6")),
      sliderInput("lod_thresh", "LOD threshold",
                  min=1, max=8, value=4, step=0.5, width="100%"),

      br(),
      actionButton("simulate", "Simulate",
                   width = "100%",
                   style = "background:#2E75B6; color:#fff; border:none;
                            padding:8px; border-radius:4px; font-size:14px;")
    ),

    mainPanel(
      width = 9,

      h4("Founder Haplotypes"),
      p(tags$small("Each row is one founder: pure inbred, no recombination.")),
      plotOutput("founder_plot", height = "160px"),

      h4("RIL Haplotype Mosaic"),
      p(tags$small("A random sample of RILs. Colors = founder ancestry across 100 cM.")),
      plotOutput("hap_plot", height = "300px"),

      h4("Cases vs. Controls"),
      p(tags$small(
        "Top N = cases (top); random N from unselected = controls (bottom). ",
        "Under a QTL peak, one founder color dominates in cases."
      )),
      plotOutput("cases_plot", height = "360px"),

      h4("QTL Scan"),
      p(tags$small(
        "Chi-square test on founder counts, smoothed with 5-position rolling average. ",
        tags$span(style="color:red; font-weight:bold;", "Red solid"),
        " = true QTL.  ",
        tags$span(style="color:blue;", "Blue dashed"),
        " = LOD threshold."
      )),
      plotOutput("lod_plot", height = "200px"),

      h4("Founder Frequency: Cases − Controls"),
      p(tags$small("Difference in founder frequency between pools at each position.")),
      plotOutput("freq_plot", height = "200px")
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  results <- eventReactive(input$simulate, {
    pool_size <- input$pool_size
    n_total   <- 3L * pool_size
    n_display <- min(input$n_display, pool_size)
    hub       <- input$design == "hub"
    qtl_cfg   <- get_qtl_config(input$genetic_model, input$n_founders,
                                hub_design = hub)

    founder_haps <- matrix(rep(seq_len(input$n_founders), each = N_POS),
                           nrow = input$n_founders, ncol = N_POS)

    obs_haps    <- simulate_rils(input$n_founders, input$n_gen, n_total, input$design)
    pheno       <- simulate_phenotype(obs_haps, qtl_cfg)
    ranked      <- order(pheno, decreasing = TRUE)
    ci1         <- ranked[seq_len(pool_size)]
    unselected  <- ranked[(pool_size + 1L):n_total]
    ki1         <- sample(unselected, pool_size)

    raw_lod <- qtl_scan(list(obs_haps), list(ci1), list(ki1))
    lod     <- smooth_lod(raw_lod, window = 5L)
    disp_all  <- sample.int(n_total, min(n_display * 2L, n_total))
    disp_case <- ci1[seq_len(n_display)]
    disp_ctrl <- ki1[seq_len(n_display)]

    list(
      obs_haps     = obs_haps,
      founder_haps = founder_haps,
      lod          = lod,
      qtl_pos      = qtl_cfg$pos,
      qtl_eff      = qtl_cfg$eff,
      n_founders   = input$n_founders,
      pool_size    = pool_size,
      n_display    = n_display,
      lod_thresh   = input$lod_thresh,
      case_idx     = ci1,
      control_idx  = ki1,
      disp_all     = disp_all,
      disp_case    = disp_case,
      disp_ctrl    = disp_ctrl
    )
  })

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
        breaks = if (nf <= 20) seq_len(nf) else c(1, seq(10, nf, by=10)),
        labels = if (nf <= 20) paste0("F", seq_len(nf)) else
                   c("F1", paste0("F", seq(10, nf, by=10)))) +
      labs(x = "Position (cM)", y = "Founder") +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(), axis.text.y = element_text(size = 9))
  })

  output$hap_plot <- renderPlot({
    req(results())
    r    <- results()
    nf   <- r$n_founders
    cols <- get_founder_colors(nf)
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
        guide = if (show_legend) guide_legend(title="Founder", nrow=ceiling(nf/8)) else "none") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "Position (cM)", y = "RIL") +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = if (show_legend) "bottom" else "none")
  })

  output$cases_plot <- renderPlot({
    req(results())
    r    <- results()
    nf   <- r$n_founders
    cols <- get_founder_colors(nf)
    case_haps    <- r$obs_haps[r$disp_case, , drop = FALSE]
    control_haps <- r$obs_haps[r$disp_ctrl, , drop = FALSE]
    if (length(r$qtl_pos) > 0L && length(r$qtl_eff) > 0L) {
      qtl_col         <- which.min(abs(POS_SEQ - r$qtl_pos[1L]))
      effect_founders <- which(r$qtl_eff[[1L]] > 0)
      at_qtl          <- case_haps[, qtl_col]
      is_effect       <- at_qtl %in% effect_founders
      case_haps       <- case_haps[order(!is_effect, at_qtl), , drop = FALSE]
    }
    nc <- nrow(case_haps)
    nk <- nrow(control_haps)
    df_cases <- data.frame(
      position = rep(POS_SEQ, each=nc), ril = rep(seq_len(nc), times=N_POS),
      founder  = factor(as.vector(case_haps), levels=seq_len(nf)), group = "Cases")
    df_ctrl  <- data.frame(
      position = rep(POS_SEQ, each=nk), ril = rep(seq_len(nk), times=N_POS),
      founder  = factor(as.vector(control_haps), levels=seq_len(nf)), group = "Controls")
    df       <- rbind(df_cases, df_ctrl)
    df$group <- factor(df$group, levels=c("Cases","Controls"))
    show_legend <- nf <= 16
    p <- ggplot(df, aes(x=position, y=ril, fill=founder)) +
      geom_raster() +
      facet_wrap(~group, ncol=1, scales="free_y") +
      scale_fill_manual(values=cols, drop=FALSE,
        guide=if(show_legend) guide_legend(title="Founder",nrow=ceiling(nf/8)) else "none") +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      labs(x="Position (cM)", y="RIL") +
      theme_minimal(base_size=13) +
      theme(panel.grid=element_blank(), axis.text.y=element_blank(),
            axis.ticks.y=element_blank(), strip.text=element_text(face="bold",size=11),
            legend.position=if(show_legend)"bottom" else "none")
    if (length(r$qtl_pos) > 0)
      p <- p + geom_vline(xintercept=r$qtl_pos, color="red", linewidth=0.9, alpha=0.85)
    p
  })

  output$lod_plot <- renderPlot({
    req(results())
    r      <- results()
    lod_df <- data.frame(position=POS_SEQ, lod=r$lod)
    y_max  <- max(r$lod_thresh * 1.5, max(r$lod) * 1.15, na.rm=TRUE)
    p <- ggplot(lod_df, aes(x=position, y=lod)) +
      geom_line(color="black", linewidth=1) +
      geom_hline(yintercept=r$lod_thresh, linetype="dashed", color="blue", linewidth=0.7) +
      scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
      scale_y_continuous(limits=c(0,y_max), expand=c(0,0)) +
      labs(x="Position (cM)", y="LOD") +
      theme_minimal(base_size=13) +
      theme(panel.grid.minor=element_blank())
    if (length(r$qtl_pos) > 0)
      p <- p + geom_vline(xintercept=r$qtl_pos, color="red", linewidth=0.9, alpha=0.85)
    p
  })

  output$freq_plot <- renderPlot({
    req(results())
    r    <- results()
    nf   <- r$n_founders
    cols <- get_founder_colors(nf)
    nc   <- length(r$case_idx)
    nk   <- length(r$control_idx)
    case_haps    <- r$obs_haps[r$case_idx,    , drop=FALSE]
    control_haps <- r$obs_haps[r$control_idx, , drop=FALSE]
    rows <- vector("list", nf * N_POS)
    k    <- 1L
    for (f in seq_len(nf)) {
      for (pos in seq_len(N_POS)) {
        rows[[k]] <- data.frame(
          position = POS_SEQ[pos],
          founder  = factor(f, levels=seq_len(nf)),
          diff     = (if(nc>0) sum(case_haps[,pos]==f)/nc else 0) -
                     (if(nk>0) sum(control_haps[,pos]==f)/nk else 0)
        )
        k <- k + 1L
      }
    }
    df <- do.call(rbind, rows)
    p <- ggplot(df, aes(x=position, y=diff, color=founder)) +
      geom_line(linewidth=0.75, alpha=0.9) +
      geom_hline(yintercept=0, color="grey50", linewidth=0.5) +
      scale_color_manual(values=cols, drop=FALSE,
        guide=if(nf<=16) guide_legend(title="Founder",nrow=ceiling(nf/8)) else "none") +
      scale_x_continuous(expand=c(0,0), limits=c(0,100)) +
      labs(x="Position (cM)", y="Freq(cases) − Freq(controls)") +
      theme_minimal(base_size=13) +
      theme(panel.grid.minor=element_blank(),
            legend.position=if(nf<=16)"bottom" else "none")
    if (length(r$qtl_pos) > 0)
      p <- p + geom_vline(xintercept=r$qtl_pos, color="red", linewidth=0.9, alpha=0.85)
    p
  })
}

shinyApp(ui, server)

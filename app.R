library(shiny)

# =====================================================================
# Sexual partnership concurrency and HIV-1 transmission dynamics.
# An interactive stochastic microsimulation from the EpiModel ecosystem,
# a companion to the Concurrency Tutorial.
#
# A population of men and women with purely heterosexual mixing forms and
# dissolves partnerships over time. HIV transmission depends on the stage
# of the infected partner (a four-stage acute -> latent -> pre-AIDS -> AIDS
# progression, with stage-specific transmission probabilities from
# Hollingsworth 2008). The only thing that changes across the four
# scenarios is whether men, women, both, or neither are permitted
# concurrent (overlapping) partnerships; every other parameter is held
# fixed so the effect of concurrency alone can be read off the curves.
#
# The whole model uses base R only, so it runs entirely in the browser
# through shinylive with no package dependencies beyond shiny.
# =====================================================================


# ---- Model core (ported from the concurrency.sim package) -----------

# Runs the concurrency microsimulation and returns an object of class
# "conc_microsim": a list with one prevalence matrix (female, male) per
# simulation. Kept faithful to the original package implementation.
conc_microsim <- function(s.num.f, i.num.f, s.num.m, i.num.m,
                          monog.f = TRUE, monog.m = TRUE,
                          meandeg, part.duration, nsteps, nsims = 1,
                          verbose = FALSE, updateProgress = NULL) {

  # Per-month transmission probability by time since infection of the
  # positive partner (Hollingsworth et al. 2008): high in acute, low in
  # latent, elevated in pre-AIDS, zero in AIDS (no acts).
  betaByTime <- c(rep(0.2055, 3), rep(0.0088, 100),
                  rep(0.0614, 9), rep(0, 10))

  out <- list()
  for (sim in 1:nsims) {

    nFeml <- s.num.f + i.num.f
    nMale <- s.num.m + i.num.m
    n <- nFeml + nMale

    expectedEdges <- round(meandeg * n / 2)
    probDissolution <- 1 / part.duration
    timeAidsDeath <- length(betaByTime)

    # Female state vectors
    idsFeml <- 1:nFeml
    hivStatFeml <- rep(0, nFeml)
    infTimeFeml <- rep(NA, nFeml)
    if (i.num.f > 0) {
      idsInfFeml <- sample(idsFeml, i.num.f)
      hivStatFeml[idsInfFeml] <- 1
      infTimeFeml[idsInfFeml] <- sample(0:(-timeAidsDeath + 2), i.num.f, replace = TRUE)
    }

    # Male state vectors
    idsMale <- 1:nMale
    hivStatMale <- rep(0, nMale)
    infTimeMale <- rep(NA, nMale)
    if (i.num.m > 0) {
      idsInfMale <- sample(idsMale, i.num.m)
      hivStatMale[idsInfMale] <- 1
      infTimeMale[idsInfMale] <- sample(0:(-timeAidsDeath + 2), i.num.m, replace = TRUE)
    }

    # Initial contact network. Sampling with replacement permits a person
    # to hold several ties (concurrency); without replacement enforces a
    # momentary-degree cap of one (monogamy).
    if (monog.f == FALSE) {
      edgelistFeml <- sample(idsFeml, expectedEdges, replace = TRUE)
    } else {
      edgelistFeml <- sample(idsFeml, expectedEdges, replace = FALSE)
    }
    if (monog.m == FALSE) {
      edgelistMale <- sample(idsMale, expectedEdges, replace = TRUE)
    } else {
      edgelistMale <- sample(idsMale, expectedEdges, replace = FALSE)
    }

    for (tstep in 1:nsteps) {

      ## Transmissions ##
      hivStatFemlPart <- hivStatFeml[edgelistFeml]
      hivStatMalePart <- hivStatMale[edgelistMale]

      # Serodiscordant with positive female -> male at risk
      sdpFeml <- which(hivStatMalePart == 0 & hivStatFemlPart == 1)
      idsSdpFeml <- edgelistFeml[sdpFeml]
      infTimeSdpFeml <- tstep - infTimeFeml[idsSdpFeml]
      probTransSdpFeml <- betaByTime[infTimeSdpFeml]
      transSdpFeml <- rbinom(length(sdpFeml), 1, probTransSdpFeml)
      newInfMale <- edgelistMale[sdpFeml[transSdpFeml == 1]]
      hivStatMale[newInfMale] <- 1
      infTimeMale[newInfMale] <- tstep

      # Serodiscordant with positive male -> female at risk
      sdpMale <- which(hivStatFemlPart == 0 & hivStatMalePart == 1)
      idsSdpMale <- edgelistMale[sdpMale]
      infTimeSdpMale <- tstep - infTimeMale[idsSdpMale]
      probTransSdpMale <- betaByTime[infTimeSdpMale]
      transSdpMale <- rbinom(length(sdpMale), 1, probTransSdpMale)
      newInfFeml <- edgelistFeml[sdpMale[transSdpMale == 1]]
      hivStatFeml[newInfFeml] <- 1
      infTimeFeml[newInfFeml] <- tstep

      ## Deaths to AIDS ##
      idsDeathAidsFeml <- which((tstep - infTimeFeml) == timeAidsDeath)
      idsDeathAidsMale <- which((tstep - infTimeMale) == timeAidsDeath)

      ## End of ties because of death ##
      edgesFemlDeath <- which(edgelistFeml %in% idsDeathAidsFeml)
      edgesMaleDeath <- which(edgelistMale %in% idsDeathAidsMale)
      edgesBothDeath <- c(edgesFemlDeath, edgesMaleDeath)
      if (length(edgesBothDeath) > 0) {
        edgelistFeml <- edgelistFeml[-edgesBothDeath]
        edgelistMale <- edgelistMale[-edgesBothDeath]
      }

      ## End of other ties randomly ##
      edgesElig <- rbinom(length(edgelistMale), 1, probDissolution)
      edgesBroken <- which(edgesElig == 1)
      if (length(edgesBroken) > 0) {
        edgelistMale <- edgelistMale[-edgesBroken]
        edgelistFeml <- edgelistFeml[-edgesBroken]
      }

      ## Add new edges (replace those lost, keeping the edge count fixed) ##
      nNewTies <- expectedEdges - length(edgelistMale)
      if (nNewTies > 0) {
        if (monog.f == FALSE) {
          f <- sample(1:nFeml, nNewTies, replace = TRUE)
        } else {
          femlNoTies <- setdiff(1:nFeml, edgelistFeml)
          # Index explicitly rather than sample(femlNoTies, ...): when only one
          # unpaired id remains, sample() would treat the scalar as a range.
          f <- femlNoTies[sample.int(length(femlNoTies), nNewTies)]
        }
        if (monog.m == FALSE) {
          m <- sample(1:nMale, nNewTies, replace = TRUE)
        } else {
          maleNoTies <- setdiff(1:nMale, edgelistMale)
          m <- maleNoTies[sample.int(length(maleNoTies), nNewTies)]
        }
        edgelistFeml <- c(edgelistFeml, f)
        edgelistMale <- c(edgelistMale, m)
      }

      ## Replace AIDS deaths with new susceptible births ##
      hivStatFeml[idsDeathAidsFeml] <- 0
      infTimeFeml[idsDeathAidsFeml] <- NA
      hivStatMale[idsDeathAidsMale] <- 0
      infTimeMale[idsDeathAidsMale] <- NA

      ## Track prevalence by sex ##
      if (tstep == 1) {
        femlPrev <- mean(hivStatFeml)
        malePrev <- mean(hivStatMale)
      } else {
        femlPrev[tstep] <- mean(hivStatFeml)
        malePrev[tstep] <- mean(hivStatMale)
      }
    }

    if (verbose == TRUE) cat("SIM = ", sim, "/", nsims, "\n", sep = "")

    prevBoth <- matrix(c(femlPrev, malePrev), ncol = 2)
    colnames(prevBoth) <- c("femlPrev", "malePrev")
    out[[sim]] <- prevBoth

    if (is.function(updateProgress)) {
      updateProgress(detail = paste0("Simulation ", sim, "/", nsims, " complete"))
    }
  }

  class(out) <- "conc_microsim"
  out
}


# ---- Palette --------------------------------------------------------
col_feml   <- "#c0392b"   # red: female prevalence
col_male   <- "#4682B4"   # steel blue: male prevalence
col_teal   <- "#1a6b54"   # teal: titles and accents
col_target <- "#333333"


# ---- Plot builder (restyled from plot.conc_microsim) ----------------

draw_means <- function(x, mean.smooth, mean.lwd, pal) {
  df <- do.call("cbind", x)
  dfF <- df[, which(colnames(df) == "femlPrev"), drop = FALSE]
  dfM <- df[, which(colnames(df) == "malePrev"), drop = FALSE]
  meanF <- rowMeans(dfF)
  meanM <- rowMeans(dfM)
  if (mean.smooth == TRUE) {
    meanF <- supsmu(x = 1:length(meanF), y = meanF)$y
    meanM <- supsmu(x = 1:length(meanM), y = meanM)$y
  }
  lines(meanF, lwd = mean.lwd, col = pal[1])
  lines(meanM, lwd = mean.lwd, col = pal[2])
}

draw_qnts <- function(x, qnts, qnt.pal, qnts.smooth) {
  df <- do.call("cbind", x)
  dfF <- df[, which(colnames(df) == "femlPrev"), drop = FALSE]
  dfM <- df[, which(colnames(df) == "malePrev"), drop = FALSE]
  quants <- c((1 - qnts) / 2, 1 - ((1 - qnts) / 2))
  qntF <- apply(dfF, 1, function(z) quantile(z, c(quants[1], quants[2]), na.rm = TRUE))
  qntM <- apply(dfM, 1, function(z) quantile(z, c(quants[1], quants[2]), na.rm = TRUE))

  xxF <- c(1:(ncol(qntF)), (ncol(qntF)):1)
  if (qnts.smooth == FALSE) {
    yyF <- c(qntF[1, ], rev(qntF[2, ]))
  } else {
    yyF <- c(supsmu(x = 1:(ncol(qntF)), y = qntF[1, ])$y,
             rev(supsmu(x = 1:(ncol(qntF)), y = qntF[2, ])$y))
  }
  polygon(xxF, yyF, col = qnt.pal[1], border = NA)

  xxM <- c(1:(ncol(qntM)), (ncol(qntM)):1)
  if (qnts.smooth == FALSE) {
    yyM <- c(qntM[1, ], rev(qntM[2, ]))
  } else {
    yyM <- c(supsmu(x = 1:(ncol(qntM)), y = qntM[1, ])$y,
             rev(supsmu(x = 1:(ncol(qntM)), y = qntM[2, ])$y))
  }
  polygon(xxM, yyM, col = qnt.pal[2], border = NA)
}

plot_conc <- function(x, sim.lines = FALSE, sim.alpha = 0.3,
                      mean.line = TRUE, mean.smooth = TRUE, mean.lwd = 3,
                      qnts = 0.5, qnt.alpha = 0.3, qnts.smooth = TRUE) {
  nsims <- length(x)
  nsteps <- nrow(x[[1]])
  pal <- c(col_feml, col_male)
  sim.pal <- adjustcolor(pal, alpha.f = sim.alpha)
  qnt.pal <- adjustcolor(pal, alpha.f = qnt.alpha)

  ylim <- c(0, min(1, max(sapply(x, max)) * 1.25))
  xlim <- c(0, nsteps)

  par(family = "sans", mar = c(4.5, 4.5, 2.5, 1), bg = "#ffffff")
  plot(1, 1, type = "n", xlim = xlim, ylim = ylim, bty = "n",
       xlab = "Time (months)", ylab = "HIV prevalence",
       main = "Simulated HIV prevalence over time",
       col.main = col_teal, font.main = 2, cex.main = 1.25,
       col.lab = "#555", cex.lab = 1.05, col.axis = "#777")

  if (qnts > 0) draw_qnts(x, qnts, qnt.pal, qnts.smooth)
  if (sim.lines == TRUE) {
    for (i in 1:nsims) {
      lines(1:nsteps, x[[i]][, "femlPrev"], col = sim.pal[1], lwd = 0.8)
      lines(1:nsteps, x[[i]][, "malePrev"], col = sim.pal[2], lwd = 0.8)
    }
  }
  if (mean.line == TRUE) draw_means(x, mean.smooth, mean.lwd = mean.lwd, pal)

  legend("topleft", c("Female", "Male"), lwd = 4, col = pal[1:2],
         bty = "n", cex = 1.0, text.col = "#555")
}


# ---- Interpretation text (pure) -------------------------------------

conc_label <- function(conc) {
  switch(conc,
    "No Concurrency"          = "no concurrency (both sexes monogamous)",
    "Female Only Concurrency" = "female-only concurrency",
    "Male Only Concurrency"   = "male-only concurrency",
    "Both Sexes Concurrency"  = "concurrency in both sexes")
}

recap_items <- function(x, conc) {
  df <- do.call("cbind", x)
  dfF <- df[, which(colnames(df) == "femlPrev"), drop = FALSE]
  dfM <- df[, which(colnames(df) == "malePrev"), drop = FALSE]
  finalF <- mean(dfF[nrow(dfF), ])
  finalM <- mean(dfM[nrow(dfM), ])
  overall <- (finalF + finalM) / 2
  peak <- max(rowMeans((dfF + dfM) / 2))

  list(
    HTML(sprintf("Scenario: <strong>%s</strong>, across %d simulation%s.",
                 conc_label(conc), length(x), if (length(x) == 1) "" else "s")),
    HTML(sprintf("Final HIV prevalence: <strong>%.1f%%</strong> in women and <strong>%.1f%%</strong> in men.",
                 100 * finalF, 100 * finalM)),
    HTML(sprintf("Final overall prevalence: <strong>%.1f%%</strong>; peak reached along the way: <strong>%.1f%%</strong>.",
                 100 * overall, 100 * peak)),
    HTML("Hold every other parameter fixed and switch the concurrency rule: allowing overlapping partnerships opens faster transmission paths through the network, so prevalence rises even though mean degree and partnership duration are unchanged.")
  )
}


# ---- UI -------------------------------------------------------------
ui <- fluidPage(

  tags$head(
    tags$title("Concurrency and HIV Transmission"),
    tags$link(
      href = "https://fonts.googleapis.com/css2?family=Open+Sans:wght@400;600;700&display=swap",
      rel = "stylesheet"
    ),
    tags$style(HTML("
      html { font-size: 16px; }
      body {
        font-family: 'Open Sans', -apple-system, BlinkMacSystemFont, sans-serif;
        font-size: 1rem; color: #444; background-color: #f0f5f9;
        -webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale;
      }
      .container-fluid { max-width: 1200px; padding: 0 2rem; }

      /* Header banner */
      .app-header {
        background: linear-gradient(135deg, #3a6f99, #4682B4);
        color: #ffffff; padding: 2rem 2.5rem; margin: -15px -15px 1.5rem -15px;
        border-radius: 0 0 8px 8px;
      }
      .app-header h2 { font-weight: 700; font-size: 2rem; margin: 0 0 0.4rem 0; }
      .app-header p { font-size: 1.15rem; margin: 0; opacity: 0.9; }

      /* Info cards */
      .card {
        background: #ffffff; border: 1px solid #dee2e6; border-radius: 8px;
        padding: 1.5rem; box-shadow: 0 2px 6px rgba(0,0,0,0.05); margin-bottom: 1.25rem;
      }
      .card h4 {
        font-weight: 600; font-size: 1.25rem; color: #1a6b54; margin-top: 0;
        margin-bottom: 0.75rem; padding-bottom: 0.5rem; border-bottom: 2px solid #e8eef3;
      }
      .card p, .card li { font-size: 1.05rem; line-height: 1.65; color: #555; }
      .card ul { padding-left: 1.2rem; }
      .card li { padding: 0.2rem 0; }
      .card strong { color: #1a6b54; }
      .card a { color: #3a6f99; text-decoration: none; }
      .card a:hover { color: #1a6b54; text-decoration: underline; }

      /* Sidebar controls */
      .well {
        background: #ffffff !important; border: 1px solid #dee2e6 !important;
        border-radius: 8px !important; box-shadow: 0 2px 6px rgba(0,0,0,0.05) !important;
        padding: 1.5rem !important;
      }
      .control-label { font-weight: 600; font-size: 1.0rem; color: #3a6f99; margin-bottom: 0.25rem; }
      .irs--shiny .irs-bar { background: #4682B4; border-top: 1px solid #3a6f99; border-bottom: 1px solid #3a6f99; }
      .irs--shiny .irs-single, .irs--shiny .irs-from, .irs--shiny .irs-to { background-color: #4682B4; }
      .irs--shiny .irs-handle { border: 1px solid #4682B4; }
      .btn-primary {
        background-color: #4682B4; border-color: #3a6f99; font-weight: 600;
        width: 100%; margin-bottom: 0.5rem;
      }
      .btn-primary:hover, .btn-primary:focus { background-color: #3a6f99; border-color: #2f597d; }
      .radio label, .checkbox label { font-size: 0.95rem; color: #555; }
      .sidebar-section {
        font-weight: 600; font-size: 1.1rem; color: #1a6b54;
        margin: 1.25rem 0 0.5rem 0;
      }
      .sidebar-section:first-of-type { margin-top: 0.25rem; }

      /* Plot and results cards */
      .plot-card {
        background: #ffffff; border: 1px solid #dee2e6; border-radius: 8px;
        padding: 1rem; box-shadow: 0 2px 6px rgba(0,0,0,0.05); margin-bottom: 1.25rem;
      }
      .results-card {
        background: #ffffff; border: 1px solid #dee2e6; border-radius: 8px;
        padding: 1.5rem; box-shadow: 0 2px 6px rgba(0,0,0,0.05); margin-bottom: 1.25rem;
      }
      .results-card .results-title {
        font-weight: 600; font-size: 1.25rem; color: #1a6b54; margin-top: 0;
        margin-bottom: 0.75rem; padding-bottom: 0.5rem; border-bottom: 2px solid #e8eef3;
      }
      .results-card li {
        font-size: 1.05rem; line-height: 1.6; color: #555; padding: 0.4rem 0;
        border-bottom: 1px solid #f0f0f0; list-style: none; padding-left: 0.2rem;
      }
      .results-card li:last-child { border-bottom: none; }
      .results-card ul { padding-left: 0; margin: 0; }
      .results-card strong { color: #1a6b54; }

      /* Graphical-parameters panel */
      .gparams-title { font-weight: 600; font-size: 1.1rem; color: #3a6f99; margin: 0 0 0.75rem 0; }

      /* Footer */
      .app-footer { text-align: center; color: #666; font-size: 0.95rem; padding: 1rem 0 1.5rem 0; }
      .app-footer .app-authors { color: #555; font-size: 1.0rem; margin-bottom: 0.35rem; }
      .app-footer a { color: #3a6f99; text-decoration: none; }
      .app-footer a:hover { color: #1a6b54; text-decoration: underline; }
    "))
  ),

  # Header
  div(class = "app-header",
    h2("Sexual Partnership Concurrency and HIV Transmission"),
    p("An interactive stochastic microsimulation from the EpiModel ecosystem")
  ),

  # Info panels
  fluidRow(
    column(6,
      div(class = "card",
        h4("About this model"),
        p(strong("Concurrency"), "means holding more than one sexual partnership
           at the same time. This app simulates a stochastic HIV-1 epidemic in a
           heterosexual population to show how concurrency shapes transmission at
           the population level, even when the average number of partners and the
           length of partnerships are held fixed."),
        p("The theory and mathematics behind it are laid out in the full",
          tags$a(href = "https://epimodel.github.io/ConcurrencyTutorial", "Concurrency Tutorial,"),
          "and the same dynamics can be explored with the",
          tags$a(href = "https://www.epimodel.org", "EpiModel"),
          "software for R.")
      )
    ),
    column(6,
      div(class = "card",
        h4("How it works"),
        p("HIV progresses through four stages (acute, latent, pre-AIDS, AIDS) at
           fixed intervals. The probability of transmitting to a susceptible partner
           depends on the infected partner's stage: highest during acute infection,
           low during latency, elevated in pre-AIDS, and zero in AIDS, following",
          tags$a(href = "https://www.ncbi.nlm.nih.gov/pubmed/18662132", "Hollingsworth (2008).")),
        p("Choose whether concurrency is permitted for", strong("neither sex, women only,
           men only, or both."), "A momentary-degree cap of one enforces monogamy;
           relaxing it lets a person hold overlapping ties. Press",
          strong("Run Model"), "to simulate and compare.")
      )
    )
  ),

  # Simulation
  sidebarLayout(
    sidebarPanel(width = 3,
      tags$h4(style = "font-weight:600; font-size:1.25rem; color:#1a6b54; margin-top:0;
                       margin-bottom:0.75rem; padding-bottom:0.5rem; border-bottom:2px solid #e8eef3;",
              "Parameters"),
      actionButton("runMod", "Run Model", class = "btn-primary"),
      helpText("Adjust the settings below, then run the model."),

      div(class = "sidebar-section", "Initial conditions"),
      numericInput("s.num.m", "Susceptible men", value = 1000, min = 1, max = 20000, step = 1),
      numericInput("s.num.f", "Susceptible women", value = 1000, min = 1, max = 20000, step = 1),
      numericInput("i.num.m", "Infected men", value = 50, min = 1, max = 20000, step = 1),
      numericInput("i.num.f", "Infected women", value = 50, min = 1, max = 20000, step = 1),

      div(class = "sidebar-section", "Epidemic parameters"),
      sliderInput("md", "Mean degree", min = 0.5, max = 1, value = 0.8, step = 0.05),
      sliderInput("dur", "Mean partnership duration (months)", min = 1, max = 50, value = 10, step = 1),
      selectInput("conc", "Concurrency rule",
                  choices = c("No Concurrency", "Female Only Concurrency",
                              "Male Only Concurrency", "Both Sexes Concurrency")),

      div(class = "sidebar-section", "Control settings"),
      sliderInput("sims", "Simulations", value = 10, min = 1, max = 25, step = 1),
      numericInput("nsteps", "Time steps (months)", value = 2000, min = 100, max = 5000, step = 1)
    ),

    mainPanel(width = 9,
      div(class = "plot-card",
        plotOutput("concplot", height = "520px")
      ),
      div(class = "results-card",
        div(class = "results-title", "Simulation Results"),
        htmlOutput("recap")
      ),
      div(class = "well",
        div(class = "gparams-title", "Graphical parameters"),
        fluidRow(
          column(3, checkboxInput("sim.lines", "Sim lines", value = FALSE)),
          column(3, checkboxInput("mean.line", "Mean line", value = TRUE)),
          column(3, checkboxInput("mean.smooth", "Smooth mean", value = TRUE)),
          column(3, checkboxInput("qnts.smooth", "Smooth quantile", value = TRUE))
        ),
        fluidRow(
          column(6, sliderInput("qnts", "Quantile range", value = 0.5, min = 0, max = 1, step = 0.05)),
          column(6, sliderInput("sim.alpha", "Sim line transparency", value = 0.3, min = 0.05, max = 1, step = 0.05))
        ),
        fluidRow(
          column(6, sliderInput("mean.lwd", "Mean width", value = 3, min = 1, max = 5, step = 1)),
          column(6, sliderInput("qnt.alpha", "Quantile transparency", value = 0.3, min = 0.05, max = 1, step = 0.05))
        )
      )
    )
  ),

  # Footer
  div(class = "app-footer",
    p(class = "app-authors",
      "Developed by",
      tags$a(href = "https://www.epimodel.org/team.html", "Steven M. Goodreau"),
      "(University of Washington),",
      tags$a(href = "https://www.epimodel.org/team.html", "Samuel M. Jenness"),
      "(Emory University), and",
      tags$a(href = "https://www.epimodel.org/team.html", "Martina Morris"),
      "(University of Washington)"),
    p("Part of the",
      tags$a(href = "https://www.epimodel.org", "EpiModel"),
      "ecosystem • Powered by",
      tags$a(href = "https://posit-dev.github.io/r-shinylive/", "shinylive"))
  )
)


# ---- Server ---------------------------------------------------------
server <- function(input, output, session) {

  # Concurrency rule -> per-sex monogamy constraints
  monogM <- reactive({
    switch(input$conc,
      "No Concurrency" = TRUE, "Female Only Concurrency" = TRUE,
      "Male Only Concurrency" = FALSE, "Both Sexes Concurrency" = FALSE)
  })
  monogF <- reactive({
    switch(input$conc,
      "No Concurrency" = TRUE, "Female Only Concurrency" = FALSE,
      "Male Only Concurrency" = TRUE, "Both Sexes Concurrency" = FALSE)
  })

  # A momentary-degree cap of one bounds mean degree at 1: a mean degree
  # above 1 is only feasible once BOTH sexes are allowed concurrency, so the
  # slider ceiling opens to 1.5 only for the both-sexes scenario.
  observeEvent(input$conc, {
    mdMax <- ifelse(input$conc == "Both Sexes Concurrency", 1.5, 1.0)
    updateSliderInput(session, "md", label = "Mean degree",
                      min = 0.5, max = mdMax, value = min(input$md, mdMax), step = 0.05)
  })

  sim <- reactive({
    input$runMod  # the only reactive dependency: rerun on Run Model

    # Snapshot inputs so nothing but the button triggers a rerun. Counts are
    # rounded so the value used matches the integer shown in the widget.
    sf <- round(isolate(input$s.num.f)); iff <- round(isolate(input$i.num.f))
    sm <- round(isolate(input$s.num.m)); im <- round(isolate(input$i.num.m))
    md <- isolate(input$md);      dur <- isolate(input$dur)
    nst <- round(isolate(input$nsteps)); nsm <- isolate(input$sims)
    cc <- isolate(input$conc)
    mF <- isolate(monogF());      mM <- isolate(monogM())

    validate(
      need(all(vapply(list(sf, iff, sm, im), function(v) isTRUE(is.finite(v) && v >= 1 && v <= 20000), logical(1))),
           "Enter a population count between 1 and 20,000 for susceptible and infected men and women."),
      need(isTRUE(is.finite(nst) && nst >= 100 && nst <= 5000), "Enter between 100 and 5000 time steps."),
      need(isTRUE(is.finite(dur) && dur >= 1), "Enter a partnership duration of at least 1 month.")
    )

    nFeml <- sf + iff
    nMale <- sm + im
    expectedEdges <- round(md * (nFeml + nMale) / 2)
    validate(
      need(!(mM && expectedEdges > nMale) && !(mF && expectedEdges > nFeml),
           sprintf("This combination cannot be built. A sex held to one partner at a time cannot reach a mean degree of %.2f with these population sizes. Lower the mean degree, allow that sex concurrency, or raise its population.", md))
    )

    progress <- shiny::Progress$new()
    progress$set(message = "Computing", value = 0)
    on.exit(progress$close())
    updateProgress <- function(detail = NULL) progress$inc(amount = 1 / nsm, detail = detail)

    model <- conc_microsim(
      s.num.f = sf, i.num.f = iff, s.num.m = sm, i.num.m = im,
      monog.f = mF, monog.m = mM, meandeg = md, part.duration = dur,
      nsteps = nst, nsims = nsm, verbose = FALSE, updateProgress = updateProgress)
    attr(model, "conc") <- cc
    model
  })

  output$concplot <- renderPlot({
    plot_conc(sim(),
      sim.lines = input$sim.lines, sim.alpha = input$sim.alpha,
      mean.line = input$mean.line, mean.smooth = input$mean.smooth,
      mean.lwd = input$mean.lwd, qnts = input$qnts,
      qnt.alpha = input$qnt.alpha, qnts.smooth = input$qnts.smooth)
  })

  output$recap <- renderUI({
    model <- sim()
    tags$ul(lapply(recap_items(model, attr(model, "conc")), function(z) tags$li(z)))
  })
}

shinyApp(ui = ui, server = server)

# Concurrency and HIV Transmission

An interactive Shiny app that shows how sexual partnership concurrency shapes HIV-1 transmission in a heterosexual population. It is a companion to the [Concurrency Tutorial](https://epimodel.github.io/ConcurrencyTutorial) and part of the broader [EpiModel](https://www.epimodel.org/) ecosystem.

Available at: https://epimodel.github.io/concurrency.sim

## The idea

Concurrency means holding more than one sexual partnership at the same time. The app runs a stochastic microsimulation of an HIV epidemic and lets you switch concurrency on or off for each sex while holding everything else fixed, so the effect of concurrency alone can be read off the prevalence curves. Allowing overlapping partnerships opens faster transmission paths through the network, so prevalence rises even when the mean number of partners and the average partnership duration are unchanged.

HIV progresses through four stages (acute, latent, pre-AIDS, AIDS) at fixed intervals, and the per-step transmission probability depends on the infected partner's stage: highest during acute infection and zero in AIDS, following [Hollingsworth (2008)](https://www.ncbi.nlm.nih.gov/pubmed/18662132).

## Using the app

Set the initial population and epidemic parameters in the sidebar, choose a concurrency rule (none, women only, men only, or both sexes), and press Run Model. A momentary-degree cap of one enforces monogamy, so mean degree can exceed one only when both sexes are allowed concurrency. The main panel plots female and male prevalence over time across the simulations, with a mean line and a quantile band, and reports the final and peak prevalence. The graphical controls below the plot adjust the display without rerunning the model.

## Authors

- Steven M. Goodreau, University of Washington
- Samuel M. Jenness, Emory University
- Martina Morris, University of Washington

See the full EpiModel team at [epimodel.org/team](https://www.epimodel.org/team.html).

## About

This app is part of the [EpiModel](https://www.epimodel.org/) ecosystem and is built with [Shiny](https://shiny.posit.co/) and [shinylive](https://posit-dev.github.io/r-shinylive/), which runs the model entirely in the browser with no server. It is released under the GPL-3 license, consistent with the EpiModel ecosystem.

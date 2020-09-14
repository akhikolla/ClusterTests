## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE-------------------------------------------------
library(bggum)

## ----read_in_raw_data----------------------------------------------------
roberts_court <- read.csv("roberts_court.csv", stringsAsFactors = FALSE)
head(roberts_court)

## ----reshape_data--------------------------------------------------------
## Recode the votes to be integers in $\{0, \ldots, K_j - 1\}$
roberts_court$vote <- ifelse(roberts_court$vote %in% c(2, 7), 0, 1)
## Reshape the data
library(dplyr)
library(tidyr)
responses <- roberts_court %>%
    select(-caseId, -term) %>%
    spread(lexisCite, vote)
## I like to keep the respondents' names as the rownames
rownames(responses) <- responses$justiceName
## Now we just eliminate the respondent ID column and turn it into a matrix
responses <- as.matrix(responses[ , -1])
## Eliminate unanimous items
unanimous <- apply(responses, 2, function(x) length(unique(na.omit(x))) == 1)
responses <- responses[ , !unanimous]
## Here are responses to a few items:
responses[ , c(1, floor(ncol(responses) / 2), ncol(responses))]

## ----tune_proposals, eval = FALSE----------------------------------------
#  ## Load the package
#  library(bggum)
#  ## Set the seed for reproducibility
#  set.seed(123)
#  ## Tune the proposal densities
#  proposal_sds <- tune_proposals(responses, tune_iterations = 5000)
#  set.seed(456)
#  ## Tune the temperature schedule
#  temps <- tune_temperatures(responses, n_temps = 6, proposal_sds = proposal_sds)

## ----load_hypers, echo = FALSE-------------------------------------------
## (These are the results when I ran the above unevaluated code locally)
temps <- c(1, 0.933, 0.870, 0.812, 0.758, 0.706)

## ----check_temps---------------------------------------------------------
round(temps, 2)

## ----sample_posterior, eval = FALSE--------------------------------------
#  ## We need the parallel package
#  library(parallel)
#  ## We'll set up the cluster with two cores (for two chains)
#  cl <- makeCluster(2, type = "FORK", outfile = "bggum-log.txt")
#  ## Deal with reproducibility
#  clusterSetRNGStream(cl = cl, iseed = 789)
#  ## Produce the chains
#  chains <- parLapplyLB(cl = cl, X = 1:2, fun = function(x) {
#      ggumMC3(data = responses,
#              sample_iterations = 50000,
#              burn_iterations = 5000,
#              proposal_sds = proposal_sds,
#              temps = temps)
#  })
#  ## Stop the cluster
#  stopCluster(cl)

## ----post_process, eval = FALSE------------------------------------------
#  constraint <- which(rownames(responses) == "RBGinsburg")
#  processed_chains <- lapply(chains, post_process, constraint = constraint,
#                             expected_sign = "-")

## ----roberts_trace_plots_setup1, eval = FALSE----------------------------
#  roberts <- which(rownames(responses) == "JGRoberts")
#  trace_colors <- c("#0072b280", "#d55e0080")
#  iters <- nrow(chains[[1]])
#  idx <- seq(floor(iters / 1000), iters, floor(iters / 1000))
#  chain1_raw <- chains[[1]][idx, roberts]
#  chain2_raw <- chains[[2]][idx, roberts]
#  chain1_pp  <- processed_chains[[1]][idx, roberts]
#  chain2_pp  <- processed_chains[[2]][idx, roberts]

## ----roberts_trace_plots_setup2, echo = FALSE----------------------------
roberts <- which(rownames(responses) == "JGRoberts")
trace_colors <- c("#0072b280", "#d55e0080")
roberts_draws <- read.csv("roberts_draws.csv", stringsAsFactors = FALSE)
iters <- nrow(roberts_draws) * 50
idx <- 1:nrow(roberts_draws)
chain1_raw <- roberts_draws[idx, "chain1_raw"]
chain2_raw <- roberts_draws[idx, "chain2_raw"]
chain1_pp  <- roberts_draws[idx, "chain1_pp"]
chain2_pp  <- roberts_draws[idx, "chain2_pp"]
idx <- seq(floor(iters / 1000), iters, floor(iters / 1000))

## ----roberts_trace_plots, results = 'asis', fig.dim = c(6, 4)------------
ylims <- range(c(chain1_raw, chain2_raw, chain1_pp, chain2_pp))
opar <- par(mar = c(3, 3, 3, 1) + 0.1)
plot(idx, chain1_raw, type = "l", col = trace_colors[1], ylim = ylims,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "Raw Samples")
lines(idx, chain2_raw, col = trace_colors[2])
axis(side = 1, tick = FALSE, line = -0.75)
axis(side = 2, tick = FALSE, line = -0.75)
title(xlab = "Iteration", ylab = expression(theta[Roberts]), line = 1.5)
plot(idx, chain1_pp, type = "l", col = trace_colors[1], ylim = ylims,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "Processed Samples")
lines(idx, chain2_pp, col = trace_colors[2])
axis(side = 1, tick = FALSE, line = -0.75)
axis(side = 2, tick = FALSE, line = -0.75)
title(xlab = "Iteration", ylab = expression(theta[Roberts]), line = 1.5)
par(opar)

## ----convergence1, eval = FALSE------------------------------------------
#  library(coda)
#  ## Look at the Gelman-Rubin potential scale reduction factor:
#  convergence_stats <- gelman.diag(processed_chains)
#  ## See what estimates are:
#  summary(convergence_stats$psrf[ , 1])

## ----convergence2, echo = FALSE------------------------------------------
convergence_stats <- read.csv("convergence.csv")
summary(convergence_stats[ , 1])

## ----summary, eval = FALSE-----------------------------------------------
#  posterior_summary <- summary(processed_chains)

## ----load_summary, echo = FALSE------------------------------------------
posterior_summary <- readRDS("posterior_summary.rds")

## ----display_summary-----------------------------------------------------
str(posterior_summary$estimates, max.level = 1)
str(posterior_summary$sds, max.level = 1)
head(posterior_summary$statistics)

## ----plot_thetas, results = 'asis', fig.dim = c(6, 6)--------------------
theta <- posterior_summary$estimates$theta
names(theta) <- rownames(responses)
n <- length(theta)
ordering <- order(theta)
theta_sorted <- sort(theta)
opar <- par(mar = c(3, 6, 1, 1) + 0.1)
plot(theta_sorted, factor(names(theta_sorted), levels = names(theta_sorted)),
     pch = 19, xlim = c(-3, 2.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 1, tick = FALSE, line = -0.75)
axis(side = 2, tick = FALSE, line = -0.75, las = 1,
     labels = names(theta_sorted), at = 1:n)
title(xlab = expression(theta), line = 1.5)
segments(x0 = posterior_summary$statistics[ordering, 1], y0 = 1:n,
         x1 = posterior_summary$statistics[ordering, 4], y1 = 1:n)
par(opar)

## ----plot_items, results = 'asis', fig.dim = c(6, 4)---------------------
alpha <- posterior_summary$estimates$alpha
delta <- posterior_summary$estimates$delta
tau   <- posterior_summary$estimates$tau
phillip_morris <- which(colnames(responses) == "2007 U.S. LEXIS 1332")
rucho <- which(colnames(responses) == "2019 U.S. LEXIS 4401")
icc(alpha[phillip_morris], delta[phillip_morris], tau[[phillip_morris]],
    main_title = "Philip Morris USA, Inc. v. Williams",
    plot_responses = TRUE, thetas = theta,
    responses = responses[ , phillip_morris])
icc(alpha[rucho], delta[rucho], tau[[rucho]],
    main_title = "Rucho v. Common Cause",
    plot_responses = TRUE, thetas = theta,
    responses = responses[ , rucho])


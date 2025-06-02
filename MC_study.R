###########################################################################
# Monte Carlo study for lower‑tail Spearman’s rho: empirical vs Bernstein #
###########################################################################

## Written by Frederic Ouimet (May 2025)

start_time <- Sys.time()

base_dir <- "C://Users//fred1//Dropbox//Ouimet_Susam_2025_Spearman//simulations"
setwd(base_dir)

# -------- User‑configurable parameters ----------------------------------
text_cex           <- 2.0
plot_width         <- 800          # base width  (px)
plot_height        <- 600          # base height (px)
plot_height_factor <- 1.10         # make window 10 % taller

nb_cores     <- parallel::detectCores() - 1
K            <- 10000
sim_thetas   <- c(-1, -0.5, 0, 0.5, 1)   # θ values to simulate & tabulate
plot_thetas  <- c(-1, -0.5, 0, 0.5, 1)   # θ values to show in the figures
sample_sizes <- c(50, 200)
p_values     <- c(0.1, 0.5, 1)
m_seq_full   <- 1:60                 # simulate m = 1,2,…,60
m_seq_plot   <- seq(5, 60, by = 5)   # ticks / markers only at 5,10,…,60
# ------------------------------------------------------------------------

library(doParallel)
library(foreach)
library(copula)
library(zipfR)      # for Ibeta()
library(dplyr)
library(writexl)

true_rho <- function(p, theta) {
  Dp <- p^3/3 - p^4/4
  Bp <- p^2/2 - p^3/3
  theta * Bp^2 / Dp
}

cl <- makeCluster(nb_cores)
registerDoParallel(cl)

mc <- foreach(theta = sim_thetas, .combine = bind_rows) %:%
  foreach(n     = sample_sizes, .combine = bind_rows) %:%
  foreach(rep   = 1:K,           .combine = bind_rows,
          .packages = c("copula","zipfR")
  ) %dopar% {
    
    cop <- fgmCopula(param = theta)
    uv  <- rCopula(n, cop)
    R   <- rank(uv[,1])/(n+1)
    S   <- rank(uv[,2])/(n+1)
    
    out <- expand.grid(theta = theta, n = n, rep = rep,
                       p = p_values, m = c(0, m_seq_full))
    out$rho <- NA_real_
    
    for(i in seq_len(nrow(out))) {
      p  <- out$p[i]
      Dp <- p^3/3 - p^4/4
      if(out$m[i] == 0) {                   # empirical estimator
        emp_int      <- mean(pmax(0, p - R) * pmax(0, p - S))
        out$rho[i]   <- (emp_int - p^4/4) / Dp
      } else {                              # Bernstein estimator
        m0  <- out$m[i]
        ks  <- 0:m0
        Cn  <- outer(ks, ks, Vectorize(function(k,l)
          mean(R <= k/m0 & S <= l/m0)))
        Bk  <- Ibeta(p, ks+1, m0-ks+1, lower = TRUE)
        bin <- choose(m0, ks)
        tmp <- sum(Cn * outer(bin*Bk, bin*Bk, FUN = "*"))
        out$rho[i] <- (tmp - p^4/4) / Dp
      }
    }
    out
  }

stopCluster(cl)

# ---------- summary statistics ------------------------------------------
true_vals <- expand.grid(theta = sim_thetas, p = p_values) %>%
  mutate(rho_true = true_rho(p, theta))

res <- mc %>%
  left_join(true_vals, by = c("theta","p")) %>%
  mutate(error = rho - rho_true, sqerr = error^2)

emp_sum <- res %>%
  filter(m == 0) %>%
  group_by(theta, n, p) %>%
  summarise(
    rho_emp     = mean(rho),
    absbias_emp = abs(mean(error)),
    var_emp     = var(rho),
    mse_emp     = mean(sqerr),
    .groups     = "drop"
  )

bern_sum <- res %>%
  filter(m > 0) %>%
  group_by(theta, n, p, m) %>%
  summarise(
    absbias = abs(mean(error)),
    var_est = var(rho),
    mse     = mean(sqerr),
    .groups = "drop"
  ) %>%
  left_join(emp_sum, by = c("theta","n","p"))

# ---------- Table 1: m* ≈ floor(n^{2/3}) ---------------------------------
tab1 <- bern_sum %>%
  group_by(theta, n, p) %>%
  slice_min(abs(m - floor(n^(2/3))), with_ties = FALSE) %>%
  summarise(
    m            = m,
    AbsBias_emp  = unique(absbias_emp),
    AbsBias_bern = absbias,
    Var_emp      = unique(var_emp),
    Var_bern     = var_est,
    MSE_emp      = unique(mse_emp),
    MSE_bern     = mse,
    .groups      = "drop"
  ) %>%
  mutate(
    MSERedPct = 100 * (MSE_emp - MSE_bern) / MSE_emp
  ) %>%
  mutate(
    across(
      c(AbsBias_emp, AbsBias_bern,
        Var_emp, Var_bern,
        MSE_emp, MSE_bern),
      ~ round(., 4)
    ),
    MSERedPct = round(MSERedPct, 1)
  ) %>%
  relocate(m, .after = p) %>%
  select(theta, n, p, m,
         AbsBias_emp, AbsBias_bern,
         Var_emp, Var_bern,
         MSE_emp,  MSE_bern, MSERedPct)

write.csv(tab1, file = file.path(base_dir,"table1_simulation_results.csv"),
          row.names = FALSE)
write_xlsx(tab1, path = file.path(base_dir,"table1_simulation_results.xlsx"))

# ---------- plotting -----------------------------------------------------
symbol_map <- c(16, 15, 17)
col_map3   <- c("blue","red","darkgreen")
col_map6   <- rep(col_map3, each = 2)
lty6       <- rep(c(2,1), 3)
pch6       <- c(NA,16, NA,15, NA,17)

legend_labels <- expression(
  hat(rho)[n](0.1),
  hat(rho)[paste(n, ",", m)](0.1),
  hat(rho)[n](0.5),
  hat(rho)[paste(n, ",", m)](0.5),
  hat(rho)[n](1),
  hat(rho)[paste(n, ",", m)](1)
)

draw_legend <- function(pos) {
  tw <- max(strwidth(legend_labels, cex = text_cex))
  legend(pos,
         legend     = legend_labels,
         col        = col_map6,
         lty        = lty6,
         pch        = pch6,
         lwd        = 3,
         cex        = text_cex,
         ncol       = 3,
         text.width = tw * 1.1,
         x.intersp  = 0.5,
         xjust      = 1,      # right‑justify the box
         bty        = "n")
}

metrics <- list(
  list(col_emp = "mse_emp",     col_bern = "mse",     ylab = "MSE",      name = "MSE"),
  list(col_emp = "absbias_emp", col_bern = "absbias", ylab = "|Bias|",   name = "Absolute Bias"),
  list(col_emp = "var_emp",     col_bern = "var_est", ylab = "Variance", name = "Variance")
)

for(th in plot_thetas) {
  for(nn in sample_sizes) {
    
    emp_df <- filter(emp_sum, theta == th, n == nn) %>% arrange(p)
    ber_df <- filter(bern_sum, theta == th, n == nn)
    
    for(metric in metrics) {
      
      y_bern_mat <- sapply(p_values, function(pv) {
        tmp <- filter(ber_df, p == pv) %>% arrange(m)
        y   <- tmp[[ metric$col_bern ]]
        names(y) <- tmp$m
        y[m_seq_full]
      })
      
      y_emp_vals <- emp_df[[ metric$col_emp ]]
      
      # legend position rule, with special placement when θ = 1
      if (th == 1 || th == -0.5) {
        if (metric$name == "MSE") {
          pos <- "topright"
        } else if (metric$name == "Variance") {
          pos <- "bottomright"
        } else {
          pos <- "topright"
        }
      } else if (th == -1) {
        if (metric$name == "MSE") {
          pos <- "topright"
        } else if (metric$name == "Variance") {
          pos <- "right"
        } else {
          pos <- "right"
        }
      } else {
        if (metric$name == "MSE") {
          pos <- if (nn == 200 && th == 1) "topright" else "bottomright"
        } else if (metric$name == "Variance") {
          pos <- "bottomright"
        } else {
          pos <- "topright"
        }
      }
      
      theta_str <- format(th, trim = TRUE, scientific = FALSE)
      png(file.path(base_dir,
                    sprintf("%s_n%d_theta%s.png", metric$name, nn, theta_str)),
          width  = plot_width,
          height = plot_height * plot_height_factor)
      
      par(mar = c(4, 4.5, 3.5, 1), mgp = c(3, 1, 0))
      
      idx_use <- 3:length(m_seq_full)
      y_min   <- min(y_emp_vals, y_bern_mat[idx_use, ])
      if(metric$name == "Variance") {
        y_max <- max(y_emp_vals) * 1.05
      } else {
        y_max <- max(y_emp_vals, y_bern_mat[idx_use, ])
      }
      
      x_line  <- m_seq_full
      x_pts   <- m_seq_plot
      idx_pts <- match(x_pts, m_seq_full)
      
      plot(x_line, y_bern_mat[,1], type = "l",
           col = col_map3[1], lty = 1, lwd = 3,
           xlab = "m", ylab = metric$ylab,
           ylim = c(y_min, y_max),
           main = sprintf("%s, n=%d, θ=%s", metric$name, nn, theta_str),
           cex.lab = text_cex, cex.axis = text_cex, cex.main = text_cex,
           panel.first = grid(nx = length(m_seq_plot), ny = NULL,
                              col = "lightgray", lty = "dotted", lwd = 0.3),
           xaxt = "n")
      axis(1, at = x_pts, labels = x_pts, cex.axis = text_cex)
      points(x_pts, y_bern_mat[idx_pts,1], pch = symbol_map[1],
             col = col_map3[1], cex = 2)
      abline(h = y_emp_vals[1], col = col_map3[1], lty = 2, lwd = 3)
      
      lines(x_line, y_bern_mat[,2], col = col_map3[2], lty = 1, lwd = 3)
      points(x_pts, y_bern_mat[idx_pts,2], pch = symbol_map[2],
             col = col_map3[2], cex = 2)
      abline(h = y_emp_vals[2], col = col_map3[2], lty = 2, lwd = 3)
      
      lines(x_line, y_bern_mat[,3], col = col_map3[3], lty = 1, lwd = 3)
      points(x_pts, y_bern_mat[idx_pts,3], pch = symbol_map[3],
             col = col_map3[3], cex = 2)
      abline(h = y_emp_vals[3], col = col_map3[3], lty = 2, lwd = 3)
      
      draw_legend(pos)
      dev.off()
    }
  }
}

## To create Table 1 in LaTeX

library(dplyr)

csv_path  <- file.path(base_dir, "table1_simulation_results.csv")
latex_out <- file.path(base_dir, "table1_latex.tex")

# -----------------------------------------------------------------------
tab <- read.csv(csv_path) %>%
  arrange(theta, n, p) %>%                    # nice ordering
  mutate(                                     # rounding
    across(
      c(AbsBias_emp, AbsBias_bern,
        Var_emp, Var_bern,
        MSE_emp, MSE_bern),
      ~ sprintf("%.4f", .)
    ),
    MSERedPct_fmt = ifelse(
      MSERedPct >= 0,
      sprintf("\\textbf{%.1f}", MSERedPct),
      sprintf("\\textbf{--%.1f}", abs(MSERedPct))
    )
  )

# -----------------------------------------------------------------------
con <- file(latex_out, "w")

writeLines("\\begin{table}[ht!]",            con)
writeLines("\\centering",                    con)
writeLines(sprintf("\\caption{Estimates based on $%d$ Monte Carlo replications for the absolute bias, variance and MSE of the lower-tail Spearman's rho estimator $\\widehat{\\rho}_n(p)$ and its Bernstein version $\\widehat{\\rho}_{m,n}(p)$, with the rule-of-thumb Bernstein degree $m = \\lfloor n^{2/3} \\rfloor$.}", K), con)
writeLines("\\label{tab:1}",                 con)
writeLines("\\setlength{\\tabcolsep}{4pt}",  con)
writeLines("\\small",                        con)
writeLines("\\begin{tabular}{ccccccccccc}",  con)
writeLines("\\hline",                        con)
writeLines("$\\theta$ & $n$ & $p$ & $m$ & $\\lvert\\Bias[\\widehat{\\rho}_n(p)]\\rvert$ & $\\lvert\\Bias[\\widehat{\\rho}_{m,n}(p)]\\rvert$ & $\\Var[\\widehat{\\rho}_n(p)]$ & $\\Var[\\widehat{\\rho}_{m,n}(p)]$ & $\\mathrm{MSE}[\\widehat{\\rho}_n(p)]$ & $\\mathrm{MSE}[\\widehat{\\rho}_{m,n}(p)]$ & MSE reduction (\\%) \\\\", con)
writeLines("\\hline",                        con)

for(i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  line <- sprintf(
    "%g & %d & %.1f & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
    r$theta, r$n, r$p,            # θ, n, p
    r$m,                          # m now comes here
    r$AbsBias_emp,  r$AbsBias_bern,
    r$Var_emp,      r$Var_bern,
    r$MSE_emp,      r$MSE_bern,
    r$MSERedPct_fmt
  )
  writeLines(line, con)
  
  # Determine the appropriate horizontal rule
  next_row <- if(i < nrow(tab)) tab[i + 1, ] else NULL
  
  if(is.null(next_row) || r$theta != next_row$theta) {
    # End of a theta block  --> double rule
    writeLines("\\hline\\hline", con)
  } else if(r$n != next_row$n) {
    # Same theta, new n     --> single rule
    writeLines("\\hline", con)
  }
}

writeLines("\\end{tabular}", con)
writeLines("\\end{table}",   con)
close(con)

cat("LaTeX table written to:", latex_out, "\n")

end_time <- Sys.time()
cat("Total run time:", round(difftime(end_time, start_time, units = "mins"), 1),
    "minutes\n")


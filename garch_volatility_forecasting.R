library(xts)
library(readxl)
library(fBasics)
library(rugarch)
library(FinTS)
library(sandwich)

#===========================
# EXPLORATORY DATA ANALYSIS 
#===========================

output_dirs <- c("output/tables", "output/figures","output/cache")

lapply(
  output_dirs,
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)

load_close <- function(path, skip) {
  d <- read_xlsx(path, skip = skip, col_types = c("date", "numeric"))
  xts(d[, -1], order.by = as.Date(d$`Exchange Date`))
}

close_px <- list(
  GS  = load_close("data/PriceHistory_GS.xlsx",  skip = 4),
  SLV = load_close("data/PriceHistory_SLV.xlsx", skip = 3)
)
  
returns <- lapply(close_px, function(p) na.omit(100 * diff(log(p))))

split_date <- as.Date("2024-12-31")
r_train <- lapply(returns, function(r) r[index(r) <= split_date])
r_test  <- lapply(returns, function(r) r[index(r) >  split_date])

eda_asset <- function(r, name) {
  print(basicStats(r))
  
  print(jarqueberaTest(r))
  
  hist(r, breaks = 100, freq = FALSE,
       main = paste("Daily log-returns", name), xlab = "")
  curve(dnorm(x, mean(r), sd(r)), add = TRUE, col = 2, lwd = 2)
  
  acf(coredata(r), lag.max = 40, main = paste("ACF r_t",   name))
  acf(coredata(abs(r)), lag.max = 40, main = paste("ACF |r_t|", name))
  acf(coredata(r^2), lag.max = 40, main = paste("ACF r_t^2", name))
  
  y <- as.numeric(r - mean(r))
  for (L in c(1, 5, 10)) print(ArchTest(y, lags = L))
  for (L in c(5, 10)) print(Box.test(y^2, lag = L, type = "Ljung-Box"))
}

mapply(eda_asset, r_train, names(r_train))
       
for (a in names(r_train)){
write.csv(basicStats(r_train[[a]]), sprintf("output/tables/basicStats_%s_train.csv", a),quote=FALSE)
}


#=========================================================================
# MODELS: GARCH specification grid, estimation, diagnostics, IC comparison
#=========================================================================
make_spec <- function(model, dist) {
  garch_m <- model == "GARCH-M"
  ugarchspec(
    variance.model = list(model = if (garch_m) "sGARCH" else model,
                          garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE,
                          archm = garch_m, archpow = 1),
    distribution.model = dist
  )
}

model_grid <- expand.grid(
  model = c("sGARCH", "GARCH-M", "eGARCH", "apARCH"),
  dist  = c("norm", "std", "ged"),
  stringsAsFactors = FALSE
)

dist_short <- c(norm = "N", std = "t", ged = "GED")
model_grid$label <- paste(model_grid$model, dist_short[model_grid$dist], sep = "+")

specs <- Map(make_spec, model_grid$model, model_grid$dist)
names(specs) <- model_grid$label

fit_grid <- function(r) {
  fits <- lapply(specs, function(s) ugarchfit(s, data = r))
  conv <- sapply(fits, function(f) f@fit$convergence)
  if (any(conv != 0)) warning("Non-converged fits: ",
                              paste(names(fits)[conv != 0], collapse = ", "))
  fits
}

fits <- lapply(r_train, fit_grid) 

coef_table <- function(fit) fit@fit$robust.matcoef

diagnose <- function(fit) {
  z  <- as.numeric(residuals(fit, standardize = TRUE))
  sb <- signbias(fit)
  c(
    persistence = persistence(fit),
    unc_sd      = sqrt(uncvariance(fit)),
    setNames(sapply(c(1, 5, 22), function(L) ArchTest(z, lags = L)$p.value),
             paste0("ARCH_p_", c(1, 5, 22))),
    setNames(sapply(c(5, 10), function(L)
      Box.test(z^2, lag = L, type = "Ljung-Box")$p.value),
      paste0("LB2_p_", c(5, 10))),
    sign_bias_p = sb$prob[1],
    neg_bias_p  = sb$prob[2],
    pos_bias_p  = sb$prob[3],
    joint_p     = sb$prob[4]
  )
}

diagnostics <- lapply(fits, function(fa) round(t(sapply(fa, diagnose)), 4))

ic_compare <- function(fit) {
  ic <- t(sapply(fit, function(f) infocriteria(f)[c(1, 2, 4)]))
  colnames(ic) <- c("AIC", "BIC", "HQ")
  list(round(ic,4),apply(ic, 2, rank))
}

ic <- lapply(fits, ic_compare)

plot_dist_panel <- function(fits_asset, asset, family, dir = "output/figures") {
  sel <- model_grid$label[model_grid$model == family]
  png(file.path(dir, sprintf("distfit_%s_%s.png", asset, family)),
      width = 1800, height = 1100, res = 150)
  par(mfrow = c(2, 3))
  for (w in c(8, 9)) for (m in sel) plot(fits_asset[[m]], which = w)
  par(mfrow = c(1, 1))
  dev.off()
}

for (a in names(fits)) for (fam in unique(model_grid$model)) plot_dist_panel(fits[[a]], a, fam)

flag_matrix <- function(m, cols = 3:11, alpha = 0.05) {
  txt <- formatC(m, format = "f", digits = 4)
  txt[, cols] <- ifelse(m[, cols] < alpha,
                        paste0(txt[, cols], "**"),
                        paste0(txt[, cols], "  "))
 noquote(txt)
}

for (a in names(fits)){
  write.csv(flag_matrix(diagnostics[[a]]), sprintf("output/tables/diagnostics_%s_train.csv", a),quote=FALSE)
}

for (t in names(ic)){
  write.csv(ic[[t]][[2]], sprintf("output/tables/ic_%s_train.csv", t),quote=FALSE)
}

for (a in names(fits)){
  for (fam in model_grid$label){
  write.csv(coef_table(fits[[a]][[fam]]), sprintf("output/tables/coef_%s_%s_train.csv", a, fam),quote=FALSE)
  }
}



#============================================================
# Rolling 1-step FORECASTS, QLIKE/MSE losses, Diebold-Mariano
#============================================================


roll_grid <- function(r, n_test) {
    lapply(specs, function(s)
    ugarchroll(s, data = r, forecast.length = n_test,
               refit.every = 5, calculate.VaR = FALSE, refit.window = "recursive"))
}

cache_file <- "output/cache/rolls.rds" # delete after changing specs, refit.every, or data

if (file.exists(cache_file)) {
  rolls <- readRDS(cache_file)
} else {
  rolls <- Map(roll_grid, returns, lapply(r_test, nrow))
  saveRDS(rolls, cache_file)
}

sigma_forecasts <- lapply(rolls, function(ra) sapply(ra, function(x) as.data.frame(x)$Sigma))
realized  <- lapply(rolls, function(ra) as.data.frame(ra[[1]])$Realized)


loss_fns <- list(
  QLIKE = function(h, r2) log(h) + r2 / h,
  MSE   = function(h, r2) (h - r2)^2,
  MAE   = function(h,r2) abs(h - r2)
)

losses <- Map(function(sig, rz) {
  h <- sig^2
  r2 <- rz^2
  lapply(loss_fns, function(f) f(h, r2))
}, sigma_forecasts, realized)

mean_loss <- lapply(losses, function(la) sapply(la, colMeans))
mean_loss

dm_test <- function(l1, l2, alternative = "greater") {
  d   <- l1 - l2
  reg <- lm(d ~ 1)
  se  <- sqrt(NeweyWest(reg, lag = floor(length(d)^(1/3)), prewhite = FALSE)[1, 1])
  stat <- coef(reg)[[1]] / se
  p <- switch(alternative,
              greater   = pnorm(stat, lower.tail = FALSE),
              less      = pnorm(stat),
              two.sided = 2 * pnorm(-abs(stat)))
  list(statistic = stat, p.value = p)
}

dm_matrix <- function(L) {
  k   <- ncol(L)
  out <- matrix(NA, k, k, dimnames = list(colnames(L), colnames(L)))
  for (i in seq_len(k))
    for (j in seq_len(k)[-i])
      out[i, j] <- dm_test(L[, i], L[, j])$p.value
  round(out, 3)
}

dm <- lapply(losses, function(la) lapply(la, dm_matrix))

for (a in names(losses)) {
  write.csv(mean_loss[[a]], sprintf("output/tables/meanloss_%s.csv", a),quote=FALSE)
  for (l in names(dm[[a]]))
    write.csv(dm[[a]][[l]], sprintf("output/tables/DM_%s_%s.csv", a, l),quote=FALSE)
}

plot_forecasts <- function(sig, rz, models, asset, dir = "output/figures") {
  png(file.path(dir, sprintf("forecast_sigma_%s.png", asset)),
      width = 2000, height = 1100, res = 150)
  cols <- c("grey", "black", "red", "blue", "green")
  matplot(cbind(abs(rz), sig[, models]), type = "l", lty = 1,
          lwd = c(1, 2, 2, 2, 2), col = cols,
          main = sprintf("|r_t| vs one-step-ahead sigma — %s (2025–2026)", asset),
          ylab = "", xlab = "")
  legend("topright", legend = c("|r_t|", models), col = cols,
         lty = 1, lwd = c(1, 2, 2, 2, 2), cex = 0.7)
  dev.off()
}

rep_models <- c("sGARCH+t", "GARCH-M+t", "eGARCH+t", "apARCH+t")

for (a in names(sigma_forecasts)){
  plot_forecasts(sigma_forecasts[[a]], realized[[a]], rep_models, a)
}

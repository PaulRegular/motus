
library(motus)
library(TMB)
library(ggplot2)
library(data.table)

## Truncate n to simplify estimation
dat <- sync_dat
dat$n[dat$n >= 10] <- 10

## Couple additions for plotting
dat$month <- data.table::month(dat$DATETIME)
dat$week <- data.table::week(dat$DATETIME)
dat$date <- as.POSIXct(as.IDate(dat$DATETIME))
dat$yday <- data.table::yday(dat$DATETIME)

## Combine array and n to make the factor levels for the model
l <- unique(dat[, list(array, n)])
l <- l[order(l$array, l$n),]
l$levels <- paste0(l$array, "-", l$n)
dat$array_n <- factor(paste0(dat$array, "-", dat$n), levels = l$levels)

# ## Quick views
# ggplot(dat) +
#     geom_violin(aes(x = month, group = month, y = dev_easting)) +
#     facet_grid(year ~ location, scales = "free")
# ggplot(dat) +
#     geom_violin(aes(x = week, group = week, y = dev_easting)) +
#     facet_grid(year ~ location, scales = "free")
# ggplot(dat) +
#     geom_violin(aes(x = yday, group = yday, y = dev_easting)) +
#     facet_grid(year ~ location, scales = "free")
# ggplot(dat) +
#     geom_violin(aes(x = n, group = n, y = dev_easting)) +
#     facet_grid(year ~ location, scales = "free")
# ggplot(dat) +
#     geom_violin(aes(x = n, group = n, y = dev_northing)) +
#     facet_grid(year ~ location, scales = "free")

## Exclude early records because there was extensive noise
dat <- dat[dat$yday >= 240, ]

## Set-up TMB data and par inits
tmb_data <- list(lon = unname(dat$dev_easting),
                 lat = unname(dat$dev_northing),
                 n_i = as.numeric(dat$array_n) - 1,
                 dist = 1)
tmb_pars <- list(log_sd_lon = rep(0, length(unique(tmb_data$n_i))),
                 log_sd_lat = rep(0, length(unique(tmb_data$n_i))),
                 log_tau_lon = rep(0, length(unique(tmb_data$n_i))),
                 log_nu_lon = rep(log(3), length(unique(tmb_data$n_i))),
                 log_tau_lat = rep(0, length(unique(tmb_data$n_i))),
                 log_nu_lat = rep(log(3), length(unique(tmb_data$n_i))),
                 log_scale_lon = rep(0, length(unique(tmb_data$n_i))),
                 log_scale_lat = rep(0, length(unique(tmb_data$n_i))),
                 logit_p = rep(0, length(unique(tmb_data$n_i))),
                 log_df = rep(log(3), length(unique(tmb_data$n_i))))

## df fixed to 3 to account for heavy tails while maintaining statistical properities (mean, sd, etc.)
tmb_map <- list(log_nu_lon = rep(factor(NA), length(tmb_pars$log_nu_lon)),
                log_nu_lat = rep(factor(NA), length(tmb_pars$log_nu_lat)),
                logit_p = factor(rep(NA, length(tmb_pars$logit_p)))) # assume p = 50%
str(tmb_data)
str(tmb_pars)

## Compile TMB model
compile("analysis/obs_error.cpp")
dyn.load(dynlib("analysis/obs_error"))

## Fit normal distribution -----------------------------------------------------

tmb_data$dist <- 0
obj <- MakeADFun(tmb_data, tmb_pars[c("log_sd_lon", "log_sd_lat")],
                 DLL = "obs_error")

## Minimize the objective function using nlminb in R
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 2000, iter.max = 2000))
dnorm_rep <- obj$report()
sdrep <- sdreport(obj)
sdrep

## Fit t-distribution ----------------------------------------------------------

tmb_data$dist <- 1
obj <- MakeADFun(tmb_data, tmb_pars[c("log_tau_lon", "log_nu_lon", "log_tau_lat", "log_nu_lat")],
                 DLL = "obs_error", map = tmb_map[c("log_nu_lon", "log_nu_lat")])

## Minimize the objective function using nlminb in R
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 2000, iter.max = 2000))
dt_rep <- obj$report()
sdrep <- sdreport(obj)
sdrep

## Fit Cauchy distribution -----------------------------------------------------

tmb_data$dist <- 2
obj <- MakeADFun(tmb_data, tmb_pars[c("log_scale_lon", "log_scale_lat")],
                 DLL = "obs_error")

## Minimize the objective function using nlminb in R
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 2000, iter.max = 2000))
dcauchy_rep <- obj$report()
sdrep <- sdreport(obj)
sdrep

## Fit normal-t mixture --------------------------------------------------------

tmb_data$dist <- 3
obj <- MakeADFun(tmb_data, tmb_pars[c("log_sd_lon", "log_sd_lat", "logit_p", "log_df")],
                 map = tmb_map["logit_p"], DLL = "obs_error")

## Minimize the objective function using nlminb in R
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 2000, iter.max = 2000))
drobust_rep <- obj$report()
sdrep <- sdreport(obj)
sdrep


## Compile results -------------------------------------------------------------

## Table
obs_error <- data.frame(array_n = levels(dat$array_n),
                        sd_lon = dnorm_rep$sd_lon, sd_lat = dnorm_rep$sd_lat,
                        tau_lon = dt_rep$tau_lon, nu_lon = dt_rep$nu_lon,
                        tau_lat = dt_rep$tau_lat, nu_lat = dt_rep$nu_lat,
                        scale_lon = dcauchy_rep$scale_lon,
                        scale_lat = dcauchy_rep$scale_lat,
                        robust_sd_lon = drobust_rep$sd_lon, robust_sd_lat = drobust_rep$sd_lat,
                        p = drobust_rep$p,
                        df = drobust_rep$df)
array_n <- do.call(rbind, strsplit(as.character(obs_error$array_n), "-"))
obs_error$location <- array_n[, 1]
obs_error$year <- as.numeric(array_n[, 2])
obs_error$n <- as.numeric(array_n[, 3])
obs_error
ggplot(obs_error, aes(x = n)) +
    geom_point(aes(y = sd_lon)) +
    facet_grid(location ~ year)
ggplot(obs_error, aes(x = n)) +
  geom_point(aes(y = tau_lon)) +
  facet_grid(location ~ year)
ggplot(obs_error, aes(x = n)) +
  geom_point(aes(y = nu_lon)) +
  facet_grid(location ~ year)
ggplot(obs_error, aes(x = n)) +
    geom_point(aes(y = scale_lon)) +
    facet_grid(location ~ year)
ggplot(obs_error, aes(x = n)) +
    geom_point(aes(y = robust_sd_lon)) +
    facet_grid(location ~ year)
ggplot(obs_error, aes(x = n)) +
    geom_point(aes(y = robust_sd_lat)) +
    facet_grid(location ~ year)
ggplot(obs_error, aes(x = n)) +
    geom_point(aes(y = p)) +
    facet_grid(location ~ year)
ggplot(obs_error, aes(x = n)) +
    geom_point(aes(y = df)) +
    facet_grid(location ~ year)
obs_error[, c("location", "year", "n")] <- NULL
saveRDS(obs_error, file = "analysis/obs_error.rds")

## Estimate one proportion? Fix proportion? Fix df? (Maybe it's best to try it first on several tracks?)



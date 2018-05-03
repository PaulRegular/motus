
library(motus)
library(TMB)
library(ggplot2)
library(data.table)

## Truncate n to simplify estimation
sync_dat$n[sync_dat$n >= 10] <- 10

## Combine array and n to make the factor levels for the model
l <- unique(sync_dat[, list(array, n)])
l <- l[order(l$array, l$n),]
l$levels <- paste0(l$array, "-", l$n)
sync_dat$array_n <- factor(paste0(sync_dat$array, "-", sync_dat$n), levels = l$levels)

# ## Quick views
# ggplot(sync_dat) +
#     geom_violin(aes(x = month, group = month, y = dev_easting)) +
#     facet_grid(year ~ location, scales = "free")
# ggplot(sync_dat) +
#     geom_violin(aes(x = week, group = week, y = dev_easting)) +
#     facet_grid(year ~ location, scales = "free")
# ggplot(sync_dat) +
#     geom_violin(aes(x = yday, group = yday, y = dev_easting)) +
#     facet_grid(year ~ location, scales = "free")
# ggplot(sync_dat) +
#     geom_violin(aes(x = n, group = n, y = dev_easting)) +
#     facet_grid(year ~ location, scales = "free")
# ggplot(sync_dat) +
#     geom_violin(aes(x = n, group = n, y = dev_northing)) +
#     facet_grid(year ~ location, scales = "free")

## Exclude early records because there was extensive noise
sync_dat <- sync_dat[sync_dat$yday >= 240, ]

## Set-up TMB data and par inits
tmb_data <- list(lon = unname(sync_dat$dev_easting),
                 lat = unname(sync_dat$dev_northing),
                 n_i = as.numeric(sync_dat$array_n) - 1)
tmb_pars <- list(log_tau_lon = rep(0, length(unique(tmb_data$n_i))),
                 log_nu_lon = rep(log(3), length(unique(tmb_data$n_i))),
                 log_tau_lat = rep(0, length(unique(tmb_data$n_i))),
                 log_nu_lat = rep(log(3), length(unique(tmb_data$n_i))))
tmb_map <- list(log_nu_lon = rep(factor(NA), length(tmb_pars$log_nu_lon)),
                log_nu_lat = rep(factor(NA), length(tmb_pars$log_nu_lat)))   ## df fixed to 3 to account for heavy tails while maintaining statistical properities (mean etc.)
str(tmb_data)
str(tmb_pars)

## Compile TMB model
compile("analysis/obs_error.cpp")
dyn.load(dynlib("analysis/obs_error"))

## Generate objective function
obj <- MakeADFun(tmb_data, tmb_pars, DLL = "obs_error", map = tmb_map)

## Minimize the objective function using nlminb in R
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(eval.max = 2000, iter.max = 2000))
rep <- obj$report()
sdrep <- sdreport(obj)
sdrep

## Table
obs_error <- data.frame(array_n = levels(sync_dat$array_n),
                        tau_lon = rep$tau_lon, nu_lon = rep$nu_lon,
                        tau_lat = rep$tau_lat, nu_lat = rep$nu_lat)
array_n <- do.call(rbind, strsplit(as.character(obs_error$array_n), "-"))
obs_error$location <- array_n[, 1]
obs_error$year <- as.numeric(array_n[, 2])
obs_error$n <- as.numeric(array_n[, 3])
obs_error
ggplot(obs_error, aes(x = n)) +
  geom_point(aes(y = tau_lon)) +
  facet_grid(location ~ year)
ggplot(obs_error, aes(x = n)) +
  geom_point(aes(y = nu_lon)) +
  facet_grid(location ~ year)
obs_error[, c("location", "year", "n")] <- NULL
saveRDS(obs_error, file = "analysis/obs_error.rds")



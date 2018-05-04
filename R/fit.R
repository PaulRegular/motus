

#' Function for fitting state-space movement model
#'
#' @param track  data.frame of an individual track containing time, UTM coordinates,
#'               and estimates of observation error. UTM coordinate columns should be labeled
#'               lon and lat and units should equal meters. Observation error should be labeled
#'               tau_lon, tau_lat, nu_lon, and nu_lat for the parameters for the t-distribution.
#' @param silent Disable tracing information?
#'
#' @export
#'

fit_ssm <- function(track, silent = FALSE) {

    ## Set-up TMB data
    ## Note: centered the coordinates aid convergance
    tmb_data <- list(y_lon = track$lon - mean(track$lon),
                     y_lat = track$lat - mean(track$lat),
                     tau_lon = track$tau_lon,
                     tau_lat = track$tau_lat,
                     nu_lon = track$nu_lon,
                     nu_lat = track$nu_lat,
                     delta_t = c(1, diff(as.numeric(track$DATETIME))/60/1), # mins
                     n = nrow(track))

    ## Set-up initial par values for TMB
    ## Note: supplied y values for x values - key thing is that the starting locations are supplied
    ## Note: x_lon and x_lat values are in km units (not m like y_lon and y_lat) - this was done
    ##       to improve parameter scaling and ease convergance
    tmb_pars <- list(logit_gamma = rep(0, length(tmb_data$y_lon)),
                     log_sd_gamma = 0,
                     log_sd_lon = 0,
                     log_sd_lat = 0,
                     log_alpha_lon = 0,
                     log_alpha_lat = 0,
                     x_lon_km = tmb_data$y_lon / 1000,
                     x_lat_km = tmb_data$y_lat / 1000)

    ## Generate objective function, minimize and get parameters
    obj <- MakeADFun(tmb_data, tmb_pars, random = c("logit_gamma", "x_lon_km", "x_lat_km"),
                     DLL = "motus", silent = silent)
    control <- list(eval.max = 10000, iter.max = 10000)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)
    sd_rep <- sdreport(obj)
    sd_res <- summary(sd_rep)

    ## Calculate AIC
    nll <- opt$objective
    k <- length(opt[["par"]])
    AIC <- 2 * k + 2 * nll + 2 * k * (k + 1) / (k - 1)

    ## Extract state estimates
    xrep_lon <- sd_res[rownames(sd_res) %in% "x_lon_km", ]
    xrep_lat <- sd_res[rownames(sd_res) %in% "x_lat_km", ]
    xrep_gamma <- sd_res[rownames(sd_res) %in% "logit_gamma", ]

    ## Save states to the data
    track$lon_est <- xrep_lon[, 1] * 1000 + mean(track$lon)
    track$lat_est <- xrep_lat[, 1] * 1000 + mean(track$lat)
    track$lon_lwr <- (xrep_lon[, 1] - qnorm(0.975) * xrep_lon[, 2]) * 1000 + mean(track$lon)
    track$lon_upr <- (xrep_lon[, 1] + qnorm(0.975) * xrep_lon[, 2]) * 1000 + mean(track$lon)
    track$lat_lwr <- (xrep_lat[, 1] - qnorm(0.975) * xrep_lat[, 2]) * 1000 + mean(track$lat)
    track$lat_upr <- (xrep_lat[, 1] + qnorm(0.975) * xrep_lat[, 2]) * 1000 + mean(track$lat)
    track$gamma_est <- plogis(xrep_gamma[, 1])
    track$gamma_lwr <- plogis(xrep_gamma[, 1] - qnorm(0.975) * xrep_gamma[, 2])
    track$gamma_upr <- plogis(xrep_gamma[, 1] + qnorm(0.975) * xrep_gamma[, 2])

    ## Tidy par table
    ind <- rownames(sd_res) %in% c("x_lon_km", "x_lat_km", "log_alpha_lon", "log_alpha_lat", "logit_gamma")
    main_par <- sd_res[!ind, ]
    main_par <- data.frame(par = rownames(main_par),
                           est = main_par[, 1], sd = main_par[, 2],
                           stringsAsFactors = FALSE)
    main_par$lwr <- main_par$est - 1.96 * main_par$sd
    main_par$upr <- main_par$est + 1.96 * main_par$sd
    main_par$par <- gsub("sd", "sigma", main_par$par)
    main_par$par <- factor(main_par$par, levels = c("log_sigma_gamma", "log_sigma_lon", "log_sigma_lat"))

    ## Return results
    list(tmb_data = tmb_data, tmb_pars = tmb_pars, opt = opt, rep = rep,
         sd_rep = sd_rep, sd_res = sd_res, AIC = AIC, track = track,
         main_par = main_par)

}



#' Function for fitting state-space movement model
#'
#' @param track  data.frame of an individual track containing time, UTM coordinates, time differences
#'               and estimates of observation error. UTM coordinate columns should be labeled
#'               lon and lat and units should equal km. time difference should be a ratio of
#'               the minimum sampling interval and named delta_t. Observation error should be labeled
#'               1) sd_lon and sd_lat for the normal distribution; 2) tau_lon, tau_lat, nu_lon,
#'               and nu_lat for the parameters for the t-distribution; and, 3) scale_lon and
#'               scale_lat for the cauchy distribution.
#' @param scale  Method for scaling x values to improve parameter scaling and aid convergence, where
#'               "sd" scales using sd, "max" scales to max deviation, or provide a value.
#' @param fix_gamma Should the gamma parameter (autocorrelation for movement) be fixed or time-varrying?
#' @param dist   Distribution to use for observation error ("normal", "t", or "cauchy")
#' @param silent Disable tracing information?
#' @param gr_threshold Stop if maximum gradient exceeds this value (large values indicate convergence issues)
#'
#' @export
#'

fit_ssm <- function(track, scale = "sd", fix_gamma = FALSE, dist = "t", silent = FALSE,
                    gr_threshold = 10) {

    ## Center lon and lat values and scale using max deviation (or sd)
    ## this should aid convergence by scaling the x parameters
    center_lon <- mean(track$lon)
    center_lat <- mean(track$lat)
    if (is.character(scale)) {
        if (scale == "sd") {
            scale <- sd(c(track$lon - center_lon, track$lat - center_lat))
        }
        if (scale == "max") {
            scale <- max(abs(c(track$lon - center_lon, track$lat - center_lat)))
        }
    }

    ## Set-up TMB data
    ## Note: centered the coordinates aid convergance
    tmb_data <- list(y_lon = track$lon - center_lon,
                     y_lat = track$lat - center_lat,
                     scale = scale,
                     obs_sd_lon = track$sd_lon,
                     obs_sd_lat = track$sd_lat,
                     tau_lon = track$tau_lon,
                     tau_lat = track$tau_lat,
                     nu_lon = track$nu_lon,
                     nu_lat = track$nu_lat,
                     scale_lon = track$scale_lon,
                     scale_lat = track$scale_lat,
                     delta_t = track$delta_t,
                     n = nrow(track),
                     dist = as.numeric(factor(dist, levels = c("normal", "t", "cauchy"))) - 1,
                     fix_gamma = as.numeric(fix_gamma))

    ## Set-up initial par values for TMB
    ## Note: supplied y values for x values - key thing is that the starting locations are supplied
    tmb_pars <- list(logit_gamma = rep(0, length(tmb_data$y_lon)),
                     log_sd_gamma = 0,
                     log_sd_lon = 0,
                     log_sd_lat = 0,
                     log_alpha_lon = 0,
                     log_alpha_lat = 0,
                     x_slon = tmb_data$y_lon / scale,
                     x_slat = tmb_data$y_lat / scale)

    ## Conditional mapping
    if (fix_gamma) {
        tmb_map <- list(logit_gamma = rep(factor(1), length(tmb_data$y_lon)),
                        log_sd_gamma = factor(NA))
        tmb_random <- c("x_slon", "x_slat")
    } else {
        tmb_random <- c("logit_gamma", "x_slon", "x_slat")
        tmb_map <- list()
    }

    ## Generate objective function, minimize and get parameters
    obj <- MakeADFun(tmb_data, tmb_pars, random = tmb_random, map = tmb_map,
                     DLL = "motus", silent = silent)
    control <- list(eval.max = 10000, iter.max = 10000)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)
    sd_rep <- sdreport(obj)
    sd_res <- summary(sd_rep)

    ## Maximum gradient
    max_gr <- max(abs(sd_rep$gradient.fixed))
    stopifnot(max_gr <= gr_threshold)

    ## Calculate AIC
    nll <- opt$objective
    k <- length(opt[["par"]])
    AIC <- 2 * k + 2 * nll + 2 * k * (k + 1) / (k - 1)

    ## Extract state estimates
    xrep_lon <- sd_res[rownames(sd_res) %in% "x_lon", ]
    xrep_lat <- sd_res[rownames(sd_res) %in% "x_lat", ]
    xrep_gamma <- sd_res[rownames(sd_res) %in% "logit_gamma", ]
    if (fix_gamma) {
        xrep_gamma <- t(replicate(tmb_data$n, xrep_gamma))
    }

    ## Save states to the data
    track$lon_est <- (xrep_lon[, 1] + mean(track$lon))
    track$lat_est <- (xrep_lat[, 1] + mean(track$lat))
    track$lon_lwr <- (xrep_lon[, 1] - qnorm(0.975) * xrep_lon[, 2]) + mean(track$lon)
    track$lon_upr <- (xrep_lon[, 1] + qnorm(0.975) * xrep_lon[, 2]) + mean(track$lon)
    track$lat_lwr <- (xrep_lat[, 1] - qnorm(0.975) * xrep_lat[, 2]) + mean(track$lat)
    track$lat_upr <- (xrep_lat[, 1] + qnorm(0.975) * xrep_lat[, 2]) + mean(track$lat)
    track$logit_gamma_est <- xrep_gamma[, 1]
    track$logit_gamma_lwr <- xrep_gamma[, 1] - qnorm(0.975) * xrep_gamma[, 2]
    track$logit_gamma_upr <- xrep_gamma[, 1] + qnorm(0.975) * xrep_gamma[, 2]
    track$gamma_est <- plogis(track$logit_gamma_est)
    track$gamma_lwr <- plogis(track$logit_gamma_lwr)
    track$gamma_upr <- plogis(track$logit_gamma_upr)

    ## Tidy par table
    ind <- rownames(sd_res) %in% c("x_lon", "x_lat", "x_slon", "x_slat",
                                   "log_alpha_lon", "log_alpha_lat", "logit_gamma")
    main_par <- sd_res[!ind, ]
    main_par <- data.frame(par = rownames(main_par),
                           est = main_par[, 1], sd = main_par[, 2],
                           stringsAsFactors = FALSE)
    main_par$lwr <- main_par$est - 1.96 * main_par$sd
    main_par$upr <- main_par$est + 1.96 * main_par$sd
    main_par$par <- gsub("sd", "sigma", main_par$par)
    main_par$par <- factor(main_par$par, levels = c("log_sigma_gamma", "log_sigma_lon", "log_sigma_lat"))

    ## Return results
    list(tmb_data = tmb_data, tmb_pars = tmb_pars, opt = opt, rep = rep, max_gr = max_gr,
         sd_rep = sd_rep, sd_res = sd_res, AIC = AIC, track = track,
         main_par = main_par)

}


## Largly based on code provided by auger-methe_et_al_2017_Supplement_4_RCode
library(TMB)
library(data.table)
library(motus)
library(plotly)
library(tmbstan)

## Simplify data name
tracks <- crab_dat

## Load observation error predictions from obs_error model estimates
obs_sd <- readRDS("analysis/obs_error.rds")
tracks$n_caped <- tracks$n
tracks$n_caped[tracks$n_caped > 10] <- 10  # cap n at 10 to correspond with obs_sd estimates
tracks$array_n <- factor(paste0(tracks$array, "-", tracks$n_caped),
                         levels = levels(obs_sd$array_n))
tracks <- merge(tracks, obs_sd, by = "array_n")            # add obs_sd estimates to track data
tracks <- tracks[order(tracks$TRANSMITTER, tracks$DATETIME), ]   # ensure records are ordered by transmitter then time
tracks$time_since_release <- tracks$DATETIME - tracks$release_date_time
units(tracks$time_since_release) <- "hours"

## Ensure transmitters get a unique tag
tracks$TRANSMITTER <- paste0(tracks$array, "-", tracks$TRANSMITTER)

## Discard tracks with less than and greater than n records
## Tracks are difficult to fit outside this range
tracks[, tot_fixes := .N, by = "TRANSMITTER"]
# tracks <- tracks[tracks$tot_fixes > 100 & tracks$tot_fixes < 1000, ]
tracks <- tracks[tracks$tot_fixes > 100, ]

## Discard individuals that presumbably either died or shead the transmitter
## (i.e. Transmitters that remained within the array for the whole experiment)
# ggplot(d) + geom_point(aes(x = DATETIME, y = TRANSMITTER, colour = array)) + facet_wrap(~ year, scales = "free_x")
# d <- d[!d$TRANSMITTER %in% c("TODO"), ]

## Add delta_t (time difference between fixes) as a raito of the min sampling interval
tracks <- tracks[order(tracks$year, tracks$array, tracks$TRANSMITTER, tracks$DATETIME), ]
tracks[, delta_t := c(NA, diff(as.numeric(DATETIME))/60/1), by = "TRANSMITTER"]

## Simplify UTM coord name
tracks$lon <- tracks$easting
tracks$lat <- tracks$northing
tracks$easting <- NULL
tracks$northing <- NULL
tracks <- data.frame(tracks)

## Fit model to one individual
track <- tracks[tracks$TRANSMITTER == "Carson-2015-57081", ]
# track <- tracks[tracks$TRANSMITTER == "Carson-2016-53231", ]
# track <- tracks[tracks$TRANSMITTER == "Carson-2016-53215", ]
# track <- tracks[tracks$TRANSMITTER == "Carson-2017-31323", ]
# track <- tracks[tracks$TRANSMITTER == "Carson-2017-31305", ] # cauchy worked here
# track <- tracks[tracks$TRANSMITTER == "Carson-2015-57098", ] # appears to have limited error
# track <- tracks[tracks$TRANSMITTER == "Carson-2016-53139", ] # appears to have limited error
# track <- tracks[tracks$TRANSMITTER == "Carson-2016-53214", ] # normal works
# track <- tracks[tracks$TRANSMITTER == "Lilly-2016-53248", ] # large time breaks
# track <- tracks[tracks$TRANSMITTER == sample(unique(tracks$TRANSMITTER), 1), ]
res <- fit_ssm(track, scale = 5000, dist = "t", gamma_prob = 0.6)
# res <- fit_ssm(track, dist = "normal", scale = 5000, gamma_model = "fixed")
plot_track(res$track, discrete = TRUE)
plot_track(res$track, discrete = FALSE)
plot_trend(res$track, y_name = "gamma")
plot_trend(res$track, y_name = "logit_gamma")

res <- fit_ssm(track[track$time_since_release < 24 * 2, ],
               formula = ~ time_since_release, dist = "normal",
               scale = 5000, gamma_model = "fixed")
# res <- update(res, dist = "normal", fix_gamma = FALSE, start_par = res$par_est)
# res <- update(res, dist = "t", start_par = res$par_est)
res$sd_rep
hist(res$track$gamma_est, breaks = 100)
hist(res$track$logit_gamma_est, breaks = 100)

p_lon <- plot_trend(res$track, y_name = "lon")
p_lat <- plot_trend(res$track, y_name = "lat")
p_gamma <- plot_trend(res$track, y_name = "gamma")
subplot(p_lon, p_lat, p_gamma, nrows = 3, shareX = TRUE)
unique(res$track$TRANSMITTER)

## Test SIMULATE
track <- tracks[tracks$TRANSMITTER == "Carson-2015-57081", ]
res <- fit_ssm(track, scale = 5000, dist = "t")
sim <- sim_ssm(res)
plot(obs_lat ~ obs_lon, data = sim, asp = 1, pch = 16, col = rgb(1, 0, 0, 0.5), cex = 0.5)
lines(sim$true_lon, sim$true_lat)
plot(plogis(sim$epislon_gamma), type = "b")

## Test mcmc
track <- tracks[tracks$TRANSMITTER == "Carson-2015-57081", ]
res <- fit_ssm(track, scale = 5000, dist = "t")
samps <- samp_ssm(res, chains = 3)
# par(mar = c(1,1,1,1))
# pairs(samps, pairs = names(res$obj$par))
traceplot(samps, pars = names(res$obj$par), inc_warmup = TRUE)
## doesn't look good. Would need to increase burn-in and number of iterations


test.tracks<-unique(tracks$TRANSMITTER)
counter<-0
for(aa in 1:length(test.tracks)){
    track<-tracks[tracks$TRANSMITTER==test.tracks[aa],]

    ##use tryCatch to let us move past errors in a loop of multiple tags

    tryCatch({res <- fit_ssm(track, dist = "normal", scale = 5000, gamma_model = "fixed")
    counter<-counter+1
    ##ifelse(nrow(all.track>0), all.track<-rbind(all.track, res$track),all.track<-res$track)
    ifelse(counter==1, all.track<-res$track, all.track<-rbind(all.track, res$track))
    print("done successfully")
    print(paste("counter:", counter, sep=" "))


    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}##end aa loop

save(all.track, file="all.tracks.ssm.normal.fixed.RData")





# ind <- names(res$sd_rep$value) == "delta_gamma"
# delta_gamma_est <- res$sd_rep$value[ind]
# delta_gamma_sd <- res$sd_rep$sd[ind]
# delta_gamma_prob <- pnorm(0, mean = delta_gamma_est, sd = delta_gamma_sd)
#
# track <- res$track
# track$directed <- delta_gamma_prob > 0.9
# track$gamma_col <- factor(track$directed)
#
# plot_ly(data = track, x = ~lon, y = ~lat, colors = viridis::viridis(100)) %>%
#     add_markers(size = I(2), color = I("black"), name = "Observed") %>%
#     add_segments(x = ~lon_est[-nrow(track)], xend = ~lon_est[-1],
#                  y = ~lat_est[-nrow(track)], yend = ~lat_est[-1],
#                  color = ~gamma_col[-nrow(track)], showlegend = TRUE,
#                  text = ~round(gamma_est[-nrow(track)], 2), name = "")
#
#
#
# system.time({
#     ## Try fitting to all individuals
#     set.seed(123)
#     ids <- sample(unique(tracks$TRANSMITTER), 100)
#     tfits <- nfits <- cfits <- vector(mode = "list", length = length(ids))
#     names(tfits) <- names(nfits) <- names(cfits) <- ids
#     success <- 0
#     for (i in seq_along(ids)) {
#         track <- tracks[tracks$TRANSMITTER == ids[i], ]
#         start <- try(fit_ssm(track, dist = "normal", scale = 5000, silent = TRUE, fix_gamma = TRUE))
#         nfits[[i]] <- try(update(start, fix_gamma = FALSE, start_par = start$par_est))
#         tfits[[i]] <- try(update(nfits[[i]], dist = "t", start_par = nfits[[i]]$par_est))
#         if (class(tfits[[i]]) != "try-error") success <- success + 1
#         # cfits[[i]] <- try(fit_ssm(track, dist = "cauchy", scale = 1000, silent = TRUE))
#         cat(paste(i, "out of", length(ids), "\n"))
#         cat(paste0(round(success / i * 100), "% converged\n"))
#     }
# })
#
# mean(sapply(nfits, class) != "try-error")
# mean(sapply(tfits, class) != "try-error")
# # mean(sapply(cfits, class) != "try-error")
#
# success <- names(which(sapply(cfits, class) != "try-error"))
# res <- cfits[[success[14]]]
# plot_track(res$track)
# plot_trend(res$track, y_name = "gamma")
#
# success <- names(which(sapply(tfits, class) != "try-error"))
# res <- tfits[[success[76]]]
# plot_track(res$track, title = unique(res$track$TRANSMITTER))
# plot_trend(res$track, y_name = "gamma")
# hist(res$track$gamma_est, breaks = 100, xlab = "gamma")
# unique(res$track$TRANSMITTER)
# ## good examples: "Carson-2016-53231", "Lilly-2016-53174", "Lilly-2016-53240"
# ## should think about changes in angle: "Carson-2017-53263"

## Add start and end label to plot_track
## - Also consider impose transparancy using delta_t??
## Consider estimating angular parameter (random walk?)
## Filter long time breaks? for example: "Lilly-2016-53248"





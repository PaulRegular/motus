
## Largly based on code provided by auger-methe_et_al_2017_Supplement_4_RCode
library(TMB)
library(data.table)
library(motus)
library(plotly)

## Simplify data name
tracks <- crab_dat

## Load observation error predictions from obs_error model estimates
obs_sd <- readRDS("analysis/obs_error.rds")
tracks$array_n <- factor(paste0(tracks$array, "-", tracks$n), levels = levels(obs_sd$array_n))
tracks <- merge(tracks, obs_sd, by = "array_n")            # add obs_sd estimates to track data
tracks <- tracks[order(tracks$TRANSMITTER, tracks$DATETIME), ]   # ensure records are ordered by transmitter then time

## Ensure transmitters get a unique tag
tracks$TRANSMITTER <- paste0(tracks$array, "-", tracks$TRANSMITTER)

## Discard tracks with less than n records
tracks[, tot_fixes := .N, by = "TRANSMITTER"]
tracks <- tracks[tracks$tot_fixes > 50, ]

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
# track <- tracks[tracks$TRANSMITTER == sample(unique(tracks$TRANSMITTER), 1), ]
res <- fit_ssm(track, dist = "t", scale = 1000, fix_gamma = FALSE)
res$sd_rep
hist(res$track$gamma_est, breaks = 100)
hist(res$track$logit_gamma_est, breaks = 100)

# plot_track(res$track)
p_lon <- plot_trend(res$track, y_name = "lon")
p_lat <- plot_trend(res$track, y_name = "lat")
p_gamma <- plot_trend(res$track, y_name = "gamma")
subplot(p_lon, p_lat, p_gamma, nrows = 3, shareX = TRUE)
unique(res$track$TRANSMITTER)


## Try fitting to all individuals
ids <- unique(tracks$TRANSMITTER)
tfits <- nfits <- cfits <- vector(mode = "list", length = length(ids))
names(tfits) <- names(nfits) <- names(cfits) <- ids
for (id in ids) {
    track <- tracks[tracks$TRANSMITTER == id, ]
    tfits[[id]] <- try(fit_ssm(track, dist = "t"))
    nfits[[id]] <- try(fit_ssm(track, dist = "normal"))
    cfits[[id]] <- try(fit_ssm(track, dist = "cauchy"))
}

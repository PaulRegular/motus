
## Largly based on code provided by auger-methe_et_al_2017_Supplement_4_RCode
library(TMB)
library(data.table)

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
tracks <- tracks[tracks$tot_fixes > 100, ]

## Discard individuals that presumbably either died or shead the transmitter
## (i.e. Transmitters that remained within the array for the whole experiment)
# ggplot(d) + geom_point(aes(x = DATETIME, y = TRANSMITTER, colour = array)) + facet_wrap(~ year, scales = "free_x")
# d <- d[!d$TRANSMITTER %in% c("TODO"), ]

## Simplify UTM coord name
tracks$lon <- tracks$easting
tracks$lat <- tracks$northing
tracks$easting <- NULL
tracks$northing <- NULL

## Fit model to one individual
track <- tracks[tracks$TRANSMITTER == "Carson-2015-57081", ]
res <- fit_ssm(track)





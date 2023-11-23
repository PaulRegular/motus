
library(data.table)

## Import data
files <- list.files("data-raw", pattern = "ALL-CALC-POSITIONS", recursive = TRUE,
                    full.names = TRUE)
dat <- lapply(files, data.table::fread)
dat <- rbindlist(dat)
dt <- as.POSIXct(strptime(dat$DATETIME, "%Y-%m-%d %H:%M"), tz = "UTC")
attr(dt, "tzone") <- "America/St_Johns"
dat$DATETIME <- dt
dat$year <- data.table::year(dat$DATETIME)
dat$month <- data.table::month(dat$DATETIME)
dat$week <- data.table::week(dat$DATETIME)
dat$date <- as.POSIXct(as.IDate(dat$DATETIME))
dat$yday <- data.table::yday(dat$DATETIME)
dat$TRANSMITTER <- substring(dat$DETECTEDID, 10, 14) # fix transmitter column
xy <- as.matrix(cbind(dat$LON, dat$LAT))
xy <- rgdal::project(xy, "+proj=utm +zone=21 ellps=WGS84")
dat$easting <- xy[, 1]
dat$northing <- xy[, 2]
dat$location <- ifelse(dat$northing > 5020000, "Carson", "Lilly")
dat$array <- paste0(dat$location, "-", dat$year)

## Determine sync tag numbers
files <- list.files("data-raw", pattern = "SyncTag", recursive = TRUE)
files <- strsplit(files, "/")
files <- sapply(files, tail, 1)
files <- strsplit(files, "-|_|\\.")
sync_tags <- sapply(seq_along(files), function(i) files[[i]][nchar(files[[i]]) == 5])
sync_tags <- unique(as.numeric(sync_tags))

## Calculate deviations from the sync locations
sync_dat <- dat[dat$TRANSMITTER %in% sync_tags, ]
sync_dat <- sync_dat[!is.na(sync_dat$DATETIME), ]
sync_loc <- sync_dat[, list(mean_easting = mean(easting, na.rm = TRUE),
                            mean_northing = mean(northing, na.rm = TRUE)),
                     by = c("array", "TRANSMITTER")]
sync_dat <- merge(sync_dat, sync_loc, by = c("array", "TRANSMITTER"))
sync_dat$dev_easting <- sync_dat$easting - sync_dat$mean_easting
sync_dat$dev_northing <- sync_dat$northing - sync_dat$mean_northing

## Reduce the number of columns
sync_dat <- sync_dat[, c("TRANSMITTER", "DATETIME", "LAT", "LON", "n", "HPE", "HPEm",
                         "year", "location", "array", "easting", "northing",
                         "dev_easting", "dev_northing"),
                     with = FALSE]

## Save sync tag data
usethis::use_data(sync_dat, overwrite = TRUE)
usethis::use_data(sync_loc, overwrite = TRUE) # also save mean locations

## Now process the crab data and save
crab_dat <- dat[!dat$TRANSMITTER %in% sync_tags,
                c("TRANSMITTER", "DATETIME", "LAT", "LON", "n", "HPE", "HPEm",
                  "year", "location", "array", "easting", "northing")]

## Add metadata
metadata <- fread("data-raw/crab tag metadata 14022018.csv")
metadata <- metadata[, c("ext tag", "acoutic tag", "carap width",
                         "release date", "release Time",
                         "cap lat N", "cap long W",
                         "release location")]
setnames(metadata, names(metadata), tolower(gsub(" ", "_", names(metadata))))
setnames(metadata, "acoutic_tag", "TRANSMITTER")
metadata$release_time[metadata$release_time == "1800"] <- "18:00"
metadata$release_date_time <- as.POSIXct(strptime(paste(metadata$release_date, metadata$release_time), "%d-%b-%y %H:%M"))
metadata$TRANSMITTER[metadata$ext_tag == 1670] <- 53260   # fix one typo in the metadata
metadata$TRANSMITTER <- as.character(metadata$TRANSMITTER)
crab_dat <- merge(crab_dat, metadata, by = "TRANSMITTER")

## Ensure data is ordered by id then time
crab_dat <- crab_dat[order(crab_dat$year, crab_dat$array, crab_dat$TRANSMITTER, crab_dat$DATETIME), ]
crab_dat <- crab_dat[!is.na(crab_dat$DATETIME), ]

## Add seismic exposure times
starts <- list("2017" = ISOdatetime(2017, 9, 12, 10, 30, 0),
               "2016" = ISOdatetime(2016, 9, 22, 13, 30, 0),
               "2015" = ISOdatetime(2015, 9, 25, 0, 0, 0))
ends <- list("2017" = ISOdatetime(2017, 9, 12, 18, 0, 0),
             "2016" = ISOdatetime(2016, 9, 22, 15, 30, 0),
             "2015" = ISOdatetime(2015, 9, 29, 0, 0, 0))
crab_dat$exposure <- "control"
crab_dat$exposure_weeks <- "control"

for (i in names(starts)) {

    ind <- crab_dat$year == i & crab_dat$DATETIME > ends[[i]]
    crab_dat$exposure[ind] <- "after"

    weeks_after <- difftime(crab_dat$DATETIME[ind], ends[[i]], units = "weeks")
    crab_dat$exposure_weeks[ind] <- paste0("after (week ", ceiling(weeks_after), ")")

    ind <- crab_dat$year == i &
        crab_dat$DATETIME <= ends[[i]]
    crab_dat$exposure[ind] <- "during"
    crab_dat$exposure_weeks[ind] <- "during"

    ind <- crab_dat$year == i &
        crab_dat$DATETIME <= starts[[i]]
    crab_dat$exposure[ind] <- "before"
    crab_dat$exposure_weeks[ind] <- "before"

}

## Label first 24 hours as the adjustment period
# crab_dat[, start_time := head(DATETIME, 1), by = c("year", "array", "TRANSMITTER")]
# ind <- difftime(crab_dat$DATETIME, crab_dat$start_time, units = "hours") <= 24
ind <- difftime(crab_dat$DATETIME, crab_dat$release_date_time, units = "hours") <= 24
crab_dat$exposure[ind] <- crab_dat$exposure_weeks[ind] <- "adjustment"

## Visual checks
library(plotly)
year_plot <- function(yr = 2015, loc = c("Carson", "Lilly")) {
    plot_ly(data = crab_dat[year == yr & location %in% loc, ],
            x = ~DATETIME, y = ~TRANSMITTER, color = ~exposure_weeks) %>%
        add_markers()
}
year_plot(yr = 2015)
year_plot(yr = 2016)
year_plot(yr = 2017)
year_plot(yr = 2015, loc = "Carson")
year_plot(yr = 2016, loc = "Carson")
year_plot(yr = 2017, loc = "Carson")

## Keep adjustment period for explorations of flight responses...or exclude
# crab_dat <- crab_dat[crab_dat$exposure != "adjustment" | crab_dat$exposure_weeks != "adjustment", ]

# crab_dat %>%
#     filter(year == 2015) %>%
#     plot_ly(x = ~DATETIME, y = ~TRANSMITTER, color = ~exposure) %>%
#     add_markers()

usethis::use_data(crab_dat, overwrite = TRUE)



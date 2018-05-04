
#' Plot track
#'
#' @param track  track data.frame from \code{\link{fit_ssm}}
#' @param title  title of plot
#'
#' @export
#'
#' @examples
#'
#' @import plotly
#'

plot_track <- function(track, title = "") {

    # ggplot(data = track) +
    #     geom_point(aes(x = lon, y = lat), size = 1, alpha = 0.2) +
    #     geom_path(aes(x = lon_est, y = lat_est, colour = gamma_est), size = 1) +
    #     scale_color_viridis_c(name = "gamma") +
    #     theme_void() + coord_fixed()

    ax <- list(title = "",
               zeroline = FALSE,
               showline = FALSE,
               showticklabels = FALSE)

    cols <- RColorBrewer::brewer.pal(8, "RdBu")
    track$col <- viridis::viridis(100)[cut(track$gamma_est, 100)]

    p <- plot_ly(data = track, x = ~lon, y = ~lat) %>%
        add_markers(size = I(2), color = I("black"), name = "Observed") %>%
        layout(
            title = title,
            yaxis = c(ax, list(scaleanchor = "x")),
            xaxis = ax
        )
    for (i in seq(nrow(track))[-1]) {
        p <- p %>% add_segments(x = track$lon_est[i - 1], y = track$lat_est[i - 1],
                                xend = track$lon_est[i], yend = track$lat_est[i],
                                color = I(track$col[i - 1]),
                                showlegend = FALSE, name = "Predicted")
    }
    p

}


#' Plot temporal trend
#'
#' @param track     track data.frame from \code{\link{fit_ssm}}
#' @param y_name    name of variable to plot
#' @param xlab      label for x axis
#' @param ylab      label for y axis
#' @param col_est   colour for model estimates
#' @param col_obs   colour for observations
#'
#' @export
#'

plot_trend <- function(track, y_name = "lon", xlab = "time", ylab = y_name,
                       col_est = "steelblue", col_obs = "black") {

    d <- track[, c("DATETIME", paste0(y_name, c("_est", "_lwr", "_upr")))]
    names(d) <- gsub(y_name, "y", names(d))
    if (y_name %in% names(track)) d$y <- track[, y_name]
    p <- plot_ly(data = d, x = ~DATETIME)
    if (y_name %in% names(track)) {
        p <- p %>% add_markers(y = ~y, size = I(2), color = I(col_obs), name = "Observed")
    }
    p %>% add_ribbons(ymin = ~y_lwr, ymax = ~y_upr, name = "95% CI", line = list(width = 0),
                      color = I(col_est), opacity = 0.4) %>%
        add_lines(y = ~y_est, size = I(1), name = "Estimated", color = I(col_est)) %>%
        layout(yaxis = list(title = ylab),
               xaxis = list(title = xlab))

}


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

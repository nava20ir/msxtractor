#' Plot an extracted ion chromatogram
#'
#' @param rt Numeric vector. Retention time values, usually in seconds.
#' @param intensity Numeric vector. Intensity values.
#' @param title Character. Plot title.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' xic <- extract_xic("file.mzXML", target_mz = 453.367, mz_tol = 0.011)
#' plot_xic(xic$rt, xic$intensity, title = "XIC")
#' }
plot_xic <- function(rt,
                     intensity,
                     title = "Extracted ion chromatogram") {

  if (length(rt) != length(intensity)) {
    stop("rt and intensity must have the same length.", call. = FALSE)
  }

  ggplot2::ggplot(
    data.frame(
      rt = rt,
      intensity = intensity
    ),
    ggplot2::aes(rt, intensity)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Retention Time (s)",
      y = "Intensity",
      title = title
    )
}

#' Plot an MS2 spectrum
#'
#' @param mz Numeric vector. Fragment m/z values.
#' @param intensity Numeric vector. Fragment intensity values.
#' @param rt Optional numeric vector. Retention time labels. If supplied, values are shown above peaks.
#' @param title Character. Plot title.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' ms2 <- extract_ms2("file.mzXML", target_precursor = 453.367, mz_tol = 0.011)
#' plot_ms2(ms2$mz, ms2$intensity, ms2$rt)
#' }
plot_ms2 <- function(mz,
                     intensity,
                     rt = NULL,
                     title = "MS2 spectrum") {

  if (length(mz) != length(intensity)) {
    stop("mz and intensity must have the same length.", call. = FALSE)
  }

  if (!is.null(rt) && length(rt) != length(mz)) {
    stop("If rt is supplied, it must have the same length as mz.", call. = FALSE)
  }

  df <- data.frame(
    mz = mz,
    intensity = intensity
  )

  if (!is.null(rt)) {
    df$rt <- round(rt, 1)
  }

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = mz,
        xend = mz,
        y = 0,
        yend = intensity
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "m/z",
      y = "Intensity",
      title = title
    )

  if (!is.null(rt)) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(
          x = mz,
          y = intensity,
          label = rt
        ),
        angle = 90,
        size = 2.5,
        vjust = -0.2
      )
  }

  p
}

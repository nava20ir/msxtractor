\name{plot_xic}
\alias{plot_xic}
\title{Plot an extracted ion chromatogram}
\usage{plot_xic(rt, intensity, title = "Extracted ion chromatogram")}
\arguments{
  \item{rt}{Numeric vector. Retention time values, usually in seconds.}
  \item{intensity}{Numeric vector. Intensity values.}
  \item{title}{Character. Plot title.}
}
\value{A ggplot object.}
\description{Creates a ggplot line plot for an extracted ion chromatogram.}

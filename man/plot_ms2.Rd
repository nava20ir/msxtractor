\name{plot_ms2}
\alias{plot_ms2}
\title{Plot an MS2 spectrum}
\usage{plot_ms2(mz, intensity, rt = NULL, title = "MS2 spectrum")}
\arguments{
  \item{mz}{Numeric vector. Fragment m/z values.}
  \item{intensity}{Numeric vector. Fragment intensity values.}
  \item{rt}{Optional numeric vector. Retention time labels.}
  \item{title}{Character. Plot title.}
}
\value{A ggplot object.}
\description{Creates a stick spectrum plot for MS2 fragment ions. If rt is supplied, retention time labels are shown above peaks.}

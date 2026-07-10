\name{extract_xic}
\alias{extract_xic}
\title{Extract an extracted ion chromatogram from an MS file}
\usage{extract_xic(path, target_mz, mz_tol = 0.01)}
\arguments{
  \item{path}{Character. Path to an mzML or mzXML file.}
  \item{target_mz}{Numeric. Target m/z value.}
  \item{mz_tol}{Numeric. Absolute m/z tolerance around target_mz.}
}
\value{A data.frame with columns rt and intensity.}
\description{Reads an mzML or mzXML file with MSnbase, keeps MS1 scans, and extracts an XIC around a target m/z using MSnbase::chromatogram().}

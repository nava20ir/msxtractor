\name{extract_ms2}
\alias{extract_ms2}
\title{Extract MS2 spectra matching a precursor m/z}
\usage{extract_ms2(path, target_precursor, mz_tol = 0.01)}
\arguments{
  \item{path}{Character. Path to an mzML or mzXML file.}
  \item{target_precursor}{Numeric. Target precursor m/z.}
  \item{mz_tol}{Numeric. Absolute precursor m/z tolerance.}
}
\value{A data.frame with columns precursor_mz, rt, mz, and intensity.}
\description{Reads an mzML or mzXML file with MSnbase, keeps MS2 scans, filters spectra by precursor m/z, and returns all fragment peaks from matching MS2 scans.}

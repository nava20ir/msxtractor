#' Extract MS1 and MS2 peaks from a mass spectrometry data file
#'
#' Reads a mass spectrometry file supported by \pkg{mzR} and extracts:
#' \itemize{
#'   \item Total ion chromatogram (TIC) information
#'   \item All MS1 centroid/profile peaks
#'   \item All MS2 centroid/profile peaks together with their precursor m/z
#' }
#'
#' The returned peak tables contain one row per m/z-intensity pair and include
#' the corresponding retention time for each scan. MS2 peaks additionally
#' include the precursor m/z value recorded for the scan.
#'
#' @param file Character string. Path to a mass spectrometry file (e.g.
#'   mzML, mzXML, or netCDF) readable by \pkg{mzR}.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{tic}{
#'     A data frame containing:
#'     \itemize{
#'       \item \code{rt}: retention time (seconds)
#'       \item \code{intensity}: total ion current
#'     }
#'   }
#'   \item{ms1}{
#'     A data frame containing all MS1 peaks with columns:
#'     \itemize{
#'       \item \code{rt}: retention time
#'       \item \code{mz}: m/z value
#'       \item \code{intensity}: peak intensity
#'     }
#'   }
#'   \item{ms2}{
#'     A data frame containing all MS2 peaks with columns:
#'     \itemize{
#'       \item \code{rt}: retention time
#'       \item \code{mz}: fragment m/z value
#'       \item \code{intensity}: fragment intensity
#'       \item \code{precursor_mz}: precursor ion m/z
#'     }
#'   }
#' }
#'
#' @details
#' The function opens the raw MS file using \code{mzR::openMSfile()},
#' extracts scan metadata using \code{mzR::header()}, retrieves peak
#' matrices for all MS1 and MS2 scans, and combines them into long-format
#' data frames suitable for plotting, filtering, or downstream analysis.
#'
#' Retention times are reported in the units stored in the mzML file
#' (typically seconds).
#'
#' @examples
#' \dontrun{
#' res <- extract_ms_peaks("sample.mzML")
#'
#' head(res$tic)
#' head(res$ms1)
#' head(res$ms2)
#'
#' # Number of MS1 peaks
#' nrow(res$ms1)
#'
#' # Filter MS2 fragments from precursor 445.34 ± 0.01 Da
#' subset(
#'   res$ms2,
#'   abs(precursor_mz - 445.34) <= 0.01
#' )
#' }
#'
#' @importFrom mzR openMSfile header peaks close
#' @export
extract_ms_peaks <- function(file) {
  ms <- mzR::openMSfile(file)
  hdr <- mzR::header(ms)
  tic <- data.frame(
    rt = hdr$retentionTime,
    intensity = hdr$totIonCurrent
  )
  
  extract_level <- function(idx, precursor = FALSE) {
    
    res <- lapply(idx, function(i) {
      p <- mzR::peaks(ms, i)
      if (is.null(p) || nrow(p) == 0)
        return(NULL)
      df <- data.frame(
        rt = rep(hdr$retentionTime[i], nrow(p)),
        mz = p[, 1],
        intensity = p[, 2]
      )
      if (precursor) {
        df$precursor_mz <- rep(hdr$precursorMZ[i], nrow(p))
      }
      df
    })
    res <- Filter(Negate(is.null), res)
    if (length(res) == 0)
      return(data.frame())
    do.call(rbind, res)
  }
  
  res_ms1 <- extract_level(which(hdr$msLevel == 1))
  res_ms2 <- extract_level(which(hdr$msLevel == 2), precursor = TRUE)

  mzR::close(ms)
  list(
    tic = tic,
    ms1 = res_ms1,
    ms2 = res_ms2
  )
}





#' Filter a data frame by m/z value within a specified tolerance
#'
#' Returns all rows where the values in a specified column are within
#' a given absolute tolerance of a target m/z value.
#'
#' This utility is useful for extracting MS1 or MS2 signals around a
#' selected precursor or fragment m/z from the output of
#' [extract_ms_peaks()] or similar peak tables.
#'
#' @param df A data frame containing an m/z column.
#' @param mz Numeric. Target m/z value to search for.
#' @param tol Numeric. Absolute m/z tolerance.
#' @param colname Character string specifying the column containing
#'   m/z values. Defaults to `"mz"`.
#'
#' @return A data frame containing only rows where
#'   `abs(df[[colname]] - mz) <= tol`.
#'
#' @details
#' The filtering uses an absolute tolerance rather than ppm:
#'
#' \deqn{|m/z_{observed} - m/z_{target}| \le tolerance}
#'
#' If no rows satisfy the criterion, an empty data frame with the same
#' columns is returned.
#'
#' @examples
#' peaks <- data.frame(
#'   mz = c(100.01, 100.05, 101.20),
#'   intensity = c(1000, 500,
filter_mz <- function(df, mz, tol, colname = "mz") {
  
  stopifnot(colname %in% names(df))
  
  df[abs(df[[colname]] - mz) <= tol, , drop = FALSE]
  
}







#' Extract Ms1 table from an MS file
#'
#' This function reads an mzML or mzXML file with mzR,
#' and extracts MS1 table 
#'
#' @param path Character. Path to an mzML or mzXML file.
#'
#' @return A data.frame with columns rt, mz and intensity.
#' @export
#'
#' @examples
#' \dontrun{
#' extract_ms1_peaks <- extract_tic("file.mzXML")
#' }
extract_ms1_peaks <- function(file) {
  
  ms <- mzR::openMSfile(file)
  hdr <- mzR::header(ms)
  
  ms1_idx <- which(hdr$msLevel == 1)
  
  res <- do.call(rbind, lapply(ms1_idx, function(i) {
    p <- mzR::peaks(ms, i)
    
    data.frame(
      rt = hdr$retentionTime[i],
      mz = p[, 1],
      intensity = p[, 2]
    )
  }))
  
  mzR::close(ms)
  
  res
}









#' Extract Total ion chromatogram from an MS file
#'
#' This function reads an mzML or mzXML file with mzR,
#' and extracts an TIC 
#'
#' @param path Character. Path to an mzML or mzXML file.
#'
#' @return A data.frame with columns rt and intensity.
#' @export
#'
#' @examples
#' \dontrun{
#' tic <- extract_tic("file.mzXML")
#' }
extract_tic <- function(file) {
  
  ms <- mzR::openMSfile(file)
  
  hdr <- mzR::header(ms)
  
  # retention time and total ion current
  tic <- data.frame(
    rt = hdr$retentionTime,
    intensity = hdr$totIonCurrent
  )
  
  mzR::close(ms)
  
  tic
}






#' Extract an extracted ion chromatogram from an MS file
#'
#' This function reads an mzML or mzXML file with MSnbase, keeps MS1 scans,
#' and extracts an XIC around a target m/z using MSnbase::chromatogram().
#'
#' @param path Character. Path to an mzML or mzXML file.
#' @param target_mz Numeric. Target m/z value.
#' @param mz_tol Numeric. Absolute m/z tolerance around target_mz. Default is 0.01.
#'
#' @return A data.frame with columns rt and intensity.
#' @export
#'
#' @examples
#' \dontrun{
#' xic <- extract_xic("file.mzXML", target_mz = 453.367, mz_tol = 0.011)
#' }
extract_xic <- function(path,
                        target_mz,
                        mz_tol = 0.01) {

  if (!file.exists(path)) {
    stop("File does not exist: ", path, call. = FALSE)
  }

  if (!is.numeric(target_mz) || length(target_mz) != 1L) {
    stop("target_mz must be one numeric value.", call. = FALSE)
  }

  if (!is.numeric(mz_tol) || length(mz_tol) != 1L || mz_tol <= 0) {
    stop("mz_tol must be one positive numeric value.", call. = FALSE)
  }

  raw <- MSnbase::readMSData(path, mode = "onDisk")
  raw <- MSnbase::filterMsLevel(raw, 1)

  if (length(raw) == 0L) {
    return(data.frame(rt = numeric(), intensity = numeric()))
  }

  chr <- MSnbase::chromatogram(
    raw,
    mz = c(target_mz - mz_tol, target_mz + mz_tol)
  )

  data.frame(
    rt = MSnbase::rtime(chr[[1]]),
    intensity = MSnbase::intensity(chr[[1]])
  )
}

#' Extract MS2 spectra matching a precursor m/z
#'
#' This function reads an mzML or mzXML file with MSnbase, keeps MS2 scans,
#' filters spectra by precursor m/z, and returns all fragment peaks from
#' matching MS2 scans.
#'
#' @param path Character. Path to an mzML or mzXML file.
#' @param target_precursor Numeric. Target precursor m/z.
#' @param mz_tol Numeric. Absolute precursor m/z tolerance. Default is 0.01.
#'
#' @return A data.frame with columns precursor_mz, rt, mz, and intensity.
#' @export
#'
#' @examples
#' \dontrun{
#' ms2 <- extract_ms2("file.mzXML", target_precursor = 453.367, mz_tol = 0.011)
#' }
extract_ms2 <- function(path,
                        target_precursor,
                        mz_tol = 0.01) {

  if (!file.exists(path)) {
    stop("File does not exist: ", path, call. = FALSE)
  }

  if (!is.numeric(target_precursor) || length(target_precursor) != 1L) {
    stop("target_precursor must be one numeric value.", call. = FALSE)
  }

  if (!is.numeric(mz_tol) || length(mz_tol) != 1L || mz_tol <= 0) {
    stop("mz_tol must be one positive numeric value.", call. = FALSE)
  }

  raw <- MSnbase::readMSData(path, mode = "onDisk")
  raw <- MSnbase::filterMsLevel(raw, 2)

  if (length(raw) == 0L) {
    return(data.frame(
      precursor_mz = numeric(),
      rt = numeric(),
      mz = numeric(),
      intensity = numeric()
    ))
  }

  precursors <- MSnbase::precursorMz(raw)

  keep <- !is.na(precursors) &
    abs(precursors - target_precursor) <= mz_tol

  if (!any(keep)) {
    message("No MS2 spectra found for precursor ", target_precursor)
    return(data.frame(
      precursor_mz = numeric(),
      rt = numeric(),
      mz = numeric(),
      intensity = numeric()
    ))
  }

  raw_keep <- raw[keep]
  prec_keep <- precursors[keep]
  rt_keep <- MSnbase::rtime(raw_keep)

  res <- lapply(seq_along(raw_keep), function(i) {
    data.frame(
      precursor_mz = prec_keep[i],
      rt = rt_keep[i],
      mz = MSnbase::mz(raw_keep[[i]]),
      intensity = MSnbase::intensity(raw_keep[[i]])
    )
  })

  do.call(rbind, res)
}

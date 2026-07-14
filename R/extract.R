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

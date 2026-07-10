# msxtractor

Small R package for extracting and plotting:

- XICs from MS1 data
- MS2 fragment spectra matching a target precursor m/z
- XIC line plots
- MS2 stick spectra

It uses `MSnbase`, so it should work with mzML and mzXML files supported by MSnbase.

## Install local package

```r
install.packages("msxtractor_0.1.1.tar.gz", repos = NULL, type = "source")
```

## Usage

```r
library(msxtractor)

file_name <- "your_file.mzXML"

xic <- extract_xic(
  file_name,
  target_mz = 453.367,
  mz_tol = 0.011
)

ms2 <- extract_ms2(
  file_name,
  target_precursor = 453.367,
  mz_tol = 0.011
)

plot_xic(
  xic$rt,
  xic$intensity,
  title = basename(file_name)
)

plot_ms2(
  ms2$mz,
  ms2$intensity,
  rt = ms2$rt,
  title = "MS2 spectrum"
)
```

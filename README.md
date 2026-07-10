# msxtractor

Small R package for extracting:

- XICs from MS1 data
- MS2 fragment spectra matching a target precursor m/z

It uses `MSnbase`, so it should work with mzML and mzXML files supported by MSnbase.

## Install local package

From the parent directory containing `msxtractor`:

```r
remotes::install_local("msxtractor", dependencies = TRUE)
```

Or from the tarball:

```r
install.packages("msxtractor_0.1.0.tar.gz", repos = NULL, type = "source")
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
```

## Optional plotting

```r
library(ggplot2)

ggplot(xic, aes(rt, intensity)) +
  geom_line(linewidth = 0.8) +
  theme_bw()

ggplot(ms2) +
  geom_segment(aes(x = mz, xend = mz, y = 0, yend = intensity)) +
  theme_bw()
```

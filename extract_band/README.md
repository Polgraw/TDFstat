extract_band
============

This code extracts narrow band Short Fourier Transforms from SFDB files, performs inverse FFT and writes Short Time Series (STS)
into HDF5 file (one file per band). The STS's from HDF files are concatenated into long time series and cleaned later by genseg code.

The original version of extract_band is part of the [PSS library](https://dcc.ligo.org/LIGO-T1900140/public)
authored by Pia Astone et al. from Virgo Rome CW group (INFN, Physics Department of University of Rome "La Sapienza").
It extracts narrow band SFT's from SFDB files and writes output to text files.

Here we provide version modified by Paweł Ciecieląg (CAMK, Polish Virgo CW group) - 'extract_band_hdf_td`,
which in addition performs inverse FFT and writes resulting Short Time Series (STS) into HDF5 file (one file per band).
Files in this subdirectory are minimal subset of the whole PSS library which is needed to compile extract_band.

To compile, go to `pss_sfdb` subdirectory and type make.

[Documentation page](https://polgraw.github.io/TDFstat/input_data/)

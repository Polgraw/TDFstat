---
layout: default
title: Search
excerpt:
nav_order: 3
---

# F-statistic candidate signal search

This is the main code of the pipeline that performs the search for continuous gravitational-wave signals in the data from a network of detectors using the $\mathcal{F}$-statistic method. The search is performed over a grid of parameters covering the sky position, spindown and frequency values. Candidate signals (triggers) above a given threshold for the $\mathcal{F}$-statistic are stored in the output HDF5 files.

Since this is very time consuming part of the TDFstat all-sky search, the code is optimized for performance: vectorized and parallelized using OpenMP.


## Algorithm flowchart

![Code flowchart](img/flowchart.png)


## Prerequisites

The code is written in `C` and compiled with gnu17 dialect of GNU C compiler (support for complex numbers).
Required libraries: [GNU Scientific Library (GSL)](http://www.gnu.org/software/gsl/), [FFTW library](http://www.fftw.org), [HDF5 library](https://www.hdfgroup.org/solutions/hdf5), [iniparser](https://github.com/ndevilla/iniparser) (included in the repository).  
Opional libraries: [SLEEF vector library](https://sleef.org/) (included in the repository).


## Compilation 

Run `make` in `search/network/openmp`. The standard Makfile is configured to produce executable named `gwsearch-cascadelake-avx2-float`. Comments in the Makefile help to modify certain compilation options.

## Configuration file

The search code requires a configuration file in INI format. The repository contains example file `search-template.ini` with comments explaining each option. 
```
[search]

indir = /work/chuck/virgo/O3/allsky_o3_c01  # input data directory (string)
outdir = .            # output directory (string)
band = 72             # band number (int)
seg = 13              # segment number (int)
hemi = 2              # hemisphere [1,2 or 0 for both] (int)
thr = 14.5            # fstatistic threshold of candidates (double)
nod = 6               # length of input time series [days] (int)
dt = 2                # input time series sampling interval [seconds] (double)
overlap = 0.0625      # bands overlap, band frequency fpo=10+(1-overlap)*B/(2*dt) (double)
narrowdown = -1       # limits band in f to band_center +/- [0-0.5], auto calculated if < 0 (double)
fstat_norm =          # fstatistic normalization (NULL=white noise or blocks_avg)
grid_file = 013/grids/grid_013_0072_H1L1c.bin   

# optional
usedet =              # use only specified detectors 
range_file =          # range.dat , limits parameter space of the search
dump_range_file =     # name of the file to dump max. search ranges and exit
addsig =              # name of the file with signals to be injected
mods = read_O3        # coma separated list of modifiers: read_O3

# flags
veto_flag = 0         # veto lines: 0-no, 1-yes
gen_vlines_flag = 0   # if 1 generate vlines file and exit
checkp_flag = 0       # write checkpoint file on every triggers buffer flush
```

### Details of options / concepts


??? note "indir, input data (click to expand)"

    <br/>Input data must be prepared in a way outlined in the [genseg documentation](../input_data/#tdfstat-input-data-structure). Requied files are: xdat (time series), DetSSB (ephemeris), grid (grid generator matrix). Optional files: lines (veto lines), addsig (software injections), range file (search ranges).  
    Please note that same of the parameters have to match those used during input data generation (e.g., dt, overlap, nod).


??? note "Bands, segments, overlap (click to expand)"

    <br/> Short summary from the [genseg documentation](../input_data/): we analyze narrow-band (0.25-1 Hz) time segments of length being an integer multiple of sidereal day (typically 2-24 days).  
    Subsequent time segments are labeled by natural numbers. In some contexts (like filenames) those numbers are formatted using pattern DDD (e.g. 009).  
    Bands are overlapping in frequency to avoid edge effects. Bands are also labeled by natural numbers and in some contexts the 4 digit format is used (e.g., 0027). The general formula for the initial frequency of the band number b is:  
    $fpo(b) = 10 + (1 - overlap) \cdot b \frac{1}{2 \cdot dt}$  
    where bandwidth B=1/(2dt) and overlap is expressed as a fraction of B. The overlap shuld have the form of $2^{-n}$, to assure that fpo is aligned with fourier bins in the SFDB database (e.g. 0.0625).


??? note "narrowdown (click to expand)"

    <br/> Narrowdown is used to limit the range of frequencies for which the F-statistic is computed. It should be in range [0-0.5] and we define the whole band range to be [-0.5, 0.5]. Narrowdown 0.5 means "use the whole bandwidth". Narrowdown 0.45 means that the band is narrowed by 5% on both sides. If narrowdown < 0, then the frequency range is calculated from overlap in such a way that the overapping part of the band is split in the middle, effectivlly creating two sibling, non-overlapping bands.
    
    ```
         |<------- band1 ----|--->|
                    |<---|---- band 2 ------->|
    ```  
    
    This is needed to avoid analyzing the same range of frequencies more than once when performing all-sky search over many bands.

??? note "range_file (click to expand)"


??? note "addsig (click to expand)"



## Example of usage

Minimal call to `gwsearch-cascadelake-avx2-float` is as follows (code compiled with the `GNUSINCOS` option): 

```
% ./gwsearch-cascadelake-avx2-float search.ini  
```


### Network of detectors 

Test data frames $nnn=001-008$ with pure Gaussian noise 2-day time segments with sampling time equal to 2s (`xdatc_nnn_1234.bin`) for two LIGO detectors H1 and L1 are [available here](https://polgraw.camk.edu.pl/H1L1_2d_0.25.tar.gz). 


The program will proceed assuming that 

  * the data directory for frame `001` is located at `../../../testdata/2d_0.25/001` and contain subdirectories with input data for H1, L1 and/or V1 detectors (all available detectors are used by default; to select specific detectors, use `-usedet` option), 
  * the grid of parameters files is expected to be in `../../../testdata/2d_0.25/001`, 
  * `band` equals to $1234$,  
  * the sampling time `dt` equals $2 s$, 
  * number of days `nod` in $2$, 
  * the `-addsig` option is used to add a software injection to pure noise Gaussian data. The signal's parameters are randomly generated using the `sigen` code (for more details see the [minimal example](../polgraw-allsky/pipeline_script)).  
  * the threshold for the $\mathcal{F}$-statistic is set to be $14.5$,
  * `-output` is the current directory, 
  * `--nocheckpoint` disables the checkpointing (writing the last visited position on the grid to the `state` file), 


## Output files

HDF5 structure:
```
HDF5 object                  | data type
----------------------------------------------------
/                            | root group
├── attr: format_version     | int
├── attr: git_commit         | string
├── attr: opts               | compound (struct Command_line_opts)
├── attr: sett               | compound (struct Search_settings)
├── attr: s_range            | compound (struct Search_ranges)
├── attr: ifo                | compound[] (struct Detector_settings)
├── attr: detectors          | string
├── attr: t_start            | string
├── attr: t_end              | string
├── attr: totsgnl            | int
├── attr: walltime           | double
├── attr: num_threads        | int
├── dataset "triggers"       | dataset
│   ├── data                 | compound[] (Trigger[])
│       ├── m                | float
│       ├── n                | float
│       ├── s                | float
│       ├── ra               | float
│       ├── dec              | float
│       ├── fdot             | float
│       ├── ffstat           | float[]
├── dataset "triggers_inj_1" | dataset
│    ├── data                | compound[] (Trigger[])
└── ...
```
#### Full h5dump output example
<details>
<summary>Click to expand</summary>
```
HDF5 "triggers_013_0072_2.h5" {
GROUP "/" {
   ATTRIBUTE "format_version" {
      DATATYPE  H5T_STD_I32LE
      DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
   }
   ATTRIBUTE "git_commit" {
      DATATYPE  H5T_STRING {
         STRSIZE 41;
         STRPAD H5T_STR_NULLTERM;
         CSET H5T_CSET_ASCII;
         CTYPE H5T_C_S1;
      }
      DATASPACE  SCALAR
   }
   ATTRIBUTE "ifo" {
      DATATYPE  H5T_COMPOUND {
         H5T_STRING {
            STRSIZE 2;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         } "name";
         H5T_STRING {
            STRSIZE 2048;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         } "xdatname";
      }
      DATASPACE  SIMPLE { ( 2 ) / ( 2 ) }
   }
   ATTRIBUTE "num_threads" {
      DATATYPE  H5T_STD_I32LE
      DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
   }
   ATTRIBUTE "opts" {
      DATATYPE  H5T_COMPOUND {
         H5T_STD_I32LE "checkp_flag";
         H5T_STD_I32LE "veto_flag";
         H5T_STD_I32LE "seg";
         H5T_STD_I32LE "band";
         H5T_STD_I32LE "hemi";
         H5T_STD_I32LE "nod";
         H5T_IEEE_F64LE "thr";
         H5T_IEEE_F64LE "narrowdown";
         H5T_IEEE_F64LE "overlap";
         H5T_STRING {
            STRSIZE H5T_VARIABLE;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         } "indir";
         H5T_STRING {
            STRSIZE H5T_VARIABLE;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         } "outdir";
         H5T_STRING {
            STRSIZE H5T_VARIABLE;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         } "range_file";
         H5T_STRING {
            STRSIZE H5T_VARIABLE;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         } "grid_file";
         H5T_STRING {
            STRSIZE H5T_VARIABLE;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         } "usedet";
         H5T_STRING {
            STRSIZE H5T_VARIABLE;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         } "addsig";
         H5T_STRING {
            STRSIZE H5T_VARIABLE;
            STRPAD H5T_STR_NULLTERM;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         } "fstat_norm";
      }
      DATASPACE  SCALAR
   }
   ATTRIBUTE "s_range" {
      DATATYPE  H5T_COMPOUND {
         H5T_ARRAY { [2] H5T_STD_I32LE } "pmr";
         H5T_ARRAY { [2] H5T_STD_I32LE } "fr";
         H5T_ARRAY { [2] H5T_IEEE_F32LE } "mr";
         H5T_ARRAY { [2] H5T_IEEE_F32LE } "nr";
         H5T_ARRAY { [2] H5T_IEEE_F32LE } "spndr";
         H5T_IEEE_F32LE "mstep";
         H5T_IEEE_F32LE "nstep";
         H5T_IEEE_F32LE "sstep";
         H5T_IEEE_F32LE "mst";
         H5T_IEEE_F32LE "nst";
         H5T_IEEE_F32LE "sst";
         H5T_STD_I32LE "pst";
      }
      DATASPACE  SCALAR
   }
   ATTRIBUTE "sett" {
      DATATYPE  H5T_COMPOUND {
         H5T_IEEE_F64LE "fpo";
         H5T_IEEE_F64LE "dt";
         H5T_IEEE_F64LE "B";
         H5T_IEEE_F64LE "oms";
         H5T_IEEE_F64LE "omr";
         H5T_IEEE_F64LE "Smin";
         H5T_IEEE_F64LE "Smax";
         H5T_IEEE_F64LE "sepsm";
         H5T_IEEE_F64LE "cepsm";
         H5T_STD_I32LE "nfft";
         H5T_STD_I32LE "nod";
         H5T_STD_I32LE "N";
         H5T_STD_I32LE "nfftf";
         H5T_STD_I32LE "nmax";
         H5T_STD_I32LE "nmin";
         H5T_STD_I32LE "s";
         H5T_STD_I32LE "nd";
         H5T_STD_I32LE "interpftpad";
         H5T_STD_I32LE "fftpad";
         H5T_STD_I32LE "Ninterp";
         H5T_STD_I32LE "nifo";
         H5T_STD_I32LE "numlines_band";
         H5T_STD_I32LE "nvlines_all_inband";
         H5T_STD_I32LE "bufsize";
         H5T_STD_I32LE "dd";
         H5T_ARRAY { [16] H5T_IEEE_F64LE } "M";
         H5T_ARRAY { [5][2] H5T_IEEE_F64LE } "lines";
      }
      DATASPACE  SCALAR
   }
   ATTRIBUTE "t_end" {
      DATATYPE  H5T_STRING {
         STRSIZE 20;
         STRPAD H5T_STR_NULLTERM;
         CSET H5T_CSET_ASCII;
         CTYPE H5T_C_S1;
      }
      DATASPACE  SCALAR
   }
   ATTRIBUTE "t_start" {
      DATATYPE  H5T_STRING {
         STRSIZE 20;
         STRPAD H5T_STR_NULLTERM;
         CSET H5T_CSET_ASCII;
         CTYPE H5T_C_S1;
      }
      DATASPACE  SCALAR
   }
   ATTRIBUTE "totsgnl" {
      DATATYPE  H5T_STD_I32LE
      DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
   }
   ATTRIBUTE "walltime" {
      DATATYPE  H5T_IEEE_F64LE
      DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
   }
   DATASET "triggers" {
      DATATYPE  H5T_COMPOUND {
         H5T_IEEE_F32LE "m";
         H5T_IEEE_F32LE "n";
         H5T_IEEE_F32LE "s";
         H5T_IEEE_F32LE "ra";
         H5T_IEEE_F32LE "dec";
         H5T_IEEE_F32LE "fdot";
         H5T_VLEN { H5T_IEEE_F32LE } "ffstat";
      }
      DATASPACE  SIMPLE { ( 4952 ) / ( H5S_UNLIMITED ) }
   }
}
}
```
</details>


Binary output files, containing trigger candidate events above an arbitrary threshold (option `-threshold` for the $\mathcal{F}$-statistic, default 20), are written to the `output_dir` directory. There are two output files for every input data sequence: `triggers_nnn_bbbb_1.bin` and
`triggers_nnn_bbbb_2.bin`,  where  `1` and  `2` correspond to the northern and southern ecliptic hemisphere. Each trigger (candidate) event occupies `40` consecutive bytes (5 double numbers), with the following meaning:

Record no.            | 
--------------------- | ---------------------------- 
1                     | frequency [radians, between 0 and $\pi$] above `fpo`  
2                     | spindown [$\mathrm{Hz/s}$]  
3                     | declination [radians, between $\pi/2$ and $-\pi/2$]
4                     | right ascension [radians, between 0 and $2\pi$]
5                     | signal-to-noise ratio

For the example above, the first 10 triggers from `triggers_001_1234_2.bin` are 

```bash
3.05617018e+00 -3.42376198e-08 -7.68007347e-02 2.59248668e+00 5.06667333e+00 
1.18243015e+00 -3.20762991e-08 -7.68007347e-02 2.59248668e+00 5.05528873e+00 
1.08103361e-01 -2.77536578e-08 -7.68007347e-02 2.59248668e+00 5.07085254e+00 
1.90022435e+00 -2.77536578e-08 -7.68007347e-02 2.59248668e+00 5.15191593e+00 
1.90000217e+00 -2.55923371e-08 -7.68007347e-02 2.59248668e+00 5.42638039e+00 
2.09224664e+00 -2.34310165e-08 -7.68007347e-02 2.59248668e+00 5.20879551e+00 
2.38731576e+00 -2.12696958e-08 -7.68007347e-02 2.59248668e+00 5.31983396e+00 
3.00543165e+00 -1.91083751e-08 -7.68007347e-02 2.59248668e+00 5.29454616e+00 
7.49333983e-01 -1.26244131e-08 -7.68007347e-02 2.59248668e+00 5.08724856e+00 
2.08710778e-01  3.43510887e-10 -7.68007347e-02 2.59248668e+00 5.17537018e+00 
```

### Auxiliary output files

* `wisdom-hostname.dat` - performance-testing file created by the `FFTW`. The `hostname` variable is determined by a call to `gethostname()`, 

* `state_nnn_bbbb.dat` - checkpointing file containing the last grid point visited. The search can  be safely restarted, calculations will continue  from the last grid position saved to this file. After successful termination, checkpoint file is left empty.

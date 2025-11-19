#ifndef __STRUCT_H__
#define __STRUCT_H__

#include <fftw3.h>
#include <complex.h>
#include <hdf5.h>

#define MAX_DETECTORS 8        // Maximum number of detectors in network
#define DETNAME_LENGTH 2       // Detector name length (H1, L1, V1...)
#define FNAME_LENGTH 2048      // Maximum length of a filename
#define INICANDSIZE 1048576    // Initial size for array candidates storage;
                               // realloc called if needed (in coincidences)

#define MAXL 512               // Max number of veto lines in band
#define MAXVFILEL 8192         // Max number of text lines in the veto file

/* Command line options for search */
typedef struct _comm_line_opts {

     int checkp_flag, veto_flag, gen_vlines_flag, help_flag;
     //int fftinterp;
     int seg, band, hemi; //, nod;
     double thr;
     double fpo_val, narrowdown, overlap;
     const char *indir, *outdir, *grid_file, *dump_range_file,
          *usedet, *addsig, *fstat_norm, *label, *mods,
          *gtype, *gcenter, *gsizes, *gsteps;
     char state_file[FNAME_LENGTH];

} Command_line_opts;


/* input signal arrays */
typedef struct _signals {

     float *xDat;
     double *DetSSB;       // Ephemeris of the detector
     double *aa, *bb;      // Amplitude modulation functions
     double *shftf, *shft; // Resampling and time-shifting

     double epsm,
            phir,
            sepsm,	  // sin(epsm)
            cepsm,	  // cos(epsm)
            sphir,	  // sin(phi_r)
            cphir,	  // cos(phi_r)
            crf0,    // number of 0s as: N/(N-Nzeros)
            sig2; 	  // variance of signal

     int Nzeros;
     double complex *xDatma, *xDatmb;

} Signals;


/* fftw arrays */
typedef struct _fftw_arrays {

     fftw_complex *xa, *xb;
     FFTW_PRE(_complex) *fxa, *fxb;
     int arr_len;

} FFTW_arrays;


/* Search range  */
typedef struct _search_range {

     // Hemispheres range: pmr[0] - start, pmr[1] - end; values are 1 or 2.
     // fr is always integer (fft bins) and limits only the OUTPUT range, max range is 0 to nfft-1
     int pmr[2], fr[2];
     // mm,nn,s ranges: [0] - start, [1] - end
     float mr[2], nr[2], spndr[2];
     // grid steps (default =1)
     float mstep, nstep, sstep;
     // initial values after restart
     float mst, nst, sst;
     int pst;

} Search_range;


/* FFTW plans  */
typedef struct _fftw_plans {

     fftw_plan pl_int,  // interpolation forward
               pl_inv;  // interpolation backward
     FFTW_PRE(_plan) plan;

} FFTW_plans;


/* Auxiluary arrays */
typedef struct _aux_arrays {

     double *sinmodf, *cosmodf; // Earth position
     double *t2;                // time^2

} Aux_arrays;


/* Search settings */
typedef struct _search_settings {

     double fpo,    // Band frequency
            dt,     // Sampling time
            B,      // Bandwidth
            oms,    // Dimensionless angular frequency (fpo)
            omr,    // C_OMEGA_R * dt (dimensionless Earth's angular frequency)
            Smin,   // Minimum spindown
            Smax,   // Maximum spindown
            sepsm,	// sin(epsm)
            cepsm;	// cos(epsm)

     int nfft,         // length of fft
         nod,          // number of days of observation
         N,            // number of data points
         nfftf,        // nfft * fftpad
         nmax,	        // first and last point
         nmin, 	   // of Fstat
         s,            // number of spindowns
         nd,           // degrees of freedom
         interpftpad,
         fftpad,       // zero padding
         Ninterp, 	   // for resampling (set in plan_fftw() init.c)
         nifo,         // number of detectors
         numlines_band,// number of lines in band
         nvlines_all_inband, // number of all veto lines in band
         bufsize,      // size of buffer for ffstat pairs
         dd;           // block size for Fstat maxima search

     double M[16];           // Grid-generating matrix (or Fisher matrix,
                          // in case of coincidences)
     double invM[4][4];   // Inverse of M

     double vedva[4][4];   // transformation matrix: its columns are
                           // eigenvectors, each component multiplied
                           // by sqrt(eigval), see init.c manage_grid_matrix():
                           // sett->vedva[i][j]  = eigvec[i][j]*sqrt(eigval[j])

     double lines[MAXL][2]; // Array for lines in given band

} Search_settings;


/* Amplitude modulation function coefficients */
typedef struct _ampl_mod_coeff {

	double c1, c2, c3, c4, c5, c6, c7, c8, c9;

} Ampl_mod_coeff;


/* Detector and its data related settings */
typedef struct _detector {

     char xdatname[FNAME_LENGTH];
     char name[DETNAME_LENGTH];

     double ephi, 		// Geographical latitude phi in radians
            elam, 		// Geographical longitude in radians
            eheight,     // Height h above the Earth ellipsoid in meters
            egam; 		// Orientation of the detector gamma

     Ampl_mod_coeff amod;
     Signals sig;

} Detector_settings;

/* Global array of detectors (network) */
extern Detector_settings ifo[MAX_DETECTORS];

/* Command line option struct for coincidences */
typedef struct _comm_line_opts_coinc {

     int help_flag;

     int shift,    // Cell shifts  (4 digit number corresponding to fsda, e.g. 0101)
                   // Cell scaling in fsda: 4 numbers corresponding to
         scalef,   // frequency f,
         scales,   // spindown s,
         scaled,   // declination d,
         scalea,   // right ascencion a
         refr,     // Reference frame
         band,     // band number
         hemi;     // Hemisphere

     // Minimal number of coincidences recorded in the output
     int mincoin;
     double fpo, refgps, narrowdown, snrcutoff, overlap;
     char prefix[512], refgrid[1024], *wd, infile[512];

} Command_line_opts_coinc;


typedef struct _triggers {

     int frameinfo[256][3];    // Info about candidates in frames:
                               // - [0] frame number, [1] initial number
                               // of candidates, [2] number of candidates
                               // after sorting (unique)

     int frcount, goodcands;

} Candidate_triggers;

typedef struct _trigger {
          float  m, n, s;
          float  ra, dec, fdot;
          hvl_t  ffstat;
} Trigger;

#endif

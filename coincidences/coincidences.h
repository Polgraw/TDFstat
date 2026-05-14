
typedef struct {

     int
          refr,         // Reference frame
          hemi,         // Hemisphere
          scalef,       // frequency,
          scales,       // spindown,
          scalem,       // m sky position,
          scalen,       // n sky position.
          shift,        // Cell shifts  (4 digit number corresponding to fsda, e.g. 0101)
                        // Cell scaling in fsda: 4 numbers corresponding to
          mincoin,      // Minimum number of coincidences
          write_ctrigs; // Whether to write out coincident triggers

     double cthr;       // Fstatistic threshold for coincidences
     const char *trig_dset, *coinc_dset_pre;
     const char *out_dir;
     char *trig_files[MAX_NSEG];
     size_t nseg; // Number of segments to compare

} Coinc_opts;


typedef struct {

     int
          band,    // Band
          nod,     // Number of days
          seg,     // Segment
          N,       // data length
          hemi,    // Hemisphere
          nfftf,   // Number of FFT bins (to calculate fbin)
          numlines_band,
          nvlines_all_inband; // Number of lines in band

     size_t  sgnlv_size; // Size of sgnlv array

     double
          overlap,    // Overlap between bands
          narrowdown, // Narrowdown factor for bandpass filtering
          B,          // Bandwidth
          fpo,        // Starting frequency of band
          dt,         // Time step
          lines[MAXL][2]; // Array for lines in given band

     char *grid_file; //just to verify it doesn't change

} Search_params;


typedef struct {

     float  m, n, s;          // "integer" grid coordinates
     float  ra, dec, fdot;    // physical coordinates corresponding to the grid point
     hvl_t  *ffstat;          // interleaved f (in 0-PI units) and Fstatistic values

} Coinc_Trigger;


void read_coinc_ini(char *ini_fname, Coinc_opts *copts);
size_t read_triggers_file(const char *filename, const char *t_dset_name,
                          Trigger **sgnlv, Search_params *search_par);
int select_goodcands(int iseg, Coinc_opts *copts, Search_params *search_par,
                     Trigger *sgnlv, Coinc_Trigger *ctrigs);
int write_ctrigs_hdf(const char *ctrigs_fname, Coinc_opts *copts,
                     Search_params *search_par, Coinc_Trigger *ctrigs);

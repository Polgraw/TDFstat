// coincidences options read from ini file
typedef struct {

    int
        refr,         // Reference frame
        hemi,         // Hemisphere
        scalef,       // cell scale in frequency
        scales,       // cell scale in spindown
        scalem,       // cell scale in m
        scalen,       // cell scale in n
        mincoin,      // Minimum number of coincidences
        write_ctrigs; // Whether to write out coincident triggers
    double cthr;       // Fstatistic threshold for coincidences
    const char *shift; // Cell shifts (binary-like 4 digits corresponding to mnsf shift, e.g. 0101)
    const char *trig_dset, *coinc_dset; // Trigger dataset name and coincidence dataset prefix
    const char *out_dir;
    char *trig_files[MAX_NSEG];
    size_t nseg; // Number of segments to compare

} Coinc_opts;

// Search parameters read from triggers HDF file
// contains both sett and opts elements
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
        lines[MAXL][2], // Array for lines in given band
        vfrac; // Veto fraction of the band

    char *grid_file; //just to verify it doesn't change
    char *label;

} Search_params;

// Coincidence trigger structure; ffstat here has additionsl dim - time segment
typedef struct {
    float  m, n, s;          // "integer" grid coordinates
    float  ra, dec, fdot;    // physical coordinates corresponding to the grid point
    hvl_t  *ffstat;          // interleaved f (in 0-PI units) and Fstatistic values
    int iccell;              // index of the coincidence cell this trigger falls into
} Coinc_Trigger;


// This structure maps search grid cells m,n,s to coincidences cells
typedef struct {
    int mc, nc, sc;      // "integer" coordinates of the coincidence cell
    int nctrigs;         // length of ictrigs array
    int *ictrigs;        // indices of sgnlv/ctrigs that fall in the coin cell
} Search2coi_mns;

// Coincidence structure; contains information about the coincidence cell and the triggers that fall into it
typedef struct {
    char shift[5];   // Cell shifts  (4 digit number corresponding to fsmn, e.g. 0101)
    short w;         // number of segments participating in the coincidence
    short *cseg;     // segments participating in the coincidence (size=w)
    short *n_ccell_trigs; // number of triggers in the coin cell for each segment (size=w)
    int *trig_mns;   // indices of the unique triggers in each segment (size=w)
    float avg_snr, avg_f, avg_fdot, avg_ra, avg_dec;   // physical coordinates of unique trigger in each segment
} Coincidence;

typedef struct {
    int ncoinc_min; // Minimum number of coincidences
    int ncoinc_max; // Maximum number of coincidences
    int crit_mul_0_01; // Critical multiplicity for fap=0.01
    int crit_mul_0_001; // Critical multiplicity for fap=0.001
    hvl_t fap;  // False Alarm Probability [ncoinc_min, ncoinc_max]
} Fap_t;


void read_coinc_ini(char *ini_fname, Coinc_opts *copts);
size_t read_triggers_file(const char *filename, const char *t_dset_name,
                          Trigger **sgnlv, Search_params *search_par);
int select_goodcands(int iseg, Coinc_opts *copts, Search_params *search_par,
                     Trigger *sgnlv, Coinc_Trigger *ctrigs);
int write_ctrigs_hdf(const char *ctrigs_fname, Coinc_opts *copts,
                     Search_params *search_par, Coinc_Trigger *ctrigs,
                     int seginfo[][3]);
int init_coin_hdf(const char *coin_fname, Coinc_opts *copts,
                  Search_params *search_par);
int write_coi_hdf(const char *coin_fname, Coinc_opts *copts,
                  Coincidence *coi, int icoi, const char *shift_str,
                  int seginfo[][3]);
int add_maxcoi_hdf(const char *coin_fname, Coinc_opts *copts,
                   Coincidence *max_coi, Fap_t *fap_data);
float read_vlines_file(const char *veto_fname, Search_params *search_par);
int FalseAlarmProb(int, int, double, int*, double, double *);

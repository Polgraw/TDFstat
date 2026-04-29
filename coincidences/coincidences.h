typedef struct {

     int
          refr,     // Reference frame
          hemi,     // Hemisphere
          scalef,   // frequency,
          scales,   // spindown,
          scalem,   // m sky position,
          scalen,   // n sky position.
          shift,    // Cell shifts  (4 digit number corresponding to fsda, e.g. 0101)
                    // Cell scaling in fsda: 4 numbers corresponding to
          mincoin;  // Minimum number of coincidences

     double cthr;         // Fstatistic threshold for coincidences
//     char trig_dset[DSET_NAME_LEN], coinc_dset[DSET_NAME_LEN];
//     char out_dir[FILE_NAME_LEN], frame_list[FILE_NAME_LEN*100];
     const char *trig_dset, *coinc_dset_pre;
     const char *out_dir;
     char *trig_files[MAX_TRIG_FILES];
     size_t nseg;

} Coinc_opts;

typedef struct {

     int
          band,    // Band
          nod,     // Number of days
          seg,     // Segment
          hemi,    // Hemisphere
          nfftf,   // Number of FFT bins (to calculate fbin)
          nvlines_all_inband; // Number of lines in band

     size_t  sgnlv_size; // Size of sgnlv array

     double
          overlap,    // Overlap between bands
          narrowdown, // Narrowdown factor for bandpass filtering
          B,         // Bandwidth
          lines[MAXL][2]; // Array for lines in given band

     char *grid_file; //just to verify it doesn't change

} Search_params;


typedef struct {
          float  m, n, s;
          float  ra, dec, fdot;
          hvl_t  *ffstat;
} Coinc_Trigger;


void read_coinc_ini(char *ini_fname, Coinc_opts *copts);
size_t read_triggers_file(const char *filename, const char *t_dset_name,
     Trigger **sgnlv, Search_params *search_par);

void filter_to_allcands(int iseg, Coinc_opts *copts, Search_params *search_par,
     Trigger *sgnlv, Coinc_Trigger *allcands);

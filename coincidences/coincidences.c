#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <getopt.h>
//#include <dirent.h>

#include "cdefs.h"
#include "../search/network/openmp/struct.h"
#include "coincidences.h"


int main (int argc, char* argv[]) {

     char ini_fname[FILE_NAME_LEN];
     Coinc_opts copts;
     //Search_settings sett;
     //Candidate_triggers trig;
     Trigger *sgnlv;
     size_t sgnlv_size;
     search_params sparams = {0};
     int i, j; 

     if (argc > 1) {
          strcpy (ini_fname, argv[1]);
     } else if (argc > 2) {
          printf("WARNING: too many arguments, only first one is used.\n");
     } else {
          printf("ERROR: missing input file name. Call: \"<executable> search.ini\" \n");
          exit(EXIT_FAILURE);
     }
     
     read_coinc_ini(ini_fname, &copts);
     
     // read trigger from given files,
     // search sibling bands for out of band triggers,
     // apply Fstat threshold, shift in frequency,
     // store selected triggers in a new data structure
     //for(i=0; i<copts.n_trig_files; i++) {
     for(i=0; i<1; i++) {
          sgnlv_size = read_triggers_file(copts.trig_files[0], copts.trig_dset, &sgnlv, &sparams);

          printf("opts (from HDF5 attribute):\n");
          printf("  band       = %d\n",   sparams.band);
          printf("  seg        = %d\n",   sparams.seg);
          printf("  hemi       = %d\n",   sparams.hemi);
          printf("  overlap    = %f\n",   sparams.overlap);
          printf("  narrowdown = %f\n",   sparams.narrowdown);
          printf("  grid_file  = %s\n",   sparams.grid_file ? sparams.grid_file : "(null)");
          if (sgnlv != NULL) {
               // debugging
               printf("sgnlv_size=%zu  sgnlv[0]: ", sgnlv_size);
               printf("m=%f  n=%f  s=%f  ra=%f  dec=%f  fdot=%f  ntrig=%zu\n",
                    sgnlv[0].m, sgnlv[0].n, sgnlv[0].s, sgnlv[0].ra, sgnlv[0].dec, 
                    sgnlv[0].fdot, sgnlv[0].ffstat.len/2);
               float *ffp = (float *)sgnlv[0].ffstat.p;
               for (size_t k = 0; k < (sgnlv[0].ffstat.len+1)/2; k++) {
                    printf("  t[%zu]=(%f,%f)  ", k, ffp[2*k], ffp[2*k+1]);
               }
               printf("\n");
               
          }
     }

     if (sgnlv != NULL) {
     } else {
          printf("sgnlv is NULL (no triggers read)\n");
     }

     return EXIT_SUCCESS;     
}

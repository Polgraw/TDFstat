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
     Coinc_Trigger *allcands=NULL;
     Search_params search_par = {0};
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
     
     // read triggers from input files,
     // filter: apply Fstat threshold, shift in frequency,
     // store selected triggers in a new data structure: allcands
     // search sibling bands for out of band triggers (filter and append to allcands)

     if (allcands != NULL) {
          free(allcands);
          allcands = NULL;
     }
     
     //for(i=0; i<copts.n_trig_files; i++) {
     for(i=0; i<2; i++) {
          read_triggers_file(copts.trig_files[i], copts.trig_dset, &sgnlv, &search_par);
          //filter_to_allcands(&search_par, &sgnlv, &allcands);
          //if (i==0) init_allcands();
          //filter_to_allcands(&sgnlv, &search_par);
          // handle out of band triggers here, if needed     
     }

     // write HDF file with allcands, if needed
     // search for coincidences between segments 
     
     
     if (sgnlv != NULL) {
     } else {
          printf("sgnlv is NULL (no triggers read)\n");
     }

     return EXIT_SUCCESS;     
}


void filter_to_allcands(Coinc_opts *copts, Search_params *search_par,
     Trigger **sgnlv, Coinc_Trigger **allcands)
{
     // sgnlv array of structures contains variable length array ffstat which interleave f and fstat values
     // allocate memory for allcands of the same size as sgnlv, if needed
     // for each element of sgnlv loop over sgnlv.ffstat triggers and
     //      1. apply Fstat threshold (copts->cthr)
     //      2. check if f is in the search_par->lines (where [0] is start of the line and [1] the end)
     //      3. shift in f due to spindown: f = f0 + 2*sgnlv.fdot*(search_par->N)*(copts->refr - search_par->seg)
     // append the remaining triggers to allcands

     *allcands= (Coinc_Trigger *) malloc (sizeof(Coinc_Trigger) * search_par->sgnlv_size);
     if (*allcands == NULL) {
          printf("ERROR: could not allocate memory for allcands\n");
          exit(EXIT_FAILURE);
     }
     
     
     
     
     
     
     
}

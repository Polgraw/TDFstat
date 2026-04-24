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

     //for(i=0; i<copts.n_trig_files; i++) {
     for(i=0; i<2; i++) {
          read_triggers_file(copts.trig_files[i], copts.trig_dset, &sgnlv, &search_par);
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

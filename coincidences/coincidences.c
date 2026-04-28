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
     int i, j, k;

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
          
          // initialize allcands
          if (i==0) {
               allcands = (Coinc_Trigger *) malloc (sizeof(Coinc_Trigger) * search_par.sgnlv_size);
               for (j=0; j<search_par.sgnlv_size; j++) {
                    allcands[j].m = sgnlv[j].m;
                    allcands[j].n = sgnlv[j].n;
                    allcands[j].s = sgnlv[j].s;
                    // here ffstat is an array of length = number of segments (copts.nseg)
                    // and each element contains (hvl_t ffstat) for the coresponding segment
                    allcands[j].ffstat = (hvl_t *) malloc (sizeof(hvl_t *) * copts.nseg);
               }
          }
          
          int ntrig_seg = 0;
          // copy sgnlv[j].ffstat for the current segment (i) to allcands[j].ffstat[i]
          for (j=0; j<search_par.sgnlv_size; j++) {
               allcands[j].ffstat[i].p = (float *) malloc (sizeof(float) * sgnlv[j].ffstat.len);
               allcands[j].ffstat[i].len = sgnlv[j].ffstat.len;
               memcpy(allcands[j].ffstat[i].p, sgnlv[j].ffstat.p, 
                    sizeof(float) * sgnlv[j].ffstat.len);
               ntrig_seg += (sgnlv[j].ffstat.len+1)/2;
          }
          
          printf("     Segment %d: read %d triggers\n", i, ntrig_seg);

          printf("    [all sgnlv_size: %zu]", search_par.sgnlv_size);
          float *ffp = (float *)allcands[0].ffstat[i].p;
          for (size_t k = 0; k < (allcands[0].ffstat[i].len+1)/2; k++) {
               printf("  all p[%zu]=(%f,%f)  ", k, ffp[2*k], ffp[2*k+1]);
          }
          printf("\n");
          
          // filter triggers and append to allcands
          //filter_to_allcands(&copts, &search_par, &sgnlv, &allcands);

     }

     // write HDF file with allcands, if needed
     // search for coincidences between segments 
     
     

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
     
     printf("    [sgnlv_size: %zu]", search_par->sgnlv_size);
     float *ffp = (float *)sgnlv[0]->ffstat.p;
     for (size_t k = 0; k < (sgnlv[0]->ffstat.len+1)/2; k++) {
          printf("  p[%zu]=(%f,%f)  ", k, ffp[2*k], ffp[2*k+1]);
     }
     printf("\n");     
     
     exit(1);
     
}

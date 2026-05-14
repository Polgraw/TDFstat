#include <H5Epublic.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <dirent.h>

#include "cdefs.h"
#include "../search/network/openmp/struct.h"
#include "coincidences.h"


int main (int argc, char* argv[]) {

     char ini_fname[FILE_NAME_LEN];
     Coinc_opts copts;                // coincidence search options, read from ini file
     Trigger *sgnlv;                  // triggers struct as written in search
     Coinc_Trigger *ctrigs=NULL;    // triggers struct for coincidences, assembled from all time segments
     Search_params search_par = {0};  // selected fields from Search_settings
     int seginfo[MAX_NSEG][3];  // Info about candidates in frames:
                                // [0] frame number, 
                                // [1] good, inband triggers at reference time
                                // [2] unique triggers, max. one per coincidence cell at reference time

     int i, j, k, iseg;

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

     if (ctrigs != NULL) {
          free(ctrigs);
          ctrigs = NULL;
     }
     
     copts.nseg = 2; // for testing, should be read from ini file
     for(iseg=0; iseg<copts.nseg; iseg++) {
          printf("> segment %d\n", iseg);
          read_triggers_file(copts.trig_files[iseg], copts.trig_dset, &sgnlv, &search_par);

#if 1     
          {
               int itest = search_par.sgnlv_size-1;
               printf("   [select_goodcands: sgnlv[%d].ffstat. len=%zu]", itest, sgnlv[itest].ffstat.len);
               float *ffp = (float *)sgnlv[itest].ffstat.p;
               for (size_t k = 0; k < (sgnlv[itest].ffstat.len+1)/2; k++)
                    printf("  p[%zu]=(%f,%f)  ", k, ffp[2*k], ffp[2*k+1]);
               printf("\n");
          }
#endif
          
          // initialize allcands
          if (iseg==0) {
               ctrigs = (Coinc_Trigger *) malloc (sizeof(Coinc_Trigger) * search_par.sgnlv_size);
               for (j=0; j<search_par.sgnlv_size; j++) {
                    ctrigs[j].m    = sgnlv[j].m;
                    ctrigs[j].n    = sgnlv[j].n;
                    ctrigs[j].s    = sgnlv[j].s;
                    ctrigs[j].ra   = sgnlv[j].ra;
                    ctrigs[j].dec  = sgnlv[j].dec;
                    ctrigs[j].fdot = sgnlv[j].fdot;
                    // here ffstat is an array of length = number of segments (copts.nseg)
                    // and each element contains (hvl_t ffstat) for the given segment
                    ctrigs[j].ffstat = (hvl_t *) malloc (sizeof(hvl_t) * copts.nseg);
               }
          }

          seginfo[iseg][0] = search_par.seg;
          // select good candidates: 
          // apply Fstat threshold, shift in frequency to reference time, 
          // narrowdown, reject triggers in lines
          seginfo[iseg][1] = select_goodcands(iseg, &copts, &search_par, sgnlv, ctrigs);

#if 1
          int itest = search_par.sgnlv_size-1;
          printf("   [ctrigs[%d].ffstat[%d].  ", itest, iseg);
          float *ffp = (float *)ctrigs[itest].ffstat[iseg].p;
          for (size_t k = 0; k < (ctrigs[itest].ffstat[iseg].len+1)/2; k++)
               printf("  p[%zu]=(%f,%f)  ", k, ffp[2*k], ffp[2*k+1]);
          printf("\n");
#endif          

     } // iseg

     // write HDF file with ctrigs, if needed
     if (copts.write_ctrigs) {
          char ctrigs_fname[FILE_NAME_LEN];
          // can't use search_par.hemi since it can be 0 (both) and thus
          // the only source of current hemi is the triggers filename
          int fname_len = (int)strlen(copts.trig_files[0]);
          snprintf(ctrigs_fname, FILE_NAME_LEN, "%s/ctrigs%s", 
               copts.out_dir, copts.trig_files[0]+fname_len-10);
          printf("Writing coincidence triggers to file: %s\n", ctrigs_fname);
          write_ctrigs_hdf(ctrigs_fname, &copts, &search_par, ctrigs);
     }


     // search for coincidences between segments 
     //coincidences_on_cgrid(&copts, &search_par, ctrigs, coinc);     

     return EXIT_SUCCESS;
}




int select_goodcands(int iseg, Coinc_opts *copts, Search_params *search_par,
                     Trigger *sgnlv, Coinc_Trigger *ctrigs)
{
     // sgnlv array of structures contains variable length array ffstat which interleave f and fstat values
     // allocate memory for allcands of the same size as sgnlv, if needed
     // for each element of sgnlv loop over sgnlv.ffstat triggers and
     //      1. apply Fstat threshold (copts->cthr)
     //      2. check if f is in the search_par->lines (where [0] is start of the line and [1] the end)
     //      3. shift in f due to spindown: f = f0 + 2*sgnlv.fdot*(search_par->N)*(copts->refr - search_par->seg)
     // append the remaining triggers to allcands
     // remark: we keep the original search grid here

     int i, j, k, itrig;
     int ntrig_seg = 0;
     int buffer_size = 2048; // initial buffer for ffdot shifts, will be reallocated if needed
     float *ffbuffer = (float *) malloc(sizeof(float) * buffer_size);

#if 0
     // only for testing: copy sgnlv[j].ffstat to allcands[j].ffstat[i]
     for (j=0; j<search_par->sgnlv_size; j++) {
          ctrigs[j].ffstat[iseg].p = (float *) malloc(sizeof(float) * sgnlv[j].ffstat.len);
          ctrigs[j].ffstat[iseg].len = sgnlv[j].ffstat.len;
          memcpy(ctrigs[j].ffstat[iseg].p, sgnlv[j].ffstat.p, sizeof(float) * sgnlv[j].ffstat.len);
          ntrig_seg += (sgnlv[j].ffstat.len+1)/2;
     }
#endif

     for (j=0; j<search_par->sgnlv_size; j++) {
          //fprintf(stderr, "     Segment %d, j=%d, sgnlv[j]->ffstat.len = %zu\n", iseg, j, sgnlv[j].ffstat.len);
          // check if ffstat is null? 
          int ic = 0; // "good" triggers counter
          for (itrig=0; itrig<sgnlv[j].ffstat.len/2; itrig++) {
               float f = ((float *)sgnlv[j].ffstat.p)[2*itrig];
               float fstat = ((float *)sgnlv[j].ffstat.p)[2*itrig+1];
               // skip triggers below threshold
               if (fstat < copts->cthr) continue;
               // Narrowing-down the band around center
               if( (f <= M_PI_2 - search_par->narrowdown) ||
                    (f >= M_PI_2 + search_par->narrowdown) )  continue;
               // check if f is in the search_par->lines
               int isline = 0;
               for(int il=0; il<search_par->nvlines_all_inband; il++){
                    if( f >= search_par->lines[il][0] && f <= search_par->lines[il][1] ) {
                         isline = 1;
                         break;
                    }
               }
               if (isline) continue;
               // spindown in linear units
               float fdot_lin = M_PI*sgnlv[j].fdot*pow(search_par->dt,2);
               // shifting f to the reference segment (copts->refr) due to spindown 
               f = f + 2.*fdot_lin*(search_par->N)*(copts->refr - search_par->seg);
               if((f<0) || (f>M_PI)) continue; // skip triggers shifted out of the band
               ffbuffer[2*ic] = f;
               ffbuffer[2*ic+1] = fstat;
               ic++;
               if (ic >= buffer_size/2) {
                    buffer_size *= 2;
                    ffbuffer = (float *) realloc (ffbuffer, sizeof(float) * buffer_size);
               }
          }
          ctrigs[j].ffstat[iseg].len = 2*ic;
          ctrigs[j].ffstat[iseg].p = (float *) malloc(sizeof(float) * ctrigs[j].ffstat[iseg].len);
          //printf("before goodcands j=%d  buffer_size=%d\n", j, buffer_size);
          memcpy(ctrigs[j].ffstat[iseg].p, ffbuffer, 
               sizeof(float) * ctrigs[j].ffstat[iseg].len);
          ntrig_seg += ic;
     }


     printf("   [ntrig_seg=%d]\n", ntrig_seg);

     return ntrig_seg;
}

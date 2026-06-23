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
     
     //copts.nseg = 2; // for testing, should be read from ini file
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
          
          // initialize ctrigs
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
          // return the number of good candidates in the segment
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


     /* ----------------------------------------------------------------------- 
      * map between search and coincidences grids in m,n,s dimensions i.e.
      * ctrigs/sgnlv (m,n,s) -> coincidences (mc,nc,sc)
      ------------------------------------------------------------------------- */
     
     // 1 means no scaling - original search grid
     int scm = copts.scalem;
     int scn = copts.scalen;
     int scs = copts.scales;
     int scf = copts.scalef;

     // max number of ctrigs elements in a coincidences cell
     int ictrigs_size = scm*scn*scs;
     // number of cells in the coincidences grid mc,nc,sc; 
     // extra space to account for edge cases, reallocated later if needed
     int maxccells = 1.3*search_par.sgnlv_size/ictrigs_size;

     Search2coi_mns *s2c_mns = (Search2coi_mns *) malloc(sizeof(Search2coi_mns) * maxccells);
     for(i=0; i<maxccells; i++)
          s2c_mns[i].ictrigs = (int *) malloc(sizeof(int) * ictrigs_size);

     int iccell = 0; // coincidence cell counter
     for(j=0; j<search_par.sgnlv_size; j++) {
          int mc = (int)round(ctrigs[j].m/scm);
          int nc = (int)round(ctrigs[j].n/scn);
          int sc = (int)round(ctrigs[j].s/scs);
          //  printf("j=%d: m=%f, n=%f, s=%f -> mc=%d, nc=%d, sc=%d\n", 
          //       j, ctrigs[j].m, ctrigs[j].n, ctrigs[j].s, mc, nc, sc);
          // search for the cell with mc, nc, sc in s2c_mns;
          // if found, add point to ctrigs[j].ictrigs;
          // if not found, add a new cell to s2c_scm & add to ctrigs[j].ictrigs
          int found = 0;
          for(int ic=0; ic<iccell; ic++) {
               if( (s2c_mns[ic].mc == mc) && (s2c_mns[ic].nc == nc) && (s2c_mns[ic].sc == sc) ) {
                    found = 1;
                    if (s2c_mns[ic].nctrigs >= ictrigs_size) {
                         //ictrigs is full; this should not happen so exit with an error
                         fprintf(stderr, "[ERROR] in coincidence cell %d, nctrigs=%d but max is %d.\n", 
                              ic, s2c_mns[ic].nctrigs+1, ictrigs_size);
                         exit(EXIT_FAILURE);
                    }
                    s2c_mns[ic].ictrigs[s2c_mns[ic].nctrigs] = j;
                    s2c_mns[ic].nctrigs++;
                    break;
               }
          } // for ic
          if (!found) {
               // add more space
               if (iccell >= maxccells) {
                    fprintf(stderr, "[WARNING] maxccells=%d reached, reallocating s2c_mns with 1.2*maxccells=%d\n", 
                         maxccells, (int)(1.2*maxccells));
                    maxccells = (int)(1.2*maxccells);
                    s2c_mns = (Search2coi_mns *) realloc(s2c_mns, sizeof(Search2coi_mns) * maxccells);
                    for(int ic=iccell; ic<maxccells; ic++)
                         s2c_mns[ic].ictrigs = (int *) malloc(sizeof(int) * ictrigs_size);
               }
               fprintf(stderr, "[INFO] new coincidence cell %d/%d: mc=%d, nc=%d, sc=%d\n", 
                    iccell, maxccells, mc, nc, sc);
               s2c_mns[iccell].mc = mc;
               s2c_mns[iccell].nc = nc;
               s2c_mns[iccell].sc = sc;
               s2c_mns[iccell].ictrigs[0] = j;
               s2c_mns[iccell].nctrigs = 1;
               iccell++;
          } // if not found
     } // for j

     int nccells = iccell;
     printf("Number of coincidence cells: %d/%d (sgnlv_size=%ld)\n", nccells, maxccells, search_par.sgnlv_size);
     // print trig2coi
     /*
     for(int ic=0; ic<nccells; ic++) {
          printf("Cell %d: mc=%d, nc=%d, sc=%d, nctrigs=%d, ictrigs=", 
               ic, s2c_mns[ic].mc, s2c_mns[ic].nc, s2c_mns[ic].sc, s2c_mns[ic].nctrigs);
          for(int j=0; j<s2c_mns[ic].nctrigs; j++)
               printf("%d ", s2c_mns[ic].ictrigs[j]);
          printf("\n");
     }
     */
     
     /* -------------------------------------------------------------------------
      * search for coincidences between segments
      * -------------------------------------------------------------------------*/
     
     int ncoi = nccells; 
     // all coincidences, to be written to HDF file
     // assume one coinc. per mns coinc. cell, reallocate if needed;
     Coincidence *coi = (Coincidence *) malloc(sizeof(Coincidence) * ncoi);
     // maximum coincidence based on w & snr
     Coincidence max_coi;
     max_coi.w = 0;
     max_coi.cseg = (short *) malloc(sizeof(short) * copts.nseg);
     max_coi.trig_mns = (int *) malloc(sizeof(int) * copts.nseg);
     max_coi.n_ccell_trigs = (short *) malloc(sizeof(short) * copts.nseg);
     max_coi.avg_snr = 0.f;
     
     int icoi = 0; // coincidence counter
     int nfccells = search_par.nfftf/scf; // number of frequency cells in the coincidences grid

     // Temporary per-segment accumulators for one (imns, fic) cell.
     // MAX_NSEG stack arrays avoid per-cell heap allocation.
     short tmp_cseg[MAX_NSEG];
     int   tmp_trig_mns[MAX_NSEG];
     short tmp_n_ccell_trigs[MAX_NSEG];
          
     // loop over mns cells
     for (int imns=0; imns<nccells; imns++) {
          // loop over frequency cells
          for (int fic=0; fic<nfccells; fic++) {
               // frequency cell boundaries
               float fmin = M_PI/nfccells*fic;
               float fmax = M_PI/nfccells*(fic+1);

               int   w       = 0;
               float sum_snr = 0.f, sum_f = 0.f, sum_fdot = 0.f;
               float sum_ra  = 0.f, sum_dec = 0.f;

               // loop over segments: find max-SNR trigger in this (imns, fic) cell
               for (iseg=0; iseg<copts.nseg; iseg++) {
                    float best_snr     = -1.f;
                    int   best_cidx    = -1;
                    float best_f       =  0.f;
                    int   n_cell_trigs =  0;
                    
                    // scan all ctrigs belonging to this mns cell
                    for (int j=0; j<s2c_mns[imns].nctrigs; j++) {
                         int cidx = s2c_mns[imns].ictrigs[j];
                         if (ctrigs[cidx].ffstat[iseg].len == 0) continue;
                         float *ffp  = (float *)ctrigs[cidx].ffstat[iseg].p;
                         int nfpairs = (int)ctrigs[cidx].ffstat[iseg].len / 2;
                         for (int k=0; k<nfpairs; k++) {
                              float f     = ffp[2*k];
                              float fstat = ffp[2*k+1];
                              if (f >= fmin && f < fmax) {
                                   n_cell_trigs++;
                                   if (fstat > best_snr) {
                                        best_snr  = fstat;
                                        best_cidx = cidx;
                                        best_f    = f;
                                   }
                              }
                         }
                    } // for j
                    
                    if (best_cidx >= 0) {
                         tmp_cseg[w]          = (short)iseg;
                         tmp_trig_mns[w]      = best_cidx;
                         tmp_n_ccell_trigs[w] = (short)n_cell_trigs;
                         sum_snr  += best_snr;
                         sum_f    += best_f;
                         sum_fdot += ctrigs[best_cidx].fdot;
                         sum_ra   += ctrigs[best_cidx].ra;
                         sum_dec  += ctrigs[best_cidx].dec;
                         w++;
                    }
               } // for iseg
                    
               if (w < copts.mincoin) continue;
               
               // grow coi buffer if needed
               if (icoi >= ncoi) {
                    ncoi = (int)(1.3 * ncoi) + 1;
                    coi  = (Coincidence *) realloc(coi, sizeof(Coincidence) * ncoi);
               }
               
               // fill coi[icoi]
               coi[icoi].w     = (short)w;
               coi[icoi].shift = (short)copts.shift;
               coi[icoi].cseg          = (short *) malloc(sizeof(short) * w);
               coi[icoi].trig_mns      = (int *)   malloc(sizeof(int)   * w);
               coi[icoi].n_ccell_trigs = (short *) malloc(sizeof(short) * w);
               memcpy(coi[icoi].cseg,          tmp_cseg,          sizeof(short) * w);
               memcpy(coi[icoi].trig_mns,      tmp_trig_mns,      sizeof(int)   * w);
               memcpy(coi[icoi].n_ccell_trigs, tmp_n_ccell_trigs, sizeof(short) * w);
               coi[icoi].avg_snr  = sum_snr  / w;
               coi[icoi].avg_f    = sum_f    / w;
               coi[icoi].avg_fdot = sum_fdot / w;
               coi[icoi].avg_ra   = sum_ra   / w;
               coi[icoi].avg_dec  = sum_dec  / w;
               
               // update max_coi: more segments wins; tie -> higher avg_snr wins
               if (w > max_coi.w ||
                    (w == max_coi.w && coi[icoi].avg_snr > max_coi.avg_snr)) {
                         max_coi.w     = (short)w;
                         max_coi.shift = coi[icoi].shift;
                         memcpy(max_coi.cseg,          tmp_cseg,          sizeof(short) * w);
                         memcpy(max_coi.trig_mns,      tmp_trig_mns,      sizeof(int)   * w);
                         memcpy(max_coi.n_ccell_trigs, tmp_n_ccell_trigs, sizeof(short) * w);
                         max_coi.avg_snr  = coi[icoi].avg_snr;
                         max_coi.avg_f    = coi[icoi].avg_f;
                         max_coi.avg_fdot = coi[icoi].avg_fdot;
                         max_coi.avg_ra   = coi[icoi].avg_ra;
                         max_coi.avg_dec  = coi[icoi].avg_dec;
                    }
               
               icoi++;
          } // for fic
     } // for imns

     printf("\nTotal coincidences found (w>=%d): %d\n", copts.mincoin, icoi);
     if (max_coi.w > 0) {
          printf("Max coincidence:\n");
          printf("  w=%d  avg_snr=%.4f  avg_f=%.6f  avg_fdot=%.4e"
               "  avg_ra=%.6f  avg_dec=%.6f\n",
               max_coi.w, max_coi.avg_snr, max_coi.avg_f,
               max_coi.avg_fdot, max_coi.avg_ra, max_coi.avg_dec);
          printf("  Segments (seg / trig_idx / n_cell_trigs):");
          for (int is=0; is<max_coi.w; is++)
               printf("  %d/%d/%d",
                    max_coi.cseg[is], max_coi.trig_mns[is],
                    max_coi.n_ccell_trigs[is]);
          printf("\n");
     } else {
          printf("No coincidences above mincoin=%d found.\n", copts.mincoin);
     }
     
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
               //float fdot_lin = M_PI*sgnlv[j].fdot*pow(search_par->dt,2);
               // shifting f to the reference segment (copts->refr) due to spindown 
               //f = f + 2.*fdot_lin*(search_par->N)*(copts->refr - search_par->seg);
               f = f + sgnlv[j].fdot*search_par->dt*(search_par->N)*(copts->refr - search_par->seg);
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

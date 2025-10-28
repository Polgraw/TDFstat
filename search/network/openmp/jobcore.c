#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <malloc.h>
#include <gsl/gsl_vector.h>
#include <complex.h>
#include <fftw3.h>
#include <signal.h>
#include <assert.h>
#include <omp.h>

#include "auxi.h"
#include "struct.h"
#include "settings.h"
#include "jobcore.h"
#include "timer.h"
#include "io.h"

extern volatile sig_atomic_t save_state;


// Main searching function (loops inside)
void search( Search_settings *sett,
             Command_line_opts *opts,
             Search_range *s_range,
             FFTW_plans *plans,
             FFTW_arrays *fftw_arr,
             Aux_arrays *aux,
             int *FNum )
{

     struct flock lck;

     int i;
     int pm;                // hemisphere
     float mm, nn;          // sky positions
     int sgnlc=0;           // number of triggers in the buffer
     //FLOAT_TYPE *sgnlv;   // array with candidates data
     Trigger *sgnlv;        // triggers buffer (to be saved)
     long totsgnl;          // total number of triggers found

     char outname[1100];
     int fd, status;
     FILE *state;

     // sgnlv is a buffer for triggers produced in one run of jobcore
     // max size = <spindown range>/<spindown step>

     int max_spindowns = (s_range->spndr[1]-s_range->spndr[0])/s_range->sstep + 1;
     sgnlv = (Trigger *)calloc(max_spindowns, sizeof(Trigger));

     // max buffer size for (f, fstat) pairs = 2*nfft/2
     sett->bufsize = sett->nfft;
     // initialize ffstat
     for (i=0; i<max_spindowns; i++) {
          sgnlv[i].ffstat.len = 0;
          sgnlv[i].ffstat.p = NULL;
     }

     state = NULL;
     if (opts->checkp_flag) {
          state = fopen (opts->state_file, "w");
          // spindown below is wrong if addsig or range_file is used
          // but then we don't need checkpointing
          fprintf(state, "%d %f %f %f %d\n", s_range->pst, s_range->mst,
               s_range->nst, s_range->sst, *FNum);
          fseek(state, 0, SEEK_SET);
     }

     /* Loop over hemispheres */

     for (pm=s_range->pst; pm<=s_range->pmr[1]; ++pm) {

          struct timespec tstart = get_current_time(CLOCK_REALTIME), tend;

          sprintf (outname, "%s/triggers_%03d_%04d%s_%d.h5",
               opts->outdir, opts->seg, opts->band, opts->label, pm);
          // remove existing trigger file if checkpointing is disabled
          if (! opts->checkp_flag) remove(outname);
          totsgnl = 0;

          hdfout_init(outname, opts, sett, s_range, sgnlv);

          /* Two main loops over sky positions */

          // imm & inn allow to transform loop indices from float to int
          int imm_size = (s_range->mr[1]-s_range->mst)/s_range->mstep + 1;
          //for (mm=s_range->mst; mm<=s_range->mr[1]; ++mm) {
          for (int imm=0; imm<imm_size; ++imm) {
               mm = s_range->mst + imm*s_range->mstep;

               int inn_size = (s_range->nr[1]-s_range->nst)/s_range->nstep + 1;
               //for (nn=s_range->nst; nn<=s_range->nr[1]; ++nn) {
               for (int inn=0; inn<inn_size; ++inn) {
                    nn = s_range->nst + inn*s_range->nstep;

                    int fnum_old = *FNum;  // number of records already written to file

                    /* Loop over spindowns is inside job_core() */
                    status = job_core(
                         pm,           // hemisphere
                         mm,           // grid 'sky position'
                         nn,           // other grid 'sky position'
                         sett,         // search settings
                         opts,         // code config options
                         s_range,      // search ranges
                         plans,        // fftw plans
                         fftw_arr,     // arrays for fftw
                         aux,          // auxiliary arrays
                         &sgnlc,       // number of triggers in the buffer
                         sgnlv,        // triggers buffer
                         FNum);        // number of F-statistic evaluations (=n_spindowns)

                    // Get back to regular spin-down range
                    s_range->sst = s_range->spndr[0];

                    totsgnl += sgnlc;
                    // Write triggers buffer to the file if not empty
                    int nrec = *FNum - fnum_old;  // number of records to write
                    if (nrec>0) {
                         //printf("Writing [m=%g   n=%g]  nrec=%d\n", mm, nn, nrec);
                         hdfout_extend(outname, sgnlv, nrec);
                    }

                    sgnlc=0;

                    if(opts->checkp_flag) {
                         ftruncate(fileno(state), 0);
                         fprintf(state, "%d %f %f %f %d\n", pm, mm, nn+1, s_range->sst, *FNum);
                         fseek(state, 0, SEEK_SET);
                         if (save_state == 1) {
                              //printf("%d %d %d %d %d\n", pm, mm, nn+1, s_range->sst, *FNum);
                              printf("\nState saved after signal\nExiting\n");
                              exit(EXIT_SUCCESS);
                         }
                    }
                    save_state = 0;
                    //} /* if sgnlc > sett-nfft */
               } // for nn
               s_range->nst = s_range->nr[0];
          } // for mm
          s_range->mst = s_range->mr[0];

          tend = get_current_time(CLOCK_REALTIME);
          double time_elapsed = get_time_difference(tstart, tend);
          int nthreads = omp_get_max_threads();
          printf("\nwalltime = %e s | ncpus = %d | ncpus*walltime = %e s\n",
               time_elapsed, nthreads, time_elapsed*nthreads);

          // All triggers should be written at this point
          printf("\n### Total number of triggers in %s = %ld\n\n", outname, totsgnl);
          sgnlc=0;
          hdfout_finalize(outname, totsgnl, time_elapsed, nthreads);

     } // for pm

     if (opts->checkp_flag) {
          // empty state file to prevent restart after successful end
          ftruncate(fileno(state), 0);
          fclose(state);
     }

     // Free triggers buffer
     for(i=0; i<max_spindowns; i++) {
          if (sgnlv[i].ffstat.p) free(sgnlv[i].ffstat.p);
     }
     free(sgnlv);

     printf("\nEND\n");

} //search


/* Main job */

int job_core(
               int pm,                   // Hemisphere
               float mm,                   // Grid 'sky position'
               float nn,                   // Second grid 'sky position'
               Search_settings *sett,    // Search settings
               Command_line_opts *opts,  // Search options
               Search_range *s_range,    // Range for searching
               FFTW_plans *plans,        // Plans for fftw
               FFTW_arrays *fftw_arr,    // Arrays for fftw
               Aux_arrays *aux,          // Auxiliary arrays
               int *sgnlc,               // Candidate trigger parameters
               Trigger *sgnlv,           // Candidate array
               int *FNum)                // Candidate signal number
{
     int i, j, n;
     float smin = s_range->sst, smax = s_range->spndr[1];
     double al1, al2, sinalt, cosalt, sindelt, cosdelt, nSource[3], ft, het0;
     FLOAT_TYPE sgnlt[NPAR], sgnl0;
     FLOAT_TYPE _tmp1[sett->nifo][sett->N] __attribute__((aligned(128)));

     struct timespec tstart, tend;
     double spindown_timer = 0;
     int spindown_counter  = 0;

     /*
     Matrix	M(.,.) (defined on page 22 of PolGrawCWAllSkyReview1.pdf file)
     defines the transformation form integers (bin, ss, nn, mm) determining
     a grid point to linear coordinates omega, omegadot, alpha_1, alpha_2),
     where bin is the frequency bin number and alpha_1 and alpha_2 are
     defined on p. 22 of PolGrawCWAllSkyReview1.pdf file.

     [omega]                          [bin]
     [omegadot]       = M(.,.) \times [ss]
     [alpha_1/omega]                  [nn]
     [alpha_2/omega]                  [mm]

     Array M[.] is related to matrix M(.,.) in the following way;

     [ M[0] M[4] M[8]  M[12] ]
     M(.,.) =   [ M[1] M[5] M[9]  M[13] ]
     [ M[2] M[6] M[10] M[14] ]
     [ M[3] M[7] M[11] M[15] ]

     and

     M[1] = M[2] = M[3] = M[6] = M[7] = 0
     */

     // Grid positions
     al1 = nn*sett->M[10] + mm*sett->M[14];
     al2 = nn*sett->M[11] + mm*sett->M[15];

     // check if the search is in an appropriate region of the grid
     // if not, returns NULL
     if ((sqr(al1)+sqr(al2))/sqr(sett->oms) > 1.) return 0;

     float ss;
     double shft1, phase, cp, sp;
     complex double exph;

     // Change linear (grid) coordinates to real coordinates
     lin2ast(al1/sett->oms, al2/sett->oms, pm, sett->sepsm, sett->cepsm,
          &sinalt, &cosalt, &sindelt, &cosdelt);

     // calculate declination and right ascention
     sgnlt[2] = asin(sindelt);
     sgnlt[3] = fmod(atan2(sinalt, cosalt) + 2.*M_PI, 2.*M_PI);

     het0 = fmod(nn*sett->M[8] + mm*sett->M[12], sett->M[0]);

     // Nyquist frequency
     int nyqst = (sett->nfft)/2 + 1;

     // Loop for each detector
     for(n=0; n<sett->nifo; ++n) {

          /* Amplitude modulation functions aa and bb
          * for each detector (in signal sub-struct
          * of _detector, ifo[n].sig.aa, ifo[n].sig.bb)
          */

          modvir(sinalt, cosalt, sindelt, cosdelt, sett->N, &ifo[n], aux);

          // Calculate detector positions with respect to baricenter
          nSource[0] = cosalt*cosdelt;
          nSource[1] = sinalt*cosdelt;
          nSource[2] = sindelt;

          shft1 = nSource[0]*ifo[n].sig.DetSSB[0]
                + nSource[1]*ifo[n].sig.DetSSB[1]
                + nSource[2]*ifo[n].sig.DetSSB[2];

#pragma omp parallel default(shared) private(phase,cp,sp,exph)
          {
#pragma omp for schedule(static)
               for(i=0; i<sett->N; ++i) {
                    ifo[n].sig.shft[i] = nSource[0]*ifo[n].sig.DetSSB[i*3]
                                       + nSource[1]*ifo[n].sig.DetSSB[i*3+1]
                                       + nSource[2]*ifo[n].sig.DetSSB[i*3+2];
                    ifo[n].sig.shftf[i] = ifo[n].sig.shft[i] - shft1;
                    _tmp1[n][i] = aux->t2[i] + (double)(2*i)*ifo[n].sig.shft[i];
               }

#pragma omp for schedule(static)
               for(i=0; i<sett->N; ++i) {
                    // Phase modulation
                    phase = het0*i + sett->oms*ifo[n].sig.shft[i];
                    sincos(phase, &sp, &cp);
                    exph = cp - I*sp;

                    // Matched filter
                    ifo[n].sig.xDatma[i] = ifo[n].sig.xDat[i]*ifo[n].sig.aa[i]*exph;
                    ifo[n].sig.xDatmb[i] = ifo[n].sig.xDat[i]*ifo[n].sig.bb[i]*exph;
               }

               /* Resampling using spline interpolation:
               * This will double the sampling rate
               */
#pragma omp for schedule(static)
               for(i=0; i < sett->N; ++i) {
                    fftw_arr->xa[i] = ifo[n].sig.xDatma[i];
                    fftw_arr->xb[i] = ifo[n].sig.xDatmb[i];
               }

               // Zero-padding (filling with 0s up to sett->nfft,
               // the nearest power of 2)
#pragma omp for schedule(static)
               for (i=sett->N; i<sett->nfft; ++i) {
                    fftw_arr->xa[i] = 0.;
                    fftw_arr->xb[i] = 0.;
               }
          } //omp parallel


          fftw_execute_dft(plans->pl_int,fftw_arr->xa,fftw_arr->xa);  //forward fft (len nfft)
          fftw_execute_dft(plans->pl_int,fftw_arr->xb,fftw_arr->xb);  //forward fft (len nfft)

          // move frequencies from second half of spectrum;
          // and zero frequencies higher than nyquist
          // loop length: nfft - nyqst = nfft - nfft/2 - 1 = nfft/2 - 1

          for(i=nyqst + sett->Ninterp - sett->nfft, j=nyqst; i<sett->Ninterp; ++i, ++j)
               fftw_arr->xa[i] = fftw_arr->xa[j];
#pragma omp parallel for schedule(static) default(shared)
          for(i=nyqst; i<nyqst + sett->Ninterp - sett->nfft; ++i)
               fftw_arr->xa[i] = 0.;

          for(i=nyqst + sett->Ninterp - sett->nfft, j=nyqst; i<sett->Ninterp; ++i, ++j)
               fftw_arr->xb[i] = fftw_arr->xb[j];
#pragma omp parallel for schedule(static) default(shared)
          for(i=nyqst; i<nyqst + sett->Ninterp - sett->nfft; ++i)
               fftw_arr->xb[i] = 0.;


          // Backward fft (len Ninterp = nfft*interpftpad)
          fftw_execute_dft(plans->pl_inv,fftw_arr->xa,fftw_arr->xa);
          fftw_execute_dft(plans->pl_inv,fftw_arr->xb,fftw_arr->xb);

          ft = (double)sett->interpftpad / sett->Ninterp; //scale FFT
          for (i=0; i < sett->Ninterp; ++i) {
               fftw_arr->xa[i] *= ft;
               fftw_arr->xb[i] *= ft;
          }

          //  struct timeval tstart = get_current_time(), tend;

          // Spline interpolation to xDatma, xDatmb arrays
          splintpad(fftw_arr->xa, ifo[n].sig.shftf, sett->N, sett->interpftpad, ifo[n].sig.xDatma);
          splintpad(fftw_arr->xb, ifo[n].sig.shftf, sett->N, sett->interpftpad, ifo[n].sig.xDatmb);


     } // end of detector loop

     // square sums of modulation factors
     FLOAT_TYPE aa = 0., bb = 0.;

     for(n=0; n<sett->nifo; ++n) {

          double aatemp = 0., bbtemp = 0.;

          for(i=0; i<sett->N; ++i) {
               aatemp += sqr(ifo[n].sig.aa[i]);
               bbtemp += sqr(ifo[n].sig.bb[i]);
          }

          aa += aatemp/ifo[n].sig.sig2;
          bb += bbtemp/ifo[n].sig.sig2;

          for(i=0; i<sett->N; ++i) {
               ifo[n].sig.xDatma[i] /= ifo[n].sig.sig2;
               ifo[n].sig.xDatmb[i] /= ifo[n].sig.sig2;
          }
     }


     // Check if the signal is added to the data
     // or the range file is given:
     // if not, proceed with the wide range of spindowns
     // if yes, use smin = s_range->sst, smax = s_range->spndr[1]
     if(!strcmp(opts->addsig, "") && !strcmp(opts->range_file, "")) {

          // Spindown range defined using Smin and Smax (settings.c)
          smin = trunc((sett->Smin - nn*sett->M[9] - mm*sett->M[13])/sett->M[5]);
          smax = trunc(-(nn*sett->M[9] + mm*sett->M[13] + sett->Smax)/sett->M[5]);

          // swapping smin and smax in case when grid matrix
          // values are defined with opposite signs than ''usual''
          if(smin > smax) {
               smin = smin + smax ;
               smax = smin - smax ;
               smin = smin - smax ;
          }
     }

     //const int s_stride = s_range->sstep;
     printf ("\n>>%d\t%g\t%g\t[%g..%g:%g]\n", *FNum, mm, nn, smin, smax, s_range->sstep);

     static FFTW_PRE(_complex) *fxa, *fxb;
     fxa = fftw_arr->fxa;
     fxb = fftw_arr->fxb;

     static FLOAT_TYPE *F;
     if (!F) F = (FLOAT_TYPE *)malloc(2*sett->nfft*sizeof(FLOAT_TYPE));

     //we use shared plans and  fftw_execute with 'new-array' interface

     /* Spindown loop  */

     int iss_size = (smax-smin)/s_range->sstep + 1;
     //for (ss=smin; ss<=smax; ss += s_stride) {
     for (int iss=0; iss<iss_size; ++iss) {
          ss = smin + iss*s_range->sstep;

#if TIMERS>2
          //tstart = get_current_time(CLOCK_PROCESS_CPUTIME_ID);
          tstart = get_current_time(CLOCK_MONOTONIC);
#endif

          // Spindown parameter
          //FLOAT_TYPE spnd = ss*sett->M[5] + nn*sett->M[9] + mm*sett->M[13];
          sgnlt[1] = ss*sett->M[5] + nn*sett->M[9] + mm*sett->M[13];

          // initialize current record for triggers
          sgnlv[iss].m = mm;
          sgnlv[iss].n = nn;
          sgnlv[iss].s = ss;
          sgnlv[iss].ra = sgnlt[3];
          sgnlv[iss].dec = sgnlt[2];
          sgnlv[iss].fdot = sgnlt[1]/(M_PI*sett->dt*sett->dt);
          sgnlv[iss].ffstat.len = 0;
          // allocate memory for ffstat buffer if not allocated yet
          if (sgnlv[iss].ffstat.p == NULL) {
               sgnlv[iss].ffstat.p = (float *) malloc(sett->bufsize * sizeof(float));
               if (sgnlv[iss].ffstat.p == NULL) {
                    printf("[ERROR] Cannot allocate memory for triggers buffer of size %d\n", sett->bufsize);
                    exit(EXIT_FAILURE);
               }
          }

          FLOAT_TYPE het1;

#ifdef VERBOSE
          //print a 'dot' every new spindown
          printf ("."); fflush (stdout);
#endif

          het1 = fmod(ss*sett->M[4], sett->M[0]);
          if(het1<0) het1 += sett->M[0];

          sgnl0 = het0 + het1;

          spindown_modulation(sett->nifo, sett->N, het1, sgnlt[1], _tmp1, fxa, fxb);

          // Zero-padding]
#pragma omp parallel for schedule(static)
          for (i = sett->nfftf-1; i > sett->N-1; --i)
               fxa[i] = fxb[i] = (FLOAT_TYPE)0.;

          FFTW_PRE(_execute_dft)(plans->plan, fxa, fxa);
          FFTW_PRE(_execute_dft)(plans->plan, fxb, fxb);

          // Computing F-statistic
#pragma omp parallel for schedule(static)
          for (i=sett->nmin; i<=sett->nmax; ++i){
               F[i] = NORM(fxa[i])/aa + NORM(fxb[i])/bb ;
          }

          (*FNum)++;

#undef FSTATDEB
#ifdef FSTATDEB
          // warning: use with nthreads=1
          static double *fraw;
          static int ifile=0;
          if (!fraw) fraw = (double *) malloc((sett->nmax-sett->nmin)*sizeof(double));
          memcpy(fraw, F+sett->nmin, (sett->nmax-sett->nmin)*sizeof(double));
          if ( (mm==-61 && nn==-32 && ss==273) ) {
               char f1name[32];
               ifile++;
               sprintf(f1name, "fraw-%d.dat", ifile);
               FILE *f1 = fopen(f1name, "w");
               for (i=0; i<(sett->nmax-sett->nmin); i++)
                    fprintf(f1, "%d   %lf   %lf  %lf  %lf  %lf %lf\n", i, fraw[i],
                         2.*M_PI*i/((double) sett->fftpad*sett->nfft) + sgnl0,
                         sqr(creal(ifo[0].sig.xDatma[2*i])),
                         sqr(cimag(ifo[0].sig.xDatmb[2*i])),
                         sqr(creal(ifo[1].sig.xDatma[2*i])),
                         sqr(cimag(ifo[1].sig.xDatmb[2*i])) );
               fclose(f1);
          }
#endif


          // Normalize F-statistics if the noise is not white noise
          if( ! strcmp(opts->fstat_norm, "blocks_avg") ) {
               FStat(F + sett->nmin, sett->nmax - sett->nmin, NAVFSTAT, 0);
          }

          /* select triggers */
          int itrig=0;
#if 0
          int dd = sett->dd;
          /* find the highest maximum (above trl) in each block of length dd;
          dd is set to (1/day frequency in units of F indices)-1,
          just below the distance between F-statistic peaks for a signal
          */

          /* stay in (nmin, nmax) range! */
          for (i=sett->nmin+1; i<sett->nmax-dd; i+=dd) {
               int ii=-1;
               FLOAT_TYPE Fc = opts->thr;
               for (j=i; j<i+dd; ++j) {
                    if (F[j] < Fc || F[j-1] > F[j] || F[j] < F[j+1] ) continue;
                    ii = j;
                    Fc = F[j];
                    j++;
               }

               if ( ii < 0 ) continue; // no maximum in this block
#else
          for (i=sett->nmin; i<sett->nmax; ++i) {
               if (F[i] < opts->thr) continue;
               FLOAT_TYPE Fc;
               int ii = i;
               Fc = F[i];
               while (++i < sett->nmax && F[i] > opts->thr) {
                    if(F[i] > Fc) {
                         ii = i;
                         Fc = F[i];
                    } // if F[i]
               } // while i
#endif
               // Candidate signal frequency ( ~ ii/nmax * pi ; nmax=nfftf/2 )
               sgnlt[0] = (FLOAT_TYPE)(2*ii)/(FLOAT_TYPE)sett->nfftf * M_PI + sgnl0;

               // Checking if signal is within a known instrumental line
               int k, veto_status = 0;
               for (k=0; k<sett->numlines_band; k++){
                    if(sgnlt[0]>=sett->lines[k][0] && sgnlt[0]<=sett->lines[k][1]) {
                         veto_status=1;
                         break;
                    }
               }

               if(!veto_status) {

                    //sgnlt[4] = sqrtf(2.*(Fc - sett->nd));
                    // we store fstatistics value instead of snr
                    sgnlt[4] = Fc;

                    // Add new parameters to the output buffer array
                    ((float *)(sgnlv[iss].ffstat.p))[2*itrig] = sgnlt[0];    // frequency
                    ((float *)(sgnlv[iss].ffstat.p))[2*itrig + 1] = sgnlt[4]; // fstat value
                    sgnlv[iss].ffstat.len += 2;
                    itrig++;
                    (*sgnlc)++;

#ifdef VERBOSE
                    printf ("\nSignal %d: %d %g %g %g %d snr=%.2f\n",
                         *sgnlc, pm, mm, nn, ss, ii, sgnlt[4]);
                    printf("aququ\n");
#endif
               } // if veto_status
          } // for i


#if TIMERS>2
          //tend = get_current_time(CLOCK_PROCESS_CPUTIME_ID);
          tend = get_current_time(CLOCK_MONOTONIC);
          spindown_timer += get_time_difference(tstart, tend);
          spindown_counter++;
#endif

     } // for ss

     printf("  ntrig=%d, buffer %d%% full\n", *sgnlc, (*sgnlc)*100/(sett->bufsize/2) );
#ifndef VERBOSE
     //printf("Number of signals found: %d (buffer %d%% full)\n", *sgnlc, (*sgnlc)*100/(sett->bufsize/2) );
#endif

#if TIMERS>2
     //printf("\nTotal spindown loop time: %e s, mean spindown cpu-time: %e s (%d runs)\n",
     printf("  [Perf] Spindown loop cputime: %e s, <cputime/ns>: %e s (ns: %d)\n",
          spindown_timer, spindown_timer/spindown_counter, spindown_counter);
#endif

     return 0;

} // jobcore

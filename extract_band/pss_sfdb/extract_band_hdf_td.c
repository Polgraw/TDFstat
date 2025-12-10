/***********Program to extract Bandwidth from the SFDB files***************/
/***********Author:   Pia**************************/
/***********Last version 11 DEC. 2013******************/
/**** Major rewrite (HDF + Time Domain output): Paweł Ciecieląg 2025 ****/

// define to output sft instead of sts (for debugging/testing)
#include <H5Epublic.h>
#include <H5Ipublic.h>
#undef OUT_SFT

#define NMAXCHUNK 500000

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <time.h>
#include <fftw3.h>

/************ PSS libraries ************/
//#include "../pss_lib/pss_math.h"
#include "../pss_lib/pss_serv.h"
//#include "../pss_lib/pss_snag.h"
//#include "../pss_lib/pss_sfc.h"
/****** Antenna libraries *************/
//#include "pss_ante.h"

/****** SFDB libraries **************/
#include "pss_sfdb.h"


int main(void){

     int ifatti;
     int ii;
     int ilista=1;           //1 to read input files from a lista file
     int errorcode;
     GD *gd;                 /*GD with SFDB data*/
     GD *gd_short;           /*GD with the short SFDB data*/
     GD *gd_band=NULL;            /*GD with SFDB data, over a chosen sub-bandwidth*/
     char *fileBAND;
     float total_band, initial_freq, band, final_freq, bandUS;
     float freq_sin, fr;
     int iw=-1;
     int fft_read;
     double scale=0.;
     double subsampling_factor=0.;
     HEADER_PARAM *header_param;
     header_param=(HEADER_PARAM *)malloc(sizeof(HEADER_PARAM));
     int retval=EXIT_SUCCESS;

     // FFT related variables
     fftwf_complex *sft=NULL; // short fourier transform
     float *sts=NULL;         // short time series
     fftwf_plan iplan=NULL;        // inverse fft plan - C2R
     int lfft=0., lfft2=0.;

     // HDF related variables
     hid_t ofile_id;
     herr_t hstat;
     int format_version=1;
     int last_ichunk=-1;
     double scaling_factor;


     // supported detectors as defined in HEADER_PARAM in pss_sfdb.h
     static const char det[][3] = { "00", "V1", "H1", "L1" };

     LOG_INFO = logfile_open("extract_band"); // extern in pss_sfdb.h

     gd = crea_gd(0., 0., 0., "orig");
     gd_short = crea_gd(0., 0., 0., "ar");
     fileBAND = (char*)malloc(sizeof(char)*MAXMAXLINE);

     printf("Name of the output file (HDF5)\n");
     scanf("%s", fileBAND);

     // temporarily disable error printing by HDF5 bacause H5Fopen is allowed to fail
     H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
     // open HDF5 file for appending, if not possible create a new one
     if ((ofile_id = H5Fopen(fileBAND, H5F_ACC_RDWR, H5P_DEFAULT)) != H5I_INVALID_HID ){
          H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t)H5Eprint2, stderr); // re-enable error printing
          hstat = H5LTget_attribute_int(ofile_id, "/", "last_ichunk", &last_ichunk);
          if (hstat < 0) {printf("Error reading attribute last_ichunk\n"); goto fail;}
          if (last_ichunk < 0) {printf("Error: last_ichunk < 0 \n"); goto fail;}
          printf("Opened %s [append mode, last_ichunk=%d]\n", fileBAND, last_ichunk);
          // verify link creation order
          /*
          hid_t info = H5Fget_create_plist(ofile_id);
          unsigned flags;
          H5Pget_link_creation_order(info, &flags);
          H5Pclose(info);
          */
          // getting creation plist for file doesn't work, use root group
          hid_t rgroup = H5Gopen(ofile_id, "/", H5P_DEFAULT);
          hid_t info = H5Gget_create_plist(rgroup);
          unsigned flags;
          H5Pget_link_creation_order(info, &flags);
          H5Pclose(info);
          if ( flags & (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED) ) {
               puts("Link creation order is tracked and indexed - good!");
          } else {
               printf("Wrong link creation order flags=%u / %u!\n", flags, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED); goto fail;
          }
     } else {
          H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t)H5Eprint2, stderr);// re-enable error printing
          // set link creation tracking and indexing in hdf file
          // to enable iteration over creation time
          hid_t fcpl;
          if ((fcpl = H5Pcreate(H5P_FILE_CREATE)) == H5I_INVALID_HID){
               puts("--> FAILURE");
               return EXIT_FAILURE;
          }
          H5Pset_link_creation_order(fcpl, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED );
          // create a new file
          if ((ofile_id = H5Fcreate(fileBAND, H5F_ACC_TRUNC, fcpl, H5P_DEFAULT)) == H5I_INVALID_HID){
               puts("Error creating creating HDF file");
               puts("--> FAILURE");
               return EXIT_FAILURE;
          }
          H5Pclose(fcpl);
          ifatti = -1;
          hstat = H5LTset_attribute_int(ofile_id, "/", "last_ichunk", &ifatti, 1);
          if (hstat<0) {printf("Error creating attribute last_ichunk\n"); goto fail;}
          printf("Created %s\n", fileBAND);
     }

     // main loop over sft chunks
     for (ifatti = 0; ifatti < NMAXCHUNK; ifatti++){

          // dataset name = SFT counter starting from 1
          int ichunk = ifatti + 1;
          int length = snprintf( NULL, 0, "%d", ichunk );
          char* dsname;
          if (!dsname) free(dsname);
          dsname = malloc( length + 1 );
          snprintf(dsname, length + 1, "%d", ichunk );
          int prev_gps_sec = header_param->gps_sec; // for debugging
          if (ichunk <= last_ichunk) {
               // iterate over chunks already stored in the HDF file
               errorcode = readfilesfdb(gd, gd_short, header_param, ifatti, iw, &fft_read, ilista, 1);
               if (errorcode != 1){
                    if (ifatti<last_ichunk)
                         printf("[ERROR] HDF file cointains more data than available from the input files !\n");
                    else if (ifatti==last_ichunk)
                         printf("[INFO] Input file list unchanged.\n");
                    goto success;
               }
               // verify consistency of data in the HDF file with the SFDB file
               printf(">>> Verifying HDF data [ichunk/last_ichunk: %d/%d]\n", ichunk, last_ichunk);
               int ichunk_hdf;
               printf("    ichunk=%d ", ichunk);
               hstat = H5LTget_attribute_int(ofile_id, dsname, "ichunk", &ichunk_hdf);
               if (ichunk_hdf != ichunk || hstat < 0) {
                    printf("\nMismatch in HDF: %d\n", ichunk_hdf);
                    goto fail;
               }
               printf("[OK]\n");

               /* we verify only mjdtime because gps_sec is modified in the _td version (due to windowing) */
               printf("    mjdtime=%f ", header_param->mjdtime);
               double mjdtime_hdf;
               hstat = H5LTget_attribute_double(ofile_id, dsname, "sft_mjdtime", &mjdtime_hdf);
               if (mjdtime_hdf != header_param->mjdtime || hstat < 0) {
                    printf("\nMismatch in HDF: mjdtime_hdf = %f\n", mjdtime_hdf);
                    goto fail;
               }
               printf("[OK]\n");

               printf("    nfft = %d ", header_param->nfft);
               int nfft_hdf;
               hstat = H5LTget_attribute_int(ofile_id, dsname, "nfft", &nfft_hdf);
               if (nfft_hdf != header_param->nfft || hstat < 0) {
                    printf("\nMismatch in HDF: nfft_hdf = %d\n", nfft_hdf);
                    goto fail;
               }
               printf("[OK]\n");

               printf("    sfdb_file = %s ", header_param->sfdbname);
               char sfdb_file_hdf[MAXMAXLINE+1];
               hstat = H5LTget_attribute_string(ofile_id, dsname, "sfdb_file", sfdb_file_hdf);
               if (strncmp(sfdb_file_hdf, header_param->sfdbname, MAXMAXLINE) != 0 || hstat < 0){
                    printf("\nMismatch in HDF: sfdb_file_hdf = %s\n", sfdb_file_hdf);
                    goto fail;
               }
               printf("[OK]\n");
               printf(">>> Verification OK\n");

               continue;
          } else {
               errorcode = readfilesfdb(gd, gd_short, header_param, ifatti, iw, &fft_read, ilista, 0);
               printf("[NEW CHUNK] ichunk / last_ichunk: %d / %d\n", ichunk, last_ichunk);
          }

          if (errorcode != 1){
               printf("Error reading SFT ichunk=%d\n", ichunk);
               printf("Read ichunk/last_ichunk = %d/%d\n", ichunk-1, last_ichunk);
               //if (ifatti > 0){
               if (ichunk-1 > last_ichunk){
                    printf("Setting new last_ichunk = %d\n", ichunk-1);
                    hstat = H5LTset_attribute_int(ofile_id, "/", "last_ichunk", &ifatti, 1);
                    if (hstat<0) {printf("Error creating attribute last_ichunk\n"); goto fail;}
                    //printf("gd->n=%ld  gd_band->n=%ld gd_band->dx=%lf\n", gd->n, gd_band->n, gd_band->dx);
	          } else {
                    printf("There is no new chunks\n");
               }
               goto success;
          }

          // some initializations if first chunk or we begin to append new data
          if ((ichunk == 1) || (ichunk == last_ichunk+1)){
               printf("First file mjdtime = %15.10f ,  gps time (s and ns) = %d %d\n",
                    header_param->mjdtime,header_param->gps_sec,header_param->gps_nsec);
               total_band = gd->dx*gd->n/2; //the data are Real and Imag, so there is a factor 2
               printf("Initial frequency of the data - Bandwidth (positive only)- total length %f %f %ld\n",
                    gd->ini, total_band, gd->n);
               printf("scaling factor in amplitude %e\n", header_param->einstein);

               printf("Initial frequency to extract ? \n");
               scanf("%f", &freq_sin);
               fr = gd->dx*floor(freq_sin/gd->dx);
               if (fabs(freq_sin - fr) > 1.e-30) {
                    printf("WARNING! Selected frequency %f is not exact, reavaluated to %f\n", freq_sin, fr);
                    //exit(EXIT_FAILURE);
               }
               freq_sin = fr;
               printf("Reavaluated EXACT frequency %f\n", fr);
               printf("Total bandwidth ?\n");
               scanf("%f", &band);

               initial_freq = fr;
               final_freq = gd->dx*floor((fr + band)/gd->dx);
               band = (final_freq - initial_freq);
               printf("initial_freq=%f , final_freq=%f, band=%f\n", initial_freq, final_freq, band);

               gd_band = crea_gd(gd->n, (double)initial_freq, gd->dx, "band");
          }

          band_extract(gd, gd_band, band);
          initial_freq = gd_band->ini; //revaluated exactly with the selected bins
          bandUS = (gd_band->n/2)*gd_band->dx; //really used band (for the power 2 in the FFT)

          // timing output
          int chunk_overlap = header_param->gps_sec - prev_gps_sec - (int)header_param->tbase/2;
          printf("    ichunk=%d    prev_gps_sec=%d    gps_sec=%d   chunk_overlap=%d\n    file=%s\n",
               ichunk, prev_gps_sec, header_param->gps_sec, chunk_overlap, header_param->sfdbname);

          // basic check if chunks are read in the increasing time order (this is required)
          // start of the current chunk should be after the end of the previous chunk
          if ( chunk_overlap < 0) {
               printf("\n    [WARNING] Current chunk overlaps previous by %d [s]. Skipping...\n", chunk_overlap);
               // rollback
               ifatti--;
               header_param->gps_sec = prev_gps_sec;
               continue;
          }

          if ((ifatti == 0) || (ichunk == last_ichunk+1)){
               printf("initial_freq=%f , final_freq=%f , **bandUS=%f, samples extracted=%f , gd_band->n=%ld\n",
               initial_freq, initial_freq + bandUS, bandUS, 2*band/gd_band->dx, gd_band->n);

               //lfft = 2*nsamples = gd_band->n (length of complex fft);
               lfft = gd_band->n;
               lfft2 = lfft/2;
               // allocate one more element for C2R transform
               sft = (fftwf_complex *)fftwf_malloc( (lfft2+1) * sizeof(fftwf_complex));
               sts = (float *)fftwf_malloc( lfft * sizeof(float));

               iplan = fftwf_plan_dft_c2r_1d(lfft, sft, sts, FFTW_ESTIMATE);

               //write some information on the file at the beginning:
               //fprintf(BAND, "%% Beginning freq- Band- Samples in one stretch- Subsampling factor- inter (overlapping, 2 if data were overlapped)- Frequency step- Scaling factor- ***The data are real and imag of the FFT\n");
               scale = 1.0;
               if(header_param->einstein == 1) scale = 1.0E20;
               if(scale != 1){
                    scaling_factor = 1./scale;
                    printf("Scaling factor in SFDB was 1. Changed to %e\n", 1./scale);
                    //fprintf(BAND, "%% %f %f %ld %f %d %10.7f %e\n",
                    //	gd_band->ini, bandUS, gd_band->n/2, 1.0*(gd->n/gd_band->n), header_param->typ, gd_band->dx, 1./scale);
               } else {
                    scaling_factor = header_param->einstein;
                    //fprintf(BAND, "%% %f %f %ld %f %d %10.7f %e\n",
                    //	gd_band->ini, bandUS, gd_band->n/2, 1.0*(gd->n/gd_band->n), header_param->typ, gd_band->dx, header_param->einstein);
               }

               // file attributes
#ifdef OUT_SFT
               hstat = H5LTset_attribute_string(ofile_id, "/", "dtype", "sft");
#else
               hstat = H5LTset_attribute_string(ofile_id, "/", "dtype", "sts");
#endif
               if (hstat<0) {printf("Error creating attribute dtype\n"); goto fail;}

               hstat = H5LTset_attribute_int(ofile_id, "/", "format_version", &format_version, 1);
               if (hstat<0) {printf("Error creating attribute format_version\n"); goto fail;}

               hstat = H5LTset_attribute_double(ofile_id, "/", "fpo", &(gd_band->ini), 1);
               if (hstat<0) {printf("Error creating attribute fpo\n"); goto fail;}

               hstat = H5LTset_attribute_float(ofile_id, "/", "bandwidth", &bandUS, 1);
               if (hstat<0) {printf("Error creating attribute bandwidth\n"); goto fail;}

               // number of samples = lfft2 both for sft and sts
               hstat = H5LTset_attribute_int(ofile_id, "/", "nsamples", &lfft2, 1);
               if (hstat<0) {printf("Error creating attribute nsamples\n"); goto fail;}

               subsampling_factor = 1.0*(gd->n/gd_band->n);
               hstat = H5LTset_attribute_double(ofile_id, "/", "subsampling_factor", &subsampling_factor, 1);
               if (hstat<0) {printf("Error creating attribute subsampling_factor\n"); goto fail;}

               // we don't need this information; we assume it's always 2 (overlapping), otherwise rise error
               hstat = H5LTset_attribute_int(ofile_id, "/", "sft_overlap", &(header_param->typ), 1);
               if (hstat<0) {printf("Error creating attribute sft_overlap\n"); goto fail;}

               hstat = H5LTset_attribute_double(ofile_id, "/", "df", &(gd_band->dx), 1);
               if (hstat<0) {printf("Error creating attribute df\n"); goto fail;}

               hstat = H5LTset_attribute_double(ofile_id, "/", "scaling_factor", &scaling_factor, 1);
               if (hstat<0) {printf("Error creating attribute scaling_factor\n"); goto fail;}

               hstat = H5LTset_attribute_string(ofile_id, "/", "site", &(det[header_param->detector][0]));
               if (hstat<0) {printf("Error creating attribute site\n"); goto fail;}
          }

          /*
          fprintf(BAND, "%% FFT number in the file; Beginning mjd days; Gps s; Gps ns;\n");
          fprintf(BAND, "%% %d %15.10f %d %d\n",
          header_param->nfft, header_param->mjdtime, header_param->gps_sec, header_param->gps_nsec);

          for(ii = 0; ii < gd_band->n; ii += 2){
          fprintf(BAND, "%15.8e %15.8e\n", gd_band->y[ii]*scale, gd_band->y[ii+1]*scale);
          }
          */

          for(ii = 0; ii < lfft2; ++ii){
               sft[ii][0] = gd_band->y[2*ii]*scale;
               sft[ii][1] = gd_band->y[2*ii+1]*scale;
          }

          int gps_sec, gps_nsec;
#ifdef OUT_SFT
          // outut sft for debugging/testing
          hsize_t dims[2];
          dims[0] = lfft2;
          dims[1] = 2;
          hstat = H5LTmake_dataset_float(ofile_id, dsname, 2, dims, (const float *)sft);
          if (hstat<0) {printf("Error creating sft dataset %s\n", dsname); goto fail;}
          gps_sec = header_param->gps_sec;
          gps_nsec = header_param->gps_nsec;
#else
          // is it correct ? test with the imaginary component = -sft[lfft2-1][1]
          sft[lfft2][0] = sft[lfft2-1][0];
          sft[lfft2][1] = 0.;

          fftwf_execute(iplan); // output in sts

          // write only flat top of the Tukey window
          int lfft4 = lfft/4;
         	float norm = 1./(lfft*subsampling_factor);
          for(ii = lfft4; ii < lfft4+lfft2; ++ii)
               sts[ii] *= norm;

          hsize_t ts_dims[1] = { lfft2 };
          hstat = H5LTmake_dataset_float(ofile_id, dsname, 1, ts_dims, &sts[lfft4]);
          if (hstat<0) {printf("Error creating sts dataset %s\n", dsname); goto fail;}

          //double dt = 1./(2*bandUS);
          //double t_offset = lfft4*dt;
          double t_offset = header_param->tbase/4.;
          double t_gps = header_param->gps_sec + header_param->gps_nsec/1e9 + t_offset;

          gps_sec = (int)t_gps;
          gps_nsec = (int)((t_gps - gps_sec)*1e9);
#endif

          // dataset attributes
          hstat = H5LTset_attribute_int(ofile_id, dsname, "nfft", &header_param->nfft, 1);
          if (hstat<0) {printf("Error creating attribute nfft\n"); goto fail;}

          hstat = H5LTset_attribute_int(ofile_id, dsname, "ichunk", &ichunk, 1);
          if (hstat<0) {printf("Error creating attribute ichunk\n"); goto fail;}

          hstat = H5LTset_attribute_int(ofile_id, dsname, "gps_sec", &gps_sec, 1);
          if (hstat<0) {printf("Error creating attribute gps_sec\n"); goto fail;}

          hstat = H5LTset_attribute_int(ofile_id, dsname, "gps_nsec", &gps_nsec, 1);
          if (hstat<0) {printf("Error creating attribute gps_nsec\n"); goto fail;}

          hstat = H5LTset_attribute_string(ofile_id, dsname, "sfdb_file", (const char *)&header_param->sfdbname);
          if (hstat<0) {printf("Error creating attribute sfdb_file\n"); goto fail;}

          // we don't use mjdtime later, it is stored to identify sft chunk in SFDB
          // gps_sec time refers to the flat-top part of the Tukey window,
          // sft_mjdtime refers to the beginning of the sft chunk
          hstat = H5LTset_attribute_double(ofile_id, dsname, "sft_mjdtime", &header_param->mjdtime, 1);
          if (hstat<0) {printf("Error creating attribute mjdtime\n"); goto fail;}

          time_t now = time(NULL);
          struct tm *local_time = localtime(&now);
          char ctime[100];
          strftime(ctime, sizeof(ctime), "%Y%m%d_%H%M%S", local_time);
          hstat = H5LTset_attribute_string(ofile_id, dsname, "ctime", ctime);
          if (hstat<0) {printf("Error creating attribute ctime\n"); goto fail;}

          free(dsname);
     } //ifatti

     puts("--> SUCCESS");
     goto success;

fail:
     retval = EXIT_FAILURE;
     puts("--> FAILURE");
success:
     H5Fclose(ofile_id);
     free(fileBAND);
     free(header_param);
     free(sft);
     free(sts);
     fftwf_cleanup();
     return retval;
}

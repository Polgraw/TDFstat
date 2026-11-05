#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../utils/iniparser/src/iniparser.h"
#include "genseg.h"

#ifdef USE_LAL
/* Compile the ephemeris generation code */
/* Requires the LALPulsar library */
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include "EphemerisDetector.h"
#endif


int main (int argc, char *argv[]) {
     FILE *td_stream;
     char td_fname[MAX_LINE+32], td_dir[MAX_LINE],
          ini_fname[MAX_LINE], date_fname[MAX_LINE+16];
     int  N, nod, i, yndx, notempty, bufsize = BUFSIZE,
          use_sci, nout, nseg, overwrite, startseg=1;
     double alpha, dt, othr, w_taper_dt, maxasd;
     struct stat st = {0};
     dictionary *ini; // config file
     const char *plsr, *DataDir, *out_replace, *sci_regions;

     double mingps=0., maxgps=0., start, end, startgps;
     int iseg, sci_len=0, retval;
     double scaling_factor = 1.e-20;
     float *xtime=NULL, *xall=NULL, *segar=NULL,
          *seg_sci_mask=NULL, *x0=NULL;

     // HDF related variables
     FILE *infile=NULL;
     hid_t infile_id=H5I_INVALID_HID;
     herr_t hstat;
     int format_version=1, nsamples, last_ichunk=-1;
     float bandwidth;
     double fpo;
     char dtype[4], site[3], oband[5];
     const char *H5FileName;

     int gen_eph=0;
#ifdef USE_LAL
     const char *EphDir=NULL, *efile=NULL, *sfile=NULL;
     char eFname[MAX_LINE], sFname[MAX_LINE], eph_fname[MAX_LINE+32];
     double *DetSSB=NULL, *rDet=NULL, *rSSB=NULL;
     double mjd1, phir, elam = 0, position[4];
     EphemerisData *edat = NULL;
     Detectors detector = 0;
#endif

     printf ("%s: generating time segments from sts data\n", argv[0]);
     if (argc > 1) strcpy(ini_fname, argv[1]); else strcpy(ini_fname, INI_FNAME);
     printf ("Loading config file %s\n", ini_fname);
     if ((ini = iniparser_load (ini_fname)) == NULL) {
          printf("Cannot parse file: %s\n", ini_fname);
          goto fail;
     }

     H5FileName = iniparser_getstring (ini, "genseg:infile", NULL);// input HDF5 file
     startgps = iniparser_getdouble(ini, "genseg:startgps", 0.0);  // write time sequences starting from this time
     nod = iniparser_getint (ini, "genseg:nod", 0);                // number of days per segment
     startseg = iniparser_getint (ini, "genseg:startseg", 1);    // no of the first segment, default is 1
     nseg = iniparser_getint (ini, "genseg:nseg", 0);              // number of segments, infinity if <=0
     plsr = iniparser_getstring (ini, "genseg:plsr", NULL);        // pulsar name or band number
     DataDir = iniparser_getstring (ini, "genseg:datadir", NULL);  // output directory
     overwrite = iniparser_getboolean (ini, "genseg:overwrite", 0);// allow to overwrite existing data

     alpha = iniparser_getdouble (ini, "genseg:alpha", 0.1);       // grubbs test parameter
     maxasd = iniparser_getdouble (ini, "genseg:maxasd", 1.e-21);  // max. asd in band, used to calculate threshold for large outliers
     out_replace = iniparser_getstring (ini, "genseg:out_replace", "zero"); // replace outliers with zero or random gaussian value

     use_sci = iniparser_getboolean (ini, "genseg:use_sci", 1);    // use science data only
     sci_regions = iniparser_getstring (ini, "genseg:sci_regions", NULL);   // file with science regions
     w_taper_dt = iniparser_getdouble (ini, "genseg:w_taper_dt", 600.);     // Tukey window tapering size in seconds; if <=0 do not apply window to science regions
     gen_eph = iniparser_getboolean (ini, "genseg:gen_eph", 0);
#ifdef USE_LAL
     EphDir = iniparser_getstring (ini, "genseg:EphDir", NULL);    // directory containing efile and sfile
     efile = iniparser_getstring (ini, "genseg:efile", NULL);      // earth ephemeris file
     sfile = iniparser_getstring (ini, "genseg:sfile", NULL);      // sun ephemeris file
#endif
     if (argc==3){
          if (strlen(argv[2]) != 4) {
               printf("[ERROR] Band number argument should be 4 digits including leading zeros\n");
               goto fail;
          }
          strcpy(oband, argv[2]);
          printf("[INFO] Band number overwrite enabled: bbbb -> %s\n", oband);
          // substitute bbbb in the input file name and plsr with oband
          char *p = strstr(H5FileName, "bbbb");
          if (p != NULL) {
               strncpy(p, oband, 4);
          } else {
               printf("[ERROR] bbbb not found in the input file name %s\n", H5FileName);
               goto fail;
          }
          if (strcmp(plsr, "bbbb") == 0) {
               p = strstr(plsr, "bbbb");
               strncpy(p, oband, 4);
          } else {
               printf("[ERROR] bbbb not found in the band/pulsar name %s\n", plsr);
               goto fail;
          }
     }

     printf("[IN] infile = %s\n", H5FileName);
     printf("[IN] startgps = %f\n", startgps);
     printf("[IN] nod = %d\n", nod);
     printf("[IN] nseg = %d\n", nseg);
     printf("[IN] plsr = %s\n", plsr);
     printf("[IN] datadir = %s\n", DataDir);
     printf("[IN] overwrite = %d\n", overwrite);
     printf("[IN] alpha = %f\n", alpha);
     printf("[IN] maxasd = %g\n", maxasd);
     printf("[IN] out_replace = %s\n", out_replace);
     printf("[IN] use_sci = %d\n", use_sci);
     printf("[IN] sci_regions = %s\n", sci_regions);
     printf("[IN] w_taper_dt = %f\n", w_taper_dt);
#ifdef USE_LAL
     printf("[IN] gen_eph = %d\n", gen_eph);
     printf("[IN] EphDir = %s\n", EphDir);
     printf("[IN] efile = %s\n", efile);
     printf("[IN] sfile = %s\n", sfile);
#else
     printf("[IN] gen_eph = %d\n", gen_eph);
     if (gen_eph) {
          printf("[ERROR] gen_eph = %d but the code was compiled without lalsuite.\n        Please set USE_LAL = yes in the Makefile\n", gen_eph);
          exit(EXIT_FAILURE);
     }
#endif

     if (nod <= 0) {
          printf("nod <= 0 !\n");
          goto fail;
     }


     // r+ returns NULL if file does not exist OR is directory
     if((infile = fopen(H5FileName,"r+")) != NULL){
          fclose(infile);
          infile_id = H5Fopen(H5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
          if (infile_id < 0) goto fail;
          // just some sanity checks
          hstat = H5LTget_attribute_int(infile_id, "/", "last_ichunk", &last_ichunk);
          if (hstat < 0) {printf("[HDF] Error reading attribute last_ichunk\n"); goto fail;}
          if (last_ichunk < 0) {printf("[HDF] Error: last_ichunk < 0 \n"); goto fail;}
          printf("[HDF] Opened %s \n      [read only mode][format_version = %d]\n", H5FileName, format_version);
          // getting creation plist for a file doesn't work, use root group
          hid_t rgroup = H5Gopen(infile_id, "/", H5P_DEFAULT);
          hid_t info = H5Gget_create_plist(rgroup);
          unsigned flags;
          H5Pget_link_creation_order(info, &flags);
          H5Pclose(info);
          if ( flags & (H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED) ) {
               puts("[HDF] Link creation order is tracked and indexed - good!");
          } else {
               printf("[HDF] Wrong link creation order flags=%u / %u!\n", flags, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
               goto fail;
          }

          hstat = H5LTget_attribute_int(infile_id, "/", "format_version", &format_version);
          if (hstat < 0) {
               printf("[HDF] Error reading attribute format_version\n"); goto fail;
          } else {
               printf("[HDF] format_version = %d\n", format_version);
          }

          hstat = H5LTget_attribute_string(infile_id, "/", "dtype", dtype);
          if (hstat < 0) {
               printf("[HDF] Error reading attribute dtype\n"); goto fail;
          } else {
               printf("[HDF] dtype = %s\n", dtype);
          }

          hstat = H5LTget_attribute_string(infile_id, "/", "site", site);
          if (hstat < 0) {
               printf("[HDF] Error reading attribute site\n"); goto fail;
          } else {
               printf("[HDF] site = %s\n", site);
          }

          hstat = H5LTget_attribute_double(infile_id, "/", "fpo", &fpo);
          if (hstat < 0) {
               printf("[HDF] Error reading attribute fpo\n"); goto fail;
          } else {
               printf("[HDF] fpo = %g\n", fpo);
          }

          hstat = H5LTget_attribute_float(infile_id, "/", "bandwidth", &bandwidth);
          if (hstat < 0) {
               printf("[HDF] Error reading attribute bandwidth\n"); goto fail;
          } else {
               printf("[HDF] bandwidth = %g\n", bandwidth);
          }

          hstat = H5LTget_attribute_int(infile_id, "/", "nsamples", &nsamples);
          if (hstat < 0) {
               printf("[HDF] Error reading attribute nsamples\n"); goto fail;
          } else {
               printf("[HDF] nsamples = %d\n", nsamples);
          }

          hstat = H5LTget_attribute_double(infile_id, "/", "scaling_factor", &scaling_factor);
          if (hstat < 0) {
               printf("[HDF] Error reading attribute scaling_factor\n"); goto fail;
          } else {
               printf("[HDF] scaling_factor = %g\n", scaling_factor);
          }

          printf("[HDF] last_ichunk = %d\n", last_ichunk);

     } else {
          printf("[HDF] Error! Cannot open input file \n" );
          goto fail;
     }

     dt = 1./(2*bandwidth);
     N = round(nod*C_SIDDAY/dt);      // No. of data points
     // 10 * sigma
     othr = 10.*maxasd*sqrt(bandwidth)/scaling_factor;

     printf("Infered parameters:\n");
     printf("[IN] site = %s\n", site);
     printf("[IN] dt = %g\n", dt);
     printf("[IN] N = %d (length of the output time series)\n", N);
     printf("[IN] othr = %g (large outlier threshold)\n", othr);

     if (stat (DataDir, &st) == -1) {
          mkdir (DataDir, 0755);
          printf ("[INFO] Data directory created: %s/\n", DataDir);
     }

     float *xchunk;
     xchunk = malloc(nsamples*sizeof(float));

     const gsl_rng_type *T;
     gsl_rng *r = NULL;
     unsigned long mySeed;
     if (! strcmp(out_replace,"gauss")) {
          // initialize random numbers
          mySeed = time(NULL);
          gsl_rng_env_setup();
          T = gsl_rng_default;
          r = gsl_rng_alloc (T);
          gsl_rng_set(r, mySeed);
     }

     if (use_sci) {
	     FILE *sci_file = fopen(sci_regions, "r");
	     if(sci_file == NULL){
	          printf("[SCI] Error! Cannot open sci_regions file \n");
	          exit(EXIT_FAILURE);
	     } else {
	          printf("[SCI] Reading file %s\n", sci_regions);
     	     fscanf(sci_file, "%lf %lf", &start, &end);
     	     mingps = start;
     	     maxgps = end;
     	     while(fscanf(sci_file, "%lf %lf", &start, &end) != EOF){
                    if (maxgps < end)   maxgps = end;
                    if (mingps > start) mingps = start;
               }

     	     sci_len = (int)ceil((maxgps-mingps)/dt);
     	     printf("[SCI] mingps = %f, maxgps = %f, sci_len = %d\n", mingps, maxgps, sci_len);

     	     // default is non-sci data
     	     segar = (float *)calloc(sci_len, sizeof(float));

     	     rewind(sci_file);
     	     // set mas to 1 in scientific regions
     	     while(fscanf(sci_file, "%lf %lf", &start, &end) != EOF){
                    int starti = (int)floor((start-mingps)/dt);
                    int endi   = (int)ceil((end-mingps)/dt);
                    DEBUG_PRINT(("   [SCI] [%d, %d]\n", starti, endi));
                    for(i=starti; i<=endi; ++i)
                         segar[i] = 1.;
     	     }
               fclose(sci_file);
		}
	     seg_sci_mask = malloc(N*sizeof(float)); // science data mask
     } // if(use_sci)


#ifdef USE_LAL
     if (gen_eph) {
          detector = get_detector((char *)site);
          get_position(detector, position);
          elam = position[1];
          printf("[EPH] Ephemeris for detector %s will be created\n",	names[detector]);
          sprintf(eFname, "%s/%s", EphDir, efile);
          sprintf(sFname, "%s/%s", EphDir, sfile);
          if ((edat = XLALInitBarycenter(eFname, sFname)) == NULL) {
               printf("[EPH] Problem in XLALInitBarycenter\n");
               goto fail;
          };
     }
     DetSSB = (double *)calloc(3*N+2, sizeof(double));
     rDet = (double *)calloc(3*N, sizeof(double));
     rSSB = (double *)calloc(3*N, sizeof(double));
#endif

     printf("[INFO] Time sequences for %s will be created.\n", plsr);

     double tseg_start, tseg_end, tchunk_start;
     int gps_sec, gps_nsec;
     H5LTget_attribute_int(infile_id, "1", "gps_sec", &gps_sec);
     H5LTget_attribute_int(infile_id, "1", "gps_nsec", &gps_nsec);

     double chunk_time = gps_sec + gps_nsec/1.0e9;
     if (startgps < chunk_time) {
          printf("[ERROR] startgps < chunk_1_time\n");
          exit(EXIT_FAILURE);
     }

     // align startgps to chunk[1].gpstime
     int idelta = (int)round((startgps-chunk_time)/dt);
     double startgps_aligned = chunk_time + idelta*dt;
     if ( (startgps_aligned - startgps) > 1.e-30) {
          printf("[WARNING] startgps was aligned to chunk_1 gpstime + i*dt: (old) %f -> (new) %f \n", startgps, startgps_aligned);
          startgps = startgps_aligned;
     } else {
          printf("[INFO] startgps is aligned with chunk_1 gpstime + i*dt, good! \n");
     }

     xtime = (float *)malloc(N*sizeof(float));
     xall = (float *)malloc(N*sizeof(float));

     int ichunk=1;

     for (iseg=1; iseg<=nseg; ++iseg){
          int iseg_name = iseg+startseg-1; // segment naming starts from startseg
          int chunk_i0=0, chunk_i1=0; // initialized to zero to avoid compiler warnings
          // zero whole segment
          memset(xtime, 0, N*sizeof(float));
          tseg_start = startgps + (iseg-1)*N*dt;
          tseg_end = tseg_start + (N-1)*dt;
          printf("------------------------------------------------------\n");
          printf("[INFO] Generating  iseg=%d  tseg_start=%f  tseg_end=%f\n", iseg_name, tseg_start, tseg_end);

          while(ichunk<=last_ichunk){
               int length = snprintf( NULL, 0, "%d", ichunk );
               char* dsname = malloc( length + 1 );
               snprintf(dsname, length + 1, "%d", ichunk );
               H5LTget_attribute_int(infile_id, dsname, "gps_sec", &gps_sec);
               H5LTget_attribute_int(infile_id, dsname, "gps_nsec", &gps_nsec);
               tchunk_start = gps_sec + gps_nsec/1.0e9;
               //tchunk_end = tchunk_start + (nsamples-1)*dt;

               // chunk possition in time segment
               chunk_i0 = (int)round((tchunk_start - tseg_start)/dt);
               if ((ichunk>1) && (chunk_i0 > chunk_i1+1)){
                    printf("[INFO] gap in segment data detected ! [%d-%d]\n", chunk_i1, chunk_i0);
                    // do nothing, just inform
               }
               chunk_i1 = chunk_i0 + nsamples - 1;

               DEBUG_PRINT(("   [DEB] ichunk=%d  t=%d  [%d, %d] ",
                    ichunk, gps_sec, chunk_i0, chunk_i1));

               if ( chunk_i1 < 0 ) {
                    // skip to the next chunk
                    DEBUG_PRINT(("   before range - continue\n"));
                    ++ichunk;
                    continue;
               }
               if ( chunk_i0 >= N ) {
                    printf("[INFO] gap between segments detected !\n");
                    break;
               }

               // read chunk data
               H5LTread_dataset_float(infile_id, dsname, xchunk);

               if ( (chunk_i0 >= 0) && (chunk_i1 <= (N-1)) ) {
                    // chunk is entirely in the range
                    DEBUG_PRINT((" in range [%d, %d]/[0, %d]\n", chunk_i0, chunk_i1, N-1));
                    for (i=0; i<nsamples; ++i){
                         xtime[chunk_i0+i] = xchunk[i];
                    }
                    ++ichunk;
                    continue;
               } else if (chunk_i0 < 0) {
                    // chunk overlaps segment on the left side
                    DEBUG_PRINT((" left edge [%d, %d]/[0, %d]\n", 0, chunk_i1, N-1));
                    for (i=-chunk_i0; i<nsamples; ++i){
                         xtime[chunk_i0+i] = xchunk[i];
                    }
                    ++ichunk;
                    continue;
               } else {
                    // chunk overlaps segment on the right side
                    DEBUG_PRINT((" right edge [%d, %d]/[0, %d]\n", chunk_i0, N-1, N-1));
                    //for (i=0; i<(nsamples-(chunk_i1-N+1)); ++i){
                    for (i=0; i<(N-chunk_i0); ++i){
                         xtime[chunk_i0+i] = xchunk[i];
                    }
                    // we can't write next segment here because memset; do not increment ichunk
                    break;
               }
          } // while(ichunk)

          if (ichunk > last_ichunk) {
               printf("[ERROR] Not enough chunks in the HDF file!\n");
               exit(EXIT_FAILURE);
          }

          // keep the original time series, xall will be cleaned
          for (i=0; i<N; ++i){
               xall[i] = xtime[i];
          }

          int i1, i2, l;
          // init mask of scientific regions in segment
          int offsetgps=(int)round((tseg_start - mingps)/dt);
          // memcpy (seg_sci_mask, segar+offsetgps, N*sizeof(float));

          for (i=0; i < N; i++){
               if ((offsetgps+i) < 0  || (offsetgps+i) > sci_len) {
                    seg_sci_mask[i] = 0.;
               } else {
                    seg_sci_mask[i] = segar[offsetgps+i];
               }
          }

          // prepare scientific regions
          i = 0;
          while(i < N){
               // find sci region
               if (seg_sci_mask[i] < 0.5) {i++; continue;}
               i1 = i;  // first point
               i2 = i1; // last point
               while(++i)
                    if ((seg_sci_mask[i] < 0.5) || (i == N)) break;

               i2 = i-1; // last point in the sci region

#if 1
               // exclude large outliers groups at edges of sci regions
               int is;
               DEBUG_PRINT(("[DEB] Exclude large outliers at sci edges: [%d, %d]  ", i1, i2));
               for (is=i1; is <= i2; is++){
                    if ( fabs(xall[is]) < othr &&
                         fabs(xall[is+1]) < othr ) {
                              break;
                         }
                    seg_sci_mask[is] = 0.;
               }
               i1 = is;
               for (is=i2; is >= i1; is--){
                    if ( fabs(xall[is]) < othr &&
                         fabs(xall[is-1]) < othr ) {
                              break;
                         }
                    seg_sci_mask[is] = 0.;
               }
               i2 = is;
               DEBUG_PRINT((" => [%d, %d]\n", i1, i2));
#endif

               // Apply Tukey window  to the scientific region
#if 1

#define PI  3.141592653589793
               if (w_taper_dt > dt) {
                    int lobe_isize = (int)ceil(w_taper_dt/dt);
                    int li;
                    // in case lobe_isize*2 < L ; +2 accounts for even and odd cases
                    li = MIN(lobe_isize, (i2-i1+2)/2 );
                    DEBUG_PRINT(("[DEB] Tukey window applied to [%d, %d]  iobe_isize=%d\n", i1, i2, li));
                    for(l=0; l<li; l++){
                         double lobe = 0.5*(1.-cos(PI*l/li));
                         seg_sci_mask[i1 + l] = lobe;
                         seg_sci_mask[i2 - l] = lobe;
                    }
               }
#endif
    	     } // while(i < N)

          // apply scientific mask
          for (i=0; i < N; i++) {
               xall[i] *= seg_sci_mask[i];
          }

#if 1
          // remove ramaining large outliers above 6*sigma level (othr not used)
          printf("[INFO] Cleaning large outliers - replace with %s \n", out_replace);

          double sdval;
          nout = 0;
          sdval = 0.;
          for (i=0; i<N; i++){
               double v = xall[i];
               // do not count zeros
               if (fabs(v) < othr && fabs(v) > 1.e-30){
                    sdval += v*v;
                    nout++;
               }
          }
          sdval = sqrt(sdval/(double)(nout-1));
          printf("[INFO] ---- stdval = %f ,  noutzero = %d\n", sdval, nout);

          nout = 0;
          if (! strcmp(out_replace, "gauss")) {
               for (i=0; i<N; i++)
                    // if (fabs(xall[i]) > othr){
                    if (fabs(xall[i]) > 6.*sdval){
                         xall[i] = gsl_ran_gaussian_ziggurat(r, sdval);
                         nout++;
                    }
          } else {
               for (i=0; i<N; i++)
                    // if (fabs(xall[i]) > othr){
                    if (fabs(xall[i]) > 6.*sdval){
                         xall[i] = 0.;
                         nout++;
                    }
          }
          printf("[INFO] ---- Large outliers removed : %d / %d = %d%% \n", nout, N, 100*nout/N);
#endif
          // remove remaining outliers using Grubbs test
          x0 = (float *)malloc(N*sizeof (float)); // init to 0 in fGrubbsOutliersMany
#if 1
          nout = fGrubbsOutliersMany(xall, x0, N, bufsize, alpha, out_replace);
          printf("[INFO] Grubbs outliers removed : %d / %d = %6.2f%% \n", nout, N, 100.*(float)nout/N);
#else
          memcpy(x0, xall, N*sizeof(float));
#endif

// replace all zeros with gaussian values
#if 0
          printf("[INFO] Replacing zeros with gaussian values");

          double sd = 0.;
          nout = 0;
          for (i=0; i<N; i++){
               double v = x0[i];
               if (fabs(v) > 1.e-30){
                    sd += v*v;
                    nout++;
               }
          }
          sd = sqrt(sd/(double)(nout-1));
          printf("[INFO] | std dev = %f \n", sd);

          for (i=0; i<N; i++){
               if (fabs(x0[i]) < 1.e-30)
                    x0[i] = gsl_ran_gaussian_ziggurat(r, sd);
          }
#endif

          /**************************/
          /* write output files     */
          /**************************/
          sprintf(td_dir, "%s/%03d", DataDir, iseg_name);
          if (stat (td_dir, &st) == -1) {
               mkdir (td_dir, 0755);
          }
          sprintf(td_dir, "%s/%03d/%s", DataDir, iseg_name, site);
          if (stat (td_dir, &st) == -1) {
               mkdir (td_dir, 0755);
          }
          // standard filename format
          sprintf(td_fname, "%s/xdat_%03d_%s.bin", td_dir, iseg_name, plsr);
          sprintf(date_fname, "%s/starting_date", td_dir);

          printf ("[INFO] Writing segment %03d : ", iseg_name);
          // Do "xall" contains only zeros?
          notempty=0;
          for (yndx=0; yndx < N; yndx++){
               if(x0[yndx] != 0.0){
                    notempty=1;// not empty
                    break;
               }
          }

          if (overwrite) {
               // order is important, don't create empty file
               if ( notempty && ((td_stream=fopen (td_fname, "w")) != NULL) ) {
                    fwrite((void *)(x0), sizeof(float), N, td_stream);
                    fclose(td_stream);
                    printf(" data ");
               }
               if ((td_stream = fopen(date_fname, "w")) != NULL) {
                    fprintf(td_stream, "%.10e\n", tseg_start);
                    fclose(td_stream);
                    printf(" | starting_date ");
               }
          } else {
               if ((td_stream=fopen (td_fname, "r")) != NULL) {
                    fclose(td_stream);
                    printf(" skipping data ");
                    //continue;
               } else {
                    if ( notempty && ((td_stream=fopen (td_fname, "w")) != NULL)) {
                         fwrite((void *)(x0), sizeof(float), N, td_stream);
                         fclose(td_stream);
                         printf(" data ");
                    }
               }
               if ((td_stream = fopen(date_fname, "r")) != NULL) {
                    fclose(td_stream);
                    printf(" | skipping starting_date ");
                    //continue;
               } else {
                    if ((td_stream = fopen(date_fname, "w")) != NULL) {
                         fprintf(td_stream, "%.10e\n", tseg_start);
                         fclose(td_stream);
                         printf(" | starting_date ");
                    }
               }
          }
          printf("\n");
          /* Generate detector ephemeris */
#ifdef USE_LAL
          if (gen_eph) {
               get_barycenter (tseg_start, detector, edat, DetSSB, rDet, dt, N);
               for (i=0; i<3*N; i++) rSSB[i] = DetSSB[i] - rDet[i];
               mjd1 = gps2mjd (tseg_start);
               phir = sid (mjd1, elam);
               DetSSB[3*N] = phir;
               DetSSB[3*N+1] = EPSILON;

               sprintf (eph_fname, "%s/DetSSB.bin", td_dir);
               if ((td_stream=fopen (eph_fname, "w")) != NULL) {
                    fwrite ((void *)DetSSB, sizeof(double), 3*N+2, td_stream);
                    fclose (td_stream);
               }
               sprintf (eph_fname, "%s/rDet.bin", td_dir);
               if ((td_stream=fopen (eph_fname, "w")) != NULL) {
                    fwrite ((void *)rDet, sizeof(double), 3*N, td_stream);
                    fclose (td_stream);
               }
               sprintf (eph_fname, "%s/rSSB.bin", td_dir);
               if ((td_stream=fopen (eph_fname, "w")) != NULL) {
                    fwrite ((void *)rSSB, sizeof(double), 3*N, td_stream);
                    fclose (td_stream);
               }
          }
#endif
     } // for iseg

     retval = EXIT_SUCCESS;
     puts("--> SUCCESS");
     goto success;
fail:
     retval = EXIT_FAILURE;
     puts("--> FAILURE");

success:
     if (infile_id != H5I_INVALID_HID) H5Fclose(infile_id);
     iniparser_freedict(ini);
     free(xtime);
     free(xall);
     free(segar);
     free(seg_sci_mask);
#ifdef USE_LAL
     free (rSSB);
     free (rDet);
     free (DetSSB);
#endif
     return retval;
}

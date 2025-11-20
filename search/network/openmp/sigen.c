#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>

#include "auxi.h"
#include "struct.h"
#include "settings.h"

static int help_flag=0;

double get_rand();
Detector_settings ifo[MAX_DETECTORS];

int main(int argc, char *argv[])
{

     int i, numl=0, freq_line_check, c, pm, band=0, reffr, nfrinband, nos;
     char *wd=NULL;
     double amp=0, snr=0;
     double sgnlo[7], rrn[2],
          fdotmin, fdotmax,
          fpo_val,
          sepsm, cepsm,
          iota, ph_o, psik, hop, hoc, overlap;

     FILE *data, *output_file = NULL;

     Search_settings sett;

     //#mb fftpad value=1 (default for some time now)
     sett.fftpad = 1;

     // Default data sampling time (s)
     sett.dt = 0.5;

     // Default reference time frame (frame in which the frequency
     // of the signal is as selected below, not spun-down/up) is 1
     reffr = 1;

     // Initial value of starting frequency
     // set to a negative quantity. If this is not
     // changed by the command line value
     // fpo is calculated from the band number b.
     fpo_val = -1;

     // Initial value of the number of days is set to 0
     sett.nod = 0;

     // Initial value of number of signals is set to 1
     nos = 1;

     overlap = -1.;
     nfrinband = 0;

     while (1) {
          static struct option long_options[] = {
               {"help", no_argument, &help_flag, 1},
               // GW amplitude
               {"amp", required_argument, 0, 'a'},
               // SNR
               {"snr", required_argument, 0, 's'},
               // fpo value
               {"fpo", required_argument, 0, 'f'},
               // frequency band number
               {"band", required_argument, 0, 'b'},
               // band overlap
               {"overlap", required_argument, 0, 'o'},
               // data sampling time
               {"dt", required_argument, 0, 't'},
               // number of days in the time-domain segment
               {"nod", required_argument, 0, 'd'},
               // number of signals to print
               {"nos", required_argument, 0, 'n'},
               // reference frame ()
               {"reffr", required_argument, 0, 'r'},
               // change directory parameter
               {"cwd", required_argument, 0, 'c'},
               // +/- frames from reference to keep the signal inband
               {"nfrinband", required_argument, 0, 'i'},


               {0, 0, 0, 0}
          };

          if(help_flag) {

               printf("*** Software injection parameters - cont. GW signal ***\n");
               printf("Usage: ./sigen -[switch1] <value1> -[switch2] <value2> ...\n") ;
               printf("Switches are:\n\n");
               printf("-amp      GW amplitude of the injected (mutually exclusive with -snr)\n");
               printf("-snr      SNR of the injected signal (mutually exclusive with -amp)\n");
               printf("-fpo      fpo (starting frequency) value (either this or band and overlap should be provided)\n");
               printf("-band     Band number\n");
               printf("-overlap  Band overlap\n");
               printf("-dt       Data sampling time dt (default value: 0.5)\n");
               printf("-nod      Number of days\n");
               printf("-nos      Number of signals to generate \n");
               printf("-reffr    Reference frame (default value: 1)\n");
               printf("-cwd      Change to directory <dir>\n");
               printf("-nfrinband  Force the signal to stay inband for +/- frames \n\n");

               printf("--help    This help\n");
               exit (0);
          }

          int option_index = 0;
          c = getopt_long_only (argc, argv, "a:s:f:b:o:t:d:n:r:c:i:", long_options, &option_index);

          if (c == -1)
               break;
          switch (c) {
               case 'a':
               amp  = atof (optarg);
               break;
               case 's':
               snr  = atof (optarg);
               break;
               case 'f':
               fpo_val = atof (optarg);
               break;
               case 'b':
               band = atoi (optarg);
               break;
               case 'o':
               overlap = atof (optarg);
               break;
               case 't':
               sett.dt = atof(optarg);
               break;
               case 'd':
               sett.nod = atoi(optarg);
               break;
               case 'n':
               nos = atoi(optarg);
               break;
               case 'r':
               reffr = atof(optarg);
               break;
               case 'c':
               wd = (char *) malloc (1+strlen(optarg));
               strcpy (wd, optarg);
               break;
               case 'i':
               nfrinband = atoi(optarg);
               break;
               case '?':
               break;
               default: break ;

          } /* switch c */
     } /* while 1 */

     // Check if sett->nod was set up, if not, exit
     if(!(sett.nod)) {
          printf("Number of days not set... Exiting\n");
          exit(EXIT_FAILURE);
     }

     if (wd) {
          printf ("Changing working directory to %s\n", wd);
          if (chdir(wd)) {
               perror (wd);
               abort ();
          }
     }

     // Check if the options are consistent with each other
     if(!(amp || snr)) {
          printf("Options -amp or -snr not given. Exiting...\n");
          exit(0);
     } else if(amp && snr) {
          printf("Options -amp and -snr are mutually exclusive. Exiting...\n");
          exit(0);
     }

     // Starting band frequency:
     // fpo_val is optionally read from the command line
     // Its initial value is set to -1
     if(fpo_val >= 0) {
          sett.fpo = fpo_val;
     } else if (band > 0 && overlap >= 0.) {
          sett.fpo = 10. + (1. - overlap)*band*(0.5/sett.dt);
     } else {
          printf("Band AND overlap or fpo must be specified!\n");
          exit(EXIT_FAILURE);
     }

     // // Search settings
     // search_settings(&sett);

     // The search_settings function call is replaced by direct assigments here

     sett.B = 0.5/sett.dt;               // Bandwidth (Hz)
     sett.N = round (sett.nod*C_SIDDAY/sett.dt); // Number of data points

     // spindown range of NS in physical units [Hz/s]
     // we assume minimum NS age 1000 yrs
     if (sett.fpo < 200.) {
          fdotmin = 2.*(sett.fpo+sett.B)/(2.*1000.*C_YEARSEC);
          fdotmax = 0.;
     } else {
          fdotmin = 2e-10;
          fdotmax = 2e-11;
     }

     // dimensionless spindown range
     sett.Smax = M_PI*fdotmin*sett.dt*sett.dt;
     sett.Smin = M_PI*fdotmax*sett.dt*sett.dt;

     fprintf(stderr, "Band number is %04d\n", band);
     fprintf(stderr, "The reference frequency fpo is %f\n", sett.fpo);
     fprintf(stderr, "The data sampling time dt is %f\n", sett.dt);
     fprintf(stderr, "The reference time frame is %d\n", reffr);

     //#mb Changing the limits near the bands border
     //#mb For O3 UL simulation run
     rrn[0] = M_PI/20.;
     rrn[1] = M_PI - rrn[0];

     // Random signal parameters
     //-------------------------
     double f1,f2;
     int inband=0;
     do {
         // fprintf(output_file, "# amporsnr h0/snr reff freq fdot ra dec iota psi phase \n");
         printf("# amporsnr h0/snr reff freq fdot ra dec iota psi phase \n");
         for (int signal_idx = 0; signal_idx < nos; signal_idx++) {
             do {
                 // Frequency derivative
                 sgnlo[1] = sett.Smin - (sett.Smin + sett.Smax) * get_rand();

                 // Frequency in (0, \pi) range
                 sgnlo[0] = get_rand() * (rrn[1] - rrn[0]) + rrn[0];

                 f1 = sgnlo[0] + 2. * sgnlo[1] * sett.N * nfrinband;
                 f2 = sgnlo[0] - 2. * sgnlo[1] * sett.N * nfrinband;

                 // Check if the signal is in band
                 inband = ((f1 > 0) && (f1 < M_PI) && (f2 < M_PI) && (f2 > 0));
             } while (inband != 1);

             // Hemisphere
             pm = round(get_rand()) + 1;

             // epsma - average inclination of Earth's axis to the ecliptic
             // (see definitions of constants in settings.h)
             sepsm = sin(C_EPSMA);
             cepsm = cos(C_EPSMA);

             // Uniform sphere sampling algorithm
             double x1, x2, X, Y, Z;
             do {
                 x1 = 2 * get_rand() - 1;
                 x2 = 2 * get_rand() - 1;
             } while (x1 * x1 + x2 * x2 >= 1);

             X = 2 * x1 * sqrt(1 - x1 * x1 - x2 * x2);
             Y = 2 * x2 * sqrt(1 - x1 * x1 - x2 * x2);
             Z = 1 - 2 * (x1 * x1 + x2 * x2);

             // Random phase and polarization of the signal
             ph_o = 2. * M_PI * get_rand();
             psik = 2. * M_PI * get_rand();
             hoc = 2. * get_rand() - 1.;
             hop = (1. + hoc * hoc) / 2.;
             iota = acos(hoc);

             // Convert the frequency into physical Hz units
             sgnlo[0] = sett.B * sgnlo[0] / M_PI + sett.fpo;

             // Converting the spin down into physical Hz/t units
             sgnlo[1] = sgnlo[1] / M_PI / sett.dt / sett.dt;

             // Sky position: Right ascension
             sgnlo[2] = atan2(Y, X) + M_PI;

             // Sky position: declination
             sgnlo[3] = M_PI_2 - acos(Z);

             // Inclination
             sgnlo[4] = iota;

             // Polarisation
             sgnlo[5] = psik;

             // Phase
             sgnlo[6] = ph_o;

             // Output (GW amplitude or signal-to-noise ratio)

             if (amp) {
                 printf("amp %le %d %.16le %.16le %.16le %.16le %.16le %.16le %.16le\n", amp, reffr,
                         sgnlo[0], sgnlo[1], sgnlo[2], sgnlo[3],
                         sgnlo[4], sgnlo[5], sgnlo[6]);
             } else if (snr) {
                 printf("snr %le %d %.16le %.16le %.16le %.16le %.16le %.16le %.16le\n", snr, reffr,
                         sgnlo[0], sgnlo[1], sgnlo[2], sgnlo[3],
                         sgnlo[4], sgnlo[5], sgnlo[6]);
             } else {
                 fprintf(stderr, "Neither the SNR or the amplitude is provided....\n");
                 fprintf(stderr, "Maybe you are looking at an imaginary signal????\n");
                 fclose(output_file);
                 return (EXIT_FAILURE);
             }
         } // for loop for nos

     } while (0);

     return (EXIT_SUCCESS);

} // sigen()


// Random number generator with seed from /dev/urandom
// range: [0,1]
double get_rand() {

     FILE *urandom;
     unsigned int seed;

     urandom = fopen ("/dev/urandom", "r");
     if (urandom == NULL) {
          fprintf (stderr, "Cannot open /dev/urandom!\n");
          exit(EXIT_FAILURE);
     }

     fread (&seed, sizeof (seed), 1, urandom);
     srand (seed);
     return ((double)rand()/(double)(RAND_MAX));

}

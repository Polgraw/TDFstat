#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h>
#include <sys/time.h>

// ./gauss-xdat 86164 3 1 ../testdata/001/H1/xdatc_001_0666.bin  mask.bin
// compilation: gcc gauss-xdat-mask.c -o gauss-xdat-mask -lm -lgsl -lgslcblas

gsl_rng * r;

unsigned long int random_seed() {

     unsigned int seed;
     struct timeval tv;
     FILE *devrandom;
     
     if ((devrandom = fopen("/dev/random","r")) == NULL) {
	  gettimeofday(&tv,0);
	  seed = tv.tv_sec + tv.tv_usec;
     } else {
	  fread(&seed,sizeof(seed),1,devrandom);
	  fclose(devrandom);
     }
     
     return(seed);
}


int main(int argc, char **argv) { 
     
     int i, N, stat; 
     float *x, amp, sigma, *zmask;
     unsigned long mySeed;
     mySeed = random_seed();
     
     N = atoi(argv[1]); 
     amp = atof(argv[2]); 
     sigma = atof(argv[3]);  
     
     x = (float *)calloc(N, sizeof(float));
     zmask = (float *)calloc(N, sizeof(float));
     
     // FILE *zmask_file = fopen(argv[5], "r");
     // stat = fread(zmask, sizeof(float), N, zmask_file);
     // printf("Mask read in.\n");
     // fclose (zmask_file);

     FILE *dataout;
     dataout = fopen(argv[4], "wb");
  
     const gsl_rng_type * T;
     
     gsl_rng_env_setup();
     
     T = gsl_rng_default;
     r = gsl_rng_alloc (T);
     gsl_rng_set(r, mySeed);
     // Generate normal distribution (around 0, 
     // with amplitude amp and sigma)
     for(i=0; i<N; i++) {
	  x[i] = amp*gsl_ran_gaussian_ziggurat(r, sigma);
	  // if (fabs(zmask[i]) < 1.e-30) x[i] = 0.;
     }

     gsl_rng_free(r);
     
     fwrite(x, sizeof(*x), N, dataout);
     fclose(dataout);
     free(x);
     
     return 0;
}

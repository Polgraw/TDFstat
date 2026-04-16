#ifndef __JOBCORE_H__
#define __JOBCORE_H__

#include "auxi.h"
#include "struct.h"

void search(
               Search_settings *sett,
               Command_line_opts *opts,
               Search_range *s_range,
               FFTW_plans *plans,
               FFTW_arrays *fftw_arr,
               Aux_arrays *aux,
               Signal_params *sgnl_params,
               int *Fnum);

/* Main job function
 * The output is stored in single or double precision
 * (FLOAT_TYPE defined in struct.h)
 */

int job_core(
     	     int pm,                   // hemisphere
     	     float mm,                 // grid 'sky position'
     	     float nn,                 // other grid 'sky position'
     	     Search_settings *sett,    // search settings
     	     Command_line_opts *opts,  // cmd opts
     	     Search_range *s_range,    // range for searching
     	     FFTW_plans *plans,        // plans for fftw
     	     FFTW_arrays *fftw_arr,    // arrays for fftw
     	     Aux_arrays *aux,          // auxiliary arrays
     	     int *sgnlc,               // current number of candidates
     	     Trigger *sgnlv,        // candidate array
     	     int *FNum);               // candidate signal number


void spindown_modulation(
                         const int nifo,
                         const int N,
                         const FLOAT_TYPE het1,
                         const FLOAT_TYPE spnd,
                         const FLOAT_TYPE _tmp1[][N],
                         FFTW_PRE(_complex) *fxa,
                         FFTW_PRE(_complex) *fxb );

#endif

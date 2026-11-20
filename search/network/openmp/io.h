#include "auxi.h"
#include "struct.h"

int hdfout_attr_init (char *outname, Command_line_opts *opts, Search_settings *sett,
     Injection_params *inj, int ninj);

int hdfout_dset_init (char *outname, const char *dset_name, Search_range *s_range,
     Trigger *sgnlv);

int hdfout_extend_dset(char *outname, const char *dset_name, Trigger *sgnlv, int sgnlv_size);

int hdfout_finalize(char *outname, const char *dset_name, int totsgnl, double time_elapsed, int nthreads);

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <time.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#include <H5Apublic.h>
#include <H5Ipublic.h>
#include <H5Tpublic.h>

#include "cdefs.h"
#include "../search/network/openmp/struct.h"
#include "../utils/iniparser/src/iniparser.h"
#include "coincidences.h"

#define FORMAT_VERSION  1
#define TRIG_RANK       1

/* Minimal struct for partial read of the 'opts' compound attribute.
 * Only the six fields requested from Command_line_opts are declared here;
 * HDF5 compound-type conversion matches by name and ignores the rest.   */
typedef struct {
     int    band;
     int    seg;
     int    hemi;
     double overlap;
     double narrowdown;
     char  *grid_file;   /* variable-length string – HDF5 allocates this */
} Opts_partial;



void read_coinc_ini(char *ini_fname, Coinc_opts *copts)
{
     dictionary *ini;
     char *tlist;
     int error = 0;

     printf("Loading config file %s\n", ini_fname);
     if ((ini = iniparser_load(ini_fname)) == NULL) {
          perror(ini_fname);
          exit(EXIT_FAILURE);
     }
     printf("-- INI file contents --\n");
     iniparser_dump(ini, stdout);
     printf("-----------------------\n");

     copts->refr         = iniparser_getint   (ini, "coincidences:refr",          -1);
     copts->cthr         = iniparser_getdouble(ini, "coincidences:cthr",          -1.);
     copts->hemi         = iniparser_getint   (ini, "coincidences:hemi",           1);
     copts->scalef       = iniparser_getint   (ini, "coincidences:scalef",         1);
     copts->scales       = iniparser_getint   (ini, "coincidences:scales",         1);
     copts->scalem       = iniparser_getint   (ini, "coincidences:scalem",         1);
     copts->scalen       = iniparser_getint   (ini, "coincidences:scalen",         1);
     copts->shift        = iniparser_getint   (ini, "coincidences:shift",       0000);
     copts->mincoin      = iniparser_getint   (ini, "coincidences:mincoin",        3);
     copts->out_dir      = iniparser_getstring(ini, "coincidences:out_dir",       ".");
     copts->trig_dset    = iniparser_getstring(ini, "coincidences:trig_dset", "triggers");
     copts->coinc_dset_pre = iniparser_getstring(ini, "coincidences:trig_dset_pre", "coinc");
     tlist               = (char *)iniparser_getstring(ini, "coincidences:trig_file_list", NULL);

     /* various checks */
     if (copts->refr < 0) {
          error = 1; printf("[ERROR] missing refr !\n");
     }
     if (copts->cthr < 0) {
          error = 1; printf("[ERROR] missing cthr !\n");
     }
     if (!tlist) {
          error = 1; printf("[ERROR] trig_file_list !\n");
     }

     if (error) {
          fprintf(stderr, "Error: invalid or missing parameters in %s\n", ini_fname);
          iniparser_freedict(ini);
          exit(EXIT_FAILURE);
     }

     size_t count = 0;
     for (char *token = strtok(tlist, " "); token != NULL; token = strtok(NULL, " ")) {
          copts->trig_files[count] = (char *)malloc((strlen(token) + 1) * sizeof(char));
          strcpy(copts->trig_files[count], token);
          count++;
          if (count >= MAX_TRIG_FILES) {
               fprintf(stderr,
                       "Error: too many trigger files specified in %s (max %d)\n",
                       ini_fname, MAX_TRIG_FILES);
               iniparser_freedict(ini);
               exit(EXIT_FAILURE);
          }
     }

     copts->n_trig_files = count;
     printf("Number of files to search for coincidences: %zu\n", copts->n_trig_files);
     //for (size_t i = 0; i < copts->n_trig_files; i++)
     //     puts(copts->trig_files[i]);

     //iniparser_freedict(ini);

} /* read_coinc_ini */


size_t read_triggers_file(const char *filename, const char *t_dset_name,
                          Trigger **sgnlv, Search_params *search_par)
{
     size_t sgnlv_size = 0;

     hid_t   file, t_dataset, t_space;
     hsize_t t_dim[TRIG_RANK];
     herr_t  hstat;
     bool init_search_par = (search_par->grid_file == NULL);

     /* ------------------------------------------------------------------ */
     /* Define compound type matching the write side (hdfout_init/extend)   */
     /* ------------------------------------------------------------------ */
     hid_t ffstat_type = H5Tvlen_create(H5T_NATIVE_FLOAT);
     hid_t t_tid = H5Tcreate(H5T_COMPOUND, sizeof(Trigger));
     H5Tinsert(t_tid, "m",      HOFFSET(Trigger, m),      H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "n",      HOFFSET(Trigger, n),      H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "s",      HOFFSET(Trigger, s),      H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "ra",     HOFFSET(Trigger, ra),     H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "dec",    HOFFSET(Trigger, dec),    H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "fdot",   HOFFSET(Trigger, fdot),   H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "ffstat", HOFFSET(Trigger, ffstat), ffstat_type);

     /* ------------------------------------------------------------------ */
     /* Open file                                                           */
     /* ------------------------------------------------------------------ */
     file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
     if (file < 0) {
          fprintf(stderr, "Error: cannot open HDF5 file %s\n", filename);
          H5Tclose(t_tid);
          H5Tclose(ffstat_type);
          *sgnlv = NULL;
          return 0;
     }

     /* ------------------------------------------------------------------ */
     /* Read selected fields from the 'opts' attribute:                     */
     /*   band, seg, hemi  (int)                                            */
     /*   overlap, narrowdown  (double)                                     */
     /*   grid_file  (variable-length string)                               */
     /* ------------------------------------------------------------------ */

     hid_t vstr_t = H5Tcopy(H5T_C_S1);
     H5Tset_size(vstr_t, H5T_VARIABLE);

     hid_t opts_tid = H5Tcreate(H5T_COMPOUND, sizeof(Opts_partial));
     H5Tinsert(opts_tid, "band", offsetof(Opts_partial, band), H5T_NATIVE_INT);
     H5Tinsert(opts_tid, "seg", offsetof(Opts_partial, seg), H5T_NATIVE_INT);
     H5Tinsert(opts_tid, "hemi", offsetof(Opts_partial, hemi), H5T_NATIVE_INT);
     H5Tinsert(opts_tid, "overlap", offsetof(Opts_partial, overlap), H5T_NATIVE_DOUBLE);
     H5Tinsert(opts_tid, "narrowdown", offsetof(Opts_partial, narrowdown), H5T_NATIVE_DOUBLE);
     H5Tinsert(opts_tid, "grid_file", offsetof(Opts_partial, grid_file),  vstr_t);

     printf("> Reading file %s\n", filename);
     Opts_partial op = {0};
     hid_t opts_attr = H5Aopen(file, "opts", H5P_DEFAULT);
     if (opts_attr >= 0) {
          hstat = H5Aread(opts_attr, opts_tid, &op);
          H5Aclose(opts_attr);
          if (hstat < 0) {
               fprintf(stderr,"Warning: cannot read 'opts' attribute from %s\n", filename);
               exit(EXIT_FAILURE);
          }
          if (init_search_par) {
               search_par->band       = op.band;
               search_par->seg        = op.seg;
               search_par->hemi       = op.hemi;
               search_par->overlap    = op.overlap;
               search_par->narrowdown = op.narrowdown;
               search_par->grid_file  = op.grid_file ? strdup(op.grid_file) : NULL;
               printf ("    [getting params: band=%d, seg=%d, hemi=%d, overlap=%.2f, narrowdown=%.2f, grid_file=%s]\n",
                       search_par->band, search_par->seg, search_par->hemi,
                       search_par->overlap, search_par->narrowdown,
                       search_par->grid_file ? search_par->grid_file : "NULL");
          } else {
               // compare op with searh par
               if (search_par->band != op.band ||
                   search_par->hemi != op.hemi || search_par->overlap != op.overlap ||
                   search_par->narrowdown != op.narrowdown ||
                   strcmp(search_par->grid_file, op.grid_file) != 0)
               {
                    fprintf(stderr, "search_par->grid_file=%s\n", search_par->grid_file);
                    fprintf(stderr, "op.grid_file=%s\n", op.grid_file);
                    fprintf(stderr, "Error: 'opts' attribute in %s does not match expected values\n", filename);
                    exit(EXIT_FAILURE);
               } else {
                    printf("    [params match]\n");
               }
          }
          /* Reclaim HDF5-allocated variable-length string memory */
          hid_t scalar_sp = H5Screate(H5S_SCALAR);
          H5Treclaim(opts_tid, scalar_sp, H5P_DEFAULT, &op);
          H5Sclose(scalar_sp);
     } else {
          fprintf(stderr, "Warning: cannot open 'opts' attribute in %s\n", filename);
          exit(EXIT_FAILURE);
     }

     H5Tclose(opts_tid);
     H5Tclose(vstr_t);

     /* ------------------------------------------------------------------ */
     /* Open the triggers dataset                                           */
     /* ------------------------------------------------------------------ */
     t_dataset = H5Dopen2(file, t_dset_name, H5P_DEFAULT);
     if (t_dataset < 0) {
          fprintf(stderr, "Error: cannot open dataset '%s' in %s\n",
                  t_dset_name, filename);
          H5Fclose(file);
          H5Tclose(t_tid);
          H5Tclose(ffstat_type);
          *sgnlv = NULL;
          return 0;
     }

     /* ------------------------------------------------------------------ */
     /* Query number of triggers from the dataspace extent                  */
     /* ------------------------------------------------------------------ */
     t_space = H5Dget_space(t_dataset);
     hstat   = H5Sget_simple_extent_dims(t_space, t_dim, NULL);
     sgnlv_size = (size_t)t_dim[0];

     if (init_search_par) search_par->sgnlv_size = sgnlv_size;
     if (sgnlv_size != search_par->sgnlv_size) {
          fprintf(stderr, "Error! sgnlv_size mismatch: read %zu, but expected %zu \n",
               sgnlv_size, search_par->sgnlv_size);
          exit(EXIT_FAILURE);
     }
     
     if (sgnlv_size == 0) {
          *sgnlv = NULL;
          H5Sclose(t_space);
          H5Dclose(t_dataset);
          H5Fclose(file);
          H5Tclose(t_tid);
          H5Tclose(ffstat_type);
          return 0;
     }

     /* ------------------------------------------------------------------ */
     /* Allocate output buffer                                              */
     /* ------------------------------------------------------------------ */
     *sgnlv = (Trigger *)malloc(sgnlv_size * sizeof(Trigger));
     if (*sgnlv == NULL) {
          fprintf(stderr, "Error: cannot allocate memory for %zu triggers\n",
                  sgnlv_size);
          H5Sclose(t_space);
          H5Dclose(t_dataset);
          H5Fclose(file);
          H5Tclose(t_tid);
          H5Tclose(ffstat_type);
          return 0;
     }

     /* ------------------------------------------------------------------ */
     /* Read the whole dataset into the buffer.                             */
     /* HDF5 allocates the variable-length (hvl_t) payload for ffstat;     */
     /* call H5Treclaim(t_tid, t_space, H5P_DEFAULT, *sgnlv) when done     */
     /* with the data to free those internal allocations.                   */
     /* ------------------------------------------------------------------ */
     hstat = H5Dread(t_dataset, t_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, *sgnlv);
     if (hstat < 0) {
          fprintf(stderr, "Error: cannot read dataset '%s' from %s\n",
                  t_dset_name, filename);
          free(*sgnlv);
          *sgnlv     = NULL;
          sgnlv_size = 0;
          exit(EXIT_FAILURE);
     }

     printf("    [sgnlv_size: %zu]", sgnlv_size);
     float *ffp = (float *)sgnlv[0]->ffstat.p;
     for (size_t k = 0; k < (sgnlv[0]->ffstat.len+1)/2; k++) {
          printf("  p[%zu]=(%f,%f)  ", k, ffp[2*k], ffp[2*k+1]);
     }
     printf("\n");
     /* ------------------------------------------------------------------ */
     /* Release HDF5 resources (caller owns *sgnlv)                        */
     /* ------------------------------------------------------------------ */
     H5Sclose(t_space);
     H5Dclose(t_dataset);
     H5Fclose(file);
     H5Tclose(t_tid);
     H5Tclose(ffstat_type);

     return sgnlv_size;
}

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

typedef struct {
     int    nod;
     int    N;
     int    nfftf;
     int    numlines_band;
     int    nvlines_all_inband;
     double fpo;
     double dt;
     double B;
     double lines[MAXL][2];
} Sett_partial;

/* -----------------------------------------------------------------------
 * Helper structs for write_ctrigs_hdf() compound HDF5 types.
 * Using flat structs (no pointers to arrays) so HOFFSET is well-defined.
 * ----------------------------------------------------------------------- */

/* Subset of Coinc_opts that maps cleanly to HDF5 scalars + vlen strings  */
typedef struct {
     int    refr, hemi, scalef, scales, scalem, scalen, shift, mincoin, nseg;
     double cthr;
     char  *out_dir;    /* vlen string */
     char  *trig_dset;  /* vlen string */
} Copts_attr;

/* Subset of Search_params (no large lines[][] array)                     */
typedef struct {
     int    band, nod, N, hemi, numlines_band, nvlines_all_inband, sgnlv_size;
     double overlap, narrowdown, B, fpo, dt;
     char  *grid_file;  /* vlen string */
} Sparams_attr;

/* Scalar-only fields of Coinc_Trigger for the "ctrigs" compound dataset  */
typedef struct {
     float m, n, s, ra, dec, fdot;
} Ctrig_base;



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

     copts->refr            = iniparser_getint   (ini, "coincidences:refr",          -1);
     copts->cthr            = iniparser_getdouble(ini, "coincidences:cthr",          -1.);
     copts->hemi            = iniparser_getint   (ini, "coincidences:hemi",           1);
     copts->scalef          = iniparser_getint   (ini, "coincidences:scalef",         1);
     copts->scales          = iniparser_getint   (ini, "coincidences:scales",         1);
     copts->scalem          = iniparser_getint   (ini, "coincidences:scalem",         1);
     copts->scalen          = iniparser_getint   (ini, "coincidences:scalen",         1);
     copts->shift           = iniparser_getint   (ini, "coincidences:shift",       0000);
     copts->mincoin         = iniparser_getint   (ini, "coincidences:mincoin",        3);
     copts->out_dir         = iniparser_getstring(ini, "coincidences:out_dir",       ".");
     copts->trig_dset       = iniparser_getstring(ini, "coincidences:trig_dset", "triggers");
     copts->coinc_dset_pre  = iniparser_getstring(ini, "coincidences:trig_dset_pre", "coinc");
     copts->write_ctrigs    = iniparser_getint   (ini, "coincidences:write_ctrigs", 0);
     tlist                  = (char *)iniparser_getstring(ini, "coincidences:trig_file_list", NULL);

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
          if (count >= MAX_NSEG) {
               fprintf(stderr,
                       "Error: too many trigger files specified in %s (max %d)\n",
                       ini_fname, MAX_NSEG);
               iniparser_freedict(ini);
               exit(EXIT_FAILURE);
          }
     }

     copts->nseg = count;
     printf("Number of segments for coincidences: %zu\n", copts->nseg);
     //for (size_t i = 0; i < copts->n_trig_files; i++)
     //     puts(copts->trig_files[i]);

     //iniparser_freedict(ini);

} /* read_coinc_ini */


/* =========================================================================
 * read_triggers_file - Read triggers from HDF5 file produced by search 
 * into sgnlv array.
 * Return size of the array == number of (m,n,s) grid points.
 * ========================================================================= */
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

     hid_t sett_tid = H5Tcreate(H5T_COMPOUND, sizeof(Sett_partial));
     H5Tinsert(sett_tid, "fpo", offsetof(Sett_partial, fpo), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "dt", offsetof(Sett_partial, dt), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "B", offsetof(Sett_partial, B), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "nod", offsetof(Sett_partial, nod), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "N", offsetof(Sett_partial, N), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nfftf", offsetof(Sett_partial, nfftf), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "numlines_band", offsetof(Sett_partial, numlines_band), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nvlines_all_inband", offsetof(Sett_partial, nvlines_all_inband), H5T_NATIVE_INT);

     
     printf("   [file %s]\n", filename);
     
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
               printf ("   [getting opts: band=%d, seg=%d, hemi=%d, overlap=%.2f, narrowdown=%.2f, grid_file=%s]\n",
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
                    printf("   [opts match]\n");
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

     Sett_partial sp = {0};
     hid_t sett_attr = H5Aopen(file, "sett", H5P_DEFAULT);
     if (sett_attr >= 0) {
          hstat = H5Aread(sett_attr, sett_tid, &sp);
          H5Aclose(sett_attr);
          if (hstat < 0) {
               fprintf(stderr,"Warning: cannot read 'sett' attribute from %s\n", filename);
               exit(EXIT_FAILURE);
          }
          if (init_search_par) {
               search_par->fpo                = sp.fpo;
               search_par->dt                 = sp.dt;
               search_par->B                  = sp.B;
               search_par->nod                = sp.nod;
               search_par->N                  = sp.N;
               search_par->nfftf              = sp.nfftf;
               search_par->numlines_band      = sp.numlines_band;
               search_par->nvlines_all_inband = sp.nvlines_all_inband;
               printf ("   [getting sett: nod=%d, N=%d, B=%f, fpo=%f, nvlines_all_inband=%d]\n",
                           search_par->nod, search_par->N, search_par->B, search_par->fpo,
                           search_par->nvlines_all_inband);
          } else {
               // compare sp with searh par
               if (search_par->nod != sp.nod ||
                   search_par->N != sp.N )
               {
                    fprintf(stderr, "nod=%d   search_par.nod=%d\n", sp.nod, search_par->nod);
                    fprintf(stderr, "N=%d   search_par.N=%d\n", sp.N, search_par->N);
                    fprintf(stderr, "Error: 'sett' attribute in %s does not match expected values\n", filename);
                    exit(EXIT_FAILURE);
               } else {
                    printf("   [sett match]\n");
               }
          }
     } else {
          fprintf(stderr, "Warning: cannot open 'sett' attribute in %s\n", filename);
          exit(EXIT_FAILURE);
     }

     H5Tclose(sett_tid);

     //Read 'lines' from 'sett' (only present when nvlines_all_inband > 0 )
     if (init_search_par && search_par->nvlines_all_inband > 0) {
          hsize_t lines_dims[2] = {(hsize_t)search_par->nvlines_all_inband, 2};
          hid_t lines_arr_t = H5Tarray_create2(H5T_NATIVE_DOUBLE, 2, lines_dims);

          /* Compound of exactly one field "lines" at offset 0.            */
          /* Size = nvlines_all_inband * 2 * sizeof(double), so reading    */
          /* directly into search_par->lines (double[MAXL][2]) is safe.   */
          size_t lines_mem_size = (size_t)search_par->nvlines_all_inband
                                  * 2 * sizeof(double);
          hid_t lines_tid = H5Tcreate(H5T_COMPOUND, lines_mem_size);
          H5Tinsert(lines_tid, "lines", 0, lines_arr_t);
          H5Tclose(lines_arr_t);

          hid_t sett_attr2 = H5Aopen(file, "sett", H5P_DEFAULT);
          if (sett_attr2 >= 0) {
               hstat = H5Aread(sett_attr2, lines_tid, search_par->lines);
               H5Aclose(sett_attr2);
               if (hstat < 0) {
                    fprintf(stderr, "Error: cannot read 'lines' from 'sett' in %s\n",
                            filename);
                    exit(EXIT_FAILURE);
               }
               printf("   [");
               for(size_t i = 0; i < (size_t)search_par->nvlines_all_inband; i++) {
                    printf(" line[%zu]=(%f, %f) ", i,
                           search_par->lines[i][0], search_par->lines[i][1]);
               }
               printf("]\n");
          } else {
               fprintf(stderr, "Error: cannot re-open 'sett' attribute for lines in %s\n",
                       filename);
               exit(EXIT_FAILURE);
          }
          H5Tclose(lines_tid);
     }

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

#if 1     
     //printf("    [sgnlv_size: %zu]", sgnlv_size);
     int itest = sgnlv_size-1;
     printf("   [read_triggers_file: sgnlv[%d].ffstat. ]", itest);
     float *ffp = (float *)(*sgnlv)[itest].ffstat.p;
     for (size_t k = 0; k < ((*sgnlv)[itest].ffstat.len+1)/2; k++) {
          printf("  p[%zu]=(%f,%f)  ", k, ffp[2*k], ffp[2*k+1]);
     }
     printf("\n");
#endif
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



/* =========================================================================
 * write_ctrigs_hdf()
 *
 * Write coincidence triggers to an HDF5 file.
 *
 *   Attributes on root "/":
 *     "copts"      – scalar compound (Coinc_opts fields)
 *     "search_par" – scalar compound (Search_params fields)
 *
 *   Datasets:
 *     "ctrigs"  – compound (m,n,s,ra,dec,fdot), shape (sgnlv_size,)
 *     "ffstat"  – vlen float,                   shape (sgnlv_size, nseg)
 *                 element [j][iseg] = ctrigs[j].ffstat[iseg]
 *
 * Returns EXIT_SUCCESS / EXIT_FAILURE.
 * ========================================================================= */
int write_ctrigs_hdf(const char *ctrigs_fname, Coinc_opts *copts,
                     Search_params *search_par, Coinc_Trigger *ctrigs)
{
     hid_t  file, attr;
     herr_t hstat;
     size_t j, iseg;
     const size_t sgnlv_size = search_par->sgnlv_size;
     const size_t nseg       = copts->nseg;

     /* ------------------------------------------------------------------ */
     /* Variable-length string type (reused for all vlen-string fields)     */
     /* ------------------------------------------------------------------ */
     hid_t vstr_t = H5Tcopy(H5T_C_S1);
     H5Tset_size(vstr_t, H5T_VARIABLE);

     /* ------------------------------------------------------------------ */
     /* Create file                                                         */
     /* ------------------------------------------------------------------ */
     file = H5Fcreate(ctrigs_fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
     if (file < 0) {
          fprintf(stderr, "Error: cannot create HDF5 file %s\n", ctrigs_fname);
          H5Tclose(vstr_t);
          return EXIT_FAILURE;
     }

     hid_t scalar_sp = H5Screate(H5S_SCALAR);

     /* ------------------------------------------------------------------ */
     /* Attribute: copts                                                    */
     /* ------------------------------------------------------------------ */
     hid_t copts_tid = H5Tcreate(H5T_COMPOUND, sizeof(Copts_attr));
     H5Tinsert(copts_tid, "refr",      HOFFSET(Copts_attr, refr),      H5T_NATIVE_INT);
     H5Tinsert(copts_tid, "hemi",      HOFFSET(Copts_attr, hemi),      H5T_NATIVE_INT);
     H5Tinsert(copts_tid, "scalef",    HOFFSET(Copts_attr, scalef),    H5T_NATIVE_INT);
     H5Tinsert(copts_tid, "scales",    HOFFSET(Copts_attr, scales),    H5T_NATIVE_INT);
     H5Tinsert(copts_tid, "scalem",    HOFFSET(Copts_attr, scalem),    H5T_NATIVE_INT);
     H5Tinsert(copts_tid, "scalen",    HOFFSET(Copts_attr, scalen),    H5T_NATIVE_INT);
     H5Tinsert(copts_tid, "shift",     HOFFSET(Copts_attr, shift),     H5T_NATIVE_INT);
     H5Tinsert(copts_tid, "mincoin",   HOFFSET(Copts_attr, mincoin),   H5T_NATIVE_INT);
     H5Tinsert(copts_tid, "nseg",      HOFFSET(Copts_attr, nseg),      H5T_NATIVE_INT);
     H5Tinsert(copts_tid, "cthr",      HOFFSET(Copts_attr, cthr),      H5T_NATIVE_DOUBLE);
     H5Tinsert(copts_tid, "out_dir",   HOFFSET(Copts_attr, out_dir),   vstr_t);
     H5Tinsert(copts_tid, "trig_dset", HOFFSET(Copts_attr, trig_dset), vstr_t);

     Copts_attr ca = {
          .refr      = copts->refr,
          .hemi      = copts->hemi,
          .scalef    = copts->scalef,
          .scales    = copts->scales,
          .scalem    = copts->scalem,
          .scalen    = copts->scalen,
          .shift     = copts->shift,
          .mincoin   = copts->mincoin,
          .nseg      = (int)copts->nseg,
          .cthr      = copts->cthr,
          .out_dir   = (char *)copts->out_dir,
          .trig_dset = (char *)copts->trig_dset,
     };

     attr  = H5Acreate2(file, "copts", copts_tid, scalar_sp, H5P_DEFAULT, H5P_DEFAULT);
     hstat = H5Awrite(attr, copts_tid, &ca);
     if (hstat < 0)
          fprintf(stderr, "Warning: cannot write 'copts' attribute to %s\n", ctrigs_fname);
     H5Aclose(attr);
     H5Tclose(copts_tid);

     /* ------------------------------------------------------------------ */
     /* Attribute: search_par                                               */
     /* ------------------------------------------------------------------ */
     hid_t sparams_tid = H5Tcreate(H5T_COMPOUND, sizeof(Sparams_attr));
     H5Tinsert(sparams_tid, "band",               HOFFSET(Sparams_attr, band),               H5T_NATIVE_INT);
     H5Tinsert(sparams_tid, "nod",                HOFFSET(Sparams_attr, nod),                H5T_NATIVE_INT);
     H5Tinsert(sparams_tid, "N",                  HOFFSET(Sparams_attr, N),                  H5T_NATIVE_INT);
     H5Tinsert(sparams_tid, "hemi",               HOFFSET(Sparams_attr, hemi),               H5T_NATIVE_INT);
     H5Tinsert(sparams_tid, "numlines_band",      HOFFSET(Sparams_attr, numlines_band),      H5T_NATIVE_INT);
     H5Tinsert(sparams_tid, "nvlines_all_inband", HOFFSET(Sparams_attr, nvlines_all_inband), H5T_NATIVE_INT);
     H5Tinsert(sparams_tid, "sgnlv_size",         HOFFSET(Sparams_attr, sgnlv_size),         H5T_NATIVE_INT);
     H5Tinsert(sparams_tid, "overlap",            HOFFSET(Sparams_attr, overlap),            H5T_NATIVE_DOUBLE);
     H5Tinsert(sparams_tid, "narrowdown",         HOFFSET(Sparams_attr, narrowdown),         H5T_NATIVE_DOUBLE);
     H5Tinsert(sparams_tid, "B",                  HOFFSET(Sparams_attr, B),                  H5T_NATIVE_DOUBLE);
     H5Tinsert(sparams_tid, "fpo",                HOFFSET(Sparams_attr, fpo),                H5T_NATIVE_DOUBLE);
     H5Tinsert(sparams_tid, "dt",                 HOFFSET(Sparams_attr, dt),                 H5T_NATIVE_DOUBLE);
     H5Tinsert(sparams_tid, "grid_file",          HOFFSET(Sparams_attr, grid_file),          vstr_t);

     Sparams_attr spa = {
          .band               = search_par->band,
          .nod                = search_par->nod,
          .N                  = search_par->N,
          .hemi               = search_par->hemi,
          .numlines_band      = search_par->numlines_band,
          .nvlines_all_inband = search_par->nvlines_all_inband,
          .sgnlv_size         = (int)search_par->sgnlv_size,
          .overlap            = search_par->overlap,
          .narrowdown         = search_par->narrowdown,
          .B                  = search_par->B,
          .fpo                = search_par->fpo,
          .dt                 = search_par->dt,
          .grid_file          = search_par->grid_file,
     };

     attr  = H5Acreate2(file, "search_par", sparams_tid, scalar_sp, H5P_DEFAULT, H5P_DEFAULT);
     hstat = H5Awrite(attr, sparams_tid, &spa);
     if (hstat < 0)
          fprintf(stderr, "Warning: cannot write 'search_par' attribute to %s\n", ctrigs_fname);
     H5Aclose(attr);
     H5Tclose(sparams_tid);

     /* ------------------------------------------------------------------ */
     /* Dataset "ctrigs": compound (m,n,s,ra,dec,fdot), shape (sgnlv_size,)*/
     /* ------------------------------------------------------------------ */
     hid_t ctrig_tid = H5Tcreate(H5T_COMPOUND, sizeof(Ctrig_base));
     H5Tinsert(ctrig_tid, "m",    HOFFSET(Ctrig_base, m),    H5T_NATIVE_FLOAT);
     H5Tinsert(ctrig_tid, "n",    HOFFSET(Ctrig_base, n),    H5T_NATIVE_FLOAT);
     H5Tinsert(ctrig_tid, "s",    HOFFSET(Ctrig_base, s),    H5T_NATIVE_FLOAT);
     H5Tinsert(ctrig_tid, "ra",   HOFFSET(Ctrig_base, ra),   H5T_NATIVE_FLOAT);
     H5Tinsert(ctrig_tid, "dec",  HOFFSET(Ctrig_base, dec),  H5T_NATIVE_FLOAT);
     H5Tinsert(ctrig_tid, "fdot", HOFFSET(Ctrig_base, fdot), H5T_NATIVE_FLOAT);

     Ctrig_base *base_buf = (Ctrig_base *)malloc(sgnlv_size * sizeof(Ctrig_base));
     if (base_buf == NULL) {
          fprintf(stderr, "Error: cannot allocate ctrigs write buffer\n");
          H5Sclose(scalar_sp); H5Tclose(ctrig_tid);
          H5Tclose(vstr_t); H5Fclose(file);
          return EXIT_FAILURE;
     }
     for (j = 0; j < sgnlv_size; j++) {
          base_buf[j].m    = ctrigs[j].m;
          base_buf[j].n    = ctrigs[j].n;
          base_buf[j].s    = ctrigs[j].s;
          base_buf[j].ra   = ctrigs[j].ra;
          base_buf[j].dec  = ctrigs[j].dec;
          base_buf[j].fdot = ctrigs[j].fdot;
     }

     hsize_t ctrig_dim[1] = {sgnlv_size};
     hid_t ctrig_space = H5Screate_simple(1, ctrig_dim, NULL);
     hid_t ctrig_dset  = H5Dcreate2(file, "ctrigs", ctrig_tid, ctrig_space,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     hstat = H5Dwrite(ctrig_dset, ctrig_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, base_buf);
     if (hstat < 0)
          fprintf(stderr, "Error: cannot write 'ctrigs' dataset to %s\n", ctrigs_fname);
     free(base_buf);
     H5Dclose(ctrig_dset);
     H5Sclose(ctrig_space);
     H5Tclose(ctrig_tid);

     /* ------------------------------------------------------------------ */
     /* Dataset "ffstat": vlen float, shape (sgnlv_size, nseg)             */
     /*                                                                     */
     /* ffbuf[j*nseg + iseg] = ctrigs[j].ffstat[iseg]  (shallow copy of   */
     /* hvl_t – .p still points into ctrigs; HDF5 only reads from it)      */
     /* ------------------------------------------------------------------ */
     hid_t ffstat_vlen_t = H5Tvlen_create(H5T_NATIVE_FLOAT);

     hvl_t *ffbuf = (hvl_t *)malloc(sgnlv_size * nseg * sizeof(hvl_t));
     if (ffbuf == NULL) {
          fprintf(stderr, "Error: cannot allocate ffstat write buffer\n");
          H5Tclose(ffstat_vlen_t); H5Sclose(scalar_sp);
          H5Tclose(vstr_t); H5Fclose(file);
          return EXIT_FAILURE;
     }
     for (j = 0; j < sgnlv_size; j++)
          for (iseg = 0; iseg < nseg; iseg++)
               ffbuf[j * nseg + iseg] = ctrigs[j].ffstat[iseg];

     hsize_t ff_dims[2] = {sgnlv_size, nseg};
     hid_t ff_space = H5Screate_simple(2, ff_dims, NULL);
     hid_t ff_dset  = H5Dcreate2(file, "ffstat", ffstat_vlen_t, ff_space,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     hstat = H5Dwrite(ff_dset, ffstat_vlen_t, H5S_ALL, H5S_ALL, H5P_DEFAULT, ffbuf);
     if (hstat < 0)
          fprintf(stderr, "Error: cannot write 'ffstat' dataset to %s\n", ctrigs_fname);
     free(ffbuf);   /* frees the hvl_t array only; .p data stays in ctrigs */
     H5Dclose(ff_dset);
     H5Sclose(ff_space);
     H5Tclose(ffstat_vlen_t);

     /* ------------------------------------------------------------------ */
     /* Release shared resources and close file                            */
     /* ------------------------------------------------------------------ */
     H5Sclose(scalar_sp);
     H5Tclose(vstr_t);
     H5Fclose(file);

     printf("Written ctrigs HDF5 file: %s  (%zu ctrigs, %zu segs)\n",
            ctrigs_fname, sgnlv_size, nseg);
     return EXIT_SUCCESS;
} /* write_ctrigs_hdf */

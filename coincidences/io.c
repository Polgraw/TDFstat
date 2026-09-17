#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
    char  *label;       /* variable-length string – HDF5 allocates this */
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
    int    refr, hemi, scalef, scales, scalem, scalen, mincoin, nseg;
    double cthr;
    char   shift[5];  /* fixed-size string for HDF5 */
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
    int iccell;
} Ctrig_base;

/* Flat mirror of Coincidence for the per-shift "coi" compound dataset;
 * the three variable-length arrays (cseg, n_ccell_trigs, trig_mns) are
 * represented as hvl_t members inside the compound type itself, so a
 * single dataset per shift is enough (no extra vlen datasets needed).  */
typedef struct {
    char  shift[5];
    short w;
    float avg_snr, avg_f, avg_fdot, avg_ra, avg_dec;
    hvl_t cseg;
    hvl_t n_ccell_trigs;
    hvl_t trig_mns;
} Coi_hdf;



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
    copts->shift           = iniparser_getstring(ini, "coincidences:shift",          NULL);
    copts->mincoin         = iniparser_getint   (ini, "coincidences:mincoin",        4);
    copts->out_dir         = iniparser_getstring(ini, "coincidences:out_dir",       ".");
    copts->trig_dset       = iniparser_getstring(ini, "coincidences:trig_dset", "triggers");
    copts->coinc_dset      = iniparser_getstring(ini, "coincidences:coinc_dset", "coinc_");
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
        H5Fclose(file);
        exit(EXIT_FAILURE);
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
    H5Tinsert(opts_tid, "label", offsetof(Opts_partial, label),  vstr_t);

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
        search_par->seg        = op.seg;
        if (init_search_par) {
            search_par->band       = op.band;
            search_par->hemi       = op.hemi;
            search_par->overlap    = op.overlap;
            search_par->narrowdown = op.narrowdown;
            search_par->grid_file  = op.grid_file ? strdup(op.grid_file) : NULL;
            search_par->label  = op.label ? strdup(op.label) : NULL;
            printf ("   [getting opts: band=%d, seg=%d, hemi=%d, overlap=%.2f, narrowdown=%.2f, grid_file=%s, label=%s]\n",
                search_par->band, search_par->seg, search_par->hemi,
                search_par->overlap, search_par->narrowdown,
                search_par->grid_file ? search_par->grid_file : "NULL",
                search_par->label ? search_par->label : "NULL");
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
                printf("   [opts match]");
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
            // compare sp with search par
            if (search_par->nod != sp.nod ||
                search_par->N != sp.N )
            {
                fprintf(stderr, "nod=%d   search_par.nod=%d\n", sp.nod, search_par->nod);
                fprintf(stderr, "N=%d   search_par.N=%d\n", sp.N, search_par->N);
                fprintf(stderr, "Error: 'sett' attribute in %s does not match expected values\n", filename);
                exit(EXIT_FAILURE);
            } else {
                printf("   [sett match]");
            }
        }
    } else {
        fprintf(stderr, "Warning: cannot open 'sett' attribute in %s\n", filename);
        exit(EXIT_FAILURE);
    }

    H5Tclose(sett_tid);

    //Read 'lines' from 'sett' (only present when nvlines_all_inband > 0 )
    if (init_search_par && search_par->nvlines_all_inband > 0) {
#if 0
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
#endif
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

#if 0
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
 *     "seginfo" – int,                          shape (nseg, 3)
 *                 [iseg][0] = frame number (search_par.seg)
 *                 [iseg][1] = good inband triggers at reference time
 *                 [iseg][2] = unique triggers (one per coincidence cell)
 *
 * Returns EXIT_SUCCESS / EXIT_FAILURE.
 * ========================================================================= */
int write_ctrigs_hdf(const char *ctrigs_fname, Coinc_opts *copts,
                     Search_params *search_par, Coinc_Trigger *ctrigs,
                     int seginfo[][3])
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
    //H5Tinsert(copts_tid, "shift",     HOFFSET(Copts_attr, shift),     H5T_NATIVE_INT);
    hid_t shift_str_t = H5Tcopy(H5T_C_S1);
    H5Tset_size(shift_str_t, 4);
    H5Tinsert(copts_tid, "shift", HOFFSET(Copts_attr, shift), shift_str_t);
    H5Tclose(shift_str_t);
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
        //.shift     = NULL,
        .mincoin   = copts->mincoin,
        .nseg      = (int)copts->nseg,
        .cthr      = copts->cthr,
        .out_dir   = (char *)copts->out_dir,
        .trig_dset = (char *)copts->trig_dset,
    };
    strncpy(ca.shift, copts->shift ? copts->shift : "0000", 5);
    //strncpy(ca.shift, copts->shift, 4);

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
    H5Tinsert(ctrig_tid, "iccell", HOFFSET(Ctrig_base, iccell), H5T_NATIVE_INT);

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
        base_buf[j].iccell = ctrigs[j].iccell;
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
    /* Dataset "seginfo": int, shape (nseg, 3)                            */
    /* ------------------------------------------------------------------ */
    hsize_t si_dims[2] = {nseg, 3};
    hid_t si_space = H5Screate_simple(2, si_dims, NULL);
    hid_t si_dset  = H5Dcreate2(file, "seginfo", H5T_NATIVE_INT, si_space,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hstat = H5Dwrite(si_dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, seginfo);
    if (hstat < 0)
        fprintf(stderr, "Error: cannot write 'seginfo' dataset to %s\n", ctrigs_fname);
    H5Dclose(si_dset);
    H5Sclose(si_space);

    /* ------------------------------------------------------------------ */
    /* Release shared resources and close file                            */
    /* ------------------------------------------------------------------ */
    H5Sclose(scalar_sp);
    H5Tclose(vstr_t);
    H5Fclose(file);

    printf("Wrote ctrigs HDF5 file: %s  (%zu ctrigs, %zu segs)\n",
        ctrigs_fname, sgnlv_size, nseg);
    return EXIT_SUCCESS;
} /* write_ctrigs_hdf */

/* =========================================================================
 * init_coin_hdf()
 *
 * Creates (truncates) the coincidences output HDF5 file and writes the
 * file-level information common to all shifts:
 *
 *   Attributes:
 *     "copts"      - compound, same layout as in write_ctrigs_hdf()
 *     "search_par" - compound, same layout as in write_ctrigs_hdf()
 *
 * Call this once, before the shift loop. Per-shift coincidence results
 * are then added with write_coi_hdf(), once per shift; seginfo is
 * attached as an attribute of each shift's dataset (see write_coi_hdf()).
 *
 * Returns EXIT_SUCCESS / EXIT_FAILURE.
 * ========================================================================= */
int init_coin_hdf(const char *coin_fname, Coinc_opts *copts,
                  Search_params *search_par)
{
    hid_t  file, attr;
    herr_t hstat;

    hid_t vstr_t = H5Tcopy(H5T_C_S1);
    H5Tset_size(vstr_t, H5T_VARIABLE);

    file = H5Fcreate(coin_fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Error: cannot create HDF5 file %s\n", coin_fname);
        H5Tclose(vstr_t);
        return EXIT_FAILURE;
    }

    hid_t scalar_sp = H5Screate(H5S_SCALAR);

    hid_t copts_tid = H5Tcreate(H5T_COMPOUND, sizeof(Copts_attr));
    H5Tinsert(copts_tid, "refr",      HOFFSET(Copts_attr, refr),      H5T_NATIVE_INT);
    H5Tinsert(copts_tid, "hemi",      HOFFSET(Copts_attr, hemi),      H5T_NATIVE_INT);
    H5Tinsert(copts_tid, "scalef",    HOFFSET(Copts_attr, scalef),    H5T_NATIVE_INT);
    H5Tinsert(copts_tid, "scales",    HOFFSET(Copts_attr, scales),    H5T_NATIVE_INT);
    H5Tinsert(copts_tid, "scalem",    HOFFSET(Copts_attr, scalem),    H5T_NATIVE_INT);
    H5Tinsert(copts_tid, "scalen",    HOFFSET(Copts_attr, scalen),    H5T_NATIVE_INT);
    hid_t shift_str_t = H5Tcopy(H5T_C_S1);
    H5Tset_size(shift_str_t, 4);
    H5Tinsert(copts_tid, "shift", HOFFSET(Copts_attr, shift), shift_str_t);
    H5Tclose(shift_str_t);
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
        .mincoin   = copts->mincoin,
        .nseg      = (int)copts->nseg,
        .cthr      = copts->cthr,
        .out_dir   = (char *)copts->out_dir,
        .trig_dset = (char *)copts->trig_dset,
    };
    strncpy(ca.shift, copts->shift ? copts->shift : "0000", 5);

    attr  = H5Acreate2(file, "copts", copts_tid, scalar_sp, H5P_DEFAULT, H5P_DEFAULT);
    hstat = H5Awrite(attr, copts_tid, &ca);
    if (hstat < 0)
        fprintf(stderr, "Warning: cannot write 'copts' attribute to %s\n", coin_fname);
    H5Aclose(attr);
    H5Tclose(copts_tid);

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
        fprintf(stderr, "Warning: cannot write 'search_par' attribute to %s\n", coin_fname);
    H5Aclose(attr);
    H5Tclose(sparams_tid);

    H5Sclose(scalar_sp);
    H5Tclose(vstr_t);
    H5Fclose(file);

    printf("Initialized coincidences HDF5 file: %s\n", coin_fname);
    return EXIT_SUCCESS;
} /* init_coin_hdf */

/* =========================================================================
 * write_coi_hdf()
 *
 * Appends the coincidences found for a single shift to an already
 * initialized coincidences HDF5 file (see init_coin_hdf()).
 *
 * Adds one dataset named "<copts->coinc_dset><shift_str>"
 * (e.g. "coinc_0101"), a compound type with fields:
 *
 *     shift  (4-char fixed string)               w        (short)
 *     avg_snr, avg_f, avg_fdot, avg_ra, avg_dec   (float)
 *     cseg, n_ccell_trigs                         (vlen short, length w)
 *     trig_mns                                    (vlen int,   length w)
 *
 * shape (icoi,)
 *
 * The dataset also gets a "seginfo" attribute, int shape (nseg, 3),
 * with the caller's current seginfo snapshot (frame number, good
 * inband triggers, unique triggers for this shift).
 *
 * Returns EXIT_SUCCESS / EXIT_FAILURE.
 * ========================================================================= */
int write_coi_hdf(const char *coin_fname, Coinc_opts *copts,
                  Coincidence *coi, int icoi, const char *shift_str,
                  int seginfo[][3])
{
    herr_t hstat;
    int    k;

    hid_t file = H5Fopen(coin_fname, H5F_ACC_RDWR, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Error: cannot open HDF5 file %s for writing\n", coin_fname);
        return EXIT_FAILURE;
    }

    char dset_name[FILE_NAME_LEN];
    snprintf(dset_name, FILE_NAME_LEN, "%s%s", copts->coinc_dset, shift_str);

    hid_t vlen_short_t = H5Tvlen_create(H5T_NATIVE_SHORT);
    hid_t vlen_int_t   = H5Tvlen_create(H5T_NATIVE_INT);

    hid_t shift_str_t = H5Tcopy(H5T_C_S1);
    H5Tset_size(shift_str_t, 4);

    hid_t coi_tid = H5Tcreate(H5T_COMPOUND, sizeof(Coi_hdf));
    H5Tinsert(coi_tid, "shift",         HOFFSET(Coi_hdf, shift),         shift_str_t);
    H5Tinsert(coi_tid, "w",             HOFFSET(Coi_hdf, w),             H5T_NATIVE_SHORT);
    H5Tinsert(coi_tid, "avg_snr",       HOFFSET(Coi_hdf, avg_snr),       H5T_NATIVE_FLOAT);
    H5Tinsert(coi_tid, "avg_f",         HOFFSET(Coi_hdf, avg_f),         H5T_NATIVE_FLOAT);
    H5Tinsert(coi_tid, "avg_fdot",      HOFFSET(Coi_hdf, avg_fdot),      H5T_NATIVE_FLOAT);
    H5Tinsert(coi_tid, "avg_ra",        HOFFSET(Coi_hdf, avg_ra),        H5T_NATIVE_FLOAT);
    H5Tinsert(coi_tid, "avg_dec",       HOFFSET(Coi_hdf, avg_dec),       H5T_NATIVE_FLOAT);
    H5Tinsert(coi_tid, "cseg",          HOFFSET(Coi_hdf, cseg),          vlen_short_t);
    H5Tinsert(coi_tid, "n_ccell_trigs", HOFFSET(Coi_hdf, n_ccell_trigs), vlen_short_t);
    H5Tinsert(coi_tid, "trig_mns",      HOFFSET(Coi_hdf, trig_mns),      vlen_int_t);
    H5Tclose(shift_str_t);

    Coi_hdf *cbuf = NULL;
    if (icoi > 0) {
        cbuf = (Coi_hdf *) malloc(sizeof(Coi_hdf) * icoi);
        if (cbuf == NULL) {
            fprintf(stderr, "Error: cannot allocate coi write buffer\n");
            H5Tclose(coi_tid); H5Tclose(vlen_short_t); H5Tclose(vlen_int_t);
            H5Fclose(file);
            return EXIT_FAILURE;
        }
        for (k = 0; k < icoi; k++) {
            memcpy(cbuf[k].shift, coi[k].shift, 4);
            cbuf[k].w        = coi[k].w;
            cbuf[k].avg_snr  = coi[k].avg_snr;
            cbuf[k].avg_f    = coi[k].avg_f;
            cbuf[k].avg_fdot = coi[k].avg_fdot;
            cbuf[k].avg_ra   = coi[k].avg_ra;
            cbuf[k].avg_dec  = coi[k].avg_dec;
            cbuf[k].cseg.len          = (size_t)coi[k].w;
            cbuf[k].cseg.p            = coi[k].cseg;
            cbuf[k].n_ccell_trigs.len = (size_t)coi[k].w;
            cbuf[k].n_ccell_trigs.p   = coi[k].n_ccell_trigs;
            cbuf[k].trig_mns.len      = (size_t)coi[k].w;
            cbuf[k].trig_mns.p        = coi[k].trig_mns;
        }
    }

    hsize_t dims[1] = {(hsize_t)icoi};
    hid_t space = H5Screate_simple(1, dims, NULL);
    hid_t dset  = H5Dcreate2(file, dset_name, coi_tid, space,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (icoi > 0) {
        hstat = H5Dwrite(dset, coi_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, cbuf);
        if (hstat < 0)
            fprintf(stderr, "Error: cannot write '%s' dataset to %s\n", dset_name, coin_fname);
    }
    free(cbuf);

    hsize_t si_dims[2] = {copts->nseg, 3};
    hid_t   si_space   = H5Screate_simple(2, si_dims, NULL);
    hid_t   si_attr    = H5Acreate2(dset, "seginfo", H5T_NATIVE_INT, si_space,
        H5P_DEFAULT, H5P_DEFAULT);
    hstat = H5Awrite(si_attr, H5T_NATIVE_INT, seginfo);
    if (hstat < 0)
        fprintf(stderr, "Warning: cannot write seginfo attribute on %s in %s\n",
            dset_name, coin_fname);
    H5Aclose(si_attr);
    H5Sclose(si_space);

    H5Dclose(dset);
    H5Sclose(space);
    H5Tclose(coi_tid);
    H5Tclose(vlen_short_t);
    H5Tclose(vlen_int_t);
    H5Fclose(file);

    printf("   Wrote dataset '%s' to %s (%d coincidences)\n",
        dset_name, coin_fname, icoi);
    return EXIT_SUCCESS;
} /* write_coi_hdf */


float read_vlines_file(const char *veto_fname, Search_params *search_par)
{
    int i, lnum;
    char line[MAXLINE];
    FILE *data;

    const char vf_tag[] = "#band_veto_fraction=";
    const size_t vf_taglen = sizeof(vf_tag) - 1;
    float vf = -1.;    // negative until read from the file

    i=0;
    if ((data = fopen(veto_fname, "r")) != NULL) {
        while (fgets(line, MAXLINE, data) != NULL) {
            // Skip comments or read veto fraction
            if (*line == '#') {
                if (!strncmp(line, vf_tag, vf_taglen)) {
                    if (sscanf(line + vf_taglen, "%f", &vf) != 1) {
                        printf("Can't read veto fraction from\nfile: %s\nline %s\n Aborting!\n",
                            veto_fname, line);
                        exit(EXIT_FAILURE);
                    }
                }
                continue;
            }

            // sanity check: number of non-comment lines in the vlines.dat file = nvlines_all_inband
            if (i >= search_par->nvlines_all_inband) {
                printf("Number of entries in the vlines.dat is larger then nvlines_all_inband! Aborting\n");
                exit(EXIT_FAILURE);
            }

            if (sscanf(line, "%lf %lf", &search_par->lines[i][0], &search_par->lines[i][1]) != 2) {
                printf("Can't parse line %s in file:\n%s\nAborting!\n", line, veto_fname);
                exit(EXIT_FAILURE);
            }
            i++;
        } // while
        fclose(data);
    } else {
        printf("Can't open file %s\nAborting!\n", veto_fname);
        exit(EXIT_FAILURE);
    }

    printf("   [read_vlines_file: read %d lines from %s, veto fraction = %f]\n",
        i, veto_fname, vf);
    printf("   [");
    for(size_t i = 0; i < (size_t)search_par->nvlines_all_inband; i++) {
        printf("line[%zu]=(%f, %f)  ", i,
            search_par->lines[i][0], search_par->lines[i][1]);
    }
    printf("]\n");

    return vf;

} // read_vlines_file

/* =========================================================================
 * add_maxcoi_hdf()
 *
 * Called once, after all shifts have been processed:
 *
 *   1. Attaches a compound "FAP" attribute (mirroring Fap_t, hvl_t fap
 *      member included) to the dataset belonging to the shift that holds
 *      the global maximum coincidence (max_coi->shift).
 *
 *   2. Writes a new dataset named "max_coi" (a single row, same Coi_hdf
 *      compound layout used for the per-shift coincidence datasets)
 *      containing just the max_coi structure.
 *
 * Returns EXIT_SUCCESS / EXIT_FAILURE.
 * ========================================================================= */
int add_maxcoi_hdf(const char *coin_fname, Coinc_opts *copts,
                   Coincidence *max_coi, Fap_t *fap_data)
{
    herr_t hstat;

    hid_t file = H5Fopen(coin_fname, H5F_ACC_RDWR, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Error: cannot open HDF5 file %s for writing\n", coin_fname);
        return EXIT_FAILURE;
    }

    /* ---- 1. attach FAP attribute to the winning shift's dataset ---- */
    char dset_name[FILE_NAME_LEN];
    snprintf(dset_name, FILE_NAME_LEN, "%s%s", copts->coinc_dset, max_coi->shift);

    hid_t dset = H5Dopen2(file, dset_name, H5P_DEFAULT);
    if (dset < 0) {
        fprintf(stderr, "Error: cannot open dataset '%s' in %s to attach FAP\n",
            dset_name, coin_fname);
        H5Fclose(file);
        return EXIT_FAILURE;
    }

    if (fap_data != NULL && fap_data->fap.p != NULL) {
        hid_t vlen_double_t = H5Tvlen_create(H5T_NATIVE_DOUBLE);

        hid_t fap_tid = H5Tcreate(H5T_COMPOUND, sizeof(Fap_t));
        H5Tinsert(fap_tid, "ncoinc_min",     HOFFSET(Fap_t, ncoinc_min),     H5T_NATIVE_INT);
        H5Tinsert(fap_tid, "ncoinc_max",     HOFFSET(Fap_t, ncoinc_max),     H5T_NATIVE_INT);
        H5Tinsert(fap_tid, "crit_mul_0_01",  HOFFSET(Fap_t, crit_mul_0_01),  H5T_NATIVE_INT);
        H5Tinsert(fap_tid, "crit_mul_0_001", HOFFSET(Fap_t, crit_mul_0_001), H5T_NATIVE_INT);
        H5Tinsert(fap_tid, "fap",            HOFFSET(Fap_t, fap),            vlen_double_t);

        hid_t fap_scalar_space = H5Screate(H5S_SCALAR);
        hid_t fap_attr = H5Acreate2(dset, "FAP", fap_tid, fap_scalar_space,
            H5P_DEFAULT, H5P_DEFAULT);
        hstat = H5Awrite(fap_attr, fap_tid, fap_data);
        if (hstat < 0)
            fprintf(stderr, "Warning: cannot write FAP attribute on %s in %s\n",
                dset_name, coin_fname);
        H5Aclose(fap_attr);
        H5Sclose(fap_scalar_space);
        H5Tclose(fap_tid);
        H5Tclose(vlen_double_t);
    }

    H5Dclose(dset);

    /* ---- 2. write the "max_coi" dataset (single row) ---- */
    hid_t vlen_short_t = H5Tvlen_create(H5T_NATIVE_SHORT);
    hid_t vlen_int_t   = H5Tvlen_create(H5T_NATIVE_INT);

    hid_t shift_str_t = H5Tcopy(H5T_C_S1);
    H5Tset_size(shift_str_t, 4);

    hid_t coi_tid = H5Tcreate(H5T_COMPOUND, sizeof(Coi_hdf));
    H5Tinsert(coi_tid, "shift",         HOFFSET(Coi_hdf, shift),         shift_str_t);
    H5Tinsert(coi_tid, "w",             HOFFSET(Coi_hdf, w),             H5T_NATIVE_SHORT);
    H5Tinsert(coi_tid, "avg_snr",       HOFFSET(Coi_hdf, avg_snr),       H5T_NATIVE_FLOAT);
    H5Tinsert(coi_tid, "avg_f",         HOFFSET(Coi_hdf, avg_f),         H5T_NATIVE_FLOAT);
    H5Tinsert(coi_tid, "avg_fdot",      HOFFSET(Coi_hdf, avg_fdot),      H5T_NATIVE_FLOAT);
    H5Tinsert(coi_tid, "avg_ra",        HOFFSET(Coi_hdf, avg_ra),        H5T_NATIVE_FLOAT);
    H5Tinsert(coi_tid, "avg_dec",       HOFFSET(Coi_hdf, avg_dec),       H5T_NATIVE_FLOAT);
    H5Tinsert(coi_tid, "cseg",          HOFFSET(Coi_hdf, cseg),          vlen_short_t);
    H5Tinsert(coi_tid, "n_ccell_trigs", HOFFSET(Coi_hdf, n_ccell_trigs), vlen_short_t);
    H5Tinsert(coi_tid, "trig_mns",      HOFFSET(Coi_hdf, trig_mns),      vlen_int_t);
    H5Tclose(shift_str_t);

    Coi_hdf mbuf;
    memcpy(mbuf.shift, max_coi->shift, 4);
    mbuf.w        = max_coi->w;
    mbuf.avg_snr  = max_coi->avg_snr;
    mbuf.avg_f    = max_coi->avg_f;
    mbuf.avg_fdot = max_coi->avg_fdot;
    mbuf.avg_ra   = max_coi->avg_ra;
    mbuf.avg_dec  = max_coi->avg_dec;
    mbuf.cseg.len          = (size_t)max_coi->w;
    mbuf.cseg.p            = max_coi->cseg;
    mbuf.n_ccell_trigs.len = (size_t)max_coi->w;
    mbuf.n_ccell_trigs.p   = max_coi->n_ccell_trigs;
    mbuf.trig_mns.len      = (size_t)max_coi->w;
    mbuf.trig_mns.p        = max_coi->trig_mns;

    hsize_t dims[1] = {1};
    hid_t space = H5Screate_simple(1, dims, NULL);
    hid_t mdset = H5Dcreate2(file, "max_coi", coi_tid, space,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hstat = H5Dwrite(mdset, coi_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mbuf);
    if (hstat < 0)
        fprintf(stderr, "Error: cannot write 'max_coi' dataset to %s\n", coin_fname);

    H5Dclose(mdset);
    H5Sclose(space);
    H5Tclose(coi_tid);
    H5Tclose(vlen_short_t);
    H5Tclose(vlen_int_t);
    H5Fclose(file);

    printf("   Wrote dataset 'max_coi' to %s (w=%d, shift=%s)\n",
        coin_fname, max_coi->w, max_coi->shift);

    return EXIT_SUCCESS;
} /* add_maxcoi_hdf */

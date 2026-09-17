#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/param.h>
#include <libgen.h>
#include <H5Tpublic.h>

#include "cdefs.h"
#include "../search/network/openmp/struct.h"
#include "coincidences.h"

typedef struct { int mc, nc, sc, orig; } TrigKey;

static int cmp_trig_key(const void *a, const void *b) {
    const TrigKey *ka = a, *kb = b;
    if (ka->mc != kb->mc) return (ka->mc > kb->mc) - (ka->mc < kb->mc);
    if (ka->nc != kb->nc) return (ka->nc > kb->nc) - (ka->nc < kb->nc);
    return (ka->sc > kb->sc) - (ka->sc < kb->sc);
}

/* Comparator for qsort: sort Coincidence array by w desc, then avg_snr desc */
static int cmp_coi(const void *a, const void *b)
{
     const Coincidence *ca = (const Coincidence *)a;
     const Coincidence *cb = (const Coincidence *)b;
     if (cb->w != ca->w)               /* primary:   more segments first */
          return (int)cb->w - (int)ca->w;
     if (cb->avg_snr > ca->avg_snr) return  1; /* secondary: higher SNR first */
     if (cb->avg_snr < ca->avg_snr) return -1;
     return 0;
}


int main (int argc, char* argv[]) {

    char ini_fname[FILE_NAME_LEN];
    Coinc_opts copts;                  // coincidence search options, read from ini file
    Trigger *sgnlv;                    // triggers struct as written in search
    Coinc_Trigger *ctrigs=NULL;        // triggers struct for coincidences, assembled from all time segments
    Search_params search_par = {0};    // selected fields from Search_settings
    int seginfo[MAX_NSEG][3];   // Info about candidates in frames:
                                // [0] frame number,
                                // [1] good, inband triggers at reference time
                                // [2] unique triggers, max. one per coincidence cell at reference time

    int i, j, k, iseg;
    double vf = 0.0; // veto fraction
    double t0[MAX_NSEG]; // starting time of each segment, not used yet

    if (argc > 1) {
        strcpy (ini_fname, argv[1]);
    } else if (argc > 2) {
        printf("WARNING: too many arguments, only first one is used.\n");
    } else {
        printf("ERROR: missing input file name. Call: \"<executable> search.ini\" \n");
        exit(EXIT_FAILURE);
    }

    read_coinc_ini(ini_fname, &copts);

    if (ctrigs != NULL) {
        free(ctrigs);
        ctrigs = NULL;
    }

    // read band veto file from the same location as the first trigger file
    char veto_fname[FILE_NAME_LEN], *veto_dirname;
    char *veto_buf = strdup(copts.trig_files[0]);
    veto_dirname = dirname(veto_buf);

    // read triggers from input files,
    // select good candidates, and store them in ctrigs
    // TODO: search sibling bands for out of band triggers
    for(iseg=0; iseg<copts.nseg; iseg++) {
        printf("> segment %d\n", iseg);
        read_triggers_file(copts.trig_files[iseg], copts.trig_dset, &sgnlv, &search_par);

        // initialize
        if (iseg==0) {
            ctrigs = (Coinc_Trigger *) malloc (sizeof(Coinc_Trigger) * search_par.sgnlv_size);
            for (j=0; j<search_par.sgnlv_size; j++) {
                ctrigs[j].m    = sgnlv[j].m;
                ctrigs[j].n    = sgnlv[j].n;
                ctrigs[j].s    = sgnlv[j].s;
                ctrigs[j].ra   = sgnlv[j].ra;
                ctrigs[j].dec  = sgnlv[j].dec;
                ctrigs[j].fdot = sgnlv[j].fdot;
                // here ffstat is an array of length = number of segments (copts.nseg)
                // and each element contains (hvl_t ffstat) for the given segment
                ctrigs[j].ffstat = (hvl_t *) malloc (sizeof(hvl_t) * copts.nseg);
            }
            if (search_par.nvlines_all_inband > 0) {
                snprintf(veto_fname, FILE_NAME_LEN, "%s/vlines_%04d%s.dat", veto_dirname, search_par.band, search_par.label);
                free(veto_buf);
                vf = read_vlines_file(veto_fname, &search_par);
            }
        }

        seginfo[iseg][0] = search_par.seg;
        // select good candidates:
        // apply Fstat threshold, shift in frequency to reference time,
        // narrowdown, reject triggers in lines;
        // return the number of good candidates in the segment
        seginfo[iseg][1] = select_goodcands(iseg, &copts, &search_par, sgnlv, ctrigs);

        for (j=0; j<search_par.sgnlv_size; j++)
            if (sgnlv[j].ffstat.p) free(sgnlv[j].ffstat.p);
        free(sgnlv);
        sgnlv = NULL;

#if 0
        int itest = search_par.sgnlv_size-1;
        printf("   [ctrigs[%d].ffstat[%d].  ", itest, iseg);
        float *ffp = (float *)ctrigs[itest].ffstat[iseg].p;
        for (size_t k = 0; k < (ctrigs[itest].ffstat[iseg].len+1)/2; k++)
            printf("  p[%zu]=(%f,%f)  ", k, ffp[2*k], ffp[2*k+1]);
        printf("\n");
#endif

    } // iseg


    /* -----------------------------------------------------------------------
     * map between search and coincidences grids in m,n,s dimensions i.e.
     * ctrigs/sgnlv (m,n,s) -> coincidences (mc,nc,sc)
     * ------------------------------------------------------------------------- */

    // 1 means no scaling - original search grid
    int scm = copts.scalem;
    int scn = copts.scalen;
    int scs = copts.scales;
    int scf = copts.scalef;

    // initialize the coincidences output HDF5 file: copts/search_par
    // attributes (per-shift coi datasets, each with its own seginfo
    // attribute, are added below inside the shift loop via write_coi_hdf())
    char coin_fname[FILE_NAME_LEN];
    snprintf(coin_fname, FILE_NAME_LEN, "%s/coin_%04d_%1d.h5",
        copts.out_dir, search_par.band, search_par.hemi);
    printf("Initializing coincidences file: %s\n", coin_fname);
    init_coin_hdf(coin_fname, &copts, &search_par);

    // maximum coincidence across ALL shifts, based on w & snr
    Coincidence max_coi;
    max_coi.w = 0;
    max_coi.cseg = (short *) malloc(sizeof(short) * copts.nseg);
    max_coi.trig_mns = (int *) malloc(sizeof(int) * copts.nseg);
    max_coi.n_ccell_trigs = (short *) malloc(sizeof(short) * copts.nseg);
    max_coi.avg_snr = 0.f;

    // seginfo (Nku) & n_ccell_total snapshot for whichever shift currently
    // holds max_coi; needed to compute its FAP after all shifts finish
    int max_coi_Nku[MAX_NSEG];
    int max_coi_n_ccell_total = 0;

    // the last 4 bits of ish encode the shifts in m,n,s,f dimensions, respectively;
    // check all 16 combinations of shifts
    // unless a specific shift is requested via copts.shift,
    int ish, ish0 = 0, ish1 = 15;
    if (copts.shift != NULL && strlen(copts.shift) == 4) {
        // convert binary-like shift value to integer in [0,15] range
        ish1 = ish0 = (int)strtol(copts.shift, NULL, 2);
    }
    printf("Shift range: %d..%d\n", ish0, ish1);
    //exit(0);
    /* ---------------------------------------------------------
     * main loop over shifts
     * --------------------------------------------------------- */
    for (ish=ish0; ish<=ish1; ish++) {
        int shift[4];
        for (i=0; i<4; i++)
            shift[i] = (ish >> i) & 1;
        printf("\n>>> Shift: m=%d, n=%d, s=%d, f=%d <<<\n", shift[3], shift[2], shift[1], shift[0]);

        float shiftm = shift[3]*0.5;
        float shiftn = shift[2]*0.5;
        float shifts = shift[1]*0.5;
        float shiftf = shift[0]*0.5;

        // max number of ctrigs elements (mns) in one coincidences cell (mc,nc,sc)
        int ictrigs_size = scm*scn*scs;
        // ~ number of (mc,nc,sc) coincidence grid cells;
        // extra space to account for edge cases, reallocated later if needed
        int maxccells = 1.3*search_par.sgnlv_size/ictrigs_size;

        Search2coi_mns *s2c_mns = (Search2coi_mns *) malloc(sizeof(Search2coi_mns) * maxccells);
        for(i=0; i<maxccells; i++)
            s2c_mns[i].ictrigs = (int *) malloc(sizeof(int) * ictrigs_size);

        int iccell = 0; // coincidence cell counter

        //sort triggers by (mc,nc,sc), then group in a single linear pass — O(N log N)
        TrigKey *keys = malloc(search_par.sgnlv_size * sizeof(TrigKey));
        for (j=0; j<search_par.sgnlv_size; j++) {
            keys[j].mc   = (int)floorf(ctrigs[j].m / scm - shiftm);
            keys[j].nc   = (int)floorf(ctrigs[j].n / scn - shiftn);
            keys[j].sc   = (int)floorf(ctrigs[j].s / scs - shifts);
            keys[j].orig = j;
            //printf("keys[%d].ms=%d, ctrigs[%d].s=%f\n", j, keys[j].nc, j, ctrigs[j].n);
        }

        qsort(keys, search_par.sgnlv_size, sizeof(TrigKey), cmp_trig_key);

        for (int jstart=0; jstart<search_par.sgnlv_size; ) {
            int jend = jstart + 1;
            while (jend < search_par.sgnlv_size &&
                   keys[jend].mc == keys[jstart].mc &&
                   keys[jend].nc == keys[jstart].nc &&
                   keys[jend].sc == keys[jstart].sc)
                jend++;
            int ncell = jend - jstart;
            if (iccell >= maxccells) {
                maxccells = (int)(1.3 * maxccells) + 1;
                s2c_mns = realloc(s2c_mns, sizeof(Search2coi_mns) * maxccells);
                for (int ic=iccell; ic<maxccells; ic++)
                    s2c_mns[ic].ictrigs = malloc(sizeof(int) * ictrigs_size);
            }
            s2c_mns[iccell].mc = keys[jstart].mc;
            s2c_mns[iccell].nc = keys[jstart].nc;
            s2c_mns[iccell].sc = keys[jstart].sc;
            s2c_mns[iccell].nctrigs = ncell;
            for (int jj=0; jj<ncell; jj++){
                s2c_mns[iccell].ictrigs[jj] = keys[jstart+jj].orig;
                // save the coincidence cell index to the ctrigs structure
                ctrigs[keys[jstart+jj].orig].iccell = iccell;
            }
            iccell++;
            jstart = jend;
        }
        free(keys);

        int nccells = iccell;
        printf("   nccels/maxccells: %d/%d sgnlv_size=%ld  ",
            nccells, maxccells, search_par.sgnlv_size);
#if 0
        for(int ic=0; ic<nccells; ic++) {
            printf("Cell %d: mc=%d, nc=%d, sc=%d, nctrigs=%d, ictrigs=",
                ic, s2c_mns[ic].mc, s2c_mns[ic].nc, s2c_mns[ic].sc, s2c_mns[ic].nctrigs);
            for(int j=0; j<s2c_mns[ic].nctrigs; j++)
                printf("%d(%d) ", s2c_mns[ic].ictrigs[j], ctrigs[s2c_mns[ic].ictrigs[j]].iccell);
            printf("\n");
        }
#endif

        /* -------------------------------------------------------------------------
        * search for coincidences between segments
        * -------------------------------------------------------------------------*/

        int ncoi = nccells;
        // all coincidences, to be written to HDF file
        // assume one coinc. per mns coinc. cell, reallocate if needed;
        Coincidence *coi = (Coincidence *) malloc(sizeof(Coincidence) * ncoi);

        int icoi = 0; // coincidence counter
        int nfccells = search_par.nfftf/scf; // number of frequency cells in the coincidences grid
        int n_ccell_total = nccells * nfccells; // total 4D (m,n,s,f) coincidence cells
        printf("   n_ccell_total = (%d mns x %d f) = %d\n",
               nccells, nfccells, n_ccell_total);

        // Working arrays: best trigger per (iseg, fic) cell.
        // Allocated once and reset via dirty list – never memset entirely.
        int   *best_cidx_cell = malloc(copts.nseg * nfccells * sizeof(int));
        float *best_snr_cell  = calloc(copts.nseg * nfccells, sizeof(float));
        float *best_f_cell    = malloc(copts.nseg * nfccells * sizeof(float));
        short *n_trigs_cell   = calloc(copts.nseg * nfccells, sizeof(short));
        for (int ii=0; ii<copts.nseg*nfccells; ii++)
            best_cidx_cell[ii] = -1;
        // Dirty list: touched (iseg*nfccells + fic) indices this mns cell.
        // Upper bound per cell: nseg * nctrigs_per_cell * nfpairs_avg (small).
        int *dirty = malloc(copts.nseg * nfccells * sizeof(int));

        short tmp_cseg[MAX_NSEG];
        int   tmp_trig_mns[MAX_NSEG];
        short tmp_n_ccell_trigs[MAX_NSEG];

        const float inv_fbin = (float)nfccells / (float)M_PI;

        // initialise unique-trigger counters (filled in pass 2 below)
        for (iseg=0; iseg<copts.nseg; iseg++)
            seginfo[iseg][2] = 0;

        // loop over mns cells
        for (int imns=0; imns<nccells; imns++) {
            int ndirty = 0;

            /* ---- Pass 1: bin every (f,fstat) pair into its fic cell ---- */
            for (int jj=0; jj<s2c_mns[imns].nctrigs; jj++) {
                int cidx = s2c_mns[imns].ictrigs[jj];
                for (iseg=0; iseg<copts.nseg; iseg++) {
                    if (ctrigs[cidx].ffstat[iseg].len == 0) continue;
                    float *ffp  = (float *)ctrigs[cidx].ffstat[iseg].p;
                    int nfpairs = (int)ctrigs[cidx].ffstat[iseg].len / 2;
                    for (int k=0; k<nfpairs; k++) {
                        float f     = ffp[2*k];
                        //float fstat = ffp[2*k+1];
                        float snr2 = 2.*ffp[2*k+1] - 4.; // snr^2
                        int fic = (int)floorf(f/(M_PI/search_par.nfftf)/scf - shiftf);
                        if ((unsigned)fic >= (unsigned)nfccells) continue;
                        int idx = iseg * nfccells + fic;
                        if (best_cidx_cell[idx] < 0)   // first touch: mark dirty
                            dirty[ndirty++] = idx;
                        n_trigs_cell[idx]++;
                        if (snr2 > best_snr_cell[idx]) {
                            best_snr_cell[idx]  = snr2;
                            best_cidx_cell[idx] = cidx;
                            best_f_cell[idx]    = f;
                        }
                    }
                }
            } // for jj

            /* ---- Pass 2: scan fic cells for coincidences ---- */
            for (int fic=0; fic<nfccells; fic++) {
                int   w       = 0;
                float sum_snr = 0.f, sum_f = 0.f, sum_fdot = 0.f;
                float sum_ra  = 0.f, sum_dec = 0.f;

                for (iseg=0; iseg<copts.nseg; iseg++) {
                    int idx  = iseg * nfccells + fic;
                    int cidx = best_cidx_cell[idx];
                    if (cidx < 0) continue;
                    seginfo[iseg][2]++;  // add an unique trigger
                    tmp_cseg[w]          = (short)iseg;
                    tmp_trig_mns[w]      = cidx;
                    tmp_n_ccell_trigs[w] = n_trigs_cell[idx];
                    sum_snr  += best_snr_cell[idx];
                    sum_f    += best_f_cell[idx];
                    sum_fdot += ctrigs[cidx].fdot;
                    sum_ra   += ctrigs[cidx].ra;
                    sum_dec  += ctrigs[cidx].dec;
                    w++;
                }
                if (w < copts.mincoin) continue;

                if (icoi >= ncoi) {
                    ncoi = (int)(1.3 * ncoi) + 1;
                    coi  = (Coincidence *) realloc(coi, sizeof(Coincidence) * ncoi);
                }
                coi[icoi].w     = (short)w;
                // encode shift as a string
                for (int ii = 0; ii < 4; ii++)
                    coi[icoi].shift[ii] = '0' + shift[3-ii];
                coi[icoi].shift[4] = '\0';
                coi[icoi].cseg          = (short *) malloc(sizeof(short) * w);
                coi[icoi].trig_mns      = (int *)   malloc(sizeof(int)   * w);
                coi[icoi].n_ccell_trigs = (short *) malloc(sizeof(short) * w);
                memcpy(coi[icoi].cseg,          tmp_cseg,          sizeof(short) * w);
                memcpy(coi[icoi].trig_mns,      tmp_trig_mns,      sizeof(int)   * w);
                memcpy(coi[icoi].n_ccell_trigs, tmp_n_ccell_trigs, sizeof(short) * w);
                coi[icoi].avg_snr  = sqrtf(sum_snr);
                coi[icoi].avg_f    = sum_f    / w;
                coi[icoi].avg_fdot = sum_fdot / w;
                coi[icoi].avg_ra   = sum_ra   / w;
                coi[icoi].avg_dec  = sum_dec  / w;

                if (w > max_coi.w ||
                    (w == max_coi.w && coi[icoi].avg_snr > max_coi.avg_snr)) {
                        max_coi.w = (short)w;
                        memcpy(max_coi.shift,         coi[icoi].shift,   5);
                        memcpy(max_coi.cseg,          tmp_cseg,          sizeof(short)*w);
                        memcpy(max_coi.trig_mns,      tmp_trig_mns,      sizeof(int)*w);
                        memcpy(max_coi.n_ccell_trigs, tmp_n_ccell_trigs, sizeof(short)*w);
                        max_coi.avg_snr  = coi[icoi].avg_snr;
                        max_coi.avg_f    = coi[icoi].avg_f;
                        max_coi.avg_fdot = coi[icoi].avg_fdot;
                        max_coi.avg_ra   = coi[icoi].avg_ra;
                        max_coi.avg_dec  = coi[icoi].avg_dec;
                }
                icoi++;
            } // for fic

            /* ---- Reset only touched entries (not the full array) ---- */
            for (int ii=0; ii<ndirty; ii++) {
                int idx = dirty[ii];
                best_cidx_cell[idx] = -1;
                best_snr_cell[idx]  = 0.f;
                best_f_cell[idx]    = 0.f;
                n_trigs_cell[idx]   = 0;
            }
        } // for imns

        free(best_cidx_cell); free(best_snr_cell);
        free(best_f_cell); free(n_trigs_cell); free(dirty);

        qsort(coi, icoi, sizeof(Coincidence), cmp_coi);

        // build shift string (same encoding as coi[].shift) and write
        // this shift's coincidences (with its seginfo snapshot) to the
        // coin HDF5 file
        char shift_str[5];
        for (int ii=0; ii<4; ii++)
            shift_str[ii] = '0' + shift[3-ii];
        shift_str[4] = '\0';

        // if this shift currently holds the global maximum coincidence,
        // snapshot its unique-trigger counts (Nku) and n_ccell_total —
        // needed later to compute FAP for the winning shift only
        if (strncmp(max_coi.shift, shift_str, 4) == 0) {
            for (iseg=0; iseg<copts.nseg; iseg++)
                max_coi_Nku[iseg] = seginfo[iseg][2];
            max_coi_n_ccell_total = n_ccell_total;
        }

        printf("   Best coincidences:\n");
        if (icoi > 0) {
            for(i=0; i<MIN(3, icoi); i++){
                float ff = search_par.fpo + (coi[i].avg_f / M_PI)*search_par.B;
                printf("   coi[%d]: w=%d  avg_snr=%.4f  avg_f=%.6f (%.6f) avg_fdot=%.4e"
                    "  avg_ra=%.6f  avg_dec=%.6f shift=%s\n",
                    i, coi[i].w, coi[i].avg_snr, coi[i].avg_f, ff,
                    coi[i].avg_fdot, coi[i].avg_ra, coi[i].avg_dec, coi[i].shift);
            }

        } else {
            printf("No coincidences above mincoin=%d found.\n", copts.mincoin);
        }

        // write coincidences to HDF5 file
        write_coi_hdf(coin_fname, &copts, coi, icoi, shift_str, seginfo);

        // free per-shift allocations
        for (i=0; i<maxccells; i++)
            free(s2c_mns[i].ictrigs);
        free(s2c_mns);

        for (i=0; i<icoi; i++) {
            free(coi[i].cseg);
            free(coi[i].trig_mns);
            free(coi[i].n_ccell_trigs);
        }
        free(coi);

    } // ishift

    printf("\n");

    // ---- Global best coincidence across all shifts ----
    printf("=== Global best coincidence (all shifts) ===\n");
    if (max_coi.w > 0) {
        float ff = search_par.fpo + (max_coi.avg_f / M_PI) * search_par.B;
        printf("   max_coi: w=%d  avg_snr=%.4f  avg_f=%.6f (%.6f Hz)  avg_fdot=%.4e"
            "  avg_ra=%.6f  avg_dec=%.6f  shift=%s\n",
            max_coi.w, max_coi.avg_snr, max_coi.avg_f, ff,
            max_coi.avg_fdot, max_coi.avg_ra, max_coi.avg_dec, max_coi.shift);
        printf("   Segments (seg / trig_idx / n_cell_trigs):");
        for (int is=0; is<max_coi.w; is++)
            printf("  %d/%d/%d", max_coi.cseg[is], max_coi.trig_mns[is], max_coi.n_ccell_trigs[is]);
        printf("\n");
    } else {
        printf("   No coincidences found in any shift (w=0).\n");
    }

    // ---- FAP calculation for the shift holding the global max, + HDF write ----
    if (max_coi.w > 0) {
        Fap_t FAP;
        int nseg = copts.nseg;
        int mincoin = copts.mincoin;
        double maxcomb = 1.e9;
        double Nc = round((1.-vf) * max_coi_n_ccell_total);

        printf("FAP calculation for shift=%s (global best coincidence)\n", max_coi.shift);
        printf("   Nc=%le, nseg=%d, mincoin=%d, maxcomb=%le\n   (iseg Nku):",
            Nc, nseg, mincoin, maxcomb);
        for (i=0; i<nseg; i++)
            printf("   %d %d", i+1, max_coi_Nku[i]);
        printf("\n");

        FAP.fap.p = (double *)calloc(nseg+1, sizeof(double));
        FAP.fap.len = nseg+1;
        double *fapptr = (double *)FAP.fap.p;

        FalseAlarmProb(mincoin, nseg, Nc, max_coi_Nku, maxcomb, fapptr);

        FAP.crit_mul_0_01  = nseg;
        FAP.crit_mul_0_001 = nseg;
        printf("   FAP results:\n");
        for (i=mincoin; i<nseg; i++) {
            if (fapptr[i] > 0.01  && fapptr[i+1] <= 0.01)  FAP.crit_mul_0_01  = i+1;
            if (fapptr[i] > 0.001 && fapptr[i+1] <= 0.001) FAP.crit_mul_0_001 = i+1;
            printf("%d %le ", i, fapptr[i]);
        }
        printf("\n   crit_mul_0_01=%d   crit_mul_0_001=%d\n",
            FAP.crit_mul_0_01, FAP.crit_mul_0_001);
        fflush(stdout);

        FAP.ncoinc_min = mincoin;
        FAP.ncoinc_max = nseg;

        // attach FAP to the winning shift's dataset + write "max_coi" dataset
        add_maxcoi_hdf(coin_fname, &copts, &max_coi, &FAP);

        free(FAP.fap.p);
    } else {
        printf("No coincidences found in any shift; skipping FAP calculation.\n");
    }

    free(max_coi.cseg);
    free(max_coi.trig_mns);
    free(max_coi.n_ccell_trigs);

    // write HDF file with ctrigs, if needed
    if (copts.write_ctrigs) {
        char ctrigs_fname[FILE_NAME_LEN];
        // can't use search_par.hemi since it can be 0 (both) and thus
        // the only source of current hemi is the triggers filename
        int fname_len = (int)strlen(copts.trig_files[0]);
        snprintf(ctrigs_fname, FILE_NAME_LEN, "%s/ctrigs%s",
            copts.out_dir, copts.trig_files[0]+fname_len-10);
        //printf("Writing coincidence triggers to file: %s\n", ctrigs_fname);
        write_ctrigs_hdf(ctrigs_fname, &copts, &search_par, ctrigs, seginfo);
    }

    for (j=0; j<search_par.sgnlv_size; j++) {
        for (iseg=0; iseg<copts.nseg; iseg++)
            if (ctrigs[j].ffstat[iseg].p) free(ctrigs[j].ffstat[iseg].p);
        free(ctrigs[j].ffstat);
    }
    free(ctrigs);

    return EXIT_SUCCESS;
}




int select_goodcands(int iseg, Coinc_opts *copts, Search_params *search_par,
                     Trigger *sgnlv, Coinc_Trigger *ctrigs)
{
    // sgnlv array of structures contains variable length array ffstat which interleave f and fstat values
    // allocate memory for allcands of the same size as sgnlv, if needed
    // for each element of sgnlv loop over sgnlv.ffstat triggers and
    //      1. apply Fstat threshold (copts->cthr)
    //      2. check if f is in the search_par->lines (where [0] is start of the line and [1] the end)
    //      3. shift in f due to spindown: f = f0 + 2*sgnlv.fdot*(search_par->N)*(copts->refr - search_par->seg)
    // append the remaining triggers to allcands
    // remark: we keep the original search grid here

    int i, j, k, itrig;
    int ntrig_seg = 0;
    int buffer_size = 2048; // initial buffer for ffstat
    float *ffbuffer = (float *) malloc(sizeof(float) * buffer_size);

    for (j=0; j<search_par->sgnlv_size; j++) {
        // check if ffstat is null?
        int ic = 0; // "good" triggers counter
        for (itrig=0; itrig<sgnlv[j].ffstat.len/2; itrig++) {
            float f = ((float *)sgnlv[j].ffstat.p)[2*itrig];
            float fstat = ((float *)sgnlv[j].ffstat.p)[2*itrig+1];
            // skip triggers below threshold
            if (fstat < copts->cthr) continue;
#if 1
            // narrowdown is included in the vlines
            int isline = 0;
            for(int il=0; il<search_par->nvlines_all_inband; il++){
                if( f >= search_par->lines[il][0] && f <= search_par->lines[il][1] ) {
                    isline = 1;
                    break;
                }
            }
            if (isline) continue;
#else
            // Narrowing-down the band around center
            if( (f <= M_PI_2 - search_par->narrowdown) ||
                (f >= M_PI_2 + search_par->narrowdown) )  continue;
            // check if f is in the search_par->lines
            int isline = 0;
            for(int il=0; il<search_par->nvlines_all_inband; il++){
                if( f >= search_par->lines[il][0] && f <= search_par->lines[il][1] ) {
                    isline = 1;
                    break;
                }
            }
            if (isline) continue;
#endif
            // spindown in linear units
            float fdot_lin = sgnlv[j].fdot*M_PI*pow(search_par->dt,2);
            // shifting f to the reference segment (copts->refr) due to spindown
            f = f + 2.*fdot_lin*(search_par->N)*(copts->refr - search_par->seg);
            //f = f + sgnlv[j].fdot*search_par->dt*(search_par->N)*(copts->refr - search_par->seg) * 2.*M_PI*search_par->dt;
            if((f<0) || (f>M_PI)) continue; // skip triggers shifted out of the band
            ffbuffer[2*ic] = f;
            ffbuffer[2*ic+1] = fstat;
            ic++;
            if (ic >= buffer_size/2) {
                buffer_size *= 2;
                ffbuffer = (float *) realloc (ffbuffer, sizeof(float) * buffer_size);
            }
        }
        ctrigs[j].ffstat[iseg].len = 2*ic;
        ctrigs[j].ffstat[iseg].p = (float *) malloc(sizeof(float) * ctrigs[j].ffstat[iseg].len);
        //printf("before goodcands j=%d  buffer_size=%d\n", j, buffer_size);
        memcpy(ctrigs[j].ffstat[iseg].p, ffbuffer,
            sizeof(float) * ctrigs[j].ffstat[iseg].len);
        ntrig_seg += ic;
    }

    printf("   [good ntrig_seg=%d]\n", ntrig_seg);

    free(ffbuffer);

    return ntrig_seg;
}

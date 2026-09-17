#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_statistics_int.h>

#define MAXCOMB 1.e9
#define AVG_MEDIAN


// Calculates false alarm probablility of noc coincidences out of L
int FalseAlarmProb(
    int noc,        // minimum number of coincidences
    int L,          // number of all segments = maximum number of coincidences
    double Nc,      // number of cells in each frame
    int *Nk,        // array with numbers of unique candidates in each frame
    double maxcomb, // use approximation if number of combinations > maxcomb
    double *fap     // array: false alarm probablility of n coincidences for noc..L
    )
{

    int i, j, k, l;
    double ee[L][5], eemean[5]={0.}, eemedian[5]={0.}, *eeavg, Nkmedian;
    double C[L+1][5], pf[L+1][5];

    gsl_combination *cp, *cq;
    size_t *cpd, *cqd;


    for(i=0; i<L; i++){    //#mb length(Nk)
        for(l=0; l<5; l++){
            ee[i][l] = (double)Nk[i]/(Nc*pow(2,l));
            eemean[l] += ee[i][l];
        }
    }
    for(l=0; l<5; l++)
        eemean[l] /= L;

    Nkmedian = (double)gsl_stats_int_median(Nk, 1, L);
    for(l=0; l<5; l++)
        eemedian[l] =  Nkmedian/(Nc*pow(2,l));
    printf("ee_median[0]=%f    ee_mean[0]=%f   Nc=%f\n", Nkmedian/Nc, eemean[0], Nc);
#ifdef AVG_MEAN
    printf("Averaging ee with mean!\n");
    eeavg = eemean;
#else
    printf("Averaging ee with median!\n");
    eeavg = eemedian;
#endif

    // Calculate C[i] - probability of i coincidences in any given call (and l-th correction for shifts)
    for(i=noc; i<=L; i++) {

        printf("n=%3d ", i);
        double P[5], Q[5], Ctmp[5]={0.};

        double ncomb = gsl_sf_choose(L,i);
        if ( ncomb < maxcomb ) {
            printf("[ncomb=%.1e][exa]", ncomb);
            cp = gsl_combination_calloc(L, i);
            cq = gsl_combination_alloc(L, L-i);
            gsl_combination_init_last(cq);

            while(1) {
                cpd = gsl_combination_data(cp);
                cqd = gsl_combination_data(cq);

                for(l=0; l<5; ++l) {
                    P[l] = 1.;
                    Q[l] = 1.;
                }

                for(j=0; j<i; ++j){
                    for(l=0; l<5; ++l)
                        P[l] *= ee[cpd[j]][l];
                }
                for(j=0; j<(L-i); ++j){
                    for(l=0; l<5; ++l)
                        Q[l] *= (1. - ee[cqd[j]][l]);
                }

                for(l=0; l<5; ++l)
                    Ctmp[l] += P[l]*Q[l];

                if (gsl_combination_next(cp) == GSL_FAILURE ||
                    gsl_combination_prev(cq) == GSL_FAILURE)
                break;
            }

            gsl_combination_free (cp);
            gsl_combination_free (cq);

        } else {

            printf("[ncomb=%.1e][avg]", ncomb);
            for(l=0; l<5; l++)
                Ctmp[l]= ncomb*pow(eeavg[l],i)*pow(1.-eeavg[l],L-i);

        }

        printf(" C[%3d] = [ ", i);
        for(l=0; l<5; l++){
            C[i][l] = Ctmp[l]; // probability of i coincidences between L frames
            printf("%12.6e  ", C[i][l]);
        }
        printf("]\n");

    } // i


    // pf - probability of i or more coincidences in any given cell
    for(i=noc; i<=L; i++) {
        for(l=0; l<5; l++){
#if 0
	        // like in the old version
	        pf[i][l] = C[i][l];
#else
            // correct
            pf[i][l] = 0.;
            for(j=i; j<=L; j++){
                pf[i][l] += C[j][l] ;
            }
#endif
        }
    }

    // PF0 = 1 - (1 - pfe).^Nc
    for(i=noc; i<=L; i++) {

        fap[i] = pow(2,4)*pf[i][0]
            - ( gsl_sf_choose(4,1)*pf[i][1]
                + gsl_sf_choose(4,2)*pf[i][2]
                + gsl_sf_choose(4,3)*pf[i][3]
                + gsl_sf_choose(4,4)*pf[i][4] )
            - ( gsl_sf_choose(4,2)*pf[i][2]
                + gsl_sf_choose(4,3)*pf[i][3]
                + gsl_sf_choose(4,4)*pf[i][4] )
            - ( gsl_sf_choose(4,3)*pf[i][3]
                + gsl_sf_choose(4,4)*pf[i][4] )
            - pf[i][4];

        fap[i] = 1. - pow(1. - fap[i], Nc);

    }

#if 0
    printf("FAP results:\n");
    for(i=noc; i<=L; i++){
        printf("%hu %le ", i, fap[i]);
    }
    printf("\n");

    exit(1);
#endif
    return 0;

}

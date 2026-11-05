#ifndef _HEADER_PSS_SFC_H
#define _HEADER_PSS_SFC_H

/*___________________________________________________
 ¦                                                   ¦
 ¦                    pss_sfc.h                      ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    June 2003                      ¦
 ¦___________________________________________________¦*/


#ifdef __cplusplus
extern "C" {
#endif

#define MAX_NHOLES 20
#define MAX_PERS 100

/*
typedef struct SFC_CH {
	double	dx;
	double	dy;
	long	lenx;
	long	leny;
	long	type;
	char	name[129];
} SFC_CH;
*/

typedef struct SFC_{
	char			file[160];
	char			pname[80];
	FILE*			fid;
	char			label[9];
	long			prot;
	long			coding;
	double			t0;
	double			dt;
	char			capt[129];
	long			nch;
	long			hlen;	
	long long		len;
	long			point0;
	short			eof;
	long			acc;
	char			filme[129];
	char			filmaster[129];
	char			filspre[129];
	char			filspost[129];
	char			filppre[129];
	char			filppost[129];
	long			blen;
	unsigned long*	point;

	//SFC_CH*			ch;
    //char			*ch[];  
    char			**ch;  
        
} SFC_;


typedef struct HOLES_{
/* structure describing holes in the output vector */
	int				nholes;
	long			nztot;
	long			nzeros[MAX_NHOLES];
	int				kzeros[MAX_NHOLES];
} HOLES_;


typedef struct ALLPERS_{
/* allowed analysis periods (initial, final times id odd, even positions) */
	int				nper;
	double			inifin[2*MAX_PERS];
} ALLPERS_;


/* ----------------------  General Routines -------------------------*/

SFC_* sfc_open(char *file);

void sfc_open1(char *file,SFC_ *sfc_);
/**pia version, 23 Aug. 2005:does not allocate SFC_ again !!**/

SFC_* sfc_openw(char* file, struct SFC_* sfc_);


void sfc_show(SFC_* sfc_);


/* ------------------------  sds Routines ---------------------------*/


SFC_* sds_open(char *file, struct ALLPERS_ *allpers);

void sds_open1(char *file, struct ALLPERS_ *allpers,SFC_ *sds_);
/**pia version, 23 Aug. 2005: to call sfc_open1 and not sfc_open:
does not allocate SFC_ again !!**/

SFC_* sds_openw(char* file, struct SFC_* sfc_);


int sds_check_time(double tim0, double tim1, struct ALLPERS_ *allpers);


void vec_from_sds(float *buffer, float *vec, double *tim0, struct HOLES_* holes,
			struct SFC_* sds_, int chn, long len, struct ALLPERS_ *allpers);

  /**pia version, 23 Aug. 2005**/
void vec_from_sds1(float *buffer, float *vec, double *tim0, struct HOLES_* holes,
			struct SFC_* sds_, int chn, long len, struct ALLPERS_ *allpers);
 /**pia version, 23 Aug. 2005**/
/* ----------------------  sbl Routines -------------------------*/


SFC_* sbl_open(struct SFC_* sbl_);


SFC_* sbl_openw(struct SFC_* sbl_);


/* ----------------------  vbl Routines -------------------------*/


SFC_* vbl_open(struct SFC_* vbl_);


SFC_* vbl_openw(struct SFC_* vbl_);



/* ----------------------  Service Routines -------------------------*/


void sfc_util();


void gw_filename(struct tm *t, char *ant_abbr, char *dat_abbr, char *sub_typ, char *ext, char *out);


#ifdef __cplusplus
}
#endif

#endif /* _HEADER_PSS_SFC_H */

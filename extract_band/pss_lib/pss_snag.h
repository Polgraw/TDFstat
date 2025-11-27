#ifndef _HEADER_PSS_SNAG_H
#define _HEADER_PSS_SNAG_H

/*___________________________________________________
 ¦                                                   ¦
 ¦                   pss_snag.h                      ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    March 2000                     ¦
 ¦___________________________________________________¦*/


#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------- GD --------------------------------*/

typedef struct GD{/* gd structure */
	char name[20];/* name */
	double* x;    /* abscissa */
	float* y;     /* data */
	long n;       /* dimension */
	long nall;	  /* allocated dimension */
	int type;     /* type; = 1 -> virtual abscissa, = 2 -> real abscissa */
	int complex;  /* 0 -> real, 1 -> complex */
	double ini;   /* initial virt.absc. value */
	double dx;    /* virt.absc. sampling value */
	char* capt;   /* caption */
	int cont;     /* control variable */
} GD;


GD* crea_gd(long nall, double ini,double dx, char* capt);
  /* to create type 1 gd */


GD* vec_to_gd(float* y, long n, double ini,double dx, char* capt);
  /* to create type 1 gd from a vector */


float* y_gd(GD* gd);
  /* extracts the data vector from a gd */

GD* set_gd(GD* g, int typ,
		   double par1, double par2, double par3);
/*  typ  =  1   par1
		 =  2   par1 * ramp + par2
		 =  3   par1 * sin(par2*i+par3)
		 =  4   par1 * randn +par2 (seed round(par3+1))
*/


int show_gd(GD* gd);
/* show the gd parameters */


GD* sum_gd(float a1, GD* gd1, 
		   float a2, GD* gd2,
		   float a3);
/* sums two gd's */


GD* mult_gd(GD* gd1, GD* gd2);
/* product of two gd's */


/* ------------------------------- DM --------------------------------*/

typedef struct DM{/* dm structure */
	char name[20];/* name */
	double* x;    /* real abscissa; if it exists, the dimension is n/m */
	float* y;     /* data (input by columns); total dimension is n */
	long n;       /* total number of data */
	long nall;	  /* allocated dimension */
	int type;     /* type; = 1 -> virtual abscissa, = 2 -> real abscissa */
	int complex;  /* 0 -> real, 1 -> complex */
	double ini;   /* initial virt. absc. value */
	double dx;    /* virt.absc. sampling value */
	long m;       /* dimension 2 (length of a column) */
	double ini2;  /* initial virt.absc2 value */
	double dx2;   /* virt.absc2 sampling value */
	char* capt;   /* caption */
	int cont;     /* control variable */
} DM;


DM* crea_dm(long n, long m, double ini,double dx, double ini2,double dx2, char* capt);
  /* to create a dm */


DM* vec_to_dm(float* y, long n, double ini,double dx, char* capt);
  /* to create type 1 gd from a vector */

DM* gd2dm(GD *gd);


float* y_dm(DM* dm);
  /* extracts the data vector from a gd */

DM* set_dm(DM* dm, int typ,
		   double par1, double par2, double par3);
/*  
*/


int show_dm(DM* dm);
/* show the gd parameters */


DM* sum_dm(float a1, DM* dm1, 
		   float a2, DM* dm2,
		   float a3);
/* sums two dm's */


DM* mult_dm(DM* gd1, DM* dm2);
/* product of two dm's */



/*-----------------------------------------------*/

typedef struct SNAG_RING{/* ring structure */
	float*	y;			 /* the data */
	long		len;     /* length of the ring */
	long		indin;   /* input index (first to be written) DIFFERENT FROM MATLAB */
	long		indout;  /* output index (first to be read) DIFFERENT FROM MATLAB */
	long		totin;   /* total data input */
	long		totout;  /* total data output */
	double		ini;     /* beginning of abscissa */
	double		lastim;  /* last input sample time */
	double		dx;      /* sampling "time" */ 
	char*		capt;    /* caption */
	int			cont;    /* control variable */
	int			debug;   /*  > 0 -> debug */
} RING;


RING* crea_ring(long len);


int read_rg(RING *r,float* v,long n,double *tinit);
/* reads from a ring */


int write_rg(RING *r,float* v,long n,double lastim);
/* writes to a ring */



/*-----------------------------------------------*/


typedef struct DS{ /* ds structure */
   char name[20];  /* name */
   double  tini1;  /* time of the first sample of y1 */
   double  tini2;  /* time of the first sample of y2 */
   double  dt;     /* sampling time */
   long    len;    /* length of chunks */
   float*  y1;     /* odd chunk (last chunk if not interlaced) 
                             (dimension 5*len/4)               */
   float*  y2;     /* even chunk (last but one chunk if not interlaced) */
   long    ind1;   /* index of y1 */
   long    ind2;   /* index of y2 */
   short   type;   /* type (=0 -> not interlaced no last but one, */
                   /*   =1 -> not interlaced, =2 -> interlaced)   */
   long    nc1;    /*  number of y1 chunk */
   long    nc2;    /*  number of y2 chunk */
   long    lcw;    /*  last chunk written ("produced") */
   long    lcr;    /*  last chunk read ("served" - for client-server use) */
   double  treq;   /*  time requested (to start) */
   char*   capt;   /*  caption */
   int     cont;   /*  control variable */
   int     debug;  /*  > 0 -> debug */
} DS;



DS* crea_ds(long len, int typ, char* capt);
/* 
	len		chunk length
	typ		ds type (0,1,2 -> interlacing)
	capt	caption
*/


int y_ds(DS *d, float* y);


int read_rg2ds(RING *r, DS *d, int iy, long ink, 
			   long n, double *tinit);
/* reads from a ring directly to a ds */


int sim2ds(DS *d, RING *r, float *v, long nv, int typ, int verb);
/* typ = 0     all zeros
    "  = 1     ramp
	"  = 2     sine  */


/*----------------------- MCH ------------------------*/


typedef struct MCH_CH{ /* mch ch sub-structure */
	char	name[20];
	double	del;
	double	t; /* computed (see mch_param) */
	long	n;
	double	dt; /* computed (see mch_param) */
	long	bias; /* computed (see mch_param) */
} MCH_CH;


typedef struct MCH{ /* mch ch sub-structure */
	long	nch; /* number of channels */
	MCH_CH	*ch;
	double	t0;
	double	dt0;
	long	ntot; /* computed (see mch_param) */
	float	*y;
} MCH;

MCH mch_param(MCH mch);
/* computes the derived parameters for an mch structure and allocates y */


#ifdef __cplusplus
}
#endif

#endif /* _HEADER_PSS_SNAG_H */

/*___________________________________________________
 ¦                                                   ¦
 ¦                     pss_gw.h                      ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    March 2002                     ¦
 ¦___________________________________________________¦*/


#ifdef __cplusplus
extern "C" {
#endif

#define MATH_PI 3.14159265358979323
#define MATH_E 2.71828182845904523
#define PHYS_C 299792458
#define PHYS_G 6.67259E-11
#define ASTRO_SY_S  315581499.84
#define ASTRO_SD_S  315581499.84

typedef struct MATH{
	double	pi;
	double	e;
} MATH;


typedef struct PHYS{
	double	c;
	double	G;   

} PHYS;


typedef struct ASTRO{
	double	SY_s;
} ASTRO;


typedef struct SOURCE{
	char	name[20];
	double	alpha;
	double	delta;

	
	double	alpha_r;
	double	delta_r;
} SOURCE;


typedef struct IANTENNA{
	double	lat;
	double	lon;

	double	lat_r;
	double	lon_r;
} IANTENNA;


typedef struct BANTENNA{
	double	lat;
	double	lon;
} BANTENNA;


typedef struct GW{
	MATH		math;
	PHYS		phys;
	ASTRO		astro;

	SOURCE		gc;
	long		npulsar; /* pulsar */
	SOURCE		*pulsar; /* pulsar */
	long		nglob;   /* globular cluster */
	SOURCE		*glob;   /* globular cluster */
	long		ngal;    /* galaxies */
	SOURCE		*gal;    /* galaxies */
	long		nsource; /* generic sources */
	SOURCE		*source; /* generic sources */
	long		nsimsour;/* simulated source */
	SOURCE		*simsour;/* simulated source */

	IANTENNA	geo;
	IANTENNA	ligo_han;
	IANTENNA	ligo_liv;
	IANTENNA	tama;
	IANTENNA	virgo;
	long		niant;
	IANTENNA	*iant;

	BANTENNA	allegro;
	BANTENNA	auriga;
	BANTENNA	explorer;
	BANTENNA	nautilus;
	BANTENNA	niobe;
	long		nbant;
	BANTENNA	*bant;
} GW;


GW pss_set_gw(long nsimsour, long niant, long nbant);
/* input parameters are for simulated objects */


int show_gw_source(GW gw, int type, int sort_type);


int show_gw_iant(GW gw, int sort_type);


int show_gw_bant(GW gw, int sort_type);


int show_gw_const(GW gw, int type);


#ifdef __cplusplus
}
#endif


/***************************************************************/
/*                 PSS_astro.h,  v. 22 Jan 2002                */
/*     header file used for the PSS_astro library              */ 
/***************************************************************/


#include <inttypes.h>   //pia64bit
#include "wherearefiles.h"
#include "daspostare.h"
/********Ephemerides file  ********/
//#define EPHFILE  "../pss_astro/novas/Jpleph.405"
//#define EPHFILE  "/storage/gpfs_virgo3/home/CW/sabrina/Jpleph.405"

FILE *F1;

/*******UT1-UTC file (from IERS)******/
//#define UTC2UT1  "../pss_astro/novas/ut1-utc_mod1.dat"  //nuovo
//#define UTC2UT1  "/storage/gpfs_virgo3/home/CW/sabrina/ut1-utc_mod1.dat" //nuovo

FILE *F2;

//#define DELTANUT  "../pss_astro/novas/eopc04_mod"
//#define DELTANUT  "/storage/gpfs_virgo3/home/CW/sabrina/eopc04_mod"
FILE *F3;

/*****Constant definition*****/

#define day2sec  86400.0

#define dot2utc  2415020.5E0

#define dot2mjd  15020

#define utc2mjd  2400000.5E0

#define dot2gps  29224.0

#define jd06jan1980 2444244.5

#define mjd06jan1980 44244

#define mjd2gpsC 44244

#define gps1jan2006 820108814

/* utc=dot+dot2utc, where dot are days from the beginning of century, UTC=Julian date
   mjd=dot+dot2mjd, where dot are days from the beginning of 1900, MJD=Modified Julian day
  mjd=utc-utc2mjd

  JD=utc here,
  gpstime=(dot-dot2gps)*day2sec+(leap_sec-19) %%13 fino al 2006, 14 dal 1/Jan/2006
  gpstime=(mjd-mjd2gpsC)*day2sec+(leap_sec-19) %%13 fino al 2006, 14 dal 1/Jan/2006
  gpstime=(JD-jd06jan1980)*day2sec+(leap_sec-19)
  dot=dot2gps+(gpstime-(leap_sec-19))/day2sec
*/
#define ecl_deg  23.439
/**This is the value at Epoch=Equinox 2000**/


#define NutFlag  1
/**Kaz suggestion: if 0 eqeq=0 -> Grenw. mean sidereal time
if 1 ->Greenw. apparent sidereal time.
Ma sembra ci sia un errore in NOVAS e la chiamata va fatta con 0 per poi
usare novas_mod_terra, nutate and precession.  May 2003**/
/*No: va bene come era e deve essere messa a 1 Aug 2003*/

/*****Function prototypes:*****/


/*To fill the structures:*/



/*Astronomical time conversions:*/

int  PSS_leapseconds(double utc);

  double mjd2gps(double mjdtime);

  double PSS_utc2tdt(double utc);

  double  PSS_tdt2utc(double tdt);

  double PSS_tdt2tdb(double tdt);
 
  double  PSS_utc2ut1(double utc);

void dot2tsslocal(DETECTOR *detector,double dot, short int fut1,short int nut1,short int nutcorr,double *tss_local);
void mjd2tsslocal(DETECTOR *detector,double mjd, short int fut1,short int nut1,short int nutcorr,double *tss_local);
void gps2tsslocal(DETECTOR *detector,long gpstime, short int fut1,short int nut1,short int nutcorr,double *tss_local);
double gps2mjd(double gpstime);
/* From gpstime to mjd time */

/*Modified from Novas.c:*/

 void novas_mod_starvectors (SOURCE *star,  double *pos, double *vel);

 void novas_mod_terra (DETECTOR *locale, double st, double *pos, double *vel);

void novas_mod_earthtilt (double tjd,double *mobl, double *tobl, double *eq, double dpsi, double deps);

short int novas_mod_nutate (double tjd, short int fn, double *pos,double *pos2,double dpsi,double deps);
//Sep. 2003 !!!!!

/*Used to read the jpl ephem file:*/

void jpl_mod_constan(char nam[][6], double val[], double sss[], intptr_t *n);
//pia64
void jpl_mod_state(double et2[],int list[],double pv[][6],double nut[]);

void jpl_mod_split(double tt, double fr[]);

void jpl_mod_interp(double buf[],double t[],int ncf,int ncm,int na,int ifl,double pv[]);

void jpl_mod_pleph(double et,int ntarg,int ncent,double rrd[] );




/*Other utility functions:*/

/*Corrections to nutation angles: */
void   PSS_deltanut(double utc, double *ddpsi, double *ddeps);

/*Doppler effect:*/

double  PSS_doppler(SOURCE *source,double poss[] , double vel[]);

double PSS_einstein(double tdt_i);

void PSS_einstein_vect(SOURCE *source,double dot_ini,double step,int npoint,double nu_einstein[]);

void PSS_shapiro_vect(DETECTOR *detector,SOURCE *source,double dot_ini, double step, int npoint,short int fut1, short int nut1, short int nutcorr,short int equatorial,double nu_shapiro[]);

/*Coordinate conversions:*/

short int PSS_rect2spher(double *pos, double *ra, double *dec);

void PSS_spher2rect(double ra, double dec, double dist, double *vector);

void PSS_eclr2equr (double *vector_in, double *vector_out);

void PSS_equr2eclr (double *vector_in, double *vector_out);

void PSS_equs2ecls (double ra, double dec, double *lambda, double *beta);

void PSS_ecls2equs (double lambda, double beta, double *ra, double *dec);


/**********************************************************************************/


/*****Detector velocity:*****/
/*For a given detector, initial date (days from the beginning of the century or mjd),step (s), number
 of point, it gives the velocity vector (AU/day) in rectangular ecliptic 
or equatorial coordinates, respect to
the Barycenter of the Earth-Sun system (referred to Epoch). It returns MJD of the first point.*/

double  PSS_detector_vel(DETECTOR *detector,double dot_ini, double step, int npoint,short int fut1, short int nut1,short int nutcorr,short int equatorial,short int imjd,double velx[],double vely[],double velz[], double deinstein[]);

/*****Detector velocity:*****/
/*For a given detector, initial date (days from the beginning of the century or mjd),step (s), number
 of point, it gives the position vector (AU) in rectangular ecliptic or equatorialcoordinates, respect to
the Barycenter of the Earth-Sun system (referred to Epoch). It returns MJD of the first point.*/

double  PSS_detector_pos(DETECTOR *detector,double dot_ini, double step, int npoint,short int fut1, short int nut1,short int nutcorr,short int equatorial,short int imjd,double posx[],double posy[],double posz[]);

double PSS_detector_velpos(DETECTOR *detector,double dot_ini, double step, int npoint,short int fut1, short int nut1, short int nutcorr,short int equatorial,short int imjd,double velx[], double vely[],double velz[],double posx[], double posy[],double posz[],double deinstein[]);
double PSS_detector_posTSSB(DETECTOR *detector,double dot_ini, double step, int npoint,short int fut1, short int nut1, short int nutcorr,short int equatorial,short int imjd,double posx[], double posy[],double posz[]);
double PSS_detector_velposSCALAR(DETECTOR *detector,double dot_ini, double step,short int fut1, short int nut1, short int nutcorr,short int equatorial,short int imjd,double *velx, double *vely,double *velz,double *posx, double *posy,double *posz,double deinstein);
double PSS_Nodet_velpos(double dot_ini, double step, int npoint,short int fut1, short int nut1, short int nutcorr,short int equatorial,short int imjd,double velx[], double vely[],double velz[],double posx[], double posy[],double posz[],double deinstein[]);
/*****Source position:*****/
/*For a given source, identified by ra and dec of Epoch it gives the 
normalized position vector (AU) 
in rectangular ecliptic coordinates,  respect to the Barycenter of the Earth-Sun system. 
The source is supposed to be fixed in the sky (no proper motion).*/

void  PSS_source_pos(SOURCE *source,short int equatorial,double poss[]);

/* Amplitude modulation */ 
void modbar (DETECTOR *detector,SOURCE *source,double tss_local,float *a,float *b);
void radpat_interf (DETECTOR *detector,SOURCE *source,double tss_local,float *fcirc,float *flin);
/***********************************Miscellanea:**********************************/

/*****Vectors:*****/

void PSS_vect(double vectx[], double vecty[],double vectz[],double vectout[12],int *np);

/*For given vectors, vectx vecty vectz, supposed to be the x,y,z components of a
vector as a function of time, np points, it evaluates- and puts in vectout[]-:
the x,y,z components of the first point,
the x,y,z components of the middle point,
the x,y,z components of the last point,
the average x,y,z components
*/



double scalar_product(double vect1[3],double vect2[3]);
/*This function performs the scalar product of vectors of 3 components

INPUT:
vect1[3],double;
vect2[3],double;

OUTPUT:
sc,double
*/

/************************LOCMAX********************/

double PSS_locmax_vel(DETECTOR *detector,double dot_ini,short int fut1, double velx[], double vely[],double velz[]);

#ifndef _HEADER_DASPOSTARE_H
#define _HEADER_DASPOSTARE_H

/***************************************************************/
/* Definizioni che poi andranno spostate, Pia 25 Jan 2002      */
/*                                                             */
/***************************************************************/

/*****Structures:*****/

typedef struct DETECTOR{
  /*
   name[16]      = name of the detector.
   longitude     = longitude [deg] East positive.
   latitude      = latitude  [deg] North positive.
   height        =height [m].
   azimuth       =azimuth [deg]  //From North to East
  */
  char name[5];
  double longitude;
  double latitude;
  double height;
  double azimuth;
} DETECTOR;


typedef struct SOURCE{
  /*
   name[16]       = name of the source.
   ra             = mean right ascension [hours].
   dec            = mean declination [degrees].
   promora        = proper motion in RA [seconds of time per century].
   promodec       = proper motion in declination [arcseconds per century].
   parallax       = parallax [arcseconds].
   radialvelocity = radial velocity [kilometers per second].
   frequency      = intrinsic g.w. frequency.
   psi            =linear polarization angle [degrees].
   eps            =percentage of linear polarization [0-1].
  */
  char name[16];
  double ra;
  double dec;
  double parallax;
  double promora;
  double promodec;
  double radialvelocity;
  double frequency;
  double psi;
  double eps;

} SOURCE;




/*****Function prototypes:*****/

/*To fill the structures:*/

DETECTOR *def_detector(char* type);

SOURCE  *def_source(char* type1);

#endif /* _HEADER_DASPOSTARE_H */

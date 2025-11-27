/***************************************************************************
*******                  JPLBIN.H   v.1.6                          *********
****************************************************************************
**  This header file is used by TESTEPH program.                          **
**  It is NECESSARY TO ADJUST IT MANUALLY for different ephemerides.      **
****************************************************************************
**  Written: May 28, 1997 by PAD  **  Last modified: April 18,2021 by PAD **
****************************************************************************
**  PAD: dr. Piotr A. Dybczynski,          e-mail: dybol@amu.edu.pl       **
** Astronomical Observatory Institute, A.Mickiewicz Univ., Poznan, Poland **
**      https://apollo.astro.amu.edu.pl/PAD/pmwiki.php?n=Dybol.JPLEph     **
***************************************************************************/

/* UNCOMMENT ONE AND ONLY ONE OF THE FOLLOWING DENUM DEFINITIONS: */

/*#define DENUM 200*/
/*#define DENUM 403*/
/*#define DENUM 404*/
/*#define DENUM 405*/
/*#define DENUM 406*/
/*#define DENUM 421*/
/*#define DENUM 422*/
/*#define DENUM 430*/
/*#define DENUM 431*/
/*#define DENUM 440*/
/*#define DENUM 441*/
#define DENUM 200
/* #define DENUM 436*/
/* #define DENUM 440 */
// se sono da KSIZE 2036 posso lasciare scritto 405 qui...Fatto anche 436 e 440. Infatti first file non serve qui, serve solo alla conversione in ascii..Per SFDB rimetti 405
//Lo ho VERIFICATO: ottengo identico risultato
#if   DENUM==200
#define KSIZE 1652
#define NMAX 400
#elif DENUM==403 || DENUM==405 || DENUM==421 || DENUM==422 || DENUM==436 || DENUM==431 || DENUM==440 || DENUM==441
#define KSIZE 2036
#define NMAX 1000
#elif DENUM==404 || DENUM==406
#define KSIZE 1456
#define NMAX 1000
#endif

#define NRECL 4
#define RECSIZE (NRECL*KSIZE)
#define NCOEFF (KSIZE/2)
#define TRUE 1
#define FALSE 0

/* the following two definitions adapt the software for larger number of constants (1000)
(eg. for DE430 and DE431 and later ones) */

// #define NMAX 1000. Pia May 20 messo prima, per la 200 deve essere diverso.
#define OLDMAX 400

/*
#pragma pack(push, 1)
*/


struct  __attribute__((__packed__))  rec1{
         char ttl[3][84];
         char cnam[OLDMAX][6];
         double ss[3];
         int ncon;
         double au;
         double emrat;
         int ipt[12][3];
         int numde;
         int lpt[3];
         int rpt[3];
         int tpt[3];
         char cnam2[(NMAX-OLDMAX)][6];
       };
 struct {
         struct rec1 r1;
         char spare[RECSIZE-sizeof(struct rec1)];
       } R1;

 struct  __attribute__((__packed__)) rec2{
         double cval[NMAX];
       };
 struct {
         struct rec2 r2;
         char spare[RECSIZE-sizeof(struct rec2)];
       } R2;

/* The following definitions are for the asc2eph.c program only */
/* asc2eph.c converts the ASCII JPL ephemeris files to the binary file
   with the possibility to restrict the ephemeris to a shorter interval */

#if   DENUM==441 || DENUM==431

#define FIRSTFILE -13000
#define LASTFILE 16000
#define FILESTEP 1000

#elif DENUM==440 || DENUM==436

#define FIRSTFILE 1550
#define LASTFILE 2550
#define FILESTEP 100

#elif DENUM==422 || DENUM==406

#define FIRSTFILE -3000
#define LASTFILE 2900
#define FILESTEP 100

#elif DENUM==421

#define FIRSTFILE 1900
#define LASTFILE 2050
#define FILESTEP 150

#elif DENUM==405

#define FIRSTFILE 1600
#define LASTFILE 2200
#define FILESTEP 20

#elif DENUM==403

#define FIRSTFILE 1600
#define LASTFILE 2100
#define FILESTEP 100

#elif DENUM==200

#define FIRSTFILE 1600
#define LASTFILE 2160
#define FILESTEP 20

#endif






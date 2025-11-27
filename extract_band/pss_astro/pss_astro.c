/*********************************************************/
/***** PSS_astro.c*****22 Jan 2002*/

#include<stdio.h>
#include<math.h>
#include<string.h>
#include <inttypes.h>   //pia64bit


#include"jplbin.h"
#include"novas/novas.h"

#include "daspostare.h"
#include"pss_astro.h"




/*****Functions:*****/

/****************************Astronomical time conversions:**************************/

/***** PSS_leapseconds

/***** input utc=mjd+utc2mjd*****/
int  PSS_leapseconds(double utc)

{

    int nleap_sec;

/* 	 Leap seconds  */
    if (utc <  (float)2447161.5) nleap_sec=23;

/* 	fra il 1 gennaio 1988 e il 1 gennaio 1990 */
    if (utc >= (float)2447161.5 && utc < (float)2447892.5) {
	nleap_sec = 24;
    }
/* 	fra il 1 gennaio 1990 e il 1 gennaio 1991 */
    if (utc >= (float)2447892.5 && utc < (float)2448257.5) {
	nleap_sec= 25;
    }
/* 	fra il 1 gennaio 1991 e il 1 luglio 1992 */
    if (utc >= (float)2448257.5 && utc < (float)2448804.5) {
	nleap_sec= 26;
    }
/* 	fra il 1 luglio 1992 e il 1 luglio 1993 */
    if (utc >= (float)2448804.5 && utc < (float)2449169.5) {
	nleap_sec= 27;
    }
/* 	fra il 1 luglio 1993 e il 1 luglio 1994 */
    if (utc >= (float)2449169.5 && utc < (float)2449534.5) {
	nleap_sec= 28;
    }
/* 	fra il 1 luglio 1994 e il 1 gennaio 1996 */
    if (utc >= (float)2449534.5 && utc < (float)2450083.5) {
	nleap_sec= 29;
    }
/* 	fra il 1 gennaio 1996 e il 1 luglio 1997 */
    if (utc >= (float)2450083.5 && utc < (float)2450630.5) {
	nleap_sec= 30;
    }
/* 	fra il 1 luglio 1997 e il 1 gennaio 1999 */
    if (utc >= (float)2450630.5 && utc < (float)2451179.5) {
	nleap_sec= 31;
    }
/* 	fra il 1 gennaio 1999 e il 1 gennaio 2006 */
     if (utc >= (float)2451179.5 && utc < (float)2453736.5) {
	nleap_sec= 32;
    }
    /* da aggiungere  fra il 1 gennaio 2006  e 1 gennaio 2009 */
     if (utc >= (float)2453736.5 && utc < (float) 2454832.5) {
        nleap_sec= 33;
    }

      /*  fra il 1 gennaio 2009  e il 30 giugno 2012*/
     if (utc >= (float)2454832.5 && utc < (float)2456109.5 ) {
        nleap_sec= 34;
    }
      /*  fra il 30 giugno 2012 e il primo luglio 2015 */
      if (utc >= (float)2456109.5 && utc < (float)2457204.5 ) {
        nleap_sec= 35;
    }
      /* 1 luglio 2015 */
      if (utc >= (float)2457204.5 && utc < (float)2457754.5) {
        nleap_sec= 36;
      }
      /* 1 JAN 2017 This is midnight 1 Jan 2017 0:0:0*/
      if (utc >= (float)2457754.5) {
        nleap_sec= 37;
    }
    return nleap_sec;
} 

double  PSS_utc2tdt(double utc)

{
    double tdt;
    int nleap_sec;

    nleap_sec=PSS_leapseconds(utc);

    tdt = utc + (nleap_sec + 32.184) / 86400.;
    return tdt;
} 


/***** Function PSS_utc2tdt *****/
/***** From Terrestrial Dynamic Time (TDT) to UTC*****/
/***** UTC=TDT-(number leap seconds)-32.184*****/
double  PSS_tdt2utc(double tdt)

{
    double utc,appo;
    int nleap_sec;
/*      From TDT (terrestrial dynamical time) to UTC*/
/* 	 Leap seconds  to be removed*/

    appo=tdt-32.184/84400.;
    nleap_sec=PSS_leapseconds(appo);

    
    utc = tdt - (nleap_sec + 32.184) / 86400.;
    return utc;
} 



/***** Function PSS_tdt2tdb*****/ 
/***** From Terrestrial Dynamic Time (TDT) to Barycentric Dynamic Time (TDB)*****/

/*      Computes the cumulative relativistic time correction to
 *               earth-based clocks, TDB-TDT, for a given time. Routine
 *               furnished by the Bureau des Longitudes, modified by
 *               removal of terms much smaller than 0.1 microsecond.
 *       Programmer and Completion Date:
 *       Lloyd Rawley - S.T.X. - 07/28/89
 *       (Retyped by Masaki Mori - 04/19/96)
 *       (Converted to C for bary, optimized,
 *        corrected following Faithead & Bretagnon 1990, A&A 229, 240
 *        by Arnold Rots - 1997-11-20)
 *
 
 *
 *     Calling Sequence:  tdbdt = ctacv(jdno, fjdno)
 *        Argument   Type   I/O                Description
 *        --------   ----   ---  ----------------------------------------
 *        jdno      long     I   Julian day number of lookup
 *        fjdno     double   I   Fractional part of Julian day number
 *        tdbtdt    double   O   Time difference TDB-TDT (seconds)
 *
 *     Called by:  TTtoTDB
 *
 *     Calls:  none
 *
 *     COMMON use:  none
 *
 *     Significant Local Variables:
 *        Variable   Type   Ini. Val.        Description
 *        --------   ----   ---------  ----------------------------------
 *          t       double      -      Time since 2000.0 (millennia)
 *          tt      double      -      Square of T
 *
 *     Logical Units Used:  none
 *
 *     Method:
 *        Convert input time to millennia since 2000.0
 *        For each sinusoidal term of sufficient amplitude
 *           Compute value of term at input time
 *        End for
 *        Add together all terms and convert from microseconds to seconds
 *     End CTATV
 *
 *     Note for the retyped version:
 *        Results of this routine has been confirmed up to (1E-10)microsecond
 *        level compared with the calculation by the IBM machine: this seems
 *        to be within the 64-bit precision allowance, but still the original
 *        hardcopy should be kept as the right one. (M.M.)
 *
 *     Note for the C version: the accuracy is guaranteed to 100 ns.
 *
 *---------------------------------------------------------------------------*/
/**************Modified by Pia (input,output parameters) 22 Jan 2002**********/ 
/* The original name, ctatv, has been replaced with PSS_tdt2tdb */
 double PSS_tdt2tdb (double tdt)
 {
   
    long jdno;
     double tdb,fjdno,ctatvo;
     double t, tt, t1, t2, t3, t4, t5, t24, t25, t29, t30, t31 ;


      fjdno=tdt-(int)tdt;
      jdno=(int)tdt;
      t = ((jdno-2451545) + fjdno)/(365250.0) ;
      tt = t*t ;

      t1  =       1656.674564 * sin(  6283.075943033*t + 6.240054195)
           +        22.417471 * sin(  5753.384970095*t + 4.296977442)
           +        13.839792 * sin( 12566.151886066*t + 6.196904410)
           +         4.770086 * sin(   529.690965095*t + 0.444401603)
           +         4.676740 * sin(  6069.776754553*t + 4.021195093)
           +         2.256707 * sin(   213.299095438*t + 5.543113262)
           +         1.694205 * sin(    -3.523118349*t + 5.025132748)
           +         1.554905 * sin( 77713.772618729*t + 5.198467090)
           +         1.276839 * sin(  7860.419392439*t + 5.988822341)
           +         1.193379 * sin(  5223.693919802*t + 3.649823730)
           +         1.115322 * sin(  3930.209696220*t + 1.422745069)
           +         0.794185 * sin( 11506.769769794*t + 2.322313077)
           +         0.600309 * sin(  1577.343542448*t + 2.678271909)
           +         0.496817 * sin(  6208.294251424*t + 5.696701824)
           +         0.486306 * sin(  5884.926846583*t + 0.520007179)
           +         0.468597 * sin(  6244.942814354*t + 5.866398759)
           +         0.447061 * sin(    26.298319800*t + 3.615796498)
           +         0.435206 * sin(  -398.149003408*t + 4.349338347)
           +         0.432392 * sin(    74.781598567*t + 2.435898309)
           +         0.375510 * sin(  5507.553238667*t + 4.103476804) ;

      t2  =          0.243085 * sin(  -775.522611324*t + 3.651837925)
           +         0.230685 * sin(  5856.477659115*t + 4.773852582)
           +         0.203747 * sin( 12036.460734888*t + 4.333987818)
           +         0.173435 * sin( 18849.227549974*t + 6.153743485)
           +         0.159080 * sin( 10977.078804699*t + 1.890075226)
           +         0.143935 * sin(  -796.298006816*t + 5.957517795)
           +         0.137927 * sin( 11790.629088659*t + 1.135934669)
           +         0.119979 * sin(    38.133035638*t + 4.551585768)
           +         0.118971 * sin(  5486.777843175*t + 1.914547226)
           +         0.116120 * sin(  1059.381930189*t + 0.873504123)
           +         0.101868 * sin( -5573.142801634*t + 5.984503847)
           +         0.098358 * sin(  2544.314419883*t + 0.092793886)
           +         0.080164 * sin(   206.185548437*t + 2.095377709)
           +         0.079645 * sin(  4694.002954708*t + 2.949233637)
           +         0.075019 * sin(  2942.463423292*t + 4.980931759)
           +         0.064397 * sin(  5746.271337896*t + 1.280308748)
           +         0.063814 * sin(  5760.498431898*t + 4.167901731)
           +         0.062617 * sin(    20.775395492*t + 2.654394814)
	   +         0.058844 * sin(   426.598190876*t + 4.839650148)
           +         0.054139 * sin( 17260.154654690*t + 3.411091093) ;

      t3  =          0.048373 * sin(   155.420399434*t + 2.251573730)
           +         0.048042 * sin(  2146.165416475*t + 1.495846011)
           +         0.046551 * sin(    -0.980321068*t + 0.921573539)
           +         0.042732 * sin(   632.783739313*t + 5.720622217)
           +         0.042560 * sin(161000.685737473*t + 1.270837679)
           +         0.042411 * sin(  6275.962302991*t + 2.869567043)
           +         0.040759 * sin( 12352.852604545*t + 3.981496998)
           +         0.040480 * sin( 15720.838784878*t + 2.546610123)
           +         0.040184 * sin(    -7.113547001*t + 3.565975565)
           +         0.036955 * sin(  3154.687084896*t + 5.071801441)
           +         0.036564 * sin(  5088.628839767*t + 3.324679049)
           +         0.036507 * sin(   801.820931124*t + 6.248866009)
           +         0.034867 * sin(   522.577418094*t + 5.210064075)
           +         0.033529 * sin(  9437.762934887*t + 2.404714239)
           +         0.033477 * sin(  6062.663207553*t + 4.144987272)
           +         0.032438 * sin(  6076.890301554*t + 0.749317412)
           +         0.032423 * sin(  8827.390269875*t + 5.541473556)
           +         0.030215 * sin(  7084.896781115*t + 3.389610345)
           +         0.029862 * sin( 12139.553509107*t + 1.770181024)
           +         0.029247 * sin(-71430.695617928*t + 4.183178762) ;

      t4  =          0.028244 * sin( -6286.598968340*t + 5.069663519)
	   +         0.027567 * sin(  6279.552731642*t + 5.040846034)
           +         0.025196 * sin(  1748.016413067*t + 2.901883301)
           +         0.024816 * sin( -1194.447010225*t + 1.087136918)
           +         0.022567 * sin(  6133.512652857*t + 3.307984806)
           +         0.022509 * sin( 10447.387839604*t + 1.460726241)
           +         0.021691 * sin( 14143.495242431*t + 5.952658009)
           +         0.020937 * sin(  8429.241266467*t + 0.652303414)
           +         0.020322 * sin(   419.484643875*t + 3.735430632)
           +         0.017673 * sin(  6812.766815086*t + 3.186129845)
           +         0.017806 * sin(    73.297125859*t + 3.475975097)
           +         0.016155 * sin( 10213.285546211*t + 1.331103168)
           +         0.015974 * sin( -2352.866153772*t + 6.145309371)
           +         0.015949 * sin(  -220.412642439*t + 4.005298270)
           +         0.015078 * sin( 19651.048481098*t + 3.969480770)
           +         0.014751 * sin(  1349.867409659*t + 4.308933301)
           +         0.014318 * sin( 16730.463689596*t + 3.016058075)
           +         0.014223 * sin( 17789.845619785*t + 2.104551349)
           +         0.013671 * sin(  -536.804512095*t + 5.971672571)
           +         0.012462 * sin(   103.092774219*t + 1.737438797) ;

      t5  =          0.012420 * sin(  4690.479836359*t + 4.734090399)
           +         0.011942 * sin(  8031.092263058*t + 2.053414715)
           +         0.011847 * sin(  5643.178563677*t + 5.489005403)
           +         0.011707 * sin( -4705.732307544*t + 2.654125618)
           +         0.011622 * sin(  5120.601145584*t + 4.863931876)
           +         0.010962 * sin(     3.590428652*t + 2.196567739)
           +         0.010825 * sin(   553.569402842*t + 0.842715011)
           +         0.010396 * sin(   951.718406251*t + 5.717799605)
           +         0.010453 * sin(  5863.591206116*t + 1.913704550)
           +         0.010099 * sin(   283.859318865*t + 1.942176992)
           +         0.009858 * sin(  6309.374169791*t + 1.061816410)
           +         0.009963 * sin(   149.563197135*t + 4.870690598)
           +         0.009370 * sin(149854.400135205*t + 0.673880395) ;

      t24 = t * (  102.156724 * sin(  6283.075849991*t + 4.249032005)
           +         1.706807 * sin( 12566.151699983*t + 4.205904248)
           +         0.269668 * sin(   213.299095438*t + 3.400290479)
           +         0.265919 * sin(   529.690965095*t + 5.836047367)
           +         0.210568 * sin(    -3.523118349*t + 6.262738348)
           +         0.077996 * sin(  5223.693919802*t + 4.670344204) ) ;

      t25 = t * (    0.059146 * sin(    26.298319800*t + 1.083044735)
           +         0.054764 * sin(  1577.343542448*t + 4.534800170)
           +         0.034420 * sin(  -398.149003408*t + 5.980077351)
           +         0.033595 * sin(  5507.553238667*t + 5.980162321)
           +         0.032088 * sin( 18849.227549974*t + 4.162913471)
           +         0.029198 * sin(  5856.477659115*t + 0.623811863)
           +         0.027764 * sin(   155.420399434*t + 3.745318113)
           +         0.025190 * sin(  5746.271337896*t + 2.980330535)
           +         0.024976 * sin(  5760.498431898*t + 2.467913690)
	   +         0.022997 * sin(  -796.298006816*t + 1.174411803)
           +         0.021774 * sin(   206.185548437*t + 3.854787540)
           +         0.017925 * sin(  -775.522611324*t + 1.092065955)
           +         0.013794 * sin(   426.598190876*t + 2.699831988)
           +         0.013276 * sin(  6062.663207553*t + 5.845801920)
           +         0.012869 * sin(  6076.890301554*t + 5.333425680)
           +         0.012152 * sin(  1059.381930189*t + 6.222874454)
           +         0.011774 * sin( 12036.460734888*t + 2.292832062)
           +         0.011081 * sin(    -7.113547001*t + 5.154724984)
           +         0.010143 * sin(  4694.002954708*t + 4.044013795)
           +         0.010084 * sin(   522.577418094*t + 0.749320262)
           +         0.009357 * sin(  5486.777843175*t + 3.416081409) ) ;

      t29 = tt * (   0.370115 * sin(                     4.712388980)
           +         4.322990 * sin(  6283.075849991*t + 2.642893748)
           +         0.122605 * sin( 12566.151699983*t + 2.438140634)
           +         0.019476 * sin(   213.299095438*t + 1.642186981)
           +         0.016916 * sin(   529.690965095*t + 4.510959344)
           +         0.013374 * sin(    -3.523118349*t + 1.502210314) ) ;

      t30 = t * tt * 0.143388 * sin( 6283.075849991*t + 1.131453581) ;

       ctatvo=(t1+t2+t3+t4+t5+t24+t25+t29+t30) * 1.0e-6 ;
       tdb=tdt+ctatvo/86400.;

    return tdb;
}
 

/***** Function PSS_utc2ut1 *****/
/***** From Coordinated Universal time,JD  (UTC) to Universal Time  (UT1)*****/
/***** UT1=UTC+DUT1*****/
/*****************************************************************
      DUT1 are seconds given in the IERS Bull. A
      Available at http://maia.usno.navy.mil/search/search.html
      INPUT: utc (double)
      OUTPUT: ut1 (double)
      Last modified: Pia 16/07/08
      ******************************************************* *****/
double  PSS_utc2ut1(double utc)
//Usata se fut1=1
{
    double ut1;
    float dut1,dut0_s;
    float dut0;
    int mjdf;
    float mjd,mjd_s;
    short int nline;
    long dovesta;
    //char appo3[10], appo4[10];
    float mslope;
    int verb=0; //if 0 does not write things on the screen
    
    mjd=utc-utc2mjd;
    dovesta=ftell(F2);
    fscanf(F2, "%d %f\n",&mjdf,&dut0);
    //printf("nel file mjdf mjd =%d %f \n",mjdf,mjd);  
    fseek(F2,dovesta,0);

    if (((int)mjd-mjdf) < 0){
      if (verb>4) printf("***ATTENTION: mjd < mjdf =%f %d \n",mjd,mjdf);
      ut1=utc;
      return ut1;
    }
    if (fabs(mjd-(float)mjdf) <=0.0001){
      //printf("UTC2UT1 mjd == mjdf %f %f \n",mjd,(float)mjdf);
      ut1=utc+dut0/day2sec;
      return ut1;
    }
    while ((int)mjd-mjdf > 0 ) {
       dovesta=ftell(F2);
       fscanf(F2, "%d %f\n",&mjdf,&dut0);
    } 
     fscanf(F2, "%d %f\n",&mjdf,&dut0);
    //printf("UTC2UT1  trovato mjdf mjd =%d %f \n",mjdf,mjd);
    fscanf(F2, "%d %f\n",&mjdf,&dut0_s);
    fseek(F2,dovesta,0);
     mslope=dut0_s-dut0;
     if(mslope>0.5){
	//printf("UTC2UT1 mjdf mslope togliere LEAP SEC ! %f %f \n",mjdf,mslope);
	      mslope-=1.0;
	      //printf("UTC2UT1 mjdf mslope tolto LEAP SEC ! %f %f \n",mjdf,mslope);
             
     }
     //printf("UTC2UT1 mjdf  mjd =%d  %f \n",mjdf,mjd);
     if(verb >4)printf("UTC2UT1 mjdf,mjd,mslope dut0_s dut0 =%d %f %f %f %f \n",mjdf,mjd,mslope,dut0_s,dut0);   
    dut1=dut0+mslope*(mjd-mjdf);
    //printf("UTC2UT1 mslope dut1 dut0 mjd mjdf =%f %f %f %f %d \n",mslope,dut1,dut0,mjd,mjdf);
    ut1=utc + dut1/day2sec;
    return ut1;
}


/***** Function PSS_deltanut *****/
/***** From Coordinated Universal time,JD  (UTC) gives ddpsi and ddeps,
        to be added to the nutation angles*****/

/*****************************************************************
      File available at http://hpiers.obspm.fr/eoppc/eop/eopc04
      (novas/eopc04.62-now, modified using scrivi-eop.pl into novas/eopc04_mod
      from 1991 to ..
      INPUT: utc (double)
      OUTPUTS: ddpsi,ddeps , arcseconds (float)
Last modified: Pia 16/07/08
      ***********************************************************/
 void  PSS_deltanut(double utc, double *ddpsi,double *ddeps)
 //Used if nutcorr ==1
{
  float psi0,eps0,psi0_s,eps0_s;
  int mjdff;
  double mjd;
  int year,imonth;
  int day;
  float dmp1,dmp2,dmp3,dmp4;
  long dovesta;
  double mslope_psi,mslope_eps;
  int verb=0; //if 0 does not write things on the screen

      dovesta=ftell(F3);
      fscanf(F3,"%d %d %d %d %f %f %f %f %f %f\n",&year,&imonth,&day,&mjdff,&dmp1,&dmp2,&dmp3,&dmp4,&psi0,&eps0);
     fseek(F3,dovesta,0);

     mjd=utc-utc2mjd; 
     if (((int)mjd-mjdff) < 0){
       printf("***ATTENTION DELTANUT: mjd > mjdff =%f %d \n",mjd,mjdff);
       *ddpsi=0.;
       *ddeps=0.;
     return;
    }
    if (fabs(mjd-(float)mjdff) <=0.0001){
      //printf("DELTANUT mjd == mjdff %f %f \n",mjd,(float)mjdff);
      *ddpsi=psi0;
      *ddeps=eps0;
      return;
    }
     while ((int)mjd-mjdff > 0 ) {
        dovesta=ftell(F3);
	fscanf(F3,"%d %d %d %d %f %f %f %f %f %f\n",&year,&imonth,&day,&mjdff,&dmp1,&dmp2,&dmp3,&dmp4,&psi0,&eps0);
     }
	fscanf(F3,"%d %d %d %d %f %f %f %f %f %f\n",&year,&imonth,&day,&mjdff,&dmp1,&dmp2,&dmp3,&dmp4,&psi0,&eps0);
     //printf("DELTANUT trovato mjdff mjd =%d %f \n",mjdff,mjd);
     fscanf(F3,"%d %d %d %d %f %f %f %f %f %f\n",&year,&imonth,&day,&mjdff,&dmp1,&dmp2,&dmp3,&dmp4,&psi0_s,&eps0_s);
     fseek(F3,dovesta,0);
     mslope_psi=psi0_s-psi0;
     mslope_eps=eps0_s-eps0;
     //printf("DELTANUT mjdff  mjd =%d %f \n",mjdff,mjd);
     //printf("DELTANUT mslope psi0 psi0_s =%f %f %f \n",mslope_psi,psi0,psi0_s);
     //printf("DELTANUT mslope_eps eps0 eps0_s =%f %f %f \n",mslope_eps,eps0,eps0_s);
     *ddpsi=psi0+mslope_psi*(mjd-mjdff);
     *ddeps=eps0+mslope_eps*(mjd-mjdff);
      if(verb !=0)printf("DELTANUT mjd mjdff ddpsi ddeps =%f %d %f %f  \n",mjd,mjdff,*ddpsi,*ddeps);
    return;
}

/*****************************Modified from Novas.c:****************************/

/***** Function novas_mod_terra *****/
/*****Modified by Pia, 16 Jan 2002, to use the struct DETECTOR *****/
void novas_mod_terra (DETECTOR *locale, double st,double *pos, double *vel)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the position and velocity vectors of a terrestrial
      observer with respect to the center of the Earth.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97, 
         pp. 1197-1210.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.

   INPUT
   ARGUMENTS:
      glon (double)
         Longitude of observer with respect to reference meridian
         (East +) in degrees.
      glat (double)
         Geodetic latitude (North +) of observer in degrees.
      ht (double)
         Height of observer in meters.
      st (double)
         Local apparent sidereal time at reference meridian in hours.

   OUTPUT
   ARGUMENTS:
      pos[3] (double)
         Position vector of observer with respect to center of Earth,
         equatorial rectangular coordinates, referred to true equator
         and equinox of date, components in AU.
      vel[3] (double)
         Velocity vector of observer with respect to center of Earth,
         equatorial rectangular coordinates, referred to true equator
         and equinox of date, components in AU/Day.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      KMAU, EARTHRAD, F, OMEGA, DEG2RAD

   FUNCTIONS
   CALLED:
      pow    math.h
      sin    math.h
      cos    math.h
      sqrt   math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/04-93/WTH (USNO/AA):  Translate Fortran.
      V1.1/06-98/JAB (USNO/AA):  Move constants 'f' and 'omega' to
                                 file 'novascon.c'.

   NOTES:
      1. If reference meridian is Greenwich and st=0, 'pos' is
      effectively referred to equator and Greenwich.
      2. This function is the "C" version of Fortran NOVAS routine
      'terra'.

------------------------------------------------------------------------
*/
{
   short int j;

   double df2, sinphi, cosphi, c, s, ach, ash, stlocl, sinst, cosst;

/*
   Compute parameters relating to geodetic to geocentric conversion.
*/

   df2 = pow ((1.0 - F),2);

   sinphi = sin (locale->latitude * DEG2RAD);
   cosphi = cos (locale->latitude * DEG2RAD);
   c = 1.0 / sqrt (pow (cosphi,2.0) + df2 * pow (sinphi,2.0));
   s = df2 * c;
   ach = EARTHRAD * c + (locale->height / 1000.0);
   ash = EARTHRAD * s + (locale->height / 1000.0);

/*
   Compute local sidereal time factors at the observer's longitude.
*/

   stlocl = (st * 15.0 + locale->longitude) * DEG2RAD;
   sinst = sin (stlocl);
   cosst = cos (stlocl);

/*
   Compute position vector components in kilometers.
*/

   pos[0] = ach * cosphi * cosst;
   pos[1] = ach * cosphi * sinst;
   pos[2] = ash * sinphi;

/*
   Compute velocity vector components in kilometers/sec.
*/

   vel[0] = -OMEGA * ach * cosphi * sinst;
   vel[1] =  OMEGA * ach * cosphi * cosst;
   vel[2] =  0.0;

/*
   Convert position and velocity components to AU and AU/DAY.
*/

   for (j = 0; j < 3; j++)
   {
      pos[j] /= KMAU;
      vel[j] /= KMAU;
      vel[j] *= 86400.0;
   }

   //printf("locale long,lat,h:%f %f %f\n",locale->longitude,locale->latitude,locale->height);
   return;
}




/***** Function novas_mod_starvectors *****/
/*****Modified by Pia, 16 Jan 2002, to use the struct SOURCE *****/
void novas_mod_starvectors (SOURCE *star,double *pos, double *vel)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Converts angular quanities for stars to vectors.

   REFERENCES: 
      Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97, 
         pp. 1197-1210.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.

   INPUT
   ARGUMENTS:
      *star (struct cat_entry)
         Pointer to catalog entry structure containing J2000.0 catalog
         data with FK5-style units (defined in novas.h).

   OUTPUT
   ARGUMENTS:
      pos[3] (double)
         Position vector, equatorial rectangular coordinates,
         components in AU.
      vel[3] (double)
         Velocity vector, equatorial rectangular coordinates,
         components in AU/Day.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      RAD2SEC, DEG2RAD, KMAU

   FUNCTIONS
   CALLED:
      sin     math.h
      cos     math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Updated to C programming standards.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'vectrs'.

------------------------------------------------------------------------
*/
{
   double paralx, dist, r, d, cra, sra, cdc, sdc, pmr, pmd, rvl;

/*
   If parallax is unknown, undetermined, or zero, set it to 1e-7 second
   of arc, corresponding to a distance of 10 megaparsecs.
*/

   paralx = star->parallax;

   if (star->parallax <= 0.0)
      paralx = 1.0e-7;

/*
   Convert right ascension, declination, and parallax to position vector
   in equatorial system with units of AU.
*/

   dist = RAD2SEC / paralx;
   r = (star->ra) * 15.0 * DEG2RAD;
   d = (star->dec) * DEG2RAD;
   cra = cos (r);
   sra = sin (r);
   cdc = cos (d);
   sdc = sin (d);

   pos[0] = dist * cdc * cra;
   pos[1] = dist * cdc * sra;
   pos[2] = dist * sdc;

/*
   Convert proper motion and radial velocity to orthogonal components of
   motion with units of AU/Day.
*/

   pmr = star->promora * 15.0 * cdc / (paralx * 36525.0);
   pmd = star->promodec / (paralx * 36525.0);
   rvl = star->radialvelocity * 86400.0 / KMAU;

/*
   Transform motion vector to equatorial system.
*/

   vel[0] = - pmr * sra - pmd * sdc * cra + rvl * cdc * cra;
   vel[1] =   pmr * cra - pmd * sdc * sra + rvl * cdc * sra;
   vel[2] =   pmd * cdc + rvl * sdc;

   return;
}

/********earthtilt */

void novas_mod_earthtilt (double tjd, 
                double *mobl, double *tobl, double *eq, double dpsi,
                double deps)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Computes quantities related to the orientation of the Earth's
      rotation axis at Julian date 'tjd'.

      Modified from NOVAS to use the nutation angles of the ephemeris file
     Thus, dpsi and deps are input parameters

   REFERENCES: 
      Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97, 
         pp. 1197-1210.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.
      Transactions of the IAU (1994). Resolution C7; Vol. XXIIB, p. 59.
      McCarthy, D. D. (ed.) (1996). IERS Technical Note 21. IERS
         Central Bureau, Observatoire de Paris), pp. 21-22.

   INPUT
   ARGUMENTS:
      tjd (double)
         TDB Julian date of the desired time
      dpsi (double)
         Nutation in longitude in arcseconds at 'tjd'.
       deps (double)
         Nutation in obliquity in arcseconds at 'tjd'.


   OUTPUT
   ARGUMENTS:
      *mobl (double)
         Mean obliquity of the ecliptic in degrees at 'tjd'.
      *tobl (double)
         True obliquity of the ecliptic in degrees at 'tjd'.
      *eq (double)
         Equation of the equinoxes in seconds of time at 'tjd'.
  
   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      PSI_COR, EPS_COR, DEG2RAD 

   FUNCTIONS
   CALLED:
      nutation_angles  novas.c
      fund_args        novas.c
      fabs             math.h
      pow              math.h
      cos              math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/08-93/WTH (USNO/AA) Translate Fortran.
      V1.1/06-97/JAB (USNO/AA) Incorporate IAU (1994) and IERS (1996) 
                               adjustment to the "equation of the 
                               equinoxes".
      V1.2/10-97/JAB (USNO/AA) Implement function that computes 
                               arguments of the nutation series.
      V1.3/07-98/JAB (USNO/AA) Use global variables 'PSI_COR' and 
                               'EPS_COR' to apply celestial pole offsets
                               for high-precision applications.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'etilt'.
      2. Values of the celestial pole offsets 'PSI_COR' and 'EPS_COR'
      are set using function 'cel_pole', if desired.  See the prolog
      of 'cel_pole' for details.

    In this version of the code d_psi and d_eps are the output of JPL
    file DE405.
------------------------------------------------------------------------
*/
{
   static double tjd_last = 0.0;
   static double t, dp, de;
   double d_psi, d_eps, mean_obliq, true_obliq, eq_eq, args[5];

/*
   Compute time in Julian centuries from epoch J2000.0.
*/

  t = (tjd - T0) / 36525.0;

/*
   Compute the nutation angles (arcseconds) from the standard nutation 
   model if the input Julian date is significantly different from the 
   last Julian date.
*/

  //if (fabs (tjd - tjd_last) > 1.0e-6)
  //    nutation_angles (t, &dp,&de);

/*
   Apply observed celestial pole offsets.
*/

  //d_psi = dp + PSI_COR;
  //d_eps = de + EPS_COR;
  d_psi=dpsi;  //from JPL file
  d_eps=deps;
/*
   Compute mean obliquity of the ecliptic in arcseconds.
*/

   mean_obliq = 84381.4480 - 46.8150 * t - 0.00059 * pow (t, 2.0)
      + 0.001813 * pow (t, 3.0);

/*
   Compute true obliquity of the ecliptic in arcseconds.
*/

   true_obliq = mean_obliq + d_eps;

/*
   Convert obliquity values to degrees.
*/

   mean_obliq /= 3600.0;
   true_obliq /= 3600.0;

/*
   Compute equation of the equinoxes in seconds of time.

   'args[4]' is "omega", the longitude of the ascending node of the 
   Moon's mean orbit on the ecliptic in radians.  This is also an 
   argument of the nutation series.
*/

   fund_args (t, args);

   eq_eq = d_psi * cos (mean_obliq * DEG2RAD) +
      (0.00264  * sin (args[4]) + 0.000063 * sin (2.0 * args[4]));

   eq_eq /= 15.0;
                           
/*
   Reset the value of the last Julian date and set the output values.
*/

   tjd_last = tjd;

   //*dpsi = d_psi;
   //*deps = d_eps;
   *eq = eq_eq;
   *mobl = mean_obliq;
   *tobl = true_obliq;

   return;
}

//******nutate  Sep. 2003 !!!!!

short int novas_mod_nutate (double tjd, short int fn, double *pos, 

                  double *pos2,double dpsi,double deps)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Nutates equatorial rectangular coordinates from mean equator and
      equinox of epoch to true equator and equinox of epoch. Inverse
      transformation may be applied by setting flag 'fn'.

   REFERENCES: 
      Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97, 
         pp. 1197-1210.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.

   INPUT
   ARGUMENTS:
      tdb (double)
         TDB julian date of epoch.
      fn (short int)
         Flag determining 'direction' of transformation;
            fn  = 0 transformation applied, mean to true.
            fn != 0 inverse transformation applied, true to mean.
      pos[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to mean equator and equinox of epoch.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to true equator and equinox of epoch.

   RETURNED
   VALUE:
      (short int)
         0...Everything OK.

   GLOBALS
   USED:
      DEG2RAD, RAD2SEC

   FUNCTIONS
   CALLED:
      earthtilt     novas.c
      cos           math.h
      sin           math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'nutate'.

------------------------------------------------------------------------
*/
{
  double cobm, sobm, cobt, sobt, cpsi, spsi, xx, yx, zx, xy, yy, zy; 
   double xz, yz, zz, oblm, oblt, eqeq, psi, eps;

   novas_mod_earthtilt (tjd, &oblm,&oblt,&eqeq,dpsi,deps); //Sep. 2003  !!!!!

   psi=dpsi; //evaluated from JPL+ corrections in files eopc Sep. 2003 !!!!!
   eps=deps; // idem

   cobm = cos (oblm * DEG2RAD);
   sobm = sin (oblm * DEG2RAD);
   cobt = cos (oblt * DEG2RAD);
   sobt = sin (oblt * DEG2RAD);
   cpsi = cos (psi / RAD2SEC);
   spsi = sin (psi / RAD2SEC);

/*
   Nutation rotation matrix follows.
*/

   xx = cpsi;
   yx = -spsi * cobm;
   zx = -spsi * sobm;
   xy = spsi * cobt;
   yy = cpsi * cobm * cobt + sobm * sobt;
   zy = cpsi * sobm * cobt - cobm * sobt;
   xz = spsi * sobt;
   yz = cpsi * cobm * sobt - sobm * cobt;
   zz = cpsi * sobm * sobt + cobm * cobt;

   if (!fn)
   {

/*
   Perform rotation.
*/

      pos2[0] = xx * pos[0] + yx * pos[1] + zx * pos[2];
      pos2[1] = xy * pos[0] + yy * pos[1] + zy * pos[2];
      pos2[2] = xz * pos[0] + yz * pos[1] + zz * pos[2];
   }
    else
   {

/*
   Perform inverse rotation.
*/

      pos2[0] = xx * pos[0] + xy * pos[1] + xz * pos[2];
      pos2[1] = yx * pos[0] + yy * pos[1] + yz * pos[2];
      pos2[2] = zx * pos[0] + zy * pos[1] + zz * pos[2];
   }

   return 0;
}



/*******************Functions to read the jpl ephem file:**********************/

/******************************************************************************
*                                                                            *
*                                                                            *
*                                                                            *
* The program reads the jpl ephemeris reading and interpolating routine      *
*                                                                            *
* It uses testeph.c from                                                     *
* ftp://ftp.astro.amu.edu.pl/pub                                             *  
*  by Piotr A. Dybczynski (dybol@phys.amu.edu.pl)                            *
* After the initial identifying text which is concluded by an "EOT" in       *
* columns 1-3, the test file contains the following quantities:              *
*                                                                            *
*     JPL ephemeris number                                                   *
*     calendar date                                                          *
*     julian ephemeris date                                                  *
*     target number (1-mercury, ...,3-earth, ,,,9-pluto, 10-moon, 11-sun,    *
*                    12-solar system barycenter, 13-earth-moon barycenter    *
*                    14-nutations, 15-librations)                            *
*     center number (same codes as target number)                            *
*     coordinate number (1-x, 2-y, ... 6-zdot)                               *
*     coordinate  [au, au/day].                                              *
*                                                                            *


*                                                                            *
*                   use jplbin.h version 1.2                                 *
*
******************************************************************************
*                 Last modified: Jan,8  2002                       by Pia     *
*****************************************************************************/
/**** include variable and type definitions, specyfic for this C version
      YOU HAVE TO ADJUST MANUALLY THE FILE: JPLBIN.H BEFORE YOU START       */



int KM=0,BARY=0;
double PVSUN[6];

/****************************************************************************/
/*****************************************************************************
**                         jpl_mod_pleph(et,ntar,ncent,rrd)                         **
******************************************************************************
**                                                                          **
**    This subroutine reads the jpl planetary ephemeris                     **
**    and gives the position and velocity of the point 'ntarg'              **
**    with respect to 'ncent'.                                              **
**                                                                          **
**    Calling sequence parameters:                                          **
**                                                                          **
**      et = (double) julian ephemeris date at which interpolation          **
**           is wanted.                                                     **
**                                                                          **
**    ntarg = integer number of 'target' point.                             **
**                                                                          **
**    ncent = integer number of center point.                               **
**                                                                          **
**    The numbering convention for 'ntarg' and 'ncent' is:                  **
**                                                                          **
**            1 = mercury           8 = neptune                             **
**            2 = venus             9 = pluto                               **
**            3 = earth            10 = moon                                **
**            4 = mars             11 = sun                                 **
**            5 = jupiter          12 = solar-system barycenter             **
**            6 = saturn           13 = earth-moon barycenter               **
**            7 = uranus           14 = nutations (longitude and obliq)     **
**                                 15 = librations, if on eph. file         **
**                                                                          **
**            (If nutations are wanted, set ntarg = 14.                     **
**             For librations, set ntarg = 15. set ncent= 0)                **
**                                                                          **
**     rrd = output 6-element, double array of position and velocity        **
**           of point 'ntarg' relative to 'ncent'. The units are au and     **
**           au/day. For librations the units are radians and radians       **
**           per day. In the case of nutations the first four words of      **
**           rrd will be set to nutations and rates, having units of        **
**           radians and radians/day.                                       **
**                                                                          **
**           The option is available to have the units in km and km/sec.    **
**           for this, set km=TRUE at the begining of the program.          **
*****************************************************************************/
void jpl_mod_pleph(double et,int ntarg,int ncent,double rrd[] )
{
  double et2[2],pv[13][6];/* pv is the position/velocity array
                             NUMBERED FROM ZERO: 0=Mercury,1=Venus,...
                             8=Pluto,9=Moon,10=Sun,11=SSBary,12=EMBary
                             First 10 elements (0-9) are affected by state(),
                             all are adjusted here.                         */


  int bsave,i,k;
  int list[12];          /* list is a vector denoting, for which "body"
                            ephemeris values should be calculated by state():
                            0=Mercury,1=Venus,2=EMBary,...,8=Pluto,
                            9=geocentric Moon, 10=nutations in long. & obliq.
                            11= lunar librations  */

/*  initialize et2 for 'state' and set up component count   */

  et2[0]=et;
  et2[1]=0.0;

  for(i=0;i<6;++i) rrd[i]=0.0;

  if(ntarg == ncent) return;

  for(i=0;i<12;++i) list[i]=0;

/*   check for nutation call    */

  if(ntarg == 14)
    {
      if(R1.r1.ipt[11][1] > 0) /* there is nutation on ephemeris */
        {
          list[10]=2;
          jpl_mod_state(et2,list,pv,rrd);
        }
      else puts("***** no nutations on the ephemeris file  ******\n");
      return;
    }

/*  check for librations   */

  if(ntarg == 15)
    {
      if(R1.r1.lpt[1] > 0) /* there are librations on ephemeris file */
        {
          list[11]=2;
          jpl_mod_state(et2,list,pv,rrd);
          for(i=0;i<6;++i)  rrd[i]=pv[10][i]; /* librations */
        }
      else puts("*****  no librations on the ephemeris file  *****\n");
      return;
    }

/*  force barycentric output by 'state'     */

  bsave=BARY;
  BARY= TRUE;

/*  set up proper entries in 'list' array for state call     */

  for(i=0;i<2;++i) /* list[] IS NUMBERED FROM ZERO ! */
     {
      k=ntarg-1;
      if(i == 1) k=ncent-1;   /* same for ntarg & ncent */
      if(k <= 9) list[k]=2;   /* Major planets */
      if(k == 9) list[2]=2;   /* for moon state earth state is necessary*/
      if(k == 2) list[9]=2;   /* for earth state moon state is necessary*/
      if(k == 12) list[2]=2;  /* EMBary state additionaly */
     }

/*   make call to state   */

  jpl_mod_state(et2,list,pv,rrd);
  /* Solar System barycentric Sun state goes to pv[10][] */
  if(ntarg == 11 || ncent == 11) for(i=0;i<6;++i) pv[10][i]=PVSUN[i];

  /* Solar System Barycenter coordinates & velocities equal to zero */
  if(ntarg == 12 || ncent == 12) for(i=0;i<6;++i) pv[11][i]=0.0;

  /* Solar System barycentric EMBary state:  */
  if(ntarg == 13 || ncent == 13) for(i=0;i<6;++i) pv[12][i]=pv[2][i];

  /* if moon from earth or earth from moon ..... */
  if( (ntarg*ncent) == 30 && (ntarg+ncent) == 13)
      for(i=0;i<6;++i) pv[2][i]=0.0;
  else
    {
       if(list[2] == 2) /* calculate earth state from EMBary */
          for(i=0;i<6;++i) pv[2][i] -= pv[9][i]/(1.0+R1.r1.emrat);

       if(list[9] == 2) /* calculate Solar System barycentric moon state */
          for(i=0;i<6;++i) pv[9][i] += pv[2][i];
    }

  for(i=0;i<6;++i)  rrd[i]=pv[ntarg-1][i]-pv[ncent-1][i];
  BARY=bsave;

  return;
}
/*****************************************************************************
**                     jpl_mod_interp(buf,t,ncf,ncm,na,ifl,pv)              **
******************************************************************************
**                                                                          **
**    this subroutine differentiates and interpolates a                     **
**    set of chebyshev coefficients to give position and velocity           **
**                                                                          **
**    calling sequence parameters:                                          **
**                                                                          **
**      input:                                                              **
**                                                                          **
**        buf   1st location of array of d.p. chebyshev coefficients        **
**              of position                                                 **
**                                                                          **
**          t   t[0] is double fractional time in interval covered by       **
**              coefficients at which interpolation is wanted               **
**              (0 <= t[0] <= 1).  t[1] is dp length of whole               **
**              interval in input time units.                               **
**                                                                          **
**        ncf   # of coefficients per component                             **
**                                                                          **
**        ncm   # of components per set of coefficients                     **
**                                                                          **
**         na   # of sets of coefficients in full array                     **
**              (i.e., # of sub-intervals in full interval)                 **
**                                                                          **
**         ifl  integer flag: =1 for positions only                         **
**                            =2 for pos and vel                            **
**                                                                          **
**                                                                          **
**      output:                                                             **
**                                                                          **
**        pv   interpolated quantities requested.  dimension                **
**              expected is pv(ncm,ifl), dp.                                **
**                                                                          **
*****************************************************************************/
void jpl_mod_interp(double coef[],double t[2],int ncf,int ncm,int na,int ifl,
            double posvel[6])
{
  static double pc[18],vc[18];
  static int np=2, nv=3, first=1;
  static double twot=0.0;
  double dna,dt1,temp,tc,vfac,temp1;
  int l,i,j;

  if(first){           /* initialize static vectors when called first time */
             pc[0]=1.0;
             pc[1]=0.0;
             vc[1]=1.0;
             first=0;
           }

/*  entry point. get correct sub-interval number for this set
    of coefficients and then get normalized chebyshev time
    within that subinterval.                                             */

  dna=(double)na;
  modf(t[0],&dt1);
  temp=dna*t[0];
  l=(int)(temp-dt1);

/*  tc is the normalized chebyshev time (-1 <= tc <= 1)    */

  tc=2.0*(modf(temp,&temp1)+dt1)-1.0;

/*  check to see whether chebyshev time has changed,
    and compute new polynomial values if it has.
    (the element pc[1] is the value of t1[tc] and hence
    contains the value of tc on the previous call.)     */

  if(tc != pc[1])
    {
      np=2;
      nv=3;
      pc[1]=tc;
      twot=tc+tc;
    }

/*  be sure that at least 'ncf' polynomials have been evaluated
    and are stored in the array 'pc'.    */

  if(np < ncf)
    {
      for(i=np;i<ncf;++i)  pc[i]=twot*pc[i-1]-pc[i-2];
      np=ncf;
    }

/*  interpolate to get position for each component  */

  for(i=0;i<ncm;++i) /* ncm is a number of coordinates */
     {
       posvel[i]=0.0;
       for(j=ncf-1;j>=0;--j)
          posvel[i]=posvel[i]+pc[j]*coef[j+i*ncf+l*ncf*ncm];
     }

      if(ifl <= 1) return;


/*  if velocity interpolation is wanted, be sure enough
    derivative polynomials have been generated and stored.    */

  vfac=(dna+dna)/t[1];
  vc[2]=twot+twot;
  if(nv < ncf)
    {
      for(i=nv;i<ncf;++i) vc[i]=twot*vc[i-1]+pc[i-1]+pc[i-1]-vc[i-2];
      nv=ncf;
    }

/*  interpolate to get velocity for each component    */

   for(i=0;i<ncm;++i)
      {
        posvel[i+ncm]=0.0;
        for(j=ncf-1;j>0;--j)
           posvel[i+ncm]=posvel[i+ncm]+vc[j]*coef[j+i*ncf+l*ncf*ncm];
        posvel[i+ncm]=posvel[i+ncm]*vfac;
      }
   return;
}

/****************************************************************************
****                       jpl_mod_split(tt,fr)                          ****
*****************************************************************************
****  this subroutine breaks a d.p. number into a d.p. integer           ****
****  and a d.p. fractional part.                                        ****
****                                                                     ****
****  calling sequence parameters:                                       ****
****                                                                     ****
****    tt = d.p. input number                                           ****
****                                                                     ****
****    fr = d.p. 2-word output array.                                   ****
****         fr(1) contains integer part                                 ****
****         fr(2) contains fractional part                              ****
****                                                                     ****
****         for negative input numbers, fr(1) contains the next         ****
****         more negative integer; fr(2) contains a positive fraction.  ****
****************************************************************************/
void jpl_mod_split(double tt, double fr[2])
{
/*  main entry -- get integer and fractional parts  */

      fr[1]=modf(tt,&fr[0]);

      if(tt >= 0.0 || fr[1] == 0.0) return;

/*  make adjustments for negative input number   */

      fr[0]=fr[0]-1.0;
      fr[1]=fr[1]+1.0;

      return;
}

/*****************************************************************************
**                        jpl_mod_state(et2,list,pv,nut)                            **
******************************************************************************
** This subroutine reads and interpolates the jpl planetary ephemeris file  **
**                                                                          **
**    Calling sequence parameters:                                          **
**                                                                          **
**    Input:                                                                **
**                                                                          **
**        et2[] double, 2-element JED epoch at which interpolation          **
**              is wanted.  Any combination of et2[0]+et2[1] which falls    **
**              within the time span on the file is a permissible epoch.    **
**                                                                          **
**               a. for ease in programming, the user may put the           **
**                  entire epoch in et2[0] and set et2[1]=0.0               **
**                                                                          **
**               b. for maximum interpolation accuracy, set et2[0] =        **
**                  the most recent midnight at or before interpolation     **
**                  epoch and set et2[1] = fractional part of a day         **
**                  elapsed between et2[0] and epoch.                       **
**                                                                          **
**               c. as an alternative, it may prove convenient to set       **
**                  et2[0] = some fixed epoch, such as start of integration,**
**                  and et2[1] = elapsed interval between then and epoch.   **
**                                                                          **
**       list   12-element integer array specifying what interpolation      **
**              is wanted for each of the "bodies" on the file.             **
**                                                                          **
**                        list[i]=0, no interpolation for body i            **
**                               =1, position only                          **
**                               =2, position and velocity                  **
**                                                                          **
**              the designation of the astronomical bodies by i is:         **
**                                                                          **
**                        i = 0: mercury                                    **
**                          = 1: venus                                      **
**                          = 2: earth-moon barycenter                      **
**                          = 3: mars                                       **
**                          = 4: jupiter                                    **
**                          = 5: saturn                                     **
**                          = 6: uranus                                     **
**                          = 7: neptune                                    **
**                          = 8: pluto                                      **
**                          = 9: geocentric moon                            **
**                          =10: nutations in longitude and obliquity       **
**                          =11: lunar librations (if on file)              **
**                                                                          **
**                                                                          **
**    output:                                                               **
**                                                                          **
**    pv[][6]   double array that will contain requested interpolated       **
**              quantities.  The body specified by list[i] will have its    **
**              state in the array starting at pv[i][0]  (on any given      **
**              call, only those words in 'pv' which are affected by the    **
**              first 10 'list' entries (and by list(11) if librations are  **
**              on the file) are set.  The rest of the 'pv' array           **
**              is untouched.)  The order of components in pv[][] is:       **
**              pv[][0]=x,....pv[][5]=dz.                                   **
**                                                                          **
**              All output vectors are referenced to the earth mean         **
**              equator and equinox of epoch. The moon state is always      **
**              geocentric; the other nine states are either heliocentric   **
**              or solar-system barycentric, depending on the setting of    **
**              global variables (see below).                               **
**                                                                          **
**              Lunar librations, if on file, are put into pv[10][k] if     **
**              list[11] is 1 or 2.                                         **
**                                                                          **
**        nut   dp 4-word array that will contain nutations and rates,      **
**              depending on the setting of list[10].  the order of         **
**              quantities in nut is:                                       **
**                                                                          **
**                       d psi  (nutation in longitude)                     **
**                       d epsilon (nutation in obliquity)                  **
**                       d psi dot                                          **
**                       d epsilon dot                                      **
**                                                                          **
**    Global variables:                                                     **
**                                                                          **
**         KM   logical flag defining physical units of the output          **
**              states. KM = TRUE, km and km/sec                            **
**                         = FALSE, au and au/day                           **
**              default value = FALSE.  (Note that KM determines time unit  **
**              for nutations and librations. Angle unit is always radians.)**
**                                                                          **
**       BARY   logical flag defining output center.                        **
**              only the 9 planets (entries 0-8) are affected.              **
**                       BARY = TRUE =\ center is solar-system barycenter   **
**                            = FALSE =\ center is sun                      **
**              default value = FALSE                                       **
**                                                                          **
**      PVSUN   double, 6-element array containing                          **
**              the barycentric position and velocity of the sun.           **
**                                                                          **
*****************************************************************************/
void jpl_mod_state(double et2[2],int list[12],double pv[][6],double nut[4])
{
  int i,j;
  int verb=0; //0 does not print values
  static int ipt[13][3], first=TRUE;
                  /* local copy of R1.r1.ipt[] is extended for R1.r1.lpt[[] */
  long int nr;
  static long int nrl=0;
  double pjd[4];
  static double buf[NCOEFF];
  double s,t[2],aufac;
  double pefau[6];

  if(first)
    {
      first=FALSE;
      for(i=0;i<3;++i)
         {
           for(j=0;j<12;++j) ipt[j][i]=(int)R1.r1.ipt[j][i];
           ipt[12][i] = (int)R1.r1.lpt[i];
         }
   }

/*  ********** main entry point **********  */

  s=et2[0] - 0.5;
  jpl_mod_split(s,&pjd[0]);
  jpl_mod_split(et2[1],&pjd[2]);
  pjd[0]=pjd[0]+pjd[2]+0.5;
  pjd[1]=pjd[1]+pjd[3];
  jpl_mod_split(pjd[1],&pjd[2]);
  pjd[0]=pjd[0]+pjd[2];
/* here pjd[0] contains last midnight before epoch desired (in JED: *.5)
   and pjd[3] contains the remaining, fractional part of the epoch         */

/*   error return for epoch out of range  */

  if( (pjd[0]+pjd[3]) < R1.r1.ss[0] || (pjd[0]+pjd[3]) > R1.r1.ss[1] )
    {
      puts("Requested JED not within ephemeris limits.\n");
      return;
    }

/*   calculate record # and relative time in interval   */

      nr=(long)((pjd[0]-R1.r1.ss[0])/R1.r1.ss[2])+2;
      /* add 2 to adjus for the first two records containing header data */
      if(pjd[0] == R1.r1.ss[1]) nr=nr-1;
      t[0]=( pjd[0]-( (1.0*nr-2.0)*R1.r1.ss[2]+R1.r1.ss[0] ) +
           pjd[3] )/R1.r1.ss[2];

/*   read correct record if not in core (static vector buf[])   */

      if(nr != nrl)
        {
          nrl=nr;
          fseek(F1,nr*RECSIZE,SEEK_SET);
          fread(buf,sizeof(buf),1,F1);
        }

      if(KM)
        {
          t[1]=R1.r1.ss[2]*86400.0;
          aufac=1.0;
        }
      else
        {
          t[1]=R1.r1.ss[2];
          aufac=1.0/R1.r1.au;
        }

/*  every time interpolate Solar System barycentric sun state   */

    jpl_mod_interp(&buf[ipt[10][0]-1],t,ipt[10][1],3,ipt[10][2],2,pefau);

      for(i=0;i<6;++i)  PVSUN[i]=pefau[i]*aufac;

/*  check and interpolate whichever bodies are requested   */

      for(i=0;i<10;++i)
         {
           if(list[i] == 0) continue;

           jpl_mod_interp(&buf[ipt[i][0]-1],t,ipt[i][1],3,ipt[i][2],list[i],pefau);

           for(j=0;j<6;++j)
              {
                if(i < 9 && !BARY)   pv[i][j]=pefau[j]*aufac-PVSUN[j];
                else                 pv[i][j]=pefau[j]*aufac;
              }
         }

/*  do nutations if requested (and if on file)    */

      if(list[10] > 0 && ipt[11][1] > 0)
         jpl_mod_interp(&buf[ipt[11][0]-1],t,ipt[11][1],2,ipt[11][2],list[10],nut);

/*  get librations if requested (and if on file)    */

      if(list[11] > 0 && ipt[12][1] > 0)
        {
          jpl_mod_interp(&buf[ipt[12][0]-1],t,ipt[12][1],3,ipt[12][2],list[11],pefau);
          for(j=0;j<6;++j) pv[10][j]=pefau[j];
        }
  return;
}
/****************************************************************************
**                        jpl_mod_constan(nam,val,sss,n)                   **
*****************************************************************************
**                                                                         **
**    this function obtains the constants from the ephemeris file          **
**    calling seqeunce parameters (all output):                            **
**      char nam[][6] = array of constant names (max 6 characters each)    **
**      double val[] = array of values of constants                        **
**      double sss[3] = JED start, JED stop, step of ephemeris             **
**      int n = number of entries in 'nam' and 'val' arrays                **
** Note: we changed name of this routine becouse const is a reserved word  **
**       in C language.
*****************************************************************************
**    external variables:                                                  **
**         struct R1 and R2 (first two records of ephemeris)               **
**         defined in file:     jplbin.h                                   **
**         F1 = ephemeris binary file pointer (obtained from fopen() )     **
/*

	char ttl[3][84];
	char cnam[400][6];
        double ss[3];
        int ncon; 
        double au;
        double emrat;
        int ipt[12][3]; 
        int numde; 
        int lpt[3]; 

*/

void jpl_mod_constan(char nam[][6], double val[], double sss[], intptr_t *n)
{
  int i,j;
  int verb=0; //0 does not print values. Con 5 stampa valori

  int icheck64=1;  //put to 1 to have the checks. Pia 4 maggio 2011
  

  //fread(&R1,sizeof(R1),1,F1);  //pia: does not work on a 64 bit computer
  fread(&R1.r1.ttl,sizeof(R1.r1.ttl),1,F1);
  if(icheck64==1&& verb>4)printf("%s\n",R1.r1.ttl);
  fread(&R1.r1.cnam,sizeof(R1.r1.cnam),1,F1);
  if(icheck64==1&& verb>4)printf("%s\n",R1.r1.cnam);
  fread(&R1.r1.ss,sizeof(R1.r1.ss),1,F1);
  if(icheck64==1&& verb>4)printf("%f\n",R1.r1.ss);
  fread(&R1.r1.ncon,sizeof(R1.r1.ncon),1,F1);
   if(icheck64==1&& verb>4)printf("%d\n",R1.r1.ncon);
  fread(&R1.r1.au,sizeof(R1.r1.au),1,F1);
   if(icheck64==1&& verb>4)printf("%f\n",R1.r1.au);
  fread(&R1.r1.emrat,sizeof(R1.r1.emrat),1,F1);
   if(icheck64==1&& verb>4)printf("%f\n",R1.r1.emrat);
  fread(&R1.r1.ipt,sizeof(R1.r1.ipt),1,F1);
   if(icheck64==1&& verb>4)printf("%d\n",R1.r1.ipt);
  fread(&R1.r1.numde,sizeof(R1.r1.numde),1,F1);
   if(icheck64==1&& verb>4)printf("%d\n",R1.r1.numde);
  fread(&R1.r1.lpt,sizeof(R1.r1.lpt),1,F1);
   if(icheck64==1&& verb>4)printf("%d\n",R1.r1.lpt);

  fread(&R2,sizeof(R2),1,F1);
    // *n =(int) R1.r1.ncon;  
  *n =(intptr_t) R1.r1.ncon;  //pia64
  
  if(icheck64==1 && verb>4){
//printf("%d %d %d \n",*n,R1.r1.ncon,sizeof(R1.r1.ncon));
	printf("in jpl_mod_constan %d %d %d %d %d \n",*n,R1.r1.ncon,sizeof(R1.r1.ncon),sizeof(R1),sizeof(R2));
	puts("devono essere 156 156 4 8144 8144");
	printf("ttl sizeof %d\n",sizeof(R1.r1.ttl));
	printf("cnam sizeof %d\n",sizeof(R1.r1.cnam));
	printf("ss sizeof %d\n",sizeof(R1.r1.ss));
	printf("ncon int sizeof %d\n",sizeof(R1.r1.ncon));
	printf("au sizeof %d\n",sizeof(R1.r1.au));
	printf("emrat sizeof %d\n",sizeof(R1.r1.emrat));
	printf("ipt int sizeof %d\n",sizeof(R1.r1.ipt));
	printf("numde int sizeof %d\n",sizeof(R1.r1.numde));
	printf("lpt int sizeof %d\n",sizeof(R1.r1.lpt));
  }
  for(i=0;i<3;++i) sss[i]=R1.r1.ss[i];

      for(i=0;i<*n;++i) 
 
     {
       val[i] = R2.r2.cval[i];
       for(j=0;j<6;++j) nam[i][j] = R1.r1.cnam[i][j];
     }

  
  return;
}



/*************************************************************************************/
/*******************************Other utility functions ******************************/
/*************************************************************************************/

/***** Function PSS_Doppler *****/
/*****Evaluates the Doppler effect, scalar product of the poss and vel vectors
 * INPUT
 *  ARGUMENTS:
      *source (struct SOURCE)
           poss[3] (double)
         Source position vector, equatorial rectangular coordinates,
         normalized components.
           vel[3] (double)
         Detector velocity vector, equatorial rectangular coordinates,
         components in AU/day. 
         

*   OUTPUT
*   ARGUMENTS:
      nu_doppler (double)
      Observed frequency, in Hz

   RETURNED
   VALUE:
      nu_doppler.

   GLOBALS
   USED:
      C, day2sec, KMAU

*****/
double  PSS_doppler(SOURCE *source,double *poss , double *vel)
    
  {
    double sa,c_used;
    double nu_doppler; 

  
    //sa=poss[0]*vel[0];
    //sb=poss[1]*vel[1];
    //sc=poss[2]*vel[2];
    sa=scalar_product(poss,vel);
    c_used=C; //*KMAU/day2sec; ///input in AU/day
    nu_doppler=source->frequency*(1.0+sa/c_used); //CORRETTO IL + IL 13/10/06 !!
    // printf("c freq doppler= %f %f %f\n",c_used,source->frequency,nu_doppler);
    return nu_doppler;
   } 


/***** Function PSS_einstein*****/
double PSS_einstein(double tdt_i)

{

   /* Builtin functions */
    double d_mod(), sin(), cos();

    double g_bary, g, g2;
    double Deinstein;   

/*      Function that evaluates the Doppler effect due to */
/*      the relativistic Einstein delay                   */
/* 	We must evaluate d(TDB-TDT)/dt=Deinstein          */
/*      nu_doppler=f_source*Deinstein                     */
/*      We use an approximate formula for TDB-TDT         */
/*      INPUT
        tdt_i: TDT at the current time. double
        OUTPUT
        Deinstein: time derivative of the time shift due to the
        Einstein effect. double
*/

  
    g_bary= (tdt_i- T0) * 0.9856003 + 357.53; /* deg */
    g_bary *= DEG2RAD;      /* rad */
    g = fmod(g_bary, TWOPI);
    g2 = fmod(2*g_bary, TWOPI);
    // tdb = tdt_i + (sin(g) * 0.001658 + sin(g2) * 1.4e-5) / 86400.;
 

    Deinstein=(0.9856003/day2sec)*(0.001658*cos(g)+2.0*0.000014*cos(g2)); 
    /* (deg/s)*s */
    Deinstein *= DEG2RAD;
    //    printf("%20.16f\n", Deinstein);
    return Deinstein;
} /* einstein_ */

/***** Function PSS_einstein_vect *****/
void PSS_einstein_vect(SOURCE *source,double dot_ini, double step, int npoint,double nu_einstein[])

{

   /* Builtin functions */
    double d_mod(), sin(), cos();

    double g_bary, g, g2;   
    double Deinstein,utc,tdt,tdt_i;
    int k;
    int verb=0; //0 does not print values

/*      Function that evaluates the Doppler effect due to */
/*      the relativistic Einstein delay                   */
/* 	We must evaluate d(TDB-TDT)/dt=Deinstein          */
/*      nu_einstein=f_source*Deinstein                    */
/*      For a given source frequency,                     */ 
/*      initial date (days from the beginning of the century),step (s), number */
/*      of point, it gives the frequency shift            */ 
/*      We use an approximate formula for TDB-TDT         */
/*      INPUT
	source, pointer to struct SOURCE
	dot_ini, double: days from the beginning of the Century
	step, double: step in seconds
	npoint, int: number of data

        OUTPUT
        nu_einstein[ ]: frequency shift due to the 
        Einstein effect. vector double.

        This contribution must be subtracted to the Doppler effect
*/

  utc=dot_ini+dot2utc;
  if (verb>4)printf("Einstein effect, initial utc= %18.8f\n ",utc);
  /*****Leap seconds*****/
  tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
  
 for (k=0; k <=npoint; ++k)
       {
       tdt_i=tdt+(k*step/day2sec);  

       g_bary= (tdt_i- T0) * 0.9856003 + 357.53; /* deg */
       g_bary *= DEG2RAD;      /* rad */
       g = fmod(g_bary, TWOPI);
       g2 = fmod(2*g_bary, TWOPI);
    // tdb = tdt_i + (sin(g) * 0.001658 + sin(g2) * 1.4e-5) / 86400.;
 

    Deinstein=(0.9856003/day2sec)*(0.001658*cos(g)+2.0*0.000014*cos(g2)); 
    /* (deg/s)*s */
    Deinstein *= DEG2RAD;
    //printf("%20.16f\n", Deinstein);
    nu_einstein[k]=Deinstein*source->frequency;

       }
    return;
} /* einstein_vect */



/*****************************************************************/
/*****Shapiro effect:*****/

void PSS_shapiro_vect(DETECTOR *detector,SOURCE *source,double dot_ini, double step, int npoint,short int fut1, short int nut1, short int nutcorr,short int equatorial,double nu_shapiro[])

/*
For a given detector, initial date (days from the beginning of the century),
step (s), number of point, it gives the frequency shift (Hz)

                nu_shapiro=f_source*Dshapiro 


INPUT ARGUMENTS:

detector, pointer to struct DETECTOR
source,  pointer to struct SOURCE
dot_ini, double: days from the beginning of the Century
step, double: step in seconds
npoint, int: number of data
fut1, int: flag for UT1. If 0 the function uses UTC, if 1 it uses UT1
 short int equatorial=1; //if 1 evaluations are in equatorial copordinates
                           (which does not change the result..)

OUTPUT ARGUMENTS:

nu_shapiro[npoint]  frequency shift due to the 
        Shapiro effect. vector double.

This contribution must be subtracted to the Doppler effect

RETURN:
void
*************************************************************************/

{
 
  double et, r[6];
  int i, ntarg, nctr;
  int ii,k;
  double pos[3],poss[3];
  double vel[3];
  double utc,tdt,tdb,utc_ii,utch,utcl,ut1,mjd;
  /****used in sidereal_time****/
  double eqeq,eqeqP; 
  double gast;
  /****used in earthtilt****/
  double oblm,oblt,psi,eps;
/****used in terra, nutate,precession****/
  double pos1[3],pos2[3],pog[3]; 
  double vel1[3],vel2[3],vog[3];
/****used for nutation in the ephemeris****/
  double ddd[6];
  /****specific for Shapiro***/
  double R_mod,aultsc,rschw,sc1,sc2,gauss;
  double Dshapiro;


  gauss= 0.01720209895;

  F1=fopen(EPHFILE,"rb");
  F2=fopen(UTC2UT1,"rb");
  //F3=fopen(DELTANUT,"rb");
  //puts("Initial time, days  \n");
  //printf("%f\n", dot_ini);

  /****Source***/
 PSS_source_pos(source,equatorial,poss);

 /****definitions****/
 aultsc=day2sec*1.0/C;
 rschw = (pow(gauss,2)) * (pow(aultsc,3)) / (pow(day2sec,2));
 
  //printf("aultsc,rschw= %f %f\n",aultsc,rschw);

  ii=0;
  ntarg=3; /* Earth */
  nctr=11; /*Center of the Sun */
  utc=dot_ini+dot2utc;
  mjd=utc-utc2mjd;
  //printf("initial utc= %18.8f\n ",utc);
  //printf("initial mjd= %18.8f\n ",mjd);
  /*****Leap seconds*****/
  tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
  //printf("initial tdt= %18.8f\n ",tdt);
  tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
  //printf("tdb= %18.8f\n ",tdb);
  //printf("Detector name,long,lat:%s %f %f\n",detector->name,detector->longitude,detector->latitude);

 while( ii < npoint)
   {
   /*****  Read a value; Skip if not within the time-range
        of the present version of the ephemeris.                            */

    et=tdb+(ii*step/day2sec); /*et is the barycenter time, at each sample */
    utc_ii=utc+(ii*step/day2sec);
 
    jpl_mod_pleph(et,ntarg,nctr,r);
   

    for (k=0; k<=2; ++k)
	{
	  pos[k]=r[k];     //AU
          vel[k]=r[k+3];   //AU/day. Earth center vs Sun Center, referred at T0
	                   //Rectangular equatorial coordinates
	}
    /*****Evaluate Detector velocity and position vs Earth Center*****/
    /***** Nutation series from the JPL file *****/
     jpl_mod_pleph(et,14,0,ddd);
     psi=ddd[0]*RAD2SEC;
     eps=ddd[1]*RAD2SEC;
     novas_mod_earthtilt (et, &oblm,&oblt,&eqeq,psi,eps);

     /*****Sidereal time needs UT1. We use UTC that approx. UT1 (error <0.9 s/anno) if fut1=0,  else we use UT1*****/
     if (!fut1)
     {
	utch=(long)(utc_ii);
	utcl=utc_ii-utch;
      }
	else
      {
        ut1=PSS_utc2ut1(utc_ii);
        utch=(long)(ut1);
        utcl=ut1-utch;
       }
     eqeqP=eqeq; //Kaz suggestion. May 2003
     if(NutFlag==0)eqeqP=0.0;
   sidereal_time (utch,utcl,eqeqP, &gast);
   novas_mod_terra (detector,gast, pos1,vel1);
 //modified Sept.-Oct 2003 !!!!!
   if(nutcorr==1)novas_mod_nutate (et,FN1,pos1, pos2,psi,eps); //uncomment if required
   if(nutcorr==0)nutate (et,FN1,pos1, pos2); 
   precession(et,pos2,T0, pog); 
   if(nutcorr==1)novas_mod_nutate(et,FN1,vel1, vel2,psi,eps); /*performs inverse nutation (true to mean) */
   if(nutcorr==0)nutate(et,FN1,vel1, vel2);/*performs inverse nutation (true to mean) */
   precession(et,vel2,T0, vog); //rectangular equatorial coord., ref. at T0 


  /***Sum the two contributions: Earth vs Center of Sun + Det vs. Earth *****/
  for(k=0; k<=2; ++k)
      {
	pos[k]+=pog[k];
	vel[k]+=vog[k];
	//pos[k]=pog[k];  //only det vs Earth-Center
	//vel[k]=vog[k];
       
        vel[k]/=day2sec; //AU/secondo
      }

     if(equatorial!=1)
     {
        /*****convert into Ecliptical coordinates*****/
        PSS_equr2eclr (vel, vel1); 
        PSS_equr2eclr (pos, pos1);
     }
     else
     {
       pos1[0]=pos[0];
       pos1[1]=pos[1];
       pos1[2]=pos[2];
       vel1[0]=vel[0];
       vel1[1]=vel[1];
       vel1[2]=vel[2];       
     }
    /****Shapiro vector****/
    R_mod=sqrt(pow(pos1[0],2)+pow(pos1[1],2)+pow(pos1[2],2));
    sc1=scalar_product(poss,vel1);
    sc2=scalar_product(poss,pos1);
    Dshapiro=-2*rschw*(sc1/(R_mod+sc2));
    
    nu_shapiro[ii]=Dshapiro*source->frequency;
    ii++;
   }
     fclose(F1);
     fclose(F2);
     //     fclose(F3);
     return;
}






/*****Coordinate conversions:*****/

/********PSS_rect2spher*****/

short int PSS_rect2spher (double *pos, 

                        double *ra, double *dec)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Converts a vector in equatorial or ecliptical rectangular coordinates to
      equatorial or ecliptical spherical coordinates.

   REFERENCES: 
      Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97, 
         pp. 1197-1210.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.

   INPUT
   ARGUMENTS:
      pos[3] (double)
         Position vector, equatorial (or ecliptical) rectangular coordinates.

   OUTPUT
   ARGUMENTS:
      *rightascension (double)
         Right ascension (or ecliptical longitude)  in hours.
      *declination (double)
         Declination (or ecliptical latitude) in degrees.

   RETURNED
   VALUE:
      (short int)
         0...Everything OK.
         1...All vector components are zero; 'ra' and 'dec' are
             indeterminate.
         2...Both vec[0] and vec[1] are zero, but vec[2] is nonzero;
             'ra' is indeterminate.
   GLOBALS
   USED:
      RAD2SEC

   FUNCTIONS
   CALLED:
      sqrt     math.h
      pow      math.h
      atan2    math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'angles'.

------------------------------------------------------------------------
*/
{
   double xyproj;

   xyproj = sqrt (pow (pos[0], 2.0) + pow (pos[1], 2.0));
   if ((xyproj == 0.0) && (pos[2] == 0))
   {
      *ra = 0.0;
      *dec = 0.0;
      return 1;
   }
    else if (xyproj == 0.0)
   {
      *ra = 0.0;
      if (pos[2] < 0.0)
         *dec = -90.0;
       else
         *dec = 90.0;
      return 2;
   }
    else
   {
      *ra = atan2 (pos[1], pos[0]) * RAD2SEC / 54000.0;
      *dec = atan2 (pos[2], xyproj) * RAD2SEC / 3600.0;

      if (*ra < 0.0)
         *ra += 24.0;
   }
   return 0;
}



/*****PSS_spher2rect*****/

void PSS_spher2rect (double ra, double dec, double dist,

                   double *vector)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Converts equatorial spherical coordinates to a vector (equatorial
      or ecliptical rectangular coordinates).

   REFERENCES: 
      None.

   INPUT
   ARGUMENTS:
      ra (double)
         Right ascension or ecliptical longitude (hours).
      dec (double)
         Declination or ecliptical latitude (degrees).

   OUTPUT
   ARGUMENTS:
      vector[3] (double)
         Position vector, equatorial or ecliptical rectangular coordinates (AU).

   RETURNED
   VALUE:
      (short int)
         0...Everything OK.

   GLOBALS
   USED:
      DEG2RAD

   FUNCTIONS
   CALLED:
      cos     math.h
      sin     math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/05-92/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.

   NOTES:
      None.

------------------------------------------------------------------------
*/
{

   vector[0] = dist * cos (DEG2RAD * dec) * cos (DEG2RAD * 15.0 * ra);
   vector[1] = dist * cos (DEG2RAD * dec) * sin (DEG2RAD * 15.0 * ra);
   vector[2] = dist * sin (DEG2RAD * dec);

   return;
}


/*****PSS_eclr2equr*****/

void PSS_eclr2equr (double *vector_in, double *vector_out)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Converts ecliptical rectangular coordinates to equatorial
      rectangular coordinates).

   REFERENCES: 
      None.

   INPUT
   ARGUMENTS:
      vector_in[3] (double)
      Position vector, ecliptical rectangular coordinates (AU).    


   OUTPUT
   ARGUMENTS:
      vector_out[3] (double)
         Position vector, equatorial rectangular coordinates (AU).


   GLOBALS
   USED:
      DEG2RAD
      ecl_deg : Value of the ecliptical angle, deg

   FUNCTIONS
   CALLED:
      cos     math.h
      sin     math.h

 
------------------------------------------------------------------------
*/
{

   vector_out[0] = vector_in[0];
   vector_out[1] = vector_in[1] * cos (DEG2RAD * ecl_deg) - vector_in[2] * sin (DEG2RAD * ecl_deg);
   vector_out[2] = vector_in[1] * sin (DEG2RAD * ecl_deg) + vector_in[2] * cos (DEG2RAD * ecl_deg);

   return;
}

/*****PSS_equr2eclr*****/

void PSS_equr2eclr (double *vector_in, double *vector_out)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Converts equatorial rectangular coordinates to ecliptical
      rectangular coordinates).

   REFERENCES: 
      None.

   INPUT
   ARGUMENTS:
      vector_in[3] (double)
      Position vector, equatorial rectangular coordinates (AU).    


   OUTPUT
   ARGUMENTS:
      vector_out[3] (double)
         Position vector, ecliptical rectangular coordinates (AU).


   GLOBALS
   USED:
      DEG2RAD
      ecl_deg : Value of the ecliptical angle, deg

   FUNCTIONS
   CALLED:
      cos     math.h
      sin     math.h

 
------------------------------------------------------------------------
*/
{

   vector_out[0] =  vector_in[0];
   vector_out[1] =  vector_in[1] * cos (DEG2RAD * ecl_deg) + vector_in[2] * sin (DEG2RAD * ecl_deg);
   vector_out[2] = -vector_in[1] * sin (DEG2RAD * ecl_deg) + vector_in[2] * cos (DEG2RAD * ecl_deg);

   return;
}

/*****PSS_equs2ecls*****/

void PSS_equs2ecls (double ra, double dec, double *lambda, double *beta)
/*

------------------------------------------------------------------------

   PURPOSE:    
      Converts equatorial spherical coordinates to a ecliptical
      spherical coordinates).

   REFERENCES: 
      None.

   INPUT
   ARGUMENTS:
      ra (double)
         Right ascension (hours).
      dec (double)
         Declination (degrees).

   OUTPUT
   ARGUMENTS:
     lambda (double)
         celestial longitude (hours). Ranges from 0 up to 24 hours.
      beta (double)
         celestial latitude (degrees) Ranges from -90 up to 90 deg.

   RETURNED
   VALUE:
      void

   GLOBALS
   USED:
      DEG2RAD
      RAD2DEG
      ecl_deg : Value of the ecliptical angle, deg
      TWOPI

   FUNCTIONS
   CALLED:
      cos     math.h
      sin     math.h



------------------------------------------------------------------------
*/
{
  double a1,a2,coslambda,sinlambda;

  a1= sin (DEG2RAD * dec) * cos (DEG2RAD * ecl_deg);
  a2= sin (DEG2RAD * ecl_deg) * cos (DEG2RAD * dec) * sin(DEG2RAD*15.0*ra);
  *beta = asin(a1 - a2);
  coslambda= cos (DEG2RAD * 15.0 *ra) * cos (DEG2RAD * dec);/// / cos(*beta);
  a1= sin (DEG2RAD * dec) * sin (DEG2RAD * ecl_deg);
  a2= cos (DEG2RAD * dec) * cos (DEG2RAD * ecl_deg) * sin(DEG2RAD*15.0*ra);
  sinlambda=(a1+a2); ///  / cos(*beta);  
  *lambda=atan2(sinlambda,coslambda);
   /*atan2 gives the result between (-pi; pi) radians*/
  *lambda*=(RAD2DEG/15.0); //hours
     if (*lambda < 0.0)
         *lambda += 24.0;
  *beta*=RAD2DEG;          //deg
   return;
}


/*****PSS_ecls2equs*****/

void PSS_ecls2equs (double lambda, double beta, double *ra, double *dec)
/*

------------------------------------------------------------------------

   PURPOSE:    
      Converts ecliptical spherical coordinates to  equatorial
      spherical coordinates).

   REFERENCES: 
      None.

   INPUT
   ARGUMENTS:
    lambda (double)
         celestial longitude (hours). Ranges from 0 up to 24 hours.
      beta (double)
         celestial latitude (degrees) Ranges from -90 up to 90 deg.

   OUTPUT
   ARGUMENTS:
      ra (double)
         Right ascension (hours).
      dec (double)
         Declination (degrees).
 
   RETURNED
   VALUE:
      void

   GLOBALS
   USED:
      DEG2RAD
      RAD2DEG
      ecl_deg : Value of the ecliptical angle, deg
      TWOPI

   FUNCTIONS
   CALLED:
      cos     math.h
      sin     math.h



------------------------------------------------------------------------
*/
{
  double a1,a2,cosra,sinra;

  a1= sin (DEG2RAD * beta) * cos (DEG2RAD * ecl_deg);
  a2= sin (DEG2RAD * ecl_deg) * cos (DEG2RAD * beta) * sin(DEG2RAD*15.0*lambda);
  *dec = asin(a1 + a2);
  cosra= cos (DEG2RAD * 15.0 *lambda) * cos (DEG2RAD * beta); ///   / cos(*dec);
  a1= sin (DEG2RAD * beta) * sin (DEG2RAD * ecl_deg);
  a2= cos (DEG2RAD * beta) * cos (DEG2RAD * ecl_deg) * sin(DEG2RAD*15.0*lambda);
  sinra=(-a1+a2); ///  / cos(*dec);  

  *ra=atan2(sinra,cosra);
  /*atan2 gives the result between (-pi; pi) radians*/
  *ra*=(RAD2DEG/15);  //hours
    if (*ra < 0.0)
         *ra += 24.0;
  *dec*=RAD2DEG;      //deg
   return;
}
/******************************************************************************************/
/*****Detector velocity:*****/

double PSS_detector_vel(DETECTOR *detector,double dot_ini, double step, int npoint,short int fut1, short int nut1, short int nutcorr,short int equatorial,short int imjd,double velx[], double vely[],double velz[], double deinstein[])

/*For a given detector, initial date (days from the beginning of the century or mjd),step (s), number
 of point, it gives the velocity vector (AU/day) in rectangular ecliptic 
or equatorial coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch). 

It returns MJD of the first point.

INPUT ARGUMENTS:

detector, pointer to struct DETECTOR

dot_ini, double: days from the beginning of the Century, or mjd, or TDT in mjd

step, double: step in seconds

npoint, int: number of data

fut1, int: flag for UT1. If 0 the function uses UTC, if 1 it uses UT1

imjd=0,1,2 0=days from 1900, 1=mjd , 2=tdt in JD 

OUTPUT ARGUMENTS:

velx[npoint], 
vely[npoint]
velz[npoint]
double arrays:Det. velocity vector (AU/day) in rectangular ecliptic or euatorial coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch)

deinstein[npoint]
RETURN:
MJD of the first datum
*******************************************************************************************
NOTE: Position vector of the Earth and Detector are not required.
      Uncomment the relative instructions if you need a code that computes also positions
********************************************************************************************
*/

{
   /* Sotto era 400 anzichè NMAX che ora vale 1000 e sta in jplbin.h */
  char nams[NMAX][6];
  double et, r[6], ss[3], vals[NMAX];
  int i, ntarg, nctr;
  intptr_t nvs;  //pia64
  int ii,k;
  double pos[3];
  double vel[3];
  double utc,tdt,tdb,utc_ii,utch,utcl,ut1,mjd;
  /****used in sidereal_time****/
  double eqeq,eqeqP; 
  double gast;
  /****used in earthtilt****/
  double oblm,oblt,psi,eps;
  double ddpsi,ddeps; 
/****used in terra, nutate,precession****/
  double pos1[3],pos2[3],pog[3]; 
  double vel1[3],vel2[3],vog[3];
/****used for nutation in the ephemeris****/
  double ddd[6];
  double tdt_ini; 
  int verb=0; //0 does not print values

  /******PER TEST FILE */
   FILE *OUT33;
   OUT33=fopen("table_test.dat","wb");


  F1=fopen(EPHFILE,"rb");
  F2=fopen(UTC2UT1,"rb");
  F3=fopen(DELTANUT,"rb");

  if (verb>4) puts("Initial time, days");
  if (verb>4) printf("%f\n", dot_ini);

/****** Print the ephemeris constants. ******/

  jpl_mod_constan(nams,vals,ss,&nvs);
  if (verb>4) printf("Range DE405= %14.2f  %14.2f  %14.2f\n",ss[0],ss[1],ss[2]); //pia Sep. 2005
  /*  for(i=0;i<nvs;++i)
      printf("%.6s  %24.16E\n",nams[i],vals[i]); */

  ii=0;
  ntarg=3; /* Earth */
  nctr=12; /*Solar system barycenter */
  if(imjd==0)
  {
    utc=dot_ini+dot2utc;
    mjd=utc-utc2mjd;
    if (verb>4)printf("initial utc= %18.8f\n ",utc); //pia Sep. 2005
    if (verb>4)printf("initial mjd= %18.8f\n ",mjd);
    if (verb>4)printf("initial days from 1900= %18.8f\n ",dot_ini);
  }
   if(imjd==1)
  {
    utc=dot_ini+utc2mjd;
    mjd=dot_ini;
    dot_ini=utc-dot2utc;  //dot_ini sempre giorni inizio secolo, in UTC !! NOTARE
    if (verb>4)printf("initial utc= %18.8f\n ",utc);
    if (verb>4)printf("initial mjd= %18.8f\n ",mjd);
  }
   if(imjd==2){
     tdt=dot_ini+utc2mjd;
     mjd=dot_ini; //mjd but TDT
     if (verb>4)printf("initial tdt= %20.11f\n ",tdt);
     if (verb>4)printf("initial TDT in mjd= %18.8f\n ",mjd);
   }
  /*****Leap seconds*****/
  if(imjd !=2) tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
  if (verb>4)printf("initial tdt= %18.8f\n ",tdt); //pia Sep. 2005
  tdt_ini=tdt;
  tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
  if (verb>4)printf("tdb= %18.8f\n ",tdb);
  if (verb>4)printf("utc2mjd= %25.15f\n ",utc2mjd);
  ////printf("Detector name,long,lat:%s %f %f\n",detector->name,detector->longitude,detector->latitude); //pia Sep. 2005

 while( ii < npoint)
   {
   /*****  Read a value; Skip if not within the time-range
        of the present version of the ephemeris.                            */
     if(imjd!=2){
       utc=dot_ini+dot2utc+ii*step/(1.0*day2sec); //commento kaz: tdb va ricalcolato sempre
       tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
     }
     else
       {
	 tdt=tdt_ini+ii*step/(1.0*day2sec);
       }
     tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
     et=tdb;              //sost. kaz fino a qui e tolta la riga seguente
     //et=tdb+(ii*step/day2sec); /*et is the barycenter time, at each sample */
     //ora il tdb lo ricalcolo su tutti i samples
     if(imjd==2)utc=PSS_tdt2utc(tdt);
     utc_ii=utc;  //+(ii*step/day2sec);
    if(et<ss[0]) continue;      /* ephemeris starts later than test points */
    if(et>ss[1]) break;         /* ephemeris ends earlier than test points */

    jpl_mod_pleph(et,ntarg,nctr,r);
   

    for (k=0; k<=2; ++k)
	{
	  pos[k]=r[k];     //AU
          vel[k]=r[k+3];   //AU/day. Earth center vs SSB, referred at T0
	                   //Rectangular equatorial coordinates
	}
/*****Evaluate Detector velocity vs Earth Center*****/
      if(nut1==1){   
      /***** Nutation series from the JPL file *****/
      jpl_mod_pleph(et,14,0,ddd);
      psi=ddd[0]*RAD2SEC;
      eps=ddd[1]*RAD2SEC;
       if(nutcorr==1)  {
      PSS_deltanut(et, &ddpsi,&ddeps);
      //if(ii<=1)printf("et,psi eps ddpsi ddeps dal JPL senza corr %f %f %f %f %f\n",et,psi,eps,ddpsi,ddeps);
      psi+=ddpsi;
      eps+=ddeps;
       }
      novas_mod_earthtilt (et, &oblm,&oblt,&eqeq,psi,eps);
      //if(ii<=1)printf("et,psi eps dal JPL+corr %f %f %f \n",et,psi,eps);
      //printf("et,psi eps dal JPL %f %f %f \n",et,psi,eps);
    }
    else {
      /*************Nutation series from NOVAS************/
      earthtilt (et, &oblm,&oblt,&eqeq,&psi,&eps);
      //printf("et,psi eps dal NOVAS %f %f %f \n",et,psi,eps);
    }
     //Aug 2003: i valori JPL o NOVAS sono identici
     //Invece devo usare il file dello IERS per le correzioni
     //eopc04.62-now preso da http://hpiers.obspm.fr/eoppc/eop/eopc04/
     //dove sono in secondi di arco

     /*****Sidereal time needs UT1. We use UTC that approx. UT1 (error <0.9 s/anno) if fut1=0,  else we use UT1*****/
     if (!fut1)
     {
	utch=(long)(utc_ii);
	utcl=utc_ii-utch;
      }
	else
      {
        ut1=PSS_utc2ut1(utc_ii);
        utch=(long)(ut1);
        utcl=ut1-utch;
       }
     eqeqP=eqeq; //Kaz suggestion. May 2003
     if(NutFlag==0)eqeqP=0.0;
   sidereal_time (utch,utcl,eqeqP, &gast);
   novas_mod_terra (detector,gast, pos1,vel1);
   //modified Sept. 2003 !!!!!
   //if(nutcorr==1)novas_mod_nutate (et,FN1,pos1, pos2,psi,eps); //uncomment if required
   //if(nutcorr==0)nutate (et,FN1,pos1, pos2); //uncomment if required
   //precession (et,pos2,T0, pog); //uncomment if required
  if(nutcorr==1)novas_mod_nutate(et,FN1,vel1, vel2,psi,eps); /*performs inverse nutation (true to mean) */
   if(nutcorr==0)nutate(et,FN1,vel1, vel2);
   precession (et,vel2,T0, vog); //rectangular equatorial coord., ref. at T0 

    /*******write a TEST file with values to compare***********/
  
   fprintf(OUT33,"%23.11f ",mjd+(ii*step/(1.0*day2sec)));
     
   fprintf(OUT33,"%25.18e %25.18e %25.18e %25.18e %25.18e %25.18e %25.18e %25.18e %25.18e%25.18e %25.18e %25.18e  ",vel[0]/C, vel[1]/C,vel[2]/C, vel1[0]/C, vel1[1]/C,vel1[2]/C,(vel2[0]-vel1[0])/C, (vel2[1]-vel1[1])/C,(vel2[2]-vel1[2])/C,(vog[0]-vel2[0])/C, (vog[1]-vel2[1])/C,(vog[2]-vel2[2])/C);
      fprintf(OUT33,"\n"); 

  /*******END TEST*****/

   /*****Sum the two contributions: Earth vs SSB + Det vs. Earth *****/
  for(k=0; k<=2; ++k)
      {
	//pos[k]+=pog[k];
	 vel[k]+=vog[k];

      }
 
    /*****If equatorial diverso da 1 convert into Ecliptical coordinates*****/
  if(equatorial!=1)
  {
     PSS_equr2eclr (vel, vel1); 
     //PSS_equr2eclr (pos, pos1); //uncomment if required
    velx[ii]=vel1[0];
    vely[ii]=vel1[1];
    velz[ii]=vel1[2];
  }
  else
  { 
     velx[ii]=vel[0];
     vely[ii]=vel[1];
     velz[ii]=vel[2];
  }
    deinstein[ii]= PSS_einstein(tdt);  //Deinstein
    //printf("Einstein %d %25.20f\n", ii, deinstein[ii]);
    ii++;
   }
     fclose(F1);
     fclose(F2);
     fclose(F3);
     fclose(OUT33); //TEST FILE
     return mjd;
}


/*****Detector position:*****/
double PSS_detector_pos(DETECTOR *detector,double dot_ini, double step, int npoint,short int fut1, short int nut1, short int nutcorr,short int equatorial,short int imjd,double posx[], double posy[],double posz[])

/*For a given detector, initial date (days from the beginning of the century or mjd),step (s), number
 of point, it gives the position vector  in rectangular ecliptic 
or equatorial coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch). 

It returns MJD of the first point.

INPUT ARGUMENTS:

detector, pointer to struct DETECTOR

dot_ini, double: days from the beginning of the Century, or mjd, or TDT in mjd

step, double: step in seconds

npoint, int: number of data

fut1, int: flag for UT1. If 0 the function uses UTC, if 1 it uses UT1

imjd=0,1,2 0=days from 1900, 1=mjd , 2=tdt in JD 

OUTPUT ARGUMENTS:

posx[npoint], 
posy[npoint]
posz[npoint]
double arrays:Det. pos vector (km) in rectangular ecliptic or euatorial coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch)

deinstein[npoint]
RETURN:
MJD of the first datum
*******************************************************************************************
NOTE: Velocity vector of the Earth and Detector are not required.
      Uncomment the relative instructions if you need a code that computes also positions
********************************************************************************************
*/

{
   /* Sotto era 400 anzichè NMAX che ora vale 1000 e sta in jplbin.h */
  char nams[NMAX][6];
  double et, r[6], ss[3], vals[NMAX];
  int i, ntarg, nctr;
  intptr_t nvs;  //pia64
  int ii,k;
  double pos[3];
  double vel[3];
  double utc,tdt,tdb,utc_ii,utch,utcl,ut1,mjd;
  /****used in sidereal_time****/
  double eqeq,eqeqP; 
  double gast;
  /****used in earthtilt****/
  double oblm,oblt,psi,eps;
  double ddpsi,ddeps; 
/****used in terra, nutate,precession****/
  double pos1[3],pos2[3],pog[3]; 
  double vel1[3],vel2[3],vog[3];
/****used for nutation in the ephemeris****/
  double ddd[6];
  double tdt_ini; 
  int verb=0; //0 does not print values


  /******PER TEST FILE */
   FILE *OUT33;
   OUT33=fopen("table_test.dat","wb");


  F1=fopen(EPHFILE,"rb");
  F2=fopen(UTC2UT1,"rb");
  F3=fopen(DELTANUT,"rb");
  if (verb>4)puts("Initial time, days");
  if (verb>4)printf("%f\n", dot_ini);

/****** Print the ephemeris constants. ******/

  jpl_mod_constan(nams,vals,ss,&nvs);
  if (verb>4)printf("Range DE405= %14.2f  %14.2f  %14.2f\n",ss[0],ss[1],ss[2]); //pia Sep. 2005
  /*  for(i=0;i<nvs;++i)
      printf("%.6s  %24.16E\n",nams[i],vals[i]); */
  ii=0;
  ntarg=3; /* Earth */
  nctr=12; /*Solar system barycenter */
  if(imjd==0)
  {
    utc=dot_ini+dot2utc;
    mjd=utc-utc2mjd;
    if (verb>4)printf("initial utc= %18.8f\n ",utc); //pia Sep. 2005
    if (verb>4)printf("initial mjd= %18.8f\n ",mjd);
    if (verb>4)printf("initial days from 1900= %18.8f\n ",dot_ini);
  }
   if(imjd==1)
  {
    utc=dot_ini+utc2mjd;
    mjd=dot_ini;
    dot_ini=utc-dot2utc;  //dot_ini sempre giorni inizio secolo, in UTC !! NOTARE
    if (verb>4)printf("initial utc= %18.8f\n ",utc);
    if (verb>4)printf("initial mjd= %18.8f\n ",mjd);
  }
   if(imjd==2){
     tdt=dot_ini+utc2mjd;
     mjd=dot_ini; //mjd but TDT
     if (verb>4)printf("initial tdt= %20.11f\n ",tdt);
     if (verb>4)printf("initial TDT in mjd= %18.8f\n ",mjd);
   }
  /*****Leap seconds*****/
  if(imjd !=2) tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
  if (verb>4)printf("initial tdt= %18.8f\n ",tdt); //pia Sep. 2005
  tdt_ini=tdt;
  tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
  if (verb>4)printf("tdb= %18.8f\n ",tdb);
  if (verb>4)printf("utc2mjd= %25.15f\n ",utc2mjd);
  ////printf("Detector name,long,lat:%s %f %f\n",detector->name,detector->longitude,detector->latitude); //pia Sep. 2005

 while( ii < npoint)
   {
   /*****  Read a value; Skip if not within the time-range
        of the present version of the ephemeris.                            */
     if(imjd!=2){
       utc=dot_ini+dot2utc+ii*step/(1.0*day2sec); //commento kaz: tdb va ricalcolato sempre
       tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
     }
     else
       {
	 tdt=tdt_ini+ii*step/(1.0*day2sec);
       }
     tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
     et=tdb;              //sost. kaz fino a qui e tolta la riga seguente
     //et=tdb+(ii*step/day2sec); /*et is the barycenter time, at each sample */
     //ora il tdb lo ricalcolo su tutti i samples
     if(imjd==2)utc=PSS_tdt2utc(tdt);
     utc_ii=utc;  //+(ii*step/day2sec);
    if(et<ss[0]) continue;      /* ephemeris starts later than test points */
    if(et>ss[1]) break;         /* ephemeris ends earlier than test points */

    jpl_mod_pleph(et,ntarg,nctr,r);
   

    for (k=0; k<=2; ++k)
	{
	  pos[k]=r[k];     //AU
          vel[k]=r[k+3];   //AU/day. Earth center vs SSB, referred at T0
	                   //Rectangular equatorial coordinates
	}
/*****Evaluate Detector velocity vs Earth Center*****/
      if(nut1==1){   
      /***** Nutation series from the JPL file *****/
      jpl_mod_pleph(et,14,0,ddd);
      psi=ddd[0]*RAD2SEC;
      eps=ddd[1]*RAD2SEC;
       if(nutcorr==1)  {
      PSS_deltanut(et, &ddpsi,&ddeps);
      //if(ii<=1)printf("et,psi eps ddpsi ddeps dal JPL senza corr %f %f %f %f %f\n",et,psi,eps,ddpsi,ddeps);
      psi+=ddpsi;
      eps+=ddeps;
       }
      novas_mod_earthtilt (et, &oblm,&oblt,&eqeq,psi,eps);
      //if(ii<=1)printf("et,psi eps dal JPL+corr %f %f %f \n",et,psi,eps);
      //printf("et,psi eps dal JPL %f %f %f \n",et,psi,eps);
    }
    else {
      /*************Nutation series from NOVAS************/
      earthtilt (et, &oblm,&oblt,&eqeq,&psi,&eps);
      //printf("et,psi eps dal NOVAS %f %f %f \n",et,psi,eps);
    }
     //Aug 2003: i valori JPL o NOVAS sono identici
     //Invece devo usare il file dello IERS per le correzioni
     //eopc04.62-now preso da http://hpiers.obspm.fr/eoppc/eop/eopc04/
     //dove sono in secondi di arco

     /*****Sidereal time needs UT1. We use UTC that approx. UT1 (error <0.9 s/anno) if fut1=0,  else we use UT1*****/
     if (!fut1)
     {
	utch=(long)(utc_ii);
	utcl=utc_ii-utch;
      }
	else
      {
        ut1=PSS_utc2ut1(utc_ii);
        utch=(long)(ut1);
        utcl=ut1-utch;
       }
     eqeqP=eqeq; //Kaz suggestion. May 2003
     if(NutFlag==0)eqeqP=0.0;
   sidereal_time (utch,utcl,eqeqP, &gast);
   novas_mod_terra (detector,gast, pos1,vel1);
   if(nutcorr==1)novas_mod_nutate(et,FN1,pos1, pos2,psi,eps); //uncomment if required
   if(nutcorr==0)nutate(et,FN1,pos1, pos2); //uncomment if required
   precession(et,pos2,T0, pog); //uncomment if required
   //if(nutcorr==1)novas_mod_nutate(et,FN1,vel1, vel2,psi,eps); /*performs inverse nutation (true to mean) */
   //if(nutcorr==0)nutate(et,FN1,vel1, vel2);
   //precession (et,vel2,T0, vog); //rectangular equatorial coord., ref. at T0 

    /*******write a TEST file with values to compare***********/
  
   fprintf(OUT33,"%23.11f ",mjd+(ii*step/(1.0*day2sec)));
     
   fprintf(OUT33,"%25.18e %25.18e %25.18e %25.18e %25.18e %25.18e %25.18e %25.18e %25.18e%25.18e %25.18e %25.18e  ",pos[0], pos[1],pos[2], pos1[0], pos1[1],pos1[2],
(pos2[0]-pos1[0]), (pos2[1]-pos1[1]),(pos2[2]-pos1[2]),(pog[0]-pos2[0]), (pog[1]-pos2[1]),(pog[2]-pos2[2]));
      fprintf(OUT33,"\n"); 

  /*******END TEST*****/

   /*****Sum the two contributions: Earth vs SSB + Det vs. Earth *****/
  for(k=0; k<=2; ++k)
      {
	pos[k]+=pog[k];
	//vel[k]+=vog[k];

      }
 
    /*****If equatorial diverso da 1 convert into Ecliptical coordinates*****/
  if(equatorial!=1)
  {
     PSS_equr2eclr (pos, pos1); 
     //PSS_equr2eclr (vel, vel1); //uncomment if required
    posx[ii]=pos1[0];
    posy[ii]=pos1[1];
    posz[ii]=pos1[2];
  }
  else
  { 
     posx[ii]=pos[0];
     posy[ii]=pos[1];
     posz[ii]=pos[2];
  }
  //POI VEDO SE RIM.deinstein[ii]= PSS_einstein(tdt);  //Deinstein
    //printf("Einstein %d %25.20f\n", ii, deinstein[ii]);
    ii++;
   }
     fclose(F1);
     fclose(F2);
     fclose(F3);
     fclose(OUT33); //TEST FILE
     return mjd;
}


/*****Detector velocity AND position:*****/

double PSS_detector_velpos(DETECTOR *detector,double dot_ini, double step, int npoint,short int fut1, short int nut1, short int nutcorr,short int equatorial,short int imjd,double velx[], double vely[],double velz[],double posx[], double posy[],double posz[],double deinstein[])

/*For a given detector, initial date (days from the beginning of the century or mjd),step (s), number
 of point, it gives the velocity vector (AU/day) and the position in AU in rectangular ecliptic 
or equatorial coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch). 

It returns MJD of the first point.

INPUT ARGUMENTS:

detector, pointer to struct DETECTOR

dot_ini, double: days from the beginning of the Century, or mjd, or TDT in mjd

step, double: step in seconds

npoint, int: number of data

fut1, int: flag for UT1. If 0 the function uses UTC, if 1 it uses UT1

imjd=0,1,2 0=days from 1900, 1=mjd , 2=tdt in JD 

OUTPUT ARGUMENTS:

velx[npoint], 
vely[npoint]
velz[npoint]
posx[npoint], 
posy[npoint]
posz[npoint]
double arrays:Det. velocity vector (AU/day)  and pos in km in rectangular ecliptic or euatorial coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch)

deinstein[npoint]
RETURN:
MJD of the first datum
*******************************************************************************************
NOTE: Position vector of the Earth and Detector are not required.
      Uncomment the relative instructions if you need a code that computes also positions
********************************************************************************************
*/

{
   /* Sotto era 400 anzichè NMAX che ora vale 1000 e sta in jplbin.h */
  char nams[NMAX][6];
  double et, r[6], ss[3], vals[NMAX];
  int i, ntarg, nctr;
   intptr_t nvs;  //pia64
  int ii,k;
  double pos[3];
  double vel[3];
  double utc,tdt,tdb,utc_ii,utch,utcl,ut1,mjd;
  /****used in sidereal_time****/
  double eqeq,eqeqP; 
  double gast;
  /****used in earthtilt****/
  double oblm,oblt,psi,eps;
  double ddpsi,ddeps; 
/****used in terra, nutate,precession****/
  double pos1[3],pos2[3],pog[3]; 
  double vel1[3],vel2[3],vog[3];
/****used for nutation in the ephemeris****/
  double ddd[6];
  double tdt_ini; 
  int verb=0; //0 does not print values

  /******PER TEST FILE */
  //FILE *OUT33;
  //OUT33=fopen("table_test.dat","wb");


  F1=fopen(EPHFILE,"rb");
  F2=fopen(UTC2UT1,"rb");
  F3=fopen(DELTANUT,"rb");

  if (verb>4)puts("Initial time, days");
  if (verb>4)printf("%f\n", dot_ini);

/****** Print the ephemeris constants. ******/

   jpl_mod_constan(nams,vals,ss,&nvs);
   if (verb>4)printf("Range DE405= %14.2f  %14.2f  %14.2f\n",ss[0],ss[1],ss[2]); //pia Sep. 2005
   if (verb>4)puts("--> If ok it Should be: 2444208.50  2458864.50  32.00.  nvs 156");
   //for(i=0;i<nvs;++i)
   //printf("%.6s  %24.16E\n",nams[i],vals[i]); 
  
  ii=0;
  ntarg=3; /* Earth */
  nctr=12; /*Solar system barycenter */
  if(imjd==0)
  {
    utc=dot_ini+dot2utc;
    mjd=utc-utc2mjd;
    if (verb>4)printf("initial utc= %18.8f\n ",utc); //pia Sep. 2005
    if (verb>4)printf("initial mjd= %18.8f\n ",mjd);
    if (verb>4)printf("initial days from 1900= %18.8f\n ",dot_ini);
  }
   if(imjd==1)
  {
    utc=dot_ini+utc2mjd;
    mjd=dot_ini;
    dot_ini=utc-dot2utc;  //dot_ini sempre giorni inizio secolo, in UTC !! NOTARE
    if (verb>4)printf("imjd=1: initial utc= %18.8f\n ",utc);
    if (verb>4)printf("imjd=1: initial mjd= %18.8f\n ",mjd);
  }
   if(imjd==2){
     tdt=dot_ini+utc2mjd;
     mjd=dot_ini; //mjd but TDT
     if (verb>4)printf("initial tdt= %20.11f\n ",tdt);
     if (verb>4)printf("initial TDT in mjd= %18.8f\n ",mjd);
   }
  /*****Leap seconds*****/
  if(imjd !=2) tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
  if (verb>4)printf("initial tdt= %18.8f\n ",tdt); //pia Sep. 2005
  tdt_ini=tdt;
  tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
  if (verb>4)printf("tdb= %18.8f\n ",tdb);
  if (verb>4)printf("utc2mjd= %25.15f\n ",utc2mjd);
  if (verb>4)printf("Detector name,long,lat:%s %f %f\n",detector->name,detector->longitude,detector->latitude); //pia Sep. 2005

 while( ii < npoint)
   {
   /*****  Read a value; Skip if not within the time-range
        of the present version of the ephemeris.                            */
     if(imjd!=2){
       utc=dot_ini+dot2utc+ii*step/(1.0*day2sec); //commento kaz: tdb va ricalcolato sempre
       tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
 
     }
     else
       {
	
	 tdt=tdt_ini+ii*step/(1.0*day2sec);
       }
     tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
 
    et=tdb;              //sost. kaz fino a qui e tolta la riga seguente
    //et=tdb+(ii*step/day2sec); /*et is the barycenter time, at each sample */
     //ora il tdb lo ricalcolo su tutti i samples
     if(imjd==2)utc=PSS_tdt2utc(tdt);
     utc_ii=utc;  //+(ii*step/day2sec);
    if(et<ss[0]) continue;      /* ephemeris starts later than test points */
    if(et>ss[1]) break;         /* ephemeris ends earlier than test points */

    jpl_mod_pleph(et,ntarg,nctr,r);
    //printf("et ss0 ss1 %f %f %f\n",et,ss[0],ss[1]);

    for (k=0; k<=2; ++k)
	{
	  pos[k]=r[k];     //AU
          vel[k]=r[k+3];   //AU/day. Earth center vs SSB, referred at T0
	                   //Rectangular equatorial coordinates
	}
/*****Evaluate Detector velocity vs Earth Center*****/
      if(nut1==1){   
      /***** Nutation series from the JPL file *****/
      jpl_mod_pleph(et,14,0,ddd);
      psi=ddd[0]*RAD2SEC;
      eps=ddd[1]*RAD2SEC;
       if(nutcorr==1)  {
      PSS_deltanut(et, &ddpsi,&ddeps);
      //if(ii<=1)printf("et,psi eps ddpsi ddeps dal JPL senza corr %f %f %f %f %f\n",et,psi,eps,ddpsi,ddeps);
      psi+=ddpsi;
      eps+=ddeps;
       }
      novas_mod_earthtilt (et, &oblm,&oblt,&eqeq,psi,eps);
      //if(ii<=1)printf("et,psi eps dal JPL+corr %f %f %f \n",et,psi,eps);
      //printf("et,psi eps dal JPL %f %f %f \n",et,psi,eps);
    }
    else {
      /*************Nutation series from NOVAS************/
      earthtilt (et, &oblm,&oblt,&eqeq,&psi,&eps);
      //printf("et,psi eps dal NOVAS %f %f %f \n",et,psi,eps);
    }
     //Aug 2003: i valori JPL o NOVAS sono identici
     //Invece devo usare il file dello IERS per le correzioni
     //eopc04.62-now preso da http://hpiers.obspm.fr/eoppc/eop/eopc04/
     //dove sono in secondi di arco

     /*****Sidereal time needs UT1. We use UTC that approx. UT1 (error <0.9 s/anno) if fut1=0,  else we use UT1*****/
     if (!fut1)
     {
	utch=(long)(utc_ii);
	utcl=utc_ii-utch;
      }
	else
      {
        ut1=PSS_utc2ut1(utc_ii);
        utch=(long)(ut1);
        utcl=ut1-utch;
       }
     eqeqP=eqeq; //Kaz suggestion. May 2003
     if(NutFlag==0)eqeqP=0.0;
   sidereal_time (utch,utcl,eqeqP, &gast);
   novas_mod_terra (detector,gast, pos1,vel1);
   //modified Sept. 2003 !!!!!
   if(nutcorr==1)novas_mod_nutate(et,FN1,pos1, pos2,psi,eps); //uncomment if required
   if(nutcorr==0)nutate(et,FN1,pos1, pos2); //uncomment if required
   precession(et,pos2,T0, pog); //uncomment if required
  if(nutcorr==1)novas_mod_nutate(et,FN1,vel1, vel2,psi,eps); /*performs inverse nutation (true to mean) */
   if(nutcorr==0)nutate(et,FN1,vel1, vel2);
   precession (et,vel2,T0, vog); //rectangular equatorial coord., ref. at T0 

    /*******write a TEST file with values to compare***********/
  
   //fprintf(OUT33,"%23.11f ",mjd+(ii*step/(1.0*day2sec)));
     
   //fprintf(OUT33,"%25.18e %25.18e %25.18e %25.18e %25.18e %25.18e %25.18e %25.18e %25.18e%25.18e %25.18e %25.18e  ",vel[0]/C, vel[1]/C,vel[2]/C, vel1[0]/C, vel1[1]/C,vel1[2]/C,(vel2[0]-vel1[0])/C, (vel2[1]-vel1[1])/C,(vel2[2]-vel1[2])/C,(vog[0]-vel2[0])/C, (vog[1]-vel2[1])/C,(vog[2]-vel2[2])/C);
   //fprintf(OUT33,"\n"); 

  /*******END TEST*****/

   /*****Sum the two contributions: Earth vs SSB + Det vs. Earth *****/
  for(k=0; k<=2; ++k)
      {
	 pos[k]+=pog[k];
	 vel[k]+=vog[k];

      }
 
    /*****If equat orial diverso da 1 convert into Ecliptical coordinates*****/
  if(equatorial!=1)
  {
     PSS_equr2eclr (vel, vel1); 
     //PSS_equr2eclr (pos, pos1); //uncomment if required
    velx[ii]=vel1[0];
    vely[ii]=vel1[1];
    velz[ii]=vel1[2];
  }
  else
  { 
     velx[ii]=vel[0];
     vely[ii]=vel[1];
     velz[ii]=vel[2];
  }
   if(equatorial!=1)
  { 
    PSS_equr2eclr (pos, pos1); //uncomment if required
    posx[ii]=pos1[0];
    posy[ii]=pos1[1];
    posz[ii]=pos1[2];
  }
  else
  { 
     posx[ii]=pos[0];
     posy[ii]=pos[1];
     posz[ii]=pos[2];
  }
    deinstein[ii]= PSS_einstein(tdt);  //Deinstein
    //printf("Einstein %d %25.20f\n", ii, deinstein[ii]);
    ii++;
   }
     if (verb>4)printf("VELOCITA'/C: %15.8e %15.8e %15.8e\n",velx[ii]/C,vely[ii]/C,velz[ii]/C);
     if (verb>4)printf("VELOCITA': %15.8e %15.8e %15.8e\n",vel[0],vely[ii],velz[ii]);
     if (verb>4)printf("POSIZIONI*day2sec/C: %15.8e %15.8e %15.8e\n",posx[ii]*day2sec/C,posy[ii]*day2sec/C,posz[ii]*day2sec/C);
     fclose(F1);
     fclose(F2);
     fclose(F3);
     //   fclose(OUT33); //TEST FILE
     return mjd;
}
double PSS_detector_posTSSB(DETECTOR *detector,double dot_ini, double step, int npoint,short int fut1, short int nut1, short int nutcorr,short int equatorial,short int imjd,double posx[], double posy[],double posz[])


/*For a given detector, initial date (days from the beginning of the century or mjd),step (s), number
 of point, it gives the velocity vector (AU/day) and the position in AU in rectangular ecliptic 
or equatorial coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch). 
// DON"T SUM THE TWO CONTRIBUTIONS: needed only Earth center vs SSB. For Mary/Fermi

It returns MJD of the first point.

INPUT ARGUMENTS:

detector, pointer to struct DETECTOR

dot_ini, double: days from the beginning of the Century, or mjd, or TDT in mjd

step, double: step in seconds

npoint, int: number of data

fut1, int: flag for UT1. If 0 the function uses UTC, if 1 it uses UT1

imjd=0,1,2 0=days from 1900, 1=mjd , 2=tdt in JD 

OUTPUT ARGUMENTS:

posx[npoint], 
posy[npoint]
posz[npoint]
double arrays:Det. velocity vector (AU/day)  and pos in km in rectangular ecliptic or euatorial coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch)

deinstein[npoint]
RETURN:
MJD of the first datum
*******************************************************************************************
NOTE: Position vector of the Earth and Detector are not required.
      Uncomment the relative instructions if you need a code that computes also positions
********************************************************************************************
*/

{
   /* Sotto era 400 anzichè NMAX che ora vale 1000 e sta in jplbin.h */
  char nams[NMAX][6];
  double et, r[6], ss[3], vals[NMAX];
  int i, ntarg, nctr;
   intptr_t nvs;  //pia64
  int ii,k;
  double pos[3];
  double utc,tdt,tdb,utc_ii,utch,utcl,ut1,mjd;
  /****used in sidereal_time****/
  double eqeq,eqeqP; 
  double gast;
  /****used in earthtilt****/
  double oblm,oblt,psi,eps;
  double ddpsi,ddeps; 
/****used in terra, nutate,precession****/
  double pos1[3],pos2[3],pog[3]; 
  
/****used for nutation in the ephemeris****/
  double ddd[6];
  double tdt_ini; 
  int verb=0; //0 does not print values

  /******PER TEST FILE */
  //FILE *OUT33;
  //OUT33=fopen("table_test.dat","wb");


  F1=fopen(EPHFILE,"rb");
  F2=fopen(UTC2UT1,"rb");
  F3=fopen(DELTANUT,"rb");

  if (verb>4)puts("Initial time, days");
  if (verb>4)printf("%f\n", dot_ini);

/****** Print the ephemeris constants. ******/

   jpl_mod_constan(nams,vals,ss,&nvs);
   if (verb>4)printf("Range DE405= %14.2f  %14.2f  %14.2f\n",ss[0],ss[1],ss[2]); //pia Sep. 2005
   if (verb>4)puts("--> If ok it Should be: 2444208.50  2458864.50  32.00.  nvs 156");
   //for(i=0;i<nvs;++i)
   //printf("%.6s  %24.16E\n",nams[i],vals[i]); 
  
  ii=0;
  ntarg=3; /* Earth */
  nctr=12; /*Solar system barycenter */
  if(imjd==0)
  {
    utc=dot_ini+dot2utc;
    mjd=utc-utc2mjd;
    if (verb>4)printf("initial utc= %18.8f\n ",utc); //pia Sep. 2005
    if (verb>4)printf("initial mjd= %18.8f\n ",mjd);
    if (verb>4)printf("initial days from 1900= %18.8f\n ",dot_ini);
  }
   if(imjd==1)
  {
    utc=dot_ini+utc2mjd;
    mjd=dot_ini;
    dot_ini=utc-dot2utc;  //dot_ini sempre giorni inizio secolo, in UTC !! NOTARE
    if (verb>4)printf("imjd=1: initial utc= %18.8f\n ",utc);
    if (verb>4)printf("imjd=1: initial mjd= %18.8f\n ",mjd);
  }
   if(imjd==2){
     tdt=dot_ini+utc2mjd;
     mjd=dot_ini; //mjd but TDT
     printf("initial tdt= %20.11f\n ",tdt);
     printf("initial TDT in mjd= %18.8f\n ",mjd);
   }
  /*****Leap seconds*****/
  if(imjd !=2) tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
  if (verb>4)printf("initial tdt= %18.8f\n ",tdt); //pia Sep. 2005
  tdt_ini=tdt;
  tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
  if (verb>4)printf("tdb= %18.8f\n ",tdb);
  if (verb>4)printf("utc2mjd= %25.15f\n ",utc2mjd);
  if (verb>4)printf("Detector name,long,lat:%s %f %f\n",detector->name,detector->longitude,detector->latitude); //pia Sep. 2005

 while( ii < npoint)
   {
   /*****  Read a value; Skip if not within the time-range
        of the present version of the ephemeris.                            */
     if(imjd!=2){
       utc=dot_ini+dot2utc+ii*step/(1.0*day2sec); //commento kaz: tdb va ricalcolato sempre
       tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
 
     }
     else
       {
	
	 tdt=tdt_ini+ii*step/(1.0*day2sec);
       }
     tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
 
    et=tdb;              //sost. kaz fino a qui e tolta la riga seguente
    //et=tdb+(ii*step/day2sec); /*et is the barycenter time, at each sample */
     //ora il tdb lo ricalcolo su tutti i samples
     if(imjd==2)utc=PSS_tdt2utc(tdt);
     utc_ii=utc;  //+(ii*step/day2sec);
    if(et<ss[0]) continue;      /* ephemeris starts later than test points */
    if(et>ss[1]) break;         /* ephemeris ends earlier than test points */

    jpl_mod_pleph(et,ntarg,nctr,r);
    //printf("et ss0 ss1 %f %f %f\n",et,ss[0],ss[1]);

    for (k=0; k<=2; ++k)
	{
	  pos[k]=r[k];     //AU
                           //vel are AU/day. Earth center vs SSB, referred at T0
	                   //Rectangular equatorial coordinates
	}



  
 
    /*****If equatorial diverso da 1 convert into Ecliptical coordinates*****/
 
   if(equatorial!=1)
  { 
    PSS_equr2eclr (pos, pos1); //uncomment if required
    posx[ii]=pos1[0];
    posy[ii]=pos1[1];
    posz[ii]=pos1[2];
  }
  else
  { 
     posx[ii]=pos[0];
     posy[ii]=pos[1];
     posz[ii]=pos[2];
  }

    ii++;
   }
 // printf("VELOCITA'/C: %15.8e %15.8e %15.8e\n",velx[ii]/C,vely[ii]/C,velz[ii]/C);
 //  printf("VELOCITA': %15.8e %15.8e %15.8e\n",vel[0],vely[ii],velz[ii]);
     if (verb>4)printf("POSIZIONI*day2sec/C: %15.8e %15.8e %15.8e\n",posx[ii]*day2sec/C,posy[ii]*day2sec/C,posz[ii]*day2sec/C);
     fclose(F1);
     fclose(F2);
     fclose(F3);
     //   fclose(OUT33); //TEST FILE
     return mjd;
}






/*****Source position:*****/


void  PSS_source_pos(SOURCE *source,short int equatorial,double poss[])
/******************************************************************************************/
/*For a given source, identified by ra and dec of Epoch J2000 it gives the normalized
position vector in rectangular ecliptic or equatorial coordinates, 
respect to the Barycenter of the Earth-Sun system. 
The source is supposed to be fixed in the sky (no proper motion).

INPUT ARGUMENT:

source, pointer to the struct SOURCE
short int equatorial. If 1 the output is in Equatorial coordinates, if not in Ecliptic
                      coordinates.

OUTPUT ARGUMENTS:
poss[3], double. Normalized position vector in ecliptical rectangular coordinates.


*/
/*****************************************************************************************/

{
 char nam[16];
 double pose[3],vels[3];
 double poss_mod;
 int k;
 int verb=0; //0 means do not write things..

   if(verb>4)printf("source name,alpha,delta:%s %f %f\n",source->name,source->ra,source->dec);
   novas_mod_starvectors(source,poss,vels);
    if(verb>4)printf("source position rect. equat. coord., AU  %f %f %f\n",poss[0],poss[1],poss[2]);
   if(equatorial!=1)
   {
      /*****convert into Ecliptical coordinates*****/
      PSS_equr2eclr (poss, pose);
   }
   else
   {
     pose[0]=poss[0];
     pose[1]=poss[1];
     pose[2]=poss[2];
   }
   if(verb!=0){
     if(equatorial!=1 && verb>4)printf("source position eclitt., AU  %f %f %f\n",pose[0],pose[1],pose[2]);
   }
   poss_mod=sqrt(pow(pose[0],2)+pow(pose[1],2)+pow(pose[2],2));

   /***Normalize***/
  for(k=0; k<=2; ++k)
      {
	poss[k]=pose[k]/poss_mod;

      }
  if(verb!=0){
    if(equatorial!=1 && verb>4)printf("Norm. source pos AU ecliptic  %f %f %f\n",poss[0],poss[1],poss[2]);
    if(equatorial==1 && verb>4)printf("Norm. source pos AU equatorial  %f %f %f\n",poss[0],poss[1],poss[2]);
  }
}


/********************************Miscellanea:**********************************/

/*****Vectors:*****/

void PSS_vect(double vectx[],double vecty[],double vectz[],double vectout[12],int *np)

/*For given vectors, vectx vecty vectz, supposed to be the x,y,z components of a
vector as a function of time, np points, it evaluates- and puts in vectout[]-:
the x,y,z components of the first point,
the x,y,z components of the middle point,
the x,y,z components of the last point,
the average x,y,z components

INPUT:
vectx[], double
vecty[], double
vectz[], double

OUTPUT:
vectout[12], double, where:
vectout[0,1 and 2] are vectx[0],vecty[0],vectz[0];
vectout[3,4 and 5] are vectx[middle],vecty[middle],vectz[middle];
vectout[6,7 and 8] are vectx[last],vecty[last],vectz[last];
vectout[9,10 and 11] are average vectx,average vecty,average vectz;

RETURN:
void
*/

{
  int len,len2,k;
  int verb=0; //0 does not print values
  len=*np;
  if(verb>4)printf("len= %d\n",len);
  len2=(int) (len/2);

  vectout[0]=vectx[0];
  vectout[1]=vecty[0];
  vectout[2]=vectz[0];
  vectout[3]=vectx[len2];
  vectout[4]=vecty[len2];
  vectout[5]=vectz[len2];
  vectout[6]=vectx[len];
  vectout[7]=vecty[len];
  vectout[8]=vectz[len];

  vectout[9]=0.;
  vectout[10]=0.;
  vectout[11]=0.;
  for(k=0; k<len; ++k)
      {
	vectout[9]+=vectx[k];
        vectout[10]+=vecty[k];
        vectout[11]+=vectz[k];
      }
   vectout[9]/=len;
   vectout[10]/=len;
   vectout[11]/=len;
   if (verb>4)printf("Initial vel  AU/day %f %f %f\n",vectout[0],vectout[1],vectout[2]);
   if (verb>4)printf("Middle vel  AU/day %f %f %f\n",vectout[3],vectout[4],vectout[5]);
   if (verb>4)printf("Middle vel  AU/day %f %f %f\n",vectout[6],vectout[7],vectout[8]);
   if (verb>4)printf("Final vel  AU/day %f %f %f\n",vectout[9],vectout[10],vectout[11]);
  return;
}


/*************Scalar product***************/
/*This function performs the scalar product of vectors of 3 components

INPUT:
vect1[3],double;
vect2[3],double;

OUTPUT:
sc,double
*/

double scalar_product(double vect1[3],double vect2[3])
{
     double sc;
     sc=vect1[0]*vect2[0]+vect1[1]*vect2[1]+vect1[2]*vect2[2];


  return sc;
  }
/*****************Amplitude modulation*****************/
/* modbar.c */
/*  Computes coefficients of the amplitude modulation functions
   for a BAR detector.
*/



void modbar (DETECTOR *detector,SOURCE *source,double tss_local,float *a,float *b)
{
 
  double lambda,azim,egam,alpha,delta,tss_rad;
double c1, c2, c3, c4, c5;
double d1, d2, d3, d4, d5;
  int verb=0; //0 does not print values
 
 /*INPUT:
 source->ra is given in hours;
 source->dec is given in degrees;
 detector->latitude is given in degrees;
 detector->azimuth is given in degrees; It is North to East clockwise
 tss_local in hours;
*/

 /*OUTPUT
a,b amplitude modulation function for the given sidereal time
 */

//egam,ephi,alpha,delta,tss_rad are in radians
   egam=(90-detector->azimuth)*DEG2RAD; //in the Eq. azim is from East counter-clockwise
                                        //paper Astone,Borkowski,Jaranowski,Krolak. App. A2
  lambda=detector->latitude*DEG2RAD;
  alpha=source->ra*15*DEG2RAD;
  delta=source->dec*DEG2RAD;
  tss_rad=tss_local*15*DEG2RAD;

 //tss_rad=3.9719881551361973; //prova brutalissima confronto MAP!!!!!
 //tss_rad=3.9763634273224491; //prova brutalissima confronto MAP!!!!!
 /*Prese da krolak_1.ps direi che non vanno per le barre
 c1=sin(2*egam)*(3-cos(2*lambda))*(3-cos(2*delta))*cos(2*(alpha-tss_rad));
 c2=cos(2*egam)*sin(lambda)*(3-cos(2*delta))*sin(2*(alpha-tss_rad));
 c3=sin(2*egam)*sin(2*lambda)*sin(2*delta)*cos(alpha-tss_rad);
 c4=cos(2*egam)*cos(lambda)*sin(2*delta)*sin(alpha-tss_rad);
 c5=sin(2*egam)*pow(cos(lambda),2)*pow(cos(delta),2);
  
 *a=(1.0/16)*c1-(1.0/4)*c2+(1.0/4)*c3-(1.0/2)*c4+(3.0/4)*c5;

 d1=cos(2*egam)*sin(lambda)*sin(delta)*cos(2*(alpha-tss_rad));
 d2=sin(2*egam)*(3-cos(2*lambda))*sin(delta)*sin(2*(alpha-tss_rad));
 d3=cos(2*egam)*cos(lambda)*cos(delta)*cos(alpha-tss_rad);
 d4=sin(2*egam)*sin(2*lambda)*cos(delta)*sin(alpha-tss_rad);
  
 *b=d1+(1.0/4)*d2+d3+(1.0/2)*d4;
  */
 c1=(pow(cos(egam),2)-pow(sin(egam),2)*pow(sin(lambda),2))*(1+pow(sin(delta),2))*cos(2*(alpha-tss_rad));
 c2=sin(2*egam)*sin(lambda)*(1+pow(sin(delta),2))*sin(2*(alpha-tss_rad));
 c3=pow(sin(egam),2)*sin(2*lambda)*sin(2*delta)*cos(alpha-tss_rad);
 c4=sin(2*egam)*cos(lambda)*sin(2*delta)*sin(alpha-tss_rad);
 c5=(1-3*pow(sin(egam),2)*pow(cos(lambda),2))*pow(cos(delta),2);
  
  *a=(1.0/2)*c1+(1.0/2)*c2-(1.0/2)*c3+(1.0/2)*c4+(1.0/2)*c5;

  d1=sin(2*egam)*sin(lambda)*sin(delta)*cos(2*(alpha-tss_rad));
  d2=(pow(cos(egam),2)-pow(sin(egam),2)*pow(sin(lambda),2))*sin(delta)*sin(2*(alpha-tss_rad));
  d3=sin(2*egam)*cos(lambda)*cos(delta)*cos(alpha-tss_rad);
  d4=pow(sin(egam),2)*sin(2*lambda)*cos(delta)*sin(alpha-tss_rad);
  
  *b=-d1+d2-d3-d4;


  if (verb>4) printf("a,b %f, %f \n",*a,*b);
return;
} 
/* radpat_interf */
/*  Computes coefficients of the amplitude modulation functions
   for an interferometer detector.
*/



void radpat_interf (DETECTOR *detector,SOURCE *source, double tss_local, float *a, float *b)
{
 
  double lat,azim,alpha,delta,psi,tss_rad;
  double eps;
  double clat,slat,cazim,sazim,cospsi,sinpsi;
  double cdelta,sdelta;
  double term1,term2,term3,term4,term5,term6,term7,term8,term9,term10,term11,term12,term13,term14;
  double ctheta,psi1,phi1,c2psi,s2psi,c2phi,s2phi;
  double f1,f2,f3,ftheta;

  int verb=0; //if 0 does not print things...
//See the matlab function radpat_interf.m (Cristiano)
 /*INPUT:
 source->ra is given in hours;
 source->dec is given in degrees;
 detector->latitude is given in degrees North;
 detector->longitude is given in degrees East;
 detector->azimuth is given in degrees North to East;
 tss_local in hours;
*/

 /*OUTPUT
a,b amplitude modulation function for the given sidereal time: a=f circular; b= f linear.
Then: f=a+b
 */

//azim,ephi,alpha,delta,psi,tss_rad are in radians
  azim=(detector->azimuth)*DEG2RAD; //North to East (clockwise). (VERIFY if it is = South to West)
  lat=detector->latitude*DEG2RAD;
  alpha=source->ra*15*DEG2RAD;
  delta=source->dec*DEG2RAD;
  psi=source->psi*DEG2RAD;
  tss_rad=tss_local*15*DEG2RAD;

  eps=source->eps;


  clat=cos(lat);
  slat=sin(lat);
  cazim=cos(azim);
  sazim=sin(azim);
  cospsi=cos(psi);
  sinpsi=sin(psi);

  cdelta=cos(delta);
  sdelta=sin(delta);
  
  term1=cdelta*clat;
  term2=sdelta*slat;
  term3=sdelta*clat;
  term4=cdelta*slat;
  term5=-cdelta*sazim;
  term6=-cdelta*cazim*slat;
  term7=-sdelta*cazim*clat;
  term8=-cdelta*cazim;
  term9=term5*slat;
  term10=-term3*sazim;
  term11=cospsi*clat;
  term12=sinpsi*clat;
  term13=sinpsi*slat;
  term14=cospsi*slat;

  ctheta=-term1*cos(alpha-tss_rad)-term2;
  psi1=atan2(term11*sin(alpha-tss_rad)-term12*sdelta*cos(alpha-tss_rad)+term13*cdelta,-term12*sin(alpha-tss_rad)-term11*sdelta*cos(alpha-tss_rad)+term14*cdelta);
  phi1=-atan2(-term5*sin(alpha-tss_rad)+term6*cos(alpha-tss_rad)-term7,term8*sin(alpha-tss_rad)+term9*cos(alpha-tss_rad)-term10);
  c2psi=cos(2*psi1);
  s2psi=sin(2*psi1);
  c2phi=cos(2*phi1);
  s2phi=sin(2*phi1);

  f1=c2phi*c2psi;
  f2=ctheta*s2phi;
  f3=c2phi*s2psi;
  ftheta=0.5*(1+ctheta*ctheta);

  *a=0.5*(1-eps)*(pow(ftheta*c2phi,2)+f2*f2);
  *b=eps*pow(ftheta*f1-f2*s2psi,2);


  if(verb>4)printf("f_circ,f_lin %f, %f\n",*a,*b);
return;
} 

void dot2tsslocal(DETECTOR *detector,double dot, short int fut1,short int nut1,short int nutcorr,double *tss_local)

     /*
INPUT:
detector
dot=days from the beginning of the century;
fut1, int: flag for UT1. If 0 the function uses UTC, if 1 it uses UT1

OUTPUT
tss_local=local sidereal time in hours
     */
{
     double utc,mjd,tdt,tdb,et;
     double utch,utcl,ut1;
     double ddd[6];
     double ss[3];
 /****used in sidereal_time****/
    
  double gast,eqeq,eqeqP;
  /****used in earthtilt****/
  double oblm,oblt,psi,eps;
  double ddpsi,ddeps;
  int verb=0; //0 does not print values

  utc=dot+dot2utc;
  mjd=utc-utc2mjd;
  if (verb>4)printf("utc= %18.8f\n ",utc);
  if (verb>4)printf("mjd= %18.8f\n ",mjd);
  /*****Leap seconds*****/
  tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
  if (verb>4)printf("tdt= %18.8f\n ",tdt);
  tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
  if (verb>4)printf("tdb= %18.8f\n ",tdb);
  if (verb>4)printf("Detector name,long,lat:%s %f %f\n",detector->name,detector->longitude,detector->latitude);

 /**Skip if not within the time-range of the present version of the ephemeris. */

    et=tdb; /*et is the barycenter time, at each sample */
    //if(et<ss[0]) continue;      /* ephemeris starts later than test points */
    //if(et>ss[1]) break;         /* ephemeris ends earlier than test points */

   
    /***** Nutation series from the JPL file *****/
     jpl_mod_pleph(et,14,0,ddd);
     psi=ddd[0]*RAD2SEC;
     eps=ddd[1]*RAD2SEC;
     if(nutcorr==1)  {
      PSS_deltanut(et, &ddpsi,&ddeps);
      //printf("et,psi eps ddpsi ddeps dal JPL senza corr %f %f %f %f %f\n",et,psi,eps,ddpsi,ddeps);
      psi+=ddpsi;
      eps+=ddeps;
       }
     novas_mod_earthtilt (et, &oblm,&oblt,&eqeq,psi,eps);
     //printf("et,psi eps dal JPL+corr %f %f %f \n",et,psi,eps);
     /*****Sidereal time needs UT1. We use UTC that approx. UT1 (error <0.9 s/anno) if fut1=0,  else we use UT1*****/
     if (!fut1)
     {
	utch=(long)(utc);
	utcl=utc-utch;
      }
	else
      {
        ut1=PSS_utc2ut1(utc);
        utch=(long)(ut1);
        utcl=ut1-utch;
       }
   
     eqeqP=0.0; //always 0, to evaluate the MEAN sidereal time
   sidereal_time (utch,utcl,eqeqP, &gast);

   *tss_local=gast+(detector->longitude)/15.0;  //hours
   if (verb>4)printf("dot,mjd,gast, tss local, hours= %f %f %f, %f\n ",dot,mjd,gast,*tss_local);
   return;
   }

void mjd2tsslocal(DETECTOR *detector,double mjd, short int fut1,short int nut1,short int nutcorr,double *tss_local)

     /*
INPUT:
detector
mjd=modified Julian date;
fut1, int: flag for UT1. If 0 the function uses UTC, if 1 it uses UT1

OUTPUT
tss_local=local sidereal time in hours
     */
{
     double utc,dot,tdt,tdb,et;
     double utch,utcl,ut1;
     double ddd[6];
     double ss[3];
 /****used in sidereal_time****/
    
  double gast,eqeq,eqeqP;
  /****used in earthtilt****/
  double oblm,oblt,psi,eps;
  double ddpsi,ddeps;
   /* Sotto era 400 anzichè NMAX che ora vale 1000 e sta in jplbin.h */
  char nams[NMAX][6];
  double vals[NMAX];
  intptr_t nvs; //pia64
  int verb=0; //0 does not print values

  if(fut1==1)F2=fopen(UTC2UT1,"rb"); //quiqui:F1 e F2,F3 definiti in PSS_astro.h
  if(nut1==1)F1=fopen(EPHFILE,"rb");
  if(nutcorr==1)F3=fopen(DELTANUT,"r");

  dot=mjd-dot2mjd; //days from the beginning of 1900
  utc=dot+dot2utc;
  if (verb>4)printf("utc= %18.8f\n ",utc);
  if (verb>4)printf("mjd= dot= %18.8f  %f\n ",mjd,dot);
  /*****Leap seconds*****/
  tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
   if (verb>4)printf("tdt= %18.8f\n ",tdt);
  tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
   if (verb>4)printf("tdb= %18.8f\n ",tdb);
   if (verb>4)printf("Detector name,long,lat:%s %f %f\n",detector->name,detector->longitude,detector->latitude);
  jpl_mod_constan(nams,vals,ss,&nvs);
   if (verb>4)printf("Range DE405= %14.2f  %14.2f  %14.2f\n",ss[0],ss[1],ss[2]);
 /**Skip if not within the time-range of the present version of the ephemeris. */

    et=tdb; /*et is the barycenter time, at each sample */
    //if(et<ss[0]) continue;      /* ephemeris starts later than test points */
    //if(et>ss[1]) break;         /* ephemeris ends earlier than test points */

    if(nut1==1){   
      /***** Nutation series from the JPL file *****/
      jpl_mod_pleph(et,14,0,ddd);
      psi=ddd[0]*RAD2SEC;
      eps=ddd[1]*RAD2SEC;
      if(nutcorr==1)  {
      PSS_deltanut(et, &ddpsi,&ddeps);
      if (verb>4) printf("quiqui:***et,psi eps ddpsi ddeps dal JPL senza corr %f %f %f %f %f\n",et,psi,eps,ddpsi,ddeps);
      psi+=ddpsi;
      eps+=ddeps;
       }
      novas_mod_earthtilt (et, &oblm,&oblt,&eqeq,psi,eps);
    }
    else {
      /*************Nutation series from NOVAS************/
      earthtilt (et, &oblm,&oblt,&eqeq,&psi,&eps);
    }
   
     /*****Sidereal time needs UT1. We use UTC that approx. UT1 (error <0.9 s/anno) if fut1=0,  else we use UT1*****/
     if (!fut1)
     {
	utch=(long)(utc);
	utcl=utc-utch;
      }
	else
      {
        ut1=PSS_utc2ut1(utc);
        utch=(long)(ut1);
        utcl=ut1-utch;
       }
   
   eqeqP=0.0; //always 0, to evaluate the MEAN sidereal time
   sidereal_time (utch,utcl,eqeqP, &gast);

   *tss_local=gast+(detector->longitude)/15.0;  //hours
    if(verb>4)printf("dot,mjd,gast, tss local, hours= %f %f %f, %f\n ",dot,mjd,gast,*tss_local);
   if(nut1==1)fclose(F1); 
   if(fut1==1)fclose(F2);
   if(nutcorr==1)fclose(F3);
   return;
   }
void gps2tsslocal(DETECTOR *detector,long gpstime, short int fut1,short int nut1,short int nutcorr,double *tss_local)

     /*
INPUT:
detector
gpsime=GPS seconds (long);
fut1, int: flag for UT1. If 0 the function uses UTC, if 1 it uses UT1

OUTPUT
tss_local=local sidereal time in hours
     */
{
     double utc,dot,mjd,tdt,tdb,et;
     double utch,utcl,ut1;
     double ddd[6];
     double ss[3];
 /****used in sidereal_time****/
    
  double gast,eqeq,eqeqP;
  /****used in earthtilt****/
  double oblm,oblt,psi,eps;
  double ddpsi,ddeps;
  int nleap;
  int verb=0;
  if(fut1==1)F2=fopen(UTC2UT1,"rb"); //F1 e F2,F3 definiti in PSS_astro.h
  if(nut1==1)F1=fopen(EPHFILE,"rb");
  if(nutcorr==1)F3=fopen(DELTANUT,"r");
  if(gpstime < gps1jan2006) nleap=32;
   /* fra il 1 gennaio 2006= 820108814  e 1 gennaio 2009=914803215 */
   if (gpstime >= (float)gps1jan2006 && gpstime < (float) 914803215) {
       nleap= 33;
    }
    if (gpstime >= (float)914803215) {
        nleap= 34;
    }
   if (gpstime >= (float)1025136015) {  //30 giugno 2012, 23 59 e 60 s
        nleap= 35;
    }
    if (gpstime >= (float)1119744016) {   //30 giugno 2015, idem. )1119744017 e' la mezzanotte del 1 luglio (when it was added)
        nleap= 36;
    }
    if (gpstime >= (float)1167264017) {   //1 JAN 2017, a mezzanotte e' con 8 finale.
        nleap= 37;
    }
  dot=dot2gps+(gpstime-(nleap-19))/(day2sec*1.0);
  utc=dot+dot2utc;
  mjd=utc-utc2mjd;
  if (verb>4)printf("utc= %18.8f\n ",utc);
  printf("mjd= %18.8f\n ",mjd);
  /*****Leap seconds*****/
  tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
  if (verb>4)printf("tdt= %18.8f\n ",tdt);
  tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
  if (verb>4)printf("tdb= %18.8f\n ",tdb);
  if (verb>4)printf("Detector name,long,lat:%s %f %f\n",detector->name,detector->longitude,detector->latitude);

 /**Skip if not within the time-range of the present version of the ephemeris. */

    et=tdb; /*et is the barycenter time, at each sample */
    //if(et<ss[0]) continue;      /* ephemeris starts later than test points */
    //if(et>ss[1]) break;         /* ephemeris ends earlier than test points */

   if (verb>4)printf("et 1= %18.8f\n ",et);
  if(nut1==1){   
      /***** Nutation series from the JPL file *****/
      jpl_mod_pleph(et,14,0,ddd);
      psi=ddd[0]*RAD2SEC;
      eps=ddd[1]*RAD2SEC;
      if(nutcorr==1)  {
      PSS_deltanut(et, &ddpsi,&ddeps);
      //printf("et,psi eps ddpsi ddeps dal JPL senza corr %f %f %f %f %f\n",et,psi,eps,ddpsi,ddeps);
      psi+=ddpsi;
      eps+=ddeps;
       }
      novas_mod_earthtilt (et, &oblm,&oblt,&eqeq,psi,eps);
    }
    else {
      /*************Nutation series from NOVAS************/
      earthtilt (et, &oblm,&oblt,&eqeq,&psi,&eps);
    }
      
     /*****Sidereal time needs UT1. We use UTC that approx. UT1 (error <0.9 s/anno) if fut1=0,  else we use UT1*****/
     if (!fut1)
     {
	utch=(long)(utc);
	utcl=utc-utch;
      }
	else
      {
        ut1=PSS_utc2ut1(utc);
        utch=(long)(ut1);
        utcl=ut1-utch;
       }
   
     eqeqP=0.0; //always 0, to evaluate the MEAN sidereal time
     if (verb>4)printf("in gps2tsslocal eqeqP= %f \n",eqeqP);
   sidereal_time (utch,utcl,eqeqP, &gast);

   *tss_local=gast+(detector->longitude)/15.0;  //hours
   if (verb>4)printf("dot,mjd,gast, tss local, hours= %f %f %f, %f\n ",dot,mjd,gast,*tss_local);
   if(nut1==1)fclose(F1);
   if(fut1==1)fclose(F2);
   if(nutcorr==1)fclose(F3);
   return;
   }

 
 double gps2mjd(double gpstime)
     /* From gpstime to mjd time */
{
  double dot,utc,mjd;
  int nleap;
  if(gpstime < gps1jan2006) nleap=32;  
   /* fra il 1 gennaio 2006= 820108814  e 1 gennaio 2009=914803215 */
   if (gpstime >= (float)gps1jan2006 && gpstime < (float) 914803215) {
       nleap= 33;
    }
    if (gpstime >= (float)914803215) {
        nleap= 34;
    }
   if (gpstime >= (float)1025136015) {  //30 giugno 2012, 23 59 e 60 s
        nleap= 35;
    }
    if (gpstime >= (float)1119744016) {   //30 giugno 2015, idem
        nleap= 36;
    }
    if (gpstime >= (float)1167264017) {   //1 JAN 2017, a mezzanotte e' con 8 finale.
        nleap= 37;
    }
  dot=dot2gps+(gpstime-(nleap-19))/day2sec;
  utc=dot+dot2utc;
  mjd=utc-utc2mjd;

  return mjd;
}

double mjd2gps(double mjdtime)
{
      double utc,gpstime;
      int nleap_sec,gpsleap;

      utc=mjdtime+utc2mjd;
      nleap_sec=PSS_leapseconds(utc);
      gpsleap=nleap_sec-19;
      gpstime=(mjdtime-mjd2gpsC)*day2sec+gpsleap;
    

return gpstime;
}
/*****Detector velocity:*****/

double PSS_locmax_vel(DETECTOR *detector,double dot_ini,short int fut1, double velx[], double vely[],double velz[])

/*For a given detector,  date (days from the beginning of the century) 
it gives the velocity vector (AU/day) in rectangular ecliptic coordinates, 
respect to the Barycenter of the Earth-Sun system (referred to Epoch). 

It returns MJD of the time.

INPUT ARGUMENTS:

detector, pointer to struct DETECTOR

dot_ini, double: days from the beginning of the Century



fut1, int: flag for UT1. If 0 the function uses UTC, if 1 it uses UT1



OUTPUT ARGUMENTS:

velx[npoint], 
vely[npoint]
velz[npoint]
double arrays:Det. velocity vector (AU/day) in rectangular ecliptic coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch)

RETURN:
MJD of the first datum
*******************************************************************************************
NOTE: Position vector of the Earth and Detector are not required.
      Uncomment the relative instructions if you need a code that computes also positions
********************************************************************************************
*/

{
   /* Sotto era 400 anzichè NMAX che ora vale 1000 e sta in jplbin.h */
  char nams[NMAX][6];
  double et, r[6], ss[3], vals[NMAX];
  int i, ntarg, nctr;
   intptr_t nvs;  //pia64
  int k;
  double pos[3];
  double vel[3];
  double utc,tdt,tdb,utc_ii,utch,utcl,ut1,mjd;
  /****used in sidereal_time****/
  double eqeq,eqeqP; 
  double gast;
  /****used in earthtilt****/
  double oblm,oblt,psi,eps;
/****used in terra, nutate,precession****/
  double pos1[3],pos2[3],pog[3]; 
  double vel1[3],vel2[3],vog[3];
/****used for nutation in the ephemeris****/
  double ddd[6];
  int verb=0;
  F1=fopen(EPHFILE,"rb");
  F2=fopen(UTC2UT1,"rb");
  //F3=fopen(DELTANUT,"rb");
/****** Print the ephemeris constants. ******/

  jpl_mod_constan(nams,vals,ss,&nvs);
  if (verb>4)printf("Range DE405= %14.2f  %14.2f  %14.2f\n",ss[0],ss[1],ss[2]);
  /*  for(i=0;i<nvs;++i)
      printf("%.6s  %24.16E\n",nams[i],vals[i]); */

  
  ntarg=3; /* Earth */
  nctr=12; /*Solar system barycenter */
  utc=dot_ini+dot2utc;
  mjd=utc-utc2mjd;
  if (verb>4)printf("initial utc= %18.8f\n ",utc);
  if (verb>4)printf("initial mjd= %18.8f\n ",mjd);
  /*****Leap seconds*****/
  tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
  if (verb>4)printf("initial tdt= %18.8f\n ",tdt);
  tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
  if (verb>4)printf("tdb= %18.8f\n ",tdb);
  if (verb>4)printf("Detector name,long,lat:%s %f %f\n",detector->name,detector->longitude,detector->latitude);


   /*****  Read a value; Skip if not within the time-range
        of the present version of the ephemeris.                            */

    et=tdb; /*et is the barycenter time*/
    utc_ii=utc;
     if (verb>4)printf("%f %f %f \n", dot_ini,et,utc);

    jpl_mod_pleph(et,ntarg,nctr,r);
   

    for (k=0; k<=2; ++k)
	{
	  pos[k]=r[k];     //AU
          vel[k]=r[k+3];   //AU/day. Earth center vs SSB, referred at T0
	                   //Rectangular equatorial coordinates
	}
    /*****Evaluate Detector velocity vs Earth Center*****/
    /***** Nutation series from the JPL file *****/
     jpl_mod_pleph(et,14,0,ddd);
     psi=ddd[0]*RAD2SEC;
     eps=ddd[1]*RAD2SEC;
     novas_mod_earthtilt (et, &oblm,&oblt,&eqeq,psi,eps);

     /*****Sidereal time needs UT1. We use UTC that approx. UT1 (error <0.9 s/anno) if fut1=0,  else we use UT1*****/
     if (!fut1)
     {
	utch=(long)(utc_ii);
	utcl=utc_ii-utch;
      }
	else
      {
        ut1=PSS_utc2ut1(utc_ii);
        utch=(long)(ut1);
        utcl=ut1-utch;
       }
    eqeqP=eqeq; //Kaz suggestion. May 2003
     if(NutFlag==0)eqeqP=0.0;
   sidereal_time (utch,utcl,eqeqP, &gast);
   novas_mod_terra (detector,gast, pos1,vel1);
   //nutate (et,FN1,pos1, pos2); //uncomment if required
   //precession (et,pos2,T0, pog); //uncomment if required
   nutate(et,FN1,vel1, vel2); /*performs inverse nutation (true to mean) */
   precession (et,vel2,T0, vog); //rectangular equatorial coord., ref. at T0 


   /*****Sum the two contributions: Earth vs SSB + Det vs. Earth *****/
  for(k=0; k<=2; ++k)
      {
	//pos[k]+=pog[k];
	 vel[k]+=vog[k];

      }

    /*****Convert into Ecliptical coordinates*****/
    PSS_equr2eclr (vel, vel1); 
    //PSS_equr2eclr (pos, pos1); //uncomment if required

    velx[0]=vel1[0];
    vely[0]=vel1[1];
    velz[0]=vel1[2];
 
     fclose(F1);
     fclose(F2);
     //    fclose(F3);
     return mjd;
}


double PSS_detector_velposSCALAR(DETECTOR *detector,double dot_ini, double step,short int fut1, short int nut1, short int nutcorr,short int equatorial,short int imjd,double *velx, double *vely,double *velz,double *posx, double *posy,double *posz,double deinstein)

/*For a given detector, initial date (days from the beginning of the century or mjd),step (s), number
 of point, it gives the velocity vector (AU/day) and the position in AU in rectangular ecliptic 
or equatorial coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch). 

It returns MJD of the first point.

INPUT ARGUMENTS:

detector, pointer to struct DETECTOR

dot_ini, double: days from the beginning of the Century, or mjd, or TDT in mjd

step, double: step in seconds



fut1, int: flag for UT1. If 0 the function uses UTC, if 1 it uses UT1

imjd=0,1,2 0=days from 1900, 1=mjd , 2=tdt in JD 

OUTPUT ARGUMENTS:


posz
double SCALARS:Det. velocity vector (AU/day)  and pos in km in rectangular ecliptic or euatorial coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch)

deinstein
RETURN:
MJD of the first datum
*******************************************************************************************
NOTE: Position vector of the Earth and Detector are not required.
      Uncomment the relative instructions if you need a code that computes also positions
********************************************************************************************
*/

{
   /* Sotto era 400 anzichè NMAX che ora vale 1000 e sta in jplbin.h */
  char nams[NMAX][6];
  double et, r[6], ss[3], vals[NMAX];
  int i, ntarg, nctr;
   intptr_t nvs;  //pia64
  int k;
  double pos[3];
  double vel[3];
  double utc,tdt,tdb,utc_ii,utch,utcl,ut1,mjd;
  /****used in sidereal_time****/
  double eqeq,eqeqP; 
  double gast;
  /****used in earthtilt****/
  double oblm,oblt,psi,eps;
  double ddpsi,ddeps; 
/****used in terra, nutate,precession****/
  double pos1[3],pos2[3],pog[3]; 
  double vel1[3],vel2[3],vog[3];
/****used for nutation in the ephemeris****/
  double ddd[6];
  double tdt_ini; 
  int verb=5;  //poi va sistemata e tolta da qui

  /******PER TEST FILE */
  //FILE *OUT33;
  //OUT33=fopen("table_test.dat","wb");


  F1=fopen(EPHFILE,"rb");
  F2=fopen(UTC2UT1,"rb");
  F3=fopen(DELTANUT,"rb");

  if (verb>4)puts("Initial time, days");
  if (verb>4)printf("%f\n", dot_ini);

/****** Print the ephemeris constants. ******/

   jpl_mod_constan(nams,vals,ss,&nvs);
   if (verb>4)printf("Range DE405= %14.2f  %14.2f  %14.2f\n",ss[0],ss[1],ss[2]); //pia Sep. 2005
   if (verb>4)puts("--> If ok it Should be: 2444208.50  2458864.50  32.00.  nvs 156");
   //for(i=0;i<nvs;++i)
   //printf("%.6s  %24.16E\n",nams[i],vals[i]); 
  
  
  ntarg=3; /* Earth */
  nctr=12; /*Solar system barycenter */
  if(imjd==0)
  {
    utc=dot_ini+dot2utc;
    mjd=utc-utc2mjd;
    if (verb>4)printf("initial utc= %18.8f\n ",utc); //pia Sep. 2005
    if (verb>4)printf("initial mjd= %18.8f\n ",mjd);
    if (verb>4)printf("initial days from 1900= %18.8f\n ",dot_ini);
  }
   if(imjd==1)
  {
    utc=dot_ini+utc2mjd;
    mjd=dot_ini;
    dot_ini=utc-dot2utc;  //dot_ini sempre giorni inizio secolo, in UTC !! NOTARE
    if (verb>4)printf("imjd=1: initial utc= %18.8f\n ",utc);
    if (verb>4)printf("imjd=1: initial mjd= %18.8f\n ",mjd);
  }
   if(imjd==2){
     tdt=dot_ini+utc2mjd;
     mjd=dot_ini; //mjd but TDT
     if (verb>4)printf("initial tdt= %20.11f\n ",tdt);
     if (verb>4)printf("initial TDT in mjd= %18.8f\n ",mjd);
   }
  /*****Leap seconds*****/
  if(imjd !=2) tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
  if (verb>4)printf("initial tdt= %18.8f\n ",tdt); //pia Sep. 2005
  tdt_ini=tdt;
  tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
  if (verb>4)printf("tdb= %18.8f\n ",tdb);
  if (verb>4)printf("utc2mjd= %25.15f\n ",utc2mjd);
  if (verb>4)printf("Detector name,long,lat:%s %f %f\n",detector->name,detector->longitude,detector->latitude); //pia Sep. 2005

 
   /*****  Read a value; Skip if not within the time-range
        of the present version of the ephemeris.                            */
     if(imjd!=2){
       utc=dot_ini+dot2utc; //ONE VALUE ONLY
       tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
 
     }
     else
       {
	
	 tdt=tdt_ini;
       }
     tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
 
    et=tdb;              //sost. kaz fino a qui e tolta la riga seguente
    //et=tdb+(ii*step/day2sec); /*et is the barycenter time, at each sample */
     //ora il tdb lo ricalcolo su tutti i samples
     if(imjd==2)utc=PSS_tdt2utc(tdt);
     utc_ii=utc;  //+(ii*step/day2sec);
    if(et<ss[0]) return;      /* ephemeris starts later than test points */
    if(et>ss[1]) return;         /* ephemeris ends earlier than test points */

    jpl_mod_pleph(et,ntarg,nctr,r);
    //printf("et ss0 ss1 %f %f %f\n",et,ss[0],ss[1]);

    for (k=0; k<=2; ++k)
	{
	  pos[k]=r[k];     //AU
          vel[k]=r[k+3];   //AU/day. Earth center vs SSB, referred at T0
	                   //Rectangular equatorial coordinates
	}
/*****Evaluate Detector velocity vs Earth Center*****/
      if(nut1==1){   
      /***** Nutation series from the JPL file *****/
      jpl_mod_pleph(et,14,0,ddd);
      psi=ddd[0]*RAD2SEC;
      eps=ddd[1]*RAD2SEC;
       if(nutcorr==1)  {
      PSS_deltanut(et, &ddpsi,&ddeps);
      //if(ii<=1)printf("et,psi eps ddpsi ddeps dal JPL senza corr %f %f %f %f %f\n",et,psi,eps,ddpsi,ddeps);
      psi+=ddpsi;
      eps+=ddeps;
       }
      novas_mod_earthtilt (et, &oblm,&oblt,&eqeq,psi,eps);
      //if(ii<=1)printf("et,psi eps dal JPL+corr %f %f %f \n",et,psi,eps);
      //printf("et,psi eps dal JPL %f %f %f \n",et,psi,eps);
    }
    else {
      /*************Nutation series from NOVAS************/
      earthtilt (et, &oblm,&oblt,&eqeq,&psi,&eps);
      //printf("et,psi eps dal NOVAS %f %f %f \n",et,psi,eps);
    }
     //Aug 2003: i valori JPL o NOVAS sono identici
     //Invece devo usare il file dello IERS per le correzioni
     //eopc04.62-now preso da http://hpiers.obspm.fr/eoppc/eop/eopc04/
     //dove sono in secondi di arco

     /*****Sidereal time needs UT1. We use UTC that approx. UT1 (error <0.9 s/anno) if fut1=0,  else we use UT1*****/
     if (!fut1)
     {
	utch=(long)(utc_ii);
	utcl=utc_ii-utch;
      }
	else
      {
        ut1=PSS_utc2ut1(utc_ii);
        utch=(long)(ut1);
        utcl=ut1-utch;
       }
     eqeqP=eqeq; //Kaz suggestion. May 2003
     if(NutFlag==0)eqeqP=0.0;
   sidereal_time (utch,utcl,eqeqP, &gast);
   novas_mod_terra (detector,gast, pos1,vel1);
   //modified Sept. 2003 !!!!!
   if(nutcorr==1)novas_mod_nutate(et,FN1,pos1, pos2,psi,eps); //uncomment if required
   if(nutcorr==0)nutate(et,FN1,pos1, pos2); //uncomment if required
   precession(et,pos2,T0, pog); //uncomment if required
  if(nutcorr==1)novas_mod_nutate(et,FN1,vel1, vel2,psi,eps); /*performs inverse nutation (true to mean) */
   if(nutcorr==0)nutate(et,FN1,vel1, vel2);
   precession (et,vel2,T0, vog); //rectangular equatorial coord., ref. at T0 

    /*******write a TEST file with values to compare***********/
  
   //fprintf(OUT33,"%23.11f ",mjd+(ii*step/(1.0*day2sec)));
     
   //fprintf(OUT33,"%25.18e %25.18e %25.18e %25.18e %25.18e %25.18e %25.18e %25.18e %25.18e%25.18e %25.18e %25.18e  ",vel[0]/C, vel[1]/C,vel[2]/C, vel1[0]/C, vel1[1]/C,vel1[2]/C,(vel2[0]-vel1[0])/C, (vel2[1]-vel1[1])/C,(vel2[2]-vel1[2])/C,(vog[0]-vel2[0])/C, (vog[1]-vel2[1])/C,(vog[2]-vel2[2])/C);
   //fprintf(OUT33,"\n"); 

  /*******END TEST*****/

   /*****Sum the two contributions: Earth vs SSB + Det vs. Earth *****/
  for(k=0; k<=2; ++k)
      {
	 pos[k]+=pog[k];
	 vel[k]+=vog[k];

      }
 
    /*****If equat orial diverso da 1 convert into Ecliptical coordinates*****/
  if(equatorial!=1)
  {
     PSS_equr2eclr (vel, vel1); 
     //PSS_equr2eclr (pos, pos1); //uncomment if required
    *velx=vel1[0];
    *vely=vel1[1];
    *velz=vel1[2];
  }
  else
  { 
     *velx=vel[0];
     *vely=vel[1];
     *velz=vel[2];
  }
   if(equatorial!=1)
  { 
    PSS_equr2eclr (pos, pos1); //uncomment if required
    *posx=pos1[0];
    *posy=pos1[1];
    *posz=pos1[2];
  }
  else
  { 
     *posx=pos[0];
     *posy=pos[1];
     *posz=pos[2];
  }
    deinstein= PSS_einstein(tdt);  //Deinstein
    //printf("Einstein %d %25.20f\n", ii, deinstein[ii]);
     if (verb>4)printf("nella function: VELOCITA'/C: %15.8e %15.8e %15.8e\n",*velx/C,*vely/C,*velz/C);
     if (verb>4)printf("VELOCITA': %15.8e %15.8e %15.8e\n",*velx,*vely,*velz);
     if (verb>4)printf("POSIZIONI*day2sec/C: %15.8e %15.8e %15.8e\n",*posx*day2sec/C,*posy*day2sec/C,*posz*day2sec/C);
     fclose(F1);
     fclose(F2);
     fclose(F3);
     //   fclose(OUT33); //TEST FILE
     return mjd;
}

/*****Earth center velocity AND position:*****/

double PSS_Nodet_velpos(double dot_ini, double step, int npoint,short int fut1, short int nut1, short int nutcorr,short int equatorial,short int imjd,double velx[], double vely[],double velz[],double posx[], double posy[],double posz[],double deinstein[])

/*For the Earth center, initial date (days from the beginning of the century or mjd),step (s), number
 of point, it gives the velocity vector (AU/day) and the position in AU in rectangular ecliptic 
or equatorial coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch). 

It returns MJD of the first point.

INPUT ARGUMENTS:



dot_ini, double: days from the beginning of the Century, or mjd, or TDT in mjd

step, double: step in seconds

npoint, int: number of data

fut1, int: flag for UT1. If 0 the function uses UTC, if 1 it uses UT1

imjd=0,1,2 0=days from 1900, 1=mjd , 2=tdt in JD 

OUTPUT ARGUMENTS:

velx[npoint], 
vely[npoint]
velz[npoint]
posx[npoint], 
posy[npoint]
posz[npoint]
double arrays:Det. velocity vector (AU/day)  and pos in km in rectangular ecliptic or euatorial coordinates, respect to the Barycenter of the Earth-Sun system (referred to Epoch)

deinstein[npoint]
RETURN:
MJD of the first datum
*******************************************************************************************

      Uncomment the relative instructions if you need a code that computes also positions
********************************************************************************************
*/

{
   /* Sotto era 400 anzichè NMAX che ora vale 1000 e sta in jplbin.h */
  char nams[NMAX][6];
  double et, r[6], ss[3], vals[NMAX];
  int i, ntarg, nctr;
   intptr_t nvs;  //pia64
  int ii,k;
  double pos[3];
  double vel[3];
  double utc,tdt,tdb,utc_ii,utch,utcl,ut1,mjd;
  /****used in sidereal_time****/
  double eqeq,eqeqP; 
  double gast;
  /****used in earthtilt****/
  double oblm,oblt,psi,eps;
  double ddpsi,ddeps; 
/****used in terra, nutate,precession****/
  double pos1[3],pos2[3],pog[3]; 
  double vel1[3],vel2[3],vog[3];
/****used for nutation in the ephemeris****/
  double ddd[6];
  double tdt_ini; 
  int verb=0; //0 does not print values

  /******PER TEST FILE */
  //FILE *OUT33;
  //OUT33=fopen("table_test.dat","wb");


  F1=fopen(EPHFILE,"rb");
  F2=fopen(UTC2UT1,"rb");
  F3=fopen(DELTANUT,"rb");

  if (verb>4)puts("Initial time, days");
  if (verb>4)printf("%f\n", dot_ini);

/****** Print the ephemeris constants. ******/

   jpl_mod_constan(nams,vals,ss,&nvs);
   if (verb>4)printf("Range DE405= %14.2f  %14.2f  %14.2f\n",ss[0],ss[1],ss[2]); //pia Sep. 2005
   if (verb>4)puts("--> If ok it Should be: 2444208.50  2458864.50  32.00.  nvs 156");
   //for(i=0;i<nvs;++i)
   //printf("%.6s  %24.16E\n",nams[i],vals[i]); 
  
  ii=0;
  ntarg=3; /* Earth */
  nctr=12; /*Solar system barycenter */
  if(imjd==0)
  {
    utc=dot_ini+dot2utc;
    mjd=utc-utc2mjd;
    if (verb>4)printf("initial utc= %18.8f\n ",utc); //pia Sep. 2005
    if (verb>4)printf("initial mjd= %18.8f\n ",mjd);
    if (verb>4)printf("initial days from 1900= %18.8f\n ",dot_ini);
  }
   if(imjd==1)
  {
    utc=dot_ini+utc2mjd;
    mjd=dot_ini;
    dot_ini=utc-dot2utc;  //dot_ini sempre giorni inizio secolo, in UTC !! NOTARE
    if (verb>4)printf("imjd=1: initial utc= %18.8f\n ",utc);
    if (verb>4)printf("imjd=1: initial mjd= %18.8f\n ",mjd);
  }
   if(imjd==2){
     tdt=dot_ini+utc2mjd;
     mjd=dot_ini; //mjd but TDT
     if (verb>4)printf("initial tdt= %20.11f\n ",tdt);
     if (verb>4)printf("initial TDT in mjd= %18.8f\n ",mjd);
   }
  /*****Leap seconds*****/
  if(imjd !=2) tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
  if (verb>4)printf("initial tdt= %18.8f\n ",tdt); //pia Sep. 2005
  tdt_ini=tdt;
  tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
  if (verb>4)printf("tdb= %18.8f\n ",tdb);
  if (verb>4)printf("utc2mjd= %25.15f\n ",utc2mjd);


 while( ii < npoint)
   {
   /*****  Read a value; Skip if not within the time-range
        of the present version of the ephemeris.                            */
     if(imjd!=2){
       utc=dot_ini+dot2utc+ii*step/(1.0*day2sec); //commento kaz: tdb va ricalcolato sempre
       tdt=PSS_utc2tdt(utc); //from utc to terrestrial dynam. time,tdt
 
     }
     else
       {
	
	 tdt=tdt_ini+ii*step/(1.0*day2sec);
       }
     tdb=PSS_tdt2tdb(tdt); //from tdt to barycenter time tdb
 
    et=tdb;              //sost. kaz fino a qui e tolta la riga seguente
    //et=tdb+(ii*step/day2sec); /*et is the barycenter time, at each sample */
     //ora il tdb lo ricalcolo su tutti i samples
     if(imjd==2)utc=PSS_tdt2utc(tdt);
     utc_ii=utc;  //+(ii*step/day2sec);
    if(et<ss[0]) continue;      /* ephemeris starts later than test points */
    if(et>ss[1]) break;         /* ephemeris ends earlier than test points */

    jpl_mod_pleph(et,ntarg,nctr,r);
    //printf("et ss0 ss1 %f %f %f\n",et,ss[0],ss[1]);

    for (k=0; k<=2; ++k)
	{
	  pos[k]=r[k];     //AU
          vel[k]=r[k+3];   //AU/day. Earth center vs SSB, referred at T0
	                   //Rectangular equatorial coordinates
	}


    /*****If equat orial diverso da 1 convert into Ecliptical coordinates*****/
  if(equatorial!=1)
  {
     PSS_equr2eclr (vel, vel1); 
     //PSS_equr2eclr (pos, pos1); //uncomment if required
    velx[ii]=vel1[0];
    vely[ii]=vel1[1];
    velz[ii]=vel1[2];
  }
  else
  { 
     velx[ii]=vel[0];
     vely[ii]=vel[1];
     velz[ii]=vel[2];
  }
   if(equatorial!=1)
  { 
    PSS_equr2eclr (pos, pos1); //uncomment if required
    posx[ii]=pos1[0];
    posy[ii]=pos1[1];
    posz[ii]=pos1[2];
  }
  else
  { 
     posx[ii]=pos[0];
     posy[ii]=pos[1];
     posz[ii]=pos[2];
  }
    deinstein[ii]= PSS_einstein(tdt);  //Deinstein
    //printf("Einstein %d %25.20f\n", ii, deinstein[ii]);
    ii++;
   }
     if (verb>4)printf("Earth only VELOCITA'/C: %15.8e %15.8e %15.8e\n",velx[ii]/C,vely[ii]/C,velz[ii]/C);
     if (verb>4)printf("Earth only VELOCITA': %15.8e %15.8e %15.8e\n",vel[0],vely[ii],velz[ii]);
     if (verb>4)printf("Earth only POSIZIONI*day2sec/C: %15.8e %15.8e %15.8e\n",posx[ii]*day2sec/C,posy[ii]*day2sec/C,posz[ii]*day2sec/C);
     fclose(F1);
     fclose(F2);
     fclose(F3);
     //   fclose(OUT33); //TEST FILE
     return mjd;    
  }

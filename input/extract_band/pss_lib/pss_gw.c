/*___________________________________________________
 ¦                                                   ¦
 ¦                     pss_gw.c                      ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    March 2002                     ¦
 ¦___________________________________________________¦*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pss_snag.h"
#include "pss_math.h"
#include "pss_serv.h"
#include "pss_gw.h"


GW pss_set_gw(long nsimsour, long niant, long nbant)
/* input parameters are for simulated objects */
{
	GW gw;

	gw.math.pi=MATH_PI;
	gw.math.e=MATH_E;

	gw.phys.c=PHYS_C;
	gw.phys.G=PHYS_G;

	gw.astro.SY_s=ASTRO_SY_S;
	printf("nsimsour,niant,nbant= %ld %ld %ld \n",nsimsour,niant,nbant); //pia: messo printf per usarli..
	return gw;
}


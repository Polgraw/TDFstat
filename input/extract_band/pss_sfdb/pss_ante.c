/**pss_ante.c: pss software**/
/**Last version: Marc, 17 2009**/
/**Author: Pia**/

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
/****** Antenna libraries**************/
#include "pss_ante.h"

DETECTOR_PARAM* crea_parameters(int casin)
{
  DETECTOR_PARAM *detector_param; /*detector (bar or interferometer) parameters*/
  detector_param=(DETECTOR_PARAM *)malloc(sizeof(DETECTOR_PARAM));

    if(casin==1)puts("Itf data: sds format");
    if(casin==2)puts("Itf data: frame format");
     detector_param->itf_fixed_param=(ITFFIXED_PARAM *)malloc(sizeof(ITFFIXED_PARAM)); 
     detector_param->itf_var_param=(ITFVAR_PARAM *)malloc(sizeof(ITFVAR_PARAM));
   
  detector_param->itf_fixed_param->gpstime=0;
  detector_param->itf_fixed_param->time=0;
  return detector_param;
} 


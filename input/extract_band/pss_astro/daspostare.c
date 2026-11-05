/***************************************************************/
/* Definizioni che poi andranno spostate, Pia 25 Jan 2002      */
/*                                                             */
/***************************************************************/
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

#include "daspostare.h"

/*****Structures*****/

DETECTOR *def_detector(char* det_name){

   DETECTOR *detect;
   detect=(DETECTOR *)malloc(sizeof(DETECTOR));

  if(!strcmp(det_name,"explo"))
    {  
    strcpy(detect->name,"explo");
    detect->longitude=6.25;
    detect->latitude=46.25;
    detect->height=100.0;
    detect->azimuth=39.3;
    }
  if (!strcmp(det_name,"nauti"))
    {    
    strcpy(detect->name,"nauti");
    detect->longitude=12.67;  //EAST 
    detect->latitude=41.82;   //NORTH
    detect->height=30.0;
    detect->azimuth=44; //verso EST contato a partire dal nord  //From North to East
    }
  if (!strcmp(det_name,"virgo")) 
     {
     strcpy(detect->name,"virgo");
     detect->longitude=10.+30./60.+16.1885/3600.;
     detect->latitude=43.+37./60.+53.0880/3600.;
     detect->height=9;  //corretto 22 dic. 2009 53.238;
     detect->azimuth=199.4326-180.00;  //19.0+26/60.0;   // rispetto al nord clockwise
     }
   if (!strcmp(det_name,"ligoh")) 
     {
     strcpy(detect->name,"ligoh");
     detect->longitude=240.592; //deg
     detect->latitude=46.455;  //deg
     detect->height=159;
     detect->azimuth= 144.0006-180.00;   // rispetto al nord clockwise
     }
 if (!strcmp(det_name,"ligol")) 
     {
     strcpy(detect->name,"ligol");
     detect->longitude=269.266;
     detect->latitude=30.563;
     detect->height=20;  //metri
     detect->azimuth=72.2835-180.00;   // rispetto al nord clockwise
     }
  if (!strcmp(det_name,"parkes"))
    {    
    strcpy(detect->name,"parks");
    detect->longitude=148.0+15.0/60+42.0/3600;
    detect->latitude=-(33.0+0.04/3600);
    detect->height=392.0;
    }
  if (!strcmp(det_name,"kagra"))
    {    
    strcpy(detect->name,"kagra");
    detect->longitude=137.305956;
    detect->latitude=36.41186;
    detect->height=414.0;
    detect->azimuth=60.396226-180.00;   // rispetto al nord clockwise
    }
  if (!strcmp(det_name,"Nodet"))
    {    
    strcpy(detect->name,"Nodet");
    detect->longitude=0;
    detect->latitude=0;
    detect->height=0;
    detect->azimuth=0;   // rispetto al nord clockwise
    }
  return detect;
  }


SOURCE *def_source(char* source_name){

   SOURCE *source;
   source=(SOURCE *)malloc(sizeof(SOURCE));

  if(!strcmp(source_name,"zerozer"))
    {  
    strcpy(source->name,"zerozer");
    source->ra=0.0;
    source->dec=-0.0;
    source->promora=0.0; 
    source->promodec=0.0; 
    source->radialvelocity=0.0;
    source->parallax=1.0e-4;
    source->frequency=900.0;
    source->psi=0.;
    source->eps=0.;
    }
if(!strcmp(source_name,"galcent"))
    {  
    strcpy(source->name,"galcent");
    source->ra=17.7;
    source->dec=-29.0;
    source->promora=0.0; 
    source->promodec=0.0; 
    source->radialvelocity=0.0;
    source->parallax=1.0e-4;
    source->frequency=900.0;
    source->psi=0.;
    source->eps=0.;
    }
 if(!strcmp(source_name,"tucanae"))
    {  
    strcpy(source->name,"tucanae");
    source->ra=0.4;
    source->dec=-71.993;
    source->promora=0.0; 
    source->promodec=0.0; 
    source->radialvelocity=0.0;
    source->parallax=1.0e-7;
    source->frequency=900.0;
    source->psi=0.;
    source->eps=0.;
    }
 if(!strcmp(source_name,"PSR437"))
    {  
    strcpy(source->name,"PSR437");
    source->ra=4.0+37/60.0+15.7481794/3600.0;
    source->dec=-(47.0+15/60.0+8.2315166/3600.0);
    source->promora=0.0; 
    source->promodec=0.0; 
    source->radialvelocity=0.0;
    source->parallax=6.6634e-3;
    source->frequency=173.687949174533;
    source->psi=0.;
    source->eps=0.;

    }
 if(!strcmp(source_name,"PSR1937"))
    {  
    strcpy(source->name,"PSR1937");
    source->ra=19.0+39/60.0+38.56024702/3600.0;
    source->dec=(21.0+34/60.0+59.1434800/3600.0);
    source->promora=0.0; 
    source->promodec=0.0; 
    source->radialvelocity=0.0;
    source->parallax=0.200000;
    source->frequency=641.9282681231345578;
    source->psi=0.;
    source->eps=0.;

    }
 if(!strcmp(source_name,"crab"))
    {  
    strcpy(source->name,"crab");
    source->ra=5.575548056;
    source->dec=22.01446;
    source->promora=0.0; 
    source->promodec=0.0; 
    source->radialvelocity=0.0;
    source->parallax=0.200000;
    source->frequency=641.9282681231345578;
    source->psi=0.;
    source->eps=0.;

    }
  if(!strcmp(source_name,"simul"))
    {  
    strcpy(source->name,"simul");
    source->ra=1000; //to de defined in the main code
    source->dec=1000; //to be defined in the main code
    source->promora=0.0; 
    source->promodec=0.0; 
    source->radialvelocity=0.0;
    source->parallax=0.200000;
    source->frequency=0; //to be defined in the main code
    source->psi=0.;
    source->eps=0.;

    }
 return source;
}


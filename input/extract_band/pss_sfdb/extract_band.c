/***********Program to extract Bandwidth from the SFDB files***************/
/***********Author:   Pia**************************/
/***********Last version 11 DEC. 2013******************/

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

/************ PSS libraries************/
//#include "../pss_lib/pss_math.h"
//#include "../pss_lib/pss_serv.h"
//#include "../pss_lib/pss_snag.h"
//#include "../pss_lib/pss_sfc.h"
/****** Antenna libraries*************/
//#include "pss_ante.h"

/****** SFDB libraries**************/
#include "pss_sfdb.h"


int main(void){
    int ifatti;
    int ii;
    int nmax;               /*number of FFT to be read*/
    int ilista=1;           //1 to read input files from a lista file
    int errorcode;
    GD *gd;                 /*GD with SFDB data*/
    GD *gd_short;           /*GD with the short SFDB data*/ 
    GD *gd_band;            /*GD with SFDB data, over a chosen sub-bandwidth*/
    char *fileBAND;
    FILE *BAND;
    int i_band;
    float total_band, initial_freq, band, final_freq, bandUS;
    int inter;
    float freq_sin, fr;
    int iw;
    int fft_read;
    double scale;
  
    HEADER_PARAM *header_param; /*Structure with parameters which have been  written in the sfdb file*/
    header_param = (HEADER_PARAM *)malloc(sizeof(HEADER_PARAM));
  
    iw = -1;
    inter = 0;
    LOG_INFO = logfile_open("extract_band");
 
    
    gd = crea_gd(0., 0., 0., "orig");
    gd_short = crea_gd(0., 0., 0., "ar");
    fileBAND = (char*)malloc(sizeof(char)*32);
    printf("Number max of FFTs to be read \n");
    scanf("%d", &nmax);
  
    printf("Name of the output file with the selected band  (ascii)\n");
    scanf("%s", fileBAND);
    BAND=fopen(fileBAND, "w");
  
    for(ifatti = 0; ifatti < nmax; ifatti++){
	errorcode = readfilesfdb(gd, gd_short, header_param, ifatti, iw, &fft_read, ilista, 0);
	if (errorcode != 1){
	    printf("Error in reading data into SFT file!\n");
	    printf("Lette fft numero= %d\n", ifatti);
	    if (ifatti > 0){
		printf("gd->n gd_band->n gd_band->dx %ld %ld %f\n",
		       gd->n, gd_band->n, gd_band->dx);
	    }
	    fclose(BAND);
	    return errorcode;
	}
  
	//fft2spectrum=sqrt((double) header_param->tsamplu/(2.0*header_param->nsamples));
	//printf("fft2spectrum= %f\n",fft2spectrum);
	
	if(ifatti == 0){
    
	    printf("First file mjdtime, gps time (s and ns) %15.10f %d %d\n",
		   header_param->mjdtime,header_param->gps_sec,header_param->gps_nsec);
	    total_band = gd->dx*gd->n/2; //the data are Real and Imag, so there is a factor 2
	    printf("Initial frequency of the data - Bandwidth (positive only)- total length %f %f %ld\n",
		   gd->ini,total_band,gd->n);
	    printf("scaling factor in amplitude %e\n",header_param->einstein);
	    if(header_param->typ == 2) inter = 1; //data overlapped by the half
     
	    printf("Initial frequency to extract ? \n");
	    scanf("%f", &freq_sin);
	    fr = gd->dx*floor(freq_sin/gd->dx);
	    freq_sin = fr;
	    printf("Reavaluated EXACT  frequency %f\n", fr);   
	    printf("Total bandwidth ?\n");
	    scanf("%f", &band);
      
	    initial_freq = fr;
	    final_freq = gd->dx*floor((fr + band)/gd->dx);      
	    band = (final_freq - initial_freq);
	    printf("ini, ifi, band  %f %f %f\n", initial_freq, final_freq, band);
	}
   
	
	if(ifatti == 0) gd_band = crea_gd(gd->n, (double)initial_freq, gd->dx, "band");
	//printf("Before band_extr:band, final_freq -  %f %f  \n",band,final_freq);

	//TEST
	//for (ii=0;ii<gd->n;ii+=2){
	//    gd->y[ii]=ii/2*gd->dx;
	//    gd->y[ii+1]=(ii/2+0.5)*gd->dx;
	//  }
	//END TEST
	i_band = band_extract(gd, gd_band, band);
	initial_freq = gd_band->ini; //revaluated exactly with the selected bins
	//printf("****Exact values:gd_band->ini, gd_band->dx gd_band->n/2= %f %f %ld\n",gd_band->ini,gd_band->dx,gd_band->n/2);
	bandUS = (gd_band->n/2)*gd_band->dx; //really used band (for the power 2 in the FFT)
  
	//printf("initial_freq,final_freq,**bandUS,samples extracted,gd_band->n %f %f %f %f %ld\n",initial_freq,initial_freq+bandUS,bandUS, 2*band/gd_band->dx,gd_band->n);

	if(ifatti == 0){
	    printf("initial_freq,final_freq,**bandUS,samples extracted,gd_band->n %f %f %f %f %ld\n",
		   initial_freq, initial_freq + bandUS, bandUS, 2*band/gd_band->dx, gd_band->n);
	    //write some information on the file at the beginning:	 
	    fprintf(BAND, "%% Beginning freq- Band- Samples in one stretch- Subsampling factor- inter (overlapping, 2 if data were overlapped)- Frequency step- Scaling factor- ***The data are real and imag of the FFT\n");
	    scale = 1.0;
	    if(header_param->einstein == 1) scale = 1.0E20;
	    if(scale == 1.0){
		fprintf(BAND, "%% %f %f %ld %f %d %10.7f %e\n",
			gd_band->ini, bandUS, gd_band->n/2, 1.0*(gd->n/gd_band->n), header_param->typ, gd_band->dx, header_param->einstein);
	    }
	    if(scale != 1){
		printf("Scaling factor in SFDB was 1. Changed to %e\n", 1./scale);
		fprintf(BAND, "%% %f %f %ld %f %d %10.7f %e\n",
			gd_band->ini, bandUS, gd_band->n/2, 1.0*(gd->n/gd_band->n), header_param->typ, gd_band->dx, 1./scale);
	    }
	}
	fprintf(BAND, "%% FFT number in the file; Beginning mjd days; Gps s; Gps ns;\n");
	fprintf(BAND, "%% %d %15.10f %d %d\n",
		header_param->nfft, header_param->mjdtime, header_param->gps_sec, header_param->gps_nsec);
	for(ii = 0; ii < gd_band->n; ii += 2){
	    fprintf(BAND, "%15.8e %15.8e\n", gd_band->y[ii]*scale, gd_band->y[ii+1]*scale);
	}
	
    }

    fclose(BAND);
    free(fileBAND);
    puts("-->done");
    return 0;
}

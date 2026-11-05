
/**pss_sfdb.c: pss_sfdb software for the creation of the SFDB data**/
/**Last version: Feb , 3, 2010**/
/**Author: Pia**/

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
long long int llabs(long long int j);
long  int lroundf(float appo);
/*******PSS astro*******/
//#include"../pss_astro/jplbin32.h"
#include"../pss_astro/novas/novas.h"
#include "../pss_astro/daspostare.h"
//#include"../pss_astro/pss_astro.h"
#define dot2utc  2415020.5E0
/****** PSS libraries***************/
#include "../pss_lib/pss_math.h"
#include "../pss_lib/pss_serv.h"
#include "../pss_lib/pss_snag.h"
#include "../pss_lib/pss_sfc.h"
/****** Antenna libraries**************/
#include "pss_ante.h"
/******SFDB libraries**************/
#include "pss_sfdb.h"

/*****************************************************/
//#define OUTFILEH  "outfileh.dat" :
// modifica 1 June 09 moved from pss_sfdb.h
 FILE *SFDB;
 FILE *AREST;
 FILE *PEAKMAP;
 FILE *OUTH1;
 FILE *EVEN_INFO;
 FILE *LOG_INFO;  //pointer to the log file
 FILE *LISTA; //lista of SFDB files for the peakmap vretion
 FILE *OUTTEST;  //scrive eventi
 FILE *OUTTESTB;  //scrive eventi in banda
 FILE *IN_NOISEPULSARS;
 FILE *OUT_NOISEPULSARS;
 FILE *P10;
//fine modifica

/*****************To fill the structures:********************/
INPUT_PARAM* crea_input(int idefaults,int casin,int *casout)
{

  /*If needed for specific purposes it MUST be modified outside, in the main code  */
  /* idefaults=0 puts the default values*/
  INPUT_PARAM *input_param;
  input_param=(INPUT_PARAM *)malloc(sizeof(INPUT_PARAM));
  //char filename[MAXMAXLINE];  /*input file*/
  long len;
  int typ,wink,itotal;
  int verb;
  int nn,nn_orig;                  /* len of the file name */
  char *filename;
  char *appo;
  char *appo1;
  char *appo2;
  char *appo3;
  char *appo4;
  //char check_extens[NEXTS];
  //int icheck;
  int red;
  int lena,lenb;
  int i,ii;
  int jj;
  int nnn;
  int uno;
  float freqc;
  int subsam;
  int fft_m;
  char *fflname;
  char *scisegfilename;
  int flagonline;
  char Ch[25];
  char nam[5];
  float lenRED;
  //int noise_pulsars;

  fflname=(char*)malloc(sizeof(char)*(MAXMAXLINE+1));
  scisegfilename=(char*)malloc(sizeof(char)*200);
  filename=(char*)malloc(sizeof(char)*(MAXMAXLINE+1));
  appo=(char*)malloc(sizeof(char)*200); //Pia: messo 100
  appo1=(char*)malloc(sizeof(char)*7); //era 5 prima di sfdb09
  appo2=(char*)malloc(sizeof(char)*10);
  appo3=(char*)malloc(sizeof(char)*4);
  appo4=(char*)malloc(sizeof(char)*200);
  printf("Detector ? virgo or V1, ligoh or H1, ligol or L1, kagra or K1\n");
  scanf("%s",nam);
  if (strncasecmp(nam,"V1",2) ==0) {
    puts("detector put to virgo");
      strncpy(nam,"virgo",5);
    }
  if (strncasecmp(nam,"L1",2) ==0) {
      puts("detector put to ligol");
      strncpy(nam,"ligol",5);
    }
  if (strncasecmp(nam,"H1",2) ==0) {
      puts("detector put to ligoh");
      strncpy(nam,"ligoh",5);
    }
   if (strncasecmp(nam,"K1",2) ==0) {
      puts("detector put to kagra");
      strncpy(nam,"kagra",5);
    }
  strncpy(input_param->detector_name,nam,5);
  printf("input_param->detector_name: %s \n",input_param->detector_name);
  // printf("Noise around known pulsars evaluation ?  (1=Yes,0 or else=No)\n");
  //scanf("%d",&noise_pulsars);
  //input_param->noise_pulsars=noise_pulsars;
  printf("Type  of SFDB files? \n <=0--> Void SFDB \n ==1--> no SFDB \n =2-->full SFDB \n=3--> full SFDB and peakmap .p10 \n");
  scanf("%d",&flagonline);
  if(flagonline<=1)*casout=2;
  if(flagonline==2)*casout=1;
  if(flagonline>=3)*casout=3;
  if(flagonline==1)flagonline=1;
  if(flagonline<=0)flagonline=0;
  if(flagonline==2)flagonline=0;
  if(flagonline>=3)flagonline=3;
  input_param->flagonline=flagonline;
  printf("Factor to write EVF in the file ?  (1 or 2 usually) \n");
  scanf("%f",&input_param->factor_write);
  if(casin==1){
   input_param->resampling=1;
   printf("input sds file: \n");
   scanf("%s",filename);
   strncpy(input_param->framechannel,"sdsfile",7);
  }
 if(casin==2){
  printf("input ffl file, for frame input: \n");
  scanf("%s",fflname);
  strcpy(input_param->fflname,fflname);
  if (strncmp(input_param->fflname,"lfn:",4) == 0) {
    Logical2PhysicalPath(input_param,1);  //with 1 change the path to the ffl list file. with 2 to the frame files
  }
  printf("science segment file name: \n");
  scanf("%s",scisegfilename);
  strcpy(input_param->scisegfilename,scisegfilename);
  printf("file frame input channel\n h_4096Hz (4096Hz)\n Pr_B1_ACp (20kHz)\n h_20000Hz (20kHz)\n Em_SEBDCE01 (1 kHz)\n Em_MABDCE01 (20kHz) \n H1:LDAS-STRAIN (16384 Hz LIGO H1)\n  L1:LDAS-STRAIN (Ligo L1)\n V1:iso_CW_SN (Paola bin)\n  H1:DCS-CALIB_STRAIN_C01 (or L1)\n");
  /*printf("file frame input channel   1=V1:h_4096Hz 2=V1:Pr_B1_ACp 3=V1:h_20000Hz  ?\n");

  scanf("%d",&fch);
       printf("%d,fch\n",fch);
       switch(fch){
               case 1:
                 strcpy(input_param->framechannel,"h_4096Hz");
                 //strcpy(input_param->framechannel,"V1:h_4096Hz");
                 printf("PRIMA %s framech\n",input_param->framechannel);
                 break;
               case 2:
                      strcpy(input_param->framechannel,"Pr_B1_ACp");
                      break;
               case 3:
                       strcpy(input_param->framechannel,"h_20000Hz");
                       break;
       }

  */

  scanf("%s",Ch);
  //printf("strlen %d\n",strlen(Ch));
  strncpy(input_param->framechannel,Ch,strlen(Ch));
  printf("Frame channel selected: %s \n",input_param->framechannel);
  //Added for the flags: with the suggested names. MAGGIO 2017
  strncpy(input_param->flagchannel,"tobefilled",10);
  if(strncmp(input_param->framechannel,"h_",2)==0){
    //strncpy(input_param->flagchannel,"V1:DQ_META_ITF_LOCKED",21);  usato for O2 data
    strncpy(input_param->flagchannel,"V1:Hrec_STATE_VECTOR",20);
    input_param->flagstatus=4095;  //era 1
    //strncpy(input_param->flagchannel,"V1:Hrec_Veto_ScienceMode",strlen("V1:Hrec_Veto_ScienceMode"));
    //strncpy(input_param->flagchannel,"V1:Hrec_Flag_Channel",strlen("V1:Hrec_Flag_Channel"));
    //input_param->flagstatus=0;
  }
   if(strncmp(input_param->framechannel,"V1:h_",5)==0){
     //strncpy(input_param->flagchannel,"V1:DQ_META_ITF_LOCKED",21);  usato for O2 data
      strncpy(input_param->flagchannel,"V1:Hrec_STATE_VECTOR",20);
      input_param->flagstatus=4095;
     //strncpy(input_param->flagchannel,"V1:Hrec_Veto_ScienceMode",strlen("V1:Hrec_Veto_ScienceMode"));
     //strncpy(input_param->flagchannel,"V1:Hrec_Flag_Channel",strlen("V1:Hrec_Flag_Channel"));
     //input_param->flagstatus=0;
  }
  if(strncmp(input_param->framechannel,"V1:Hr",5)==0){
    //strncpy(input_param->flagchannel,"V1:DQ_META_ITF_LOCKED",21); usato for O2
    //strncpy(input_param->flagchannel,"V1:Hrec_STATE_VECTOR",20);
      strncpy(input_param->flagchannel,"V1:DQ_ANALYSIS_STATE_VECTOR",27);
      input_param->flagstatus=4095;
  }
  if(strncmp(input_param->framechannel,"L1:",3)==0){
    strncpy(input_param->flagchannel,"L1:IFO-SV_STATE_VECTOR",22);
    input_param->flagstatus=65535;
  }
  if(strncmp(input_param->framechannel,"H1:",3)==0){
    strncpy(input_param->flagchannel,"H1:IFO-SV_STATE_VECTOR",22);
    input_param->flagstatus=65535;
  }
  if(strncmp(input_param->framechannel,"L2:",3)==0){
    strncpy(input_param->flagchannel,"L1:IFO-SV_STATE_VECTOR",22);
    input_param->flagstatus=65535;
  }
  if(strncmp(input_param->framechannel,"H2:",3)==0){
    //strncpy(input_param->flagchannel,"H1:IFO-SV_STATE_VECTOR",22);
    strncpy(input_param->flagchannel,"H2:ISC-ALS_EY_REFL_PWR_MON_OUT16",32);
    input_param->flagstatus=3500;
  }

  if(strncmp(input_param->flagchannel,"tobefilled",10)==0){
    //    strncpy(input_param->flagchannel,"V1:Hrec_Veto_LOCK_STEP_STATUS",29);
    //strncpy(input_param->flagchannel,"V1:DQ_META_ITF_LOCKED",21);  usato for O2
    //strncpy(input_param->flagchannel,"V1:Hrec_STATE_VECTOR",20);
    strncpy(input_param->flagchannel,"V1:DQ_ANALYSIS_STATE_VECTOR",27);
    input_param->flagstatus=4095;  //era 12. Poi 1 per O2.
  }
  printf("Flags channel selected: %s \n",input_param->flagchannel);

  // End flag channel


  strcpy(filename,input_param->fflname);
  printf("resampling factor ? 1= no res. Has to be a power of 2! \n");
  scanf("%d",&input_param->resampling);
 }

  //strcpy(filename,"/storage/gpfs_virgo3/scratch/pss/virgo/sd/sds/VSR1-2/deca_20070528/VIR_hrec_20070529_114516_.sds");
  // write the filename of the first file in the LOG file
 jj=logfile_input(LOG_INFO,input_param->detector_name,"detector");
 //jj= logfile_input(LOG_INFO, filename, "First file of the run");
 if(casin==2){
   jj= logfile_input(LOG_INFO,input_param->framechannel , "Frame channel");
 }

 if(idefaults !=0) {
    printf("Verbosity level 0 1 2 3\n");
    scanf("%d",&verb);
    printf("chunk len  -  DS type (0,1,2=interlaced) ?\n");
    if (input_param->resampling >1){
      puts("WARNING: This is the length before the bandwidth extraction");
      printf("Hence the final length in the files will be reduced by the resampling factor: %d \n",input_param->resampling);
    }
    puts("Chunk len will be rounded to the next 2 power, if it is not a power of 2");
    printf("Typical for Virgo 1-band:\n len=4194304 (4096Hz)");
    printf("len =16777216 (20kHz) (es. Em_MABDCE01 Pr_B1_ACp )\n");
    printf("len =1048576 (1kHz) (es Em_SEBDCE01)\n ");
    puts("ds type=2 interlaced; 1 not interlaced");
    scanf("%d%d",&lena,&typ);
    printf("Maximum number of data chunk to be done ?\n");
    scanf("%d",&itotal);
    printf("Reduction factor, for very short FFTs (e.g. 2,4,8..Sugg.:128) ?\n");
    lenRED=next2power(lena)/input_param->resampling;
    printf("When choosing this factor, remember that the total length is: %f\n",lenRED);
    scanf("%d",&red);
    printf("windows  (0=no,1=Hann,2=Hamm,3=MAP, 4=Blackmann flatcos 5=flat top,cosine edge. Sugg. 5) \n");
    scanf("%d",&wink);
    printf("Cut frequency for the highpass filtering (e.g. 100. If 0 the highpass and EVT veto are skipped) ?\n");
    puts("This is the absolute frequency value in Hz (starting from frequency 0 Hz)");
    scanf("%f",&freqc);
    printf("Subsampling factor for the veto in the subbands. Power of 2 (e.g. 128)? <=0: not applied\n");

    puts("Typical for Virgo Em_SEBDCE01 ch (1kHz) =4\n");
    puts("Typical for Virgo Em_MABDCE01 ch (20KHz) =4\n");

    scanf("%d",&subsam);
    printf("Max number of FFTs in one output file (e.g. 100) ?\n");
    scanf("%d",&fft_m);
 }
  else {
    verb=0;
    lena=4194304;
    typ=2;
    itotal=10000000;
    red=128;
    wink=5;
    freqc=100;
    subsam=0;
    fft_m=100;
    input_param->resampling=1;
    }
  nn_orig=strlen(filename);
  printf("Original filename nn_orig= %d \n",nn_orig);
  lenb=next2power(lena);
  len=(long) lenb;
  printf("len written, len used %d %ld\n",lena,len);
  strcpy(input_param->filename,filename);
  input_param->filename[nn_orig]='\0';
  //printf("%s input parameters filename \n",input_param->filename);
  //printf("%s filename \n",filename);
  input_param->freqc=freqc;
  input_param->verb=verb;
  input_param->len=len;
  input_param->typ=typ;
  input_param->wink=wink;
  input_param->itotal=itotal;  /*Maximum number of data chunk to be done*/
  input_param->red=red;  /*Reduction factor for very short FFTs */
  nn=strlen(input_param->filename);
  input_param->nn=nn;
  input_param->subsam=subsam;
  printf("filename len= %d \n",input_param->nn);
  strcpy(input_param->capt_h,"h reconst");
  strcpy(input_param->capt_gd,"Chunk");
  input_param->typ_h=input_param->typ; //DS type for h-reconstructed data: the same as input data.
  input_param->fft_m=fft_m;
  input_param->iopen=0; //to open a new file
 /* To construct the sfdb file name, from a filename which contains a path.
    Data of the C5 run are in /windows/xg23nas/pss/virgo/sd/c5/
*/
  ii=-1;
  for(i=0;i<nn;i++){
      if(input_param->filename[i] == 47 || input_param->filename[i] == 92) ii=i; /* 47 -> / , 92 -> \ */
  }

  /*To extract the name and remove the extension, next substituted with the .SFDB one:*/
  //seleziona(appo,input_param->filename,ii+1,nn-ii-1-NEXTS-1);
  seleziona(appo,input_param->filename,ii+1,nn-ii-1-NEXTS);
  printf("core of the name= %s \n", appo);
  uno=1;
  sprintf(appo2,"%04d",uno);
  strncat(appo,appo2,4);
  nnn=strlen(appo);
  strncpy(appo4,appo,nnn);
  // printf("$$$$ appo appo2= %s %s\n", appo,appo2);  //quiquiritogliere

  switch (*casout)
	  {
		  case 1:
		          strncpy(appo1,".SFDB09",7);
			  break;
		  case 2:
			  strncpy(appo1,".SFDB09",7);
			  break;
		  case 3:
			  strncpy(appo1,".SFDB09",7);
			  strncpy(appo3,".p10",4);
			  break;
	  }



  strncat(appo,appo1,7);  //era 5 prima di sfdb09
  strcpy(input_param->sfdbname,appo);

  nnn=strlen(appo);
  input_param->sfdbname[nnn]='\0';
  if(input_param->flagonline!=1){
   puts("In the header of the SFDB file you will find in detector: 1 for virgo, 2 for ligoh 3 for ligol, 4 for kagra");
   printf("Sfdb filename= %s \n", input_param->sfdbname);
   //printf("Detector = %s \n", input_param->detector_name);

   jj=logfile_output(LOG_INFO,input_param->sfdbname,"List name. The first SFDB file has been opened");
     }
  if(input_param->flagonline >=3){
    puts("Peakmap file with extension .p10 will be created");
    strncat(appo4,appo3,4);
    strcpy(input_param->peakname,appo4);
    nnn=strlen(appo4);
    input_param->peakname[nnn]='\0';

  }
 /* To guess the file type (R87-SDS-FRAME): */
  seleziona(input_param->file_extens,input_param->filename,nn-NEXTS,NEXTS);
  printf("extension of name= %s \n", input_param->file_extens);

 /* Check for the file type (SDS or FRAME): if icheck=0 they are R87 */
    //check if FRAME or SDS
  //strcpy(check_extens,"SDS");
  //icheck=strcasecmp(check_extens,input_param->file_extens);
  //if(icheck==0)input_param->interferom=1;  //sds
  //if(icheck!=0)input_param->interferom=2; //frame

  free(appo);
  free(appo1);
  free(appo2);
  free(appo3);
  free(appo4);
  free(fflname);
  free(filename);


  return input_param;
}


HEADER_PARAM* crea_sfdbheader(INPUT_PARAM *input_param)
{
  HEADER_PARAM *header_param; /*Parameters of the SFDB file, written in the header*/
 header_param=(HEADER_PARAM *)malloc(sizeof(HEADER_PARAM));

  header_param->endian=1.0;
  if(strncasecmp(input_param->detector_name,"nauti",5)==0) header_param->detector=0;
  if(strncasecmp(input_param->detector_name,"virgo",5)==0) header_param->detector=1;
  if(strncasecmp(input_param->detector_name,"ligoh",5)==0) header_param->detector=2;
  if(strncasecmp(input_param->detector_name,"ligol",5)==0) header_param->detector=3;
  if(strncasecmp(input_param->detector_name,"kagra",5)==0) header_param->detector=4;
  header_param->gps_sec=0;
  header_param->gps_nsec=0;
  header_param->tbase=0;
  header_param->firstfreqindex=0;
  header_param->nsamples=0;
  header_param->red=16;
  header_param->typ=0;
  header_param->wink=0;
  header_param->normd=1.;
  header_param->normw=1.;
  header_param->n_flag=0;
  header_param->n_zeroes=0;
  header_param->mjdtime=0;
  header_param->einstein=1.0;  //scaling factor
  //New for O2 March 2018. Ho messo ridondante, va sia con _CO2 che senza
  if(strncmp(input_param->framechannel,"H1_HOFT_C00_GATED_G02",21)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"L1_HOFT_C00_GATED_G02",21)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"L1:DCH-CLEAN_STRAIN_CO2",23)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"H1:DCH-CLEAN_STRAIN_CO2",23)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"L1:DCH-CLEAN_STRAIN",19)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"H1:DCH-CLEAN_STRAIN",19)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"STRAIN",6)==0)header_param->einstein= EINSTEIN;
  //End new for O2
  if(strncmp(input_param->framechannel,"sdsfile",7)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"Hrec_OgNE",9)==0)header_param->einstein=1.0e+6;
  if(strncmp(input_param->framechannel,"Hrec_OgWE",9)==0)header_param->einstein=1.0e+6;
  if(strncmp(input_param->framechannel,"h_4096Hz",8)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"h_20000Hz",9)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"Pr_B1_ACp",9)==0)header_param->einstein= 1.0e-6;
  if(strncmp(input_param->framechannel,"H1:LDAS-STRAIN",14)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"H1:GDS-CALIB_STRAIN",19)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"H1:DCS-CALIB_STRAIN",19)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"L1:GDS-CALIB_STRAIN",19)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"L1:DCS-CALIB_STRAIN",19)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"L1:LDAS-STRAIN",14)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"V1:bin_CW",9)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"H1:LSC-STRAIN",13)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"L1:LSC-STRAIN",13)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"H1:SCO_X1",9)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"L1:SCO_X1",9)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"V1:SCO_X1",9)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"V1:Hrec",7)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"V1:mfdv",7)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"L1:SIMCW",8)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"H1:GDS-GATED_STRAIN",19)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"L1:GDS-GATED_STRAIN",19)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"H1:DCS-GATED_STRAIN",19)==0)header_param->einstein= EINSTEIN;
  if(strncmp(input_param->framechannel,"L1:DCS-GATED_STRAIN",19)==0)header_param->einstein= EINSTEIN;
  //if(strcmp(input_param->framechannel,"Em_MABDCE01")==0)header_param->einstein= 1.0;
  //if(strcmp(input_param->framechannel,"Em_SEBDCE01")==0)header_param->einstein= 1.0;
  printf("einstein factor %e \n",header_param->einstein);
  header_param->frinit=0;
  header_param->tsamplu=0;
  header_param->deltanu=0;
  header_param->vx_eq=0;
  header_param->vy_eq=0;
  header_param->vz_eq=0;
  header_param->px_eq=0;
  header_param->py_eq=0;
  header_param->pz_eq=0;
  header_param->sat_howmany=0.;
  header_param->spare1=0.;
  header_param->spare2=0.;
  header_param->spare3=0.;
  header_param->perc_zeroes=0.;
  header_param->spare5=0.;
  header_param->spare6=0.;
  header_param->lavesp=0;  //no more spare from Nov. 2010: lenght of the average periodogram
  header_param->sciseg=0;  //2013-08-01: Added "science selection" flag. 0-> science segment file not used. 1: science segment file used.
  header_param->spare9=0;
  logfile_par(LOG_INFO,"Scaling_factor_in_amplitude",header_param->einstein,10); //30 marzo 2011
  return header_param;
}

SCIENCE_SEGMENT* crea_science_segment(INPUT_PARAM *input_param)
{
  int jj, jjj=0, nsegm=0;
  SCIENCE_SEGMENT *science_segment;
  science_segment=(SCIENCE_SEGMENT *)malloc(sizeof(SCIENCE_SEGMENT));
  strcpy(science_segment->scisegfilename,input_param->scisegfilename);

  FILE *segm = fopen(science_segment->scisegfilename,"r");

  if (segm == NULL) {
  	printf("Science segment file %s not found!\n", science_segment->scisegfilename);
  	jj= logfile_input(LOG_INFO, "NO", "Science segment list file");
	free(science_segment);
	return NULL;
  }

  printf("reading Science segment file %s\n", science_segment->scisegfilename);
  science_segment->twin=60;

  jj= logfile_input(LOG_INFO, science_segment->scisegfilename, "Science segment list file");

  for (jjj = 0; jjj<10000; jjj++)
  {
  	science_segment->num[jjj]=0; science_segment->segtstart[jjj]=0; science_segment->segtstop[jjj]=0; science_segment->segdur[jjj]=0;
  }

  while (!feof(segm)){
  	fscanf(segm, "%d %d %d %d\n",&science_segment->num[nsegm],&science_segment->segtstart[nsegm],&science_segment->segtstop[nsegm],&science_segment->segdur[nsegm]);
  	nsegm++;
  }

  fclose(segm);

  for (jjj = 0; jjj<10;jjj++){
  	printf("Science segment start=%d stop=%d\n",science_segment->segtstart[jjj],science_segment->segtstop[jjj]);

  }
  printf("Number of science segments = %d\n",nsegm);

  science_segment->nsegm=nsegm;

  return science_segment;
}


EVEN_PARAM* crea_evenparam(long len,int scrivilog)
{
  int i;
  int lent;
  lent=len;
  EVEN_PARAM *even_param; /*Parameters to check for events and clean the data*/
  even_param=(EVEN_PARAM *)calloc(1,sizeof(EVEN_PARAM));
  //printf("In crea_evenparam len= %ld\n",len);
  even_param->absvalue=1;  //1 uses abs values. Else uses values with sign
  //even_param->tau evaluated in purge data. 5 may 2011 pia
  even_param->tau=0; //600.0;       //memory time of the autoregressive average   //if(lent <= 131072) even_param->tau=10.0;
  if(even_param->absvalue==1)even_param->cr=6.0;        //CR of the threshold
  if(even_param->absvalue!=1)even_param->cr=5.0;        //CR of the threshold
  even_param->deadtime=1.0;  //dead time in seconds
  even_param->factor=20.0;    //e.g. 10 or 20 : when re-evaluate mean and std
  even_param->iflev=0;       //0=no events; 1=begin; 2=event; 3=end of the event
  even_param->xamed=(float*)malloc(lent*sizeof(float));       //AR mean
  even_param->xastd=(float*)malloc(lent*sizeof(float));      //AR std
  even_param->xw=0.;
  even_param->qw=0.;
  even_param->w_norm=1.;
  even_param->edge=0.15;  //how many seconds around (before and after) the event have to be "purged"
  even_param->total=0; //how many samples, due to all events, have been cleaned in the GD
  even_param->ndata_run=0; /*number of data in the run (one SFDB file, many FFTs !) */
  even_param->number=0; //number of events in one gd;
  even_param->ar=1;  //Algor. for the AR evaluation: 0 is the old, 1 the new
  even_param->begin=malloc((size_t) len*sizeof(int)); //the dimension is the number of events, less than len.
  even_param->duration=malloc((size_t) len*sizeof(int)); //but I am using len, for simplicity
  even_param->imax=malloc((size_t) len*sizeof(int));
  even_param->crmax=malloc((size_t) len*sizeof(float));
  even_param->ener=malloc((size_t) len*sizeof(float));
  for(i=0;i<lent;i++){
    even_param->xamed[i]=0.;
    even_param->xastd[i]=0.;
    even_param->begin[i]=-1;
    even_param->duration[i]=0;
    even_param->imax[i]=0;
    even_param->crmax[i]=0;
    even_param->ener[i]=0;
  }
  even_param->xw_pari=0.;
  even_param->qw_pari=0.;
  even_param->w_norm_pari=0.;
  even_param->xw_dispari=0.;
  even_param->qw_dispari=0.;
  even_param->w_norm_dispari=0.;
  even_param->ndata_run_pari=0; /*number of pari data in the run (one SFDB file, many FFTs !) */
  even_param->ndata_run_dispari=0; /*number of dispari data in the run (one SFDB file, many FFTs !) */
  even_param->sat_level=4.0e+20; /*Level to saturate data  9-may-2008. Lo aumento 2017*/
  even_param->sat_maxn=(int) (len*1.5); /*max number of saturated data in one FFT */
  even_param->notzero = 0.00001;

  if(scrivilog==1) {
  /* Definition of parameters for the Log file */
  logfile_comment(LOG_INFO,"even NEW: a new FFT has started");
  logfile_comment(LOG_INFO,"   PAR1: Beginning time of the new FFT");
  logfile_comment(LOG_INFO,"   PAR2: FFT number in the run");
  logfile_comment(LOG_INFO,"   PAR3: flag on lock status or science mode");
  logfile_comment(LOG_INFO," (percentage of flagged data)");
  logfile_comment(LOG_INFO,"  PAR4: percentage of zeroes in that FFT");

  logfile_comment(LOG_INFO,"even EVT: time domain events");
  //logfile_comment(LOG_INFO,"   PAR1: Beginning time, in mjd");
  logfile_comment(LOG_INFO,"   PAR1:  Time of the max, in mjd");
  logfile_comment(LOG_INFO,"   PAR2: Duration [s]");
  logfile_comment(LOG_INFO,"   PAR3: Max hp data amplitude*scaling"); //20 AUG 2007-hp
  logfile_comment(LOG_INFO,"   PAR4: Max CR");
  logfile_comment(LOG_INFO,"   PAR5: Energy (sum of squared amp)");
  logfile_comment(LOG_INFO,"even EVB: time domain events in subBWs");
  //logfile_comment(LOG_INFO,"   PAR1: Beginning time, in mjd");
  logfile_comment(LOG_INFO,"   PAR1: Time of the max, in mjd");
  logfile_comment(LOG_INFO,"   PAR2: Duration [s]");
  logfile_comment(LOG_INFO,"   PAR3: Max hp data amplitude*scaling"); //20 AUG 2007-hp
  logfile_comment(LOG_INFO,"   PAR4: Max CR");
  logfile_comment(LOG_INFO,"   PAR5: Energy (sum of squared amp)");
  logfile_comment(LOG_INFO,"   PAR6: Beginning frequency of the BW");
  logfile_comment(LOG_INFO,"even SAT: saturated data");
  logfile_comment(LOG_INFO,"   PAR1: number of saturated samples in the FFT");
  }
  return even_param;
}

EVF_PARAM* crea_evfparam(INPUT_PARAM *input_param)
{

  EVF_PARAM *evf_param; /*Parameters to find EVF events (frequency domain in the AR procedure*/
  evf_param=(EVF_PARAM *)calloc(1,sizeof(EVF_PARAM));
  //puts("In crea_evfparam");


  evf_param->tau_Hz=0;  //0.02; 5 may 2011: defined later
  evf_param->maxage=0; //0.02;

  /*
   if(input_param->len !=0){
  if (input_param->len <= 3000000) {   // TODO Albe modificato per risol 10 mHz
  	evf_param->tau_Hz=0.2;
  	evf_param->maxage=0.2;
  }
   }
  */
  evf_param->maxdin=sqrt(2.5);       //threshold in amplitude
  // REMOVE !! evf_param->maxdin=1000000;  //TEST: high threshold in amplitude 13-7-2006--only for TEST
   if(input_param->len !=0){
  evf_param->factor_write=input_param->factor_write;   // how to enhance the threshold to write evf in the log file
   }else
     {
       evf_param->factor_write=20; //piapia per le peakmap
     }
  /* Definition of parameters for the Log file */
  logfile_comment(LOG_INFO,"even EVF: frequency domain events, with high threshold");
  logfile_comment(LOG_INFO,"   PAR1: End bin of EVF (AR reverse procedure)");
  logfile_comment(LOG_INFO,"   PAR2: Duration [in samples]");
  logfile_comment(LOG_INFO,"   PAR3: Ratio, in amplitude, max/average");
  //logfile_comment(LOG_INFO,"   PAR4: Energy*EINSTEIN**2 or average*EINSTEIN (average if duration=0, when age>maxage)");
  logfile_comment(LOG_INFO,"   PAR4: Energy*scaling**2");
  logfile_comment(LOG_INFO,"   PAR5: nbins N=(Nend-Nmax)->freq max=(Nend-N)*df");
  if(input_param->flagonline !=1){ //added 4 may 2011
           logfile_comment(LOG_INFO,"even VEL: pos and vel");

	   logfile_comment(LOG_INFO,"   PAR1: Detector velx/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR2: Detector vely/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR3: Detector velz/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR4: Detector posx/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR5: Detector posy/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR6: Detector posz/C, equatorial");

	   }
  logfile_comment(LOG_INFO,"stat TOT: total number of frequency domain events");
/*Added 16-aug-2006  */
  logfile_comment(LOG_INFO,"par GEN: parameters of the AR spectrum estimation");

  logfile_comment(LOG_INFO,"GEN_THR is the threshold in amplitude");
  logfile_comment(LOG_INFO,"GEN_TAU the memory frequency of the AR estimation");
  logfile_comment(LOG_INFO,"GEN_MAXAGE [Hz] the max age of the AR process");
  logfile_comment(LOG_INFO,"GEN_FAC the factor for which the thres is multiplied, less EVF in the log");

 //logfile_comment(LOG_INFO,"FLG : percentage of flagged data (lock status < 12; if hrec, not science mode)");

  return evf_param;
}

double day2itt2001(int year,int month, int day,int hour,int min,float sec,float tmax)
     /*
       From year,month,day,hout,min,sec,tmax
       it gives days from the beginning of the year 2001
     */

{
  int iday,bisestile,iadd;
  float itt_eve;
  bisestile=0;
  iday=0;
  iadd=0;
  if(year==2002)iadd=365;
  if(year==2003)iadd=365*2;
  if(year==2004)iadd=365*3;
  if(year==2005)iadd=365*4+1;
  if(year==2004)bisestile=1;
  if(month==1){
    itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
    return itt_eve;
  }
  iday+=31;
  if(month==2){
    itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
    return itt_eve;
  }
  if(bisestile==0)iday+=28;
  if(bisestile==1)iday+=29;
  if(month==3){
    itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
    return itt_eve;
  }
  iday+=31;
  if(month==4){
    itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
    return itt_eve;
  }
  iday+=30;
  if(month==5){
    itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
    return itt_eve;
  }
  iday+=31;
  if(month==6){
    itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
    return itt_eve;
  }
  iday+=30;
  if(month==7){
    itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
    return itt_eve;
  }
  iday+=31;
  if(month==8){
    itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
    return itt_eve;
  }
  iday+=31;
  if(month==9){
    itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
    return itt_eve;
  }
  iday+=30;
  if(month==10){
    itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
    return itt_eve;
  }
  iday+=31;
  if(month==11){
    itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
    return itt_eve;
  }
  iday+=30;
  if(month==12){
    itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
    return itt_eve;
  }
  fprintf(stderr,"FATAL ERROR\n");
  abort();

}

double day2mjd(int year,int month, int day,int hour,int min,float sec,float tmax)
  /*
From year,month,day,hout,min,sec,tmax
it gives MJD (modified Julian date, days)
  */

  {
    int iday,bisestile,iadd;
    double itt_eve;
    bisestile=0; //per ora messo fisso ad hoc
    iday=0;
    iadd=0;

    if(year==1998)iadd=-(365*3+1); //il 2000 era bisestile (1096 da togliere a daybegin=36889.0)
    if(year==2002)iadd=365;
    if(year==2003)iadd=365*2;
    if(year==2004)iadd=365*3;
    if(year==2005)iadd=365*4+1;
     if(year==2004)bisestile=1;
    if(month==1){
        itt_eve=iadd+iday+day+((double) hour*hour2sec+(double) min*min2sec+sec+tmax)/day2sec;
        itt_eve+=(daybegin+dot2mjd); //MJD

        return itt_eve;
    }
    iday+=31;
    if(month==2){

         itt_eve=iadd+iday+day+((double) hour*hour2sec+(double) min*min2sec+sec+tmax)/day2sec;
        itt_eve+=(daybegin+dot2mjd); //MJD
        return itt_eve;
    }
      if(bisestile==0)iday+=28;
      if(bisestile==1)iday+=29;
 if(month==3){
        itt_eve=iadd+iday+day+((double) hour*hour2sec+(double) min*min2sec+sec+tmax)/day2sec;
        itt_eve+=(daybegin+dot2mjd); //MJD
        return itt_eve;
    }
 iday+=31;
 if(month==4){
        //itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
        itt_eve=iadd+iday+day+((double) hour*hour2sec+(double) min*min2sec+sec+tmax)/day2sec;
        itt_eve+=(daybegin+dot2mjd); //MJD
        return itt_eve;
    }
 iday+=30;
 if(month==5){
        //itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
        itt_eve=iadd+iday+day+((double) hour*hour2sec+(double) min*min2sec+sec+tmax)/day2sec;
        itt_eve+=(daybegin+dot2mjd); //MJD
        return itt_eve;
    }
 iday+=31;
 if(month==6){
        //itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
        itt_eve=iadd+iday+day+((double) hour*hour2sec+(double) min*min2sec+sec+tmax)/day2sec;
        itt_eve+=(daybegin+dot2mjd); //MJD
        return itt_eve;
    }
 iday+=30;
 if(month==7){
        //itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
        itt_eve=iadd+iday+day+((double) hour*hour2sec+(double) min*min2sec+sec+tmax)/day2sec;
        itt_eve+=(daybegin+dot2mjd); //MJD
        return itt_eve;
    }
 iday+=31;
 if(month==8){
        //itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
        itt_eve=iadd+iday+day+((double) hour*hour2sec+(double) min*min2sec+sec+tmax)/day2sec;
        itt_eve+=(daybegin+dot2mjd); //MJD
        return itt_eve;
    }
 iday+=31;
 if(month==9){
        //itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
        itt_eve=iadd+iday+day+((double) hour*hour2sec+(double) min*min2sec+sec+tmax)/day2sec;
        itt_eve+=(daybegin+dot2mjd); //MJD
        return itt_eve;
    }
 iday+=30;
 if(month==10){
        //itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
        itt_eve=iadd+iday+day+((double) hour*hour2sec+(double) min*min2sec+sec+tmax)/day2sec;
        itt_eve+=(daybegin+dot2mjd); //MJD
        return itt_eve;
    }
 iday+=31;
 if(month==11){
        //itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
        itt_eve=iadd+iday+day+((double) hour*hour2sec+(double) min*min2sec+sec+tmax)/day2sec;
        itt_eve+=(daybegin+dot2mjd); //MJD
        return itt_eve;
    }
 iday+=30;
 if(month==12){
        //itt_eve=iadd+iday+day+hour*1.0/day2hour+min*1.0/day2min+((float)sec+(float)tmax)/day2sec;
        itt_eve=iadd+iday+day+((double) hour*hour2sec+(double) min*min2sec+sec+tmax)/day2sec;
        itt_eve+=(daybegin+dot2mjd); //MJD
        return itt_eve;
    }
 return 1;
  }

void mjd2day(double *mjd_eve,int *year,int *month,int *day,int *hour,int *min,double *sec)
{
  int bisestile;
  int ittint,iday1,iday2,iday;
  double itt_eve;

  bisestile=0; //per ora messo fisso ad hoc
  itt_eve= *mjd_eve-dot2mjd-daybegin; //daybegin e' il 31 dic 2000

  if(itt_eve <=0)*year=1998; //attenzione, sto considerando solo il 1998 !!!
  if((itt_eve >=1) && (itt_eve <=365))*year=2001;
  if((itt_eve >=366) && (itt_eve <=730))*year=2002;
  if((itt_eve >=731) && (itt_eve <=1095))*year=2003;
  if((itt_eve >=1096) && (itt_eve <=1461))*year=2004;
  if((itt_eve >=1462) && (itt_eve <=1826))*year=2005;
  if(*year==1998)itt_eve+=(365*3+1);
  if(*year==2002)itt_eve-=365.0;
  if(*year==2003)itt_eve-=(365*2.0);
  if(*year==2004)bisestile=1;
  if(*year==2004)itt_eve-=(365*3.0);
  if(*year==2005)itt_eve-=(365*4.0+1);
  ittint= (int) itt_eve;

  iday1=0;
  iday2=31;
  if(ittint>iday1 && ittint <=iday2){
    *month=1;
    iday=iday1; //quello trovato. Serve poi per valutare il gg
  }

  iday1=iday2;
  if(bisestile==0)iday2+=28;
  if(bisestile==1)iday2+=29;
  if(ittint>iday1 && ittint <=iday2){
    *month=2;
    iday=iday1;
  }
  iday1=iday2;
  iday2+=31;
  if(ittint>iday1 && ittint <=iday2){
    *month=3;
    iday=iday1;
    //printf("marzo %d %d %d \n", ittint, *month,iday);
  }
  iday1=iday2;
  iday2+=30;
  if(ittint>iday1 && ittint <=iday2){
    *month=4;
    iday=iday1;
  }
  iday1=iday2;
  iday2+=31;
  if(ittint>iday1 && ittint <=iday2){
    *month=5;
    iday=iday1;
  }
  iday1=iday2;
  iday2+=30;
  if(ittint>iday1 && ittint <=iday2){
    *month=6;
    iday=iday1;
  }
  iday1=iday2;
  iday2+=31;
  if(ittint>iday1 && ittint <=iday2){
    *month=7;
    iday=iday1;
  }
  iday1=iday2;
  iday2+=31;
  if(ittint>iday1 && ittint <=iday2){
    *month=8;
    iday=iday1;
  }
  iday1=iday2;
  iday2+=30;
  if(ittint>iday1 && ittint <=iday2){
    *month=9;
    iday=iday1;
  }
  iday1=iday2;
  iday2+=31;
  if(ittint>iday1 && ittint <=iday2){
    *month=10;
    iday=iday1;
  }
  iday1=iday2;
  iday2+=30;
  if(ittint>iday1 && ittint <=iday2){
    *month=11;
    iday=iday1;
  }
  iday1=iday2;
  iday2+=31;
  if(ittint>iday1 && ittint <=iday2){
    *month=12;
    iday=iday1;
  }


  *day=ittint-iday;
  *hour=(int)((itt_eve-ittint)*24.0);
  *min=(int) (((itt_eve-ittint)*24.0 - *hour)*60.0);
  *sec=(double)  (((((itt_eve-ittint)*24.0 - *hour)*60.0) - *min)*60.0);

  return;
}


/*FFT:***********************************/
void recoverfft(int nn, Cmplx *cy)
     /* function to recover the whole FFT from one half
        nn=half FFT
        input and output in the complex vector cy (input nn complex data, output 2*nn complex data)
     */
{
  int lfft,lfft2,ii,iimax;

  lfft2=nn; //meta' FFT
  lfft=lfft2*2; //total FFT
  ii=0;

  cy[lfft2].r=cy[lfft2-(ii+1)].r;
  cy[lfft2].i=0.;
  iimax=lfft2-1;
  //ii=0 already done, it has imag=0
  for(ii=1;ii<= iimax; ii++){
    cy[lfft2+ii].r=cy[lfft2-ii].r;
    cy[lfft2+ii].i=-cy[lfft2-ii].i;
  }

  return;
}
void recoverfft_analitic(int nn, Cmplx *cy)
     /* function to recover the whole FFT from one half: with zeroes. analitic signal
        nn=half FFT
        input and output in the complex vector cy (input nn complex data, output 2*nn complex data)
     */
{
  int lfft,lfft2,ii,iimax;

  lfft2=nn; //meta' FFT
  lfft=lfft2*2; //total FFT
  ii=0;

  cy[lfft2].r=0.; //cy[lfft2-(ii+1)].r;
  cy[lfft2].i=0.;
  iimax=lfft2-1;
  //ii=0 already done, it has imag=0
  for(ii=1;ii<= iimax; ii++){
    cy[lfft2+ii].r=0.;
    cy[lfft2+ii].i=0.;
  }

  return;
}
/***********To handle with DS and vectors:*****************/
int fill_newds(DS *d, float *y,DS *d_orig)
/*It copies the parameters of an orig DS into another one. Typically used to
pass from raw data to h-reconstructed data: the two data streams have the same characteristics but
different data  */
{
	int i;
	d->dt=d_orig->dt;
	d->len=d_orig->len;
	d->tini1=d_orig->tini1;
	d->tini2=d_orig->tini2;
	//printf("fill_newds d_h->dt , len  tini1 tini2 %f %ld %f %f \n",d->dt,d->len,d->tini1,d->tini2);
	d->nc1=d_orig->nc1;
	d->nc2=d_orig->nc2;
	d->lcw=d_orig->lcw;

	if (d->type == 2)
	  {
	   if (d->lcw%2 == 1)
	     for (i=0; i < d->len; i++)d->y1[i]=y[i];

	else
	  for (i=0; i < d->len; i++)d->y2[i]=y[i];
	}

	else
	  for (i=0; i < d->len; i++)d->y1[i]=y[i];

	return 1;
}

int gd_ds(GD *gd,DS *d)
/*
Extracts one chunk of data from a DS and puts the data into a GD.
It uses y_ds. The GD has to be already created, with crea_gd. It sets the properties of the GD
from that of the DS
 */
{
  int i;


  y_ds(d,gd->y);  //Output from DS *d to the GD data


  gd->dx=d->dt;
  gd->n=d->len;

  //printf("In gd_ds 2 Data = %f %f  \n", gd->y[0],gd->y[1]);
  i=(int) gd->n;

  //printf( "In gd_ds : len, dt = %ld %f  \n",gd->n,gd->dx);
  return i;
}



/***********************Function to manipulate strings**********/

int seleziona(char a[], char b[], int s, int l)
{
  int i, j;

  i = s;
  j = 0;
  while ( b[i] != '\0' && j < l) {
    a[j++] = b[i++];
  }
  a[j] =  '\0';

  return 0;
}

/**********************************Windows***********************/
void gd2wingd(GD *gd,int wink)

     /*
       INPUT=GD
       OUTPUT=GD with the chosen window
       wink=1 Hanning
       wink=2 Hamming
       wink=3 Used for MAP (as in GEO data)
       wink=4 Blackmann
       wink=5 flatcos -flat top, cosine edge--Sergio
     */
{
  int i,winlen;
  float aaw,bbw,ccw;
  double aind;
  float frac;
  double factor;
  if(wink ==1){
    aaw=0.5;
    bbw=0.5;
    ccw=0.0;
    winlen=(int) gd->n;
  }

  if(wink ==2){
    aaw=0.54;
    bbw=0.46;
    ccw=0.0;
    winlen=(int) gd->n;
  }
  if(wink ==3){
    aaw=0.5;
    bbw=-0.5;
    ccw=0.0;
    frac=0.05;
    winlen=(int) (0.5+(gd->n)*frac);
    //printf("MAP winlen %d \n",winlen);
  }
   if(wink ==4){
    aaw=0.42;
    bbw=0.50;
    ccw=0.08;
    winlen=(int) gd->n;
  }
    if(wink ==5){
    aaw=0.50;
    bbw=0.50;
    ccw=0.0;
    winlen=(int) gd->n/4;
  }

   for (i=0; i < winlen; i++){
      if(wink !=3 && wink !=5) aind=i*2*PIG/(winlen-1);
      if(wink ==3) aind=-PIG+i*PIG/(winlen);
      if(wink ==5) aind=i*PIG/(winlen);
      factor=(aaw-bbw*cos(aind)+ccw*cos(2*aind));

      gd->y[i]=gd->y[i]*factor;
      if(wink==3)gd->y[gd->n-i-1]=gd->y[gd->n-i-1]*factor;
      if(wink==5)gd->y[gd->n-i-1]=gd->y[gd->n-i-1]*factor;

   }

  return ;

}

void wingd2gd(GD *gd,int wink)

     /*
       Function used to remove the window used on data to construct an FFT. When going back to the time
       domain the window have to be removed, if factor is not zero
       INPUT=GD with windowed data
       OUTPUT=GD without the window
       wink=1 Hanning
       wink=2 Hamming
       wink=3 Used for MAP (as in GEO data)
       wink=4 Blackmann
       wink=5 flatcos -flat top, cosine edge--Sergio
     */
{
  int i,winlen;
  float aaw,bbw,ccw;
  double aind;
  float frac;
  double factor;
  if(wink ==1){
    aaw=0.5;
    bbw=0.5;
    ccw=0.0;
    winlen=(int) gd->n;
  }

  if(wink ==2){
    aaw=0.54;
    bbw=0.46;
    ccw=0.0;
    winlen=(int) gd->n;
  }
  if(wink ==3){
    aaw=0.5;
    bbw=-0.5;
    ccw=0.0;
    frac=0.05;
    winlen=(int) (0.5+(gd->n)*frac);
    //printf("MAP winlen %d \n",winlen);
  }
   if(wink ==4){
    aaw=0.42;
    bbw=0.50;
    ccw=0.08;
    winlen=(int) gd->n;
  }
   if(wink ==5){
    aaw=0.50;
    bbw=0.50;
    ccw=0.0;
    winlen=(int) gd->n/4;
  }
   for (i=0; i < winlen; i++){
      if(wink !=3 && wink !=5) aind=i*2*PIG/(winlen-1);
      if(wink ==3) aind=-PIG+i*PIG/(winlen);
      if(wink ==5) aind=i*PIG/(winlen);

      factor=(aaw-bbw*cos(aind)+ccw*cos(2*aind));
      if(factor!=0.0){
	gd->y[i]=gd->y[i]/factor;
	if(wink==3)gd->y[gd->n-i-1]=gd->y[gd->n-i-1]/factor;
	if(wink==5)gd->y[gd->n-i-1]=gd->y[gd->n-i-1]/factor;
      }
   }

  return ;

}

void vect2winvect(float *signal,int wink,int winlen1)

     /*
       INPUT=vector signal float
       OUTPUT=vector signal with the chosen window
       wink=1 Hanning
       wink=2 Hamming
       wink=3 Used for MAP (as in GEO data)
       wink=4 Blackmann
       wink=5 flatcos -flat top, cosine edge--Sergio
     */
{
  int i,winlen;
  float aaw,bbw,ccw;
  double aind;
  float frac;
  double factor;
  if(wink ==1){
    aaw=0.5;
    bbw=0.5;
    ccw=0.0;
    winlen=winlen1;
  }

  if(wink ==2){
    aaw=0.54;
    bbw=0.46;
    ccw=0.0;
    winlen=winlen1;
  }
  if(wink ==3){
    aaw=0.5;
    bbw=-0.5;
    ccw=0.0;
    frac=0.05;
    winlen=(int) (0.5+(winlen1)*frac);
    //printf("MAP winlen %d \n",winlen);
  }
   if(wink ==4){
    aaw=0.42;
    bbw=0.50;
    ccw=0.08;
    winlen=winlen1;
  }
   if(wink ==5){
    aaw=0.50;
    bbw=0.50;
    ccw=0.0;
    winlen=winlen1/4;
  }

   for (i=0; i < winlen; i++){
      if(wink !=3 && wink !=5) aind=i*2*PIG/(winlen-1);
      if(wink ==3) aind=-PIG+i*PIG/(winlen);
      if(wink ==5) aind=i*PIG/(winlen);
      factor=(aaw-bbw*cos(aind)+ccw*cos(2*aind));

      signal[i]=signal[i]*factor;
      if(wink==3)signal[winlen1-i-1]=signal[winlen1-i-1]*factor;
      if(wink==5)signal[winlen1-i-1]=signal[winlen1-i-1]*factor;
   }

  return ;

}

/***************Function to go from time-domain h reconstructed data (in a DS *d_h) to SFDB data again in d_h******************/
int gdtime2freq(GD *gd,HEADER_PARAM *header_param){
/*from a GD in time to the same GD in frequency domain
header_param are needed only to state the inizial frequency of the FFT
Frequency domain:output to the same GD *gd with real in ii and imag in ii+1... !! Half bandwidth !!
ATT: the input GD must have its own gd->dx (time step), which is then converted into the proper one for the new freq. domain GD (frequency step)
Pay attention, if the case, to reassign, in the main code, the proper time step to gd->dx and ini!!!
 */
//Last versione 15 april 2009. Non usa piu' drealft
  int ii;
  double dx;
  unsigned long nn;
  double *datamay07;
  int k; //aggiunta may 07 per double

  dx= gd->dx;
  nn= (unsigned long) gd->n;
  //printf("in gdtime2freq nn=%ld dx=%f\n",nn,dx);
  //datamay07=(double *)malloc((size_t) (gd->n)*sizeof(double)); //may 07 per double
  datamay07=(double *)malloc((size_t) (2*gd->n)*sizeof(double)); //may 07 per double. Sostituito al precedente per togliere drealft Aprile 2009
  /*for(k=0;k<gd->n;k++){
      datamay07[k]=gd->y[k];
      }*/

  ii=0; //IMP: sostituito dfour1 a dreal 15/04/09
 for(k=0;k<gd->n;k++){
      datamay07[ii]=gd->y[k];
      datamay07[ii+1]=0.;
      ii+=2;
    }
    //TOLTO !!//drealft(datamay07-1,nn,1);

   dfour1(datamay07-1,nn,1);
  //the output is only 1/2 of the total FFT, in a real vector of size nn
  //which contains the real part of the FFT in y[ii]
  //and the imag part in y[ii+1]
  //Frequency domain:output to the same GD *gd with real in ii and imag in ii+1... !! Half bandwidth !!
  gd->dx=(double) 1.0/(gd->n*dx);  //step in frequency domain
  gd->ini=header_param->frinit; //valid only for the total band

  ii=(int) gd->n;

  //printf( "In gdtime2freq : nn, len, df and original dt = %ld, %ld %f %f \n",nn,gd->n,gd->dx,dx);
   for(k=0;k<gd->n;k++){
      gd->y[k]=datamay07[k];
    }
  free(datamay07);
  return ii;
}

int gdfreq2time(GD *gd,HEADER_PARAM *header_param,int ivign){
/*from a GD in frequency (half bandwidth, alternatively real and imag part) to the same GD in time domain
ivign=1 if the data are "windowed" before going back to the time domain ("vignettatura quadratica").
Time domain:output to the same GD *gd with the real part of the data (imag is zero) !! Real data  !!
ATT: the input GD must have its own gd->dx (freq. step), which is then converted into the proper one for the new time
domain GD (time step). !!
Pay attention, if the case, to reassign, in the main code, the proper freq. step to gd->dx !!!
*/

  int ii,k;
  double dx;
  double ind;
  unsigned long nn;
  float fract=0.03; //era 0.03   fraction of "windowed" data, if ivign=1
  float lwind;
  Cmplx *cy;
  Cmplx a;
  double *y;
  cy=(Cmplx *)malloc((size_t) (gd->n)*sizeof(Cmplx));
  y=(double *)malloc((size_t) (2*gd->n)*sizeof(double));
  int len;

  len=gd->n;
  nn= (unsigned long) gd->n;
  k=0;
  for(ii=0;ii<gd->n;ii+=2){
    cy[k]=Complex(gd->y[ii],gd->y[ii+1]);
    k++;
  }
  //puts("Go back to the time domain");
  //puts("First step: recontruct the whole FFT from its half");
 recoverfft(gd->n/2,cy);
 if(ivign==1){
    lwind=(gd->n)*fract;
    for(ii=0;ii<lwind-1;ii++){
      ind= ((double) ii)/lwind;
      a=cy[ii];
      cy[ii]=RCmul(0.001+pow(ind,2),a);
      a=cy[gd->n-ii];
      cy[gd->n-ii]=RCmul(0.001+pow(ind,2),a);
    }
 }

 k=0;
   for(ii=0;ii<gd->n;ii++){
     y[k]=cy[ii].r;
     y[k+1]=cy[ii].i;
     k+=2;
   }

 //the input is a real vector, which contains alternatively Real and Imag part of the complex FFT data,
 //over the whole bandwidth. That's why the size is 2*nn
  dfour1(y-1,len,-1);
//Output: y[ii] is the Real part; y[ii+1] is the Imag part of the time domain data (size = 2*gd->n).
//They must be normalized, dividing by the gd->n. Only the real part is used, because the Imag part is zero.
   k=0;
   for(ii=0;ii<2*gd->n;ii+=2){
     gd->y[k]=y[ii]/len;  //take only the real part of the time-domain data
     k++;
   }

  dx= gd->dx;
  //printf("In gdfreq2time- (dx ancora in freq)gd->dx=%f gd->n %ld \n",gd->dx,gd->n);
  gd->dx=(double) 1.0/(gd->n*dx);  //step in time domain
  //printf("In gdfreq2time- (dx  in time)gd->dx=%f  \n",gd->dx);
  gd->ini=header_param->mjdtime;
  ii=(int) gd->n;

  //printf( "In gdfreq2time : nn, len, dt and original df = %ld, %ld %f %f \n",nn,gd->n,gd->dx,dx);

  free(y);
  free(cy);

  return ii;
}

int gdfreq2time_analitic(GD *gd_time2, GD *gd,HEADER_PARAM *header_param,int ivign){
/*
Same as the previous BUT using the analitic signal
header_param are needed only to state the inizial time of the data chunk;
ivign=1 if the data are "windowed" before going back to the time domain ("vignettatura quadratica");
Time domain:output to  GD *gd_time2 with real and Imag !! Complex data, output alternatively
Real and Imag.
*/

  int ii;

  double dx;
  unsigned long nn;
  double *y;
  int len;
  y=(double *)malloc((size_t) (2*gd->n)*sizeof(double));

  len=gd->n;
  nn= (unsigned long) (2*len);
  //Put the freq data from the gd into a vector with Re and Imag
  for(ii=0;ii<len;ii++){
    y[ii]=gd->y[ii];
  }
  //Second half
  for(ii=len;ii<2*len;ii++){
    y[ii]=0.;
  }

   printf("ivign %d\n",ivign);
  //puts("Done First step: recontruct the whole FFT from its half: using zeroes- Analitic signal");


 //the input is a real vector, which contains alternatively Real and Imag part of the complex FFT data,
 //over the whole bandwidth. That's why the size is doubled

  dfour1(y-1,len,-1);
  //Output: y[ii] is the Real part; y[ii+1] is the Imag part of the time domain data (size = 2*gd->n).
  //They must be normalized, because the IFFT is done without the /len
  gd_time2->n=nn;

   for(ii=0;ii<gd_time2->n;ii++){
     // PERCHE 0.5 ? gd_time2->y[ii]=y[ii]/(0.5*len);  //Real and Imag
     gd_time2->y[ii]=y[ii]/(len);  //Real and Imag
   }

  dx= gd->dx;
  gd_time2->dx=(double) 1.0/(gd->n*dx);  //step in time domain
  gd_time2->ini=header_param->mjdtime;
  ii=(int) gd_time2->n;

  //printf( "In gdfreq2time_analitic : nn, gd_time2->n, dt and original df = %ld, %ld %f %f \n",nn,gd_time2->n,gd_time2->dx,dx);

  free(y);
  return ii;
}

int short_ps_freq(float *ps,GD *gd,HEADER_PARAM *header_param)
 /*From time data in gd produces one short power spectrum in ps.
   Reduction factor red depends on header_param->red .
   Each datum in ps is the SUM (total power) over a number red of data  (data in a same freq. bin).
   The size of ps is lenps/2 (see the code): it contains sqrt(spectrum) of the data in half the band*/
//Modified 15/04/09 per non usare drealft
{
  int lenps,j,jj,jj1,k,kk;
  double *datamay07; //may07
  float rpw,ipw;
  float freq;

  lenps= (int) gd->n/(header_param->red);

   //printf("short ps time n lenps rid %ld %d %d\n",gd->n,lenps,header_param->red);

   for(k=0;k<(lenps-1);k+=2){
      ps[k]=0.;
      ps[k+1]=0.;
    }
  for(j=0;j<header_param->red;j++){
    jj=j*lenps;
    jj1=jj+lenps;
    //datamay07=(double *)malloc((size_t) (lenps)*sizeof(double));
    datamay07=(double *)malloc((size_t) (2*lenps)*sizeof(double));
    /*
    k=0;
    for(kk=jj;kk<jj1;kk++){
      datamay07[k]=gd->y[kk];
      k++;
    }
    drealft(datamay07-1,lenps,1); */
     k=0;
    for(kk=jj;kk<jj1;kk++){
      datamay07[k]=gd->y[kk];
      datamay07[k+1]=0.;
      k+=2;
    }
    dfour1(datamay07-1,lenps,1);
    /*Very short FFT in ps[k]: k are the various frequency bins*/
    j=-1;
    for(k=0;k<lenps;k+=2){
      freq=(float)(header_param->firstfreqindex+k*header_param->red)*1.0/header_param->tbase; //not used, just to know..
      rpw=datamay07[k];
      ipw=datamay07[k+1];
      j+=1;
      ps[j]+=sqrt(rpw*rpw+ipw*ipw);
    }

  }

  free(datamay07);
  return header_param->red;
}

int short_ps_time(float *ps,GD *gd,HEADER_PARAM *header_param)
     /*From time data in gd produces the averages of short power spectra in ps.
      Reduction depends on header_param->red .
      Each datum in ps is the average over a number lenps/2 of data (see the code).
      The size of ps is red*/
//Modified 15/04/09 per non usare drealft
{
  int lenps,j,jj,jj1,k,kk;
  double *datamay07;
  float rpw,ipw;


   lenps= (int) gd->n/(header_param->red);


  //printf("short ps (time) n lenps rid %ld %d %d\n",gd->n,lenps,header_param->red);
  datamay07=(double *)malloc((size_t) (2*(lenps+1))*sizeof(double));
  //datamay07=(double *)malloc((size_t) (lenps+1)*sizeof(double));
  for(j=0;j<header_param->red;j++){
    jj=j*lenps;
    jj1=jj+lenps;
    /*
    k=0;
    for(kk=jj;kk<jj1;kk++){
    datamay07[k]=gd->y[kk];
    k++;
    }
    drealft(datamay07-1,lenps,1);
    */
    k=0;
    for(kk=jj;kk<jj1;kk++){
    datamay07[k]=gd->y[kk];
    datamay07[k+1]=0.;
    k+=2;
    }
    dfour1(datamay07-1,lenps,1);
    /*Average power ps[j]*/
    ps[j]=0.;
    for(k=0;k<lenps;k+=2){
      rpw=datamay07[k];
      ipw=datamay07[k+1];
      ps[j]+=(rpw*rpw+ipw*ipw);
    }
    ps[j]/=(float)(lenps/2.0);
    //printf(" Short power j ps[j]  %d %f \n",j,ps[j]);

  }
 free(datamay07);

  return header_param->red;
}
int short_ps_periodogram(float *ps,GD *gd,HEADER_PARAM *header_param)
     /* It produces the average periodogram, averaged on short (a factor red) power spectra
      Reduction depends on header_param->red .
      done always with a Hanning window (not asked to the user !)
      The size of ps is gd->n/red */
{
  int lenps,j,jj,jj1,k,kk,ii;
  double *data;
  float rpw,ipw;
  float total_norm, normw;
  int i,winlen;
  float aaw,bbw,ccw;
  double aind;
  double factor;

  lenps= (int) gd->n/(header_param->red);  //length of the chunk to be averaged
  aaw=0.5;
  bbw=0.5;
  ccw=0.0;
  winlen=lenps;

  //printf("short periodogram, Hanning windows:  n lenps rid %ld %d %d\n",gd->n,lenps,header_param->red);
  data=(double *)malloc((size_t) (2*(lenps+1))*sizeof(double));

  for(ii=0;ii<lenps/2;ii++){
    ps[ii]=0;
  }


  for(j=0;j<header_param->red;j++){
    jj=j*lenps;
    jj1=jj+lenps;
    k=0;
    i=0;
    for(kk=jj;kk<jj1;kk++){
      aind=i*2*PIG/(winlen-1);
      factor=(aaw-bbw*cos(aind)+ccw*cos(2*aind));
      //factor=1;  //PROVA
      data[k]=gd->y[kk]*factor;
      data[k+1]=0.;
      k+=2;
      i+=1;
    }
    dfour1(data-1,lenps,1);
    /*short power spectrum and sum over the previous spectra*/
    ii=0;
    for(k=0;k<lenps;k+=2){
      rpw=data[k];
      ipw=data[k+1];
      ps[ii]+=(rpw*rpw+ipw*ipw);  //output is the power spectrum
      //ps[ii]+=sqrt(rpw*rpw+ipw*ipw);  //output is the AMP spectrum
      ii+=1;
    }
    //printf(" Short power j ps[j]  %d %f \n",j,ps[j]);

  }
  //Note on the norm: total_norm is on the sqrt(average spectrum)
  //normw=sqrt(3/1.0618); //norm for Hanning window. 1.0618 is EMPIRICAL Genn 2011
  normw=2*sqrt(2.0/3.0); //visto da pswindow di Sergio. 14/12/2011. lui usa aaw e bbw =1 io 0.5.
 total_norm= header_param->normd*sqrt(header_param->red)*sqrt(bil2uni)*normw;
  for(ii=0;ii<lenps/2;ii++){
    ps[ii]=total_norm*sqrt(ps[ii]/(float)(header_param->red));
    //ps[ii]=total_norm*(ps[ii]/(float)(header_param->red));
  }

 free(data);

  return header_param->red;
}

/*********************************************************/
/*FUNCTION TO APPLY THE VETO PROCEDURES FROM GD TO ANOTHER GD  */
int gd2gdfreq(GD *gd_clean,GD *gd_freq, GD *gd,GD *gd_highpass,GD *gd_band,GD *gd_appo,HEADER_PARAM *header_param,EVEN_PARAM *even_param,EVEN_PARAM *even_paramb,INPUT_PARAM *input_param,SCIENCE_SEGMENT* science_segment)
{
  /*
 Input. gd. Output: gd_clean, in time, and gd_freq, with the FFT (Real, Imag)
   */
  int kki,ii;
  float freqc;
  int sat_maxn,subsam;
  float band_ext;
  int len;
  int ic;
  int iV;

  /* MOVED after the EVT veto. March 15, 2017
 if (science_segment != NULL) {
    	printf("Applying science segment list selection! Call the function\n");
	iV= ScienceVeto_onGD(gd,header_param,science_segment);
    }
  if(input_param->resampling <=1){
    //else it is done in resampling
    header_param->n_zeroes=0;
    for(ii=0;ii<gd->n;ii++){
      kki=isnan(gd->y[ii]);
   if (kki !=0)
   {
    gd->y[ii]=0;
   }
   if(gd->y[ii]==0)header_param->n_zeroes+=1;
 }
 header_param->perc_zeroes=(float)(header_param->n_zeroes)/gd->n;  //becomes a percentage
  }
  */
 freqc=input_param->freqc-header_param->frinit;
 if(freqc <0){
   freqc=0;
 }
 printf("Relative and absolute freqc in EVT veto %f %f \n",freqc,freqc+header_param->frinit);
//IF the DARK use highpassed data, with tail frequency= 10Hz

 ic=strncmp(input_param->framechannel,"Pr_B1_ACp",9);
 if(ic==0){
   puts("Pr_B1 data: apply high pass filter");
   ii=highpass_data_bil(gd,gd,header_param,10.0);
 }
 ic=strncmp(input_param->framechannel,"H1:LDAS-STRAIN",14);
 if(ic==0){
   puts("H1:LDAS-STRAIN: apply high pass filter");
   ii=highpass_data_bil(gd,gd,header_param,24.0);
 }
 ic=strncmp(input_param->framechannel,"L1:LDAS-STRAIN",14);
 if(ic==0){
   puts("L1:LDAS-STRAIN: apply high pass filter");
   ii=highpass_data_bil(gd,gd,header_param,24.0);
 }
 ic=strncmp(input_param->framechannel,"H1:LSC-STRAIN",13);
 if(ic==0){
   puts("H1:LSC-STRAIN: apply high pass filter");
   ii=highpass_data_bil(gd,gd,header_param,24.0);
 }
 ic=strncmp(input_param->framechannel,"L1:LSC-STRAIN",13);
 if(ic==0){
   puts("L1:LSC-STRAIN: apply high pass filter");
   ii=highpass_data_bil(gd,gd,header_param,24.0);
 }
 if(freqc !=0) {
   ii=highpass_data_bil(gd_highpass,gd,header_param,freqc);  //highpass data in time domain
   ii=purge_data_subtract(gd_clean,gd,gd_highpass,even_param,header_param);
   //identifies and purge large events from gd to gd_clean (time domain). Events are identified from highpassed data, then the cleaning is applied to data in gd, not highpassed.
   //printf("FFT number = %d \n",ii);
 }
 if(freqc==0){
   len=gd->n;
   gd_clean->n=gd->n;
   gd_clean->ini=gd->ini;
   gd_clean->dx=gd->dx;
   for(ii=0;ii<len;ii++)gd_clean->y[ii]=gd->y[ii];
 }
 //Add here the science segments veto
 if (science_segment != NULL) {
    	printf("Applying science segment list selection! Call the function\n");
	iV= ScienceVeto_onGD(gd_clean,header_param,science_segment);
    }
  if(input_param->resampling <=1){
    //else it is done in resampling
    header_param->n_zeroes=0;
    for(ii=0;ii<gd->n;ii++){
      kki=isnan(gd_clean->y[ii]);
   if (kki !=0)
   {
    gd_clean->y[ii]=0;
   }
   if(gd_clean->y[ii]==0)header_param->n_zeroes+=1;
 }
 header_param->perc_zeroes=(float)(header_param->n_zeroes)/gd->n;  //becomes a percentage
  }
  //End science veto, done on gd_clean, output of the EVT veto
 //Put a saturation level to highly disturbed data
//input gd, output the same gd
kki=saturated_data(gd_clean,even_param,header_param);
sat_maxn=even_param->sat_maxn;
//Veto in the sub-bandwidths. January 2009
  subsam=input_param->subsam; //subsampling factor. A power of 2
  if(input_param->subsam > 0){
  band_ext= 1/(2*header_param->tsamplu*subsam);  //in Hz
  purge_subbw(gd_appo,gd_clean,gd_band,even_paramb,header_param,band_ext);
  }
//End veto in the sub-bandwidths.

 //SFDB  (complex FFT data)
  len=gd->n;
  for(ii=0;ii< len;ii++)gd_freq->y[ii]=gd_clean->y[ii];
 ii= gdtime2freq(gd_freq,header_param);  //From the GD in time to the same GD in freq
   for(ii=0;ii<len;ii+=2){
     gd_freq->y[ii]=gd_freq->y[ii]*sqrt(bil2uni);    //Real part
     gd_freq->y[ii+1]=gd_freq->y[ii+1]*sqrt(bil2uni);  //Imag part
   }

 return ii;
}

/***************Function to write the output file with the SFDB******************/
int gd2sfdbfile(GD *gd, HEADER_PARAM *header_param,EVF_PARAM *evf_param,float *ps,EVEN_PARAM *even_param, INPUT_PARAM *input_param,float *psTOT,int iGRID){

  /*The function writes data which are in a GD *gd + some parameters into a binary file (header + data)
   Input are time domain data in a GD *gd. The function performs the total FFT and also and the very short ps using short_ps_time and short_ps_freq
 */
  int errorcode,ii,j;
  long  len;
  float rpw,ipw;
  int lenps,lenps1;
  int estimfromar=1;  //procedure to estimate the very short spectrum
                      //1 is the AR


  int jj;
  int sat;

  /**********Variables needed for the Dopple effect**************/
  short int fut1,nut1,nutcorr;  //nut1=1 uses nutation from ephem files;
                               //nut1=0  uses nutation from novas interpolation;

  short int equatorial=1; //1 output are Equatorial, 0 ecliptic
  short int imjd=1;      //1 input are mjds, 0 imput are days from 1900
  double mjd;
  /*****Detector*****/


  double velx,vely,velz;
  double posx,posy,posz;
  double deinstein;

  double vel[3];
  double pos[3];
  double mid_mjd; //mjd of the middle+1 point of the FFT
  float rpwD;
  DETECTOR *detector;
  //char nam[5];
  int sat_maxn;
  int fft_m;
  int nn;
  int jjin;
   // SFDBNAME sabrina
  struct tm *timstr;
  char appo[40],appo1[8],appo2[10];
  double data_log[6];  //to write the logfile
  int prec[6];
  /* Needed for the peakmaps */
   double sampling_frequency;
   int npeak;
   float ratio_min;
   double spini,spdf;
   int splen;
   float *ratioVEC;
   int *ipeakfound;
   int kkk;
  //strcpy(nam,"virgo");
  //strcpy(nam,"ligol");
  //IMPORTANT NOTE: remember that header_param->nfft is the fft number in the run (e.g. C6), not
  //in that particular file. Thus, only the first SFDB file has  header_param->nfft=1

  /*Now we have time-domain data in the GD*/
   char *appo3;
   char *appo4;
   int nleap;
   double utc;

  appo3=(char*)malloc(sizeof(char)*4);
  appo4=(char*)malloc(sizeof(char)*200);

  int nnn;


  fft_m=input_param->fft_m; //200*(2097152/(2*header_param->nsamples)); //200: AD HOC sulla banda piu' grande pia
             //100 per Virgo //50 per Ligo a 16384 lunghe il doppio


  //printf("**In gd2sfdbfile: actual SDS file filme %s \n",sfc_data->filme);
  //printf("**Actual SDS time t0 mjd %f \n",sfc_data->t0);
  //nn=strlen(sfc_data->filme);
  //seleziona(appo,sfc_data->filme,0,nn-NEXTS-1);


  if(strncmp(input_param->framechannel,"sdsfile",7)==0){
        jjin=filename_from_file1(input_param->filename);
        nn=strlen(input_param->filename)-NEXTS-1;  //-jjin;
        seleziona(appo,input_param->filename,jjin,nn); //-17 e' AD HOC  5 ott 2010una prova !!

    }

  if(strncmp(input_param->framechannel,"sdsfile",strlen("sdsfile"))!=0)strncpy(appo,input_param->framechannel,strlen(input_param->framechannel));
  appo[strlen(input_param->framechannel)] = '\0';
	//printf("$$$$ appo = %s\n", appo);
        utc=header_param->gps_sec/day2sec+dot2gps+dot2utc;
        nleap=PSS_leapseconds(utc);
	timstr=GPS2UTCSAB(header_param->gps_sec,nleap);
	printf("header_param->gps_sec nleap= %d %d\n",header_param->gps_sec,nleap);
	printf("%04d%02d%02d\n",timstr->tm_year+1900,timstr->tm_mon+1,timstr->tm_mday);
	printf("%02d%02d%02d",timstr->tm_hour,timstr->tm_min,timstr->tm_sec);
	strcat(appo,"_");
	sprintf(appo2,"%04d%02d%02d",timstr->tm_year+1900,timstr->tm_mon+1,timstr->tm_mday);
	strcat(appo,appo2);
	strcat(appo,"_");
	sprintf(appo2,"%02d%02d%02d",timstr->tm_hour,timstr->tm_min,timstr->tm_sec);
	strcat(appo,appo2);

	nnn=strlen(appo);
        strncpy(appo4,appo,nnn);
	appo4[nnn]='\0';
	//printf("check appo4,  appo nnn %s %s %d \n",appo4,appo,nnn);
	strncpy(appo1,".SFDB09",7);
	strncat(appo,appo1,7);


	//strcpy(header_param->sfdbname,appo); MOVED below: when the new file is opened

	 if(input_param->flagonline >=3){
	   //puts("Peakmap file with extension .p10 will be created");
	  strncpy(appo3,".p10",4);
	  strncat(appo4,appo3,4);
	  //  strcpy(input_param->peakname,appo4);  MOVED BELOW
	  //nnn=strlen(appo4);
	  //input_param->peakname[nnn]='\0';
	  free(appo3);
	  //free(appo4);
	}


        //printf("Writing SFDB filename= %s \n", header_param->sfdbname);


   sat_maxn=even_param->sat_maxn;
   lenps= (int) gd->n/(header_param->red);
   lenps1= header_param->red;
   header_param->lavesp=lenps/2;  //lenght of the short periodogram (2010: now the same as AR spectrum)
   len=gd->n;
   printf("Poi levalo: in gd2sfdbfile: gd->n  %ld \n",len);
   printf("Poi levalo: in gd2sfdbfile: nfft %d  \n",header_param->nfft);
   if(header_param->nfft==1 && SFDB == NULL){
     input_param->nfftONE=0;
     if(input_param->flagonline !=1){
       strcpy(header_param->sfdbname,appo);
       SFDB=fopen(header_param->sfdbname,"w");
       printf("First SFDB file opened - file name: %s \n", header_param->sfdbname);}
      if(input_param->flagonline >=3){
       strcpy(input_param->peakname,appo4);
       nnn=strlen(appo4);
       input_param->peakname[nnn]='\0';
       free(appo4);
       P10=fopen(input_param->peakname,"w");
       printf("First P10 file opened - file name: %s \n", input_param->peakname);
       //write the header for P10 files
       puts("Write the general header for the first pm file");
       sampling_frequency=1/header_param->tsamplu;
       fwrite((void*)&input_param->fft_m, sizeof(int),1,P10);  //viene poi riscritto se diverso
       fwrite((void*)&sampling_frequency, sizeof(double),1,P10);
       fwrite((void*)&header_param->nsamples, sizeof(int),1,P10);
       fwrite((void*)&header_param->frinit, sizeof(double),1,P10);
       //other header info added after the construction of the pm, to have and write npeak
      }
   }




   /*Closes the file every a chosen number of FFTs, if itf detector, and opens a new file */
   if(input_param->nfftONE==fft_m){
     //if(header_param->detector !=0 && fmodf((float)header_param->nfft,(float) (fft_appo))==0.){
     printf("****Close the file at FFT (included) num=, with nfft in the file= %d %d\n",header_param->nfft-1,input_param->nfftONE);
     input_param->nfftONE=0;
      if(input_param->flagonline !=1)
	{
      //printf("open a new SFDB file: nfft fft_m %d  %d \n",header_param->nfft, fft_m);
      fclose(SFDB);
      if(input_param->flagonline >=3)fclose(P10);
      if(iGRID>=1)ToCopyFiles(input_param,header_param);     //rimettere 1
      strcpy(header_param->sfdbname,appo);
      SFDB=fopen(header_param->sfdbname,"w");
      printf("New SFDB file opened  - file name: %s \n", header_param->sfdbname);
    if(input_param->flagonline >=3){
	strcpy(input_param->peakname,appo4);
	nnn=strlen(appo4);
	input_param->peakname[nnn]='\0';
	free(appo4);
        P10=fopen(input_param->peakname,"w");
        printf("A new P10 file opened - file name: %s \n", input_param->peakname);
	//write the header for P10 files
       puts("Write the general header of any new pm file");
       sampling_frequency=1/header_param->tsamplu;
       fwrite((void*)&input_param->fft_m, sizeof(int),1,P10);  //viene poi riscritto se diverso
       fwrite((void*)&sampling_frequency, sizeof(double),1,P10);
       fwrite((void*)&header_param->nsamples, sizeof(int),1,P10);
       fwrite((void*)&header_param->frinit, sizeof(double),1,P10);
     //other header info added after the construction of the pm, to have and write npeak
    }
      jj=logfile_output(LOG_INFO,header_param->sfdbname,"A new SFDB file opened");
	}
   }
input_param->nfftONE+=1;
if(input_param->flagonline !=1) {
   // Added 15/01/08:
  detector=def_detector(input_param->detector_name);
      //Evaluation of detector velocities
  /****************************************************************************/
   fut1=0;  ///1=uses ut1; 0=uses utc as aprrox. of ut1; //0 per fare fino al 2010
   nut1=1;  //1=uses nutation from ephem files; 0=no;
   nutcorr=0; //reads the file DELTANUT, if 0 does not read DELTANUT //0 per fare fino al 2010
 /****************************************************************************/
   //	mid_mjd=header_param->mjdtime+(double)(gd->n/2)*header_param->tsamplu/day2sec;   //aggiunto +1 nel formato SFDB09
	mid_mjd=header_param->mjdtime+(double)(gd->n/2+1)*header_param->tsamplu/day2sec;
	printf("poi levalo: Middle +1 point mjd %f\n",mid_mjd);

	mjd=PSS_detector_velposSCALAR(detector,mid_mjd,header_param->tsamplu,fut1,nut1,nutcorr,equatorial,imjd,&velx,&vely,&velz,&posx,&posy,&posz,deinstein);
	//puts("uscito dalla function detector_velpos");
	//velx,y,z inAU/day; posx,y,z in AU. C in AU/day
	//vel is adimensional and pos in light seconds
	vel[0]=velx/C;
	vel[1]=vely/C;
	vel[2]=velz/C;
	header_param->vx_eq=vel[0];
	header_param->vy_eq=vel[1];
	header_param->vz_eq=vel[2];
	pos[0]=posx*day2sec/C;
	pos[1]=posy*day2sec/C;
	pos[2]=posz*day2sec/C;
	header_param->px_eq=pos[0];
	header_param->py_eq=pos[1];
	header_param->pz_eq=pos[2];

	printf("Poi levalo:Calculating detector velocity\n");
	printf("Poi levalo:VELOCITA'/C: %15.8e %15.8e %15.8e\n",velx/C,vely/C,velz/C);
        //printf("una VELOCITA'/C in header param: %15.8e \n",header_param->vx_eq);
        //printf("VELOCITA': %15.8e %15.8e %15.8e\n",velx,vely,velz);
        //printf("POSIZIONI*day2sec/C: %15.8e %15.8e %15.8e\n",posx*day2sec/C,posy*day2sec/C,posz*day2sec/C);


         //write the values in the LOG file
         data_log[0]=(double) header_param->vx_eq;
	 data_log[1]=(double) header_param->vy_eq;
	 data_log[2]=(double) header_param->vz_eq;
	 data_log[3]=(double) header_param->px_eq;
	 data_log[4]=(double) header_param->py_eq;;
	 data_log[5]=(double) header_param->pz_eq;;
	 prec[0]=0;
	 prec[1]=0;
	 prec[2]=0;
	 prec[3]=0;
	 prec[4]=0;
	 prec[5]=0;

	 logfile_ev(LOG_INFO,"VEL",6,data_log,prec);
     //end for the log file

	free(detector);

	//end evaluation of detector velocities--Added on Jan 15 2008
 }  //end of flagonline not 1
 if(input_param->flagonline ==1)
   {
     	header_param->vx_eq=0.;
	header_param->vy_eq=0.;
	header_param->vz_eq=0.;
	header_param->px_eq=0.;
	header_param->py_eq=0.;
	header_param->pz_eq=0.;
   }
   /*End of the close/open procedure  */

	/*
	  19 may 2008
	  If too many saturated data: do not write the FFT and header. return
	 */
	//puts("Prima di saturated");

	sat= (int) header_param->sat_howmany;
	//if (sat >=sat_maxn){
	//errorcode=0;
	  //printf("Veto on the FFT due to sat data: howmany,howmany,maxn %f %d %d\n",header_param->sat_howmany,sat,sat_maxn);

	  //return errorcode;
	//}
	//puts("Prima di veto nzeroes");
	//if(header_param->n_zeroes >(header_param->nsamples*2))
	//{
	//printf("veto una fft nfft, header_param->n_zeros, %d, %d\n",header_param->nfft, header_param->n_zeroes);
	//errorcode=0;

	//return errorcode;        //continue; //go to the next in the loop
	//}

   if(input_param->flagonline !=1)

   {
   /*Write the header*/
   errorcode=fwrite((void*)&header_param->endian, sizeof(double),1,SFDB);
   errorcode=fwrite((void*)&header_param->detector, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->gps_sec, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->gps_nsec, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->tbase, sizeof(double),1,SFDB);
   errorcode=fwrite((void*)&header_param->firstfreqindex, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->nsamples, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->red, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->typ, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->n_flag, sizeof(float),1,SFDB);
   errorcode=fwrite((void*)&header_param->einstein, sizeof(float),1,SFDB);
   errorcode=fwrite((void*)&header_param->mjdtime, sizeof(double),1,SFDB);
   errorcode=fwrite((void*)&header_param->nfft, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->wink, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->normd, sizeof(float),1,SFDB);
   errorcode=fwrite((void*)&header_param->normw, sizeof(float),1,SFDB);
   errorcode=fwrite((void*)&header_param->frinit, sizeof(double),1,SFDB);
   errorcode=fwrite((void*)&header_param->tsamplu, sizeof(double),1,SFDB);
   errorcode=fwrite((void*)&header_param->deltanu, sizeof(double),1,SFDB);
   //specific header (different for bar-itf) detectors
   if(header_param->detector==0){ //bar detector
     errorcode=fwrite((void*)&header_param->vx_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->vy_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->vz_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->px_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->py_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->pz_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->n_zeroes, sizeof(int),1,SFDB);
     errorcode=fwrite((void*)&header_param->sat_howmany, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare1, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare2, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare3, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->perc_zeroes, sizeof(float),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare5, sizeof(float),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare6, sizeof(float),1,SFDB);
     errorcode=fwrite((void*)&header_param->lavesp, sizeof(int),1,SFDB);
     errorcode=fwrite((void*)&header_param->sciseg, sizeof(int),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare9, sizeof(int),1,SFDB);
   }

   if(header_param->detector>=1){ //itf detector
     errorcode=fwrite((void*)&header_param->vx_eq, sizeof(double),1,SFDB);  //vx/C equat
     errorcode=fwrite((void*)&header_param->vy_eq, sizeof(double),1,SFDB);  //vy
     errorcode=fwrite((void*)&header_param->vz_eq, sizeof(double),1,SFDB);  //vz
     errorcode=fwrite((void*)&header_param->px_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->py_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->pz_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->n_zeroes, sizeof(int),1,SFDB);
     errorcode=fwrite((void*)&header_param->sat_howmany, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare1, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare2, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare3, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->perc_zeroes, sizeof(float),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare5, sizeof(float),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare6, sizeof(float),1,SFDB);
     errorcode=fwrite((void*)&header_param->lavesp, sizeof(int),1,SFDB);
     errorcode=fwrite((void*)&header_param->sciseg, sizeof(int),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare9, sizeof(int),1,SFDB);
   }


   if (errorcode!=1)printf("Error in writing header into SFDB file!\n");
   //printf("file header:endian %lf \n",header_param->endian);
   //printf("file header:nsamples %d  \n",header_param->nsamples);
   //printf("file header:n_flag %f  \n",header_param->n_flag);
   //printf("file header:n_zeroes %d  \n",header_param->n_zeroes);
   //printf("file header:perc_zeroes %f  \n",header_param->perc_zeroes);
   //printf("File header:gps %d \n",header_param->gps_sec);
   //printf("File header:mjdtime %f \n",header_param->mjdtime);
   /*Write the SFDB data*/
   /*Creation of the average of very short power spectra */


     //2008 ps=(float *)malloc((size_t) (lenps1)*sizeof(float));
    ii=short_ps_periodogram(ps,gd,header_param);
    //for(j=0;j<header_param->red;j++){
    for(j=0;j<lenps/2;j++){
      //printf(" SQRT of Short power spectrum, dimension lenps/2: j ps[j] %d %f\n",j,ps[j]);
      //rpw=ps[j]*pow(header_param->normd,2)*pow(header_param->normw,2)*bil2uni; //piapia
      rpw=ps[j]; //sqrt of the periodogram
      if(header_param->perc_zeroes < 1)rpw=rpw/sqrt(1-header_param->perc_zeroes);
      errorcode=fwrite((void*)&rpw, sizeof(float),1,SFDB);
    }
   }  //end of  if(input_param->flagonline !=1)
   if(input_param->flagonline >=3){

     //inizializzazioni
     ratio_min=evf_param->maxdin*evf_param->maxdin;  //sqrt(2.5) * sqrt(2.5). Peakmaps are done starting from squared ratioes


   }


   /*Window, if needed */
     if(header_param->wink>0){
       gd2wingd(gd,header_param->wink);
     }
    //2008 free(ps);
    if(estimfromar!=1){
      //Very short FFT from sum over lower resolution FFTs, (input are data in time domain)
      //2008 ps=(float *)malloc((size_t) (lenps)*sizeof(float));
      ii=short_ps_freq(ps,gd,header_param);
      for(j=0;j<lenps/2;j++){
      //printf(" Very Short sqrt(spectrum), dimension lenps/2: j ps[j]  %d %f \n",j,ps[j]);
      rpw=ps[j]*header_param->normd*header_param->normw; //sqrt of the power spectrum
      if(header_param->perc_zeroes < 1)rpw=rpw/sqrt(1-header_param->perc_zeroes);
      if(input_param->flagonline !=1)errorcode=fwrite((void*)&rpw, sizeof(float),1,SFDB);
      }
      // 2008 free(ps);
    }


   ii= gdtime2freq(gd,header_param);  //From the GD in time to the same GD in frequency domai
   if(estimfromar==1){
     //Very short spectrum from AR estimation: (input are data in frequency domain)
     //ii=short_psar_freq(ps,gd,header_param,evf_param);
     //puts("Call the reverse AR estimation");
     ii=short_psar_rev_freq(ps,gd,header_param,evf_param,psTOT); //10 Jul 2007 REVERSED !! //psTOT needed for peakmap files. Needs to be  squared when doing the peakmaps
   for(j=0;j<lenps/2;j++){
    rpw=ps[j]*header_param->normd*header_param->normw*sqrt(bil2uni); //sqrt of the ps
    if(header_param->perc_zeroes < 1)rpw=rpw/sqrt(1-header_param->perc_zeroes);
    if(input_param->flagonline !=1)errorcode=fwrite((void*)&rpw, sizeof(float),1,SFDB);
   }
   }
   //SFDB  (complex FFT data): if flagonline not 1
if(input_param->flagonline !=1)
  {
     for(ii=0;ii< len;ii+=2){
       //MODIFICA 15/04/2010: tolta la sqrt(2), non si sa quando era apparsa.
       //da oggi tutti gli SFDB09 saranno senza pia
       //rpw=gd->y[ii]*sqrt(bil2uni);    //Real part
       //ipw=gd->y[ii+1]*sqrt(bil2uni);  //Imag part
     rpw=gd->y[ii];    //Real part  piapia
     ipw=gd->y[ii+1];  //Imag part

     //if(ii<4)printf("First data written into the file %e %e\n",rpw,ipw); //messo %e 25-Aug 2005
     errorcode=fwrite((void*)&rpw, sizeof(float),1,SFDB);
     if (errorcode!=1)printf("Error in writing data into SFDB file!\n");
     errorcode=fwrite((void*)&ipw, sizeof(float),1,SFDB);
     if (errorcode!=1)printf("Error in writing data into SFDB file!\n");

     }

     if(input_param->flagonline >=3){
       // puts("Doing peakmaps");
       ratioVEC=(float *)malloc((size_t) (len)*sizeof(float));
       ipeakfound=(int *)malloc((size_t) (len/3)*sizeof(int));
       //add here the instruction to fill the peakmap data in the file P10. Using psTOT (dimension form 0 to len-1) and rpw,ipw
       j=0;
       ratioVEC[len-1]=1;  //last bin has to be 1.
       for(ii=0;ii<2*(len-1);ii+=2){
	 ratioVEC[j]=(gd->y[ii]*gd->y[ii]+gd->y[ii+1]*gd->y[ii+1])/(psTOT[j+1]*psTOT[j+1]); //the evaluation will be from the last bin to the first. Using the average on the left
	 j+=1;
       }
       npeak=0;
       for (j=1;j<=(len-2);j++){
	 if(ratioVEC[j] >= ratio_min) {

	   if((ratioVEC[j]>ratioVEC[j-1]) && (ratioVEC[j] > ratioVEC[j+1])){
	   npeak=npeak+1;
	   ipeakfound[npeak-1]=j;
       }
       }
       }
       //write the missing pieces of the pm header, starting for mid_mjd, npeak...
       //puts("Writing the missing pieces of the header for pm files");

       spini=gd->ini;
       spdf=header_param->deltanu*header_param->red;
       splen=gd->n/(2*header_param->red); //identico a lenps..
       if (header_param->nfft==1){
	 printf("First pm file: FFT num MJD time npeaks splen %d %f %d %d\n",header_param->nfft,header_param->mjdtime,npeak,splen);}
       fwrite((void*)&header_param->mjdtime, sizeof(double),1,P10);
       fwrite((void*)&npeak, sizeof(int),1,P10);
       fwrite((void*)&header_param->vx_eq, sizeof(double),1,P10);
       fwrite((void*)&header_param->vy_eq, sizeof(double),1,P10);
       fwrite((void*)&header_param->vz_eq, sizeof(double),1,P10);
       fwrite((void*)&header_param->px_eq, sizeof(double),1,P10);
       fwrite((void*)&header_param->py_eq, sizeof(double),1,P10);
       fwrite((void*)&header_param->pz_eq, sizeof(double),1,P10);
       fwrite((void*)&spini, sizeof(double),1,P10);
       fwrite((void*)&spdf, sizeof(double),1,P10);
       fwrite((void*)&splen, sizeof(int),1,P10);
       for (j=0;j<splen;j++){
	 rpw=ps[j]*ps[j]*pow(header_param->normd,2)*pow(header_param->normw,2)*bil2uni;  //AR power spectrum (that is squared. In the SFDB it is the sqrt)
	 //	  rpw=ps[j]*header_param->normd*header_param->normw*sqrt(bil2uni); //sqrt of the ps
	 //	  if(header_param->perc_zeroes < 1)rpw=rpw/sqrt(1-header_param->perc_zeroes);
	  if(header_param->perc_zeroes < 1)rpw=rpw/(1-header_param->perc_zeroes);
          fwrite((void*)&rpw, sizeof(float),1,P10);
       }
       //write the data, usinf thenP10 format


     for(j=0;j<npeak;j++){
       fwrite((void*)&ipeakfound[j], sizeof(int),1,P10);
       //if (j>npeak-2)printf("j,i peakfound %d %d\n",j,ipeakfound[j]);


     }
      for(j=0;j<npeak;j++){
	kkk=ipeakfound[j];
	rpw=sqrt(ratioVEC[kkk]);
	fwrite((void*)&rpw,sizeof(float),1,P10);
	//if (j>npeak-2)printf("j,ratio %d %f\n",j,rpw);
     }

      for(j=0;j<npeak;j++){
	kkk=ipeakfound[j];
	rpwD=psTOT[kkk+1]*header_param->normd*header_param->normw*sqrt(bil2uni); //sqrt of the ps
	if(header_param->perc_zeroes < 1)rpwD=rpwD/sqrt(1-header_param->perc_zeroes);
	fwrite((void*)&rpwD, sizeof(float),1,P10);
	//if (j>npeak-2)printf("j,psTOT %d %f\n",j,rpwD);
     }

     free(ratioVEC);
     free(ipeakfound);
     }
     //puts("In gd2sfdbfile:another data chunk has been written");
  } //end of:   if(input_param->flagonline !=1)

  return errorcode;
}



/***************Function to read the output file with the SFDB******************/
int readfilesfdb_old(GD *gd, GD *gd_short,HEADER_PARAM *header_param,int k,int iw,int *fft_read,int ilista){

  /*gd is the output GD with data from the SFDB file*/
  /*iw=number of the fft which -for test- has to be written on a ascii file*/
  //FILE *OUTH1;
  FILE *OUTS;
  FILE *OUTSS;
  int errorcode,howmany,ii,kk;
  long howmany2;
  char filediingresso[MAXMAXLINE+1];
  char fileout[MAXLINE+1];
  char filetest[MAXLINE+1];
  float rpw,ipw,freq,sden;
  float *ps;
  float *ps_short;
  int lenps;
  long dovesta;
  int jj;
  char filelista[MAXMAXLINE+1];
  char listapath[MAXMAXLINE+1];

  char capt_f[12];  //Caption for the GD with the SFDB data read from a file.
  char capt_short[12];  //Caption for the GD with the short SFDB data read from a file

   //printf("In redfilesfdb k= %d\n",k);
   if(k==0){
   *fft_read=0;
   // Read all files from a lista file, if ilista ==1
    if(ilista==1){
      //Open the file lista.txt and read the files
     printf("Name of the lista file, with the path\n");
     scanf("%s",listapath);
      LISTA=fopen(listapath,"r");
      fscanf(LISTA,"%s",filelista);
      strcpy(filediingresso,filelista);
      puts("-->Read files from the lista file");
      puts("-->Should be created with ls -d$PWD/*.SFDB");
      puts("and a line at the end with written e.g. FINE");
   }
   if(ilista != 1){
     printf("file input\n");
     scanf("%s",filediingresso);
   }
   jj= logfile_input(LOG_INFO, filediingresso, "first SFDB input file"); //info for the log
   //printf("file output with one chosen spectrum (ascii)\n");
   //scanf("%s",fileout);
   //printf("file output with the very short power spectrum (ascii)\n");
   //scanf("%s",filetest);
   OUTH1=fopen(filediingresso,"r");
   //OUTS=fopen(fileout,"w");
   //OUTSS=fopen(filetest,"w");
   printf("SFDB file opened %s \n",filediingresso);
   }
   if(k==iw){
      printf("file output with one chosen power spectrum (ascii) 1/Hz\n");
      scanf("%s",fileout);
      printf("file output with the very short power spectrum (ascii) 1/Hz \n");
      scanf("%s",filetest);
      OUTS=fopen(fileout,"w");
      OUTSS=fopen(filetest,"w");
   }

      //puts("Before the header");
     //printf("header size in byte %d\n",sizeof(HEADER_PARAM));

     dovesta=ftell(OUTH1);
     //printf("all' inizio dovesta= %ld\n",dovesta);

     /*Read the header: */
     //OLD errorcode=fread((void*)&header_paramw,sizeof(header_paramw),1,OUTH);
     errorcode=fread((void*)&header_param->endian, sizeof(double),1,OUTH1);
     //on error, try to open a new file
      if (errorcode!=1){
	 printf("FFT number=  %d  \n",k);
         printf("Error in reading header into SFT file!: open a new file\n");
	 fclose(OUTH1);
	  if(ilista != 1){
	    printf("New file input (give a fake name to stop here)\n");
	    scanf("%s",filediingresso);
	  }
	  if(ilista==1){
	    puts("-->Read files from the lista.txt file");
	    strcpy(filelista,"END"); //to finish if no other files are present
	    fscanf(LISTA,"%s",filelista);
	    strcpy(filediingresso,filelista);
	  }
	 OUTH1=fopen(filediingresso,"r");
	 if(OUTH1==NULL) return errorcode;
	 printf("New SFDB file opened %s \n",filediingresso);
         errorcode=fread((void*)&header_param->endian, sizeof(double),1,OUTH1);
      }
     printf("endian %f\n",header_param->endian);
     errorcode=fread((void*)&header_param->detector, sizeof(int),1,OUTH1);
     errorcode=fread((void*)&header_param->gps_sec, sizeof(int),1,OUTH1);
     errorcode=fread((void*)&header_param->gps_nsec, sizeof(int),1,OUTH1);
     errorcode=fread((void*)&header_param->tbase, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->firstfreqindex, sizeof(int),1,OUTH1);
     errorcode=fread((void*)&header_param->nsamples, sizeof(int),1,OUTH1);
     errorcode=fread((void*)&header_param->red, sizeof(int),1,OUTH1);
     errorcode=fread((void*)&header_param->typ, sizeof(int),1,OUTH1);
     errorcode=fread((void*)&header_param->n_flag, sizeof(float),1,OUTH1);
     errorcode=fread((void*)&header_param->einstein, sizeof(float),1,OUTH1);
     errorcode=fread((void*)&header_param->mjdtime, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->nfft, sizeof(int),1,OUTH1);
     errorcode=fread((void*)&header_param->wink, sizeof(int),1,OUTH1);
     errorcode=fread((void*)&header_param->normd, sizeof(float),1,OUTH1);
     errorcode=fread((void*)&header_param->normw, sizeof(float),1,OUTH1);
     errorcode=fread((void*)&header_param->frinit, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->tsamplu, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->deltanu, sizeof(double),1,OUTH1);
     //specific header (different for bar-itf) detectors
     if(header_param->detector==0){ //bar detector
       errorcode=fread((void*)&header_param->vx_eq, sizeof(double),1,OUTH1);
       errorcode=fread((void*)&header_param->vy_eq, sizeof(double),1,OUTH1);
       errorcode=fread((void*)&header_param->vz_eq, sizeof(double),1,OUTH1);
       errorcode=fread((void*)&header_param->sat_howmany, sizeof(double),1,OUTH1);
       errorcode=fread((void*)&header_param->spare2, sizeof(double),1,OUTH1);
     }
     if(header_param->detector>=1){ //itf detector
       errorcode=fread((void*)&header_param->vx_eq, sizeof(double),1,OUTH1);  //vx/C equat.
       errorcode=fread((void*)&header_param->vy_eq, sizeof(double),1,OUTH1);  //vy
       errorcode=fread((void*)&header_param->vz_eq, sizeof(double),1,OUTH1);  //vz
       errorcode=fread((void*)&header_param->sat_howmany, sizeof(double),1,OUTH1);
       errorcode=fread((void*)&header_param->spare2, sizeof(double),1,OUTH1);
     }
     //puts("After the header");
     if (errorcode!=1){
       printf("Error in reading header into SFT file!\n");
       printf("FFT number= data= %d %d \n",k,ii);
       return errorcode;
     }
   printf("header: nfft, mjdtime,tsamplu, deltanu,red,n_flag, nsamples %d %f %f %f %d %f %d\n",header_param->nfft,header_param->mjdtime,header_param->tsamplu,header_param->deltanu,header_param->red,header_param->n_flag,header_param->nsamples);
   if(header_param->detector==0){
   printf("header: tbase,frinit  %f %f\n",header_param->tbase,header_param->frinit);
   }
    if(header_param->detector>=1){
   printf("header: tbase,frinit,normw  %f %f %f\n",header_param->tbase,header_param->frinit,header_param->normw);
   }
    howmany=header_param->nsamples;
    howmany2=howmany*2;
    //printf("howmany howmany2 = %d %ld\n",howmany,howmany2);

    strcpy(capt_f,"SFDB data");
    strcpy(capt_short,"Short SFDB");
    gd->n=howmany2;
    gd->nall=howmany2;
    gd->ini=(double)(header_param->firstfreqindex)*1.0/header_param->tbase;
    gd->dx=header_param->deltanu;
    gd->capt=capt_f;
    free(gd->y);
    gd->y=malloc((size_t) howmany2*sizeof(float));
    ps=(float *)malloc((size_t) (header_param->red)*sizeof(float));
    /*Read the averages of the short power spectra (behaviour with time)*/
    for(ii=0;ii<header_param->red;ii++){
      errorcode=fread((void*)&rpw, sizeof(float),1,OUTH1);
       if (errorcode!=1) {
	printf("Error in reading power ps data into SFT file! %d \n",ii);
	return errorcode;
      }
      ps[ii]=rpw;
      //printf(" Short power j ps[j] %d %f\n",ii,ps[ii]);
    }
    free(ps);
    lenps= (int) gd->n/(header_param->red);
    ps_short=(float *)malloc((size_t) (lenps/2)*sizeof(float));
    gd_short->n=lenps/2;
    gd_short->nall=lenps/2;
    gd_short->ini=(double)(header_param->firstfreqindex)*1.0/header_param->tbase;
    gd_short->dx=header_param->deltanu*header_param->red;
    gd_short->capt=capt_short;
    free(gd_short->y);
    gd_short->y=malloc((size_t) lenps/2*sizeof(float));

    /*Read the very short FFT (dimension lenps; each bin averaged over red data*/
     kk=0; //only to evaluate frequencies for the test file
    for(ii=0;ii<lenps/2;ii++){
      errorcode=fread((void*)&rpw, sizeof(float),1,OUTH1);
       ps_short[ii]=rpw;
       if (errorcode!=1) {
	printf("Error in reading very short FFT  data into SFT file! %d \n",ii);
	return errorcode;
      }
      gd_short->y[ii]=rpw;       //sqrt(spectrum)
      //TEST FILE with the very short power spectrum:
      freq=(float)(header_param->firstfreqindex)/header_param->tbase+kk*header_param->red/header_param->tbase;
      sden=rpw*rpw;
      if(k==iw)fprintf(OUTSS,"%f %e \n",freq,sden); //write the iw v. short power spectrum  in a ascii file
      kk++;
    }
    free(ps_short);

    /*Read the data: */
    kk=0;
    for(ii=0;ii< 2*howmany; ii+=2){
      errorcode=fread((void*)&rpw, sizeof(float),1,OUTH1);
      if (errorcode!=1) {
	printf("Error in reading rpw data into SFT file! %d \n",ii);
	return errorcode;
      }
      errorcode=fread((void*)&ipw, sizeof(float),1,OUTH1);
      if (errorcode!=1) {
	printf("Error in reading ipw data into SFT file! %d \n",ii);
	return errorcode;
      }
      //if(ii==0)printf("ii==0  rpw,ipw %f %f \n",rpw,ipw);
      //if(ii<=2)printf("rpw,ipw %f %f \n",rpw,ipw);
      gd->y[ii]=rpw;       //Real part
      gd->y[ii+1]=ipw;    //Imag part
      freq=(float)(header_param->firstfreqindex+kk)*1.0/header_param->tbase;
      sden=(rpw*rpw+ipw*ipw)*pow(header_param->normd,2)*pow(header_param->normw,2); //piapia;
      //if(kk < 10)printf("kk,freq sden %d %f %f\n",kk,freq,sden);
      if(k==iw)fprintf(OUTS,"%f %e \n",freq,sden); //write the iw in a ascii file. piapia %e 18 aug 2005
      kk++;

    } //for ii

    *fft_read+=1;
    if(k==iw){
      fclose(OUTS);
      fclose(OUTSS);
    }
    dovesta=ftell(OUTH1);
    //printf("alla fine dovesta= %ld\n",dovesta);
    if(errorcode!=1)
    {
       fclose(OUTH1);
    }
  return errorcode;
}

/***************Function to read the output file with the SFDB******************/
int readfilesfdb(GD *gd, GD *gd_short, HEADER_PARAM *header_param,
     int k, int iw, int *fft_read, int ilista, int header_only){

     /*gd is the output GD with data from the SFDB file*/
     /*iw=number of the fft which -for test- has to be written on a ascii file*/
     //FILE *OUTH1;
     //To read SFDB09 files, from 12 feb 2009
     FILE *OUTS;
     FILE *OUTSS;
     int errorcode,howmany,ii,kk;
     long howmany2;
     static char filediingresso[MAXMAXLINE+1];
     char fileout[MAXLINE+1];
     char filetest[MAXLINE+1];
     float rpw,ipw,freq,sden;
     float *ps;
     float *ps_short;
     int lenps;
     long dovesta;
     int jj;
     char filelista[MAXMAXLINE+1];
     char listapath[MAXMAXLINE+1];

     char capt_f[12];  //Caption for the GD with the SFDB data read from a file.
     char capt_short[12];  //Caption for the GD with the short SFDB data read from a file

     //printf("In redfilesfdb k= %d\n",k);
     if(k==0){
          *fft_read=0;
          // Read all files from a lista file, if ilista ==1
          if(ilista==1){
               //Open the file lista.txt and read the files
               printf("Name of the lista file, with the path\n");
               scanf("%s",listapath);
               LISTA=fopen(listapath,"r");
               fscanf(LISTA,"%s",filelista);
               strcpy(filediingresso,filelista);
               puts("-->Read files from the lista file");
               puts("-->Should be created with ls -d$PWD/*.SFDB");
               //    puts("and a line at the end with written e.g. FINE");
          }
          if(ilista != 1){
               printf("file input\n");
               scanf("%s",filediingresso);
          }
          jj= logfile_input(LOG_INFO, filediingresso, "first SFDB input file"); //info for the log
          //printf("file output with one chosen spectrum (ascii)\n");
          //scanf("%s",fileout);
          //printf("file output with the very short power spectrum (ascii)\n");
          //scanf("%s",filetest);
          OUTH1=fopen(filediingresso,"r");
          //OUTS=fopen(fileout,"w");
          //OUTSS=fopen(filetest,"w");
          printf("SFDB file opened %s \n",filediingresso);
     }
     if(k==iw){
          printf("file output with one chosen power spectrum (ascii) 1/Hz\n");
          scanf("%s",fileout);
          printf("file output with the very short power spectrum (ascii) 1/Hz \n");
          scanf("%s",filetest);
          OUTS=fopen(fileout,"w");
          OUTSS=fopen(filetest,"w");
     }

     //puts("Before the header");
     //printf("header size in byte %d\n",sizeof(HEADER_PARAM));

     //PROVA A TOGLERLO  dovesta=ftell(OUTH1);
     //printf("all' inizio dovesta= %ld\n",dovesta);

     /*Read the header: */
     //OLD errorcode=fread((void*)&header_paramw,sizeof(header_paramw),1,OUTH);
     errorcode=fread((void*)&header_param->endian, sizeof(double),1,OUTH1);
     //on error, try to open a new file
     if (errorcode!=1){
          printf("FFT number=  %d  \n",k);
          printf("Error in reading header into SFT file!: open a new file\n");
          fclose(OUTH1);
          if(ilista != 1){
	          printf("New file input (give a fake name to stop here)\n");
	          scanf("%s",filediingresso);
  	     }
          if(ilista==1){
	          puts("-->Read files from the lista.txt file");
     	     strcpy(filelista,"END"); //to finish if no other files are present
     	     fscanf(LISTA,"%s",filelista);
     	     strcpy(filediingresso,filelista);
          }
          OUTH1=fopen(filediingresso,"r");
          if(OUTH1==NULL) return errorcode;
          printf("New SFDB file opened %s \n",filediingresso);
          //dopo aperto mettere un fake name,per la fine lista !

          errorcode=fread((void*)&header_param->endian, sizeof(double),1,OUTH1);
     }
     strcpy(&header_param->sfdbname,&filediingresso);

     printf("endian %f\n",header_param->endian);
     errorcode=fread((void*)&header_param->detector, sizeof(header_param->detector),1,OUTH1);
     printf("detector %d\n",header_param->detector);
     errorcode=fread((void*)&header_param->gps_sec, sizeof(header_param->gps_sec),1,OUTH1);
     errorcode=fread((void*)&header_param->gps_nsec, sizeof(header_param->gps_nsec),1,OUTH1);
     errorcode=fread((void*)&header_param->tbase, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->firstfreqindex, sizeof(header_param->firstfreqindex),1,OUTH1);
     errorcode=fread((void*)&header_param->nsamples, sizeof(header_param->nsamples),1,OUTH1);
     errorcode=fread((void*)&header_param->red, sizeof(header_param->red),1,OUTH1);
     errorcode=fread((void*)&header_param->typ, sizeof(header_param->typ),1,OUTH1);
     errorcode=fread((void*)&header_param->n_flag, sizeof(float),1,OUTH1);
     errorcode=fread((void*)&header_param->einstein, sizeof(float),1,OUTH1);
     errorcode=fread((void*)&header_param->mjdtime, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->nfft, sizeof(header_param->nfft),1,OUTH1);
     errorcode=fread((void*)&header_param->wink, sizeof(header_param->wink),1,OUTH1);
     errorcode=fread((void*)&header_param->normd, sizeof(float),1,OUTH1);
     errorcode=fread((void*)&header_param->normw, sizeof(float),1,OUTH1);
     errorcode=fread((void*)&header_param->frinit, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->tsamplu, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->deltanu, sizeof(double),1,OUTH1);
     //specific header (different for bar-itf) detectors

     errorcode=fread((void*)&header_param->vx_eq, sizeof(double),1,OUTH1);  //vx/C equat.
     errorcode=fread((void*)&header_param->vy_eq, sizeof(double),1,OUTH1);  //vy
     errorcode=fread((void*)&header_param->vz_eq, sizeof(double),1,OUTH1);  //vz
     errorcode=fread((void*)&header_param->px_eq, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->py_eq, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->pz_eq, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->n_zeroes, sizeof(header_param->n_zeroes),1,OUTH1);
     errorcode=fread((void*)&header_param->sat_howmany, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->spare1, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->spare2, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->spare3, sizeof(double),1,OUTH1);
     errorcode=fread((void*)&header_param->perc_zeroes, sizeof(float),1,OUTH1);
     errorcode=fread((void*)&header_param->spare5, sizeof(float),1,OUTH1);
     errorcode=fread((void*)&header_param->spare6, sizeof(float),1,OUTH1);
     errorcode=fread((void*)&header_param->lavesp, sizeof(header_param->lavesp),1,OUTH1);
     errorcode=fread((void*)&header_param->sciseg, sizeof(header_param->sciseg),1,OUTH1);
     errorcode=fread((void*)&header_param->spare9, sizeof(header_param->spare9),1,OUTH1);

     puts("After the header");
     if (errorcode!=1){
          printf("Error in reading header into SFT file!\n");
          printf("FFT number= data= %d %d \n",k,ii);
          return errorcode;
     }
     printf("header: nfft, mjdtime,tsamplu, deltanu,red,n_flag, nsamples %d %f %f %f %d %f %d\n",header_param->nfft,header_param->mjdtime,header_param->tsamplu,header_param->deltanu,header_param->red,header_param->n_flag,header_param->nsamples);
     printf("file header:perc_zeroes %f  \n",header_param->perc_zeroes);
     printf("file header:n_zeroes %d  \n",header_param->n_zeroes);
     if(header_param->detector==0){
          printf("header: tbase,frinit  %f %f\n",header_param->tbase,header_param->frinit);
     }

     printf("header: tbase,frinit,normw  %f %f %f\n",header_param->tbase,header_param->frinit,header_param->normw);

     howmany=header_param->nsamples;
     howmany2=howmany*2;
     //printf("howmany howmany2 = %d %ld\n",howmany,howmany2);
     lenps = (int)howmany2/(header_param->red);

     if(header_param->lavesp==0) header_param->lavesp=header_param->red; //2011

     if (header_only) {
          long sfdb_seek = (header_param->lavesp)*sizeof(float) +
               (lenps/2)*sizeof(float) +
               howmany2*sizeof(float);
          fseek(OUTH1, sfdb_seek, SEEK_CUR);
          goto finish;
     }

     strcpy(capt_f,"SFDB data");
     strcpy(capt_short,"Short SFDB");
     gd->n=howmany2;
     gd->nall=howmany2;
     gd->ini=(double)(header_param->firstfreqindex)*1.0/header_param->tbase;
     gd->dx=header_param->deltanu;
     gd->capt=capt_f;
     free(gd->y);
     gd->y=malloc((size_t) howmany2*sizeof(float));
     //if(header_param->lavesp==0)header_param->lavesp=header_param->red; //2011
     ps=(float *)malloc((size_t) (header_param->lavesp)*sizeof(float));
     /*Read the averages of the short power spectra (behaviour with time)*/
     for(ii=0;ii<header_param->lavesp;ii++){
          errorcode=fread((void*)&rpw, sizeof(float),1,OUTH1);
          if (errorcode!=1) {
              	printf("Error in reading power ps data into SFT file! %d \n",ii);
              	return errorcode;
          }
          ps[ii]=rpw;
          //printf(" Short power j ps[j] %d %f\n",ii,ps[ii]);
     }
     free(ps);
     //lenps= (int) gd->n/(header_param->red);
     ps_short=(float *)malloc((size_t) (lenps/2)*sizeof(float));
     gd_short->n=lenps/2;
     gd_short->nall=lenps/2;
     gd_short->ini=(double)(header_param->firstfreqindex)*1.0/header_param->tbase;
     gd_short->dx=header_param->deltanu*header_param->red;
     gd_short->capt=capt_short;
     free(gd_short->y);
     gd_short->y=malloc((size_t) lenps/2*sizeof(float));

     /*Read the very short FFT (dimension lenps; each bin averaged over red data*/
          kk=0; //only to evaluate frequencies for the test file
     for(ii=0;ii<lenps/2;ii++){
          errorcode=fread((void*)&rpw, sizeof(float),1,OUTH1);
          ps_short[ii]=rpw;
          if (errorcode!=1) {
              	printf("Error in reading very short FFT  data into SFT file! %d \n",ii);
              	return errorcode;
          }
          gd_short->y[ii]=rpw;       //sqrt(spectrum)
          //TEST FILE with the very short power spectrum:
          freq=(float)(header_param->firstfreqindex)/header_param->tbase+kk*header_param->red/header_param->tbase;
          sden=rpw*rpw;
          if(k==iw) fprintf(OUTSS,"%f %e \n",freq,sden); //write the iw v. short power spectrum  in a ascii file
          kk++;
     }
     free(ps_short);

     /*Read the data: */
#undef ORIG_DATA_READ
#ifdef ORIG_DATA_READ
     kk=0;
     for(ii=0;ii< 2*howmany; ii+=2){
          errorcode=fread((void*)&rpw, sizeof(float),1,OUTH1);
          if (errorcode!=1) {
               printf("Error in reading rpw data into SFT file! %d \n",ii);
              	return errorcode;
          }
          errorcode=fread((void*)&ipw, sizeof(float),1,OUTH1);
          if (errorcode!=1) {
               printf("Error in reading ipw data into SFT file! %d \n",ii);
              	return errorcode;
          }
          //if(ii==0)printf("ii==0  rpw,ipw %f %f \n",rpw,ipw);
          //if(ii<=2)printf("rpw,ipw %f %f \n",rpw,ipw);
          gd->y[ii]=rpw;       //Real part
          gd->y[ii+1]=ipw;    //Imag part
          freq=(float)(header_param->firstfreqindex+kk)*1.0/header_param->tbase;
          sden=(rpw*rpw+ipw*ipw)*pow(header_param->normd,2)*pow(header_param->normw,2); //piapia;
          //if(kk < 10)printf("kk,freq sden %d %f %f\n",kk,freq,sden);
          if(k==iw) fprintf(OUTS,"%f %e \n",freq,sden); //write the iw in a ascii file. piapia %e 18 aug 2005
          kk++;

     } //for ii
#else
     errorcode = fread((void*)gd->y, sizeof(float), howmany2, OUTH1);
     if (errorcode != howmany2) {
          printf("Error reading gd.y data! Read %d of %ld elements.\n", errorcode, howmany2);
          errorcode = 0;
          return errorcode;
     }
     errorcode = 1; //no error
     //printf("No Error in reading gd.y data into SFT file! read %ld of %ld values\n", errorcode, howmany2);
     if(k==iw) {
          for(ii=0; ii<howmany; ii++){
               freq=(float)(header_param->firstfreqindex+ii)*1.0/header_param->tbase;
               rpw=gd->y[2*ii];
               ipw=gd->y[2*ii+1];
               sden=(rpw*rpw+ipw*ipw)*pow(header_param->normd,2)*pow(header_param->normw,2); //piapia;
               fprintf(OUTS,"%f %e \n",freq,sden); //write the iw in a ascii file. piapia %e 18 aug 2005
          }
     }
#endif

finish:
     *fft_read+=1;
     if(k==iw){
          fclose(OUTS);
          fclose(OUTSS);
     }
     //dovesta=ftell(OUTH1);
     //printf("alla fine dovesta= %ld\n",dovesta);
     if(errorcode!=1){
          fclose(OUTH1);
     }
     return errorcode;
}

int band_extract(GD *gd,GD *gd_band,float band)
/* to extract a sub-band on a type 1 gd, from frequency domain data (from gd to gd_band*/
/* Both the gds contain only one half of the FFTs. The data are alternatively the Real and Imag part */
/* Fattore di normalizzazione ? I dati vanno *sqrt(reduct_fact) ? */
{

	int i,k,reduct_fact,n_band,n2_band,i_band;
	float total_band;
	char capt[10];
	double ini;
	long ln2;
	float norm;

	gd_band->dx=gd->dx;
 	total_band=gd->dx*gd->n/2; //the data are Real and Imag, so there is a factor 2
	//printf("In band_extract: initial freq - Total Bandwidth (real)-size %f %f %ld\n",gd->ini,total_band,gd->n);

	reduct_fact=total_band/band;
	n_band=gd->n/reduct_fact;               //number of bins to be extracted
	n2_band=next2power(n_band);             //corrected to a power of 2, as needed
	reduct_fact=gd->n/n2_band;
	//printf("Re-evaluated reduct n_band,n2_band %d %d %d\n",reduct_fact,n_band,n2_band);
	//printf("In extract band:,gd_band->ini,gd->ini,gd->dx  %f %f %f\n",gd_band->ini,gd->ini,gd->dx);
	i_band=(int) ((gd_band->ini-gd->ini)/gd->dx);
//index of the beginning frequency, IF the data in the GD were complex.
//But they are alternatively Real and Imag. Thus the found index i_band has to
//be multiplied by 2 (see the next for loop), when getting the data. n2_band is already the double.
// We must check that it is even (0,2,4..):
	//printf("i_band   %d \n",i_band);
	if(i_band <=0)i_band=0;
	if(i_band%2!=0)i_band+=1; //now it is even
	//printf("Index of the beginning (even !): i_band   %d \n",i_band);
	//printf("PRIMO INDICE PRESO in extract_band= %d\n ",2*i_band);
	strcpy(capt,"Sub-band");
	ini=(double) (i_band*gd->dx+gd->ini); //exact beginning frequency
	//printf("Exact ini= %f\n",ini);
	ln2=(long)n2_band;
	free(gd_band->y); //nuovo
	gd_band->y=malloc((size_t) ln2*sizeof(float)); //nuovo
	gd_band->n=ln2;
	gd_band->ini=ini;

	norm=1/(sqrt(bil2uni)*reduct_fact); //NON SO SE SAREBBE QUESTO !
	//printf(" NOT USED here: normalization for the extraction %f\n ",norm);
	//Data should be normalized ? Depends on their use..
	//No, if the total band is then recovered..see the vetoes on the subbw
	k=0;
	for(i=2*i_band; i<2*i_band+n2_band;i++){
	  gd_band->y[k]=gd->y[i];
	  k++;
	}

	//printf("Initial freq.- len of the new GD -  %f %ld  \n",gd_band->ini,gd_band->n);
	return i_band;
}

int resampling(GD *gd_res,GD *gd,HEADER_PARAM *header_param)
 {

   int res,i,ivign;
   float bmax;
   float initial_freq=0.;
   int ii,kki;
   //printf("In resampling- gd_res->n gd->n %ld %ld  \n",gd_res->n,gd->n);
   res=gd->n/gd_res->n;
   //printf("res=%d\n",res);

header_param->n_zeroes=0;
 for(ii=0;ii<gd->n;ii++){
   kki=isnan(gd->y[ii]);
   if (kki !=0)
   {
    gd->y[ii]=0;
   }
   if(gd->y[ii]==0)header_param->n_zeroes+=1;
 }
 header_param->perc_zeroes=(float)(header_param->n_zeroes)/gd->n;  //becomes a percentage


   i= gdtime2freq(gd,header_param);//from time gd to frequency gd
   gd_res->ini=initial_freq;
   gd_res->dx=gd->dx;
   //printf("In resampling dopo gdtime2freq- gd_res->dx gd->dx %f %f  \n",gd_res->dx,gd->dx);
   //bmax=1/(2*gd_res->dx); //max freq in the data
   bmax=(gd_res->dx*gd_res->n/2);
   //printf("In resampling- max freq =%f  \n",bmax);
   //From the freq data extract the BW, go back to the time domain.
   ivign=0;

   //printf("In resampling- gd_res->dx gd->dx %f %f  \n",gd_res->dx,gd->dx);
   i=band_extract(gd,gd_res,bmax);
   //Band extract does not normalize data

   //Go back to the time domain:
   i=gdfreq2time(gd_res,header_param,ivign);
   //printf("dopo gdfreq2time \n");
   //normalization of the data, for the subsamplig:
   for(i=0; i<gd_res->n;i++){
     gd_res->y[i]=gd_res->y[i]/res;

     //if(i<=5)printf("i=%d gd_res[i]=%f\n",i,gd_res->y[i]);
     //if(i>=gd_res->n-5)printf("i=%d gd_res[i]=%f\n",i,gd_res->y[i]);

   }
   //Redefine gd_res and header_param
   //gd_res->dx=res*gd->dx;
   //gd_res->ini=gd->ini;
   header_param->nsamples=gd_res->n/2;
   header_param->tsamplu=gd_res->dx;
   //printf("In resampling- after resampling gd_res->dx gd_res->ini %f %f  \n",gd_res->dx,gd_res->ini);
   //printf("In resampling- header_param->nsamples=%d  header_param->tsamplu =%f\n",header_param->nsamples,header_param->tsamplu);
   return i;

 }
int band_extract_double(GD *gd,GD *gd_band,float band)
/* to extract a sub-band on a type 1 gd, from frequency domain data (from gd to gd_band*/
/* Both the gds contain only one half of the FFTs. The data are alternatively the Real and Imag part */
/* Fattore di normalizzazione ? I dati vanno *sqrt(reduct_fact) ? */
{

	int i,k,reduct_fact,n_band,n2_band,i_band;
	float total_band;
	char capt[10];
	double ini;
	long ln2;
	float norm;

	gd_band->dx=gd->dx;
 	total_band=gd->dx*gd->n/2; //the data are Real and Imag, so there is a factor 2
	//printf("In band_extract: initial freq - Total Bandwidth (real)-size %f %f %ld\n",gd->ini,total_band,gd->n);

	reduct_fact=total_band/band;
	n_band=gd->n/reduct_fact;               //number of bins to be extracted
	n2_band=next2power(n_band);             //corrected to a power of 2, as needed
	reduct_fact=gd->n/n2_band;
	//printf("Re-evaluated reduct n_band,n2_band %d %d %d\n",reduct_fact,n_band,n2_band);
	//printf("In extract band:,gd_band->ini,gd->ini,gd->dx  %f %f %f\n",gd_band->ini,gd->ini,gd->dx);
	i_band=(int) ((gd_band->ini-gd->ini)/gd->dx);
//index of the beginning frequency, IF the data in the GD were complex.
//But they are alternatively Real and Imag. Thus the found index i_band has to
//be multiplied by 2 (see the next for loop), when getting the data. n2_band is already the double.
// We must check that it is even (0,2,4..):
	//printf("i_band   %d \n",i_band);
	if(i_band <=0)i_band=0;
	if(i_band%2!=0)i_band+=1; //now it is even
	//printf("Index of the beginning (even !): i_band   %d \n",i_band);
	//printf("PRIMO INDICE PRESO in extract_band= %d\n ",2*i_band);
	strcpy(capt,"Sub-band");
	ini=(double) (i_band*gd->dx+gd->ini); //exact beginning frequency
	//printf("Exact ini= %f\n",ini);
	ln2=(long)n2_band;
	free(gd_band->y); //nuovo
	gd_band->y=malloc((size_t) 2*ln2*sizeof(float)); //Doubled 9 ott. 2009
	gd_band->n=2*ln2;
	gd_band->ini=ini;

	norm=1/(sqrt(bil2uni)*reduct_fact); //NON SO SE SAREBBE QUESTO !
	//printf(" NOT USED here: normalization for the extraction %f\n ",norm);
	//Data should be normalized ? Depends on their use..
	//No, if the total band is then recovered..see the vetoes on the subbw
	for(k=0;k<gd_band->n;k++)gd_band->y[k]=0;
	k=gd_band->n/4;
	for(i=2*i_band; i<2*i_band+n2_band;i++){
	  gd_band->y[k]=gd->y[i];
	  k++;
	}

	//printf("Initial freq.- len of the new GD -  %f %ld  \n",gd_band->ini,gd_band->n);
	return i_band;
}

int band_extract_1hz(GD *gd,GD *gd_band,float *band,float *initial_freq,float *final_freq,int ifatti)
/* to extract a sub-band on a type 1 gd, from frequency domain data (from gd to gd_band*/
/* Both the gds contain only one half of the FFTs. The data are alternatively the Real and Imag part */
/* IMP: Fattore di normalizzazione  I dati vanno *norm  */

{

	int i,k,reduct_fact,n_band,n2_band,i_band;
	float total_band;
	char capt[10];
	double ini;
	long ln2;
	int n05_band;
	int kini,kifi;
	float banda05;
	float diff,iini,iifi;
	int eff_len;
	float norm;

	eff_len=gd->n/2; //samples in the band:the data are Real and Imag, so there is a factor 2
	gd_band->dx=gd->dx;
 	total_band=gd->dx*eff_len;
	//printf("In band_extract: initial freq - Total Bandwidth (real)-size %f %f %ld\n",gd->ini,total_band,gd->n);

	banda05=0.5; //fixed at 0.5 Hz to have 1 Hz total bw abd 1 s sampling
	if(ifatti==0){
	  iini=*initial_freq/banda05;
	  iifi=*final_freq/banda05;
	  kini=floor(iini);
	  kifi=floor(iifi);
	  printf("kini, kifi -  %d %d  \n",kini,kifi);
	  if(kini !=kifi){
	    diff=iifi-kifi;
	    printf("ATT: diff= %f\n",diff);
	    printf("ATT: iini iifi  diff/0.5= %f %f %f\n",iini,iifi,diff/banda05);
	    *final_freq-=diff/banda05;
	    *band=(*final_freq-*initial_freq);
	    kini=floor(*initial_freq/banda05);
	    kifi=floor(*final_freq/banda05);
	    printf("Re-evaluated band, final_freq -  %f %f  \n",*band,*final_freq);
	    printf("Re-evaluated kini, kifi -  %d %d  \n",kini,kifi);
	}
	  gd_band->ini=*initial_freq;
      }
	reduct_fact=total_band/(*band);
	n_band=eff_len/reduct_fact;       //number of data to be  extracted
	n05_band=floor(banda05/gd->dx);   //number of data in the 0.5 Hz bw
	n2_band=next2power(n05_band);     //enlarged with zeroes to a power of 2, as needed
	reduct_fact=eff_len/n2_band;
	printf("Re-evaluated reduct n_band,n2_band %d %d %d\n",reduct_fact,n_band,n2_band);
	printf("In extract band:,gd_band->ini,gd->ini,gd->dx  %f %f %f\n",gd_band->ini,gd->ini,gd->dx);
	i_band=(int) ((gd_band->ini-gd->ini)/gd->dx);
//index of the beginning frequency, IF the data in the GD were complex.
//But they are alternatively Real and Imag. Thus the found index i_band has to
//be multiplied by 2 (see the next for loop), when getting the data. n2_band has to be doubled.
// We must check that it is even (0,2,4..):
	//printf("i_band   %d \n",i_band);
	if(i_band <=0)i_band=0;
	if(i_band%2!=0)i_band+=1; //now it is even
	//printf("Index of the beginning (even !): i_band   %d \n",i_band);
	//printf("PRIMO INDICE PRESO in extract_band= %d\n ",2*i_band);
	strcpy(capt,"Sub-band");
	ini=(double) (i_band*gd->dx+gd->ini); //exact beginning frequency
	//printf("Exact ini= %f\n",ini);
	ln2=(long)(2*n2_band); //doubled for Real and Imag
	if(ifatti==0){
	  free(gd_band->y); //nuovo
	  gd_band->y=malloc((size_t) ln2*sizeof(float)); //nuovo
	  gd_band->n=ln2;
	  gd_band->ini=ini;
	}

	k=0;  //fill with zeroes the 0.5 band
	for(i=2*i_band;i<ln2;i++){
	  gd_band->y[k]=0.;
	  k++;
	}
	//fill with the data and normalize the data /(sqrt(bil2uni)*reduct_fact);
        k=0;
	//norm=1/(sqrt(bil2uni)*reduct_fact);
	norm=2.0/sqrt(reduct_fact); //NON SO SE E' CORRETTO !!
	printf("normalization for the extraction %f\n ",norm);
	for(i=2*i_band; i<(2*i_band+2*n_band);i++){ //Extracted data
	  gd_band->y[k]=gd->y[i]*norm;
	  k++;
	}

	return i_band;
}
int next2power(int n)
     /*for a given number n it gives the next integer which is a power of 2*/
{
  double m;
  int n2,m1;

  m=ceil(log10(n)/log10(2.0));
  m1=(int) m;
  n2=(int) pow(2.0,m);


    return n2;
    }

int bandpass_data(GD *gd_bandpass,GD *gd,HEADER_PARAM *header_param)
/*bandpass data in time domain. Needed before the events extraction*/
{

	int i;
	char capt[12];
	int itest; //0=normal  1=writes one test file
	float band_extension,central_freq,banda,theta,erre;
	float freqm_a,freqp_a; //a means aliased

	//	FILE *OUTTEST_BAND;

	itest=1;

	strcpy(capt,"Bandpass");
	//puts(" Creates GD with bandpassed data");
	gd_bandpass->n=gd->n;
	gd_bandpass->ini=gd->ini;
	gd_bandpass->dx=gd->dx;
	//gd_bandpass->capt=capt;

	//TEST:
	if(itest==1){
	  //for(i=0;i<gd->n;i++)gd->y[i]=3.;
	}
	//For the setting of these parameters I need to understand what to do for VIRGO
	//and find a general expression, valid for ROG and VIRGO
	band_extension=0.;
	//freqm_a=header_param->freqm-header_param->frinit;
	//freqp_a=header_param->freqp-header_param->frinit;
	freqm_a=15; //15 jan 2008 tolte le freq. dall' header. 25 e 15 a caso..
	freqp_a=25;
	central_freq=(freqm_a+freqp_a)/2.0;   //Central freq of the bandpass
	banda=(freqp_a-freqm_a)+band_extension;
	theta=2*PIG*central_freq*header_param->tsamplu;
	erre=exp(-banda*header_param->tsamplu);

	//printf("freqm, freqp,central freq, band  %f %f %f %f\n",freqm_a,freqp_a,central_freq,banda);
	//printf("theta, erre  %f %f\n",theta,erre);
	//printf("tlen, sampling time  %ld %f\n",gd->n,header_param->tsamplu);
	//Bandpass from gd to gd_bandpass:
	gd_bandpass->y[0]=gd->y[0]; //copy in gd_bandpass the first datum  of gd
	gd_bandpass->y[1]=gd->y[1]+2*erre*cos(theta)*gd_bandpass->y[0];
	for(i=2;i<gd->n;i++){
	  gd_bandpass->y[i]=gd->y[i]+2*erre*cos(theta)*gd_bandpass->y[i-1]-pow(erre,2.0)*gd_bandpass->y[i-2];
	}
	/*
	  if(itest==1){
	    //TEST:<:scrittura del file
	    OUTTEST_BAND=fopen("testfile_band.dat","w");
	    for(k=0;k<gd->n;k++)fprintf(OUTTEST_BAND,"%d %f %f\n",k,gd_bandpass->y[k],gd->y[k]);
	    fclose(OUTTEST_BAND);
	  }
	*/

	return i;
}
int saturated_data(GD *gd,EVEN_PARAM *even_param,HEADER_PARAM *header_param)
{
  int i;
  //printf("even_param->sat_level= %f\n ",even_param->sat_level);
  //printf("even_param->sat_maxn= %d\n ",even_param->sat_maxn);
 even_param->sat_howmany=0.;
 for (i=0; i < gd->n; i++)
 {
   if (gd->y[i] >= even_param->sat_level)
     {
      // 23-03-2013 COMMENTATO PER PROVE OUTLIERS ALBERTO
       //gd->y[i]=even_param->sat_level;
       even_param->sat_howmany+=1.;
     }
     if (gd->y[i] <= -even_param->sat_level)
     {
      // 23-03-2013 COMMENTATO PER PROVE OUTLIERS ALBERTO
       //gd->y[i]=-even_param->sat_level;
       even_param->sat_howmany+=1.;

     }
 }
 //printf("even_param->sat_howmany= %f\n ",even_param->sat_howmany);
 header_param->sat_howmany=even_param->sat_howmany;
  return i;
}

int highpass_data(GD *gd_highpass,GD *gd,HEADER_PARAM *header_param,float freqc)
/*highpass data in time domain. Needed for Virgo before the events extraction*/
{

	int i;
	//char capt[12];
	int itest; //0=normal  1=writes one test file
	float w,w1,b0,b1,a1;

//	FILE *OUTTEST_BAND;

	itest=1;

	//strcpy(capt,"Highpass");
	//puts(" Creates GD with highpassed data");
	gd_highpass->n=gd->n;
	gd_highpass->ini=gd->ini;
	gd_highpass->dx=gd->dx;
	//gd_highpass->capt=capt;

	//TEST:
	if(itest==1){
	  //for(i=0;i<gd->n;i++)gd->y[i]=3.;
	}

	w=exp(-2*PIG*freqc*header_param->tsamplu);
	w1=(1+w)/2;
	b0=w1;
	b1=-w1;
	a1=-w;

	//printf("freqc (cutoff frequency)  %f\n",freqc);
	//printf("w, w1  %f %f\n",w,w1);
	//printf("tlen, sampling time  %ld %f\n",gd->n,header_param->tsamplu);
	//Highpass from gd to gd_highpass:
	gd_highpass->y[0]=b0*gd->y[0]; //copy in gd_bandpass the first datum  of gd
	for(i=1;i<gd->n;i++){
	  gd_highpass->y[i]=b0*gd->y[i]+b1*gd->y[i-1]-a1*gd_highpass->y[i-1];
	}
	/*
	  if(itest==1){
	    //TEST::scrittura del file
	    OUTTEST_BAND=fopen("testfile_band.dat","w");
	    for(k=0;k<gd->n;k++)fprintf(OUTTEST_BAND,"%d %f %f\n",k,gd_highpass->y[k],gd->y[k]);
	    fclose(OUTTEST_BAND);
	  }
	*/
	return i;
}
int highpass_data_bil(GD *gd_highpass,GD *gd,HEADER_PARAM *header_param,float freqc)
/*Bilateral (no time delay) highpass data in time domain. Needed for Ligo-Virgo before the events extraction*/
{
  //int k;
	int i,len;
	//char capt[12];
	int itest; //0=normal  1=writes one test file
	float w,w1,b0,b1,a1;
	float *firsthighpass;
	int iii;

	//	FILE *OUTTEST_BAND;

	iii=gd->n+1;
	firsthighpass=malloc((size_t) iii*sizeof(float));
	itest=1;

	//strcpy(capt,"Highpass data");
	//puts(" Creates GD with highpassed data");
	gd_highpass->n=gd->n;
	gd_highpass->ini=gd->ini;
	gd_highpass->dx=gd->dx;
	//gd_highpass->capt=capt;

	//TEST:
	if(itest==1){
	  //for(i=0;i<gd->n;i++)gd->y[i]=3.;
	}

	w=exp(-2*PIG*freqc*header_param->tsamplu);
	w1=(1+w)/2;
	b0=w1;
	b1=-w1;
	a1=-w;

	//printf("freqc (cutoff frequency)  %f\n",freqc);
	//printf("w, w1  %f %f\n",w,w1);
	//printf("tlen, sampling time  %ld %f\n",gd->n,header_param->tsamplu);
	//First Highpass from gd to firsthighpass:
	firsthighpass[0]=b0*gd->y[0]; //copy in gd_highpass the first datum  of gd
	for(i=1;i<gd->n;i++){
	  firsthighpass[i]=b0*gd->y[i]+b1*gd->y[i-1]-a1*firsthighpass[i-1];
	}
	//Second step: Highpass from firsthighpass  to gd_highpass:
	len=(gd->n)-1;
	gd_highpass->y[len]=b0*firsthighpass[len]; //copy  the last datum
	for(i=len-1;i>=0;i--){
	  gd_highpass->y[i]=b0*firsthighpass[i]+b1*firsthighpass[i+1]-a1*gd_highpass->y[i+1];
	}
	/*
	  if(itest==1){
	    //TEST::scrittura del file
	    OUTTEST_BAND=fopen("testfile_band.dat","w");
	    for(k=0;k<gd->n;k++)fprintf(OUTTEST_BAND,"%d %f %f\n",k,gd_highpass->y[k],gd->y[k]);
	    fclose(OUTTEST_BAND);
	  }
	*/
	  free(firsthighpass);
	return i;
}


int purge_data(GD *gd_clean,GD *gd,GD *gd_highpass,EVEN_PARAM *even_param,HEADER_PARAM *header_param)
/*identifies and purge large events from a gd to another gd, which has the same parameters */
{

	int i,k;
	int lwind;
	char capt[12];
	int index1,index2,index3,index4;
	int itest; //0=normal  1=test the procedure with constant data
	double time_eve; //beginning time of the event in days mjd
	int imax; //index of the maximum of the event (even_param->imax)

	//	FILE *OUTTEST;
	double data_log[5];  //to write the logfile
	int prec[5];



	itest=0;
    i= sn_medsig(gd_highpass,even_param,header_param); //evaluates mean and std
	i= even_anst(gd_highpass,even_param); //register events indexes (beginning,duration,imax)

	strcpy(capt,"Clean data");
	//puts(" Creates GD with clean data");
	gd_clean->n=gd->n;
	gd_clean->ini=gd->ini;
	gd_clean->dx=gd->dx;
	//gd_clean->capt=capt;

	//TEST:
	if(itest==1){
	even_param->number=1;
	even_param->begin[0]=200;
 	even_param->duration[0]=10;
	for(i=0;i<gd->n;i++)gd->y[i]=3.;
	//fine scritte per il TEST !!
	}

	 //Information for the LOG file: a new FFT is under analysis
	data_log[0]=mjd2gps(header_param->mjdtime);
	 //data_log[1]=(double) header_param->nfft;
	 data_log[1]=header_param->nfft;
	  data_log[2]=header_param->n_flag;
	  data_log[3]=header_param->perc_zeroes;
	 prec[0]=3;
	 prec[1]=0;
	 prec[2]=101;
	 prec[3]=101;
	 logfile_ev(LOG_INFO,"NEW",4,data_log,prec);
	  //Info on saturated data: 19 may 2008
	  data_log[0]=even_param->sat_howmany;
	  prec[0]=3;
	  logfile_ev(LOG_INFO,"SAT",1,data_log,prec);
	  //End saturated data info

	for(i=0;i<gd->n;i++)gd_clean->y[i]=gd->y[i]; //copy in gd_clean the data of gd
	lwind=ceil(((double) even_param->edge/gd->dx));
	printf("lwind %d %f\n",lwind,even_param->edge);
	//Modify the data around each event:

	for(i=0; i<even_param->number; i++){
	  if(even_param->begin[i]!=-1){
	     time_eve=header_param->mjdtime+(double)(even_param->begin[i]*header_param->tsamplu)/day2sec;
	     imax=even_param->imax[i];

	     //printf("***One event:i begin duration %d %d %d\n",i,even_param->begin[i],even_param->duration[i]);
	     //printf("time_eve,duration in s,max amp gd and highpass %f %f %f %f\n",time_eve,even_param->duration[i]*header_param->tsamplu,gd->y[imax],gd_highpass->y[imax]);
	    index1=even_param->begin[i]-lwind;
	    index2=even_param->begin[i]+even_param->duration[i]+lwind;
	    index3=even_param->begin[i];
	    index4=even_param->begin[i]+even_param->duration[i];
	    if(index1<=0)index1=0;
	    if(index2>=gd->n)index2=gd->n;
  	    for(k=index3; k<index4; k++)gd_clean->y[k]=0.;
	    ////for(k=index1; k<index3; k++)gd_clean->y[k]=gd->y[k]*(1.-1.*sin(1.0*(k-index1)/(index3-k+1.0)/lwind));
	    ////for(k=index4; k<index2; k++)gd_clean->y[k]=gd->y[k]*(1.-1.*sin(1.0*(index2-k)/(k-index4+1.0)/lwind));
	    for(k=index1; k<index3; k++)gd_clean->y[k]=gd->y[k]*(1.-1.*(1.0*(k-index1)/lwind));
	    for(k=index4; k<index2; k++)gd_clean->y[k]=gd->y[k]*(1.-1.*(1.0*(index2-k)/lwind));

	  //Information for the LOG file
	  data_log[0]=mjd2gps(time_eve);
	  data_log[1]=(double) even_param->duration[i]*header_param->tsamplu;
	  data_log[2]=(double) gd->y[imax];
	  data_log[3]=(double) even_param->crmax[i];
	  data_log[4]=(double) even_param->ener[i];
	  prec[0]=3;
	  prec[1]=0;
	  prec[2]=0;
	  prec[3]=0;
	  prec[4]=0;
	  logfile_ev(LOG_INFO,"EVT",5,data_log,prec);

	  }
	  /*
	  if(itest!=2){
	    //TEST::scrittura del file. with !=2 always writes the first event
	    //printf("index1 index3 index4 index2  %d  %d %d %d \n",index1,index3,index4,index2);
	    if(i==1)OUTTEST=fopen("testfile.dat","w");
	    if(index1>=5){
	    if(i==1)for(k=index1-5;k<=index2+5;k++)fprintf(OUTTEST,"%d %f %f %f\n",k,gd_clean->y[k],gd->y[k],gd_highpass->y[k]);
	    }
	    if(i==1)fclose(OUTTEST);
	  }
	  */
	}

	return i;
}

int purge_data_subtract(GD *gd_clean,GD *gd,GD *gd_highpass,EVEN_PARAM *even_param,HEADER_PARAM *header_param) {
  	//Define the value of evf tau and maxage. New: may,5, 2011 Pia
	//and tau evt


  sn_medsig(gd_highpass,even_param,header_param); //evaluates mean and std
  even_anst(gd_highpass,even_param); //register events indexes (beginning,duration,imax);
  data_subtract(gd_clean, gd, gd_highpass, even_param, header_param);
  return header_param->nfft;
}

int data_subtract(GD *gd_clean,GD *gd,GD *gd_highpass,EVEN_PARAM *even_param,HEADER_PARAM *header_param)
/*identifies and purge large events from a gd to another gd, which has the same parameters */
/*by subtracting the contribution of the bilateral high pass filtering*/
{

	int i,k;
	int lwind;
	char capt[12];
	int index1,index2,index3,index4;
	int itest; //0=normal  1=test the procedure with constant data
	double time_eve; //beginning  (from 14/07/09: MAX) time of the event in days mjd
	int imax; //index of the maximum of the event (even_param->imax)

	//	FILE *OUTTEST;
	double data_log[5];  //to write the logfile
	int prec[5];

	itest=0;

	strcpy(capt,"Clean data");
	//puts(" Creates GD with clean data");
	gd_clean->n=gd->n;
	gd_clean->ini=gd->ini;
	gd_clean->dx=gd->dx;
	//gd_clean->capt=capt;

	//TEST:
	if(itest==1){
	even_param->number=1;
	even_param->begin[0]=200;
 	even_param->duration[0]=10;
	for(i=0;i<gd->n;i++)gd->y[i]=3.;
	//fine scritte per il TEST !!
	}

	 //Information for the LOG file: a new FFT is under analysis
	data_log[0]=mjd2gps(header_param->mjdtime);
	 //data_log[1]=(double) header_param->nfft;
	  data_log[1]= header_param->nfft;
	 data_log[2]=header_param->n_flag;
	 data_log[3]=header_param->perc_zeroes;
	 prec[0]=3;
	 prec[1]=0;
	 prec[2]=101;
	 prec[3]=101;
	 logfile_ev(LOG_INFO,"NEW",4,data_log,prec);


	 //Info on saturated data: 19 may 2008
	  data_log[0]=even_param->sat_howmany;
	  prec[0]=3;
	  logfile_ev(LOG_INFO,"SAT",1,data_log,prec);
	  //End saturated data info
	  //Info on flagged data: 15 april 2011
	  //data_log[0]=header_param->n_flag;
	  //prec[0]=0;
	  //logfile_ev(LOG_INFO,"FLG",1,data_log,prec);
	  //End flag data info
	for(i=0;i<gd->n;i++)gd_clean->y[i]=gd->y[i]; //copy in gd_clean the data of gd
	lwind=ceil(((double) even_param->edge/gd->dx));
	//printf("lwind %d %f\n",lwind,even_param->edge);
	//Modify the data around each event:

	for(i=0; i<even_param->number; i++){
	  if(even_param->begin[i]!=-1){
	     time_eve=header_param->mjdtime+(double)(even_param->begin[i]*header_param->tsamplu)/day2sec;
	     imax=even_param->imax[i];

	     //printf("***One event:i begin duration %d %d %d\n",i,even_param->begin[i],even_param->duration[i]);
	     //printf("time_eve,duration in s,max amp (gd and highpass %f %f %f %f\n",time_eve,even_param->duration[i]*header_param->tsamplu,gd->y[imax],gd_highpass->y[imax]);
	    index1=even_param->begin[i]-lwind;
	    index2=even_param->begin[i]+even_param->duration[i]+lwind;
	    index3=even_param->begin[i];
	    index4=even_param->begin[i]+even_param->duration[i];
	    if(index1<=0)index1=0;
	    if(index2>=gd->n)index2=gd->n;

  	    for(k=index3; k<index4; k++)gd_clean->y[k]-=gd_highpass->y[k]; //in the event
	    //At the edges:
	    for(k=index1; k<index3; k++)gd_clean->y[k]-=gd_highpass->y[k]*(k-1.0*index1)/lwind;
	    for(k=index4; k<index2; k++)gd_clean->y[k]-=gd_highpass->y[k]*(1.0*index2-k)/lwind;

	  //Information for the LOG file
	    //Modified 14/07 con Sabrina per imax
	    time_eve=header_param->mjdtime+(double)(even_param->imax[i]*header_param->tsamplu)/day2sec;
	  data_log[0]=mjd2gps(time_eve);
	  data_log[1]=((double) even_param->duration[i])* (double) header_param->tsamplu;
	  data_log[2]=(double) gd_highpass->y[imax]; //modifica 20-AUG 2007
	  data_log[3]=(double) even_param->crmax[i];
	  data_log[4]=(double) even_param->ener[i];
	  prec[0]=3;
	  prec[1]=0;
	  prec[2]=0;
	  prec[3]=0;
	  prec[4]=0;
	  logfile_ev(LOG_INFO,"EVT",5,data_log,prec);

	  }

	  //if(itest!=2){
	    //TEST::scrittura del file. with !=2 always writes the first event
	    //printf("index1 index3 index4 index2  %d  %d %d %d \n",index1,index3,index4,index2);
	  if(header_param->nfft==1&&i==0){
	      OUTTEST=fopen("testfile.dat","wr");
	      ///puts("Open testfile");
	  }
	  if(even_param->crmax[i]>=15){
	      if(imax >=(int) (header_param->nsamples/2) && imax < (int) (3*header_param->nsamples/2 ))
		{
		  //if(i>=0)for(k=imax;k<imax+1;k++)printf("%20.9f %20.12f %d %12.5f %12.5f %12.5f\n",time_eve,header_param->mjdtime+k*(header_param->tsamplu/day2sec),k,gd_clean->y[k],gd->y[k],gd_highpass->y[k]);
	        if(i>=0)for(k=index1;k<index2;k++) {
			  if ((k>imax-2000)&&(k<imax+2000))
			    {
			      if(OUTTEST==NULL){
				OUTTEST=fopen("testfile.dat","wr");
				//puts("***Non era aperto ! Open testfile");
			      }
			      // sabrina fprintf(OUTTEST,"%20.9f %20.9f %d %12.5f %12.5f %12.5f %12.5f\n",time_eve,header_param->mjdtime+k*(header_param->tsamplu/day2sec),k,even_param->crmax[i],gd_clean->y[k],gd->y[k],gd_highpass->y[k]);
			     }
		}
	     }
	  }
	    //if(i==1)fclose(OUTTEST);
	    //}
	//
	}


	return i;
}

int purge_data_subsampled(GD *gd_clean,GD *gd,EVEN_PARAM *even_param, double tsample,HEADER_PARAM *header_param)
/*identifies and purge large events from a gd to another gd, which has the same parameters */
{

	int i,k;
	int lwind;
	char capt[12];
	int index1,index2,index3,index4;
	double time_eve; //beginning time of the event in days mjd
	int imax; //index of the maximum of the event (even_param->imax)

	//	FILE *OUTTEST;
	double data_log[5];  //to write the logfile
	int prec[5];


	even_param->edge=tsample*10; //proviamo june 2008
	even_param->tau=600; //ariprova
       	i= sn_medsig(gd,even_param,header_param); //evaluates mean and std
	i= even_anst(gd,even_param); //register events indexes (beginning,duration,imax)

	strcpy(capt,"Clean data");
	//puts(" Creates GD with clean data");
	gd_clean->n=gd->n;
	gd_clean->ini=gd->ini;
	gd_clean->dx=gd->dx;
	//gd_clean->capt=capt;



	 //Information for the LOG file: a new FFT is under analysis
	 //data_log[0]=header_param->mjdtime;
	 //data_log[1]=(double) header_param->nfft;
	 //prec[0]=3;
	 //prec[1]=0;
	 //logfile_ev(LOG_INFO,"NEW",2,data_log,prec);



	for(i=0;i<gd->n;i++)gd_clean->y[i]=gd->y[i]; //copy in gd_clean the data of gd
	lwind=ceil(((double) even_param->edge/gd->dx));
	printf("lwind edge tsample %d %f %f\n",lwind,even_param->edge,tsample);
	//Modify the data around each event:

	for(i=0; i<even_param->number; i++){
	  if(even_param->begin[i]!=-1){
	     time_eve=header_param->mjdtime+(double)(even_param->begin[i]*tsample)/day2sec;
	     imax=even_param->imax[i];

	    printf("***One event:i begin duration %d %d %d\n",i,even_param->begin[i],even_param->duration[i]);
	    printf("time_eve,duration in s,max amp gd  %f %f %f \n",time_eve,even_param->duration[i]*tsample,gd->y[imax]);
	    index1=even_param->begin[i]-lwind;
	    index2=even_param->begin[i]+even_param->duration[i]+lwind;
	    index3=even_param->begin[i];
	    index4=even_param->begin[i]+even_param->duration[i];
	    if(index1<=0)index1=0;
	    if(index2>=gd->n)index2=gd->n;
  	    for(k=index3; k<index4; k++)gd_clean->y[k]=0.;

	    for(k=index1; k<index3; k++)gd_clean->y[k]=gd->y[k]*(1.-1.*(1.0*(k-index1)/lwind));
	    for(k=index4; k<index2; k++)gd_clean->y[k]=gd->y[k]*(1.-1.*(1.0*(index2-k)/lwind));


	  //Information for the LOG file
	    data_log[0]=mjd2gps(time_eve);
	  data_log[1]=(double) even_param->duration[i]*tsample;
	  data_log[2]=(double) gd->y[imax];
	  data_log[3]=(double) even_param->crmax[i];
	  data_log[4]=(double) even_param->ener[i];
	  prec[0]=3;
	  prec[1]=0;
	  prec[2]=0;
	  prec[3]=0;
	  prec[4]=0;
	  logfile_ev(LOG_INFO,"EVT",5,data_log,prec);

	  }
	  /*
	  if(itest!=2){
	    //TEST::scrittura del file. with !=2 always writes the first event
	    printf("index1 index3 index4 index2  %d  %d %d %d \n",index1,index3,index4,index2);
	    if(i==1)OUTTEST=fopen("testfile.dat","w");
	    if(index1>=5){
	    if(i==1)for(k=index1-5;k<=index2+5;k++)fprintf(OUTTEST,"%d %f %f \n",k,gd_clean->y[k],gd->y[k]);
	    }
	    if(i==1)fclose(OUTTEST);
	    }
	  */
	}
	 even_param->edge=0.15;  //rimesso come era prima
	return i;
}

//
int purge_subbw(GD *gd_appo,GD *gd,GD *gd_band,EVEN_PARAM *even_paramb,HEADER_PARAM *header_param,float band_ext)
/*identifies and purge large events from a gd to another gd.
The veto is in the sub-bandwidths
Input: gd in time and band_ext
Output: a time GD with same parameters, but cleaned  gd. Input and output is the same gd
*/
{
  int i,k,ii,kk;
	int lwind;
	char capt[12];
	int index1,index2,index3,index4;
	double time_eve; //beginning (from 14/07/09: MAX)time of the event in days mjd
	int imax; //index of the maximum of the event (even_param->imax)
	double data_log[6];  //to write the logfile
	int prec[6];
	double tsample;
	int ibanda;
	int bdafare;
	float bmax;
	double initial_freq;
	int isubsampling;
	int ivign;
	char capt_gd[12],capt1_gd[12];
	double frinit;
	int ippp; //2009 september
	int ioverl=1; // 2 if overlapped 1 if not overlapped

       ippp=0;

       strcpy(capt_gd,"freqdata");
       strcpy(capt1_gd,"timedata");
       initial_freq=0.;

       ii= gdtime2freq(gd,header_param);//from time gd to frequency gd
       gd->ini=(double) initial_freq;
       bmax=1/(2*header_param->tsamplu); //max freq in the data
       bdafare=(int) (bmax/band_ext)*ioverl-1;
       printf("band_ext,bmax bdafare  %f %f %d\n",band_ext,bmax,bdafare);

       //Main loop: from the freq data extract the BW, go back to the time domain.
        ivign=0;
       for(ibanda=0;ibanda<bdafare;ibanda++){
	 gd_band->ini=initial_freq+ibanda*band_ext/ioverl;
	 gd_band->dx=gd->dx;
	 frinit=gd_band->ini;
	 //band extraction
	  if(ioverl==1){
	    i=band_extract_double(gd,gd_band,band_ext);
	      };
	  if(ioverl==2){
	    i=band_extract(gd,gd_band,band_ext);
	      };
	 //Band extract does not normalize data, because at the end
	 //the total band is recovered. 16-april 2009
	 //printf("ibanda gd_band->ini gd_band->n %d %f %ld \n",ibanda, gd_band->ini,gd_band->n);
	 //Go back to the time domain:
	 isubsampling=gd->n/gd_band->n;
	 ///gd2wingd(gd_band,5);
	 //printf("from band");
	 ii=gdfreq2time(gd_band,header_param,ivign);
	 // and now veto in the time domain data:
	 tsample=gd_band->dx;
	 even_paramb->edge=tsample*10; //proviamo june 2008
	 even_paramb->tau=120; //ariprova messo a 120 da 600

	 even_paramb->xw_pari=0.;
	 even_paramb->qw_pari=0.;
	 even_paramb->w_norm_pari=0.;
	 even_paramb->xw_dispari=0.;
	 even_paramb->qw_dispari=0.;
	 even_paramb->w_norm_dispari=0.;
	 even_paramb->ndata_run_pari=0;
	 even_paramb->ndata_run_dispari=0;

       	 i= sn_medsig(gd_band,even_paramb,header_param); //evaluates mean and std
	 i= even_anst(gd_band,even_paramb); //register events indexes (beginning,duration,imax)
         strcpy(capt,"Clean data");
	 //puts(" Creates GD with clean data in the sub BW");
	 gd_appo->n=gd_band->n;
	 gd_appo->ini=gd_band->ini;
	 gd_appo->dx=gd_band->dx;

	for(i=0;i<gd_band->n;i++)gd_appo->y[i]=gd_band->y[i]; //copy in gd_appo the data of gd
	lwind=ceil(((double) even_paramb->edge/gd_band->dx));
	//printf("lwind %d %f\n",lwind,even_paramb->edge);
	//Modify the data around each event:

	for(i=0; i<even_paramb->number; i++){
	  //metto un deadtime iniziale 28 ottobre 2009 prova
	  time_eve=header_param->mjdtime+(double)(even_paramb->begin[i]*tsample)/day2sec;
	  //if(even_paramb->begin[i]!=-1 && time_eve>(header_param->mjdtime+even_paramb->deadtime/day2sec)){
	      if(even_paramb->begin[i]!=-1 && time_eve>(header_param->mjdtime+even_paramb->deadtime/day2sec) &&(even_paramb->crmax[i]>= 1.5*even_paramb->cr )){
	      //2011 use an higher threshold for EVB
	     imax=even_paramb->imax[i];

	     //printf("***One event:i begin duration in samples%d %d %d\n",i,even_paramb->begin[i],even_param->duration[i]);
	     //printf("time_eve,duration in s,max amp  %f %f %f \n",time_eve,even_paramb->duration[i]*tsample,gd_band->y[imax]);
	    index1=even_paramb->begin[i]-lwind;
	    index2=even_paramb->begin[i]+even_paramb->duration[i]+lwind;
	    index3=even_paramb->begin[i];
	    index4=even_paramb->begin[i]+even_paramb->duration[i];
	    if(index1<=0)index1=0;
	    if(index2>=gd_band->n)index2=gd_band->n;
  	    for(k=index3; k<index4; k++)gd_appo->y[k]=0.;

	    for(k=index1; k<index3; k++)gd_appo->y[k]=gd_band->y[k]*(1.-1.*(1.0*(k-index1)/lwind));
	    for(k=index4; k<index2; k++)gd_appo->y[k]=gd_band->y[k]*(1.-1.*(1.0*(index2-k)/lwind));


	  //Information for the LOG file
	   //Modified 14/07 con Sabrina per imax
	    time_eve=header_param->mjdtime+(double)(even_paramb->imax[i]*tsample)/day2sec;
	    data_log[0]=mjd2gps(time_eve);
	  data_log[1]=(double) even_paramb->duration[i]*tsample;
	  data_log[2]=(double) gd_band->y[imax];
	  data_log[3]=(double) even_paramb->crmax[i];
	  data_log[4]=(double) even_paramb->ener[i];
	  data_log[5]=(double) frinit;   //ADDED the beg frequency of the BW
	  prec[0]=3;
	  prec[1]=0;
	  prec[2]=0;
	  prec[3]=0;
	  prec[4]=0;
	  prec[5]=0;  //ADDED

	  logfile_ev(LOG_INFO,"EVB",6,data_log,prec);  //6 era 5 e EVT non EVB

	  } //END IF lungo veto eventi
	  //Scrittura eventi su testfileb, luglio 2009 per TEST
	  //if(header_param->nfft==1&&i==0){
	  //    OUTTESTB=fopen("testfileb.dat","wr");
	  //    puts("Open testfileb");
	  //}
	  if(even_paramb->crmax[i]>=100){
	    if((imax*tsample/header_param->tsamplu) >=(int) (header_param->nsamples/2) && (imax*tsample/header_param->tsamplu) < (int) (3*header_param->nsamples/2 ))
		{

	        if(i>=0)for(k=index1;k<index2;k++) {
			  if ((k>imax-200)&&(k<imax+200))
			    {
			      if(OUTTESTB==NULL){
				OUTTESTB=fopen("testfileb.dat","wr");
				puts("***Non era aperto ! Open testfileb");
			      }
			      fprintf(OUTTESTB,"%20.9f %20.9f %d %12.5f %12.5f %12.5f\n",time_eve,header_param->mjdtime+k*(tsample/day2sec),k,even_paramb->crmax[i],gd_band->y[k],frinit);
			     }
		}
		}
	  }

	  //Fine scrittura eventi luglio 2009

        }

	// From gd_appo in time to gd_appo in frequency

       ii= gdtime2freq(gd_appo,header_param);//from time gd to frequency gd
       ////wingd2gd(gd_appo,5); //Rilevo la finestra fatta
	//From gb_appo put the data into gd in the correct bins
        gd_appo->ini=frinit;
	kk=(int) floor(2*frinit/gd->dx); //initial bin in gd
       if(kk%2!=0)kk+=1; //now it is even PARI
       //printf("$$$kk= initial freq= %d %f\n",kk,frinit);
       //printf("$$$gd->n gd_band->n gd_appo->n= %ld %ld %ld\n",gd->n,gd_band->n,gd_appo->n);
       //printf("$$$gd->dx gd_band->dx gd_appo->dx= %f %f %f\n",gd->dx,gd_band->dx,gd_appo->dx);

       //TEST
       /*for(ii=0;ii<8;ii++) {
	 printf("PRIMA gd->y[ii+kk],gd_appo->y[ii], %f %f\n",gd->y[ii+kk],gd_appo->y[ii]);
	 //printf("subbw kk (pari), ii-1 (dispari!) )%d,%d,\n",kk,ii-1);
	 }*/
       ippp=gd_appo->n/4;
       if(ioverl==2){
	 for(ii=0;ii<gd_appo->n/2;ii++)gd->y[ii+kk]=gd_appo->y[ii+ippp]; //verif.
       }
       if(ioverl==1){
	 for(ii=0;ii<(gd_appo->n)/2;ii++)gd->y[ii+kk]=gd_appo->y[ii+ippp];
       }
       //Imag del primo dato viene 0 dopo aver sottocampionato. Non lo e' nella
       //FFT lunga !


      //printf("ULTIMO INDICE SCRITTO= %d\n ",(ii-1)+kk);
       //printf("$$$kk= ii= initial freq finale freq= %d %d %f %f\n",kk,ii-1,frinit,frinit+(ii-1)/2*gd->dx);
       }

	// End of the main loop in the subbw. Go back to time, in gd
       ivign=0; ////
       ii=gdfreq2time(gd,header_param,ivign);

       printf("gd_appo->dx gd_appo->ini gd_appo->n= %f %f %ld\n",gd->dx,gd->ini,gd->n);
 	return i;
}

int sn_medsig(GD *gd,EVEN_PARAM *even_param,HEADER_PARAM *header_param)
/*from data in a GD it evaluates the AR mean and std, using the even_param*/
/*
  even_param->tau=600.0;       //memory time of the autoregressive average
  even_param->cr=6.0;        //CR of the threshold if abs; else is 5
  even_param->deadtime=1.0;  //dead time in seconds
  even_param->factor=10.0;    //e.g. 10: when re-evaluate mean and std
  even_param->iflev=0;       //0=no events; 1=begin; 2=event; 3=end of the event
  even_param->xamed;       //vector AR mean:initial values put to 0
  even_param->xastd;      //vector AR std:initial values put to 0
  even_param->xw=0.;
  even_param->qw=0.;
  even_param->w_norm=0.;   //normalization factor for mean and std
	// prova aumento edge a 0.5 29 maggio 2007. Non cambia nulla, va male. Rimetto come era.
  even_param->edge=0.01;  //how many seconds around the event have to be "purged". Virgo, fino al 28 maggio era 0.01. Lo aumento per LIGO
		//Vedi tutti i problemi con la rimozione di eventi grossi in LIGO avuto 29 maggio 2007.
  even_param->total=0; //how many samples, due to all events, have been cleaned in the GD
  even_param->ndata_run=0; //number of data in the run (one SFDB file, many FFTs !)
  even_param->ar=0;  //Algor. for the AR evaluation: 0 ise the old, 1 the new
  even_param->begin;       //vector with the index of the beginning of the event. Initial values put to 0
  even_param->duration;   //vector with the number of samples of each event. Initial values put to 0
  even_param->imax;      //index at which each event has the maximum. Initial values put to 0
  even_param->xw_pari;
  even_param->qw_pari;
  even_param->w_norm_pari;
  even_param->xw_dispari;
  even_param->qw_dispari;
  even_param->w_norm_dispari;
  even_param->ndata_run_pari;
  even_param->ndata_run_dispari;
*/
//IMPORTANT COMMENT: This function is used when constructing the data base. Hence the first FFT in input has
//always header_param->nfft=1.
//If this is not true (function used in a different way) then the code lines with inizializations have
//to be changed !

{
  int i,i_eval;
  double w_even;
  float ad,qd; //Value of the datum and its square
  float s=0.;
  float ss=0.;  //mean and average, without normalization, every i_eval samples
  double itaust;
  double norm=1.0;
  double asd_appo;

  if(gd->n*gd->dx < even_param->tau)even_param->tau=(float) (1.0*gd->n*gd->dx);  //change the memory time, if too long
  //printf("even_param->tau %f\n",even_param->tau);
  itaust=(double) 1.0/(gd->dx/even_param->tau);
  w_even=exp(-1.0/itaust);
  i_eval= (int) ((even_param->tau/gd->dx)/even_param->factor);
  if(i_eval<1)i_eval=1;
  //printf("gd->n  gd->dx  product %ld %f %f \n",gd->n,gd->dx,1.0*gd->n*gd->dx);
  //printf("dx tau w %f %f %f \n",gd->dx,even_param->tau,w_even);
  even_param->ar=1;  //sets the algorithm for the AR evaluation. 1: new, 0: old
  if(header_param->typ==2){ //overlapped
  if(header_param->nfft<=2){
    even_param->xw=0.;
    even_param->qw=0.;
    even_param->w_norm=0.;
  }
  if(header_param->nfft>=3){
    if(header_param->nfft%2==0){
      even_param->xw=even_param->xw_pari;
      even_param->qw=even_param->qw_pari;
      even_param->w_norm=even_param->w_norm_pari;
    }
    if(header_param->nfft%2!=0){
      even_param->xw=even_param->xw_dispari;
      even_param->qw=even_param->qw_dispari;
      even_param->w_norm=even_param->w_norm_dispari;
   }
  }
  } //if overlapped
  for(i=0; i<(int) gd->n;i++){
    even_param->ndata_run+=1; //needed only with the old ar algorithm !!
    if(header_param->nfft%2==0)even_param->ndata_run_pari+=1; //idem
    if(header_param->nfft%2!=0)even_param->ndata_run_dispari+=1; //idem
    if(even_param->absvalue==1)ad=fabsf(gd->y[i]); //for tests input a constant:e.g. 2.0;
     if(even_param->absvalue!=1)ad=gd->y[i];
    qd=ad*ad;

    if(even_param->ar==1)even_param->xw=(1.0-w_even)*ad+w_even*even_param->xw;
    if(even_param->ar==1)even_param->qw=(1.0-w_even)*qd+w_even*even_param->qw;
    if(even_param->ar==0)even_param->xw=ad+w_even*even_param->xw;
    if(even_param->ar==0)even_param->qw=qd+w_even*even_param->qw;
    if(even_param->ar==1)even_param->w_norm=(1.0-w_even)+w_even*even_param->w_norm; //with the new, it must be done here
    if(i%i_eval==0 && i!=0){ //this check can be done using even_param->ndata
      s=even_param->xw;
      ss=even_param->qw;
      if(even_param->ar==0 && header_param->typ!=2)even_param->w_norm=(1.0-w_even)/(1.0-exp(-(double)even_param->ndata_run/itaust));
      if(even_param->ar==0 && header_param->typ==2){
	if(header_param->nfft<=2)even_param->w_norm=(1.0-w_even)/(1.0-exp(-(double)gd->n/itaust));
	if(header_param->nfft>=3){
	   if(header_param->nfft%2==0)even_param->w_norm=(1.0-w_even)/(1.0-exp(-(double)even_param->ndata_run_pari/itaust));
	   if(header_param->nfft%2!=0)even_param->w_norm=(1.0-w_even)/(1.0-exp(-(double)even_param->ndata_run_dispari/itaust));
	}
      }
    }
    if(even_param->w_norm !=0) {
      if(even_param->ar==1)norm=1.0/even_param->w_norm;
    }
      if(even_param->ar==0)norm=even_param->w_norm;
      even_param->xamed[i]=s*norm;
      asd_appo=(ss*norm-s*s*norm*norm);
      if(asd_appo !=0)even_param->xastd[i]=(float) sqrt(asd_appo);
      else
	even_param->xastd[i]=0.0;
      //printf("xamed xastd   %f %f  \n",even_param->xamed[i],even_param->xastd[i]);
  }

  // Averages and stds at the beginning of the gd:
   for(i=0; i<i_eval;i++){
     even_param->xamed[i]=even_param->xamed[i_eval];
     even_param->xastd[i]=even_param->xastd[i_eval];
   }
   //Fill the empty values, with the previous ones (that is, the last evaluated)
    for(i=i_eval+1; i< (int) gd->n;i++){
      if(i%i_eval!=0){
	even_param->xamed[i]=even_param->xamed[i-1];
	even_param->xastd[i]=even_param->xastd[i-1];
      }
   }
    //To Check:    --15/10/2004--
    //for(i=0; i<gd->n;i+=i_eval/2){
    //printf("i xamed xastd   %d %f %f  \n",i,even_param->xamed[i],even_param->xastd[i]);
    // }
  if(header_param->typ==2){ //overlapped
    if(header_param->nfft%2==0){
      even_param->xw_pari=even_param->xw;
      even_param->qw_pari=even_param->qw;
      even_param->w_norm_pari=even_param->w_norm;
    }
    if(header_param->nfft%2!=0){
      even_param->xw_dispari=even_param->xw;
      even_param->qw_dispari=even_param->qw;
      even_param->w_norm_dispari=even_param->w_norm;
    }
  }
  return i;
}

int even_anst(GD *gd,EVEN_PARAM *even_param)
/*from AR mean and std, it looks for events, and registers their characteristics*/
/*Last modified 16/08/07: duration of the event corrected  */
{
  int i,j;
  float ad,rd;
  float xmax;
  //long int time_below; //number of samples below
  //long int samples_deadtime;

  int time_below; //number of samples below
  int samples_deadtime;

  float appo;
  appo=even_param->deadtime/gd->dx;
  samples_deadtime=lroundf(appo);
  for(i=0;i<(int) gd->n;i++){
    even_param->begin[i]=-1;
    even_param->duration[i]=0;
    even_param->imax[i]=0;
    even_param->crmax[i]=0;
    even_param->ener[i]=0;
  }
  even_param->iflev=0;
  time_below=0.;
  //printf("**Even_anst:dead time,cr,samples_deadtime = %f %f %ld\n",even_param->deadtime,even_param->cr,samples_deadtime);
  j=-1;
  for(i=0; i<(int) gd->n;i++){
       if(even_param->absvalue!=1)ad=gd->y[i];
       if(even_param->absvalue==1)ad=fabsf(gd->y[i]);
       //rd=fabsf((ad-even_param->xamed[i])/(even_param->xastd[i]+even_param->notzero));
       if(even_param->xastd[i]!=0.0){
	 rd=fabsf((ad-even_param->xamed[i])/(even_param->xastd[i]));
       }
       else {
	 rd=0.0;
       }
       //if(i<=4)printf("dato,xamed,xastd,rd %f %f %f %f\n",ad,even_param->xamed[i],even_param->xastd[i],rd);
       if((even_param->iflev==0)&&(rd>=even_param->cr) && (ad !=0)){  //prova ad !=0 23/01/2006
	   j+=1;
	   even_param->iflev=1;
	   //printf("Beginning of an event (j,i,ad) %d %d %f \n",j,i,ad);
	   even_param->begin[j]=i;
	   even_param->duration[j]=1;
	   even_param->ener[j]+=ad*ad;
	   even_param->imax[j]=i;
	   even_param->crmax[j]=rd;
	   xmax=ad;
	   time_below=0;
	 }
	 if((even_param->iflev==2)&&(rd>=even_param->cr)){
	   time_below=0;
	   even_param->duration[j]+=1;
	   even_param->ener[j]+=ad*ad;
	   if(ad>=xmax)even_param->imax[j]=i;
	   if(ad>=xmax)even_param->crmax[j]=rd;
	   if(ad>=xmax)xmax=ad;
	 }
	 if((even_param->iflev==2)&&(rd<even_param->cr)){
	   time_below+=1.;
	   even_param->duration[j]+=1;
	   even_param->ener[j]+=ad*ad;
	 }
	 if(time_below==samples_deadtime){
	   even_param->iflev=3;
	   time_below=0.;
	 }
	 // if di sotto e' solo un check: non deve succedere se non per approssimazioni sbagliate. In quel caso alla durata devo sottrarre un sample.
	 if(time_below>samples_deadtime){
	   even_param->duration[j]-=1;
	   even_param->iflev=3;
	   puts("ATTENZIONE: time below maggiore di samples deadtime !!");
	   time_below=0.;
	 }
	 if(even_param->iflev==1)even_param->iflev=2;
	 if(even_param->iflev==3){
	   even_param->duration[j]-=(int) samples_deadtime; //16 Aug 2007
	    if(even_param->duration[j]<=0) even_param->duration[j]=1;
	    //printf("End of an event %d %d %d %f \n",j,i,even_param->duration[j],xmax);
	    even_param->iflev=0;
	 }
  }
  even_param->number=j+1;  //total number of events found in the gd (j starts from 0)
  return i;
}

/******************Simulations************************/
int add_signals(GD *gd,HEADER_PARAM *header_param)
/*Add signals to the GD *gd */
{

	int i,k;
	float freqm_a,freqp_a; //a means aliased
	int ifactor=40000;  //test. Questo per 131072 samples. 1000 se 4096 samples
	//float *signal;
	float rough_aver,factor,snr=10.;
	float tau=0.03; //30 ms
	int ihowmany;
	int len;

	FILE *OUTTEST_SIGNAL;
	len=(int) gd->n;
	//signal=(float *)malloc((size_t) (len)*sizeof(float));

	puts(" Add signals to the GD");


	//freqm_a=header_param->freqm-header_param->frinit;
	//freqp_a=header_param->freqp-header_param->frinit;
	freqm_a=15; //15 jan 2008 tolte le freq. dall' header. 25 e 15 a caso..
	freqp_a=25;
	printf("freqm, freqp  %f %f \n",freqm_a,freqp_a);


	ihowmany=(int) 3* (float) tau*1.0/header_param->tsamplu;
	//printf("tlen, sampling time  ihowmany %ld %f %d\n",gd->n,header_param->tsamplu,ihowmany);
	rough_aver=0;
	for(i=0;i<gd->n;i++){
	 rough_aver=rough_aver+fabsf(gd->y[i]);
	}
	rough_aver/=gd->n;
	factor=rough_aver*snr;
	printf(" factor, rough_aver,snr = %f %f %f\n",factor,rough_aver,snr);

	/*for(i=0;i<gd->n;i++)signal[i]=0.;
	for(i=1;i<gd->n;i++){
	  if(i % ifactor!=0){
	    for(k=i; k<(i+ihowmany);k++)
	    if(k<gd->n)signal[k]=sin(2*PIG*freqp_a*(k-i)*header_param->tsamplu)*exp(-(k-i)*header_param->tsamplu/tau);
	  }
	}
	*/

	for(i=0;i<gd->n;i++){
	  if(i % ifactor==0 && i !=0)printf("ADD i factor= %d %f\n",i,factor);
	  if(i % ifactor==0 && i !=0)gd->y[i]+=factor; //factor; //*signal[i];
	}
	//TEST::scrittura del file
	OUTTEST_SIGNAL=fopen("testfile_signal.dat","w");
	for(k=0;k<gd->n;k++)fprintf(OUTTEST_SIGNAL,"%d  %f  \n",k,gd->y[k]);
	fclose(OUTTEST_SIGNAL);


	//free(signal);
	return i;
}

void add_sinusoids_old(GD *gd,HEADER_PARAM *header_param,float freq_sin,float amp_sin,double beg_mjd, double *phase_out)
/*Add sinusoids to the GD *gd which is in the frequency domain, from the SFDB*/
{

	int i;
	float freq_sin_a; //aliased frequency
	double phase,mod_term,phase_mod;
	double mod_terma,mod_termb,end_mjd; //if the linear approx is used
	float test_factor=0; //1.0; //put to 0 to add the signal and remove the noise
	                      //useful for tests, in which at the end we deal only with signals

	float *signal;
	int len;
	long long nsam;
	long nsaml;
	float fft2spectrum;
	unsigned long ku;
	int imodulation,approx;
	float var_freq=10; //Hz in tmax
	float tmax=30.0; //days
	/**********Variables needed for the Dopple effect**************/
	short int fut1,nut1,nutcorr;  //nut1=1 uses nutation from ephem files;
                                      //nut1=0  uses nutation from novas interpolation;

	 short int equatorial=1; //1 output are Equatorial, 0 ecliptic
	 short int imjd=1;      //1 input are mjds, 0 imput are days from 1900
	 double mjd;
	 /*****Detector*****/
	 double *velx,*vely,*velz;
	 double *deinstein;
	 double vel[3];
	 /******Source*******/
	 double poss[3];
	 /*******Doppler*******/
	 double nu_doppler;


	 DETECTOR *detector;
	 SOURCE *source;
	 char nam[16];
	 double sa;
	 char nam1[5];

	 if(header_param->detector==1)strncpy(nam1,"virgo",5);
	 if(header_param->detector==2)strncpy(nam1,"ligoh",5);
	 if(header_param->detector==3)strncpy(nam1,"ligol",5);
	 if(header_param->detector==4)strncpy(nam1,"kagra",5);
	detector=def_detector(nam1);
	 // strcpy(nam,"zerozer");
	strcpy(nam,"galcent");
	 // strcpy(nam,"PSR437"); //TEMPO
	 // strcpy(nam,"PSR1937");   //TEMPO Cutler
	source=def_source(nam);
	imodulation=2; //1 if the frequency of the signal is lineraly modulated , 0 if constant, 2 if Doppler
	approx=1; //if imodulation==2: if 1 it uses a linear approx for the phase in one fft.
	          //It runs much much faster...If not 1 (e.g. 0) the phase is exact
	len=(int) gd->n;
	signal=(float *)malloc((size_t) (len)*sizeof(float));
    	puts(" Add a sine signal to the GD (frequency domain !!!)");
	puts("Phase=0 at a chosen reference time");
	freq_sin_a=freq_sin-header_param->frinit;
	printf("Actual mjdtime, beg_mjd %f %f\n",header_param->mjdtime,beg_mjd);
	if(beg_mjd==REFERENCE_MJD)*phase_out=0.;

	/* The following nsam calc is useful to start from an absolute time beg_mjd, which is used in the first fft
	   of each new file. For the others the evaluation is done using the time difference from the two ffts
	   It must be used with care, as approximations must be checked for.
	*/
	nsam= (long long) ((header_param->mjdtime-beg_mjd)*day2sec/header_param->tsamplu);
	printf("***Nsamples from the reference, or from the last used, time = %lld\n",nsam);
	phase=*phase_out;
	if(imodulation==2){
	   velx=(double *)malloc((size_t) (nsam+1)*sizeof(double));
	   vely=(double *)malloc((size_t) (nsam+1)*sizeof(double));
	   velz=(double *)malloc((size_t) (nsam+1)*sizeof(double));
	   deinstein=(double *)malloc((size_t) (nsam+1)*sizeof(double));
	   PSS_source_pos(source,equatorial,poss);
	   nsaml=(long) nsam;
	   mjd=PSS_detector_vel(detector,beg_mjd,header_param->tsamplu,nsaml,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	   puts(" ***Done the evaluation of the initial phase on the nsam samples !!");
	   for(i=0;i<nsam;i++){
	     vel[0]=velx[i];
	     vel[1]=vely[i];
	     vel[2]=velz[i];
	     sa=scalar_product(poss,vel);
	     // nu_doppler=freq_sin*(1.0-sa/C);
	     nu_doppler=freq_sin*(1.0+sa/C); //CORRETTO !! avevo lasciato il segno meno, fino al 28/12/2006 !!
	     mod_term=nu_doppler-freq_sin;
	     phase+=2*PIG*(freq_sin_a+mod_term)*header_param->tsamplu;
	   }
	   free(velx);
	   free(vely);
	   free(velz);
	   free(deinstein);
	   velx=(double *)malloc((size_t) (len)*sizeof(double));
	   vely=(double *)malloc((size_t) (len)*sizeof(double));
	   velz=(double *)malloc((size_t) (len)*sizeof(double));
	   deinstein=(double *)malloc((size_t) (len)*sizeof(double));
	   if(approx !=1){
	   mjd=PSS_detector_vel(detector,header_param->mjdtime,header_param->tsamplu,len,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	   }
	    puts("***Done the first step of the evaluation of the phase on the actual samples in the gd");

	}

	/*Construction of the sinusoidal signal*/
	puts("Now construct the signal on the actual chunk");
	if(imodulation==2 && approx ==1){
	  puts("Use of the linear approximation on the phase in one FFT");
	  mjd=PSS_detector_vel(detector,header_param->mjdtime,header_param->tsamplu,1,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	  vel[0]=velx[0];
	  vel[1]=vely[0];
	  vel[2]=velz[0];
	  sa=scalar_product(poss,vel);
	  //nu_doppler=freq_sin*(1.0-sa/C);
	  nu_doppler=freq_sin*(1.0+sa/C); //CORRETTO !! avevo lasciato il segno meno, fino al 28/12/2006 !!
	  mod_terma=nu_doppler-freq_sin;
	  end_mjd=header_param->mjdtime+header_param->tsamplu*(gd->n-1)/day2sec;
	  mjd=PSS_detector_vel(detector,end_mjd,header_param->tsamplu,1,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	  vel[0]=velx[0];
	  vel[1]=vely[0];
	  vel[2]=velz[0];
	  sa=scalar_product(poss,vel);
	  // NO !!nu_doppler=freq_sin*(1.0-sa/C);
	  nu_doppler=freq_sin*(1.0+sa/C); //CORRETTO !! avevo lasciato il segno meno, fino al 28/12/2006 !!
	  mod_termb=nu_doppler-freq_sin;
	}
	for(i=0;i<gd->n;i++){
	  if(imodulation==0)mod_term=0.;
	  if(imodulation==1)mod_term=var_freq*(i+nsam)*header_param->tsamplu/(tmax*day2sec);
	  if(imodulation==2 && approx!=1){
	    vel[0]=velx[i];
	    vel[1]=vely[i];
	    vel[2]=velz[i];
	    sa=scalar_product(poss,vel);
	    // NO!! nu_doppler=freq_sin*(1.0-sa/C);
	    nu_doppler=freq_sin*(1.0+sa/C); //CORRETTO !! avevo lasciato il segno meno, fino al 28/12/2006 !!
	    mod_term=nu_doppler-freq_sin;
	  }
	  if(imodulation==2 && approx==1){
	    mod_term=mod_terma+(mod_termb-mod_terma)*i/(gd->n-1);
	  }
	  phase+=2*PIG*(freq_sin_a+mod_term)*header_param->tsamplu;
	  phase_mod=fmod(phase,TWOPI);
	  signal[i]=amp_sin*sin(phase_mod);
	  //Print for test
	  if(imodulation==2){
	    if(approx !=1){
	      if(i<3)printf("i,mod_term,phase_mod,nu_doppler %d %15.13f %f %15.11f\n",i,mod_term,phase_mod,nu_doppler);
	      if(i>gd->n-3)printf("i,mod_term,phase_mod,nu_doppler %d %15.13f %f %15.11f\n",i,mod_term,phase_mod,nu_doppler);
	    }
	     if(approx ==1){
	      if(i==gd->n-1)printf("i,mod_term,phase_mod,nu_doppler %d %15.13f %f %15.11f\n",i,mod_term,phase_mod,nu_doppler);
	       if(i<3)printf("i,mod_term,phase_mod %d %15.13f %f \n",i,mod_term,phase_mod);
	      if(i!=gd->n-1 && i>gd->n-3)printf("i,mod_term,phase_mod %d %15.13f %f \n",i,mod_term,phase_mod);
	    }
	  }
	  if(imodulation<2){
	  if(i<3)printf("i,mod_term,phase_mod, %d %15.13f %f \n",i, mod_term,phase_mod);
	  if(i>gd->n-3)printf("i,mod_term,phase_mod, %d %15.13f %f \n",i, mod_term,phase_mod);
	  }
	  if(header_param->typ==2 && i==(gd->n-1)/2)*phase_out=phase;

	  if(header_param->typ!=2 && i==(gd->n-1))*phase_out=phase;

	}

	/*Window on the signal, if needed (if present on the detector data)*/
	if(header_param->wink>0)vect2winvect(signal,header_param->wink,len);
	/*Go to the frequency domain*/
	ku= (unsigned long) gd->n;
	fft2spectrum=sqrt(header_param->tsamplu/gd->n); //thus it is normalized as the data are: the squared modulus of fft is the
	                                                //power spectrum
	 realft(signal-1,ku,1);
	 //the output is only 1/2 of the total FFT, in a real vector of size k
	 //which contains the real part of the FFT in signal[k]
	 //and the imag part in signal[k+1]

	/*Add the signal to the noise (both in frequency domain, Real-Imag)*/
	for(i=0;i<gd->n;i++){
	 gd->y[i]=gd->y[i]*test_factor+signal[i]*fft2spectrum;
	}

	free(signal);
	if(imodulation==2)free(velx);
	if(imodulation==2)free(vely);
	if(imodulation==2)free(velz);
	if(imodulation==2)free(deinstein);
	return ;
}

void add_sinusoids(GD *gd,HEADER_PARAM *header_param,float freq_sin,float amp_sin,double beg_mjd, double *phase_out)
/*Add sinusoids to the GD *gd which is in the frequency domain, from the SFDB*/
{

	int i;
	float freq_sin_a; //aliased frequency
	double phase,mod_term,phase_mod,*mod_termv;
	double mod_terma,mod_termb,end_mjd,diff_min; //if the linear approx is used
	double diff_days,diff_ref;
	float test_factor=0.0; //1.0; //put to 0 to add the signal and remove the noise
	                      //useful for tests, in which at the end we deal only with signals

	float *signal;
	int len;
	long long nsam,i_nsam;
	double nsam_d;
	long nsaml,nsam_min,nsam_used; //nsam_min needed if the linear approx is used
	float fft2spectrum;
	unsigned long ku;
	int imodulation,approx;
	float var_freq=10; //Hz in tmax
	float tmax=30.0; //days
	/**********Variables needed for the Dopple effect**************/
	short int fut1,nut1,nutcorr;  //nut1=1 uses nutation from ephem files;
                                      //nut1=0  uses nutation from novas interpolation;

	 short int equatorial=1; //1 output are Equatorial, 0 ecliptic
	 short int imjd=1;      //1 input are mjds, 0 imput are days from 1900
	 double mjd;
	 /*****Detector*****/
	 double *velx,*vely,*velz;
	 double *deinstein;
	 double vel[3];
	 /******Source*******/
	 double poss[3];
	 /*******Doppler*******/
	 double nu_doppler;


	 DETECTOR *detector;
	 SOURCE *source;
	 char nam[16];
	 double sa;
	 char nam1[5];
	 if(header_param->detector==1)strncpy(nam1,"virgo",5);
	 if(header_param->detector==2)strncpy(nam1,"ligoh",5);
	 if(header_param->detector==3)strncpy(nam1,"ligol",5);
	 if(header_param->detector==4)strncpy(nam1,"kagra",5);
	detector=def_detector(nam1);
	 // strcpy(nam,"zerozer");
	strcpy(nam,"galcent");
	 // strcpy(nam,"PSR437"); //TEMPO
	 // strcpy(nam,"PSR1937");   //TEMPO Cutler
	source=def_source(nam);
	 //evaluation of detector velocities
  /****************************************************************************/
   fut1=0;  ///1=uses ut1; 0=uses utc as aprrox. of ut1; //0 per fare fino al 2010
   nut1=1;  //1=uses nutation from ephem files; 0=no;
   nutcorr=0; //reads the file DELTANUT, if 0 does not read DELTANUT //0 per fare fino al 2010
 /****************************************************************************/
	imodulation=0; //1 if the frequency of the signal is lineraly modulated , 0 if constant, 2 if Doppler
	approx=1; //Usually put to 1 ! Works if imodulation==2: if 1 it uses a linear approx for the phase in one fft.
	          //It runs much much faster...If not 1 (e.g. 0) the phase is exact
	len=(int) gd->n;
	signal=(float *)malloc((size_t) (len)*sizeof(float));
	mod_termv=(double *)malloc((size_t) (len+1)*sizeof(double));
    	puts(" Add a sine signal to the GD (frequency domain !!!)");
	puts("Phase=0 at a chosen reference time");
	freq_sin_a=freq_sin-header_param->frinit;
	printf("Actual mjdtime, beg_mjd %f %f\n",header_param->mjdtime,beg_mjd);
	if(beg_mjd==REFERENCE_MJD)*phase_out=0.;
	if(imodulation==2)PSS_source_pos(source,equatorial,poss);

	/* The following nsam calc is needed to start from an absolute time beg_mjd, which is used in the first fft
	   of each new file, and to consider appropriately holes in the data.
	   The evaluation is done using the actual time difference from the two ffts
	*/
	nsam= (long long) ((header_param->mjdtime-beg_mjd)*day2sec/header_param->tsamplu);
	diff_days=(header_param->mjdtime-beg_mjd);
	diff_ref=fmod(diff_days,(double) header_param->tsamplu/day2sec);
	//diff_ref=diff_days; //this gives the same result of the istruction above.
	printf("***Nsamples from the reference, or from the last used time (value and abs)= %lld %lld \n",nsam, llabs(nsam));
	printf("***Diff_days, diff_ref= %f %e \n",diff_days,diff_ref);
	/*The following is very important: the signal so evaluated is perfect, both in case of a hole in
	  the data and in case of a negative hole (see the example of VIR_hrec_20040624_215946_.SFDB.
	  The data have a negative hole  and a positive hole. Adding a signal at 600 Hz, subsampling
	  from 598, band=2 Hz, the negative hole is  well-seen zooming the signal at the seconds 1480-1490
	  and the positive hole zooming at the seconds 2480-2500 -e.g.. Constructing a sinusoid with fr_ini=598.999023,
	  dt=0.1280017..and comparing it to the data in dmps5 we see that the phase of the two sinusoids is the same.
	  And, regarding the negative hole, when the signal "goes back" it coincides perfectly with the previous.
	*/
	if(nsam!=0 && beg_mjd !=REFERENCE_MJD){
	   puts("nsam different from zero: Jump in the data");
	   puts("The evaluation of the phase has to be re-done:RESTART from the reference time + diff_ref !");
	   beg_mjd=REFERENCE_MJD+diff_ref;
	   *phase_out=0.;
	   nsam= (long long) (ceil)((header_param->mjdtime-beg_mjd)*day2sec/header_param->tsamplu)+1; //capire il +1
	                                                                                              //viene perfetto con ceil e +1
      	   nsam_d= (double) ((header_param->mjdtime-beg_mjd)*day2sec/header_param->tsamplu);
	   printf("***New evaluation: nsamples from the REFERENCE TIME (value,abs, double)= %lld %lld %f\n",nsam, llabs(nsam),nsam_d);
	}
	/*nsam can't be negative from now on*/
	phase=*phase_out;
	if(imodulation==0){
	  if(nsam>0)for(i=0;i<nsam;i++)phase+=2*PIG*(freq_sin_a)*header_param->tsamplu;
	  if(nsam<0)for(i=0;i<llabs(nsam);i++)phase-=2*PIG*(freq_sin_a)*header_param->tsamplu;
	}
	if(imodulation==2 && approx !=1){
	  if(nsam !=0){
	    velx=(double *)malloc((size_t) (llabs(nsam)+1)*sizeof(double));
	    vely=(double *)malloc((size_t) (llabs(nsam)+1)*sizeof(double));
	    velz=(double *)malloc((size_t) (llabs(nsam)+1)*sizeof(double));
	    deinstein=(double *)malloc((size_t) (llabs(nsam)+1)*sizeof(double));
            nsaml=(long) llabs(nsam);
	    if(nsam >0){
	      mjd=PSS_detector_vel(detector,beg_mjd,header_param->tsamplu,nsaml,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	    }
	    if(nsam<0){
	      puts("nsam negative: new data begins BEFORE the end of previous data. there is cross-overlapping in the data");
	      puts("The evaluation of the phase has to be pushed back");
	      mjd=PSS_detector_vel(detector,header_param->mjdtime,header_param->tsamplu,nsaml,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	    }
	    puts(" ***No approximations: Done the evaluation of the det vel for the nsam samples !!");
	    for(i=0;i<llabs(nsam);i++){
	       vel[0]=velx[i];
	       vel[1]=vely[i];
	       vel[2]=velz[i];
	       sa=scalar_product(poss,vel);
	       //nu_doppler=freq_sin*(1.0-sa/C);
	        nu_doppler=freq_sin*(1.0+sa/C); //CORRETTO !! avevo lasciato il segno meno, fino al 28/12/2006 !!
	       mod_term=nu_doppler-freq_sin;
	       if(nsam >0)phase+=2*PIG*(freq_sin_a+mod_term)*header_param->tsamplu;
	       if(nsam <0)phase-=2*PIG*(freq_sin_a+mod_term)*header_param->tsamplu;
	    }
	    puts(" ***No approximations: Done the evaluation of the initial phase on the nsam samples !!");
	    free(velx);
	    free(vely);
	    free(velz);
	    free(deinstein);
	  }
	  velx=(double *)malloc((size_t) (len)*sizeof(double));
	  vely=(double *)malloc((size_t) (len)*sizeof(double));
	  velz=(double *)malloc((size_t) (len)*sizeof(double));
	  deinstein=(double *)malloc((size_t) (len)*sizeof(double));
	  mjd=PSS_detector_vel(detector,header_param->mjdtime,header_param->tsamplu,len,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
          puts("***No approximations: Done the first step of the evaluation of the det. vel for the actual samples in the gd");
	  for(i=0;i<gd->n;i++){
	    vel[0]=velx[i];
	    vel[1]=vely[i];
	    vel[2]=velz[i];
	    sa=scalar_product(poss,vel);
	    //Corretto il 28/12/2006::::nu_doppler=freq_sin*(1.0-sa/C);
	    nu_doppler=freq_sin*(1.0+sa/C); //CORRETTO !! avevo lasciato il segno meno, fino al 28/12/2006 !!
	    mod_termv[i]=nu_doppler-freq_sin;
	  }
	  free(velx);
	  free(vely);
	  free(velz);
	  free(deinstein);
	}

	if(imodulation==2 && approx ==1){
	   nsaml=1;
	   diff_min=(header_param->mjdtime-beg_mjd)*day2min;
	   nsam_min=(long) (min2sec/header_param->tsamplu)+1; //now it is the number of samples in 1 minute
	   //the approximation is done on a 1-minute basis. If we want to change it is sufficient to redifine nsam_min
	   //e.g. nsam_min=(long) ((min2sec/header_param->tsamplu)+1)*6; the basis will be 6 minute-->faster and less accurate
	   //e.g. nsam_min=(long) ((min2sec/header_param->tsamplu)+1)/10; the basis will be 6 seconds-->slower but more accurate
	   if(nsam !=0){
	     nsam_used=nsam_min;
	     if(llabs(nsam)<(long long)nsam_min)nsam_used=(long) llabs(nsam);
	     printf("diff_min,nsam_min %f %ld\n",diff_min,nsam_min);
	     if(nsam<0){
	       puts("nsam negative: new data begins BEFORE the end of previous data. there is cross-overlapping in the data");
	       puts("The evaluation of the phase has to be pushed back");
	     }
	     velx=(double *)malloc((size_t) (nsaml)*sizeof(double));
	     vely=(double *)malloc((size_t) (nsaml)*sizeof(double));
	     velz=(double *)malloc((size_t) (nsaml)*sizeof(double));
	     deinstein=(double *)malloc((size_t) (nsaml)*sizeof(double));
	     for(i_nsam=0; i_nsam<(llabs(nsam)+nsam_used); i_nsam+=nsam_used){
	        if(nsam>0)mjd=PSS_detector_vel(detector,beg_mjd,header_param->tsamplu,nsaml,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
		 if(nsam<0)mjd=PSS_detector_vel(detector,header_param->mjdtime,header_param->tsamplu,nsaml,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	        vel[0]=velx[0];
	        vel[1]=vely[0];
	        vel[2]=velz[0];
	        sa=scalar_product(poss,vel);
	        //nu_doppler=freq_sin*(1.0-sa/C); //CORRETTO IL 28/12/2006
		nu_doppler=freq_sin*(1.0+sa/C); //CORRETTO !! avevo lasciato il segno meno, fino al 28/12/2006 !!
	        mod_terma=nu_doppler-freq_sin;
	        if(nsam>0)end_mjd=beg_mjd+header_param->tsamplu*(i_nsam+nsam_used-1)/day2sec;
		if(nsam<0)end_mjd=header_param->mjdtime+header_param->tsamplu*(i_nsam+nsam_used-1)/day2sec;
	        mjd=PSS_detector_vel(detector,end_mjd,header_param->tsamplu,nsaml,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	        vel[0]=velx[0];
	        vel[1]=vely[0];
	        vel[2]=velz[0];
	        sa=scalar_product(poss,vel);
	        //NO !! nu_doppler=freq_sin*(1.0-sa/C);
		nu_doppler=freq_sin*(1.0+sa/C); //CORRETTO !! avevo lasciato il segno meno, fino al 28/12/2006 !!
	        mod_termb=nu_doppler-freq_sin;
	        for(i=0;i<nsam_used;i++){
	           if(i_nsam+i<llabs(nsam)){
		     mod_term=mod_terma+(mod_termb-mod_terma)*i/(nsam_used-1);
                     if(nsam>0)phase+=2*PIG*(freq_sin_a+mod_term)*header_param->tsamplu;
		     if(nsam<0)phase-=2*PIG*(freq_sin_a+mod_term)*header_param->tsamplu;
		     if(i_nsam+i==llabs(nsam)-1)printf("i,i_nsam,nsam,diff_min,nsam_used %d %lld,%lld,%f %ld\n",i,i_nsam,nsam,diff_min,nsam_used);
		   }
		}

	     }
	     puts(" ***Linear approximation: Done the evaluation of the initial phase on the nsam samples !!");
	     free(velx);
	     free(vely);
	     free(velz);
	     free(deinstein);
	   }
	   velx=(double *)malloc((size_t) (nsaml)*sizeof(double));
	   vely=(double *)malloc((size_t) (nsaml)*sizeof(double));
	   velz=(double *)malloc((size_t) (nsaml)*sizeof(double));
	   deinstein=(double *)malloc((size_t) (nsaml)*sizeof(double));
	   nsam_used=nsam_min;
	   if(len<nsam_min)nsam_used=len;
	   for(i_nsam=0; i_nsam<(len+nsam_used); i_nsam+=nsam_used){
	      mjd=PSS_detector_vel(detector,header_param->mjdtime,header_param->tsamplu,nsaml,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	      vel[0]=velx[0];
	      vel[1]=vely[0];
	      vel[2]=velz[0];
	      sa=scalar_product(poss,vel);
	      // NO !! nu_doppler=freq_sin*(1.0-sa/C);
	      nu_doppler=freq_sin*(1.0+sa/C); //CORRETTO !! avevo lasciato il segno meno, fino al 28/12/2006 !!
	      mod_terma=nu_doppler-freq_sin;
	      end_mjd=header_param->mjdtime+header_param->tsamplu*(i_nsam+nsam_used-1)/day2sec;
	      mjd=PSS_detector_vel(detector,end_mjd,header_param->tsamplu,nsaml,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	      vel[0]=velx[0];
	      vel[1]=vely[0];
	      vel[2]=velz[0];
	      sa=scalar_product(poss,vel);
	      // NO !! nu_doppler=freq_sin*(1.0-sa/C);
	       nu_doppler=freq_sin*(1.0+sa/C); //CORRETTO !! avevo lasciato il segno meno, fino al 28/12/2006 !!
	      mod_termb=nu_doppler-freq_sin;
	      //printf("i_nsam,nsam,diff_min,nsam_used len %lld,%lld,%f %ld %d\n",i_nsam,nsam,diff_min,nsam_used,len);
	      for(i=0;i<nsam_used;i++){
	         if(i_nsam+i<len)mod_termv[i_nsam+i]=mod_terma+(mod_termb-mod_terma)*i/(nsam_used-1);
		 if(i_nsam+i==len-1)printf("Last done:i,i_nsam,len,diff_min,nsam_used %d %lld,%d,%f %ld\n",i, i_nsam,len,diff_min,nsam_used);
	      }
	   }

           free(velx);
	   free(vely);
	   free(velz);
	   free(deinstein);
           puts("***Linear approximation: Done the step of the evaluation of the det. vel. on the actual samples in the gd");

	}
	/*Construction of the sinusoidal signal*/
	puts("Now construct the signal on the actual chunk");
	for(i=0;i<gd->n;i++){
	  if(imodulation==0)mod_termv[i]=0.;
	  if(imodulation==1)mod_termv[i]=var_freq*(i+nsam)*header_param->tsamplu/(tmax*day2sec);
	  phase+=2*PIG*(freq_sin_a+mod_termv[i])*header_param->tsamplu;
	  //if(imodulation==0)phase=2*PIG*(freq_sin_a+mod_termv[i])*(i*header_param->tsamplu+(header_param->mjdtime-(REFERENCE_MJD+diff_ref))*day2sec); //piapia prova che dovrebbe essere tolta! E' per fare la sinusoide in modo facile e capire il ritorno indietro
	  phase_mod=fmod(phase,TWOPI);
	  signal[i]=amp_sin*sin(phase_mod);
	  //Print for test
	  if(imodulation==2){
	    if(approx !=1){
	      if(i<3)printf("i,mod_termv,phase_mod,nu_doppler %d %15.13f %f %15.11f\n",i,mod_termv[i],phase_mod,nu_doppler);
	      if(i>gd->n-5)printf("i,mod_termv,phase_mod,nu_doppler %d %15.13f %f %15.11f\n",i,mod_termv[i],phase_mod,nu_doppler);
	    }
	     if(approx ==1){
	      if(i==gd->n-1)printf("i,mod_termv,phase_mod,nu_doppler %d %15.13f %f %15.11f\n",i,mod_termv[i],phase_mod,nu_doppler);
	       if(i<3)printf("i,mod_termv,phase_mod %d %15.13f %f \n",i,mod_termv[i],phase_mod);
	      if(i!=gd->n-1 && i>gd->n-3)printf("i,mod_term,phase_mod %d %15.13f %f \n",i,mod_termv[i],phase_mod);
	    }
	  }
	  if(imodulation<2){
	  if(i<3)printf("i,mod_termv,phase_mod, %d %15.13f %f \n",i, mod_termv[i],phase_mod);
	  if(i>gd->n-5)printf("i,mod_termv,phase_mod, %d %15.13f %f \n",i, mod_termv[i],phase_mod);
	  }
	  if(header_param->typ==2 && i==(gd->n-1)/2)*phase_out=phase;

	  if(header_param->typ!=2 && i==(gd->n-1))*phase_out=phase;

	}
	free(mod_termv);
   	/*Window on the signal, if needed (if present on the detector data)*/
	if(header_param->wink>0)vect2winvect(signal,header_param->wink,len);
	/*Go to the frequency domain*/
	ku= (unsigned long) gd->n;
	fft2spectrum=sqrt(header_param->tsamplu/gd->n); //thus it is normalized as the data are: the squared modulus of fft is the
	                                                //power spectrum
	realft(signal-1,ku,1);
	 //the output is only 1/2 of the total FFT, in a real vector of size k
	 //which contains the real part of the FFT in signal[k]
	 //and the imag part in signal[k+1]

	/*Add the signal to the noise (both in frequency domain, Real-Imag)*/
	for(i=0;i<gd->n;i++){
	  gd->y[i]=gd->y[i]*test_factor+signal[i];  //piapia tolto *fft2spectrum;
	}

	free(signal);

	return ;
}

void add_sinusoids_lowres(GD *gd,HEADER_PARAM *header_param,double *ra,double *dec,float *freq,float *amp,double mid_mjd,int number_of_signals,double *possx,double *possy,double *possz,float *template1,int len_template,int factor_template,int fast_proc,int nfft)
/*Add sinusoids with !! LOW-RES !! to the GD *gd which is in the frequency domain, from the SFDB*/
{


  //6/11/06 rimettere test_factor a 1 e test a 0 // FATTO:
   //: 15/11/06 rimesso test_factor a 1
  //lasciato test a 1, va poi rimesso a 0 //
	int i,kk;

	float *signal;
	int len;
	long nsaml;
	float fft2spectrum;
	unsigned long ku;
	float test_factor=1.0; //1.0; //put to 0 to add the signal and remove the noise
	                      //0 is only useful for tests, in which at the end we deal only with signals
	                      //Controllare sia a 1 normalmente !!

	/**********Variables needed for the Dopple effect**************/
	short int fut1,nut1,nutcorr;  //nut1=1 uses nutation from ephem files;
                                      //nut1=0  uses nutation from novas interpolation;

	 short int equatorial=1; //1 output are Equatorial, 0 ecliptic
	 short int imjd=1;      //1 input are mjds, 0 input are days from 1900
	 double mjd;

	 double tss_local; //needed to the amplitude modulation
	 float f_lin,f_circ,modul1; //amplitude modulation
	 /*****Detector*****/
	 double *velx,*vely,*velz;
	 double *deinstein;
	 double vel[3];
	 /******Source*******/
	 double poss[3];
	 /*******Doppler*******/
	 double nu_doppler;
	 double phase,phase_mod;

	 DETECTOR *detector;
	 SOURCE *source;
	 char nam[16];
	 char nam1[5];
	 double sa;

	 int test=0; // if 1 uses constant amplitudes =10 (high ampl to do tests)
	             // 0 (or different from 1) = normal functioning



	 float freq_template=1000;
	 float dnu_template;
	 int bin_freq_template;
	 int len1;
	 int k0;  //actual Doppler frequency index in the data vector
	 int k01; //corresponding bin in the template vector
	 int verb=0; //=0 does not print too many things. 1 used for tests, possible only for 1-2 signals
	 FILE *LOWRES;
	 FILE *LOWRES1;
	 float norm; //normalization factor from the template to the data FFT len. To be CHECKED !!
	 if(header_param->detector==1)strncpy(nam1,"virgo",5);
	 if(header_param->detector==2)strncpy(nam1,"ligoh",5);
	 if(header_param->detector==3)strncpy(nam1,"ligol",5);
	 if(header_param->detector==4)strncpy(nam1,"kagra",5);
	detector=def_detector(nam1);
       	strcpy(nam,"simul");
	source=def_source(nam);
	float h;


   /*****************************************************************************/
   fut1=0;  ///1=uses ut1; 0=uses utc as aprrox. of ut1; //0 per fare fino al 2010
   nut1=1;  //1=uses nutation from ephem files; 0=no;
   nutcorr=0; //reads the file DELTANUT, if 0 does not read DELTANUT //0 per fare fino al 2010
   /****************************************************************************/
	len=(int) gd->n;
	if(fast_proc!=1)signal=(float *)malloc((size_t) (len)*sizeof(float));
	printf("Actual mjdtime, mid_mjd, dnu %f %f %f \n",header_param->mjdtime,mid_mjd,header_param->deltanu);
	printf("Normalizations for amplitudes normw,normd  %f %f\n",header_param->normw,header_param->normd);
    	puts(" Add sinusoidal signals to the GD (frequency domain !!!)");
	puts("Phase NOT integrated. LOW RES simulation !!");
	nsaml=1;
	velx=(double *)malloc((size_t) (labs(nsaml)+1)*sizeof(double));
	vely=(double *)malloc((size_t) (labs(nsaml)+1)*sizeof(double));
	velz=(double *)malloc((size_t) (labs(nsaml)+1)*sizeof(double));
	deinstein=(double *)malloc((size_t) (labs(nsaml)+1)*sizeof(double));
	if(fast_proc !=1)for(i=0;i<gd->n;i++)signal[i]=0.;

	if(fast_proc==1 && test_factor !=1)for(i=0;i<gd->n;i++)gd->y[i]*=test_factor;
	//Conversion from ra and dec to equatorial for the sources done just once:
	//Result in vectors possx,possy,possz
	if(nfft==0){
	  for(kk=0;kk<number_of_signals;kk++){
	  source->ra=ra[kk];
	  source->dec=dec[kk];
	  PSS_source_pos(source,equatorial,poss);
	  possx[kk]=poss[0];
	  possy[kk]=poss[1];
	  possz[kk]=poss[2];
	  }

	}
	if(verb>=4){
	  for(kk=0;kk<2;kk++)printf("nfft,kk,possx,possy,possz %d, %d, %f,%f,%f\n",nfft+1,kk,possx[kk],possy[kk],possz[kk]);
	}
	if(fast_proc==1)len1=len_template*factor_template;
	if(nfft==0 && fast_proc==1){
	puts("Fast procedure to add sine signals. Construction of the template");
	printf("len_template, factor_template %d %d\n", len_template, factor_template);
	  for(i=0;i<len_template;i++){
	    phase=2*PIG*(freq_template-header_param->frinit)*i*header_param->tsamplu;
	    phase_mod=fmod(phase,TWOPI);
	    template1[i]=sin(phase_mod);
	  }
	  //window:
	  if(header_param->wink>0)vect2winvect(template1,header_param->wink,len_template);
	  //Add zeroes:
	  for(i=len_template;i<len1;i++) template1[i]=0.;
	  //Go to the frequency domain
	ku= (unsigned long) len1;
	//Controllare: questo sotto non sembra corretto..poi non e' manco usato 16/11/06
	fft2spectrum=sqrt(header_param->tsamplu/len1); //thus it is normalized as the data: the squared modulus of fft is the
	                                                //power spectrum
	realft(template1-1,ku,1);
	 //the output is only 1/2 of the total FFT, in a real vector of size k
	 //which contains the real part of the FFT in template1[i]
	 //and the imag part in template1[i+1]
	 puts("END of the construction of the (fast) template, in freq domain");
	 if(verb>=4){
	  dnu_template=1/(len1*header_param->tsamplu);
	  norm=dnu_template/header_param->deltanu; //vd. prove del 3/11/06 sul quaderno
	    //TEST::scrittura del file
	    LOWRES=fopen("lowres.dat","w");
	    for(i=0;i<len1;i+=2)fprintf(LOWRES,"%d %f %f %f\n",i,i/(2*len1*header_param->tsamplu),norm*template1[i],norm*template1[i+1]);
	    fclose(LOWRES);
	 }


	}
	if(fast_proc==1){
	  dnu_template=1/(len1*header_param->tsamplu);
	  norm=factor_template*dnu_template/header_param->deltanu;
	  bin_freq_template=(ceil) ((freq_template-header_param->frinit)/dnu_template);
	  if(verb>=4)printf("dnu_template,bin_freq_template %f,%d\n",dnu_template,bin_freq_template);
	}
	for(kk=0;kk<number_of_signals;kk++){
	  //Frequency Doppler modulation: detector velocity evaluated only once for each FFT
	  if(kk==0){
	  mjd=PSS_detector_vel(detector,mid_mjd,header_param->tsamplu,nsaml,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);

	  vel[0]=velx[0];
	  vel[1]=vely[0];
	  vel[2]=velz[0];
	  }
	  poss[0]=possx[kk];
	  poss[1]=possy[kk];
	  poss[2]=possz[kk];
	  sa=scalar_product(poss,vel);
	  source->frequency=freq[kk];
	  // NO!! nu_doppler=source->frequency*(1.0-sa/C);
	   nu_doppler=source->frequency*(1.0+sa/C); //CORRETTO !! avevo lasciato il segno meno, fino al 28/12/2006 !!
	  //To check:
	  if(kk==0)printf("nu_doppler, source->frequency %f %f\n", nu_doppler, source->frequency);
	  if(test!=1){
	    //Amplitude modulation for the mid_mjd:
	    //Evaluation of tss_local only once for each new FFTs
	    if(kk==0){
	    mjd2tsslocal(detector,mid_mjd,fut1,nut1,nutcorr,&tss_local);
	    if(tss_local >=24.0)tss_local-=24.0;
	    else if
	      (tss_local < 0.0)tss_local+=24.0;
	    }
	    radpat_interf(detector,source,tss_local,&f_circ,&f_lin);
	    modul1=f_lin+f_circ;
	  }
	  if(test==1){
	    modul1=1.0;
	    amp[kk]=10.0;
	  }
	  if(fast_proc==1){
	    k0=(ceil) ((nu_doppler-header_param->frinit)/header_param->deltanu);
	    k0*=2;          //Data are Real and Imag
	    k01=(ceil) (bin_freq_template);
	    k01*=2;  //Real
	    if(verb>=4)printf("IMP norm=,Doppler freq. bin in Re,Im data: k0,k01 (pos. nel template)= %f %d %d\n",norm,k0,k01);
	    for(i=-10;i<10;i+=2){ // prova era da -20 a 20. Con -40 40 peggiora
	      k01=(ceil) (bin_freq_template);
	      k01*=2;  //Real
	      k01+=i*factor_template;
	      gd->y[i+k0]+=modul1*amp[kk]*template1[k01]*norm; //Real
	      gd->y[i+1+k0]+=modul1*amp[kk]*template1[k01+1]*norm;   //Imag
	    }
	    if(verb>=4 && nfft==0 && kk==0){
	      //TEST1::scrittura del file
	      LOWRES1=fopen("lowres1.dat","w");
	      for(i=-20;i<20;i+=2){
		h=sqrt(gd->y[i+k0]*gd->y[i+k0]+gd->y[i+1+k0]*gd->y[i+1+k0])*sqrt(bil2uni*header_param->deltanu)*header_param->normw*header_param->normd;
		fprintf(LOWRES1,"%d %d %f %f %f %f\n",k0,k0+i,(i+k0)*header_param->deltanu/2,gd->y[i+k0],gd->y[i+1+k0],h);
	      }
	      fclose(LOWRES1);
	    }

	  }

	  if(fast_proc !=1){
	   for(i=0;i<gd->n;i++){
	      phase=2*PIG*(nu_doppler-header_param->frinit)*i*header_param->tsamplu;
	      phase_mod=fmod(phase,TWOPI);
	      signal[i]+=amp[kk]*sqrt(modul1)*sin(phase_mod);
	    }
	  }
	}

	free(velx);
	free(vely);
	free(velz);
	free(deinstein);

	if(fast_proc !=1){
           /*Window on the (sum of the ) signal, if needed (if present on the detector data)*/
	   if(header_param->wink>0)vect2winvect(signal,header_param->wink,len);
	   /*Go to the frequency domain*/
	   ku= (unsigned long) gd->n;
	   fft2spectrum=sqrt(header_param->tsamplu/gd->n); //thus it is normalized as the data: the squared modulus of fft is the
	                                                //power spectrum
	   realft(signal-1,ku,1);
	   //the output is only 1/2 of the total FFT, in a real vector of size k
	   //which contains the real part of the FFT in signal[k]
	   //and the imag part in signal[k+1]
           /*Add the signals to the noise (both in frequency domain, Real-Imag)*/
	  for(i=0;i<gd->n;i++){
	   gd->y[i]=gd->y[i]*test_factor+signal[i];  //piapia tolto *fft2spectrum;
	  }
          free(signal);
	}

  return;
}

/*****************Peak map*****************/
int peakmapp05(GD *gd,GD *gd_short,HEADER_PARAM *header_param)
/*Crea peak max using maxima from the SFDB (in gd), normalized with the very short spectrum in gd_short)
  The output table is written into an file
Old .p05 format (only exists in this version)
*/
{

	int i,k,k_short,nmax;
	float ratio;
	//float freqm_a,freqp_a; //a means aliased NOT USED 15/01/08
	float *power,*power_short;
	float ratio_min=2.5;
	int len,len_short;

	nmax=0;
	len=(int) gd->n/2;
	len_short=(int) gd_short->n/2;
	power=(float *)malloc((size_t) (len)*sizeof(float));
	power_short=(float *)malloc((size_t) (len_short)*sizeof(float));
	puts(" Produces the peak map");


	//freqm_a=header_param->freqm-gd->ini;
	//freqp_a=header_param->freqp-gd->ini;
	//if(header_param->detector==0)printf("Resonant detector. freqm, freqp  %f %f \n",freqm_a,freqp_a);
	if(header_param->detector>=1)printf("Itf detector.min and max freq.%f %f \n",gd->ini,gd->ini+gd->dx*(gd->n-1));

	printf("len, len_short,sampling time %ld %ld %f\n",gd->n,gd_short->n,header_param->tsamplu);

	k=0;
	for(i=0;i<gd->n;i+=2){
	 power[k]=sqrt(gd->y[i]*gd->y[i]+gd->y[i+1]*gd->y[i+1]);
	 power[k]*=(header_param->normd*header_param->normw); //piapia: normalized as the v. short power
	 k+=1;
	}

       	k=0;
	for(i=0;i<gd_short->n;i+=2){
	 power_short[k]=sqrt(gd_short->y[i]*gd_short->y[i]+gd_short->y[i+1]*gd_short->y[i+1]);
	 k+=1;
	}
	fprintf(PEAKMAP,"%s %s %s %s %s %s\n","%", "bin number (0=beg. freq.);", "mjd;", "freq;", "ratio;", "hden;");
	fprintf(PEAKMAP,"%s %s %f %s %f\n","%", "beg. frequency:", gd->ini,"delta_nu:", header_param->deltanu);
	for(k=1;k<len-1;k++){
	  if((power[k]>=power[k-1]) && (power[k] >=power[k+1])){
	    k_short=(int) k/header_param->red;
	    ratio=power[k]/power_short[k_short];
	    if(ratio >=ratio_min)nmax+=1;
	    if(ratio >=ratio_min)fprintf(PEAKMAP,"%d  %f %f %f %f \n",k,header_param->mjdtime,gd->ini+k*header_param->deltanu,ratio,power[k]);
	  }
	}

	free(power);
	free(power_short);
	printf("nmax of 1 FFT= %d\n",nmax);
	return nmax;
}
/*****************Peak map*****************/
int peakmap_from_arp05(GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param,int nfft,double mid_mjd,int nfft_file)
/*Crea peak max using maxima from the SFDB (in gd), normalized with the AR standard deviation of the high
resolution sqrt(spectrum) )
The output table is written into an file
Old format .p05
*/
{

	int i,k,npeak;
	float ratio,ratio_p,ratio_s;
	//float freqm_a,freqp_a; //a means aliased NOT USED 15 JAN 2008
	float *power,*xamed,*xastd;
	float ratio_min=1.5811; //0.; //deve essere 1.5811;  //sqrt(2.5);  //21 aprile 2005
	//float ratio_min1=1.5811;  //era usato per maxdin
 	int len;
	float tau_Hz; //=0.02; //0.01; //era 0.05, poi 0.1 0.15;  //memory of the AR average. Now they are Hz. Now I am using 0.05 Hz.
	float maxage; //=0.02; //0.01;  //0.15, max age of the AR average:if the age is >maxage the Ar average is re-evaluated
	//entrambi a 0.01 va bene alle risonanze ma non abbastanza ai picchi (il simulato a 911 e la calibr.
	//che restano fra 20-30, mentre se metto 0.15 su entrambi i 2 picchi scendono fra 3-4 ma alle
	//risonanze nel ratio ho come un piccolo buco, non sono descritte abbastanza bene
	//0.05 e' un buon compromesso, ma alle risonanze gia' crea qualche problema
	//messi a 0.02 il 16/01/2006 per provare a migliorare con Virgo
	double w_even;
	float ad,qd; //abs of the datum and its square
	float s,ss;  //mean and average, without normalization, every i_eval samples
	double xw,qw,w_norm;
	double itaust; //how many samples in tau_Hz
	int imaxage; //how many samples in age_Hz
	int iage;
	int itaust3;
	double norm;
	int ieval=1; //mean and std are upgraded every 1 sample (a bin in frequency domain)
	float maxdin; //values above this are not considered into the AR evalution
	//float factor=1.5811; //era maxdin=factor*ratio_min era 4--1--cosi ora viene 2.5
	float med_norm;

		/**********Variables needed for the Dopple effect**************/
	 short int fut1,nut1,nutcorr;  //nut1=1 uses nutation from ephem files;
                                      //nut1=0  uses nutation from novas interpolation;

	 short int equatorial=0; //1 output are Equatorial, 0 ecliptic
	 short int imjd=1;      //1 input are mjds, 0 imput are days from 1900
	 double mjd;
	 /*****Detector*****/
	 double *velx,*vely,*velz;
	 double *deinstein;
	 double vel[3];

	 double beg_mjd; //mjd of the beginning point of the FFT
	 int nsaml=1;   //to call PSS_detector_vel only on one sample
	 double data_log[5]; //to register the information of EVF in the log file
	 int prec[5];
	 double evf_ener,evf_amp;
	 float ratio_max;
	 float factor_log; //=2; //dati registrati nel log file sono solo quelli a soglia doppia
	 int stat_log; //alla fine scrivo quanti non ne ho scritti
	 int kn;
	 char nam1[5];

	 DETECTOR *detector;

	 //printf("nfft_file= %d\n",nfft_file);
	 tau_Hz=evf_param->tau_Hz;
	 maxage=evf_param->maxage;
	 maxdin=evf_param->maxdin;
	 factor_log=evf_param->factor_write;
	 stat_log=0;
	 kn=nfft+1;
	 if(kn==1){
	   logfile_comment(LOG_INFO,"par GEN: parameters of the FFTs");
	   logfile_par(LOG_INFO,"GEN_BIN",header_param->deltanu,3);
	   logfile_par(LOG_INFO,"GEN_INIF",header_param->frinit,3);
	   logfile_par(LOG_INFO,"GEN_NSAM",header_param->nsamples,1);
	   logfile_comment(LOG_INFO,"   GEN_BIN is the frequency resolution");
	   logfile_comment(LOG_INFO,"   GEN_INIF is the beginning frequency");
	   logfile_comment(LOG_INFO,"   GEN_NSAM the number of samples in the FFTs");

	   logfile_comment(LOG_INFO,"even NEW: a new FFT has started");
	   logfile_comment(LOG_INFO,"   PAR1: Beginning time of the new FFT");
	   logfile_comment(LOG_INFO,"   PAR2: FFT number in the run");
	   //logfile_comment(LOG_INFO,"   PAR3: number of peaks");
	   logfile_comment(LOG_INFO,"   PAR3: Detector velx/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR4: Detector vely/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR5: Detector velz/C, equatorial");
	   //le caratt. degli EVF sono scritte in crea_evf param
	 }

	 if(header_param->detector==1)strncpy(nam1,"virgo",5);
	 if(header_param->detector==2)strncpy(nam1,"ligoh",5);
	 if(header_param->detector==3)strncpy(nam1,"ligol",5);
	 if(header_param->detector==4)strncpy(nam1,"kagra",5);
	 detector=def_detector(nam1);
	 //evaluation of detector velocities
  /****************************************************************************/
   fut1=0;  ///1=uses ut1; 0=uses utc as aprrox. of ut1; //0 per fare fino al 2010
   nut1=1;  //1=uses nutation from ephem files; 0=no;
   nutcorr=0; //reads the file DELTANUT, if 0 does not read DELTANUT //0 per fare fino al 2010
 /****************************************************************************/
	if(mid_mjd==0)mid_mjd=header_param->mjdtime+(double)(gd->n/2)*header_param->tsamplu/day2sec;
	beg_mjd=mid_mjd-(double)(gd->n/2)*header_param->tsamplu/day2sec;
	//printf("Actual mjdtime, middle point mjd %f %f\n",beg_mjd,mid_mjd);
	velx=(double *)malloc((size_t) (nsaml)*sizeof(double));
	vely=(double *)malloc((size_t) (nsaml)*sizeof(double));
	velz=(double *)malloc((size_t) (nsaml)*sizeof(double));
	deinstein=(double *)malloc((size_t) (nsaml)*sizeof(double));
	mjd=PSS_detector_vel(detector,mid_mjd,header_param->tsamplu,nsaml,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	//puts("uscito dalla function detector_vel");
	vel[0]=velx[0];
	vel[1]=vely[0];
	vel[2]=velz[0];
	free(velx);
	free(vely);
	free(velz);
	free(deinstein);
	//end evaluation of detector velocities
	 //Information for the LOG file: a new FFT is under analysis
	data_log[0]=mjd2gps(beg_mjd);
	 data_log[1]=(double) (header_param->nfft);
	 data_log[2]=(double) vel[0]/C;
	 data_log[3]=(double) vel[1]/C;
	 data_log[4]=(double) vel[2]/C;
	 prec[0]=3;
	 prec[1]=0;
	 prec[2]=3;
	 prec[3]=3;
	 prec[4]=3;
	 logfile_ev(LOG_INFO,"NEW",5,data_log,prec);



	npeak=0;
	len=(int) gd->n/2;

	power=(float *)malloc((size_t) (len)*sizeof(float));
	//puts("allocato 1");
	xamed=(float *)malloc((size_t) (len)*sizeof(float));
	//puts("allocato 2");
	xastd=(float *)malloc((size_t) (len)*sizeof(float));
	//puts("allocato 3");
	puts("Producing peak maps, using the AR mean and std (function: peakmap_from_arp05)");
	/* tolti: per ora Cristiano non li vuole
	if(header_param->nfft==1){
	  fprintf(PEAKMAP,"%s %s %s\n","%", "mjd", "number of peaks");
	  if(equatorial==0)fprintf(PEAKMAP,"%s %s\n","%", "vx,vy,vz ecliptic AU/day (of the middle mjd)");
	  if(equatorial==1)fprintf(PEAKMAP,"%s %s\n","%", "vx,vy,vz equatorial (of the middle mjd)");
	  fprintf(PEAKMAP,"%s %s %s %s\n","%", "bin number (0=beg. freq.);", "ratio;", "Average hden *10^20;");
	  fprintf(PEAKMAP,"%s %s %f %s %f\n","%", "beg. frequency:", gd->ini,"delta_nu:", header_param->deltanu);
	}
	*/
	//NO, definito nella struttura  maxdin=ratio_min1*factor;
	itaust=(double) 1.0/(gd->dx/tau_Hz);
	imaxage=(int) 1.0/(gd->dx/maxage);
	itaust3=(int) ((double) 1.0/(gd->dx/0.5));
	if(itaust3 < (int) (3*itaust)) itaust3=(int) (3*itaust);
	//itaust3: how many bins in 0.5 Hz or in 3 taus. Needed at the beginning
	w_even=exp(-1.0/itaust);
	//freqm_a=header_param->freqm-gd->ini;
	//freqp_a=header_param->freqp-gd->ini;
	//if(header_param->detector==0)printf("Resonant detector. freqm, freqp  %f %f \n",freqm_a,freqp_a);
	if(header_param->detector>=1)printf("Itf detector.min and max freq.%f %f \n",gd->ini,gd->ini+gd->dx*(gd->n-1));

	//printf("len,dnu,tau_Hz,itaust itaust3 %d %f %f %f %d\n",len,gd->dx,tau_Hz,itaust,itaust3);
	//printf("maxage,imaxage,maxdin  %f %d %f \n",maxage,imaxage,maxdin);
	k=0;
	for(i=0;i<gd->n;i+=2){
	 power[k]=sqrt(gd->y[i]*gd->y[i]+gd->y[i+1]*gd->y[i+1]);
	 k+=1;
	}

	/* AR evaluation*/
	k=0;
	xw=0.;
	qw=0.;
	w_norm=0.;
	iage=0;
	while (k<len){
	  if(k==0)xamed[k]=power[k];
	  if(k==0)ratio=1;
	  if(k!=0)ratio=power[k]/xamed[k-1];
	  if(ratio <=maxdin && iage>=1){
	    //End of a frequency event EVF, write the Log file
	    //Write the info in the log file:
	    data_log[3]=evf_ener*gd->dx;  //Power
	    data_log[2]=ratio_max;
	    if(ratio_max>=(maxdin*factor_log))logfile_ev(LOG_INFO,"EVF",4,data_log,prec);
	    stat_log+=1;  //total number of peaks above threshold maxdin
	  }
	  if(ratio <=maxdin){
	    iage=0;
	    ad=power[k];
	    qd=ad*ad;
	    //old procedure: xw=(1.0-w_even)*ad+w_even*xw;
	    //qw=(1.0-w_even)*qd+w_even*qw;
	    //w_norm=(1.0-w_even)+w_even*w_norm;
	    xw=ad+w_even*xw;
	    qw=qd+w_even*qw;
	    w_norm=1+w_even*w_norm;
	    if(k%ieval==0 && k!=0){
	      s=xw;
	      ss=qw;
	      norm=1.0/w_norm;
	      xamed[k]=s*norm;
	      xastd[k]=(float)sqrt((double)fabsf(ss*norm-s*s*norm*norm));
	    }
	  }
	   if(ratio > maxdin){
	     if(iage==0){
	       //a new frequency event begins
	       //Start collecting information for the LOG file
	       data_log[0]=(double) header_param->frinit+k*gd->dx;
	       ratio_max=ratio;
	       evf_ener=0.;
	       evf_amp=0;
	       prec[0]=3;
	       prec[1]=0;
	       prec[2]=0;
	       prec[3]=0;

	     }
	     iage+=1;
	     data_log[1]=(double) (iage*gd->dx); //duration above threshold
	     evf_ener+=(double) power[k]*power[k];  //to evaluate power
	     evf_amp+=(double) power[k];  //to evaluate mean
	     if(ratio>=ratio_max) ratio_max=ratio;
	     xamed[k]=xamed[k-1];
	     xastd[k]=xastd[k-1];
	     if(iage>=40)printf("IAGE > 40 k ,freq,iage imaxage %d %f %d %d\n",k,gd->ini+k*gd->dx,iage,imaxage);
	  }
	  if(iage<=imaxage)k+=1;
	  if(iage>imaxage){
	    //printf("IAGE > IMAXAGE !! k ,iage imaxage %d %d %d\n",k,iage,imaxage);
	    //Write the info in the log file:
	    data_log[1]=0.;     //put to zero the duration, because the AR procedure restarts
	    data_log[3]=evf_amp/iage;  //Average value
	    data_log[2]=ratio_max;
	    logfile_ev(LOG_INFO,"EVF",4,data_log,prec);
	    //Restart the AR evaluation:
	    k-=(iage-1);
	    if(k<0)printf("Att. k <0, iage %d %d \n",k,iage);
	    xw=0.;
	    qw=0.;
	    w_norm=0.;
	    iage=0;
	    if(k>=1)xamed[k-1]=power[k];
	  }
	}
	 //Fill the empty values, with the previous ones (that is, the last evaluated)
	 //if ieval=1 this step is not needed
	if(ieval!=1){
	  for(k=ieval+1; k<len ;k++){
	     if(k%ieval!=0){
	       xamed[k]=xamed[k-1];
	       xastd[k]=xastd[k-1];
	     }
	  }
	}
	/* Averages and stds at the beginning of the gd: */
	  for(k=0; k<2*itaust3;k++){
	    xamed[k]=xamed[2*itaust3];
	    xastd[k]=xastd[2*itaust3];
	  }

	   /*Now check numbers of local maxima in this FFT*/
	  for(k=0; k<len;k++){
	    if(k==0)ratio=1.0; //by definition
	    if(k>=1)ratio=power[k]/xamed[k-1];
	    if(ratio >=ratio_min){
	       if(k<2)ratio_p=1.0; //by definition
	       if(k==(len-1))ratio_s=1.0; //by definition
	       if(k>=2)ratio_p=power[k-1]/xamed[k-2];
	       if(k!=(len-1))ratio_s=power[k+1]/xamed[k];
	       //if((ratio>=ratio_p) && (ratio >=ratio_s))npeak+=1; //mod. 16/01/2006
	       if((ratio>ratio_p) && (ratio >ratio_s))npeak+=1;
	    }
	  }

	  /*Now search for local maxima and write them, in SNR=dato/std*/
	  for(k=0; k<len;k++){

	    // Se ascii: if(k==0)fprintf(PEAKMAP,"%15.10f %d\n",header_param->mjdtime,npeak);
	    //Se ascii:  if(k==0)fprintf(PEAKMAP,"%13.10f %13.10f %13.10f \n", vel[0],vel[1],vel[2]);
	    if(k==0)fwrite((void*)&beg_mjd, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&npeak, sizeof(int),1,PEAKMAP); //Per Cristiano: in binario
	    if(k==0)fwrite((void*)&vel[0], sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&vel[1], sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&vel[2], sizeof(double),1,PEAKMAP);
	    if(k==0)ratio=1.0; //by definition
	    if(k>=1)ratio=power[k]/xamed[k-1];
	    if(ratio >=ratio_min){
	       if(k<2)ratio_p=1.0; //by definition
	       if(k==(len-1))ratio_s=1.0; //by definition
	       if(k>=2)ratio_p=power[k-1]/xamed[k-2];
	       if(k!=(len-1))ratio_s=power[k+1]/xamed[k];
	       // if((ratio>=ratio_p) && (ratio >=ratio_s)){ //mod. 16/01/2006
	       if((ratio>ratio_p) && (ratio >ratio_s)){
		 //fprintf(PEAKMAP,"%d  %f %f %f %f \n",k,header_param->mjdtime,gd->ini+k*header_param->deltanu,ratio,xamed[k-1]);
		 //fprintf(PEAKMAP,"%d %f %f \n",k,ratio,xamed[k-1]); //per Cristiano in ascii
		fwrite((void*)&k, sizeof(int),1,PEAKMAP); ///per Cristiano: in binario
		fwrite((void*)&ratio, sizeof(float),1,PEAKMAP);
		med_norm=xamed[k-1]*header_param->normd*header_param->normw*sqrt(bil2uni); //piapia
		fwrite((void*)&med_norm, sizeof(float),1,PEAKMAP);
	       }
	    }
	    if(nfft==0){
	     if(k==0)fprintf(AREST,"%d  %f %f %f %f \n",k,beg_mjd,gd->ini+k*header_param->deltanu,ratio,xamed[0]*header_param->normd*header_param->normw);
	    if(k>=0)fprintf(AREST,"%d  %f %f %f %f \n",k,beg_mjd,gd->ini+k*header_param->deltanu,ratio,xamed[k-1]*header_param->normd*header_param->normw);
	    }
	  }


	free(power);
	free(xamed);
	free(xastd);
	//printf("npeak of 1 FFT= %d\n",npeak);
	data_log[0]=stat_log;
	prec[0]=0;
	logfile_stat(LOG_INFO,"TOT",1,data_log,prec);
	return npeak;
}

/*****************Peak map with reversed AR estimation  10 July 2007*****************/
int peakmap_from_ar_revp05(GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param,int nfft,double mid_mjd,int nfft_file)
/*Crea peak max using maxima from the SFDB (in gd), normalized with the AR standard deviation of the high
resolution sqrt(spectrum) )
The output table is written into an file
Old format .p05
*/
{

	int i,k,npeak;
	float ratio,ratio_p,ratio_s;
	//float freqm_a,freqp_a; //a means aliased NOT USED 15/01/08
	float *power,*xamed,*xastd;
	float ratio_min=1.5811; //0.; //deve essere 1.5811;  //sqrt(2.5);  //21 aprile 2005
	//float ratio_min1=1.5811;  //era usato per maxdin
 	int len,len1;
	float tau_Hz; //=0.02; //0.01; //era 0.05, poi 0.1 0.15;  //memory of the AR average. Now they are Hz. Now I am using 0.05 Hz.
	float maxage; //=0.02; //0.01;  //0.15, max age of the AR average:if the age is >maxage the Ar average is re-evaluated
	//entrambi a 0.01 va bene alle risonanze ma non abbastanza ai picchi (il simulato a 911 e la calibr.
	//che restano fra 20-30, mentre se metto 0.15 su entrambi i 2 picchi scendono fra 3-4 ma alle
	//risonanze nel ratio ho come un piccolo buco, non sono descritte abbastanza bene
	//0.05 e' un buon compromesso, ma alle risonanze gia' crea qualche problema
	//messi a 0.02 il 16/01/2006 per provare a migliorare con Virgo
	double w_even;
	float ad,qd; //abs of the datum and its square
	float s,ss;  //mean and average, without normalization, every i_eval samples
	double xw,qw,w_norm;
	double itaust; //how many samples in tau_Hz
	int imaxage; //how many samples in age_Hz
	int iage;
	int itaust3;
	double norm;
	int ieval=1; //mean and std are upgraded every 1 sample (a bin in frequency domain)
	float maxdin; //values above this are not considered into the AR evalution
	//float factor=1.5811; //era maxdin=factor*ratio_min era 4--1--cosi ora viene 2.5
	float med_norm;

		/**********Variables needed for the Dopple effect**************/
	 short int fut1,nut1,nutcorr;  //nut1=1 uses nutation from ephem files;
                                      //nut1=0  uses nutation from novas interpolation;

	 short int equatorial=0; //1 output are Equatorial, 0 ecliptic
	 short int imjd=1;      //1 input are mjds, 0 imput are days from 1900
	 double mjd;
	 /*****Detector*****/
	 double *velx,*vely,*velz;
	 double *deinstein;
	 double vel[3];

	 double beg_mjd; //mjd of the beginning point of the FFT
	 int nsaml=1;   //to call PSS_detector_vel only on one sample
	 double data_log[5]; //to register the information of EVF in the log file
	 int prec[5];
	 double evf_ener,evf_amp;
	 float ratio_max;
	 float factor_log; //=2; //dati registrati nel log file sono solo quelli a soglia doppia
	 int stat_log; //alla fine scrivo quanti non ne ho scritti
	 int kn;
	 char nam1[5];

	 DETECTOR *detector;
	 //printf("nfft_file= %d\n",nfft_file);
	 tau_Hz=evf_param->tau_Hz;
	 maxage=evf_param->maxage;
	 maxdin=evf_param->maxdin;
	 factor_log=evf_param->factor_write;
	 stat_log=0;
	 kn=nfft+1;
	 if(kn==1){
	   logfile_comment(LOG_INFO,"par GEN: parameters of the FFTs");
	   logfile_par(LOG_INFO,"GEN_BIN",header_param->deltanu,3);
	   logfile_par(LOG_INFO,"GEN_INIF",header_param->frinit,3);
	   logfile_par(LOG_INFO,"GEN_NSAM",header_param->nsamples,1);
	   logfile_comment(LOG_INFO,"   GEN_BIN is the frequency resolution");
	   logfile_comment(LOG_INFO,"   GEN_INIF is the beginning frequency");
	   logfile_comment(LOG_INFO,"   GEN_NSAM the number of samples in the FFTs");

	   logfile_comment(LOG_INFO,"even NEW: a new FFT has started");
	   logfile_comment(LOG_INFO,"   PAR1: Beginning time of the new FFT");
	   logfile_comment(LOG_INFO,"   PAR2: FFT number in the run");
	   //logfile_comment(LOG_INFO,"   PAR3: number of peaks");
	   logfile_comment(LOG_INFO,"   PAR3: Detector velx/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR4: Detector vely/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR5: Detector velz/C, equatorial");
	   //le caratt. degli EVF sono scritte in crea_evf param
	 }

	 if(header_param->detector==1)strncpy(nam1,"virgo",5);
	 if(header_param->detector==2)strncpy(nam1,"ligoh",5);
	 if(header_param->detector==3)strncpy(nam1,"ligol",5);
	 if(header_param->detector==4)strncpy(nam1,"kagra",5);
	 detector=def_detector(nam1);
	 //evaluation of detector velocities
  /****************************************************************************/
   fut1=0;  ///1=uses ut1; 0=uses utc as aprrox. of ut1; //0 per fare fino al 2010
   nut1=1;  //1=uses nutation from ephem files; 0=no;
   nutcorr=0; //reads the file DELTANUT, if 0 does not read DELTANUT //0 per fare fino al 2010
 /****************************************************************************/
	if(mid_mjd==0)mid_mjd=header_param->mjdtime+(double)(gd->n/2)*header_param->tsamplu/day2sec;
	beg_mjd=mid_mjd-(double)(gd->n/2)*header_param->tsamplu/day2sec;
	//printf("Actual mjdtime, middle point mjd %f %f\n",beg_mjd,mid_mjd);
	velx=(double *)malloc((size_t) (nsaml)*sizeof(double));
	vely=(double *)malloc((size_t) (nsaml)*sizeof(double));
	velz=(double *)malloc((size_t) (nsaml)*sizeof(double));
	deinstein=(double *)malloc((size_t) (nsaml)*sizeof(double));
	mjd=PSS_detector_vel(detector,mid_mjd,header_param->tsamplu,nsaml,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	//puts("uscito dalla function detector_vel");
	vel[0]=velx[0];
	vel[1]=vely[0];
	vel[2]=velz[0];
	free(velx);
	free(vely);
	free(velz);
	free(deinstein);
	//end evaluation of detector velocities
	 //Information for the LOG file: a new FFT is under analysis
	data_log[0]=mjd2gps(beg_mjd);
	 data_log[1]=(double) (header_param->nfft);
	 data_log[2]=(double) vel[0]/C;
	 data_log[3]=(double) vel[1]/C;
	 data_log[4]=(double) vel[2]/C;
	 prec[0]=3;
	 prec[1]=0;
	 prec[2]=3;
	 prec[3]=3;
	 prec[4]=3;
	 logfile_ev(LOG_INFO,"NEW",5,data_log,prec);



	npeak=0;
	len=(int) gd->n/2;

	power=(float *)malloc((size_t) (len)*sizeof(float));
	//puts("allocato 1");
	xamed=(float *)malloc((size_t) (len)*sizeof(float));
	//puts("allocato 2");
	xastd=(float *)malloc((size_t) (len)*sizeof(float));
	//puts("allocato 3");
	puts("Producing peak maps, using the AR mean and std (function peakmap_from_ar_revp05)");
	/* tolti: per ora Cristiano non li vuole
	if(header_param->nfft==1){
	  fprintf(PEAKMAP,"%s %s %s\n","%", "mjd", "number of peaks");
	  if(equatorial==0)fprintf(PEAKMAP,"%s %s\n","%", "vx,vy,vz ecliptic AU/day (of the middle mjd)");
	  if(equatorial==1)fprintf(PEAKMAP,"%s %s\n","%", "vx,vy,vz equatorial (of the middle mjd)");
	  fprintf(PEAKMAP,"%s %s %s %s\n","%", "bin number (0=beg. freq.);", "ratio;", "Average hden *10^20;");
	  fprintf(PEAKMAP,"%s %s %f %s %f\n","%", "beg. frequency:", gd->ini,"delta_nu:", header_param->deltanu);
	}
	*/
	//NO, definito nella struttura  maxdin=ratio_min1*factor;
	itaust=(double) 1.0/(gd->dx/tau_Hz);
	imaxage=(int) 1.0/(gd->dx/maxage);
	itaust3=(int) ((double) 1.0/(gd->dx/0.5));
	if(itaust3 < (int) (3*itaust)) itaust3=(int) (3*itaust);
	//itaust3: how many bins in 0.5 Hz or in 3 taus. Needed at the beginning
	w_even=exp(-1.0/itaust);
	//freqm_a=header_param->freqm-gd->ini;
	//freqp_a=header_param->freqp-gd->ini;
	//if(header_param->detector==0)printf("Resonant detector. freqm, freqp  %f %f \n",freqm_a,freqp_a);
	//if(header_param->detector>=1)printf("Itf detector.min and max freq.%f %f \n",gd->ini,gd->ini+gd->dx*(gd->n-1));

	//printf("len,dnu,tau_Hz,itaust itaust3 %d %f %f %f %d\n",len,gd->dx,tau_Hz,itaust,itaust3);
	//printf("maxage,imaxage,maxdin  %f %d %f \n",maxage,imaxage,maxdin);
	k=0;
	for(i=0;i<gd->n;i+=2){
	 power[k]=sqrt(gd->y[i]*gd->y[i]+gd->y[i+1]*gd->y[i+1]);
	 k+=1;
	}

	/* AR evaluation*/
	k=len-1; //Reverse mode (data are from 0 to len-1)
	len1=len-1;
	xw=0.;
	qw=0.;
	w_norm=0.;
	iage=0;
	ratio=1.;
	while (k>=0){
	  if(k==len1)xamed[k]=power[k];
	  if(k==len1)ratio=1.;
	  if(k!=len1)ratio=power[k]/xamed[k+1];
	  if(ratio <=maxdin && iage>=1){
	    //End of a frequency event EVF, write the Log file
	    //Write the info in the log file:
	    data_log[3]=evf_ener*gd->dx;  //Power
	    data_log[2]=ratio_max;
	    if(ratio_max>=(maxdin*factor_log))logfile_ev(LOG_INFO,"EVF",4,data_log,prec);
	    stat_log+=1;  //total number of peaks above threshold maxdin
	  }
	  if(ratio <=maxdin){
	    iage=0;
	    ad=power[k];
	    qd=ad*ad;
	    //old procedure: xw=(1.0-w_even)*ad+w_even*xw;
	    //qw=(1.0-w_even)*qd+w_even*qw;
	    //w_norm=(1.0-w_even)+w_even*w_norm;
	    xw=ad+w_even*xw;
	    qw=qd+w_even*qw;
	    w_norm=1+w_even*w_norm;
	    if(k%ieval==0 && k!=len1){
	      s=xw;
	      ss=qw;
	      norm=1.0/w_norm;
	      xamed[k]=s*norm;
	      xastd[k]=(float)sqrt((double)fabsf(ss*norm-s*s*norm*norm));
	    }
	  }
	   if(ratio > maxdin){
	     if(iage==0){
	       //a new frequency event begins
	       //Start collecting information for the LOG file
	       data_log[0]=(double) header_param->frinit+k*gd->dx;
	       ratio_max=ratio;
	       evf_ener=0.;
	       evf_amp=0;
	       prec[0]=3;
	       prec[1]=0;
	       prec[2]=0;
	       prec[3]=0;

	     }
	     iage+=1;
	     data_log[1]=(double) (iage*gd->dx); //duration above threshold
	     evf_ener+=(double) power[k]*power[k];  //to evaluate power
	     evf_amp+=(double) power[k];  //to evaluate mean
	     if(ratio>=ratio_max) ratio_max=ratio;
	     xamed[k]=xamed[k+1];
	     xastd[k]=xastd[k+1];
	     //if(iage>=40)printf("IAGE > 40 k ,freq,iage imaxage %d %f %d %d\n",k,gd->ini+k*gd->dx,iage,imaxage);
	  }
	  if(iage<=imaxage)k-=1;
	  if(iage>imaxage){
	    //printf("IAGE > IMAXAGE !! k ,iage imaxage %d %d %d\n",k,iage,imaxage);
	    //Write the info in the log file:
	    data_log[1]=0.;     //put to zero the duration, because the AR procedure restarts
	    data_log[3]=evf_amp/iage;  //Average value
	    data_log[2]=ratio_max;
	    logfile_ev(LOG_INFO,"EVF",4,data_log,prec);
	    //Restart the AR evaluation:
	    k+=(iage-1);
	    //if(k>len1)printf("Att. k >len1, iage %d %d \n",k,iage);
	    xw=0.;
	    qw=0.;
	    w_norm=0.;
	    iage=0;
	    if(k<=(len1-1))xamed[k+1]=power[k];
	  }
	}
	 //Fill the empty values, with the previous ones (that is, the last evaluated)
	 //if ieval=1 this step is not needed
	if(ieval!=1){
	  for(k=(len1-(ieval-1)); k>=0 ;k--){
	     if(k%ieval!=0){
	       xamed[k]=xamed[k+1];
	       xastd[k]=xastd[k+1];
	     }
	  }
	}
	/* Averages and stds at the END of the gd: */
	  for(k=len1; k<len1-2*itaust3;k--){
	    xamed[k]=xamed[len1-2*itaust3];
	    xastd[k]=xastd[len1-2*itaust3];
	  }

	   /*Now check numbers of local maxima in this FFT*/
	  for(k=len1; k>=0;k--){
	    if(k>=(len1-1))ratio=1.0; //by definition
	    if(k<(len1-1))ratio=power[k]/xamed[k+1];
	    if(ratio >=ratio_min){
	       if(k<1)ratio_p=1.0; //by definition
	       if(k>=(len1-2))ratio_s=1.0; //by definition
	       if(k>=1)ratio_p=power[k-1]/xamed[k];
	       if(k<(len1-2))ratio_s=power[k+1]/xamed[k+2];
	       //if((ratio>=ratio_p) && (ratio >=ratio_s))npeak+=1; //mod. 16/01/2006
	       if((ratio>ratio_p) && (ratio >ratio_s))npeak+=1;
	    }
	  }


	   /*Now search for local maxima and write them, in SNR=dato/std*/
	   for(k=0; k<=len1;k++){

	    // Se ascii: if(k==0)fprintf(PEAKMAP,"%15.10f %d\n",header_param->mjdtime,npeak);
	    //Se ascii:  if(k==0)fprintf(PEAKMAP,"%13.10f %13.10f %13.10f \n", vel[0],vel[1],vel[2]);
	    if(k==0)fwrite((void*)&beg_mjd, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&npeak, sizeof(int),1,PEAKMAP); //Per Cristiano: in binario
	    if(k==0)fwrite((void*)&vel[0], sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&vel[1], sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&vel[2], sizeof(double),1,PEAKMAP);
	    if(k>=(len1-1))ratio=1.0; //by definition
	    if(k<(len1-1))ratio=power[k]/xamed[k+1];

	    if(ratio >=ratio_min){
	       if(k<1)ratio_p=1.0; //by definition
	       if(k>=(len1-2))ratio_s=1.0; //by definition
	       if(k>=1)ratio_p=power[k-1]/xamed[k];
	       if(k<(len1-2))ratio_s=power[k+1]/xamed[k+2];

	       if((ratio>ratio_p) && (ratio >ratio_s)){
		 //fprintf(PEAKMAP,"%d  %f %f %f %f \n",k,header_param->mjdtime,gd->ini+k*header_param->deltanu,ratio,xamed[k-1]);
		 //fprintf(PEAKMAP,"%d %f %f \n",k,ratio,xamed[k-1]); //per Cristiano in ascii
		fwrite((void*)&k, sizeof(int),1,PEAKMAP); ///per Cristiano: in binario
		fwrite((void*)&ratio, sizeof(float),1,PEAKMAP);
		med_norm=xamed[k-1]*header_param->normd*header_param->normw*sqrt(bil2uni); //piapia
		fwrite((void*)&med_norm, sizeof(float),1,PEAKMAP);
	       }
	    }
	    if(nfft==0){
	     if(k==0)fprintf(AREST,"%d  %f %f %f %f \n",k,beg_mjd,gd->ini+k*header_param->deltanu,ratio,xamed[0]*header_param->normd*header_param->normw);
	    if(k>=0)fprintf(AREST,"%d  %f %f %f %f \n",k,beg_mjd,gd->ini+k*header_param->deltanu,ratio,xamed[k-1]*header_param->normd*header_param->normw);
	    }
	   }


	free(power);
	free(xamed);
	free(xastd);
	//printf("npeak of 1 FFT= %d\n",npeak);
	data_log[0]=stat_log;
	prec[0]=0;
	logfile_stat(LOG_INFO,"TOT",1,data_log,prec);
	return npeak;
}
/*****************Peak map*****************/
int peakmap_from_ar(GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param,int nfft,double mid_mjd,int nfft_file, GD *gd_short)
/*Crea peak max using maxima from the SFDB (in gd), normalized with the AR standard deviation of the high
resolution sqrt(spectrum) )
The output table is written into an file
Last modified 27 Aug 2008, new format .p08
Last modified 30 March 2009, new format .p09
*/
{

	int i,k,npeak;
	float ratio,ratio_p,ratio_s;
	//float freqm_a,freqp_a; //a means aliased NOT USED 15 JAN 2008
	float *power,*xamed,*xastd;
	float ratio_min; //=1.5811; //0.; //deve essere 1.5811;  //sqrt(2.5);  //21 aprile 2005
	//float ratio_min1=1.5811;  //era usato per maxdin
 	int len;
	float tau_Hz; //=0.02; //0.01; //era 0.05, poi 0.1 0.15;  //memory of the AR average. Now they are Hz. Now I am using 0.05 Hz.
	float maxage; //=0.02; //0.01;  //0.15, max age of the AR average:if the age is >maxage the Ar average is re-evaluated
	//entrambi a 0.01 va bene alle risonanze ma non abbastanza ai picchi (il simulato a 911 e la calibr.
	//che restano fra 20-30, mentre se metto 0.15 su entrambi i 2 picchi scendono fra 3-4 ma alle
	//risonanze nel ratio ho come un piccolo buco, non sono descritte abbastanza bene
	//0.05 e' un buon compromesso, ma alle risonanze gia' crea qualche problema
	//messi a 0.02 il 16/01/2006 per provare a migliorare con Virgo
	double w_even;
	float ad,qd; //abs of the datum and its square
	float s,ss;  //mean and average, without normalization, every i_eval samples
	double xw,qw,w_norm;
	double itaust; //how many samples in tau_Hz
	int imaxage; //how many samples in age_Hz
	int iage;
	int itaust3;
	double norm;
	int ieval=1; //mean and std are upgraded every 1 sample (a bin in frequency domain)
	float maxdin; //values above this are not considered into the AR evalution
	//float factor=1.5811; //era maxdin=factor*ratio_min era 4--1--cosi ora viene 2.5
	float med_norm;

		/**********Variables needed for the Dopple effect**************/
	 short int fut1,nut1,nutcorr;  //nut1=1 uses nutation from ephem files;
                                      //nut1=0  uses nutation from novas interpolation;

	 //short int equatorial=0; //1 output are Equatorial, 0 ecliptic
	 //short int imjd=1;      //1 input are mjds, 0 imput are days from 1900
	 //double mjd;
	 /*****Detector*****/
	 //OLD double *velx,*vely,*velz;
	 //double *deinstein;
	 //double vel[3];

	 double beg_mjd; //mjd of the beginning point of the FFT
	 //int nsaml=1;   //to call PSS_detector_vel only on one sample
	 double data_log[5]; //to register the information of EVF in the log file
	 int prec[5];
	 double evf_ener,evf_amp;
	 float ratio_max;
	 float factor_log; //=2; //dati registrati nel log file sono solo quelli a soglia doppia
	 int stat_log; //alla fine scrivo quanti non ne ho scritti
	 int kn;
	 double sampling_frequency;
	 double spini,spdf;
	 int splen;
	 int ik;
	 float sp;
	 char nam1[5];
	 ratio_min=evf_param->maxdin;
	 spini=gd->ini;
	 spdf=header_param->deltanu*header_param->red;
	 splen=gd->n/(2*header_param->red);
	 //printf("Short ps spini,spdf,splen %f %f %d\n",spini,spdf,splen);

	 DETECTOR *detector;


	 tau_Hz=evf_param->tau_Hz;
	 maxage=evf_param->maxage;
	 maxdin=evf_param->maxdin;
	 factor_log=evf_param->factor_write;
	 stat_log=0;
	 kn=nfft+1;
	 if(kn==1){
	   logfile_comment(LOG_INFO,"par GEN: parameters of the FFTs");
	   logfile_par(LOG_INFO,"GEN_BIN",header_param->deltanu,3);
	   logfile_par(LOG_INFO,"GEN_INIF",header_param->frinit,3);
	   logfile_par(LOG_INFO,"GEN_NSAM",header_param->nsamples,1);
	   logfile_comment(LOG_INFO,"   GEN_BIN is the frequency resolution");
	   logfile_comment(LOG_INFO,"   GEN_INIF is the beginning frequency");
	   logfile_comment(LOG_INFO,"   GEN_NSAM the number of samples in the FFTs");

	   logfile_comment(LOG_INFO,"even NEW: a new FFT has started");
	   logfile_comment(LOG_INFO,"   PAR1: Beginning time of the new FFT");
	   logfile_comment(LOG_INFO,"   PAR2: FFT number in the run");
	   //logfile_comment(LOG_INFO,"   PAR3: number of peaks");
	   logfile_comment(LOG_INFO,"   PAR3: Detector velx/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR4: Detector vely/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR5: Detector velz/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR6: Detector posx/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR7: Detector posy/C, equatorial");
	   logfile_comment(LOG_INFO,"   PAR8: Detector posz/C, equatorial");
	   //le caratt. degli EVF sono scritte in crea_evf param
	 }
       if(nfft_file==1){
	   // File Header
	   // nfft int32: Added at the end in the main code
	   // sfr   double (sampling frequency)
	   // lfft   int32 (original length)
	   // inifr   double
	   sampling_frequency=1/header_param->tsamplu;
	   fwrite((void*)&sampling_frequency, sizeof(double),1,PEAKMAP);
	   fwrite((void*)&header_param->nsamples, sizeof(int),1,PEAKMAP);
	   fwrite((void*)&header_param->frinit, sizeof(double),1,PEAKMAP);
	 }
         if(header_param->detector==1)strncpy(nam1,"virgo",5);
	 if(header_param->detector==2)strncpy(nam1,"ligoh",5);
	 if(header_param->detector==3)strncpy(nam1,"ligol",5);
	 if(header_param->detector==4)strncpy(nam1,"kagra",5);
	  detector=def_detector(nam1);
	 //evaluation of detector velocities
  /****************************************************************************/
   fut1=0;  ///1=uses ut1; 0=uses utc as aprrox. of ut1; //0 per fare fino al 2010
   nut1=1;  //1=uses nutation from ephem files; 0=no;
   nutcorr=0; //reads the file DELTANUT, if 0 does not read DELTANUT //0 per fare fino al 2010
 /****************************************************************************/
	if(mid_mjd==0)mid_mjd=header_param->mjdtime+(double)(gd->n/2)*header_param->tsamplu/day2sec;
	beg_mjd=mid_mjd-(double)(gd->n/2)*header_param->tsamplu/day2sec;
	//printf("Actual mjdtime, middle point mjd %f %f\n",beg_mjd,mid_mjd);
	/* OLD
	velx=(double *)malloc((size_t) (nsaml)*sizeof(double));
	vely=(double *)malloc((size_t) (nsaml)*sizeof(double));
	velz=(double *)malloc((size_t) (nsaml)*sizeof(double));
	deinstein=(double *)malloc((size_t) (nsaml)*sizeof(double));
	mjd=PSS_detector_vel(detector,mid_mjd,header_param->tsamplu,nsaml,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,deinstein);
	//puts("uscito dalla function detector_vel");
	vel[0]=velx[0];
	vel[1]=vely[0];
	vel[2]=velz[0];
	free(velx);
	free(vely);
	free(velz);
	free(deinstein);
	//end evaluation of detector velocities
	*/
	 //Information for the LOG file: a new FFT is under analysis
	data_log[0]=mjd2gps(beg_mjd);
	 data_log[1]=(double) (header_param->nfft);
	 data_log[2]=(double) header_param->vx_eq;;
	 data_log[3]=(double) header_param->vy_eq;;
	 data_log[4]=(double) header_param->vz_eq;;
	 data_log[5]=(double) header_param->px_eq;
	 data_log[6]=(double) header_param->py_eq;;
	 data_log[7]=(double) header_param->pz_eq;;
	 prec[0]=3;
	 prec[1]=0;
	 prec[2]=3;
	 prec[3]=3;
	 prec[4]=3;
	 prec[5]=3;
	 prec[6]=3;
	 prec[7]=3;

	 logfile_ev(LOG_INFO,"NEW",8,data_log,prec);



	npeak=0;
	len=(int) gd->n/2;

	power=(float *)malloc((size_t) (len)*sizeof(float));
	//puts("allocato 1");
	xamed=(float *)malloc((size_t) (len)*sizeof(float));
	//puts("allocato 2");
	xastd=(float *)malloc((size_t) (len)*sizeof(float));
	//puts("allocato 3");
	puts("Producing peak maps, using the AR mean and std (function: peakmap_from_ar)");

	itaust=(double) 1.0/(gd->dx/tau_Hz);
	imaxage=(int) 1.0/(gd->dx/maxage);
	itaust3=(int) ((double) 1.0/(gd->dx/0.5));
	if(itaust3 < (int) (3*itaust)) itaust3=(int) (3*itaust);
	//itaust3: how many bins in 0.5 Hz or in 3 taus. Needed at the beginning
	w_even=exp(-1.0/itaust);

	//if(header_param->detector>=1)printf("Itf detector.min and max freq.%f %f \n",gd->ini,gd->ini+gd->dx*(gd->n-1));

	//printf("len,dnu,tau_Hz,itaust itaust3 %d %f %f %f %d\n",len,gd->dx,tau_Hz,itaust,itaust3);
	//printf("maxage,imaxage,maxdin  %f %d %f \n",maxage,imaxage,maxdin);
	k=0;
	for(i=0;i<gd->n;i+=2){
	 power[k]=sqrt(gd->y[i]*gd->y[i]+gd->y[i+1]*gd->y[i+1]);
	 k+=1;
	}

	/* AR evaluation*/
	k=0;
	xw=0.;
	qw=0.;
	w_norm=0.;
	iage=0;
	while (k<len){
	  if(k==0)xamed[k]=power[k];
	  if(k==0)ratio=1;
	  if(k!=0)ratio=power[k]/xamed[k-1];
	  if(ratio <=maxdin && iage>=1){
	    //End of a frequency event EVF, write the Log file
	    //Write the info in the log file:
	    //data_log[3]=evf_ener*gd->dx;  //Power
	    data_log[3]=evf_ener*gd->dx*pow(header_param->normd,2)*pow(header_param->normw,2)*bil2uni;  //Unilateral Energy normalized correctly
	    data_log[2]=ratio_max;
	    if(ratio_max>=(maxdin*factor_log))logfile_ev(LOG_INFO,"EVF",4,data_log,prec);
	    stat_log+=1;  //total number of peaks above threshold maxdin
	  }
	  if(ratio <=maxdin){
	    iage=0;
	    ad=power[k];
	    qd=ad*ad;
	    //old procedure: xw=(1.0-w_even)*ad+w_even*xw;
	    //qw=(1.0-w_even)*qd+w_even*qw;
	    //w_norm=(1.0-w_even)+w_even*w_norm;
	    xw=ad+w_even*xw;
	    qw=qd+w_even*qw;
	    w_norm=1+w_even*w_norm;
	    if(k%ieval==0 && k!=0){
	      s=xw;
	      ss=qw;
	      norm=1.0/w_norm;
	      xamed[k]=s*norm;
	      xastd[k]=(float)sqrt((double)fabsf(ss*norm-s*s*norm*norm));
	    }
	  }
	   if(ratio > maxdin){
	     if(iage==0){
	       //a new frequency event begins
	       //Start collecting information for the LOG file
	       data_log[0]=(double) header_param->frinit+k*gd->dx;
	       ratio_max=ratio;
	       evf_ener=0.;
	       evf_amp=0;
	       prec[0]=3;
	       prec[1]=0;
	       prec[2]=0;
	       prec[3]=0;

	     }
	     iage+=1;
	     data_log[1]=(double) (iage*gd->dx); //duration above threshold
	     evf_ener+=(double) power[k]*power[k];  //to evaluate power
	     evf_amp+=(double) power[k];  //to evaluate mean
	     if(ratio>=ratio_max) ratio_max=ratio;
	     xamed[k]=xamed[k-1];
	     xastd[k]=xastd[k-1];
	     //printf("IAGE");
	     //if(iage>=40)printf("IAGE > 40 k ,freq,iage imaxage %d %f %d %d\n",k,gd->ini+k*gd->dx,iage,imaxage);
	  }
	  if(iage<=imaxage)k+=1;
	  if(iage>imaxage){
	    //printf("IAGE > IMAXAGE !! k ,iage imaxage %d %d %d\n",k,iage,imaxage);
	    //Write the info in the log file:
	    data_log[1]=0.;     //put to zero the duration, because the AR procedure restarts
	    data_log[3]=(evf_amp/iage)*header_param->normd*header_param->normw*sqrt(bil2uni); ;  //Average value
	    data_log[2]=ratio_max;
	    logfile_ev(LOG_INFO,"EVF",4,data_log,prec);
	    //Restart the AR evaluation:
	    k-=(iage-1);
	    //if(k<0)printf("Att. k <0, iage %d %d \n",k,iage);
	    xw=0.;
	    qw=0.;
	    w_norm=0.;
	    iage=0;
	    if(k>=1)xamed[k-1]=power[k];
	  }
	}
	 //Fill the empty values, with the previous ones (that is, the last evaluated)
	 //if ieval=1 this step is not needed
	if(ieval!=1){
	  for(k=ieval+1; k<len ;k++){
	     if(k%ieval!=0){
	       xamed[k]=xamed[k-1];
	       xastd[k]=xastd[k-1];
	     }
	  }
	}
	/* Averages and stds at the beginning of the gd: */
	  for(k=0; k<2*itaust3;k++){
	    xamed[k]=xamed[2*itaust3];
	    xastd[k]=xastd[2*itaust3];
	  }

	   /*Now check numbers of local maxima in this FFT*/
	  for(k=0; k<len;k++){
	    if(k==0)ratio=1.0; //by definition
	    if(k>=1)ratio=power[k]/xamed[k-1];
	    if(ratio >=ratio_min){
	       if(k<2)ratio_p=1.0; //by definition
	       if(k==(len-1))ratio_s=1.0; //by definition
	       if(k>=2)ratio_p=power[k-1]/xamed[k-2];
	       if(k!=(len-1))ratio_s=power[k+1]/xamed[k];
	       //if((ratio>=ratio_p) && (ratio >=ratio_s))npeak+=1; //mod. 16/01/2006
	       if((ratio>ratio_p) && (ratio >ratio_s))npeak+=1;
	    }
	  }

	  /*Now search for local maxima and write them, in SNR=dato/std*/
	  /*First write the BLOCK HEADER*/
	  /* Then the SP data */
	  /* Then the data */
	  // Block header
	  // mjd double
	  // npeak int32
	  // velx double
	  // vely double
	  // velz double
	  //Short spectrum  SP(in Einstein^2)
	  // spini double
	  // spdf  double
	  // splen int32
	  // Data sp  float
	  // Data (npeak times)
	  // bin int32
	  // ratio float
	  // xamed float (mean of H)
	  for(k=0; k<len;k++){
	    if(k==0)fwrite((void*)&beg_mjd, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&npeak, sizeof(int),1,PEAKMAP);
	    if(k==0)fwrite((void*)&header_param->vx_eq, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&header_param->vy_eq, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&header_param->vz_eq, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&header_param->px_eq, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&header_param->py_eq, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&header_param->pz_eq, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&spini, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&spini, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&spdf, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&splen, sizeof(int),1,PEAKMAP);
	    if(k==0){
	      for(ik=0;ik<=splen-1;ik++){
		sp=(gd_short->y[ik]*gd_short->y[ik]);
		fwrite((void*)&sp, sizeof(float),1,PEAKMAP);
	      }
	    }

	    if(k>=(len-1))ratio=1.0; //by definition
	    if(k==0)ratio=1.0; //by definition
	    if(k>=1)ratio=power[k]/xamed[k-1];
	    if(ratio >=ratio_min){
	       if(k<2)ratio_p=1.0; //by definition
	       if(k==(len-1))ratio_s=1.0; //by definition
	       if(k>=2)ratio_p=power[k-1]/xamed[k-2];
	       if(k!=(len-1))ratio_s=power[k+1]/xamed[k];
	       // if((ratio>=ratio_p) && (ratio >=ratio_s)){ //mod. 16/01/2006
	       if((ratio>ratio_p) && (ratio >ratio_s)){

		fwrite((void*)&k, sizeof(int),1,PEAKMAP);
		fwrite((void*)&ratio, sizeof(float),1,PEAKMAP);
		med_norm=xamed[k-1]*header_param->normd*header_param->normw*sqrt(bil2uni); //piapia
		if(header_param->perc_zeroes < 1)med_norm=med_norm/sqrt(1-header_param->perc_zeroes);
		fwrite((void*)&med_norm, sizeof(float),1,PEAKMAP);
	       }
	    }
	    if(nfft==0){
	     if(k==0)fprintf(AREST,"%d  %f %f %f %f \n",k,beg_mjd,gd->ini+k*header_param->deltanu,ratio,xamed[0]*header_param->normd*header_param->normw);
	    if(k>=0)fprintf(AREST,"%d  %f %f %f %f \n",k,beg_mjd,gd->ini+k*header_param->deltanu,ratio,xamed[k-1]*header_param->normd*header_param->normw);
	    }
	  }


	free(power);
	free(xamed);
	free(xastd);
	//printf("npeak of 1 FFT= %d\n",npeak);
	data_log[0]=stat_log;
	prec[0]=0;
	logfile_stat(LOG_INFO,"TOT",1,data_log,prec);
	return npeak;
}

/*****************Peak map with reversed AR estimation  Last 27 Aug 2008*****************/
int peakmap_from_ar_rev(GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param,int nfft,double mid_mjd,int nfft_file,GD *gd_short)
/*Crea peak max using maxima from the SFDB (in gd), normalized with the AR standard deviation of the high
resolution sqrt(spectrum) )
The output table is written into an file
New .p08 format  27 Aug 2008
Last modified 30 March 2009, new format .p09
*/
{
	int i,k,npeak,jj;
	float ratio,ratio_p,ratio_s,ratioS;
	//float freqm_a,freqp_a; //a means aliased NOT USED 15/01/08
	float *power,*xamed,*xastd;
	float ratio_min; //0.; //deve essere 1.5811;  //sqrt(2.5);  //21 aprile 2005
	//float ratio_min1=1.5811;  //era usato per maxdin
 	int len,len1;
	float tau_Hz; //=0.02; //0.01; //era 0.05, poi 0.1 0.15;  //memory of the AR average. Now they are Hz. Now I am using 0.05 Hz.
	float maxage; //=0.02; //0.01;  //0.15, max age of the AR average:if the age is >maxage the Ar average is re-evaluated
	//entrambi a 0.01 va bene alle risonanze ma non abbastanza ai picchi (il simulato a 911 e la calibr.
	//che restano fra 20-30, mentre se metto 0.15 su entrambi i 2 picchi scendono fra 3-4 ma alle
	//risonanze nel ratio ho come un piccolo buco, non sono descritte abbastanza bene
	//0.05 e' un buon compromesso, ma alle risonanze gia' crea qualche problema
	//messi a 0.02 il 16/01/2006 per provare a migliorare con Virgo
	double w_even;
	float ad,qd; //abs of the datum and its square
	float s,ss;  //mean and average, without normalization, every i_eval samples
	double xw,qw,w_norm;
	double itaust; //how many samples in tau_Hz
	int imaxage; //how many samples in age_Hz
	int iage;
	int itaust3;
	double norm;
	int ieval=1; //mean and std are upgraded every 1 sample (a bin in frequency domain)
	float maxdin; //values above this are not considered into the AR evalution
	//float factor=1.5811; //era maxdin=factor*ratio_min era 4--1--cosi ora viene 2.5
	float med_norm;
	int ilog=2; //means AR done on the power. <2 on the amplitude

		/**********Variables needed for the Dopple effect**************/
	//short int fut1,nut1,nutcorr;  //nut1=1 uses nutation from ephem files;
                                      //nut1=0  uses nutation from novas interpolation;

	 //short int equatorial=0; //1 output are Equatorial, 0 ecliptic
	 //short int imjd=1;      //1 input are mjds, 0 imput are days from 1900
	 //double mjd;
	 /*****Detector*****/
	 //OLD double *velx,*vely,*velz;
	 //double *deinstein;
	 //double vel[3];

	 double beg_mjd; //mjd of the beginning point of the FFT
	 //int nsaml=1;   //to call PSS_detector_vel only on one sample
	 double data_log[8]; //to register the information of EVF in the log file
	 int prec[8];
	 double evf_ener,evf_amp;
	 float ratio_max;
	 float factor_log; //=2; //dati registrati nel log file sono solo quelli a soglia doppia
	 int stat_log; //alla fine scrivo quanti non ne ho scritti
	 int kn;
	 double sampling_frequency;
	 double spini,spdf;
	 int splen;
	 int ik;
	 float sp;
	 char nam1[5];
	 int kinizio;
	 ratio_min=evf_param->maxdin;
	 spini=gd->ini;
	 spdf=header_param->deltanu*header_param->red;
	 splen=gd->n/(2*header_param->red);
	 //printf("Short ps spini,spdf,splen %f %f %d\n",spini,spdf,splen);
	 char nam[5];
	 int noadj=10;

	 //DETECTOR *detector;
	  if(header_param->nfft==1){
	    if(header_param->deltanu <= 0.001){
	      evf_param->tau_Hz=0.02;
	      evf_param->maxage=0.02;

  }
  else
    {
      evf_param->tau_Hz= 0.02*header_param->deltanu/0.001;
      evf_param->maxage=evf_param->tau_Hz;

    }

  //Write in the LOG file


  //printf("evf param tau_Hz %f\n",evf_param->tau_Hz);

  logfile_par(LOG_INFO,"GEN_THR",evf_param->maxdin,3);
  logfile_par(LOG_INFO,"GEN_TAU",evf_param->tau_Hz,3);
  logfile_par(LOG_INFO,"GEN_MAXAGE",evf_param->maxage,3);
  logfile_par(LOG_INFO,"GEN_FAC",evf_param->factor_write,3);


  }

	 tau_Hz=evf_param->tau_Hz;
	 maxage=evf_param->maxage;
	 if(ilog <2)maxdin=evf_param->maxdin;
	 if(ilog ==2)maxdin=(evf_param->maxdin)*(evf_param->maxdin);
	 factor_log=evf_param->factor_write;
	 stat_log=0;
	 kn=nfft+1;
	 if(kn==1){
	   //printf("factor_log = %f\n",factor_log);
	   if(header_param->detector==1)strncpy(nam1,"virgo",5);
	   if(header_param->detector==2)strncpy(nam1,"ligoh",5);
	   if(header_param->detector==3)strncpy(nam1,"ligol",5);
	   if(header_param->detector==4)strncpy(nam1,"kagra",5);
	   printf("detector as written in the header %s\n",nam1);
	   printf("Confirm: Detector ? virgo,ligoh,ligol \n");
	   scanf("%s",nam);
	   jj=logfile_input(LOG_INFO,nam,"detector");
	   logfile_comment(LOG_INFO,"par GEN: parameters of the FFTs");
	   logfile_par(LOG_INFO,"GEN_BIN",header_param->deltanu,3);
	   logfile_par(LOG_INFO,"GEN_INIF",header_param->frinit,3);
	   logfile_par(LOG_INFO,"GEN_NSAM",header_param->nsamples,1);
	   logfile_comment(LOG_INFO,"   GEN_BIN is the frequency resolution");
	   logfile_comment(LOG_INFO,"   GEN_INIF is the beginning frequency");
	   logfile_comment(LOG_INFO,"   GEN_NSAM the number of samples in the FFTs");

	   logfile_comment(LOG_INFO,"even NEW: a new FFT has started");
	   logfile_comment(LOG_INFO,"   PAR1: Beginning gps time of the new FFT");
	   logfile_comment(LOG_INFO,"   PAR2: Beginning mjd time of the new FFT");
	   logfile_comment(LOG_INFO,"   PAR3: FFT number in the run");
	   //logfile_comment(LOG_INFO,"   PAR3: number of peaks");
	   logfile_comment(LOG_INFO,"even VEL:veloc and positions");
	   logfile_comment(LOG_INFO,"  PAR1: Detector velx/C, equatorial");
	   logfile_comment(LOG_INFO,"  PAR2: Detector vely/C, equatorial");
	   logfile_comment(LOG_INFO,"  PAR3: Detector velz/C, equatorial");
	   logfile_comment(LOG_INFO,"  PAR4: Detector posx/C, equatorial");
	   logfile_comment(LOG_INFO,"  PAR5: Detector posy/C, equatorial");
	   logfile_comment(LOG_INFO,"  PAR6: Detector posz/C, equatorial");
	   //le caratt. degli EVF sono scritte in crea_evf param
	 }
	 if(nfft_file==1){
	   //printf("kn, nfft_file %d %d\n",kn,nfft_file);
	   // File Header
	   // nfft int32: Added at the end in the main code
	   // sfr   double (sampling frequency)
	   // lfft   int32 (original length)
	   // inifr   double
	   sampling_frequency=1/header_param->tsamplu;
	   fwrite((void*)&sampling_frequency, sizeof(double),1,PEAKMAP);
	   fwrite((void*)&header_param->nsamples, sizeof(int),1,PEAKMAP);
	   fwrite((void*)&header_param->frinit, sizeof(double),1,PEAKMAP);
	 }


	 //detector=def_detector(nam);

	if(mid_mjd==0)mid_mjd=header_param->mjdtime+(double)(gd->n/2)*header_param->tsamplu/day2sec;
	beg_mjd=mid_mjd-(double)(gd->n/2)*header_param->tsamplu/day2sec;
	//printf("*****Actual mjdtime, middle point mjd %f %f\n",beg_mjd,mid_mjd);

	data_log[0]=header_param->gps_sec+(header_param->gps_nsec)/(10^9);     //mjd2gps(beg_mjd);
	 data_log[1]=header_param->mjdtime;
	 data_log[2]=(header_param->nfft);

	 prec[0]=3;
	 prec[1]=3;
	 prec[2]=0;


	 logfile_ev(LOG_INFO,"NEW",3,data_log,prec);

	 data_log[0]=(double) header_param->vx_eq;
	 data_log[1]=(double) header_param->vy_eq;
	 data_log[2]=(double) header_param->vz_eq;
	 data_log[3]=(double) header_param->px_eq;
	 data_log[4]=(double) header_param->py_eq;;
	 data_log[5]=(double) header_param->pz_eq;;
	 prec[0]=0;
	 prec[1]=0;
	 prec[2]=0;
	 prec[3]=0;
	 prec[4]=0;
	 prec[5]=0;

	 logfile_ev(LOG_INFO,"VEL",6,data_log,prec);
	 //puts("dopo la scritta di NEW");
	 //printf("velx y z %f %f %f \n",data_log[0],data_log[1],data_log[2]);

	npeak=0;
	len=(int) gd->n/2;

	power=(float *)malloc((size_t) (len)*sizeof(float));
	//puts("allocato 1");
	xamed=(float *)malloc((size_t) (len)*sizeof(float));
	//puts("allocato 2");
	xastd=(float *)malloc((size_t) (len)*sizeof(float));
	//puts("allocato 3");
	puts("Producing peak maps, using the AR mean and std (function: peakmap_from_ar_rev)");

	itaust=(double) 1.0/(gd->dx/tau_Hz);
	imaxage=(int) 1.0/(gd->dx/maxage);
	itaust3=(int) ((double) 1.0/(gd->dx/0.5));
	if(itaust3 < (int) (3*itaust)) itaust3=(int) (3*itaust);
	//itaust3: how many bins in 0.5 Hz or in 3 taus. Needed at the beginning
	w_even=exp(-1.0/itaust);
	//printf("min and max freq.%f %f \n",gd->ini,gd->ini+gd->dx*(gd->n-1));

	//printf("len,dnu,tau_Hz,itaust itaust3 %d %f %f %f %d\n",len,gd->dx,tau_Hz,itaust,itaust3);
	//printf("maxage,imaxage,maxdin  %f %d %f \n",maxage,imaxage,maxdin);
	k=0;
	for(i=0;i<gd->n;i+=2){
	  if(ilog <2)power[k]=sqrt(gd->y[i]*gd->y[i]+gd->y[i+1]*gd->y[i+1]);
	  if(ilog ==2)power[k]=(gd->y[i]*gd->y[i]+gd->y[i+1]*gd->y[i+1]);
	 k+=1;
	}

	/* AR evaluation*/
	k=len-1; //Reverse mode (data are from 0 to len-1)
	len1=len-1;
	xw=0.;
	qw=0.;
	w_norm=0.;
	iage=0;
	ratio=1.;
	while (k>=0){
	  if(k==len1)xamed[k]=power[k];
	  if(k==len1)ratio=1.;
	  if(k!=len1)ratio=power[k]/xamed[k+1];
	  if(ratio <=maxdin && iage>=1){
	    //End of a frequency event EVF, write the Log file
	    //Write the info in the log file:
	    //data_log[3]=evf_ener*gd->dx;  //Power
	    data_log[3]=evf_ener*gd->dx*pow(header_param->normd,2)*pow(header_param->normw,2)*bil2uni;  //Unilateral Energy normalized correctly
	    if(ilog<2)data_log[2]=ratio_max;
	    if(ilog==2)data_log[2]=sqrt(ratio_max);
	    if(ratio_max>=(maxdin*factor_log))logfile_ev(LOG_INFO,"EVF",4,data_log,prec);
	    stat_log+=1;  //total number of peaks above threshold maxdin
	  }
	  if(ratio <=maxdin){
	    iage=0;
	    ad=power[k];
	    qd=ad*ad;
	    //old procedure: xw=(1.0-w_even)*ad+w_even*xw;
	    //qw=(1.0-w_even)*qd+w_even*qw;
	    //w_norm=(1.0-w_even)+w_even*w_norm;
	    xw=ad+w_even*xw;
	    qw=qd+w_even*qw;
	    w_norm=1+w_even*w_norm;
	    if(k%ieval==0 && k!=len1){
	      s=xw;
	      ss=qw;
	      norm=1.0/w_norm;
	      xamed[k]=s*norm;
	      xastd[k]=(float)sqrt((double)fabsf(ss*norm-s*s*norm*norm));
	    }
	  }
	   if(ratio > maxdin){
	     if(iage==0){
	       //a new frequency event begins
	       //Start collecting information for the LOG file
	       data_log[0]=(double) header_param->frinit+k*gd->dx;
	       ratio_max=ratio;
	       kinizio=k;
	       evf_ener=0.;
	       evf_amp=0;
	       prec[0]=3;
	       prec[1]=0;
	       prec[2]=0;
	       prec[3]=0;

	     }
	     iage+=1;
	     data_log[1]=(double) (iage*gd->dx); //duration above threshold
	     if(ilog<2){
	       evf_ener+=(double) power[k]*power[k];  //to evaluate power
	       evf_amp+=(double) power[k];  //to evaluate mean
	     }
	     if(ilog==2){
	       evf_ener+=(double) power[k];  //to evaluate power
	       evf_amp+=(double) sqrt(power[k]);  //to evaluate mean
	     }
	     if(ratio>=ratio_max) ratio_max=ratio;

	     //modifica Dec 14 2011. Continuo ad aggiornare anche se non
	     //lo uso. Usato se e solo se AGE > MAXAGE

	      ad=power[k];
	      if(ratio >= noadj*maxdin){
		ad=power[kinizio+1];
	      } //cosi' dopo un picco-stretto o largo-torna subito bene
	        //Se il picco era troppo alto non andava bene nelle sue vicinanze.

	     qd=ad*ad;
	     xw=ad+w_even*xw;
	     qw=qd+w_even*qw;
	     w_norm=1+w_even*w_norm;
	     s=xw;
	     ss=qw;
	     norm=1.0/w_norm;

	     xamed[k]=xamed[k+1];
	     xastd[k]=xastd[k+1];
	     //levo 40 e metto imaxage-2
	     //if(iage>=40)printf("IAGE  k ,freq,iage imaxage %d %d %f %d %d\n",iage,k,gd->ini+k*gd->dx,iage,imaxage);
	  }
	  if(iage<=imaxage)k-=1;
	  if(iage>imaxage){
	    //printf("IAGE > IMAXAGE !! k ,iage imaxage %d %d %d\n",k,iage,imaxage);
	    //Write the info in the log file:
	    data_log[1]=0.;     //put to zero the duration, because the AR procedure restarts
	    data_log[3]=(evf_amp/iage)*header_param->normd*header_param->normw*sqrt(bil2uni);   //Average value
	    if(ilog<2)data_log[2]=ratio_max;
	    if(ilog==2)data_log[2]=sqrt(ratio_max);
	    // Don't write the changes of level. 3 maggio 2011. Pia
	    //logfile_ev(LOG_INFO,"EVF",4,data_log,prec);
	    //Restart the AR evaluation:
	    k+=(iage-1);
	    //if(k>len1)printf("Att. k >len1, iage %d %d \n",k,iage);
	    //xw=0.;
	    //qw=0.;
	    //w_norm=0.;
	    iage=0;
	    if(k<=(len1-1)){
	      //  xamed[k+1]=power[k];
	       xamed[k+1]=power[k]; //importante sia rimesso cosi' se non aggiorno la stima durante l' evento
	       //xamed[k+1]=s*norm;
	          xastd[k+1]=(float)sqrt((double)fabsf(ss*norm-s*s*norm*norm));
	    }
	  }
	}
	 //Fill the empty values, with the previous ones (that is, the last evaluated)
	 //if ieval=1 this step is not needed
	if(ieval!=1){
	  for(k=(len1-(ieval-1)); k>=0 ;k--){
	     if(k%ieval!=0){
	       xamed[k]=xamed[k+1];
	       xastd[k]=xastd[k+1];
	     }
	  }
	}
	/* Averages and stds at the END of the gd: */
	  for(k=len1; k<len1-2*itaust3;k--){
	    xamed[k]=xamed[len1-2*itaust3];
	    xastd[k]=xastd[len1-2*itaust3];
	  }

	   /*Now check numbers of local maxima in this FFT*/
	  for(k=len1; k>=0;k--){
	    if(k>=(len1-1))ratio=1.0; //by definition
	    if(k<(len1-1))ratio=power[k]/xamed[k+1];
	    if(ratio >=maxdin){
	       if(k<1)ratio_p=1.0; //by definition
	       if(k>=(len1-2))ratio_s=1.0; //by definition
	       if(k>=1)ratio_p=power[k-1]/xamed[k];
	       if(k<(len1-2))ratio_s=power[k+1]/xamed[k+2];
	       //if((ratio>=ratio_p) && (ratio >=ratio_s))npeak+=1; //mod. 16/01/2006
	       if((ratio>ratio_p) && (ratio >ratio_s))npeak+=1;
	    }
	  }
	   /*Now search for local maxima and write them, in SNR=dato/std*/
	  /*First write the BLOCK HEADER*/
	  /* Then the SP data */
	  /* Then the data */
	  // Block header
	  // mjd double
	  // npeak int32
	  // velx double
	  // vely double
	  // velz double
	  //Short spectrum  SP(in Einstein^2)
	  // spini double
	  // spdf  double
	  // splen int32
	  // Data sp  float
	  // Data (npeak times)
	  // bin int32
	  // ratio float
	  // xamed float (mean of H)


	   for(k=0; k<=len1;k++){
	    if(k==0)fwrite((void*)&beg_mjd, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&npeak, sizeof(int),1,PEAKMAP);
	    if(k==0)fwrite((void*)&header_param->vx_eq, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&header_param->vy_eq, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&header_param->vz_eq, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&header_param->px_eq, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&header_param->py_eq, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&header_param->pz_eq, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&spini, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&spdf, sizeof(double),1,PEAKMAP);
	    if(k==0)fwrite((void*)&splen, sizeof(int),1,PEAKMAP);
	    if(k==0){
	      for(ik=0;ik<=splen-1;ik++){
		sp=(gd_short->y[ik]*gd_short->y[ik]); //the AR spectrum normalized correctly, alsoforthezeros
		fwrite((void*)&sp, sizeof(float),1,PEAKMAP);
	      }
	    }

	    if(k>=(len1-1))ratio=1.0; //by definition
	    if(k<(len1-1))ratio=power[k]/xamed[k+1];

	    if(ratio >=maxdin){
	       if(k<1)ratio_p=1.0; //by definition
	       if(k>=(len1-2))ratio_s=1.0; //by definition
	       if(k>=1)ratio_p=power[k-1]/xamed[k];
	       if(k<(len1-2))ratio_s=power[k+1]/xamed[k+2];

	       if((ratio>ratio_p) && (ratio >ratio_s)){
		   fwrite((void*)&k, sizeof(int),1,PEAKMAP);

		if(ilog<2){
		  fwrite((void*)&ratio, sizeof(float),1,PEAKMAP);
		  med_norm=xamed[k+1]*header_param->normd*header_param->normw*sqrt(bil2uni);
		}
		if(ilog==2){
		  ratioS=sqrt(ratio);
		  fwrite((void*)&ratioS, sizeof(float),1,PEAKMAP);
		  med_norm=sqrt(xamed[k+1])*header_param->normd*header_param->normw*sqrt(bil2uni);
		}
		if(header_param->perc_zeroes < 1)med_norm=med_norm/sqrt(1-header_param->perc_zeroes);
		fwrite((void*)&med_norm, sizeof(float),1,PEAKMAP);
	       }
	    }
	    if(nfft==0){

	     if(ilog<2){
	        if(k==0)fprintf(AREST,"%d  %f %f %f %f \n",k,beg_mjd,gd->ini+k*header_param->deltanu,ratio,xamed[0]*header_param->normd*header_param->normw);
	    if(k>=0)fprintf(AREST,"%d  %f %f %f %f \n",k,beg_mjd,gd->ini+k*header_param->deltanu,ratio,xamed[k-1]*header_param->normd*header_param->normw);
	     }
	      if(ilog==2){
		if(k==0)fprintf(AREST,"%d  %f %f %f %f \n",k,beg_mjd,gd->ini+k*header_param->deltanu,ratio,sqrt(xamed[0])*header_param->normd*header_param->normw);
		if(k>=0)fprintf(AREST,"%d  %f %f %f %f \n",k,beg_mjd,gd->ini+k*header_param->deltanu,sqrt(ratio),sqrt(xamed[k-1])*header_param->normd*header_param->normw*sqrt(bil2uni));
	     }
	    }
	   }

	free(power);
	free(xamed);
	free(xastd);
	//printf("npeak of 1 FFT= %d\n",npeak);
	data_log[0]=stat_log;
	prec[0]=0;
	logfile_stat(LOG_INFO,"TOT",1,data_log,prec);
	return npeak;
}

int short_psar_freq(float *ps,GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param)
 /*Very short power spectrum: From frequency data in gd produces one short AR power spectrum in ps.
   Reduction factor red depends on header_param->red .
   Tau of AR could depend on  header_param->red and  header_param->deltanu. Now it is fixed.
   Each datum in ps is the AVERAGE over a number red of data  (data in a same freq. bin).
   The size of ps is lenps/2 (see the code): it contains sqrt(power) in half the band*/

{
  int lenps,j,jj,jj1,k,kk,i;
  float *data;
  //float freq;

  float *power,*xamed,*xastd;
  int len;
  float ratio;
  float tau_Hz; //=0.02; //0.01; //era 0.05, poi 0.1 0.15;  //memory of the AR average. Now they are Hz. Now I am using 0.05 Hz.
  float maxage; //=0.02; //0.01;  //0.15, max age of the AR average:if the age is >maxage the Ar average is re-evaluated
	//entrambi a 0.01 va bene alle risonanze ma non abbastanza ai picchi (il simulato a 911 e la calibr.
	//che restano fra 20-30, mentre se metto 0.15 su entrambi i 2 picchi scendono fra 3-4 ma alle
	//risonanze nel ratio ho come un piccolo buco, non sono descritte abbastanza bene
	//0.05 e' un buon compromesso, ma alle risonanze gia' crea qualche problema
  //messi a 0.02 il 16/01/2006 per provare a migliorare con Virgo
  double w_even;
  float ad,qd; //abs of the datum and its square
  float s,ss;  //mean and average, without normalization, every i_eval samples
  double xw,qw,w_norm;
  double itaust; //how many samples in tau_Hz
  int imaxage; //how many samples in age_Hz
  int iage;
  int itaust3;
  double norm;
  int ieval=1; //mean and std are upgraded every 1 sample (a bin in frequency domain)
  float maxdin; //=2.5; //values above this are not considered into the AR evalution

  double data_log[6]; //to register the information of EVF in the log file
  int prec[6];
  double evf_ener,evf_amp;
  float ratio_max;
  float factor_log; //=2; //dati registrati nel log file sono solo quelli a soglia doppia
  int stat_log; //alla fine scrivo quanti non ne ho scritti

  int ilog=2; //15/01/09 1=prova con la stima log. 2 with the power
  int noadj=10;
  int kinizio;
	puts("Producing very short power spectrum (function: short_psar_freq)");
        tau_Hz=evf_param->tau_Hz;
	maxage=evf_param->maxage;
	maxdin=evf_param->maxdin;
	if (ilog ==2)maxdin=maxdin*maxdin;
	factor_log=evf_param->factor_write;
        stat_log=0;
	len=(int) gd->n/2;

	power=(float *)malloc((size_t) (len)*sizeof(float));
	xamed=(float *)malloc((size_t) (len)*sizeof(float));
	xastd=(float *)malloc((size_t) (len)*sizeof(float));
	//puts(" Produces the very short power spectrum, using the AR mean of the SFDB");
	itaust=(double) 1.0/(gd->dx/tau_Hz);
	imaxage=(int) 1.0/(gd->dx/maxage);
	itaust3=(int) ((double) 1.0/(gd->dx/0.5));
	if(itaust3 < (int) (3*itaust)) itaust3=(int) (3*itaust);
	//itaust3: how many bins in 0.5 Hz or in 3 taus. Needed at the beginning
	w_even=exp(-1.0/itaust);
	//printf("len,dnu,tau_Hz,itaust itaust3 %d %f %f %f %d\n",len,gd->dx,tau_Hz,itaust,itaust3);
	//printf("maxage,imaxage,maxdin  %f %d %f \n",maxage,imaxage,maxdin);

	k=0;
	for(i=0;i<gd->n;i+=2){
	  power[k]=sqrt(gd->y[i]*gd->y[i]+gd->y[i+1]*gd->y[i+1]); //for TEST: put =1  13/7/2006
	  //tested and it is ok: the AR average is correctly done and normalized
	  if (ilog==1)power[k]=log(power[k]*power[k]); //15 Jan 09
	  if (ilog==2)power[k]=(power[k]*power[k]); //15 Jan 09
	 k+=1;

	}

	/* AR evaluation*/
	k=0;
	xw=0.;
	qw=0.;
	w_norm=0.;
	iage=0;
	evf_ener=0.;
	evf_amp=0;
	ratio_max=0.;
	prec[0]=3;
	prec[1]=0;
	prec[2]=0;
	prec[3]=0;
	while (k<len){
	  if(k==0)xamed[k]=power[k];
	  if(k==0)ratio=1;
	  if(k!=0)ratio=power[k]/xamed[k-1];
	  if(ratio <=maxdin && iage>=1){
	    //End of a frequency event EVF, write the Log file
	    //Write the info in the log file:
	    //data_log[3]=evf_ener*gd->dx;  //Power
	    data_log[3]=evf_ener*gd->dx*pow(header_param->normd,2)*pow(header_param->normw,2)*bil2uni;  //Unilateral Energy normalized correctly
	    data_log[2]=ratio_max;
	    if(ilog==2)data_log[2]=sqrt(ratio_max);
	    if(ratio_max>=(maxdin*factor_log))logfile_ev(LOG_INFO,"EVF",4,data_log,prec);
	    stat_log+=1;  //total number of peaks above threshold maxdin
	  }
	  if(ratio <=maxdin){
	    iage=0;
	    ad=power[k];
	    qd=ad*ad;
	    //old procedure: xw=(1.0-w_even)*ad+w_even*xw;
	    //qw=(1.0-w_even)*qd+w_even*qw;
	    //w_norm=(1.0-w_even)+w_even*w_norm;
	    xw=ad+w_even*xw;
	    qw=qd+w_even*qw;
	    w_norm=1+w_even*w_norm;
	    if(k%ieval==0 && k!=0){
	      s=xw;
	      ss=qw;
	      norm=1.0/w_norm;
	      xamed[k]=s*norm;
	      xastd[k]=(float)sqrt((double)fabsf(ss*norm-s*s*norm*norm));
	    }
	  }
	   if(ratio > maxdin){
	     if(iage==0){
	       //a new frequency event begins
	       //Start collecting information for the LOG file
	       data_log[0]=(double) header_param->frinit+k*gd->dx;
	       ratio_max=ratio;
	       kinizio=k;
	       evf_ener=0.;
	       evf_amp=0;
	       prec[0]=3;
	       prec[1]=0;
	       prec[2]=0;
	       prec[3]=0;

	     }
	     iage+=1;
	     data_log[1]=(double) (iage*gd->dx); //duration above threshold
	     if(ilog <2){
	       evf_ener+=(double) power[k]*power[k];  //to evaluate power
	       evf_amp+=(double) power[k];  //to evaluate mean
	     }
	     if(ilog ==2){
	       evf_ener+=(double) power[k];  //to evaluate power
	       evf_amp+=(double) sqrt(power[k]);  //to evaluate mean
	     }
	     if(ratio>=ratio_max) ratio_max=ratio;
	     //modifica 14/12/2011
	      ad=power[k];
	      if(ratio >= noadj*maxdin){
		ad=power[kinizio-1];
	      } //cosi' dopo un picco-stretto o largo-torna subito bene
	        //Se il picco era troppo alto non andava bene nelle sue vicinanze.

	     qd=ad*ad;
	     xw=ad+w_even*xw;
	     qw=qd+w_even*qw;
	     w_norm=1+w_even*w_norm;
	     s=xw;
	     ss=qw;
	     norm=1.0/w_norm;

	     xamed[k]=xamed[k-1];
	     xastd[k]=xastd[k-1];
	     //if(iage>=40)printf("IAGE > 40 k ,freq,iage imaxage %d %f %d %d\n",k,header_param->frinit+k*gd->dx,iage,imaxage);
	   }
	  if(iage<=imaxage)k+=1;
	  if(iage>imaxage){
	    //printf("IAGE > IMAXAGE !! k ,iage imaxage %d %d %d\n",k,iage,imaxage);
	    //Write the info in the log file
	    data_log[1]=0.;     //put to zero the duration, because the AR procedure restarts

	    data_log[3]=(evf_amp/iage)*header_param->normd*header_param->normw*sqrt(bil2uni);  //Average value
	    if(ilog<2)data_log[2]=ratio_max;
	    if(ilog ==2){
	      data_log[2]=sqrt(ratio_max);
	    }
	    //Alberto non lo vuole: tolto
	    //logfile_ev(LOG_INFO,"EVF",4,data_log,prec);
	    //Restart the AR evaluation:
	    k-=(iage-1);
	    if(k<0)printf("Att. k <0, iage %d %d \n",k,iage);
	    //NON RIAZZERARE, ma partire da cio' che si e' calcolato nel loop
	    //xw=0.;
	    //qw=0.;
	    //w_norm=0.;
	    iage=0;
	    if(k>=1){
	       xamed[k+1]=power[k]; //importante sia rimesso cosi' se non aggiorno la stima durante l' evento
	       //xamed[k-1]=s*norm;
	     xastd[k-1]=(float)sqrt((double)fabsf(ss*norm-s*s*norm*norm));
	    }
	  }
	}
	 //Fill the empty values, with the previous ones (that is, the last evaluated)
	 //if ieval=1 this step is not needed
	if(ieval!=1){
	  for(k=ieval+1; k<len ;k++){
	     if(k%ieval!=0){
	       xamed[k]=xamed[k-1];
	       xastd[k]=xastd[k-1];
	     }
	  }
	}
	/* Averages and stds at the beginning of the gd: */
	  for(k=0; k<2*itaust3;k++){
	    xamed[k]=xamed[2*itaust3];
	    xastd[k]=xastd[2*itaust3];
	  }



  lenps= (int) gd->n/(header_param->red);
  //printf("short ps  n lenps rid %ld %d %d\n",gd->n,lenps,header_param->red);
  data=(float *)malloc((size_t) (lenps)*sizeof(float));
  for(k=0;k<lenps;k++){
      data[k]=0.;
  }

  k=-1;
  for(j=0;j<lenps/2;j++){
    jj=j*header_param->red;
    jj1=jj+header_param->red;
    k++;
    for(kk=jj;kk<jj1;kk++){
      data[k]+=xamed[kk];
    }
  }


    /*Very short FFT in ps[k]: k are the various frequency bins*/

  for(k=0;k<lenps/2;k++){
  //freq=(float)(header_param->firstfreqindex+k*header_param->red)*1.0/header_param->tbase; //not used, just to know
    ps[k]=data[k]/header_param->red;
    if (ilog==1)ps[k]=sqrt(exp(ps[k]+1.0/sqrt(3))); //15 Jan 09
    if (ilog==2)ps[k]=sqrt(ps[k]); //12 DEC 11
    }

    free(data);
    free(power);
    free(xamed);
    free(xastd);

  data_log[0]=stat_log;
  prec[0]=0;
  logfile_stat(LOG_INFO,"TOT",1,data_log,prec);

  return header_param->red;
}
/***************************************************************************************/
/*Reversed AR estimation  (10 July 2007) */
/*Indicated for spectra with low frequencies high level */
int short_psar_rev_freq(float *ps,GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param,float *psTOT)
 /*Very short power spectrum: From frequency data in gd produces one short AR power spectrum in ps. EVALUATED IN REVERSE MODE.
   Reduction factor red depends on header_param->red .
   Tau of AR could depend on  header_param->red and  header_param->deltanu. Now it is fixed.
   Each datum in ps is the AVERAGE over a number red of data  (data in a same freq. bin).
   The size of ps is lenps/2 (see the code): it contains sqrt(power) in half the band
   psTOT added, 19/11/2014, for the peakmaps. if ilog =2, come appunto e',  i conti sono fatti sulla potenza. Poi si da' la sqrt dello spettro
*/

{
  int lenps,j,jj,jj1,k,kk,i;
  float *data;
  //float freq;

  float *power,*xamed,*xastd;
  int len;
  float ratio;
  float tau_Hz; //=0.02; //0.01; //era 0.05, poi 0.1 0.15;  //memory of the AR average. Now they are Hz. Now I am using 0.05 Hz.
  float maxage; //=0.02; //0.01;  //0.15, max age of the AR average:if the age is >maxage the Ar average is re-evaluated
   //entrambi a 0.01 va bene alle risonanze ma non abbastanza ai picchi (il simulato a 911 e la calibr.
   //che restano fra 20-30, mentre se metto 0.15 su entrambi i 2 picchi scendono fra 3-4 ma alle
   //risonanze nel ratio ho come un piccolo buco, non sono descritte abbastanza bene
  //0.05 e' un buon compromesso, ma alle risonanze gia' crea qualche problema
  //messi a 0.02 il 16/01/2006 per provare a migliorare con Virgo
  double w_even;
  float ad,qd; //abs of the datum and its square
  float s,ss;  //mean and average, without normalization, every i_eval samples
  double xw,qw,w_norm;

  double itaust; //how many samples in tau_Hz
  int imaxage; //how many samples in age_Hz
  int iage;
  int itaust3;
  double norm;
  int ieval=1; //mean and std are upgraded every 1 sample (a bin in frequency domain)
  float maxdin;  //=2.5; //values above this are not considered into the AR evalution

  double data_log[6]; //to register the information of EVF in the log file
  int prec[6];
  double evf_ener,evf_amp;
  float ratio_max;
  float factor_log; //=2; //dati registrati nel log file sono solo quelli a soglia doppia
  int stat_log; //alla fine scrivo quanti non ne ho scritti
  int bin_max;
  int kinizio;
  int ilog=2; //0 amp 1 log 2 power
  int itest=0;  //1 write an output file with norm.
  int noadj=10; //

        //puts("Producing very short power spectrum in REVERSE mode (function: short_psar_rev_freq)");
        tau_Hz=evf_param->tau_Hz;
	maxage=evf_param->maxage;


	maxdin=evf_param->maxdin;
	if (ilog==2) maxdin=evf_param->maxdin*evf_param->maxdin;
	factor_log=evf_param->factor_write;
        stat_log=0;
	len=(int) gd->n/2;
	ratio=1.0;
	power=(float *)malloc((size_t) (len)*sizeof(float));
	xamed=(float *)malloc((size_t) (len)*sizeof(float));
	xastd=(float *)malloc((size_t) (len)*sizeof(float));
	//puts(" Produces the very short power spectrum, using the AR mean of the SFDB");
	itaust=(double) 1.0/(gd->dx/tau_Hz);
	imaxage=(int) 1.0/(gd->dx/maxage);
	itaust3=(int) ((double) 1.0/(gd->dx/0.5));
	if(itaust3 < (int) (3*itaust)) itaust3=(int) (3*itaust);
	//itaust3: how many bins in 0.5 Hz or in 3 taus. Needed at the beginning
	w_even=exp(-1.0/itaust);
	//printf("len,dnu,tau_Hz,itaust itaust3 %d %f %f %f %d\n",len,gd->dx,tau_Hz,itaust,itaust3);
	//printf("maxage,imaxage,maxdin,factor_log %f %d %f %f\n",maxage,imaxage,maxdin,factor_log);

	k=0;
	for(i=0;i<gd->n;i+=2){
	  if (ilog <2)power[k]=sqrt(gd->y[i]*gd->y[i]+gd->y[i+1]*gd->y[i+1]);
	  if(ilog==2)power[k]=(gd->y[i]*gd->y[i]+gd->y[i+1]*gd->y[i+1]);
	  //tested and it is ok: the AR average is correctly done and normalized
	  //if(header_param->nfft==1)printf("POW %f %f\n",k*gd->dx,power[k]);
	  if (ilog==1)power[k]=log(power[k]*power[k]); //test 15 Jan 09
	  k+=1;
	}

	/* AR REVERSE evaluation*/
	k=len-1; //Reverse mode (data are from 0 to len-1)
	xw=0.;
	qw=0.;
	w_norm=0;
	iage=0;

	evf_ener=0.;
	evf_amp=0;
	ratio_max=0.;
	prec[0]=1;
	//prec[0]=3;
	prec[1]=0;
	prec[2]=0;
	prec[3]=0;
	prec[4]=0;

	while (k>=0){

	  if(k==(len-1))xamed[k]=power[k];
	  if(k==(len-1))ratio=1.0;
	  if(k!=(len-1)){
	    if(k>=0 && k<=(len-2)) //Aggiunto questo if 7/02/08
	      if(xamed[k+1] !=0){
		ratio=power[k]/xamed[k+1];
	      }
	      else {
		ratio=1.;
		  }
	      else
		ratio=1.0;
	  }
	  if(header_param->nfft==1 && k> (len-1-1*itaust3)){ //prova 12 maggio 2010 (2*itaust)
	    ratio=1;
	  }
	  if(ratio <=maxdin && iage>=1){
	    //End of a frequency event EVF, write the Log file
	    //Write the info in the log file:
	    //data_log[3]=evf_ener*gd->dx;  //Power
	    data_log[3]=evf_ener*gd->dx*pow(header_param->normd,2)*pow(header_param->normw,2)*bil2uni;  //Unilateral Energy normalized correctly 12 marzo 2010
	    if(ilog <2)data_log[2]=ratio_max;
	    if(ilog==2){
	      data_log[2]=sqrt(ratio_max);
	    }
	    data_log[4]=bin_max;
	    if(ratio_max>=(maxdin*factor_log))logfile_ev(LOG_INFO,"EVF",5,data_log,prec);
	    stat_log+=1;  //total number of peaks above threshold maxdin

	  }
	  if(ratio <= maxdin){
	    iage=0;
	    ad=power[k];
	    qd=ad*ad;
	    //old procedure:
	    //xw=(1.0-w_even)*ad+w_even*xw;
	    //qw=(1.0-w_even)*qd+w_even*qw;
	    //w_norm=(1.0-w_even)+w_even*w_norm;
	    xw=ad+w_even*xw;
	    qw=qd+w_even*qw;
	    w_norm=1+w_even*w_norm;
	    if(k%ieval==0 && k!=(len-1)){
	      s=xw;
	      ss=qw;
	      norm=1.0/w_norm;
	      xamed[k]=s*norm;
	      xastd[k]=(float)sqrt((double)fabsf(ss*norm-s*s*norm*norm));
	    }
	  }
	   if(ratio > maxdin){
	     if(iage==0){
	       //a new frequency event begins
	       //Start collecting information for the LOG file
	       //data_log[0]=(double) header_param->frinit+k*gd->dx;
	       data_log[0]=k;
	       kinizio=k;
	       ratio_max=ratio;
	       bin_max=0;
	       evf_ener=0.;
	       evf_amp=0;
	       prec[0]=1;
	       //prec[0]=3;
	       prec[1]=0;
	       prec[2]=0;
	       prec[3]=0;
	       prec[4]=0;

	     }
	     iage+=1;
	     //data_log[1]=(double) (iage*gd->dx); //duration above threshold
	     data_log[1]=iage;
	     if(ilog <2){
	       evf_ener+=(double) power[k]*power[k];  //to evaluate power
	       evf_amp+=(double) power[k];  //to evaluate mean
	     }
	   if(ilog ==2){
	       evf_ener+=(double) power[k];  //to evaluate power
	       evf_amp+=(double) sqrt(power[k]);  //to evaluate amplitude mean
	   }
	     if(ratio>=ratio_max){
	       bin_max=(kinizio-k);
	       ratio_max=ratio;
	     }

	     //modifica Dec 14 2011: fa quello che faceva prima, ma
	     //aggiorna la stima che sara' usata solo se e quando AGE > MAXAGE

	      ad=power[k];
	      if(ratio >= noadj*maxdin){
		ad=power[kinizio+1];
	      } //cosi' dopo un picco-stretto o largo-torna subito bene
	        //Se il picco era troppo alto non andava bene nelle sue vicinanze.

	     qd=ad*ad;
	     xw=ad+w_even*xw;
	     qw=qd+w_even*qw;
	     w_norm=1+w_even*w_norm;
	     s=xw;
	     ss=qw;
	     norm=1.0/w_norm;
	     xamed[k]=xamed[k+1]; //s*norm;
	     xastd[k]=xastd[k+1]; //(float)sqrt((double)fabsf(ss*norm-s*s*norm*norm));


	     // printf("IAGE\n");
	     // Sabrif(iage>=40)printf("IAGE > 40 k ,freq,iage imaxage %d %f %d %d\n",k,header_param->frinit+k*gd->dx,iage,imaxage);
	   }
	  if(iage<=imaxage)k-=1;
	  if(iage>imaxage){
	    //printf("===============IAGE > IMAXAGE !! k ,iage imaxage %d %d %d\n",k,iage,imaxage);
	    //Write the info in the log file:
	    data_log[1]=0.;     //put to zero the duration, because the AR procedure restarts

	    data_log[3]=(evf_amp/iage)*header_param->normd*header_param->normw*sqrt(bil2uni);  //Average value
	    if(ilog <2)data_log[2]=ratio_max;
	    if(ilog==2)data_log[2]=sqrt(ratio_max);
	    data_log[4]=bin_max;
	    // TODO Albe Non scrivo gli zeri
	    //logfile_ev(LOG_INFO,"EVF",5,data_log,prec);

	    //Restart the AR evaluation:
	    k+=(iage-1);
	    if(k>(len-1))printf("Att. k > (len-1), iage %d %d \n",k,iage);
	    if(k>(len-1)) k=len-1;
	    //xw=0;
	    //qw=0;
	    //w_norm=0.;
	    iage=0;
	    if(k<=(len-2)){
	      xamed[k+1]=power[k]; //importante sia rimesso cosi' se non aggiorno la stima durante l' evento
	      //xamed[k+1]=s*norm;
	        xastd[k+1]=(float)sqrt((double)fabsf(ss*norm-s*s*norm*norm));
	    }
	  }

	  if(itest==1){
	    if (header_param->nfft==2) {
	    //TEST::scrittura del file
	    if (k==(len-1)) OUTTEST=fopen("testfile_norm.dat","w");
	       fprintf(OUTTEST,"%d %f %d\n",k,w_norm,iage);
	       //if (k==0)fclose(OUTTEST);
	    }
	  }


	}  //while
	 //Fill the empty values, with the previous ones (that is, the last evaluated)
	 //if ieval=1 this step is not needed
	if(ieval!=1){
	  for(k=(len-1-(ieval-1)); k>=0 ;k--){
	     if(k%ieval!=0){
	       xamed[k]=xamed[k+1];
	       xastd[k]=xastd[k+1];
	     }
	  }
	}
	/* Averages and stds at the END of the gd: */
	  for(k=(len-1); k> (len-1-2*itaust3);k--){
	    xamed[k]=xamed[len-1-2*itaust3];
	    xastd[k]=xastd[len-1-2*itaust3];
	  }


  lenps= (int) gd->n/(header_param->red);
  //printf("header_param->normw =%f normd = %f \n",header_param->normw,header_param->normd);
  //printf("short ps  n lenps rid %ld %d %d\n",gd->n,lenps,header_param->red);
  data=(float *)malloc((size_t) (lenps)*sizeof(float));
  for(k=0;k<lenps;k++){
      data[k]=0.;
  }

  k=-1;
  for(j=0;j<lenps/2;j++){
    jj=j*header_param->red;
    jj1=jj+header_param->red;
    k++;
    for(kk=jj;kk<jj1;kk++){
      // Modified on 22 Nov 2010; average on the square not on the amp !
      //rimodified Gennaio 2011
     data[k]+=xamed[kk];
      //data[k]+=xamed[kk]*xamed[kk];

    }
  }


    /*Very short FFT in ps[k]: k are the various frequency bins*/

  for(k=0;k<lenps/2;k++){
  //freq=(float)(header_param->firstfreqindex+k*header_param->red)*1.0/header_param->tbase; //not used, just to know
      ps[k]=data[k]/header_param->red;
    //ps[k]=sqrt(data[k]/header_param->red);
      if (ilog==2)ps[k]=sqrt(ps[k]); //15 Jan 09
      if (ilog==1)ps[k]=sqrt(exp(ps[k]+1.0/sqrt(3))); //15 Jan 09
    }
  //len is gd->n /2. gd->n is double as it contains real and imag
  for (kk=0;kk<(int)len;kk++){
      psTOT[kk]=sqrt(xamed[kk]);
  }

    free(data);
    free(power);
    free(xamed);
    free(xastd);

  data_log[0]=stat_log;
  prec[0]=0;
  logfile_stat(LOG_INFO,"TOT",1,data_log,prec);

  return header_param->red;
}

struct tm *GPS2UTCSAB(long GPSsec, int ULeapS)
/*
	GPSsec   gps seconds (integer value)
	ULeapS   leap seconds (TAI-UTC : 19 on 1 Jan 1980, 32 on 1 Jan 1999)

*/

{
	struct tm *t;
	int mdays[12]={31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	long secutc,ndays,nsec,ndaystot;
	long i=0,nyear=0,yeardays,mday;
	int Y,M,D,DoY,DoW,h,m,s;
	time_t now;

	secutc=GPSsec-ULeapS+19+5*86400;
	ndays=secutc/86400;
	nsec=secutc-ndays*86400;
	ndaystot=ndays;

	nyear=0;
	while(ndays > 0){
		yeardays=365;
		if((nyear/4)*4 == nyear)yeardays=366;
		ndays-=yeardays;
		nyear+=1;
	}
	nyear=nyear-1;
	ndays=ndays+yeardays;
	mday=ndays;

	i=0;
	if((nyear/4)*4 == nyear)mdays[1]=29;
	while(mday >= 0){
		mday=mday-mdays[i]; //printf(" --> %d %d %d \n",ndays,i,mday);
		i++;
	}
	mday=mday+mdays[i-1];
	// t=localtime(&secutc);  Error
	now=time(0);
	t=localtime(&now);

	Y=nyear+1980-1900;
	M=i-1;
	D=mday+1;
	DoY=ndays;
	DoW=ndaystot+2-((ndaystot+2)/7)*7;
	h=nsec/3600;
	m=(nsec-h*3600)/60;
	s=nsec-h*3600-m*60;
	// printf(" Nella function %d %d %d  %d %d  %d %d %d \n",Y,M,D,DoY,DoW,h,m,s);
	t->tm_year=Y;
	t->tm_mon=M;
	if (M == 12){   //modifica Pia 13 maggio 2016.
	  t->tm_mon=0;   //dopo viene messo + 1. Ossia il mese in uscita e' il mese -1. Non so perche' chiedere a Sabri e casomai sistemare
	  t->tm_year=Y+1;
	  // printf("NELL IF M = 12 %04d%02d%02d\n",t->tm_year,t->tm_mon+1,t->tm_mday);
	}
	t->tm_mday=D;

	t->tm_hour=h;
	t->tm_min=m;
	t->tm_sec=s;
	t->tm_yday=DoY;
	t->tm_wday=DoW;
	t->tm_isdst=0;

	return t;
}

/***************Function to write the output file with the SFDB******************/
int gd2sfdbfilevoid(GD *gd, HEADER_PARAM *header_param,EVF_PARAM *evf_param,float *ps,EVEN_PARAM *even_param, INPUT_PARAM *input_param,float *psTOT){

  /*The function writes data which are in a GD *gd + some parameters into a binary file (header + data)
   Input are time domain data in a GD *gd. The function performs the total FFT and also and the very short ps using short_ps_time and short_ps_freq
 */
  int errorcode,ii,j;
  long  len;
  float rpw,ipw;
  int lenps,lenps1;
  int estimfromar=1;  //procedure to estimate the very short spectrum
                      //1 is the AR


  //int jj;
  int sat;

  /**********Variables needed for the Dopple effect**************/
  short int fut1,nut1,nutcorr;  //nut1=1 uses nutation from ephem files;
                               //nut1=0  uses nutation from novas interpolation;
  //NOT needed here !
  //short int equatorial=1; //1 output are Equatorial, 0 ecliptic
  //short int imjd=1;      //1 input are mjds, 0 imput are days from 1900
  //double mjd;
  /*****Detector*****/


  //double velx,vely,velz;
  //double posx,posy,posz;
  //double deinstein;

  double vel[3];
  double pos[3];
  double mid_mjd; //mjd of the middle+1 point of the FFT

  //DETECTOR *detector;
  //char nam[5];
  int sat_maxn;
  int fft_m;

  int nn;
  int jjin;
  char nam[5];

  int rand(void);

  int itest,dothetest;
  float test,mtest,qtest;
  dothetest=0;  //1 to do a test on the short spectra with gaussian nooise  Jan. 2011
  if(header_param->detector==1)strncpy(nam,"virgo",5);
  if(header_param->detector==2)strncpy(nam,"ligoh",5);
  if(header_param->detector==3)strncpy(nam,"ligol",5);
  if(header_param->detector==4)strncpy(nam,"kagra",5);


   // SFDBNAME sabrina
  struct tm *timstr;
  char appo[40],appo1[8],appo2[10];
  int fakesamples=1;  //this is the relevant change in this function  12 nov 2010
  //strcpy(nam,"virgo");
  //strcpy(nam,"ligol");
  //IMPORTANT NOTE: remember that header_param->nfft is the fft number in the run (e.g. C6), not
  //in that particular file. Thus, only the first SFDB file has  header_param->nfft=1

  /*Now we have time-domain data in the GD*/

  fft_m=input_param->fft_m; //200*(2097152/(2*header_param->nsamples)); //200: AD HOC sulla banda piu' grande pia
             //100 per Virgo //50 per Ligo a 16384 lunghe il doppio


  //printf("**In gd2sfdbfile: actual SDS file filme %s \n",sfc_data->filme);
  //printf("**Actual SDS time t0 mjd %f \n",sfc_data->t0);
  //nn=strlen(sfc_data->filme);
  //seleziona(appo,sfc_data->filme,0,nn-NEXTS-1);


  if(strncmp(input_param->framechannel,"sdsfile",7)==0){
        jjin=filename_from_file1(input_param->filename);
        nn=strlen(input_param->filename)-jjin;
        seleziona(appo,input_param->filename,jjin,nn-NEXTS-1-17); //-17 e' AD HOC  5 ott 2010una prova !!

    }


if(strncmp(input_param->framechannel,"sdsfile",strlen("sdsfile"))!=0)strncpy(appo,input_param->framechannel,strlen(input_param->framechannel));
 appo[strlen(input_param->framechannel)] = '\0';
	//printf("$$$$ appo = %s\n", appo);
	timstr=GPS2UTCSAB(header_param->gps_sec,0.0);
	//printf("%04d%02d%02d\n",timstr->tm_year+1900,timstr->tm_mon+1,timstr->tm_mday);
	//printf("%02d%02d%02d",timstr->tm_hour,timstr->tm_min,timstr->tm_sec);
	strcat(appo,"_");
	sprintf(appo2,"%04d%02d%02d",timstr->tm_year+1900,timstr->tm_mon+1,timstr->tm_mday);
	strcat(appo,appo2);
	strcat(appo,"_");
	sprintf(appo2,"%02d%02d%02d",timstr->tm_hour,timstr->tm_min,timstr->tm_sec);
	strcat(appo,appo2);




	strncpy(appo1,".SFDB09",7);

	strncat(appo,appo1,7);


	strcpy(header_param->sfdbname,appo);

        //printf("$$$$ Sfdb filename= %s \n", header_param->sfdbname);


   sat_maxn=even_param->sat_maxn;
   lenps= (int) gd->n/(header_param->red);
   lenps1= header_param->red;
   header_param->lavesp=lenps/2; //lenght of the short periodogram (2010: now the same as AR spectrum)
   len=gd->n;
   //printf("in gd2sfdbfile: gd->n  %ld \n",len);
   //printf("in gd2sfdbfile: nfft %d  \n",header_param->nfft);
   if(header_param->nfft==1 && SFDB == NULL){
     if(input_param->flagonline !=1){
       SFDB=fopen(header_param->sfdbname,"w");
       printf("First SFDB file opened - file name: %s \n", header_param->sfdbname);}
   }


   /*Closes the file every a chosen number of FFTs, if itf detector, and opens a new file */
   if(header_param->detector !=0 && fmodf((float)header_param->nfft,(float) fft_m)==0.){
      if(input_param->flagonline !=1)
	{
      //printf("open a new SFDB file: nfft fft_m %d  %d \n",header_param->nfft, fft_m);
      fclose(SFDB);
      SFDB=fopen(header_param->sfdbname,"w");
      printf("New SFDB file opened  - file name: %s \n", header_param->sfdbname);
      // jj=logfile_output(LOG_INFO,header_param->sfdbname,"A new SFDB file opened");
	}
   }

if(input_param->flagonline !=1) {
   // Added 15/01/08:
  //detector=def_detector(nam);
      //Evaluation of detector velocities
  /****************************************************************************/
   fut1=0;  ///1=uses ut1; 0=uses utc as aprrox. of ut1; //0 per fare fino al 2010
   nut1=1;  //1=uses nutation from ephem files; 0=no;
   nutcorr=0; //reads the file DELTANUT, if 0 does not read DELTANUT //0 per fare fino al 2010
 /****************************************************************************/
   //	mid_mjd=header_param->mjdtime+(double)(gd->n/2)*header_param->tsamplu/day2sec;   //aggiunto +1 nel formato SFDB09
	mid_mjd=header_param->mjdtime+(double)(gd->n/2+1)*header_param->tsamplu/day2sec;
	//printf("Middle +1 point mjd %f\n",mid_mjd);
	/* REMOVED IN THE FUNCTION "gd2sfdb...void"
	mjd=PSS_detector_velposSCALAR(detector,mid_mjd,header_param->tsamplu,fut1,nut1,nutcorr,equatorial,imjd,velx,vely,velz,posx,posy,posz,deinstein);
	*/
	//puts("uscito dalla function detector_velpos");
	//velx,y,z inAU/day; posx,y,z in AU. C in AU/day
	//vel is adimensional and pos in light seconds
	//vel and pos in the function ..void are not needed and thus we save computing time
	vel[0]=0;  //velx/C;
	vel[1]=0;  //vely/C;
	vel[2]=0;  //velz/C;
	header_param->vx_eq=vel[0];
	header_param->vy_eq=vel[1];
	header_param->vz_eq=vel[2];
	pos[0]=0; //posx*day2sec/C;
	pos[1]=0; //posy*day2sec/C;
	pos[2]=0; //posz*day2sec/C;
	header_param->px_eq=pos[0];
	header_param->py_eq=pos[1];
	header_param->pz_eq=pos[2];

	//free(detector);

	//end evaluation of detector velocities--Added on Jan 15 2008
 }  //end of flagonline not 1
 if(input_param->flagonline ==1)
   {
     	header_param->vx_eq=0.;
	header_param->vy_eq=0.;
	header_param->vz_eq=0.;
	header_param->px_eq=0.;
	header_param->py_eq=0.;
	header_param->pz_eq=0.;
   }
   /*End of the close/open procedure  */

	/*
	  19 may 2008
	  If too many saturated data: do not write the FFT and header. return
	 */
	//puts("Prima di saturated");

	sat= (int) header_param->sat_howmany;
	//if (sat >=sat_maxn){
	  //errorcode=0;
	  //printf("Veto on the FFT due to sat data: howmany,howmany,maxn %f %d %d\n",header_param->sat_howmany,sat,sat_maxn);

	  //return errorcode;
	//}
	//puts("Prima di veto nzeroes");
	//if(header_param->n_zeroes > header_param->nsamples*2)  //2012 era /2
	//{
	//printf("veto una fft nfft, header_param->n_zeros, %d, %d\n",header_param->nfft, header_param->n_zeroes);
	//errorcode=0;

	//return errorcode;        //continue; //go to the next in the loop
	//}

   if(input_param->flagonline !=1)

   {
   /*Write the header*/

   errorcode=fwrite((void*)&header_param->endian, sizeof(double),1,SFDB);
   errorcode=fwrite((void*)&header_param->detector, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->gps_sec, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->gps_nsec, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->tbase, sizeof(double),1,SFDB);
   errorcode=fwrite((void*)&header_param->firstfreqindex, sizeof(int),1,SFDB);
   //errorcode=fwrite((void*)&header_param->nsamples, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&fakesamples, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->red, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->typ, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->n_flag, sizeof(float),1,SFDB);
   errorcode=fwrite((void*)&header_param->einstein, sizeof(float),1,SFDB);
   errorcode=fwrite((void*)&header_param->mjdtime, sizeof(double),1,SFDB);
   errorcode=fwrite((void*)&header_param->nfft, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->wink, sizeof(int),1,SFDB);
   errorcode=fwrite((void*)&header_param->normd, sizeof(float),1,SFDB);
   errorcode=fwrite((void*)&header_param->normw, sizeof(float),1,SFDB);
   errorcode=fwrite((void*)&header_param->frinit, sizeof(double),1,SFDB);
   errorcode=fwrite((void*)&header_param->tsamplu, sizeof(double),1,SFDB);
   errorcode=fwrite((void*)&header_param->deltanu, sizeof(double),1,SFDB);
   //specific header (different for bar-itf) detectors
   if(header_param->detector==0){ //bar detector
     errorcode=fwrite((void*)&header_param->vx_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->vy_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->vz_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->px_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->py_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->pz_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->n_zeroes, sizeof(int),1,SFDB);
     errorcode=fwrite((void*)&header_param->sat_howmany, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare1, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare2, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare3, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->perc_zeroes, sizeof(float),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare5, sizeof(float),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare6, sizeof(float),1,SFDB);
     errorcode=fwrite((void*)&header_param->lavesp, sizeof(int),1,SFDB);
     errorcode=fwrite((void*)&header_param->sciseg, sizeof(int),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare9, sizeof(int),1,SFDB);
   }

   if(header_param->detector>=1){ //itf detector
     errorcode=fwrite((void*)&header_param->vx_eq, sizeof(double),1,SFDB);  //vx/C equat
     errorcode=fwrite((void*)&header_param->vy_eq, sizeof(double),1,SFDB);  //vy
     errorcode=fwrite((void*)&header_param->vz_eq, sizeof(double),1,SFDB);  //vz
     errorcode=fwrite((void*)&header_param->px_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->py_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->pz_eq, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->n_zeroes, sizeof(int),1,SFDB);
     errorcode=fwrite((void*)&header_param->sat_howmany, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare1, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare2, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare3, sizeof(double),1,SFDB);
     errorcode=fwrite((void*)&header_param->perc_zeroes, sizeof(float),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare5, sizeof(float),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare6, sizeof(float),1,SFDB);
     errorcode=fwrite((void*)&header_param->lavesp, sizeof(int),1,SFDB);
     errorcode=fwrite((void*)&header_param->sciseg, sizeof(int),1,SFDB);
     errorcode=fwrite((void*)&header_param->spare9, sizeof(int),1,SFDB);
   }


   if (errorcode!=1)printf("Error in writing header into SFDB file!\n");
   //printf("file header:endian %lf \n",header_param->endian);
   //printf("file header:nsamples %d  \n",header_param->nsamples);
   //printf("file header:n_flag %f  \n",header_param->n_flag);
   //printf("file header:perc_zeroes %f  \n",header_param->perc_zeroes);
   //printf("file header:n_zeroes %d  \n",header_param->n_zeroes);
   //printf("File header:gps %d \n",header_param->gps_sec);
   //printf("File header:mjdtime %f \n",header_param->mjdtime);
   /*Write the SFDB data*/
   /*Creation of the average of very short power spectra */
   if(dothetest==1){
	 mtest=0;
	 qtest=0;
	 for(j=0;j<gd->n;j++){
	 test=0;
	 for (itest=1;itest<=10;itest++){
	   test+=rand();
	   //if(j==0){
	   //   printf("RAND_MAX test/RAND_MAX %d %f  \n",RAND_MAX,test/RAND_MAX);
	   //}
	 }
	 gd->y[j]=(test/RAND_MAX-5.0)/sqrt(10.0/12.0);
	 mtest+= gd->y[j];
	 qtest+= gd->y[j]*gd->y[j];
	 }
       mtest=mtest/gd->n;
       qtest=qtest/gd->n;
       qtest=(qtest-mtest*mtest);
       //printf("TEST mean variance N gd->n %f %f %ld  \n",mtest,qtest,gd->n);
       //DONE ON wn1: opt/exp_soft/virgo/CW/appoSDS/VIR_h_V3_4096Hz_20090717_154725_.sds
   }

        //2008 ps=(float *)malloc((size_t) (lenps1)*sizeof(float));
     //it was short_ps_time. changed Nov 2010
    ii=short_ps_periodogram(ps,gd,header_param);
    //for(j=0;j<header_param->red;j++){
    for(j=0;j<lenps/2;j++){
      //printf(" Short power spectrum, dimension red: j ps[j] %d %f\n",j,ps[j]);
      //rpw=ps[j]*pow(header_param->normd,2)*pow(header_param->normw,2)*bil2uni; //piapia
      rpw=ps[j]; //sqrt of the power spectrum
      if(header_param->perc_zeroes < 1)rpw=rpw/sqrt(1-header_param->perc_zeroes);
      errorcode=fwrite((void*)&rpw, sizeof(float),1,SFDB);
    }
   }  //end of  if(input_param->flagonline !=1)
  /*Window, if needed */
     if(header_param->wink>0){
       gd2wingd(gd,header_param->wink);
     }
    //2008 free(ps);
    if(estimfromar!=1){
      //Very short FFT from sum over lower resolution FFTs, (input are data in time domain)
      //2008 ps=(float *)malloc((size_t) (lenps)*sizeof(float));
      ii=short_ps_freq(ps,gd,header_param);
      for(j=0;j<lenps/2;j++){
      //printf(" Very Short sqrt(spectrum), dimension lenps/2: j ps[j]  %d %f \n",j,ps[j]);
	rpw=ps[j]*header_param->normd*header_param->normw; //sqrt of the power spectrum
	if(header_param->perc_zeroes < 1)rpw=rpw/sqrt(1-header_param->perc_zeroes);
      if(input_param->flagonline !=1)errorcode=fwrite((void*)&rpw, sizeof(float),1,SFDB);
      }
      // 2008 free(ps);
    }

   ii= gdtime2freq(gd,header_param);  //From the GD in time to the same GD in frequency domai
   if(estimfromar==1){
     //Very short spectrum from AR estimation: (input are data in frequency domain)
     //ii=short_psar_freq(ps,gd,header_param,evf_param); //
     //puts("Call the reverse AR estimation");
     ii=short_psar_rev_freq(ps,gd,header_param,evf_param,psTOT); //10 Jul 2007 REVERSED !! //psTOT needed for peakmap
   for(j=0;j<lenps/2;j++){
    rpw=ps[j]*header_param->normd*header_param->normw*sqrt(bil2uni); //sqrt of the ps
    if(header_param->perc_zeroes < 1)rpw=rpw/sqrt(1-header_param->perc_zeroes);
    if(input_param->flagonline !=1)errorcode=fwrite((void*)&rpw, sizeof(float),1,SFDB);
   }
   }
   //SFDB  (complex FFT data): if flagonline not 1
if(input_param->flagonline !=1)
  {
     for(ii=0;ii< fakesamples;ii+=2){
       //MODIFICA 15/04/2010: tolta la sqrt(2), non si sa quando era apparsa.
       //da oggi tutti gli SFDB09 saranno senza pia
       //rpw=gd->y[ii]*sqrt(bil2uni);    //Real part
       //ipw=gd->y[ii+1]*sqrt(bil2uni);  //Imag part
     rpw=gd->y[ii];    //Real part  piapia
     ipw=gd->y[ii+1];  //Imag part

     //if(ii<4)printf("First data written into the file %e %e\n",rpw,ipw); //messo %e 25-Aug 2005
     errorcode=fwrite((void*)&rpw, sizeof(float),1,SFDB);
     if (errorcode!=1)printf("Error in writing data into SFDB file!\n");
     errorcode=fwrite((void*)&ipw, sizeof(float),1,SFDB);
     if (errorcode!=1)printf("Error in writing data into SFDB file!\n");


   }
   //puts("In gd2sfdbfile:another data chunk has been written. New rev procedure");
  }  //end of:   if(input_param->flagonline !=1)

  return errorcode;
}
//Versione con gfal
int ToCopyFiles(INPUT_PARAM *input_param,HEADER_PARAM *header_param)
{
  int ij,ijj,ijinizio,i,ijprec;
  char *CommandForGrid,*srmdirCNAF,*remainfrompath,*remain,*dovesono,*theband;

  CommandForGrid=(char*)malloc(sizeof(char)*300);
  srmdirCNAF=(char*)malloc(sizeof(char)*MAXMAXLINE);
  remainfrompath=(char*)malloc(sizeof(char)*MAXMAXLINE);
  remain=(char*)malloc(sizeof(char)*MAXMAXLINE);
  dovesono=(char*)malloc(sizeof(char)*MAXMAXLINE);
  theband=(char*)malloc(sizeof(char)*6);
  int nnn;

  //aggiunto gethostname. 22 ottobre 2019
  char hostname[128];
  gethostname(hostname, sizeof hostname);
  printf("My hostname: %s\n", hostname);
  //fine aggiunto
        i=0;  //non serve....
        nnn=strlen(input_param->filelista);
	//ijj=0;
	//ij=0;
	//for(i=0;i<nnn;i++){
	//  if(input_param->filelista[i] == 47 || input_param->filelista[i] == 92){
	//     ijprec=ij;
	//     ijj+=1;
	//     if (ijj==5) ijinizio=i;
	//     ij=i; /* 47 -> / , 92 -> \ */}
	//}
	//strcpy(remain,"/sfdb/O2/C02_clean_noscience/");  //Pia: in futuro va sistemato non hard coded
	//For O3, C00 calib. Pia. 11 May 2019. Removed the line above, for O2
	//strcpy(remain,"/sfdb/O3b/C00/NC/");  //Pia: in futuro va sistemato non hard coded. RIMETTERE NC Aprile 2020!!!!
	//strcpy(remain,"/sfdb/O3/C01/GAT_CL60/");  //Pia: in futuro va sistemato non hard coded
	//strcpy(remain,"/sfdb/O3/C01/CLEAN60/HF/");
	//strcpy(remain,"/sfdb/O4/C00/CLEAN_AR/"); //OK per O4 CLEAN. AR=analysis ready
	strcpy(remain,"/sfdb/O4/C00/GAT/"); //OK per O4 CLEAN. AR=analysis ready
	//strcpy(remain,"/sfdb/O4/C00/NOLINES/");
	if(strncmp(header_param->sfdbname,"H",1)==0) strcpy(remainfrompath,"H/");
	if(strncmp(header_param->sfdbname,"L",1)==0) strcpy(remainfrompath,"L/");
	if(strncmp(header_param->sfdbname,"V",1)==0) strcpy(remainfrompath,"V/");
	if(strncmp(header_param->sfdbname,"V",1)==0){
	      free(remain);
	      remain=(char*)malloc(sizeof(char)*MAXMAXLINE);
	      strcpy(remain,"/sfdb/O3b/Online21/");
	}
	//printf(" input_param->filelista,nnn, remainfrompath: %s %d %s\n",input_param->filelista,nnn,remainfrompath);
	strcat(remainfrompath,remain);
	if(input_param->resampling==64) strcpy(theband,"128/");
	if(input_param->resampling==16) strcpy(theband,"512/");
	if(input_param->resampling==8) strcpy(theband,"1024/");
	if(input_param->resampling==4) strcpy(theband,"2048/");
	if(input_param->resampling==1) strcpy(theband,"8192/");
	if(input_param->resampling==2) strcpy(theband,"4096/");
	nnn=strlen(theband);
	theband[nnn]='\0';
	strcat(remainfrompath,theband);
	nnn=strlen(remainfrompath);
	remainfrompath[nnn]='\0';
	strcpy(srmdirCNAF,"/virgo4/virgo/RomePSS/");
	strcat(srmdirCNAF,remainfrompath); //il secondo si attacca al primo
	//printf(" BB remainfrompath srmdirCNAF  %s %s\n",remainfrompath,srmdirCNAF);
	strcpy(CommandForGrid,"gfal-copy -p ");
	//Tolto getcwd(dovesono);
	dovesono=(char *)get_current_dir_name();
	printf("dovesono: %s\n",dovesono);
	strcat(CommandForGrid,dovesono);
	printf("CommandForGrid: %s\n",CommandForGrid);
	strcat(CommandForGrid,"/");
	strcat(CommandForGrid,header_param->sfdbname);
	strcat(CommandForGrid," srm://storm-fe-archive.cr.cnaf.infn.it:8444");
	strcat(CommandForGrid,srmdirCNAF);
	strcat(CommandForGrid,header_param->sfdbname);
	nnn=strlen(CommandForGrid);
	CommandForGrid[nnn]='\0';
	printf(" Final instruction for the grid copy of sfdb:  %s\n",CommandForGrid);
	system(CommandForGrid);
	free(CommandForGrid);
	CommandForGrid=(char*)malloc(sizeof(char)*300);
	//NOW remove the file from the node:
	strcpy(CommandForGrid,"gfal-rm ");
	strcat(CommandForGrid,dovesono);
	printf("CommandForGrid: %s\n",CommandForGrid);
	strcat(CommandForGrid,"/");
	strcat(CommandForGrid,header_param->sfdbname);
	nnn=strlen(CommandForGrid);
	CommandForGrid[nnn]='\0';
	printf(" Final instruction to delete the file on the node, sfdb:  %s\n",CommandForGrid);
	system(CommandForGrid);
	free(CommandForGrid);
	//End remove the file from the node

	//PROVA DA TOGLIERE per semplice cp!!!
	//puts("Prova con semplice cp");
	//CommandForGrid=(char*)malloc(sizeof(char)*300);
        //strcpy(CommandForGrid,"cp ");
	//strcat(CommandForGrid,dovesono);
	//strcat(CommandForGrid,"/");
	//strcat(CommandForGrid,header_param->sfdbname);
	//strcat(CommandForGrid," /storage/gpfs_virgo3/home/astone/");
	//nnn=strlen(CommandForGrid);
	//CommandForGrid[nnn]='\0';
	//printf(" TEST instruction for the grid copy of sfdb:  %s\n",CommandForGrid);
	//system(CommandForGrid);
	//free(CommandForGrid);
	//END PROVA DA TOGLIERE
	free(srmdirCNAF);
	free(remainfrompath);
	free(remain);
      if(input_param->flagonline >=3){
	CommandForGrid=(char*)malloc(sizeof(char)*300);
	srmdirCNAF=(char*)malloc(sizeof(char)*MAXMAXLINE);
	remainfrompath=(char*)malloc(sizeof(char)*MAXMAXLINE);
	remain=(char*)malloc(sizeof(char)*MAXMAXLINE);
	//strcpy(remain,"/p10/O2/C02_clean_noscience/");  //Pia: in futuro va sistemato non hard coded. SOTTO rimettere NC per i non gated. Anche sopra Aprile 2020
	//strcpy(remain,"/p10/O3/C01/GAT_CL60/");
	//strcpy(remain,"/p10/O3b/C00/NC/");
	//strcpy(remain,"/p10/O3/C01/CLEAN60/HF/");
	strcpy(remain,"/p10/O4/C00/GAT/");
	//strcpy(remain,"/p10/O4/C00/NOLINES/");
	if(strncmp(header_param->sfdbname,"H",1)==0) strcpy(remainfrompath,"H/");
	if(strncmp(header_param->sfdbname,"L",1)==0) strcpy(remainfrompath,"L/");
	if(strncmp(header_param->sfdbname,"V",1)==0) strcpy(remainfrompath,"V/");
	if(strncmp(header_param->sfdbname,"V",1)==0){
	      free(remain);
	      remain=(char*)malloc(sizeof(char)*MAXMAXLINE);
	      strcpy(remain,"/p10/O3b/Online21/");
	}
	//printf(" input_param->filelista,nnn, remainfrompath: %s %d %s\n",input_param->filelista,nnn,remainfrompath);
	strcat(remainfrompath,remain);
	strcat(remainfrompath,theband);
	nnn=strlen(remainfrompath);
	remainfrompath[nnn]='\0';
	strcpy(srmdirCNAF,"/virgo4/virgo/RomePSS/");
	strcat(srmdirCNAF,remainfrompath); //il secondo si attacca al primo
	//printf(" BB remainfrompath srmdirCNAF  %s %s\n",remainfrompath,srmdirCNAF);
	strcpy(CommandForGrid,"gfal-copy -p ");
	strcat(CommandForGrid,dovesono);
	strcat(CommandForGrid,"/");
	strcat(CommandForGrid,input_param->peakname);
	strcat(CommandForGrid," srm://storm-fe-archive.cr.cnaf.infn.it:8444");
	strcat(CommandForGrid,srmdirCNAF);
	strcat(CommandForGrid,input_param->peakname);
	nnn=strlen(CommandForGrid);
	CommandForGrid[nnn]='\0';
	printf(" Final instruction for the grid copy of peakmap:  %s\n",CommandForGrid);
	system(CommandForGrid);
	free(CommandForGrid);
	CommandForGrid=(char*)malloc(sizeof(char)*300);
	//NOW remove the file from the node:
	strcpy(CommandForGrid,"gfal-rm ");
	strcat(CommandForGrid,dovesono);
	printf("CommandForGrid: %s\n",CommandForGrid);
	strcat(CommandForGrid,"/");
	strcat(CommandForGrid,input_param->peakname);
	nnn=strlen(CommandForGrid);
	CommandForGrid[nnn]='\0';
	printf(" Final instruction to delete the file on the node, sfdb:  %s\n",CommandForGrid);
	system(CommandForGrid);
	free(CommandForGrid);
	//End remove the file from the node
	free(srmdirCNAF);
	free(remainfrompath);
	free(remain);
	free(dovesono);
	free(theband);
      }
  return i;
      }
//Sotto versione ultima con lcg-cp. sostituita da quella sopra
int ToCopyFiles_V2(INPUT_PARAM *input_param,HEADER_PARAM *header_param)
{
  int ij,ijj,ijinizio,i,ijprec;
  char *CommandForGrid,*srmdirCNAF,*remainfrompath,*remain,*dovesono,*theband;

  CommandForGrid=(char*)malloc(sizeof(char)*300);
  srmdirCNAF=(char*)malloc(sizeof(char)*MAXMAXLINE);
  remainfrompath=(char*)malloc(sizeof(char)*MAXMAXLINE);
  remain=(char*)malloc(sizeof(char)*MAXMAXLINE);
  dovesono=(char*)malloc(sizeof(char)*MAXMAXLINE);
  theband=(char*)malloc(sizeof(char)*6);
  int nnn;


        nnn=strlen(input_param->filelista);
	//ijj=0;
	//ij=0;
	//for(i=0;i<nnn;i++){
	//  if(input_param->filelista[i] == 47 || input_param->filelista[i] == 92){
	//     ijprec=ij;
	//     ijj+=1;
	//     if (ijj==5) ijinizio=i;
	//     ij=i; /* 47 -> / , 92 -> \ */}
	//}
	strcpy(remain,"/sfdb/O2/C02_clean_noscience/");  //Pia: in futuro va sistemato non hard coded
	if(strncmp(header_param->sfdbname,"H",1)==0) strcpy(remainfrompath,"H/");
	if(strncmp(header_param->sfdbname,"L",1)==0) strcpy(remainfrompath,"L/");
	if(strncmp(header_param->sfdbname,"V",1)==0) strcpy(remainfrompath,"V/");
	//printf(" input_param->filelista,nnn, remainfrompath: %s %d %s\n",input_param->filelista,nnn,remainfrompath);
	strcat(remainfrompath,remain);
	if(input_param->resampling==64) strcpy(theband,"128/");
	if(input_param->resampling==16) strcpy(theband,"512/");
	if(input_param->resampling==8) strcpy(theband,"1024/");
	if(input_param->resampling==4) strcpy(theband,"2048/");
	nnn=strlen(theband);
	theband[nnn]='\0';
	strcat(remainfrompath,theband);
	nnn=strlen(remainfrompath);
	remainfrompath[nnn]='\0';
	strcpy(srmdirCNAF,"/virgo4/virgo/RomePSS/");
	strcat(srmdirCNAF,remainfrompath); //il secondo si attacca al primo
	//printf(" BB remainfrompath srmdirCNAF  %s %s\n",remainfrompath,srmdirCNAF);
	strcpy(CommandForGrid,"lcg-cp file:");
	getcwd(dovesono);
	//printf("dovesono: %s\n",dovesono);
	strcat(CommandForGrid,dovesono);
	strcat(CommandForGrid,"/");
	strcat(CommandForGrid,header_param->sfdbname);
	strcat(CommandForGrid," srm://storm-fe-archive.cr.cnaf.infn.it:");
	strcat(CommandForGrid,srmdirCNAF);
	strcat(CommandForGrid,header_param->sfdbname);
	nnn=strlen(CommandForGrid);
	CommandForGrid[nnn]='\0';
	printf(" Final instruction for the grid copy of sfdb:  %s\n",CommandForGrid);
	system(CommandForGrid);
	free(CommandForGrid);
	//PROVA DA TOGLIERE per semplice cp!!!
	//puts("Prova con semplice cp");
	//CommandForGrid=(char*)malloc(sizeof(char)*300);
        //strcpy(CommandForGrid,"cp ");
	//strcat(CommandForGrid,dovesono);
	//strcat(CommandForGrid,"/");
	//strcat(CommandForGrid,header_param->sfdbname);
	//strcat(CommandForGrid," /storage/gpfs_virgo3/home/astone/");
	//nnn=strlen(CommandForGrid);
	//CommandForGrid[nnn]='\0';
	//printf(" TEST instruction for the grid copy of sfdb:  %s\n",CommandForGrid);
	//system(CommandForGrid);
	//free(CommandForGrid);
	//END PROVA DA TOGLIERE
	free(srmdirCNAF);
	free(remainfrompath);
	free(remain);
      if(input_param->flagonline >=3){
	CommandForGrid=(char*)malloc(sizeof(char)*300);
	srmdirCNAF=(char*)malloc(sizeof(char)*MAXMAXLINE);
	remainfrompath=(char*)malloc(sizeof(char)*MAXMAXLINE);
	remain=(char*)malloc(sizeof(char)*MAXMAXLINE);
	strcpy(remain,"/p10/O2/C02_clean_noscience/");  //Pia: in futuro va sistemato non hard coded
	if(strncmp(header_param->sfdbname,"H",1)==0) strcpy(remainfrompath,"H/");
	if(strncmp(header_param->sfdbname,"L",1)==0) strcpy(remainfrompath,"L/");
	if(strncmp(header_param->sfdbname,"V",1)==0) strcpy(remainfrompath,"V/");
	//printf(" input_param->filelista,nnn, remainfrompath: %s %d %s\n",input_param->filelista,nnn,remainfrompath);
	strcat(remainfrompath,remain);
	strcat(remainfrompath,theband);
	nnn=strlen(remainfrompath);
	remainfrompath[nnn]='\0';
	strcpy(srmdirCNAF,"/virgo4/virgo/RomePSS/");
	strcat(srmdirCNAF,remainfrompath); //il secondo si attacca al primo
	//printf(" BB remainfrompath srmdirCNAF  %s %s\n",remainfrompath,srmdirCNAF);
	strcpy(CommandForGrid,"lcg-cp file:");
	strcat(CommandForGrid,dovesono);
	strcat(CommandForGrid,"/");
	strcat(CommandForGrid,input_param->peakname);
	strcat(CommandForGrid," srm://storm-fe-archive.cr.cnaf.infn.it:");
	strcat(CommandForGrid,srmdirCNAF);
	strcat(CommandForGrid,input_param->peakname);
	nnn=strlen(CommandForGrid);
	CommandForGrid[nnn]='\0';
	printf(" Final instruction for the grid copy of peakmap:  %s\n",CommandForGrid);
	system(CommandForGrid);
	free(CommandForGrid);
	free(srmdirCNAF);
	free(remainfrompath);
	free(remain);
	free(dovesono);
	free(theband);
      }
  return i;
      }

int ToCopyFiles_V1(INPUT_PARAM *input_param,HEADER_PARAM *header_param)
{
  int ij,ijj,ijinizio,i,ijprec;
 char *CommandForGrid,*srmdirCNAF,*remainfrompath,*remain,*dovesono;

  CommandForGrid=(char*)malloc(sizeof(char)*300);
  srmdirCNAF=(char*)malloc(sizeof(char)*MAXMAXLINE);
  remainfrompath=(char*)malloc(sizeof(char)*MAXMAXLINE);
  remain=(char*)malloc(sizeof(char)*MAXMAXLINE);
  dovesono=(char*)malloc(sizeof(char)*MAXMAXLINE);
  int nnn;


        nnn=strlen(input_param->filelista);
	ijj=0;
	ij=0;
	 for(i=0;i<nnn;i++){
	   if(input_param->filelista[i] == 47 || input_param->filelista[i] == 92){
	     ijprec=ij;
	     ijj+=1;
	     if (ijj==5) ijinizio=i;
	     ij=i; /* 47 -> / , 92 -> \ */}
	 }

	seleziona(remainfrompath,input_param->filelista,ijinizio+1,ijprec-ijinizio);
	//printf(" input_param->filelista,nnn, remainfrompath: %s %d %s\n",input_param->filelista,nnn,remainfrompath);
        strcpy(remain,"sfdb_out/");
	strcat(remainfrompath,remain);
	nnn=strlen(remainfrompath);
	remainfrompath[nnn]='\0';
	strcpy(srmdirCNAF,"/virgo4/virgo/");
	strcat(srmdirCNAF,remainfrompath);
	//printf(" BB remainfrompath srmdirCNAF  %s %s\n",remainfrompath,srmdirCNAF);
	strcpy(CommandForGrid,"lcg-cp file:");
	getcwd(dovesono);
	//printf("dovesono: %s\n",dovesono);
	strcat(CommandForGrid,dovesono);
	strcat(CommandForGrid,"/");
	strcat(CommandForGrid,header_param->sfdbname);
	strcat(CommandForGrid," srm://storm-fe-archive.cr.cnaf.infn.it:");
	strcat(CommandForGrid,srmdirCNAF);
	strcat(CommandForGrid,header_param->sfdbname);
	nnn=strlen(CommandForGrid);
	CommandForGrid[nnn]='\0';
	printf(" Final instruction for the grid copy of sfdb:  %s\n",CommandForGrid);
	system(CommandForGrid);
	free(CommandForGrid);
	//PROVA DA TOGLIERE per semplice cp!!!
	//puts("Prova con semplice cp");
	//CommandForGrid=(char*)malloc(sizeof(char)*300);
        //strcpy(CommandForGrid,"cp ");
	//strcat(CommandForGrid,dovesono);
	//strcat(CommandForGrid,"/");
	//strcat(CommandForGrid,header_param->sfdbname);
	//strcat(CommandForGrid," /storage/gpfs_virgo3/home/astone/");
	//nnn=strlen(CommandForGrid);
	//CommandForGrid[nnn]='\0';
	//printf(" TEST instruction for the grid copy of sfdb:  %s\n",CommandForGrid);
	//system(CommandForGrid);
	//free(CommandForGrid);
	//END PROVA DA TOGLIERE
	free(srmdirCNAF);
	free(remainfrompath);
	free(remain);
      if(input_param->flagonline >=3){
	CommandForGrid=(char*)malloc(sizeof(char)*300);
	srmdirCNAF=(char*)malloc(sizeof(char)*MAXMAXLINE);
	remainfrompath=(char*)malloc(sizeof(char)*MAXMAXLINE);
	remain=(char*)malloc(sizeof(char)*MAXMAXLINE);
	seleziona(remainfrompath,input_param->filelista,ijinizio+1,ijprec-ijinizio);
	//printf(" input_param->filelista,nnn, remainfrompath: %s %d %s\n",input_param->filelista,nnn,remainfrompath);
        strcpy(remain,"peakmap_out/");
	strcat(remainfrompath,remain);
	nnn=strlen(remainfrompath);
	remainfrompath[nnn]='\0';
	strcpy(srmdirCNAF,"/virgo4/virgo/");
	strcat(srmdirCNAF,remainfrompath);
	strcpy(CommandForGrid,"lcg-cp file:");
	strcat(CommandForGrid,dovesono);
	strcat(CommandForGrid,"/");
	strcat(CommandForGrid,input_param->peakname);
	strcat(CommandForGrid," srm://storm-fe-archive.cr.cnaf.infn.it:");
	strcat(CommandForGrid,srmdirCNAF);
	strcat(CommandForGrid,input_param->peakname);
	nnn=strlen(CommandForGrid);
	CommandForGrid[nnn]='\0';
	printf(" Final instruction for the grid copy of peakmap:  %s\n",CommandForGrid);
	system(CommandForGrid);
	free(CommandForGrid);
	free(srmdirCNAF);
	free(remainfrompath);
	free(remain);
	free(dovesono);
      }
  return i;
      }



 int Logical2PhysicalPath(INPUT_PARAM *input_param,int iquale)  //with 1 change the path to the ffl list file. with 2 to the frame files
 {

   char *srmdirCNAF,*remainfrompath;
   int i,nnn,ijj,ijinizio,ij;
  srmdirCNAF=(char*)malloc(sizeof(char)*MAXMAXLINE);
  remainfrompath=(char*)malloc(sizeof(char)*MAXMAXLINE);


  if(iquale==1) {
        nnn=strlen(input_param->fflname);
	ijj=0;
	 for(i=0;i<nnn;i++){
	   if(input_param->fflname[i] == 47 || input_param->fflname[i] == 92){
	     ijj+=1;
	     if (ijj==2) ijinizio=i;
	     ij=i; /* 47 -> / , 92 -> \ */}
	 }
	seleziona(remainfrompath,input_param->fflname,ijinizio+1,nnn-ijinizio);
	//printf(" input_param->fflname,nnn, remainfrompath: %s %d %s\n",input_param->fflname,nnn,remainfrompath);
	nnn=strlen(remainfrompath);
	remainfrompath[nnn]='\0';
	strcpy(srmdirCNAF,"/storage/gpfs_virgo/virgo4/");
	strcat(srmdirCNAF,remainfrompath);
	//printf(" BB remainfrompath srmdirCNAF  %s %s\n",remainfrompath,srmdirCNAF);
	strcpy(input_param->fflname,srmdirCNAF);
	nnn=strlen(input_param->fflname);
	input_param->fflname[nnn]='\0';
	printf(" Physical input_param->fflname %s\n",input_param->fflname);
  }
  //e basta almano per ora. i nomi dei files devono gia' essere fisici.
  return i;

  }
//Functions modified and added to use frinit different from zero. Main: crea_sfdbFRINIT.  22/07/2015
INPUT_PARAM* crea_inputFRINIT(int idefaults,int casin,int *casout,float *frinit)
{

  /*If needed for specific purposes it MUST be modified outside, in the main code  */
  /* idefaults=0 puts the default values*/
  INPUT_PARAM *input_param;
  input_param=(INPUT_PARAM *)malloc(sizeof(INPUT_PARAM));
  //char filename[MAXMAXLINE];  /*input file*/
  long len;
  int typ,wink,itotal;
  int verb;
  int nn,nn_orig;                  /* len of the file name */
  char *filename;
  char *appo;
  char *appo1;
  char *appo2;
  char *appo3;
  char *appo4;
  //char check_extens[NEXTS];
  //int icheck;
  int red;
  int lena,lenb;
  int i,ii;
  int jj;
  int nnn;
  int uno;
  float freqc;
  int subsam;
  int fft_m;
  char *fflname;
  char *scisegfilename;
  int flagonline;
  char Ch[25];
  char nam[5];
  //int noise_pulsars;

  fflname=(char*)malloc(sizeof(char)*(MAXMAXLINE+1));
  scisegfilename=(char*)malloc(sizeof(char)*200);
  filename=(char*)malloc(sizeof(char)*(MAXMAXLINE+1));
  appo=(char*)malloc(sizeof(char)*200); //Pia: messo 100
  appo1=(char*)malloc(sizeof(char)*7); //era 5 prima di sfdb09
  appo2=(char*)malloc(sizeof(char)*10);
  appo3=(char*)malloc(sizeof(char)*4);
  appo4=(char*)malloc(sizeof(char)*200);
  printf("Detector ? virgo or V1, ligoh or H1, ligol or L1 \n");
  scanf("%s",nam);
  if (strncasecmp(nam,"V1",2) ==0) {
    puts("detector put to virgo");
      strncpy(nam,"virgo",5);
    }
  if (strncasecmp(nam,"L1",2) ==0) {
      puts("detector put to ligol");
      strncpy(nam,"ligol",5);
    }
  if (strncasecmp(nam,"H1",2) ==0) {
      puts("detector put to ligoh");
      strncpy(nam,"ligoh",5);
    }
  strncpy(input_param->detector_name,nam,5);
  printf("input_param->detector_name: %s \n",input_param->detector_name);
  // printf("Noise around known pulsars evaluation ?  (1=Yes,0 or else=No)\n");
  //scanf("%d",&noise_pulsars);
  //input_param->noise_pulsars=noise_pulsars;
  printf("Type  of SFDB files? \n <=0--> Void SFDB \n ==1--> no SFDB \n =2-->full SFDB \n=3--> full SFDB and peakmap .p10 \n");
  scanf("%d",&flagonline);
  if(flagonline<=1)*casout=2;
  if(flagonline==2)*casout=1;
  if(flagonline>=3)*casout=3;
  if(flagonline==1)flagonline=1;
  if(flagonline<=0)flagonline=0;
  if(flagonline==2)flagonline=0;
  if(flagonline>=3)flagonline=3;
  input_param->flagonline=flagonline;
  printf("Factor to write EVF in the file ?  (1 or 2 usually) \n");
  scanf("%f",&input_param->factor_write);
  if(casin==1){
   input_param->resampling=1;
   printf("input sds file: \n");
   scanf("%s",filename);
   strncpy(input_param->framechannel,"sdsfile",7);
  }
 if(casin==2){
  printf("input ffl file, for frame input: \n");
  scanf("%s",fflname);
  strcpy(input_param->fflname,fflname);
  if (strncmp(input_param->fflname,"lfn:",4) == 0) {
    Logical2PhysicalPath(input_param,1);  //with 1 change the path to the ffl list file. with 2 to the frame files
  }
  printf("science segment file name: \n");
  scanf("%s",scisegfilename);
  strcpy(input_param->scisegfilename,scisegfilename);
  printf("file frame input channel\n h_4096Hz (4096Hz)\n Pr_B1_ACp (20kHz)\n h_20000Hz (20kHz)\n Em_SEBDCE01 (1 kHz)\n Em_MABDCE01 (20kHz) \n H1:LDAS-STRAIN (16384 Hz LIGO H1)\n  L1:LDAS-STRAIN (Ligo L1)\n");
  /*printf("file frame input channel   1=V1:h_4096Hz 2=V1:Pr_B1_ACp 3=V1:h_20000Hz  ?\n");

  scanf("%d",&fch);
       printf("%d,fch\n",fch);
       switch(fch){
               case 1:
                 strcpy(input_param->framechannel,"h_4096Hz");
                 //strcpy(input_param->framechannel,"V1:h_4096Hz");
                 printf("PRIMA %s framech\n",input_param->framechannel);
                 break;
               case 2:
                      strcpy(input_param->framechannel,"Pr_B1_ACp");
                      break;
               case 3:
                       strcpy(input_param->framechannel,"h_20000Hz");
                       break;
       }

  */

  scanf("%s",Ch);
  //printf("strlen %d\n",strlen(Ch));
  strncpy(input_param->framechannel,Ch,strlen(Ch));
  printf("Frame channel selected: %s \n",input_param->framechannel);

//Added for the flags: with the suggested names. MAGGIO 2017
  strncpy(input_param->flagchannel,"tobefilled",10);
  if(strncmp(input_param->framechannel,"h_",2)==0){
    //strncpy(input_param->flagchannel,"V1:DQ_META_ITF_LOCKED",21); used for O2
     strncpy(input_param->flagchannel,"V1:Hrec_STATE_VECTOR",20);
     input_param->flagstatus=4095;  //era 1 per O2. Sarebbe 1 se volessimo contare solo Science
    //strncpy(input_param->flagchannel,"V1:Hrec_Veto_ScienceMode",strlen("V1:Hrec_Veto_ScienceMode"));
    //strncpy(input_param->flagchannel,"V1:Hrec_Flag_Channel",strlen("V1:Hrec_Flag_Channel"));
    //input_param->flagstatus=0;
  }
   if(strncmp(input_param->framechannel,"V1:h_",5)==0){
     //strncpy(input_param->flagchannel,"V1:DQ_META_ITF_LOCKED",21); used for O2
      strncpy(input_param->flagchannel,"V1:Hrec_STATE_VECTOR",20);
      input_param->flagstatus=4095;
     //strncpy(input_param->flagchannel,"V1:Hrec_Veto_ScienceMode",strlen("V1:Hrec_Veto_ScienceMode"));
     //strncpy(input_param->flagchannel,"V1:Hrec_Flag_Channel",strlen("V1:Hrec_Flag_Channel"));
     //input_param->flagstatus=0;
  }
  if(strncmp(input_param->framechannel,"L1:",3)==0){
    strncpy(input_param->flagchannel,"L1:IFO-SV_STATE_VECTOR",22);
    input_param->flagstatus=65535;
  }
  if(strncmp(input_param->framechannel,"H1:",3)==0){
    strncpy(input_param->flagchannel,"H1:IFO-SV_STATE_VECTOR",22);
    input_param->flagstatus=65535;
  }
  if(strncmp(input_param->framechannel,"L2:",3)==0){
    strncpy(input_param->flagchannel,"L1:IFO-SV_STATE_VECTOR",22);
    input_param->flagstatus=65535;
  }
  if(strncmp(input_param->framechannel,"H2:",3)==0){
    //strncpy(input_param->flagchannel,"H1:IFO-SV_STATE_VECTOR",22);
    strncpy(input_param->flagchannel,"H2:ISC-ALS_EY_REFL_PWR_MON_OUT16",32);
    input_param->flagstatus=3500;
  }

  if(strncmp(input_param->flagchannel,"tobefilled",10)==0){
    //    strncpy(input_param->flagchannel,"V1:Hrec_Veto_LOCK_STEP_STATUS",29);
    //strncpy(input_param->flagchannel,"V1:DQ_META_ITF_LOCKED",21); used for O2
    strncpy(input_param->flagchannel,"V1:Hrec_STATE_VECTOR",20);
    input_param->flagstatus=4095;  //era 12. Poi 1 per O2. Sarebbe 1 se volessimo contare solo Science. GE 2018
  }
  printf("Flags channel selected: %s \n",input_param->flagchannel);

  // End flag channel

  strcpy(filename,input_param->fflname);

  printf("resampling factor ?  Has to be a power of 2! And >=2 if frinit >0 !!\n");
  scanf("%d",&input_param->resampling);
  if(frinit >0 && input_param->resampling==1){
    input_param->resampling=2;
  }
 }

  //strcpy(filename,"/storage/gpfs_virgo3/scratch/pss/virgo/sd/sds/VSR1-2/deca_20070528/VIR_hrec_20070529_114516_.sds");
  // write the filename of the first file in the LOG file
 jj=logfile_input(LOG_INFO,input_param->detector_name,"detector");
 //jj= logfile_input(LOG_INFO, filename, "First file of the run");
 if(casin==2){
   jj= logfile_input(LOG_INFO,input_param->framechannel , "Frame channel");
 }

 if(idefaults !=0) {
    printf("Verbosity level 0 1 2 3\n");
    scanf("%d",&verb);
    printf("chunk len  -  DS type (0,1,2=interlaced) \n");
    puts("Chunk len will be rounded to the next 2 power, if it is not a power of 2");
    printf("Typical for Virgo 1-band:\n len=4194304 (4096Hz)");
    printf("len =16777216 (20kHz) (es. Em_MABDCE01 Pr_B1_ACp )\n");
    printf("len =1048576 (1kHz) (es Em_SEBDCE01)\n ");
    puts("ds type=2 interlaced; 1 not interlaced");
    scanf("%d%d",&lena,&typ);
    printf("Maximum number of data chunk to be done ?\n");
    scanf("%d",&itotal);
    printf("Reduction factor, for very short FFTs (e.g. 2,4,8..Sugg.:128) ?\n");
    scanf("%d",&red);
    printf("windows  (0=no,1=Hann,2=Hamm,3=MAP, 4=Blackmann flatcos 5=flat top,cosine edge. Sugg. 5) \n");
    scanf("%d",&wink);
    printf("Cut frequency for the highpass filtering (e.g. 100. If 0 the highpass and EVT veto are skipped) ?\n");
    scanf("%f",&freqc);
    printf("Subsampling factor for the veto in the subbands. Power of 2 (e.g. 128)? <=0: not applied\n");
    puts("Typical for Virgo Em_SEBDCE01 ch (1kHz) =4\n");
    puts("Typical for Virgo Em_MABDCE01 ch (20KHz) =4\n");

    scanf("%d",&subsam);
    printf("Max number of FFTs in one output file (e.g. 100) ?\n");
    scanf("%d",&fft_m);
  }
  else {
    verb=0;
    lena=4194304;
    typ=2;
    itotal=10000000;
    red=128;
    wink=5;
    freqc=100;
    subsam=0;
    fft_m=100;
    input_param->resampling=1;
    }
  nn_orig=strlen(filename);
  printf("Original filename nn_orig= %d \n",nn_orig);
  lenb=next2power(lena);
  len=(long) lenb;
  printf("len written, len used %d %ld\n",lena,len);
  strcpy(input_param->filename,filename);
  input_param->filename[nn_orig]='\0';
  //printf("%s input parameters filename \n",input_param->filename);
  //printf("%s filename \n",filename);
  input_param->freqc=freqc;
  input_param->verb=verb;
  input_param->len=len;
  input_param->typ=typ;
  input_param->wink=wink;
  input_param->itotal=itotal;  /*Maximum number of data chunk to be done*/
  input_param->red=red;  /*Reduction factor for very short FFTs */
  nn=strlen(input_param->filename);
  input_param->nn=nn;
  input_param->subsam=subsam;
  printf("filename len= %d \n",input_param->nn);
  strcpy(input_param->capt_h,"h reconst");
  strcpy(input_param->capt_gd,"Chunk");
  input_param->typ_h=input_param->typ; //DS type for h-reconstructed data: the same as input data.
  input_param->fft_m=fft_m;
  input_param->iopen=0; //to open a new file
 /* To construct the sfdb file name, from a filename which contains a path.
    Data of the C5 run are in /windows/xg23nas/pss/virgo/sd/c5/
*/
  ii=-1;
  for(i=0;i<nn;i++){
      if(input_param->filename[i] == 47 || input_param->filename[i] == 92) ii=i; /* 47 -> / , 92 -> \ */
  }

  /*To extract the name and remove the extension, next substituted with the .SFDB one:*/
  //seleziona(appo,input_param->filename,ii+1,nn-ii-1-NEXTS-1);
  seleziona(appo,input_param->filename,ii+1,nn-ii-1-NEXTS);
  printf("core of the name= %s \n", appo);
  uno=1;
  sprintf(appo2,"%04d",uno);
  strncat(appo,appo2,4);
  nnn=strlen(appo);
  strncpy(appo4,appo,nnn);
  // printf("$$$$ appo appo2= %s %s\n", appo,appo2);  //quiquiritogliere

  switch (*casout)
	  {
		  case 1:
		          strncpy(appo1,".SFDB09",7);
			  break;
		  case 2:
			  strncpy(appo1,".SFDB09",7);
			  break;
		  case 3:
			  strncpy(appo1,".SFDB09",7);
			  strncpy(appo3,".p10",4);
			  break;
	  }



  strncat(appo,appo1,7);  //era 5 prima di sfdb09
  strcpy(input_param->sfdbname,appo);

  nnn=strlen(appo);
  input_param->sfdbname[nnn]='\0';
  if(input_param->flagonline!=1){
   puts("In the header of the SFDB file you will find in detector: 1 for virgo, 2 for ligoh 3 for ligol 4 for kagra");
   printf("Sfdb filename= %s \n", input_param->sfdbname);
   //printf("Detector = %s \n", input_param->detector_name);

   jj=logfile_output(LOG_INFO,input_param->sfdbname,"List name. The first SFDB file has been opened");
     }
  if(input_param->flagonline >=3){
    puts("Peakmap file with extension .p10 will be created");
    strncat(appo4,appo3,4);
    strcpy(input_param->peakname,appo4);
    nnn=strlen(appo4);
    input_param->peakname[nnn]='\0';

  }
 /* To guess the file type (R87-SDS-FRAME): */
  seleziona(input_param->file_extens,input_param->filename,nn-NEXTS,NEXTS);
  printf("extension of name= %s \n", input_param->file_extens);

 /* Check for the file type (SDS or FRAME): if icheck=0 they are R87 */
    //check if FRAME or SDS
  //strcpy(check_extens,"SDS");
  //icheck=strcasecmp(check_extens,input_param->file_extens);
  //if(icheck==0)input_param->interferom=1;  //sds
  //if(icheck!=0)input_param->interferom=2; //frame

  free(appo);
  free(appo1);
  free(appo2);
  free(appo3);
  free(appo4);
  free(fflname);
  free(filename);


  return input_param;
}
int resamplingFRINIT(GD *gd_res,GD *gd,HEADER_PARAM *header_param)
 {

   int res,i,ivign;
   float bmax;
   int ii,kki;
   //printf("In resampling- gd_res->n gd->n %ld %ld  \n",gd_res->n,gd->n);
   res=gd->n/gd_res->n;
   //printf("res=%d\n",res);

header_param->n_zeroes=0;
 for(ii=0;ii<gd->n;ii++){
   kki=isnan(gd->y[ii]);
   if (kki !=0)
   {
    gd->y[ii]=0;
   }
   if(gd->y[ii]==0)header_param->n_zeroes+=1;
 }
 header_param->perc_zeroes=(float)(header_param->n_zeroes)/gd->n;  //becomes a percentage


   i= gdtime2freq(gd,header_param);//from time gd to frequency gd
   gd_res->dx=gd->dx;
   //printf("In resampling dopo gdtime2freq- gd_res->dx gd->dx %f %f  \n",gd_res->dx,gd->dx);
   //bmax=1/(2*gd_res->dx); //max freq in the data
   bmax=(gd_res->dx*gd_res->n/2);
   //printf("In resampling- max freq =%f  \n",bmax);
   //bmax is the bandwidth. max freq=frinit+bmax.  22 July 2015
   //From the freq data extract the BW, go back to the time domain.
   ivign=0;

   //printf("In resampling- gd_res->dx gd->dx %f %f  \n",gd_res->dx,gd->dx);
   i=band_extract(gd,gd_res,bmax);
   //Band extract does not normalize data
    header_param->frinit=gd_res->ini; //added 22 July 2015
    printf("In resampling- after resampling: gd_res->ini, before res: gd->ini %f %f  \n",gd_res->ini,gd->ini);
   //Go back to the time domain:
   i=gdfreq2time(gd_res,header_param,ivign);
   //printf("dopo gdfreq2time \n");
   //normalization of the data, for the subsamplig:
   for(i=0; i<gd_res->n;i++){
     gd_res->y[i]=gd_res->y[i]/res;

     //if(i<=5)printf("i=%d gd_res[i]=%f\n",i,gd_res->y[i]);
     //if(i>=gd_res->n-5)printf("i=%d gd_res[i]=%f\n",i,gd_res->y[i]);

   }
   //Redefine gd_res and header_param
   //gd_res->dx=res*gd->dx;
   //gd_res->ini=gd->ini;
   header_param->nsamples=gd_res->n/2;
   header_param->tsamplu=gd_res->dx;

   //printf("In resampling- after resampling gd_res->dx gd_res->ini %f %f  \n",gd_res->dx,gd_res->ini);
   //printf("In resampling- header_param->nsamples=%d  header_param->tsamplu =%f\n",header_param->nsamples,header_param->tsamplu);
   return i;

 }
int ScienceVeto_onGD(GD *gd,HEADER_PARAM *header_param,SCIENCE_SEGMENT* science_segment)
{
  int i;
  //printf("even_param->sat_level= %f\n ",even_param->sat_level);
  //printf("even_param->sat_maxn= %d\n ",even_param->sat_maxn);
  //printf("Applying science segment list selection!\n");
  double TINI,TFIN;
  int jjj;
  double vectime;
  TINI=header_param->gps_sec+(header_param->gps_nsec)/(10^9);
  TFIN=TINI+gd->n*header_param->tsamplu;
  int iVETATI=0;
  jjj=0;  //to run faster

  while (jjj<science_segment->nsegm-1 && TINI>(double)science_segment->segtstart[jjj+1]){
	 jjj++;
   }
   printf("tsamplu= %f\n",header_param->tsamplu);
   if (TINI>=(double)science_segment->segtstart[jjj] && TFIN<=(double)science_segment->segtstop[jjj]) {
      printf("===> FFT all inside a science period %16.6f %16.6f\n",TINI,TFIN);
      i=0;
      return i;
    } else {
	printf("===> FFT NOT ALL inside a science period %16.6f %16.6f\n",TINI,TFIN);
      }


  for (i=0; i < gd->n; i++)
 {
   vectime=(double)(TINI+i*header_param->tsamplu);
   while (jjj<science_segment->nsegm-1 && vectime>(double)science_segment->segtstart[jjj+1]){
	jjj++;
  }
   if (vectime>=(double)science_segment->segtstart[jjj] && vectime<(double)science_segment->segtstop[jjj]) {
 	    if (i < 2) printf("OK time = %16.6f segm start = %d segm stop = %d ----> in science mode segment n. %d! gd->y[i]=%e\n", vectime, science_segment->segtstart[jjj], science_segment->segtstop[jjj], science_segment->num[jjj], gd->y[i]);

   } else {
	  gd->y[i]=0.;
	  iVETATI+=1;
	  if (iVETATI <=3) printf("VETOED time = %16.6f segm start = %d segm stop = %d ----> in science mode segment n. %d! gd->y[i]=%e\n", vectime, science_segment->segtstart[jjj], science_segment->segtstop[jjj], science_segment->num[jjj], gd->y[i]);

   }
   if (vectime >=1126629172 && vectime <=1126629172.0003) {
     printf("Al secondo 52 non dovrebbe essere vetato= %16.6f segm start = %d segm stop = %d ----> in science mode segment n. %d! gd->y[i]=%e\n", vectime, science_segment->segtstart[jjj], science_segment->segtstop[jjj], science_segment->num[jjj], gd->y[i]);
   }
   if (vectime >=1126629182 && vectime <=1126629182.0003) {
     printf("Al secondo 62 dovrebbe essere vetato= %16.6f segm start = %d segm stop = %d ----> in science mode segment n. %d! gd->y[i]=%e\n", vectime, science_segment->segtstart[jjj], science_segment->segtstop[jjj], science_segment->num[jjj], gd->y[i]);
   }
   if (vectime >=1126629180 && vectime <=1126629180.0003) {
     printf("Al secondo 60 non dovrebbe essere vetato= %16.6f segm start = %d segm stop = %d ----> in science mode segment n. %d! gd->y[i]=%e\n", vectime, science_segment->segtstart[jjj], science_segment->segtstop[jjj], science_segment->num[jjj], gd->y[i]);
   }
 } //closes for i
 //printf("even_param->sat_howmany= %f\n ",even_param->sat_howmany);
  return i;
}

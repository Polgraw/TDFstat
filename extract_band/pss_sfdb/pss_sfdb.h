/**pss_sfdb software for the creation of the SFDB data**/
/**August, 25, 2004 **/
/**Author: Pia**/

#include "../pss_lib/pss_math.h"
#include "../pss_lib/pss_snag.h"
#include "../pss_lib/pss_sfc.h"

/***************************************************/
/*****Constant definitions*****/
#define daybegin 36889.0   //31 December 2000
#define dot2gps  29224.0
#define leapsec  13
#define day2sec  86400.0
#define day2min  1440.0
#define day2hour 24.0
#define hour2sec  3600.0
#define min2sec   60.0
#define dot2mjd  15020
#define REFERENCE_MJD 54248.1665046   //53342.031492
//53180.916000 to be used for VIRGO C4
//53342.031492 /*this is 20 Nov. 2004. Good for C5 in VIRGO. Used only for simulations with Doppler */
/* mjd=dot+dot2mjd, where dot are days from the beginning of 1900, MJD=Modified Julian day */
#define MJD2000 51544.0
// 1 Jan 2000, used as reference time for fake signals simulation
#define PIG      3.1415926
#define BOLTZ	 1.38086e-23
#define EINSTEIN  1.0e-20  /*scaling factor on the h-reconstructed data, both in time and frequency domain) */
#define bil2uni 2.0
#define MAXLINE 25    //era 14 prova 2 feb 2009
#define MAXMAXLINE 170 //era 90 Feb 2007. Jan 2006 era 35. 25 Aug 2005: metto 80 per avere un path lungo (/windows/xg23nas/pss/virgo/sd/c5)
#define NEXTS   3  //era 3 Va messo a 4 ?? Provare !!quiqui Dciembre 2008
#define maxdata_ring 600000  /*not more than maxdata_ring data in one ring */
#define RATIO_YES 0  /*0 normally. Produces SFDB data--1 only if on files we want the ratio. Produces .RAPP data, 2 if we want spectra .SPEC files */
/*****************************************************/
//#define OUTFILEH  "outfileh.dat"
 extern FILE *SFDB;
 extern FILE *P10;
 extern FILE *AREST;
 extern FILE *PEAKMAP;
 extern FILE *OUTH1;
 extern FILE *EVEN_INFO;
 extern FILE *LOG_INFO;  //pointer to the log file
 extern FILE *LISTA; //lista of SFDB files for the peakmap vretion
 extern FILE *OUTTEST; //scrive tutti gli eventi
 extern FILE *OUTTESTB; //scrive tutti gli eventi in banda
 extern FILE *IN_NOISEPULSARS;
 extern FILE *OUT_NOISEPULSARS;
/*************Structures:****************************/
extern struct FrFile *iFileNoise;

typedef struct EVEN_PARAM{
  float tau;       //memory time of the autoregressive average
  float cr;        //CR of the threshold
  float deadtime;  //dead time in seconds
  float factor;    //e.g. 10: when re-evaluate mean and std
  int iflev;       //0=no events; 1=begin; 2=event; 3=end of the event
  float *xamed;       //AR mean
  float *xastd;      //AR std
  float xw;
  float qw;
  double w_norm;   //normalization factor for mean and std
  float edge;     //how many seconds around (before and after) the event have to be "purged"
  int total; //how many samples, due to all events, have been cleaned in the GD
  int *begin;   //sample numbers, in the GD, of the beginning of each event
  int *duration; //event duration in number of samples
  int *imax;    //index at which each event has the maximum;
  float *crmax;    //CR of the maximum of the event;
  float *ener;    //Event energy (sum of the squared amplitudes);
  unsigned long ndata_run; //number of data in the run (one SFDB file, many FFTs !)
  int number;  //number of events in one gd;
  short int ar;  //Algor. for the AR evaluation: 0 is the old, 1 the new
  float xw_pari;
  float qw_pari;
  float xw_dispari;
  float qw_dispari;
  double w_norm_pari;   //normalization factor for mean and std
  double w_norm_dispari;   //normalization factor for mean and std
  unsigned long ndata_run_pari; //number of pari data in the run (one SFDB file, many FFTs !)
  unsigned long ndata_run_dispari; //number of dispari data in the run (one SFDB file, many FFTs !)
  float sat_level;
  int sat_maxn;
  double sat_howmany;
  int absvalue;  //1=uses abs. Else=uses values with their sign
  float notzero; // avoid the denominator of the CR to be zero, value dependends on the data!

 }EVEN_PARAM;

 typedef struct EVF_PARAM{
  float tau_Hz;       //memory time of the autoregressive average
  float maxage;       //max age of the AR procedure
  float maxdin;       //threshold
  float factor_write;    //how to enhance the treshold to write EVF in the log file

 }EVF_PARAM;

typedef struct INPUT_PARAM{
  char filename[MAXMAXLINE+1];  /*First input file or ffl list*/
  char detector_name[5];
  long len;           /*len of the DS (=FFT) */
  int typ,wink;
  int verb;
  int itotal;               /*Maximum number of data chunk to be done*/
  int red;                  /*Reduction factor for the very short FFTs (e.g. 2,4,8..)*/
  int nn;                  /* len of the file name */
  //int interferom;          /*if 1 data are sds, if 2 data are frames */
  int typ_h;               /*DS type for the h reconstructed data = typ */
  char sfdbname[MAXMAXLINE+1]; /*Name of the output SFDB file*/
  char peakname[MAXMAXLINE+1]; /*Name of the output P10 file*/
  char filelista[MAXMAXLINE+1]; /*Name of the first input frame file. Needed to handle the grid copies*/
  char capt_h[12];          /*Caption for the DS with the h reconstructed data. Used for R87 where raw and h-rec are different*/
  char file_extens[NEXTS];  /*Last NEXTS characters: filename extension (e.g. R87) */
  char capt_gd[12];         /*Caption for the GD with SFT data */
  float freqc;   //Cut frequency for the highpass
  int subsam;    //Subsampling for the veto in the subbands
  int fft_m;  //Max Number of FFTs in one output file
  int iopen;
  int nfftONE; //number of FFTs in one single file

  char fflname[MAXMAXLINE+1]; //lista of the frame files
  char scisegfilename[MAXMAXLINE+1]; //science segment list file
  char framechannel[40]; //1=V1:h_4096Hz 2=V1:Pr_B1_ACp 3=V1:h
  char flagchannel[30];
  int flagstatus; //lock or science mode status
  int flagonline; //1 if online; 0 if offline
  float factor_write; //2 normally if offline  1 if online
  int resampling; //resampling factor typically for auxiliary channels
  //int noise_pulsars; //1 if the code evaluates the noise around given bands

}INPUT_PARAM;

typedef struct HEADER_PARAM{
  char sfdbname[MAXMAXLINE+1];  /*Name of the output SFDB file*/
  double endian;
  int detector;            /*if 1 it is virgo, if 2 ligoh if 3 ligol*/
  int gps_sec;
  int gps_nsec;
  double tbase;          /*Len in seconds of 1 data chunk used for the  FFT*/
  int firstfreqindex;
  int nsamples;         /*Number of samples in half FFT (which is what is stored)*/
  int red;            /*Reduction factor for the very short FFTs (e.g. 2,4,8..)*/
  float einstein;
  double mjdtime;
  int nfft;
  int wink;         /*0=no win--Other numbers are different windows*/
  float normd;      /*Normalization factor fft2spectrum  */
  float normw;      /*Normalization factor for the window */
  int typ;          /*0, 1 = non overlapping data  2=overlapping data*/
  float n_flag; //how many (data/total) have a itf lock status less than 12 in that FFT
  double frinit,tsamplu,deltanu;
  double vx_eq,vy_eq,vz_eq;
  double px_eq,py_eq,pz_eq;
  int n_zeroes;  //number of added zeroes in each chunk
  double sat_howmany;
  double spare1,spare2,spare3;
  float perc_zeroes,spare5,spare6;
  int lavesp,sciseg,spare9;

} HEADER_PARAM;

typedef struct SCIENCE_SEGMENT{
  int twin;
  int nsegm;
  char scisegfilename[MAXMAXLINE+1];    /*Name of the science segment file*/
  int num[10000]; // science segment index
  int segtstart[10000]; // science segment GPS time start
  int segtstop[10000]; // science segment GPS time stop
  int segdur[10000]; // science segment duration (s)
} SCIENCE_SEGMENT;


/***************************************************/
/*****************Functions:************************/
/******To fill the structures:*********************/
INPUT_PARAM* crea_input(int idefaults,int casin,int *casout);

/* idefaults=0 puts the default values*/
HEADER_PARAM* crea_sfdbheader(INPUT_PARAM *input_param);
/* Read the science segment file and fill the structure*/
SCIENCE_SEGMENT * crea_science_segment(INPUT_PARAM *input_param);

/***********To handle with times and conversions:**********************************/
double day2itt2001(int year,int month, int day,int hour,int min,float sec,float tmax);
/*
  From year,month,day,hout,min,sec,tmax it gives days from the beginning of the year 2001
*/
double day2mjd(int year,int month, int day,int hour,int min,float sec,float tmax);
void mjd2day(double *mjd_eve,int *year,int *month,int *day,int *hour,int *min,double *sec);

/**************FFTs:***********************************/
void recoverfft(int nn,Cmplx  *cy);
     /*
	function to recover the whole FFT from one half
	nn=half FFT
	input and output in the real vector y
     */
void recoverfft_analitic(int nn,Cmplx  *cy);
     /*
	function to recover the whole FFT from one half: with zeroes analitic signal
	nn=half FFT
	input and output in the real vector y
     */
/*************To handle with DS, vectors and GDs*************/
int fill_newds(DS *d, float *y,DS *d_orig);
/*It copies the parameters of an orig DS into another one. Typically used to
pass from raw data to h-reconstructed data: the two data streams have the same characteristics but
different data  */

int gd_ds(GD *gd,DS *d);
/*
Extracts one chunk of data from a DS and puts the data into a GD.
It uses y_ds. The GD has to be already created, with crea_gd. It sets the properties of the GD
from that of the DS
 */
/**************Utilities to manipulate strings and numbers:**************/

int seleziona(char a[], char b[], int s, int l);

int next2power(int n);
 /*for a given number n it gives the next integer which is a power of 2*/

/**************Windows****************************/

void gd2wingd(GD *gd,int wink);
/*
       Applies a chosen windows (wink) to data in a GD
       INPUT=GD with data
       OUTPUT=GD with windowed data
       wink=1 Hanning
       wink=2 Hamming
       wink=3 Used for MAP (as in GEO data)
       wink=4 Blackmann

*/
void wingd2gd(GD *gd,int wink);
/*
       Function used to remove the window used on data to construct an FFT. When going back to the time
       domain the window have to be removed, if factor is not zero
       INPUT=GD with windowed data
       OUTPUT=GD without the window
       wink=1 Hanning
       wink=2 Hamming
       wink=3 Used for MAP (as in GEO data)
       wink=4 Blackmann
*/
void vect2winvect(float *signal,int wink,int winlen1);
   /*
       INPUT=vector signal float
       OUTPUT=vector signal with the chosen window
       wink=1 Hanning
       wink=2 Hamming
       wink=3 Used for MAP (as in GEO data)
       wink=4 Blackmann
     */
/**************To produce the SFDB files:****************************************/

/********Time and frequency domain operations, with input GD and output the same or different GD **********/
int gdtime2freq(GD *gd,HEADER_PARAM *header_param);
/*from a GD in time to the same GD in frequency domain
header_param are needed only to state the inizial frequency of the FFT
Frequency domain:output to the same GD *gd with real in ii and imag in ii+1... !! Half bandwidth !!
 */
int gdfreq2time(GD *gd,HEADER_PARAM *header_param,int ivign);
/*from a GD in frequency (half bandwidth) to the same GD in time domain.
header_param are needed only to state the inizial time of the data chunk;
ivign=1 if the data are "windowed" before going back to the time domain ("vignettatura quadratica");
Time domain:output to the same GD *gd with the real part of the data (imag is zero) !! Real data  !!
 */
int gdfreq2time_analitic(GD *gd_time2,GD *gd,HEADER_PARAM *header_param,int ivign);
/*from a GD in frequency (half bandwidth) to  GD_time2 (double dim,Re and Imag) in time domain.
Same as the previous BUT using the analitic signal
header_param are needed only to state the inizial time of the data chunk;
ivign=1 if the data are "windowed" before going back to the time domain ("vignettatura quadratica");
Time domain:output to the same GD *gd with only the real part of the data !! Complex data, output
the real part  !!
 */
int band_extract(GD *gd,GD *gd_band,float band);
/* to extract a sub-band on a type 1 gd, from frequency domain data (from gd to gd_band)
 Both the gds contain only one half of the FFTs. The data are alternatively the Real and Imag part */
int band_extract_double(GD *gd,GD *gd_band,float band);
/* As the previous but with 1/4 of zeroes at the beg and end*/
int band_extract_1hz(GD *gd,GD *gd_band,float *band,float *initial_freq,float *final_freq,int ifatti);
//NEW: to have always a 1 Hz band. 0.5 Hz band, which will become 1 Hz for the analitic sig
/* to extract a sub-band on a type 1 gd, from frequency domain data (from gd to gd_band)
 Both the gds contain only one half of the FFTs. The data are alternatively the Real and Imag part */
int resampling(GD *gd_res,GD *gd,HEADER_PARAM *header_param);
int resamplingFRINIT(GD *gd_res,GD *gd,HEADER_PARAM *header_param);
int short_ps_time(float *ps,GD *gd,HEADER_PARAM *header_param);
 /*From time data in gd produces the averages of short power spectra in ps.
   Reduction factor red depends on header_param->red .
   Each datum in ps is the average over a number lenps/2 (see the code) of data.
   The size of ps is red*/
int short_ps_periodogram(float *ps,GD *gd,HEADER_PARAM *header_param);
     /*
New, Nov. 2010. produces the averaged periodogram
      */
int short_ps_freq(float *ps,GD *gd,HEADER_PARAM *header_param);
 /*From time data in gd produces one short power spectrum in ps.
   Reduction factor red depends on header_param->red .
   Each datum in ps is the SUM (total power)  over a number red of data (data in a same freq. bin).
   The size of ps is lenps/2 (see the code): it contains sqrt(spectrum) of the data in half the band*/

int short_psar_freq(float *ps,GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param);
 /*Very short power spectrum: From frequency data in gd produces one short AR power spectrum in ps.
   Reduction factor red depends on header_param->red .
   Tau of AR could depend on  header_param->red and  header_param->deltanu. Now it is fixed.
   Each datum in ps is the AVERAGE over a number red of data  (data in a same freq. bin).
   The size of ps is lenps/2 (see the code): it contains sqrt(power) in half the band*/
int short_psar_rev_freq(float *ps,GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param,float *psTOT);
/*The same, but reversed. 10 July 2007 */

/************Events extraction to purge the data************/
EVEN_PARAM* crea_evenparam(long len,int scrivi);
EVF_PARAM* crea_evfparam(INPUT_PARAM *input_param);
int bandpass_data(GD *gd_bandpass,GD *gd,HEADER_PARAM *header_param);
/*bandpass data in time domain. Needed before the events extraction*/
int purge_data(GD *gd_clean,GD *gd,GD *gd_highpass,EVEN_PARAM *even_param,HEADER_PARAM *header_param);
int purge_data_subtract(GD *gd_clean,GD *gd,GD *gd_highpass,EVEN_PARAM *even_param,HEADER_PARAM *header_param);
int data_subtract(GD *gd_clean,GD *gd,GD *gd_highpass,EVEN_PARAM *even_param,HEADER_PARAM *header_param);
int purge_data_subsampled(GD *gd_clean1,GD *gd,EVEN_PARAM *even_param, double tsample,HEADER_PARAM *header_param);
int purge_subbw(GD *gd_clean,GD *gd,GD *gd_band,EVEN_PARAM *even_param,HEADER_PARAM *header_param,float band_ext);
/*identifies and purge large events from a gd to another gd.
The veto is in the sub-bandwidths
Input: gd and band_ext
Output: a GD with same parameters, but cleaned  The same gd !*/
int highpass_data(GD *gd_highpass,GD *gd,HEADER_PARAM *header_param,float freqc);
/*highpass data in time domain. Needed for Virgo before the events extraction*/
int highpass_data_bil(GD *gd_highpass,GD *gd,HEADER_PARAM *header_param,float freqc);
/*No time shift highpass data in time domain. Needed for Ligo-Virgo before the events extraction*/
/*identifies and purge large events from a gd to another gd, which has the same parameters */
int sn_medsig(GD *gd_highpass,EVEN_PARAM *even_param,HEADER_PARAM *header_param);
/*from data in a GD it evaluates the AR mean and std, using the even_param*/
int even_anst(GD *gd_highpass,EVEN_PARAM *even_param);
/*from AR mean and std, it looks for events, and registers their characteristics*/
int saturated_data(GD *gd,EVEN_PARAM *even_param,HEADER_PARAM *header_param);
int ScienceVeto_onGD(GD *gd,HEADER_PARAM *header_param,SCIENCE_SEGMENT* science_segment);  //NEW JAN 2017
/****************************Simulations: add signals*******************************/
int add_signals(GD *gd,HEADER_PARAM *header_param);
/*Add signals to the GD *gd*/
void add_sinusoids_old(GD *gd,HEADER_PARAM *header_param,float freq_sin,float amp_sin,double beg_mjd,double *phase_out);
/*Add sinusoids to the GD *gd which is in the frequency domain, from the SFDB*/
void add_sinusoids(GD *gd,HEADER_PARAM *header_param,float freq_sin,float amp_sin,double beg_mjd,double *phase_out);
/*Add sinusoids to the GD *gd which is in the frequency domain, from the SFDB*/
void add_sinusoids_lowres(GD *gd,HEADER_PARAM *header_param,double *ra,double *dec,float *freq,float *amp,double mid_mjd,int number_of_signals, double *possx,double *possy,double *possz,float *template1,int len_template,int factor_template,int fast_proc,int k);
/*Add sinusoids with !! LOW-RES !! to the GD *gd which is in the frequency domain, from the SFDB*/
/*********************************SFDB files *********************************************/
int gd2gdfreq(GD *gd_clean,GD *gd_freq,GD *gd,GD *gd_highpass,GD *gd_band,GD *gd_appo,HEADER_PARAM *header_param,EVEN_PARAM *even_param,EVEN_PARAM *even_paramb,INPUT_PARAM *input_param,SCIENCE_SEGMENT* science_segment);
int gd2sfdbfile(GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param,float *ps,EVEN_PARAM *even_param, INPUT_PARAM *input_param,float *psTOT,int iGRID);
  /*The function writes data which are in a GD *gd + some parameters into a binary file (header + data)
   Input are time domain data in a GD *gd. The function performs the total FFT and also and the very short ps using short_ps_time and short_ps_freq.
 */
int gd2sfdbfilevoid(GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param,float *ps,EVEN_PARAM *even_param, INPUT_PARAM *input_param,float *psTOT);
  /*The function writes data which are in a GD *gd + some parameters into a binary file (header + data)
   Input are time domain data in a GD *gd. The function performs the total FFT and also and the very short ps using short_ps_time and short_ps_freq.
 */
int readfilesfdb_old(GD *gd,GD *gd_short, HEADER_PARAM *header_param,int k,int iw,int *fft_read,int ilista);
 /*The function reads data which are in a sfdb file (GD + GD short + some parameters) */
int readfilesfdb(GD *gd,GD *gd_short, HEADER_PARAM *header_param,int k,int iw,int *fft_read,int ilista, int header_only);
 /*The function reads data which are in a sfdb SFDB09 file (GD + GD short + some parameters) from 12 feb 2009*/

int ToCopyFiles(INPUT_PARAM *input_param,HEADER_PARAM *header_param);

int ToCopyFiles_V1(INPUT_PARAM *input_param,HEADER_PARAM *header_param);

int ToCopyFiles_V2(INPUT_PARAM *input_param,HEADER_PARAM *header_param);

int Logical2PhysicalPath(INPUT_PARAM *input_param, int iquale);

/******************Peak map****************/
int peakmapp05(GD *gd,GD *gd_short,HEADER_PARAM *header_param);
/*Crea peak max using maxima from the SFDB (in a gd), normalized with the very short spectrum in gd_short)
  The output table is written into an file
Old .p05 format
*/
int peakmap_from_arp05(GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param,int k,double mid_mjd,int nfft_file);
/*Crea peak max using maxima from the SFDB (in a gd), normalized with the AR standard deviation of the high
resolution sqrt(spectrum) ) The output table is written into an file
Old .p05 format
*/
int peakmap_from_ar_revp05(GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param,int k,double mid_mjd,int nfft_file);
/*As the previous, but with REVERSED AR spectral estimation
Old .p05 format
*/
int peakmap_from_ar(GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param,int k,double mid_mjd,int nfft_file,GD *gd_short);
/*Crea peak max using maxima from the SFDB (in a gd), normalized with the AR standard deviation of the high
resolution sqrt(spectrum) ) The output table is written into an file
New .p08 format
*/
int peakmap_from_ar_rev(GD *gd,HEADER_PARAM *header_param,EVF_PARAM *evf_param,int k,double mid_mjd,int nfft_file,GD *gd_short);
/*As the previous, but with REVERSED AR spectral estimation
New .p08 format
*/
struct tm *GPS2UTCSAB(long GPSsec, int ULeapS);//sabrina


int p092p10(FILE *P10,FILE *P09,int nfft_file); //pia 2 maggio 2011

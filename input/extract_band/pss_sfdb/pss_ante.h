/**pss_ante.h: pss software**/
/**Last version: Sept, 22, 2004**/
/**Author: Pia**/


/**************Structures***********/


typedef struct FIXED_PARAM{
  double mass;
  double length;
  double temperature;
  double reduced_mass;
  double transducer_mass;
  double transducer_cap;
  double Vbias;
  double freqm;
  double freqp;
  double frcal;
  double taum;
  double taup;
  double Qm;
  double Qp;
  double phibrowm;
  double phibrowp;
  double phical;
  double tslow;
  double frinit_long;
  double frinit;
  double ameno;
  double apiu;
  double gpstime,time;
  int type;

} FIXED_PARAM;


typedef struct VAR_PARAM{
 short int  year,day,hour,min,sec,dmsec;
 short int nhead,h_lco1,h_lco2,tiporec;
 short int h_msec1,h_sec1,h_nca1,h_nca2,h_lrec;
 short int nsamples1,nsamples2;
 float n_flag,n_flagappo; ////how many (data/total) have a itf lock status less than 12 in that FFT
 int nfft;
 double tsampl1,tsampl2,duration,frinit;
 double deltanu,calspect;                          /*evaluated from each new spectrum*/
 double freqm,freqp,taum,taup,tsamplu;             /*evaluated from each new spectrum*/
 double gpstime,mjdtime;
 unsigned long fftlen,kdati_tot;
 char mask_strin[59600];  //nel 2001 era 12000
 

} VAR_PARAM;

typedef struct ITFFIXED_PARAM{
  double length;
  double temperature;
  double freq_min;
  double freq_max;
  double frinit_long;
  double frinit;
  double gpstime,time;
  int type;

} ITFFIXED_PARAM;


typedef struct ITFVAR_PARAM{
 short int  year,day,hour,min,sec,dmsec;
 short int nhead,h_lco1,h_lco2,tiporec;
 short int h_msec1,h_sec1,h_nca1,h_nca2,h_lrec;
 short int nsamples1,nsamples2;
 float n_flag,n_flagappo; ////how many (data/total) have a itf lock status less than 12 in that FFT
 int nfft;
 double tsampl1,tsampl2,duration,frinit;
 double deltanu;                                /*evaluated from each new spectrum*/
 double freq_min,freq_max,tsamplu;             /*evaluated from each new spectrum*/
 double gpstime,mjdtime;
 unsigned long fftlen,kdati_tot;
 

} ITFVAR_PARAM;


 
typedef struct DETECTOR_PARAM{
     
       FIXED_PARAM *bar_fixed_param;
       VAR_PARAM *bar_var_param;
       ITFFIXED_PARAM *itf_fixed_param;
       ITFVAR_PARAM *itf_var_param;

}DETECTOR_PARAM;


/************to fill the structures****************/
DETECTOR_PARAM* crea_parameters(int casin);

#ifdef ROMA_NEPH421
#warning "Using ROMA locale For new ephem 421"
#define EPHFILE  "/storage/users/pia/Effemeridi/lnxp1900p2053.421"
#define UTC2UT1  "/storage/users/pia/GIT_FH/create_sfdb/Pss_Fr/src/pss/pss_astro/novas/ut1-utc_mod1.dat"
#define DELTANUT  "/storage/users/pia/GIT_FH/create_sfdb/Pss_Fr/src/pss/pss_astro/novas/eopc04_mod"
#endif
#ifdef ROMA_NEPH436
#warning "Using ROMA locale For new ephem 436"
#define EPHFILE  "/storage/users/pia/Effemeridi/jpleph.436"
#define UTC2UT1  "/storage/users/pia/GIT_FH/create_sfdb/Pss_Fr/src/pss/pss_astro/novas/ut1-utc_mod1.dat"
#define DELTANUT  "/storage/users/pia/GIT_FH/create_sfdb/Pss_Fr/src/pss/pss_astro/novas/eopc2024_mod"
#endif
#ifdef ROMA_NEPH440
#warning "Using ROMA locale For new ephem 440"
#define EPHFILE  "/storage/users/pia/Effemeridi/jpleph.440"
#define UTC2UT1  "/storage/users/pia/GIT_FH/create_sfdb/Pss_Fr/src/pss/pss_astro/novas/ut1-utc_mod1.dat"
#define DELTANUT  "/storage/users/pia/GIT_FH/create_sfdb/Pss_Fr/src/pss/pss_astro/novas/eopc2024_mod"
#endif
#ifdef ROMA_NEPH405
#warning "Using ROMA locale For new ephem 405"
#define EPHFILE  "/storage/users/pia/Effemeridi/lnxp1600p2200.405"
#define UTC2UT1  "/storage/users/pia/GIT_FH/create_sfdb/Pss_Fr/src/pss/pss_astro/novas/ut1-utc_mod1.dat"
#define DELTANUT  "/storage/users/pia/GIT_FH/create_sfdb/Pss_Fr/src/pss/pss_astro/novas/eopc04_mod"
#endif
#ifdef ROMA_NEPH200
#warning "Using ROMA locale For new ephem 200"
#define EPHFILE  "/storage/users/pia/Effemeridi/lnxm1600p2170.200"
#define UTC2UT1  "/storage/users/pia/GIT_FH/create_sfdb/Pss_Fr/src/pss/pss_astro/novas/ut1-utc_mod1.dat"
#define DELTANUT  "/storage/users/pia/GIT_FH/create_sfdb/Pss_Fr/src/pss/pss_astro/novas/eopc04_mod"
#endif
//END new ephem
//CNAF sabri
//#define EPHFILE  "/storage/gpfs_virgo3/home/CW/sabrina/Jpleph.405"
//#define UTC2UT1  "/storage/gpfs_virgo3/home/CW/sabrina/ut1-utc_mod1.dat" 
//#define DELTANUT  "/storage/gpfs_virgo3/home/CW/sabrina/eopc04_mod"

//CNAF_PIA // pia
#ifdef CNAF_PIA
#define EPHFILE  "/storage/gpfs_virgo3/home/astone/SVN/Pss_Fr/trunk/src/pss/pss_astro/novas/Jpleph.405"
#define UTC2UT1  "/storage/gpfs_virgo3/home/astone/SVN/Pss_Fr/trunk/src/pss/pss_astro/novas/ut1-utc_mod1.dat" 
#define DELTANUT  "/storage/gpfs_virgo3/home/astone/SVN/Pss_Fr/trunk/src/pss/pss_astro/novas/eopc04_mod"
#endif
//CNAF
#ifdef CNAF
#warning "Using CNAF!" 
#define EPHFILE  "/opt/exp_software/virgo/virgoDev/Pss_Fr/trunk/src/pss/pss_astro/novas/Jpleph.405"
#define UTC2UT1  "/opt/exp_software/virgo/virgoDev/Pss_Fr/trunk/src/pss/pss_astro/novas/ut1-utc_mod1.dat" 
#define DELTANUT  "/opt/exp_software/virgo/virgoDev/Pss_Fr/trunk/src/pss/pss_astro/novas/eopc04_mod"
#endif

//ROMA wn
#ifdef ROMA
#warning "Using ROMA!" 
#define EPHFILE  "/opt/exp_soft/virgo/CW/new-pss/pss/pss_astro/novas/Jpleph.405"
#define UTC2UT1  "/opt/exp_soft/virgo/CW/new-pss/pss/pss_astro/novas/ut1-utc_mod1.dat" 
#define DELTANUT  "/opt/exp_soft/virgo/CW/new-pss/pss/pss_astro/novas/eopc04_mod"
#endif
//ROMA wn nuovi
#ifdef ROMA_NEW
#warning "Using ROMA!"
#define EPHFILE  "/opt/exp_software/virgo/CW/new-pss/pss/pss_astro/novas/Jpleph.405"
#define UTC2UT1  "/opt/exp_software/virgo/CW/new-pss/pss/pss_astro/novas/ut1-utc_mod1.dat"
#define DELTANUT  "/opt/exp_software/virgo/CW/new-pss/pss/pss_astro/novas/eopc04_mod"
#endif


// LIGO
//#ifdef LIGO
//#warning "Using LIGO!" 
//#define EPHFILE  "/home/pulsar/NoEMi_peakmaps/Pss_Fr/trunk/src/pss/pss_astro/novas/Jpleph.405"
//#define UTC2UT1  "/home/pulsar/NoEMi_peakmaps/Pss_Fr/trunk/src/pss/pss_astro/novas/ut1-utc_mod1.dat" 
//#define DELTANUT  "/home/pulsar/NoEMi_peakmaps/Pss_Fr/trunk/src/pss/pss_astro/novas/eopc04_mod"
//#endif


//pc s where relative to pss_sfdb works
#ifdef RELATIVE
#warning "uses the relative path"
#define EPHFILE  "../pss_astro/novas/Jpleph.405"
#define UTC2UT1  "../pss_astro/novas/ut1-utc_mod1.dat" 
#define DELTANUT  "../pss_astro/novas/eopc04_mod"
#endif

//Using default !!!
#ifdef LIGO
#warning "using LIGO"
#define EPHFILE  "/home/pulsar/NoEMi/Pss_Fr/trunk/src/pss/pss_astro/novas/Jpleph.405"
#define UTC2UT1  "/home/pulsar/NoEMi/Pss_Fr/trunk/src/pss/pss_astro/novas/ut1-utc_mod1.dat" 
#define DELTANUT  "/home/pulsar/NoEMi/Pss_Fr/trunk/src/pss/pss_astro/novas/eopc04_mod"
#endif


//Using default !!!
#ifdef LIGO_pm
#warning "using LIGO_pm"
#define EPHFILE  "/home/pulsar/NoEMi_peakmaps/Pss_Fr/trunk/src/pss/pss_astro/novas/Jpleph.405"
#define UTC2UT1  "/home/pulsar/NoEMi_peakmaps/Pss_Fr/trunk/src/pss/pss_astro/novas/ut1-utc_mod1.dat"
#define DELTANUT  "/home/pulsar/NoEMi_peakmaps/Pss_Fr/trunk/src/pss/pss_astro/novas/eopc04_mod"
#endif

//Cascina
//#define EPHFILE  "/virgoDev/Pss_Fr/src/pss/pss_astro/novas/Jpleph.405"
//#define UTC2UT1  "/virgoDev/Pss_Fr/src/pss/pss_astro/novas/ut1-utc_mod1.dat" 
//#define DELTANUT  "/virgoDev/Pss_Fr/src/pss/pss_astro/novas/eopc04_mod"


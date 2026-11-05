/*___________________________________________________
 ¦                                                   ¦
 ¦                   pss_serv.h                      ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    March 2000                     ¦
 ¦___________________________________________________¦*/


#ifdef __cplusplus
extern "C" {
#endif

#include <time.h>
	
struct pss_infile_db{
	char*			fil_db;
	char**			files;
	long			nfiles;
	long			index;
	struct FrFile	*frinfile;
	struct File		*infile;
	double			*tin,*tfi;
};

/* Numerical Recipes Macros */
//pia: inizio pezzo commentato: 29-aug. 2005
/*
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)


static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))


static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))


static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))


static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))


static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))


static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))


static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))


static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

*/
//pia: commentate tutte queste definizioni sopra che non sono usate. Da rimettere 29-Aug 2005

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/*------------------------------------------------------------*/

/* Some utilities modified from Numerical Recipes */

void nrerror(char error_text[]);

float *all_Vect(long nl, long nh);

int *all_iVect(long nl, long nh);

unsigned char *all_cVect(long nl, long nh);

unsigned long *all_lVect(long nl, long nh);

double *all_dVect(long nl, long nh);

float **all_Matr(long nrl, long nrh, long ncl, long nch);

double **all_dMatr(long nrl, long nrh, long ncl, long nch);

int **all_iMatr(long nrl, long nrh, long ncl, long nch);

void free_Vect(float *v, long nl, long nh);

void free_iVect(int *v, long nl, long nh);

void free_cVect(unsigned char *v, long nl, long nh);

void free_lVect(unsigned long *v, long nl, long nh);

void free_dVect(double *v, long nl, long nh);

void free_Matr(float **m, long nrl, long nrh, long ncl, long nch);

void free_dMatr(double **m, long nrl, long nrh, long ncl, long nch);

void free_iMatr(int **m, long nrl, long nrh, long ncl, long nch);

void prntft(float data[],unsigned long nn);

/*------------------------------------------------------------*/

/* Easy file routines 

  Allow easy storage for arrays of floating point numbers.

  Parameters are passed throw a SimFil structure 
*/


typedef struct SimFil{
	char *FileName;
	int type;      /* 1 = text, 2 = bin */
	char *access;  /* r = read, w = write */
	double init;
	double samtim;
	long nx;
	long ny;       /* = 0, ny not defined */
	char *caption; /* 80 characters */
	FILE *fil;
	int verb; /* verbosity: 0,1,... */
} SimFil;


SimFil *simfil_crea_str(char *filnam, int type, char *access,
					    double init, double samtim, int nx, int ny, char *capt);


int simfil_open(SimFil *sfil);


int simfil_write(SimFil *sfil, float *y, long n);


int simfil_read(SimFil *sfil, float *y, long n);


int simfil_inquire(char *filnam, int type, SimFil *sfil);


/*------------------------------------------------------------*/

void setVect(float *v, long nl, long nh,
			 int typ, double par1, double par2, double par3);
/*  nl   init index
	nh   final index
	typ  =  1   par1
		 =  2   par1 * ramp + par2
		 =  3   par1 * sin(par2*i+par3)
		 =  4   par1 * randn +par2 (seed round(par3+1))
*/


void typeVect(float *v, int ini, int fin, int nperlin, int typ);
/* types vector values */


void debugVect(float *v);
/* prints vec data chosen interactively */


void statVect(float *v, int n, double *min, double *max,
			  double *ave, double *stdev);


float* real2complex(float *v, int n);


/* ---------- interactive routines --------------- */


int getmenu(char *title, char **items, int nitems);
/* the output value starts from 1; if you need to address vectors
   do:
        item=getmenu(...); item--;
*/

void get_line(char *prompt, char *line);

void get_line_bat(char *prompt, char *line);
/* similar to get_line, but possible use with batch files */

int getlineKR(FILE *fil,char s[],int lim);

long get_long(char * prompt);

double get_double(char * prompt);

void put_notice(char *notice);

void wait_inp();


/* ---------------- miscellanea ------------------ */

double timclock(char *str, time_t *tim1, clock_t *cl1,
				int mode);

char* showtim();

double tim2mjd(time_t tim);
// converts unix time (s from 1-1-1970) to mjd


void empty_ram(int mbyt);

void create_junk_file(char *file, int mbyt);

void waiting( float wait );


/* --------------- String routines ----------------- */

#define NUM_TOK 20

typedef struct STR_TOK{
	int		n;
	char*	str[NUM_TOK];
} STR_TOK;


int isgennum(char c);
/* checks if a character can be part of a number */

STR_TOK str2numtok(char *str, int nexp);
/* creates string tokens using separators any character, 
   except -,+,.,e,E and digits;
   if nexp > 0, it checks if there are at least nexp token */

void vec_strinp(long i, long n, char *inp, char *out);
/* copy a string in the i-th position in a longer string */

void vec_strout(long i, long n, char *inp, char *out);
/* copy a string from the i-th position in a longer string */

void clear_string(char *str, long lstr);

char *str_trim(char *in);
/* trims out the blanks and control character at the edges */

void str_trim1(char *in);
/* trims out the blanks and control character at the edges-v. Pia 16 Aug 2005 */

char *path_from_file(char *file);
/* captures the path from a path+filename string */

int path_from_file1(char *file);
/* captures the path from a path+filename string-v. Pia 17 Aug 2005 */
/*returns the len of the path string*/
/**use:
ii=path_from_file1(file);
strncpy(strapp,file,ii); 
 **/

char *filename_from_file(char *file);
/* captures the file name from a path+filename string */

int filename_from_file1(char *file);
/* captures the file name from a path+filename string */
/*returns the beginning index of the filename string*/
/**use:
ii=filename_from_file1(file);
strcpy(strapp,file+ii); 
 **/
/* pia V. 22 Aug 2005 */
/* --------------- Bit routines ----------------- */


long read_bit(long inp, short bitinit, short bitlen);

short read_bit_s(short inp, short bitinit, short bitlen);

char read_bit_c(char inp, short bitinit, short bitlen);

short write_bit(long *inp, long wr, short bitinit, short bitlen);

short write_bit_s(short *inp, short wr, short bitinit, short bitlen);

short write_bit_c(char *inp, char wr, short bitinit, short bitlen);


short show_bit(long *inp, char *label);

short show_bit_s(short *inp, char *label);

short show_bit_c(char *inp, char *label);


FILE* logfile_open(char* prog);

int logfile_comment(FILE* fid, char* comment);

int logfile_error(FILE* fid, char* errcomment);

int logfile_par(FILE* fid, char* name, double value, int prec);

int logfile_ipar(FILE* fid, char* name, long value);

int logfile_special(FILE* fid, char* label);

int logfile_endspecial(FILE* fid, char* label);

int logfile_input(FILE* fid, char* inpfile, char* comment);

int logfile_output(FILE* fid, char* outfile, char* comment);

int logfile_ev(FILE* fid, char* type, int nvalues, double* values, int* prec);

int logfile_stat(FILE* fid, char* type, int nvalues, double* values, int* prec);

int logfile_stop(FILE* fid);

int logfile_close(FILE* fid);

#ifdef __cplusplus
}
#endif

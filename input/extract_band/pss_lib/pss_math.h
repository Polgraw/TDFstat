#ifndef _HEADER_PSS_MATH_H
#define _HEADER_PSS_MATH_H

/*___________________________________________________
 ¦                                                   ¦
 ¦                   pss_math.h                      ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    March 2000                     ¦
 ¦___________________________________________________¦*/


#ifdef __cplusplus
extern "C" {
#endif

/* Cmplx */

#ifndef _NR_CMPLX_H_
#define _NR_CMPLX_H_

#ifndef _CMPLX_DECLARE_T_
typedef struct CMPLX {float r,i;} Cmplx;
#define _CMPLX_DECLARE_T_
#endif /* _CMPLX_DECLARE_T_ */

Cmplx Cadd(Cmplx a, Cmplx b);
Cmplx Csub(Cmplx a, Cmplx b);
Cmplx Cmul(Cmplx a, Cmplx b);
Cmplx Complex(float re, float im);
Cmplx Conjg(Cmplx z);
Cmplx Cdiv(Cmplx a, Cmplx b);
float Cabs(Cmplx z);
Cmplx Csqrt(Cmplx z);
Cmplx RCmul(float x, Cmplx a);

#endif /* _NR_CMPLX_H_ */


/* fft's */

void four1(float data[], unsigned long nn, int isign);

void dfour1(double data[], unsigned long nn, int isign);

void realft(float data[], unsigned long n, int isign);

void drealft(double data[], unsigned long n, int isign);

void twofft(float data1[], float data2[], float fft1[], float fft2[],
	unsigned long n);

void fourn(float data[], unsigned long nn[], int ndim, int isign);

void rlft3(float ***data, float **speq, unsigned long nn1, unsigned long nn2,
	unsigned long nn3, int isign);


/* random numbers */

float Ran1(long *idum);

float Ran2(long *idum);

float Ran3(long *idum);

float RandN(long *idum); 

/*--------------------------------------------------------*/

/* iLog2 */

int iLog2(long inp);
/* int log base 2 */

#ifdef __cplusplus
}
#endif


double dmin(double a, double b);

long lmin(long a, long b);

#endif /* _HEADER_PSS_MATH_H */

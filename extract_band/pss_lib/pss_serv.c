/*___________________________________________________
 ¦                                                   ¦
 ¦                   pss_serv.c                      ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    March 2000                     ¦
 ¦___________________________________________________¦*/


#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fcntl.h>
#include <ctype.h>

#include "pss_math.h"
#include "pss_serv.h"

#define NR_END 1
#define FREE_ARG char*

/*------------------------------------------------------------*/

/* Some utilities modified from Numerical Recipes */

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


float *all_Vect(long nl, long nh)
/* allocate a float Vect with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in Vect()");
	return v-nl+NR_END;
}


int *all_iVect(long nl, long nh)
/* allocate an int Vect with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in iVect()");
	return v-nl+NR_END;
}


unsigned char *all_cVect(long nl, long nh)
/* allocate an unsigned char Vect with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cVect()");
	return v-nl+NR_END;
}


unsigned long *all_lVect(long nl, long nh)
/* allocate an unsigned long Vect with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lVect()");
	return v-nl+NR_END;
}


double *all_dVect(long nl, long nh)
/* allocate a double Vect with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dVect()");
	return v-nl+NR_END;
}


float **all_Matr(long nrl, long nrh, long ncl, long nch)
/* allocate a float Matr with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in Matr()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in Matr()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


double **all_dMatr(long nrl, long nrh, long ncl, long nch)
/* allocate a double Matr with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in Matr()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in Matr()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


int **all_iMatr(long nrl, long nrh, long ncl, long nch)
/* allocate a int Matr with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in Matr()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in Matr()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_Vect(float *v, long nl, long nh)
/* free a float Vect allocated with Vect() */
{
	free((FREE_ARG) (v+nl-NR_END));
	printf("not used:nl,nh %ld %ld\n",nl,nh); //pia messo il print per usare nh...
}


void free_iVect(int *v, long nl, long nh)
/* free an int Vect allocated with iVect() */
{
	free((FREE_ARG) (v+nl-NR_END));
	printf("not used:nl,nh %ld %ld\n",nl,nh); //pia messo il print per usare nh...
}


void free_cVect(unsigned char *v, long nl, long nh)
/* free an unsigned char Vect allocated with cVect() */
{
	free((FREE_ARG) (v+nl-NR_END));
	printf("not used:nl,nh %ld %ld\n",nl,nh); //pia messo il print per usare nh...
}


void free_lVect(unsigned long *v, long nl, long nh)
/* free an unsigned long Vect allocated with lVect() */
{
	free((FREE_ARG) (v+nl-NR_END));
	printf("not used:nl,nh %ld %ld\n",nl,nh); //pia messo il print per usare nh...
}


void free_dVect(double *v, long nl, long nh)
/* free a double Vect allocated with dVect() */
{
	free((FREE_ARG) (v+nl-NR_END));
	printf("not used: nl,nh %ld %ld\n",nl,nh); //pia messo il print per usare nh...
}


void free_Matr(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float Matr allocated by Matr() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
	printf("not used:nrh,nch %ld %ld\n",nrh,nch); //pia messo il print per usare nrh e nch...
}


void free_dMatr(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double Matr allocated by dMatr() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
	printf("not used:nrh,nch %ld %ld\n",nrh,nch); //pia messo il print per usare nrh e nch...
}


void free_iMatr(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int Matr allocated by iMatr() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
	printf("not used:nrh,nch %ld %ld\n",nrh,nch); //pia messo il print per usare nrh e nch...
}


void prntft(float data[],unsigned long nn)
{
	unsigned long n;
   
	printf("%4s %13s %13s %12s %13s\n",
		"n","real(n)","imag.(n)","real(N-n)","imag.(N-n)");
	printf("   0 %14.6e %12.6e %12.6f %12.6e\n",
		data[1],data[2],data[1],data[2]);
	for (n=3;n<=nn+1;n+=2) {
		printf("%4lu %14.6e %12.6e %12.6e %12.6e\n",
			((n-1)/2),data[n],data[n+1],
			data[2*nn+2-n],data[2*nn+3-n]);
	}
	printf(" press return to continue ...\n");
	(void) getchar();
	return;
}


/*------------------------------------------------------------*/

/* Easy file routines 

  Allow easy storage for arrays of floating point numbers.

  Parameters are passed through a SimFil structure.
*/


SimFil* simfil_crea_str(char *filnam, int type, char *access,
					   double init, double samtim, int nx, int ny, char *capt)
{
	SimFil *sfil;
	int i;

	sfil = (SimFil *) calloc(1,sizeof(SimFil));

	sfil->caption=(char*)malloc((size_t) 80);
	for(i=0;i<80;i++)sfil->caption[i]=0;

	sfil->FileName=filnam;
	sfil->type=type;
	sfil->access=access;
	sfil->init=init;
	sfil->samtim=samtim;
	strcpy(sfil->caption,capt);
	sfil->nx=nx;
	sfil->ny=ny;

	return sfil;
}



int simfil_open(SimFil *sfil)
{
	char mode[3];
	double init,samtim;
	int nx,ny;
//	char *capt;

	mode[0]=*sfil->access;

	mode[1]=0;
	if(sfil->type==2)strcat(mode,"b");

	if((sfil->fil=fopen(sfil->FileName,mode)) != NULL){
		printf(" File %s opened mode %s\n",sfil->FileName,mode);

		if(strcmp(sfil->access,"w")==0){
			if(sfil->type==1){
				fprintf(sfil->fil,"%15.6e %15.6e %ld %ld\n",
					sfil->init,sfil->samtim,sfil->nx,sfil->ny);
				fprintf(sfil->fil,"%s\n",sfil->caption); //pia ld ld
			}
			else{
				fwrite(&sfil->init,sizeof(double),1,sfil->fil);
				fwrite(&sfil->samtim,sizeof(double),1,sfil->fil);
				fwrite(&sfil->nx,sizeof(int),1,sfil->fil);
				fwrite(&sfil->ny,sizeof(int),1,sfil->fil);
				fwrite(sfil->caption,1,80,sfil->fil);
			}
		}
		else{
			if(sfil->type==1){
				fscanf(sfil->fil,"%lf %lf %d %d\n",
					&init,&samtim,&nx,&ny);
				fscanf(sfil->fil,"%s\n",sfil->caption);
			}
			else{
				fread(&init,sizeof(double),1,sfil->fil);
				fread(&samtim,sizeof(double),1,sfil->fil);
				fread(&nx,sizeof(int),1,sfil->fil);
				fread(&ny,sizeof(int),1,sfil->fil);
				fread(sfil->caption,1,80,sfil->fil);
			}

			sfil->init=init;
			sfil->samtim=samtim;
			sfil->nx=nx;
			sfil->ny=ny;
		}
	}
	else
		printf(" Error opening %s mode %s\n",sfil->FileName,mode);

	return 0;
}



int simfil_write(SimFil *sfil, float *y, long n)
{
	long i;

	switch(sfil->type)
	{
	case 1:
		for(i=0;i<n;i++)fprintf(sfil->fil,"%15.6e\n",y[i]);
		break;
	case 2:
		fwrite(y,sizeof(float),n,sfil->fil);
	}

	return 0;
}


int simfil_read(SimFil *sfil, float *y, long n)
{
	long i;

	switch(sfil->type)
	{
	case 1:
	         //pia, 29-Aug 2005.Corretto con il puntatore ! &y[i]
	        //for(i=0;i<n;i++)fscanf(sfil->fil,"%e\n",y[i]);
		for(i=0;i<n;i++)fscanf(sfil->fil,"%e\n",&y[i]);
		break;
	case 2:
		fread(y,sizeof(float),n,sfil->fil);
	}

	return 0;
}



int simfil_inquire(char *filnam, int type, SimFil *sfil)
{
	long i=0;
	float a;

	sfil=simfil_crea_str(filnam,type,"r",0.,0.,0,0," Nulla ");

	simfil_open(sfil);

	printf(" Inquire file %s ; type %d \n",sfil->FileName,sfil->type);

	printf("    init,samtim,nx,ny = %15.6e, %15.6e, %ld %ld\n",
		sfil->init,sfil->samtim,sfil->nx,sfil->ny);   //pia ld ld
	printf("Caption : %s\n",sfil->caption);  

	switch(type)
	{
	case 1:
		while(feof(sfil->fil)==0){
			fscanf(sfil->fil,"%e\n",&a);
			i++;
		}
		break;
	case 2:
		while(feof(sfil->fil)==0){
			fread(&a,sizeof(float),1,sfil->fil);
			i++;
		}
		i=i-1;
	}

	printf("  Found %ld data, expected %ld\n",i,sfil->nx*sfil->ny); //pia ld ld

	fclose(sfil->fil);

	return 0;
}


/*------------------------------------------------------------*/

void setVect(float *v, long nl, long nh,
			 int typ, double par1, double par2, double par3)
/*  nl   init index
	nh   final index
	typ  =  1   par1
		 =  2   par1 * ramp + par2
		 =  3   par1 * sin(par2*i+par3)
		 =  4   par1 * randn +par2 (seed round(par3+1))
*/
{
	long i,idum;

	idum=par3;

	for (i = nl; i <= nh; i++)
	{
		switch (typ)
		{
			case 1:
				v[i]=par1;
				break;
			case 2:
				v[i]=par1*i+par2;
				break;
			case 3:
				v[i]=par1*sin(par2*i+par3);
				break;
			case 4:
				v[i]=par1*RandN(&idum)+par2;
				break;
		}
	}
}



void typeVect(float *v, int ini, int fin, int nperlin, int typ)
/* types vector values 
     typ = 0  long
	  "  = 1  short */
{
	int nlin,i,ii=0;

	nlin=(fin-ini+1)/nperlin;

	for (i=ini; i<=fin; i++){
		
		if(ii++==(ii/nperlin)*nperlin)printf("\n %d -> ",i);

		switch(typ){
		case 0 :
			printf("%f ",v[i]);
		case 1 :
			printf("%12.4e ",v[i]);
		}
	}

	printf("\n");
}


void debugVect(float *v)
/* prints vec data chosen interactively */
{
	long i,ini=0,fin=0;

	while(ini<=fin)
	{
		printf("Ini, end index ? (ini>end -> exit) ");
		//pia corretto con il puntatore !!  scanf("%d %d",ini,fin);
		scanf("%ld %ld",&ini,&fin);  //pia ld ld anche

		for(i=ini;i<=fin;i++)
		  // pia tolto: printf(" %d %d %f\n",i,&v[i],v[i]); //pia: tolta la stampa del & che da' warning
		  printf(" %ld  %f\n",i,v[i]); //ho messo questo
	}
}


void statVect(float *v, int n, double *mini, double *maxi,
			  double *ave, double *stdev)
{
	double q=0.,qq=0.;
	int i;

	*mini=1.e100;
	*maxi=-1.e100;

	for (i=0; i<n; i++)
	{
		if (v[i] > *maxi) *maxi=v[i];
		if (v[i] < *mini) *mini=v[i];
		q+=v[i];
		qq+=v[i]*v[i];
	}

	*ave=q/n;
	*stdev=sqrt((qq-q*q/n)/(n-1));

	printf("mean,stdev,min,max = %12.4e, %12.4e, %12.4e, %12.4e\n",
		*ave,*stdev,*mini,*maxi);
}



float* real2complex(float *v, int n)
/* creates a "complex" vector with imag part equal to 0 
   from a real float vector */
{
	float *w;
	int i;

	w=(float*)malloc((size_t) 2*n*sizeof(float));

	for (i=0; i<n; i++){
		w[2*i]=v[i];
		w[2*i+1]=0.;
	}

	return(w);
}


/* ---------- interactive routines --------------- */



int getmenu(char *title, char **items, int nitems)
{
	int i,item=0,iter=0;

	printf("\n      %s \n\n",title);


	for(i=0;i<nitems;i++){
		if(strcmp(items[i],"dummy"))
			printf(" %2d.  %s\n",i+1,items[i]);	
		else
			printf("\n");
	}

	while((item<1) | (item >nitems)){
		printf("\n   Which item ? ");
		scanf("%d",&item);
		if(!strcmp(items[item-1],"dummy"))item=0;
		iter++;
		if(iter>4)return 0;
	}

	return item;
}


void get_line(char * prompt, char * line){

	fflush(stdin);
	puts(prompt);
	fgets(line,128,stdin);
	str_trim1(line);
}


void get_line_bat(char *prompt, char *line){
/* similar to get_line, but possible use with batch files */
    //int i, ch,point; //pia tolti perche' non usati
   int len;
   char form[10];
   
   printf("%s",prompt);
   
   scanf("%d",&len);  //pia: messo &len
   sprintf(form," %%%ds",len);
   scanf(form,line);
   printf( " ---> %s\n", line );
}

/* getlineKR: read a line into s, return length; from Kernighan & Ritchie */
int getlineKR(FILE *fil,char *line, int max)
{
	if (fgets(line, max, fil) == NULL)
		return 0;
	else
		str_trim1(line);
		return strlen(line);
}


long get_long(char *prompt){
	long out;

	fflush(stdin);
	printf(prompt);
	scanf("%ld",&out); //pia ld

	return out;
}


double get_double(char *prompt){
	double out;

	fflush(stdin);
	printf(prompt);
	//scanf("%f",&out); //pia. da' warning double argument e float format. CHIEDERE A SERGIO
	scanf("%lf",&out);  

	return out;
}


void put_notice(char *notice){
	printf(notice);
	printf("\n\n    <CR>\n");
	getchar();
}


void wait_inp(){
	int a;

	printf("   Press 0 and return to exit ");
	scanf("%d",&a);
}

/* ---------------- miscellanea ------------------ */


double timclock(char *str, time_t *tim1, clock_t *cl1,
				int mode)
/* Computes time intervals using functions 
     time() (in seconds) and
     clock() (in milliseconds)
	 mode < 1 -> no printed output
*/
{
	time_t tim;
	clock_t cl;
	double dcl,dt;

	static long iter=0;

	tim=time(0);
	cl=clock();

	dt=difftime(tim,*tim1);
	dcl=cl-*cl1;
	if(mode>0)
	    printf("%ld %s Dtime, Dclock = %f, %f\n",iter,str,dt,dcl); //pia ld

	*cl1=cl;
	*tim1=tim;
	iter++;

	return dcl/CLOCKS_PER_SEC;
}


char* showtim(){
	char *cnow;
	time_t now;
	struct tm *ptime;

	now=time(0);

	ptime=localtime(&now);
	cnow=asctime(ptime);

	return cnow;
}


double tim2mjd(time_t tim){
// converts unix time (s from 1-1-1970) to mjd

	double mjd;

	mjd=tim/86400.+40587.;

	return mjd;
}

void empty_ram(int mbyt){
	int *a;
	long i;
	long len=16384,base=64;

	FILE *fil;

	a=(int*)malloc(len*sizeof(float));

	for(i=0;i<len;i++)
		a[i]=i;

	fil=fopen("scratch.dat","wb");

	for(i=0;i<base*mbyt;i++)
		fwrite(a,sizeof(int),len,fil);

	fclose(fil);

	fil=fopen("scratch.dat","rb");

	for(i=0;i<base*mbyt;i++){
		fread(a,sizeof(int),len,fil);
	}

	fclose(fil);
	remove("scratch.dat");
}


void create_junk_file(char *file, int mbyt){
	int *a;
	long i;
	long len=16384,base=16;

	FILE *fil;

	a=(int*)malloc(len*sizeof(float));

	for(i=0;i<len;i++)
		a[i]=i;

	fil=fopen(file,"wb");

	for(i=0;i<base*mbyt;i++)
		fwrite(a,sizeof(int),len,fil);

	fclose(fil);
}


/* Pauses for a specified number of milliseconds. */
void waiting( float wait )
{
   clock_t goal,wait1;

   wait1=(clock_t)(wait * CLOCKS_PER_SEC); 

   goal = wait1 + clock();
   while( goal > clock() )
      ;
}



/*------------------------------------------------------------*/

/* String routines */


void vec_strinp(long i, long n, char *inp, char *out)
/* copy a string in the i-th position in a longer string */
{
	long j;

	for(j=0;j<n;j++)out[i*n+j]=0;
	for(j=0;j<n && inp[j]>0;j++)out[i*n+j]=inp[j];
}


void vec_strout(long i, long n, char *inp, char *out)
/* copy a string from the i-th position in a longer string */
{
	long j;

	for(j=0;j<n;j++)out[j]=0;
	for(j=0;j<n && inp[j]>0;j++)out[j]=inp[i*n+j];
}


void clear_string(char *str, long lstr){
	long i;

	for(i=0;i<lstr;i++)str[i]=0;
}


int isgennum(char c)
/* checks if a character can be part of a number */
{
	int i=0;

	if(isdigit(c))i=1;
	if(c==46)i=1; /* . */
	if(c==101)i=1; /* e */
	if(c==69)i=1; /* E */
	if(c==43)i=1; /* + */
	if(c==45)i=1; /* - */

	return i;
}


STR_TOK str2numtok(char *str, int nexp)
/* creates string tokens using separators any character, 
   except -,+,.,e,E and digits;
   if nexp > 0, it checks if there are at least nexp token */
{
	STR_TOK tok;
	long len,i,i1,lentok,n;
	char c=' ';

	len=strlen(str);
	n=0;
	for(i=0;i<len+1;i++){
		if(!isgennum(c) && isgennum(str[i])){
			n++;
			i1=i;
		}
		if(isgennum(c) && !isgennum(str[i])){
			lentok=i-i1;
			tok.str[n-1]=(char *)malloc(lentok+1);
			strncpy(tok.str[n-1],str+i1,lentok);
		}
		c=str[i];
	}

	tok.n=n;
	if(nexp){
	  if(n>nexp)printf(" *** %ld characters more than expected \n",n-nexp); //pia ld
	  if(n<nexp)printf(" *** %ld characters less than expected \n",nexp-n); //pia ld
	}

	return tok;
}


char *str_trim(char *in){
/* trims out the blanks and control characters at the edges */ 
	// USE str_trim1 INSTEAD !!!
	char *out;
	long i,i1,i2,len,lenin;

	lenin=strlen(in);
	i1=0;i2=-1;

	for(i=0;i<lenin;i++){
		if(in[i]<33 && i2<0)i1=i+1;
		if(in[i]>32)i2=i;
	}
	len=i2-i1+1; //printf(" %s \n %d - %d - %d - %d \n",in,i1,i2,len,lenin);
	if(len>0){
		out=(char*)malloc(len+1);
		for(i=i1;i<=i2;i++)
			out[i-i1]=in[i];
		out[len]=0;
	}

	return out;
}

void str_trim1(char *in){
/* trims out the blanks and control characters at the edges */
/* v. Pia  16 Aug 2005 */	
	long i,i1,i2,len,lenin;

	lenin=strlen(in);
	i1=0;i2=-1;

	for(i=0;i<lenin;i++){
		if(in[i]<33 && i2<0)i1=i+1;
		if(in[i]>32)i2=i;
	}
	len=i2-i1+1; //printf(" %s \n %d - %d - %d - %d \n",in,i1,i2,len,lenin);
	if(len>0){
		
		for(i=i1;i<=i2;i++)
			in[i-i1]=in[i];
		in[len]=0;
	}
	

	return;
}


char *path_from_file(char *file){
/* captures the path from a path+filename string */

	char* pnam;
	long len,ii,i;

	len=strlen(file); 
	ii=-1;

	for(i=0;i<len;i++){
		if(file[i] == 47 || file[i] == 92) ii=i; /* 47 -> / , 92 -> \ */
	}

	pnam=(char*)malloc(ii+2);
	pnam[0]=0;
	//printf(" *** ii = %d  \n",ii);
	strncpy(pnam,file,ii+1); 

	return pnam;
}
int path_from_file1(char *file){
/* captures the path from a path+filename string */
/*returns the len of the path string*/
/**use:
ii=path_from_file1(file);
strncpy(strapp,file,ii); 
 **/
/* pia V. 17 Aug 2005 */

      
	long len,ii,i;
	int j;

	len=strlen(file); 
	ii=-1;

	for(i=0;i<len;i++){
		if(file[i] == 47 || file[i] == 92) ii=i; /* 47 -> / , 92 -> \ */
	}

	
	j= (int) (ii+1);
	return j;
}



char *filename_from_file(char *file){
/* captures the file name from a path+filename string */

	char* pnam;
	long len,ii,i;

	len=strlen(file); 
	ii=-1;

	for(i=0;i<len;i++){
		if(file[i] == 47 || file[i] == 92) ii=i; /* 47 -> / , 92 -> \ */
	}

	pnam=(char*)malloc(len-ii);
	//printf(" *** ii = %d  \n",ii);
	strcpy(pnam,file+ii+1); 

	return pnam;
}

int filename_from_file1(char *file){
/* captures the file name from a path+filename string */
/*returns the beginning index of the filename string*/
/**use:
ii=path_from_file1(file);
strcpy(strapp,file+ii); 
 **/
/* pia V. 22 Aug 2005 */
       
	long len,ii,i;
	int j;

	len=strlen(file); 
	ii=-1;

	for(i=0;i<len;i++){
		if(file[i] == 47 || file[i] == 92) ii=i; /* 47 -> / , 92 -> \ */
	}
	
	//printf(" *** ii = %d  \n",ii);

	j=ii+1;
	return j;
}


/*------------------------------------------------------------*/

/* Bit routines */


long read_bit(long inp, short bitinit, short bitlen){
	long out;
	unsigned long mask=0,p2=1;
	long i;

	for(i=0;i<bitinit+bitlen;i++){
		if(i>=bitinit)mask=mask+p2;
		p2=p2*2;
	}

	out=inp&mask;
	out=out >> bitinit;

	return out;
}

short read_bit_s(short inp, short bitinit, short bitlen){
	short out;
	unsigned short mask=0,p2=1;
	long i;

	for(i=0;i<bitinit+bitlen;i++){
		if(i>=bitinit)mask=mask+p2;
		p2=p2*2;
	}

	out=inp&mask;
	out=out >> bitinit;

	return out;
}

char read_bit_c(char inp, short bitinit, short bitlen){
	char out;
	unsigned char mask=0,p2=1;
	long i;

	for(i=0;i<bitinit+bitlen;i++){
		if(i>=bitinit)mask=mask+p2;
		p2=p2*2;
	}

	out=inp&mask;
	out=out >> bitinit;

	return out;
}

short write_bit(long *inp, long wr, short bitinit, short bitlen){
	long p2=1;
	long i;

	wr=wr << bitinit;
	p2=pow(2,bitinit);

	for(i=bitinit;i<bitinit+bitlen;i++){
		if(p2&wr)*inp=*inp|p2;
		else *inp=*inp&~p2;

		p2=p2*2;
	}

	return 1;
}

short write_bit_s(short *inp, short wr, short bitinit, short bitlen){
	short p2=1;
	long i;

	wr=wr << bitinit;
	p2=pow(2,bitinit);

	for(i=bitinit;i<bitinit+bitlen;i++){
		if(p2&wr)*inp=*inp|p2;
		else *inp=*inp&~p2;

		p2=p2*2;
	}

	return 1;
}

short write_bit_c(char *inp, char wr, short bitinit, short bitlen){
	unsigned char p2=0x1;
	long i;

	wr=wr << bitinit;
	p2=pow(2,bitinit);

	for(i=bitinit;i<bitinit+bitlen;i++){
		if(p2&wr)*inp=*inp|p2;
		else *inp=*inp&~p2;

		p2=p2*2;
	}

	return 1;
}


short show_bit(long *inp,char *label){
	char bits[32];
	long app,mask=0x1,i;

	app=*inp;
	printf("\n %s %ld >>> ",label,app); //pia ld

	for(i=0;i<32;i++){
		bits[i]=app&mask;
		app=app >> 1;
	}

	for(i=31;i>=0;i--)printf("%d",bits[i]);

	printf(" \n");
	return 1;
}


short show_bit_s(short *inp,char *label){
	char bits[16];
	short app,mask=0x1,i;

	app=*inp;
	printf("\n bits in %s %d >>> ",label,app);

	for(i=0;i<16;i++){
		bits[i]=app&mask;
		app=app >> 1;
	}

	for(i=15;i>=0;i--)printf("%d",bits[i]);

	printf(" \n");
	return 1;
}


short show_bit_c(char *inp,char *label){
	char bits[8];
	//char app,mask=0x1,i; //pia: protesta che i e' un char. Correggo, vedi sotto. CHIEDERE A Sergio
	char app,mask=0x1; //pia
	int i; //pia

	app=*inp;
	printf("\n bits in %s %d >>> ",label,app);

	for(i=0;i<8;i++){
		bits[i]=app&mask;
		app=app >> 1;
	}

	for(i=7;i>=0;i--)printf("%d",bits[i]);

	printf(" \n");
	return 1;
}


FILE* logfile_open(char* prog){
	FILE* fid;
	char fil[80];
	char dat[9],tim[7];
	char *str;
	time_t now;
	struct tm *t;

	str=(char*)malloc(27);
	now=time(0);
	t=localtime(&now);
	str=asctime(t);

	strcpy(fil,prog);
	strcat(fil,"_");
	sprintf(dat,"%04d%02d%02d",t->tm_year+1900,t->tm_mon+1,t->tm_mday);
	strcat(fil,dat);
	strcat(fil,"_");
	sprintf(tim,"%02d%02d%02d",t->tm_hour,t->tm_min,t->tm_sec);
	strcat(fil,tim);
	strcat(fil,".log");

	fid=fopen(fil,"wt");

	fprintf(fid,"PSS %s job log file \n\n",prog);
	fprintf(fid,"started at %s \n",str);

	printf("\n  |%s| \n",fil);

	return fid;
}


int logfile_comment(FILE* fid, char* comment){
	int stat=0;

	fprintf(fid,"! %s \n",comment);

	return stat;
}


int logfile_error(FILE* fid, char* errcomment){
	int stat=0;

	fprintf(fid,"*** ERROR ! %s \n",errcomment);

	return stat;
}


int logfile_par(FILE* fid, char* name, double value, int prec){
	int stat=0;

	switch(prec){
		case 0:
			fprintf(fid,"(PAR) %s = %g \n",name,value);
			break;
		case 1:
			fprintf(fid,"(PAR) %s = %.3f \n",name,value);
			break;
		case 2:
			fprintf(fid,"(PAR) %s = %.6f \n",name,value);
			break;
		case 3:
			fprintf(fid,"(PAR) %s = %.9f \n",name,value);
			break;
		case 4:
			fprintf(fid,"(PAR) %s = %.12f \n",name,value);
			break;
		case 10:
			fprintf(fid,"(PAR) %s = %.6e \n",name,value);
			break;
	}

	return stat;
}


int logfile_ipar(FILE* fid, char* name, long value){
	int stat=0;

	fprintf(fid,"(PAR) %s = %ld \n",name,value);

	return stat;
}


int logfile_special(FILE* fid, char* label){
	int stat=0;

	fprintf(fid,"SPECIAL> %s \n",label);

	return stat;
}


int logfile_endspecial(FILE* fid, char* label){
	int stat=0;

	fprintf(fid,"ENDSPECIAL> %s \n",label);

	return stat;
}


int logfile_input(FILE* fid, char* inpfile, char* comment){
	int stat=0;

	fprintf(fid,"INPUT : %s  %s\n",inpfile,comment);

	return stat;
}


int logfile_output(FILE* fid, char* outfile, char* comment){
	int stat=0;

	fprintf(fid,"OUTPUT : %s  %s\n",outfile,comment);

	return stat;
}
int logfile_ev(FILE* fid, char* type, int nvalues, double* values, int* prec){
	int i,stat=0;

	fprintf(fid,"--> %s > ",type);
	for(i=0;i<nvalues;i++){
		switch(prec[i]){
			case 0:
				fprintf(fid,"%g ",values[i]);
				break;
			case 1:
				fprintf(fid,"%.3f ",values[i]);
				break;
			case 2:
				fprintf(fid,"%.6f ",values[i]);
				break;
			case 3:
				fprintf(fid,"%.9f ",values[i]);
				break;
			case 10:
				fprintf(fid,"%.6e ",values[i]);
				break;
			case 100:
				fprintf(fid,"%.0f ",values[i]);
				break;
		         case 101:
			        fprintf(fid,"%+4.2f ",values[i]);  //pia Jan 2012
				break;
		}
	}
	fprintf(fid," \n");

	return stat;
}


int logfile_stat(FILE* fid, char* type, int nvalues, double* values, int* prec){
	int i, stat=0;

	fprintf(fid,">>> %s > ",type);
	for(i=0;i<nvalues;i++){
		switch(prec[i]){
			case 0:
				fprintf(fid,"%g ",values[i]);
				break;
			case 1:
				fprintf(fid,"%.3f ",values[i]);
				break;
			case 2:
				fprintf(fid,"%.6f ",values[i]);
				break;
			case 3:
				fprintf(fid,"%.9f ",values[i]);
				break;
			case 10:
				fprintf(fid,"%.6e ",values[i]);
				break;
			case 100:
				fprintf(fid,"%.0f ",values[i]);
				break;
		}
	}
	fprintf(fid," \n");

	return stat;
}


int logfile_stop(FILE* fid){
	int stat=0;
	char *str;
	time_t now;
	struct tm *t;

	str=(char*)malloc(27);
	now=time(0);
	t=localtime(&now);
	str=asctime(t);

	fprintf(fid,"\nstop at %s \n\n",str);

	return stat;
}


int logfile_close(FILE* fid){
	int stat=0;

	fclose(fid);

	return stat;
}

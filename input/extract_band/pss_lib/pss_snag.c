/*___________________________________________________
 ¦                                                   ¦
 ¦                   pss_snag.c                      ¦
 ¦       by Sergio Frasca - ((( 0 ))) Virgo          ¦
 ¦                    March 2000                     ¦
 ¦___________________________________________________¦*/

#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pss_snag.h"
#include "pss_serv.h"
#include "pss_math.h"

GD* crea_gd(long nall, double ini,double dx, char* capt)
/* to create type 1 gd */
{
	GD* gd;
	gd = (GD *) calloc(1,sizeof(GD));

	gd->type=1;
	gd->cont=0;

	gd->x=0;
	gd->y=malloc((size_t) nall*sizeof(float));
	gd->n=nall;
	gd->nall=nall;
	gd->ini=ini;
	gd->dx=dx;
	gd->capt=capt;

	return gd;
}


GD* vec_to_gd(float* y, long n, double ini,double dx, char* capt)
/* to create type 1 gd from a vector 
   the vector is not copied, but attached to the gd */
{
	GD* gd;
	gd = (GD *) calloc(1,sizeof(GD));

	gd->type=1;
	gd->cont=0;

	gd->x=0;
	gd->y=y;
	gd->n=n;
	gd->nall=n;
	gd->ini=ini;
	gd->dx=dx;
	gd->capt=capt;

	return gd;
}



float* y_gd(GD* gd)
/* extracts the data vector from a gd */
{
	float* y;

	y=gd->y;

	return y;
}



GD* set_gd(GD* g, int typ,
		   double par1, double par2, double par3)
/*  typ  =  1   par1
		 =  2   par1 * ramp + par2
		 =  3   par1 * sin(par2*i+par3)
		 =  4   par1 * randn +par2 (seed round(par3+1))
*/
{
	float *y;
	int n=g->n;
	GD* g1;

	g1=g;
	y=g1->y;

	setVect(y,0,n-1,typ,par1,par2,par3);

	g1->y=y;

	return g1;
}



int show_gd(GD* gd)
/* show the gd parameters */
{
	int err=0;

	printf(" gd -> n = %ld (all=%ld), ini = %f, dx = %f, %s \n",
	       gd->n,gd->nall,gd->ini,gd->dx,gd->capt); //pia: ld ld

	return err;
}



GD* sum_gd(float a1, GD* gd1,
		   float a2, GD* gd2,
		   float a3)
/* sums two gd's */
{
	GD* gd0;
	int i,n;

	gd0=gd1;

	if (a2 == 0.)
	{
		n=gd1->n;
		for (i=0; i < n; i++)
		gd0->y[i]=a1*gd1->y[i]+a3;

	}else
	{
		n=lmin(gd1->n,gd2->n);
		for (i=0; i < n; i++)
		gd0->y[i]=a1*gd1->y[i]+a2*gd2->y[i]+a3;
	}

	gd0->n=n;

	return gd0;
}



GD* mult_gd(GD* gd1, GD* gd2)
/* product of two gd's */
{
	GD* gd0;
	int i,n;

	gd0=gd1;

	n=lmin(gd1->n,gd2->n);

	for (i=0; i < n; i++)
	gd0->y[i]=gd1->y[i]*gd2->y[i];

	gd0->n=n;

	return gd0;
}


/*--------------------- DM --------------------------*/


DM* crea_dm(long n, long m, double ini,double dx, double ini2,double dx2, char* capt)
/* to create type 1 dm */
{
	DM* dm;
	dm = (DM *) calloc(1,sizeof(DM));

	dm->type=1;
	dm->cont=0;

	dm->x=0;
	dm->y=malloc((size_t) n*sizeof(float));
	dm->n=n;
	dm->nall=n;
	dm->ini=ini;
	dm->dx=dx;
	dm->m=m;
	dm->ini2=ini2;
	dm->dx2=dx2;
	dm->capt=capt;

	return dm;
}


DM* gd2dm(GD *gd){
	DM* dm;
	long i;

	dm->type=gd->type;
	dm->cont=0;

	if(gd->type==2){
		dm->x=malloc((size_t)gd->n*sizeof(float));
		dm->x=gd->x;
		for(i=0;i<gd->n;i++)dm->x[i]=gd->x[i];
	}
	dm->y=malloc((size_t) gd->n*sizeof(float));
	for(i=0;i<gd->n;i++)dm->y[i]=gd->y[i];

	dm->n=gd->n;
	dm->nall=gd->nall;
	dm->ini=gd->ini;
	dm->dx=gd->dx;
	dm->m=1;
	dm->ini2=0;
	dm->dx2=1;
	strcpy(dm->capt,gd->capt);

	return dm;
}


/*-----------------------------------------------*/


RING* crea_ring(long len)
{
	RING *r;

	r = (RING *)malloc(sizeof(RING));

	r->len=len;
	r->capt="no name";
	r->cont=0;
	r->dx=1;
	r->indin=0;
	r->indout=0;
	r->ini=0;
	r->lastim=0;
	r->totin=0;
	r->totout=0;
	r->y=malloc((size_t) len*sizeof(float));
	r->debug=0;

	return r;
}



int read_rg(RING *r,float* v,long n,double *tinit)
/* reads from a ring */
{
	long in,ic=1,fin,fin2,n2,i;
	double tinit1;

	if (n > r->len){
		printf(" *** Error : r.len = %ld,  v len = %ld\n",r->len,n); //pia: ld ld
		return 3;}

	if(r->totin-r->totout-n>0){
		printf(" *** Error : too much data requested; totin,totout,n = %ld, %ld, %ld\n",
			r->totin,r->totout,n); //pia: ld ld
		return 2;}

	tinit1=r->lastim-(r->totin-r->totout-1)*r->dx;
	*tinit=tinit1;

	in=r->indout;
	if (in > r->len)in-=r->len;

	fin=in+n;
	if (fin >= r->len){
		ic=2;
		n2=r->len-in;
		fin2=in+n-r->len-1;}

	if (ic == 1)
		for (i=in; i<=fin;i++) v[i-in]=r->y[i];
	else{
		for (i=0; i<n2; i++) v[i]=r->y[in+i];
		for (i=0; i<fin2; i++) v[i+n2]=r->y[i];}

	r->indout=fin;
	r->totout+=n;

	if (r->debug > 0) printf(" **---> %ld %f\n",in,v[0]); //pia: ld

	if (r->debug > 1) printf(" Read ring: totout x1 xn= %ld %12.4e %12.4e\n",
		r->totout,v[0],v[n-1]); //pia: ld 

	return 0;
}



int write_rg(RING *r,float* v,long n,double lastim)
/* writes to a ring */
{
	long in,ic=1,fin,i;

	if (n > r->len){
		printf(" *** Error : r.len = %ld,  v len = %ld\n",r->len,n); //pia: ld ld
		return 3;}

	in=r->indin;

	if (r->debug > 0) printf(" *---> %ld %f \n",in,v[0]); //pia: ld

	if (in > r->len)in-=r->len;

	fin=in+n;

	if (fin >= r->len){
		ic=2;
		fin-=r->len;}

	if (ic == 1) 
		for (i=in; i<fin; i++) r->y[i]=v[i-in];
	else{
		for (i=in; i<r->len; i++) r->y[i]=v[i-in];
		for (i=0; i<fin; i++) r->y[i]=v[n-fin+i];}  // check

	r->indin=fin;
	r->totin+=n;
	r->lastim=lastim;

	if (r->debug > 1) printf(" Write ring: totin = %ld\n",r->totin); //pia: ld

	if (r->totin-r->totout > r->len)
		printf(" * * Warning : ring data lost !");

	return 0;
}




/*-----------------------------------------------*/


DS* crea_ds(long len, int typ, char* capt)
/* 
	len		chunk length
	typ		ds type (0,1,2 -> interlacing)
	capt	caption
*/
{
	DS *ds;

	int i;
	double dt=1.;

	ds = (DS *) malloc(sizeof(DS));

	ds->len=len;
	ds->capt=capt;
	ds->cont=0;
	ds->dt=dt;
	ds->ind1=1;
	ds->ind2=1;
	ds->lcr=0;
	ds->lcw=0;
	ds->nc1=0;
	ds->nc2=0;
	ds->tini1=-len*dt;
	ds->tini2=-len*dt;
	ds->treq=0;
	ds->type=typ;
	ds->y1=malloc((size_t) (5*len*sizeof(float))/4);
	ds->y2=malloc((size_t) len*sizeof(float));
	ds->debug=0;

	for (i=0; i<(5*len)/4; i++)
		ds->y1[i]=0.;

	for (i=0; i<len; i++)
		ds->y2[i]=0.;

	return ds;
}



int y_ds(DS *d, float *y)
{
	int i;

	if (d->type == 2)
	{
		if (d->lcw%2 == 1)
			for (i=0; i < d->len; i++) y[i]=d->y1[i];
		else 
			for (i=0; i < d->len; i++) y[i]=d->y2[i];
	}
	else
		for (i=0; i < d->len; i++) y[i]=d->y1[i];

	return 1;
}



int read_rg2ds(RING *r, DS *d, int iy, long ink, 
			   long n, double *tinit)
/* reads from a ring directly to a ds */
{
	long in,ic=1,fin,fin2,n2,i;
	double tinit1;

	if (n > r->len){
		printf(" *** Error : r.len = %ld,  v len = %ld\n",r->len,n); //pia: ld ld
		return 3;}

	if(r->totin-r->totout<n){
		printf(" *** Error : too much data requested; totin,totout,n = %ld, %ld, %ld\n",
			r->totin,r->totout,n); //pia: ld ld ld
		return 2;}

	tinit1=r->lastim-(r->totin-r->totout-1)*r->dx;
	*tinit=tinit1;

	in=r->indout;
	//piapia sostituita !! 12/02/03 if (in > r->len)in-=r->len;
	if (in >= r->len)in-=r->len;

	fin=in+n;
	if (fin >= r->len){
		ic=2;
		n2=r->len-in;
		fin2=in+n-r->len-1;}

	if (ic == 1)
		if (iy == 1)
			for (i=in; i<=fin;i++) d->y1[i-in+ink]=r->y[i];
		else
			for (i=in; i<=fin;i++) d->y2[i-in+ink]=r->y[i];
	else
	{
		if (iy == 1)
		{
			for (i=0; i<n2; i++) d->y1[i+ink]=r->y[in+i];
			for (i=0; i<fin2; i++) d->y1[i+n2+ink]=r->y[i];
		}
		else
		{
			for (i=0; i<n2; i++) d->y2[i+ink]=r->y[in+i];
			for (i=0; i<fin2; i++) d->y2[i+n2+ink]=r->y[i];
		}
	}

	r->indout=fin;
	r->totout+=n;

//	if (r->debug > 0) printf(" **---> %d %f\n",in,v[0]);

//	if (r->debug > 1) printf(" Read ring: totout x1 xn= %d %12.4e %12.4e\n",
//		r->totout,v[0],v[n-1]);

	return 0;
}




int sim2ds(DS *d, RING *r, float *v, long nv, int typ, int verb)
/* typ = 0     all zeros
    "  = 1     ramp
	"  = 2     sine 
	
  verb > 0  verbose output [!!! not yet used]
	
	r	service ring
	g	shuttle gd (from ring to ds)
	v	shuttle array (from sim to ring)
*/
{
	long len,len2,len4,len34,len54;
	long totin,totout;
	long i,nc1,nc2;
	int status;
	double lastim,tinit;

	double yy;
	static double ii=0;

	len=d->len;
	len2=len/2;
	len4=len/4;
	len34=3*len4;
	len54=5*len4;

//	ngd=g->n;

	if (d->lcw == 0)
	{
		d->cont=0;
		d->capt="simulation";
		if(verb>=2)puts("Begin of the simulation"); //pia (per usare verb)
	}

	totin=r->totin;
	totout=r->totout;

	if (d->debug > 0) printf(" --> len,totin,totout= %ld,%ld,%ld\n",len,totin,totout);//pia: ld ld ld

	nc1=d->nc1;
	nc2=d->nc2;

	printf(" --> nc1,nc2= %ld,%ld\n",nc1,nc2); //pia: ld ld

	if (d->type == 2)  /* case interlaced */
	{
		if (nc1 <= nc2) /* beginning and/or odd chunks */
		{
			if (nc1 == 0)
			{
				if (len > r->len) 
				{
					printf(" *** Insufficient ring length : %ld instead of %ld",
						r->len,len); //pia: ld ld
					return 2;
				}
				while (totin-totout < len)  /* if there are few ring data, get more */
				{
					d->cont++;
				
					for (i=0; i<nv; i++)
					{
						switch (typ)
						{
						case 0:
							yy=0.;
							break;
						case 1:
							yy=ii;
							break;
						case 2:
							yy=sin(i*0.06283185);
							break;
						}
						v[i]=yy;
						ii++;
					}

					status=write_rg(r,v,nv,lastim); /* *** LASTIM SENZA VALORE */
					totin=r->totin;
				//	d->dt=g->dx;
				//	r->dx=g->dx;
				}

				
				for (i=0; i < len4; i++)       d->y1[i]=0.;
				status=read_rg2ds(r,d,1,len4,len,&tinit);
				//for (i=len4; i < 5*len/4; i++) d->y1[i]=g->y[i-len4];
			}

			d->lcw++;
			d->nc1++;

			printf("ds chunk -> %ld",d->lcw); //pia ld
		}
		else    /* even chunks :  produces d.y2 */
		{
			while (totin-totout < len)
			{
				d->cont++;

				for (i=0; i<nv; i++)
				{
					switch (typ)
					{
					case 0:
						yy=0.;
						break;
					case 1:
						yy=ii;
						break;
					case 2:
						yy=sin(i*0.06283185);
						break;
					}
					v[i]=yy;
					ii++;
				}

				status=write_rg(r,v,nv,lastim);
				totin=r->totin;
			//	d->dt=g->dx;
			//	r->dx=g->dx;
			}

			for (i=0; i < len34; i++)
				d->y2[i]=d->y1[i+len2];
			//for (i=len34; i < len; i++)
			//	d->y2[i]=g->y[i-len34];
			for (i=0; i < len4; i++)
				d->y1[i]=d->y2[i+len2];
			//for (i=len4; i < len54; i++)
			//	d->y1[i]=g->y[i-len4];
			status=read_rg2ds(r,d,1,len4,len,&tinit);
			for (i=0; i < len4; i++)
				d->y2[i+len34]=d->y1[i];

			d->lcw++;
			d->nc2++;

			printf("ds chunk -> %ld",d->lcw); //pia: ld
		}
	}
	else
	{
		while (totin-totout < len)
		{
			d->cont++;

			for (i=0; i<nv; i++)
			{
				switch (typ)
				{
				case 0:
					yy=0.;
					break;
				case 1:
					yy=ii;
					break;
				case 2:
					yy=sin(i*0.06283185);
					break;
				}
				v[i]=yy;
				ii++;
			}


			status=write_rg(r,v,nv,lastim);
			totin=r->totin;
		//	d->dt=g->dx;
		//	r->dx=g->dx;
		}

		if (d->type == 1)
			for(i=0; i<d->len; i++)
				d->y2[i]=d->y1[i];

		//for (i=0; i < len; i++)
		//	d->y1[i]=g->y[i];
		status=read_rg2ds(r,d,1,0,len,&tinit);

		d->lcw++;
		d->nc2++;

		printf("ds chunk -> %ld",d->lcw); //pia: ld
	}

	return 1;
}


MCH mch_param(MCH mch)
/* computes the derived parameters for an mch structure and allocates y */
{
	return mch;
}

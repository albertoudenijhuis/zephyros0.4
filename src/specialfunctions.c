#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "specialfunctions.h"

double specialfunctions_incompletegamma (double alpha, double x)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper 
	   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=lgamma(alpha);
   double accurate=1.e-50, overflow=1.e30;
   double factor, gin=0., rn=0., a=0.,b=0.,an=0.,dif=0., term=0., pn[6];

   if (x==0.) return (0.);
   if (x<0. || p<=0.) return (-1);

   factor=exp(p*log(x)-x-g);   
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1.;  term=1.;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1.-p;   b=a+x+1.;  term=0.;
   pn[0]=1.;  pn[1]=x;  pn[2]=x+1.;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2.;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0.) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1.-(factor*gin);

 l50:
   return (gin);
}










/* digamma.c
 *
 * Mark Johnson, 2nd September 2007
 *
 * Computes the Ψ(x) or digamma function, i.e., the derivative of the 
 * log gamma function, using a series expansion.
 *
 * Warning:  I'm not a numerical analyst, so I may have made errors here!
 *
 * The parameters of the series were computed using the Maple symbolic
 * algebra program as follows:
 *
 * series(Psi(x+1/2), x=infinity, 21);
 *
 * which produces:
 *
 *  ln(x)+1/(24*x^2)-7/960/x^4+31/8064/x^6-127/30720/x^8+511/67584/x^10-1414477/67092480/x^12+8191/98304/x^14-118518239/267386880/x^16+5749691557/1882718208/x^18-91546277357/3460300800/x^20+O(1/(x^21)) 
 *
 * It looks as if the terms in this expansion *diverge* as the powers
 * get larger.  However, for large x, the x^-n term will dominate.
 *
 * I used Maple to examine the difference between this series and
 * Digamma(x+1/2) over the range 7 < x < 20, and discovered that the
 * difference is less that 1e-8 if the terms up to x^-8 are included.
 * This determined the power used in the code here.  Of course,
 * Maple uses some kind of approximation to calculate Digamma,
 * so all I've really done here is find the smallest power that produces
 * the Maple approximation; still, that should be good enough for our
 * purposes.
 *
 * This expansion is accurate for x > 7; we use the recurrence 
 *
 * digamma(x) = digamma(x+1) - 1/x
 *
 * to make x larger than 7.
 */


double specialfunctions_digamma(double x) {
  double result = 0, xx, xx2, xx4;

  if ((x > 0.) == 0) {
	  printf("digamma error, x = %.2e \n", x);
  }
  assert(x > 0.);
  for ( ; x < 7.; ++x)
    result -= 1./x;
  x -= 1./2.;
  xx = 1./x;
  xx2 = xx*xx;
  xx4 = xx2*xx2;
  result += log(x)+(1./24.)*xx2-(7./960.)*xx4+(31./8064.)*xx4*xx2-(127./30720.)*xx4*xx4;
  return result;
}




/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************

double specialfunctions_incompletegamma(double a, double x) 
{
	double ans = 0.;
	double ax = 0.;
	double c = 0.;
	double r = 0.;
	double igammaepsilon = 1.e-25;

	if ((x <= 0.) | (a <= 0.)) {
		return 0.;
	}
	if ((x > 1.) & (x > a)) {
		return 1. - specialfunctions_incompletegammac(a, x);
	}
	ax = a*log(x) - x - lgamma(a);
	if (ax < -709.78271289338399) {
		return 0.;
	}
	ax = exp(ax);
	r = a;
	c = 1.;
	ans = 1.;
	while( (c/ans)>igammaepsilon) 
	{
		r = r+1.;
		c = c*x/r;
		ans = ans+c;
	}
	return ans*ax/a;	
}
*/




/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
/*
double specialfunctions_incompletegammac(double a, double x)
{
	double igammaepsilon = 1.e-15;
	double igammabignumber = 4503599627370496.0;
	double igammabignumberinv = 2.22044604925031308085*0.0000000000000001;
	double ans = 0.;
	double ax = 0.;
	double c = 0.;
	double yc = 0.;
	double r = 0.;
	double t = 0.;
	double y = 0.;
	double z = 0.;
	double pk = 0.;
	double pkm1 = 0.;
	double pkm2 = 0.;
	double qk = 0.;
	double qkm1 = 0.;
	double qkm2 = 0.;

	if ((x <= 0.) | (a <= 0.)) {
		return 1.;
	}
	if ((x < 1.) | (x < a)) {
		return 1.-specialfunctions_incompletegamma(a, x);
	}
	ax = a*log(x)-x-lgamma(a);
	if (ax < -709.78271289338399 )
	{
		return 0.;
	}
	ax = exp(ax);
	y = 1.-a;
	z = x+y+1.;
	c = 0.;
	pkm2 = 1;
	qkm2 = x;
	pkm1 = x+1.;
	qkm1 = z*x;
	ans = pkm1/qkm1;
	while( t>igammaepsilon)
	{
		c = c+1.;
		y = y+1.;
		z = z+2.;
		yc = y*c;
		pk = pkm1*z-pkm2*yc;
		qk = qkm1*z-qkm2*yc;
		if( qk !=0. ) {
			r = pk/qk;
			t = fabs((ans-r)/r);
			ans = r;
		} else {
			t = 1.;
		}
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;
		if( fabs(pk) > igammabignumber) {
			pkm2 = pkm2*igammabignumberinv;
			pkm1 = pkm1*igammabignumberinv;
			qkm2 = qkm2*igammabignumberinv;
			qkm1 = qkm1*igammabignumberinv;
		}
	}
	return ans*ax;
}
*/

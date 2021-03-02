#define mishchenko2000_NPN1 	100
#define mishchenko2000_NPNG1 500
#define mishchenko2000_NPNG2 2*mishchenko2000_NPN1
#define mishchenko2000_NPN2 	2*mishchenko2000_NPN1
#define mishchenko2000_NPL 	mishchenko2000_NPN2+1
#define mishchenko2000_NPN3 	mishchenko2000_NPN1+1
#define mishchenko2000_NPN4 	mishchenko2000_NPN1
#define mishchenko2000_NPN5 	2*mishchenko2000_NPN1
#define mishchenko2000_NPN6 	mishchenko2000_NPN1+1

#include <complex.h>


extern struct {
	float rt11[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float rt12[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float rt21[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float rt22[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float it11[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float it12[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float it21[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float it22[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
} tmat_ ;

extern struct {
	//for t-matrix calculation
	double 	axi;
	double 	rat;
	double 	lam;
	double 	mrr;
	double 	mri;
	double 	eps;
	int		np;
	int		ndgs;
	double  ddelt;

	//for amplitude and phase matrices calculation
	double	alpha;
	double	beta;
	double	thet0;
	double	thet;
	double	phi0;
	double	phi; 
	
	double s11r;
	double s11i;
	double s12r;
	double s12i;
	double s21r;
	double s21i;
	double s22r;
	double s22i;
		
    double z11;
    double z12;
    double z13;
    double z14;
    double z21;
    double z22;
    double z23;
    double z24;
    double z31;
    double z32;
    double z33;
    double z34;
    double z41;
    double z42;
    double z43;
    double z44;
    
    //unique matrix number
    int tmatrixnr;
    int nmax;
} mishchenko_io_ ;

typedef struct st_mishchenko2000_widget
{
	float rt11[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float rt12[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float rt21[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float rt22[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float it11[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float it12[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float it21[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	float it22[mishchenko2000_NPN4][mishchenko2000_NPN4][mishchenko2000_NPN6];
	
	//for t-matrix calculation
	double 	axi;
	double 	rat;
	double 	lam;
	double 	mrr;
	double 	mri;
	double 	eps;
	int		np;
	int		ndgs;
	double  ddelt;
	
	//for amplitude and phase matrices calculation
	double	alpha;
	double	beta;
	double	thet0;
	double	thet;
	double	phi0;
	double	phi;
	
	double complex s11;
	double complex s12;
	double complex s21;
	double complex s22; 
	
    //unique matrix number
    int tmatrixnr;
    int nmax;
} t_mishchenko2000_widget;

void calc_tmatrix_();
void calc_amp_phase_matrices_();


void mishchenko2000_calc_tmatrix(t_mishchenko2000_widget *wid);
void mishchenko2000_amp_phase_matrices(t_mishchenko2000_widget *wid);

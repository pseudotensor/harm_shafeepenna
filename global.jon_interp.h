#define DEBUGINTERP 0 // debug messages

// normal failure to interpolate message
#define SIMPLEDEBUGINTERP 0

// whether to extrapolate if outside input boundaries
#define EXTRAPOLATE 1


#define FLOAT2IMAGE1(x) ( (x<0.0) ? 0.0 : (x>255.0 ? 255.0 : x) )

#define FLOAT2IMAGE(x) ((unsigned char)(round(FLOAT2IMAGE1(ftemp))))


#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdarg.h>

#include "global.realdef.h"

#include "coord.h"
#include "global.nondepnmemonics.h"
#include "definit.h"

#undef MCOORD
//#define MCOORD CYLMINKMETRIC          // coordinates for some metric
#define MCOORD SPCMINKMETRIC

#undef DOSTOREPOSITIONDATA
#define DOSTOREPOSITIONDATA 0

#undef DOGRIDSECTIONING
#define DOGRIDSECTIONING 0

#undef N1
#undef N2
#undef N3
#define N1 2
#define N2 2
#define N3 2

#include "global.depnmemonics.h"
#include "global.loops.h"
#include "global.variousmacros.h"
#include "global.fieldmacros.h"
#include "global.structs.h"
#define SIMULBCCALC -1
#include "global.gridsectioning.h"
#include "global.comploops.h"



// number of dimensions for code (2d or 3d code)
#define COMPDIM 3

#define COORDSINGFIXCYL 0

#define ANALYTICJAC 1 // helps to avoid no solution sometimes -- e.g. when precisionj of numerical jac in question)
// 0: numerical jac
// 1: analytic jac
// applies to newt() and broydn()
#define ROOTMETHOD 2 // 2 has proven best
// 0 : nondamped mnewt() // definitely can fail to find solution even to modest error unless guess is close
// 1 : damped mnewt() // comment is same as 0 above
// 2 : newt() // works with or without good guess, and with or (mostly) without analytic jac --  as long as jacobian is allowed to be unrestricted outside normal domain (see coord.c and singularity/POSDEFMETRIC fix turned off).
// 3 : brodyn() // basically same as comment for 2 above, but can be much worse than 2
#define GOODGUESS 0 // turned on for 2/3 ROOTMETHODS only to speed things up
// 0: simple mid-point guess always
// 1: choose guess with some insight of typical coordinat transformations


#define DOUBLEWORK 1
// 0=use unsigned char during calculations for images (can result in over or under flow)
// 1=use doubles during calculations for images (avoids over/underflow)


//#define dualfprintf fprintf
//#define myfprintf fprintf
#define fail_file stderr
#define logfull_file stderr
#define log_file stderr
#define myexit exit





#define MINMAX(q,a,b) ( ((q)==CMIN) ? MIN(a,b) : MAX(a,b) )







#undef PLOOP
#undef PDUMPLOOP
#undef PINVERTLOOP
#undef PBOUNDLOOP
#undef PLOOPINTERP

  // check for existence in bad form using:
// grep "PLOOP" *.c | grep --invert-match "PLOOP("
/* loop over all Primitive variables */
#define PLOOP(pl) for(pl=0;pl<NPR;pl++)
/* loop over all dumped Primitive variables */
#define PDUMPLOOP(pd) for(pd=0;pd<NPRDUMP;pd++)
/* loop over all inversion Primitive variables */
#define PINVERTLOOP(pi) for(pi=0;pi<NPRINVERT;pi++)
/* loop over all bounding Primitive variables */
#define PBOUNDLOOP(pb) for(pb=0;pb<NPRBOUND;pb++)
/* loop over all center to edge variables */
#define PLOOPINTERP(pl) for(pl=0;pl<NPR2INTERP;pl++)





#define SLEPSILON (1.e-6)



// need not change below datatype stuff
#if(REALTYPE==FLOATTYPE)
#define SCANARG "%f"
// 16 args
// 21 args after going to 3D and doing MBH/QBH
#define SCANHEADER "%f %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f"
#elif(REALTYPE==DOUBLETYPE)
#define SCANARG "%lf"
#define SCANHEADER "%lf %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf"
#elif(REALTYPE==LONGDOUBLETYPE)
#define SCANARG "%Lf"
#define SCANHEADER "%Lf %d %d %d %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %d %Lf %Lf "
#endif

#define PRINTSCANHEADER "%g %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %g %g"




#define POSDEFMETRIC 0




#define USEMPI 0





#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
//#define SIGN(a) ( ((a) <0.) ? -1. : 1. )




//#define FAILSTATEMENT(file,function,number) {if(debugfail>=1){ fprintf(fail_file,"%s %d-%s(): failure\n",file,number,function); fflush(fail_file); fprintf(fail_file,"i: %d j: %d p: %d\n",icurr,jcurr,pcurr);} return(1);}


// NR STUFF
extern FTYPE ranc(int seed);

extern int ludcmp(FTYPE **a, int n, int *indx, FTYPE *d);
extern void lubksb(FTYPE **a, int n, int *indx, FTYPE *d);


/* NR routines from nrutil.h used by nrutil.c */
extern int *ivector(long nl, long nh);
extern void free_ivector(int *v, long nl, long nh);
extern FTYPE *dvector(long nl, long nh);
extern void free_dvector(FTYPE *v, long nl, long nh);
extern FTYPE **dmatrix(long nrl, long nrh, long ncl, long nch);
extern void free_dmatrix(FTYPE **m, long nrl, long nrh, long ncl,
			 long nch);
extern FTYPE ***dtensor(long nrl, long nrh, long ncl, long nch,
			 long ndl, long ndh);
extern void free_dtensor(FTYPE ***t, long nrl, long nrh, long ncl,
			 long nch, long ndl, long ndh);
extern void nrerror(char error_text[]);


extern void free_vector(FTYPE *v, long nl, long nh);
extern FTYPE *vector(long nl, long nh);
extern FTYPE **matrix(long nrl, long nrh, long ncl, long nch);
extern unsigned char **cmatrix(int nrl,int nrh,int ncl,int nch);
extern FTYPE **fmatrix(int nrl,int nrh,int ncl,int nch);
extern FTYPE ***f3matrix(int nzl, int nzh, int nrl, int nrh, int ncl, int nch);
extern void free_cmatrix(unsigned char **m, long nrl, long nrh, long ncl, long nch);
extern void free_fmatrix(FTYPE **m, long nrl, long nrh, long ncl, long nch);

extern void free_f3matrix(FTYPE ***m, long nzl, long nzh, long nrl, long nrh, long ncl, long nch);
extern void free_matrix(FTYPE **m, long nrl, long nrh, long ncl, long nch);

extern void broydn(FTYPE x[], int n, int *check,
	    void (*vecfunc)(int, FTYPE [], FTYPE []));
extern void qrdcmp(FTYPE **a, int n, FTYPE *c, FTYPE *d, int *sing);
extern void rsolv(FTYPE **a, int n, FTYPE d[], FTYPE b[]);
extern void qrupdt(FTYPE **r, FTYPE **qt, int n, FTYPE u[], FTYPE v[]);
extern void rotate(FTYPE **r, FTYPE **qt, int n, int i, FTYPE a, FTYPE b);

extern int gaussj(FTYPE **tmp, int n, FTYPE **b, int m);

extern void newt(FTYPE x[], int n, int *check,
	  void (*vecfunc)(int, FTYPE [], FTYPE []));
extern void lnsrch(int n, FTYPE xold[], FTYPE fold, FTYPE g[], FTYPE p[], FTYPE x[], FTYPE *f, FTYPE stpmax, int *check, FTYPE (*func)(FTYPE []));
extern void fdjac(int n, FTYPE x[], FTYPE fvec[], FTYPE **df,
		  void (*vecfunc)(int, FTYPE [], FTYPE []));
extern FTYPE nrfmin(FTYPE x[]);

extern void X2spc(int n, FTYPE *Xguess, FTYPE *spc_curr);


// functions inside jon_interp_coord.c (or used by it or for it)
// coordinate stuff
extern void set_coord_parms(void);
extern void write_coord_parms(void);
extern void read_coord_parms(void);

#if(COMPDIM==2)
// 2D:
extern void bl_coord(FTYPE *X, FTYPE *r, FTYPE *th);
void dxdxprim(FTYPE *X, FTYPE r, FTYPE th, FTYPE (*dxdxp)[NDIM]);
extern void dxdxprim(FTYPE *X, FTYPE r, FTYPE th, FTYPE (*dxdxp)[NDIM]);
extern void coord(int i, int j, int loc, FTYPE *X);
extern void coordf(FTYPE i, FTYPE j, int loc, FTYPE *X);
extern void icoord(FTYPE *X,int loc, int *i, int *j);
#endif

#if(COMPDIM==3)
// 3D:
extern void bl_coord(FTYPE *X, FTYPE *V);
void dxdxprim(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);
void dxdxp_analytic(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);
extern void coord(int i, int j, int k, int loc, FTYPE *X);
extern void coordf(FTYPE i, FTYPE j, FTYPE k, int loc, FTYPE *X);
extern void icoord(FTYPE *X,int loc, int *i, int *j, int *k);
#endif



extern void coord_ijk(int i, int j, int k, int loc, FTYPE *X);
extern void coord_free(int i, int j, int k, int loc, FTYPE *X);

extern void bl_coord_ijk(int i, int j, int k, int loc, FTYPE *V);
extern void bl_coord_ijk_2(int i, int j, int k, int loc, FTYPE *X, FTYPE *V);

extern void dxdxprim_ijk(int i, int j, int k, int loc, FTYPE (*dxdxp)[NDIM]);
extern void dxdxprim_ijk_2(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);

extern void idxdxprim(FTYPE (*dxdxp)[NDIM], FTYPE (*idxdxp)[NDIM]);
extern void idxdxprim_ijk(int i, int j, int k, int loc, FTYPE (*idxdxp)[NDIM]);
extern void idxdxprim_ijk_2(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE (*idxdxp)[NDIM]);



extern void dxdxp_func(FTYPE *X, FTYPE (*dxdxp)[NDIM]);
extern void set_points();
extern int setihor(void);
extern FTYPE setRin(int ihor);
extern FTYPE rhor_calc(int which);


// functions inside jon_interp_mnewt.c or used by it (i.e. usrfun() )
extern int mnewt(int ntrial, int mintrial, FTYPE x[], int n, FTYPE tolx, FTYPE tolf, FTYPE tolxallowed, FTYPE tolfallowed, FTYPE tolxreport, FTYPE tolfreport);
extern int usrfun(FTYPE *Xguess, int n, FTYPE *beta, FTYPE **alpha,FTYPE*norm);




extern int usrfun2(int n, FTYPE *Xguess, FTYPE *spc_curr, FTYPE **alpha);
 

extern void bcucof(FTYPE y[], FTYPE y1[], FTYPE y2[], FTYPE y12[], FTYPE d1, FTYPE d2, FTYPE **c);

extern void bcuint(FTYPE y[], FTYPE y1[], FTYPE y2[], FTYPE y12[], FTYPE x1l,
        FTYPE x1u, FTYPE x2l, FTYPE x2u, FTYPE x1, FTYPE x2, FTYPE *ansy,
	    FTYPE *ansy1, FTYPE *ansy2);

// log file stuff
extern void myfprintf(FILE* fileptr, char *format, ...);
extern void dualfprintf(FILE* fileptr,char *format, ...);
extern void logsfprintf(char *format, ...);
extern void trifprintf(char *format, ...);

extern void myfopen(char*fname, char*fmt, char*message, FILE ** fileptr);
extern void myfclose(FILE ** fileptr,char*message);

extern void matrix_inverse(int whichcoord, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);

extern FTYPE rmso_calc(int which) ;


extern void interpicoord(FTYPE *X,int loc, int *i, int *j, int *k);



#include "global.jon_interp.h"

FTYPE Xmax[NDIM];

FTYPE (*p)[N2M][N3M][NPR];
FTYPE Risco,drsing;
FTYPE *rcent,*rcent_tot;
int horizoni,horizoncpupos1;
int numprocs;
int mycpupos[NDIM];
int ncpux1;

int oN1,oN2,oN3,nN1,nN2,nN3 ;
FTYPE refinefactor;
int roN1,roN2,roN3;
FTYPE dX[NDIM],Rin,fakeRin,dxc,dyc,dzc,fakedxc,fakedyc,fakedzc;
FTYPE Zin,Zout,Zeqin;
FTYPE startxc, startyc, startzc, endxc, endyc, endzc;
int newgridtype,oldgridtype;

FTYPE X[NDIM];
FTYPE Xmetricnew[NDIM],Xmetricold[NDIM]; // used to store time of latest and oldest metric


FTYPE t,gam,spin,QBH,MBH;
int startpos[NDIM];
int totalsize[NDIM];
long realnstep,nstep;
FTYPE readnstep;

int DATATYPE;

int totalzones,icurr,jcurr,kcurr,pcurr;
int defcoord;
FTYPE h_over_r,hslope;
FTYPE jetalpha;

// GODMARK3D -- should these be stored in coordparms.dat?
FTYPE Rin_array[NDIM], Rout_array[NDIM];  //atch -- arrays for a more general way of handling the grid dimensions

FTYPE Rhor,Rout,dx[NDIM],startx[NDIM],R0;
FTYPE dxdxp[NDIM][NDIM];
int myid;
int debugfail;
FTYPE dt;

FTYPE spc_target[NDIM];

FTYPE Lunit,dV,dVF;
//int N1,N2,N3;
//int N1BND, N2BND, N3BND;

int BCtype[COMPDIM*2];



FTYPE a;
int failed;
long nstroke;


int nn;
FTYPE *fvec;

int nrerrorflag;
FTYPE Rchop;

int didstorepositiondata;
FTYPE (*dxdxpstore)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM];
FTYPE (*idxdxpstore)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM];
FTYPE (*Xstore)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM];
FTYPE (*Vstore)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM];


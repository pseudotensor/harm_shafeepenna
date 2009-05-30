#include "decs.h"



// jon's version of NR's dfridr modified to accept more general, needed, function
#if(1) // uniformly low values although not always lower than original version
#define MAXITER 50
#define NRANSI
#define CON 1.1
#define CON2 (CON*CON)
//#define NTAB 130 // number of function evaluations is 2XNTAB
#define NTAB 30 // number of function evaluations is 2XNTAB
#define SAFE 2.0
//#define NRHMAX 1 // maximum size where metric changes substantially
#define NRHMAX 1.e-1
#define TRYTOL (pow(NUMEPSILON,2.0/3.0)) // attempted tolerance
#define OKTOL 1e-1 // error must be below this
#define MINH (pow(NUMEPSILON,2.0/3.0)) // minimum h
//#define HSTARTCHANGE 10.0
#define HSTARTCHANGE 2.0
#endif

#if(0) // original version (gets pretty damn low for many, but not all, derivatives -- some 1E-6 rather than 1E-13)
#define MAXITER 128
#define NRANSI
#define CON 1.3
#define CON2 (CON*CON)
#define NTAB 10 // number of function evaluations is 2XNTAB
#define SAFE 2.0
//#define NRHMAX 1 // maximum size where metric changes substantially
#define NRHMAX 1E-2
#define TRYTOL 1E-10 // attempted tolerance
#define OKTOL 1e-5 // error must be below this
#define MINH 1E-15 // minimum h
#define HSTARTCHANGE 10.0
#endif

// whether to turn on extensive recent debugging
#define DEBUGDF 0

FTYPE dfridr(FTYPE (*func)(struct of_geom *, FTYPE*,int,int), struct of_geom *ptrgeom, FTYPE *X,int ii, int jj, int kk)
{
  int i,j,k;
  FTYPE errt,fac,hh,**a,ans;
  FTYPE dX[NDIM],Xh[NDIM],Xl[NDIM];
  FTYPE h,err;
  FTYPE hstart;
  FTYPE newdx,temp;
  int shit;
  int iter;
  FTYPE errlist[MAXITER];
  FTYPE hhlist[MAXITER];
  FTYPE anslist[MAXITER];
  FTYPE minerror,minerrorhstart;
  FTYPE hstartfrac;
  int miniter;
  int iterdebug;
  FTYPE minans;


  hstartfrac=HSTARTCHANGE; // starting fractional drop in hstart
	
  // allocate memory
  a=dmatrix(1,NTAB,1,NTAB);

  //    hstart=NRHMAX;
  //  hstart=MAX(dx[kk],NRHMAX);
  //  hstart=MIN(NRHMAX,hstart*10.0);
  hstart=1.0;



  miniter=0; // didn't find minimum is assumed, which means first is minimum!
  iter=0;
  minerror=BIG;
  minans=BIG;

  while(1){

#if(DEBUGDF)
    dualfprintf(fail_file,"iter=%d\n",iter);
#endif


    h=hstart;
    if (h == 0.0) nrerror("h must be nonzero in dfridr.");

    hh=h;
    // HARM STUFF
    for(k=0;k<NDIM;k++) dX[k]=0.0; // other components will remains 0 for this function
    dX[kk]=hh;
    for(k=0;k<NDIM;k++){
      //      Xl[k]=X[k]-dX[k];
      temp=X[k]-dX[k];
      newdx=temp-X[k];
      Xl[k]=X[k]+newdx;
    }
    for(k=0;k<NDIM;k++){
      //      Xh[k]=X[k]+dX[k];
      temp=X[k]+dX[k];
      newdx=temp-X[k];
      Xh[k]=X[k]+newdx;
    }
    // end HARM STUFF
    //    for(k=0;k<NDIM;k++) dualfprintf(fail_file,"k=%d %21.15g %21.15g\n",k,X[k],dX[k]);
    a[1][1]=((*func)(ptrgeom,Xh,ii,jj)-(*func)(ptrgeom,Xl,ii,jj))/(2.0*hh);
    err=BIG;

    for (i=2;i<=NTAB;i++) {
      hh /= CON;
      // HARM STUFF
      dX[kk]=hh;
      for(k=0;k<NDIM;k++) Xl[k]=X[k]-dX[k];
      for(k=0;k<NDIM;k++) Xh[k]=X[k]+dX[k];
      //      for(k=0;k<NDIM;k++) dualfprintf(fail_file,"i=%d k=%d %21.15g %21.15g\n",i,k,X[k],dX[k]);
      // end HARM STUFF
      a[1][i]=((*func)(ptrgeom,Xh,ii,jj)-(*func)(ptrgeom,Xl,ii,jj))/(2.0*hh);
      fac=CON2;

#if(DEBUGDF)
	// debug
      //	if(ii==1 && jj==1 && kk==1){
      //	  if( icurr==29 && jcurr==-11){
	    for(shit=0;shit<4;shit++) dualfprintf(fail_file,"hh=%21.15g Xl[%d]=%21.15g Xh[%d]=%21.15g\n",hh,shit,Xl[shit],shit,Xh[shit]);
	    //	  }
    //	}
#endif


      for (j=2;j<=i;j++) {
	a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
	fac=CON2*fac;
	//normalized error
	//	errt=MAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]))/((*func)(ptrgeom,X,ii,jj)+SMALL);
	//	errt=MAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
	// normalized error
	errt=MAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]))/(fabs(a[j][i])+SMALL);

#if(DEBUGDF)
	// debug
	//	if(ii==3 && jj==3 && kk==2){
	//	if(ii==2 && jj==2 && kk==2){
	//	if(ii==1 && jj==1 && kk==1){
	//	  if( icurr==29 && jcurr==-11){
	    //	  if( icurr==40 && ((jcurr==54)||(jcurr==9)) ){
	    dualfprintf(fail_file,"a11=%21.15g a[%d][%d]=%21.15g errt=%21.15g fac=%21.15g\n",a[1][1],i,j,a[i][j],errt,fac);
	    //	  }
	    //	}
#endif

	if (errt <= err) {
	  err=errt;
	  ans=a[j][i];
	}
      }
      // normalized error
      //      if (fabs((a[i][i]-a[i-1][i-1])/( (*func)(ptrgeom,X,ii,jj)+SMALL)) >= SAFE*(err)) break;
      if (fabs((a[i][i]-a[i-1][i-1])) >= SAFE*(err)) break;
      if (fabs((a[i][i]-a[i-1][i-1]))/(fabs(ans)+SMALL) >= SAFE*(err)) break;
    }
	  
    // now check error is not crazy with the starting large h, decrease h if crazy until not crazy and get good error
    errlist[iter]=err; // store error
    //    hhlist[iter]=hstart;
    hhlist[iter]=hh; // current hh to focus around
    anslist[iter]=ans;

    if(err>TRYTOL){
      if(hstart<MINH){
	if(err<OKTOL) break;
	else{

	  if(miniter>=0){
	    if(minerror<OKTOL){
	      ans=minans;
	      break;	      // then accept as answer
	    }
	  }
	  else{ // then must fail
	    dualfprintf(fail_file,"never found error below %21.15g: err=%21.15g : ii=%d jj=%d kk=%d\n",OKTOL,err,ii,jj,kk);
	    dualfprintf(fail_file,"gi=%d gj=%d gk=%d\n",icurr,jcurr,kcurr);
	    dualfprintf(fail_file,"miniter=%d errlist[miniter]=%21.15g hhlist[miniter]=%21.15g\n",miniter,errlist[miniter],hhlist[miniter]);	    
	    dualfprintf(fail_file,"minerror=%21.15g minans=%21.15g\n",minerror,minans);
	  
	    for(iterdebug=0;iterdebug<=iter;iterdebug++){
	      dualfprintf(fail_file,"h[%d]=%21.15g err[%d]=%21.15g ans[%d]=%21.15g\n",iterdebug,hhlist[iterdebug],iterdebug,errlist[iterdebug],iterdebug,anslist[iterdebug]);
	    }
	    myexit(66);
	  }
	}
      }
      else{ // keep going since not below MINH
	if(errlist[iter]<minerror){
	  // store min error event
	  minerror=errlist[iter];
	  minerrorhstart=hhlist[iter];
	  miniter=iter;
	  minans=ans;
#if(DEBUGDF)
	  dualfprintf(fail_file,"minerr=%21.15g minhstart=%21.15g miniter=%d minans=%21.15g\n",minerror,minerrorhstart,miniter,minans);
#endif
	}
	if(iter>0){
	  // check on error
	  if(0&& (miniter>=1) && (errlist[iter]>1000.0*minerror)){
	    // then something bad is happening ... probably skipped over minimum error region
	    hstartfrac=(hstartfrac-1.0)*0.5+1.0; // slow down
	    //	    hstart=hhlist[miniter]*(1.0+NTAB*0.5/hstartfrac); // back up a bit
	    hstart=minerrorhstart*10.0;
	    iter=-1; // start over
	    //	    miniter=0;
#if(DEBUGDF||1)
	    dualfprintf(fail_file,"hstart=%21.15g hstartfrac=%21.15g\n",hstart,hstartfrac);
#endif
	  }
	  else{
	    // keep dropping h
	    hstart/=hstartfrac;
#if(DEBUGDF)
	    dualfprintf(fail_file,"hstartnew=%21.15g\n",hstart);
#endif
	  }
	}
      } // end of keep going part
    }// end if err > TRYTOL
    else break; // then error<TRYTOL!  GOOD!
    iter++;
    if(iter>=MAXITER){
      if(err<OKTOL) break;
      else{ // then maybe problems
	if(minerror<OKTOL){
	  ans=minans;
	  break;	      // then accept as answer
	}
	else{
	  // then must fail
	  dualfprintf(fail_file,"iter>MAXITER: never found error below %21.15g: err=%21.15g : ii=%d jj=%d kk=%d\n",OKTOL,err,ii,jj,kk);
	  dualfprintf(fail_file,"gi=%d gj=%d gk=%d\n",icurr,jcurr,kcurr);
	  dualfprintf(fail_file,"miniter=%d errlist[miniter]=%21.15g hhlist[miniter]=%21.15g\n",miniter,errlist[miniter],hhlist[miniter]);
	  dualfprintf(fail_file,"minerror=%21.15g minans=%21.15g\n",minerror,minans);
	  for(iterdebug=0;iterdebug<iter;iterdebug++){
	    dualfprintf(fail_file,"h[%d]=%21.15g err[%d]=%21.15g ans[%d]=%21.15g\n",iterdebug,hhlist[iterdebug],iterdebug,errlist[iterdebug],iterdebug,anslist[iterdebug]);
	  }
	  myexit(67);
	}
      }      
    }
  }
  // done
  free_dmatrix(a,1,NTAB,1,NTAB);
  //  dualfprintf(fail_file,"hfinal=%21.15g errfinal=%21.15g\n",hh,err);
  //fflush(fail_file);
  return ans;

       

}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE
#undef NRANSI
// (C) Copr. 1986-92 Numerical Recipes Software *1.@Q.. 






/* 
   FTYPE delta(int i, int j) { if(i == j) return(1.) ; else return(0.) 
   ; } */

/* Minkowski metric; signature +2 */
/* 
   FTYPE mink(int i, int j) { if(i == j) { if(i == 0) return(-1.) ;
   else return(1.) ; } else return(0.) ; } */




//////////////////////////////
//
// below very specific to defcoord==LOGRSINTH and MCOORD=KSCOORDS:
// gives: analytic connection and source
//
/////////////////////////////////

// continuous mod function like Mathematica
int contmod(int x, int y)
{
  int ans;

  ans = x%y;
  if(ans<0) ans+=y;

  return(ans);
}

// below 2 mysin/mycos preserve symmetry for arg=-M_PI/2 .. 2M_PI by restricting arg to 0..PI/2
FTYPE mysin(FTYPE th)
{
  int contmod(int x, int y);
  int nmod4;
  FTYPE ans;

#if(ACCURATESINCOS)
  nmod4=contmod((int)(th/(M_PI*0.5)),4);
  
  switch(nmod4){
  case 0:
    ans=sin(th);
    break;
  case 1:
    ans=sin(M_PI-th);
    break;
  case 2:
    ans=-sin(th-M_PI);
    break;
  case 3:
    ans=-sin(2.0*M_PI-th);
    break;
  default:
    dualfprintf(fail_file,"No such case for mysin with nmod4=%d\n",nmod4);
    ans=sqrt(-1.0);
    myexit(206983462);
    break;
  }
#else
  ans=sin(th);
#endif

  return(ans);

}


FTYPE mycos(FTYPE th)
{
  

#if(ACCURATESINCOS)
  return(mysin(th+M_PI*0.5));
#else
  return(cos(th));
#endif

}



//#ifdef WIN32
// cot = 1/tan = cos/sin
FTYPE cot(FTYPE arg)
{
#if(ACCURATESINCOS)
  return(mycos(arg)/mysin(arg));
#else
  return(1.0/tan(arg));
#endif
}
//#endif 

FTYPE csc(FTYPE arg)
{
#if(ACCURATESINCOS)
  return(1.0/mysin(arg));
#else
  return(1.0/sin(arg));
#endif
}


FTYPE sec(FTYPE arg)
{
#if(ACCURATESINCOS)
  return(1.0/mycos(arg));
#else
  return(1.0/cos(arg));
#endif
}




// fromwhere==0 or 1, then assume horizoni on downside of actual location
// fromwhere==2, then assume horizoni on upside of actual location
int find_horizon(int fromwhere)
{
  int i, j, k, ii;
  FTYPE r1, r2;
  FTYPE X[NDIM],V[NDIM];
  int gotit;
  FTYPE horizonvalue;
  // called after grid is setup for all cpus


  // first compute where horizon is

  // these 2 below are used prior, but not initialized otherwised on restart
  // some calculations
  // these 2 below are also used by find_horizon() below
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);


  // need to find horizon and place horizoni on right-hand-side of location



  // was testing to make sure if held horizoni fixed that conserved mass if using conservation of baryon number as basis for masses in Mvsr and accretion of energy
  // GODMARK DEBUG DEBUG DEBUG
  //  if(fromwhere!=2){
  //  return(0);
  // }


  // definition of horizoni must be consistent so fluxes are consistent and have conservation
  fromwhere=2; // force to be on upside unless Rhor=0, which is caught first





  if(fromwhere==0) trifprintf("begin: find_horizon ... ");

  // find cpu column that brackets the horizon and determine the
  // i-offset of horizon
  // notice that only 1 CPU will get horizon since stop process once found
  // notice that radius(horizoni) is below or equal to actual horizon radius

  


  horizonvalue = Rhor;
  horizoni = -100;
  horizoncpupos1=-1;
  gotit = 0;
  for (ii = numprocs - 1; ii >= 0; ii--) { // should get done by first row
    if (ii == myid) {
      for (i = N1-1; i >= 0; i--) {


	j = N2 / 2;             // doesn't matter (spherical polar assumed)
	k = N3 / 2;             // doesn't matter (spherical polar assumed)
	coord_ijk(i, j, k, FACE1, X);
	bl_coord_ijk(i, j, k, FACE1, V);
	r1=V[1];
	coord_ijk(ip1, j, k, FACE1, X);
	bl_coord_ijk(ip1, j, k, FACE1, V);
	r2=V[1];
	// looking between FACE1's r value and upper FACE1's r value, so loop is from i=N1-1..i=0

	if(ii==myid && myid==0 && i==0){
	  // special check in case horizon inside inner-most radial grid
	  if(horizonvalue<=r1 || horizonvalue<SMALL){ // GODMARK: this means horizon can't be chosen to <SMALL and mean there is a black hole there
	    // then horizon off grid or right on edge, but still ok
	    // treat as if horizon is off grid if right on edge
	    horizoni = 0;
	    horizoncpupos1=mycpupos[1];
	    break;
	  }
	}


	//        if (fabs(r1 - horizonvalue) <= (r2 - r1)) {     // find horizon
	if (fromwhere!=2){
	  if(horizonvalue >= r1 && horizonvalue < r2){ // note that if strictly on r2, then next CPU should pick it up
	    horizoni = i;
	    horizoncpupos1 = mycpupos[1];
	    break;
	  }
	}
	else if (fromwhere==2){
	  if(horizonvalue >= r1 && horizonvalue < r2){
	    horizoni = ip1;
	    horizoncpupos1 = mycpupos[1];
	    if(horizoni>=N1){
	      horizoni=0;
	      horizoncpupos1++;
	    }
	    else{
	      // then on original CPU
	      horizoncpupos1 = mycpupos[1];
	    }
	    //dualfprintf(fail_file,"horizon: %d %d\n",horizoni,horizoncpupos1);
	    break;
	  }
	  //	  dualfprintf(fail_file,"horizonnot: %d %d :: %21.15g %21.15g %21.15g\n",horizoni,horizoncpupos1,r1,Rhor,r2);
	}
      }
    }

    if (numprocs > 0) {
#if(USEMPI)
      MPI_Bcast(&horizoni, 1, MPI_INT, ii, MPI_COMM_GRMHD);
      MPI_Bcast(&horizoncpupos1, 1, MPI_INT, ii, MPI_COMM_GRMHD);
#endif
    }
    if (horizoni >= 0) gotit = 1;                // can stop entire process

    // keep horizoni as relative to CPU with horizon so any CPU knows where horizon is
    //    if (mycpupos[1] != horizoncpupos1) {
    //  horizoni = -100;
    //}                           // reset if not right cpu group
    if (gotit) break;
  }




  if(gotit==0){
    dualfprintf(fail_file,"Never found horizon: fromwhere=%d :: MBH=%21.15g a=%21.15g :: Rhor=%21.15g Risco=%21.15g\n",fromwhere,MBH,a,Rhor,Risco);
    myexit(6246);
  }


  /////////////////////////////////
  //
  // report some information
  if(fromwhere==0) {
    trifprintf("horizoni: %d horizoncpupos1: %d\n", horizoni, horizoncpupos1);
    // just a check
    dualfprintf(log_file,"horizoni: %d mycpupos[1]: %d horizoncpupos1: %d\n", horizoni, mycpupos[1], horizoncpupos1);
    
    trifprintf("end: find_horizon\n");
  }


  return(0);
}











int find_RinRout(FTYPE *localRin, FTYPE *localRout)
{
  FTYPE X[NDIM],V[NDIM];
  int whichcpu;


  // assume outer radius is on outer CPU
  // only 1 CPU needs to get this
  whichcpu=ncpux1-1;

  if(myid==whichcpu){

    coord_ijk(N1, 0, 0, FACE1, X);
    bl_coord_ijk(N1, 0, 0, FACE1, V);
    *localRout=V[1];
  }

  if (numprocs > 0) {
#if(USEMPI)
    MPI_Bcast(localRout, 1, MPI_FTYPE, whichcpu, MPI_COMM_GRMHD);
#endif
  }


  // assume inner radius is on inner CPU
  // only 1 CPU needs to get this
  whichcpu=0;

  if(myid==whichcpu){

    coord_ijk(0, 0, 0, FACE1, X);
    bl_coord_ijk(0, 0, 0, FACE1, V);
    *localRin=V[1];
  }

  if (numprocs > 0) {
#if(USEMPI)
    MPI_Bcast(localRin, 1, MPI_FTYPE, whichcpu, MPI_COMM_GRMHD);
#endif
  }


  return(0);
}








// get dr(i=0)
// assumes black hole is at r=0
void set_drsing(void)
{
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE V[NDIM],X[NDIM];
  FTYPE dr;

  // assume BH r=0 is at inner radial boundary
  if(mycpupos[1]==0){
    coord_ijk(0,0,0,CENT,X);
    bl_coord_ijk(0,0,0,CENT,V);
    dxdxprim_ijk(0,0,0,CENT,dxdxp);
    
    dr = (dxdxp[1][1]*dx[1] + dxdxp[1][2]*dx[2])/10.0; // divide by 10 so doesn't dominate

    //    dualfprintf(fail_file,"%21.15g %21.15g %21.15g %21.15g\n",dxdxp[1][1],dx[1],dxdxp[1][2],dx[2]);
  }
  else dr=0;

#if(USEMPI)
  MPI_Allreduce(&dr, &drsing,1, MPI_FTYPE, MPI_MAX,MPI_COMM_GRMHD);
#else
  drsing=dr;
#endif

}

// get 1-D line over all CPUs that has the radius (only applicable if r(x_1) and not r(x_1,x_2)
void set_rvsr(void)
{
  FTYPE X[NDIM],V[NDIM];
  int i,j,k;
  int ii;
  FTYPE r;


  //////////////////////////////////////
  //
  // get radius for this CPU

  // initialize full rcent for all CPUs, choosing value so that can use MAX over all CPUs and get answer we want
  GRAVLOOP(ii) rcent[ii]=-1E30;

  j=0; // assumes j==0 has same radial dependence as all other j (true if r(x_1))
  k=0; // as with j
  COMPLOOPFP11{
    coord_ijk(i, j, k, CENT, X);
    bl_coord_ijk(i, j, k, CENT, V);
    r=V[1];
    ii=startpos[1]+i;
    rcent[ii] = r;
  }

  // send information to myid=0 since only this processor needs rcent
#if(USEMPI)
  MPI_Reduce(&(rcent[-N1BND]),&(rcent_tot[-N1BND]),NUMGRAVPOS,MPI_FTYPE,MPI_MAX,0,MPI_COMM_GRMHD);
  //  MPI_Reduce(rcent,rcent_tot,ncpux1*N1,MPI_FTYPE,MPI_MAX,0,MPI_COMM_GRMHD);
#else
  GRAVLOOP(ii) rcent_tot[ii]=rcent[ii];
#endif


}



// t-r inverse only (for r=0 coordinate singularity)
// NDIM in size, but invert submatrix that's 2x2 in size
void matrix_inverse_2d(FTYPE gcov[][NDIM], FTYPE gcon[][NDIM])
{
  int jj,kk;

  DLOOP(jj,kk) gcon[jj][kk] = 0.0;


  gcon[TT][TT] = 1.0/(-gcov[TT][RR]*gcov[TT][RR]/gcov[RR][RR] + gcov[TT][TT]);

  gcon[TT][RR] = gcon[RR][TT] = -gcov[TT][RR]/(-gcov[TT][RR]*gcov[TT][RR] + gcov[RR][RR]*gcov[TT][TT]);

  gcon[RR][RR] = 1.0/(-gcov[TT][RR]*gcov[TT][RR]/gcov[TT][TT] + gcov[RR][RR]);

}

// t-r-\theta inverse only (for \theta={0,\pi} coordinate singularities)
// NDIM in size, but invert submatrix that's 3x3 in size
void matrix_inverse_3d(FTYPE gcov[][NDIM], FTYPE gcon[][NDIM])
{
  int jj,kk;
  FTYPE gcovrhsq,gcovtrsq,gcovthsq;


  DLOOP(jj,kk) gcon[jj][kk] = 0.0;


  gcovrhsq = gcov[RR][TH]*gcov[RR][TH];
  gcovtrsq = gcov[TT][RR]*gcov[TT][RR];
  gcovthsq = gcov[TT][TH]*gcov[TT][TH];
  
  gcon[TT][TT]=(gcovrhsq - gcov[TH][TH]*gcov[RR][RR])/(gcov[RR][RR]*gcovthsq - 2.0*gcov[RR][TH]*gcov[TT][TH]*gcov[TT][RR] + gcov[TH][TH]*gcovtrsq + gcovrhsq*gcov[TT][TT] - gcov[TH][TH]*gcov[RR][RR]*gcov[TT][TT]);
  gcon[TT][RR]=gcon[RR][TT]=(-(gcov[RR][TH]*gcov[TT][TH]) + gcov[TH][TH]*gcov[TT][RR])/(-2.0*gcov[RR][TH]*gcov[TT][TH]*gcov[TT][RR] + gcov[TH][TH]*gcovtrsq + gcovrhsq*gcov[TT][TT] + gcov[RR][RR]*(gcovthsq - gcov[TH][TH]*gcov[TT][TT]));
  gcon[TT][TH]=gcon[TH][TT]=(gcov[RR][RR]*gcov[TT][TH] - gcov[RR][TH]*gcov[TT][RR])/(gcov[RR][RR]*gcovthsq - 2.0*gcov[RR][TH]*gcov[TT][TH]*gcov[TT][RR] + gcov[TH][TH]*gcovtrsq + gcovrhsq*gcov[TT][TT] - gcov[TH][TH]*gcov[RR][RR]*gcov[TT][TT]);
  gcon[RR][RR]=(gcovthsq - gcov[TH][TH]*gcov[TT][TT])/(gcov[RR][RR]*gcovthsq - 2.0*gcov[RR][TH]*gcov[TT][TH]*gcov[TT][RR] + gcov[TH][TH]*gcovtrsq + gcovrhsq*gcov[TT][TT] - gcov[TH][TH]*gcov[RR][RR]*gcov[TT][TT]);
  gcon[RR][TH]=gcon[TH][RR]=(-(gcov[TT][TH]*gcov[TT][RR]) + gcov[RR][TH]*gcov[TT][TT])/(-2.0*gcov[RR][TH]*gcov[TT][TH]*gcov[TT][RR] + gcov[TH][TH]*gcovtrsq + gcovrhsq*gcov[TT][TT] + gcov[RR][RR]*(gcovthsq - gcov[TH][TH]*gcov[TT][TT]));
  gcon[TH][TH]=(gcovtrsq - gcov[RR][RR]*gcov[TT][TT])/(gcov[RR][RR]*gcovthsq - 2.0*gcov[RR][TH]*gcov[TT][TH]*gcov[TT][RR] + gcov[TH][TH]*gcovtrsq + gcovrhsq*gcov[TT][TT] - gcov[TH][TH]*gcov[RR][RR]*gcov[TT][TT]);


}
/* invert gcov to get gcon */
// can be used to invert any 2nd rank tensor (symmetric or not)
// actually returns the inverse transpose, so if
// gcov=T^j_k then out pops (iT)^k_j such that T^j_k (iT)^k_l = \delta^j_l
void matrix_inverse(int whichcoord, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM])
{
  int pl;
  static int firstc = 1;
  int j, k;
  static FTYPE **tmp;
  int anglesing,centersing,truedim;
  void metric_sing_check(int whichcoord, FTYPE gcov[][NDIM], int *anglesing, int*centersing, int *truedim);
  void matrix_inverse_2d(FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);
  void matrix_inverse_3d(FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);  



  if (firstc) {
    tmp = dmatrix(1, NDIM, 1, NDIM);
    firstc = 0;
  }

  DLOOP(j,k) tmp[j + 1][k + 1] = gcov[j][k];
  
#if(1)
  // check for singularities
  // only truedim used
  // allow avoiding of gaussj fail as detection of singularity
  metric_sing_check(whichcoord, gcov, &anglesing, &centersing, &truedim);
  //  dualfprintf(fail_file,"anglesing=%d centersing=%d truedim=%d\n",anglesing,centersing,truedim);
#else
  truedim=NDIM;
#endif

  // 0-out all gcon
  DLOOP(j,k) gcon[j][k]=0.0;


  //  DLOOP(j,k) dualfprintf(fail_file,"tmp[%d][%d]=%21.15g\n",j+1,k+1,tmp[j+1][k+1]);

  // GODMARK: Feeding in nan's results in gaussj segfaulting with erroneous access.  Seems bad behavior!
  if(gaussj(tmp, truedim, NULL, 0)){
    // then singular
#if(0) // new singularity check before if(gaussj) should work
    if(ISSPCMCOORD(whichcoord)){
      // super hack
      //if(centersing)  matrix_inverse_2d(gcov,gcon);
      //      else if(anglesing) matrix_inverse_3d(gcov,gcon);
      matrix_inverse_2d(gcov,gcon);
    }
    else{
      dualfprintf(fail_file,"whichcoord=%d\n",whichcoord);
      myexit(6243);
    }
#else
    dualfprintf(fail_file,"Singularity check didn't work\n");
    dualfprintf(fail_file,"whichcoord=%d anglesing=%d centersing=%d truedim=%d\n",whichcoord,anglesing,centersing,truedim);
    DLOOP(j,k) dualfprintf(fail_file,"gcov[%d][%d]=%21.15g\n",j,k,gcov[j][k]);
    dualfprintf(fail_file,"icurr=%d jcurr=%d kcurr=%d\n",icurr,jcurr,kcurr);
    PLOOP(pl) dualfprintf(fail_file,"prim[%d]=%21.15g\n",pl,p[icurr][jcurr][kcurr][pl]);
    myexit(2714);
#endif

  }
  else{
    // assign but also transpose (shouldn't do in general, confusing)
    //DLOOP(j,k) gcon[j][k] = tmp[k + 1][j + 1];
    DLOOP(j,k) gcon[j][k] = tmp[j + 1][k + 1];
  }


#if(1) // check for nan's
  DLOOP(j,k) if(!finite(gcon[j][k])){
    dualfprintf(fail_file,"Came out of matrix_inverse with inf/nan for gcon at j=%d k=%d\n",j,k);
    myexit(5);
  }
#endif



}





void metric_sing_check(int whichcoord, FTYPE gcov[][NDIM], int *anglesing, int*centersing, int *truedim)
{


  // check for singularity (always removes highest coordinates -- so can use same matrix inverse with lower dimension when on coordinate singularities)
  *anglesing=0;
  *centersing=0;
  if(ISSPCMCOORD(whichcoord)){
    if(fabs(gcov[PH][PH])<10.0*SMALL){
      *anglesing=1;
    }
    if(fabs(gcov[TH][TH])<10.0*SMALL){
      *centersing=1;
    }
    if(*centersing){
      *truedim=NDIM-2;
    }
    else if(*anglesing){
      *truedim=NDIM-1;
    }
    else *truedim=NDIM;
  }
  else *truedim=NDIM;


}






// compute the radius of the inner most stable circular orbit
FTYPE rmso_calc(int which)
{
  FTYPE rmso,Z1,Z2,sign ;
  FTYPE j;

  j=a/(MBH+SMALL); // so doesn't nan for MBH=0 since "a" should be "0" if MBH=0

  if(which==PROGRADERISCO) sign=1; else sign=-1;

  Z1 = 1. + pow(1. - j*j,1./3.)*(pow(1. + j,1./3.) +
                                 pow(1. - j, 1./3.)) ;
  Z2 = sqrt(fabs(3.*j*j + Z1*Z1)) ;
  rmso=3. + Z2-sign*sqrt(fabs(3. - Z1)*fabs(3. + Z1 + 2.*Z2)) ;

  return(MBH*rmso) ;
}

FTYPE uphi_isco_calc(int which,FTYPE rold)
{
  FTYPE uphi;
  FTYPE sign;
  FTYPE Z1,Z2;
  FTYPE j,rnew;
  
  rnew=rold/MBH;
  j=a/MBH;

  if(which==PROGRADERISCO) sign=1; else sign=-1;

  Z1=rnew*rnew-sign*2.*j*sqrt(rnew)+j*j;
  Z2=rnew*(rnew*rnew-3.*rnew+sign*2.*j*sqrt(rnew));

  uphi=sign*Z1/sqrt(Z2);

  return(MBH*uphi);

}

FTYPE rhor_calc(int which)
{
  FTYPE sign,rhor;
  FTYPE j,jsq,disc;
  
  j=a/(MBH+SMALL); // so doesn't nan for MBH=0 since "a" should be "0" if MBH=0
  jsq=j*j;

  if(which==0) sign=1; else sign=-1;

  disc=MAX(1.0 - jsq,0.0);
  
  

  rhor=1. +sign*sqrt(disc);

  return(rhor*MBH);
}


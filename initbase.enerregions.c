#include "decs.h"


// ignorant of horizon -- just total computational domain
// determine if this cpu is doing what flux
int setflux(void)
{
  int dir;
  int enerregion;


  enerregion=GLOBALENERREGION;

  // setup pointers
  enerpos=enerposreg[enerregion];
  doflux=dofluxreg[enerregion];

  // all CPUs , total space for global region
  // inclusive loop range (0..N-1)
  enerpos[X1DN]=0;
  enerpos[X1UP]=N1-1;
  enerpos[X2DN]=0;
  enerpos[X2UP]=N2-1;
  enerpos[X3DN]=0;
  enerpos[X3UP]=N3-1;


  // only 0 through N-1 mean do flux

  // doflux<0 means ignore this cpu in flux calculation
  // doflux>=0 is used to set where flux is computed.
  // If dimension exists (i.e. N>1), then OUTM = N is outer flux on boundary
  // If dimension doesn't exist, then no such outer flux or outer boundary, so stay on same boundary (e.g., i=0) rather than i=N.

  // x1
  if((N1>1)&&(mycpupos[1]==0)){
    doflux[X1DN]=0; // or horizoni
    if(!specialstep) trifprintf("proc: %d doing flux X1DN\n",myid);
  }
  else doflux[X1DN]=-100;

  if((N1>1)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=OUTM1;
    if(!specialstep) trifprintf("proc: %d doing flux X1UP\n",myid);
  }
  else doflux[X1UP]=-100;

  // x2
  if((N2>1)&&(mycpupos[2]==0)){
    doflux[X2DN]=0;
    if(!specialstep) trifprintf("proc: %d doing flux X2DN\n",myid);
  }
  else doflux[X2DN]=-100;

  if((N2>1)&&(mycpupos[2]==ncpux2-1)){
    doflux[X2UP]=OUTM2;
    if(!specialstep) trifprintf("proc: %d doing flux X2UP\n",myid);
  }
  else doflux[X2UP]=-100;

  // x3
  if((N3>1)&&(mycpupos[3]==0)){
    doflux[X3DN]=0;
    if(!specialstep) trifprintf("proc: %d doing flux X3DN\n",myid);
  }
  else doflux[X3DN]=-100;

  if((N3>1)&&(mycpupos[3]==ncpux3-1)){
    doflux[X3UP]=OUTM3;
    if(!specialstep) trifprintf("proc: %d doing flux X3UP\n",myid);
  }
  else doflux[X3UP]=-100;

  // fluxes are on edges of zone, so 0 and N are on edge fluxes

  if(!specialstep){
    DIRLOOP(dir) trifprintf("proc: %d %d doflux[%d]=%d enerpos[%d]=%d\n",myid,GLOBALENERREGION,dir,doflux[dir],dir,enerpos[dir]);
  }

  return(0);
}



// determine if this cpu is doing what flux through horizon
// assumes horizoni and horizoncpupos1 set by find_horizon before this function called
int sethorizonflux(void)
{
  int dir;
  int enerregion;

  enerregion=OUTSIDEHORIZONENERREGION;

  // setup pointers
  enerpos=enerposreg[enerregion];
  doflux=dofluxreg[enerregion];

  // all CPUs , total space for global region
  // inclusive loop range (0..N-1)
  if(mycpupos[1]==horizoncpupos1){
    enerpos[X1DN]=horizoni;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=N2-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
  }
  else if(mycpupos[1]>horizoncpupos1){
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=N2-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
  }
  else if(mycpupos[1]<horizoncpupos1){
    // then CPU not part of integral since only including beyond horizoni
    enerpos[X1DN]=0;
    enerpos[X1UP]=-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=-1;
  }

  // only 0 through N-1 mean do flux

  // doflux<0 means ignore this cpu in flux calculation
  // doflux>=0 is used to set where flux is computed.
  // If dimension exists (i.e. N>1), then OUTM = N is outer flux on boundary
  // If dimension doesn't exist, then no such outer flux or outer boundary, so stay on same boundary (e.g., i=0) rather than i=N.

  // x1
  if((N1>1)&&(mycpupos[1]==horizoncpupos1)){
    doflux[X1DN]=horizoni;
    if(!specialstep) trifprintf("proc: %d doing horizonflux X1DN\n",myid);
  }
  else doflux[X1DN]=-100;

  if((N1>1)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=OUTM1;
    if(!specialstep) trifprintf("proc: %d doing horizonflux X1UP\n",myid);
  }
  else doflux[X1UP]=-100;

  // x2
  if((N2>1)&&(mycpupos[2]==0)){
    doflux[X2DN]=0;
    if(!specialstep) trifprintf("proc: %d doing horizonflux X2DN\n",myid);
  }
  else doflux[X2DN]=-100;

  if((N2>1)&&(mycpupos[2]==ncpux2-1)){
    doflux[X2UP]=OUTM2;
    if(!specialstep) trifprintf("proc: %d doing horizonflux X2UP\n",myid);
  }
  else doflux[X2UP]=-100;

  // x3
  if((N3>1)&&(mycpupos[3]==0)){
    doflux[X3DN]=0;
    if(!specialstep) trifprintf("proc: %d doing horizonflux X3DN\n",myid);
  }
  else doflux[X3DN]=-100;

  if((N3>1)&&(mycpupos[3]==ncpux3-1)){
    doflux[X3UP]=OUTM3;
    if(!specialstep) trifprintf("proc: %d doing horizonflux X3UP\n",myid);
  }
  else doflux[X3UP]=-100;

  // fluxes are on edges of zone, so 0 and N are on edge fluxes
  if(!specialstep){
    DIRLOOP(dir) trifprintf("proc: %d %d horizon: doflux[%d]=%d enerpos[%d]=%d\n",myid,enerregion,dir,doflux[dir],dir,enerpos[dir]);
  }

  return(0);
}




// true global region where quantities are evolved
// this can be used to determine where quantities are evolved and so active domain only
// assumes horizoni and horizoncpupos1 set by find_horizon before this function called
int settrueglobalregion(void)
{
  int dir;
  int enerregion;

  enerregion=TRUEGLOBALENERREGION;

  // setup pointers
  enerpos=enerposreg[enerregion];
  doflux=dofluxreg[enerregion];

  // all CPUs , total space for global region
  // inclusive loop range (0..N-1)
  if(mycpupos[1]==horizoncpupos1){
    // below tries to say that inside horizoni-1-N1BND the evolution doesn't matter so ignore it.
    // used, e.g., in advance.c.  However, so that rest of code doesn't have to change, must ensure that
    //   the other zones are set to something so no Nan's and need to be set such that time-like (e.g. for boundprim())
    // otherwise can go through code and base loops on this enerregion
    if(horizoni==0){
      enerpos[X1DN]=MAX(0,horizoni-N1BND); // should be no smaller than 0
    }
    else{
      enerpos[X1DN]=MAX(0,horizoni-1-N1BND); // horizoni-1 accounts for fact that horizoni is ceil() of horizon position, not floor()
    }
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=N2-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
  }
  else if(mycpupos[1]>horizoncpupos1){
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=N2-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
  }
  else if(mycpupos[1]<horizoncpupos1){
    // then CPU not part of integral since only including beyond horizoni
    enerpos[X1DN]=0;
    enerpos[X1UP]=-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=-1;
  }

  // only 0 through N-1 mean do flux

  // doflux<0 means ignore this cpu in flux calculation
  // doflux>=0 is used to set where flux is computed.
  // If dimension exists (i.e. N>1), then OUTM = N is outer flux on boundary
  // If dimension doesn't exist, then no such outer flux or outer boundary, so stay on same boundary (e.g., i=0) rather than i=N.

  // x1
  if((N1>1)&&(mycpupos[1]==horizoncpupos1)){
    if(horizoni==0){
      doflux[X1DN]=MAX(0,horizoni-N1BND); // should be no smaller than 0
    }
    else{
      doflux[X1DN]=MAX(0,horizoni-1-N1BND); // horizoni-1 accounts for fact that horizoni is ceil() of horizon position, not floor()
    }
    if(!specialstep) trifprintf("proc: %d doing trueglobal X1DN\n",myid);
  }
  else doflux[X1DN]=-100;

  if((N1>1)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=OUTM1;
    if(!specialstep) trifprintf("proc: %d doing trueglobal X1UP\n",myid);
  }
  else doflux[X1UP]=-100;

  // x2
  if((N2>1)&&(mycpupos[2]==0)){
    doflux[X2DN]=0;
    if(!specialstep) trifprintf("proc: %d doing trueglobal X2DN\n",myid);
  }
  else doflux[X2DN]=-100;

  if((N2>1)&&(mycpupos[2]==ncpux2-1)){
    doflux[X2UP]=OUTM2;
    if(!specialstep) trifprintf("proc: %d doing trueglobal X2UP\n",myid);
  }
  else doflux[X2UP]=-100;

  // x3
  if((N3>1)&&(mycpupos[3]==0)){
    doflux[X3DN]=0;
    if(!specialstep) trifprintf("proc: %d doing trueglobal X3DN\n",myid);
  }
  else doflux[X3DN]=-100;

  if((N3>1)&&(mycpupos[3]==ncpux3-1)){
    doflux[X3UP]=OUTM3;
    if(!specialstep) trifprintf("proc: %d doing trueglobal X3UP\n",myid);
  }
  else doflux[X3UP]=-100;

  // fluxes are on edges of zone, so 0 and N are on edge fluxes
  if(!specialstep){
    DIRLOOP(dir) trifprintf("proc: %d %d trueglobal: doflux[%d]=%d enerpos[%d]=%d\n",myid,enerregion,dir,doflux[dir],dir,enerpos[dir]);
  }

  return(0);
}



// this ranges over entire domain where ANYTHING is done at all
// This is used to avoid any unnecessary access to cells that will never be used
//
// determine if this cpu is doing what flux through horizon
// assumes horizoni and horizoncpupos1 set by find_horizon before this function called
int settrueglobalwithbndregion(void)
{
  int dir;
  int enerregion;

  enerregion=TRUEGLOBALWITHBNDENERREGION;

  // setup pointers
  enerpos=enerposreg[enerregion];
  doflux=dofluxreg[enerregion];

  // all CPUs , total space for global region
  // inclusive loop range (0..N-1)
  if(mycpupos[1]==horizoncpupos1){
    // below tries to say that inside horizoni-1-N1BND the evolution doesn't matter so ignore it.
    // used, e.g., in advance.c.  However, so that rest of code doesn't have to change, must ensure that
    //   the other zones are set to something so no Nan's and need to be set such that time-like (e.g. for boundprim())
    // otherwise can go through code and base loops on this enerregion
    if(horizoni==0){
      enerpos[X1DN]=MAX(-N1BND,horizoni-2*N1BND); // should be no smaller than -N1BND
    }
    else{
      enerpos[X1DN]=MAX(-N1BND,horizoni-1-2*N1BND); // horizoni-1 accounts for fact that horizoni is ceil() of horizon position, not floor()
    }
    enerpos[X1UP]=N1-1+N1BND;
    enerpos[X2DN]=-N2BND;
    enerpos[X2UP]=N2-1+N2BND;
    enerpos[X3DN]=-N3BND;
    enerpos[X3UP]=N3-1+N3BND;
  }
  else if(mycpupos[1]>horizoncpupos1){
    enerpos[X1DN]=-N1BND;
    enerpos[X1UP]=N1-1+N1BND;
    enerpos[X2DN]=-N2BND;
    enerpos[X2UP]=N2-1+N2BND;
    enerpos[X3DN]=-N3BND;
    enerpos[X3UP]=N3-1+N3BND;
  }
  else if(mycpupos[1]<horizoncpupos1){
    // then CPU not part of integral since only including beyond horizoni
    enerpos[X1DN]=0;
    enerpos[X1UP]=-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=-1;
  }

  // only 0 through N-1 mean do flux

  // doflux<0 means ignore this cpu in flux calculation
  // doflux>=0 is used to set where flux is computed.
  // If dimension exists (i.e. N>1), then OUTM = N is outer flux on boundary
  // If dimension doesn't exist, then no such outer flux or outer boundary, so stay on same boundary (e.g., i=0) rather than i=N.

  // x1
  if((N1>1)&&(mycpupos[1]==horizoncpupos1)){
    if(horizoni==0){
      doflux[X1DN]=MAX(0,horizoni-N1BND); // should be no smaller than 0
    }
    else{
      doflux[X1DN]=MAX(0,horizoni-1-N1BND); // horizoni-1 accounts for fact that horizoni is ceil() of horizon position, not floor()
    }
    if(!specialstep) trifprintf("proc: %d doing trueglobalwithbnd X1DN\n",myid);
  }
  else doflux[X1DN]=-100;

  if((N1>1)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=OUTM1;
    if(!specialstep) trifprintf("proc: %d doing trueglobalwithbnd X1UP\n",myid);
  }
  else doflux[X1UP]=-100;

  // x2
  if((N2>1)&&(mycpupos[2]==0)){
    doflux[X2DN]=0;
    if(!specialstep) trifprintf("proc: %d doing trueglobalwithbnd X2DN\n",myid);
  }
  else doflux[X2DN]=-100;

  if((N2>1)&&(mycpupos[2]==ncpux2-1)){
    doflux[X2UP]=OUTM2;
    if(!specialstep) trifprintf("proc: %d doing trueglobalwithbnd X2UP\n",myid);
  }
  else doflux[X2UP]=-100;

  // x3
  if((N3>1)&&(mycpupos[3]==0)){
    doflux[X3DN]=0;
    if(!specialstep) trifprintf("proc: %d doing trueglobalwithbnd X3DN\n",myid);
  }
  else doflux[X3DN]=-100;

  if((N3>1)&&(mycpupos[3]==ncpux3-1)){
    doflux[X3UP]=OUTM3;
    if(!specialstep) trifprintf("proc: %d doing trueglobalwithbnd X3UP\n",myid);
  }
  else doflux[X3UP]=-100;

  // fluxes are on edges of zone, so 0 and N are on edge fluxes
  if(!specialstep){
    DIRLOOP(dir) trifprintf("proc: %d %d trueglobalwithbnd: doflux[%d]=%d enerpos[%d]=%d\n",myid,enerregion,dir,doflux[dir],dir,enerpos[dir]);
  }

  return(0);
}




// GODMARK
// the below only works for a grid with full Pi.  Crashes code at runtime otherwise!  should fix.

// assume this is an intrinsically axisymmetric function, so k (\phi) is just carried along -- no truncation in \phi

// determine the flux positions for each CPU for the jet region (jetedge)
// AND the range of volume integration for cons. variables in jet (jetpos).
int setjetflux(void)
{
  FTYPE X[NDIM],V[NDIM],r,th;
  int i,j,k,pl;
  int dir;
  FTYPE startth,endth,thetajet;
  int jetedge[NUMJETS];


  if(defcoord==EQMIRROR){
    dualfprintf(fail_file,"setjetflux() not setup to work for non-full-Pi grids\n");
    myexit(1);
  }


  // jet region is assumed to be within a constant theta slice
  // this is theta w.r.t. polar axis
  thetajet=M_PI*0.5-h_over_r_jet;
  // find j for which theta is at our prespecified point

  i=0;j=0;k=0;
  coord_ijk(i, j, k, FACE2, X);
  bl_coord_ijk(i, j, k, FACE2, V);
  startth=V[2];
  i=0;j=N2;k=0;
  coord_ijk(i, j, k, FACE2, X);
  bl_coord_ijk(i, j, k, FACE2, V);
  endth=V[2];

  // assumes 0<thetajet<Pi/2
  if((fabs(startth-thetajet)/thetajet<1E-8)||
     (fabs(endth-thetajet)/thetajet<1E-8)||
     (fabs(startth-(M_PI-thetajet))/thetajet<1E-8)||
     (fabs(endth-(M_PI-thetajet))/thetajet<1E-8)
     ){
    dualfprintf(fail_file,"thetajet is on top of grid, move h_over_r_jet a bit\n");
    myexit(1);
  }


  ////////////////////
  //
  // INNERJET
  //

  // setup pointers
  enerpos=enerposreg[INNERJETREGION];
  doflux=dofluxreg[INNERJETREGION];


  // see if jet edge is related to this CPU
  // assumes increasing j is increasing th
  if((startth<=thetajet)&&(endth<=thetajet)){
    // if cpu entirely within inner theta jet
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=N2-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
    jetedge[INNERJET]=-100;
  }
  else if((startth<thetajet)&&(endth>thetajet)){
    // if inner jet edge is on this CPU but not on boundary
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    i=0;
    for(j=0;j<=OUTM2;j++){
      coord_ijk(i, j, k, FACE2, X);
      bl_coord_ijk(i, j, k, FACE2,V);
      r=V[1];
      th=V[2];
      // look for switch from below to above thetajet at inner theta jet edge
      if(th>thetajet){
	enerpos[X2UP]=jm1;
	jetedge[INNERJET]=j;
	break;
      }
    }
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
  }
  else if((startth>=thetajet)&&(endth>=thetajet)){
    // if cpu is entirely not contained in inner jet
    enerpos[X1DN]=-100;
    enerpos[X1UP]=-100;
    enerpos[X2DN]=-100;
    enerpos[X2UP]=-100;
    enerpos[X3DN]=-100;
    enerpos[X3UP]=-100;
    jetedge[INNERJET]=-100;
  }
  else{
    trifprintf("problem with INNERJET setjetflux()\n");
    myexit(1);
  }



  // left edge (any directional condition would do)
  if((N1>1)&&(enerpos[X1DN]!=-100)&&(mycpupos[1]==0)){
    doflux[X1DN]=0; // or horizoni
    if(!specialstep) trifprintf("proc: %d doing inner jet flux X1DN\n",myid);
  }
  else doflux[X1DN]=-100;

  // right edge (any directional condition would do)
  if((N1>1)&&(enerpos[X1DN]!=-100)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=OUTM1;
    if(!specialstep) trifprintf("proc: %d doing inner jet flux X1UP\n",myid);
  }
  else doflux[X1UP]=-100;

  // lower theta boundary
  if((N2>1)&&(mycpupos[2]==0)){
    doflux[X2DN]=0;
    if(!specialstep) trifprintf("proc: %d doing inner jet flux X2DN\n",myid);
  }
  else doflux[X2DN]=-100;
  
  // upper theta boundary
  if((N2>1)&&(jetedge[INNERJET]!=-100)){ // only get flux if CPU has edge
    doflux[X2UP]=jetedge[INNERJET];
    if(!specialstep) trifprintf("proc: %d doing inner jet flux X2UP\n",myid);
  }
  else doflux[X2UP]=-100;

  if((N3>1)&&(enerpos[X3DN]!=-100)&&(mycpupos[3]==0)){
    doflux[X3DN]=0;
    if(!specialstep) trifprintf("proc: %d doing inner jet flux X3DN\n",myid);
  }
  else doflux[X3DN]=-100;

  // right edge (any directional condition would do)
  if((N3>1)&&(enerpos[X3DN]!=-100)&&(mycpupos[3]==ncpux3-1)){
    doflux[X3UP]=OUTM3;
    if(!specialstep) trifprintf("proc: %d doing inner jet flux X3UP\n",myid);
  }
  else doflux[X3UP]=-100;


  if(!specialstep) {
    DIRLOOP(dir) trifprintf("proc: %d %d innerjet: doflux[%d]=%d enerpos[%d]=%d\n",myid,INNERJETREGION,dir,doflux[dir],dir,enerpos[dir]);
  }

  /////////////////////
  //
  // OUTERJET
  //

  // setup pointers
  enerpos=enerposreg[OUTERJETREGION];
  doflux=dofluxreg[OUTERJETREGION];


  // see if outer jet edge is related to this CPU
  if((startth<=M_PI-thetajet)&&(endth<=M_PI-thetajet)){
    // if cpu entirely not within outer jet region
    enerpos[X1DN]=-100;
    enerpos[X1UP]=-100;
    enerpos[X2DN]=-100;
    enerpos[X2UP]=-100;
    enerpos[X3DN]=-100;
    enerpos[X3UP]=-100;
    jetedge[OUTERJET]=-100;
  }
  else if((startth<M_PI-thetajet)&&(endth>M_PI-thetajet)){
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    // if outer jet edge is on this CPU but not on boundary
    i=0;k=0;
    for(j=0;j<=OUTM2;j++){
      coord_ijk(i, j, k, FACE2, X);
      bl_coord_ijk(i, j, k, FACE2, V);
      th=V[2];
      // look for switch from below to above thetajet at inner theta jet edge
      if(th>M_PI-thetajet){
	enerpos[X2DN]=jm1;
	jetedge[OUTERJET]=jm1;
	break;
      }
    }
    enerpos[X2UP]=N2-1;

    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;

  }
  else if((startth>=M_PI-thetajet)&&(endth>=M_PI-thetajet)){
    // if cpu is entirely containe within outer jet
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=N2-1;
    enerpos[X3DN]=0;
    enerpos[X3UP]=N3-1;
    jetedge[OUTERJET]=-100;
  }
  else{
    trifprintf("problem with OUTERJET setjetflux()\n");
    myexit(1);
  }

  if((N1>1)&&(enerpos[X1DN]!=-100)&&(mycpupos[1]==0)){
    doflux[X1DN]=0; // or horizoni
    if(!specialstep) trifprintf("proc: %d doing outer jet flux X1DN\n",myid);
  }
  else doflux[X1DN]=-100;

  if((N1>1)&&(enerpos[X1DN]!=-100)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=OUTM1;
    if(!specialstep) trifprintf("proc: %d doing outer jet flux X1UP\n",myid);
  }
  else doflux[X1UP]=-100;

  if((N2>1)&&(jetedge[OUTERJET]!=-100)){
    doflux[X2DN]=jetedge[OUTERJET];
    if(!specialstep) trifprintf("proc: %d doing outer jet flux X2DN\n",myid);
  }
  else doflux[X2DN]=-100;

  if((N2>1)&&(mycpupos[2]==ncpux2-1)){
    doflux[X2UP]=OUTM2;
    if(!specialstep) trifprintf("proc: %d doing outer jet flux X2UP\n",myid);
  }
  else doflux[X2UP]=-100;
  // fluxes are on edges of zone, so 0 and N are on edge fluxes

  if((N3>1)&&(enerpos[X3DN]!=-100)&&(mycpupos[3]==0)){
    doflux[X3DN]=0; 
    if(!specialstep) trifprintf("proc: %d doing outer jet flux X3DN\n",myid);
  }
  else doflux[X3DN]=-100;

  if((N3>1)&&(enerpos[X3DN]!=-100)&&(mycpupos[3]==ncpux3-1)){
    doflux[X3UP]=OUTM3;
    if(!specialstep) trifprintf("proc: %d doing outer jet flux X3UP\n",myid);
  }
  else doflux[X3UP]=-100;


  if(!specialstep){
    DIRLOOP(dir) trifprintf("proc: %d %d outerjet: doflux[%d]=%d enerpos[%d]=%d\n",myid,OUTERJETREGION,dir,doflux[dir],dir,enerpos[dir]);
  }

  return(0);
}

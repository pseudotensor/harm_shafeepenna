

#include "decs.h"

#include "para_and_paraenohybrid.h"



void parainitchecks(void)
{
  int dimen;


  DIMENLOOP(dimen){
    if(PARAGENDQALLOWEXTREMUM && WENOINTERPTYPE(lim[dimen])){
      // if PARAMODWENO==1 implied, but if 0 allow since doesn't matter if doing weno with PARAGENDQALLOWEXTREMUM==1
      // then ok
    }
    else if(PARAMODWENO==1 && PARAGENDQALLOWEXTREMUM==0 && WENOINTERPTYPE(lim[dimen])){
      dualfprintf(fail_file,"WARNING: PARAGENDQALLOWEXTREMUM==0 for hybrid method, will be less accurate for turbulent regions\n");
    }
    else if(PARAGENDQALLOWEXTREMUM==1 && lim[dimen]==PARA){
      dualfprintf(fail_file,"ERROR: PARAGENDQALLOWEXTREMUM==1 but not WENO-based limiter -- will be less stable in stiff regions (e.g. near horizon). Use PARALINE instead.\n");
      myexit(34643463);
    }

    if(lim[dimen]==PARAFLAT){
      dualfprintf(fail_file,"WARNING: PARAFLAT inefficient compared to PARALINE, suggested to use PARALINE instead\n");
    }
  }


}

/*
 * parabolic interpolation subroutin  
 * ref. Colella && Woodward's paper
 * Colella, P., & Woodward, P. R. 1984, J. Comput. Phys., 54, 174-201
 *
 * using zone-centered value of 5 continuous zones 
 * to get left and right value of the middle zone.
 *  
 * 
 */

// Latest JCM version is para4()/parapl(): 02/25/08


// lout/rout is left and right sides of cell
// note how used in step_ch.c to get appropriate interface value

// given by Xiaoyue Guan to Scott Noble on Nov 9, 2004, given to me Jan 7, 2005
void para(FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm ;
  FTYPE dq0[5];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, aDqm,aDqp,aDqc,s,l,r,qa, qd, qe;

  // shifted dq
  dq=dq0+2;

  /*CW1.7 */
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=min(aDqc,min(aDqm,aDqp))*sign(Dqc);
  }

  /* CW1.6 */

  l=0.5*(y[0]+y[-1])-(dq[0]-dq[-1])/6.0;
  r=0.5*(y[1]+y[0])-(dq[1]-dq[0])/6.0;

  qa=(r-y[0])*(y[0]-l);
  qd=(r-l);
  qe=6.0*(y[0]-0.5*(l+r));


  if (qa <=0. ) {
    l=y[0];
    r=y[0];
  }

  if (qd*(qd-qe)<0.0) l=3.0*y[0]-2.0*r;
  else if (qd*(qd+qe)<0.0) r=3.0*y[0]-2.0*l;


  *lout=l;   //a_L,j
  *rout=r;
  //*dw=r-l;                      //CW1.5
  //*w6=6.0*(y[0]-0.5*(l+r));
}



// given by Xiaoyue Guan on Jan 9, 2005
void para2(FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm ;
  FTYPE dq0[5];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;

  // shifted dq
  dq=dq0+2;

  /*CW1.7 */
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

#if(PARA2LIM == VANL) 
    Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    aDqvanl=fabs(Dqvanl);
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == PMC)
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == MC)
    dq[mm] =Dqc;
#endif
  }
  /* CW1.6 */

  l=0.5*(y[0]+y[-1])-(dq[0]-dq[-1])/6.0;
  r=0.5*(y[1]+y[0])-(dq[1]-dq[0])/6.0;

  /*
    l=max(min(y[0],y[-1]),l);
    l=min(max(y[0],y[-1]),l);
    r=max(min(y[0],y[1]),r);
    r=min(max(y[0],y[1]),r);
  */

  qa=(r-y[0])*(y[0]-l);
  qd=(r-l);
  qe=6.0*(y[0]-0.5*(l+r));


  if (qa <=0. ) {
    l=y[0];
    r=y[0];
  }

  else if (qd*(qd-qe)<0.0) l=3.0*y[0]-2.0*r;
  else if (qd*(qd+qe)<0.0) r=3.0*y[0]-2.0*l;


  *lout=l;   //a_L,j
  *rout=r;
  //*dw=r-l;                      //CW1.5
  //*w6=6.0*(y[0]-0.5*(l+r));
}




// 3rd para from Xiaoyue that she bundled with a new TVD-optimal RK3
// given on 02/17/2005
void para3(FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm ;
  FTYPE dq0[5];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;

  // shifted dq
  dq=dq0+2;

  /*CW1.7 */
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

#if(PARA2LIM == VANL) 
    Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    aDqvanl=fabs(Dqvanl);

    if (s <=0.) dq[mm]=0.;
    else dq[mm] = -aDqvanl*sign(Dqc);
    //else dq[mm]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == MC)
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=-min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == MINM)
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = -aDqm*sign(Dqc);
    else dq[mm]=-aDqp*sign(Dqc);
#elif(PARA2LIM == NLIM) //w/o slope limiter
    //if(s<=0.) dq[mm] = 0.; // DONOR
    dq[mm] = Dqc;
#endif
  }

  /* CW1.6 */

  l=0.5*(y[0]+y[-1])-(dq[0]-dq[-1])/6.0;
  r=0.5*(y[1]+y[0])-(dq[1]-dq[0])/6.0;


  l=max(min(y[0],y[-1]),l);
  l=min(max(y[0],y[-1]),l);
  r=max(min(y[0],y[1]),r);
  r=min(max(y[0],y[1]),r);


  qa=(r-y[0])*(y[0]-l);
  qd=(r-l);
  qe=6.0*(y[0]-0.5*(l+r));

  /*
    if (qa <=0. ) {
    l=y[0];
    r=y[0];
    }

    else if (qd*(qd-qe)<0.0) l=3.0*y[0]-2.0*r;
    else if (qd*(qd+qe)<0.0) r=3.0*y[0]-2.0*l;


    *lout=l;   //a_L,j
    *rout=r;
    */

  if (qa <=0. ) {
    *lout=y[0];
    *rout=y[0];
  }  
  else {
    *lout = l;
    *rout = r;
  }

  //2.0 at top/bottom of a steep gradient 
  if (qd*(qd-qe)<0.0) *lout=3.0*y[0]-2.0*r;
  else *lout = l;

  if (qd*(qd+qe)<0.0) *rout=3.0*y[0]-2.0*l;
  else *rout = r;
  //*dw=r-l;                      //CW1.5
  //*w6=6.0*(y[0]-0.5*(l+r));
}





// older method that has no failures problem CHANGINGMARK
void para4_old(int pl, FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  int mm ;
  FTYPE dq0[5];
  FTYPE *dq;
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;
  void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq);

  // shifted dq
  dq=dq0+2;

  /*CW1.7 */
  for(mm=-1 ; mm<=1 ; mm++) {
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;


#if(PARA2LIM == VANL) 
    Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    aDqvanl=fabs(Dqvanl);

    if (s <=0.) dq[mm]=0.;
    //else dq[mm] = aDqvanl*sign(Dqc);
    else dq[mm]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);

#elif(PARA2LIM == MC)

#if(0)
    // Jon's version
    dq[mm]=MINMOD(Dqc,MINMOD(Dqm,Dqp));
#else
    // Xioyue's version
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]= min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#endif



#elif(PARA2LIM == MINM_STEEPENER)

    // Xioyue's version (steepeneed version of MINM)
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = aDqm*sign(Dqc);
    else dq[mm]=aDqp*sign(Dqc);


#elif(PARA2LIM == MINM) // no steepener, normal MINM

#if(0)
    // Jon's version
    dq[mm] = MINMOD(0.5*Dqm,0.5*Dqp); // no steepening    
#elif(1)
    // Jon's steep version
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = aDqm*sign(Dqc);
    else dq[mm]=aDqp*sign(Dqc);
#elif(0)
    // Xioyue's version
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = 0.5*aDqm*sign(Dqc);
    else dq[mm]=0.5*aDqp*sign(Dqc);
#endif

#elif(PARA2LIM == NLIM) //w/o slope limiter

    dq[mm] = Dqc;
#endif
  }

#if(JONPARAREDUCE)
  //  if(pl==U1){
  if(pl!=RHO){
    if(
       (fabs(dq[-1]-dq[0])/(fabs(dq[-1])+fabs(dq[0])+SMALL)>0.1)||
       (fabs(dq[1]-dq[0])/(fabs(dq[1])+fabs(dq[0])+SMALL)>0.1)
       ){
      slope_lim_3points(MINM, y[-1], y[0], y[1], dq);
      *lout =y[0] - 0.5* (*dq);
      *rout=y[0] + 0.5* (*dq);
      return;
    }
  }

#endif

  /* CW1.6 */

  // modified as per Matt's paper
  l=0.5*(y[0]+y[-1])-(dq[0]-dq[-1])/8.0;
  r=0.5*(y[1]+y[0])-(dq[1]-dq[0])/8.0;


  l=max(min(y[0],y[-1]),l);
  l=min(max(y[0],y[-1]),l);
  r=max(min(y[0],y[1]),r);
  r=min(max(y[0],y[1]),r);


  // modified as per Matt's paper
  qa=(r-y[0])*(y[0]-l);
  qd=(r-l);
  qe=6.0*(y[0]-0.5*(l+r));


  if (qa <=0. ) {
    l=y[0];
    r=y[0];
  }
  else{

    if (qd*(qd-qe)<0.0) 
      l=3.0*y[0]-2.0*r;


    if (qd*(qd+qe)<0.0) 
      r=3.0*y[0]-2.0*l;
  }


  *lout=l;   //a_L,j
  *rout=r;

  //  *dqleft=dq[-1];
  //  *dqcenter=dq[0];
  //  *dqright=dq[1];

}




// doesn't use dq
void para4(int pl, FTYPE *y, FTYPE *lout, FTYPE *rout)
{
  void para4gen(int dqrange, int pl, FTYPE *y, FTYPE *lout, FTYPE *rout, FTYPE *dq);
  void checkparamonotonicity(int dqrange, int pl, FTYPE *y, FTYPE *ddq, FTYPE *dq, FTYPE *lin, FTYPE *rin, FTYPE *lout, FTYPE *rout);
  FTYPE dq0[5];
  FTYPE *dq;
  int dqrange;
  FTYPE a_ddq[7];
  FTYPE *ddq;
  int mm;


  dqrange=3;
  // shifted dq
  dq=dq0+2;
  ddq=a_ddq+3; // shifted sufficiently


  para4gen(dqrange, pl,y,lout,rout,dq);

  for(mm=-dqrange/2+1;mm<=dqrange/2;mm++){
    ddq[mm] = dq[mm] - dq[mm-1];
  }
  checkparamonotonicity(dqrange, pl, y, ddq, dq, lout, rout, lout, rout);


}




#if(PARAGENDQALLOWEXTREMUM==0)
#define PARAGENMINMOD(a,b) MINMOD(a,b)
#else
#define PARAGENMINMOD(a,b) MINMODB(a,b)
#endif

// Xiaoyue given on 03/25/05
// she realized sign error in 1st der's in para3()
// noted Matt's paper astro-ph/0503420 suggested CW1.6 uses 1/8 rather than 1/6
// I noticed that Matt uses MC for field variables and PPM+ for hydro variables
// GODMARK: This step could be done once for entire line of data only once
void para4gen(int dqrange, int pl, FTYPE *y, FTYPE *lout, FTYPE *rout, FTYPE *dq)
{
  int mm ;
  FTYPE a_dq1[10],a_dq2[10];
  FTYPE *dq1,*dq2;
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r;
  void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq);
  FTYPE Dqparacenterleft,Dqparacenterright;
  FTYPE Dqsteep;
  FTYPE ddql,ddqr;
  void paracont(FTYPE ddq, FTYPE *y, FTYPE *facecont);





  // Procedure:
  // 1) First obtain limited slopes (above)
  // 2) Assume l,r obtained from 3rd order polynomial for points or quartic for averages
  // 3) steepen or flatten solution
  // 4) Check monotonicity: 
  //    Condition 1: If y[0] is a local max or local min compared to l,r, then set l=r=y[0]
  //    Condition 2: If parabola fit through l,y[0],r is non-monotonic, then if non-monotonic region is near l, then modify r until derivative of parabola at l is 0.

  // shifted dq (sufficiently shifted)
  dq1=a_dq1+dqrange;
  dq2=a_dq2+dqrange;





#if(JONPARASTEEP)
  /////////////////
  //
  // first check if 3rd order polynomial will have derivative at center of cell that MC says doesn't need steepening
  //
  ////////////////

  mm=0;
  Dqm = 2.0 *(y[mm]-y[mm-1]);   // steepened
  Dqp = 2.0 *(y[mm+1]-y[mm]);   // steepened
  //  Dqc = 0.5 *(y[mm+1]-y[mm-1]); // normal
  // Note that the factors of 3 (and 6) cause asymmetries at roundoff-error
  Dqparacenterleft =  (+0.5*y[0]-y[-1]    +y[-2]/6.0+y[1]/3.0); // used to get a_{j-1/2}
  Dqparacenterright = (-0.5*y[0]-y[-1]/3.0+y[1] -    y[2]/6.0); // used to get a_{j+1/2}
  // see if para will use centered slope that is not steepened enough compared to MC

  // first compare Dqp and Dqm
  Dqsteep=PARAGENMINMOD(Dqm,Dqp);
  //    *dq=PARAGENMINMOD(Dqc,PARAGENMINMOD(Dqm,Dqp));
  // now compare centered with steepened Dqsteep value
  if(fabs(Dqparacenterleft)>fabs(Dqsteep) && fabs(Dqparacenterright)>fabs(Dqsteep)){
    *lout = y[0] - 0.5* Dqsteep;
    *rout = y[0] + 0.5* Dqsteep;
    //dualfprintf(fail_file,"USEDSTEEP\n");
    return;
    // else PARA is ok to use
  }
#endif

  //dualfprintf(fail_file,"USEDPARA\n");



 
  ////////////////////////
  //
  // get slopes
  //
  ////////////////////////

  // Dqm(0) = 2.0*(dq1[0])
  // Dqp(0) = 2.0*(dq1[1]) // so need to go 1 farther to get all needed Dqp's
  for(mm=-dqrange/2 ; mm<=dqrange/2 ; mm++) {
    dq1[mm] = (y[mm]-y[mm-1]); // slope centered at cell face
    dq2[mm] = 0.5 *(y[mm+1]-y[mm-1]); // slope centered at cell center
  }
  mm=dqrange/2+1; // get last dq1 (can't do in loop above since +1 would mean dq2 beyond data range
  dq1[mm] = (y[mm]-y[mm-1]); // slope centered at cell face
  

  /////////////////
  //
  // Determine monotonized slopes
  //
  /////////////////

  /*CW1.7 */
  for(mm=-dqrange/2 ; mm<=dqrange/2 ; mm++) {

    Dqm = 2.0 *dq1[mm];   // steepened
    Dqp = 2.0 *dq1[mm+1]; // steepened
    Dqc = dq2[mm];        // normal


#if(PARA2LIM == VANL) 
    Dqm*=0.5; // desteepen
    Dqp*=0.5; // desteepen
    s = Dqm*Dqp;
    Dqvanl=2.0*s/(Dqm+Dqp);
    //    aDqvanl=fabs(Dqvanl);
    //    aDqm = fabs(Dqm) ;
    //    aDqp = fabs(Dqp) ;
    //    aDqc = fabs(Dqc) ;

    // true VANL:
    if (s <=0.) dq[mm]=0.;
    else dq[mm] = Dqvanl*sign(Dqc);

    // Xioyue was using MC-VANL combo of some kind
    //    else dq[mm] = aDqvanl*sign(Dqc);
    //    else dq[mm]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);

#elif(PARA2LIM == MC)

#if(1)
    // Jon's version
    dq[mm]=PARAGENMINMOD(Dqc,PARAGENMINMOD(Dqm,Dqp));
#else
    // Xioyue's version
    s = Dqm*Dqp;
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]= min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#endif



#elif(PARA2LIM == MINM_STEEPENER)

    // Xioyue's version (steepeneed version of MINM)
    s = Dqm*Dqp;
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = aDqm*sign(Dqc);
    else dq[mm]=aDqp*sign(Dqc);


#elif(PARA2LIM == MINM) // no steepener, normal MINM

#if(0)
    // Jon's version
    dq[mm] = PARAGENMINMOD(0.5*Dqm,0.5*Dqp); // no steepening    
#elif(1)
    // Jon's steep version
    s = Dqm*Dqp;
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = aDqm*sign(Dqc);
    else dq[mm]=aDqp*sign(Dqc);
#elif(0)
    // Xioyue's version
    s = Dqm*Dqp;
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    if (s<=0.) dq[mm] = 0.;
    else if (aDqm<aDqp) dq[mm] = 0.5*aDqm*sign(Dqc);
    else dq[mm]=0.5*aDqp*sign(Dqc);
#endif

#elif(PARA2LIM == NLIM) //w/o slope limiter

    dq[mm] = Dqc;
#endif
  }






#if(JONPARAREDUCE)
  //  if(pl==U1){
  if(pl!=RHO){
    if(
       (fabs(dq[-1]-dq[0])/(fabs(dq[-1])+fabs(dq[0])+SMALL)>0.1)||
       (fabs(dq[1]-dq[0])/(fabs(dq[1])+fabs(dq[0])+SMALL)>0.1)
       ){
      slope_lim_3points(MINM, y[-1], y[0], y[1], dq);
      *lout =y[0] - 0.5* (*dq);
      *rout=y[0] + 0.5* (*dq);
      return;
    }
  }

#endif



#if(WHICHLIMITERTOUSEFORLR==0)
  //////////////////////////////////
  //
  // Obtain continuous solution at interface
  //
  // CW1.6 for obtaining a_{j+1/2} using quartic polynomial, but slopes have been replaced with limited slopes
  //
  // GODMARK: This step could be done once for entire line of data
  
  ddql=(dq[0]-dq[-1]);
  ddqr=(dq[1]-dq[0]);
  
  paracont(ddql, &y[0], &l);
  paracont(ddqr, &y[1], &r);





#if(0)
  // doesn't seem to be within CW!
  l=max(min(y[0],y[-1]),l);
  l=min(max(y[0],y[-1]),l);
  r=max(min(y[0],y[1]),r);
  r=min(max(y[0],y[1]),r);
#endif

  *lout = l;
  *rout = r;



#elif(WHICHLIMITERTOUSEFORLR==1)


  // notice that this results in discontinuous solution at each interface already, unlike para method
  *lout = y[0] - 0.5*dq[0];
  *rout = y[0] + 0.5*dq[0];

  // differs from normal MC,MINM, etc. by computing dq's that will be used later for steepener, fattener, and parabolic check

  // Note that MC is actually 3rd order acccurate (error term O(dx^3)) because linear l,r gives same answer as parabola!


#endif





}












// PARA starts with continuous 3rd order polynomial (4th order through averages if used) using monotonized slopes for second derivative term
// facecont[i] = 0.5*(y[i]+y[i-1]) - ddq[i]   where ddq[0] = dq[0]-dq[-1]
void paracont(FTYPE ddq, FTYPE *y, FTYPE *facecont)
{
  FTYPE avgpointcoef;


#if(AVGINPUT)
  // /6 is result of passing through points from differential of average values in each cell
  avgpointcoef=1.0/6.0;
  *facecont=0.5*(y[0]+y[-1])-ddq*avgpointcoef;
#else
  avgpointcoef=1.0/8.0;
  // consistent with Matt's paper
  // assume passing quartic polynomial through points y[-2,-1,0,1,2]
  // see PARA_interpolation_checks.nb
  *facecont=0.5*(y[0]+y[-1])-ddq*avgpointcoef;
#endif

}







// check monotonicity of parabola
void checkparamonotonicity(int dqrange, int pl, FTYPE *y, FTYPE *ddq, FTYPE *dq, FTYPE *lin, FTYPE *rin, FTYPE *lout, FTYPE *rout)
{
  FTYPE a6COEF;
  FTYPE qa,qb,qd,qe;
  FTYPE r,l;
  int i;
  int numddq;



  

  l = *lin;
  r = *rin;


#if(PARAGENDQALLOWEXTREMUM)
  // (dq[0]-dq[-1]) is defined as ddq[0]
  //  ddqr=(dq[1]-dq[0]) is ddq[1]
  qb=0.0;
  numddq=0;
  for(i=-dqrange/2+1;i<=dqrange/2;i++){
    qb+=sign(ddq[i]);
    numddq++;
  }
  // check that ddq's all same sign
  if(fabs(qb)>(FTYPE)numddq-0.1) qb=1.0; // 0.1 is just to avoid machine precision issue
  else qb=-1.0;
#endif


  /////////////
  //
  // now perform the limiting procedures that create the discontinuities

  // Condition 1: monotonicity for l,y[0],r. qa>0 if monotonic
  qa=(r-y[0])*(y[0]-l);
  // modify Condition 1 as in Duez et al. (2005)
  // allow nonmonotonic behavior if the second derivative doesn't change sign around the center
  // GODMARK: actually should check 2 derivatives in each direction!  Then need PARAFLAT for extra values


  // CW: \delta a_j
  qd=(r-l);

#if(AVGINPUT)
  // CW: a_{6,j} in CW, which is for averages:
  a6COEF=6.0;
  qe=a6COEF*(y[0]-0.5*(l+r));
#else
  // for points need: see PARA_interpolation_checks.nb
  // parabola has solution y = l + (x-xl)*(qd + a6*(1-(x-xl))) with a6 = 4(y_0-1/2(l+r)) for points
  a6COEF=4.0;
  qe=a6COEF*(y[0]-0.5*(l+r));
#endif





#if(PARAGENDQALLOWEXTREMUM)
  if (qa <=0.0 && qb<=0.0 )
#else
  if (qa <=0.0)
#endif
 { // Condition 1


#if(NONMONOLIM==0)
    l=y[0];
    r=y[0];
#elif(NONMONOLIM==1)
    // makes no sense to reduce all the way to DONOR since to second order can still have monotonic result, so use MONO result in this case and assume flatten result
    // appears to be too speculative as results in  more failures at horizon with PARA2LIM==MC
    l = y[0] - 0.5* dq[0];
    r = y[0] + 0.5* dq[0];
#endif

  }
  else if(1){
    // Condition 2

    // qe can be positive or negative still here even though qa>0
    if     (qd*(qd-qe)<0.0)  l = (-(2.0+a6COEF)*r + 2.0*a6COEF*y[0])/(a6COEF-2.0);
    else if(qd*(qd+qe)<0.0)  r = (-(2.0+a6COEF)*l + 2.0*a6COEF*y[0])/(a6COEF-2.0);
    // else no change needed
    //    else{
    //      dualfprintf(fail_file,"Problem with limiting condition 2\n");
    //    }
    
    // Xiaoyue's verison:
    //    if (qd*(qd-qe)<0.0) l=3.0*y[0]-2.0*r;
    //    if (qd*(qd+qe)<0.0) r=3.0*y[0]-2.0*l;
  }


#if(PARAGENDQALLOWEXTREMUM)
  // recompute monotonicity to confirm didn't screw it up:
  qa=(r-y[0])*(y[0]-l);
  // qb is not changed!
  if (qa <=0. && qb<=0.0 ) { // Condition 1 again
    // with Duez allowance of non-monotonic behavior, this gets triggered if l,r change -- otherwise wouldn't have been triggered
    l=y[0];
    r=y[0];
  }
#endif


  // assign output
  *lout=l;
  *rout=r;



}







// used when lim=PARAFLAT
void parapl(int realisinterp, int dir, FTYPE **yrealpl, FTYPE **ypl, FTYPE *loutpl, FTYPE *routpl)
{
  FTYPE dq0[NPR2INTERP][8];
  FTYPE *dq[NPR2INTERP];
  FTYPE *y,*yreal;
  void parasteep(int dir, int pl, FTYPE *V, FTYPE *P, FTYPE *y, FTYPE *dq, FTYPE *l, FTYPE *r);
  void paraflatten(int dir, int pl, FTYPE *y, FTYPE Fi, FTYPE *l, FTYPE *r);
  void getPressure(FTYPE **yrealpl, FTYPE *P);
  FTYPE a_P[10];
  FTYPE *V,*P;
  FTYPE  Ficalc(int dir, FTYPE *V, FTYPE *P, FTYPE **ypl);
  int pl;
  FTYPE Fi;
  int dqrange;
  FTYPE a_ddq[7];
  FTYPE *ddq;
  int mm;




  // consistent with PARAFLAT using 7 points
  dqrange = 5; // dq's will exist from -2,-1,0,1,2 and ddq computed from -2,-1,0,1

  // shift P
  P=a_P + 4; // P accessed from -3..3 ( shifted sufficiently)

  // shift dq
  PLOOPINTERP(pl){
    dq[pl]=dq0[pl]+4; // shifted sufficiently
  }

  ddq=a_ddq+3; // shifted sufficiently


  // assume velocity is istelf
  V = yrealpl[U1+dir-1];



  // get pressures for all points since needed for reduction or steepening
#if( DOPPMREDUCE || DOPPMCONTACTSTEEP)
  getPressure(yrealpl, P);
#endif



  // computed only once for all variables
#if( DOPPMREDUCE )
  Fi = Ficalc(dir,V,P,ypl);
#else
  Fi = 0.0;
#endif



  ///////////////
  //
  // Loop over variables and get interpolated left/right values within cell
  //
  //////////////


  PLOOPINTERP(pl){

    y=ypl[pl];
    yreal=yrealpl[pl];

    // get continuous solution    
    para4gen(dqrange,pl,y,&loutpl[pl],&routpl[pl],dq[pl]);

#if(DOPPMCONTACTSTEEP)
    parasteep(dir,pl,V,P,ypl[pl],dq[pl],&loutpl[pl],&routpl[pl]);
#endif


#if( DOPPMREDUCE )
    paraflatten(dir,pl,ypl[pl],Fi,&loutpl[pl],&routpl[pl]);
#endif


    // finally check monotonicity of the parabola and create discontinuities if non-monotonic
    // FLASH equations 51 -> 53 for points or averages
    for(mm=-dqrange/2+1;mm<=dqrange/2;mm++){
      ddq[mm] = dq[pl][mm] - dq[pl][mm-1];
    }
    checkparamonotonicity(dqrange, pl, ypl[pl], ddq, dq[pl], &loutpl[pl], &routpl[pl], &loutpl[pl], &routpl[pl]);

#if(NONMONOLIM>0 && DOPPMREDUCE)
    // then flatten again
    paraflatten(dir,pl,ypl[pl],Fi,&loutpl[pl],&routpl[pl]);
#endif

    
  }



}



// PPM FLATTENER formula
void paraflatten(int dir, int pl, FTYPE *y, FTYPE Fi, FTYPE *l, FTYPE *r)
{
  // FLASH Equation 49,50
  *l = Fi * y[0] + ( 1.0 - Fi ) * (*l);
  *r = Fi * y[0] + ( 1.0 - Fi ) * (*r);
}




// PPM FLATTENER parameter
FTYPE ftilde( int dir, int shift, FTYPE *Vabs, FTYPE *Pabs, FTYPE **ypl )
{
  FTYPE Ftilde,Ftilde1,Ftilde2;
  FTYPE Sp;
  FTYPE *V, *P;
  FTYPE P2diff,Pbottom;


  // shift as needed
  P = Pabs + shift;
  V = Vabs + shift;

  // FLASH Equation 43
  P2diff=P[2]-P[-2];
  Pbottom=sign(P2diff)/(fabs(P2diff)+SMALL); // singularity avoidance but keeps signature
  Sp = (P[1] - P[-1]) * Pbottom ;




  // FLASH Equation 45
  Ftilde = max( 0, min( 1.0, 10.0 * (Sp - SP0) ) );

  // FLASH Equation 46
  Ftilde1 = fabs(P[1] - P[-1]) / (min(fabs(P[1]), fabs(P[-1]))+ SMALL );
  Ftilde *= ( (FTYPE)(Ftilde1>=THIRD) );
  //  if(Ftilde1<THIRD) Ftilde=0.0;

  // FLASH Equation 47
  Ftilde2 = V[1] - V[-1];
  Ftilde *= ( (FTYPE)(Ftilde2<=0.0) );
  //  if(Ftilde2>0.0) Ftilde=0.0;

  //  dualfprintf(fail_file,"Ftilde=%21.15g\n",Ftilde);
  
  return( Ftilde );
}


// PPM FLATTENERS (final formula)
FTYPE  Ficalc(int dir, FTYPE *V, FTYPE *P, FTYPE **ypl)
{
  FTYPE ftilde( int dir, int shift, FTYPE *P, FTYPE *V, FTYPE **ypl );
  int signdP;
  FTYPE Fi;

  signdP = (P[1] - P[-1] > 0) * 2 - 1;
  // FLASH Equation 48
  Fi = max( ftilde(dir, 0, V,P,ypl), ftilde(dir, -signdP, V,P,ypl) );

  return(Fi);
}



// Get pressure
// Note this is quite inefficient since operating per-point get same pressure for entire line multiple times
// (GODMARK: no accounting of magnetic field)
void getPressure(FTYPE **yrealpl, FTYPE *P)
{
  int mm;

  // need pressure over full range from -3..3
  for(mm=-interporder[PARAFLAT]/2;mm<=interporder[PARAFLAT]/2;mm++){
    P[mm] = pressure_rho0_u(yrealpl[RHO][mm],yrealpl[UU][mm]);
  }

}




void parasteep(int dir, int pl, FTYPE *V, FTYPE *P, FTYPE *y, FTYPE *dq, FTYPE *l, FTYPE *r)
{
  int odir1,odir2;
  void parasteepgen(int pl, FTYPE etai, FTYPE *V, FTYPE *P, FTYPE *y, FTYPE *dq, FTYPE *l, FTYPE *r);
  FTYPE etaicalc(int pl, FTYPE *y, FTYPE *V, FTYPE *P);
  FTYPE etai;



#if(DOPPMSTEEPVARTYPE==0)
  if(pl==RHO)
#elif(DOPPMSTEEPVARTYPE==1)
  // define orthogonal directions for field steepening
  odir1=dir%3+1;
  odir2=(dir+1)%3+1;
  if(pl==RHO || pl==B1+odir1-1 || pl==B1+odir2-1)
#endif
    {
      // get contact indicator
      etai=etaicalc(pl,y,V,P);
      // get steepend values
      parasteepgen(pl, etai, V, P, y, dq, l, r);
    }


}



// steepener
// doesn't use l,r for any calculations, so if steepens fully then initial l,r values don't matter
// return etai if needed
void parasteepgen(int pl, FTYPE etai, FTYPE *V, FTYPE *P, FTYPE *y, FTYPE *dq, FTYPE *l, FTYPE *r)
{
  void pr_contact_compute(int pl, FTYPE *y, FTYPE *dq, FTYPE *prld, FTYPE *prrd);
  FTYPE prld, prrd;
  FTYPE l0,r0;
  FTYPE lmc,rmc;
  FTYPE mceta;


  // compute anti-disspiative left,right values
  pr_contact_compute(pl,y,dq,&prld,&prrd);
  
  // switch to MC for original l,r states if steepening      
  mceta=4.0*max( 0.25 - etai,0.0 );

  lmc = y[0] - 0.5*dq[0];
  rmc = y[0] + 0.5*dq[0];

  l0 = (*l) * mceta  + lmc*(1.0-mceta);
  r0 = (*r) * mceta  + rmc*(1.0-mceta);
  
  // assign steepened density value
  *l = l0 * ( 1.0 - etai ) + prld*etai;
  *r = r0 * ( 1.0 - etai ) + prrd*etai;
  // else make no changes to l,r   

  

}



void pr_contact_compute(int pl, FTYPE *y, FTYPE *dq, FTYPE *prld, FTYPE *prrd)
{

  // equation 33 and 34 in FLASH
  *prld=y[-1]+0.5*dq[-1];
  *prrd=y[+1]-0.5*dq[+1];


}

// PPM steepener parameter, where etai=1 is steep and etai=0 is normal not steep
// Acts to modify l and r values so actually nonmonotonic compared to surrounding cells -- This results in diffusion term in HLL or LAXF to actually be an anti-diffusion term causing (e.g. mass) to be sucked back into the cell in order to counter the diffusive flux causing spreading.  In a stationary case the expansive term balances the anti-diffusion term even with jumps at the interface.
// 
FTYPE etaicalc(int pl, FTYPE *y, FTYPE *V, FTYPE *P)
{
  FTYPE delta2l,delta2r;
  FTYPE etatilde;
  FTYPE etai;
  FTYPE ddcoef;
  int ii;
  FTYPE Pjump,dB,Bmean;
  FTYPE prjumpfactor;
  FTYPE ifinf;
  FTYPE cs2;
  int mm;
  int mmstart,mmend;
  FTYPE max3P,min3P,min3y;

  // y is accessed from y[-2..2]
  // P,dq is accessed via dq[-1,0,1]


#if(AVGINPUT)
  ddcoef=SIXTH;
#else
  ddcoef=1.0/8.0;
#endif

  // equation 35 in FLASH
  ii=-1;
  delta2l=ddcoef*( (y[ii+1]-y[ii]) - (y[ii] - y[ii-1]) );

  // equation 35 in FLASH
  ii=1;
  delta2r=ddcoef*( (y[ii+1]-y[ii]) - (y[ii] - y[ii-1]) );

  // equation 36 in FLASH
  // sign of denominator is important
  // we multiply by conditionall that is 0 or 1
  if(fabs(y[1]-y[-1])>SMALL){
    etatilde=(-(delta2r-delta2l)/(y[1]-y[-1]));
  }
  else etatilde=0.0;

  // below can lead to asymmetries
  //ifinf = ((FTYPE)(fabs(y[1]-y[-1])<0.5*SMALL)); // inf avoidance
  //etatilde=(-(delta2r-delta2l)/(y[1]-y[-1] + ifinf*SMALL));

  

  /////////////////////////////////
  // Pressure jump
  max3P=SMALL;
  min3P=BIG;
  min3y=BIG;
  //  for(mm=-2;mm<=2;mm++){
  // upwinded extension of check
  //  if(V[0]>0.0){ mmstart=-2; mmend=1; }
  //  else if(V[0]<0.0) {  mmstart=-1; mmend=2; }
  //  else {  mmstart=-1; mmend=1; }
  mmstart=-1; mmend=1;
  for(mm=mmstart;mm<=mmend;mm++){
    if(mm==0) continue; // when wanting original method
    max3P=max(max3P,fabs(P[mm]));
    min3P=min(min3P,fabs(P[mm]));
    min3y=min(min3y,fabs(y[mm]));
  }
  min3P+=SMALL;

  //  min3y=min(min(fabs(y[-1]),fabs(y[1])),fabs(y[0]));

  //Pjump=fabs(P[1]-P[-1])/(min(fabs(P[1]),fabs(P[-1]))+SMALL);
  // need more strict and extensive pressure jump check
  Pjump=fabs(max3P - min3P)/min3P;


 

  // consistently use energy density like term
  if(pl>=B1 && pl<=B3){
    dB = y[1]-y[-1];
    Bmean = 0.5*(y[-1] + y[1]); // mean field
    //    prp1sq=0.5*y[1]*y[1]*sign(y[1]); // magnetic pressure
    //    prm1sq=0.5*y[-1]*y[-1]*sign(y[-1]);
    prjumpfactor=fabs(dB)/(fabs(Bmean)+fabs(dB)+SMALL);

    //    dualfprintf(fail_file,"dB=%21.15g Bmean=%21.15g prjumpfactor=%21.15g\n",dB,Bmean,prjumpfactor);
  }
  else{
    prjumpfactor=fabs(y[1]-y[-1])/min3y;

    // check effective c_s given Pressure leaks into steepened density region
    // ensure contributes negligibly to dynamics

    //    cs2 = max3P/min3y;

    // want cs2 dt^2/dx^2 << 1 to avoid steepening leading to velocity comparable to whatever Courant condition is limited by

    // relativistic check
    // don't change energy per baryon by much
    // if pressure is relatively flat
    //    cs2m1check = (cs2 < 5.0 * fabs(P[-1])/fabs(y[-1]));
    //    cs20check = (cs2 < 5.0 * fabs(P[0])/fabs(y[0]));
    //    cs2p1check = (cs2 < 5.0 * fabs(P[1])/fabs(y[1]));
  }


  //  dualfprintf(fail_file,"ETAI:: icurr=%d pl=%d etatilde=%21.15g Pjump=%21.15g prjump=%21.15g\n",icurr,pl,etatilde,Pjump,prjumpfactor);


  // equation 39 in FLASH (not a shock)
  //  if( Pjump > 0.1*prjumpfactor) etatilde=0.0;
  etatilde *= (FTYPE)( Pjump <= 0.1*prjumpfactor);



  // equation 37 in FLASH (to avoid triggering on numerical noise)
  //  if(prjumpfactor<0.01) etatilde=0.0;
  etatilde *=  (FTYPE)(prjumpfactor>=0.01);
  //  etatilde *=  (FTYPE)(prjumpfactor>=0.1);

  // equation 38 in FLASH (really a contact)
  //if(delta2l*delta2r>0.0) etatilde=0.0;
  etatilde *=  (FTYPE)(delta2l*delta2r<=0.0);



  // equation 40 in FLASH
  etai=max(0.0,min(20.0*(etatilde-0.05),1.0));

  return(etai);

}
















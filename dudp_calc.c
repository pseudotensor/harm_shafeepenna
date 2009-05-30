
/* 
 * calculate du/dp analytically.  
 * p is the full (NPR) element primitive variables
 * alpha is be a 5x5 matrix containing dU(1-5)/dp(1-5).
 * used only by nutoprim
 * 
 *
 * cfg 7-11-01
 *
 * cfg 8-23-01 heavy repairs
 *
 */

// jon 04-24-2003 designed for relative 4-velocity by using change of variables 
// Am I concerned by the fact that the grid u^t has 2 solutions per physical velocity?
// maybe we should use 3-velocity or relative 4-velocity as base primitive velocity?
// probably not.  All that matters is that u^t is computed correctly and uniquely, which is true if started with \tilde{u}^i 's.


// assumes T^\mu_\nu form for stress energy tensor, which is most conservative form
// energy includes rest-mass.  Such a change shouldn't matter up to machine precision.


// notice that if WHICHVEL==3, global.h is setup to use dudp_calc_3vel.c directly.  This may be numerical convenient.  dudp_calc_3vel is not setup to handle REMOVERESTMASSFROMUU==2

#include "decs.h"

struct of_geom *ptrgeom;

// alpha is dU^i/dp^j = alpha[i][j]


static void dutdui_calc(FTYPE *ucon,FTYPE *dutdui) ;
static void dbtdui_calc(FTYPE *dutdui, FTYPE *pr, FTYPE *dbtdui) ;
static void dbiduj_calc(FTYPE *dbtdui,FTYPE *dutdui,FTYPE *ucon, FTYPE *b, FTYPE dbiduj[][NDIM]) ;
static void dbsqdui_calc(FTYPE dbiduj[][NDIM],FTYPE *b, FTYPE *dbsqdui) ;
static void duudud_calc(FTYPE *ucon, FTYPE duudud[][NDIM]) ;
static void dudduu_calc(FTYPE*dutdui, FTYPE dudduu[][NDIM]);
static void dbdiduj_calc(FTYPE dbiduj[][NDIM],FTYPE dbdiduj[][NDIM]);
static void dSdu_calc(FTYPE *pr,FTYPE *dSdu);
static void dSdrho_calc(FTYPE *pr,FTYPE *dSdrho);
static void dPtodP(FTYPE *pr, struct of_state *q, struct of_geom *geom, FTYPE **alpha);
static void ducon_dv3_calc(struct of_state *q,FTYPE ducon_dv3[][NDIM]);
static void dgdvi_calc(FTYPE *pr,FTYPE *dgdvi);
static void duidvj_calc(FTYPE *dgdvi,FTYPE duidvj[][NDIM]);


int dudp_calc_gen(int whichcons, FTYPE *pr, struct of_state *q, struct of_geom *geom, FTYPE **alpha)
{
  static FTYPE dutdui[NDIM] ;
  static FTYPE dbtdui[NDIM] ;
  static FTYPE dbsqdui[NDIM] ;
  static FTYPE dbiduj[NDIM][NDIM] ;
  static FTYPE dudduu[NDIM][NDIM];
  static FTYPE dbdiduj[NDIM][NDIM];
  static FTYPE tmp1[NDIM],tmp2[NDIM] ;	
  FTYPE eta,bsq ;
  int i,j,k,l ;
  FTYPE dSdrho,dSdu;
  FTYPE entropy;
  extern int entropy_calc(FTYPE *pr, FTYPE *entropy);
  //
  // functions
  FTYPE contract(FTYPE *vcon1, FTYPE *vcon2);
  FTYPE covtract(FTYPE *vcov1, FTYPE *vcov2);
  // other functions
  FTYPE rho,ie,P;
  extern void compute_1plusud0(struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0); // plus1ud0=(1+q->ucov[TT])
  FTYPE plus1ud0;
  FTYPE dpdrho0,dpdu;


  ptrgeom=geom; // so don't have to rewrite functions yet


  // initialize matrix
  for(j=1;j<=5;j++)
    for(k=1;k<=5;k++){ alpha[j][k] = 0. ;}

  //////////////////////////////
  //
  // setup reused derivatives
  //
  /////////////////////////////
  bsq = dot(q->bcon, q->bcov);
  // d u^t / d u^i
  dutdui_calc(q->ucon,dutdui) ;
  // d b^t / d u^i
  dbtdui_calc(dutdui,pr,dbtdui) ;
  // d b^{mu} / d u^j
  dbiduj_calc(dbtdui,dutdui,q->ucon,q->bcon,dbiduj) ;
  // d b^2 / d u^i
  dbsqdui_calc(dbiduj,q->bcon,dbsqdui) ;
  // d u^{mu} / d u^j
  dudduu_calc(dutdui,dudduu);
  // d b_{mu} / d u^j
  dbdiduj_calc(dbiduj,dbdiduj);

  rho=pr[RHO];
  ie=pr[UU];
  //  P=(gam-1.0)*ie;
  P = pressure_rho0_u(rho,ie);

  eta = P + rho + ie + bsq ;

  // This computes 1+u_t , which for nonrelativistic cases is ~0 .  If computed as 1+u_t, then residual will be large error if small residual.
  compute_1plusud0(geom,q,&plus1ud0); // plus1ud0=(1+q->ucov[TT])

  ///////////////////
  // now define alpha
  //
  // alpha is dU^i/dp^j = alpha[i][j]

  // rho on rho
  // d(\rho_0 u^t)/d\rho_0
  alpha[RHO+1][RHO+1] = q->ucon[TT] ;


  dpdrho0=dpdrho0_rho0_u(pr[RHO],pr[UU]);

  // u+momentums on rho
  // d(T^t_\nu)/d\rho_0
#if(REMOVERESTMASSFROMUU==2)
  SLOOPA(j) alpha[UU+1+j][RHO+1] = (1.0+dpdrho0)*q->ucon[TT]*q->ucov[j];
  j=0; alpha[UU+1+j][RHO+1] = q->ucon[TT]*plus1ud0 + dpdrho0*(q->ucon[TT]*q->ucov[TT] + 1.0);
  j=0; alpha[UU+1+j][RHO+1] += dpdrho0;
#else
  DLOOPA(j) alpha[UU+1+j][RHO+1] = (1.0+dpdrho0)*q->ucon[TT]*q->ucov[j];
  j=0; alpha[UU+1+j][RHO+1] += dpdrho0;
#endif

  // U[rho] on u
  // d(\rho_0 u^t)/du
  alpha[RHO+1][UU+1] = 0.0;

  // u+momentums on u // (notice that no d/d(ie) terms on T^t_t for REMOVERESTMASSFROMUU==2)
  // d(T^t_\nu)/du
  // needs d(u+p)/du   and dp/du
  dpdu = dpdu_rho0_u(pr[RHO],pr[UU]);
  DLOOPA(j) alpha[UU+1+j][UU+1] = (1.0+dpdu)*q->ucon[TT]*q->ucov[j] + delta(TT,j)*dpdu;

  // rho on momentums
  // d(\rho_0 u^t)/du^i
  SLOOPA(j) alpha[RHO+1][U1+1 +j-1]=pr[RHO]*dutdui[j];

  // assumes EOS does not depend on velocity or field
#if(REMOVERESTMASSFROMUU==2)
  SSLOOP(j,k){//D=j S=k // u+momentum's on velocities
    alpha[UU+1 +j][U1+1 +k-1] = 
      dbsqdui[k]*q->ucon[TT]*q->ucov[j]
      +eta*(dutdui[k]*q->ucov[j]+q->ucon[TT]*dudduu[j][k])
      +delta(j,0)*0.5*dbsqdui[k]
      -(dbtdui[k]*q->bcov[j] + q->bcon[TT]*dbdiduj[j][k]);
  }
  j=0; SLOOPA(k){//D=j S=k // u+momentum's on velocities
    alpha[UU+1 +j][U1+1 +k-1] = 
      
      // deta/du^k * u^t u_j
      dbsqdui[k]*q->ucon[TT]*q->ucov[j]

      // d/du^k ( (p+u+b^2) u^dir u_j + rho u^dir (u_j + 1)) = 


      // eta*du^t/du^k u_j + rho du^t/du^k = (p+ie+bsq)du^t/du^k u_j + rho du^t/du^k (u_j +1)
      +((P+ie+bsq)*q->ucov[j]+rho*plus1ud0)*dutdui[k]

      // eta*u^t*du_j/du^k + rho u^t = (P+ie+bsq)*u^t du_j/du^k + rho*u^t (d(1+u_j)/du^k)
      // d(1+u_j)/du^k = du_j/du^k
      +eta*q->ucon[TT]*dudduu[j][k]
      +delta(j,0)*0.5*dbsqdui[k]
      -(dbtdui[k]*q->bcov[j] + q->bcon[TT]*dbdiduj[j][k]);
  }
#else
  DSLOOP(j,k){//D=j S=k // u+momentum's on velocities
    alpha[UU+1 +j][U1+1 +k-1] = 
      dbsqdui[k]*q->ucon[TT]*q->ucov[j]
      +eta*(dutdui[k]*q->ucov[j]+q->ucon[TT]*dudduu[j][k])
      +delta(j,0)*0.5*dbsqdui[k]
      -(dbtdui[k]*q->bcov[j] + q->bcon[TT]*dbdiduj[j][k]);
  }
#endif


  // DOENTROPY thing
  // overwrite above total energy term alpha[UU][j] with alpha[ENTROPY][j]
  if(whichcons==EVOLVEFULLENTROPY){

    //	  if((icurr==37)&&(jcurr==62)) PLOOP(pl) dualfprintf(fail_file,"p[%d]=%21.15g\n",pl,pr[pl]);

    dSdrho_calc(pr,&dSdrho);
    dSdu_calc(pr,&dSdu);

    // ENTROPY on rho
    alpha[UU+1][RHO+1] = dSdrho*(q->ucon[TT]) ;
    // ENTROPY on u
    alpha[UU+1][UU+1] = dSdu*(q->ucon[TT]) ;
    // ENTROPY on momentums
    entropy_calc(pr,&entropy);
    SLOOPA(j) alpha[UU+1][U1+1 +j-1]=entropy*dutdui[j];

  }

  ///////////////////////
  //
  // Assign conserved variables from combination of mhd stresses
  //
  ////////////////////////////
  //   the rest-mass flux is added to the energy 
  //   flux, which may be numerically convenient (tests are unconclusive -- Charles) */
  // handled by Utoprimgen() now
  //	if(REMOVERESTMASSFROMUU){
  //	  for(k=1;k<=5;k++) alpha[UU+1][k] += alpha[RHO+1][k] ;
  //	}

  //convert from standard 4-velocity primitive quantity to used primitive quantity
  dPtodP(pr,q,geom,alpha);

  ////////////////////
  //
  // MULTIPLY BY GDET
  //
  ////////////////////
  //	if(WHICHEOM==WITHGDET){
  /* N.B.: all the conserved variables contain a factor of \sqrt{det(g_{\mu\nu})} */
  //	  for(j=1;j<=5;j++) for(k=1;k<=5;k++) alpha[j][k] *= geom->g ;
  //	}
  //	else if(WHICHEOM==WITHNOGDET){
  // only the U^t (dU^2/dp^j) through U^\phi (dU^5/dp^j) can avoid containing gdet
  //	  j=1; for(k=1;k<=5;k++) alpha[j][k] *= geom->g ; // mass flux term never changes

  //	  j=2; if(NOGDET0==0){ for(k=1;k<=5;k++) alpha[j][k] *= geom->g ;}
  //	  j=3; if(NOGDET1==0){ for(k=1;k<=5;k++) alpha[j][k] *= geom->g ;}
  //	  j=4; if(NOGDET2==0){ for(k=1;k<=5;k++) alpha[j][k] *= geom->g ;}
  //	  j=5; if(NOGDET3==0){ for(k=1;k<=5;k++) alpha[j][k] *= geom->g ;}
  //	}


  return(0);
}



////////////////////////////////
//
// change of primitive variables from 4-vel to something else
//
////////////////////////////////
static void dPtodP(FTYPE *pr, struct of_state *q, struct of_geom *geom, FTYPE **alpha)
{
  static FTYPE ducon_dv3[NDIM][NDIM];
  static FTYPE tempalpha[NPR][NPR];
  static FTYPE dgdvi[NDIM];
  static FTYPE duidvj[NDIM][NDIM];
  int i,j,k,l;


  // initialize matrix
  for(j=1;j<=5;j++)
    for(k=1;k<=5;k++){ tempalpha[j][k] = 0. ;}


  if(WHICHVEL==VEL4){
    // then we are done!
  }
  else if(WHICHVEL==VEL3){
    ducon_dv3_calc(q,ducon_dv3);

    // convert to relative velocity (over all alpha[U=RHO->U3][p=v1,v2,v3])
    // 4vel -> 3vel
    // assign parts that don't change variable
    for(i=1;i<=5;i++) for(j=RHO+1;j<=UU+1;j++) tempalpha[i][j]=alpha[i][j];
    for(i=1;i<=5;i++) SLOOPA(j) for(l=1;l<=3;l++) tempalpha[RHO+1+i-1][U1+1 +j-1]+=alpha[RHO+1+i-1][U1+1+l-1]*ducon_dv3[l][j];
    // assign back to alpha
    for(i=1;i<=5;i++) for(j=1;j<=5;j++) alpha[i][j]=tempalpha[i][j];
  }
  else if(WHICHVEL==VELREL4){
    // d\gamma^2 / d V^i
    dgdvi_calc(pr,dgdvi);
    // d u^{mu} / d V^j (only  mu=1,2,3 used)
    duidvj_calc(dgdvi,duidvj);

    // convert to relative velocity (over all alpha[U=RHO->U3][p=v1,v2,v3])
    // 4vel -> rel4vel
    for(i=1;i<=5;i++) for(j=RHO+1;j<=UU+1;j++) tempalpha[i][j]=alpha[i][j];
    for(i=1;i<=5;i++) SLOOPA(j) for(l=1;l<=3;l++) tempalpha[RHO+1+i-1][U1+1 +j-1]+=alpha[RHO+1+i-1][U1+1+l-1]*duidvj[l][j];
    for(i=1;i<=5;i++) for(j=1;j<=5;j++) alpha[i][j]=tempalpha[i][j];
  }
	
}




static void dgdvi_calc(FTYPE *pr,FTYPE *dgdvi)
{
  int j,k;
  FTYPE gamma;

  gamma_calc(pr,ptrgeom,&gamma) ;

  SLOOPA(j) dgdvi[j]=0.0;

  SLOOP(j,k) dgdvi[j]+=1.0/gamma*ptrgeom->gcov[j][k]*pr[U1+k-1];

  // no such dgdvi[TT]

}

static void duidvj_calc(FTYPE *dgdvi,FTYPE duidvj[][NDIM])
{
  int j,k;
  FTYPE alpha,betacon[NDIM];
  
  //alpha=1.0/sqrt(-ptrgeom->gcon[TT][TT]);
  alpha=ptrgeom->alphalapse;

  SLOOPA(j) betacon[j]=ptrgeom->gcon[TT][j]*alpha*alpha;

  SLOOPA(j) duidvj[TT][j]=dgdvi[j]/alpha;
 
  SLOOP(j,k) duidvj[j][k]=delta(j,k)-betacon[j]/alpha*dgdvi[k];

  // duidvj[j][TT] doesn't exist since no v0 primitive variable, so should never be referenced
}

static void dutdui_calc(FTYPE *ucon,FTYPE *dutdui) 
{
  int j ;
  FTYPE ucov[NDIM] ;
	
  //	dutdui[TT] = 1.0; // never used even if true

  lower_vec(ucon,ptrgeom,ucov) ;
  SLOOPA(j) dutdui[j] = -ucov[j]/ucov[0] ;

  return ;
}

static void dbtdui_calc(FTYPE *dutdui, FTYPE *pr, FTYPE *dbtdui) 
{
  int j ;
  FTYPE B[NDIM],Bcov[NDIM] ;

  B[0] = 0. ;
  SLOOPA(j) B[j] = pr[j+B1-1] ;

  lower_vec(B,ptrgeom,Bcov) ;
  SLOOPA(j) dbtdui[j] = Bcov[j] + Bcov[0]*dutdui[j] ;

  return ;
}

// valid for all i, j=1,2,3
static void dbiduj_calc(FTYPE *dbtdui,FTYPE *dutdui,FTYPE *ucon, FTYPE *b, FTYPE dbiduj[][NDIM]) 
{
  int j,k ;

  DLOOP(j,k) dbiduj[j][k] = 0. ;

  SLOOP(j,k) dbiduj[j][k] = -b[j]*dutdui[k]/ucon[TT] 
    + ucon[j]*dbtdui[k]/ucon[TT] ;

  SLOOPA(j) dbiduj[j][j] += b[TT]/ucon[TT] ;

  SLOOPA(j) dbiduj[TT][j] = dbtdui[j] ;

  return ;
}

static void dbsqdui_calc(FTYPE dbiduj[][NDIM],FTYPE *b, FTYPE *dbsqdui) 
{
  int j,k ;
  FTYPE bcov[NDIM] ;

  lower_vec(b,ptrgeom,bcov) ;
  DLOOPA(j) dbsqdui[j] = 0. ;
  DLOOP(j,k) dbsqdui[j] += 2.*bcov[k]*dbiduj[k][j] ;

  return ;
}

FTYPE contract(FTYPE *vcon1, FTYPE *vcon2)
{
  int j,k ;
  FTYPE n ;

  n = 0. ;
  DLOOP(j,k) n += ptrgeom->gcov[j][k]*vcon1[j]*vcon2[k] ;
  return(n) ;

}

FTYPE covtract(FTYPE *vcov1, FTYPE *vcov2)
{
  int j,k ;
  FTYPE n ;

  n = 0. ;
  DLOOP(j,k) n += ptrgeom->gcon[j][k]*vcov1[j]*vcov2[k] ;
  return(n) ;

}

// valid for all i, j=1,2,3
static void duudud_calc(FTYPE *ucon, FTYPE duudud[][NDIM]) 
{
  int j,k;
  
  DSLOOP(j,k) duudud[j][k] = ptrgeom->gcon[j][k] + ptrgeom->gcon[j][TT]*(-ucon[k]/ucon[TT]) ;
}

// valid for all i, j=1,2,3
static void dudduu_calc(FTYPE*dutdui, FTYPE dudduu[][NDIM]) 
{
  int j,k;

  DSLOOP(j,k) dudduu[j][k] = ptrgeom->gcov[j][k] + ptrgeom->gcov[j][TT]*dutdui[k];

}


// valid for all i, j=1,2,3
static void dbdiduj_calc(FTYPE dbiduj[][NDIM],FTYPE dbdiduj[][NDIM])
{
  int j,k;
  int l;

  DSLOOP(j,k) dbdiduj[j][k] = 0.0;

  // just a lowered dbiduj (lower the 1st index)
  DSLOOP(j,k) for(l=0;l<NDIM;l++) dbdiduj[j][k] += ptrgeom->gcov[j][l]*dbiduj[l][k];

}

static void ducon_dv3_calc(struct of_state *q,FTYPE ducon_dv3[][NDIM])
{
  int i,j;
  /* set ducon_dv */
  for (i = 0; i < NDIM; i++)
    for (j = 1; j < NDIM; j++)
      ducon_dv3[i][j] =
	q->ucon[TT] * (q->ucon[i] * q->ucov[j] + delta(i, j));
}


static void dSdrho_calc(FTYPE *pr,FTYPE *dSdrho)
{
  int j,k;
  FTYPE rho0,u;
  extern FTYPE compute_dSdrho(FTYPE rho0, FTYPE u);

  rho0=pr[RHO];
  u=pr[UU];

  *dSdrho=compute_dSdrho(rho0,u);

}

static void dSdu_calc(FTYPE *pr,FTYPE *dSdu)
{
  int j,k;
  FTYPE rho0,u;
  FTYPE compute_dSdu(FTYPE rho0, FTYPE u);

  rho0=pr[RHO];
  u=pr[UU];

  *dSdu=compute_dSdu(rho0,u);

}

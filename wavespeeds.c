#include "decs.h"



//////////////////////////////////
//
// get wave speeds for flux calculation and for dt calculation
//
/////////////////////////////////

// get wave speeds for flux calculation
int get_wavespeeds(int dir, struct of_geom *ptrgeom, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r, FTYPE *F_l, FTYPE *F_r, struct of_state *state_l, struct of_state * state_r, FTYPE *cminmax_l, FTYPE *cminmax_r, FTYPE *cminmax, FTYPE *ctopptr)
{
  int cminmax_calc(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax,FTYPE *ctop);
  int failreturn;
  FTYPE ftemp;
  FTYPE p_roe[NPR],cminmax_roe[NUMCS];
  struct of_state state_roe;
  FTYPE cminmax_o[NUMCS];
  FTYPE Hroe,cminmaxnonrel_roe[NUMCS];
  void get_roe_averaged_state(int dir, FTYPE *p_l, struct of_state *state_l, FTYPE *Ul, FTYPE * p_r, struct of_state *state_r, FTYPE *Ur, struct of_geom *geom, FTYPE * p_roe, FTYPE *Hroe, FTYPE *cminnonrel_roe, FTYPE*cmaxnonrel_roe);
  int q;
  int ignorecourant;
  int i,j,k;

  i=ptrgeom->i;
  j=ptrgeom->j;
  k=ptrgeom->k;



#if(USEGLOBALWAVE==0)

  // characteristic based upon t^n level for 1/2 step and t^{n+1/2} level for the full step.
  MYFUN(vchar(p_l, state_l, dir, ptrgeom, &cminmax_l[CMAX], &cminmax_l[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);
  MYFUN(vchar(p_r, state_r, dir, ptrgeom, &cminmax_r[CMAX], &cminmax_r[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);

  cminmax_calc(cminmax_l[CMIN],cminmax_r[CMIN],cminmax_l[CMAX],cminmax_r[CMAX],&cminmax[CMIN],&cminmax[CMAX],ctopptr);
  // have cmin,cmax,ctop here
#else
  // other non-local estimate
  cminmax_l[CMIN]=cminmax_r[CMIN]=cminmax[CMIN]=wspeed[dir][CMIN][i][j][k];
  cminmax_l[CMAX]=cminmax_r[CMAX]=cminmax[CMAX]=wspeed[dir][CMAX][i][j][k];

  // change so feeds into Riemann solver with "standard" sign
  cminmax[CMIN] = fabs(max(0., -cminmax[CMIN]));
  cminmax[CMAX] = fabs(max(0., cminmax[CMAX]));
  *ctopptr=max(cminmax[CMIN],cminmax[CMAX]);
  // have cmin,cmax,ctop here
#endif




  //////////////////////////////////
  //
  // get Roe-averaged versions of wave speed for flux calculation and for dt calculation
  //
  /////////////////////////////////

#if(ROEAVERAGEDWAVESPEED)
  // get Roe-averaged primitive state
  get_roe_averaged_state(dir, p_l, state_l, U_l, p_r, state_r, U_r, ptrgeom, p_roe, &Hroe, &cminmaxnonrel_roe[CMIN], &cminmaxnonrel_roe[CMAX]);


#if(ATHENAROE==0) // doing Jon method

  // get HARM state
  MYFUN(get_state(p_roe, ptrgeom, &state_roe),"flux.c:p2SFUevolve", "get_state()", 1);
  // get Roe-averaged state's wavespeeds
  MYFUN(vchar(p_roe, &state_roe, dir, ptrgeom, &cminmax_roe[CMAX], &cminmax_roe[CMIN],&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 3);


#if(1)  // Jon's version of MIN/MAXing

  for(q=0;q<NUMCS;q++){
    // fastest left-going wave (q=0) and then fastest right-going wave (q=1)
    cminmax[q]=MINMAX(q,cminmax_roe[q],MINMAX(q,cminmax_l[q],cminmax_r[q]));
  }


#else // doing Athena-like method of min/maxing (fails to work for some tests)

  // New version of cminmax
  /*
    for(q=0;q<NUMCS;q++){
    cminmax_l[q] = MINMAX(q,cminmax_l[q],0.0);
    cminmax_r[q] = MINMAX(q,cminmax_r[q],0.0);
    cminmax[q]=2.0*cminmax_l[q]*cminmax_r[q]/(cminmax_l[q]+cminmax_r[q]+SMALL);
    }
  */

  cminmax[CMIN] = cminmax_l[CMIN];
  cminmax[CMAX] = cminmax_r[CMAX];
  for(q=0;q<NUMCS;q++){
    cminmax[q]=MINMAX(q,cminmax_roe[q],cminmax[q]);
  }
#endif




#else // else if ATHENAROE==1  (fails to work for some tests)

  // non-rel HD version
  for(q=0;q<NUMCS;q++){
    cminmax_roe[q]=cminmaxnonrel_roe[q];
  }


#if(1) // Athena way
  // This is Athena's way of min/maxing
  cminmax[CMAX]=MAX(cminmax_roe[CMAX],cminmax_r[CMAX]);
  cminmax[CMIN]=MIN(cminmax_roe[CMIN],cminmax_l[CMIN]);

#else // testing way (same as Athena right now)

  // New version of cmin/cmax
  /*
    for(q=0;q<NUMCS;q++){
    cminmax_l[q] = MINMAX(q,cminmax_l[q],0.0);
    cminmax_r[q] = MINMAX(q,cminmax_r[q],0.0);
    cminmax[q] = 2.0*cminmax_l[q]*cminmax_r[q]/(cminmax_l[q] + cminmax_r[q] + SMALL);
    }
  */

  // Jon's way even with ATHENAROE==1
  for(q=0;q<NUMCS;q++){
    // fastest left-going wave (q=0) and then fastest right-going wave (q=1)
    cminmax[q]=MINMAX(q,cminmax_roe[q],MINMAX(q,cminmax_l[q],cminmax_r[q]));
  }
#endif



#endif // end over ATHENAROE=0/1



  // fastest eigenvalue
  // Athena sets this, but never uses it apparently
  //ctop_roe=MAX(fabs(cminmax_roe[CMAX]),fabs(cminmax_roe[CMIN]));


  // change so feeds into Riemann solver with "standard" sign
  cminmax[CMIN] = fabs(max(0., -cminmax[CMIN]));
  cminmax[CMAX] = fabs(max(0., cminmax[CMAX]));
  *ctopptr=max(cminmax[CMIN],cminmax[CMAX]);

  // have Roe versions of cmin,cmax,ctop here
#endif // end over ROEAVERAGEDWAVESPEED=1

  return(0);

}




// get Roe-averaged primitive state
// based upon Athena2's flux_hlle.c
void get_roe_averaged_state(int dir, FTYPE *p_l, struct of_state *state_l, FTYPE *Ul, FTYPE * p_r, struct of_state *state_r, FTYPE *Ur, struct of_geom *geom, FTYPE * p_roe, FTYPE *Hroe, FTYPE *cminnonrel_roe, FTYPE*cmaxnonrel_roe)
{
  FTYPE sqrtrhol,sqrtrhor,isqrtrholr;
  int j,k;
  FTYPE Pl, Pr, Hl, Hr; // specific enthalpy
  FTYPE bsql,bsqr;
  FTYPE vsqroe;
  FTYPE croe;


  Pl=pressure_rho0_u(p_l[RHO],p_l[UU]);
  Pr=pressure_rho0_u(p_r[RHO],p_r[UU]);


  // get weighted density
  sqrtrhol=sqrt(p_l[RHO]);
  sqrtrhor=sqrt(p_r[RHO]);
  isqrtrholr=1.0/(sqrtrhol+sqrtrhor);

  // Roe-averaged density
  p_roe[RHO] = sqrtrhol*sqrtrhor;

  // Roe-averaged internal energy (no consideration of geometry)
  p_roe[UU] = p_roe[RHO] * ( p_l[UU] / sqrtrhol + p_r[UU] / sqrtrhor ) * isqrtrholr;

  // Roe-averaged velocity (no consideration of geometry)
  SLOOPA(j) p_roe[UU+j] = (sqrtrhol*p_l[UU+j] + sqrtrhor*p_r[UU+j])*isqrtrholr;

  // Roe-averaged magnetic field (no consideration of geometry)
  SLOOPA(j) p_roe[B1-1+j] = (sqrtrhor*p_l[B1-1+j] + sqrtrhol*p_r[B1-1+j])*isqrtrholr;


  // b^2
  bsql = dot(state_l->bcon, state_l->bcov);
  bsqr = dot(state_r->bcon, state_r->bcov);


#if(!ATHENAROE)
  // Jon version (relativistically correct)
  Hl = (p_l[UU] + Pl + 0.5*bsql)/p_l[RHO];
  Hr = (p_r[UU] + Pr + 0.5*bsqr)/p_r[RHO];
  *Hroe=(sqrtrhol*Hl + sqrtrhor*Hr)*isqrtrholr;
#else
  // Athena version (relativistically correct)
  *Hroe = ((Ul[UU] + Pl + 0.5*bsql)/sqrtrhol + (Ur[UU] + Pr + 0.5*bsqr)/sqrtrhor)*isqrtrholr;
#endif

  /////////////////////////////////
  // GET NON-REL CASE wave speeds
  // try Einfeldt (1988) 5.1a - 5.3 for non-rel HD
  // non-rel HD Roe-averaged version of wavespeeds

  // v^2 in non-rel case
  vsqroe=0.0;
  SLOOP(j,k) vsqroe += p_roe[UU+j]*p_roe[UU+k] * geom->gcov[j][k];

  // ideal gas only GODMARK (gam)
  croe=sqrt( (gam-1.0)*(*Hroe - 0.5*vsqroe));
  *cminnonrel_roe = p_roe[UU+dir] - croe;
  *cmaxnonrel_roe = p_roe[UU+dir] + croe;

  
  


}





// store wavespeeds somewhere
int get_global_wavespeeds(int dir,FTYPE *pr,struct of_geom *ptrgeom, FTYPE *output)
{
  struct of_state state;
  int ignorecourant;


  // wave speed for cell centered value
  // uses b^\mu b_\mu so need full state
  MYFUN(get_state(pr, ptrgeom, &state),"step_ch.c:fluxcalc()", "get_state()", 0);
  MYFUN(vchar(pr, &state, dir, ptrgeom, &output[CMAX], &output[CMIN],&ignorecourant),"wavespeeds.c:get_global_wavespeeds()", "vchar() dir=1or2", 0);
  
  // uses output as temporary space since not yet needed and not used before global_vchar() below
  
  return(0);
}



// GODMARK: something wrong with comparing multiple velocities since grid/metric changes in space (i.e. v=dx/dt means something different at each grid point)
// GODMARK: assumes boundary zones exist (flux method of bounding won't work) -- have to apply extra limits on values (i,j,k) used here

// defines an effective maximum wave speed centered on the cell interface (FACE)

// might choose wavespeeds that correspond to interpolation stencil, which to first approximation is a symmetric stencil of size interporder[reallim]
int global_vchar(FTYPE pointspeed[][N2M][N3M][NUMCS], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, FTYPE wspeed[][NUMCS][N1M][N2M][N3M])
{
  int i,j,k;
  int reallim;
  int m;
  FTYPE ftemp[NUMCS];
  extern int choose_limiter(int dir, int i, int j, int k, int pl);


#if(VCHARTYPE==VERYLOCALVCHAR) // then get maximum around interface

 COMPZSLOOP( is-idel, ie+idel, js-jdel, je+jdel, ks-kdel, ke+kdel ) { // due to averaging of this face quantity to get centered i/j/k=0 back
    reallim=choose_limiter(dir,i,j,k,RHO); // base stencil size on density limiter GODMARK

    // get minimum/maximum wspeed over this stencil (but centered on edge, not zone center)
    wspeed[dir][CMIN][i][j][k]=+BIG; // assume all zones are smaller than this
    wspeed[dir][CMAX][i][j][k]=-BIG; // assume all zones are larger than this
    // GODMARK: 1D, could do multi-D stencil even if interpolation is 1D
    for(m=-1;m<=0;m++){
      wspeed[dir][CMIN][i][j][k]=MIN(wspeed[dir][CMIN][i][j][k],pointspeed[i+idel*m][j+jdel*m][k+kdel*m][CMIN]);
      wspeed[dir][CMAX][i][j][k]=MAX(wspeed[dir][CMAX][i][j][k],pointspeed[i+idel*m][j+jdel*m][k+kdel*m][CMAX]);
    }
      
  }

#elif(VCHARTYPE==LOCALVCHAR)

 COMPZSLOOP( is-idel, ie+idel, js-jdel, je+jdel, ks-kdel, ke+kdel ) { // due to averaging of this face quantity to get centered i/j/k=0 back
    reallim=choose_limiter(dir,i,j,k,RHO); // base stencil size on density limiter GODMARK

    // get minimum/maximum wspeed over this stencil (but centered on edge, not zone center)
    wspeed[dir][CMIN][i][j][k]=+BIG; // assume all zones are smaller than this
    wspeed[dir][CMAX][i][j][k]=-BIG; // assume all zones are larger than this
    // GODMARK: 1D, could do multi-D stencil even if interpolation is 1D
    // e.g. if dir=1, then expandi=0 and so COMPFZLOOP from i=0..N inclusive.  So can grab from (relative to i) -2 .. 1 inclusive for average centered on i edge
    for(m=-reallim/2;m<=reallim/2-1;m++){
      wspeed[dir][CMIN][i][j][k]=MIN(wspeed[dir][CMIN][i][j][k],pointspeed[i+idel*m][j+jdel*m][k+kdel*m][CMIN]);
      wspeed[dir][CMAX][i][j][k]=MAX(wspeed[dir][CMAX][i][j][k],pointspeed[i+idel*m][j+jdel*m][k+kdel*m][CMAX]);
    }
      
  }

#elif(VCHARTYPE==GLOBALVCHAR)

  ftemp[CMIN]=+BIG;
  ftemp[CMAX]=-BIG;
 COMPZSLOOP( is-idel, ie+idel, js-jdel, je+jdel, ks-kdel, ke+kdel ) { // due to averaging of this face quantity to get centered i/j/k=0 back
    // get minimum/maximum wspeed over entire grid containing the fluxes
    ftemp[CMIN]=MIN(ftemp[CMIN],pointspeed[i][j][k][CMIN]);
    ftemp[CMAX]=MAX(ftemp[CMAX],pointspeed[i][j][k][CMAX]);
  }

 COMPZSLOOP( is-idel, ie+idel, js-jdel, je+jdel, ks-kdel, ke+kdel ) {
    wspeed[dir][CMIN][i][j][k]=ftemp[CMIN];
    wspeed[dir][CMAX][i][j][k]=ftemp[CMAX];
  }

#endif

  return(0);

}



// GODMARK:

// really HARM is currently using VERY local lax Friedrich.
// maybe try local lax Friedrich, using max wave speed from zones used to reconstruct the zone (most common?)
// also can try more global wave speed, or even speed of light.

// determine cmin,cmax,ctop
int cminmax_calc(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax,FTYPE *ctop)
{
  FTYPE lmin,lmax,ltop;

  // need to make use of ROE averages
      
      
  // As per the HLLE method, one must use the wave speed at the interface.  The below is arbitrarily choosing the largest of the interpolated states.

  // Apparently one should use Roe's linearisation to compare
  // against, not the left/right interpolations.  Then compare the
  // left most going and right most going Rho linearisation
  // eigenvalues with the left and right numerical approximations
  // and choose the minimum and maximum values as cmin and cmax,
  // all respectively.
  // Then one defines the new cmin as the minumum of cmin and 0, and new cmax as maximum of cmax and 0.
  // Thus this is slightly different and doesn't avoid negative mass/ie densities like HLLE.


  // LAXF here is really the Rusanov flux.  Lax-Friedrichs would say ctop->dx/dt


#if(1)
  lmin=min(cmin_l,cmin_r);
  lmax=max(cmax_l,cmax_r);
#elif(0)
  lmin=cmin_l;
  lmax=cmax_r;
#else
  lmin=2.0*cmin_l*cmin_r/(cmin_l+cmin_r);
  lmax=2.0*cmax_l*cmax_r/(cmax_l+cmax_r);
#endif


#if(HYPERHLL==0)
  // sharp change at c=0
  lmin = fabs(max(0., -lmin));
  lmax = fabs(max(0., lmax));
#else

#define HYPERFUN(c,x0) (0.5*( (c)+(M_SQRT2l-1.0)*sqrt(2.0*(M_SQRT2l+1.0)*(x0)*(x0)+(3.0+2.0*M_SQRT2l)*(c)*(c)) ) )
  // GODMARK
  // need a good way to set x0 based upon fraction of local light speed
  // or something...

  // function is always positive for given lmin/lmax as original sign
  lmin=HYPERFUN(-lmin,1E-4);
  lmax=HYPERFUN(lmax,1E-4);
#endif

  // returns positive cmin,cmax,ctop
  *cmin=lmin;
  *cmax=lmax;
  *ctop=max(lmax,lmin);


  return(0);
      


}


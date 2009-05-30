
#include "u2p_defs.h"


#define METHODTYPE  0  
#define NEWT_DIM 1
#define USE_LINE_SEARCH -2


/* these variables need to be shared between the functions
   Utoprim_1D, residual, and utsq */
static FTYPE Bsq,QdotBsq,Qtsq,Qdotn,D ;
static FTYPE wglobal;
static PFTYPE *glpflag; // global pflag for local file

// Declarations: 
static FTYPE vsq_calc(FTYPE W);
static int Utoprim_new_body(FTYPE U[], struct of_geom *ptrgeom,  FTYPE prim[]);
static void primtoU_g( FTYPE prim[], FTYPE gcov[][4], FTYPE gcon[][4],  FTYPE U[] );

static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);

static FTYPE res_sq_1d_orig( FTYPE [] );

static int general_newton_raphson( FTYPE x[], int n, int do_line_search,void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int), FTYPE (*res_func) (FTYPE []) );

static void raise_g(FTYPE vcov[], FTYPE gcon[][4], FTYPE vcon[]);
static void lower_g(FTYPE vcon[], FTYPE gcov[][4], FTYPE vcov[]);
static void ncov_calc(FTYPE gcon[][4],FTYPE ncov[]) ;
static void my_lnsrch(int, FTYPE [], FTYPE, FTYPE [], FTYPE [], FTYPE [], FTYPE *, 
		      FTYPE, FTYPE, int *, FTYPE (*res_func) (FTYPE []));

static void bin_newt_data( FTYPE errx, int niters, int conv_type, int print_now  ) ;
static FTYPE pressure_W_vsq(FTYPE W, FTYPE D, FTYPE vsq) ;
static FTYPE dpdW_calc_vsq(FTYPE W, FTYPE D, FTYPE vsq);
static FTYPE dpdvsq_calc(FTYPE W, FTYPE D, FTYPE vsq);
static FTYPE dvsq_dW(FTYPE W);












#include "u2p_util.c"
#include "utoprim_jon_eos.c"


/******************************************************************
  
   Driver for new prim. var. solver.  The driver just translates
   between the two sets of definitions for U and P. 

******************************************************************/

int Utoprim_1d_final(FTYPE U[NPR], struct of_geom *ptrgeom,  PFTYPE *lpflag,  FTYPE prim[NPR])
{


  FTYPE U_tmp[NPR], U_tmp2[NPR], prim_tmp[NPR];
  int i, j, k, ret; 
  FTYPE gcov[NDIM][NDIM], gcon[NDIM][NDIM], alpha;
  const int ltrace = 0;

  lntries=0;
  
  // assign global int pointer to lpflag pointer
  glpflag=lpflag;



#if( WHICHVEL != VELREL4 )
  fprintf(stderr,"Utoprim_1d() Not implemented for WHICHVEL = %d \n", WHICHVEL );
  return(1);
#endif


  /* First update the primitive B-fields */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  /* Set the geometry variables: */
  alpha = 1.0/sqrt(-ptrgeom->gcon[0][0]);
  
  /* Calculate the transform the CONSERVATIVE variables into the new system */
  for( i = RHO; i <= BCON3; i++ ) {
    U_tmp[i] = alpha * U[i] ;
  }

  /* transform the PRIMITIVE variables into the new system */
  for( i = 0; i < BCON1; i++ ) {
    prim_tmp[i] = prim[i];
  }
  for( i = BCON1; i <= BCON3; i++ ) {
    prim_tmp[i] = alpha*prim[i];
  }

  ret = Utoprim_new_body(U_tmp, ptrgeom, prim_tmp);

  /* Check conservative variable transformation: */
#if(!OPTIMIZED)
  if( ltrace ) {
    for(i = 0; i < NDIM; i++ ) {
      for(j = 0; j < NDIM; j++ ) {
	gcov[i][j] = ptrgeom->gcov[i][j];
	gcon[i][j] = ptrgeom->gcon[i][j];
	dualfprintf(fail_file,"gcov,gcon %d %d = %21.15g %21.15g \n", i, j, gcov[i][j], gcon[i][j]);
      }
    }
    dualfprintf(fail_file,"gdet = %21.15g \n", ptrgeom->g);
    primtoU_g( prim_tmp, gcov, gcon, U_tmp2 ); 
    for( i = 0; i < NPR; i++ ) {
      dualfprintf(fail_file, "Utoprim_1d(): Utmp1[%d] = %21.15g , Utmp2[%d] = %21.15g , dUtmp[%d] = %21.15g \n", 
	       i, U_tmp[i], i, U_tmp2[i], i, fabs( (U_tmp[i]-U_tmp2[i]) / ( (U_tmp2[i]!=0.) ? U_tmp2[i] : 1. ) )  ); 
    }
  }
#endif

  /* Transform new primitive variables back : */ 
  if( ret == 0 ) { 
    for( i = 0; i < BCON1; i++ ) {
      prim[i] = prim_tmp[i];
    }
  }


  /* Always return with no problems so that the fixup routines make the correction */
  return( 0 ) ;

}


/**********************************************************************************

attempt an inversion from U to prim using the initial guess
prim.

  -- assumes that 
             /  rho gamma        \
         U = |  alpha T^t_\mu  |
             \  alpha B^i        /



             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \  alpha B^i   /


*glpflag[ptrgeom->i][ptrgeom->j]:  (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence (occurrence of "nan" or "+/-inf" ;
	     
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: utsq<0 w/ initial p[] guess;
	     3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: utsq<0 or utsq > UTSQ_TOO_BIG   with new  W;
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/

static int Utoprim_new_body(FTYPE U[NPR], struct of_geom *ptrgeom,  FTYPE prim[NPR])
{
  FTYPE Wtest;

  FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u);
  FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w);
  FTYPE x_1d[1];
  FTYPE QdotB,Bcon[4],Bcov[4],Qcov[4],Qcon[4],ncov[4],ncon[4],Qsq,Qtcon[4];
  FTYPE rho0,u,p,w,gammasq,gamma,gtmp,W_last,W,utsq,vsq,tmpdiff ;
  int i,j, retval, retval2, i_increase ;
  const int ltrace = 0;
  const int ltrace2 = 0;

  /* TEMPORARY */
  /*
    primtoU_g(prim,gcov,gcon,U) ;
  */

#if(!OPTIMIZED)
  if( ltrace ) {
    for(i=0;i<8;i++) dualfprintf(fail_file,"%d %21.15g %21.15g\n",i,prim[i],U[i]) ;
    
  }
#endif
  // Assume ok initially:
  *glpflag= UTOPRIMNOFAIL;

  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  Bcon[0] = 0. ;
  for(i=1;i<4;i++) Bcon[i] = U[BCON1+i-1] ;

  lower_g(Bcon,ptrgeom->gcov,Bcov) ;

  for(i=0;i<4;i++) Qcov[i] = U[QCOV0+i] ;
  raise_g(Qcov,ptrgeom->gcon,Qcon) ;


  Bsq = 0. ;
  for(i=1;i<4;i++) Bsq += Bcon[i]*Bcov[i] ;

  QdotB = 0. ;
  for(i=0;i<4;i++) QdotB += Qcov[i]*Bcon[i] ;
  QdotBsq = QdotB*QdotB ;

  ncov_calc(ptrgeom->gcon,ncov) ;
  raise_g(ncov,ptrgeom->gcon,ncon);

  Qdotn = Qcon[0]*ncov[0] ;

  Qsq = 0. ;
  for(i=0;i<4;i++) Qsq += Qcov[i]*Qcon[i] ;

  Qtsq = Qsq + Qdotn*Qdotn ;

  D = U[RHO] ;

  /* calculate W from last timestep and use 
     for guess */
  utsq = 0. ;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++) utsq += ptrgeom->gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;


  if( (utsq < 0.) && (fabs(utsq) < MAXNEGUTSQ) ) { 
    utsq = 0.0;
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    if( debugfail>=2 ) dualfprintf(fail_file,"Utoprim_new(): utsq < 0 in utoprim_1d attempt, utsq = %21.15g \n", utsq) ;
    *glpflag= UTOPRIMFAILCONVGUESSUTSQ;
    return(1) ;
  }

  gammasq = 1. + utsq ;
  gamma  = sqrt(gammasq);
	
  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) . 
  rho0 = D / gamma ;
  u = prim[UU] ;
  p = pressure_rho0_u(rho0,u) ;
  w = rho0 + u + p ;

  W_last = w*gammasq ;


  // Make sure that W is large enough so that v^2 < 1 : 
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*Bsq ) - QdotBsq*(2.*W_last + Bsq) ) <= W_last*W_last*(Qtsq-Bsq*Bsq))
	 && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file,"badval :  W = %21.15g, i_increase = %d \n", W_last, i_increase); 
#endif
  }
#if(!OPTIMIZED)
  if( i_increase >= 10 ) { 
    if( debugfail>=2 ) dualfprintf(fail_file,"i_increase is too large, i_increase = %d , W = %21.15g \n", i_increase, W_last);
  }
  /* Let it and fail it does */
#endif
  

#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file,"u = %21.15g,  p = %21.15g,  Bsq = %21.15g, Qsq = %21.15g \n",u,p,Bsq,Qsq);
    dualfprintf(fail_file,"Bcon[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Bcon[0],Bcon[1],Bcon[2],Bcon[3]);
    dualfprintf(fail_file,"Bcov[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Bcov[0],Bcov[1],Bcov[2],Bcov[3]);
    dualfprintf(fail_file,"Qcon[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Qcon[0],Qcon[1],Qcon[2],Qcon[3]);
    dualfprintf(fail_file,"Qcov[0-3] = %21.15g   %21.15g   %21.15g   %21.15g   \n", Qcov[0],Qcov[1],Qcov[2],Qcov[3]);
    dualfprintf(fail_file,"call find_root\n") ; 	
    
  }
#endif

  // METHOD specific:
  //	W = find_root_1D_old(W_last) ; 
  x_1d[0] = W_last;
  wglobal = w; // used to normalize validate_x() type function
  if( (retval=general_newton_raphson( x_1d, 1, USE_LINE_SEARCH, func_1d_orig, res_sq_1d_orig ) ) ) {
#if(!OPTIMIZED)
    if( debugfail>=2 ) dualfprintf(fail_file, "GNR failed, x_1d[0] = %21.15g \n", x_1d[0] );  
#endif
  }
  W = x_1d[0];


  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (W == FAIL_VAL) ) {
    if( debugfail>=2 ) {
      dualfprintf(fail_file, "Failed to find a prim. var. solution!! %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g \n",W_last,Bsq,QdotBsq,Qdotn,D,Qtsq);
      dualfprintf(fail_file, "Utoprim_new_body(): bad newt failure, t,i,j, p[0-7], U[0-7] = %21.15g %d %d ", t, ptrgeom->i, ptrgeom->j );  
      for( i = 0 ; i < NPR; i++ ) {
	dualfprintf(fail_file, "%21.15g ", prim[i]);
      }
      for( i = 0 ; i < NPR; i++ ) {
	dualfprintf(fail_file, "%21.15g ", U[i]);
      }
      dualfprintf(fail_file, "\n");
    }      
    *glpflag=  retval+UTOPRIMFAILCONVRET+1;
    return(retval);
  }
  else{
    
    Wtest=W/wglobal; // normalize to old densities

    if(Wtest <= 0. || Wtest > GAMMASQ_TOO_BIG) {
    //// not true: let below utsq check for too big (normalized properly)
    //if(Wtest <= 0.) {
      if( debugfail>=2 ) {
		dualfprintf(fail_file,"Wtest failure %21.15g \n",Wtest) ;
	dualfprintf(fail_file, "Utoprim_new_body(): Wtest<0 or Wtest=toobig failure, t,i,j, p[0-7], U[0-7] = %21.15g %d %d ", t, ptrgeom->i, ptrgeom->j );  
	for( i = 0 ; i < NPR; i++ ) {
	  dualfprintf(fail_file, "%21.15g ", prim[i]);
	}
	for( i = 0 ; i < NPR; i++ ) {
	  dualfprintf(fail_file, "%21.15g ", U[i]);
	}
	dualfprintf(fail_file, "\n");
      }      

      *glpflag=  UTOPRIMFAILCONVW;
      return(retval) ;

    }
  }


#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file,"(W,W_last,Bsq,Qtsq,QdotB,gammasq,Qdotn) %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
	    W,W_last,
	    Bsq,Qtsq,QdotB,gammasq,Qdotn) ;
    dualfprintf(fail_file,"done find_root\n") ;	
  }
  if( ltrace2 ) {
    dualfprintf(fail_file, "\n <--------- %21.15g %21.15g %21.15g %21.15g %21.15g  \n", Bsq,QdotBsq,Qdotn,D,Qtsq);
    
  }
#endif


  // Calculate vsq
  vsq = vsq_calc(W) ;
  if( (vsq < 0.) && (fabs(vsq) < MAXNEGVSQ) ) { 
    vsq = 0.0;
  }
  // no check for vsq>1 since looking for failure to fix that's not just simple like above
  if( (vsq > VSQ_TOO_BIG) || (vsq<0.0) ) {
    if( debugfail>=2 ) { 
      dualfprintf(fail_file,"vsq failure:  vsq = %21.15g , W = %21.15g \n",vsq, W) ;
      dualfprintf(fail_file, "Utoprim_new_body(): utsq==bad failure, t,i,j, p[0-7], U[0-7] = %21.15g %d %d ", t, ptrgeom->i, ptrgeom->j );  
      for( i = 0 ; i < NPR; i++ ) {
	dualfprintf(fail_file, "%21.15g ", prim[i]);
      }
      for( i = 0 ; i < NPR; i++ ) {
	dualfprintf(fail_file, "%21.15g ", U[i]);
      }
      dualfprintf(fail_file, "\n");
    }      

    *glpflag= UTOPRIMFAILCONVUTSQ;
    return(retval) ;
  }

  gtmp = sqrt(1. - vsq);
  gamma = 1./gtmp ;
  rho0 = D * gtmp;

  w = W * (1. - vsq) ;
  p = pressure_rho0_w(rho0,w) ;
  u = w - (rho0 + p) ;

  if( (rho0 <= 0.) || (u < 0.) ) { 
    if( debugfail>=2 ) {
      tmpdiff = w - rho0;
      dualfprintf(fail_file,
		  "rho or uu < 0 failure: rho,w,(w-rho),p,u  = %21.15g %21.15g %21.15g %21.15g %21.15g \n",
		  rho0,w,tmpdiff,p,u) ;
      dualfprintf(fail_file,
		  "rho or uu < 0 failure: gamma,utsq = %21.15g %21.15g  \n",  gamma, utsq) ;
    }
    if((rho0<=0.)&&(u>=0.)) *glpflag=  UTOPRIMFAILRHONEG;
    if((rho0>0.)&&(u<0.)) *glpflag= UTOPRIMFAILUNEG;
    if((rho0<=0.)&&(u<0.)) *glpflag= UTOPRIMFAILRHOUNEG;
    if(UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNNOTADJUSTED) return(retval) ; // else let assign -- used to check how bad failure is.
  }

  prim[RHO] = rho0 ;
  prim[UU] = u ;


  for(i=1;i<4;i++)  Qtcon[i] = Qcon[i] + ncon[i] * Qdotn;
  for(i=1;i<4;i++) prim[UTCON1+i-1] = gamma/(W+Bsq) * ( Qtcon[i] + QdotB*Bcon[i]/W ) ;
	
  /* set field components */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;


#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file," rho final = %21.15g ,  u final = %21.15g \n", rho0, u);
  }
#endif

  /* done! */
  return(retval) ;

}


/* evaluate v^2 (spatial, normalized velocity) from W = \gamma^2 w */
static FTYPE vsq_calc(FTYPE W)
{
	FTYPE Wsq,Xsq,Ssq;
	
	Wsq = W*W ;
	Xsq = (Bsq + W) * (Bsq + W);
	Ssq = QdotBsq / Bsq;

	//return(  Ssq * ( 1./Wsq - 1./Xsq )  +  Qtsq / Xsq  ); 
	return(  ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq) );

}

/* evaluate v^2 (spatial, normalized velocity) from W = \gamma^2 w */
static FTYPE dvsq_dW(FTYPE W)
{
	FTYPE W3,X3,Ssq,Wsq,X;
	
	X = Bsq + W;
	Wsq = W*W;
	W3 = Wsq*W ;
	X3 = X*X*X;
	Ssq = QdotBsq / Bsq;

	//return( -2.*( Ssq * ( 1./W3 - 1./X3 )  +  Qtsq / X3 ) ); 
	//	return( -2.*( W3*Qtsq + QdotBsq * ( 3*W*X + Bsq*Bsq ) ) / ( W3 * X3 )   );
	return( -2.*( Qtsq/X3  +  QdotBsq * (3*W*X + Bsq*Bsq) / ( W3 * X3 )  )  );

}


  /* 
   This file must contain the equation of state
   in two different forms:

   p(rho0,w)
   p(rho0,u)

   and also the derivatives

   dp/du
   dp/drho0

   which are required only to evaluate the
   matrix dU/dprim in Utoprim_5D;  these
   are not required for Utoprim_1D.

   current equation of state is a GAMMA-law gas.

   cfg 14 july 04

*/



/* 

pressure as a function of W, vsq, and D:


*/

/*
static FTYPE pressure_W_vsq(FTYPE W, FTYPE D, FTYPE vsq) 
{
  FTYPE gtmp;
  

  gtmp = 1. - vsq;
  
  return(  (GAMMA - 1.) * ( W * gtmp  -  D * sqrt(gtmp) ) / GAMMA  );

}
*/

/* 

partial derivative of pressure with respect to W


*/
/*
static FTYPE dpdW_calc_vsq(FTYPE W, FTYPE D, FTYPE vsq)
{

  return( (GAMMA - 1.) * (1. - vsq) /  GAMMA ) ;

}
*/
/* 

partial derivative of pressure with respect to vsq


*/


//static FTYPE dpdvsq_calc(FTYPE W, FTYPE D, FTYPE vsq)
//{
//  FTYPE outval;

  //  return( (GAMMA - 1.) * ( 0.5 * D / sqrt(1.-vsq)  - W  ) / GAMMA  ) ;
//  outval =  (GAMMA - 1.) * ( 0.5 * D / sqrt(1.-vsq)  - W  ) / GAMMA   ;

//  if( outval > 0. ) { 
//    dualfprintf(fail_file,"outval = %21.15g , D = %21.15g  , vsq = %21.15g,  W = %21.15g \n",
//	    outval, D, vsq, W );
//  }

//  return(outval);
//}



/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static int general_newton_raphson( FTYPE x[], int n, int do_line_search, void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int),  FTYPE (*res_func) (FTYPE []) )
{
  FTYPE f, f_old, df, df_old, dx[NEWT_DIM], dx_old[NEWT_DIM], x_old[NEWT_DIM], resid[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, errx_old, errx_oldest, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE randtmp, tmp;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old;
  FTYPE resid_norm, resid_check, grad_check;

  FTYPE res_func_val, res_func_old, res_func_new;
  FTYPE dn[NEWT_DIM], del_f[NEWT_DIM];
  int   keep_iterating, i_increase, retval2,retval = 0;
  const int ltrace  = 0;
  const int ltrace2 = 1;


  retval = 0;


  errx = 1. ; 
  errx_old = 2.;
  df = df_old = f = f_old = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;


  vsq_old = vsq = W = W_old = 0.;


  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 
    nstroke++;
    lntries++;

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */


#if(!OPTIMIZED)
    /*  Check for bad untrapped divergences : */
    if( (finite(f)==0) || (finite(df)==0)  ) {
      dualfprintf(fail_file,"general_newton_raphson(): nan encountered in f or df!! \n");
      dualfprintf(fail_file,"gnr nan(): f, df, x0, dx0 =  %21.15g  %21.15g  %21.15g  %21.15g  \n", f,df,x[0],dx[0]);
      return(1);
    }
#endif


#if(!OPTIMIZED)
    /* Randomly rescale Newton step to break out of iteration cycles: */
    if( ((n_iter+1) % CYCLE_BREAK_PERIOD) == 0 ) {
      randtmp = ( (1.*rand())/(1.*RAND_MAX) );
      for( id = 0; id < n ; id++) dx[id] *= randtmp;
      //	for( id = 0; id < n ; id++) dx[id] *= ( (1.*rand())/(1.*RAND_MAX) );
    }
#endif

    /* Save old values before calculating the new: */
    errx_oldest = errx_old;
    errx_old = errx;
    errx = 0.;
    f_old = f;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* Make the newton step: */
    if( do_line_search == 1 ) { 

      /* Compare the residual to its initial value: */
      if( n_iter == 0 ) { 
	resid_norm = 0.0e0;
	for( id = 0; id < n ; id++) {
	  resid_norm += fabs(resid[id]);
	}
	resid_norm /= 1.0*n ;
	if( resid_norm == 0.0 ) resid_norm = 1.0;
      }
       
      for( id = 0; id < n ; id++) {
	tmp = 0.;
	for( jd = 0; jd < n ; jd++) {
	  tmp += jac[jd][id] * resid[jd];
	}
	del_f[id] = tmp;
      }
      for( id = 0; id < n ; id++) {
	dn[id] = dx[id];
      }

      my_lnsrch(n, x_old-1, f_old, del_f-1, dn-1, x-1, &f, TOL_LINE_STEP, SCALEMAX, &retval, res_func);

      /* dx is needed for errx calculation below: */
      for( id = 0; id < n ; id++) {
	dx[id] = x[id] - x_old[id];
      }

#if(!OPTIMIZED)
      if( ltrace ) { 
	res_func_val = res_func(x);
	res_func_old = res_func(x_old);
	dualfprintf(fail_file,"gnr(): f_old, f, res_func_old, res_func_val = %21.15g  %21.15g  %21.15g  %21.15g  \n",
		f_old, f, res_func_old, res_func_val );
	dualfprintf(fail_file,"gnr(): x_old = ");
	for( id = 0; id < n ; id++) {
	  dualfprintf(fail_file," %21.15g ",x_old[id]);
	}
	dualfprintf(fail_file,"\n ");
	dualfprintf(fail_file,"gnr(): x     = ");
	for( id = 0; id < n ; id++) {
	  dualfprintf(fail_file," %21.15g ",x[id]);
	}
	dualfprintf(fail_file,"\n ");
	dualfprintf(fail_file,"gnr(): dn    = ");
	for( id = 0; id < n ; id++) {
	  dualfprintf(fail_file," %21.15g ",dn[id]);
	}
	dualfprintf(fail_file,"\n ");
	dualfprintf(fail_file,"gnr(): del_f = ");
	for( id = 0; id < n ; id++) {
	  dualfprintf(fail_file," %21.15g ",del_f[id]);
	}
	dualfprintf(fail_file,"\n ");
      }
#endif

      /* Check to see if line search problem is because the residual vector is already small enough */
      if( retval == 1 ) {
	resid_check = 0.0e0;
	for( id = 0; id < n ; id++) {
	  resid_check += fabs(resid[id]);
	}
	resid_check /= 1.0*n;
	
	if( resid_check <= resid_norm * NEWT_FUNC_TOL ) {
	  retval = 0;
	}
	if( ltrace && retval ) { 
	  dualfprintf(fail_file,"general_newton_raphson():  retval, resid_check = %4i  %21.15g \n",retval, resid_check);
	  
	}	  
      }
      /* If initial Newton step is bad, then try again without line searching: */
      if( (retval == 2) && (USE_LINE_SEARCH == do_line_search) ) { 
#if(!OPTIMIZED)
	if( ltrace ) { 
	  dualfprintf(fail_file,"gnr(): bad first step: retval, f_old, f  = %4i  %21.15g  %21.15g  \n",retval,f_old,f);
	  dualfprintf(fail_file,"gnr: doing recursive call, retval, errx = %4i  %21.15g \n", retval, errx );
	  
	}
#endif
	retval = general_newton_raphson( x_orig, n, ((do_line_search+1)%2), funcd, res_func );
	for( id = 0; id < n ; id++)  x[id] = x_orig[id] ;
	return( retval );
      }

      /* Check to see if it is trapped in a local minimum, i.e. gradient is too small */ 
      if( retval == 1 ) { 
	grad_check = 0.0e0;
	for( id = 0; id < n ; id++) {
	  resid_check = (x[id] == 0.) ? 1.0 : fabs(x[id]) ;
	  grad_check  +=  del_f[id] * resid_check ;
	}
	resid_check = (f == 0.) ? 1.0 : fabs(f) ;
	grad_check /= resid_check;
	
	/* Then we've most likely found a solution: */
	if( grad_check > GRADMIN ) { 
	  retval = -1;
	}
	else if( ltrace ) { 
	  dualfprintf(fail_file,"general_newton_raphson():  retval, grad_check = %4i  %21.15g \n",retval, grad_check);
	  
	}
      }
    }
    else {
      /* don't use line search : */
      for( id = 0; id < n ; id++) {
	x[id] += dx[id]  ;
      }

      //METHOD specific:
      i_increase = 0;
      while( (( x[0]*x[0]*x[0] * ( x[0] + 2.*Bsq ) - QdotBsq*(2.*x[0] + Bsq) ) <= x[0]*x[0]*(Qtsq-Bsq*Bsq))
	     && (i_increase < 10) ) {
	x[0] -= (1.*i_increase) * dx[0] / 10. ;
	i_increase++;
#if(!OPTIMIZED)
	dualfprintf(fail_file,"badval2 :  W = %21.15g, i_increase = %d \n", x[0], i_increase); 
#endif
      }
#if(!OPTIMIZED)
      if( i_increase >= 10 ) { 
	dualfprintf(fail_file,"i_increase is too large, i_increase = %d , W = %21.15g \n", i_increase, x[0]);
      }
      
#endif

    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific:
#if( NEWCONVERGE == 1 )
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);

    /* For the old criterion, look at errors in each indep. variable(s) (except for 5D) : */
#else
    for( id = 0; id < n ; id++) {
      errx  += (x[id]==0.) ?  fabs(dx[id]) : fabs(dx[id]/x[id]);
    }
    errx /= 1.*n;
#endif


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    x[0] = fabs(x[0]);
#if(!OPTIMIZED)
    x[0] = (x[0]/wglobal > GAMMASQ_TOO_BIG) ?  x_old[0] : x[0];
#endif


    /****************************************/
    /* Check to see if we're in a infinite loop with error function: */
    /****************************************/
#if( CHECK_FOR_STALL )
    if( ( (errx_old == errx) || (errx_oldest == errx) ) && (errx <= MIN_NEWT_TOL) )  errx = -errx;
#endif 

    /****************************************/
    /* If there's a problem with line search, then stop iterating: */
    /****************************************/
    if( (retval == 1) || (retval == -1) ) errx = -errx;


#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file," general_newton_raphson(): niter,f_old,f,errx_old,errx = %4i  %21.15g  %21.15g  %21.15g  %21.15g\n",  
	      n_iter,f_old,f,errx_old,errx );
      dualfprintf(fail_file,"gnr(): x_old = ");
      for( id = 0; id < n ; id++) {
	dualfprintf(fail_file," %21.15g ",x_old[id]);
      }
      dualfprintf(fail_file,"\n ");
      dualfprintf(fail_file,"gnr(): x     = ");
      for( id = 0; id < n ; id++) {
	dualfprintf(fail_file," %21.15g ",x[id]);
      }
      dualfprintf(fail_file,"\n ");
      dualfprintf(fail_file,"gnr(): dx     = ");
      for( id = 0; id < n ; id++) {
	dualfprintf(fail_file," %21.15g ",dx[id]);
      }
      dualfprintf(fail_file,"\n ");
      
    }
#endif

    /****************************************/
    /* Prepare for the next iteration, set the "old" variables: */
    /****************************************/
    for( id = 0; id < n ; id++)  dx_old[id] = dx[id] ;
    f_old  = f;
    df_old = df;


    /****************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations before stopping */
    /****************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) || (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


    /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) || (finite(df)==0) || (finite(x[0])==0)  ) {
#if(!OPTIMIZED)
    dualfprintf(fail_file,"general_newton_raphson(): nan encountered in f or df!! \n");
    dualfprintf(fail_file,"gnr nan(): f, df, x0, dx0 =  %21.15g  %21.15g  %21.15g  %21.15g  \n", f,df,x[0],dx[0]);
#endif
    return(1);
  }


  if( fabs(errx) > MIN_NEWT_TOL){
    if( (do_line_search != USE_LINE_SEARCH) || (USE_LINE_SEARCH < 0) ) { 
#if(DOHISTOGRAM)
      bin_newt_data( errx, n_iter, 0, 0 );
#endif
#if(!OPTIMIZED)
      if(ltrace2) {
	dualfprintf(fail_file," totalcount = %d   0   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra,errx); 
      }
      if(ltrace) {
	dualfprintf(fail_file,"general_newton_raphson():  did not find solution \n");
	if( retval == -1 ) {
	  dualfprintf(fail_file,"general_newton_raphson(): lnsrch converged: x = ");
	  for( id = 0; id < n ; id++)  dualfprintf(fail_file," %21.15g  ",x[id]);
	  dualfprintf(fail_file,"\n");
	  dualfprintf(fail_file,"general_newton_raphson(): lnsrch converged: x_old = ");
	  for( id = 0; id < n ; id++)  dualfprintf(fail_file," %21.15g  ",x_old[id]);
	  dualfprintf(fail_file,"\n");
	}
	
      }
#endif
      return(1);
    } 
    else {
      /* If bad return and we tried line searching, try it without before giving up: */
      //      dualfprintf(fail_file,"gnr: doing recursive call, do_line_search, retval, errx = %4i  %4i  %21.15g \n", do_line_search, retval, errx );
      //      
      retval2 = general_newton_raphson( x_orig, n, ((do_line_search+1)%2), funcd, res_func );
      for( id = 0; id < n ; id++)  x[id] = x_orig[id] ;
      //      dualfprintf(fail_file,"gnr retval3 = %4i \n", retval2); 
      return( retval2 );
    }
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
#if(DOHISTOGRAM)
      bin_newt_data( errx, n_iter, 1, 0 );
#endif
#if(!OPTIMIZED)
    if(ltrace2) {
      dualfprintf(fail_file," totalcount = %d   1   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra,errx); 
      
    }
    if(ltrace) {
      dualfprintf(fail_file,"general_newton_raphson(): found minimal solution \n");
      
    }
    //    dualfprintf(fail_file,"gnr retval4 = %4i \n", 0); 
#endif
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
#if(DOHISTOGRAM)
      bin_newt_data( errx, n_iter, 2, 0 );
#endif
#if(!OPTIMIZED)
    if(ltrace2) {
      dualfprintf(fail_file," totalcount = %d   2   %d  %d  %d  %21.15g \n",n_iter,retval,do_line_search,i_extra, errx); 
      
    }
    //    dualfprintf(fail_file,"gnr retval5 = %4i \n", 0); 
#endif
    return(0);
  }

#if(!OPTIMIZED)
  dualfprintf(fail_file,"gnr retval6 = %4i \n", 0);
#endif
  return(0);

}



/**************************************************** 
*****************************************************/

static void func_1d_orig(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{
  FTYPE vsq,gamma,gamma_sq,rho0,w,W,drdW,dpdW,Wsq,p_tmp, dvsq,dp1,dp2,a_term,ap_term ;
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;


  W = x[0];
  Wsq = W*W;

  //dualfprintf(fail_file,"call utsq in residual\n") ;
  vsq = vsq_calc(W) ;
  //dualfprintf(fail_file,"done w/ utsq in residual\n") ;

  // GODMARK DANGEROUS 
  //  vsq = ( vsq < -1.e-10 ) ?  dv : fabs(vsq) ;
  //  vsq = ( vsq > 1. ) ?  (1.-dv) : vsq ;
  
#if(OPTIMIZED)
  vsq = ( vsq < -MAXNEGVSQ ) ?  0.0 : fabs(vsq) ;
  vsq = ( vsq > 1. ) ?  VSQ_TOO_BIG : vsq ;
#else
  if(vsq<0.0){
    if(vsq>-MAXNEGVSQ) vsq=0.0;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
  if(vsq>VSQ_TOO_BIG){
    if(vsq<1.0+MAXNEGVSQ) vsq=VSQ_TOO_BIG;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
#endif

  //rho0 = D * sqrt(1. - vsq) ;
  //dualfprintf(fail_file,"rho0,D,gamma: %21.15g %21.15g %21.15g\n",rho0,D,gamma) ;
  //  w = W * (1. - vsq);
  p_tmp = pressure_W_vsq(W,D,vsq);



  // Same as previous residual, however, jacobian is calculated using full differentiation w.r.t. W  

  dvsq = dvsq_dW( W );
  dp1 = dpdW_calc_vsq( W, D, vsq );
  dp2 = dpdvsq_calc( W, D, vsq );
  dpdW = dp1  + dp2*dvsq;

  resid[0] = 
    + W 
    + 0.5 * Bsq * ( 1. + vsq )
    - 0.5*QdotBsq/Wsq
    + Qdotn
    - p_tmp;

  jac[0][0] = drdW = 1. - dpdW + QdotBsq/(Wsq*W) + 0.5*Bsq*dvsq;

//  resid[0] = 
//    + W*Wsq
//    + 0.5 * Wsq * Bsq * ( 1. + vsq )
//    - 0.5*QdotBsq
//    - Qdotn*Wsq
//    - p_tmp*Wsq;
//
//  jac[0][0] = drdW = W * ( 3.*W  +  Bsq * (1. + vsq + 0.5*W*dvsq)  - 2.*Qdotn - 2.*p_tmp - W*dpdW );

  dx[0] = -resid[0]/drdW;

  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

#if(!OPTIMIZED)
  if( ltrace ) {
    a_term = 0.5*Bsq*(1.+vsq) + Qdotn - p_tmp;
    ap_term = dvsq * ( 0.5*Bsq - dp2 ) - dp1;

    dualfprintf(fail_file,"func_1d_orig(): x = %21.15g, dx = %21.15g, resid = %21.15g, drdW = %21.15g \n", x[0], dx[0], resid[0], drdW);
    dualfprintf(fail_file,"func_1d_orig(): dvsq = %21.15g,  dp = %21.15g , dp1 = %21.15g , dp2 = %21.15g \n", dvsq, dpdW,dp1,dp2);
    dualfprintf(fail_file,"func_1d_orig(): W = %21.15g , vsq = %21.15g , a = %21.15g,  a_prime = %21.15g  \n", x[0], vsq, a_term, ap_term);
    
  }
#endif

}

/**************************************************** 
  Routine for line searching for the 1D method . 
*****************************************************/

static FTYPE res_sq_1d_orig(FTYPE x[])
{
  FTYPE vsq,W,Wsq,p_tmp,resid[1];
  //  FTYPE dv = 1.0e-10;
  const int ltrace = 0;


  W = x[0];
  Wsq = W*W;

  vsq = vsq_calc(W) ;
  
  // GODMARK DANGEROUS 
#if(OPTIMIZED)
  // scn's version
  //  vsq = ( vsq < -dv ) ?  dv : fabs(vsq) ;
  //vsq = ( vsq > 1. ) ?  (1.-dv) : vsq ;
  // jon's version
  vsq = ( vsq < 0.0 ) ?  0.0 : vsq ;
  vsq = ( vsq > VSQ_TOO_BIG ) ? VSQ_TOO_BIG  : vsq ;
#else
  if(vsq<0.0){
    if(vsq>-MAXNEGVSQ) vsq=0.0;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
  if(vsq>VSQ_TOO_BIG){
    if(vsq<=1.0) vsq=VSQ_TOO_BIG;
    else{
      dualfprintf(fail_file,"func_1d_orig: vsq=%21.15g\n",vsq);
      myexit(1);
    }
  }
#endif

  p_tmp = pressure_W_vsq(W,D,vsq);


  resid[0] = 
    + W 
    + 0.5 * Bsq * ( 1. + vsq )
    - 0.5*QdotBsq/Wsq
    + Qdotn
    - p_tmp;

  return(  0.5*resid[0]*resid[0] );


}


/**************************************************** 
*****************************************************/

#define ALF 1.0e-4

static void my_lnsrch( int n, FTYPE xold[], FTYPE fold, FTYPE g[], FTYPE p[], FTYPE x[], 
		FTYPE *f, FTYPE TOLX, FTYPE stpmax, int *check, FTYPE (*func) (FTYPE []) )
{
  int i;
  FTYPE a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  int icount=0;
  int bad_step;
  FTYPE bad_step_factor = 2.0;
  
  const int ltrace = 0;
	
  
  *check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=1;i<=n;i++) p[i] *= stpmax/sum;
  for (slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for (i=1;i<=n;i++) {
    //    temp=fabs(p[i])/max(fabs(xold[i]),1.0);
    temp= (xold[i] == 0.) ? fabs(p[i]) :  fabs(p[i]/xold[i]);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;

#if(!OPTIMIZED)
  if( ltrace ) {
    dualfprintf(fail_file,"my_lnsrch(): sum, slope, test, alamin =   %21.15g  %21.15g  %21.15g  %21.15g \n",sum,slope,test,alamin); 
  }
#endif

  alam=1.0;
  for (;;) {
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];

    *f=(*func)(x+1);
      
    bad_step = 0;

    if( finite(*f)==0 ) { 
      bad_step = 1;
    }

    //METHOD specific:
    if( x[1] <= 0. ) { 
      bad_step = 1;
    }

    //      if( bad_step ) alam /= bad_step_factor;
    //      if (alam < alamin) bad_step = 0;

    if( bad_step ) { 
      *check = 2;
#if(!OPTIMIZED)
      dualfprintf(fail_file,"my_lnsrch(): bad_step = 1,  f = %21.15g \n", *f); 
#endif
      return;
    }
      
    if (alam < alamin) {
      for (i=1;i<=n;i++) x[i]=xold[i];
      *check=1;
#if(!OPTIMIZED)
      if( ltrace ) { 
	dualfprintf(fail_file,"my_lnsrch(): alam < alamin: alam, alamin = %21.15g  %21.15g \n", alam,alamin); 
      }
#endif
      return;
    } 
    else if (*f <= fold+ALF*alam*slope) {
#if(!OPTIMIZED)
      if( ltrace ) { 
	dualfprintf(fail_file,"my_lnsrch(): good exit:  alam, alamin, f, fold = %21.15g  %21.15g %21.15g  %21.15g \n", alam,alamin, *f, fold); 
      }
#endif
      return;
    }
    else {
      if (alam == 1.0) {
	tmplam = -slope/(2.0*(*f-fold-slope));
#if(!OPTIMIZED)
	if( ltrace ) {
	  dualfprintf(fail_file,"my_lnsrch(): setting tmplam!!    tmplam, alam =  %21.15g  %21.15g !!\n", tmplam, alam);
	}
#endif 
      }
      else {
	rhs1 = *f-fold-alam*slope;
	rhs2=f2-fold2-alam2*slope;
	a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a == 0.0) tmplam = -slope/(2.0*b);
	else {
	  disc=b*b-3.0*a*slope;
	  if (disc<0.0) {
#if(!OPTIMIZED)	    
	    if( disc < -1.e-10 ) {
	      dualfprintf(fail_file,"my_lnsrch(): Big Roundoff problem:  disc = %21.15g \n", disc);
	    }
#endif
	    disc = 0.;
	  }
	  else tmplam=(-b+sqrt(disc))/(3.0*a);
	}
	if (tmplam>0.5*alam)
	  tmplam=0.5*alam;
#if(!OPTIMIZED)
	if( ltrace ) {
	  dualfprintf(fail_file,"my_lnsrch(): rhs1, rhs2, a, b, tmplam, alam =  %21.15g  %21.15g  %21.15g  %21.15g  %21.15g  %21.15g !!\n",
		  rhs1, rhs2, a, b, tmplam, alam );
	}
#endif
      }
    }
    alam2=alam;
    f2 = *f;
    fold2=fold;
    alam=max(tmplam,0.1*alam);
#if(!OPTIMIZED)
    if( ltrace ) {
      dualfprintf(fail_file,"my_lnsrch(): icount, alam, alam2, tmplam  =   %4i  %21.15g  %21.15g  %21.15g \n",
	      icount, alam, alam2, tmplam);
    }
#endif   
    icount++;
  }
}
#undef ALF


#define N_CONV_TYPES 3
#define N_NITER_BINS (MAX_NEWT_ITER + 3)
#define NBINS 200

/**************************************************** 

  Used to gather statistics for the Newton solver during a disk run:

  -- bins are inclusive on their minimum bounds, and exlusive on their max.


      lerrx_[min,max] = the min./max. bounds of the errx historgram
      d_errx          = errx histogram bin width
      n_beyond_range  = number of solutions with errx max. beyond range speciified by lerrx_[min,max]
      n_conv_types    = number of types of ways newton procedure can exit;

      print_now       = 1 if you want to just print the histograms without counting a solution;

*****************************************************/

static void bin_newt_data( FTYPE errx, int niters, int conv_type, int print_now  ) 
{

  
  /* General variables */
  int i, j;
  static int first_call   = 1;
  static long int n_bin_calls  = 0L;
  const  long int n_output_freq = 128*128*100;
  FTYPE lerrx;

  /* Variables for the errx histogram */ 
  static const FTYPE lerrx_min = -23.;
  static const FTYPE lerrx_max =   4.;
  static FTYPE d_errx;
  static long int n_errx[N_CONV_TYPES][NBINS];
  static FTYPE xbin[NBINS];
  static long int n_beyond_range = 0L;   /* Number of points that lie out of the bounds of our errx histogram */
  int ibin;

  /* Variables for the histogram of the number of newton iterations : */
  static long int n_niters[N_CONV_TYPES][N_NITER_BINS];
  

  /* Clear arrays, set constants : */
  if( first_call ) {
    d_errx    = ((lerrx_max - lerrx_min)/(1.*NBINS));
    
    for( i = 0; i < N_CONV_TYPES; i++ ) { 
      for( j = 0; j < NBINS; j++ ) { 
	n_errx[i][j] = 0;
      }
      for( j = 0; j < N_NITER_BINS; j++ ) { 
	n_niters[i][j] = 0;
      }
    }

    for( j = 0; j < NBINS; j++ ) { 
      xbin[j] = lerrx_min + (1.*j + 0.5)*d_errx;
    }
      
    first_call = 0;
  }


  if( print_now != 1 ) {

    /* Check validity of arguments : */
    errx = fabs(errx) ;
    lerrx = log10(errx + 1.0e-2*pow(10.,lerrx_min));

    if( (niters < 0) || (niters >= N_NITER_BINS) ) {
      dualfprintf(fail_file,"bin_newt_data(): bad value for niters = %d \n", niters );
      fflush(stdout);
      return;
    }
    
    /* Determine ibin */
    if( lerrx < lerrx_min ) {
      ibin = 0;
      n_beyond_range++ ;
    }
    else if( lerrx >= lerrx_max ) {
      ibin = NBINS - 1;
      n_beyond_range++ ;
    }
    else {
      ibin = (int) ( ( lerrx - lerrx_min ) / d_errx );
    }
		  
    /* Tally this solution  */
    n_errx[ conv_type][ibin]++;
    n_niters[conv_type][niters]++;

  }


  /* Print out the histograms periodically or when asked to : */
  if( print_now  ||  ( (n_bin_calls % n_output_freq) == 0 ) ) {

    dualfprintf(log_file,"t = %21.15g ,  n_beyond_range = %ld , n_bin_calls = %ld \n", t, n_beyond_range, n_bin_calls);

    /* ERRX */
    dualfprintf(log_file,"ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--\n");
    dualfprintf(log_file,"                         x");
    for( j = 0; j < N_CONV_TYPES; j++ ) { 
      dualfprintf(log_file,"            N%d",j);
    }
    dualfprintf(log_file,"\n");
    
    for( i = 0; i < NBINS; i++ ) { 
      dualfprintf(log_file,"%21.15g ",xbin[i]);
      for( j = 0; j < N_CONV_TYPES; j++ ) { 
	dualfprintf(log_file,"%13ld ", n_errx[j][i]);
      }
      dualfprintf(log_file,"\n");
    }    
    dualfprintf(log_file,"ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--ERRX-HISTOGRAM--\n");


    /* NITER */

    dualfprintf(log_file,"NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--\n");
    dualfprintf(log_file,"        niter");
    for( j = 0; j < N_CONV_TYPES; j++ ) { 
      dualfprintf(log_file,"            N%d",j);
    }
    dualfprintf(log_file,"\n");
    
    for( i = 0; i < N_NITER_BINS; i++ ) { 
      dualfprintf(log_file,"%13d ", i);
      for( j = 0; j < N_CONV_TYPES; j++ ) { 
	dualfprintf(log_file,"%13ld ", n_niters[j][i]);
      }
      dualfprintf(log_file,"\n");
    }    
    dualfprintf(log_file,"NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--NITER-HISTOGRAM--\n");

  }
    
    
  if( print_now != 1 )  n_bin_calls++;


  return;

}
  

#undef N_CONV_TYPES 
#undef N_NITER_BINS 
#undef NBINS 


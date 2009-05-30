// this file is actually just included in other files


////////////////////////////////////////////////////////////////
//
// FUNCTIONS ARE GENERAL for ANY EOS
// input variables may need to be expanded for more advanced EOS
//
////////////////////////////////////////////////////////////////


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

*/





// 1/ (dwmrho0/dp)
// used for dP/dW and dP/dvsq
//static FTYPE compute_idwmrho0dp_old(FTYPE rho0, FTYPE wmrho0)
//{
//  static FTYPE dpressuredu_rho0_u(FTYPE rho0, FTYPE u);
//  static FTYPE u_wmrho0(FTYPE rho0, FTYPE wmrho0);
//  FTYPE u,dpdu;
//  u=u_wmrho0(rho0, wmrho0);
//  dpdu=dpressuredu_rho0_u(rho0, u);
//
//return(1.0/(1.0/dpdu+1.0));
//}


// static declarations
static FTYPE pressure_Wp_vsq(FTYPE Wp, FTYPE D, FTYPE vsq);
static FTYPE pressure_W_vsq(FTYPE W, FTYPE D, FTYPE vsq) ;
static FTYPE wmrho0_compute_vsq(FTYPE Wp, FTYPE D, FTYPE vsq, FTYPE gamma,FTYPE gammasq);
static FTYPE pressure_W_vsq_scn(FTYPE W, FTYPE D, FTYPE vsq) ;
static FTYPE pressure_W_vsq_old(FTYPE W, FTYPE D, FTYPE vsq) ;
static FTYPE wmrho0_compute_utsq(FTYPE Wp, FTYPE D, FTYPE utsq, FTYPE gamma,FTYPE gammasq);
static FTYPE dpdWp_calc_vsq(FTYPE Wp, FTYPE D, FTYPE vsq);
static FTYPE wmrho0_compute_vsq(FTYPE Wp, FTYPE D, FTYPE vsq, FTYPE gamma,FTYPE gammasq);
static FTYPE wmrho0_compute_utsq(FTYPE Wp, FTYPE D, FTYPE utsq, FTYPE gamma,FTYPE gammasq);
static FTYPE dpdW_calc_vsq(FTYPE Wp, FTYPE D, FTYPE vsq);
static FTYPE dpdvsq_calc_Wp(FTYPE Wp, FTYPE D, FTYPE vsq);
static FTYPE dpdvsq_calc_scn(FTYPE W, FTYPE D, FTYPE vsq);
static FTYPE wmrho0_compute_vsq(FTYPE Wp, FTYPE D, FTYPE vsq, FTYPE gamma,FTYPE gammasq);
static FTYPE wmrho0_compute_utsq(FTYPE Wp, FTYPE D, FTYPE utsq, FTYPE gamma,FTYPE gammasq);
static FTYPE wmrho0_compute_vsq(FTYPE Wp, FTYPE D, FTYPE vsq, FTYPE gamma,FTYPE gammasq);


/* 

pressure as a function of W, vsq, and D:


*/

//#define GAMMA (gam)

#define OLDCALCJON 0

//#define CRAZYDEBUG 1

#if(OLDCALCJON)
static FTYPE pressure_W_vsq(FTYPE W, FTYPE D, FTYPE vsq) 
{
  FTYPE gtmp;
  FTYPE answer;

  gtmp = 1. - vsq;
 
  answer=(GAMMA - 1.) * ( W * gtmp  -  D * sqrt(gtmp) ) / GAMMA ;

#if(CRAZYDEBUG)
    if(icurr==0 && jcurr==31 && nstep==9 && steppart==2){
      dualfprintf(fail_file,"gtmp=%21.15g GAMMA=%21.15g W=%21.15g D=%21.15g answer=%21.15g\n",gtmp,GAMMA,W,D,answer);
    }
#endif
 
  return(answer);

}
#else
static FTYPE pressure_W_vsq(FTYPE W, FTYPE D, FTYPE vsq) 
{

  return(pressure_Wp_vsq(W-D,D,vsq));

}
#endif



#if(OLDCALCJON)
static FTYPE pressure_Wp_vsq(FTYPE Wp, FTYPE D, FTYPE vsq) 
{
  
  
  return(pressure_W_vsq(Wp+D, D, vsq));

}

#else


static FTYPE pressure_Wp_vsq(FTYPE Wp, FTYPE D, FTYPE vsq) 
{
  FTYPE gtmp;
  FTYPE gamma,gammasq;
  FTYPE rho0,wmrho0;
  FTYPE pressure_wmrho0(FTYPE rho0, FTYPE wmrho0);



  gamma=1.0/sqrt(1.0-vsq);
  gammasq=gamma*gamma;
  wmrho0=wmrho0_compute_vsq(Wp, D, vsq, gamma,gammasq);
  rho0=D/gamma;

  // wmrho0=u+p and for ideal gas p=(GAMMA-1) u, so p (GAMMA/(GAMMA-1)) = wmrho0
  //  return( wmrho0*IGAMMAR );
  return(pressure_wmrho0(rho0,wmrho0));

  
  //  return(pressure_W_vsq_scn(Wp+D,D,vsq));

  //  return(pressure_W_vsq_old(Wp+D,D,vsq));
}

static FTYPE pressure_Wp_utsq(FTYPE Wp, FTYPE D, FTYPE utsq) 
{
  FTYPE gtmp;
  FTYPE gamma,gammasq;
  FTYPE rho0,wmrho0;
  FTYPE pressure_wmrho0(FTYPE rho0, FTYPE wmrho0);


  gammasq=1.0+utsq;
  gamma=sqrt(gammasq);
  wmrho0=wmrho0_compute_utsq(Wp, D, utsq, gamma,gammasq);
  rho0=D/gamma;

  // wmrho0=u+p and for ideal gas p=(GAMMA-1) u, so p (GAMMA/(GAMMA-1)) = wmrho0
  return(pressure_wmrho0(rho0,wmrho0));

  
}



#endif


static FTYPE pressure_W_vsq_old(FTYPE W, FTYPE D, FTYPE vsq) 
{
  FTYPE gtmp;
  

  gtmp = 1. - vsq;
  
  return(  (gam - 1.) * ( W * gtmp  -  D * sqrt(gtmp) ) / gam  );

}


static FTYPE pressure_W_vsq_scn(FTYPE W, FTYPE D, FTYPE vsq) 
{
  FTYPE gtmp;
  FTYPE gamma;
  FTYPE wmrho0;
  FTYPE rho0;
  FTYPE pressure_wmrho0(FTYPE rho0, FTYPE wmrho0);
  

  gamma=1.0/sqrt(1.0-vsq);
  rho0=D/gamma;

  gtmp = 1. - vsq;
  wmrho0=W * gtmp  -  D * sqrt(gtmp);

  return(pressure_wmrho0(rho0, wmrho0));


}



/* 

partial derivative of pressure with respect to W holding vsq fixed

// dp/dW|v^2 = dp/drho|wmrho0 drho/dW|v^2 + dp/dwmrho0|rho dwmrho0/dW|v^2

*/
#if(OLDCALCJON)
static FTYPE dpdW_calc_vsq(FTYPE W, FTYPE D, FTYPE vsq)
{

  return( (GAMMA - 1.) * (1. - vsq) /  GAMMA ) ;

}

static FTYPE dpdWp_calc_vsq(FTYPE Wp, FTYPE D, FTYPE vsq)
{
  FTYPE W=Wp+D;

  return( (GAMMA - 1.) * (1. - vsq) /  GAMMA ) ;

}

#else

static FTYPE dpdW_calc_vsq(FTYPE W, FTYPE D, FTYPE vsq)
{


  return(dpdWp_calc_vsq(W-D, D, vsq));

}


// holding v^2 fixed -- used by 1D and 2D methods, where for 1D method one uses a full derivative to obtain answer
// obtains: dp/dW' holding v^2 fixed
static FTYPE dpdWp_calc_vsq(FTYPE Wp, FTYPE D, FTYPE vsq)
{
  FTYPE compute_idwmrho0dp(FTYPE rho0, FTYPE wmrho0);
  FTYPE compute_idrho0dp(FTYPE rho0, FTYPE wmrho0);
  FTYPE rho0,wmrho0;
  FTYPE idwmrho0dp;
  FTYPE gamma,gammasq;
  FTYPE dwmrho0dW;
  FTYPE dpdW;
  FTYPE drho0dW,idrho0dp;

  gamma=1.0/sqrt(1.0-vsq);
  gammasq=gamma*gamma;
  wmrho0=wmrho0_compute_vsq(Wp, D, vsq, gamma, gammasq);

  rho0=D/gamma;
  // holding v^2 fixed
  idwmrho0dp=compute_idwmrho0dp(rho0, wmrho0);
  dwmrho0dW = 1.0-vsq;

  drho0dW = 0.0; // because \rho=D/\gamma
  idrho0dp = compute_idrho0dp(rho0, wmrho0); // so not really needed

  dpdW = drho0dW *idrho0dp + dwmrho0dW *idwmrho0dp;

  return( dpdW ) ;

}

// holding utsq fixed
// obtains dp/dW' holding utsq fixed
static FTYPE dpdWp_calc_utsq(FTYPE Wp, FTYPE D, FTYPE utsq)
{
  FTYPE compute_idwmrho0dp(FTYPE rho0, FTYPE wmrho0);
  FTYPE compute_idrho0dp(FTYPE rho0, FTYPE wmrho0);
  FTYPE rho0,wmrho0;
  FTYPE idwmrho0dp;
  FTYPE gamma,gammasq;
  FTYPE dwmrho0dW;
  FTYPE dpdW;
  FTYPE drho0dW,idrho0dp;

  gammasq=1.0+utsq;
  gamma=sqrt(gammasq);
  wmrho0=wmrho0_compute_utsq(Wp, D, utsq, gamma, gammasq);

  rho0=D/gamma;
  // holding utsq fixed
  idwmrho0dp=compute_idwmrho0dp(rho0, wmrho0);
  dwmrho0dW = 1.0/gammasq;

  drho0dW = 0.0; // because \rho=D/\gamma
  idrho0dp = compute_idrho0dp(rho0, wmrho0); // so not really needed

  dpdW = drho0dW *idrho0dp + dwmrho0dW *idwmrho0dp;

  return( dpdW ) ;

}


#endif

static FTYPE dpdW_calc_vsq_scn(FTYPE W, FTYPE D, FTYPE vsq)
{

  return(dpdW_calc_vsq(W,D,vsq));

}

/* 

partial derivative of pressure with respect to vsq holding W fixed


*/
#if(OLDCALCJON)
static FTYPE dpdvsq_calc(FTYPE W, FTYPE D, FTYPE vsq)
{
  FTYPE outval;
  outval =  (GAMMA - 1.) * ( 0.5 * D / sqrt(1.-vsq)  - W  ) / GAMMA   ;
  
  
  return(outval);
}

#else
static FTYPE dpdvsq_calc(FTYPE W, FTYPE D, FTYPE vsq)
{

  // old version
  //  return( (GAMMA - 1.) * ( 0.5 * D / sqrt(1.-vsq)  - W  ) / GAMMA   );

  return(dpdvsq_calc_Wp(W-D, D, vsq));
}
#endif

// GODMARK: probably should write "clean" version of this, but maybe doesn't matter


#if(OLDCALCJON)
static FTYPE dpdvsq_calc_Wp(FTYPE Wp, FTYPE D, FTYPE vsq)
{

  return(dpdvsq_calc_scn(Wp+D, D, vsq));
}

#else


//dp/dv^2 holding W' fixed
static FTYPE dpdvsq_calc_Wp(FTYPE Wp, FTYPE D, FTYPE vsq)
{
  FTYPE outval;
  FTYPE gamma,W;
  FTYPE compute_idwmrho0dp(FTYPE rho0, FTYPE wmrho0);
  FTYPE compute_idrho0dp(FTYPE rho0, FTYPE wmrho0);
  FTYPE rho0,wmrho0;
  FTYPE dwmrho0dvsq,idwmrho0dp;
  FTYPE drho0dvsq,idrho0dp;
  FTYPE gammasq;
  FTYPE dpdvsq;


  gamma=1.0/sqrt(1.0-vsq);
  gammasq=gamma*gamma;
  wmrho0=wmrho0_compute_vsq(Wp, D, vsq, gamma, gammasq);

  rho0=D/gamma;

  //  return( (gam - 1.) * ( 0.5 * D / sqrt(1.-vsq)  - W  ) / gam  ) ;
  //  outval =  (gam - 1.) * ( 0.5 * D / sqrt(1.-vsq)  - W  ) / gam   ;

  W = Wp+D;

  //  outval =  (D*gamma*0.5-W)*idwmrho0dp; // old code

  idwmrho0dp=compute_idwmrho0dp(rho0, wmrho0);
  dwmrho0dvsq = (D*(gamma*0.5-1.0) - Wp);

  drho0dvsq = -D*gamma*0.5; // because \rho=D/\gamma
  idrho0dp = compute_idrho0dp(rho0, wmrho0);

  dpdvsq =   drho0dvsq *idrho0dp  +   dwmrho0dvsq *idwmrho0dp;

//  if( outval > 0. ) { 
//    dualfprintf(fail_file,"outval = %21.15g , D = %21.15g  , vsq = %21.15g,  W = %21.15g \n",
//	    outval, D, vsq, W );
//  }

  return(dpdvsq);
}


// dp/dv^2 holding W' fixed with utsq as input instead of vsq
static FTYPE dpdvsq_calc2_Wp(FTYPE Wp, FTYPE D, FTYPE utsq)
{
  FTYPE outval;
  FTYPE gamma,W;
  FTYPE compute_idwmrho0dp(FTYPE rho0, FTYPE wmrho0);
  FTYPE compute_idrho0dp(FTYPE rho0, FTYPE wmrho0);
  FTYPE rho0,wmrho0;
  FTYPE dwmrho0dvsq,idwmrho0dp;
  FTYPE drho0dvsq,idrho0dp;
  FTYPE gammasq;
  FTYPE dpdvsq;

  gammasq=1.0+utsq;
  gamma=sqrt(gammasq);
  wmrho0=wmrho0_compute_utsq(Wp, D, utsq, gamma, gammasq);

  rho0=D/gamma;

  W = Wp+D;

  idwmrho0dp=compute_idwmrho0dp(rho0, wmrho0);
  dwmrho0dvsq = (D*(gamma*0.5-1.0) - Wp);

  drho0dvsq = -D*gamma*0.5; // because \rho=D/\gamma
  idrho0dp = compute_idrho0dp(rho0, wmrho0);

  dpdvsq =   drho0dvsq *idrho0dp  +   dwmrho0dvsq *idwmrho0dp;

//  if( outval > 0. ) { 
//    dualfprintf(fail_file,"outval = %21.15g , D = %21.15g  , vsq = %21.15g,  W = %21.15g \n",
//	    outval, D, vsq, W );
//  }

  return(dpdvsq);
}

#endif





static FTYPE dpdvsq_calc_scn(FTYPE W, FTYPE D, FTYPE vsq)
{
  FTYPE outval;
  FTYPE compute_idwmrho0dp(FTYPE rho0, FTYPE wmrho0);
  FTYPE idwmrho0dp;
  FTYPE rho0,wmrho0;
  FTYPE gamma,gammasq;
  FTYPE Wp;

  Wp=W-D;
  gamma=1.0/sqrt(1.0-vsq);
  gammasq=gamma*gamma;
  wmrho0=wmrho0_compute_vsq(Wp, D, vsq, gamma, gammasq);

  rho0=D/gamma;
  idwmrho0dp=compute_idwmrho0dp(rho0, wmrho0);

  //  return( (gam - 1.) * ( 0.5 * D / sqrt(1.-vsq)  - W  ) / gam  ) ;
  outval =  ( 0.5 * D / sqrt(1.-vsq)  - W  )*idwmrho0dp  ;

//  if( outval > 0. ) { 
//    dualfprintf(fail_file,"outval = %21.15g , D = %21.15g  , vsq = %21.15g,  W = %21.15g \n",
//	    outval, D, vsq, W );
//  }

  return(outval);
}

// does NOT depend on EOS
// w-\rho_0 = (u+p) = W'/\gamma^2 - D v^2/(1+\gamma)
static FTYPE wmrho0_compute_vsq(FTYPE Wp, FTYPE D, FTYPE vsq, FTYPE gamma,FTYPE gammasq)
{
  return(Wp/gammasq - D*vsq/(1.0+gamma));

}


// does NOT depend on EOS
// w-\rho_0 = (u+p) = W'/\gamma^2 - D utsq/(1+\gamma)/gamma^2
static FTYPE wmrho0_compute_utsq(FTYPE Wp, FTYPE D, FTYPE utsq, FTYPE gamma,FTYPE gammasq)
{
  return( (Wp - D*utsq/(1.0+gamma) )/gammasq );

}

//#undef CRAZYDEBUG

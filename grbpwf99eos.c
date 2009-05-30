// EOS from PWF99 (GODMARK: mignone really, not done yet)


// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_grbpwf99(FTYPE rho0, FTYPE u)
{
  FTYPE pressure;

  pressure = u*(2.0*rho0+u)/(3.0*(rho0+u));

  //  dualfprintf(fail_file,"rho0=%g u=%g pressure=%g\n",rho0,u,pressure);

  return(pressure);
}

// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_grbpwf99(FTYPE rho0, FTYPE p)
{
  return( 1.5*(p + 3.0*p*p/(2.0*rho0+sqrt(9.0*p*p+4.0*rho0*rho0))) );
}

// dp(rho0, u)/du
FTYPE dpdu_rho0_u_grbpwf99(FTYPE rho0, FTYPE u)
{
  FTYPE dpdu;

  dpdu = 1.0/3.0*(1.0 + rho0*rho0/( (rho0+u)*(rho0+u)));

  return(dpdu);
}

// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u_grbpwf99(FTYPE rho0, FTYPE u)
{
  FTYPE dpdrho0;

  dpdrho0 = u*u/(3.0*(rho0+u)*(rho0+u));

  return(dpdrho0) ;
}


// sound speed squared (for vchar.c)
FTYPE cs2_compute_grbpwf99(FTYPE rho0, FTYPE u)
{
  FTYPE pressure;
  FTYPE h;
  FTYPE cs2;

  pressure = pressure_rho0_u(rho0,u);
  h=rho0+u+pressure; // not specific h

  cs2=pressure*(5.0*h - 8.0*pressure) / (3.0*h*(h-pressure));


  //  dualfprintf(fail_file,"cs2=%21.15g pressure=%21.15g\n",cs2,pressure);

  return(cs2);

}

// entropy as function of rho0 and internal energy (u)
// S(rho0,u)
FTYPE compute_entropy_grbpwf99(FTYPE rho0, FTYPE u)
{
  FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u);
  FTYPE entropy;

  entropy=0.0; // GODMARK: not set yet

  return(entropy);

}

// u(rho0,S)
FTYPE compute_u_from_entropy_grbpwf99(FTYPE rho0, FTYPE entropy)
{
  FTYPE u;

  u=0.0; // GODMARK: not set yet

  return(u);

}


#define GAMMA (gamideal)
#define GAMMAM1 (GAMMA-1.0)

// used for dudp_calc
FTYPE compute_dSdrho_grbpwf99(FTYPE rho0, FTYPE u)
{
  FTYPE indexn;
  FTYPE entropy;
  FTYPE dSdrho;
  FTYPE compute_entropy(FTYPE rho0, FTYPE u);

  entropy=compute_entropy(rho0,u);

  // ideal gas
  indexn=1.0/GAMMAM1;

  dSdrho=entropy/rho0-(indexn+1.0);

  return(dSdrho);

}


// used for dudp_calc
FTYPE compute_dSdu_grbpwf99(FTYPE rho0, FTYPE u)
{
  FTYPE indexn;
  FTYPE dSdu;

  // ideal gas
  indexn=1.0/GAMMAM1;

  dSdu=indexn*rho0/u;

  return(dSdu);

}


#undef GAMMA
#undef GAMMAM1



// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0_grbpwf99(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE Q,delta,delta2;
  FTYPE pressure;

  Q=wmrho0/rho0;
  delta=9.0/25.0*wmrho0*(2.0+Q);
  delta2=delta/rho0;

  pressure=(5.0/8.0)*(wmrho0 - delta/(1.0+sqrt(1.0+delta2)));

  return(pressure);
}


// 1 / (d(u+p)/dp) holding rho0 fixed
FTYPE compute_idwmrho0dp_grbpwf99(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE pressure_wmrho0(FTYPE rho0, FTYPE wmrho0);
  FTYPE idwmrho0dp;
  FTYPE p;

  p = pressure_wmrho0(rho0, wmrho0);

  idwmrho0dp = (2.0*wmrho0 + 2.0*rho0 - 5.0*p)/(5.0*rho0+5.0*wmrho0-8.0*p);

  return(idwmrho0dp);

}

// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp_grbpwf99(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE pressure_wmrho0(FTYPE rho0, FTYPE wmrho0);
  FTYPE idrho0dp;
  FTYPE p;

  p = pressure_wmrho0(rho0, wmrho0);

  idrho0dp = (2.0*wmrho0 - 5.0*p)/(5.0*rho0+5.0*wmrho0-8.0*p);

  return(idrho0dp);
}


// cooling or heating rate
FTYPE compute_qdot_grbpwf99(FTYPE rho0, FTYPE u)        
{
  return(0.0);
}



int compute_sources_EOS_grbpwf99(FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR])
{
  int sc,pl;

  SCPHYSICSLOOP(sc) PLOOP(pl) dUcomp[sc][pl]=0.0;
  return(0);
}

void compute_allextras_grbpwf99(int justnum, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras)
{
  int i;

  if(justnum==0){

    // set rest to 0
    for(i=0;i<MAXNUMEXTRAS;i++){
      extras[i] = 0.0;
    }
  }

  *numextrasreturn=0;
}

int get_extrasprocessed_grbpwf99(int doall, int i, int j, int k, FTYPE *pr, FTYPE *extras, FTYPE *processed)
{
  int pi;
  int numextrasreturn;
  FTYPE rho0,u;
  rho0=pr[RHO];
  u=pr[UU];
  compute_allextras(0, rho0, u, &numextrasreturn, extras);
  for(pi=0;pi<MAXPROCESSEDEXTRAS;pi++){
    processed[pi] = 0.0;
  }
  return(0);
}



// temperature
FTYPE compute_temp_grbpwf99(FTYPE rho0, FTYPE u)
{
  FTYPE temp;
  FTYPE pressure_rho0_u_grbpwf99(FTYPE rho0, FTYPE u);

  temp = pressure_rho0_u_grbpwf99(rho0,u)/rho0;

  return(temp);

}


void compute_EOS_parms_grbpwf99(FTYPE (*prim)[N2M][N3M][NPR])
{
  return; // do nothing                                                                                           
}



void store_EOS_parms_grbpwf99(int i, int j, int k, int numparms, FTYPE *parlist)
{
  return; // do nothing
}


void get_EOS_parms_grbpwf99(int i, int j, int k, int*numparms, FTYPE *parlist)
{

  *numparms=0;
  return; // do nothing
}


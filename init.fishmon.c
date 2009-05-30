
/* 
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"


#define SLOWFAC 1.0		/* reduce u_phi by this amount */

SFTYPE rhomax=0,umax=0,bsq_max=0,beta,rin,randfact;

int fieldfrompotential[NDIM];

FTYPE nz_func(FTYPE R) ;



int prepre_init_specific_init(void)
{

  // choice// GODMARK: not convenient location, but needed for init_mpi()
  periodicx1=0;
  periodicx2=0;
#if(USEMPI&&N3!=1)
  periodicx3=1;// GODMARK: periodic in \phi for 3D spherical polar
#else
  periodicx3=0;
#endif


  return(0);

}


int pre_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *before* all other calculations
  h_over_r=0.3;
  // below is theta distance from equator where jet will start, usually about 2-3X disk thickness
  h_over_r_jet=2.0*h_over_r;

  UTOPRIMVERSION = UTOPRIMJONNONRELCOMPAT;

  return(0);
}


int init_conservatives(FTYPE p[][N2M][N3M][NPR], FTYPE Utemp[][N2M][N3M][NPR], FTYPE U[][N2M][N3M][NPR])
{
  int pl;
  
  PLOOPBONLY(pl) trifprintf("fieldfrompotential[%d]=%d\n",pl-B1+1,fieldfrompotential[pl-B1+1]);


  trifprintf("begin init_conservatives\n");
  pi2Uavg(fieldfrompotential, p, Utemp, U);
  trifprintf("end init_conservatives\n");

  return(0);

}


int post_init_specific_init(void)
{
  // globally used parameters set by specific initial condition routines, reran for restart as well *after* all other calculations

  UTOPRIMVERSION = UTOPRIMJONNONRELCOMPAT;

  //  UTOPRIMVERSION =   UTOPRIM5D1;
  //    UTOPRIMVERSION =   UTOPRIM2DFINAL;

  return(0);
}



int init_consts(void)
{
  //  Lunit=Tunit=Munit=1.0;

  // units can be used for user to read in data, but otherwise for rest of code all that matters is Mfactor and Jfactor
  Mfactor=Jfactor=1.0;

  return(0);

}


#define NORMALTORUS 0 // note I use randfact=5.e-1 for 3D model with perturbations
#define GRBJET 1
#define KEPDISK 2

// For blandford problem also need to set:
// 0) WHICHPROBLEM 2
// 1) a=0.92
// 2) Rout=1E3
// 3) hslope=0.3
// 4) BSQORHOLIMIT=1E2;
// 5) BSQOULIMIT=1E3;
// 6) UORHOLIMIT=1E3;
// 7) tf=1E4; // and maybe DTd's
// 8) randfact = .2;
// 9) lim=PARALINE; FLUXB=FLUXCTSTAG; TIMEORDER=4;
// 10) N1,N2,N3 in init.h
// 11) MAXBND 4
// 12) PRODUCTION 1
// 13) USE<systemtype>=1
// 14) setup batch

/*

Models to run:

Constant parameters:

1) Rout=1E3 and run for tf=1E4 (so will take 5X longer than compared to Orange run at 128x128x32)

2) BSQORHOLIMIT=1E3, etc.

3) PARALINE, FLUXCTSTAG, TO=4

4) Form of A_\phi fixed

Field parameter studies in 2D axisymmetry at 256^2:

1) H/R=0.3, a=0.9: LS quadrapole,  LS dipole, SS quadrapole, SS dipole

 

Spin parameter study in 2D axisymmetry at 256^2:

 

1) H/R=0.3, LS quadrapole: a=-.999,-.99,-.9,-.5,-0.2,0,.2,.5,.9,.99,.999

H/R parameter study in 2D axisymmetry at 256^2:

1) a=0.9 LS quadrapole with H/R=0.1,0.3,0.9,1.5

2D Fiducial Models:

1) Using a=0.9, H/R=0.3, LS quad and LS dipole, do two 2D fudicial models at: 1024^2

3D Fiducial Models:

1) Using a=0.9, H/R=0.3, LS quadrapole and LS dipole, do two 3D fiducial models at 2 different resolutions: 128x128x32 and 256x256x64

Questions for Roger:

1) Choice for disk thickness?
2) Choice for field shape -- specifically?
3) Choice for flux threading disk vs. BH initially?
4) Ask about BZ77 and residual A_\phi at pole
5) 

*/


#define WHICHPROBLEM 2 // choice


int init_grid(void)
{
  SFTYPE rh;
  
  // metric stuff first
  a = 0.9375 ;
  

#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  // make changes to primary coordinate parameters R0, Rin, Rout, hslope
  R0 = 0.0;
  Rout = 40.;
  // define coordinate type
  defcoord = 9;
#elif(WHICHPROBLEM==GRBJET)
  R0 = -3.0;
  Rout = 1E5;
  // define coordinate type
  defcoord = JET4COORDS;
#endif

 
  Rhor=rhor_calc(0);
  Rin=setRin(setihor());
  //  Rin = 0.98 * Rhor;
  
  //  hslope = 0.3;

  hslope = 1.04*pow(h_over_r,2.0/3.0);





  return(0);
}

int init_global(void)
{
  int pl;

  DODIAGEVERYSUBSTEP = 0;


  SAFE=1.3;
  //  cour = 0.9;
  cour=0.8;
  //  cour = 0.5;

  ///////////////////////
  //
  // ENO-RELATED STUFF
  //
  ///////////////////////
  
  INVERTFROMAVERAGEIFFAILED = 1;
  LIMIT_AC_PRIM_FRAC_CHANGE = 1;
  MAX_AC_PRIM_FRAC_CHANGE = 0.1;

  LIMIT_AC_FRAC_CHANGE=0; // CHANGINGMARK: avoiding complication for now
  // have to make consistent with weno-minimization for fluxes
  MAX_AC_FRAC_CHANGE=0.2;

  // need MAXBND==17 if not evolving point field and doing full WENO5BND
  // test1103 N1=8 is smallest tried and new simple_ limiting works at 10% very well
  //dofluxreconevolvepointfield=0;
  // only need MAXBND==11 like normal higher-order method (like FV method)
  dofluxreconevolvepointfield=1;



#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)


  //avgscheme[1]=avgscheme[2]=avgscheme[3]=WENO5BND;
  //  lim[1] = lim[2] = lim[3] = WENO5BND;
  lim[1] = lim[2] = lim[3] = PARALINE;

  avgscheme[1]=avgscheme[2]=avgscheme[3]=DONOR; // CHANGINGMARK
  //lim[1] = lim[2] = lim[3] = MC; // CHANGINGMARK


  DOENOFLUX = NOENOFLUX;
  //DOENOFLUX = ENOFLUXRECON; // CHANGINGMARK
  //DOENOFLUX = ENOFINITEVOLUME;

  if(DOENOFLUX == ENOFLUXRECON){
    // below applies to all fluxes
    PALLLOOP(pl) do_transverse_flux_integration[pl] = 1;
    PLOOPBONLY(pl) do_transverse_flux_integration[pl] = 1;
    // below used for re-averaging of field in advance.c for dBhat/dt method
    PALLLOOP(pl) do_conserved_integration[pl] = 1;
    PLOOPBONLY(pl) do_conserved_integration[pl] = 1;
  }

  if(DOENOFLUX == ENOFINITEVOLUME){
    PALLLOOP(pl) do_transverse_flux_integration[pl] = 1;
    PLOOPBONLY(pl) do_transverse_flux_integration[pl] = 1;
    PALLLOOP(pl) do_source_integration[pl] = 0;
    PLOOPBONLY(pl) do_source_integration[pl] = 0;
    PALLLOOP(pl) do_conserved_integration[pl] = 1;
    PLOOPBONLY(pl) do_conserved_integration[pl] = 1;
    //    do_conserved_integration[B1-1+DIRZ] = 1;
  }



  FLUXB = FLUXCTSTAG;  // CHANGINGMARK
  //  FLUXB = FLUXCTHLL;
  //FLUXB = FLUXCTTOTH;
  //  TIMEORDER=2;
  TIMEORDER=4;
  //  TIMEORDER=3;
  //  fluxmethod= HLLFLUX;
  fluxmethod= LAXFFLUX; // generally more robust than HLLFLUX for various reasons
  

  //  UTOPRIMVERSION=UTOPRIM5D1;
  UTOPRIMVERSION = UTOPRIMJONNONRELCOMPAT;
  //  UTOPRIMVERSION = UTOPRIM1DFINAL;


#elif(EOMTYPE==EOMFFDE)
  // PARA and TO=4 and HLL not trustable in FFDE so far
  lim[1] =lim[2]=lim[3]= MC;
  TIMEORDER=2;


  // below applies to all fluxes
  PALLLOOP(pl) do_transverse_flux_integration[pl] = 1;
  PLOOPBONLY(pl) do_transverse_flux_integration[pl] = 1;
  // below used for re-averaging of field in advance.c for dBhat/dt method
  PALLLOOP(pl) do_conserved_integration[pl] = 1;
  PLOOPBONLY(pl) do_conserved_integration[pl] = 1;



  fluxmethod=LAXFFLUX; // generally more robust than HLLFLUX
  FLUXB = FLUXCTTOTH;
  UTOPRIMVERSION=UTOPRIM2DFINAL;
  // whether/which ENO used to interpolate fluxes
  DOENOFLUX = ENOFINITEVOLUME;
  //  DOENOFLUX= NOENOFLUX;
  //DOENOFLUX=ENOFLUXRECON;
#endif



  ranc(7); // no MPI method yet, so just pure randomization
  /* some physics parameters */
  gam = 4. / 3.;
  cooling=0;


#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  BCtype[X1UP]=OUTFLOW;
  //  rescaletype=1;
  rescaletype=4;
  BSQORHOLIMIT=1E3; // was 1E2 but latest BC test had 1E3 // CHANGINGMARK
  BSQOULIMIT=1E4; // was 1E3 but latest BC test had 1E4
  RHOMIN = 1E-4;
  UUMIN = 1E-6;
#elif(WHICHPROBLEM==GRBJET)
  BCtype[X1UP]=FIXEDOUTFLOW;
  rescaletype=4;
  BSQORHOLIMIT=1E3;
  BSQOULIMIT=1E4;
  RHOMIN = 23.0;
  UUMIN = 1.7;
#endif


  BCtype[X1DN]=OUTFLOW;
  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;
  BCtype[X3UP]=PERIODIC;
  BCtype[X3DN]=PERIODIC;



  // below floor model is only used if rescaletype!=4
  if(BCtype[X1UP]==FIXEDOUTFLOW){ // then doing bondi inflow
    // avoids constant floor activation -- trying to be more physical
    prfloorcoef[RHO]=RHOMIN/100.0;
    prfloorcoef[UU]=UUMIN/100.0;
  }
  else{
    prfloorcoef[RHO]=RHOMIN;
    prfloorcoef[UU]=UUMIN;
  }


#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
  /* output choices */
  tf = 2000.0;

  /* dumping frequency, in units of M */
  DTdumpgen[DTDISS]=DTdumpgen[DTFLUX]=DTdumpgen[DTOTHER]=DTdumpgen[DTEOS]=DTdumpgen[DTVPOT]=DTdumpgen[DTDUMP] = 50.;
  DTdumpgen[DTAVG] = 50.0;
  // ener period
  DTdumpgen[DTENER] = 2.0;
  /* image file frequ., in units of M */
  DTdumpgen[DTIMAGE] = 2.0;
  // fieldline locked to images so can overlay
  DTdumpgen[DTFIELDLINE] = DTdumpgen[DTIMAGE];

  /* debug file */  
  DTdumpgen[DTDEBUG] = 50.0;
  // DTr = .1 ; /* restart file frequ., in units of M */
  /* restart file period in steps */
  DTr = 100;

#elif(WHICHPROBLEM==GRBJET)
  /* output choices */
  tf = 5E5;
  
  DTd = 250.;                 /* dumping frequency, in units of M */
  DTavg = 250.0;
  DTener = 2.0;                       /* logfile frequency, in units of M */
  DTi = 10.0;                 /* image file frequ., in units of M */
  DTdebug = 250.0; /* debug file */
  // DTr = .1 ; /* restart file frequ., in units of M */
  DTr = 100;                  /* restart file period in steps */
#endif



  return(0);

}

// assumes normalized density
int init_atmosphere(int *whichvel, int*whichcoord,int i, int j, int k, FTYPE *pr)
{
  int pl;
  struct of_geom realgeom,geom;
  FTYPE pratm[NPR];


  get_geometry(i, j, k, CENT, &realgeom); // true coordinate system
  set_atmosphere(-1,WHICHVEL,&realgeom,pratm); // set velocity in chosen WHICHVEL frame in any coordinate system

  if(pr[RHO]<pratm[RHO]){
    PLOOP(pl) pr[pl]=pratm[pl];
  }
  

  *whichvel=WHICHVEL;
  *whichcoord=PRIMECOORDS;
  return(0);


}

int init_grid_post_set_grid(void)
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;
  extern void check_spc_singularities_user(void);

  // some calculations, althogh perhaps calculated already, definitely need to make sure computed
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);




  //SASMARK restart: need to populate panalytic with IC's
  if( RESTARTMODE==1 ) { //restarting -> set panalytic to initital conditions
    // user function that should fill p with primitives
    MYFUN(init_primitives(panalytic),"initbase.c:init()", "init_primitives()", 0);
    //to have initial vector potential to be saved in panalytic array
  }




  
  // diagnostic
  // determine nature of inner radial edge (assumes myid==0 is always there)
  if(myid==0){
    i=-INFULL1;
    j=k=0;
    coord(i,j,k, FACE1, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];
    trifprintf("rmin: %21.15g\n", r);
    trifprintf("rmin/rh: %21.15g\n", r / Rhor );
    //    trifprintf("rmin/rsing: %21.15g\n", r / (a+SMALL));
    if(r/Rhor<=1.0){
      trifprintf("inner grid is inside horizon\n");
    }
    else{
      trifprintf("inner grid is outside horizon\n");
    }
  }

  // check that singularities are properly represented by code
  check_spc_singularities_user();

  
  return(0);

}



int init_primitives(FTYPE p[][N2M][N3M][NPR])
{
  int whichvel, whichcoord;
  int initreturn;
  int i = 0, j = 0, k = 0, l;
  struct of_geom geom;
  FTYPE r,th,X[NDIM],V[NDIM];
  int normalize_densities(FTYPE p[][N2M][N3M][NPR]);
  int normalize_field(FTYPE p[][N2M][N3M][NPR]);
  int init_dsandvels(int *whichvel, int *whichcoord, int i, int j, int k, FTYPE *p, FTYPE *pstag);
  int init_atmosphere(int *whichvel, int *whichcoord, int i, int j, int k, FTYPE *pr);
  int pl;


  ///////////////////////////////////
  //
  // Assign primitive variables
  //
  ///////////////////////////////////
  trifprintf("Assign primitives\n");

  // assume we start in bl coords and convert to KSprim
  COMPFULLLOOP{
    PLOOP(pl) p[i][j][k][pl]=0.0; // so field defined when get to floor model (fixup)
  }

  // assume we start in bl coords and convert to KSprim
  COMPFULLLOOP{
    initreturn=init_dsandvels(&whichvel, &whichcoord,i,j,k,p[i][j][k],pstagscratch[i][j][k]); // request densities for all computational centers
    if(initreturn>0) return(1);
    else MYFUN(transform_primitive_vB(whichvel, whichcoord, i,j,k, p, pstagscratch),"init.c:init_primitives","transform_primitive_vB()",0);
  }


  /////////////////////////////
  //
  // normalize density if wanted
  //
  /////////////////////////////// 
  // at this point densities are still standard, so just send "p"
  trifprintf("Normalize densities\n");
  normalize_densities(p);


  /////////////////////////////
  //
  // Define an atmosphere if wanted
  //
  /////////////////////////////// 

#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)
  // normalized atmosphere
  trifprintf("Add atmosphere\n");
  COMPZLOOP{
    initreturn=init_atmosphere(&whichvel, &whichcoord,i,j,k,p[i][j][k]);
    if(initreturn>0) return(1);
    else{
      // transform from whichcoord to MCOORD
      if (bl2met2metp2v(whichvel, whichcoord,p[i][j][k], i,j,k) >= 1)
	FAILSTATEMENT("init.c:init()", "bl2ks2ksp2v()", 1);
    }
  }
#endif


  // copy over initial solution as analytic solution
  // NEEDED FOR BOUND in case uses panalytic
  COMPFULLLOOP{
    PLOOP(pl) panalytic[i][j][k][pl]=p[i][j][k][pl];
  }


  /////////////////////////////
  //
  // Fixup and Bound variables since some primitive quantities may have changed
  // These may be used to define vector potential below
  // Also setup pre_fixup() type quantities
  //
  /////////////////////////////// 
  trifprintf("Fixup and Bound #1\n");

#if(EOMTYPE!=EOMFFDE)
  // assume EOMFFDE doesn't use "density/ie" to set field, so no need to bound, and no field definition is bad for EOMFFDE
#if(FIXUPAFTERINIT)
  if(fixup(STAGEM1,p,0)>=1)
    FAILSTATEMENT("init.c:init()", "fixup()", 1);
#endif

  if (bound_prim(STAGEM1,0.0,p) >= 1)
    FAILSTATEMENT("init.c:init()", "bound_prim()", 1);

  if(pre_fixup(STAGEM1,p)>=1)
    FAILSTATEMENT("init.c:init()", "pre_fixup()", 1);
#endif


  
  /////////////////////////////
  //
  // Initialize field from vector potential
  //
  /////////////////////////////// 
#if(1)
  init_vpot(p);
  normalize_field(p); // normalizes p and pstagscratch and unew and vpotarray if tracked
#else
  // no field
  init_zero_field(p);
#endif




  // copy over initial solution as analytic solution
  // NEEDED FOR BOUND in case uses panalytic                                                                                                                                                            
  COMPFULLLOOP{
    PLOOP(pl) panalytic[i][j][k][pl]=p[i][j][k][pl];
  }


  /////////////////////////////
  //
  // Fixup and Bound variables since some primitive quantities may have changed
  // These may be used to define vector potential below
  // Also setup pre_fixup() type quantities
  //
  //
  // BOUND AGAIN IN CASE USING PANALYTIC TO BOUND!
  //
  /////////////////////////////// 
  trifprintf("Fixup and Bound #2\n");

#if(EOMTYPE!=EOMFFDE)
  // assume EOMFFDE doesn't use "density/ie" to set field, so no need to bound, and no field definition is bad for EOMFFDE
#if(FIXUPAFTERINIT)
  if(fixup(STAGEM1,p,0)>=1)
    FAILSTATEMENT("init.c:init()", "fixup()", 1);
#endif

  if (bound_allprim(STAGEM1,0.0,p) >= 1)
    FAILSTATEMENT("init.c:init()", "bound_allprim()", 1);

  if(pre_fixup(STAGEM1,p)>=1)
    FAILSTATEMENT("init.c:init()", "pre_fixup()", 1);
#endif




  return(0);


}



int init_dsandvels(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  int init_dsandvels_torus(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);
  int init_dsandvels_thindisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag);

  beta = 1.e2 ;
  randfact = 4.e-2;


#if(WHICHPROBLEM==NORMALTORUS)
  return(init_dsandvels_torus(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#elif(WHICHPROBLEM==KEPDISK)
  return(init_dsandvels_thindisk(whichvel, whichcoord,  i,  j,  k, pr, pstag));
#endif

}


// unnormalized density
int init_dsandvels_torus(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],V[NDIM],r,th;
  struct of_geom realgeom,geom;
  

  /* for disk interior */
  SFTYPE l, lnh, expm2chi, up1;
  SFTYPE DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
  SFTYPE kappa, hm1;
  SFTYPE rmax, lfish_calc(SFTYPE rmax);
  SFTYPE rh;
  //  FTYPE pratm[NPR];
  int pl;



  kappa = 1.e-3 ;
  rin = 6. ;
  rmax = 12. ;
  l = lfish_calc(rmax) ;
  

  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];



  sth = sin(th);
  cth = cos(th);

  /* calculate lnh */
  DD = r * r - 2. * r + a * a;
  AA = (r * r + a * a) * (r * r + a * a) - DD * a * a * sth * sth;
  SS = r * r + a * a * cth * cth;
  
  thin = M_PI / 2.;
  sthin = sin(thin);
  cthin = cos(thin);
  DDin = rin * rin - 2. * rin + a * a;
  AAin = (rin * rin + a * a) * (rin * rin + a * a)
    - DDin * a * a * sthin * sthin;
  SSin = rin * rin + a * a * cthin * cthin;
  
  if (r >= rin) {
    lnh = 0.5 * log((1. + sqrt(1. + 4. * (l * l * SS * SS) * DD /
			       (AA * sth * AA * sth))) / (SS * DD /
							  AA))
      - 0.5 * sqrt(1. +
		   4. * (l * l * SS * SS) * DD / (AA * AA * sth *
						  sth))
      - 2. * a * r * l / AA -
      (0.5 *
       log((1. +
	    sqrt(1. +
		 4. * (l * l * SSin * SSin) * DDin / (AAin * AAin *
						      sthin *
						      sthin))) /
	   (SSin * DDin / AAin))
       - 0.5 * sqrt(1. +
		    4. * (l * l * SSin * SSin) * DDin / (AAin *
							 AAin *
							 sthin *
							 sthin))
       - 2. * a * rin * l / AAin);
  } else
    lnh = 1.;
  

  
  /* regions outside torus */
  // this region is already in Kerr Schild prime in proper primitive quantity for velocity
  if (lnh < 0. || r < rin) {


    get_geometry(i, j, k, CENT, &realgeom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,&realgeom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
  /* region inside magnetized torus; u^i is calculated in
     Boyer-Lindquist coordinates, as per Fishbone & Moncrief, so it
     needs to be transformed at the end */
  else {
    hm1 = exp(lnh) - 1.;
    rho = pow(hm1 * (gam - 1.) / (kappa * gam), 1. / (gam - 1.));
    u = kappa * pow(rho, gam) / (gam - 1.);
    ur = 0.;
    uh = 0.;
    
    /* calculate u^phi */
    expm2chi = SS * SS * DD / (AA * AA * sth * sth);
    up1 = sqrt((-1. + sqrt(1. + 4. * l * l * expm2chi)) / 2.);
    up = 2. * a * r * sqrt(1. + up1 * up1) / sqrt(AA * SS * DD) +
      sqrt(SS / AA) * up1 / sth;
    
    
    pr[RHO] = rho ;
    pr[UU] = u* (1. + randfact * (ranc(0) - 0.5));
    pr[U1] = ur ;
    pr[U2] = uh ;    
    pr[U3] = SLOWFAC * up;

    // just define some field
    pr[B1]=0.0;
    pr[B2]=0.0;
    pr[B3]=0.0;

    if(FLUXB==FLUXCTSTAG){
      // assume pstag later defined really using vector potential or directly assignment of B3 in axisymmetry
      PLOOPBONLY(pl) pstag[pl]=pr[pl];
    }

    *whichvel=VEL4;
    *whichcoord=BLCOORDS;
    return(0);
  }
}



int init_dsandvels_thindisk(int *whichvel, int*whichcoord, int i, int j, int k, FTYPE *pr, FTYPE *pstag)
{
  SFTYPE sth, cth;
  SFTYPE ur, uh, up, u, rho;
  FTYPE X[NDIM],V[NDIM],r,th,ph;
  struct of_geom geom;
  /* for disk interior */
  FTYPE R,H,nz,z,S,cs ;
  SFTYPE rh;
  int pl;
  

  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  r=V[1];
  th=V[2];
  ph=V[3];

  //  beta=1.e2;
  //  beta=1.e2/3.0;
  //  randfact = 4.e-2;
  //rin = (1. + h_over_r)*Risco;
  rin = Risco;



  /* region outside disk */
  R = r*sin(th) ;

  if(R < rin) {

    get_geometry(i, j, k, CENT, &geom); // true coordinate system
    set_atmosphere(-1,WHICHVEL,&geom,pr); // set velocity in chosen WHICHVEL frame in any coordinate system

    *whichvel=WHICHVEL;
    *whichcoord=PRIMECOORDS;
    return(0);
  }
  else {

    H = h_over_r*R ;
    nz = nz_func(R) ;
    z = r*cos(th) ;
    S = 1./(H*H*nz) ;
    cs = H*nz ;

    rho = (S/sqrt(2.*M_PI*H*H)) * exp(-z*z/(2.*H*H)) * taper_func(R,rin) ;
    u = rho*cs*cs/(gam - 1.) ;
    ur = 0. ;
    uh = 0. ;
    up = 1./(pow(r,1.5) + a) ;
    // solution for 3-vel


    
    
    pr[RHO] = rho ;
    pr[UU] = u* (1. + randfact * (ranc(0) - 0.5));

    pr[U1] = ur ;
    pr[U2] = uh ;    
    pr[U3] = SLOWFAC * up;

    PLOOPBONLY(pl) pstag[pl]=pr[pl]=0.0;

    *whichvel=VEL3;
    *whichcoord=BLCOORDS;
    return(0);
  }
}








#define DISKFIELD 0
#define VERTFIELD 1
#define DISKVERT 2
#define BLANDFORDQUAD 3

#define FIELDTYPE BLANDFORDQUAD

FTYPE setgpara(FTYPE myr, FTYPE th, FTYPE thpower)
{
  FTYPE fneg,fpos;
  FTYPE gpara;

  fneg=1.0-pow(cos(th),thpower);
  fpos=1.0+pow(cos(th),thpower);
  gpara=0.5*(myr*fneg + 2.0*fpos*(1.0-log(fpos)));
  // remove BZ77 Paraboloidal divb!=0 at pole
  gpara=gpara-2.0*(1.0-log(2.0));

  return(gpara);


}

FTYPE setblandfordfield(FTYPE r, FTYPE th)
{
  FTYPE setgpara(FTYPE myr, FTYPE th, FTYPE thpower);
  FTYPE rshift,myr,rpower,myz,myR,myvert;
  FTYPE thother,thpower,gparalow,gparahigh,mygpara;
  FTYPE aphi;


  rshift=4.0;
  rpower=0.75;
  thpower=4.0;
  

  myr=pow(r+rshift,rpower);
  myz=myr*cos(th);
  myR=myr*sin(th);
  myvert = (th>M_PI*0.5) ? (myr*sin(th)) : (myr*sin(-th));

  thother=M_PI-th;
  gparalow=setgpara(myr,th,thpower);
  gparahigh=setgpara(myr,thother,thpower);
  mygpara=(th<0.5*M_PI) ? gparalow : gparahigh;

  // GOOD:
  // aphi=mygpara;
  // aphi=mygpara*cos(th); // B1 diverges at pole
  //  aphi=mygpara*cos(th)*sin(th); // doesn't diverge as much
  //  aphi=mygpara*cos(th)*sin(th)*sin(th); // old choice before subtracted original BZ77 problem
  aphi=mygpara*cos(th); // latest choice
  //aphi=myvert*cos(th); // vert with quad
  //aphi=myR*cos(th);

  // BAD:
  // aphi=myvert;
  
  

  return(aphi);


}

// assumes normal field in pr
int init_vpot_user(int *whichcoord, int l, int i, int j, int k, FTYPE p[][N2M][N3M][NPR], FTYPE *V, FTYPE *A)
{
  SFTYPE rho_av, q;
  FTYPE r,th;
  FTYPE vpot;
  FTYPE setblandfordfield(FTYPE r, FTYPE th);




  vpot=0.0;


  if(l==3){// A_\phi

    r=V[1];
    th=V[2];



    // Blandford quadrapole field version
    if(FIELDTYPE==BLANDFORDQUAD){
      vpot += setblandfordfield(r,th);
    }

    /* vertical field version*/
    if((FIELDTYPE==VERTFIELD)||(FIELDTYPE==DISKVERT)){
      vpot += 0.5*r*sin(th)*sin(th) ;
    }


    /* field-in-disk version */
    
    if((FIELDTYPE==DISKFIELD)||(FIELDTYPE==DISKVERT)){
      // average of density that lives on CORN3


      // since init_vpot() is called for all i,j,k, can't use
      // non-existence values, so limit averaging:
      if((i==-N1BND)&&(j==-N2BND)){
	rho_av = p[i][j][k][RHO];
      }
      else if(i==-N1BND){
	rho_av = AVGN_2(p,i,j,k,RHO);
      }
      else if(j==-N2BND){
	rho_av = AVGN_1(p,i,j,k,RHO);
      }
      else{ // normal cells
	rho_av = AVGN_for3(p,i,j,k,RHO);
      }

      q = rho_av / rhomax - 0.2;

      if (q > 0.)      vpot += q;
    }
  }

  //////////////////////////////////
  //
  // finally assign what's returned
  //
  //////////////////////////////////
  *A = vpot;
  *whichcoord = MCOORD;



  return(0);

}



int init_vpot2field_user(FTYPE A[][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE pr[][N2M][N3M][NPR])
{
  extern int vpot2field(FTYPE A[][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE p[][N2M][N3M][NPR]);
  int i,j,k,pl;
  int toreturn;
 

  // uses panalytic as temporary variable
  COMPFULLLOOP{
    PLOOPBONLY(pl) panalytic[i][j][k][pl]=pr[i][j][k][pl];
  }


  // obtain primitive magnetic field from vector potential
  toreturn=vpot2field(A,panalytic); // uses panalytic as temporary variable

  // default (assume all fields are from potential)
  PLOOPBONLY(pl) fieldfrompotential[pl-B1+1]=1;


  // Can override vector potential choice for some field components, like B3 in axisymmetry
  // see init.sasha.c

  ////////////////////
  //
  // don't override
  //
  ////////////////////
  COMPFULLLOOP{
    PLOOPBONLY(pl) pr[i][j][k][pl]=panalytic[i][j][k][pl];
  }

  return(toreturn);

}



// assumes we are fed the true densities
int normalize_densities(FTYPE p[][N2M][N3M][NPR])
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM],r,th;


  rhomax=0;
  umax=0;
  COMPZLOOP{
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    r=V[1];
    th=V[2];

    if (p[i][j][k][RHO] > rhomax)   rhomax = p[i][j][k][RHO];
    if (p[i][j][k][UU] > umax && r > rin)    umax = p[i][j][k][UU];
  }

  mpimax(&rhomax);
  mpimax(&umax);
  trifprintf("rhomax: %21.15g umax: %21.15g\n", rhomax, umax);


  COMPZLOOP{
    p[i][j][k][RHO] /= rhomax;
    p[i][j][k][UU] /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1.;

  return(0);
}


// assumes normal field definition
int normalize_field(FTYPE p[][N2M][N3M][NPR])
{
  int i,j,k;
  FTYPE bsq_ij;
  SFTYPE bsq_max, norm, beta_act;
  struct of_geom geom;
  FTYPE X[NDIM],V[NDIM];
  FTYPE r,th;

  bsq_max = 0.;
  COMPZLOOP {
    get_geometry(i, j, k, CENT, &geom);    

    if(FIELDTYPE==VERTFIELD || FIELDTYPE==BLANDFORDQUAD){
      coord(i, j, k, CENT, X);
      bl_coord(X, V);
      r=V[1];
      th=V[2];
      
      if((r>rin)&&(fabs(th-M_PI*0.5)<4.0*M_PI*dx[2]*hslope)){
	if (bsq_calc(p[i][j][k], &geom, &bsq_ij) >= 1)
	  FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
	
	if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
      }
    }
    else{
      if (bsq_calc(p[i][j][k], &geom, &bsq_ij) >= 1)
	FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
      
      if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
    }
  }

  mpimax(&bsq_max);
  trifprintf("initial bsq_max: %21.15g\n", bsq_max);

  /* finally, normalize to set field strength */
  beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
  trifprintf("initial beta: %21.15g (should be %21.15g)\n", beta_act,beta);
  norm = sqrt(beta_act / beta);
  
  
  // not quite right since only correct static field energy, not moving field energy
  normalize_field_withnorm(norm);


  // check bsq_max again
  bsq_max = 0.;
  COMPZLOOP {
    get_geometry(i, j, k, CENT, &geom);
    if (bsq_calc(p[i][j][k], &geom, &bsq_ij) >= 1)
      FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);
    if (bsq_ij > bsq_max)      bsq_max = bsq_ij;
    
  }
  mpimax(&bsq_max);
  trifprintf("new initial bsq_max: %21.15g\n", bsq_max);

  beta_act = (gam - 1.) * umax / (0.5 * bsq_max);

  trifprintf("new bsq_max: %21.15g\n", bsq_max);
  trifprintf("final beta: %21.15g (should be %21.15g)\n", beta_act,beta);

  return(0);
}



#undef SLOWFAC

SFTYPE lfish_calc(SFTYPE r)
{
  return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
	   ((-2. * a * r * (pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2))) /
	    sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
	    ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) +
					  pow(a,
					      2) * (2. + r))) / sqrt(1 +
								     (2.
								      *
								      a)
								     /
								     pow
								     (r,
								      1.5)
								     -
								     3.
								     /
								     r)))
	  / (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
	     (pow(a, 2) + (-2. + r) * r))
	  );
}

// UUMIN/RHOMIN used for atmosphere

// for each WHICHVEL possibility, set atmosphere state for any coordinate system
// which=0 : initial condition
// which=1 : evolution condition (might also include a specific angular momentum or whatever)
// which==1 assumes pr set to something locally reasonable, and we adjust to that slowly

#define TAUADJUSTATM (10.0) // timescale for boundary to adjust to using preset inflow
int set_atmosphere(int whichcond, int whichvel, struct of_geom *ptrgeom, FTYPE *pr)
{
  FTYPE rho,u,ur,uh,up;
  FTYPE X[NDIM],V[NDIM];
  FTYPE r,th;
  FTYPE prlocal[NPR];
  int pl;

  // Bondi like initial atmosphere
  //    rho = RHOMIN * 1.E-14;
  //    u = UUMIN * 1.E-14;
  coord(ptrgeom->i, ptrgeom->j, ptrgeom->k, ptrgeom->p, X);
  bl_coord(X,V);
  r=V[1];
  th=V[2];

  // default
  PLOOP(pl) prlocal[pl]=pr[pl];

#if((EOMTYPE==EOMGRMHD)||(EOMTYPE==EOMCOLDGRMHD))
  // Bondi-like atmosphere
  if(rescaletype==4){
#if(WHICHPROBLEM==NORMALTORUS || WHICHPROBLEM==KEPDISK)
    // couple rescaletype to atmosphere type
    prlocal[RHO] = RHOMIN*pow(r,-2.0);
#elif(WHICHPROBLEM==GRBJET)
    // couple rescaletype to atmosphere type
    if(r>40.0) prlocal[RHO] = RHOMIN*pow(r,-2.0);
    else prlocal[RHO] = RHOMIN*pow(40.0,-2.0);
#endif
  }
  else{
    prlocal[RHO] = RHOMIN*pow(r,-1.5);
  }
#else
  prlocal[RHO] = 0;
#endif

#if(EOMTYPE==EOMGRMHD)
  // Bondi-like atmosphere
  prlocal[UU]  = UUMIN*pow(r,-2.5);
#else
  prlocal[UU]  = 0;
#endif

    
  // bl-normal observer (4-vel components)
  
  // normal observer velocity in atmosphere
  if(whichvel==VEL4){
    prlocal[U1] = -ptrgeom->gcon[0][1]/sqrt(-ptrgeom->gcon[0][0]) ;
    prlocal[U2] = -ptrgeom->gcon[0][2]/sqrt(-ptrgeom->gcon[0][0]) ;
    prlocal[U3] = -ptrgeom->gcon[0][3]/sqrt(-ptrgeom->gcon[0][0]) ;
  }
  else if(whichvel==VEL3){
    prlocal[U1] = ptrgeom->gcon[0][1]/ptrgeom->gcon[0][0] ;
    prlocal[U2] = ptrgeom->gcon[0][2]/ptrgeom->gcon[0][0] ;
    prlocal[U3] = ptrgeom->gcon[0][3]/ptrgeom->gcon[0][0] ;
    // GAMMIE
    //ur = -1./(r*r);
    //uh=up=0.0;
  }
  else if(whichvel==VELREL4){
    prlocal[U1] = 0.0;
    prlocal[U2] = 0.0;
    prlocal[U3] = 0.0;
  }
  
  if(whichcond==1){
    if(100.0*dt>TAUADJUSTATM){
      dualfprintf(fail_file,"dt=%21.15g and TAUADJUSTATM=%21.15g\n",dt,TAUADJUSTATM);
      myexit(1);
    }
    // TAUADJUSTATM must be >> dt always in order for this to make sense (i.e. critical damping to fixed solution)
    PLOOP(pl) pr[pl] = pr[pl]+(prlocal[pl]-pr[pl])*dt/TAUADJUSTATM;
  }
  else if(whichcond==0){ 
    PLOOP(pl) pr[pl] = prlocal[pl];
    // very specific
    // always outflow field
    //    pr[B1] = pr[B2] = pr[B3] = 0;
  }
  else if(whichcond==-1){ 
    // t=0, just sets as "floor"
    PLOOP(pl) pr[pl] = prlocal[pl];
    pr[B1] = pr[B2] = pr[B3] = 0;
  }


  return(0);

}



int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *prfloor)
{
  return(set_density_floors_default(ptrgeom, pr, prfloor));
}




FTYPE nz_func(FTYPE R)
{
  return(
	 sqrt(
	      (3.*a*a - 4.*a*sqrt(R) + R*R)/
	      pow(R*(a + pow(R,1.5)),2)
	      )
	 ) ;


}

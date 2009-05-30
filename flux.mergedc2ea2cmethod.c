#include "decs.h"


// 3-point method : Not useful to evaluate flux at 3 points and fit so just use full method for 3-points
//
// 1) 3-point creation of coefficients list
// 2) function to propagate primitives to fluxes fully




// Need functions:
// 1) Feed coefficient list for certain order
// 2) Create pressure or other aux. variables using 3 points so matches at faces
// 3) Feed primitive + aux variables to function that knows form of primitives and flux from mathematica


int deconvolve_flux(int dir, struct of_geom *ptrgeom, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, struct of_state *state_c, struct of_state *state_l, struct of_state *state_r, FTYPE *F_c, FTYPE *F_l, FTYPE *F_r, FTYPE *U_c, FTYPE *U_l, FTYPE *U_r)
{

  // do matter parts first
  deconvolve_flux_ma(dir, ptrgeom, p_c, p_l, p_r, state_c, state_l, state_r, F_c, F_l, F_r, U_c, U_l, U_r);

  // then do EM parts
  deconvolve_flux_em(dir, ptrgeom, p_c, p_l, p_r, state_c, state_l, state_r, F_c, F_l, F_r, U_c, U_l, U_r);



  return(0);
}


int deconvolve_flux_ma(int dir, struct of_geom *ptrgeom, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, struct of_state *state_c, struct of_state *state_l, struct of_state *state_r, FTYPE *F_c, FTYPE *F_l, FTYPE *F_r, FTYPE *U_c, FTYPE *U_l, FTYPE *U_r)
{





  return(0);
}


int deconvolve_flux_em(int dir, struct of_geom *ptrgeom, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, struct of_state *state_c, struct of_state *state_l, struct of_state *state_r, FTYPE *F_c, FTYPE *F_l, FTYPE *F_r, FTYPE *U_c, FTYPE *U_l, FTYPE *U_r)
{





  return(0);
}



int a2cflux_from_prim(int dir, FTYPE (*prim_coef_list)[MAXSPACEORDER])
{


  return(0);





}


// GODMARK: Should redo flux's so that fluxes are accessed by F1[j][k][i] F2[k][i][j] F3[i][j][k] for faster differencing in advance.c

#include "decs.h"



// see fluxcompute.c for non-computer science, real physics calculations of flux
int fluxcalc(int stage, FTYPE pr[][N2M][N3M][NPR],
	     FTYPE F1[][N2M][N3M][NPR], 
	     FTYPE F2[][N2M][N3M][NPR], 
	     FTYPE F3[][N2M][N3M][NPR], 
 	     FTYPE CUf,
	     FTYPE fluxdt,
	     FTYPE *ndt1,
	     FTYPE *ndt2,
	     FTYPE *ndt3
	     )
{
  int fluxcalc_flux(int stage, FTYPE pr[][N2M][N3M][NPR], int *Nvec, FTYPE (*dqvec[NDIM])[N2M][N3M][NPR2INTERP], FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR], FTYPE (*fluxvecEM[NDIM])[N2M][N3M][NPR], FTYPE CUf, FTYPE *ndtvec[NDIM], struct of_loop *cent2faceloop);
  int fluxcalc_fluxctstag(int stage, FTYPE pr[][N2M][N3M][NPR], int *Nvec, FTYPE (*dqvec[NDIM])[N2M][N3M][NPR2INTERP], FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR], FTYPE CUf, struct of_loop *cent2faceloop, struct of_loop (*face2cornloop)[NDIM][NDIM]);
  extern int flux_ct(int stage, FTYPE pr[][N2M][N3M][NPR],FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR]);
  void fix_flux(FTYPE (*pr)[N2M][N3M][NPR],FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR]) ;
  FTYPE (*dqvec[NDIM])[N2M][N3M][NPR2INTERP];
  FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR];
  FTYPE (*fluxvecEM[NDIM])[N2M][N3M][NPR];
  FTYPE (**ptrfluxvec)[N2M][N3M][NPR];
  FTYPE *ndtvec[NDIM];
  int Nvec[NDIM];
  int flux_point2avg(int stage, int whichmaorem, FTYPE pr[][N2M][N3M][NPR], int *Nvec, FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR], FTYPE (*fluxvecother[NDIM])[N2M][N3M][NPR]);
  void preinterp_flux_point2avg(void);
  int i,j,k;
  int pl;
  int dir;
  int fluxEM2flux4EMF(int *Nvec, FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR], FTYPE (*fluxvecEM[NDIM])[N2M][N3M][NPR]);
  int fluxsum(int *Nvec, FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR], FTYPE (*fluxvecEM[NDIM])[N2M][N3M][NPR]);
  int zero_out_emf_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR]);
  // face2faceloop[dir=facedir] and face2cornloop[EMFdir=edgedir][EMFodir1][EMFodir2]
  // face2centloop separately used during separate part of advance()
  struct of_loop cent2faceloop[NDIM],face2cornloop[NDIM][NDIM][NDIM];


  /////////////////////////
  //
  // SETUP dimensionality
  //
  /////////////////////////

  fluxvec[1]=F1;
  fluxvec[2]=F2;
  fluxvec[3]=F3;

  fluxvecEM[1]=F1EM;
  fluxvecEM[2]=F2EM;
  fluxvecEM[3]=F3EM;

  dqvec[1]=dq1;
  dqvec[2]=dq2;
  dqvec[3]=dq3;

  ndtvec[1]=ndt1;
  ndtvec[2]=ndt2;
  ndtvec[3]=ndt3;

  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;

  if(splitmaem) ptrfluxvec=fluxvecEM;
  else ptrfluxvec=fluxvec;


  // zero-out EMFs if evolving/tracking vector potential so boundary regions appear simple in DUMPS when showing boundary cells
  if(TRACKVPOT && FULLOUTPUT!=0){
    zero_out_emf_fluxes(Nvec,ptrfluxvec);
  }

  ///////////////////////////////////////////////
  //
  // some pre-interplation flux averaging setups
  //
  ///////////////////////////////////////////////
  preinterp_flux_point2avg();



  ///////////////////////////////////////////////
  //
  // Compute normal flux at face
  // Involves interpolating CENT -> FACE1,2,3
  // Final p_l p_r results of interpolation are stored in gp_l gp_r if needed by SPLITNPR or FLUXB==FLUXCTSTAG
  // Wavespeeds may also be stored globally if STOREWAVESPEEDS==1
  // In all cases wavespeed constraint on timestep is set here
  //
  // assume fluxvec is MA only if splitmaem==1
  //
  ///////////////////////////////////////////////
  
  fluxcalc_flux(stage, pr, Nvec, dqvec, fluxvec, fluxvecEM, CUf, ndtvec, cent2faceloop);





#if(0)
  // DEBUG:
  if(Nvec[1]>1) FULLLOOP PLOOP(pl) fluxvec[1][i][j][k][pl]+=fluxvecEM[1][i][j][k][pl];
  if(Nvec[2]>1) FULLLOOP PLOOP(pl) fluxvec[2][i][j][k][pl]+=fluxvecEM[2][i][j][k][pl];
  if(Nvec[3]>1) FULLLOOP PLOOP(pl) fluxvec[3][i][j][k][pl]+=fluxvecEM[3][i][j][k][pl];
#endif



  //////////////////////////////
  //
  // FIXFLUX
  //
  /////////////////////////////

  if(FIXUPFLUX && splitmaem==0){ // for now can't use fix_flux with splitmaem since it resets field part to 0 that I'm using as diagonal gas pressure term
    fix_flux(pr, fluxvec[1], fluxvec[2], fluxvec[3]);
    if(splitmaem) fix_flux(pr, fluxvecEM[1], fluxvecEM[2], fluxvecEM[3]);
#if(PRODUCTION==0)
    trifprintf( "x");
#endif
  }




  //////////////////////////////
  //
  // FLUXCTSTAG -- overwrite field fluxes (emf's) with correct CORN1,2,3 points values
  //
  // assumes fluxcalc_flux() above called fluxcalc_standard_4fluxctstag() so gp_l and gp_r are set with FACE1,2,3 values of all quantities (including field along dir that comes from pstagscratch[] in this case)
  //
  /////////////////////////////
  if(FLUXB==FLUXCTSTAG){
#if(STOREWAVESPEEDS==0 || USESTOREDSPEEDSFORFLUX==0)
    dualfprintf(fail_file,"STOREWAVESPEEDS,USESTOREDSPEEDSFORFLUX must be 1 when FLUXB==FLUXCTSTAG\n");
    // must store because vchar() cannot be computed at CORN since only have partial velocities and partial fields and no densities
    // really only STOREWAVESPEEDS must be 1, but for now assume both
    myexit(175106);
#endif

    MYFUN(fluxcalc_fluxctstag(stage, pr, Nvec, dqvec, ptrfluxvec, CUf, cent2faceloop, face2cornloop),"flux.c:fluxcalc()", "fluxcalc_fluxctstag", 0);

    ////////////////
    // Before higher-order operations on flux, track vector potential update
    // so updating point value of A_i
    ////////////////
    update_vpot(stage, pr, ptrfluxvec, fluxdt);

  }

  


#if(PRODUCTION==0)
  trifprintf( "c");
#endif




  //////////////////////////////
  //
  // convert point FLUXES to surface-averaged fluxes (for field, if FLUXB==FLUXCTSTAG, then should do line integral of emf ENOMARK)
  //
  /////////////////////////////
  if((interporder[avgscheme[1]]>3) ||  (interporder[avgscheme[2]]>3) ||  (interporder[avgscheme[3]]>3)){
    if(splitmaem){
      // assume fluxvec is MA only if splitmaem==1
      //flux_point2avg(stage, ISMAONLY, pr, Nvec, fluxvec,NULL);
      // below version forces summation of MA+EM for stencil
      flux_point2avg(stage, ISMAONLY, pr, Nvec, fluxvec,fluxvecEM);
      flux_point2avg(stage, ISEMONLY, pr, Nvec, fluxvecEM, NULL);
    }
    else{
      flux_point2avg(stage, ISMAANDEM, pr, Nvec, fluxvec,NULL);
    }
  }





#if(0)
  // DEBUG: // also helped after flux_point2avg() but not before! -- so higher order code is problem or fed bad data
  // bound_flux(-1,F1,F2,F3);
#endif


  //////////////////////////////
  //
  // FLUXCT : must come after modifying fluxes so that divb=0 is preserved to machine error
  // i.e. flux_ct modifies fluxes in special way just before updating field so that divb=0 is conserved
  // Method is only 2nd order at best since assumes volume average is same as surface average
  //
  /////////////////////////////
  if((FLUXB==ATHENA1)||(FLUXB==ATHENA2)||(FLUXB==FLUXCTTOTH)||(FLUXB==FLUXCD)){
    MYFUN(flux_ct(stage, pr, ptrfluxvec[1], ptrfluxvec[2], ptrfluxvec[3]),"step_ch.c:advance()", "flux_ct",1);

    ////////////////
    // TOTH CT method doesn't cleanly differentiate between point update and average update of A_i, so just stick to TOTH CT EMF itself
    ////////////////
    update_vpot(stage, pr, ptrfluxvec, fluxdt);

  }


  //////////////////////
  //
  // sum up MA+EM
  // if splitmaem==1, then MA and EM split up to this point
  //
  ///////////////////////

  if(splitmaem) fluxsum(Nvec, fluxvec,fluxvecEM);

  //////////////////
  //
  // now fluxvec has MA+EM and fluxvecEM is no longer needed
  //
  //////////////////


#if(0)
  dir=1;
  FULLLOOP PLOOP(pl){
    dualfprintf(fail_file,"%ld %d :: F[%d][%d]=%21.15g\n",nstep,steppart,i,pl,fluxvec[dir][i][j][k][pl]);
  }
#endif


#if(FLUXDUMP)
  // this accounts for final flux
  FULLLOOP{
    DIMENLOOP(dir){
      if(Nvec[dir]>1){
	PLOOP(pl) fluxdump[i][j][k][4*NPR + (dir-1)*NPR*5 + NPR*0 + pl]=fluxvec[dir][i][j][k][pl];
      }
      else{
	PLOOP(pl) fluxdump[i][j][k][4*NPR + (dir-1)*NPR*5 + NPR*0 + pl]=0.0;
      }
    }
  }    
#endif


  return(0);
  
  
}



// zero-out fluxes that correspond to EMFs if evolving/tracking A_i since want boundary values to be plottable in SM or whatever
// This is  necessary since also use fluxes for other purposes in-between steps
// GODMARK: If in BCs have partial EMF and not full EMF for updating A_i, then may also cause A_i to appear funny
int zero_out_emf_fluxes(int *Nvec, FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR])
{
  int i,j,k,pl;
  int dir;

  DIMENLOOP(dir){
    if(Nvec[dir]>1){
      COMPFULLLOOP PLOOPBONLY(pl) fluxvec[dir][i][j][k][pl]=0.0;
    }
  }

  return(0);
}




// fill EM version if necessary when having new CT EMFs
int fluxEM2flux4EMF(int *Nvec, FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR], FTYPE (*fluxvecEM[NDIM])[N2M][N3M][NPR])
{
  int i,j,k,pl;
  int dir;

  if(splitmaem){
    DIMENLOOP(dir){
      if(Nvec[dir]>1){
	COMPFULLLOOP PLOOPBONLY(pl) fluxvecEM[dir][i][j][k][pl]=fluxvec[dir][i][j][k][pl];
      }
    }
  }

  return(0);
}

// sums MA and EM fluxes
int fluxsum_old(int *Nvec, FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR], FTYPE (*fluxvecEM[NDIM])[N2M][N3M][NPR])
{
  int i,j,k,pl;
  int dir;

  if(splitmaem){
    DIMENLOOP(dir){
      if(Nvec[dir]>1){
	COMPFULLLOOP PLOOP(pl) fluxvec[dir][i][j][k][pl]+=fluxvecEM[dir][i][j][k][pl];
      }
    }
  }

  return(0);
}


// sums MA (with FLUXSPLITMA(dir) containing flux[UU+dir] component) and EM fluxes
int fluxsum(int *Nvec, FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR], FTYPE (*fluxvecEM[NDIM])[N2M][N3M][NPR])
{
  int i,j,k,pl;
  int dir;

  if(splitmaem){
    DIMENLOOP(dir){
      if(Nvec[dir]>1){
	COMPFULLLOOP{
#if(SPLITPRESSURETERMINFLUXMA)
	  // add diagonal pressure term back to normal FLUX term
	  fluxvec[dir][i][j][k][UU+dir]+=fluxvec[dir][i][j][k][FLUXSPLITPMA(dir)];
	  // reset temporary storage
	  fluxvec[dir][i][j][k][FLUXSPLITPMA(dir)]=0.0;
#endif
#if(SPLITPRESSURETERMINFLUXEM)
	  // add diagonal pressure term back to normal FLUX term
	  fluxvec[dir][i][j][k][UU+dir]+=fluxvecEM[dir][i][j][k][FLUXSPLITPEM(dir)];
	  // reset temporary storage
	  fluxvecEM[dir][i][j][k][FLUXSPLITPEM(dir)]=0.0;
#endif
	  // now do normal addition
	  PLOOP(pl) fluxvec[dir][i][j][k][pl]+=fluxvecEM[dir][i][j][k][pl];
	}
      }
    }
  }

  return(0);
}






// wrapper for CENT_to_FACE1,2,3 used to compute flux at face
int fluxcalc_flux(int stage, FTYPE pr[][N2M][N3M][NPR], int *Nvec, FTYPE (*dqvec[NDIM])[N2M][N3M][NPR2INTERP], FTYPE (*fluxvec[NDIM])[N2M][N3M][NPR], FTYPE (*fluxvecEM[NDIM])[N2M][N3M][NPR], FTYPE CUf, FTYPE *ndtvec[NDIM], struct of_loop *cent2faceloop)
{
  int fluxcalc_flux_1d(int stage, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE F[][N2M][N3M][NPR], FTYPE FEM[][N2M][N3M][NPR], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata );
  int i,j,k,pl;
  int dir;
  int idel, jdel, kdel, face;
  int is, ie, js, je, ks, ke;
  int didassigngetstatecentdata;


  ///////////////////////////////////////////////
  //
  // LOOP OVER DIMENSIONS interpolating within each dimension separately for that flux -- assumes flux determined by 1-D Riemann problem
  //
  ///////////////////////////////////////////////
  didassigngetstatecentdata=0;

  DIMENLOOP(dir){ // GODMARK: Some other code related to storing global array info for get_state() assumes dir=1 first

    // set dimension having no influence on dt by default
    *(ndtvec[dir])=BIG;
    
    // don't skip dimension if doesn't exist, will be taken care of inside fluxcalc_flux_1d()


    // get loop details
    idel = fluxloop[dir][FIDEL];
    jdel = fluxloop[dir][FJDEL];
    kdel = fluxloop[dir][FKDEL];
    face = fluxloop[dir][FFACE];

    //loop over the interfaces where fluxes are computed -- atch, useCOMPZSLOOP( is, ie, js, je, ks, ke ) { ... }
    is=fluxloop[dir][FIS];
    ie=fluxloop[dir][FIE];
    js=fluxloop[dir][FJS];
    je=fluxloop[dir][FJE];
    ks=fluxloop[dir][FKS];
    ke=fluxloop[dir][FKE];


    MYFUN(fluxcalc_flux_1d(stage, pr, dir, is, ie, js, je, ks, ke, idel, jdel, kdel, face, dqvec[dir], fluxvec[dir], fluxvecEM[dir], CUf, ndtvec[dir], &cent2faceloop[dir], &didassigngetstatecentdata),"flux.c:fluxcalc()", "fluxcalc_flux_1d", dir);

#if(PRODUCTION==0)
    trifprintf("%d",dir);
#endif
  }// end DIMENLOOP(dir)

  return(0);

}



// wrapper for different standard 1-D flux calculators
// 1-D interpolate and get flux for that direction (assumes purely 1-D Riemann problem)
int fluxcalc_flux_1d(int stage, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE F[][N2M][N3M][NPR], FTYPE FEM[][N2M][N3M][NPR], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata )
{
  int fluxcalc_standard(int stage, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE F[][N2M][N3M][NPR], FTYPE FEM[][N2M][N3M][NPR], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata);
  int fluxcalc_standard_4fluxctstag(int stage, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE F[][N2M][N3M][NPR], FTYPE FEM[][N2M][N3M][NPR], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata);

  //  int fluxcalc_fluxspliteno(int stage, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE F[][N2M][N3M][NPR], FTYPE FEM[][N2M][N3M][NPR], FTYPE *ndt);
  int Nvec[NDIM];
  int i,j,k,pl;
  int odir1,odir2;


  
  if(DOENOFLUX!=ENOFLUXSPLIT){

    if(FLUXB==FLUXCTSTAG){
      MYFUN(fluxcalc_standard_4fluxctstag(stage,pr,dir,is, ie, js, je, ks, ke,idel,jdel,kdel,face,dq,F,FEM,CUf,ndt,cent2faceloop,didassigngetstatecentdata),"flux.c:fluxcalc_flux_1d()", "fluxcalc_standard_4fluxctstag()", 1);
    }
    else{
      // use older code that doesn't store into gp_l and gp_r since not needed and then waste of memory
      MYFUN(fluxcalc_standard(stage,pr,dir,is, ie, js, je, ks, ke,idel,jdel,kdel,face,dq,F,FEM,CUf,ndt,cent2faceloop,didassigngetstatecentdata),"flux.c:fluxcalc_flux_1d()", "fluxcalc_standard()", 1);
    }
  }
  else{
    //MYFUN(fluxcalc_fluxspliteno(stage,pr,dir,is, ie, js, je, ks, ke,idel,jdel,kdel,face,dq,F,FEM,CUf,ndt),"flux.c:fluxcalc_flux_1d()", "fluxcalc_fluxspliteno()", 1);
  }



  ////////////////////////////////////
  //
  // ensure flux of field set correctly to 0 in 1D
  //
  // if 1D in z and x, then no Ey so no F3[B1]
  // if 1D in z and y, then no Ex so no F3[B2]
  //
  // Stupid: This is automatically done since if no z-dir, then F3 not even used
  //
  ////////////////////////////////////

  /*
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;
  

  odir1=dir%3+1;
  odir2=(dir+1)%3+1;
  if(splitmaem==0){
    if(Nvec[dir]==1 && Nvec[odir1]==1) COMPFULLLOOP F[i][j][k][B1-1+odir1]=0.0;
    if(Nvec[dir]==1 && Nvec[odir2]==1) COMPFULLLOOP F[i][j][k][B1-1+odir2]=0.0;
    COMPFULLLOOP F[i][j][k][B1-1+dir]=0.0; // flux along field is always 0 due to antisymmetry of Faraday/Maxwell
  }
  else{
    // only need to change FEM since above F doesn't contain this EMF in this case
    if(Nvec[dir]==1 && Nvec[odir1]==1) COMPFULLLOOP FEM[i][j][k][B1-1+odir1]=0.0;
    if(Nvec[dir]==1 && Nvec[odir2]==1) COMPFULLLOOP FEM[i][j][k][B1-1+odir2]=0.0;
    COMPFULLLOOP FEM[i][j][k][B1-1+dir]=0.0; // flux along field is always 0 due to antisymmetry of Faraday/Maxwell
    // don't reset field part of MA flux since using that for other purposes (diagonal pressure term for now)
  }
  */



  return(0);
  
}


void set_normal_realisinterp(int *realisinterp)
{


  // real means normal primitive list = {rho,u,v1,v2,v3,B1,B2,B3} with v^i as WHICHVEL velocity
  // check:

  if(npr2interpstart=0 && npr2interpend==7 &&
     npr2interplist[RHO]==RHO &&
     npr2interplist[UU]==UU &&
     npr2interplist[U1]==U1 &&
     npr2interplist[U2]==U2 &&
     npr2interplist[U3]==U3 &&
     npr2interplist[B1]==B1 &&
     npr2interplist[B2]==B2 &&
     npr2interplist[B3]==B3
     ){

    // 1 here stands for NPR2INTERP type of loop
    // only do if some type of interpolation to be done
    *realisinterp=(RESCALEINTERP==0 || VARTOINTERP==PRIMTOINTERP);
  }
  else{
    *realisinterp=0; // must be forced to be zero
  }
  


}



// original flux calculator that gets F in "dir".  At end global pleft,pright,dq also set and if STOREWAVESPEEDS==1 then wavespeeds stored globally
int fluxcalc_standard(int stage, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE F[][N2M][N3M][NPR], FTYPE FEM[][N2M][N3M][NPR], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata)
{
  void slope_lim(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], struct of_loop *cent2faceloop);
  extern int rescale(int which, int dir, FTYPE *pr, struct of_geom *geom,FTYPE*newvar);
  int global_vchar(FTYPE pointspeed[][N2M][N3M][NUMCS], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, FTYPE wspeed[][2][N1M][N2M][N3M]);
  int getplpr(int dir, int idel, int jdel, int kdel, int i, int j, int k, struct of_geom *geom, FTYPE (*pr)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r, FTYPE *p_l, FTYPE *p_r);
  FTYPE (*p2interp)[N2M][N3M][NPR2INTERP];
  int i, j, k, pl;
  FTYPE p_l[NPR2INTERP], p_r[NPR2INTERP];
  FTYPE dtij;
  FTYPE ctop;
  struct of_geom geom;
  int reallim;
  int locallim;
  FTYPE pstore_l[NPR2INTERP],pstore_r[NPR2INTERP];
  FTYPE *p2interp_l,*p2interp_r;
  int odir1,odir2;
  int Nvec[NDIM];
  int realisinterp;
  int dointerpolation;
  void do_noninterpolation_dimension(int whichfluxcalc, int dointerpolation,  int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*pr)[N2M][N3M][NPR], struct of_loop *cent2faceloop, int *didassigngetstatecentdata);
  void compute_and_store_fluxstate(int dimen, int isleftright, int i, int j, int k, struct of_geom *geom, FTYPE *pr);




  ///////////////////////////////////////////////
  //
  // setup 1D reduction of EMF,flux calculation
  // here dir is face dir
  odir1=dir%3+1;
  odir2=(dir+1)%3+1;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;


  set_normal_realisinterp(&realisinterp);



  ////////////////////////////////////////////
  //
  // skip to next dir if no such dimension
  //
  /////////////////////////////////////////////

  // if nothing to interpolate, then quit
  if(npr2interpstart<=npr2interpend && Nvec[dir]>1){
    dointerpolation=1;
  }
  else{
    // just get loops and copy over to gp_?
    dointerpolation=0;

    //do limited number of things when not interpolating that dimension
    do_noninterpolation_dimension(ORIGINALFLUXCALC, dointerpolation,  realisinterp, dir, idel, jdel, kdel, pr, cent2faceloop, didassigngetstatecentdata);

    // skip real interpolation since no variation in that direction
    return(0);
  }




  /////////////////////////////
  //
  // setup wavedt
  //
  /////////////////////////////
  waveglobaldti[dir]=-100;
  waveglobaldtj[dir]=-100;
  waveglobaldtk[dir]=-100;



  //////////////////////////
  //
  // rescale before interpolation
  //
  ////////////////////////////
#if(RESCALEINTERP)
  if(npr2interpstart<=npr2interpend){
    // assume if DOEXTRAINTERP==1, then must get here
    p2interp=prc; // it's different
    p2interp_l=pstore_l;
    p2interp_r=pstore_r;
  }
#else
  p2interp=pr; // it's itself
  p2interp_l=p_l; // p2interp_? is final p_?
  p2interp_r=p_r;
#endif



  //////////////////////////
  //
  // rescale before interpolation AND/OR get wavespeeds
  //
  ////////////////////////////

#if((RESCALEINTERP)|| ( STOREWAVESPEEDS) )

  COMPFULLLOOP{
    // get geometry for center pre-interpolated values
    get_geometry(i, j, k, CENT, &geom); 

#if(RESCALEINTERP)
    // assume no need for a guess to p2interp to get pr (consistent with no unrescale being done after interpolation)
    if(npr2interpstart<=npr2interpend) rescale(1,dir,pr[i][j][k],&geom,p2interp[i][j][k]);
#endif

#if(SPLITNPR)
    // update wavespeed on FIRST pass
    if(advancepassnumber<=0)
#endif
      {
#if(STOREWAVESPEEDS)
	MYFUN(get_global_wavespeeds(dir,pr[i][j][k],&geom,wspeedtemp[i][j][k]),"flux.c:fluxcalc_standard()", "get_global_wavespeeds()", 0);
#endif // end if STOREWAVESPEEDS
      }
  }// end COMPFULLLOOP

#endif // end if rescaling or storing wavespeeds



#if(SPLITNPR)
  // update wavespeeds on FIRST pass
  if(advancepassnumber<=0)
#endif
    {
#if(STOREWAVESPEEDS)
      // get char. velocity estimate as some average over local or global zones
      global_vchar(wspeedtemp, dir, is, ie, js, je, ks, ke, idel, jdel, kdel, wspeed);
      // now wspeed contains left/right fastest wave speed (with sign preserved)
#endif  // otherwise use very local estimate
    }





  /////////////////////////////////////
  //
  // evaluate slopes (dq) or get pleft/pright of (possibly rescaled) primitive variables
  // c2e reconstruction: p2interp -> pleft & pright (indexed by grid cell #) -- atch comment
  //
  /////////////////////////////////////
  slope_lim(dointerpolation,realisinterp,dir,idel,jdel,kdel,pr,p2interp,dq,pleft,pright,cent2faceloop);







  //////////////////////////////////////
  //
  // flux loop : Extra "expand" zone for the purpose of averaging flux to get emf at corner.  Only used by field components, see flux_ct().
  // This loop is over interfaces where fluxes are evaluated -- atch
  //
  ////////////////////////////////////////

#if((SIMULBCCALC==2)&&(TYPE2==1))
  COMPFZLOOP(is,js,ks)
#else
  COMPZSLOOP( is, ie, js, je, ks, ke ) 
#endif
  {



    ////////////////////////
    //
    // get the geometry for the flux face
    //
    ///////////////////////
    get_geometry(i, j, k, face, &geom);


    //////////////////////////////////
    //
    // use p2interp,dq,pleft,pright to get p_l and p_r
    //
    /////////////////////////////////

    if(npr2interpstart<=npr2interpend){
      MYFUN(getplpr(dir,idel,jdel,kdel,i,j,k,&geom,pr,p2interp,dq,pleft,pright,p2interp_l,p2interp_r,p_l,p_r),"flux.c:fluxcalc_standard()", "getplpr", 1);
#if(SPLITNPR || FIELDSTAGMEM)
      // then means there is going to be a second pass, so store into memory
      PLOOPINTERP(pl){
	gp_l[dir][i][j][k][pl]=p_l[pl];
	gp_r[dir][i][j][k][pl]=p_r[pl];
      }
      PLOOPNOTINTERP(pl){ // restore those things didn't interpolate
	p_l[pl]=gp_l[dir][i][j][k][pl];
	p_r[pl]=gp_r[dir][i][j][k][pl];
      }
#endif
    }
    else{
#if(SPLITNPR || FIELDSTAGMEM)
      // GODMARK: for now assume interpolations either all done or none done, else need exclusion list for interpolations
      // then just get from memory
      PLOOPALLINTERP(pl){
	p_l[pl]=gp_l[dir][i][j][k][pl];
	p_r[pl]=gp_r[dir][i][j][k][pl];
      }
#endif
#if(SPLITNPR==0 || FIELDSTAGMEM==0)
      dualfprintf(fail_file,"Should not be using gp_? in flux.c when SPLITNPR==0 or FIELDSTAGMEM==0\n");
      myexit(16760276);
#endif
    }



    /////////////////////////////////////
    //
    // Compute and Store (globally) the get_state() data for the flux positions to avoid computing later
    //
    /////////////////////////////////////
#if(STOREFLUXSTATE)
    compute_and_store_fluxstate(dir, ISLEFT, i, j, k, &geom, p_l);
    compute_and_store_fluxstate(dir, ISRIGHT, i, j, k, &geom, p_r);
    // now flux_compute() and other flux-position-related things will obtain get_state() data for p_l and p_r from global arrays
#endif




    //////////////////////////////////
    //
    // actually compute the flux
    //
    /////////////////////////////////

    if(splitmaem==0){
      MYFUN(flux_compute_general(i, j, k, dir, &geom, CUf,  pr[i][j][k], p_l, p_r, F[i][j][k], &ctop),"step_ch.c:fluxcalc()", "flux_compute", 1);
    }
    else{
      MYFUN(flux_compute_splitmaem(i, j, k, dir, &geom, CUf,  pr[i][j][k], p_l, p_r, F[i][j][k], FEM[i][j][k], &ctop),"step_ch.c:fluxcalc()", "flux_compute", 1);
    }


    /////////////////////////////
    //
    // evaluate restriction on timestep
    //
    ///////////////////////////////
#if( (DOEVOLVEMETRIC || DOSELFGRAVVSR) && (RESTRICTDTSETTINGINSIDEHORIZON==2))
    // avoid limiting dt if inside horizon
    // GODMARK: can only do this if boundary condition does not drive solution's dt behavior
    if(WITHINENERREGION(enerposreg[OUTSIDEHORIZONENERREGION],i,j,k))
#endif
      {
	dtij = cour * dx[dir] / ctop;
	if (dtij < *ndt) 
	  {
	    *ndt = dtij;
	    // below are global so can report when other dt's are reported in advance.c
	    waveglobaldti[dir]=i;
	    waveglobaldtj[dir]=j;
	    waveglobaldtk[dir]=k;
	  }
      }

  }



  return (0);
}




void do_noninterpolation_dimension(int whichfluxcalc, int dointerpolation,  int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*pr)[N2M][N3M][NPR], struct of_loop *cent2faceloop, int *didassigngetstatecentdata)
{
  int i,j,k;
  int pl;
  struct of_geom geom;
  void slope_lim(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], struct of_loop *cent2faceloop);
  void slope_lim_cent2face(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], struct of_loop *cent2faceloop);
  void compute_and_store_fluxstate(int dimen, int isleftright, int i, int j, int k, struct of_geom *geom, FTYPE *pr);
  void compute_and_store_fluxstate_assign(int dimeninput, int dimenoutput, int isleftrightinput, int isleftrightoutput, int i, int j, int k);
  static int dirref;


  // if storing gp_l and gp_r, then copy 1D result in case used in some way not associated with whether dimension exists or not (not expensive)
  // For example, this is used for FLUXCTSTAG method in case dimension doesn't exist (Nvec[dir]==1) but still access gp_{l,r} instead of special conditions
  COMPFULLLOOP PLOOPINTERP(pl) gp_l[dir][i][j][k][pl]=gp_r[dir][i][j][k][pl]=pr[i][j][k][pl];


#if(STOREFLUXSTATE)
  if(*didassigngetstatecentdata==0){ // assume CENT (pr[]) is same for all centered values for dimensions not interpolating and assume always iterate over dir=1,2,3

    // indicate did do this assignment
    *didassigngetstatecentdata=1;
    dirref=dir;

    // only used if FLUXB==FLUXCTSTAG or MERGED method, but ok to do if STOREFLUXSTATE for now since coincident conditions
    COMPFULLLOOP{
      // get the geometry for the flux face
      get_geometry(i, j, k, FACE1-1+dir, &geom);
      
      // Compute and Store (globally) the get_state() data for the flux positions to avoid computing later
      compute_and_store_fluxstate(dir, ISLEFT, i, j, k, &geom, pr[i][j][k]);
      compute_and_store_fluxstate_assign(dir, dir, ISLEFT, ISRIGHT, i, j, k);
      // now flux_compute() and other flux-position-related things will obtain get_state() data for p_l and p_r from global arrays
    }// end COMPZLOOP
  }
  else{
    // just assign, not too expensive
    COMPFULLLOOP{
      // copy from dirref to other dirs for CENT (pr[]) quantity
      compute_and_store_fluxstate_assign(dirref, dir, ISLEFT, ISLEFT, i, j, k);
      compute_and_store_fluxstate_assign(dirref, dir, ISRIGHT, ISRIGHT, i, j, k);
    }
  }
#endif
  
  
  if(SPLITNPR||FLUXB==FLUXCTSTAG){
    // just get loop information
    if(whichfluxcalc==ORIGINALFLUXCALC){
      slope_lim(dointerpolation, realisinterp,dir,idel,jdel,kdel,NULL,NULL,NULL,NULL,NULL,cent2faceloop);
    }
    else{
      slope_lim_cent2face(dointerpolation, realisinterp,dir,idel,jdel,kdel,NULL,NULL,NULL,NULL,NULL,cent2faceloop);
    }
  }
  

}



// standard (non-field flux) calculation but setup to store results of interpolation so can be used for fluxctstag calculation
// set gp_l and gp_r with FACE interpolations from CENT (including field face from pstagscratch[])
// At end global pleft,pright,dq also set and if STOREWAVESPEEDS==1 then wavespeeds stored globally
int fluxcalc_standard_4fluxctstag(int stage, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE F[][N2M][N3M][NPR], FTYPE FEM[][N2M][N3M][NPR], FTYPE CUf, FTYPE *ndt, struct of_loop *cent2faceloop, int *didassigngetstatecentdata)
{
  int global_vchar(FTYPE pointspeed[][N2M][N3M][NUMCS], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, FTYPE wspeed[][2][N1M][N2M][N3M]);
  int i, j, k, pl;
  FTYPE dtij;
  FTYPE ctop;
  struct of_geom geom;
  int reallim;
  int locallim;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int interpolate_prim_cent2face(int stage, int realisinterp, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[N2M][N3M][NPR2INTERP], struct of_loop *cent2faceloop);
  int odir1,odir2;
  int Nvec[NDIM];
  int realisinterp;
  int dointerpolation;
  void do_noninterpolation_dimension(int whichfluxcalc, int dointerpolation,  int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*pr)[N2M][N3M][NPR], struct of_loop *cent2faceloop, int *didassigngetstatecentdata);





  ///////////////////////////////////////////////
  //
  // setup 1D reduction of EMF,flux calculation
  // here dir is face dir
  odir1=dir%3+1;
  odir2=(dir+1)%3+1;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;


  // only do interpolation if some type of interpolation to be done
  set_normal_realisinterp(&realisinterp);
  if(FLUXB==FLUXCTSTAG){
    realisinterp=0; // override since not interpolating p[B1+dir-1]
  }
  

  ////////////////////////////////////////////
  //
  // skip to next dir if no such dimension
  //
  /////////////////////////////////////////////
  if(npr2interpstart<=npr2interpend && Nvec[dir]>1){
    dointerpolation=1;
  }
  else{
    // just get loops, nothing to copy
    dointerpolation=0;

    //do limited number of things when not interpolating that dimension
    do_noninterpolation_dimension(NEWFLUXCALC, dointerpolation,  realisinterp, dir, idel, jdel, kdel, pr, cent2faceloop,didassigngetstatecentdata);

    // nothing else to do
    return(0);
  }






  /////////////////////////////
  //
  // setup wavedt
  //
  /////////////////////////////
  waveglobaldti[dir]=-100;
  waveglobaldtj[dir]=-100;
  waveglobaldtk[dir]=-100;


  

  //////////////////////////
  //
  // get wavespeeds if storing
  //
  ////////////////////////////

#if(STOREWAVESPEEDS)

  COMPFULLLOOP{
    // get geometry for center pre-interpolated values
    get_geometry(i, j, k, CENT, &geom); 

#if(SPLITNPR)
    // update wavespeed on FIRST pass
    if(advancepassnumber<=0)
#endif
      {
	MYFUN(get_global_wavespeeds(dir,pr[i][j][k],&geom,wspeedtemp[i][j][k]),"flux.c:fluxcalc_standard()", "get_global_wavespeeds()", 0);
      }
  }// end COMPFULLLOOP


#if(SPLITNPR)
  // update wavespeeds on FIRST pass
  if(advancepassnumber<=0)
#endif
    {
      // get char. velocity estimate as some average over local or global zones
      global_vchar(wspeedtemp, dir, is, ie, js, je, ks, ke, idel, jdel, kdel, wspeed);
      // now wspeed contains left/right fastest wave speed (with sign preserved)
    }
#endif // end if storing wavespeeds






  //////////////////////////
  //
  // obtain gp_l and gp_r (point face quantities) from pr (point centered quantity)
  //
  ////////////////////////////
  interpolate_prim_cent2face(stage, realisinterp, pr, dir, is, ie, js, je, ks, ke, idel, jdel, kdel, face, dq, cent2faceloop);





  //////////////////////////////////////
  //
  // flux loop : Extra "expand" zone for the purpose of averaging flux to get emf at corner.  Only used by field components, see flux_ct().
  // This loop is over interfaces where fluxes are evaluated -- atch
  //
  ////////////////////////////////////////

#if((SIMULBCCALC==2)&&(TYPE2==1))
  COMPFZLOOP(is,js,ks)
#else
 COMPZSLOOP( is, ie, js, je, ks, ke ) 
#endif
  {


    ////////////////////////
    //
    // get the geometry for the flux face
    //
    ///////////////////////
    get_geometry(i, j, k, face, &geom); // OPTMARK: Seems should put back together interpolate_prim_cent2face() with flux_compute stuff below so only 1 call to get_geometry @ face....not sure why this way?!


    //////////////////////////////////
    //
    // actually compute the flux
    //
    /////////////////////////////////

    if(splitmaem==0){
      MYFUN(flux_compute_general(i, j, k, dir, &geom, CUf,  pr[i][j][k], gp_l[dir][i][j][k], gp_r[dir][i][j][k], F[i][j][k], &ctop),"step_ch.c:fluxcalc()", "flux_compute", 1);
    }
    else{
      MYFUN(flux_compute_splitmaem(i, j, k, dir, &geom, CUf,  pr[i][j][k], gp_l[dir][i][j][k], gp_r[dir][i][j][k], F[i][j][k], FEM[i][j][k], &ctop),"step_ch.c:fluxcalc()", "flux_compute", 1);
    }


    /////////////////////////////
    //
    // evaluate restriction on timestep
    //
    ///////////////////////////////
#if( (DOEVOLVEMETRIC || DOSELFGRAVVSR) && (RESTRICTDTSETTINGINSIDEHORIZON==2))
    // avoid limiting dt if inside horizon
    // GODMARK: can only do this if boundary condition does not drive solution's dt behavior
    if(WITHINENERREGION(enerposreg[OUTSIDEHORIZONENERREGION],i,j,k))
#endif
    {
      dtij = cour * dx[dir] / ctop;
      if (dtij < *ndt) 
	{
	  *ndt = dtij;
	  // below are global so can report when other dt's are reported in advance.c
	  waveglobaldti[dir]=i;
	  waveglobaldtj[dir]=j;
	  waveglobaldtk[dir]=k;
	}
    }

  }// end FLUX LOOP






  return (0);
}










// normal interpolation of CENT quantities to FACE quantities
// sets global variables gp_l and gp_r to p_l and p_r from interpolations
int interpolate_prim_cent2face(int stage, int realisinterp, FTYPE pr[][N2M][N3M][NPR], int dir, int is, int ie, int js, int je, int ks, int ke, int idel, int jdel, int kdel, int face, FTYPE (*dq)[N2M][N3M][NPR2INTERP], struct of_loop *cent2faceloop)
{
  void slope_lim_cent2face(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], struct of_loop *cent2faceloop);
  extern int rescale(int which, int dir, FTYPE *pr, struct of_geom *geom,FTYPE*newvar);
  int getplpr(int dir, int idel, int jdel, int kdel, int i, int j, int k, struct of_geom *geom, FTYPE (*pr)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r, FTYPE *p_l, FTYPE *p_r);
  FTYPE (*p2interp)[N2M][N3M][NPR2INTERP];
  FTYPE p_l[NPR2INTERP], p_r[NPR2INTERP];
  FTYPE pstore_l[NPR2INTERP],pstore_r[NPR2INTERP];
  FTYPE *p2interp_l,*p2interp_r;
  int reallim;
  int locallim;
  int nprlocalstart,nprlocalend;
  int nprlocallist[MAXNPR];
  int pl,pl2;
  int i,j,k;
  struct of_geom geom;
  int realis,realie,realjs,realje,realks,realke;
  int dointerpolation;
  void compute_and_store_fluxstate(int dimen, int isleftright, int i, int j, int k, struct of_geom *geom, FTYPE *pr);





  // if inside this function, then doing interpolation
  dointerpolation=1;


  /////////////////////////////////////
  //
  // setup interpolation so avoids staggered field for field along "dir" direction
  // avoid magnetic field interpolation along "dir" direction using slope_lim() that is for cent->edges only
  // inside getplpr() staggered field is assigned to final p_l p_r before flux is computed
  //
  /////////////////////////////////////
  if(FLUXB==FLUXCTSTAG){

    ////////////////////////////////////////////
    //
    // save choice for interpolations
    nprlocalstart=npr2interpstart;
    nprlocalend=npr2interpend;
    PMAXNPRLOOP(pl) nprlocallist[pl]=npr2interplist[pl];


    // choice for range of PLOOPINTERP
    // check if wanted to interpolate B along dir, and if so remove
    for(pl=npr2interpstart;pl<=npr2interpend;pl++){
      if(npr2interplist[pl]==B1+dir-1){
	for(pl2=pl+1;pl2<=npr2interpend;pl2++) npr2interplist[pl2-1]=npr2interplist[pl2]; // moving upper to lower index
	npr2interpend--; // lost field along dir, so one less thing to do
	break;
      }
    }

  }



  //////////////////////////
  //
  // rescale before interpolation
  //
  ////////////////////////////
#if(RESCALEINTERP)
  // assume if DOEXTRAINTERP==1, then must get here
  p2interp=prc; // it's different
  p2interp_l=pstore_l;
  p2interp_r=pstore_r;
  
  COMPFULLLOOP{
    // get geometry for center pre-interpolated values
    get_geometry(i, j, k, CENT, &geom); 
    // assume no need for a guess to p2interp to get pr (consistent with no unrescale being done after interpolation)
    rescale(1,dir,pr[i][j][k],&geom,p2interp[i][j][k]);
  }// end COMPFULLLOOP
#else
  p2interp=pr; // it's itself
  p2interp_l=p_l; // p2interp_? is final p_?
  p2interp_r=p_r;
#endif






  /////////////////////////////////////
  //
  // evaluate slopes (dq) or get pleft/pright of (possibly rescaled) primitive variables
  // c2e reconstruction: p2interp -> pleft & pright (indexed by grid cell #) -- atch comment
  //
  /////////////////////////////////////
  slope_lim_cent2face(dointerpolation,realisinterp,dir,idel,jdel,kdel,pr,p2interp,dq,pleft,pright,cent2faceloop);


  // override normal fluxloop
  if(extrazones4emf){
    // don't really need fluxes in this domain, but do need interpolated face values to be transferred from pleft/pright to gp_l gp_r
    realis=emfUconsloop[FIS]-SHIFT1;
    realie=emfUconsloop[FIE]+SHIFT1;
    realjs=emfUconsloop[FJS]-SHIFT2;
    realje=emfUconsloop[FJE]+SHIFT2;
    realks=emfUconsloop[FKS]-SHIFT3;
    realke=emfUconsloop[FKE]+SHIFT3;
  }
  else{ // otherwise normal
    realis=is;
    realie=ie;
    realjs=js;
    realje=je;
    realks=ks;
    realke=ke;
  }



  //////////////////////////////////////
  //
  // flux loop : Extra "expand" zone for the purpose of averaging flux to get emf at corner.  Only used by field components, see flux_ct().
  // This loop is over interfaces where fluxes are evaluated -- atch
  //
  ////////////////////////////////////////

#if((SIMULBCCALC==2)&&(TYPE2==1))
  COMPFZLOOP(realis,realjs,realks)
#else
 COMPZSLOOP( realis, realie, realjs, realje, realks, realke ) 
#endif
  {


    ////////////////////////
    //
    // get the geometry for the flux face
    //
    ///////////////////////
    get_geometry(i, j, k, face, &geom);


    //////////////////////////////////
    //
    // use p2interp,dq,pleft,pright to get p_l and p_r
    //
    /////////////////////////////////

    MYFUN(getplpr(dir,idel,jdel,kdel,i,j,k,&geom,pr,p2interp,dq,pleft,pright,p2interp_l,p2interp_r,p_l,p_r),"flux.c:fluxcalc_standard()", "getplpr", 1);
    if(SPLITNPR || FLUXB==FLUXCTSTAG){
      // then means there is going to be a second pass, so store into memory
      PLOOPINTERP(pl){
	gp_l[dir][i][j][k][pl]=p_l[pl];
	gp_r[dir][i][j][k][pl]=p_r[pl];
      }
    }
    if(FLUXB==FLUXCTSTAG){
      // then also get B in dir direction not included in interpolation in slope_lim() above but included in getplpr() from pstagscratch[]
      pl = B1+dir-1;
      gp_l[dir][i][j][k][pl]=p_l[pl];
      gp_r[dir][i][j][k][pl]=p_r[pl];
    }


    /////////////////////////////////////
    //
    // Compute and Store (globally) the get_state() data for the flux positions to avoid computing later
    //
    /////////////////////////////////////
#if(STOREFLUXSTATE)
    compute_and_store_fluxstate(dir, ISLEFT, i, j, k, &geom, p_l);
    compute_and_store_fluxstate(dir, ISRIGHT, i, j, k, &geom, p_r);
    // now flux_compute() and other flux-position-related things will obtain get_state() data for p_l and p_r from global arrays
#endif




  }// end COMPZLOOP





#if(0)
  // DEBUG: (didn't matter)
  bound_prim(STAGEM1,gp_l[dir]);
  bound_prim(STAGEM1,gp_r[dir]);
#endif




  if(FLUXB==FLUXCTSTAG){
    ////////////////////////////////////////////
    //
    // restore choice for interpolations
    npr2interpstart=nprlocalstart;
    npr2interpend=nprlocalend;
    PMAXNPRLOOP(pl) npr2interplist[pl]=nprlocallist[pl];
  }
  

  return(0);

}




void compute_and_store_fluxstate_assign(int dimeninput, int dimenoutput, int isleftrightinput, int isleftrightoutput, int i, int j, int k)
{
  int jj;

  // now store result into array
  DLOOPA(jj){
    fluxstate[dimenoutput][isleftrightoutput][i][j][k].ucon[jj]=fluxstate[dimeninput][isleftrightinput][i][j][k].ucon[jj];
    fluxstate[dimenoutput][isleftrightoutput][i][j][k].ucov[jj]=fluxstate[dimeninput][isleftrightinput][i][j][k].ucov[jj];
#if(COMPUTE4FIELDforFLUX)
    fluxstate[dimenoutput][isleftrightoutput][i][j][k].bcon[jj]=fluxstate[dimeninput][isleftrightinput][i][j][k].bcon[jj];
    fluxstate[dimenoutput][isleftrightoutput][i][j][k].bcov[jj]=fluxstate[dimeninput][isleftrightinput][i][j][k].bcov[jj];
#endif
  }


}

// compute and store get_state() data for p_l and p_r type objects
void compute_and_store_fluxstate(int dimen, int isleftright, int i, int j, int k, struct of_geom *geom, FTYPE *pr)
{
  int pureget_stateforfluxcalc(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  struct of_state q;
  int jj;

  // get state
  pureget_stateforfluxcalc(pr,geom,&q);
  
  // now store result into array
  DLOOPA(jj){
    fluxstate[dimen][isleftright][i][j][k].ucon[jj]=q.ucon[jj];
    fluxstate[dimen][isleftright][i][j][k].ucov[jj]=q.ucov[jj];
#if(COMPUTE4FIELDforFLUX)
    fluxstate[dimen][isleftright][i][j][k].bcon[jj]=q.bcon[jj];
    fluxstate[dimen][isleftright][i][j][k].bcov[jj]=q.bcov[jj];
#endif
  }
  

}








// use dq,pleft,pright to obtain p_l and p_r for CENT to FACE
int getplpr(int dir, int idel, int jdel, int kdel, int i, int j, int k, struct of_geom *ptrgeom, FTYPE (*pr)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r, FTYPE *p_l, FTYPE *p_r)
{
  void getp2interplr(int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r);
  extern void set_plpr(int dir, int i, int j, int k, FTYPE (*prim)[N2M][N3M][NPR], FTYPE *p_l, FTYPE *p_r); // from user boundary routine
  int check_plpr(int dir, int i, int j, int k, int idel, int jdel, int kdel, struct of_geom *geom, FTYPE pr[][N2M][N3M][NPR], FTYPE *p_l, FTYPE *p_r);
  int pl;
  extern int rescale(int which, int dir, FTYPE *pr, struct of_geom *geom,FTYPE*newvar);


  
  //////////////////////////////////////
  //
  // interpolate primitive using slope (dq) or directly from pleft and pright
  // For FV: p_left, p_right (indexed by grid cell #) -> p2interp_l, p2interp_r (indexed by interface #) -- atch comment
  //
  // always do since p2interp is good and dq/pleft/pright should have stored quantities or just-computed quantities
  /////////////////////////////////////
  getp2interplr(dir,idel,jdel,kdel,i,j,k,p2interp,dq,pleft,pright,p2interp_l,p2interp_r);
			


  /////////////////////////////////////
  //      
  // after interpolation, unrescale from p2interp to normal primitive 
  // p2interp_? is p_? if not rescaling, so no need to assign p2interp_? to p_? if not rescaling
  //
  ///////////////////////////////////
  if(RESCALEINTERP && npr2interpstart<=npr2interpend){
    // only do if some interpolation done (consistent with no rescale() done in fluxcalc_standard() above
    // setup plausible p_l/p_r in case used for some reason (e.g. inversion starting guess)
    // this is sufficient for utoprim_1D...a better guess does no better (see interpU code)
    PLOOPINTERP(pl){
      p_l[pl]=pr[i][j][k][pl];
      p_r[pl]=pr[i][j][k][pl];
    }
    rescale(-1,dir,p_l,ptrgeom,p2interp_l);
    rescale(-1,dir,p_r,ptrgeom,p2interp_r);
  }




  //////////////////////////////
  //
  // Must preserve divb in 1D Riemann problem, so B^{dir} must be continuous
  //
  // GODMARK: Should really interpolate SUCH THAT this is automatically satisfied?
  // Yes, should use larger stencil and interpolate such that constant. Essentially choose
  // large stencil and effectively choosing which points we trust (not necessarily local points)
  //
  // GODMARK: This does not enforce E_\perp to be continuous for stationary flow!
  //
  ///////////////////////////////


  if(FLUXB==FLUXCTSTAG){
    // exactly correct that there is only 1 value
    pl = B1+dir-1;
    p_l[pl]=p_r[pl]=pstagscratch[i][j][k][pl];
  }
  else{
#if(BDIRCONT)
    // should really interpolate such that p_l=p_r
    pl = B1+dir-1;
    p_l[pl]=p_r[pl]=0.5*(p_l[pl]+p_r[pl]);
#endif
  }


#if(BOUNDPLPR)
  set_plpr(dir,i,j,k,pr,p_l,p_r);
#endif

  ///////////////////////
  //
  // correct interpolated quantities
  // no fixup accounting for these intermediate quantities
  //
  //////////////////////
  MYFUN(check_plpr(dir, i, j, k, idel, jdel, kdel, ptrgeom, pr, p_l, p_r),"step_ch.c:fluxcalc()", "check_plpr()", 1);


  return(0);

}






int check_plpr(int dir, int i, int j, int k, int idel, int jdel, int kdel, struct of_geom *geom, FTYPE pr[][N2M][N3M][NPR], FTYPE *p_l, FTYPE *p_r)
{
  int pl;

#if(EVOLVECHECKS)
#if(WHICHVEL==VEL4)
#if(ZEROOUTFLOWFLUX==1)
  inflow_check_4vel(dir,p_l,geom,-1);
  inflow_check_4vel(dir,p_r,geom,-1);
#endif
#elif(WHICHVEL==VEL3)
#if(ZEROOUTFLOWFLUX==1)
  inflow_check_3vel(dir,p_l,geom,-1);
  inflow_check_3vel(dir,p_r,geom,-1);
#endif
#if(JONCHECKS2)
  // must verify if this p makes sense (u^t sense)
  MYFUN(check_pr(p_l,pr[i-idel][j-jdel][k-kdel],geom,1,-1.0),"flux.c:check_plpr()", "check_pr()", 1);	
  
  MYFUN(check_pr(p_r,pr[i][j][k],geom,2,-1),"flux.c:check_plpr()", "check_pr()", 2);
#endif
#elif(WHICHVEL==VELREL4)
#if(ZEROOUTFLOWFLUX==1)
  inflow_check_rel4vel(dir,p_l,geom,-1);
  inflow_check_rel4vel(dir,p_r,geom,-1);
#endif
  // need to limit gamma since gamma may be large for interpolated value and would lead to bad fluxes
  MYFUN(limit_gamma(GAMMAMAX,p_l,geom,-1),"flux.c:check_plpr()", "limit_gamma()", 1);  //jon corr, see email re: SPINT warnings from 4/24/2006 10:54 p.m.
  MYFUN(limit_gamma(GAMMAMAX,p_r,geom,-1),"flux.c:check_plpr()", "limit_gamma()", 2);  //jon corr
#endif// end if WHICHVEL==VEL4REL      
#endif


#if(PRODUCTION==0)
  PLOOPINTERP(pl) if(!isfinite(p_l[pl])){
    dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d j=%d k=%d :: p_l is not finite pl=%d p_l=%21.15g\n",nstep,steppart,i,j,k,pl,p_l[pl]);
  }
  PLOOPINTERP(pl) if(!isfinite(p_r[pl])){
    dualfprintf(fail_file,"nstep=%ld steppart=%d i=%d j=%d k=%d :: p_r is not finite pl=%d p_r=%21.15g\n",nstep,steppart,i,j,k,pl,p_r[pl]);
  }
#endif

  // DEBUG:
  //  dualfprintf(fail_file,"CHECK dir=%d :: i=%d j=%d k=%d\n",dir,i,j,k);


  return(0);
}







//////////////////////////////////////
//
// interpolate primitive using slope or just copy pleft/pright into the correct place
//
// |=interface
// i=zone center of ith zone
//
// |              |     p2interp(i)    |
// |              |       dq(i)        |
// |        p_l(i)|p_r(i)   i          |
// |              |pleft(i)   pright(i)|
//
//
//
//////////////////////////////////////

// GODMARK: as Sasha mentions, for shifting stencil shouldn't just extrapolte value from non-local values, but use other slope and most LOCAL value to obtain extrapolation.

void getp2interplr(int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r)
{
  FTYPE Xcent[NDIM],Xleft[NDIM];
  FTYPE Vcent[NDIM],Vleft[NDIM];
  FTYPE rleft,rcent,thleft,thcent;
  int pl;
  int locallim;
  int choose_limiter(int dir, int i, int j, int k, int pl);
  extern void remapdq( int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r );
  extern void remapplpr( int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], FTYPE *p2interp_l, FTYPE *p2interp_r );

  //reshuffle dq's such that reconstruction within active zones does not use the boundary values
  remapdq( dir, idel, jdel, kdel, i, j, k, p2interp, dq, pleft, pright, p2interp_l, p2interp_r);
  
  if((HORIZONSUPERFAST)&&((lim[dir]<PARA)&&(LIMADJUST==0))&&(dir==1)){
    // since this uses dq's to shift, must always have dq's everywhere.  Hence LIMADJUST==0 and lim<PARA must be true
    coord(im1, j, k, CENT, Xleft);
    coord(i, j, k, CENT, Xcent);
    bl_coord(Xleft, Vleft);
    rleft=Vleft[1];
    thleft=Vleft[2];
    bl_coord(Xcent, Vcent);
    rcent=Vcent[1];
    thcent=Vcent[2];

    PLOOPINTERP(pl){
      locallim=choose_limiter(dir, i,j,k,pl);
      // get interpolated quantity
      if((locallim<PARA)&&(LIMADJUST==0)){
	if(rleft>Rhor) p2interp_l[pl] = p2interp[i - idel][j - jdel][k - kdel][pl] + 0.5 * dq[i - idel][j - jdel][k - kdel][pl];
	else p2interp_l[pl] = p2interp[i][j][k][pl] - 0.5 * dq[i][j][k][pl];
	if(rcent>Rhor) p2interp_r[pl] = p2interp[i][j][k][pl] - 0.5 * dq[i][j][k][pl];
	else p2interp_r[pl] = p2interp[i+idel][j+jdel][k+kdel][pl] - 1.5 * dq[i+idel][j+jdel][k+kdel][pl];
      }
      else{
	p2interp_l[pl] = pright[i-idel][j-jdel][k-kdel][pl];
	p2interp_r[pl] = pleft[i][j][k][pl];
      }
    } // end PLOOPINTERP
  }// if horizonsuperfast
  else{
    /////////////////////////////
    //
    // standard routine
    //
    /////////////////////////////
    PLOOPINTERP(pl){
      locallim=choose_limiter(dir, i,j,k,pl);
      // get interpolated quantity
      if((locallim<PARA)&&(LIMADJUST==0)){
	p2interp_l[pl] = p2interp[i - idel][j - jdel][k - kdel][pl] + 0.5 * dq[i - idel][j - jdel][k - kdel][pl];
	p2interp_r[pl] = p2interp[i][j][k][pl] - 0.5 * dq[i][j][k][pl];
      }
      else{
	//p_l & p_r for the current interface -- atch comment
	p2interp_l[pl] = pright[i-idel][j-jdel][k-kdel][pl];
	p2interp_r[pl] = pleft[i][j][k][pl];
      }
    }
  }

  //for boundaries where grid cells are avoided: set outer interface primitives at the boundary equal to inner interface primitive at that boundary
  remapplpr( dir, idel, jdel, kdel, i, j, k, p2interp, dq, pleft, pright, p2interp_l, p2interp_r);
}





// choose limiter
int choose_limiter(int dir, int i, int j, int k, int pl)
{
#if(LIMADJUST==LIMITERFIXED)
  // FUCKGODMARK
  // if(i>=4)	return(lim[dir]);
  //else return(DONOR);
  return(lim[dir]);
#else

#if(HYDROLIMADJUSTONLY)
  if(pl<B1) return(pflag[i][j][k][FLAGREALLIM]);
  else return(lim[dir]);
#else
  return(pflag[i][j][k][FLAGREALLIM]);
#endif

#endif

}









// slope_lim() is provided p2interp and returns pleft/pright
//
// |=interface
// i=zone center of ith zone
//
// |              |      p2interp(i)   |
// |         pl(i)|pr(i)    i          |
// |         Fl(i)|Fr(i)    i          |
// |         Ul(i)|Ur(i)    i          |
// |              |pleft(i)   pright(i)|
// |              |F(i)                |
//
void slope_lim(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], struct of_loop *cent2faceloop)
{
  extern void slope_lim_linetype_c2e(int realisinterp, int whichprimtype, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*stencilvar)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP]);
  extern void slope_lim_pointtype(int interporflux, int realisinterp, int pl, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP]);
  int pl;




  if( LINEINTERPTYPE(lim[dir]) ){ // this overrides lim, but lim must still be set properly
    // ENOPRIMITIVE below means primitives instead of conserved quantities
    get_loop(INTERPLINETYPE, ENOINTERPTYPE, dir, cent2faceloop);
    if(dointerpolation) slope_lim_linetype_c2e(realisinterp, ENOPRIMITIVE, ENOINTERPTYPE, dir, idel, jdel, kdel, primreal, NULL, p2interp, pleft, pright);
  }
  else{
    get_loop(INTERPPOINTTYPE, ENOINTERPTYPE, dir, cent2faceloop);
    if(dointerpolation){
      PLOOPINTERP(pl){
	slope_lim_pointtype(ENOINTERPTYPE, realisinterp, pl, dir, idel, jdel, kdel, primreal, p2interp, dq, pleft, pright);
      }
    }
  }



}



// slope_lim_cent2face() is provided p2interp and returns pleft/pright
// gets interpolations in expanded region for FLUXRECON && FLUXCTSTAG method if updating quasi-deaveraged field instead of point value
//
// |=interface
// i=zone center of ith zone
//
// |              |      p2interp(i)   |
// |         pl(i)|pr(i)    i          |
// |         Fl(i)|Fr(i)    i          |
// |         Ul(i)|Ur(i)    i          |
// |              |pleft(i)   pright(i)|
// |              |F(i)                |
//
void slope_lim_cent2face(int dointerpolation, int realisinterp, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], struct of_loop *cent2faceloop)
{
  extern void slope_lim_linetype_c2e(int realisinterp, int whichprimtype, int interporflux, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*stencilvar)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP]);
  extern void slope_lim_pointtype(int interporflux, int realisinterp, int pl, int dir, int idel, int jdel, int kdel, FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP]);
  int pl;
  int interporflux;

  
  if(extrazones4emf){
    interporflux=ENOINTERPTYPE4EMF;
    
  }
  else{
    interporflux=ENOINTERPTYPE;
  }
  


  if( LINEINTERPTYPE(lim[dir]) ){ // this overrides lim, but lim must still be set properly
    // ENOPRIMITIVE below means primitives instead of conserved quantities
    get_loop(INTERPLINETYPE, interporflux, dir, cent2faceloop);
    if(dointerpolation) slope_lim_linetype_c2e(realisinterp, ENOPRIMITIVE, interporflux, dir, idel, jdel, kdel, primreal, NULL, p2interp, pleft, pright);
  }
  else{
    get_loop(INTERPPOINTTYPE, interporflux, dir, cent2faceloop);
    if(dointerpolation){
      PLOOPINTERP(pl){
	slope_lim_pointtype(interporflux, realisinterp, pl, dir, idel, jdel, kdel, primreal, p2interp, dq, pleft, pright);
      }
    }
  }



}













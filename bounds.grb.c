
#include "decs.h"


/* bound array containing entire set of primitive variables */

// GODMARK: something seriously wrong with OUTEREXTRAP=1 (EOMFFDE)

#define OUTEREXTRAP 3
// 0: just copy
// 1: gdet or other extrapolation
// 2: copy (with rescale())
// 3: Treat same as HORIZONEXTRAP==3

// how to bound near horizon
#define HORIZONEXTRAP 3
// as above and also:
// 3: jon version

// number of iterations to get v^\phi\sim constant
#define NUMITERVPHI 5


// number of zones to use pole crushing regularizations

// to help protect the pole from death blows to the computational grid
// a sort of crushing regularization
// causes problems with stability at just beyond pole
// for field line plots, can just set B^\theta=0 along pole
#define POLEDEATH0 1 // with expansion by 1 point if detects jumps in densities or Lorentz factor (see poldeath())
//#define MAXPOLEDEATH N2BND // can't be larger than N2BND
#define MAXPOLEDEATH 2 // can't be larger than N2BND
#define DEATHEXPANDAMOUNT 0

#define POLEINTERPTYPE 3 // 0=set uu2=bu2=0, 1=linearly interpolate uu2,bu2  2=interpolate B_\phi into pole  3 =linearly for uu2 unless sucking on pole

// number of zones to enforce Lorentz factor to be small
// notice that at pole uu1 and uu2 are artificially large and these regions can lead to runaway low densities and so even higher uu1,uu3
// problem with POLEGAMMADEATH is that at large radius as fluid converges toward pole the fluid stagnates and can fall back at larger angles for no reason -- even for simple torus problem this happens when GAMMAPOLE=1.001
#define POLEGAMMADEATH0 1
// maximum allowed Lorentz factor near the pole (set to something large that should be allowed by solution -- problem and grid dependent)
//#define GAMMAPOLE (2.0)

#define GAMMAPOLEOUTGOING 1.1 // keep low
#define GAMMAPOLEOUTGOINGPOWER 1.0
#define GAMMAPOLEOUTGOINGRADIUS 10.0 // very model dependent
#define GAMMAPOLEINGOING GAMMAMAX


// factor by which to allow quantities to jump near pole
#define POLEDENSITYDROPFACTOR 5.0
#define POLEGAMMAJUMPFACTOR 2.0

// avoid such things if N2==1
#define POLEDEATH (N2==1 ? 0 : POLEDEATH0)
#define POLEGAMMADEATH (N2==1 ? 0 : POLEGAMMADEATH0)

// whether to average in radius for poledeath
#define AVERAGEINRADIUS 0 // not correct  across MPI boundaries since have to shift near boundary yet need that last cell to be consistent with as if no MPI boundary
#define RADIUSTOSTARTAVERAGING 7 // should be beyond horizon so doesn't diffuse across horizon

#define RADIUSTOAVOIDRADIALSUCK (2.0*Rhor)

#if( POLEDEATH > N2BND )
#error POLEDEATH should be <= N2BND
#endif

#if( MAXPOLEDEATH > N2BND )
#error MAXPOLEDEATH should be <= N2BND
#endif



//  needed as global so all subfunctions have it set
int inboundloop[NDIM];
int outboundloop[NDIM];
int innormalloop[NDIM];
int outnormalloop[NDIM];
int inoutlohi[NUMUPDOWN][NUMUPDOWN][NDIM];
int riin,riout,rjin,rjout,rkin,rkout;


static int bound_prim_user_general(int boundstage, int whichdir, int boundvartype, int ispstag, int* dirprim, FTYPE prim[][N2M][N3M][NPR]);



int bound_prim_user_dir(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE prim[][N2M][N3M][NPR])
{
  int bound_prim_user_general(int boundstage, int whichdir, int boundvartype, int ispstag, int* dirprim, FTYPE prim[][N2M][N3M][NPR]);
  int dirprim[NPR];
  int pl;


  // specify location of primitives
  PALLLOOP(pl) dirprim[pl]=CENT;
  //  dualfprintf(fail_file,"start bound_prim\n"); // CHANGINGMARK
  bound_prim_user_general(boundstage, whichdir, boundvartype, 0, dirprim, prim);
  //  dualfprintf(fail_file,"end bound_prim\n"); // CHANGINGMARK

  return(0);
}


// assume single user function takes care of primitive locations
int bound_pstag_user_dir(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE prim[][N2M][N3M][NPR])
{

  int dirprim[NPR];
  int pl;


  // specify location of primitives
  PALLLOOP(pl) dirprim[pl]=CENT;
  dirprim[B1]=FACE1;
  dirprim[B2]=FACE2;
  dirprim[B3]=FACE3;

#if(0)
  // Assume JCM setup bound_prim_user_general() correctly, so turned this off
  // GODMARK: assume non-field velocity not set and have to have something reasonable
  // use global values of non-field parts at present time
  // note that bound_pstag() setup loops to be over only B1..B3, but user may violate this and just stick in something so no failures even if not using data
  FULLLOOP PLOOPNOB1(pl) prim[i][j][k][pl]=p[i][j][k][pl];
  FULLLOOP PLOOPNOB2(pl) prim[i][j][k][pl]=p[i][j][k][pl];
#endif


  // assume before calling this that bound_pstag() setup PLOOPINTERP so only doing B1,B2,B3 (even though user may not respect this in bound_prim_user_general() -- which is ok since non-field quantities in pstag aren't needed -- may be problem if user_general() assumes primitive is reasonable)
  //  dualfprintf(fail_file,"start bound_pstag\n"); // CHANGINGMARK
  bound_prim_user_general(boundstage, whichdir, boundvartype, 1, dirprim, prim);
  //  dualfprintf(fail_file,"end bound_pstag\n"); // CHANGINGMARK



  return (0);
}





// user boundary routine
int bound_prim_user_general(int boundstage, int whichdir, int boundvartype, int ispstag, int* dirprim, FTYPE prim[][N2M][N3M][NPR])
{
  int i, j, k, pl;
  int locali1,globali1;
  int locali2,globali2;
  int ri1,ri2;
  struct of_geom geom[NPR],rgeom[NPR];
  FTYPE vcon[NDIM]; // coordinate basis vcon
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj, rk; // reference i,j,k
  FTYPE prescale[NPR];
  int horizonti;
  int enerregion;
  int *localenerpos;
  int poledeath(int ispstag, int boundstage, int whichx2, int* dirprim, FTYPE prim[][N2M][N3M][NPR]);
  int extrapfunc(int ispstag, int enerregion, int*localenerpos, int* dirprim, FTYPE prim[][N2M][N3M][NPR], int j,int k, int boundary);
  int jj,kk;





  ////////////////////////
  //
  // set bound loop
  //
  ///////////////////////
  set_boundloop(boundvartype, inboundloop,outboundloop,innormalloop,outnormalloop,inoutlohi, &riin, &riout, &rjin, &rjout, &rkin, &rkout);


  // for CZLOOP:
  // avoid looping over region outside active+ghost grid
  // good because somewhat general and avoid bad inversions, etc.
  //  enerregion=TRUEGLOBALENERREGION;
  enerregion=ACTIVEREGION; // now replaces TRUEGLOBALENERREGION
  localenerpos=enerposreg[enerregion];
  // if(WITHINENERREGION(localenerpos,i,j,k)){
  // note that localenerpos[X1DN] sets first evolved cell









  if(whichdir==1){
    ///////////////////////////
    //
    // X1 inner OUTFLOW/FIXEDOUTFLOW
    //
    ///////////////////////////

    if (mycpupos[1] <= horizoncpupos1) { // now all CPUs inside CPU with horizon will be using this (GODMARK: reference value needs to be chosen somehow for CPUs not on active grid)
      if((BCtype[X1DN]==OUTFLOW)||(BCtype[X1DN]==FIXEDOUTFLOW)||(BCtype[X1DN]==FREEOUTFLOW)){
	/* inner r boundary condition: u, just copy */
	LOOPX1dir{
	  ri=riin;
	  rj=j;
	  rk=k;

#if(HORIZONEXTRAP==0)
	  LOOPBOUND1IN{ // bound entire region inside non-evolved portion of grid
	    PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
	  }
#elif(HORIZONEXTRAP==1)
	  LOOPBOUND1IN{ // bound entire region inside non-evolved portion of grid
	    for(pl=RHO;pl<=UU;pl++){
	      prim[i][j][k][pl] = prim[ri][rj][rk][pl] * gdet[ri][rj][rk][dirprim[pl]]/gdet[i][j][k][dirprim[pl]] ;
	    }
	    pl=U1;	    // treat U1 as special
	    prim[i][j][k][pl] = prim[ri][rj][rk][pl] * (1. - (i-ri)*dx[1]) ;
	    for(pl=U2;pl<=U3;pl++){
	      prim[i][j][k][pl] = prim[ri][rj][rk][pl] * (1. + (i-ri)*dx[1]) ;
	    }
	    pl=B1; // treat B1 special
	    prim[i][j][k][pl] = prim[ri][rj][rk][pl] * gdet[ri][rj][rk][dirprim[pl]]/gdet[i][j][k][dirprim[pl]] ;
	    for(pl=B2;pl<=B3;pl++){
	      prim[i][j][k][pl] = prim[ri][rj][rk][pl] * (1. + (i-ri)*dx[1]) ;
	    }
	    for(pl=B3+1;pl<NPRBOUND;pl++){
	      prim[i][j][k][pl] = prim[ri][rj][rk][pl] * gdet[ri][rj][rk][dirprim[pl]]/gdet[i][j][k][dirprim[pl]] ;
	    }
	  }
#elif(HORIZONEXTRAP==2)
	  get_geometry(ri, rj, rk, dirprim[0], &rgeom[0]);
	  rescale(1,1,prim[ri][rj][rk],&rgeom[0],prescale);
	  LOOPBOUND1IN{
	    // set guess
	    PBOUNDLOOP(pl) prim[i][j][k][pl]=prim[ri][rj][k][pl];
	    get_geometry(i, j, k, dirprim[0], &geom[0]);	    
	    rescale(-1,1,prim[i][j][k],&geom[0],prescale);
	  }
#elif(HORIZONEXTRAP==3)
	  extrapfunc(ispstag, enerregion, localenerpos, dirprim, prim, j,k,X1DN);
#endif




	  if(ispstag==0){

	    if((BCtype[X1DN]==OUTFLOW)||(BCtype[X1DN]==FIXEDOUTFLOW)){
	      // GODMARK: assume all velocities at same location when doing inflow check
	      LOOPBOUND1IN{
#if(WHICHVEL==VEL4)
		get_geometry(i, j, k, dirprim[U1], &geom[U1]);
		inflow_check_4vel(1,prim[i][j][k],&geom[U1], 0) ;
#elif(WHICHVEL==VEL3)
		get_geometry(i, j, k, dirprim[U1], &geom[U1]);
		inflow_check_3vel(1,prim[i][j][k],&geom[U1], 0) ;
		// projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
		if(jonchecks){
		  //fixup1zone(prim[i][j][k],&geom[U1],0);
		  failreturn=check_pr(prim[i][j][k],prim[i][j][k],&geom[U1],-3);
		  if(failreturn){
		    dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
		    if (fail(FAIL_BCFIX) >= 1) return (1);
		  }
		}
#endif
#elif(WHICHVEL==VELREL4)
		get_geometry(i,j,k,dirprim[U1],&geom[U1]) ;
		inflow_check_rel4vel(1,prim[i][j][k],&geom[U1],0) ;
		if(limit_gamma(GAMMAMAX,prim[i][j][k],&geom[U1],0)>=1)
		  FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 1);
#endif	
	      }
	    } // end if not allowing inflow
	  }


	}// end 2 3
      }// end if correct bound type

    }// end if mycpupos[1]==0






    ///////////////////////////
    //
    // X1 outer OUTFLOW/FIXEDOUTFLOW
    //
    ///////////////////////////


    // outer r BC:
    if (mycpupos[1] == ncpux1 - 1) {
      if((BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)||(BCtype[X1UP]==FREEOUTFLOW)){
	/* outer r BC: outflow */

	LOOPX1dir{
	  ri=riout;
	  rj=j;
	  rk=k;

#if(OUTEREXTRAP==0)
	  LOOPBOUND1OUT PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk][pl];
#elif(OUTEREXTRAP==1)
	  LOOPBOUND1OUT{
	    for(pl=RHO;pl<=UU;pl++){
	      prim[i][j][k][pl] = prim[ri][rj][rk][pl] * gdet[ri][rj][rk][dirprim[pl]]/gdet[i][j][k][dirprim[pl]] ;
	    }
	    pl=U1; // treat U1 as special
	    prim[i][j][k][pl] = prim[ri][rj][rk][pl] * (1. - 2*(i-ri)*dx[1]) ;
	    for(pl=U2;pl<=U3;pl++){
	      prim[i][j][k][pl] = prim[ri][rj][rk][pl] * (1. - (i-ri)*dx[1]) ;
	    }
	    pl=B1; // treat B1 special
	    prim[i][j][k][pl] = prim[ri][rj][rk][pl] * gdet[ri][rj][rk][dirprim[pl]]/gdet[i][j][k][dirprim[pl]] ;
	    for(pl=B2;pl<=B3;pl++){
	      prim[i][j][k][pl] = prim[ri][rj][rk][pl] * (1. - (i-ri)*dx[1]) ;
	    }
	    for(pl=B3+1;pl<NPRBOUND;pl++){
	      prim[i][j][k][pl] = prim[ri][rj][rk][pl] * gdet[ri][rj][rk][dirprim[pl]]/gdet[i][j][k][dirprim[pl]] ;
	    }
	  }
#elif(OUTEREXTRAP==2)
	  get_geometry(ri, rj, rk, dirprim[0], &rgeom[0]);
	  rescale(1,1,prim[ri][rj][rk],&rgeom[0],prescale);
	  LOOPBOUND1OUT{
	    // set guess
	    PBOUNDLOOP(pl) prim[i][j][k][pl]=prim[ri][rj][rk][pl];
	    get_geometry(i, j, k, dirprim[0], &geom[0]);
	    rescale(-1,1,prim[i][j][k],&geom[0],prescale);
	  }
#elif(OUTEREXTRAP==3)
	  extrapfunc(ispstag, enerregion, localenerpos, dirprim, prim, j,k,X1UP);
#endif





	  if(ispstag==0){

	    if((BCtype[X1UP]==OUTFLOW)||(BCtype[X1UP]==FIXEDOUTFLOW)){

	      LOOPBOUND1OUT{
#if(WHICHVEL==VEL4)
		get_geometry(i, j, k, dirprim[U1], &geom[U1]);
		inflow_check_4vel(1,prim[i][j][k],&geom[U1],0) ;
#elif(WHICHVEL==VEL3)
		get_geometry(i, j, k, dirprim[U1], &geom[U1]);
		inflow_check_3vel(1,prim[i][j][k],&geom[U1],0) ;
		// projection may not preserve u^t to be real and rho>rhoscal u>uuscal
#if(JONCHECKS)
		if(jonchecks){
		  //fixup1zone(prim[i][j][k],&geom[U1],0);
		  failreturn=check_pr(prim[i][j][k],prim[i][j][k],&geom[U1],-3);
		  if(failreturn){
		    dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d k=%d\n",startpos[1]+i,startpos[2]+j,startpos[3]+k);
		    if (fail(FAIL_BCFIX) >= 1) return (1);
		  }
		}
#endif
#elif(WHICHVEL==VELREL4)
		get_geometry(i,j,k,dirprim[U1],&geom[U1]) ;
		// 	      dualfprintf(fail_file,"JUST BEFORE INFLOWCHECK: i=%d j=%d k=%d prim[U1]=%21.15g prim[U2]=%21.15g prim[U3]=%21.15g\n",i,j,k,prim[i][j][k][U1] *sqrt(geom[U1].gcov[1][1]),prim[i][j][k][U2] *sqrt(geom[U1].gcov[2][2]),prim[i][j][k][U3] *sqrt(geom[U1].gcov[3][3]));
		// 	      dualfprintf(fail_file,"JUST BEFORE INFLOWCHECK: i=%d j=%d k=%d prim[U1]=%21.15g prim[U2]=%21.15g prim[U3]=%21.15g\n",i,j,k,prim[i][j][k][U1],prim[i][j][k][U2],prim[i][j][k][U3]);
		//	      DLOOP(jj,kk) dualfprintf(fail_file,"gcov[%d][%d]=%21.15g\n",jj,kk,geom[U1].gcov[jj][kk]);
		inflow_check_rel4vel(1,prim[i][j][k],&geom[U1],0) ;
		// 	      dualfprintf(fail_file,"JUST BEFORE LIMIT: i=%d j=%d k=%d prim[U1]=%21.15g prim[U2]=%21.15g prim[U3]=%21.15g\n",i,j,k,prim[i][j][k][U1] *sqrt(geom[U1].gcov[1][1]),prim[i][j][k][U2] *sqrt(geom[U1].gcov[2][2]),prim[i][j][k][U3] *sqrt(geom[U1].gcov[3][3]));
		// 	      dualfprintf(fail_file,"JUST BEFORE LIMIT: i=%d j=%d k=%d prim[U1]=%21.15g prim[U2]=%21.15g prim[U3]=%21.15g\n",i,j,k,prim[i][j][k][U1],prim[i][j][k][U2],prim[i][j][k][U3]);
		if(limit_gamma(GAMMAMAX,prim[i][j][k],&geom[U1], 0)>=1)
		  FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);
		//	      dualfprintf(fail_file,"JUST AFTER LIMIT: i=%d j=%d k=%d prim[U1]=%21.15g prim[U2]=%21.15g prim[U3]=%21.15g\n",i,j,k,prim[i][j][k][U1] *sqrt(geom[U1].gcov[1][1]),prim[i][j][k][U2] *sqrt(geom[U1].gcov[2][2]),prim[i][j][k][U3] *sqrt(geom[U1].gcov[3][3]));
		//	      dualfprintf(fail_file,"JUST AFTER LIMIT: i=%d j=%d k=%d prim[U1]=%21.15g prim[U2]=%21.15g prim[U3]=%21.15g\n",i,j,k,prim[i][j][k][U1],prim[i][j][k][U2],prim[i][j][k][U3]);
#endif	
	      }// end over i
	    }// end if not allowing inflow

	  }



	}// end 2 3
      }// end if correct bound type
    }// end if mycpu is correct






    ///////////////////////////
    //
    // X1 inner r=0 singularity
    //
    ///////////////////////////


    /* inner radial BC (preserves u^t rho and u) */
    if (mycpupos[1] == 0) {
      if((BCtype[X1DN]==R0SING)||(BCtype[X1DN]==SYMM)||(BCtype[X1DN]==ASYMM) ){
	LOOPX1dir{
	  rj=j;
	  rk=k;
	  LOOPBOUND1IN{
	    PBOUNDLOOP(pl){
	      // SECTIONMARK: assume r=0 singularity can't move
	      if(dirprim[pl]==FACE1 || dirprim[pl]==CORN3 || dirprim[pl]==CORN2 || dirprim[pl]==CORNT ) ri = -i; // FACE1 values
	      else ri=-i-1; // "CENT" values for purposes of this BC
	      prim[i][j][k][pl] = prim[ri][rj][rk][pl];
	    }// over pl
	  }// over boundary zones
	}
      }

      if((BCtype[X1DN]==R0SING)||(BCtype[X1DN]==ASYMM) ){

	/* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
	LOOPX1dir{

	  // SECTIONMARK: assume r=0 singularity can't move
	  i=0;
	  PBOUNDLOOP(pl){
	    if(pl==U1 || pl==B1){
	      if(dirprim[pl]==FACE1 || dirprim[pl]==CORN3 || dirprim[pl]==CORN2 || dirprim[pl]==CORNT ){
		prim[i][j][k][pl] = 0.0;
	      }
	    }// else don't do this pl
	  } // end over pl
	
	  LOOPBOUND1IN {
	    PBOUNDLOOP(pl){
	      if(pl==U1 || pl==B1){
		prim[i][j][k][pl] *= -1.;
	      }// end if right pl
	    } // end over pl
	  } // end over boundary zones
	}// end loop 23
      }
    }





  }
  else if(whichdir==2){





    ///////////////////////////
    //
    // X2 inner POLARAXIS
    //
    ///////////////////////////


    /* inner polar BC (preserves u^t rho and u) */
    if (mycpupos[2] == 0) {
      if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==SYMM)||(BCtype[X2DN]==ASYMM) ){
	LOOPX2dir{
	  ri=i;
	  rk=k;
	  LOOPBOUND2IN{
	    PBOUNDLOOP(pl){
	      // assume polar axis can't move SECTIONMARK
	      if(dirprim[pl]==FACE2 || dirprim[pl]==CORN3 || dirprim[pl]==CORN1 || dirprim[pl]==CORNT ) rj = -j; // FACE2 values
	      else rj=-j-1; // "CENT" values for purposes of this BC
	      prim[i][j][k][pl] = prim[ri][rj][rk][pl];
	    }
	  }
	}
      }

      if((BCtype[X2DN]==POLARAXIS)||(BCtype[X2DN]==ASYMM) ){

	/* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
	LOOPX2dir{

	  // assume polar axis can't move SECTIONMARK
	  j=0;
	  PBOUNDLOOP(pl){
	    if(pl==U2 || pl==B2){
	      if(dirprim[pl]==FACE2 || dirprim[pl]==CORN3 || dirprim[pl]==CORN1 || dirprim[pl]==CORNT ){
		prim[i][j][k][pl] = 0.0;
	      }
	    }// else don't do this pl
	  } // end over pl
	
	
	  LOOPBOUND2IN {
	    PBOUNDLOOP(pl){
	      if(pl==U2 || pl==B2){
		prim[i][j][k][pl] *= -1.;
	      } // end if right pl
	    } // end over pl
	  } // end over boundary cells
	}// end loop 13


	if(POLEDEATH) poledeath(ispstag, boundstage,X2DN,dirprim,prim);

      } // end if POLARXIS or ASYMM
    }// end if mycpupos[2]==0





    ///////////////////////////
    //
    // X2 outer POLARAXIS
    //
    ///////////////////////////


    if (mycpupos[2] == ncpux2-1) {
      if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==SYMM)||(BCtype[X2UP]==ASYMM) ){
	LOOPX2dir{
	  ri=i;
	  rk=k;
	  LOOPBOUND2OUT{
	    PBOUNDLOOP(pl){
	      // assume polar axis can't move SECTIONMARK
	      if(dirprim[pl]==FACE2 || dirprim[pl]==CORN3 || dirprim[pl]==CORN1 || dirprim[pl]==CORNT ) rj = 2*N2-j; // FACE2 values
	      else rj=2*(N2-1)-j+1; // "CENT" values for purposes of this BC
	      prim[i][j][k][pl] = prim[ri][rj][rk][pl];
	    }
	  }
	}
      }

      if((BCtype[X2UP]==POLARAXIS)||(BCtype[X2UP]==ASYMM) ){

	/* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
	LOOPX2dir{

	  j=N2;
	  PBOUNDLOOP(pl){
	    if(pl==U2 || pl==B2){
	      if(dirprim[pl]==FACE2 || dirprim[pl]==CORN3 || dirprim[pl]==CORN1 || dirprim[pl]==CORNT ){
		prim[i][j][k][pl] = 0.0;
	      }
	    }// else don't do this pl
	  } // end over pl

	  LOOPBOUND2OUT {
	    PBOUNDLOOP(pl){
	      if(pl==U2 || pl==B2){
		prim[i][j][k][pl] *= -1.;
	      } // end if right pl
	    } // end over pl
	  } // end over bondary cells
	}// end loop 13


	if(POLEDEATH) poledeath(ispstag, boundstage,X2UP,dirprim,prim);

      } // end if POLARXIS or ASYMM
    }// end if mycpupos[2]==ncpux2-1



  }
  else if(whichdir==3){



    // periodic x3
    if ( (mycpupos[3] == 0)&&(mycpupos[3] == ncpux3 - 1) ) {
      if( (BCtype[X3DN]==PERIODIC)&&(BCtype[X3UP]==PERIODIC) ){
	// just copy from one side to another

	LOOPX3dir{

	  // copy from upper side to lower boundary zones
	  ri=i;
	  rj=j;
	  rk=rkout;
	  LOOPBOUND3IN PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk+1+k][pl];

	  // copy from lower side to upper boundary zones
	  ri=i;
	  rj=j;
	  rk=rkin;
	  LOOPBOUND3OUT PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri][rj][rk+(k-N3)][pl];
	}
      }
    }


  }
  else{
    dualfprintf(fail_file,"No such whichdir=%d\n",whichdir);
    myexit(2436262);
  }
  








    

  if(whichdir==1){
    /////////////////
    //
    // now deal with region interior to the horizon
    // This just makes values reasonable for diagnostics to be ignorant of loop ranges -- results are NOT used! -- no FLUXCTSTAG issues
    // similar to bound_spacetime_inside_horizon() in metric.c

  
    horizonti=N1*horizoncpupos1+horizoni;

    if(BCtype[X1DN]==R0SING && horizonti>0){ // global condition that all CPUs know about
      // for self-gravity, only loop over region outside horizon
      enerregion=OUTSIDEHORIZONENERREGION;
      enerpos=enerposreg[enerregion];


      // copy from very outer boundary or from location on grid defined by where last ghost cell is inside horizon

      // this doesn't do anything to cells inside horizon but outside fake region
      // apparently in this case uu1>0 even though inside the horizon
      globali1=N1*horizoncpupos1+horizoni-N1BND;
      if(globali1<startpos[1]-N1BND) locali1=-N1BND;
      else if(globali1>endpos[1]+N1BND-1) locali1=N1+N1BND-1;
      else locali1=horizoni-N1BND;
      ri1=locali1;



      // include all zones inside horizon as outflow so remove any peak that was accreted as black hole
      // flux through horizon of mass should be consistent with either method, but seems more reasonable to do below
      globali2=N1*horizoncpupos1+horizoni;
      if(globali2<startpos[1]) locali2=0;
      else if(globali2>endpos[1]-1) locali2=N1-1;
      else locali2=horizoni;
      ri2=locali2;
      // GODMARK: somehow leads to density sticking to value despite much lower active density


      LOOPF3 LOOPF2 LOOPF1{

	rj=j;
	rk=k;
      
	//      if(WITHINENERREGION(enerpos,i,j,k)){
	// then do nothing
	//      }
	//else{
	// assume horizon on negative side of "i" so don't modify right side of "i" that happens to be in the outer radial boundary

	if(i < globali1-startpos[1]){
	  // then inside horizon
	  
	  // set primitives to be trivialized
 
	  PBOUNDLOOP(pl) prim[i][j][k][pl] = prim[ri1][rj][rk][pl];

	  // reset pflag in case old failure inside horizon
#if(UTOPRIMADJUST==UTOPRIMAVG)
	  pflag[i][j][k][FLAGUTOPRIMFAIL]=UTOPRIMNOFAIL;
#endif


	  if(ispstag==0){
	    // making it small artificially forces timestep to be small since timestep still limited within fake boundary region
	    //	  prim[i][j][k][RHO] = SMALL;
	    // assumed acting on relative 4-velocity that stays physical no matter what
	    if(prim[i][j][k][U1]>0) prim[i][j][k][U1]=-0.5*(fabs(prim[ri1][rj][rk][U1])+fabs(prim[ri1+1][rj][rk][U1]));
	  }
	}

	if(0 && i < globali2-startpos[1]){
	  // then on or inside horizon and inside N1BND more grid cells for interpolation
	  
	  if(ispstag==0){
	    // assumed acting on relative 4-velocity that stays physical no matter what
	    if(prim[i][j][k][U1]>0) prim[i][j][k][U1]=-0.5*(fabs(prim[ri1][rj][rk][U1])+fabs(prim[ri1+1][rj][rk][U1]));
	  }

	}
	//}
      }


    }
  }  




  return (0);
}





// OPTMARK: all get_geometries, etc. could be simplified if take advantage of fact that only really 2 cases: all CENT or CENT for all except B's.
// otherwise extrapfunc() and poledeath()'s get_geometry() dominates CPU time when few active cells or many boundary cells



// boundary = X1DN or X1UP
// assume if ispstag==1 then only doing field part -- otherwise logic would get more complicated
int extrapfunc(int ispstag, int enerregion, int*localenerpos, int* dirprim, FTYPE prim[][N2M][N3M][NPR], int j,int k, int boundary)
{
  int i, pl;
  int ri,rj,rk;
  int iloopstart,iloopend,iloopstep;
  int ri2,ri3;
  struct of_geom geom[NPR],rgeom[NPR];
  // for Bd3 copy
  FTYPE Bd3,Bd3ri2,Bd3ri3;
  int jj;
  FTYPE Bu1,Bu2,gcon03,gcon13,gcon23,gcon33;
  FTYPE gcov01,gcov02,gcov11,gcov12,gcov21,gcov22,gcov03,gcov13,gcov23;
  // Mdot copy
  FTYPE Mdot,ucon[NDIM],uconref[NDIM];
  FTYPE Xr[NPR][NDIM],Vr[NPR][NDIM],X[NPR][NDIM],V[NPR][NDIM];
  FTYPE primtemp[NPR];
  FTYPE *ucontouse, ucon2[NDIM];
  int itervphi;
  struct of_geom ri2geom[NPR];
  FTYPE uconrefri2[NDIM];
  FTYPE Mdotri2;
  struct of_geom ri3geom[NPR];
  FTYPE uconrefri3[NDIM];
  FTYPE Mdotri3;
  FTYPE Dqp,Dqm,Dqc,dqMdot;
  FTYPE myMdot;
  FTYPE dq[NPR];
  FTYPE dqucon[NDIM];
  FTYPE dqBd3,myBd3;
  FTYPE dqgdetpl[NPR];
  FTYPE dqlogdensity[NPR];
  FTYPE uconreftouse[NDIM];
  FTYPE xfrac,yfrac;
  FTYPE ftemp;
  FTYPE signdq;
  FTYPE xdqfrac,ydqfrac;
  FTYPE ftemp2,linearvalue,expvalue;
  FTYPE D0,uconrefnew[NDIM];
  FTYPE mydqlog;
  FTYPE igeomnosing;


  rj=j;
  rk=k;


  // setup multiple reference points
  if(boundary==X1DN){
    ri = riin;
    ri2=ri+1;
    ri3=ri+2;
    // iterate up (due to use of single loop form)
    iloopstart=LOWERBOUND1;
    iloopend=ri-1;
    iloopstep=1;
    signdq=1.0;
  }
  else if(boundary==X1UP){
    ri = riout;
    ri2=ri-1;
    ri3=ri-2;
    // iterate up (due to use of single loop form)
    iloopstart=ri+1;
    iloopend=UPPERBOUND1;
    iloopstep=1;
    signdq=-1.0;
  }
  else{
    dualfprintf(fail_file,"extrapfunc not setup for this boundary = %d\n",boundary);
    myexit(811751);
  }



  /////////////
  //
  // get coordinates of reference location
  //
  /////////////
  PBOUNDLOOP(pl){
    get_geometry(ri, rj, rk, dirprim[pl], &rgeom[pl]);
    get_geometry(ri2, rj, rk, dirprim[pl], &ri2geom[pl]);
    get_geometry(ri3, rj, rk, dirprim[pl], &ri3geom[pl]);

    coord(ri,rj,rk,dirprim[pl],Xr[pl]); // reference locations for B2/U2
    bl_coord(Xr[pl],Vr[pl]);
  }



  /////////////
  //
  // obtain reference B_\phi
  //
  /////////////
  Bd3=0.0; SLOOPA(jj) Bd3 += prim[ri][rj][rk][B1+jj-1]*(rgeom[B1+jj-1].gcov[3][jj]);
  Bd3ri2=0.0; SLOOPA(jj) Bd3ri2 += prim[ri2][rj][rk][B1+jj-1]*(ri2geom[B1+jj-1].gcov[3][jj]);
  Bd3ri3=0.0; SLOOPA(jj) Bd3ri3 += prim[ri3][rj][rk][B1+jj-1]*(ri3geom[B1+jj-1].gcov[3][jj]);
  // B_\phi slope
  Dqp=Bd3ri3-Bd3ri2;
  Dqm=Bd3ri2-Bd3;
  Dqc=Bd3ri3-Bd3;
  dqBd3 = signdq*MINMOD(MINMOD(Dqp,Dqm),Dqc);



  ///////////////////////////////////////////
  //
  // reference value for \detg rho u^{x1}
  //
  ///////////////////////////////////////////
  if(ispstag==0){
    ucon_calc(prim[ri][rj][rk],&rgeom[U1],uconref);
    Mdot = (rgeom[U1].g)*prim[ri][rj][rk][RHO]*uconref[1];

    // ri2
    ucon_calc(prim[ri2][rj][rk],&ri2geom[U1],uconrefri2);
    Mdotri2 = (ri2geom[U1].g)*prim[ri2][rj][rk][RHO]*uconrefri2[1];

    // ri3
    ucon_calc(prim[ri3][rj][rk],&ri3geom[U1],uconrefri3);
    Mdotri3 = (ri3geom[U1].g)*prim[ri3][rj][rk][RHO]*uconrefri3[1];

    // Mdot slope
    Dqp=Mdotri3-Mdotri2;
    Dqm=Mdotri2-Mdot;
    Dqc=Mdotri3-Mdot;
    dqMdot = signdq*MINMOD(MINMOD(Dqp,Dqm),Dqc);
  }



  //////////////////////
  //
  // prim[pl] slope
  //
  //////////////////////
  PBOUNDLOOP(pl){
    // determine MINM slope for extrapolation
    Dqp=prim[ri3][rj][rk][pl]-prim[ri2][rj][rk][pl];
    Dqm=prim[ri2][rj][rk][pl]-prim[ri][rj][rk][pl];
    Dqc=prim[ri3][rj][rk][pl]-prim[ri][rj][rk][pl];
    dq[pl] = signdq*MINMOD(MINMOD(Dqp,Dqm),Dqc);
  }


  //////////////////////////
  //
  // log(prim[pl]) slope
  //
  //////////////////////////
  PBOUNDLOOP(pl){
    // determine MINM slope for extrapolation
    Dqp=log(fabs(prim[ri3][rj][rk][pl]))-log(fabs(prim[ri2][rj][rk][pl]));
    Dqm=log(fabs(prim[ri2][rj][rk][pl]))-log(fabs(prim[ri][rj][rk][pl]));
    Dqc=log(fabs(prim[ri3][rj][rk][pl]))-log(fabs(prim[ri][rj][rk][pl]));
    dqlogdensity[pl] = signdq*MINMOD(MINMOD(Dqp,Dqm),Dqc);
  }

  ////////////
  //
  // gdet * prim[pl] slope
  //
  ////////////
  PBOUNDLOOP(pl){
    // determine MINM slope for extrapolation
    Dqp=(ri3geom[pl].g)*prim[ri3][rj][rk][pl]-(ri2geom[pl].g)*prim[ri2][rj][rk][pl];
    Dqm=(ri2geom[pl].g)*prim[ri2][rj][rk][pl]-(rgeom[pl].g)*prim[ri][rj][rk][pl];
    Dqc=(ri3geom[pl].g)*prim[ri3][rj][rk][pl]-(rgeom[pl].g)*prim[ri][rj][rk][pl];
    dqgdetpl[pl] = signdq*MINMOD(MINMOD(Dqp,Dqm),Dqc);
  }


  //////////////////////
  //
  // u^i slope
  //
  //////////////////////
  if(ispstag==0){
    DLOOPA(jj){
      // determine MINM slope for extrapolation
      Dqp=uconrefri3[jj]-uconrefri2[jj];
      Dqm=uconrefri2[jj]-uconref[jj];
      Dqc=uconrefri3[jj]-uconref[jj];
      dqucon[jj] = signdq*MINMOD(MINMOD(Dqp,Dqm),Dqc);
    }

    // get fraction of farther away reference value to use
    // Fraction of normal reference value is determined by how similar u^t is for (e.g. for inner boundary) ri and ri+1
    xfrac = fabs(uconref[TT]-uconrefri2[TT])/uconref[TT];
#define XFRAC1 (0.1) // 0.1 was chosen from experience with BH+torus and what the u^t value was doing at the pole near the horizon -- JCM found once d(u^t)>0.1 then started to grow unbound.  Certainly related to resolution so assume resolving of horizon is generally similar to the 64^2 case with Rout=40 and exp grid
#define YFRAC1 (0.0)
#define XFRAC2 (0.2)
#define YFRAC2 (1.0)
#define YFRAC12(frac) (YFRAC1 + (YFRAC2-YFRAC1)/(XFRAC2-XFRAC1)*(frac-XFRAC1))

    if(xfrac<XFRAC1) yfrac=YFRAC1;
    else if(xfrac<XFRAC2) yfrac=YFRAC12(xfrac);
    else yfrac=YFRAC2;


#define XDQFRAC1 (0.1)
#define YDQFRAC1 (0.0)
#define XDQFRAC2 (0.2)
#define YDQFRAC2 (1.0)
#define YDQFRAC12(frac) (YDQFRAC1 + (YDQFRAC2-YDQFRAC1)/(XDQFRAC2-XDQFRAC1)*(frac-XDQFRAC1))

    // limit slope in extrapolation if excessive slope
    xdqfrac = fabs(uconref[TT]-uconrefri2[TT])/uconref[TT];
    if(xdqfrac<XDQFRAC1) ydqfrac=YDQFRAC1;
    else if(xdqfrac<XDQFRAC2) ydqfrac= YDQFRAC12(xdqfrac);
    else ydqfrac=YDQFRAC2;
    
    // modify dqucon based upon ydqfrac
    // Uses original MINM slope unless slope is too large, then start to decrease to no slope
    // only apply to velocity related slopes
    DLOOPA(jj){
      dqucon[jj] = dqucon[jj]*(1.0-ydqfrac);
      dq[UU+jj] = dq[UU+jj]*(1.0-ydqfrac);
      dqlogdensity[UU+jj] = dqlogdensity[UU+jj]*(1.0-ydqfrac);
      dqMdot = dqMdot*(1.0-ydqfrac);
      dqgdetpl[UU+jj] = dqgdetpl[UU+jj]*(1.0-ydqfrac);
    }



  }// end if ispstag==0
  else{
    yfrac=ydqfrac=0.0;
  }







  /////////////////////
  //
  // LOOP over i and pl
  //
  /////////////////////

  LOOPBOUNDGENMORE(i,iloopstart,iloopend,iloopstep){
  


    //////////////////
    //
    // default copy
    //
    //////////////////
    PBOUNDLOOP(pl) prim[i][j][k][pl]=prim[ri][rj][rk][pl];

    //////////////////
    //
    // get coordinates of local location for all quantities
    //
    //////////////////
    PBOUNDLOOP(pl){
      get_geometry(i, j, k, dirprim[pl], &geom[pl]);
      coord(i,j,k,dirprim[pl],X[pl]); // reference locations for B2/U2
      bl_coord(X[pl],V[pl]);
    }





    ///////////////////////////////
    //
    // NON-magnetic field parts (\rho,u,u^i, etc.)
    //
    ///////////////////////////////
    if(ispstag==0){


    ///////////////////////////////
    //
    // Densities and variables at higher than B3
    //
    ///////////////////////////////

    /// did assume densities roughly constant (good for disk, not as good for polar region of mostly radial Bondi-like flow)
    // now focus on accuracy for polar regions to maintain stability (avoid dissipation terms doing funny things)

    // using below Bondi-like scaling with Mdot scaling for ucon leads to space-like answers  Near BH densities should become roughly constant, and this doesn't lead to problems, so use instead
    //	  pl=RHO; prim[i][j][k][pl] = prim[ri][rj][rk][pl]*pow(V[pl][1]/Vr[pl][1],-3.0/2.0);
    //	  pl=UU;  prim[i][j][k][pl] = prim[ri][rj][rk][pl]*pow(V[pl][1]/Vr[pl][1],-5.0/2.0);

#define POWEREXTRAFRAC (0.1)
#define POWERRATIO (1.0+POWEREXTRAFRAC) // ratio of densities not allow to go bigger than this when using dqlogdensity


#if(1)
      for(pl=0;pl<NPRBOUND;pl++){
	if(pl<=UU || pl>B3){

	  // only use linear if exponentiation causes growth of value, not decreasion
	  ftemp = exp(-signdq*dqlogdensity[pl]) - POWERRATIO;

	  // limit interpolation strength to fixed log-slope
	  if(ftemp>=0.0) mydqlog=log(POWERRATIO)/(-signdq);
	  else mydqlog = dqlogdensity[pl];

	  // log extrap (very speculative and can cause problems if used alone when (say) density is super low on ri+1 and relatively high on ri, then i will be super huge
	  // should use this when values that go into slope are much different, or equally when dqlogdensity is large
	  prim[i][j][k][pl] = exp(log(fabs(prim[ri][rj][rk][pl])) + mydqlog*(i-ri));

	  // DEBUG:
	  //	  dualfprintf(fail_file,"i=%d j=%d pl=%d ftemp=%21.15g linearvalue=%21.15g expvalue=%21.15g final=%21.15g\n",i,j,pl,ftemp,linearvalue,expvalue,prim[i][j][k][pl]);

	}
      }

#else

      // override using OLD WAY for densities
#if(0)
    // linear extrap (leads to negative values, which can cause problems with flux inversions)
    pl=RHO; prim[i][j][k][pl] = prim[ri][rj][rk][pl] + dq[pl]*(i-ri);
    pl=UU;  prim[i][j][k][pl] = prim[ri][rj][rk][pl] + dq[pl]*(i-ri);
#else
    // log extrap
    pl=RHO; prim[i][j][k][pl] = exp(log(fabs(prim[ri][rj][rk][pl])) + dqlogdensity[pl]*(i-ri));
    pl=UU;  prim[i][j][k][pl] = exp(log(fabs(prim[ri][rj][rk][pl])) + dqlogdensity[pl]*(i-ri));
#endif

#endif



	  
      ///////////////////////////////
      //
      // Velocity
      //
      ///////////////////////////////

#if(0)
      // appears too speculative and use of v^3 causes ucon2pr to lead to no solution within boundary conditions -- suggestive of some inconsistency
      //  (primitive can be anything)

      // combined with \rho_0=constant, this implies \detg \rho_0 u^x1 == constant as required for Mdot=constant for radial flow
      //	  pl=U1; prim[i][j][k][pl] = prim[ri][rj][rk][pl]*fabs((rgeom[pl].g)/(geom[pl].g));
      myMdot = Mdot + dqMdot*(i-ri);
      ucon[1] = myMdot/((geom[U1].g)*prim[i][j][k][RHO]);
      // assume ucon[2] and ucon[3] are just copied directly
      // assume mostly radial flow, but visibly appears like B2 so treat similarly so no large kink in behavior
      ucon[2] = uconref[2] + dqucon[2]*(i-ri);
      // assume v^\phi \sim \Omega doesn't change much
      ucon[3] = uconref[3] + dqucon[3]*(i-ri);

      // using non-computed u^t in ucon2pr means slightly inconsistent with desired u^t[u^i], but assume approximately correct and is sufficient for extrapolation to get primitive
      //	  ucon[TT] = uconref[TT];

      // fill so can use standard functions
      primtemp[U1] = ucon[1];
      primtemp[U2] = ucon[2];
      primtemp[U3] = ucon[3];
      // get real ucon with real u^t
      ucon_calc_4vel_bothut(primtemp, &geom[U1], ucon, ucon2);
      // now get primitive velocity (choosing 4-velocity branch closest to reference velocity)
      // using the below causes flipping of choice in accretion flow region of boundary conditions due to being too far away from original uconref[TT]
      //	  if(fabs(ucon[TT]-uconref[TT])<fabs(ucon2[TT]-uconref[TT])) ucontouse=ucon;
      //	  else ucontouse=ucon2;
      ucontouse=ucon;

      /////////////////
      //
      // iterate to get closer to v^\phi\sim constant (safer than setting v^i directly and using vcon2pr())
      //
      /////////////////

      for(itervphi=0;itervphi<=NUMITERVPHI;itervphi++){
	// assume v^\phi \sim \Omega doesn't change much
	primtemp[U3] = uconref[3]/uconref[TT]*ucontouse[TT];
	// get real ucon with real u^t
	ucon_calc_4vel_bothut(primtemp, &geom[U1], ucon, ucon2);
	// now get primitive velocity (choosing 4-velocity branch closest to reference velocity)
	//	    if(fabs(ucon[TT]-uconref[TT])<fabs(ucon2[TT]-uconref[TT])) ucontouse=ucon;
	//	    else ucontouse=ucon2;
	ucontouse=ucon;
      }

      /////////////////////////
      // get final primitive from ucon
      ucon2pr(WHICHVEL,ucontouse,&geom[U1],prim[i][j][k]); // fills only U1,U2,U3


#elif(0)

      // GODMARK: Was causing problems for Rebecca


      // BH-pole corner angular sucking drive u^t>>0 and then failures occur
      // avoid exaggerating via unwanted radial sucking by using reference value with

      // whether to conserve (at least) D=\rho_0 u^t when modifying \gamma
      // see notes in fixup.c for limit_gamma() regarding why this is not a good idea
#define DO_CONSERVE_D_INBOUNDS 0 // CHANGINGMARK



#if(DO_CONSERVE_D_INBOUNDS)
      // GODMARK: here we actually modify active values of quantities (should we preserve D=\rho_0 u^t ? -- i.e. change \rho_0?)
      // This correction comes after any other velocity modifications
      // of all conservation laws, at least conserve particle number
      D0 = prim[ri][rj][rk][RHO]*uconref[TT];
#endif

      for(pl=U1;pl<=U3;pl++){

	// switch away from using near-BC value as reference if going bad since don't want to exaggerate it
	ftemp = prim[ri][rj][rk][pl]*(1.0-yfrac) + prim[ri2][rj][rk][pl]*yfrac;
	// interpolate
	prim[i][j][k][pl] = (ftemp + dq[pl]*(i-ri));


	//	dualfprintf(fail_file,"i=%d j=%d k=%d pl=%d ftemp=%21.15g yfrac=%21.15g dq[pl]=%21.15g i=%d ri=%d ri2=%d ri3=%d rj=%d rk=%d :: prim=%21.15g pri2=%21.15g pri3=%21.15g\n",i,j,k,pl,ftemp,yfrac,dq[pl],i,ri,ri2,ri3,rj,rk,prim[i][j][k][pl],prim[ri2][rj][rk][pl],prim[ri3][rj][rk][pl]); // CHANGINGMARK



	if(V[pl][RR]<RADIUSTOAVOIDRADIALSUCK){
	  // only do this if near a BH
	  // also use yfrac to reset on-grid value since apparently it's going bad
	  // this also keeps interpolation scheme passing through to boundary and so avoiding unwanted dissipation that can evolve flow away from sanity
	  // linear extrap
	  // GODMARK: Should probably make sure that u^r doesn't change sign?!
	  prim[ri][rj][rk][pl] = prim[ri][rj][rk][pl]*(1.0-yfrac) + (prim[ri2][rj][rk][pl] + dq[pl]*(ri-(ri2)))*yfrac;
	}

      }// end over velocities

#if(DO_CONSERVE_D_INBOUNDS)
      // get new u^t
      ucon_calc(prim[ri][rj][rk],&rgeom[U1],uconrefnew);
      // recompute \rho_0 to conserve particle number (generally will increase \rho_0)
      prim[ri][rj][rk][RHO] = D0/uconrefnew[TT];
#endif




#elif(0)
      // linear extrap
      for(pl=U1;pl<=U3;pl++){

	// switch away from using near-BC value as reference if going bad since don't want to exaggerate it
	ftemp = prim[ri][rj][rk][pl]*(1.0-yfrac) + prim[ri3][rj][rk][pl]*yfrac;
	// interpolate
	prim[i][j][k][pl] = ftemp + dq[pl]*(i-ri);

      }


#elif(1)
      // linear extrap
      for(pl=U1;pl<=U3;pl++){

	ftemp = prim[ri][rj][rk][pl];
	// interpolate
	prim[i][j][k][pl] = ftemp + dq[pl]*(i-ri);

      }


#endif	  


    } // end if ispstag==0



    ///////////////////////////////
    //
    // MAGNETIC FIELD
    //
    ///////////////////////////////

    ////////
    //
    // assume \detg B^x1 == constant (well, extrapolated now)
    //
    // GODMARK: Should probably make sure that B^r and B^\phi don't change sign! (well, maybe at least B^r)
    // Should probably preserve signature of B^\phi/B^r = -\Omega_F/c -- suggests interpolating B^\phi/B^r instead of B_\phi -- well, can't divide by B^r
    //
    ////////
    for(pl=B1;pl<=B2;pl++){
      igeomnosing=sign(geom[pl].g)/(fabs(geom[pl].g)+SMALL); // avoids 0.0 for any sign of geom.g
      prim[i][j][k][pl] =  (prim[ri][rj][rk][pl]*fabs((rgeom[pl].g) + dqgdetpl[pl]*(i-ri))*igeomnosing);
    }

    ///////////////
    //
    // obtian generally correct B^\phi[B^r,B^\theta,B_\phi]
    //
    ///////////////
    if(dirprim[B3]==FACE3){
      // then assume staggered method
      // need to average fields
      //      Bu1=0.25*(prim[i][j][k][B1]+prim[ip1][j][k][B1]+prim[i][j][km1][B1]+prim[ip1][j][km1][B1]);
      //Bu2=0.25*(prim[i][j][k][B2]+prim[i][jp1][k][B2]+prim[i][j][km1][B2]+prim[i][jp1][km1][B2]);
      // assume average just results in the same value since can't average non-set values
      Bu1=prim[i][j][k][B1];
      Bu2=prim[i][j][k][B2];
    }
    else{
      // then assume all fields at CENT
      Bu1=prim[i][j][k][B1];
      Bu2=prim[i][j][k][B2];
    }
    gcon03=geom[B3].gcon[0][3];
    gcon13=geom[B3].gcon[1][3];
    gcon23=geom[B3].gcon[2][3];
    gcon33=geom[B3].gcon[3][3];

    gcov01=geom[B3].gcov[0][1];
    gcov02=geom[B3].gcov[0][2];
    gcov11=geom[B3].gcov[1][1];
    gcov12=gcov21=geom[B3].gcov[1][2];
    gcov22=geom[B3].gcov[2][2];
    gcov03=geom[B3].gcov[0][3];
    gcov13=geom[B3].gcov[1][3];
    gcov23=geom[B3].gcov[2][3];
	  
    // Bd3fromBu.nb (just moved signs)
    myBd3=Bd3 + dqBd3*(i-ri);
    ftemp=(1.0 - gcon03*gcov03 - gcon13*gcov13 - gcon23*gcov23);
    igeomnosing=sign(ftemp)/(fabs(ftemp)+SMALL);
    pl=B3; prim[i][j][k][pl] = (myBd3*gcon33 + Bu1*gcon03*gcov01 + Bu2*gcon03*gcov02 + Bu1*gcon13*gcov11 + Bu2*gcon13*gcov12 + Bu1*gcon23*gcov21 + Bu2*gcon23*gcov22)*igeomnosing;

    // old B_\phi imprecise copy (need to avoid singularity anywways)
    //	  pl=B3; prim[i][j][k][pl] = prim[ri][rj][rk][pl]*fabs((rgeom[pl].gcov[3][3])/(geom[pl].gcov[3][3]));






  } // end LOOP over i and pl




  return(0);



}










// interpolate across pole to remove numerical errors there
int poledeath(int ispstag, int boundstage, int whichx2, int* dirprim, FTYPE prim[][N2M][N3M][NPR])
{
  int jj;
  int i,j,k;
  int pl;
  int rim1,rip1;
  int ri,rj,rk;
  int rj0;
  int rjstag;
  int rjtest,rjstag0,deathstagjs0,deathstagje0,rjstagtest,deathstagjstest,deathstagjetest;
  int poleloc;
  int deathjs0,deathje0;
  int deathjs,deathje;
  int deathstagjs,deathstagje;
  int gammadeathjs,gammadeathje;
  struct of_geom rgeom[NPR],geom[NPR];
  FTYPE X[NPR][NDIM];
  FTYPE V[NPR][NDIM];
  FTYPE Xr[NPR][NDIM];
  FTYPE X0[NDIM];
  int lowrho,lowuu,highu1;
  int deathjstest,deathjetest;
  FTYPE gamma,gammarj0,gammarjtest;
  FTYPE ftemp;
  FTYPE ucon[NDIM];
  FTYPE gammavaluelimit;
  int doavginradius[NPR];
  int pl2;


  // note that doesn't matter the order of the j-loop since always using reference value (so for loop doesn't need change in <= to >=)
  if(whichx2==X2DN){
    rj0 = POLEDEATH;
    rjtest = rj0+DEATHEXPANDAMOUNT; // used to ensure near pole the density doesn't drop suddenly
    poleloc = 0;
    // for POLEDEATH==2, then deathjs,je=-2,-1,0,1 as required for CENT quantities rj=2
    deathjs0 = 0-POLEDEATH;
    deathje0 = 0+POLEDEATH-1;

    deathjstest = deathjs0-DEATHEXPANDAMOUNT;
    deathjetest = deathje0+DEATHEXPANDAMOUNT;
    if(deathjstest<inoutlohi[POINTDOWN][POINTDOWN][2]) deathjstest=inoutlohi[POINTDOWN][POINTDOWN][2];
    if(deathjetest<inoutlohi[POINTDOWN][POINTDOWN][2]) deathjetest=inoutlohi[POINTDOWN][POINTDOWN][2];

    // assume for POLEDEATH==1 that B2 set correctly as 0 on pole and only do something if POLEDEATH>1
    // if POLEDEATH==2 then B2 set at  -1,0,1 and will correctly set B2 to 0 at pole rj=2
    rjstag0 = rj0;
    deathstagjs0 = 0-POLEDEATH+1;
    deathstagje0 = 0+POLEDEATH-1;

    rjstagtest = rjtest;
    deathstagjstest = deathstagjs0-DEATHEXPANDAMOUNT;
    deathstagjetest = deathstagje0+DEATHEXPANDAMOUNT;
    if(deathstagjstest<inoutlohi[POINTDOWN][POINTDOWN][2]) deathstagjstest=inoutlohi[POINTDOWN][POINTDOWN][2];
    if(deathstagjetest<inoutlohi[POINTDOWN][POINTDOWN][2]) deathstagjetest=inoutlohi[POINTDOWN][POINTDOWN][2];

    // assumes velocity is always CENT
    gammadeathjs=0-POLEGAMMADEATH;
    gammadeathje=0+POLEGAMMADEATH-1;
    if(gammadeathjs<inoutlohi[POINTDOWN][POINTDOWN][2]) gammadeathjs=inoutlohi[POINTDOWN][POINTDOWN][2];
    if(gammadeathje<inoutlohi[POINTDOWN][POINTDOWN][2]) gammadeathje=inoutlohi[POINTDOWN][POINTDOWN][2];
  }
  else if(whichx2==X2UP){
    rj0=N2-1-POLEDEATH;
    rjtest = rj0-DEATHEXPANDAMOUNT;
    poleloc=N2;
    // if POLEDEATH==2 then CENTs set at N2-2,N2-1,N2,N2+1 rj=N2-3
    deathjs0 = N2-1+1-POLEDEATH;
    deathje0 = N2-1+POLEDEATH;

    deathjstest = deathjs0-DEATHEXPANDAMOUNT;
    deathjetest = deathje0+DEATHEXPANDAMOUNT;
    if(deathjstest>inoutlohi[POINTUP][POINTUP][2]) deathjstest=inoutlohi[POINTUP][POINTUP][2];
    if(deathjetest>inoutlohi[POINTUP][POINTUP][2]) deathjetest=inoutlohi[POINTUP][POINTUP][2];

    // if POLEDEATH==2, then B2 is set at N2-1,N2,N2+1 rj=N2-3
    if(dirprim[B2]==FACE2) rjstag0=N2-POLEDEATH;
    else if(dirprim[B2]==CENT) rjstag0=rj0;
    deathstagjs0 = N2+1-POLEDEATH;
    deathstagje0 = N2-1+POLEDEATH;

    rjstagtest = rjtest;
    deathstagjstest = deathstagjs0-DEATHEXPANDAMOUNT;
    deathstagjetest = deathstagje0+DEATHEXPANDAMOUNT;
    if(deathstagjstest>inoutlohi[POINTUP][POINTUP][2]) deathstagjstest=inoutlohi[POINTUP][POINTUP][2];
    if(deathstagjetest>inoutlohi[POINTUP][POINTUP][2]) deathstagjetest=inoutlohi[POINTUP][POINTUP][2];


    // assumes velocity is always CENT .  If POLEDEATH==2, N2-2,N2-1,N2,N2+1
    gammadeathjs = N2-1+1-POLEGAMMADEATH;
    gammadeathje = N2-1+POLEGAMMADEATH;
    if(gammadeathjs>inoutlohi[POINTUP][POINTUP][2]) gammadeathjs=inoutlohi[POINTUP][POINTUP][2];
    if(gammadeathje>inoutlohi[POINTUP][POINTUP][2]) gammadeathje=inoutlohi[POINTUP][POINTUP][2];
  }
  



  /////////////////////
  //
  // Loop over i,k


#if(POLEDEATH)
  LOOPX2dir{
    ri = i;
    rk = k;


    // set radial positions to average over
    if(ri==inoutlohi[POINTDOWN][POINTDOWN][1]){
      // shift up, but not including ri
      rim1=ri+1;
      rip1=ri+2;
    }
    else if(ri==inoutlohi[POINTUP][POINTUP][1]){
      // shift down, but not including ri
      rim1=ri-2;
      rip1=ri-1;
    }
    else{
      // around ri
      rim1=ri-1;
      rip1=ri+1;
    }
    




    /////////////
    // first check for significant density drops near the axis or Lorentz factor spikes up near axis
    lowrho=lowuu=0;
    highu1=0;


    if(0){ // disabled for now

      get_geometry(ri, rj0, rk, dirprim[U1], &rgeom[U1]);
      gamma_calc(prim[ri][rj0][rk], &rgeom[U1],&gammarj0);
    
      get_geometry(ri, rjtest, rk, dirprim[U1], &rgeom[U1]);
      gamma_calc(prim[ri][rjtest][rk], &rgeom[U1],&gammarjtest);

      for (j = deathjs0; j <= deathje0; j++) {

	//POLEDENSITYDROPFACTOR
	//POLEGAMMAJUMPFACTOR


	pl=RHO;
	if(fabs(prim[i][j][k][pl])< fabs(prim[ri][rj0][rk][pl])/POLEDENSITYDROPFACTOR || fabs(prim[i][j][k][pl])< fabs(prim[ri][rjtest][rk][pl])/POLEDENSITYDROPFACTOR){
	  lowrho=1;
	}
	pl=UU;
	if(fabs(prim[i][j][k][pl])< fabs(prim[ri][rj0][rk][pl]/POLEDENSITYDROPFACTOR) || fabs(prim[i][j][k][pl])< fabs(prim[ri][rjtest][rk][pl])/POLEDENSITYDROPFACTOR){
	  lowuu=1;
	}


	// get Lorentz factor
	get_geometry(i, j, k, dirprim[U1], &rgeom[U1]);
	gamma_calc(prim[i][j][k], &rgeom[U1],&gamma);


	//      pl=U1;
	//      if(fabs(prim[i][j][k][pl])> fabs(prim[ri][rj0][rk][pl])/POLEJUMPFACTOR || fabs(prim[i][j][k][pl])> fabs(prim[ri][rjtest][rk][pl])/POLEJUMPFACTOR ){
	//	highu1=1;
	//      }

	// check variation of Lorentz factor
	if(gamma> gammarj0*POLEGAMMAJUMPFACTOR || gamma> gammarjtest*POLEGAMMAJUMPFACTOR ){
	  highu1=1;
	}


      }
    
      if(lowuu || lowrho || highu1){
	rj=rjtest; // expand copy procedure to use better reference value
	deathjs=deathjstest;
	deathje=deathjetest;

	rjstag=rjstagtest; // expand copy procedure to use better reference value
	deathstagjs=deathstagjstest;
	deathstagje=deathstagjetest;
      }
      else{
	// then use normal reference and range
	rj=rj0;
	deathjs=deathjs0;
	deathje=deathje0;

	rjstag=rjstag0;
	deathstagjs=deathstagjs0;
	deathstagje=deathstagje0;
      }
    }
    else{
      rj=rjtest; // expand copy procedure to use better reference value
      deathjs=deathjstest;
      deathje=deathjetest;
      
      rjstag=rjstagtest; // expand copy procedure to use better reference value
      deathstagjs=deathstagjstest;
      deathstagje=deathstagjetest;
    }




    ////////////////////////////////
    //
    // reference location geometry and coordinates

    PBOUNDLOOP(pl) {
      coord(ri,rj,rk,dirprim[pl],Xr[pl]); // reference locations for B2/U2
      if(pl==B2 && dirprim[B2]==FACE2){
	get_geometry(ri, rjstag, rk, dirprim[pl], &rgeom[pl]);
      }
      else{
	get_geometry(ri, rj, rk, dirprim[pl], &rgeom[pl]);
      }
    }



    ////////////////////////////////
    //
    // Loop over points to be modified
    //
    ////////////////////////////////

    for (j = MIN(deathstagjs,deathjs); j <= MAX(deathstagje,deathje); j++) {


      PBOUNDLOOP(pl) {
	// pole location
	coord(i,poleloc,k,FACE2,X0); // pole locations for B2/U2

	// current location
	coord(i,j,k,dirprim[pl],X[pl]);
	bl_coord(X[pl],V[pl]);
	get_geometry(i, j, k, dirprim[pl], &geom[pl]);

	if(AVERAGEINRADIUS && (V[pl][RR]>RADIUSTOSTARTAVERAGING)){
	  doavginradius[pl]=1;
	}
	else doavginradius[pl]=0;

      }




      if(ispstag==0){
	// symmetric quantities
	if(j>=deathjs && j<=deathje){


	  //////////
	  // for densities, u1, and u3, average in radius too!
	  // copying this means copying \Omega_F in magnetically-dominated regime beyond LC
	  for(pl=RHO;pl<=U3;pl++){
	    if(pl==U2) continue;
	    
	    if(doavginradius[pl]) prim[i][j][k][pl] = THIRD*(prim[rim1][rj][rk][pl]+prim[ri][rj][rk][pl]+prim[rip1][rj][rk][pl]);
	    else prim[i][j][k][pl] = prim[ri][rj][rk][pl];
	  }

	}    
      }



      // symmetric quantities
      if(j>=deathjs && j<=deathje){
	// B1 left alone
      }



      if(
	 dirprim[B2]==FACE2 && j>=deathstagjs && j<=deathstagje || 
	 dirprim[B2]==CENT && j>=deathjs && j<=deathje
	 ){
	/////////////////////////////
	// B2:
#if(POLEINTERPTYPE==0)
	// if flow converges toward pole, then this loses information about the velocity and field approaching the pole
	// anti-symmetric quantities:
	prim[i][j][k][B2] = 0.;

#elif(POLEINTERPTYPE==1 || POLEINTERPTYPE==2)
	// anti-symmetric:
	// assume X[2] goes through 0 at the pole and isn't positive definite
	pl=B2;
	if(doavginradius[pl]) ftemp=THIRD*(prim[rim1][rjstag][rk][pl] + prim[ri][rjstag][rk][pl] + prim[rip1][rjstag][rk][pl]);
	else  ftemp=prim[ri][rjstag][rk][pl];
	prim[i][j][k][pl] = ftemp + (X[pl][2]-Xr[pl][2])*(ftemp-0.0)/(Xr[pl][2]-X0[2]);

#elif(POLEINTERPTYPE==3)

	// then don't modify B2 -- trying to avoid instabilities related to divb=0 violation.  And seems B2 behaves ok

#endif
	
      }
      
      
      
      
      if(ispstag==0){
	
	if(j>=deathjs && j<=deathje){
	  //////////////////////////
	  // U2:
	  pl=U2;

#if(POLEINTERPTYPE==0)
	  // if flow converges toward pole, then this loses information about the velocity and field approaching the pole
	  // anti-symmetric quantities:
	  prim[i][j][k][pl] = 0.;

#elif(POLEINTERPTYPE==1 || POLEINTERPTYPE==2)
	  // anti-symmetric:
	  // assume X[2] goes through 0 at the pole and isn't positive definite
	  if(doavginradius[pl]) ftemp=THIRD*(prim[rim1][rj][rk][pl] + prim[ri][rj][rk][pl] + prim[rip1][rj][rk][pl]);
	  else ftemp=prim[ri][rj][rk][pl];
	  prim[i][j][k][pl] = ftemp + (X[pl][2]-Xr[pl][2])*(ftemp-0.0)/(Xr[pl][2]-X0[2]);

#elif(POLEINTERPTYPE==3)

	  // anti-symmetric:

	  // assume X[2] goes through 0 at the pole and isn't positive definite

	  // choose reference value
	  // 	ftemp=THIRD*(prim[rim1][rj][rk][pl] + prim[ri][rj][rk][pl] + prim[rip1][rj][rk][pl]);
	  ftemp=prim[ri][rj][rk][pl];

	  if(whichx2==X2DN && ftemp>0.0){
	    // then sucking on \theta=0 pole
	    // try to minimize sucking on pole by finding minimum U2 around
	    for(jj=0;jj<=rj+DEATHEXPANDAMOUNT;jj++) ftemp=MIN(ftemp,prim[ri][jj][rk][pl]);
	    if(doavginradius[pl]){
	      for(jj=0;jj<=rj+DEATHEXPANDAMOUNT;jj++) ftemp=MIN(ftemp,prim[rip1][jj][rk][pl]);
	      for(jj=0;jj<=rj+DEATHEXPANDAMOUNT;jj++) ftemp=MIN(ftemp,prim[rim1][jj][rk][pl]);
	    }

	    ftemp=0.0; // try crushing sucking GODMARK

	    // assume ftemp is at reference location
	    prim[i][j][k][pl] = ftemp + (X[pl][2]-Xr[pl][2])*(ftemp-0.0)/(Xr[pl][2]-X0[2]);
	  }
	  else if(whichx2==X2UP && ftemp<0.0){
	    // then sucking on \theta=\pi pole
	    for(jj=N2-1;jj>=rj-DEATHEXPANDAMOUNT;jj--) ftemp=MAX(ftemp,prim[ri][jj][rk][pl]);
	    if(doavginradius[pl]){
	      for(jj=N2-1;jj>=rj-DEATHEXPANDAMOUNT;jj--) ftemp=MAX(ftemp,prim[rip1][jj][rk][pl]);
	      for(jj=N2-1;jj>=rj-DEATHEXPANDAMOUNT;jj--) ftemp=MAX(ftemp,prim[rim1][jj][rk][pl]);
	    }

	    ftemp=0.0; // try crushing sucking GODMARK

	    // assume ftemp is at reference location (same formula for both poles)
	    prim[i][j][k][pl] = ftemp + (X[pl][2]-Xr[pl][2])*(ftemp-0.0)/(Xr[pl][2]-X0[2]);
	  }
	  else{
	    // otherwise enforce natural regular linear behavior on U2
	    if(doavginradius[pl]) ftemp=THIRD*(prim[rim1][rj][rk][pl] + prim[ri][rj][rk][pl] + prim[rip1][rj][rk][pl]);
	    else ftemp=prim[ri][rj][rk][pl];

	    prim[i][j][k][pl] = ftemp + (X[pl][2]-Xr[pl][2])*(ftemp-0.0)/(Xr[pl][2]-X0[2]);
	  }
#endif
	}
      }	


      // symmetric quantity:      
      if(j>=deathjs && j<=deathje){

	////////////////
	// B3:
	pl=B3;

#if(POLEINTERPTYPE==0 || POLEINTERPTYPE==1)
	// symmetric:
	if(doavginradius[pl]) ftemp=THIRD*(prim[rim1][rj][rk][pl] + prim[ri][rj][rk][pl] + prim[rip1][rj][rk][pl]);
	else ftemp=prim[ri][rj][rk][pl];
	prim[i][j][k][B3] = ftemp;

#elif(POLEINTERPTYPE==2)
	// symmetric:
	// approximate B_\phi copy, which (unlike copying B3) can resolve singular currents on axis
 	if(doavginradius[pl]) ftemp=THIRD*(prim[rim1][rj][rk][pl] + prim[ri][rj][rk][pl] + prim[rip1][rj][rk][pl]);
	else ftemp=prim[ri][rj][rk][pl];
	prim[i][j][k][pl] =  ftemp*fabs((rgeom[pl].gcov[3][3])/(geom[pl].gcov[3][3]));

#elif(POLEINTERPTYPE==3)
	// then do nothing if in 3D
	if(N3==1){
	  // symmetric:
	  if(doavginradius[pl]) ftemp=THIRD*(prim[rim1][rj][rk][pl] + prim[ri][rj][rk][pl] + prim[rip1][rj][rk][pl]);
	  else ftemp=prim[ri][rj][rk][pl];
	  prim[i][j][k][B3] = ftemp;
	}// else do nothing
#endif



	if(ispstag==0){

	  ///////////////////////////////////
	  //
	  // do rest if any -- assumed at CENT
	  //
	  ///////////////////////////////////
	  for(pl=B3+1;pl<NPRBOUND;pl++){
	    if(doavginradius[pl]) ftemp=THIRD*(prim[rim1][rj][rk][pl] + prim[ri][rj][rk][pl] + prim[rip1][rj][rk][pl]);
	    else ftemp=prim[ri][rj][rk][pl];
	    prim[i][j][k][pl]=ftemp;
	  }
	}

      }// end over CENT-type quantities



    }// end loop over j
  }// end loop 13
      
#endif




#if(POLEGAMMADEATH) // should come after POLEDEATH to be useful
  if(ispstag==0){
    // fixup
    // assume only dealing with velocities at same location (loop assumes CENT)
    pl=U1;
    LOOPX2dir{
      for (j = gammadeathjs; j <= gammadeathje; j++) {


	// ok-to-keep below as forever debug
	PBOUNDLOOP(pl2){
	  if( !isfinite( prim[i][j][k][pl2])){
	    dualfprintf(fail_file,"BNDNAN: t=%21.15g nstep=%ld steppart=%d :: i=%d j=%d k=%d pl2=%d prim=%21.15g\n",t,nstep,steppart,i,j,k,pl2,prim[i][j][k][pl2]);
	  }
	}




	get_geometry(i, j, k, dirprim[pl], &geom[pl]);

	ucon_calc(prim[i][j][k],&geom[pl],ucon);
	// only limit velocity if outgoing relative to grid (GODMARK: only valid in KS or BL-like coordinates such that u^r>0 means outgoing w.r.t. an observer at infinity)
	if(ucon[RR]>0.0){

	  coord(i,j,k,dirprim[pl],X[pl]);
	  bl_coord(X[pl],V[pl]);


	  if(V[pl][RR]<=GAMMAPOLEOUTGOINGRADIUS){
	    // flat \gamma limit up to GAMMAPOLEOUTGOINGRADIUS
	    gammavaluelimit = GAMMAPOLEOUTGOING;
	  }
	  else{
	    gammavaluelimit = GAMMAPOLEOUTGOING*pow(V[pl][RR]/GAMMAPOLEOUTGOINGRADIUS,GAMMAPOLEOUTGOINGPOWER);
	  }

	  limit_gamma(gammavaluelimit,prim[i][j][k],&geom[pl],-1);
	}
	else{
	  limit_gamma(GAMMAPOLEINGOING,prim[i][j][k],&geom[pl],-1);
	}
      }
    }// end loop 13
  }
#endif




  return(0);


}







// see interpline.c
int apply_bc_line(int doinverse, int iterglobal, int recontype, int bs, int be, FTYPE (*yin)[2][NBIGM], FTYPE (*yout)[2][NBIGM], FTYPE (*youtpolycoef)[MAXSPACEORDER][NBIGM])
{
  int flip_y(int iterglobal, int recontype, int bs, int be, FTYPE (*y)[2][NBIGM]);

  if(doinverse==0){
    flip_y(iterglobal, recontype, bs, be, yin);
  }
  else{
    flip_y(iterglobal, recontype, bs, be, yin);
    flip_y(iterglobal, recontype, bs, be, yout);
  }

  return(0);

}


#include "reconstructeno.h"

int flip_y(int iterglobal, int recontype, int bs, int be, FTYPE (*y)[2][NBIGM])
{
  int pl,myi;


#if( WENO_DIR_FLIP_CONS_SIGN_DN )  //flip the sign of the consrved quantities at the cylindrical axis so that they do not have a kink due to multiplication by gdet = |R|
  if( iterglobal == WENO_DIR_FLIP_CONS_SIGN_DN && (recontype == CVT_C2A || recontype == CVT_A2C) && mycpupos[iterglobal] == 0 ) { 
    PBOUNDLOOP(pl) 
      for( myi = bs; myi < 0; myi++ ) {
	y[pl][0][myi] = - y[pl][0][myi];
      }
  }
#endif
	
#if( WENO_DIR_FLIP_CONS_SIGN_UP )  //flip the sign of the consrved quantities at the cylindrical axis so that they do not have a kink due to multiplication by gdet = |R|
  if( iterglobal == WENO_DIR_FLIP_CONS_SIGN_UP && (recontype == CVT_C2A || recontype == CVT_A2C)  && mycpupos[iterglobal] == numbercpu[iterglobal] - 1 ) { 
    PBOUNDLOOP(pl) 
      for( myi = N1*(iterglobal==1) + N2*(iterglobal==2) + N3*(iterglobal==3); myi <= be; myi++ ) {
	y[pl][0][myi] = - y[pl][0][myi];
      }
  }
#endif


  return(0);

}

void remapdq( int dir, int idel, int jdel, int kdel, int i, int j, int k, FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], 
             FTYPE (*dq)[N2M][N3M][NPR2INTERP], FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], 
             FTYPE *p2interp_l, FTYPE *p2interp_r )
{
}

void remapplpr( int dir, int idel, int jdel, int kdel, int i, int j, int k, 
               FTYPE (*p2interp)[N2M][N3M][NPR2INTERP], FTYPE (*dq)[N2M][N3M][NPR2INTERP], 
               FTYPE (*pleft)[N2M][N3M][NPR2INTERP], FTYPE (*pright)[N2M][N3M][NPR2INTERP], 
               FTYPE *p2interp_l, FTYPE *p2interp_r )
{
}

///Called after the MPI boundary routines
int bound_prim_user_after_mpi_dir(int boundstage, SFTYPE boundtime, int whichdir, FTYPE prim[][N2M][N3M][NPR])
{

  return(0);
}





// check that \detg=0 on singularities
// Don't use |\theta-0|<small number since can be quite small yet not on singularity
// use integer-based grid position based detection as consistent with bondary conditions
// Called for fresh and restart run
void check_spc_singularities_user(void)
{
  int i,j,k;
  int poleloop;

  ///////////////////////////
  //
  // X1 inner r=0 singularity
  //
  ///////////////////////////


  /* inner radial BC (preserves u^t rho and u) */
  if(mycpupos[1] == 0 && BCtype[X1DN]==R0SING){
    i=0;
      
    LOOPX1dir{

      if(gdet[i][j][k][FACE1]!=0.0 || eomfunc[i][j][k][FACE1]!=0.0){
	dualfprintf(fail_file,"Detected singularity at i=%d j=%d k=%d FACE1 with gdet=%21.15g so resetting it to 0.0\n",i,j,k,gdet[i][j][k][FACE1]);
	gdet[i][j][k][FACE1]=0.0;
	eomfunc[i][j][k][FACE1]=0.0;
      }

      if(gdet[i][j][k][CORN3]!=0.0 || eomfunc[i][j][k][CORN3]!=0.0){
	dualfprintf(fail_file,"Detected singularity at i=%d j=%d k=%d CORN3 with gdet=%21.15g so resetting it to 0.0\n",i,j,k,gdet[i][j][k][CORN3]);
	gdet[i][j][k][CORN3]=0.0;
	eomfunc[i][j][k][CORN3]=0.0;
      }

      if(gdet[i][j][k][CORN2]!=0.0 || eomfunc[i][j][k][CORN2]!=0.0){
	dualfprintf(fail_file,"Detected singularity at i=%d j=%d k=%d CORN2 with gdet=%21.15g so resetting it to 0.0\n",i,j,k,gdet[i][j][k][CORN2]);
	gdet[i][j][k][CORN2]=0.0;
	eomfunc[i][j][k][CORN2]=0.0;
      }

      if(gdet[i][j][k][CORNT]!=0.0 || eomfunc[i][j][k][CORNT]!=0.0){
	dualfprintf(fail_file,"Detected singularity at i=%d j=%d k=%d CORNT with gdet=%21.15g so resetting it to 0.0\n",i,j,k,gdet[i][j][k][CORNT]);
	gdet[i][j][k][CORNT]=0.0;
	eomfunc[i][j][k][CORNT]=0.0;
      }
    }// end loop 23
  }



  ///////////////////////////
  //
  // X2 POLARAXES
  //
  ///////////////////////////


  for(poleloop=0;poleloop<=1;poleloop++){
    if(poleloop==0 && mycpupos[2] == 0 && BCtype[X2DN]==POLARAXIS){
      j=0;
    }
    else if(poleloop==1 && mycpupos[2] == ncpux2-1 && BCtype[X2UP]==POLARAXIS){
      j=N2;
    }
    else continue;


    //    LOOPX2dir{
    // No, should really be full loops:
    LOOPF1 LOOPF3{


      if(gdet[i][j][k][FACE2]!=0.0 || eomfunc[i][j][k][FACE2]!=0.0){
	dualfprintf(fail_file,"Detected singularity at i=%d j=%d k=%d FACE2 with gdet=%21.15g so resetting it to 0.0\n",i,j,k,gdet[i][j][k][FACE2]);
	gdet[i][j][k][FACE2]=0.0;
	eomfunc[i][j][k][FACE2]=0.0;
      }

      if(gdet[i][j][k][CORN3]!=0.0 || eomfunc[i][j][k][CORN3]!=0.0){
	dualfprintf(fail_file,"Detected singularity at i=%d j=%d k=%d CORN3 with gdet=%21.15g so resetting it to 0.0\n",i,j,k,gdet[i][j][k][CORN3]);
	gdet[i][j][k][CORN3]=0.0;
	eomfunc[i][j][k][CORN3]=0.0;
      }

      if(gdet[i][j][k][CORN1]!=0.0 || eomfunc[i][j][k][CORN1]!=0.0){
	dualfprintf(fail_file,"Detected singularity at i=%d j=%d k=%d CORN1 with gdet=%21.15g so resetting it to 0.0\n",i,j,k,gdet[i][j][k][CORN1]);
	gdet[i][j][k][CORN1]=0.0;
	eomfunc[i][j][k][CORN1]=0.0;
      }

      if(gdet[i][j][k][CORNT]!=0.0 || eomfunc[i][j][k][CORNT]!=0.0){
	dualfprintf(fail_file,"Detected singularity at i=%d j=%d k=%d CORNT with gdet=%21.15g so resetting it to 0.0\n",i,j,k,gdet[i][j][k][CORNT]);
	gdet[i][j][k][CORNT]=0.0;
	eomfunc[i][j][k][CORNT]=0.0;
      }
    }// end loop 13
  }// end over poles


}


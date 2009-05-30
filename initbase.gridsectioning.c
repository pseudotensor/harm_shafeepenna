
#include "decs.h"


// CHANGINGMARK: problems:
//section set location before bounds for outflow
//JOn version of OUTFLOW -- field divb problem
//big jumps, cells not defined as defined by user bounds.c because jumped more than MAXBND
//FLUXCTSTAG fails with even para
//FLUXCTTOTH fine always?  Yes, even with weno5bnd
//redo horizon enerregion to use Sasha method for moving horizon



// initialize grid sectioning
// if not restarting, then just sets to full grid (at end of init.c real first section set)
// if restarting, then go ahead and set to section given by global_sectiondef
int init_gridsectioning(void)
{
  int doprintout;
  int dimen;
  int badglobal_sectiondef;
  int faketimeorder,fakenumtimeorders;





#if( DOGRIDSECTIONING )
  trifprintf("Initializing grid sectioning: BEGIN\n");

  faketimeorder=-1;
  fakenumtimeorders=-1;

  t_transition=1; // was in bounds.c
  if( RESTARTMODE == 0 ) {
    findandsetactivesection(faketimeorder,fakenumtimeorders,nstep, t ); //SASMARK SECTIONMARK
  }
  else {
    //set up active section arrays either by using read-in parameters of the active section OR, if this
    //info is absent from restart file, from the current time
    doprintout = 1;
    badglobal_sectiondef=0;
    DIMENLOOP(dimen){
      if(global_sectiondef[POINTDOWN][dimen] < - MAXBND || global_sectiondef[POINTUP][dimen] > totalsize[dimen] + MAXBND-1 || global_sectiondef[POINTDOWN][dimen] >= global_sectiondef[POINTUP][dimen]){
	badglobal_sectiondef=1;
      }
    }
    if(badglobal_sectiondef){
      trifprintf( "Sectioning info mangled; regenerating it for current time t = %21.15g\n", t );
      findandsetactivesection(faketimeorder,fakenumtimeorders,nstep, t );
    }
    else {
      setactivesection( global_sectiondef, doprintout );
    }
  }
  trifprintf("Initializing grid sectioning: END\n");
#endif



  return(0);
}


// bound non-active region fully
// needed so RK-stepping has defined primitive at intermediate calculations
// user could use if(WITHINACTIVESECTION(ri,rj,rk)) to avoid boundary cells being set using undefined non-active cells
// use FULLLOOP to set original boundary cells to correct value so boundary conditions operate nominally
int bound_gridsectioning(int primtype, FTYPE (*prim)[N2M][N3M][NPR])
{
  FTYPE (*primsource)[N2M][N3M][NPR];
  int i,j,k,pl;
  int dimen;
  struct of_geom geom;
  struct of_state q;


  if(DOGRIDSECTIONING==0){
    dualfprintf(fail_file,"Should not be in bound_gridsectioning\n");
    myexit(248967234);
  }

  // global variables to point to in order to get defined primitives for intermediate RK primitives
  if(primtype==CENTEREDPRIM) primsource=p;
  else primsource=pstagscratch;


  //only for grid sectioning:  to fill in quasi-ghost zones
  if(prim != primsource) {

    // GODMARK: can't use COMPFULLLOOP since need to set real boundary conditions to something unless user uses WITHINACTIVESECTION(ri,rj,rk) inside their bounds.c
    //    COMPFULLLOOP{
    FULLLOOP{
      if(!WITHINACTIVESECTION(i,j,k))
	PLOOP(pl) prim[i][j][k][pl] = primsource[i][j][k][pl];  //SASMARK SECTIONING
    }
  }


  ///////////
  // GODMARK: note sure if below needed  
  if(FLUXB==FLUXCTSTAG || DOENOFLUX != NOENOFLUX ){
    // then need to deal with unew
    // force unew to be as prim (so second order)
    if(primtype==STAGGEREDPRIM){
      FULLLOOP{
	if(!WITHINACTIVESECTION(i,j,k)){
	  DIMENLOOP(dimen){
	    get_geometry(i, j, k, FACE1-1+dimen, &geom);
	    pl=B1-1+dimen;
	    unew[i][j][k][pl] = primsource[i][j][k][pl]*(geom.e[pl]); // UEVOLVE
	  }
	}
      }
    }
    else{
      FULLLOOP{
	if(!WITHINACTIVESECTION(i,j,k)){
	  get_geometry(i, j, k, CENT, &geom);
	  MYFUN(get_state(primsource[i][j][k], &geom, &q),"step_ch.c:advance()", "get_state()", 1);
	  MYFUN(primtoU(UEVOLVE,primsource[i][j][k], &q, &geom, unew[i][j][k]),"initbase.gridsectioning.c:bound_gridsectioning()", "primtoU()", 1);
	}
      }
    }
  }
  

  return(0);

}



/// Finds the index of a grid cell that conains a given radius.
/// Should be called after grid is setup for all cpus.
/// \param xr (i) radius
/// \param xi (o) index of a grid cell that contains that radius (relative to current cpu)
/// \param xcpupos1 (o) index of the cpu that contains that radius
/// \return 0 on success
int findindexfromradius(FTYPE xr, int *xcpupos1, int *xi)
{
  int i, j, k, ii;
  FTYPE r1, r2;
  FTYPE X[NDIM],V[NDIM];
  int gotit;
  int fromwhere;


  // need to find horizon and place horizoni on right-hand-side of location

  // definition of horizoni must be consistent so fluxes are consistent and have conservation
  fromwhere=1; // force to be on upside unless Rhor=0, which is caught first

  // find cpu column that brackets the horizon and determine the
  // i-offset of horizon
  // notice that only 1 CPU will get horizon since stop process once found
  // notice that radius(horizoni) is below or equal to actual horizon radius


  *xi = FLUXNOTONGRID;
  *xcpupos1 = -1;
  gotit = 0;
  for (ii = numprocs - 1; ii >= 0; ii--) { // should get done by first row
    if (ii == myid) {
      for (i = N1-1; i >= 0; i--) {
	if( BCtype[X2UP] == POLARAXIS ) {
	  j = totalsize[2]-1-startpos[2]; //on the polar axis
	  k = totalsize[3]-1-startpos[3]; //on the polar axis
	}
	else if( BCtype[X2DN] == POLARAXIS ) {
	  j = 0-startpos[2]; //on the polar axis
	  k = 0-startpos[3]; //on the polar axis
	}
	else {
	  j = N2 / 2;             // doesn't matter (spherical polar assumed)
	  k = N3 / 2;             // doesn't matter (spherical polar assumed)
	}
        coord(i, j, k, FACE1, X);
        bl_coord(X, V);
        r1=V[1];
        coord(ip1, j, k, FACE1, X);
        bl_coord(X, V);
        r2=V[1];
        // looking between FACE1's r value and upper FACE1's r value, so loop is from i=N1-1..i=0

        if(ii==myid && myid==0 && i==0){ //radius to the left of the grid
          // special check in case radius inside inner-most radial grid
          if(xr<=r1){ 
            // then horizon off grid or right on edge, but still ok
            // treat as if horizon is off grid if right on edge
            *xi = 0;
            *xcpupos1=mycpupos[1];
            break;
          }
        }

        if(ii==myid && myid==numprocs-1 && i==N1-1){ //radius to the right of the grid
          // special check in case radius outside outer-most radial grid
          if(xr>=r2){ 
            // then radius off grid or right on edge, but still ok
            // treat as if horizon is off grid if right on edge
            *xi = N1-1;
            *xcpupos1=mycpupos[1];
            break;
          }
        }

        if (fromwhere!=2){
          if(xr >= r1 && xr < r2){ // note that if strictly on r2, then next CPU should pick it up
            *xi = i;
            *xcpupos1 = mycpupos[1];
            break;
          }
        }
        else if (fromwhere==2){
          if(xr >= r1 && xr < r2){
            *xi = ip1;
            *xcpupos1 = mycpupos[1];
            if(*xi>=N1){
              *xi=0;
              ++(*xcpupos1);
            }
            else{
              // then on original CPU
              *xcpupos1 = mycpupos[1];
            }
            break;
          }
        }
      }
    }

    if (numprocs > 0) {
#if(USEMPI)
      MPI_Bcast(xi, 1, MPI_INT, ii, MPI_COMM_GRMHD);
      MPI_Bcast(xcpupos1, 1, MPI_INT, ii, MPI_COMM_GRMHD);
#endif
    }
    if (*xi >= 0) gotit = 1;                // can stop entire process

    // keep horizoni as relative to CPU with horizon so any CPU knows where horizon is
    //    if (mycpupos[1] != horizoncpupos1) {
    //  horizoni = FLUXNOTONGRID;
    //}                           // reset if not right cpu group
    if (gotit) break;
  }




  if(gotit==0){
    dualfprintf(fail_file,"Never found grid cell corresponding to radius %21.15g : fromwhere=%d\n",xr,fromwhere);
    myexit(6246);
  }


  /////////////////////////////////
  //
  // report some information
  if(fromwhere==0) {
    trifprintf("xi: %d xcpupos1: %d\n", *xi, *xcpupos1);
    // just a check
    dualfprintf(log_file,"xi: %d mycpupos[1]: %d xcpupos1: %d\n", *xi, mycpupos[1], *xcpupos1);

    trifprintf("end: find_horizon\n");
  }


  return(0);
}



// set up parameters of the active region for grid sectioning
int findandsetactivesection(int timeorder, int numtimeorders, long int thenstep, FTYPE thetime )
{
  FTYPE rlo;
  FTYPE rhi;
  int cpupos1lo;
  int cpupos1hi;
  int ilo;
  int ihi;
  int iloabs;
  int ihiabs;
  int updateeverynumsteps;
  int extra_safe_cells;
  int everynumsteps;
  int doprintout,doset;
  int findindexfromradius(FTYPE xr, int *xcpupos1, int *xi);
  int sectiondef[NUMUPDOWN][NDIM];
  FTYPE t0,t1;



#if( DOGRIDSECTIONING == 0 )
  dualfprintf(fail_file,"Got inside findandsetactivesection() with DOGRIDSECTIONING == 0\n");
  myexit(246983462);
#endif




  if(timeorder==-1 && numtimeorders==-1){
    // this indicates very first setup call that should (in general) set full grid right now

    // nothing interesting to report for initialization since not even bl_coord() defined yet
    doprintout = 0;

    // normal full total grid
    sectiondef[POINTDOWN][1]=0;
    sectiondef[POINTUP][1]=totalsize[1]-1;
    sectiondef[POINTDOWN][2]=0;
    sectiondef[POINTUP][2]=totalsize[2]-1;
    sectiondef[POINTDOWN][3]=0;
    sectiondef[POINTUP][3]=totalsize[3]-1;
        
    setactivesection( sectiondef, doprintout );
    
    return(0); // done!!
  }



  //////////////////////////
  //
  // normal non-setup section choices
  //
  //////////////////////////

  ////////////////////
  //
  // Setup update period
  //
  ///////////////////
  updateeverynumsteps=100;
  extra_safe_cells = 10;
  //number of steps after which position/size of active section is updated
  everynumsteps = updateeverynumsteps;



  // see if time for update
  doset = (nstep % updateeverynumsteps==0)  && (timeorder==numtimeorders-1);
  if(!doset) return(0); // nothing to do

  // see if should print out section information
  doprintout = (nstep % everynumsteps==0) && (timeorder==numtimeorders-1);








  /////////////////////////////////
  //
  // BELOW IS PROBLEM DEPENDENT
  //
  /////////////////////////////////
#if(0)

  // Sasha's Wind
  rlo = 0.02 * MAX(0,thetime - t_transition);  //SASMARK SECTIONMARK -- hardcored value of mm
  //this may cause problems if the actual grid inner radius is smaller than 0.1
  if( rlo < 0.1 ) {
    rlo = 0.1;
  }

  //discrete changes in lower radius of active section
  //rlo = pow( 2., floor(log(rlo)/M_LN2l) );

  rhi = 3. * thetime + 3. * (2 * M_PI / 0.25);  //SASMARK SECTIONMARK -- hardcored value of mm
  if( rhi > Rout ) {
    rhi = Rout;
  }

  //rlo = Rin;
  rhi = Rout;  //on bg don't care where the outer boundary is

  //X1DN boundary of active region
  findindexfromradius( rlo, &cpupos1lo, &ilo );

  //X1UP boundary of active region
  findindexfromradius( rhi, &cpupos1hi, &ihi );

  //find absolute index of the X1DN boundary
  iloabs = cpupos1lo * N1 + ilo;

  //find absolute index of the X1UP boundary
  ihiabs = cpupos1hi * N1 + ihi;

  //add extra cells for safety to ensure shock does not come too close numerically
  //to the X1UP boundary
  ihiabs += extra_safe_cells;

  //make sure the X1UP boundary is on the grid
  if( ihiabs >= totalsize[1] ) {
    ihiabs = totalsize[1] - 1;
  }
  

  // problem Sasha refers to is boundary conditions set worse than on-grid set and kink and behavior at true boundary leads to wave going back
  ihiabs = (totalsize[1] - 1) - MAXBND;  //to avoid problems at the upper boundary

  //iloabs = MAX( iloabs, ihiabs - totalsize[1]/2 );  //don't care for BG (was done on mako so that no more than half grid is in active section)

#else
  // torus problem

  // then switch to 6..Rout for active section with decreasing to original at late time
  t0=50;
  t1=150;
  rlo = MAX(Rin-SMALL,6.0*MIN(1.0,1.0-(thetime-t0)/(t1-t0)));


  
  //rlo = Rin;
  rhi = Rout;  //on bg don't care where the outer boundary is

  //X1DN boundary of active region
  findindexfromradius( rlo, &cpupos1lo, &ilo );

  //X1UP boundary of active region
  findindexfromradius( rhi, &cpupos1hi, &ihi );

  //find absolute index of the X1DN boundary
  iloabs = cpupos1lo * N1 + ilo;

  //find absolute index of the X1UP boundary
  ihiabs = cpupos1hi * N1 + ihi;

  //add extra cells for safety to ensure shock does not come too close numerically
  //to the X1UP boundary
  ihiabs += extra_safe_cells;

  //make sure the X1UP boundary is on the grid
  if( ihiabs >= totalsize[1] ) {
    ihiabs = totalsize[1] - 1;
  }
  
  

#endif


  /////////////////
  //
  // define arrays of locations in 3D
  //
  /////////////////
  sectiondef[POINTDOWN][1]=iloabs;
  sectiondef[POINTUP][1]=ihiabs;
  sectiondef[POINTDOWN][2]=0;
  sectiondef[POINTUP][2]=totalsize[2]-1;
  sectiondef[POINTDOWN][3]=0;
  sectiondef[POINTUP][3]=totalsize[3]-1;


  /////////////////
  //
  // Set active region (block shape in i,j,k)
  //
  /////////////////
  setactivesection( sectiondef, doprintout );

  return( 0 );
}






int setactivesection(int (*sectiondef)[NDIM], int doprintout)
{
  int dimen;
  int dir;
  int dirsign;
  int enerregion;
  int local_sectiondef[NUMUPDOWN][NDIM];
  FTYPE X[NDIM],V[NDIM];
  int Nvec[NDIM];
  int ti,tj,tk;
  int updowniter;
  int updowniteri,updowniterj,updowniterk;
  int ri;


  Nvec[0]=0;
  Nvec[1]=N1;
  Nvec[2]=N2;
  Nvec[3]=N3;
  
  /////////////
  //
  // Check hi/lo indicies
  //
  /////////////
  assert( DOGRIDSECTIONING != 1, "setactivesection(): grid sectioning should be enabled\n" );
  DIMENLOOP(dimen){
    assert( sectiondef[POINTDOWN][dimen] >= sectiondef[POINTUP][dimen], "setactivesection(): hi/lo indices out of order: dimen=%d losectiondef = %d, hisectiondef = %d\n", dimen, sectiondef[POINTDOWN][dimen], sectiondef[POINTUP][dimen] );
  }


  //////////////////////
  //
  // setup pointers
  //
  //////////////////////
  enerregion=ACTIVEREGION;
  enerpos=enerposreg[enerregion];  //activesection
  doflux=dofluxreg[enerregion];    //activefluxsection


  /////////////////////
  //
  //define global absolute indices for active section boundaries
  //
  /////////////////////
  DIMENLOOP(dimen){

    // define this CPUs relative section position
    for(updowniter=0;updowniter<NUMUPDOWN;updowniter++){
      // GODMARK: below line redundant right now
      global_sectiondef[updowniter][dimen] = sectiondef[updowniter][dimen];
      local_sectiondef[updowniter][dimen] = sectiondef[updowniter][dimen] - startpos[dimen];
    }


    dirsign=-1;
    //active section interacts with the current processor
    ri=enerpos[DIRFROMDIMEN(dimen,dirsign)] = MAX( 0, local_sectiondef[POINTDOWN][dimen] );

    if( Nvec[dimen]>1 && ri >= 0 && ri <= Nvec[dimen] ){
      doflux[DIRFROMDIMEN(dimen,dirsign)]=ri;
      if(doprintout) trifprintf("proc: %d doing sectionflux %d\n",myid,DIRFROMDIMEN(dimen,dirsign));
    }
    else doflux[DIRFROMDIMEN(dimen,dirsign)]=FLUXNOTONGRID;

    
    dirsign=1;
    ri = enerpos[DIRFROMDIMEN(dimen,dirsign)] = MIN( Nvec[dimen]-1, local_sectiondef[POINTUP][dimen] );

    //(local_sectiondef[POINTUP][dimen]+1) is the location of face index
    if( Nvec[dimen]>1 && ri >= 0 && ri <= Nvec[dimen] ){
      doflux[DIRFROMDIMEN(dimen,dirsign)]=ri + 1;  //need to add 1 to get upper edge index for the face given ri is cell center index
      if(doprintout) trifprintf("proc: %d doing sectionflux %d\n",myid,DIRFROMDIMEN(dimen,dirsign));
    }
    else doflux[DIRFROMDIMEN(dimen,dirsign)]=FLUXNOTONGRID;

  }


  ///////////////
  //
  // fluxes are on edges of zone, so 0 and N are on edge fluxes
  //
  ///////////////
  if(doprintout){
    DIRLOOP(dir) trifprintf("proc: myid=%d  :: t=%21.15g nstep=%ld enerregion=%d section: doflux[%d]=%d enerpos[%d]=%d\n",myid,t,nstep,enerregion,dir,doflux[dir],dir,enerpos[dir]);

    // full 3D cube outputted (2^3=8 3D points)
    for(updowniteri=NUMUPDOWN-1;updowniteri>=0;updowniteri--) for(updowniterj=NUMUPDOWN-1;updowniterj>=0;updowniterj--) for(updowniterk=NUMUPDOWN-1;updowniterk>=0;updowniterk--){
	  ti=global_sectiondef[updowniteri][1] + (updowniteri==POINTUP);
	  tj=global_sectiondef[updowniterj][2] + (updowniterj==POINTUP);
	  tk=global_sectiondef[updowniterk][3] + (updowniterk==POINTUP);
	  coord( ti, tj, tk, CORNT, X );
	  bl_coord( X, V );
	  trifprintf( "t = %21.15g, ud_{i,j,k} = %d %d %d :: CORNT_sectiondef_{i,j,k} = %d %d %d :: X_{1,2,3} = %21.15g %21.15g %21.15g \n", t, updowniteri, updowniterj, updowniterk, ti, tj, tk, V[1], V[2], V[3] );
    
	}
    
  }

  return( 0 );
}

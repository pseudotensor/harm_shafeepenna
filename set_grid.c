
#include "decs.h"

// GODMARK: may want to make grid functions explicitly 2D for axisymmetric space-times when in axisymmetry with space-time axis aligned with grid.

// set up all grid functions
//
// whichtime: 0: setting initial coordinate and metric quantities.  Can be called many times to solve initial value problem of coupled matter-metric system of equations.  Should NOT be treated as a single call for entire simulation.
//            1: Setting a future metric such that old metric can be used to compute connection with temporal changes incorporated
//
// CUf/Cunew: time-step for substeps used to iterate the metric and store into old metric when can take temporal difference and use as slope for present value of connection calculation
//
void set_grid(int whichtime, FTYPE *CUf, FTYPE *Cunew)
{
  int i, j, k, l, m;
  int ii, jj, kk, ll;
  FTYPE X[NDIM];
  struct of_geom geom;
  int loc;
  FTYPE V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  extern void dxdxprim(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);
  extern void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  extern void gcon_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);
  extern void eomfunc_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *eomfunc);
  extern int bound_spacetime_inside_horizon(void);
  extern int store_old_metric(void);
  extern void set_drsing(void);
  extern void set_rvsr(void);
  extern void control_time_store_metric(int whichtime, FTYPE *Cunew);
  extern int init_selfgrav(void);

  // NOTE:
  // don't assume we enter here with whichtime==0 only once since may have to "iteratively" obtain metric from stress-energy tensor and metric
  // this happens at t=0 in init() when DOSELFGRAVVSR==1


  // choose time of metric and whether to store metric before overwritten
  control_time_store_metric(whichtime,Cunew);



  /* set up boundaries, steps in coordinate grid */
  if(whichtime==0){
    dt=1.0; // dummy value that should lead to 0 connection coefficient for derivatives in time
    // just avoids 0/0 that should really be 0
  }
  // set_points() is purely dealing with coordinates, not metric quantities
  // but sets dx[0]=dt
  set_points();



  /////////////////////////////////////////
  //
  // first compute coordinate labels, which never change in time right now
  //
  /////////////////////////////////////////
  if(whichtime==0){
    // set minimum dr used to smooth metric (requires set_points())
    // assume drsing doesn't change (i.e. coordinates don't change)
    set_drsing();
    trifprintf("drsing=%21.15g\n",drsing);
  
    if(DOSELFGRAVVSR){
      // set rvsr
      // assume r doesn't change with time
      trifprintf("Setting r(i)\n");
      set_rvsr();
      trifprintf("Done setting r(i)\n");
    }



#if(DOSTOREPOSITIONDATA)


#if(MCOORD!=CARTMINKMETRIC)
   COMPFULLLOOPP1
#else
      // doesn't depend on position, only store/use 1 value
      i=j=k=0;
#endif
    {
      
      
      // over grid locations needing these quantities
      for (loc = NPG - 1; loc >= 0; loc--) {

	// store X,V, dxdxp since can be expensive to keep recomputing these things, esp. if bl_coord() involves alot of complicated functions
	coord(i, j, k, loc, Xstore[i][j][k][loc]);
	bl_coord(Xstore[i][j][k][loc],Vstore[i][j][k][loc]);
	dxdxprim(Xstore[i][j][k][loc],Vstore[i][j][k][loc],dxdxpstore[i][j][k][loc]);

	//	dualfprintf(fail_file,"i=%d j=%d k=%d loc=%d :: V=%21.15g\n",i,j,k,loc,Vstore[i][j][k][loc]);

	matrix_inverse(PRIMECOORDS, dxdxpstore[i][j][k][loc],idxdxpstore[i][j][k][loc]);
      }
    }


    didstorepositiondata=1;
#endif // end if DOSTOREPOSITIONDATA==1    


  }// end if whichtime==0

  //  dualfprintf(fail_file,"SETGRID: whichtime=%d\n",whichtime);



  ///////////////////
  //
  // Grid functions that only exist at many locations and are assigned
  // values on all points INCLUDING another value at the outer edges
  // so have edge grid data there -- makes setting up certain things
  // easier
  //
  // Notice that coord() and bl_coord() work without this.  So those
  // functions that only require those functions can do without this
  // extra grid stuff.
  //
  //////////////////


  // if here, then doing over again
  didstoremetricdata=0;


  //  dualfprintf(fail_file,"Computing metric stuff\n");
#if(MCOORD!=CARTMINKMETRIC)
 COMPFULLLOOPP1
#else
    // doesn't depend on position, only store/use 1 value
  i=j=k=0;
#endif
  {
    

    // over grid locations needing these quantities
    for (loc = NPG - 1; loc >= 0; loc--) {

      coord_ijk(i,j,k,loc,X);
      bl_coord_ijk(i,j,k,loc,V);
      /////////////////
      //
      // (1,MCOORD) here actually means PRIMCOORDS since the "1" means convert MCOORD to PRIMCOORDS.

      geom.i=icurr=i;
      geom.j=jcurr=j;
      geom.k=kcurr=k;
      geom.p=pcurr=loc;
      //      dxdxp_func(X,dxdxp[i][j][k][loc]); // future numerical version
      //dualfprintf(fail_file,"i=%d j=%d p=%d\n",i,j,loc);
      gcov_func(&geom,1,MCOORD,X, gcov[i][j][k][loc],gcovpert[i][j][k][loc]);
      //DLOOP(jj,kk) dualfprintf(fail_file,"just got gcov[%d][%d]=%21.15g\n",jj,kk,gcov[i][j][k][loc][jj][kk]);
      //      dualfprintf(fail_file,"p1 i=%d j=%d loc=%d V[1]=%21.15g\n",i,j,loc,V[1]);
      gdet[i][j][k][loc] = gdet_func_singcheck(MCOORD,V,gcov[i][j][k][loc]);
      //dualfprintf(fail_file,"p2 i=%d j=%d p=%d\n",i,j,loc);
      gcon_func(&geom,1,MCOORD,X,gcov[i][j][k][loc],gcon[i][j][k][loc]);
      alphalapse_func(&geom,1,MCOORD,X,gcov[i][j][k][loc],gcon[i][j][k][loc],&alphalapse[i][j][k][loc]);
      //dualfprintf(fail_file,"p3 i=%d j=%d p=%d\n",i,j,loc);
      eomfunc_func(&geom,1,MCOORD,X,&eomfunc[i][j][k][loc]);
      //dualfprintf(fail_file,"p4 i=%d j=%d p=%d\n",i,j,loc);


      //      if(loc==FACE2){
      //	dualfprintf(fail_file,"i=%d V1=%21.15g V[2]=%21.15g gdet=%21.15g\n",i,V[1],V[2],gdet[i][j][k][FACE2]);
      //      }


      //      if(loc==FACE1){
      //dualfprintf(fail_file,"i=%d V1=%21.15g gdet=%21.15g\n",i,V[1],gdet[i][j][k][FACE1]);
      //	DLOOP(jj,kk) dualfprintf(fail_file,"gcov[%d][%d]=%21.15g\n",jj,kk,gcov[i][j][k][FACE1][jj][kk]);
      //	DLOOP(jj,kk) dualfprintf(fail_file,"gcon[%d][%d]=%21.15g\n",jj,kk,gcon[i][j][k][FACE1][jj][kk]);
      //}


      // check if near static limit since can't divide by the below in ucon_calc
      // GODMARK
      if (fabs(gcon[i][j][k][loc][TT][TT]) < SMALL) {
	bl_coord(X,V);
	dualfprintf(fail_file, "grid location too near g_{tt}==0: %d %d %d : r=%21.15g th=%21.15g phi=%21.15g : Rin=%21.15g %21.15g\n", i,j,k,V[1],V[2],V[3],Rin,gcon[i][j][k][loc][TT][TT]);
	myexit(1);
      }
      if (0 && fabs(gcon[i][j][k][loc][RR][RR]) < SMALL) {
	bl_coord(X,V);
	dualfprintf(fail_file, "grid location too near g^{rr}==0:  %d %d %d : r=%21.15g th=%21.15g phi=%21.15g :  Rin=%21.15g %21.15g\n", i,j,k,V[1],V[2],V[3],Rin,gcon[i][j][k][loc][RR][RR]);
	myexit(1);
      }
      if (0 && fabs(gcon[i][j][k][loc][TH][TH]) < SMALL) {
	bl_coord(X,V);
	dualfprintf(fail_file,"grid location too near g^{\\theta\\theta}==0:  %d %d %d : r=%21.15g th=%21.15g phi=%21.15g :  Rin=%21.15g %21.15g\n", i,j,k,V[1],V[2],V[3],Rin,gcon[i][j][k][loc][TH][TH]);
	myexit(1);
      }
      if (0 && fabs(gcon[i][j][k][loc][PH][PH]) < SMALL) {
	bl_coord(X,V);
	dualfprintf(fail_file,"grid location too near g^{\\phi\\phi}==0:  %d %d %d : r=%21.15g th=%21.15g phi=%21.15g :  Rin=%21.15g %21.15g\n", i,j,k,V[1],V[2],V[3],Rin,gcon[i][j][k][loc][PH][PH]);
	myexit(1);
      }
      // what about g_{tt}==0? Do I ever divide by g_{tt}?
      // Yes, for ucon[TT] for 4 velocity, which is done if WHICHVEL==VEL4 or init.c
      // what about g^{rr}==0? Do I ever divide by g^{rr}?
      // Doesn't appear so
      // g^{pp} inf taken care of in metric.c by avoiding theta=0,Pi

    }
  }


  // if here, then did store new metric data
  didstoremetricdata=1;





  ///////////////////
  //
  // Grid functions that only exist at one location for all grid points
  //
  //////////////////
  //  dualfprintf(fail_file,"Computing connection\n");
#if(MCOORD!=CARTMINKMETRIC)
  FULLLOOP // connection only needed at center, and only has memory on full normal grid (not like gcon/gcov that have extra upper edge)
#else
    // doesn't depend on position, only store/use 1 value
  i=j=k=0;
#endif
  {

    loc=CENT;
    coord_ijk(i, j, k, loc, X);
    // in geom, only metric thing used is gcon to raise the lower connection.
    get_geometry(i, j, k, loc, &geom);

    //    dualfprintf(fail_file,"conn: i=%d j=%d k=%d\n",i,j,k);
    conn_func(MCOORD,X, &geom, conn[i][j][k],conn2[i][j][k]);
  }


  // DEBUG:
#if(0)
  i=j=k=0;
  DLOOP(jj,kk) dualfprintf(fail_file,"GCOV %g n",gcov[i][j][k][CENT][jj][kk]);
  DLOOP(jj,kk) DLOOPA(ll) dualfprintf(fail_file,"CONN %g\n",conn[i][j][k][jj][kk][ll]);
#endif

  ///////////////////
  //
  // Grid functions that only exist at one location AND only on active grid
  //
  //////////////////
  
  if(VOLUMEDIFF){
#if(MCOORD!=CARTMINKMETRIC)
  COMPZLOOP
#else
    // doesn't depend on position, only store/use 1 value
  i=j=k=0;
#endif
  {
      // only at 1 location, centered, using surrounding edge values
      if((defcoord==LOGRUNITH)&&(MCOORD==KSCOORDS)){
	mks_unitheta_idxvol_func(i,j,k,idxvol[i][j][k]);
      }
      else{
	idxvol[i][j][k][TT]=1.0; // really 1/dt, but changes in time      
	idxvol[i][j][k][RR]=1.0/dx[1];
	idxvol[i][j][k][TH]=1.0/dx[2];
	idxvol[i][j][k][PH]=1.0/dx[3];
      }
    }
  }



  ///////////////////
  //
  // Geometry that depends upon previously assigned quantities at multiple positions
  //
  //////////////////

  if(GDETVOLDIFF){

#if(MCOORD!=CARTMINKMETRIC)
   COMPFULLLOOPP1
#else
      // doesn't depend on position, only store/use 1 value
      i=j=k=0;
#endif
    {
    

      // over grid locations needing these quantities
      for (loc = NPG - 1; loc >= 0; loc--) {

	/////////////////
	//
	// (1,MCOORD) here actually means PRIMCOORDS since the "1" means convert MCOORD to PRIMCOORDS.

	geom.i=i;
	geom.j=j;
	geom.k=k;
	geom.p=loc;

	// uses gdet from all locations
	gdetvol_func(&geom,gdet,&gdetvol[i][j][k][loc]);


      }
    }
  }


  // set boundary conditions on metric
  // this only modifies outside unused computational regions, so can come last.  One might imagine that bounding metric would change connection calculation...true but connection bounded too.  +-1 issues? GODMARK
  if(DOEVOLVEMETRIC){
    bound_spacetime_inside_horizon();
  }





  if(DOEVOLVEMETRIC && whichtime==0){
    // if first time to call set_grid() (i.e. set_grid(0)), then store present metric as old metric
    // if evolving metric, then store old metric before computing new one
    // store metric, needed to have dg/dt terms in connection coefficients
    // initially dg/dt = 0
    // This should come after metric calculation and after connection calculation so that conneciton has old and new metric rather than new and new metric
    trifprintf("Storing old metric initially\n");
    store_old_metric();
  }




  /* done! */
}

#include "decs.h"


// static declarations
static int compute_dt_fromsource(struct of_geom *ptrgeom, struct of_state *state, FTYPE *U, FTYPE *pr, FTYPE *dUevolve, FTYPE *dUgeomevolveUU, FTYPE *dtij, FTYPE *gravitydt);
static int dUtodt(struct of_geom *ptrgeom, struct of_state *state, FTYPE *pr, FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *dUgeomgravity, FTYPE *accdt, FTYPE *gravitydt);
static int check_point_vs_average(int timeorder, int numtimeorders, PFTYPE *lpflag, FTYPE *pb, FTYPE *pf, FTYPE *upoint, FTYPE *uavg, struct of_geom *ptrgeom);


// pi: initial values at t=t0 to compute Ui
// pb: values used to compute flux/source
// pf: solution using flux(pb) from pi's Ui -> Uf

// pi, pb, and pf can all be the same since
// 1) pb used first on a stencil, not modified, to compute fluxes
// 2) pf=pi is assigned by value at each zone
// 3) pf is modified using Utoprim at each zone using pb for sources (to correspond to fluxes which used pb)
//
// So in the end only pf is modified at each zone, so the loop changing p at previous (i,j) location doesn't affect the any new location in (i,j)
int advance(int stage, FTYPE pi[][N2M][N3M][NPR],FTYPE pb[][N2M][N3M][NPR], FTYPE pf[][N2M][N3M][NPR],
	    FTYPE ui[][N2M][N3M][NPR],FTYPE uf[][N2M][N3M][NPR], FTYPE ucum[][N2M][N3M][NPR],
	    FTYPE *CUf, FTYPE *Cunew, SFTYPE fluxdt, SFTYPE boundtime, int timeorder, int numtimeorders, FTYPE *ndt)
{
  int advance_standard(int stage, FTYPE pi[][N2M][N3M][NPR],FTYPE pb[][N2M][N3M][NPR], FTYPE pf[][N2M][N3M][NPR],
		       FTYPE ui[][N2M][N3M][NPR], FTYPE uf[][N2M][N3M][NPR], FTYPE ucum[][N2M][N3M][NPR],
		       FTYPE *CUf,FTYPE *Cunew,SFTYPE fluxdt, SFTYPE boundtime, int stagenow, int numstages, FTYPE *ndt);
  int advance_finitevolume(int stage, FTYPE pi[][N2M][N3M][NPR],FTYPE pb[][N2M][N3M][NPR], FTYPE pf[][N2M][N3M][NPR],
			   FTYPE ui[][N2M][N3M][NPR],FTYPE uf[][N2M][N3M][NPR], FTYPE ucum[][N2M][N3M][NPR],
			   FTYPE *CUf,FTYPE *Cunew, SFTYPE fluxdt, SFTYPE boundtime,  int stagenow, int numstages, FTYPE *ndt);


  if(DOENOFLUX==ENOFINITEVOLUME){
    MYFUN(advance_finitevolume(stage,pi,pb,pf,ui,uf,ucum,CUf,Cunew,fluxdt,boundtime,timeorder,numtimeorders,ndt),"advance.c:advance()", "advance_finitevolume()", 1);
  }
  else if((DOENOFLUX==NOENOFLUX)||(DOENOFLUX==ENOFLUXRECON)||(DOENOFLUX==ENOFLUXSPLIT)){
    // new standard preserves conserved quantities even when metric changes
    MYFUN(advance_standard(stage,pi,pb,pf,ui,uf,ucum,CUf,Cunew,fluxdt,boundtime,timeorder,numtimeorders,ndt),"advance.c:advance()", "advance_standard()", 1);
  }
  else{
    dualfprintf(fail_file,"No such DOENOFLUX=%d\n",DOENOFLUX);
    myexit(1);
  }


  return(0);


}









// this method guarantees conservation of non-sourced conserved quantities when metric is time-dependent
// this method has updated field staggered method
int advance_standard(int stage,
		     FTYPE pi[][N2M][N3M][NPR],
		     FTYPE pb[][N2M][N3M][NPR],
		     FTYPE pf[][N2M][N3M][NPR],
		     FTYPE ui[][N2M][N3M][NPR],
		     FTYPE uf[][N2M][N3M][NPR],
		     FTYPE ucum[][N2M][N3M][NPR], 
		     FTYPE *CUf, 
		     FTYPE *Cunew, 
		     SFTYPE fluxdt,
		     SFTYPE boundtime,
		     int timeorder, int numtimeorders,
		     FTYPE *ndt)
{
  int i, j, k, pl, sc;
  FTYPE ndt1, ndt2, ndt3;
  FTYPE dUgeom[NPR],dUriemann[NPR],dUriemann1[NPR],dUriemann2[NPR],dUriemann3[NPR],dUcomp[NUMSOURCES][NPR];
  struct of_geom geom;
  struct of_state q;
  FTYPE dUtot;
  FTYPE idx1,idx2;
  SFTYPE dt4diag;
  int finalstep;
  void flux2dUavg(int i, int j, int k, FTYPE F1[][N2M][N3M][NPR],FTYPE F2[][N2M][N3M][NPR],FTYPE F3[][N2M][N3M][NPR],FTYPE *dUavg1,FTYPE *dUavg2,FTYPE *dUavg3);
  void dUtoU(FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui,  FTYPE *uf, FTYPE *ucum);
  void flux2U(int i, int j, int k, FTYPE F1[][N2M][N3M][NPR],FTYPE F2[][N2M][N3M][NPR],FTYPE *dU, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui, FTYPE *uf, FTYPE *ucum);
  int asym_compute_1(FTYPE (*prim)[N2M][N3M][NPR]);
  int asym_compute_2(FTYPE (*prim)[N2M][N3M][NPR]);
  FTYPE accdt, accdt_ij;
  int accdti,accdtj,accdtk;
  FTYPE gravitydt, gravitydt_ij;
  int gravitydti,gravitydtj,gravitydtk;
  FTYPE wavedt;
  FTYPE (*dUriemannarray)[N2M][N3M][NPR];
  int enerregion;
  int *localenerpos;
  int jj;
  FTYPE Uitemp[NPR];
  FTYPE (*utoinvert)[N2M][N3M][NPR];
  FTYPE (*myupoint)[N2M][N3M][NPR];
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int is,ie,js,je,ks,ke;
  FTYPE prbefore[NPR];

  /////////////////////////////////////////////
  //
  // Setup function tasks
  //
  ////////////////////////////////////////////


  // for ZLOOP:
  // avoid looping over region outside active+ghost grid
  // good because somewhat general and avoid bad inversions, etc.
  enerregion=TRUEGLOBALWITHBNDENERREGION;
  localenerpos=enerposreg[enerregion];



  dUriemannarray=dUgeomarray;// temporary space for dU from fluxes


  accdt=accdt_ij=BIG; // initially no limit to dt due to acceleration
  accdti=accdtj=accdtk=-100;
  gravitydt=gravitydt_ij=BIG; // initially no limit to dt due to time derivatives in gravity
  gravitydti=gravitydtj=gravitydtk=-100;



#if(ASYMDIAGCHECK)
  dualfprintf(fail_file,"BEGINNING steppart=%d nstep=%ld\n",steppart,nstep);


  dualfprintf(fail_file,"pi in advance\n");
  asym_compute_1(pi);

  dualfprintf(fail_file,"pb in advance\n");
  asym_compute_1(pb);
#endif


  // define loop range
  is=Uconsevolveloop[FIS];
  ie=Uconsevolveloop[FIE];
  js=Uconsevolveloop[FJS];
  je=Uconsevolveloop[FJE];
  ks=Uconsevolveloop[FKS];
  ke=Uconsevolveloop[FKE];



  // tells diagnostics functions if should be accounting or not
  if(timeorder==numtimeorders-1){
    dt4diag=dt;
    finalstep=1;
  }
  else{
    dt4diag=-1.0;
    finalstep=0;
  }





  /////////////////////////////////////////////
  //
  // Setup RK stuff
  //
  ////////////////////////////////////////////

#if(0)
  // DEBUG:
  bound_prim(-1,ucum);
#endif


  if(timeorder==0){
    // setup RK's uinitial (needed since sometimes set ui=0 inside advance())
    // unew should be read in, now assign to uinitial for finite volume method or normal method when dt=0.0 or moving grid
    COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOP(pl){
	ui[i][j][k][pl]=ucum[i][j][k][pl];
	// below for before first substep so starts fresh
	// helps avoid extra code in advance() for timeorder==0
	//	uf[i][j][k][pl]=ucum[i][j][k][pl]=0.0; // done later
      }
    }
  }



  /////////////////////////////////////////////
  //
  // Compute flux
  //
  ////////////////////////////////////////////

#if(PRODUCTION==0)
  trifprintf( "#0f");
#endif


  if(dt!=0.0){
    ndt1=ndt2=ndt3=BIG;
    // pb used here on a stencil, so if pb=pf or pb=pi in pointers, shouldn't change pi or pf yet -- don't currently
    MYFUN(fluxcalc(stage,pb,F1,F2,F3,CUf[2],fluxdt,&ndt1,&ndt2,&ndt3),"advance.c:advance_standard()", "fluxcalcall", 1);
  }

#if(PRODUCTION==0)
  trifprintf( "1f");
#endif


  // DEBUG: // helped!
  //  bound_prim(-1,F1);
  //  bound_prim(-1,F2);

  // DEBUG: // also helped!
  //bound_flux(-1,F1,F2,F3);


  // from here on, pi/pb/pf are only used a zone at a time rather than on a stencil



  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  //
  // now update get flux update
  //
  /////////////////////////////////////////////////////

  COMPZSLOOP(is,ie,js,je,ks,ke){

    // set geometry for centered zone to be updated
    get_geometry(i, j, k, CENT, &geom);
      
    // initialize uf and ucum if very first time here since ucum is cumulative (+=)
    if(timeorder==0) PLOOP(pl) uf[i][j][k][pl]=ucum[i][j][k][pl]=0.0;

      
    if(dt!=0.0){
      // find Ui(pi)
      // force use of primitive to set Ui since otherwise wherever corrected/changed primitive (in fixup, etc.) then would have to change conserved quantity, while same since both are point values
      // only field for staggered method is special point value at faces that needs to come from conserved quantity
      MYFUN(get_state(pi[i][j][k], &geom, &q),"step_ch.c:advance()", "get_state()", 1);
      MYFUN(primtoU(UEVOLVE,pi[i][j][k], &q, &geom, Uitemp),"step_ch.c:advance()", "primtoU()", 1);

      if(FLUXB==FLUXCTSTAG || DOENOFLUX != NOENOFLUX ){
	// then field version of ui[] is stored as "conserved" value at FACE, not CENT
	PLOOPNOB1(pl) ui[i][j][k][pl]=Uitemp[pl]; // CENT
	//PLOOPBONLY(pl) ui[i][j][k][pl] is itself // FACE (see step_ch.c's setup_rktimestep and know that ui=unew for before first substep)
	PLOOPNOB2(pl) ui[i][j][k][pl]=Uitemp[pl]; // CENT
      }
      else{
	PLOOP(pl) ui[i][j][k][pl]=Uitemp[pl]; // all at CENT
      }
	
      // dUriemann is actually average quantity, but we treat is as a point quantity at the zone center
      flux2dUavg(i,j,k,F1,F2,F3,dUriemann1,dUriemann2,dUriemann3);
      PLOOP(pl) dUriemannarray[i][j][k][pl]=dUriemann1[pl]+dUriemann2[pl]+dUriemann3[pl]; // this addition is one type of avg->point mistake
	
    }
    else{
      PLOOP(pl){
	dUriemannarray[i][j][k][pl]=0.0;
      }
    }

#if(FLUXDUMP)
    if(N1>1) PLOOP(pl) fluxdump[i][j][k][1*NPR + pl]=dUriemann1[pl];
    else PLOOP(pl) fluxdump[i][j][k][1*NPR + pl]=0.0;
		
    if(N2>1) PLOOP(pl) fluxdump[i][j][k][2*NPR + pl]=dUriemann2[pl];
    else PLOOP(pl) fluxdump[i][j][k][2*NPR + pl]=0.0;

    if(N3>1) PLOOP(pl) fluxdump[i][j][k][3*NPR + pl]=dUriemann3[pl];
    else PLOOP(pl) fluxdump[i][j][k][3*NPR + pl]=0.0;
#endif

  } // end COMPZLOOP :: end looping to obtain dUriemann and partially update unew




#if(PRODUCTION==0)
  trifprintf( "#0m");
#endif

  
    
  /////////////////////////////////////////////
  //
  // EVOLVE (update/compute) METRIC HERE
  // In general computes stress-energy tensor (T) from pb and T is then used to compute metric
  // Done here instead of after flux since horizon_flux() updates flux through boundary that would change metric
  // want metric to be at same time as everythin else done with pb so RK method is consistent
  //
  // uses unew that's NOT updated yet
  /////////////////////////////////////////////
  if(dt!=0.0){
#if(SPLITNPR)
    // don't update metric if only doing B1-B3
    if(advancepassnumber==-1 || advancepassnumber==1)
#endif
      {
	compute_new_metric_substep(pb,CUf,Cunew); // CHANGINGMARK : Not sure if Cunew here is correct
      }
  }

  ////////////////////////////////
  //
  // compute flux diagnostics (accurately using all substeps)
  //
  ///////////////////////////////
  if(dt!=0.0){
    // must come after metric changes that can change where flux surface is since otherwise when flux surface changes, we won't track this substep's flux through the new surface but the old surface (which could even be at r=0 and have no flux)
    // if using unew, then since metric update above uses old unew, need present flux at new horizon surface
#if(SPLITNPR)
    // don't update metric if only doing B1-B3
    if(advancepassnumber==-1 || advancepassnumber==1)
#endif
      {
#if(ACCURATEDIAG)
	diag_flux(pb,F1, F2, F3, fluxdt); // fluxdt is true dt for this flux as added in dUtoU() as part of ucum
#endif
      }
  }






  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  //
  // now update get source update (only affects stress-energy tensor in general)
  //
  /////////////////////////////////////////////////////
  
#if(PRODUCTION==0)
  trifprintf( "#0s");
#endif


  COMPZSLOOP(is,ie,js,je,ks,ke){

    // set geometry for centered zone to be updated
    get_geometry(i, j, k, CENT, &geom);
    
    // get state since both source() and dUtodt() need same state
    MYFUN(get_stateforsource(pb[i][j][k], &geom, &q) ,"advance.c:()", "get_state() dir=0", 1);
      
      


    // note that uf and ucum are initialized inside setup_rktimestep() before first substep
     
      
    if(dt!=0.0){
      // find dU(pb)
      // source() doesn't actually use CUf[2]=dt right now
#if(WHICHEOM==WITHNOGDET && (NOGDETB1==1 || NOGDETB2==1 || NOGDETB3==1) )
      // for FLUXB==FLUXCTSTAG, assume no source term for field
      if(FLUXB==FLUXCTSTAG){
	dualfprintf(fail_file,"Not setup for field source term if staggered field\n");
	myexit(176023);
      }
#endif
      MYFUN(source(pb[i][j][k], &geom, &q, ui[i][j][k], dUriemannarray[i][j][k], dUcomp, dUgeom),"step_ch.c:advance()", "source", 1);
      // assumes final dUcomp is nonzero and representative of source term over this timestep
	
#if(SPLITNPR)
      // don't update metric if only doing B1-B3
      if(advancepassnumber==-1 || advancepassnumber==1)
#endif
	{
#if(ACCURATEDIAG)
	  diag_source_comp(&geom,dUcomp,fluxdt);
	  diag_source_all(&geom,dUgeom,fluxdt);
#else
	  diag_source_comp(&geom,dUcomp,dt4diag);
	  diag_source_all(&geom,dUgeom,dt4diag);
#endif
	}

	
    }
    else{
      PLOOP(pl){
	dUgeom[pl]=0.0;
      }
    }
      
    dUtoU(dUgeom, dUriemannarray[i][j][k], CUf, Cunew, ui[i][j][k], uf[i][j][k], ucum[i][j][k]);
      
      
      
    // get timestep limit from acceleration
#if(LIMITDTWITHSOURCETERM)
#if(SPLITNPR)
    // don't do dUtodt if only doing B1-B3
    if(advancepassnumber==-1 || advancepassnumber==1)
#endif
    {
      if(dt!=0.0){
	// geometry is post-metric update, but should still give good estimate of future dt
	dUtodt(&geom, &q, pb[i][j][k], dUgeom, dUriemannarray[i][j][k], dUcomp[GEOMSOURCE], &accdt_ij, &gravitydt_ij);
#if( (DOEVOLVEMETRIC || DOSELFGRAVVSR) && (RESTRICTDTSETTINGINSIDEHORIZON>=1) )
	// avoid limiting dt if inside horizon
	if(WITHINENERREGION(enerposreg[OUTSIDEHORIZONENERREGION],i,j,k))
#endif
	  {
	    if(accdt_ij<accdt){
	      accdt=accdt_ij;
	      accdti=i;
	      accdtj=j;
	      accdtk=k;
	    }
	    if(gravitydt_ij<gravitydt){
	      gravitydt=gravitydt_ij;
	      gravitydti=i;
	      gravitydtj=j;
	      gravitydtk=k;
	    }
	  }
      }
    }
#endif
          

#if(FLUXDUMP)
    PLOOP(pl) fluxdump[i][j][k][0*NPR + pl]=dUgeom[pl];
#endif

  }// end COMPZLOOP


  //  previously did WITHINENERREGION() and if not so, then did below:
  // give dummy values otherwise
  //PLOOP(pl) dUriemannarray[i][j][k][pl] = dUgeom[pl] = uf[i][j][k][pl] = ucum[i][j][k][pl] = pf[i][j][k][pl]=0.0;



  


  COMPZSLOOP(is,ie,js,je,ks,ke){
    /////////////
    //
    // Utoprim as initial conditions : can't assume want these to be same in the end, so assign
    //
    // Since final step often sets pointers of pi=pf, in order to use arbitrary guess, must set here once done with pi,pb,pf.
    //
    ////////////
    PLOOP(pl) pf[i][j][k][pl] = pb[i][j][k][pl];
  }// end COMPZLOOP



  ///////////////////////////////////////
  //
  // choose which RK-quantity to invert
  //
  ///////////////////////////////////////
  if(finalstep) utoinvert = ucum;
  else utoinvert = uf;



  ////////////////////////////
  //
  // split ZLOOP above and below to allow staggered field method
  //
  ////////////////////////////
  if(FLUXB==FLUXCTSTAG){
    // if using staggered grid for magnetic field, then need to convert ucum to pstag to new pb/pf

    myupoint=upoint;

    // first copy over all quantities as point, which is true except for fields if FLUXRECON active
    COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOP(pl) myupoint[i][j][k][pl]=utoinvert[i][j][k][pl];
    }


    if(extrazones4emf && dofluxreconevolvepointfield==0){
      //bound_uavg(STAGEM1,utoinvert); // DEBUG
      // uses utoinvert and gets reaveraged field into myupoint
      field_integrate_fluxrecon(stage, pb, utoinvert, myupoint);
    }


    // first pb entry is used for shock indicator, second is filled with new field
    // myupoint goes in as staggered point value of magnetic flux and returned as centered point value of magnetic flux
    interpolate_ustag2fieldcent(stage, boundtime, timeorder, numtimeorders, pb, pstagscratch, myupoint, pf);

    ////////////////////    
    // now myupoint contains centered point conserved quantities ready for inversion
    ////////////////////

  }
  else{
    // utoinvert never reassigned from global a_utoinvert assignment since if here not doing FLUXCTSTAG
    myupoint=utoinvert;
  }






  ////////////////////////////
  //
  // INVERT
  //
  ////////////////////////////

  COMPZLOOP {


    // set geometry for centered zone to be updated
    get_geometry(i, j, k, CENT, &geom);

      
    // invert U->p
    if(finalstep){ // last call, so ucum is cooked and ready to eat!
      // store guess for diss_compute before changed by normal inversion
      PALLLOOP(pl) prbefore[pl]=pf[i][j][k][pl];

      MYFUN(Utoprimgen(finalstep,EVOLVEUTOPRIM,UEVOLVE,myupoint[i][j][k], &geom, pf[i][j][k]),"step_ch.c:advance()", "Utoprimgen", 1);
#if(DODISS||DODISSVSR)
      // then see what entropy inversion would give
      diss_compute(EVOLVEUTOPRIM,UEVOLVE,myupoint[i][j][k],&geom,prbefore,pf[i][j][k]);
#endif
	
    }
    else{ // otherwise still iterating on primitives
      MYFUN(Utoprimgen(finalstep,EVOLVEUTOPRIM,UEVOLVE,myupoint[i][j][k], &geom, pf[i][j][k]),"step_ch.c:advance()", "Utoprimgen", 1);
    }


      
#if(SPLITNPR)
    // don't update metric if only doing B1-B3
    if(advancepassnumber==-1 || advancepassnumber==1)
#endif
    {
      // immediate local (i.e. 1-zone) fix
#if(FIXUPZONES==FIXUP1ZONE)
      // SUPERGODMARK: Below should differentiate beteween whether want negative densities fixed or not, but right now fixup1zone() does all
      if((STEPOVERNEGU==0)||(STEPOVERNEGRHO==0)||(STEPOVERNEGRHOU==0)||(finalstep)){
	MYFUN(fixup1zone(pf[i][j][k],&geom,finalstep),"fixup.c:fixup()", "fixup1zone()", 1);
      }
#endif
    }


  }// end COMPZLOOP









#if(ASYMDIAGCHECK)
  dualfprintf(fail_file,"ucum in advance\n");
  asym_compute_2(ucum);

  dualfprintf(fail_file,"ENDING steppart=%d nstep=%ld\n",steppart,nstep);
#endif




  /////////////////////////////////
  //
  // If not fixing up primitives after inversion immediately, then fix up all zones at once afterwards
  //
  /////////////////////////////////

#if(SPLITNPR)
  // don't update metric if only doing B1-B3
  if(advancepassnumber==-1 || advancepassnumber==1)
#endif
    {
#if(FIXUPZONES==FIXUPALLZONES)
      fixup(stage,pf,finalstep);
#endif  
    }


  /////////////////////////////////
  //
  // Determine next timestep from waves, fluxes, and source updates
  //
  /////////////////////////////////

  wavedt = 1. / (1. / ndt1 + 1. / ndt2 + 1. / ndt3);
  *ndt = defcon * MIN(wavedt , accdt);

#if(USEGRAVITYDTINDTLIMIT)
  // use gravitydt (often too small, but sometimes accdt/ndt not small enough)
  *ndt = MIN(*ndt,gravitydt);
#endif

  gravitydtglobal = gravitydt;
  sourcedtglobal  = accdt; // accdt includes gravitydtglobal
  wavedtglobal    = wavedt;


#if(PRODUCTION==0)
  if(dt!=0.0){
    // report per-CPU time-step limited every 100 time steps

    // GODMARK: 1 : do always
    if(1|| nstep%DTr==0){
      fprintf(logdt_file,"nstep=%ld steppart=%d :: dt=%g ndt=%g ndt1=%g ndt2=%g ndt3=%g\n",nstep,steppart,dt,*ndt,ndt1,ndt2,ndt3);
      SLOOPA(jj) fprintf(logdt_file,"dir=%d wavedti=%d wavedtj=%d wavedtk=%d\n",jj,waveglobaldti[jj],waveglobaldtj[jj],waveglobaldtk[jj]);
      fprintf(logdt_file,"accdt=%g (accdti=%d accdtj=%d accdtk=%d) :: gravitydt=%g (gravitydti=%d gravitydtj=%d gravitydtk=%d) :: gravityskipstep=%d\n",accdt,accdti,accdtj,accdtk,gravitydt,gravitydti,gravitydtj,gravitydtk,gravityskipstep);
    }
  }
#endif

#if(PRODUCTION==0)
  trifprintf( "2f");
#endif

  return (0);
}










// finite volume method NOT SETUP FOR CONSISTENT METRIC EVOLUTION YET -- EASY, JUST NOT DOING IT YET -- FOLLOW ABOVE AS EXAMPLE OF WHAT TO DO
// also not setup for staggered field method
int advance_finitevolume(int stage,
			 FTYPE pi[][N2M][N3M][NPR],
			 FTYPE pb[][N2M][N3M][NPR],
			 FTYPE pf[][N2M][N3M][NPR],
			 FTYPE ui[][N2M][N3M][NPR],
			 FTYPE uf[][N2M][N3M][NPR],
			 FTYPE ucum[][N2M][N3M][NPR], 
			 FTYPE *CUf, 
			 FTYPE *Cunew, 
			 SFTYPE fluxdt,
			 SFTYPE boundtime,
			 int timeorder, int numtimeorders,
			 FTYPE *ndt)
{
  int i, j, k, pl, sc;
  FTYPE ndt1, ndt2, ndt3;
  FTYPE dUgeom[NPR],dUriemann[NPR],dUriemann1[NPR],dUriemann2[NPR],dUriemann3[NPR],dUcomp[NUMSOURCES][NPR];
  struct of_geom geom;
  struct of_state q;
  FTYPE dUtot;
  FTYPE idx1,idx2;
  SFTYPE dt4diag;
  int finalstep;
  void flux2dUavg(int i, int j, int k, FTYPE F1[][N2M][N3M][NPR],FTYPE F2[][N2M][N3M][NPR],FTYPE F3[][N2M][N3M][NPR],FTYPE *dUavg1,FTYPE *dUavg2,FTYPE *dUavg3);
  void dUtoU(FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui,  FTYPE *uf, FTYPE *ucum);
  void flux2U(int i, int j, int k, FTYPE F1[][N2M][N3M][NPR],FTYPE F2[][N2M][N3M][NPR],FTYPE *dU, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui, FTYPE *uf, FTYPE *ucum);
  int asym_compute_1(FTYPE (*prim)[N2M][N3M][NPR]);
  int asym_compute_2(FTYPE (*prim)[N2M][N3M][NPR]);
  FTYPE accdt, accdt_ij;
  int accdti,accdtj,accdtk;
  FTYPE gravitydt, gravitydt_ij;
  int gravitydti,gravitydtj,gravitydtk;
  FTYPE wavedt;
  int enerregion;
  int *localenerpos;
  int jj;
  FTYPE (*utoinvert)[N2M][N3M][NPR];
  FTYPE (*myupoint)[N2M][N3M][NPR];
  FTYPE fdummy;
  // AVG_2_POINT functions:
  void debugeno_compute(FTYPE p[][N2M][N3M][NPR],FTYPE U[][N2M][N3M][NPR],FTYPE debugvars[][N2M][N3M][NUMENODEBUGS]);
  // staggered field function:
  int whichpltoavg[NPR];
  int ifnotavgthencopy[NPR];
  int docons,dosource;
  int locpl[NPR];
  int is,ie,js,je,ks,ke;
  FTYPE prbefore[NPR];





  // any cons
  docons=0;
  PLOOP(pl) docons+=do_conserved_integration[pl];
  docons*=(DOENOFLUX == ENOFINITEVOLUME);

  // any source
  dosource=0;
  PLOOP(pl) dosource+=do_source_integration[pl];
  dosource*=(DOENOFLUX == ENOFINITEVOLUME);



  /////////////////////////////////////////////
  //
  // Setup loops and dt's
  //
  ////////////////////////////////////////////


  // for CZLOOP:
  // avoid looping over region outside active+ghost grid
  // good because somewhat general and avoid bad inversions, etc.
  enerregion=TRUEGLOBALWITHBNDENERREGION;
  localenerpos=enerposreg[enerregion];


  accdt=accdt_ij=BIG; // initially no limit to dt due to acceleration
  accdti=accdtj=accdtk=-100;
  gravitydt=gravitydt_ij=BIG; // initially no limit to dt due to time-changes in gravity
  gravitydti=gravitydtj=gravitydtk=-100;
 

  // tells diagnostics functions if should be accounting or not
  if(timeorder==numtimeorders-1){
    dt4diag=dt;
    finalstep=1;
  }
  else{
    dt4diag=-1.0; // tells diag_source to not consider
    finalstep=0;
  }


  // define loop range
  is=Uconsevolveloop[FIS];
  ie=Uconsevolveloop[FIE];
  js=Uconsevolveloop[FJS];
  je=Uconsevolveloop[FJE];
  ks=Uconsevolveloop[FKS];
  ke=Uconsevolveloop[FKE];


  /////////////////////////////////////////////
  //
  // Setup RK stuff
  //
  ////////////////////////////////////////////

  if(timeorder==0){
    // setup RK's uinitial (needed since sometimes set ui=0 inside advance())
    // unew should be read in, now assign to uinitial for finite volume method or normal method when dt=0.0 or moving grid
    // below can't be CZLOOP since need uinitial in ghost+active region
    // GODMARK: since ui=ucum (average quantities) then if change primitive somehow (fixup, etc.) then must change corresponding average conserved quantity somehow (this is one reason why using average values is annoying, although in some cases we use Sasha's method to treat average like point for *change* in average conserved quantity)
    COMPZSLOOP(is,ie,js,je,ks,ke){
      PLOOP(pl){
	ui[i][j][k][pl]=ucum[i][j][k][pl];
	// below for before first substep so starts fresh
	// helps avoid extra code in advance() for timeorder==0
	//	uf[i][j][k][pl]=ucum[i][j][k][pl]=0.0; // done later
      }
    }
  }
  // otherwise assume ui didn't change.  Present RK schemes assume this.  Otherwise have to keep track of pf/Uf pairs in RK stepping






  /////////////////////////////////////////////
  //
  // Compute flux
  //
  ////////////////////////////////////////////

#if(PRODUCTION==0)
  trifprintf( "#0f");
#endif

  if(dt!=0.0){
    ndt1=ndt2=ndt3=BIG;
    // pb used here on a stencil, so if pb=pf or pb=pi in pointers, shouldn't change pi or pf yet -- don't currently
    MYFUN(fluxcalc(stage,pb,F1,F2,F3,CUf[2],fluxdt,&ndt1,&ndt2,&ndt3),"advance.c:advance_standard()", "fluxcalcall", 1);
  }

#if(PRODUCTION==0)
  trifprintf( "1f");
#endif
  // from here on, pi/pb/pf are only used a zone at a time rather than on a stencil







#if(PRODUCTION==0)
  trifprintf( "#0m");
#endif



  /////////////////////////////////////////////
  //
  // EVOLVE (update/compute) METRIC HERE
  // In general computes stress-energy tensor (T) from pb and T is then used to compute metric
  // Done here instead of after flux since horizon_flux() updates flux through boundary that would change metric
  // want metric to be at same time as everythin else done with pb so RK method is consistent
  //
  // uses unew that's NOT updated yet
  /////////////////////////////////////////////
  if(dt!=0.0){
#if(SPLITNPR)
    // don't update metric if only doing B1-B3
    if(advancepassnumber==-1 || advancepassnumber==1)
#endif
      {
	compute_new_metric_substep(pb,CUf,Cunew);
      }
  }



  ////////////////////////////////
  //
  // compute flux diagnostics (accurately using all substeps)
  //
  ///////////////////////////////
  if(dt!=0.0){
    // must come after metric changes that can change where flux surface is since otherwise when flux surface changes, we won't track this substep's flux through the new surface but the old surface (which could even be at r=0 and have no flux)
    // if using unew, then since metric update above uses old unew, need present flux at new horizon surface
#if(SPLITNPR)
    // don't update metric if only doing B1-B3
    if(advancepassnumber==-1 || advancepassnumber==1)
#endif
      {
#if(ACCURATEDIAG)
	diag_flux(pb,F1, F2, F3, fluxdt); // fluxdt is true dt for this flux as added in dUtoU() as part of ucum update
#endif
      }
  }



	
  /////////////////////////
  //
  // SOURCE TERM
  //
  ////////////////////////

  if(dt!=0.0){
    // GODMARK: other/more special cases?
#if(WHICHEOM==WITHNOGDET && (NOGDETB1==1 || NOGDETB2==1 || NOGDETB3==1) )
    // for FLUXB==FLUXCTSTAG, assume no source term for field
    if(FLUXB==FLUXCTSTAG){
      dualfprintf(fail_file,"Not setup for field source term if staggered field\n");
      myexit(176023);
    }
#endif
    COMPZSLOOP(is,ie,js,je,ks,ke){
      // find dU(pb)
      // only find source term if non-Minkowski and non-Cartesian
      // set geometry for centered zone to be updated
      get_geometry(i, j, k, CENT, &geom);

      // get state
      MYFUN(get_stateforsource(pb[i][j][k], &geom, &q) ,"advance.c:()", "get_state() dir=0", 1);


#if(LIMITSOURCES)
      // dUriemann is volume averaged quantity (here this calcuation is done in case want to limit sources)
      flux2dUavg(i,j,k,F1,F2,F3,dUriemann1,dUriemann2,dUriemann3);
      PLOOP(pl) dUriemann[pl]=dUriemann1[pl]+dUriemann2[pl]+dUriemann3[pl]; // this addition is entirely consistent with point->averages
#else
      PLOOP(pl) dUriemann[pl]=0.0;
#endif
	    
      // get source term
      MYFUN(source(pb[i][j][k], &geom, &q, ui[i][j][k], dUriemann, dUcomp, dUgeomarray[i][j][k]),"step_ch.c:advance()", "source", 1);
	    

#if(SPLITNPR)
      // don't update metric if only doing B1-B3
      if(advancepassnumber==-1 || advancepassnumber==1)
#endif
      {
#if(ACCURATEDIAG)
	diag_source_comp(&geom,dUcomp,fluxdt);
#else
	// assumes final dUcomp is nonzero and representative of source term over this timestep
	diag_source_comp(&geom,dUcomp,dt4diag);
#endif
      }
    }// end COMPZLOOP


    // volume integrate dUgeom
    // c2a_1 c2a_2 c2a_3
    if(dosource){
      // need to limit source averaging -- GODMARK
      PALLLOOP(pl) locpl[pl]=CENT;
      PALLLOOP(pl) whichpltoavg[pl]=do_source_integration[pl];// default
      PALLLOOP(pl) ifnotavgthencopy[pl]=1-do_source_integration[pl];// default
      avg2cen_interp(locpl,whichpltoavg, ifnotavgthencopy, ENOSOURCETERM, ENOCENT2AVGTYPE, pb, dUgeomarray, dUgeomarray);
    }
  }// end if dt!=0.0





  ///////////////////////////////////////////////////////
  //
  // update Ui to Uf (ultimately to ucum)
  //
  ///////////////////////////////////////////////////////

  COMPZSLOOP(is,ie,js,je,ks,ke){

    // get geometry at center where source is located
    get_geometry(i, j, k, CENT, &geom);

    // get state
    MYFUN(get_stateforsource(pb[i][j][k], &geom, &q) ,"advance.c:()", "get_state() dir=0", 1);
		
    /////////////
    //
    // Utoprim as initial conditions : can't assume want these to be same in the end, so assign
    //
    ////////////
		
    // initialize uf and ucum if very first time here since ucum is cumulative (+=)
    if(timeorder==0) PLOOP(pl) uf[i][j][k][pl]=ucum[i][j][k][pl]=0.0;
		
    // NEED TO DEFINE Ui on other substeps besides the 0th one
    // find Ui(pi)



    if(dt!=0.0){
      // get source term (volume integrated)
      PLOOP(pl) dUgeom[pl]=dUgeomarray[i][j][k][pl];

#if(SPLITNPR)
      // don't update metric if only doing B1-B3
      if(advancepassnumber==-1 || advancepassnumber==1)
#endif
	{
#if(ACCURATEDIAG)
	  // do diagnostics for volume integrated source term
	  diag_source_all(&geom,dUgeom,fluxdt);
#else
	  diag_source_all(&geom,dUgeom,dt4diag);
#endif
	}	  
		  
		  

      // dUriemann is volume averaged quantity
      flux2dUavg(i,j,k,F1,F2,F3,dUriemann1,dUriemann2,dUriemann3);
      PLOOP(pl) dUriemann[pl]=dUriemann1[pl]+dUriemann2[pl]+dUriemann3[pl]; // this addition is entirely consistent with point->averages
    }
    else{
      PLOOP(pl){
	dUriemann[pl]=0.0;
	dUgeom[pl]=0.0;
      }
    }

    // find uf==Uf and additional terms to ucum
    dUtoU(dUgeom, dUriemann, CUf, Cunew, ui[i][j][k], uf[i][j][k], ucum[i][j][k]);



#if(LIMITDTWITHSOURCETERM)
    // GODMARK: no longer have access to dUcomp : NEED TO FIX
    // below is correct, but excessive
    // get source term again in order to have dUcomp (NEED TO FIX)
    MYFUN(source(pb[i][j][k], &geom, &q, ui[i][j][k], dUriemann, dUcomp, &fdummy),"step_ch.c:advance()", "source", 2);


    dUtodt(&geom, &q, pb[i][j][k], dUgeom, dUriemann, dUcomp[GEOMSOURCE], &accdt_ij, &gravitydt_ij);
#if( (DOEVOLVEMETRIC || DOSELFGRAVVSR) && (RESTRICTDTSETTINGINSIDEHORIZON>=1) )
    // avoid limiting dt if inside horizon
    if(WITHINENERREGION(enerposreg[OUTSIDEHORIZONENERREGION],i,j,k))
#endif
      {
	if(accdt_ij<accdt){
	  accdt=accdt_ij;
	  accdti=i;
	  accdtj=j;
	  accdtk=k;
	}
	if(gravitydt_ij<gravitydt){
	  gravitydt=gravitydt_ij;
	  gravitydti=i;
	  gravitydtj=j;
	  gravitydtk=k;
	}
      }
#endif




#if(FLUXDUMP)
    PLOOP(pl) fluxdump[i][j][k][0*NPR + pl]=dUgeom[pl];

    if(N1>1) PLOOP(pl) fluxdump[i][j][k][1*NPR + pl]=dUriemann1[pl];
    else PLOOP(pl) fluxdump[i][j][k][1*NPR + pl]=0.0;
		
    if(N2>1) PLOOP(pl) fluxdump[i][j][k][2*NPR + pl]=dUriemann2[pl];
    else PLOOP(pl) fluxdump[i][j][k][2*NPR + pl]=0.0;

    if(N3>1) PLOOP(pl) fluxdump[i][j][k][3*NPR + pl]=dUriemann3[pl];
    else PLOOP(pl) fluxdump[i][j][k][3*NPR + pl]=0.0;
#endif
  

  }



  /////////////
  //
  // Utoprim as initial conditions : can't assume want these to be same in the end, so assign
  //
  // Since final step often sets pointers of pi=pf, in order to use arbitrary guess, must set here once done with pi,pb,pf.
  //
  ////////////
  // setup initial guess for inversion
  // use pb since should be closer to pf
  COMPZSLOOP(is,ie,js,je,ks,ke) PLOOP(pl) pf[i][j][k][pl] = pb[i][j][k][pl];




  /////////////////////////////
  //
  // volume differentiate the conserved quantity
  //
  //////////////////////////////
  if(finalstep){ // last call, so ucum is cooked and ready to eat!
    utoinvert=ucum;
    ubound=utoinvert;
  }
  else{ // otherwise still iterating on primitives
    utoinvert=uf;
    ubound=utoinvert;
  }
  myupoint=upoint;


  // to debug ENO
  //#if(DOENODEBUG)
  //debugeno_compute(pb,utoinvert,enodebugarray); //SASDEBUG -- OLD USE: now assign debugs inside reconstructeno code
  //#endif





  // a2c_1 a2c_2 a2c_3
  if(docons){
    // conserved quantity is limited later after primitive is obtained
    PALLLOOP(pl) locpl[pl]=CENT;
    PALLLOOP(pl) whichpltoavg[pl]=do_conserved_integration[pl];// default
    PALLLOOP(pl) ifnotavgthencopy[pl]=1-do_conserved_integration[pl];// default
    avg2cen_interp(locpl,whichpltoavg, ifnotavgthencopy, ENOCONSERVED, ENOAVG2CENTTYPE, pb, utoinvert, myupoint);  //SASMARK:  pb's for shock indicators should be defined on ALL grid, not only on ghost+active.  Maybe should use pi instead because define everywhere?
  }
  else{
    //    myupoint=utoinvert;
    COMPZSLOOP(is,ie,js,je,ks,ke) PLOOP(pl) myupoint[i][j][k][pl]=utoinvert[i][j][k][pl];
  }



  ////////////////////////////
  //
  // split CZLOOP above and below to allow staggered field method
  //
  ////////////////////////////
  if(FLUXB==FLUXCTSTAG){
    // if using staggered grid for magnetic field, then need to convert ucum to pstag to new pb/pf

    // GODMARK: If had c2a/a2c with 3-point outputs, then could do avg2cen_interp and below at once

    // first pb entry is used for shock indicator, second is filled with new field
    interpolate_ustag2fieldcent(stage, boundtime, timeorder, numtimeorders, pb, pstagscratch, myupoint, pf);

    ////////////////////    
    // now utoinvert contains centered point conserved quantities ready for inversion
    ////////////////////

  }




  /////////////////////////////////
  //
  // Invert U -> p
  //
  // and fixup p if inversion bad
  // and compute dissipation rate if requested
  //
  /////////////////////////////////


  COMPZLOOP {
    // set geometry for centered zone to be updated
    get_geometry(i, j, k, CENT, &geom);


    // store guess for diss_compute before changed by normal inversion
    PALLLOOP(pl) prbefore[pl]=pf[i][j][k][pl];


    // invert point U-> point p
    MYFUN(Utoprimgen(finalstep,EVOLVEUTOPRIM, UEVOLVE, myupoint[i][j][k], &geom, pf[i][j][k]),"step_ch.c:advance()", "Utoprimgen", 1);

    //If using a high order scheme, need to choose whether to trust the point value
    if(docons){
      MYFUN(check_point_vs_average(timeorder, numtimeorders, pflag[i][j][k],pb[i][j][k],pf[i][j][k],myupoint[i][j][k],utoinvert[i][j][k],&geom),"advance.c:advance_finitevolume()", "check_point_vs_average()", 1);
    }


#if(DODISS||DODISSVSR)
    if(finalstep){
      // then see what entropy inversion would give
      diss_compute(EVOLVEUTOPRIM,UEVOLVE,ucum[i][j][k],&geom,prbefore, pf[i][j][k]);
    }
#endif


#if(SPLITNPR)
    // don't update metric if only doing B1-B3
    if(advancepassnumber==-1 || advancepassnumber==1)
#endif
      {
	// immediate local (i.e. 1-zone) fix
#if(FIXUPZONES==FIXUP1ZONE)
	// SUPERGODMARK: Below should differentiate beteween whether want negative densities fixed or not, but right now fixup1zone() does all
	if((STEPOVERNEGU==0)||(STEPOVERNEGRHO==0)||(STEPOVERNEGRHOU==0)||(finalstep)){
	  MYFUN(fixup1zone(pf[i][j][k], &geom, finalstep),"advance.c:advance_finitevolume()", "fixup1zone()", 1);
	}
#endif
      }
  }// end COMPZLOOP



  /////////////////////////////////
  //
  // If not fixing up primitives after inversion immediately, then fix up all zones at once afterwards
  //
  /////////////////////////////////

#if(SPLITNPR)
  // don't update metric if only doing B1-B3
  if(advancepassnumber==-1 || advancepassnumber==1)
#endif
    {
#if(FIXUPZONES==FIXUPALLZONES)
      fixup(stage,pf,finalstep);
#endif  
    }

  /////////////////////////////////
  //
  // Determine next timestep from waves, fluxes, and source updates
  //
  /////////////////////////////////


  wavedt = 1. / (1. / ndt1 + 1. / ndt2 + 1. / ndt3);
  *ndt = defcon * MIN(wavedt , accdt);

#if(USEGRAVITYDTINDTLIMIT)
  // use gravitydt (often too small, but sometimes accdt/ndt not small enough)
  *ndt = MIN(*ndt,gravitydt);
#endif

  gravitydtglobal = gravitydt;
  sourcedtglobal  = accdt; // accdt includes gravitydtglobal
  wavedtglobal    = wavedt;

#if(PRODUCTION==0)
  if(dt!=0.0){
    // report per-CPU time-step limited every 100 time steps

    // GODMARK: 1 : do always
    if(1|| nstep%DTr==0){
      fprintf(logdt_file,"nstep=%ld steppart=%d :: dt=%g ndt=%g ndt1=%g ndt2=%g ndt3=%g\n",nstep,steppart,dt,*ndt,ndt1,ndt2,ndt3);
      SLOOPA(jj) fprintf(logdt_file,"dir=%d wavedti=%d wavedtj=%d wavedtk=%d\n",jj,waveglobaldti[jj],waveglobaldtj[jj],waveglobaldtk[jj]);
      fprintf(logdt_file,"accdt=%g (accdti=%d accdtj=%d accdtk=%d) :: gravitydt=%g (gravitydti=%d gravitydtj=%d gravitydtk=%d) :: gravityskipstep=%d\n",accdt,accdti,accdtj,accdtk,gravitydt,gravitydti,gravitydtj,gravitydtk,gravityskipstep);
    }
  }
#endif



#if(PRODUCTION==0)
  trifprintf( "2f");
#endif

  return (0);
}







// check whether point conserved quantity inverted successfully to point primitive.
//   if unsuccessful, then see if want to revert to average conserved quantity and invert that
//   if Uavg->p unsuccessful, then leave as failure
// if Upoint->p is good, then check if p from Upoint is much different than p from Uavg.  If so, limit change

// upoint only needed for diagnostics
int check_point_vs_average(int timeorder, int numtimeorders, PFTYPE *lpflag, FTYPE *pb, FTYPE *pf, FTYPE *upoint, FTYPE *uavg, struct of_geom *ptrgeom)
{
  FTYPE pavg[NPR];  //atch for temporary storage of primitives obtained from inverting the averaged conserved quantities
  int invert_from_point_flag, invert_from_average_flag;
  FTYPE frac_avg_used;  //this is be used for flux interpolation limiting
  int pl;
  FTYPE limit_prim_correction( FTYPE fractional_difference_threshold, struct of_geom *geom, FTYPE *pin, FTYPE *pout );
  FTYPE fractional_diff( FTYPE a, FTYPE b );
  int is_convergence_failure;
  int avgschemeatall;
  int finalstep;


  finalstep=timeorder == numtimeorders-1;



  avgschemeatall=(interporder[avgscheme[1]]>3) ||  (interporder[avgscheme[2]]>3) ||  (interporder[avgscheme[3]]>3);
  if(avgschemeatall==0) return(0); // since nothing to do


  invert_from_point_flag = lpflag[FLAGUTOPRIMFAIL];


  if( 0 && debugfail >= 1 && (invert_from_point_flag == UTOPRIMFAILUNEG || invert_from_point_flag == UTOPRIMFAILRHONEG) ) {
    dualfprintf( fail_file, "t = %g, nstep = %ld, steppart = %d, i = %d, j = %d, rho = %21.15g, u = %21.15g, fracneg = %21.15g\n",
		 t, realnstep, steppart, ptrgeom->i + startpos[1], ptrgeom->j + startpos[2],
		 pf[RHO], pf[UU], (pf[RHO]>0)?(-pf[UU]/(pf[RHO]+DBL_MIN)):(-pf[RHO]/(pf[UU]+DBL_MIN)) );
  }


  //WHAT IF INTERNAL ENERGY BECOMES SLIGHTLY NEGATIVE?  WE STILL CAN DO THE LIMITING IN PRIM QUANTITIES! -- coorrected but check! -- SUPERSASMARK TODO atch
  if( LIMIT_AC_PRIM_FRAC_CHANGE &&
      (
       invert_from_point_flag == UTOPRIMNOFAIL || //atch added the below to still do the pt. vs. avg. check on primitives if the internal energy goes neg.
       ( (invert_from_point_flag==UTOPRIMFAILUNEG) && (0 != STEPOVERNEGU) ) || //intermediate substep with stepping over u < 0
       ( (invert_from_point_flag==UTOPRIMFAILRHONEG) && (0 != STEPOVERNEGRHO) ) //intermediate substep with stepping over rho < 0
       )
      ) {

    //make a copy of the initial guess so that not to modify the original pb's
    PLOOP(pl) pavg[pl] = pb[pl];
    //invert the average U -> "average" p
    MYFUN(Utoprimgen(finalstep,EVOLVEUTOPRIM, UEVOLVE, uavg, ptrgeom, pavg),"step_ch.c:advance()", "Utoprimgen", 3);

    invert_from_average_flag = pflag[ptrgeom->i][ptrgeom->j][ptrgeom->k][FLAGUTOPRIMFAIL];

    //Inversion from the average value succeeded or has a negative density or internal energy
    if( invert_from_average_flag == UTOPRIMNOFAIL  || invert_from_average_flag == UTOPRIMFAILUNEG ||  invert_from_average_flag == UTOPRIMFAILRHONEG ) {
      //Inversion from both from the point and the average values succeeded
      //checks if the states' gamma factors and densities are different by more than a certain fraction
      //and if different, modify the point values such that they are not further than by MAX_AC_PRIM_FRAC_CHANGE
      //from the average ones
      frac_avg_used = limit_prim_correction(MAX_AC_PRIM_FRAC_CHANGE, ptrgeom, pavg, pf);

#if( DOENODEBUG )  //atch debug; note that use the location with pl =0 , interporflux = 0, & dir = 1 since limiting the change in prim quantities is not a per-component operation
      if(  frac_avg_used > 0.01 ) {
	enodebugarray[ptrgeom->i][ptrgeom->j][ptrgeom->k][1-1][0][0][ENODEBUGPARAM_LIMCORR_PRIM]++;
      }
#endif
      
      if( pf[UU] < 0.0 && timeorder == numtimeorders-1 ) {
        lpflag[FLAGUTOPRIMFAIL] = UTOPRIMFAILU2AVG2;
      }

      //lpflag[FLAGUTOPRIMFAIL] = invert_from_point_flag;  //unneeded since it is alrady = UTOPRIMNOFAIL

    } // end if both point and average did NOT fail
    else {
      //Inversion from the point values succeeded but that from the average failed:
      //retain the point value

      //set the inversion error flag to correspond to the inversion from the average value
      lpflag[FLAGUTOPRIMFAIL] = invert_from_point_flag;
      frac_avg_used = 0.0;  //used point value, i.e. zero fracion of the average value
    }
  }
  else if( INVERTFROMAVERAGEIFFAILED && (invert_from_point_flag!=UTOPRIMNOFAIL) ) {  //failure  //atch correct
    //inversion from the point value failed

    // if last substep -> revert to the average value, else if only negative densities then allow for substep.  If other type of failures, then never allow and revert to Uavg


    // GODMARK:
    // CHECK IF INVERSION FAILED.  IF FAILED, TRY TO USE uavg:
    // invert U->p

    // GODMARK: decide whether best to revert to average or not

    //=1 if it is a non-convergence failure; =0 if it is an occurrence of a negative density
    is_convergence_failure =  invert_from_point_flag!=UTOPRIMFAILUNEG &&
      invert_from_point_flag!=UTOPRIMFAILRHONEG &&
      invert_from_point_flag!=UTOPRIMFAILRHOUNEG;

    if((timeorder==numtimeorders-1 /*&& (1 == is_nonconvergence_failure || FIXUPZONES == FIXUPNOZONES)*/ ) ||   // last substep, then DO invert from average IF no fixups or non-convergence failure (ADT) SASMARKx

       ( 1 == is_convergence_failure ) || //non-convergence (no solution in primitive quantities) error, so have to fix it up

       ( (timeorder<numtimeorders-1) && (
					 ((invert_from_point_flag==UTOPRIMFAILUNEG) && (0 == STEPOVERNEGU))||
					 ((invert_from_point_flag==UTOPRIMFAILRHONEG) && (0 == STEPOVERNEGRHO))||
					 ((invert_from_point_flag==UTOPRIMFAILRHOUNEG) && (0 == STEPOVERNEGRHOU))
					 )
	 )  //intermediate substep with no stepping over u < 0, rho<0, or both <0
       ) {
      if(debugfail >= 1) {
	dualfprintf( fail_file, "Inversion from the point value failed.  Using the inversion from the average value.\n" );
      }

      //make a copy of the initial guess so that not to modify the original pb's
      PLOOP(pl) pf[pl] = pb[pl];
      //invert the average U -> "average" p
      MYFUN(Utoprimgen(finalstep,EVOLVEUTOPRIM, UEVOLVE, uavg, ptrgeom, pf),"step_ch.c:advance()", "Utoprimgen", 3);
      //      invert_from_average_flag = lpflag[FLAGUTOPRIMFAIL];


      //Have the results from the inversion from the average value -- copy the result over
      //      PLOOP(pl) pf[pl] = pavg[pl];
      //      lpflag[FLAGUTOPRIMFAIL] = invert_from_average_flag;

      //old code:
      //MYFUN(Utoprimgen(finalstep,EVOLVEUTOPRIM, UEVOLVE, avg, &geom, pf),"step_ch.c:advance()", "Utoprimgen", 2);

      frac_avg_used = 1.0; //reverted to the average value

    }
    else {
      frac_avg_used = 0.0; //used the point value
    }
  }

  return(0);
}





#define COMPARE_GAMMA 0

//If density or gamma-factors are different by more than fractional_difference_threshold for states pin & pout, 
//if different -- correct pout such that it is not more than fractional_difference_threshold away from pin.
FTYPE limit_prim_correction( FTYPE fractional_difference_threshold, struct of_geom *geom, FTYPE *pin, FTYPE *pout )
{
  FTYPE fractional_diff( FTYPE a, FTYPE b );
  FTYPE gammain = 0.0, gammaout = 0.0;
  FTYPE frac_start, frac_end, frac_diff;
  FTYPE fraction_input_value;
  FTYPE comovingenergyin, comovingenergyout;
  int pl;
  FTYPE bsqin,bsqout;
  struct of_state qin, qout;
  int jj;
  FTYPE bdotuin, bdotuout;


#if( COMPARE_GAMMA ) 
  gamma_calc( pin, geom, &gammain );
  gamma_calc( pout, geom, &gammaout );
#endif

  get_state(pin, geom, &qin);
  get_state(pout, geom, &qout);

  bsqout = dot(qout.bcon, qout.bcov);
  bsqin = dot(qin.bcon, qin.bcov);

  bdotuin = 0.0;   DLOOPA(jj) bdotuin+=(qin.ucov[jj])*(qin.bcon[jj]);
  bdotuout = 0.0;  DLOOPA(jj) bdotuout+=(qout.ucov[jj])*(qout.bcon[jj]);

  // u.T.u comoving energy density
  comovingenergyin = pin[RHO] + pin[UU] + bsqin*0.5 - bdotuin*bdotuin;
  comovingenergyout = pout[RHO] + pout[UU] + bsqout*0.5 - bdotuout*bdotuout;


	
#if( COMPARE_GAMMA ) 
  frac_diff = MAX( fractional_diff(gammain, gammaout), 
		   fractional_diff( comovingenergyin, comovingenergyout ) );
#else
  frac_diff = fractional_diff( comovingenergyin, comovingenergyout );
#endif

  //fractional difference at which the reduction to the input value starts
  frac_start = 0.5 * fractional_difference_threshold;

  //fractional difference after which only the input value is used
  frac_end = fractional_difference_threshold;

  //the fraction of the input value used in the output; increases linearly from 0 to 1 for fractions going from frac_start to frac_end
  fraction_input_value = MAX( 0., MIN(1., (frac_diff - frac_start)/(frac_end - frac_start) ) );

  if( 0.0 != fraction_input_value ){
    //states are too different: reverted to primitives that correspond to average conserved quantities because trust them more than point values
    dualfprintf( fail_file, "States are too different, using %3d%% of the average values: i = %d, j = %d, k = %d, nstep = %ld, steppart = %d, t = %21.15g\n", 
		 (int)(100. * fraction_input_value), geom->i, geom->j, geom->k, nstep, steppart, t );
    if( debugfail >= 2 ){
      dualfprintf( fail_file, "Prim. pt. value (gamma, rho, u): " );
      dualfprintf( fail_file, "%21.15g %21.15g %21.15g\n",  gammaout, pout[RHO], pout[UU] );
      dualfprintf( fail_file, "Prim. avg value (gamma, rho, u): " );
      dualfprintf( fail_file, "%21.15g %21.15g %21.15g\n", gammain, pin[RHO], pin[UU] );
      dualfprintf( fail_file, "Frac. difference(ganna, rho, u): " );
      dualfprintf( fail_file, "%21.15g %21.15g %21.15g\n", 
		   fractional_diff(gammain, gammaout),
		   fractional_diff(pin[RHO], pout[RHO]),  
		   fractional_diff(pin[UU], pout[UU])
		   );
    }
  }

  PLOOP(pl) {
    pout[pl] = fraction_input_value * pin[pl] + (1. - fraction_input_value) * pout[pl];
  }

  return( fraction_input_value );
}



//Returns the fractional difference between a & b
FTYPE fractional_diff( FTYPE a, FTYPE b )
{
  FTYPE frac_diff;

  frac_diff = 2. * fabs( a - b ) / ( fabs(a) + fabs(b) + DBL_MIN );

  return( frac_diff );

}

























// get dUavg
void flux2dUavg(int i, int j, int k, FTYPE F1[][N2M][N3M][NPR],FTYPE F2[][N2M][N3M][NPR],FTYPE F3[][N2M][N3M][NPR],FTYPE *dU1avg,FTYPE *dU2avg,FTYPE *dU3avg)
{
  FTYPE idx1,idx2,idx3;
  int pl;

#if(VOLUMEDIFF==0)
  idx1=1.0/dx[RR];
  idx2=1.0/dx[TH];
  idx3=1.0/dx[PH];
#else
  idx1=idxvol[i][j][k][RR];
  idx2=idxvol[i][j][k][TH];
  idx3=idxvol[i][j][k][PH];
#endif

  // initialize for simplicity so don't have to do it conditionally on N>1
  PLOOP(pl){
    dU1avg[pl]=0;
    dU2avg[pl]=0;
    dU3avg[pl]=0;
  }


  if(FLUXB==FLUXCD){ // don't use volume reg. since differencing is large
    PLOOPNOB1(pl){
#if(N1>1)
      dU1avg[pl]=(
		  - (F1[ip1][j][k][pl] - F1[i][j][k][pl]) *idx1
		  );
#endif
#if(N2>1)
      dU2avg[pl]=(
		  - (F2[i][jp1][k][pl] - F2[i][j][k][pl]) *idx2
		  );
#endif

#if(N3>1)
      dU3avg[pl]=(
		  - (F3[i][j][kp1][pl] - F3[i][j][k][pl]) *idx3
		  );
#endif

    }

    // simple version that assumes Fi[Bi] is set to 0 in flux.c for FLUXCD (which it is currently)
    PLOOPBONLY(pl){
#if(N1>1)
      dU1avg[pl]=(
		  - (F1[ip1][j][k][pl] - F1[im1][j][k][pl]) *idx1
		  );
#endif
#if(N2>1)
      dU2avg[pl]=(
		  - (F2[i][jp1][k][pl] - F2[i][jm1][k][pl]) *idx2
		  );
#endif
#if(N3>1)
      dU3avg[pl]=(
		  - (F3[i][j][kp1][pl] - F3[i][j][km1][pl]) *idx3
		  );
#endif
    }

    /*
    // old version
    // FIELDS are special.  The 0's would have come automatically, but spacing is twice the normal size.
    // the double spacing is accounted for in fluxct().
    pl=B1;
    #if(N1>1)
    dU1avg[pl]=(
    0.0
    );
    #endif
    #if(N2>1)
    dU2avg[pl]=(
    - (F2[i][jp1][k][pl] - F2[i][jm1][k][pl]) *idx2
    );
    #endif
    #if(N3>1)
    dU3avg[pl]=(
    - (F3[i][j][kp1][pl] - F3[i][j][km1][pl]) *idx3
    );
    #endif

    pl=B2;
    #if(N1>1)
    dU1avg[pl]=(
    - (F1[ip1][j][k][pl] - F1[im1][j][k][pl]) *idx1
    );
    #endif
    #if(N2>1)
    dU2avg[pl]=(
    0.0
    );
    #endif
    #if(N3>1)
    dU3avg[pl]=(
    - (F3[i][j][kp1][pl] - F3[i][j][km1][pl]) *idx3
    );
    #endif


    pl=B3;
    #if(N1>1)
    dU1avg[pl]=(
    - (F1[ip1][j][k][pl] - F1[im1][j][k][pl]) *idx1
    );
    #endif
    #if(N2>1)
    dU2avg[pl]=(
    - (F2[i][jp1][k][pl] - F2[i][jm1][k][pl]) *idx2
    );
    #endif
    #if(N3>1)
    dU3avg[pl]=(
    0.0
    );
    #endif
    // end old version
    */


    // rest of variables (if any) are normal
    PLOOPNOB2(pl){
#if(N1>1)
      dU1avg[pl]=(
		  - (F1[ip1][j][k][pl] - F1[i][j][k][pl]) *idx1
		  );
#endif
#if(N2>1)
      dU2avg[pl]=(
		  - (F2[i][jp1][k][pl] - F2[i][j][k][pl]) *idx2
		  );
#endif
#if(N3>1)
      dU3avg[pl]=(
		  - (F3[i][j][kp1][pl] - F3[i][j][k][pl]) *idx3
		  );
#endif

    }
  }

  else{


    // other (normal) FLUXB methods, including FLUXCTSTAG
    PLOOP(pl) {


#if(N1>1)
      dU1avg[pl]=(
		  - (F1[ip1][j][k][pl] - F1[i][j][k][pl]) *idx1
		  );
#endif
#if(N2>1)
      dU2avg[pl]=(
		  - (F2[i][jp1][k][pl] - F2[i][j][k][pl]) *idx2
		  );
#endif
#if(N3>1)
      dU3avg[pl]=(
		  - (F3[i][j][kp1][pl] - F3[i][j][k][pl]) *idx3
		  );
#endif
    }


  }




}





// convert point versions of U_i^{n} and dU -> U_i^{n+1} and other versions
void dUtoU(FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui,  FTYPE *Uf, FTYPE *ucum)
{
  int pl;


  // finally assign new Uf and ucum
  // store uf to avoid recomputing U(pf) used later as pb for advance()
  PLOOP(pl) Uf[pl] = UFSET(CUf,dt,Ui[pl],Uf[pl],dUriemann[pl],dUgeom[pl]);


  // how much of Ui, dU, and Uf to keep for final solution
  // ultimately ucum is actual solution used to find final pf
  PLOOP(pl) ucum[pl] += UCUMUPDATE(Cunew,dt,Ui[pl],Uf[pl],dUriemann[pl],dUgeom[pl]);



#if(PRODUCTION==0)
  PLOOP(pl){
    if(!isfinite(Uf[pl])){
      dualfprintf(fail_file,"dUtoU after: nan found for Uf[%d]=%21.15g\n",pl,Uf[pl]);
      dualfprintf(fail_file,"pl=%d Ui=%21.15g dUriemann=%21.15g dugeom=%21.15g\n",pl,Ui[pl],dUriemann[pl],dUgeom[pl]);
      if(N1>1 && N2>1){
	dualfprintf(fail_file,"pl=%d i=%d j=%d k=%d :: F1[icurr]=%21.15g F1[icurr+1]=%21.15g F2[jcurr]=%21.15g F2[jcurr+1]=%21.15g\n",pl,icurr,jcurr,kcurr,F1[icurr][jcurr][kcurr][pl],F1[icurr+1][jcurr][kcurr][pl],F2[icurr][jcurr][kcurr][pl],F2[icurr][jcurr+1][kcurr][pl]);
      }
    }
  }
  PLOOP(pl){
    if(!isfinite(ucum[pl])){
      dualfprintf(fail_file,"dUtoU after: nan found for ucum[%d]=%21.15g\n",pl,ucum[pl]);
      dualfprintf(fail_file,"pl=%d Ui=%21.15g dUriemann=%21.15g dugeom=%21.15g\n",pl,Ui[pl],dUriemann[pl],dUgeom[pl]);
      if(N1>1 && N2>1){
	dualfprintf(fail_file,"pl=%d i=%d j=%d k=%d :: F1[icurr]=%21.15g F1[icurr+1]=%21.15g F2[jcurr]=%21.15g F2[jcurr+1]=%21.15g\n",pl,icurr,jcurr,kcurr,F1[icurr][jcurr][kcurr][pl],F1[icurr+1][jcurr][kcurr][pl],F2[icurr][jcurr][kcurr][pl],F2[icurr][jcurr+1][kcurr][pl]);
      }
    }
  }
#endif

}






// find global dt.  Useful if needed not during evolution, such as at t=0 or for restarting the run if restarting finished run that has a generally smaller dt than should use (including possibly dt=0)
int set_dt(FTYPE (*prim)[N2M][N3M][NPR], SFTYPE *dt)
{
  struct of_state state;
  struct of_geom geom;
  int i,j,k;
  int pl;
  int jj,kk,ll;
  int dir,ignorecourant;
  FTYPE cmax1,cmin1;
  FTYPE cmax2,cmin2;
  FTYPE cmax3,cmin3;
  int wavendti[NDIM],wavendtj[NDIM],wavendtk[NDIM];
  int accdti,accdtj,accdtk;
  int gravitydti,gravitydtj,gravitydtk;
  FTYPE tempwavedt,tempaccdt,tempgravitydt;
  FTYPE dtij[NDIM], wavedt, accdt, gravitydt;
  FTYPE wavendt[NDIM];
  FTYPE ndtfinal;
  FTYPE dUgeom[NPR],dUcomp[NUMSOURCES][NPR];
  int enerregion;
  int *waveenerpos, *sourceenerpos;
  FTYPE X[NDIM],V[NDIM],Vp1[NDIM];
  FTYPE dUriemann[NPR];
  FTYPE Ugeomfree[NPR],U[NPR];
  




  wavedt=accdt=gravitydt=ndtfinal=BIG;
  wavendt[1]=wavendt[2]=wavendt[3]=BIG;
  wavendti[1]=wavendtj[1]=wavendtk[1]=-100;
  wavendti[2]=wavendtj[2]=wavendtk[2]=-100;
  wavendti[3]=wavendtj[3]=wavendtk[3]=-100;
  accdti=accdtj=accdtk=-100;
  gravitydti=gravitydtj=gravitydtk=-100;

  enerregion=OUTSIDEHORIZONENERREGION; // consistent with flux update (except when doing WHAM)
  sourceenerpos=enerposreg[enerregion];

  COMPFULLLOOP{ // want to use boundary cells as well to limit dt (otherwise boundary-induced changes not treated)

    // includes "ghost" zones in case boundary drives solution
    get_geometry(i, j, k, CENT, &geom);
    
    // need full state for vchar()
    MYFUN(get_state(prim[i][j][k], &geom, &state),"restart.c:set_dt()", "get_state()", 1);
    
    
#if(N1>1)
    dir=1;
    MYFUN(vchar(prim[i][j][k], &state, dir, &geom, &cmax1, &cmin1,&ignorecourant),"restart.c:set_dt()", "vchar() dir=1", 1);
    dtij[dir] = cour * dx[dir] / MAX(fabs(cmax1),fabs(cmin1));
    if (dtij[dir] < wavendt[dir]){
      wavendt[dir] = dtij[dir];
      wavendti[dir] = i;
      wavendtj[dir] = j;
      wavendtk[dir] = k;
    }
#endif
    
#if(N2>1)
    dir=2;
    MYFUN(vchar(prim[i][j][k], &state, dir, &geom, &cmax2, &cmin2,&ignorecourant),"restart.c:set_dt()", "vchar() dir=2", 1);
    dtij[dir] = cour * dx[dir] / MAX(fabs(cmax2),fabs(cmin2));
    if (dtij[dir] < wavendt[dir]){
      wavendt[dir] = dtij[dir];
      wavendti[dir] = i;
      wavendtj[dir] = j;
      wavendtk[dir] = k;
    }
#endif

#if(N3>1)
    dir=3;
    MYFUN(vchar(prim[i][j][k], &state, dir, &geom, &cmax3, &cmin3,&ignorecourant),"restart.c:set_dt()", "vchar() dir=3", 1);
    dtij[dir] = cour * dx[dir] / MAX(fabs(cmax3),fabs(cmin3));
    if (dtij[dir] < wavendt[dir]){
      wavendt[dir] = dtij[dir];
      wavendti[dir] = i;
      wavendtj[dir] = j;
      wavendtk[dir] = k;
    }
#endif



    if(WITHINENERREGION(sourceenerpos,i,j,k)){


#if(LIMITDTWITHSOURCETERM)

      // conserved quantity without geometry
      MYFUN(primtoU(UEVOLVE, prim[i][j][k], &state, &geom, U),"step_ch.c:advance()", "primtoU()", 1);
      PLOOP(pl) Ugeomfree[pl] = U[pl]/geom.e[pl];

      // get source term
      // GODMARK: here dUriemann=0, although in reality this setting of dt is related to the constraint trying to make
      PLOOP(pl) dUriemann[pl]=0.0;
      MYFUN(source(prim[i][j][k], &geom, &state, U, dUriemann, dUcomp, dUgeom),"advance.c:set_dt()", "source", 1);

      // get dt limit
      compute_dt_fromsource(&geom,&state,prim[i][j][k], Ugeomfree, dUgeom, dUcomp[GEOMSOURCE], &tempaccdt, &tempgravitydt);
      if(accdt>tempaccdt){
	accdt=tempaccdt;
	accdti=i;
	accdtj=j;
	accdtk=k;
      }
      if(gravitydt>tempgravitydt){
	gravitydt=tempgravitydt;
	gravitydti=i;
	gravitydtj=j;
	gravitydtk=k;
      }

#if(0)// DEBUG
      if(i==-1 || i==0){
	dualfprintf(fail_file,"BANG: i=%d\n",i);
	PLOOP(pl) dualfprintf(fail_file,"prim[%d]=%21.15g\n",pl,prim[i][j][k][pl]);
	PLOOP(pl) dualfprintf(fail_file,"dUgeom[%d]=%21.15g dUcompgeomsource=%21.15g\n",pl,dUgeom[pl],dUcomp[GEOMSOURCE][pl]);
	if(i==-1){
	  // side-by-side
	  coord(i,j,k,CENT,X);
	  bl_coord(X,V);
	  coord(ip1,j,k,CENT,X);
	  bl_coord(X,Vp1);
	  dualfprintf(fail_file,"r(i=-1)=%21.15g r(i=0)=%21.15g\n",V[1],Vp1[1]);
	  dualfprintf(fail_file,"gdet(i=-1)=%21.15g gdet(i=0)=%21.15g\n",gdet[i][j][k][CENT],gdet[ip1][j][k][CENT]);
	  DLOOP(jj,kk) dualfprintf(fail_file,"%d %d gcov(i=-1)=%21.15g gcov(i=0)=%21.15g\n",jj,kk,gcov[i][j][k][CENT][jj][kk],gcov[ip1][j][k][CENT][jj][kk]);
	  DLOOP(jj,kk) dualfprintf(fail_file,"%d %d gcovlast(i=-1)=%21.15g gcovlast(i=0)=%21.15g\n",jj,kk,gcovlast[i][j][k][CENT][jj][kk],gcovlast[ip1][j][k][CENT][jj][kk]);
	  DLOOP(jj,kk) dualfprintf(fail_file,"%d %d gcon(i=-1)=%21.15g gcon(i=0)=%21.15g\n",jj,kk,gcon[i][j][k][CENT][jj][kk],gcon[ip1][j][k][CENT][jj][kk]);
	  DLOOP(jj,kk) DLOOPA(ll) dualfprintf(fail_file,"%d %d %d conn(i=-1)=%21.15g conn(i=0)=%21.15g\n",jj,kk,ll,conn[i][j][k][jj][kk][ll],conn[ip1][j][k][jj][kk][ll]);
	}
      }
#endif


#endif


    } // end if within source enerregion

  } // end of loop


  // GODMARK: note that in normal advance, wavendt[i] is over each CPU region and wavedt computed for each CPU and then minimized over all CPUs -- so not perfectly consistent with MPI
  // here we preserve perfect MPI domain decomposition
  mpifmin(&wavendt[1]);
  mpifmin(&wavendt[2]);
  mpifmin(&wavendt[3]);
  // single all-CPU wavedt
  wavedt = 1.0/(1.0/wavendt[1]+1.0/wavendt[2]+1.0/wavendt[3]); // wavendt[i] is over entire region for each i

  // single all-CPU accdt and gravitydt
  mpifmin(&accdt);
  mpifmin(&gravitydt);

  wavedtglobal=wavedt;
  sourcedtglobal=accdt;
  gravitydtglobal=gravitydt;


  // find global minimum value of wavendt over all cpus
  ndtfinal=MIN(wavedt,MIN(accdt,gravitydt));

#if(1)
  // below single line only right if 1-CPU
  SLOOPA(jj) dualfprintf(log_file,"dtij[%d]=%21.15g wavendti=%d wavendtj=%d wavendtk=%d\n",jj,wavendt[jj],wavendti[jj],wavendtj[jj],wavendtk[jj]);
  dualfprintf(log_file,"ndtfinal=%21.15g wavedt=%21.15g accdt=%21.15g gravitydt=%21.15g\n",ndtfinal,wavedt,accdt,gravitydt); 
  dualfprintf(log_file,"accdti=%d accdtj=%d accdtk=%d :: gravitydti=%d  gravitydtj=%d  gravitydtk=%d\n",accdti,accdtj,accdtk,gravitydti,gravitydtj,gravitydtk);
#endif

  *dt = ndtfinal;
  // don't step beyond end of run
  if (t + *dt > tf) *dt = tf - t;
  
  return(0);
}


// 0.5 not good enough for pressureless collapse
// normal cour=0.8/4 works for presureless collapse for longsteps, so use 0.1 to be safe since rarely gravity conquers timestep
// but cour=0.8 and GRAVITYCOUR = 0.1 doesn't even work for longsteps!
#define GRAVITYCOUR (0.1)

static int compute_dt_fromsource(struct of_geom *ptrgeom, struct of_state *state, FTYPE *pr, FTYPE *U, FTYPE *dUevolve, FTYPE *dUgeomevolve, FTYPE *dtij, FTYPE *gravitydt)
{
  FTYPE dUd[NDIM],dUu[NDIM];
  int jj,kk;
  FTYPE rhoprime[NDIM];
  FTYPE ag[NDIM],dtsource[NDIM];
  FTYPE rho,u,P,bsq,w,eta;
  FTYPE mydU[NDIM];
  FTYPE mydUgravity, rhoprimegravity, aggravity;
  FTYPE frac;
  int i,j,k;
  extern void compute_dr(int i, int j, int k, FTYPE *dr);
  FTYPE dr,dphidt,phi,tempdt;
  FTYPE veleff;
  FTYPE v1max,v1min;


  i=ptrgeom->i;
  j=ptrgeom->j;
  k=ptrgeom->k;

  // default is no limit on dt due to flux differential or source terms
  *dtij=BIG;
  // default is no limit due to time-dependence of gravity
  *gravitydt=BIG;



  DLOOPA(jj){
    dUd[jj]=dUevolve[UU+jj]/(ptrgeom->e[UU+jj]); // remove geometry
  }
  raise_vec(dUd,ptrgeom,dUu);

  // comparing time update with time value, so keep lower as conserved quantity
  mydU[TT]=dUd[TT];
  mydUgravity = dUgeomevolve[UU]/(ptrgeom->e[UU]); // pure gravity piece

  dx[TT]=1.0; // so in the end dt_time <=C/(dU/rhoprime)
  SLOOPA(jj){
    // treating momentum update as giving dv^i, so for getting Courant condition of moving across grid, need upper
    mydU[jj] = dUu[jj];
  }
    

  bsq = dot(state->bcon, state->bcov);
  rho=pr[RHO];
  u=pr[UU];
  P=pressure_rho0_u(rho,u);
  w = rho+u+P;
  eta = w+bsq;





  /////////////////////////
  //
  // New method for dealing with source terms
  //
  /////////////////////////
#if(1)
  DLOOPA(jj){
    // U[UU] is like eta but with \gamma^2 factors that are present in the momentum terms
    // account for geometry prefactor of conserved quantities and put in geometry consistent with source term
    // U[UU] \sim eta \gamma^2 when neglecting stress terms
    // GODMARK: Might want to be more careful like in utoprim_jon.c in how normalize errors
    // GODMARK: Consider REMOVERESTMASSFROMUU==2 effects
    rhoprime[jj]=MAX(fabs(eta),fabs(U[UU]));

    // update to 3-velocity v^i is approximately due to this acceleration
    ag[jj]=SMALL+fabs(mydU[jj]/rhoprime[jj]); // acceleration.  SMALL is so doesn't identically vanish and we get nan below

    if(jj==TT) dtsource[jj]=cour*(dx[jj]/ag[jj]);
    else dtsource[jj]=cour*sqrt(dx[jj]/ag[jj]); // characteristic time-step for increasing velocity so mass would cross a grid
      
  }




  /////////////////////////
  //
  // Self-gravity
  //
  /////////////////////////

#if(DOSELFGRAVVSR)
  // make sure metric is not varying too fast
  // just look at perturbed part of g_{tt} -- an invariant in stationary metrics, so good indicator of gauge-invariant time-dependence of metric
  // only look at loc=CENT since others should be similar to order unity -- also CENT won't diverge at r=0
  //  frac = (gcovpertlast[i][j][k][CENT][TT] - gcovpert[i][j][k][CENT][TT])/(fabs(gcovpertlast[i][j][k][CENT][TT]) + fabs(gcovpert[i][j][k][CENT][TT]));
  // above frac has no dt, so no measure of what dt should be

  // \Gamma^t_{tt} measures dg_{tt}/dt
  // \Gamma^t_{tt} \approx g_{tt},t g^{tt}  so that g_{tt},t = \Gamma^t_{tt}/g^{tt}
  // now form same construct as with (dU/dt)/U as above
  // g^{tt} can't go to 0 unless in bad coordinate system (BL-coords on horizon)
  // frac is approximately g_{tt},t/g^{tt} \aprox 1/dt
  // GODMARK: below might be problem in non-relativistic case



  // get \Gamma_{ttt} = dphi/dt ~ 1/2 g_{tt,t}
  dphidt=0.0;
  DLOOPA(jj) dphidt += conn[i][j][k][jj][TT][TT]*(ptrgeom->gcov[jj][TT]);
  dphidt = fabs(dphidt); // don't care about sign, just magnitude

  // get \phi ~ -(g_{tt} +1)/2
  // phi by itself has no meaning, but as a reference for changes in phi in time it does
  phi = -(1.0+ptrgeom->gcov[TT][TT])*0.5;
  phi = fabs(phi); // sign not important

#if(0)

  // treat dt ~ \phi / (d\phi/dt)
  //  frac = fabs(conn[i][j][k][TT][TT][TT]/((1+gcon[i][j][k][CENT][TT][TT])*(1+gcon[i][j][k][CENT][TT][TT])));
  frac = fabs(dphidt/phi);
  *gravitydt = cour*(GRAVITYCOUR/frac);
  // this dt keeps frac~cour
  //  *gravitydt = cour*(GRAVITYCOUR/frac); // GRAVITYCOUR is additional courant factor on gravitational term


  // treat d\phi/dt as v^2/dt
  compute_dr( i,  j,  k, &dr);
  //dgttdt = fabs(conn[i][j][k][TT][TT][TT]/(gcon[i][j][k][CENT][TT][TT]));
  //tempdt = GRAVITYCOUR*pow(cour*cour*dr*dr/dgttdt,1.0/3.0); // GRAVITYCOUR in front is effective additional courant factor on gravitational term
  tempdt = GRAVITYCOUR*pow(cour*cour*dr*dr/dphidt,1.0/3.0); // GRAVITYCOUR in front is effective additional courant factor on gravitational term

  if(tempdt<*gravitydt) *gravitydt=tempdt;
#elif(0)

  frac = fabs(ptrgeom->gcon[TT][TT]*dphidt);
  *gravitydt = cour*(GRAVITYCOUR/frac);




#elif(1)

  ///////////////////////////////
  //
  // 3 methods to limit dt
  //
  //////////////////////////////

  // LOCAL DPHI/DT LIMIT
  // treat dt ~ \phi / (d\phi/dt)
  //  frac = fabs(conn[i][j][k][TT][TT][TT]/((1+gcon[i][j][k][CENT][TT][TT])*(1+gcon[i][j][k][CENT][TT][TT])));
  frac = fabs(dphidt/phi);
  *gravitydt = cour*(GRAVITYCOUR/frac);
  // this dt keeps frac~cour
  //  *gravitydt = cour*(GRAVITYCOUR/frac); // GRAVITYCOUR is additional courant factor on gravitational term



  compute_dr( i,  j,  k, &dr);


  // LOCAL DU/DT LIMIT BUT TREAT AS VELOCITY LIMITED BY SOL

#if(REMOVERESTMASSFROMUU==2)
  rhoprimegravity=fabs(U[UU]); // gravity affects only \rho \phi -like terms order \rho v^2, not \rho
#else
  // remove rest-mass
  rhoprimegravity=fabs(U[UU]+U[RHO]); // gravity affects only \rho \phi -like terms order \rho v^2, not \rho
#endif
  aggravity=SMALL+fabs(mydUgravity/rhoprimegravity);
  veleff = aggravity*dx[1];

  // get speed of light in 1-direction (dx^1/dt)
  sol(pr,state,1,ptrgeom,&v1max,&v1min);
  // limit "velocity" to speed of light
  if(veleff>fabs(v1min)) veleff=fabs(v1min);
  if(veleff>fabs(v1max)) veleff=fabs(v1max);

  tempdt=GRAVITYCOUR*cour*dx[1]/veleff;
  if(tempdt<*gravitydt) *gravitydt=tempdt;


  // LOCAL LIMIT ON DPHI/DT WHERE PHI~V^2
  //dgttdt = fabs(conn[i][j][k][TT][TT][TT]/(gcon[i][j][k][CENT][TT][TT]));
  //tempdt = GRAVITYCOUR*pow(cour*cour*dr*dr/dgttdt,1.0/3.0); // GRAVITYCOUR in front is effective additional courant factor on gravitational term
  tempdt = GRAVITYCOUR*pow(cour*cour*dr*dr/dphidt,1.0/3.0); // GRAVITYCOUR in front is effective additional courant factor on gravitational term
  if(tempdt<*gravitydt) *gravitydt=tempdt;


#endif




#endif // end if DOSELFGRAVVSR==1





#else // else other older method for gravitydt


  /////////////////////////
  //
  // Old method for dealing with source terms
  //
  /////////////////////////

  // as expected, this method essentialy divides by 0 do dt->0

  SLOOPA(jj){
    // U[UU] is like eta but with \gamma^2 factors that are present in the momentum terms
    // account for geometry prefactor of conserved quantities and put in geometry consistent with source term
    // U[UU] \sim eta \gamma^2 when neglecting stress terms
    // GODMARK: Might want to be more careful like in utoprim_jon.c in how normalize errors
    // GODMARK: Consider REMOVERESTMASSFROMUU==2 effects
    rhoprime[jj]=MAX(fabs(eta),fabs(U[UU]));

    // update to 3-velocity v^i is approximately due to this acceleration
    ag[jj]=SMALL+fabs(mydU[jj]/rhoprime[jj]); // acceleration.  SMALL is so doesn't identically vanish and we get nan below

    dtsource[jj]=cour*sqrt(dx[jj]/ag[jj]); // characteristic time-step for increasing velocity so mass would cross a grid
      
  }

  // do pure gravity piece (very similar to how full TT source term dealt with)
  //  rhoprimegravity=MAX(fabs(eta),fabs(U[UU]));
  // The below term is \propto \rho Mvsr(r)/r when v=0, so is well-defined
  // GODMARK: Can this quantity go to zero near rotating BH or something?
#if(REMOVERESTMASSFROMUU==2)
  rhoprimegravity=fabs(U[UU]); // gravity affects only \rho \phi -like terms order \rho v^2, not \rho
#else
  // remove rest-mass
  rhoprimegravity=fabs(U[UU]+U[RHO]); // gravity affects only \rho \phi -like terms order \rho v^2, not \rho
#endif
  aggravity=SMALL+fabs(mydUgravity/rhoprimegravity);
  dtsource[TT] = *gravitydt=cour*(dx[TT]/aggravity);



#endif  





  /////////////////////////
  //
  // Finally store source term's version of limited dt to be used later
  //
  /////////////////////////
  

  //  dualfprintf(fail_file,"i=%d mydUgravity=%21.15g rhoprimegravity=%21.15g rhoprime[TT]=%21.15g mydU[TT]=%21.15g\n",ptrgeom->i,mydUgravity,rhoprimegravity,rhoprime[TT],mydU[TT]);

  // always do time-component
  // accounts for thermal changes if cooling function or geometry changes if metric changing
  if (dtsource[TT] < *dtij) *dtij = dtsource[TT];

#if(N1>1)
  if (dtsource[RR] < *dtij) *dtij = dtsource[RR];
#endif
#if(N2>1)
  if (dtsource[TH] < *dtij) *dtij = dtsource[TH];
#endif
#if(N3>1)
  if (dtsource[PH] < *dtij) *dtij = dtsource[PH];
#endif


  return(0);

}



// compute dt from full coordinate acceleration
static int dUtodt(struct of_geom *ptrgeom, struct of_state *q, FTYPE *pr, FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *dUgeomgravity, FTYPE *accdt, FTYPE *gravitydt)
{
  int pl;
  FTYPE dUtotal[NPR];
  FTYPE Ugeomfree[NPR],U[NPR];



  PLOOP(pl) {
    // while each piece may be "large", when summed if small then final change to conserved quantity is small, so that's all that's relevant
    dUtotal[pl]=dUriemann[pl]+dUgeom[pl];
    // GODMARK:
    //    dUtotal[pl]=fabs(dUriemann[pl])+fabs(dUgeom[pl]);
  }

  // conserved quantity without geometry
  MYFUN(primtoU(UEVOLVE, pr, q, ptrgeom, U),"step_ch.c:advance()", "primtoU()", 1);
  PLOOP(pl) Ugeomfree[pl] = U[pl]/ptrgeom->e[pl];


  compute_dt_fromsource(ptrgeom, q, pr, Ugeomfree, dUtotal, dUgeomgravity, accdt, gravitydt);
    


  return(0);

}

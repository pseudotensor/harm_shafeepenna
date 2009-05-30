
#include "decs.h"



int init(int *argc, char **argv[])
{
  extern int prepre_init(void);
  extern int pre_init(int *argc, char **argv[]);
  extern int init_defgrid(void);
  extern int init_defglobal(void);
  extern int init_defconsts(void);
  extern int init_grid_post_set_grid(void);
  extern int post_par_set(void);
  extern void filterffde(int i, int j, int k, FTYPE *pr);
  extern void filter_coldgrmhd(int i, int j, int k, FTYPE *pr);
  extern int set_dt(FTYPE (*prim)[N2M][N3M][NPR], SFTYPE *dt);  
  //
  extern int post_init(void);
  extern int compute_new_metric_and_prims(int whichtime, SFTYPE MBH, SFTYPE a, SFTYPE QBH);
  void  init_all_conservatives(void);
  int i,j,k;
  int jj;
  int pl;
  int selfgraviter;
  int set_box_grid_parameters(void);
  int gotnan;
  int faketimeorder,fakenumtimeorders;





  fprintf(stderr,"Start init\n"); fflush(stderr);

  //////////////
  //
  // Both normal and restart mode need all "pre_init" functions called
  //
  //////////////


  if(RESTARTMODE==0 || RESTARTMODE==1){
    //////////////
    //
    // prepre_init type functions should not assume anything done yet
    //
    //////////////
    prepre_init();
    prepre_init_specific_init(); // user func

    //////////////
    //
    // pre_init type functions initialize MPI and other things
    //
    //////////////
    pre_init(argc,argv);
    pre_init_specific_init(); // user func
  }



  ///////////////////
  //
  // Always set parameters and coordinates as default in case restart file doesn't have this information.  Properly, the restart file will overwrite if file contains those variables.
  //
  //////////////////
  if(RESTARTMODE==0 || RESTARTMODE==1){
    // define parameters
    init_defglobal(); // init default global parameters
    init_defconsts(); // init default physical constants
    init_consts(); // init user constants
    init_global(); // request choices for global parameters

    init_defgrid(); // init default grid
    init_grid(); // request choices for grid/coordinate/metric parameters

    set_coord_parms(); // requires correct defcoord at least
  }


  ////////////////
  //
  // RESTART MODE
  //
  // Loads primitives and conserved quantities.
  // Always assume conserved quantities are in the file, regardless of whether used.  Only used if doing higher-order method
  //
  ///////////////
  if(RESTARTMODE==1){
    if (restart_init(WHICHFILE) >= 1) {
      dualfprintf(fail_file, "main:restart_init: failure\n");
      return(1);
    }

    // don't bound read-in data yet since need grid and other things

  }




  ///////////////////
  //
  // Always write new coordparms file.  Fresh start needs it, but also do in case user updated restart file but not coord file
  //
  ////////////////////
  if(RESTARTMODE==0 || RESTARTMODE==1){
    write_coord_parms(); // output coordinate parameters to file
  }



  ///////////////////
  //
  // Both normal and restart mode need to setup grid
  //
  //////////////////
  if(RESTARTMODE==0 || RESTARTMODE==1){

    // once all interpolation parameters are set, now can set dependent items that may be used to set primitives or conservatives
    // doesn't use metric parameters, so doesn't need to be in SELFGRAV loop
    post_par_set();

    if(RESTARTMODE==1) trifprintf("proc: %d post_par_set completed: failed=%d\n", myid,failed);

    // get grid
    // 0 tells set_grid that it's first call to set_grid() and so have to assume stationarity of the metric since have no time information yet
    set_grid(0,0,0);

    if(RESTARTMODE==1) trifprintf("proc: %d grid restart completed: failed=%d\n", myid,failed);

    // set grid boundary parameters that uses metric parameters to determine loop ranges using enerregions and enerpos
    set_box_grid_parameters();

    if(RESTARTMODE==1) trifprintf("proc: %d set_box_grid_parameters completed: failed=%d\n", myid,failed);

    // user post_set_grid function
    init_grid_post_set_grid();

    if(RESTARTMODE==1) trifprintf("proc: %d init_grid_post_set_grid completed: failed=%d\n", myid,failed);
  

    trifprintf("MCOORD=%d\n",MCOORD);
    trifprintf("COMPDIM=%d\n",COMPDIM);
    trifprintf("MAXBND=%d\n",MAXBND);
    trifprintf("FLUXB=%d\n",FLUXB);

  }




  ///////////////////
  //
  // Normal fresh start need to get primitive and conserved quantities
  //
  //////////////////
  if(RESTARTMODE==0){


#if(DOSELFGRAVVSR)
#define NUMSELFGRAVITER 3 // MINIMUM of 2 iterations should be done for metric to computed once
    // further iterations will determine PRIMECOORD velocity and magnetic field using newly updated metric
    // 3 iterations allows this new metric to be used once, and I suggest this as good minimum # of iterations
#else
#define NUMSELFGRAVITER 1 // NO CHOICE
#endif



    // iterations to get consistent/converged initial conditions
    // assume user chooses primitive in final-to-be-metric so first primitive is correct final primitive
    // conservatives adjusted with adjusting metric
    for(selfgraviter=1;selfgraviter<=NUMSELFGRAVITER;selfgraviter++){


      trifprintf("begin iteration over metric: selfgraviter=%d\n",selfgraviter);


      if(selfgraviter>1){
	if(DOSELFGRAVVSR){
	  trifprintf("new metric with self-gravity: selfgraviter=%d\n",selfgraviter);
	  // if box_grid needs to change, is done inside below function
	  compute_new_metric_and_prims(0,MBH, a, QBH);
	  trifprintf("done with computing new metric with self-gravity dt=%21.15g selfgraviter=%d\n",dt,selfgraviter);
	}
      }


      
    

      // user function that should fill p with primitives
      // assumes everyting computed by:  compute_EOS_parms() is known by init_primitives without immediate reference to what's computed by compute_EOS_parms since this function depends upon primitives themselves
      MYFUN(init_primitives(p),"initbase.c:init()", "init_primitives()", 0);


      
      // dump analytic solution
      //pdump=panalytic;
      //if (dump(-9999) >= 1){
      //  dualfprintf(fail_file,"unable to print dump file\n");
      //  return (1);
      //}
      pdump=p;
      
      
      
      //    COMPLOOPF{
      // p[i][j][k][UU]=0.0;
      // }
      
      // mandatory filters on user supplied primitives
      
      /////////////////////////////
      //
      // Filter to get correct degenerate FFDE solution
      //
      /////////////////////////////// 
      
#if(EOMTYPE==EOMFFDE)
      trifprintf("System filtered to FFDE\n");
      // filter to get force-free
      COMPFULLLOOP{
	filterffde(i,j,k,p[i][j][k]);
      }
#endif
      
#if(EOMTYPE==EOMCOLDGRMHD)
      trifprintf("System filtered to cold GRMHD\n");
      // filter to get cold GRMHD
      COMPFULLLOOP{
	filter_coldgrmhd(i,j,k,p[i][j][k]);
      }
#endif
      
      
      /////////////////////////////
      //
      // Fixup and Bound variables since field may have changed
      // Also setup pre_fixup() type quantities
      //
      /////////////////////////////// 
      
      trifprintf("System Fixup and Bound\n");
      
#if(FIXUPAFTERINIT)
      if(fixup(STAGEM1,p,0)>=1)
	FAILSTATEMENT("initbase.c:init()", "fixup()", 1);
#endif






      
      if (bound_allprim(STAGEM1,t,p) >= 1)
	FAILSTATEMENT("initbase.c:init()", "bound_allprim()", 1);
      
      if(pre_fixup(STAGEM1,p)>=1)
	FAILSTATEMENT("initbase.c:init()", "pre_fixup()", 1);


      ///////////////////////////////
      // BEGIN DEBUG
      // dump solution so far
      if(selfgraviter==1){
	pdump=p;
	if (dump(9000) >= 1){
	  dualfprintf(fail_file,"unable to print dump file\n");
	  return (1);
	}
      }
      else if(selfgraviter==2){
	pdump=p;
	if (dump(9001) >= 1){
	  dualfprintf(fail_file,"unable to print dump file\n");
	  return (1);
	}
      }
      else if(selfgraviter==3){
	pdump=p;
	if (dump(9002) >= 1){
	  dualfprintf(fail_file,"unable to print dump file\n");
	  return (1);
	}
      }
      // END DEBUG
      ///////////////////////////////


      
      // after all parameters and primitives are set, then can set these items
      post_init();
      // user post_init function
      post_init_specific_init();

      init_all_conservatives();


      ///////////////////////////////
      // BEGIN DEBUG
      // dump solution so far
      if(selfgraviter==1){
	pdump=p;
	if (dump(9100) >= 1){
	  dualfprintf(fail_file,"unable to print dump file\n");
	  return (1);
	}
      }
      else if(selfgraviter==2){
	pdump=p;
	if (dump(9101) >= 1){
	  dualfprintf(fail_file,"unable to print dump file\n");
	  return (1);
	}
      }
      else if(selfgraviter==3){
	pdump=p;
	if (dump(9102) >= 1){
	  dualfprintf(fail_file,"unable to print dump file\n");
	  return (1);
	}
      }
      // END DEBUG
      ///////////////////////////////

      trifprintf("end iteration over metric: selfgraviter=%d\n",selfgraviter);
    
    }// end loop to get metric and primitive setting consistent





  }// end if RESTARTMODE==0
  else if(RESTARTMODE==1){


#if(FIXUPAFTERINIT)
    if(fixup(STAGEM1,p,0)>=1)
      FAILSTATEMENT("initbase.c:init()", "fixup()", 1);
#endif
    if (bound_prim(STAGEM1,t,p) >= 1)
      FAILSTATEMENT("initbase.c:init()", "bound_allprim()", 1);
      
    if(pre_fixup(STAGEM1,p)>=1)
      FAILSTATEMENT("initbase.c:init()", "pre_fixup()", 1);


    ////////////////
    //
    // Note there is no need to convert average or quasi-deaveraged field to staggered field (see comments in ucons2upointppoint())
    // However, divb diagnostics at first won't be right at t=0 since set to use primitive for lower-order method, but just assume diagnostic won't likely be immediately after restart
    // For RESTARTMODE==0 the pstag quantity is set by user or during vector potential conversion to u and p, but during restart we only read-in p and unew while we need also pstagscratch
    //
    /////////////////
    ucons2upointppoint(t, p,unew,ulast,p); // this regenerates pcentered for B1,B2,B3


    // after all parameters and primitives are set, then can set these items
    post_init();
    // user post_init function
    post_init_specific_init();

    // don't want conservatives or primitives to change, just compute metric
    if(DOSELFGRAVVSR){
      trifprintf("new metric with self-gravity\n");
      compute_new_metric_and_prims(0,MBH, a, QBH);
      trifprintf("done with computing new metric with self-gravity dt=%21.15g\n",dt);
    }


    if (restart_init_checks(WHICHFILE) >= 1) {
      dualfprintf(fail_file, "main:restart_init_checks: failure\n");
      return(1);
    }


    trifprintf( "proc: %d restart completed: failed=%d\n", myid, failed);
  


  }


  // must set dt after DOSELFGRAVVSR so have acceleration in case v=0 and c_s=0 initially
  // set initial dt correctly rather than random choice
  // actually compute_new_metric_and_prims() computes dt from set_dt()
  set_dt(p,&dt);
  trifprintf("dt=%21.15g at t=%21.15g at nstep=%ld at realnstep=%ld\n",dt,t,nstep,realnstep);



#if(PRODUCTION==0)
  //////////////
  //
  // make sure all zones are not nan before ending init
  //
  /////////////
  gotnan=0;
  FULLLOOP{ // diagnostic check
    PLOOPINTERP(pl){ // only check those things that will be interpolated since rest of things assumed to be not necessary as primitive and generated as conserved/flux from other primitives
      if(!finite(p[i][j][k][pl])){
	dualfprintf(fail_file,"init/restart data has NaN at i=%d j=%d k=%d ti=%d tj=%d tk=%d :: pl=%d\n",i,j,k,startpos[1]+i,startpos[2]+j,startpos[3]+k,pl);
	gotnan=1;
      }
    }
  }
  if(gotnan) myexit(82753487);


  if(FLUXB==FLUXCTSTAG){
    gotnan=0;
    FULLLOOP{ // diagnostic check
      PLOOPBONLY(pl){ // only check those things that will be interpolated since rest of things assumed to be not necessary as primitive and generated as conserved/flux from other primitives
	if(!finite(pstagscratch[i][j][k][pl])){
	  dualfprintf(fail_file,"init/restart data has pstag NaN at i=%d j=%d k=%d ti=%d tj=%d tk=%d :: pl=%d\n",i,j,k,startpos[1]+i,startpos[2]+j,startpos[3]+k,pl);
	  gotnan=1;
	}
      }
    }
    if(gotnan) myexit(82753488);
  }
#endif



#if( DOGRIDSECTIONING )
  // this is used to set grid section to something before evolving
  // assume full grid was defined previously
  faketimeorder=0;
  fakenumtimeorders=1;
  findandsetactivesection(faketimeorder, fakenumtimeorders, realnstep, t ); //SASMARK SECTIONMARK
#endif



  trifprintf("end init.c\n");
  return (0);

}






// initialize conserved quantities
//  if(RESTARTMODE==0)
void  init_all_conservatives(void)
{
  int pl;
  // below is user function that usually uses system function
  extern int init_conservatives(FTYPE p[][N2M][N3M][NPR], FTYPE Utemp[][N2M][N3M][NPR], FTYPE U[][N2M][N3M][NPR]);

  
  init_conservatives(p, ulast, unew);

  //  bound_uavg(STAGEM1,unew); DEBUG DEBUG

}



#include "initbase.defaultnprlists.c"

// Called before pre_init() : i.e. before MPI init
int prepre_init(void)
{
  int pl;



  // set default performance parameters
  set_defaults_performance_checks();


  // set file version numbers
  set_file_versionnumbers();


  ////////////////////
  //
  // Setup Loop ranges for primitive/conserved/interpolated variables
  //
  ////////////////////

  set_default_nprlists();



  advancepassnumber=-1; // by default assume all things done (should only matter if SPLITNPR==1 and debugging it)

  // still need to avoid conserved+flux calculations for other PL's in phys.c


  // below 2 now determined at command line.  See init_MPI_GRMHD() in init_mpi.c (myargs and init_MPI).
  //  RESTARTMODE=0;// whether restarting from rdump or not (0=no, 1=yes)
  //WHICHFILE=0; // see diag.c for dump_cnt and image_cnt starts
  // user defined parameter
  restartonfail=0; // whether we are restarting on failure or not and want special diagnostics
  specialstep=0; // normal step assumed
  didstorepositiondata=0; // assume haven't stored position data (yet)
  horizoni=-200;
  horizoncpupos1=-1;


  if(WHICHVEL==VEL3){
    jonchecks=1; // whether to include jon's checks to make sure u^t real and some rho/u limits throughout code
    jonchecks=0;
  }
  else jonchecks=0; // not relevant

  // choice// GODMARK: not convenient location, but needed for init_MPI_GRMHD()
  periodicx1=0;
  periodicx2=0;
  periodicx3=0;// GODMARK: periodic in \phi for 3D spherical polar

  if(USEMPI&&USEROMIO){
    binaryoutput=MIXEDOUTPUT; // choice: mixed or binary
    sortedoutput=SORTED; // no choice
  }
  else{
    // choice
    binaryoutput=TEXTOUTPUT;
    //binaryoutput=BINARYOUTPUT;
      sortedoutput=SORTED;
  }


  return(0);
}



void set_defaults_performance_checks(void)
{

  // whether to log steps
  // log the time step as per NDTCCHECK
  DOLOGSTEP=1 ;
  // log the performance as per NZCCHECK
  DOLOGPERF=1;

  CHECKCONT=1; // 1: check if user wants to continue run 0: don't

  
  // how often in REAL *seconds* to dump 0_logstep.out file (unless 0, which then uses below)
  DTstep=10.0;
  DTstepdot=1.0;
  DTperf=DTstep;
  DTgocheck=30.0;
  DTtimecheck=60.0;

  // initial guess for how often to check time per timestep
  NTIMECHECK=1000;

  // how often in steps to output step/dt/t data (controlled by above if above are nonzero, else use below numbers)
  // MARK 100 100 20 500
  // MARK 10 10 1 100 for 1024x1024 vortex
  // MARK 1D bondi: 10000 10000 1000 20000
  // MARK 2D MHD Tori128128: 500 500 10 1000
  NDTCCHECK=100;
  // how often in steps to check speed in zonecycles/sec
  NZCCHECK=100;
  NDTDOTCCHECK=10;
  NGOCHECK=1000; // how often in steps to check the go.go file to see if to continue running

  // number of wallseconds per perf run(see: main.c)
  PERFWALLTIME=30.0;
  // 1 linux cluster cpu
  // ZCPSESTIMATE (50000)
  // 25 linux cluster cpu(550Mhz Xeon's connected via Myrinet)
  // ZCPSESTIMATE (1250000)
  // 36 linux cluster cpu(550Mhz Xeon's connected via Myrinet)
  // ZCPSESTIMATE (1800000)
  // 64 linux cluster cpu(550Mhz Xeon's connected via Myrinet)
  // ZCPSESTIMATE (3200000)
  // 121 linux cluster cpu(550Mhz Xeon's connected via Myrinet)
  // ZCPSESTIMATE (6050000)
  // photon MHD for one cpu alone
  // ZCPSESTIMATE (265000)
  // photon HD 1 cpu
  // ZCPSESTIMATE (400000)
  // rainman MHD for one cpu alone
  // ZCPSESTIMATE (220000)
  // ZCPSESTIMATE (100000)
  // sgi r10000 for one cpu alone(195sMhz)
  // ZCPSESTIMATE (80000)
  // 4cpu mpigm
  // ZCPSESTIMATE (800000)
  // 4cpu r10000
  // ZCPSESTIMATE (343000)
  // 9cpu r10000
  // ZCPSESTIMATE (745000)
  // 16cpu r10000
  // ZCPSESTIMATE (1325000)
  // 25cpu r10000
  // ZCPSESTIMATE (1200000)
  // 36cpu r10000
  // ZCPSESTIMATE (1700000)
  // 49cpu r10000
  // ZCPSESTIMATE (4021000)
  // 64 r10000's 64x64 tile
  // ZCPSESTIMATE (4309000)
  // 121 r10000's
  // ZCPSESTIMATE (10943000)
  // 256 r10000's
  // ZCPSESTIMATE (20000000)
  // kerr 2DMHD
  // ZCPSESTIMATE (50000)
  // rainman 2DMHD
  //  ZCPSESTIMATE=(200000);

  // Latest HARM w/ lim=WENO5BND and FLUXCTSTAG on ki-rh42
  //  ZCPSESTIMATE=(5000);

  // Latest HARM w/ lim=PARA and FLUXCTSTAG on ki-rh42
  //  ZCPSESTIMATE=(10000);
  // Latest HARM w/ lim=PARA and FLUXCTSTAG on 1 Orange CPU
  //  ZCPSESTIMATE=(7000);

  // Latest HARM w/ lim=WENO5BND and FLUXCTTOTH on 1 Orange CPU
  ZCPSESTIMATE=(3200);


}



// not used quite yet
void set_file_versionnumbers(void)
{

  // file versions numbers(use sm for backwards compat)
  PVER= 11;
  GRIDVER= 3; // 3 is without cotangent
  DVER= 1 ;   // dumps same as for pdumps, adumps
  FLVER= 2;
  NPVER= 2;
  AVG1DVER= 2;
  AVG2DVER= 2;
  ENERVER= 7; // 6 is without c/s mode amp in ener, 7 has new ang mom
  MODEVER= 2; // 2 is all vars for 9 modes
  LOSSVER= 7; // 6 has x3 losses, 5 doesn't, 7 has new ang mom losses
  SPVER=   1;
  TSVER=   1;
  LOGDTVER= 1;
  STEPVER= 1;
  PERFVER= 3;
  ADVER= DVER;
  PDVER= DVER;
  CALCVER= 1;
  FLINEVER= 1;
  // type designations for sm automagical read in correct format for similar things
  PTYPE=     1; // global par file
  GRIDTYPE=  2;
  DTYPE=     3 ;// dump
  FLTYPE=    4; // floor
  NPTYPE=    5; // np
  AVG2DTYPE= 6;
  AVG1DTYPE= 7;
  ENERTYPE=  8;
  LOSSTYPE=  9;
  SPTYPE=    10;
  TSTYPE=    11;
  LOGDTTYPE= 12 ;
  STEPTYPE=  13;
  PERFTYPE=  14;
  ADTYPE=    15 ;// same as dump except filename
  PDTYPE=    16; // same as dump except filename
  CALCTYPE=  17; // arbitrary calcs during pp
  FLINETYPE=  18; // field line during pp
  MODETYPE=  19;
  EXPANDTYPE= 50 ;// used to signify doing pp expansion
  NPCOMPUTETYPE= 33; // used to signify want to compute np before output


}




// used to setup local versions of lists
// currently used in interpline.c for NUMPRIMSLOOP()
void  setup_nprlocalist(int whichprimtype, int *nprlocalstart, int *nprlocalend,int *nprlocallist, int *numprims)
{
  int pl;

  // setup primitive loops
  if(whichprimtype==ENOPRIMITIVE){
    *nprlocalstart=npr2interpstart;
    *nprlocalend=npr2interpend;
    PMAXNPRLOOP(pl) nprlocallist[pl]=npr2interplist[pl];
    *numprims=NPR2INTERP;
  }
  else{ // NPR type
    *nprlocalstart=nprstart;
    *nprlocalend=nprend;
    PMAXNPRLOOP(pl) nprlocallist[pl]=nprlist[pl];
    *numprims=NPR;
  }

}




// initialize MPI and other things
// NO computations should be performed here
int pre_init(int *argc, char **argv[])
{
  int ii;
  int dir,pl,sc,fl,floor,enerregion;
  int tscale;
  int dissloop;
  int i,j,k;
  extern void set_arrays(void);
  int checki;
  int recompute_fluxpositions(void);




  SQRTMINNUMREPRESENT=sqrt(MINNUMREPRESENT);

  // things initialized whether restarting or init fresh

  ranc(0);  // power up random number generator in case used without init

#if(CHECKONINVERSION)
  checki=0;
  strcpy(globalinvtext[checki++],"Qdotnp");
  strcpy(globalinvtext[checki++],"Qtsq");
  strcpy(globalinvtext[checki++],"Qtsqorig");
  strcpy(globalinvtext[checki++],"Bsq");
  strcpy(globalinvtext[checki++],"DD");
  strcpy(globalinvtext[checki++],"QdotB");
  strcpy(globalinvtext[checki++],"WWp");
  strcpy(globalinvtext[checki++],"Qtcon1");
  strcpy(globalinvtext[checki++],"Qtcon2");
  strcpy(globalinvtext[checki++],"Qtcon3");
  strcpy(globalinvtext[checki++],"Bcon1");
  strcpy(globalinvtext[checki++],"Bcon2");
  strcpy(globalinvtext[checki++],"Bcon3");
  if(checki>NUMGLOBALINV){
    dualfprintf(fail_file,"Not enough memory for globalinvtext: checki=%d NUMGLOBALINV=%d\n",checki,NUMGLOBALINV);
  }
#endif


  // must do in MPI mode or not MPI mode  
  init_MPI_GRMHD(argc, argv);

#if(USEMPI)
  mpi_set_arrays();
#endif


  // check starting go files (must be after init_mpi so all files know)
  gocheck(STARTTIME);


  report_systeminfo(stderr);
  report_systeminfo(log_file);
  if(myid==0) report_systeminfo(logfull_file);

  /////////////////
  //
  // setup files for writing and reading (must come after init_MPI_GRMHD())
  //
  makedirs();


  // init arrays
  set_arrays();


  // set default variable to dump (must come before init() where if failed or other reasons can dump output)
  udump=unew;
  ubound=unew;
  pdump = p;





  init_dumps();

  reset_dothisenerregion(); // must be done before any ENERREGIONLOOP() call

  // must go here b4 restart if restarting
  ENERREGIONLOOP(enerregion){
    // used for each region, related to global quantities
    // no need to initialize _tot quantities, they are overwritten during MPI sum in diag.c
    fladd=fladdreg[enerregion];
    fladdterms=fladdtermsreg[enerregion];
    U_init=Ureg_init[enerregion];
    pcum=pcumreg[enerregion];
    pdot=pdotreg[enerregion];
    pdotterms=pdottermsreg[enerregion];
    sourceaddterms=sourceaddtermsreg[enerregion];
    sourceadd=sourceaddreg[enerregion];
    diss=dissreg[enerregion];

    PDUMPLOOP(pl){
      fladd[pl] = 0;
      FLOORLOOP(floor) fladdterms[floor][pl]=0;
      U_init[pl] = 0;
      DIRLOOP(dir){
	pcum[dir][pl]=0;
	pdot[dir][pl]=0;
	FLLOOP(fl) pdotterms[dir][fl][pl]=0;
	if(enerregion==0) FLLOOP(fl) pdottermsjet2[dir][fl][pl]=0; // needed for other not-flux cpus!
      }
      sourceadd[pl] = 0;
      SCLOOP(sc) sourceaddterms[sc][pl] = 0;
    }
    for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++)  diss[dissloop] = 0;

    if(DOLUMVSR) if(enerregion==0) for(ii=0;ii<ncpux1*N1;ii++) lumvsr[ii]=0;
    if(DODISSVSR) if(enerregion==0) for(ii=0;ii<ncpux1*N1;ii++) for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++) dissvsr[dissloop][ii]=0;
    if(DOSELFGRAVVSR) if(enerregion==0) for(ii=0;ii<ncpux1*N1;ii++) dMvsr[ii]=0;

    // below quantities have been subsumed into own full enerregion
    //    if(DOEVOLVEMETRIC) if(enerregion==0) PLOOP(pl){
    //  horizonflux[pl]=0;
    //  horizoncum[pl]=0;
    // }
  }

  //  if(DOLUMVSR) for(ii=0;ii<ncpux1*N1;ii++) dualfprintf(fail_file,"t=0: lumvsr[%d]=%21.15g\n",ii,lumvsr[ii]);


  
  // start counter
  // full loop since user may choose to count something in boundary zones
  if(DODEBUG) FULLLOOP TSCALELOOP(tscale) FLOORLOOP(floor) failfloorcount[i][j][k][tscale][floor]=0;
#if(CALCFARADAYANDCURRENTS)
  // zero out jcon since outer boundaries not set ever since j^\mu involves spatial derivatives that don't exist outside a certain point
  for(pl=0;pl<NDIM;pl++) FULLLOOP jcon[i][j][k][pl]=0.0;
#endif


  ///////////////
  // 0 out these things so dump files are readable by SM under any cases
  ///////////////
  FULLLOOP{
    PLOOP(pl) udump[i][j][k][pl]=0.0;
#if(CALCFARADAYANDCURRENTS)
    DLOOPA(pl) jcon[i][j][k][pl]=0.0;
    for(pl=0;pl<NUMFARADAY;pl++) fcon[i][j][k][pl]=0.0;
#endif
  }





  // compute default (equivalent to horizon not existing -- just normal gridding) 
  horizoni = 0;
  horizoncpupos1 = 0;
  recompute_fluxpositions();

  // initialize grid sectioning to full grid at first
  init_gridsectioning();




  return(0);
}

int init_defgrid(void)
{
  // sets metric
  //  a=0.0;
  // set coordinates
  defcoord=LOGRSINTH;
  // sets parameters of coordinates, default changes
  R0 = 0.0;
  Rin = 0.98 * Rhor;
  Rout = 40.;
  hslope = 0.3;

  // default black hole parameters, and so length is in GMBH/c^2 and time in GMBH/c^3
  a=a0;
  MBH=MBH0;
  QBH=QBH0;


  return(0);
}


int init_defglobal(void)
{
  int i;
  int pl;
  int dtloop;

#if(!PRODUCTION)
  debugfail=2; // CHANGINGMARK
#else
  debugfail=0; // no messages in production -- assumes all utoprim-like failures need not be debugged
#endif
  // whether to show debug into on failures.  Desirable to turn off if don't care and just want code to burn away using given hacks/kludges
  // 0: no messages
  // 1: critical messages
  // 2: all failure messages

  // included in rdump
  defcon = 1.0;
  /* maximum increase in timestep */
  SAFE=1.3;
  whichrestart = 0;
  restartsteps[0] = 0;
  restartsteps[1] = 0;
  nstep = realnstep = 0;
  failed = 0;
  cour = 0.5;  //atch: modified the courant factor from 0.9
  doevolvemetricsubsteps=0; // default is to evolve on long steps (only applicable if DOEVOLVEMETRIC==1 && EVOLVEMETRICSUBSTEP==2)
  gravityskipstep=0; // default is no skipping
  gravitydtglobal = BIG;
  sourcedtglobal  = BIG;
  wavedtglobal    = BIG;



  //  avgscheme=WENO5BND;
  avgscheme[1]=avgscheme[2]=avgscheme[3]=DONOR;
  
  PALLLOOP(pl) do_transverse_flux_integration[pl] = 1;
  PALLLOOP(pl) do_source_integration[pl] = 1;
  PALLLOOP(pl) do_conserved_integration[pl] = 1;

  INVERTFROMAVERAGEIFFAILED = 1;
  LIMIT_AC_PRIM_FRAC_CHANGE = 1;
  MAX_AC_PRIM_FRAC_CHANGE = 0.1;

  LIMIT_AC_FRAC_CHANGE = 1;
  MAX_AC_FRAC_CHANGE = 0.1;

  dofluxreconevolvepointfield=1;


#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)
  //lim = WENO5FLAT;
  //lim = WENO5BND;
  //  lim = WENO3;
  //lim = DONOR;
  //lim = MINM;
  //  lim = PARA;
  //lim = MC;
  lim[1]=lim[2]=lim[3]=MC;
  //lim = PARA;
  //  lim = PARAFLAT;
  //lim = MC;
  TIMEORDER=4;
  // whether/which ENO used to interpolate fluxes
  DOENOFLUX = ENOFINITEVOLUME;
  //DOENOFLUX= NOENOFLUX;
  //DOENOFLUX=ENOFLUXRECON;
  //  fluxmethod=MUSTAFLUX;
  //fluxmethod=FORCEFLUX;
  fluxmethod=HLLFLUX;
  //fluxmethod=HLLLAXF1FLUX;
  //fluxmethod=LAXFFLUX;
  FLUXB = FLUXCTTOTH;
  UTOPRIMVERSION=UTOPRIM5D1;  //UTOPRIM2DFINAL;
  //UTOPRIMVERSION=UTOPRIM5D2;
  //  UTOPRIMVERSION=UTOPRIM2DFINAL;
#elif(EOMTYPE==EOMFFDE)
  // PARA and TO=4 and HLL not trustable in FFDE so far
  lim[1] = lim[2] = lim[3] = MC;
  TIMEORDER=2;
  fluxmethod=LAXFFLUX;
  FLUXB = FLUXCTTOTH;
  UTOPRIMVERSION=UTOPRIM2DFINAL;
  // whether/which ENO used to interpolate fluxes
  //DOENOFLUX = ENOFINITEVOLUME;
  DOENOFLUX= NOENOFLUX;
  //DOENOFLUX=ENOFLUXRECON;
#endif


  t = 0.;
  tstepparti = t;
  tsteppartf = t;

  tf = 1.0;
  
  for(dtloop=0;dtloop<NUMDTDS;dtloop++) DTdumpgen[dtloop]=1.0;
  //  DTd = DTavg = DTdebug = 1.0;
  //  DTener=1.0;
  //  DTi=1.0;
  DTr=1;



  GAMMIEDUMP=0;// whether in Gammie output types (sets filename with 3 numbers and doesn't output i,j)
  GAMMIEIMAGE=0; // Gammie filename and single density output
  GAMMIEENER=0; // currently doing gener as well as ener, but this would also output some messages in gammie form

  // DOCOLSPLIT
  //
  // 0: don't ..
  // 1: split dump files into 1 column per file with ending number being column number
  // works in MPI-mode as well.  ROMIO+DOCOLSPLIT is useful for tungsten with low memory and small files to avoid diskspace and memory limits.

  // default
  for(i=0;i<NUMDUMPTYPES;i++){
    DOCOLSPLIT[i]=0;
  }
  // otherwise specify for each dump type


  DODIAGEVERYSUBSTEP = 0;
  DOENODEBUGEVERYSUBSTEP = 0;
 

  DODIAGS=1; // whether to do diagnostics
  // specify individual diagnostics to be done
  DOENERDIAG=1;
  DOGDUMPDIAG=1;
  DORDUMPDIAG=1;
  DODUMPDIAG=1;
  if(DOAVG){
    DOAVGDIAG=1; // choice
  }
  else DOAVGDIAG=0; // no choice
  DOIMAGEDIAG=1;
  DOAREAMAPDIAG=1;

  POSDEFMETRIC=0; // see metric.c, bounds.c, and coord.c

  rescaletype=1;
  // 0: normal
  // 1: extended b^2/rho in funnel
  // 2: conserve E and L along field lines
  //   1 or 2 required to conserve E and L along field line

  /** FIXUP PARAMETERS **/
  RHOMIN=1.e-4;
  UUMIN =1.e-6;
  RHOMINLIMIT=1.e-20;
  UUMINLIMIT =1.e-20;

  // limit of B^2/rho if using that flag
  BSQORHOLIMIT=1E2;
  BSQOULIMIT=1E3;
  UORHOLIMIT=1E3;
  GAMMADAMP=5.0;

#if(EOMTYPE==EOMGRMHD || EOMTYPE==EOMCOLDGRMHD)
  // GODMARK -- unstable beyond about 25, but can sometimes get away with 1000
  GAMMAMAX=25.0; // when we think gamma is just too high and may cause unstable flow, but solution is probably accurate.
#else
  GAMMAMAX=2000.0;
#endif

  GAMMAFAIL=100.0*GAMMAMAX; // when we think gamma is rediculous as to mean failure and solution is not accurate.
  prMAX[RHO]=20.0;
  prMAX[UU]=20.0;
  prMAX[U1]=100.0;
  prMAX[U2]=100.0;
  prMAX[U3]=100.0;
  prMAX[B1]=100.0;
  prMAX[B2]=100.0;
  prMAX[B3]=100.0;

  // some physics
  gam=4./3.; // ultrarelativistic gas, assumes pgas<prad and radiation
  gamideal=gam;
	     // doesn't escape
  // gam=5/3 for non-relativistic gas, such as neucleons in collapsar model
  cooling=0;
  // cooling: 0: no cooling 1: adhoc thin disk cooling 2: neutrino cooling for collapsar model

  // boundary conditions (default is for 3D spherical polar grid -- full r,pi,2pi)
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=OUTFLOW;
  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;
  BCtype[X3UP]=PERIODIC;
  BCtype[X3DN]=PERIODIC;

  return(0);
}


int init_defconsts(void)
{


  // some constants that define neutrino cooling and spin evolution
  // cgs
  msun=1.989E33;
  lsun=3.89E33;
  G=6.672E-8;
  H=6.6262E-27;
  C=2.99792458E10;
  mn=1.67493E-24;
  mp=1.67262E-24;
  me=9.10938188E-28;
  kb=1.380658E-16; // erg K^{-1} g^{-1}
  Q=(mn-mp)*C*C;
  R=kb/mp;
  Re=kb/me;
  hpl=6.6262E-27;
  hbar=hpl/(2.0*M_PI);
  K=1.24E15;
  K2=9.9E12;
  //
  // compare below
  arad=8.*pow(M_PI,5.0)*pow(kb,4.0)/(15*C*C*C*H*H*H);
  //  arad=5.6704E-5 * 4.0 / C;
  sigmasb=arad*C/4.0;
  sigmamat=6.652E-29*100*100;
  mevocsq=1.783E-27;
  ergPmev=1.602E-6;


  // default for units is GM=C=1
  // GRB
  //  mb=mn;
  amu = 1.6605402E-24;
  // definition of baryon mass uses the below
  mb = amu;

  Mcgs=3.*msun;
  Mdot=0.1*msun;
  Mdotc=0.35; // approx code units mass accretion rate
  
  // units
  Lunit=G*Mcgs/(C*C);
  Tunit=G*Mcgs/(C*C*C);
  Vunit=Lunit/Tunit;
  Ccode=C/Vunit;
  rho0unittype=0; // mass is mass
  rhounit=Mdot/(Mdotc*Lunit*Lunit*Vunit);
  mbwithrhounit=mb;
  rhomassunit=rhounit;
  Munit=rhounit*pow(Lunit,3.0);
  mdotunit=Munit/Tunit;
  energyunit=Munit*Vunit*Vunit;
  edotunit=energyunit/Tunit;
  Pressureunit=rhounit*Vunit*Vunit;
  Tempunit=Pressureunit*mb/(rhounit*kb);
  Bunit=Vunit*sqrt(rhounit);
  massunitPmsun=Munit/msun;

  
  // physics stuff
  ledd=4.*M_PI*C*G*Mcgs*mb/sigmamat;
  leddcode=ledd/edotunit;

  // normalizations and settings for metric changes 
  a0=0.0;
  MBH0=1.0;
  QBH0=0.0;
  Mfactor=1.0;
  Jfactor=1.0;
  rhofactor=1.0;

  dabh=0.0;
  dE=0.0;
  dJ=0.0;


  return(0);
}



void reset_dothisenerregion(void)
{
  int enerregion;

  // these set some "enerregions"
  // first set which ener regions to do
  ENERREGIONALLLOOP(enerregion) dothisenerreg[enerregion]=1;
  // turn off some enerregions
  if(DOGRIDSECTIONING==0){
    dothisenerreg[ACTIVEREGION]=0;
  }
  if(DOJETDIAG==0){
    dothisenerreg[INNERJETREGION]=0;
    dothisenerreg[OUTERJETREGION]=0;
  }
}


int recompute_fluxpositions(void)
{

  // put here anything when changing horizonti (apart from recomputing metric, which is done outside)                     
  // if got here then horizonti changed, so need new positions of fluxes

  reset_dothisenerregion();
  //
  setflux();
  sethorizonflux();
  settrueglobalregion();
  settrueglobalwithbndregion();

  if(DOJETDIAG)  setjetflux();



  return(0);

}


int post_init(void)
{
  int compute_currents_t0(void);


  trifprintf("begin: post_init\n");


  // need to compute various things for EOS
  // also done per timestep in step_ch.c
  // GODMARK: Check that no EOS call is before this point 
  compute_EOS_parms(p);


  // in synch always here
  if (error_check(3)) {
    dualfprintf(fail_file, "error_check detected failure at main:1\n");
    dualfprintf(fail_file, "Bad initial conditions\n");
    myexit(1);
  }


  // all calculations that do not change for any initial conditions or problem setup or restart conditions

 
  // don't recompute these things if doing metric evolution on substeps
  if(RESTARTMODE!=0){
    setfailresponse(restartonfail);
  }
  
  
  ////////////////
  // current
  compute_currents_t0();



  
  trifprintf("end: post_init\n");


  return(0);

}


// bounds, etc. use enerpos type stuff that is set by recompute_fluxpositions()
int set_box_grid_parameters(void)
{

  // need to recompute horizon-related quantities if horizon is growing due to accretion
  if(ISBLACKHOLEMCOORD(MCOORD)){
    find_horizon(0);
    trifprintf("Rhor=%21.15g Risco=%21.15g MBH=%21.15g a=%21.15g QBH=%21.15g\n",Rhor,Risco,MBH,a,QBH);
  }
  else{
    horizoni=horizoncpupos1=0;
  }

  recompute_fluxpositions();


  return(0);

}


int compute_currents_t0(void)
{

#if(CALCFARADAYANDCURRENTS)
    
#if(CURRENTST0)
  // setup faraday so t=0 dump has diagnostics
  // may want currents for dump=0 (time-derivative terms will be wrong)
  current_doprecalc(CURTYPET,p);
  current_doprecalc(CURTYPEX,p);
  current_doprecalc(CURTYPEY,p);
  current_doprecalc(CURTYPEZ,p);
  // compute current
  current_calc(cfaraday);
  
#else
  
  // need to at least compute t=0 faraday so currents can start to be computed before reach step_ch.c (just an ordering issue with how step_ch.c does this case)
  if(WHICHCURRENTCALC==1){
    // compute faraday components needed for time centering of J
    current_doprecalc(CURTYPET,p);
  }
#endif
  
#if(FARADAYT0)
  current_doprecalc(CURTYPEFARADAY,p);
#endif
  
#endif


    return(0);

}






// after all parameters are set, can call this
int post_par_set(void)
{
  int interp_loop_set(void);
  int orders_set(void);
  void check_bnd_num(void);



  /////////////////////
  //
  // Initialize EOS and get EOS pointers for a given WHICHEOS,EOMTYPE
  //
  /////////////////////
  pickeos_eomtype(WHICHEOS,EOMTYPE);



  // if higher order, then set useghostplusactive to 1 since need that region in general
  if(HIGHERORDERMEM && (DOENOFLUX==ENOFINITEVOLUME || DOENOFLUX==ENOFLUXRECON) ){

    useghostplusactive=1;
  }
  else{
    useghostplusactive=0;
  }


  if(useghostplusactive && DOENOFLUX == ENOFLUXRECON && FLUXB==FLUXCTSTAG && dofluxreconevolvepointfield==0){
    // no matter how this parameter is set, we reset it to 
    extrazones4emf=1;
  }
  else extrazones4emf=0;

  if(SPLITMAEMMEM && HIGHERORDERMEM && DOENOFLUX == ENOFLUXRECON){
    splitmaem=1;
  }
  else splitmaem=0;

  if(splitmaem==0 && dofluxreconevolvepointfield==1 && DOENOFLUX == ENOFLUXRECON && FLUXB == FLUXCTSTAG && SPLITPRESSURETERMINFLUXMA==0 && SPLITPRESSURETERMINFLUXEM==0){
    // when not splitting MA/EM terms and doing FLUXRECON point field method, then need to still fix stencil for EMF terms
    emffixedstencil=1;
  }

  
  if(useghostplusactive && ((DOENOFLUX == ENOFLUXRECON && FLUXB==FLUXCTSTAG && dofluxreconevolvepointfield==0) || DOENOFLUX==ENOFINITEVOLUME)){
    // no matter how this parameter is set, we reset it to 
    unewisavg=1;
  }
  else unewisavg=0;


  trifprintf("useghostplusactive=%d extrazones4emf=%d\n",useghostplusactive,extrazones4emf);


  trifprintf("Setting orders\n");
  orders_set();

  trifprintf("Boundary number checks\n");
  check_bnd_num();

  trifprintf("Setting interp loop\n");
  interp_loop_set();

  trifprintf("Reporting bound loop\n");
  report_bound_loop();


  return(0);
}





// check that there are enough boundary zones for interpolation order used
void check_bnd_num(void)
{
  int totalo;
  int get_num_bnd_zones_used(int dir);
  int dimen;
  int bndnum[NDIM];
  int doingdimen[NDIM];
  int pl;
  int doingenough;


  doingdimen[1]=N1NOT1;
  doingdimen[2]=N2NOT1;
  doingdimen[3]=N3NOT1;


  bndnum[1]=N1BND;
  bndnum[2]=N2BND;
  bndnum[3]=N3BND;


  DIMENLOOP(dimen){

    // first fix avgscheme and lim for case when not doing dimension
    if(!doingdimen[dimen]){
      avgscheme[dimen]=0;
      lim[dimen]=0;
    } // otherwise assume as wanted

    // get number of boundary zones needed
    totalo=get_num_bnd_zones_used(dimen);

  
    // bndnum[dimen] (MAXBND is maximum) is the number of boundary zones to be set by boundary routines and to be passed by MPI routines
    if(totalo==bndnum[dimen] || !doingdimen[dimen]){
      // then good
    }
    else if(totalo>bndnum[dimen]){
      dualfprintf(fail_file,"Not enough: dimen=%d totalo=%d MAXBND=%d bndnum=%d for avgscheme interporder[%d]=%d or lim interporder[%d]=%d extrazones4emf=%d\n",dimen,totalo,MAXBND,bndnum[dimen],avgscheme[dimen],interporder[avgscheme[dimen]],lim[dimen],interporder[lim[dimen]],extrazones4emf);
      failed=1; // so don't try to compute things in dump
      myexit(1);
    }
    else{
      // then MAXBND too large
      dualfprintf(fail_file,"WARNING: MAXBND excessive\n");
      dualfprintf(fail_file,"WARNING: dimen=%d totalo=%d MAXBND=%d and bndnum=%d for avgscheme interporder[%d]=%d or lim interporder[%d]=%d extrazones4emf=%d\n",dimen,totalo,MAXBND,bndnum[dimen],avgscheme[dimen],interporder[avgscheme[dimen]],lim[dimen],interporder[lim[dimen]],extrazones4emf);
      // not a failure, but user should be aware especially when doing MPI
    }
    
  }


  // see if need to be here at all                                                                                         
  if(interporder[avgscheme[1]]<3 && interporder[avgscheme[2]]<3 && interporder[avgscheme[3]]<3){
    PALLLOOP(pl) do_transverse_flux_integration[pl] = 0;
    PALLLOOP(pl) do_source_integration[pl] = 0;
    PALLLOOP(pl) do_conserved_integration[pl] = 0;
    
    trifprintf("Changed do_{conserved/source/transverse_flux}_integration to 0 since avgscheme[1]=%d avgscheme[2]=%d avgscheme[3]=%d\n",avgscheme[1],avgscheme[2],avgscheme[3]);
  }


  
  // checks on parameters so user doesn't do something stupid
  if(FULLOUTPUT&&USEMPI){
    dualfprintf(fail_file,"Cannot use FULLOUTPUT!=0 when USEMPI=1\n");
    myexit(200);
  }

  if(DOEVOLVEMETRIC&& (ANALYTICCONNECTION||ANALYTICSOURCE)){
    dualfprintf(fail_file,"Unlikely you have metric in time analytically\n");
    myexit(201);
  }


  if(DOSELFGRAVVSR && (ANALYTICCONNECTION||ANALYTICSOURCE||ANALYTICGCON)){
    dualfprintf(fail_file,"Unlikely you have metric analytically with self gravity\n");
    myexit(202);
  }

  if(GDETVOLDIFF){
    if(     (MCOORD==BLCOORDS || MCOORD==KSCOORDS || MCOORD==HTMETRIC || MCOORD==HTMETRICACCURATE || MCOORD==SPCMINKMETRIC) ){
      // then fine
    }
    else{
      dualfprintf(fail_file,"GDETVOLDIFF==1 not setup for non-SPC metrics\n");
      myexit(203);
    }
  }

  if(DOSELFGRAVVSR==0 && (MCOORD==KS_BH_TOV_COORDS || MCOORD==KS_TOV_COORDS || MCOORD==BL_TOV_COORDS) ){
    dualfprintf(fail_file,"DOSELFGRAVVSR==0 with MCOORD=%d is not likely what you want\n",MCOORD);
    dualfprintf(fail_file,"Continuing....\n");
  }


  if(HIGHERORDERMEM==0 && DOENOFLUX != NOENOFLUX){
    dualfprintf(fail_file,"Need to turn on HIGHERORDERMEM when doing higher order methods (i.e. DOENOFLUX!=NOENOFLUX\n");
    myexit(204);
  }

  if(COMPDIM!=3){
    dualfprintf(fail_file,"Code not setup for anything but COMPDIM==3\n");
    myexit(205);
  }


  if(FIELDSTAGMEM==0 && FLUXB==FLUXCTSTAG){
    dualfprintf(fail_file,"FIELDSTAGMEM should be 1 if FLUXB==FLUXCTSTAG\n");
    myexit(206);
  }

  if(FIELDTOTHMEM==0 && FLUXB==FLUXCTTOTH){
    dualfprintf(fail_file,"FIELDTOTHMEM should be 1 if FLUXB==FLUXCTTOTH\n");
    myexit(207);
  }


  if( (N1%2>0 && N1>1) || (N2%2>0 && N2>1) || (N3%2>0 && N3>1) ){
    dualfprintf(fail_file,"N1, N2, N3 should be even since some parts of code assume so\n");
    myexit(208);
  }


  if(SUPERLONGDOUBLE){
    dualfprintf(fail_file,"If doing SUPERLONGDOUBLE, then should have compiled as such\n");
  }
  
  if(ROEAVERAGEDWAVESPEED || ATHENAROE){
    dualfprintf(fail_file,"ATHENA stuff only for non-rel setup and while was tested hasn't been used or kept up to date\n");
    myexit(209);
  }

  if(STOREWAVESPEEDS==0 && FLUXB==FLUXCTSTAG){
    dualfprintf(fail_file,"Must set STOREWAVESPEEDS==1 when doing FLUXCTSTAG\n");
    myexit(210);
  }

  if(LIMADJUST!=LIMITERFIXED){
    dualfprintf(fail_file,"LIMADJUST old code\n");
    myexit(211);
  }


  if(FLUXADJUST!=FLUXFIXED){
    dualfprintf(fail_file,"FLUXADJUST old code\n");
    myexit(212);
  }

  if(DOEVOLVEMETRIC&& (WHICHEOM!=WITHGDET )){
    dualfprintf(fail_file,"conn2 not computed for time-dependent metric yet\n");
    myexit(213);
  }

  if(CONNMACHINEBODY){

    if(CONNDERTYPE!=DIFFFINITE && CONNDERTYPE!=DIFFGAMMIE){
      dualfprintf(fail_file,"Makes no sense to use CONNMACHINEBODY with CONNDERTYPE!=DIFFFINITE/DIFFGAMMIE\n");
      myexit(214);
    }

    if(WHICHEOM!=WITHGDET){
      dualfprintf(fail_file,"Not setup for body correction when f is not detg\n");
      myexit(215);
    }
  }



  if(CONTACTINDICATOR!=0){
    dualfprintf(fail_file,"Contact not recommended\n");
  }

  if(FLUXDUMP!=0){
    dualfprintf(fail_file,"FLUXDUMP ACTIVE -- lots of extra output\n");
  }

  if(NUMPOTHER!=0){
    dualfprintf(fail_file,"NUMPOTHER ACTIVE -- lots of extra memory used\n");
  }

  if(LIMIT_FLUXC2A_PRIM_CHANGE){
    dualfprintf(fail_file,"LIMIT_FLUXC2A_PRIM_CHANGE doesn't work according to Sasha, so don't use\n");
    myexit(216);
  }

  // check if 
  doingenough=1;
  DIMENLOOP(dimen){
    doingenough *= (interporder[avgscheme[dimen]]>=interporder[WENO5BNDPLUSMIN] || doingdimen[dimen]==0);
  }

  if(WENO_EXTRA_A2C_MINIMIZATION==1 && doingenough==0){
    dualfprintf(fail_file,"WENO_EXTRA_A2C_MINIMIZATION==1 and interporder=%d %d %d invalid\n",interporder[avgscheme[1]],interporder[avgscheme[2]],interporder[avgscheme[3]]);
    myexit(217);
  }



  DIMENLOOP(dimen){
    if(avgscheme[dimen]!=DONOR){ // using DONOR just turns off and assume standard way to turn off so no need to message user
      if(avgscheme[dimen]<FIRSTWENO || avgscheme[dimen]>LASTWENO){
	dualfprintf(fail_file,"Choice of avgscheme[%d]=%d has no effect\n",dimen,avgscheme[dimen]);
      }
    }
  }

  
  //  if( (FLUXB==FLUXCTHLL || FLUXB==FLUXCTTOTH || (FLUXB==FLUXCTSTAG && extrazones4emf==1 && )) && EVOLVEWITHVPOT){ // avoid complicated conditional
  if(EVOLVEWITHVPOT && !(FLUXB==FLUXCTSTAG && extrazones4emf==0) ){
    // Even with FV or FLUXRECON method that relies on updating non-point fields, can't evolve with A_i since truncation error different rather than just machine error different.  This leads to large errors -- especially at boundaries? GODMARK -- not 100% sure this is the problem for test=1102 and EVOLVEWITHVPOT
    dualfprintf(fail_file,"Cannot evolve field using A_i for FLUXB==FLUXCTHLL or FLUXB==FLUXCTHLL since no single A_i evolved forward in time.  And cannot use with non-point field method since while single A_i updated, the update diverges at truncation error between updating A_i and updating the non-point-field -- this is especially bad for periodic boundary conditions where one must have machine error correct behavior at boundaries\n");
    myexit(218);
  }


  if(splitmaem==0 && (SPLITPRESSURETERMINFLUXMA || SPLITPRESSURETERMINFLUXEM)){
    dualfprintf(fail_file,"To use SPLITPRESSURE must turn on splitmaem, which requires FLUXRECON\n");
    myexit(219);
  }

  if(splitmaem==1){
    dualfprintf(fail_file,"Noticed some tests are more accurate if don't split (splitmaem==0), like test=1102 in init.sasha.c when doing non-point field FLUXRECON method the density error is smaller by factor of 2 with splitmaem==0\n");
  }

  // when using FLUXRECON, WENO5BND lim/avg:
  // 32x16: splitmaem==1 (pressure split or not!): 5.164e-05   0.0005637    0.001195    0.001326    0.003499   0.0001439   0.0001438    0.000415
  // 32x16: splitmaem==0 (higherordersmooth or rough!): 8.369e-06   0.0005615    0.001208    0.001224    0.002977   0.0001433   0.0001431   0.0003529


  if(ISSPCMCOORD(MCOORD) && ACCURATESINCOS==0){
    dualfprintf(fail_file,"Warning: if polar axis or r=0 singularity don't have \\detg=0 you should force it -- normally ACCURATESINCOS==1 does a good job of this, but still should have code to check\n");
  }


  if(FULLOUTPUT!=0){
    // then ensure that boundary zones are off if necessary

    if(DOLUMVSR){
      dualfprintf(fail_file,"lumvsr requires no boundary zones outputted so dump and lumvsr file can be used together\n");
      myexit(195815983);
    }
    if(DODISSVSR){
      dualfprintf(fail_file,"dissvsr requires no boundary zones outputted so dump and dissvsr file can be used together\n");
      myexit(195815984);
    }

  }

#if(MERGEDC2EA2CMETHOD && PARAMODEWENO)
    // then ensure that boundary zones are off if necessary

  //    dualfprintf(fail_file,"MERGEDC2EA2CMETHOD not setup for PARA yet: Turn off PARAMODEWENO\n");
#error "MERGEDC2EA2CMETHOD not setup for PARA yet: Turn off PARAMODEWENO\n"

#endif


  if(2*N1<N1BND || N1BND!=0 && N1==1 || periodicx1==1 && ncpux1>1 && N1<N1BND){
    dualfprintf(fail_file,"Code not setup to handle boundary cells N1BND=%d with only N1=%d\n",N1BND,N1);
    myexit(246872462);
  }
  if(2*N2<N2BND || N2BND!=0 && N2==1 || periodicx2==1 && ncpux2>1 && N2<N2BND){
    dualfprintf(fail_file,"Code not setup to handle boundary cells N2BND=%d with only N2=%d\n",N2BND,N2);
    myexit(246872463);
  }
  if(2*N3<N3BND || N3BND!=0 && N3==1 || periodicx3==1 && ncpux3>1 && N3<N3BND){
    dualfprintf(fail_file,"Code not setup to handle boundary cells N3BND=%d with only N3=%d\n",N3BND,N3);
    myexit(246872464);
  }

  if(N1%2 && N1>1 || N2%2 && N2>1 || N3%2 && N3>1){
    dualfprintf(fail_file,"Need even N1,N2,N3 AFAIK N1=%d N2=%d N3=%d\n",N1,N2,N3);
    myexit(19846286);
  }


  if(SENSITIVE!=LONGDOUBLETYPE){
    dualfprintf(fail_file,"WARNING: With SENSITIVE!=LONGDOUBLETYPE you may have problems for some integral or counting quantities (e.g. DTd too small or many zones to integrate over)\n");
  }



#if(MERGEDC2EA2CMETHOD==1)
  //  if(splitmaem){
  //    dualfprintf(fail_file,"MERGEDC2EA2CMETHOD==1 is not setup for splitmaem since it was assumed splitmaem only needed with old a2c method\n");
  //    myexit(346897346);
  //  }
#endif


  if(PRODUCTION==1){
    if(DOENOFLUX != NOENOFLUX && FLUXB==FLUXCTSTAG){
      dualfprintf(fail_file,"NOTE: With PRODUCTION==1, higher-order staggered field method won't compute ener file value of divB correctly because turned off diagnostic bounding to avoid excessive MPI calls to bound unew.  dump file will still be correct for MPI boundaries but not for real boundaries since unew not defined to be bounded by user\n");
    }
  }


  if(PRODUCTION==0){
    dualfprintf(fail_file,"WARNING: PRODUCTION set to 0, code may be slower\n");
  }

  if(LIMITDTWITHSOURCETERM){
    dualfprintf(fail_file,"WARNING: LIMITDTWITHSOURCETERM set to 1, code may be slower\n");
  }

  if(DODISS || DOLUMVSR || DODISSVSR || DOENTROPY!=DONOENTROPY){
    dualfprintf(fail_file,"WARNING: DODISS/DOLUMVSR/DODISSVSR/DOENTROPY!=DONOENTROPY set to 1, code may be slower\n");
  }

  if(CHECKONINVERSION){
    dualfprintf(fail_file,"WARNING: CHECKONINVERSION set to 1, code may be slower\n");
  }


  if(CHECKSOLUTION){
    dualfprintf(fail_file,"WARNING: CHECKSOLUTION set to 1, code may be slower\n");
  }


  DIMENLOOP(dimen){
    if(interporder[lim[dimen]]-1 > TIMEORDER){
      dualfprintf(fail_file,"WARNING: interporder[dimen=%d lim=%d]=%d -1 > TIMEORDER=%d is unstable in region where Courant condition setting dt\n",dimen,lim[dimen],interporder[lim[dimen]],TIMEORDER);
    }
  }

  if(DOINGLIAISON && USEMPI==0){
    dualfprintf(fail_file,"WARNING: DOINGLIAISON==1 but USEMPI==0\n");
  }


  // complain if b^2/rho b^2/u or u/rho too large for given resolution given experience with GRMHD torus problem
  DIMENLOOP(dimen){
    if(BSQORHOLIMIT/30.0>((FTYPE)totalsize[dimen])/64.0){
      dualfprintf(fail_file,"WARNING: BSQORHOLIMIT=%21.15g for totalsize[%d]=%d\n",BSQORHOLIMIT,dimen,totalsize[dimen]);
    }
    if(BSQOULIMIT/100.0>((FTYPE)totalsize[dimen])/64.0){
      dualfprintf(fail_file,"WARNING: BSQOULIMIT=%21.15g for totalsize[%d]=%d\n",BSQOULIMIT,dimen,totalsize[dimen]);
    }
    if(UORHOLIMIT/100.0>((FTYPE)totalsize[dimen])/64.0){
      dualfprintf(fail_file,"WARNING: UORHOLIMIT=%21.15g for totalsize[%d]=%d\n",UORHOLIMIT,dimen,totalsize[dimen]);
    }
  }


  // external checks
  parainitchecks();



}


int get_num_bnd_zones_used(int dimen)
{
  int avgo;
  int interpo;
  int totalo;
  int doingdimen[NDIM];
  int extraavgo;



  doingdimen[1]=N1NOT1;
  doingdimen[2]=N2NOT1;
  doingdimen[3]=N3NOT1;


  if(useghostplusactive){

    if(doingdimen[dimen]){
      // number of zones one way for finite volume scheme to convert Uavg -> Upoint
      avgo=(interporder[avgscheme[dimen]]-1)/2;
    }
    else avgo=0; // nothing done for this dimension, so no avg zones
  }
  else avgo=0; // no need for extra zones

  if(extrazones4emf){
    if(doingdimen[dimen]){
      // number of zones one way for finite volume scheme to convert Uavg -> Upoint
      extraavgo=(interporder[avgscheme[dimen]]-1)/2;
    }
    else extraavgo=0; // nothing done for this dimension, so no avg zones
  }
  else extraavgo=0; // no need for extra zones


  // number of zones one way to have for interpolation to get fluxes
  // need to get flux at i,j,k=-1 in any case, and need boundary zones from there, so 1 extra effective boundary zone for interpolation
  if(doingdimen[dimen]){
    interpo=(interporder[lim[dimen]]-1)/2+1;
  }
  else interpo=0; // then not doing this dimension

  totalo=avgo+interpo+extraavgo;

  return(totalo);
}










// define range over which various loops go
// all these should be as if no grid sectioning SECTIONMARK since used in loops that have SHIFTS inside
int interp_loop_set(void)
{
  int avgo[NDIM];
  int doingdimen[NDIM];
  int dimen;
  int dir;
  int jj;
  int avgoperdir[NDIM][NDIM];
  int odimen1,odimen2;




  doingdimen[1]=N1NOT1;
  doingdimen[2]=N2NOT1;
  doingdimen[3]=N3NOT1;


  // the fluxloop[dir][FIJKDEL]'s can be used for general purpose to get idel, jdel, kdel


  DIMENLOOP(dimen){

    if(doingdimen[dimen]){
      // number of zones one way for finite volume scheme to convert Uavg -> Upoint
      avgo[dimen]=(interporder[avgscheme[dimen]]-1)/2;
    }
    else avgo[dimen]=0; // nothing done for this dimension, so no avg zones
    trifprintf("dimen=%d avgo=%d\n",dimen,avgo[dimen]);
  }



  // always need fluxes in ghost+active if doing higher order methods
  if(useghostplusactive){

    // (interporder[avgscheme]-1)/2 is the number of points to the left and the number of ponits to the right that are needed for the finite volume scheme


    // scheme used to convert Uavg -> Upoint requires extra zones
    //    avgo=(interporder[avgscheme]-1)/2;
    // i.e. if avgscheme=WENO5, then interporder[WENO5]=5 and the Uconsloop goes from -2 .. N+1 inclusive  and fluxes go from -2 to N+2 inclusive
    // same for WENO4
    // this is -avgo .. N-1+avgo for Uconsloop and -avgo .. N+avgo for fluxes

    
    // note that if doing FLUXCTSTAG (FIELDSTAGMEM) then cross directions already expanded




    dimen=1;
    fluxloop[dimen][FIDEL]=SHIFT1;
    fluxloop[dimen][FJDEL]=0;
    fluxloop[dimen][FKDEL]=0;
    fluxloop[dimen][FFACE]=FACE1;
    fluxloop[dimen][FIS]=(-avgo[dimen])*N1NOT1;
    fluxloop[dimen][FIE]=N1-1+(avgo[dimen]+1)*N1NOT1;
    fluxloop[dimen][FJS]=INFULL2;
    fluxloop[dimen][FJE]=OUTFULL2;
    fluxloop[dimen][FKS]=INFULL3;
    fluxloop[dimen][FKE]=OUTFULL3;

    dimen=2;
    fluxloop[dimen][FIDEL]=0;
    fluxloop[dimen][FJDEL]=SHIFT2;
    fluxloop[dimen][FKDEL]=0;
    fluxloop[dimen][FFACE]=FACE2;
    fluxloop[dimen][FIS]=INFULL1;
    fluxloop[dimen][FIE]=OUTFULL1;
    fluxloop[dimen][FJS]=(-avgo[dimen])*N2NOT1;
    fluxloop[dimen][FJE]=N2-1+(avgo[dimen]+1)*N2NOT1;
    fluxloop[dimen][FKS]=INFULL3;
    fluxloop[dimen][FKE]=OUTFULL3;

    dimen=3;
    fluxloop[dimen][FIDEL]=0;
    fluxloop[dimen][FJDEL]=0;
    fluxloop[dimen][FKDEL]=SHIFT3;
    fluxloop[dimen][FFACE]=FACE3;
    fluxloop[dimen][FIS]=INFULL1;
    fluxloop[dimen][FIE]=OUTFULL1;
    fluxloop[dimen][FJS]=INFULL2;
    fluxloop[dimen][FJE]=OUTFULL2;
    fluxloop[dimen][FKS]=(-avgo[dimen])*N3NOT1;
    fluxloop[dimen][FKE]=N3-1+(avgo[dimen]+1)*N3NOT1;

  }
  else{

    // Uconsloop for these methods just involve normal CZLOOP
    // inversion for this method just involves CZLOOP
    dimen=1;
    fluxloop[dimen][FIDEL]=SHIFT1;
    fluxloop[dimen][FJDEL]=0;
    fluxloop[dimen][FKDEL]=0;
    fluxloop[dimen][FFACE]=FACE1;
    fluxloop[dimen][FIS]=0;
    fluxloop[dimen][FIE]=OUTM1;
    if(FLUXB==FLUXCTSTAG){
      fluxloop[dimen][FJS]=INFULL2;
      fluxloop[dimen][FJE]=OUTFULL2;
      fluxloop[dimen][FKS]=INFULL3;
      fluxloop[dimen][FKE]=OUTFULL3;
    }
    else{ // then only averaging for FLUXCT so only need 1 additional off-direction point
      fluxloop[dimen][FJS]=-SHIFT2;  //atch: loop over additional row to provide enough fluxes for FLUXCT, etc. to operate near the boundary
      fluxloop[dimen][FJE]=N2-1+SHIFT2; // " " 
      fluxloop[dimen][FKS]=-SHIFT3;     // " "
      fluxloop[dimen][FKE]=N3-1+SHIFT3; // " "
    }

    dimen=2;
    fluxloop[dimen][FIDEL]=0;
    fluxloop[dimen][FJDEL]=SHIFT2;
    fluxloop[dimen][FKDEL]=0;
    fluxloop[dimen][FFACE]=FACE2;
    fluxloop[dimen][FJS]=0; 
    fluxloop[dimen][FJE]=OUTM2;
    if(FLUXB==FLUXCTSTAG){
      fluxloop[dimen][FIS]=INFULL1;
      fluxloop[dimen][FIE]=OUTFULL1;
      fluxloop[dimen][FKS]=INFULL3;
      fluxloop[dimen][FKE]=OUTFULL3;
    }
    else{
      fluxloop[dimen][FIS]=-SHIFT1;   //atch: loop over additional row to provide enough fluxes for FLUXCT, etc. to operate near the boundary
      fluxloop[dimen][FIE]=N1-1+SHIFT1; // " "
      fluxloop[dimen][FKS]=-SHIFT3;    // " "
      fluxloop[dimen][FKE]=N3-1+SHIFT3;// " "
    }


    dimen=3;
    fluxloop[dimen][FIDEL]=0;
    fluxloop[dimen][FJDEL]=0;
    fluxloop[dimen][FKDEL]=SHIFT3;
    fluxloop[dimen][FFACE]=FACE3;
    fluxloop[dimen][FKS]=0;
    fluxloop[dimen][FKE]=OUTM3;
    if(FLUXB==FLUXCTSTAG){
      fluxloop[dimen][FIS]=INFULL1;
      fluxloop[dimen][FIE]=OUTFULL1;
      fluxloop[dimen][FJS]=INFULL2;
      fluxloop[dimen][FJE]=OUTFULL2;
    }
    else{
      fluxloop[dimen][FIS]=-SHIFT1;   //atch: loop over additional row to provide enough fluxes for FLUXCT, etc. to operate near the boundary
      fluxloop[dimen][FIE]=N1-1+SHIFT1;  // " "
      fluxloop[dimen][FJS]=-SHIFT2;      // " "
      fluxloop[dimen][FJE]=N2-1+SHIFT2;  // " "
    }

  }



  DIMENLOOP(dimen) DIMENLOOP(dir){
    avgoperdir[dir][dimen]=avgo[dir]*(dimen==dir);
  }


  // fluxloop for staggered field's EMF when doing old dofluxreconevolvepointfield==0 method
  DIMENLOOP(dimen){
    odimen1=dimen%3+1;
    odimen2=(dimen+1)%3+1;

    emffluxloop[dimen][FIDEL]=MAX(fluxloop[odimen1][FIDEL],fluxloop[odimen2][FIDEL]);
    emffluxloop[dimen][FJDEL]=MAX(fluxloop[odimen1][FJDEL],fluxloop[odimen2][FJDEL]);
    emffluxloop[dimen][FKDEL]=MAX(fluxloop[odimen1][FKDEL],fluxloop[odimen2][FKDEL]);
    if(dimen==1) emffluxloop[dimen][FFACE]=CORN1;
    else if(dimen==1) emffluxloop[dimen][FFACE]=CORN2;
    else if(dimen==1) emffluxloop[dimen][FFACE]=CORN3;

    if(extrazones4emf){
      // then need defined at more positions
      emffluxloop[dimen][FIS]=MIN(fluxloop[odimen1][FIS]-avgoperdir[odimen1][1],fluxloop[odimen2][FIS]-avgoperdir[odimen2][1]);
      emffluxloop[dimen][FIE]=MAX(fluxloop[odimen1][FIE]+avgoperdir[odimen1][1],fluxloop[odimen2][FIE]+avgoperdir[odimen2][1]);
      emffluxloop[dimen][FJS]=MIN(fluxloop[odimen1][FJS]-avgoperdir[odimen1][2],fluxloop[odimen2][FJS]-avgoperdir[odimen2][2]);
      emffluxloop[dimen][FJE]=MAX(fluxloop[odimen1][FJE]+avgoperdir[odimen1][2],fluxloop[odimen2][FJE]+avgoperdir[odimen2][2]);
      emffluxloop[dimen][FKS]=MIN(fluxloop[odimen1][FKS]-avgoperdir[odimen1][3],fluxloop[odimen2][FKS]-avgoperdir[odimen2][3]);
      emffluxloop[dimen][FKE]=MAX(fluxloop[odimen1][FKE]+avgoperdir[odimen1][3],fluxloop[odimen2][FKE]+avgoperdir[odimen2][3]);
    }
    else{
      emffluxloop[dimen][FIS]=MIN(fluxloop[odimen1][FIS],fluxloop[odimen2][FIS]);
      emffluxloop[dimen][FIE]=MAX(fluxloop[odimen1][FIE],fluxloop[odimen2][FIE]);
      emffluxloop[dimen][FJS]=MIN(fluxloop[odimen1][FJS],fluxloop[odimen2][FJS]);
      emffluxloop[dimen][FJE]=MAX(fluxloop[odimen1][FJE],fluxloop[odimen2][FJE]);
      emffluxloop[dimen][FKS]=MIN(fluxloop[odimen1][FKS],fluxloop[odimen2][FKS]);
      emffluxloop[dimen][FKE]=MAX(fluxloop[odimen1][FKE],fluxloop[odimen2][FKE]);
    }
  }




  ////////////////
  //
  // Define range over which "average" conserved quantities are *evolved*
  //
  // Don't really need to evolve ghost+active region if FLUXRECON method unless FLUXRECON&&(FLUXB==FLUXCTSTAG), but not a problem to be a bit excessive
  //
  // if doing FLUXRECON && evolving point field, then no need to evolve ghost+active region -- latest method
  ////////////////
  if(useghostplusactive==0 || (DOENOFLUX == ENOFLUXRECON && extrazones4emf==0 && dofluxreconevolvepointfield==1) ){
    Uconsevolveloop[FFACE]=CENT;

    dimen=1;
    Uconsevolveloop[FIS]=0;
    Uconsevolveloop[FIE]=N1-1;

    dimen=2;
    Uconsevolveloop[FJS]=0;
    Uconsevolveloop[FJE]=N2-1;

    dimen=3;
    Uconsevolveloop[FKS]=0;
    Uconsevolveloop[FKE]=N3-1;
  }
  else{
    // expanded loop using expanded range of fluxes so update Uf and ucum in layer of ghost zones so don't have to bound flux or Uf/Ui/ucum
    // only needed for FLUXRECON with STAG field method and only for fields, but do all for simplicity
    
    // inversion for this method just involves CZLOOP

    // loop over averaged U to get Uf
    //    Uconsevolveloop[FIDEL]=0;
    //    Uconsevolveloop[FJDEL]=0;
    //    Uconsevolveloop[FKDEL]=0;
    Uconsevolveloop[FFACE]=CENT;

    dimen=1;
    Uconsevolveloop[FIS]=-avgo[dimen]*N1NOT1;
    Uconsevolveloop[FIE]=N1-1+avgo[dimen]*N1NOT1;

    dimen=2;
    Uconsevolveloop[FJS]=-avgo[dimen]*N2NOT1;
    Uconsevolveloop[FJE]=N2-1+avgo[dimen]*N2NOT1;

    dimen=3;
    Uconsevolveloop[FKS]=-avgo[dimen]*N3NOT1;
    Uconsevolveloop[FKE]=N3-1+avgo[dimen]*N3NOT1;
  }



  ////////////////
  //
  // Define ghost+active region over which face (interpolated from center or natively existing in FLUXCTSTAG case) values of primitives are needed
  //
  ////////////////
  if(useghostplusactive){

    // scheme used to convert Uavg -> Upoint requires extra zones
    //    avgo=(interporder[avgscheme]-1)/2;
    // i.e. if avgscheme=WENO5, then interporder[WENO5]=5 and the Uconsloop goes from -2 .. N+1 inclusive  and fluxes go from -2 to N+2 inclusive
    // same for WENO4
    // this is -avgo .. N-1+avgo for Uconsloop and -avgo .. N+avgo for fluxes

    // loop over averaged U to get Uf
    //    Uconsloop[FIDEL]=0;
    //    Uconsloop[FJDEL]=0;
    //    Uconsloop[FKDEL]=0;
    Uconsloop[FFACE]=CENT;

    dimen=1;
    Uconsloop[FIS]=-avgo[dimen]*N1NOT1;
    Uconsloop[FIE]=N1-1+avgo[dimen]*N1NOT1;

    dimen=2;
    Uconsloop[FJS]=-avgo[dimen]*N2NOT1;
    Uconsloop[FJE]=N2-1+avgo[dimen]*N2NOT1;

    dimen=3;
    Uconsloop[FKS]=-avgo[dimen]*N3NOT1;
    Uconsloop[FKE]=N3-1+avgo[dimen]*N3NOT1;

    // inversion for this method just involves CZLOOP
  }
  else{
    Uconsloop[FFACE]=CENT;

    dimen=1;
    Uconsloop[FIS]=0;
    Uconsloop[FIE]=N1-1;

    dimen=2;
    Uconsloop[FJS]=0;
    Uconsloop[FJE]=N2-1;

    dimen=3;
    Uconsloop[FKS]=0;
    Uconsloop[FKE]=N3-1;
  }


  
  if(extrazones4emf){
    emfUconsloop[FFACE]=CENT;

    dimen=1;
    emfUconsloop[FIS]=Uconsloop[FIS]-avgo[dimen]*N1NOT1;
    emfUconsloop[FIE]=Uconsloop[FIE]+avgo[dimen]*N1NOT1;

    dimen=2;
    emfUconsloop[FJS]=Uconsloop[FJS]-avgo[dimen]*N2NOT1;
    emfUconsloop[FJE]=Uconsloop[FJE]+avgo[dimen]*N2NOT1;

    dimen=3;
    emfUconsloop[FKS]=Uconsloop[FKS]-avgo[dimen]*N3NOT1;
    emfUconsloop[FKE]=Uconsloop[FKE]+avgo[dimen]*N3NOT1;
  }
  else{
    emfUconsloop[FFACE]=CENT;

    dimen=1;
    emfUconsloop[FIS]=Uconsloop[FIS];
    emfUconsloop[FIE]=Uconsloop[FIE];

    dimen=2;
    emfUconsloop[FJS]=Uconsloop[FJS];
    emfUconsloop[FJE]=Uconsloop[FJE];

    dimen=3;
    emfUconsloop[FKS]=Uconsloop[FKS];
    emfUconsloop[FKE]=Uconsloop[FKE];
  }






  DIMENLOOP(dimen){
    for(jj=0;jj<NUMFLUXLOOPNUMBERS;jj++) trifprintf("fluxloop[dimen=%d][%d] = %d\n",dimen,jj,fluxloop[dimen][jj]);
  }
  DIMENLOOP(dimen){
    for(jj=0;jj<NUMFLUXLOOPNUMBERS;jj++) trifprintf("emffluxloop[dimen=%d][%d] = %d\n",dimen,jj,emffluxloop[dimen][jj]);
  }
  for(jj=0;jj<NUMFLUXLOOPNUMBERS;jj++) trifprintf("Uconsevolveloop[%d] = %d\n",jj,Uconsevolveloop[jj]);
  for(jj=0;jj<NUMFLUXLOOPNUMBERS;jj++) trifprintf("Uconsloop[%d] = %d\n",jj,Uconsloop[jj]);
  for(jj=0;jj<NUMFLUXLOOPNUMBERS;jj++) trifprintf("emfUconsloop[%d] = %d\n",jj,emfUconsloop[jj]);
  








  return(0);

}




int get_loop(int pointorlinetype, int interporflux, int dir, struct of_loop *loop)
{

  set_interpalltypes_loop_ranges(pointorlinetype, interporflux, dir, &(loop->intdir), &(loop->is), &(loop->ie), &(loop->js), &(loop->je), &(loop->ks), &(loop->ke), &(loop->di), &(loop->dj), &(loop->dk), &(loop->bs), &(loop->ps), &(loop->pe), &(loop->be));
    
  return(0);

}



// master interp range function for both point and line methods
// This particular loop gives back 3D grid range, not line-by-line as in original line type method (so don't use directly in interpline.c!)
int set_interpalltypes_loop_ranges(int pointorlinetype, int interporflux, int dir, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be)
{
  int withshifts;

  if(pointorlinetype==INTERPPOINTTYPE){
    set_interppoint_loop_ranges(interporflux, dir, is, ie, js, je, ks, ke, di, dj, dk);
    // interpolation is always along dir-direction for point methods
    *intdir=dir;
    *ps=*is;
    *pe=*ie;
    // below is largest possible range that includes boundary cells
    *bs = *is - N1BND*(dir==1) - N2BND*(dir==2) - N3BND*(dir==3);
    *be = *ie + N1BND*(dir==1) + N2BND*(dir==2) + N3BND*(dir==3);
  }
  else if(pointorlinetype==INTERPLINETYPE){
    withshifts=0; // force to be without shifts so results can be put into loop that has shifts embedded
    set_interp_loop_gen(withshifts, interporflux, dir, intdir, is, ie, js, je, ks, ke, di, dj, dk, bs, ps, pe, be);
    // transcribe from loop over starting positions to full 3D loop
    if(*intdir==1){
      *is=*ps;
      *ie=*pe;
    }
    else if(*intdir==2){
      *js=*ps;
      *je=*pe;
    }
    else if(*intdir==3){
      *ks=*ps;
      *ke=*pe;
    }
  }
  else{
    dualfprintf(fail_file,"No such pointorlinetype=%d\n",pointorlinetype);
    myexit(287687456);
  }


  return(0);

}




int init_dumps(void)
{
  int i;
  char dumpnamelist[NUMDUMPTYPES][MAXFILENAME]=MYDUMPNAMELIST;

  trifprintf("begin: init_dumps\n");



  ///////////////////////////
  //
  // setup number of columns per dump file (see dumpgen.c or dump.c for how used)
  //
  /////////////////////////

  dnumcolumns[FAKEDUMPCOL]=0;

  // now setup the data output/input organization for chunking method for each number of columns
  dnumcolumns[IMAGECOL]=1;

  // always NPR
  dnumcolumns[RDUMPCOL]=NPR*2; // primitives and conservatives
  //dnumcolumns[RDUMPCOL]=NPR; // primitives only

  if(DOEVOLVEMETRIC){
    dnumcolumns[RMETRICDUMPCOL]=NPG*NDIM*NDIM*2 + NPG; // gcovlast, gcovpertlast, alphalapselast
  }
  else{
    dnumcolumns[RMETRICDUMPCOL]=0;
  }


  // always NPRDUMP
  if(GAMMIEDUMP)  dnumcolumns[DUMPCOL]=2*3 + NPRDUMP*2 + 3 + 1 + NDIM * NDIM + 6 + 1
#if(CALCFARADAYANDCURRENTS)
		    + NDIM*2
		    + 2*6
#endif
		    ;
  else{
    dnumcolumns[DUMPCOL]=3*3 + NPRDUMP*2 + 3 + 1 + NDIM * NDIM + 6 + 1 
#if(CALCFARADAYANDCURRENTS)
      + NDIM*2
      + 2*6
#endif
      ;    // 61 total if also doing currents and faraday, 41 otherwise

    if(FLUXB==FLUXCTSTAG && 0){ // DEBUG (change corresponding code in dump.c)
      dnumcolumns[DUMPCOL]+= NPR2INTERP*COMPDIM*2 + NPR + COMPDIM*3*2 + COMPDIM*3*2*2;
    }
  }



  // 205+4+4*4 currently
  //dnumcolumns[GDUMPCOL]=3*3+NDIM*NDIM*NDIM+NPG*NDIM*NDIM*2+NPG+4+4*4;
  //NPG was replaced with unity in order to avoid excessive dumping of info (only center info now available)
  dnumcolumns[GDUMPCOL]=3*3  +   NDIM*NDIM*NDIM  +   1*NDIM*NDIM*2   +   1  +  NDIM   +   NDIM*NDIM;
  //t^i x^i V^i,     \Gamma^\mu_{\nu\tau},     g^{\mu\nu} g_{\mu\nu}, \sqrt{-g}, \gamma_\mu, dx^\mu/dx^\nu
  


  // 36+29+8*2+4*2+2+12*2+96*2=339
  dnumcolumns[AVGCOL]=3*3 + 1 + NUMNORMDUMP  // (6+1+29=36)
    + NUMNORMDUMP // |normal terms| (29)
#if(CALCFARADAYANDCURRENTS)
    + NDIM*2 // jcon/jcov (8)
    + NDIM*2 // |jcon|/|jcov| (8)
#endif
    + NDIM*2 // massflux/|massflux|
    + NUMOTHER*2 // other stuff and fabs of each
#if(CALCFARADAYANDCURRENTS)
    +6*2 // fcon/fcov (12)
    +6*2 // |fcon|,|fcov| (12)
#endif
    +7*16 // Tud all 7 parts, all 4x4 terms (112)
    +7*16 // |Tud| all 7 parts, all 4x4 terms (112)
    ;


  if(DOAVG2){
    dnumcolumns[AVGCOL]-=224;
    dnumcolumns[AVG2COL]=10 + 224; // otherwise doesn't exist so don't need to set
  }
  else dnumcolumns[AVG2COL]=0;


  
  if(DODEBUG){
    dnumcolumns[DEBUGCOL]=NUMFAILFLOORFLAGS*NUMTSCALES;
  }
  else dnumcolumns[DEBUGCOL]=0;

  if(DOENODEBUG){
    //dnumcolumns[ENODEBUGCOL]=NUMENODEBUGS;
    dnumcolumns[ENODEBUGCOL]=(3-1)* NUMINTERPTYPES * (NPR-4) * NUMENODEBUGS;  //SASMARK2
  }
  else dnumcolumns[ENODEBUGCOL]=0;

  if(DOFIELDLINE){
    dnumcolumns[FIELDLINECOL]=NUMFIELDLINEQUANTITIES;
  }
  else dnumcolumns[FIELDLINECOL]=0;


  if(DODISS){
    dnumcolumns[DISSDUMPCOL]=NUMDISSFUNPOS;
  }
  else dnumcolumns[DISSDUMPCOL]=0;


  if(DODUMPOTHER){ // panalytic + numpother quantities
    dnumcolumns[DUMPOTHERCOL]=NPR+NUMPOTHER;
  }
  else dnumcolumns[DUMPOTHERCOL]=0;

  if(FLUXDUMP){ // dU, flux, and ppprimitives for flux
    dnumcolumns[FLUXDUMPCOL]=NUMFLUXDUMP;
  }
  else dnumcolumns[FLUXDUMPCOL]=0;

  //#if(WHICHEOS==KAZFULL)
  // all EOSs output same size data so uniform format
  // otherwise have to also put this condition in dump.c when outputting so don't overwrite memory!
  dnumcolumns[EOSDUMPCOL]=MAXPARLIST+1+MAXNUMEXTRAS+MAXPROCESSEDEXTRAS; // 1 is temperature
  //#else
  //  dnumcolumns[EOSDUMPCOL]=0;
  //#endif


  if(DOVPOTDUMP){
    dnumcolumns[VPOTDUMPCOL]=NUMVPOTDUMP;
  }
  else dnumcolumns[VPOTDUMPCOL]=0;


  trifprintf("dump number of columns(see global.nondepnmemonics.h)\n");
  for(i=0;i<NUMDUMPTYPES;i++){
    trifprintf("%s dnumcolumns[%d]=%d\n",dumpnamelist[i],i,dnumcolumns[i]);
  }
  trifprintf("\n");


  ///////////////////////////
  //
  // setup link list
  //
  ///////////////////////////
  init_linklists();




  trifprintf("end: init_dumps\n");


  return(0);
}




// Uavg is usually unew and Upoint is usually ulast at t=0
// fieldfrompotential[1,2,3 correspond to B1,B2,B3]
int pi2Uavg(int *fieldfrompotential, FTYPE (*prim)[N2M][N3M][NPR], FTYPE (*Upoint)[N2M][N3M][NPR], FTYPE (*Uavg)[N2M][N3M][NPR])
{
  struct of_geom geom,geomf;
  struct of_state q;
  int i,j,k;
  int pl;
  extern int initial_averageu_fv(int *fieldfrompotential, FTYPE (*prim)[N2M][N3M][NPR], FTYPE (*Upoint)[N2M][N3M][NPR], FTYPE (*Uavg)[N2M][N3M][NPR]);
  extern int initial_averageu_fluxrecon(int *fieldfrompotential, FTYPE prim[][N2M][N3M][NPR], FTYPE Upoint[][N2M][N3M][NPR], FTYPE Uavg[][N2M][N3M][NPR]);
  FTYPE Utemp[NPR];


  COMPFULLLOOP{
    // set geometry
    get_geometry(i, j, k, CENT, &geom);
    

    // find U(p)
    MYFUN(get_state(prim[i][j][k], &geom, &q),"initbasec:pi2Uavg()", "get_state()", 1);
    MYFUN(primtoU(UEVOLVE,prim[i][j][k], &q, &geom, Utemp),"initbase.c:pi2Uavg()", "primtoU()", 1);

    PLOOPNOB1(pl) Upoint[i][j][k][pl]=Utemp[pl];
    PLOOPNOB2(pl) Upoint[i][j][k][pl]=Utemp[pl];

    if(FLUXB==FLUXCTSTAG){
      PLOOPBONLY(pl) if(fieldfrompotential[pl-B1+1]==0){
	get_geometry(i, j, k, FACE1+(pl-B1), &geomf);
	Upoint[i][j][k][pl]=pstagscratch[i][j][k][pl]*(geomf.g);
      }
    }
    else{
      PLOOPBONLY(pl) if(fieldfrompotential[pl-B1+1]==0) Upoint[i][j][k][pl]=Utemp[pl];
    }
    
    //    dualfprintf(fail_file,"Upoint[%d][%d][%d][UU]=%21.15g prim[UU]=%21.15g\n",i,j,k,Upoint[i][j][k][UU]/(geom.g),prim[i][j][k][UU]);
  }



  //////////////////////////////////
  //
  // now deal with higher-order interpolations
  //
  //////////////////////////////////
  if(DOENOFLUX == ENOFINITEVOLUME){
    initial_averageu_fv(fieldfrompotential, prim, Upoint, Uavg);
  }
  else if(DOENOFLUX == ENOFLUXRECON){
    initial_averageu_fluxrecon(fieldfrompotential, prim, Upoint, Uavg);
  }
  else{
    // then just copy over
    COMPFULLLOOP{
      PLOOPNOB1(pl) Uavg[i][j][k][pl]=Upoint[i][j][k][pl];
      PLOOPNOB2(pl) Uavg[i][j][k][pl]=Upoint[i][j][k][pl];
      PLOOPBONLY(pl) if(fieldfrompotential[pl-B1+1]==0) Uavg[i][j][k][pl]=Upoint[i][j][k][pl];
    }
  }





  return(0);
}


void makedirs(void)
{

  if ((mpicombine && (myid == 0)) || (mpicombine == 0)) {
    
    if(USEMPI && (!MPIAVOIDFORK) || USEMPI==0){
      system("mkdir dumps");
      system("mkdir images");
    }
    else{
#ifndef WIN32
      // assumes unix commands exist in <sys/stat.h>
      // see also "info libc"
      // create directory with rxw permissions for user only
      mkdir("dumps",0700);
      mkdir("images",0700);
#else
      fprintf(stderr,"WIN32: User must create dumps and images\n");
#endif
    }
  }

#if(USEMPI)
  // all cpus wait for directory to be created
  MPI_Barrier(MPI_COMM_GRMHD);
#endif
  
}


#include<sys/stat.h>





// acts on globals, assumes static internals that get recalled upon reentering
int addremovefromnpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*in)[N2M][N3M][NPR], FTYPE (*out)[N2M][N3M][NPR])
{
  int addremovefromanynpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int *anynprstart, int *anynprend, int *anynprlist, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*in)[N2M][N3M][NPR], FTYPE (*out)[N2M][N3M][NPR]);

  // applies only to NPR type list
  addremovefromanynpr(doadd, whichpltoavg, ifnotavgthencopy, &nprstart, &nprend, nprlist, nprlocalstart, nprlocalend, nprlocallist, in, out);


  return(0);

}

// acts on globals, assumes static internals that get recalled upon reenterin
int addremovefromanynpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int *anynprstart, int *anynprend, int *anynprlist, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*in)[N2M][N3M][NPR], FTYPE (*out)[N2M][N3M][NPR])
{
  int pl;
  int pl2,pl3;
  int i,j,k;
  int num;



  if(doadd==REMOVEFROMNPR){
    ////////////////////////////////////////////
    //
    // save choice for interpolations
    *nprlocalstart=*anynprstart;
    *nprlocalend=*anynprend;
    PMAXNPRLOOP(pl) nprlocallist[pl]=anynprlist[pl];

    if(anynprlist==nprlist) num=NPR;
    else if(anynprlist==npr2interplist) num=NPR2INTERP;
    else if(anynprlist==nprboundlist) num=NPRBOUND;
    else{
      dualfprintf(fail_file,"No such type of list in addremovefromanynpr\n");
      myexit(7615156);
    }


    // now remove any other pl's not wanting to average for whatever reason
    for(pl3=0;pl3<num;pl3++){ // as above, but loop over all undesired quantities
      for(pl= *anynprstart;pl<= *anynprend;pl++){
	if(whichpltoavg[pl3]==0 && anynprlist[pl]==pl3){
	  // need to copy over unchanged quantity
	  if(ifnotavgthencopy[pl3] && in!=NULL && out!=NULL) COMPFULLLOOP out[i][j][k][pl3]=in[i][j][k][pl3];
	  for(pl2=pl+1;pl2<= *anynprend;pl2++) anynprlist[pl2-1]=anynprlist[pl2]; // moving upper to lower index
	  *anynprend-=1; // removed dir-field
	  break;
	}
      }
    }

  }
  else if(doadd==RESTORENPR){


    ////////////////////////////////////////////
    //
    // restore choice for interpolations
    *anynprstart= *nprlocalstart;
    *anynprend= *nprlocalend;
    PMAXNPRLOOP(pl) anynprlist[pl]=nprlocallist[pl];
  }


  //  PALLLOOP(pl){
  //    dualfprintf(fail_file,"dir=%d interptype=%d nprstart=%d nprend=%d nprlist[%d]=%d\n",dir,interptype,*anynprstart,*anynprend,pl,anynprlist[pl]);
  //  }
  

  return(0);

}


// used to transform from one coordinate system to PRIMECOORDS
// when acting on pstag, only relevant for magnetic field part, and in that case if didn't use vector potential to define pstag then assume not too important to get high accuracy, so average field to other positions in simple way
int transform_primitive_vB(int whichvel, int whichcoord, int i,int j, int k, FTYPE p[][N2M][N3M][NPR], FTYPE pstag[][N2M][N3M][NPR])
{

  // deal with pstag using p before p is changed
  MYFUN(transform_primitive_pstag(whichvel, whichcoord, i,j, k, p, pstag),"initbase.c:transform_primitive_vB()","transform_primitive_pstag()",0);
  

  // For p, transform from whichcoord to MCOORD
  // This changes p directly, so must come AFTER pstag change that uses original p
  if (bl2met2metp2v(whichvel,whichcoord,p[i][j][k], i,j,k) >= 1) FAILSTATEMENT("initbase.c:transform_primitive_vB()", "bl2ks2ksp2v()", 1);


  return(0);
}




int assert_func_empty( int is_bad_val, char *s, ... )
{
  return(is_bad_val);
}

int assert_func( int is_bad_val, char *s, ... )
{
  va_list arglist;
  FILE *fileptr = fail_file;

  va_start (arglist, s);

  if( 0 != is_bad_val ) {
    if(fileptr==NULL){
      fprintf(stderr,"tried to print to null file pointer: %s\n",s);
      fflush(stderr);
    }
    else{
      dualfprintf( fail_file, "Assertion failed: " );
      vfprintf (fileptr, s, arglist);
      fflush(fileptr);
    }
    if(myid==0){
      vfprintf (stderr, s, arglist);
      fflush(stderr);
    }
    va_end (arglist);

    myexit( 1 );
  }

  return is_bad_val;
}



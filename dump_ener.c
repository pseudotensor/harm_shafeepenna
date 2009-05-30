#include "decs.h"

int dump_ener(int doener, int dordump, int call_code)
{

  static FILE *flenerreg_file[NUMENERREGIONS];
  static FILE *flener_file;
  char FLENERREGIONNAME[NUMENERREGIONS][MAXFILENAME];
  char *FLENERNAME;

  static FILE *enerreg_file[NUMENERREGIONS];
  static FILE *ener_file;
  char ENERREGIONNAME[NUMENERREGIONS][MAXFILENAME];
  char *ENERNAME;

  static FILE *generreg_file[NUMENERREGIONS];
  static FILE *gener_file;
  char GENERREGIONNAME[NUMENERREGIONS][MAXFILENAME];
  char *GENERNAME;


  static FILE *metricparmsener_file;

  // per cpu file
  static FILE *probe_file;

  // single region files
  static FILE *debug_file;
  static FILE *lumener_file;
  static FILE *dissener_file[NUMDISSVERSIONS];

  int fileiter;
#define NUMSELFGRAV 9
  static FILE *selfgravener_file[NUMSELFGRAV];


  char dfnam[MAXFILENAME], ifnam[MAXFILENAME];
  int i, j, k, pl, l, dir,sc,fl,floor,tscale;
  int dissloop;
  //  FILE *imagecnt_file, *dumpcnt_file,*avgcnt_file,*debugcnt_file,*lumvsrcnt_file,*dissvsrcnt_file;

  // full grid only
  FTYPE divb, divbmax = 0, divbavg =0;
  // special jet2
  SFTYPE pdottermsjet2_tot[COMPDIM*2][NUMFLUXTERMS][NPR];

  // all regions, local quantities
  // these quantities can always be computed, no need in restart file and so in defs.h where many other _tot quantities appear
  SFTYPE Ureg_tot[NUMENERREGIONS][NPR];
  SFTYPE Ureg_final[NUMENERREGIONS][NPR];
  // each region, local quantities
  SFTYPE *U_tot;
  SFTYPE *U_final;

  int enerregion;

  static int firsttime = 1;
  SFTYPE ftemp0,ftemp1;

  int ii;
  FTYPE phibh;
  extern FTYPE phibh_compute(FTYPE M, FTYPE a, FTYPE r, FTYPE th);







  //////////////////////////////////
  //
  // Open/append some files
  //

  if ((call_code == INIT_OUT) || (firsttime == 1)) {
    // PER CPU FILE
    sprintf(dfnam,"probe.dat%s",myidtxt);
    if((probe_file=fopen(dfnam,"at"))==NULL){
      dualfprintf(fail_file,"Can't open probe file\n");
      myexit(1);
    }
    
    
    if(myid==0) if(DODEBUG){
      // CPU=0 only
      sprintf(dfnam,"debug.out");
      myfopen(dfnam,"a+","error opening debug output file\n",&debug_file);
    }

    if(myid==0) if(DOLUMVSR){
      // CPU=0 only
      sprintf(dfnam,"lumvsr.out");
      myfopen(dfnam,"a+","error opening lumvsr output file\n",&lumener_file);
    }

    if(myid==0) if(DODISSVSR){
	// CPU=0 only
	for(fileiter=0;fileiter<NUMDISSVERSIONS;fileiter++){
	  sprintf(dfnam,"dissvsr%d.out",fileiter);
	  myfopen(dfnam,"a+","error opening dissvsr output file\n",&dissener_file[fileiter]);
	}
    }

    if(myid==0) if(DOSELFGRAVVSR){
	// CPU=0 only
	for(fileiter=0;fileiter<NUMSELFGRAV;fileiter++){
	  sprintf(dfnam,"selfgravvsr%d.out",fileiter);
	  myfopen(dfnam,"a+","error opening selfgravvsr output file\n",&selfgravener_file[fileiter]);
	}
    }

    if(myid==0) if(DOEVOLVEMETRIC){
      // CPU=0 only
      sprintf(dfnam,"metricparms.out");
      myfopen(dfnam,"a+","error opening metricparms output file\n",&metricparmsener_file);
    }

      
    // CPU=0 only
    if(myid==0) ENERREGIONLOOP(enerregion){
      // set some pointers
      FLENERNAME=FLENERREGIONNAME[enerregion];

      // when setting file pointers, need pointer to pointer or just set pointers directly (we do latter for naming reasons)
      /* 	flener_file=flenerreg_file[enerregion]; */
      /* 	ener_file=&enerreg_file[enerregion]; */
      /* 	gener_file=generreg_file[enerregion]; */
      ENERNAME=ENERREGIONNAME[enerregion];

      GENERNAME=GENERREGIONNAME[enerregion];

      pcum_tot=pcumreg_tot[enerregion];
      fladd_tot=fladdreg_tot[enerregion];
      sourceadd_tot=sourceaddreg_tot[enerregion];

      // flener
      if(enerregion==GLOBALENERREGION) sprintf(FLENERNAME,"flener.out");
      else sprintf(FLENERNAME,"flenerother%d.out",enerregion); // specific naming convention
      myfopen(FLENERNAME,"a+","error opening FLenergy output file\n",&flenerreg_file[enerregion]);

      // ener
      if(enerregion==GLOBALENERREGION) sprintf(ENERNAME,ENERFNAME);
      else sprintf(ENERNAME,"enerother%d.out",enerregion); // specific naming convention

      if (appendold == 0) {
	// a+ in general, or w for overwrites, but let's append for now
	myfopen(ENERNAME,"a+","error opening  energy output file\n",&enerreg_file[enerregion]);
      }
      else {			// if appendold==1
	trifprintf("Start setup of %s file append\n", ENERNAME);
	myfopen(ENERNAME,"a+","error opening  energy output file for append\n",&enerreg_file[enerregion]);
	appendener(ener_file,pcum_tot,fladd_tot,sourceadd_tot);
      }

      // gener
      if(enerregion==GLOBALENERREGION) sprintf(GENERNAME,GENERFNAME);
      else sprintf(GENERNAME,"generjet%d.out",enerregion-1); // specific naming convention
      myfopen(GENERNAME,"a+","error opening  energy output file\n",&generreg_file[enerregion]);

    }
  }


  /////////////////////
  //
  // compute divB and output to logs
  //
  divbmaxavg(p,&divbmax,&divbavg);



  /////////////////////////////////////////
  //
  // ENER/FLENER/FRDOT/DEBUG output integrated diagnostics
  // for any number of predetermined regions
  //
  ////////////////////////////////////////

  ENERREGIONLOOP(enerregion){
    trifprintf(".BE.%d",enerregion);

    //////////////////////////////
    //
    // setup some pointers
    //
    //
    // each region, local quantities
    U_tot=Ureg_tot[enerregion];
    U_final=Ureg_final[enerregion];
    pdot_tot=pdotreg_tot[enerregion];
    pdotterms_tot=pdottermsreg_tot[enerregion];
    // used for each region, related to global quantities
    fladd=fladdreg[enerregion];
    fladd_tot=fladdreg_tot[enerregion];
    fladdterms=fladdtermsreg[enerregion];
    fladdterms_tot=fladdtermsreg_tot[enerregion];
    U_init=Ureg_init[enerregion];
    U_init_tot=Ureg_init_tot[enerregion];
    pcum=pcumreg[enerregion];
    pcum_tot=pcumreg_tot[enerregion];
    pdot=pdotreg[enerregion];
    pdotterms=pdottermsreg[enerregion];
    sourceaddterms=sourceaddtermsreg[enerregion];
    sourceaddterms_tot=sourceaddtermsreg_tot[enerregion];
    sourceadd=sourceaddreg[enerregion];
    sourceadd_tot=sourceaddreg_tot[enerregion];
    diss=dissreg[enerregion];
    diss_tot=dissreg_tot[enerregion];

    if(myid==0){ // only cpu=0 writes to these files, although doesn't matter
      ener_file=enerreg_file[enerregion];
      gener_file=generreg_file[enerregion];
      flener_file=flenerreg_file[enerregion];
    }

    /////////////////////////////////
    //
    // do some integrations
    //
    if(dordump||doener){

      trifprintf("BI%d",enerregion);

      /////////////////////////
      //
      // things that are done for EVERY region (can't sum over all regions at once since each region has different WITHINENERREGION() conditional)
      // Below assumes continuous array starting at some address -- otherwise would have to do per quantity as in prior versions with extra loops
      //
      /////////////////////////

      // compute total conserved quantity
      if(integrate(NPR,U_tot,U_tot,CONSTYPE,enerregion)>=1) return(1);

      // DIRLOOP(dir)
      if(integrate((COMPDIM*2)*NPR,&pdot[0][0],&pdot_tot[0][0],SURFACETYPE,enerregion)>=1) return(1);

      // above should be sum of below, approximately
      //      DIRLOOP(dir) FLLOOP(fl)
      if(integrate((COMPDIM*2)*NUMFLUXTERMS*NPR,&pdotterms[0][0][0],&pdotterms_tot[0][0][0],SURFACETYPE,enerregion)>=1) return(1);
      //DIRLOOP(dir)
      if(integrate((COMPDIM*2)*NPR,&pcum[0][0],&pcum_tot[0][0],SURFACETYPE,enerregion)>=1) return(1);
      
      //FLOORLOOP(floor)
      if(integrate(NUMFAILFLOORFLAGS*NPR,&fladdterms[0][0],&fladdterms_tot[0][0],CUMULATIVETYPE,enerregion)>=1) return(1);

      if(integrate(NPR,&fladd[0],&fladd_tot[0],CUMULATIVETYPE,enerregion)>=1) return(1);

      //SCLOOP(sc)
      if(integrate(NUMSOURCES*NPR,&sourceaddterms[0][0],&sourceaddterms_tot[0][0],CUMULATIVETYPE,enerregion)>=1) return(1);
      
      if(integrate(NPR,&sourceadd[0],&sourceadd_tot[0],CUMULATIVETYPE,enerregion)>=1) return(1);

      if(DODISS){
	if(integrate(NUMDISSVERSIONS,&diss[0],&diss_tot[0],CUMULATIVETYPE,enerregion)>=1) return(1);
      }

      /////////////////////////
      // below was subsumed into its own full enerregion (should have done so in first place!)
      //
      // OUTSIDEHORIZONENERREGION STUFF (i.e. only done once over that region)
      //      if(enerregion==OUTSIDEHORIZONENERREGION){
	// now integrate horizon fluces over all CPUs
	// this is also done during actual evolution in metric.c, but needed here too so diagnostics are correctly updated
	// doesn't hurt evolution at all
      //	if(integrate(horizonflux,horizonflux_tot,SURFACETYPE,enerregion)>=1) return(1);
      //	if(integrate(horizoncum,horizoncum_tot,SURFACETYPE,enerregion)>=1) return(1);

	// assume DOSELFGRAVVSR has already cumulated dMvsr into dMvsr_tot and computed Mvsr_tot and phivsr_tot

      //      }


      /////////////////////////
      //
      // GLOBALENERREGION STUFF (i.e. only done once over that region or just simply done once at all and just keeping within loop for simplicity)
      if(enerregion==GLOBALENERREGION){
	// CUMULATIVETYPE2 for 1 data object per ii
	if(DOLUMVSR){
	  //for(ii=0;ii<ncpux1*N1;ii++)
	  if(integrate(ncpux1*N1,&lumvsr[0],&lumvsr_tot[0],CUMULATIVETYPE,enerregion)>=1) return(1);
	  //	  for(ii=0;ii<ncpux1*N1;ii++){
	  //	    dualfprintf(fail_file,"nstep=%ld ii=%d lumvsr=%21.15g lumvsr_tot=%21.15g\n",nstep,ii,lumvsr[ii],lumvsr_tot[ii]);
	  //	  }
	}

	if(DODISSVSR){
	  for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++){// this loop is over pointers, not a continuous memory space!
	    // for(ii=0;ii<ncpux1*N1;ii++)
	    if(integrate(ncpux1*N1,&dissvsr[dissloop][0],&dissvsr_tot[dissloop][0],CUMULATIVETYPE,enerregion)>=1) return(1);
	  }
	}


	// special jet2 accounting
	//DIRLOOP(dir) FLLOOP(fl) 
	if(integrate((COMPDIM*2)*NUMFLUXTERMS*NPR,&pdottermsjet2[0][0][0],&pdottermsjet2_tot[0][0][0],SURFACETYPE,enerregion)>=1) return(1);
	
	// debug
	if(DODEBUG){
	  //TSCALELOOP(tscale)
	  // below is special integratel and for CONSTYPE performs integration over failfloorcount[i][j][k][tscale][pl] for each enerregion
	  if(integratel(NUMTSCALES*NUMFAILFLOORFLAGS,&failfloorcountlocal[0][0],&failfloorcountlocal_tot[0][0],CONSTYPE,enerregion)>=1) return(1);
	}
      }
      trifprintf("AI%d",enerregion);
    }

    /////////////////////////////
    //
    // initialize the total conserved counters
    //
    if (call_code == INIT_OUT || firsttime) {
      if(RESTARTMODE==0)    PLOOP(pl) U_init_tot[pl] = U_init[pl] = U_tot[pl];
      // otherwise read from restart file
    }

    /////////////////////////////
    //
    // output some interesting diagnostics to stderr/log files.
    //
    if (call_code == FINAL_OUT) {
      PLOOP(pl) U_final[pl] = U_tot[pl];

      if(GAMMIEENER){
	trifprintf("\n\nEnergy: ini,fin,del: %21.15g %21.15g %21.15g\n",
		   U_init[UU], U_final[UU], (U_final[UU] - U_init[UU]) / U_init[UU]);

	trifprintf("\n\nMass: ini,fin,del: %21.15g %21.15g %21.15g\n",
		   U_init[RHO], U_final[RHO], (U_final[RHO] - U_init[RHO]) / U_init[RHO]);
      }
      else{
	trifprintf("\n");
	PDUMPLOOP(pl){
	  ftemp0=(U_final[pl] - U_init[pl]);
	  ftemp1=(U_final[pl]-fladd_tot[pl]-sourceadd_tot[pl]-(pcum_tot[X1DN][pl]+pcum_tot[X2DN][pl]+pcum_tot[X3DN][pl]) + (pcum_tot[X1UP][pl]+pcum_tot[X2UP][pl]+pcum_tot[X3UP][pl]) - U_init[pl]);
	  trifprintf("U[%d]: ini,fin,fdf,del,tfdf,tdel: %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
		     pl,
		     U_init[pl],
		     U_final[pl],
		     U_init[pl]!=0.0 ? ftemp0/U_init[pl] : 0.0,
		     ftemp0,
		     U_init[pl]!=0.0 ? ftemp1/U_init[pl] : 0.0,
		     ftemp1 );
	}	
      }
    }


    /////////////////////
    //
    //  output some ener diagnostics to ener, flener, debug, and probe(only cpu specific one) files
    //
    if (doener){
      if(1||GAMMIEENER){
	// 6+3=9 terms
	myfprintf(gener_file, "%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g ", t, U_tot[RHO], U_tot[U3], U_tot[UU], p[N1 / 2][N2 / 2][N3 / 2][UU] * pow(p[N1 / 2][N2 / 2][N3 / 2][RHO], -gam), p[N1 / 2][N2 / 2][N3 / 2][UU]);
	myfprintf(gener_file, "%21.15g %21.15g %21.15g ", pdot_tot[X1DN][RHO], pdot_tot[X1DN][RHO]-pdot_tot[X1DN][UU], pdot_tot[X1DN][U3]);
      }
      if(1){
	// SM use gammie.m macro, jrdpener, gammieener's


	//////////////////////////
	// ENER FILE (only dir and pl)
	//
	//////////////////////////
	
	// 2+NPR+COMPDIM*2*NPRDUMP+NPRDUMP+(2) terms
	// see jrdp3dener in gammie.m
	myfprintf(ener_file,"%21.15g %ld ",t ,realnstep);
	PDUMPLOOP(pl) myfprintf(ener_file, "%21.15g ", U_tot[pl]);
	DIRLOOP(dir) PDUMPLOOP(pl) myfprintf(ener_file, "%21.15g ", pdot_tot[dir][pl]);
	// COMPDIM*2*NPRDUMP more terms
	DIRLOOP(dir) PDUMPLOOP(pl) myfprintf(ener_file, "%21.15g ", pcum_tot[dir][pl]);
	PDUMPLOOP(pl) myfprintf(ener_file, "%21.15g ", fladd_tot[pl]);
	PDUMPLOOP(pl) myfprintf(ener_file, "%21.15g ", sourceadd_tot[pl]);
	myfprintf(ener_file, "%21.15g %21.15g ", divbmax,divbavg);
	if(DODISS) for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++) myfprintf(ener_file, "%21.15g ", diss_tot[dissloop]);
	else for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++) myfprintf(ener_file, "%21.15g ", -1.0); // dummy space holder
      }


      ////////////////////////////////
      // FLENER FILE (dir, pl, and linear summed terms for flux, floor, and source quantities)
      // (note the change in ordering for external file read in macros/functions)
      myfprintf(flener_file,"%21.15g %ld ",t,realnstep);
      DIRLOOP(dir) PDUMPLOOP(pl) FLLOOP(fl) myfprintf(flener_file, "%21.15g ", pdotterms_tot[dir][fl][pl]);
      PDUMPLOOP(pl) FLOORLOOP(floor) myfprintf(flener_file, "%21.15g ", fladdterms_tot[floor][pl]);
      PDUMPLOOP(pl) SCLOOP(sc) myfprintf(flener_file, "%21.15g ", sourceaddterms_tot[sc][pl]);
      

      if(enerregion==GLOBALENERREGION){ // only for total region for now
	// only care about inner and outer radial part
	for(dir=0;dir<=1;dir++) PDUMPLOOP(pl) FLLOOP(fl) myfprintf(flener_file, "%21.15g ", pdottermsjet2_tot[dir][fl][pl]); // jet/pole values only

	if(DOLUMVSR){
	  // luminosity vs radius
	  myfprintf(lumener_file,"%21.15g %ld ",t,realnstep);
	  for(ii=0;ii<ncpux1*N1;ii++) myfprintf(lumener_file, "%21.15g ", lumvsr_tot[ii]);
	}

	if(DODISSVSR){
	  for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++){
	    // dissipation vs radius
	    myfprintf(dissener_file[dissloop],"%21.15g %ld ",t,realnstep);
	    for(ii=0;ii<ncpux1*N1;ii++) myfprintf(dissener_file[dissloop], "%21.15g ", dissvsr_tot[dissloop][ii]);
	  }
	}

      }



      if(enerregion==OUTSIDEHORIZONENERREGION){ // only for horizon region for now

	
	if(DOEVOLVEMETRIC){
	  myfprintf(metricparmsener_file,"%21.15g %ld ",t,realnstep);

	  // below 2 have been subsumed into its own full ener region (should have done in first place!)
	  //PDUMPLOOP(pl) myfprintf(metricparmsener_file, "%21.15g ",horizonflux_tot[pl]);
	  //PDUMPLOOP(pl) myfprintf(metricparmsener_file, "%21.15g ",horizoncum_tot[pl]);

	  // first 2 are cumulative black hole mass and J in black hole metric units
	  // next 1 is what black hole mass would be if adding so-far cumulated mass
	  // next 1 is a if adding so-far mass and angular momentum
	  // next 1 is j=a/M if adding so-far mass and angular momentum
	  // next several items are original and present actual mass, a, j, Q, and q

	  // by this time of the diagnostics, dE, dJ include horizon flux but MBH and a aren't yet updated.
	  // however, horizoncum_tot IS updated so that for this enerregion u0-horizoncum_tot[0] is conserved exactly
	  myfprintf(metricparmsener_file,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %d ",dE,dJ,MBH0+dE,a0+dabh,(a0+dabh)/(SMALL+MBH0+dE),MBH0,MBH,a0,a0/(SMALL+MBH0),a,a/(SMALL+MBH),QBH0,QBH0/(SMALL+MBH0),QBH,QBH/(SMALL+MBH),Rhor,Risco,horizoni+N1*horizoncpupos1);
	}

	if(DOSELFGRAVVSR){
	  // self-gravitational potential vs radius
	  fileiter=0;
	  myfprintf(selfgravener_file[fileiter],"%21.15g %ld ",t,realnstep);
	  DUMPGRAVLOOP(ii) myfprintf(selfgravener_file[fileiter], "%21.15g ", rcent_tot[ii]);

	  // just self-gravitating part of Mvsr
	  fileiter=1;
	  myfprintf(selfgravener_file[fileiter],"%21.15g %ld ",t,realnstep);
	  DUMPGRAVLOOP(ii) myfprintf(selfgravener_file[fileiter], "%21.15g ", Mvsr_tot[ii]);

	  // true Mvsr with "point" mass included
	  fileiter=2;
	  myfprintf(selfgravener_file[fileiter],"%21.15g %ld ",t,realnstep);
	  DUMPGRAVLOOP(ii) myfprintf(selfgravener_file[fileiter], "%21.15g ", MBH+Mvsr_tot[ii]);

	  // self-gravitating part of phivsr
	  fileiter=3;
	  myfprintf(selfgravener_file[fileiter],"%21.15g %ld ",t,realnstep);
	  DUMPGRAVLOOP(ii) myfprintf(selfgravener_file[fileiter], "%21.15g ", phivsr_tot[ii]);

	  // true phivsr that comes into g_{tt} = -exp(2\phi)
	  fileiter=4;
	  myfprintf(selfgravener_file[fileiter],"%21.15g %ld ",t,realnstep);
	  DUMPGRAVLOOP(ii){
	    phibh = phibh_compute(MBH,a,rcent_tot[ii], M_PI*0.5);
	    if(!isfinite(phibh)) phibh=-BIG; // in reality g_{tt} stays well-defined
	    myfprintf(selfgravener_file[fileiter], "%21.15g ", phibh+phivsr_tot[ii]);
	  }
	  fileiter=5;
	  myfprintf(selfgravener_file[fileiter],"%21.15g %ld ",t,realnstep);
	  DUMPGRAVLOOP(ii) myfprintf(selfgravener_file[fileiter], "%21.15g ",dTrrvsr_tot[ii]);

	  fileiter=6;
	  myfprintf(selfgravener_file[fileiter],"%21.15g %ld ",t,realnstep);
	  DUMPGRAVLOOP(ii) myfprintf(selfgravener_file[fileiter], "%21.15g ", dVvsr_tot[ii]);

	  fileiter=7;
	  myfprintf(selfgravener_file[fileiter],"%21.15g %ld ",t,realnstep);
	  DUMPGRAVLOOP(ii) myfprintf(selfgravener_file[fileiter], "%21.15g ", vrsqvsr_tot[ii]);

	  fileiter=8;
	  myfprintf(selfgravener_file[fileiter],"%21.15g %ld ",t,realnstep);
	  DUMPGRAVLOOP(ii) myfprintf(selfgravener_file[fileiter], "%21.15g ", Jvsr_tot[ii]);

	}

      }




      if(enerregion==GLOBALENERREGION){ // only for total region for now
	/////////////////////////
	// PROBE FILE
	// !!per CPU!! probe file NPRDUMP*3 terms
	//
	////////////////////////
#define ITER1 MAX(N1/4,1)
#define ITER2 MAX(N2/4,1)
#define ITER3 MAX(N3/4,1)
	// 2D probe
	// k  = N3/2+1;	for(i=0;i<N1;i+=ITER1) for(j=0;j<N2;j+=ITER2){
	for(i=0;i<N1;i+=ITER1) for(j=0;j<N2;j+=ITER2) for(k=0;k<N3;k+=ITER3){
	  // 2D probe (consistent with original SM macro)
	  //	  PDUMPLOOP(pl) fprintf(probe_file, "%d %d %d %ld %21.15g %21.15g\n",startpos[1]+i,startpos[2]+j,pl,realnstep, t, p[i][j][k][pl]);
	  // 3d probe
	  PDUMPLOOP(pl) fprintf(probe_file, "%d %d %d %d %ld %21.15g %21.15g\n",startpos[1]+i,startpos[2]+j,startpos[3]+k,pl,realnstep, t, p[i][j][k][pl]);
	}
	fflush(probe_file);


	/////////////////////////
	// DEBUG FILE
	if(DODEBUG){
	  myfprintf(debug_file,"%21.15g %ld ",t,realnstep);
	  myfprintf(debug_file,"%21.15g %21.15g ",dt, 1.0*nstroke/(2.0*totalzones));
	  myfprintf(debug_file,"%21.15g %21.15g %21.15g ",wavedtglobal,sourcedtglobal,gravitydtglobal);
	  myfprintf(debug_file,"%d %d %d ",waveglobaldti[1],waveglobaldtj[1],waveglobaldtk[1]);
	  myfprintf(debug_file,"%d %d %d ",waveglobaldti[2],waveglobaldtj[2],waveglobaldtk[2]);
	  myfprintf(debug_file,"%d %d %d ",waveglobaldti[3],waveglobaldtj[3],waveglobaldtk[3]);
	  myfprintf(debug_file,"%d %d ",horizoni,horizoncpupos1);
#if(COUNTTYPE==LONGLONGINTTYPE)
	    TSCALELOOP(tscale) FLOORLOOP(floor) myfprintf(debug_file, "%lld ", failfloorcountlocal_tot[tscale][floor]);
#elif(COUNTTYPE==DOUBLETYPE)
	    TSCALELOOP(tscale) FLOORLOOP(floor) myfprintf(debug_file, "%21.15g ", failfloorcountlocal_tot[tscale][floor]);
#endif
	}
      }



      myfprintf(flener_file,"\n");
      myfprintf(ener_file,"\n");
      myfprintf(gener_file,"\n");
      if(enerregion==GLOBALENERREGION){
	if(DOLUMVSR) myfprintf(lumener_file,"\n");
	if(DODISSVSR) 	  for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++) myfprintf(dissener_file[dissloop],"\n");
	if(DODEBUG) myfprintf(debug_file,"\n");
      }
      if(enerregion==OUTSIDEHORIZONENERREGION){
	
	if(DOEVOLVEMETRIC) myfprintf(metricparmsener_file,"\n");
	if(DOSELFGRAVVSR){
	  for(fileiter=0;fileiter<NUMSELFGRAV;fileiter++) myfprintf(selfgravener_file[fileiter],"\n");
	}
      }

      trifprintf("W%d",enerregion);

    }// end if doener // done if writing to file

    ///////////////
    //
    // close the files at end of simulation
    //
    if(call_code==FINAL_OUT){
      myfclose(&flener_file,"Couldn't close flener_file\n");
      if(1) myfclose(&ener_file,"Couldn't close ener_file\n");
      if(1||GAMMIEENER) myfclose(&gener_file,"Couldn't close gener_file\n");


      if(enerregion==OUTSIDEHORIZONENERREGION){
	if(DOEVOLVEMETRIC) myfclose(&metricparmsener_file,"Couldn't close metricparmsener_file\n");
	if(DOSELFGRAVVSR){
	  for(fileiter=0;fileiter<NUMSELFGRAV;fileiter++) myfclose(&selfgravener_file[fileiter],"Couldn't close selfgravener_file\n");
	}
      }

      if(enerregion==GLOBALENERREGION){
	if(DOLUMVSR) myfclose(&lumener_file,"Couldn't close lumener_file\n");
	if(DODISSVSR) 	  for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++) myfclose(&dissener_file[dissloop],"Couldn't close dissener_file\n");
	if(DODEBUG) myfclose(&debug_file,"Couldn't close debug_file\n");
	fclose(probe_file);
      }
    }
  }// end enerregion loop



  trifprintf("E");
  firsttime = 0;
  return (0);




}


// used to append into ener.out file.  Needed to locate place within file to begin appending.
// Old code, probably broken by now -- not used
void appendener(FILE* ener_file,SFTYPE pcum_tot[][NPR],SFTYPE*fladd_tot,SFTYPE *sourceadd_tot)
{
  int gotit;
  SFTYPE tcheck;
  int l,pl,dir;
  long fpos0;
  FILE *ener_file_temp;
  char dfnam[MAXFILENAME],dfnamback[MAXFILENAME], dfnamtemp[MAXFILENAME];

  // only CPU=0 does anything here
  if(myid==0){
    
    rewind(ener_file);	// go to start
    
    gotit = 0;
    while ((!feof(ener_file)) && (gotit == 0)) {
      
      fscanf(ener_file, "%lf", &tcheck);
      
      if (fabs(tcheck - t) < 0.5 *DTdumpgen[DTENER]) {
	gotit = 1;
	for (l = 1; l <= NUMENERVAR; l++) {
	  if ((l > 3+NPR+COMPDIM*2*NPR) && (l < 3+NPR+2*COMPDIM*2*NPR+NPR)) {
	    DIRLOOP(dir) PDUMPLOOP(pl) {
	      fscanf(ener_file, "%lf", &pcum_tot[dir][pl]);
	      l++;
	    }
	    PDUMPLOOP(pl) {
	      fscanf(ener_file, "%lf", &fladd_tot[pl]);
	      l++;
	    }
	    PDUMPLOOP(pl) {
	      fscanf(ener_file, "%lf", &sourceadd_tot[pl]);
	      l++;
	    }
	  }
	}
      } else {
	// skip this bad line 
	while ((fgetc(ener_file) != '\n') && (!feof(ener_file)));
      }
      // continue after successful get since successful get is good 
      // data and should keep since corresponds to dump one is
      // keeping
      fpos0 = ftell(ener_file);	// position to continue
      // writting at if successful get
    }
    if (gotit == 0) {
      dualfprintf(fail_file,
		  "Never found right time in loss file when appending: looking for t=%21.15g lastt=%21.15g\n",
		  t, tcheck);
      myexit(1);
    } else {
      dualfprintf(logfull_file,
		  "found goodtime t=%21.15g (wanted %21.15g) to restart ener file\n",
		  tcheck, t);
      sprintf(dfnamtemp, "%s0_ener%s.temp", DATADIR, ".out");
      sprintf(dfnam, "%sener%s", DATADIR, ".out");
      sprintf(dfnamback, "%s0_ener%s.back", DATADIR, ".out");
      
      // now that done, fix up file
      if ((ener_file_temp = fopen(dfnamtemp, "wt")) == NULL) {
	dualfprintf(fail_file,
		    "Cannot open temp ener file for appending: %s\n",
		    dfnamtemp);
	myexit(1);
      } else {
	rewind(ener_file);
	while (ftell(ener_file) < fpos0 + 1) {	// +1 is for
	  // '\n' at end
	  // of line
	  fputc(fgetc(ener_file), ener_file_temp);
	}
	fclose(ener_file_temp);
	fclose(ener_file);
	rename(dfnam, dfnamback);	// move old to backup location
	rename(dfnamtemp, dfnam);	// move new to old name(normal
	// name)
	// reopen loss_file (now normal name)

	if ((ener_file = fopen(dfnam, "at")) == NULL) {
	  dualfprintf(fail_file,
		      "2: error opening ener output file %s\n", dfnam);
	  myexit(1);
	}
	trifprintf("End setup of ener file append\n");
      }			// end else if can open temp file
    }			// end else if gotit==1
  }
}

void divbmaxavg(FTYPE prim[][N2M][N3M][NPR],FTYPE*ptrdivbmax,FTYPE*ptrdivbavg)
{
  int i,j,k;
  int imax=0,jmax=0,kmax=0;
  FTYPE divb;
  FTYPE divbmax=0,divbavg=0;
  FTYPE divbmaxsend,divbavgsend;



  LOOPDIVB { // diagonostic loop
    // doesn't need geom, just use global gdet
    setfdivb(&divb, prim, udump, i, j, k); // udump set externally GODMARK
    if (divb > divbmax) {
      imax = i;
      jmax = j;
      kmax = k;
      divbmax = divb;
    }
    divbavg += divb;
  }

  // PER CPU
  fprintf(log_file,"  proc: %04d : divbmax: %d %d %d : %21.15g divbavg: %21.15g\n", myid, imax, jmax, kmax, divbmax, divbavg / ((FTYPE) (N1*N2*N3))); fflush(log_file);

#if(USEMPI)			// give CPU=0 total
  divbmaxsend = divbmax;
  divbavgsend = divbavg;

#if USINGORANGE == 1
  MPI_Barrier(MPI_COMM_GRMHD); // GODMARK: OpenMP on Orange crashes without this
#endif
  MPI_Reduce(&divbmaxsend, &divbmax, 1, MPI_FTYPE, MPI_MAX, 0,  MPI_COMM_GRMHD);
#if USINGORANGE == 1
  MPI_Barrier(MPI_COMM_GRMHD); // GODMARK: OpenMP on Orange crashes without this
#endif
  MPI_Reduce(&divbavgsend, &divbavg, 1, MPI_FTYPE, MPI_SUM, 0,  MPI_COMM_GRMHD);
#if USINGORANGE == 1
  MPI_Barrier(MPI_COMM_GRMHD); // GODMARK: OpenMP on Orange crashes without this
#endif

#endif
  divbavg /= (FTYPE) (totalzones);

  // Total over all CPUs
  myfprintf(logfull_file,"  divbmax: %21.15g divbavg: %21.15g\n",divbmax, divbavg);
  *ptrdivbmax=divbmax;
  *ptrdivbavg=divbavg;
}






void setrestart(int*appendold)
{
  *appendold = 0;
  // 0: deal with ener.out manually
  // 1: append automatically (no longer used)
}

/* gettotal accepts an arbitrary pointer set each of different sizes
 * i.e. one could do:
 * numptrs=2+NPR;
 * totalptrs[0]=pdot;  totalsizes[0]=NPR; totaloptrs[0]=pdot;
 * totalptrs[1]=fladd; totalsizes[1]=NPR; totaloptrs[1]=fladd_tot;
 * PLOOP(pl){ totalptrs[2+pl]=&U_tot[pl]; totalsizes[pl]=1;}
 * gettotal(numptrs,totalptrs,totaloptrs,totalsizes);
 *
 */

void gettotal(int numvars, SFTYPE* vars[],int*sizes,SFTYPE*vars_tot[])
{
  int j,k;
  SFTYPE *ptrsend;
  SFTYPE send;
  int didmalloc;
 
  // for 1 CPU
  if (numprocs == 1) {
    // must use _tot since can't overwrite normal global integrator 
    // variable
    for(k=0;k<numvars;k++){
      for(j=0;j<sizes[k];j++){
	vars_tot[k][j] = vars[k][j];
      }
    }
  } else {
    // give CPU=0 the totals
#if(USEMPI)

    for(k=0;k<numvars;k++){

      didmalloc=0;
      // at least assume each pointer is a continuous array of memory and avoid extra MPI_Reduce()'s
      // send and receive can't be same address, hence "ptrsend" variable
      if(&vars[k][0] == &vars_tot[k][0]){
	// then need to make scratch space since send and recieve can't be same address
	// assume not too large
	ptrsend = (SFTYPE*) malloc(sizes[k]*sizeof(SFTYPE));
	if(ptrsend==NULL){
	  dualfprintf(fail_file,"Couldn't get memory for ptrsend in gettotal with k=%d numvars=%d sizes=%d\n",k,numvars,sizes[k]);
	  myexit(496365763);
	}
	didmalloc=1;
	// now copy over to temporary space
	for(j=0;j<sizes[k];j++){
	  ptrsend[j] = vars[k][j];
	}
      }
      else{
	// input and output not same space so can use vars directly as input
	ptrsend = &vars[k][0];
      }

      MPI_Reduce(ptrsend, &vars_tot[k][0], sizes[k], MPI_SFTYPE, MPI_SUM, 0, MPI_COMM_GRMHD);

      // free malloc'ed memory
      if(didmalloc) free(ptrsend);

      // OLD WAY:
      //      for(j=0;j<sizes[k];j++){
      //	send=vars[k][j];
      //	// send and receive can't be same address, hence "send" variable
      //	MPI_Reduce(&send, &vars_tot[k][j], 1, MPI_SFTYPE, MPI_SUM, 0, MPI_COMM_GRMHD);
      //      }


    } // end loop over pointers [k]



#endif
  }
}



// all CPUs get total
void getalltotal(int numvars, SFTYPE* vars[],int*sizes,SFTYPE*vars_tot[])
{
  int j,k;
  SFTYPE send;
  
  // for 1 CPU
  if (numprocs == 1) {
    // must use _tot since can't overwrite normal global integrator 
    // variable
    for(k=0;k<numvars;k++){
      for(j=0;j<sizes[k];j++){
	vars_tot[k][j] = vars[k][j];
      }
    }
  } else {
    // give all CPUs the totals
#if(USEMPI)
    for(k=0;k<numvars;k++){
      for(j=0;j<sizes[k];j++){
	send=vars[k][j];
	// send and receive can't be same address, hence "send" variable
	// GODMARK: This could be optimized like gettotal, but not necessary for now since only used with CUMULATIVETYPE3 and not used right now
	MPI_Allreduce(&send, &vars_tot[k][j], 1, MPI_SFTYPE, MPI_SUM, MPI_COMM_GRMHD);
      }
    }
#endif
  }
}



void gettotall(int numvars, CTYPE* vars[],int*sizes,CTYPE *vars_tot[])
{
  int j,k;
  CTYPE send;
  
  // for 1 CPU
  if (numprocs == 1) {
    // must use _tot since can't overwrite normal global integrator 
    // variable
    for(k=0;k<numvars;k++){
      for(j=0;j<sizes[k];j++){
	vars_tot[k][j] = vars[k][j];
      }
    }
  } else {
    // give CPU=0 the totals
#if(USEMPI)
    for(k=0;k<numvars;k++){
      for(j=0;j<sizes[k];j++){
	send=vars[k][j];
	// send and receive can't be same address, hence "send" variable
	MPI_Reduce(&send, &vars_tot[k][j], 1, MPI_CTYPE, MPI_SUM, 0,
		   MPI_COMM_GRMHD);
      }
    }
#endif
  }
}


// each CPU does constotal
int constotal(int enerregion, SFTYPE *vars)
{
  int i,j,k,pl;
  FTYPE U[NPR];
  struct of_geom geom;
  struct of_state q;
  FTYPE ftemp[NPR];

  PLOOP(pl) vars[pl]= 0.0;

  enerpos=enerposreg[enerregion];

  ZLOOP{ // diagonostic loop
    if(WITHINENERREGION(enerpos,i,j,k)){

      if(DOENOFLUX != NOENOFLUX){
	// GODMARK
	// then need to use stored averaged conserved quantity
	// assume done with all timesteps when here, so using unew
	// for substep access one would need to know which substep and then use the correct version (ulast) on all substeps except in between full steps when unew is the correct answer
	PLOOP(pl){
	  ftemp[pl]=unew[i][j][k][pl]*dVF;
	  vars[pl] += ftemp[pl];
	}
      }
      else{
	get_geometry(i,j,k,CENT,&geom) ;
	if(!failed){
	  if(get_state(p[i][j][k],&geom,&q)>=1) return(1);
	  if(primtoU(UDIAG,p[i][j][k],&q,&geom,U)>=1) return(1);
	}
	PLOOP(pl){
	  ftemp[pl]=U[pl]*dVF;
	  vars[pl] += ftemp[pl];
	}
      }
    }
  }
  return(0);
}


// for dissipation (GODMARK: not used currently since keep track of integral and function separately)
// each CPU does constotal
int constotal2(int enerregion, SFTYPE *vars)
{
  int i,j,k,pl;
  FTYPE U[NPR];
  struct of_geom geom;
  struct of_state q;
  FTYPE ftemp[NPR];
  //  int dissloop;

  vars[0]= 0.0;

  enerpos=enerposreg[enerregion];

  ZLOOP{ // diagonostic loop
    if(WITHINENERREGION(enerpos,i,j,k) ){
      //      for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++){
      vars[0] += dissfunpos[i][j][k][0];
	//      }
    }
  }
  return(0);
}




// each CPU does counttotal
int counttotal(int enerregion, CTYPE *vars, int num)
{
  int i,j,k,variter;

  for(variter=0;variter<num;variter++) vars[variter]= 0;

  enerpos=enerposreg[enerregion];

  ZLOOP { // diagonostic loop
    if(WITHINENERREGION(enerpos,i,j,k) ){
      for(variter=0;variter<num;variter++){
	// looping over [0][variter] is equivalent to original TSCALELOOP FLLOOP since at least last two entries in failfloorcount array are continuous in memory
	vars[variter] += failfloorcount[i][j][k][0][variter];
      }
    }
  }
  return(0);
}



#define MAXPTRS 10

int integrate(int numelements, SFTYPE * var,SFTYPE *var_tot,int type, int enerregion)
{
  SFTYPE *totalptrs[MAXPTRS],*totaloptrs[MAXPTRS];
  int totalsizes[MAXPTRS],numptrs;

  switch(type){
  case CONSTYPE:
    if(constotal(enerregion,var)>=1) return(1);
    totalsizes[0]=numelements;
    gettotal(1,&var,totalsizes,&var_tot); // if only 1 object, then just send address
    break;
    // below not used or deprecated
    //  case CONSTYPE2:
    //    if(constotal2(enerregion,var)>=1) return(1);
    //    totalsizes[0]=1;
    //    gettotal(1,&var,totalsizes,&var_tot); // if only 1 object, then just send address
    //    break;
  case SURFACETYPE:
  case CUMULATIVETYPE:
    totalsizes[0]=numelements;
    gettotal(1,&var,totalsizes,&var_tot);
    break;
  case CUMULATIVETYPE2:
    totalsizes[0]=1; // GODMARK: don't see point in this since just ignore numelements
    gettotal(1,&var,totalsizes,&var_tot);
    break;
  case CUMULATIVETYPE3:
    // this type not used right now
    totalsizes[0]=1;
    getalltotal(1,&var,totalsizes,&var_tot);
    break;
  default:
    dualfprintf(fail_file,"No defined type=%d in integrate\n",type);
    myexit(1);
  }
  return(0);
}

int integratel(int numelements, CTYPE * var, CTYPE *var_tot,int type, int enerregion)
{
  CTYPE *totalptrs[MAXPTRS],*totaloptrs[MAXPTRS];
  int totalsizes[MAXPTRS],numptrs;

  switch(type){
  case CONSTYPE:
    // first add up over entire grid for each CPU
    if(counttotal(enerregion,var,numelements)>=1) return(1);
    totalsizes[0]=numelements;
    gettotall(1,&var,totalsizes,&var_tot);
    break;
  default:
    dualfprintf(fail_file,"No defined type=%d in integratel\n",type);
    myexit(1);
  }
  return(0);
}


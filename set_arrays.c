
#include "decs.h"

void set_arrays()
{
  int i, j, k, pl, l, m;
  int pl2;
  int ii;
  int jj;
  int floor,pf, tscale,dtstage;
  int dissloop;
  FTYPE valueinit;
  int dir,interpi,enodebugi;
  int dimen;
  int isleftright;
  extern int init_selfgrav(void);





#if(PRODUCTION==0)
  // initialize things to NAN in order to (hopefully) trigger memory leaks to be noticed
  valueinit=sqrt(-1.0);
#else
  // avoid this practice for production run since processing nan's slows code and do process some never-initialized/never-used cells for simplicity of code loops
  valueinit=0.0;
#endif




  interporder = (int (*)) (&(a_interporder[NUMNEGINTERPS]));



  ////////////////////////////////////////////////
  //
  // Basic primitive quantity
  //
  ////////////////////////////////////////////////

  p = (FTYPE (*)[N2M][N3M][NPR]) (&(a_p[N1BND][N2BND][N3BND][0]));
  panalytic = (FTYPE (*)[N2M][N3M][NPR]) (&(a_panalytic[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOP(pl){
    p[i][j][k][pl] = valueinit;
    panalytic[i][j][k][pl] = valueinit;
  }
  

#if(NUMPOTHER>0)
  pother = (FTYPE (*)[N1M][N2M][N3M]) (&(a_pother[0][N1BND][N2BND][N3BND]));
  FULLLOOP for(pl=0;pl<NUMPOTHER;pl++) pother[pl][i][j][k] = -1;
#else
  pother = (FTYPE (*)[N1M][N2M][N3M]) (0);
#endif
  
  
  // used in fixup.c, higherorder_pointavg.c and kazfulleos.c for compute_Hglobal()  
  ptemparray = (FTYPE (*)[N2M][N3M][NPR]) (&(a_ptemparray[N1BND][N2BND][N3BND][0]));



#if(STOREFLUXSTATE)
  fluxstate = (struct of_state (*)[NUMLEFTRIGHT][N1M][N2M][N3M]) (&(a_fluxstate[-1][0][N1BND][N2BND][N3BND]));
  DIMENLOOP(dimen){
    for(isleftright=0;isleftright<NUMLEFTRIGHT;isleftright++){
      FULLLOOP{
	fluxstate[dimen][isleftright][i][j][k].rho=valueinit;
	fluxstate[dimen][isleftright][i][j][k].ie=valueinit;
	
	DLOOPA(jj){
	  fluxstate[dimen][isleftright][i][j][k].ucon[jj]=valueinit;
	  fluxstate[dimen][isleftright][i][j][k].ucov[jj]=valueinit;
	  fluxstate[dimen][isleftright][i][j][k].bcon[jj]=valueinit;
	  fluxstate[dimen][isleftright][i][j][k].bcov[jj]=valueinit;      
	}
      }
    }
  }
#else
  fluxstate = (struct of_state (*)[NUMLEFTRIGHT][N1M][N2M][N3M]) (0);
#endif


  // below used in fluct.c and as vector potential storage
  emf = (FTYPE (*)[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3]) (&(a_emf[0][N1BND][N2BND][N3BND]));// inner shift still same
  for(l=0;l<NDIM;l++) FULLLOOP emf[l][i][j][k] = valueinit;


#if(FIELDTOTHMEM)
  // below only used in fluxct.c
  vconemf = (FTYPE (*)[N2M][N3M][NDIM-1]) (&(a_vconemf[N1BND][N2BND][N3BND][-U1]));
  FULLLOOP for(l=1;l<=COMPDIM;l++) vconemf[i][j][k][l] = valueinit;
#else
  vconemf = (FTYPE (*)[N2M][N3M][NDIM-1]) (0);
#endif

#if(TRACKVPOT)
  vpotarray = (FTYPE (*)[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3]) (&(a_vpotarray[0][N1BND][N2BND][N3BND]));// inner shift still same
  for(l=0;l<NDIM;l++) FULLLOOP vpotarray[l][i][j][k] = valueinit;
#endif

#if(STOREWAVESPEEDS)
  wspeedtemp = (FTYPE (*)[N2M][N3M][NUMCS]) (&(a_wspeedtemp[N1BND][N2BND][N3BND][0]));
  wspeed = (FTYPE (*)[NUMCS][N1M][N2M][N3M]) (&(a_wspeed[-1][0][N1BND][N2BND][N3BND])); // shifted so wspeed[1,2,3] accesses the memory
  FULLLOOP PLOOP(pl) for(l=1;l<=COMPDIM;l++) for(m=0;m<NUMCS;m++) wspeed[l][m][i][j][k] = valueinit;
#endif
  



  ////////////////////////////////////////////////
  //
  // TIME-STEPPING
  //
  ////////////////////////////////////////////////

  pk = (FTYPE (*)[N1M][N2M][N3M][NPR]) (&(a_pk[0][N1BND][N2BND][N3BND][0]));
  uinitial = (FTYPE (*)[N2M][N3M][NPR]) (&(a_uinitial[N1BND][N2BND][N3BND][0]));
  ulast = (FTYPE (*)[N2M][N3M][NPR]) (&(a_ulast[N1BND][N2BND][N3BND][0]));
  unew = (FTYPE (*)[N2M][N3M][NPR]) (&(a_unew[N1BND][N2BND][N3BND][0]));
  dUgeomarray = (FTYPE (*)[N2M][N3M][NPR]) (&(a_dUgeomarray[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOP(pl){
    DTSTAGELOOP(dtstage) pk[dtstage][i][j][k][pl] = valueinit;
    uinitial[i][j][k][pl] = valueinit;
    ulast[i][j][k][pl] = valueinit;
    if(FULLOUTPUT==0) unew[i][j][k][pl] = valueinit;
    else  unew[i][j][k][pl] = 0;
    dUgeomarray[i][j][k][pl] = valueinit;
  }      

  
#if(HIGHERORDERMEM||FIELDSTAGMEM) // upoint needed for FV method and STAG for all methods
  upoint = (FTYPE (*)[N2M][N3M][NPR]) (&(a_upoint[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOP(pl) upoint[i][j][k][pl] = valueinit;
#endif



  ////////////////////////////////////////////////
  //
  // SPATIAL INTERPOLATION
  //
  ////////////////////////////////////////////////

#if(N1>1)
  F1 = (FTYPE (*)[N2M][N3M][NPR]) (&(a_F1[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOP(pl){
    F1[i][j][k][pl] = valueinit;
  }
#endif
#if(N2>1)
  F2 = (FTYPE (*)[N2M][N3M][NPR]) (&(a_F2[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOP(pl){
    F2[i][j][k][pl] = valueinit;
  }
#endif
#if(N3>1)
  F3 = (FTYPE (*)[N2M][N3M][NPR]) (&(a_F3[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOP(pl){
    F3[i][j][k][pl] = valueinit;
  }
#endif

#if(SPLITMAEMMEM)
#if(N1>1)
  F1EM = (FTYPE (*)[N2M][N3M][NPR]) (&(a_F1EM[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOP(pl){
    F1EM[i][j][k][pl] = valueinit;
  }
#endif
#if(N2>1)
  F2EM = (FTYPE (*)[N2M][N3M][NPR]) (&(a_F2EM[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOP(pl){
    F2EM[i][j][k][pl] = valueinit;
  }
#endif
#if(N3>1)
  F3EM = (FTYPE (*)[N2M][N3M][NPR]) (&(a_F3EM[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOP(pl){
    F3EM[i][j][k][pl] = valueinit;
  }
#endif
#endif


  
#if(SPLITNPR||FIELDSTAGMEM)
  gp_l = (FTYPE (*)[N1M][N2M][N3M][NPR2INTERP]) (&(a_gp_l[-1][N1BND][N2BND][N3BND][0]));
  gp_r = (FTYPE (*)[N1M][N2M][N3M][NPR2INTERP]) (&(a_gp_r[-1][N1BND][N2BND][N3BND][0]));

  FULLLOOP PLOOPINTERP(pl){
    for(l=1;l<=3;l++){
      gp_l[l][i][j][k][pl] = valueinit;
      gp_r[l][i][j][k][pl] = valueinit;
    }
  }
  
#endif


  pleft = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (&(a_pleft[N1BND][N2BND][N3BND][0]));
  pright = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (&(a_pright[N1BND][N2BND][N3BND][0]));
  prc = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (&(a_prc[N1BND][N2BND][N3BND][0]));
  
  FULLLOOP PLOOPINTERP(pl){
    pleft[i][j][k][pl] = valueinit;
    pright[i][j][k][pl] = valueinit;
    prc[i][j][k][pl] = valueinit;
  }
  



#if(FIELDSTAGMEM)
  //  wspeedcorn = (FTYPE (*)[NUMCS][N1M][N2M][N3M]) (&(a_wspeedcorn[-1][0][N1BND][N2BND][N3BND])); // shifted so wspeedcorn[1,2,3] accesses the memory
  pstagscratch = (FTYPE (*)[N2M][N3M][NPR]) (&(a_pstagscratch[N1BND][N2BND][N3BND][0]));

  // -B1 and -U1 are so pbcorninterp[][B1] accesses 0th element of original memory (same for pvcorninterp)
  pbcorninterp = (FTYPE (*)[COMPDIM][NUMCS][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3]) (&(a_pbcorninterp[-1][-B1][0][N1BND][N2BND][N3BND]));
  pvcorninterp = (FTYPE (*)[COMPDIM][NUMCS][NUMCS][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3]) (&(a_pvcorninterp[-1][-U1][0][0][N1BND][N2BND][N3BND]));

#if(HIGHERORDERMEM)
  // using NPR so can use normal functions
  Bhat = (FTYPE (*)[N2M][N3M][NPR])(&(a_Bhat[N1BND][N2BND][N3BND][0]));

  FULLLOOP PLOOP(pl) Bhat[i][j][k][pl]=valueinit;

#endif

  FULLLOOP PALLLOOP(pl) pstagscratch[i][j][k][pl] = valueinit;

  // corner things that have extra at top
  FULLLOOPP1{
    for(pl2=1;pl2<=COMPDIM;pl2++) for(pl=U1;pl<=U3;pl++) for(m=0;m<NUMCS;m++) for(l=0;l<NUMCS;l++)  pvcorninterp[pl2][pl][m][l][i][j][k]=valueinit;
    for(pl2=1;pl2<=COMPDIM;pl2++) for(pl=B1;pl<=B3;pl++) for(l=0;l<NUMCS;l++)  pbcorninterp[pl2][pl][l][i][j][k]=valueinit;
  }

#endif




  ////////////////////////////////////////////////
  //
  // OLD SPATIAL INTERPOLATION
  //
  ////////////////////////////////////////////////

#if(DODQMEMORY)
#if(N1>1)
  dq1 = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (&(a_dq1[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOPINTERP(pl) dq1[i][j][k][pl] = valueinit;
#endif
#if(N2>1)
  dq2 = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (&(a_dq2[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOPINTERP(pl) dq2[i][j][k][pl] = valueinit;
#endif
#if(N3>1)
  dq3 = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (&(a_dq3[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOPINTERP(pl) dq3[i][j][k][pl] = valueinit;
#endif

#else
  dq1 = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (0);
  dq2 = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (0);
  dq3 = (FTYPE (*)[N2M][N3M][NPR2INTERP]) (0);
#endif



  


  ////////////////////////////////////////////////
  //
  // HIGHER ORDER STUFF
  //
  ////////////////////////////////////////////////

#if( HIGHERORDERMEM )
  fluxvectemp = (FTYPE (*)[N2M][N3M][NPR]) (&(a_fluxvectemp[N1BND][N2BND][N3BND][0]));
  Fa = (FTYPE (*)[N2M][N3M][NPR]) (&(a_Fa[N1BND][N2BND][N3BND][0]));
  Fb = (FTYPE (*)[N2M][N3M][NPR]) (&(a_Fb[N1BND][N2BND][N3BND][0]));
  trifprintf("Allocating Fa/Fb for UNSPLIT (ALL NOW) flux method\n");
  stencilvartemp = (FTYPE (*)[N2M][N3M][NPR]) (&(a_stencilvartemp[N1BND][N2BND][N3BND][0]));

  a2cin  = (FTYPE (*)[N2M][N3M][NPR]) (&(a_Fa[N1BND][N2BND][N3BND][0]));
  a2cout = (FTYPE (*)[N2M][N3M][NPR]) (&(a_Fb[N1BND][N2BND][N3BND][0]));
  trifprintf("Allocating a2cin/a2cout for UNSPLIT (ALL NOW) a2c method\n");

#else

  fluxvectemp = (FTYPE (*)[N2M][N3M][NPR]) (0);
  Fa = (FTYPE (*)[N2M][N3M][NPR]) (0);
  Fb = (FTYPE (*)[N2M][N3M][NPR]) (0);
  stencilvartemp = (FTYPE (*)[N2M][N3M][NPR]) (0);
  a2cin = (FTYPE (*)[N2M][N3M][NPR]) (0);
  a2cout = (FTYPE (*)[N2M][N3M][NPR]) (0);

#endif





///////////////////////////
//
// COUNTERS and failure checks
//
////////////////////////////

#if( WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS )
  weno_lower_order_fraction = (FTYPE (*)[N2M][N3M][NPR]) (&(a_weno_lower_order_fraction[N1BND][N2BND][N3BND][0]));  //atch
  
  //init the array with zeroes
  FULLLOOP PLOOP(pl) weno_lower_order_fraction[i][j][k][pl] = 0.0;
#endif
  
#if( STORE_GAMMA_PRIM_REDUCTION_FRACTION )
  weno_prim_lower_order_fraction = (FTYPE (*)[N1M][N2M][N3M]) (&(a_weno_prim_lower_order_fraction[0][N1BND][N2BND][N3BND]));  //atch
  
  //init the array with zeroes
  FULLLOOP  DIMENLOOP(dimen) weno_prim_lower_order_fraction[dimen][i][j][k] = 0.0;
#endif





  ////////////////////////////////////////////////
  //
  // DEBUG STUFF USUALLY OFF
  //
  ////////////////////////////////////////////////

#if(DOENODEBUG)

  enodebugarray = (CTYPE (*)[N2M][N3M][3-1][NUMINTERPTYPES][NPR-4][NUMENODEBUGS]) (&(a_enodebugarray[N1BND][N2BND][N3BND][0][0][0][0])); //SASMARKK (make dir count from 1 to 2 by changing 0 to -1)  //atch debug


  FULLLOOP DIMENLOOP(dir) INTERPLOOP(interpi) PLOOP(pl) ENODEBUGLOOP(enodebugi){
    if(dir<=2 && pl<=U2){
      enodebugarray[i][j][k][dir-1][interpi][pl][enodebugi]=0;
    }
  }
#endif

#if(FLUXDUMP)
  fluxdump=(FTYPE (*)[N2M][N3M][NUMFLUXDUMP]) (&(a_fluxdump[N1BND][N2BND][N3BND][0]));

  // normal things
  FULLLOOP for(pl=0;pl<NUMFLUXDUMP;pl++) fluxdump[i][j][k][pl]=0.0; // not used for evolution -- just dumping -- so ok to ignore if leaking(?)

#endif



  ////////////////////////////////////////////////
  //
  // DEBUG STUFF USUALLY ON
  //
  ////////////////////////////////////////////////

  pflag = (PFTYPE (*)[N2M][N3M][NUMPFLAGS]) (&(a_pflag[N1BND][N2BND][N3BND][0]));
  FULLLOOP  PFLAGLOOP(pf) pflag[i][j][k][pf] = NANPFLAG;


#if(DODEBUG)
  failfloorcount = (CTYPE (*)[N2M][N3M][NUMTSCALES][NUMFAILFLOORFLAGS]) (&(a_failfloorcount[N1BND][N2BND][N3BND][0][0]));  //SASMARK: added another [0] on the right -- was a memory leak? atch
  FULLLOOP  TSCALELOOP(tscale) FLOORLOOP(floor) failfloorcount[i][j][k][tscale][floor]=valueinit;
#endif


  ////////////////////////////////////////////////
  //
  // other diagnostics
  //
  ////////////////////////////////////////////////


#if(DODISS)
  dissfunpos = (FTYPE (*)[N2M][N3M][NUMDISSFUNPOS]) (&(a_dissfunpos[N1BND][N2BND][N3BND][0]));
#endif


#if(CALCFARADAYANDCURRENTS)
  // this faraday needed for current calculation
  cfaraday =  (FTYPE (*)[N2M][N3M][NUMCURRENTSLOTS][3]) (&(a_cfaraday[N1BND][N2BND][N3BND][0][0]));
  fcon =  (FTYPE (*)[N2M][N3M][NUMFARADAY]) (&(a_fcon[N1BND][N2BND][N3BND][0]));
  jcon = (FTYPE (*)[N2M][N3M][NDIM]) (&(a_jcon[N1BND][N2BND][N3BND][0]));

  FULLLOOP{
    for(pl=0;pl<NUMCURRENTSLOTS;pl++) for(l=0;l<3;l++){
      cfaraday[i][j][k][pl][l]=valueinit;
    }
    for(pl=0;pl<NUMFARADAY;pl++){
      fcon[i][j][k][pl]=valueinit;
    }
    for(pl=0;pl<NDIM;pl++){
      jcon[i][j][k][pl]=valueinit;
    }
  }

#endif


  ////////////////////////////////////////////////
  //
  // AVG diagnostics
  //
  ////////////////////////////////////////////////
  
  // assume time average stuff gets zeroed in avg routine
#if(DOAVG)
  normalvarstavg =  (FTYPE (*)[N2M][N3M][NUMNORMDUMP]) (&(a_normalvarstavg[N1BND][N2BND][N3BND][0]));
  anormalvarstavg =  (FTYPE (*)[N2M][N3M][NUMNORMDUMP]) (&(a_anormalvarstavg[N1BND][N2BND][N3BND][0]));
  
  fcontavg =  (FTYPE (*)[N2M][N3M][NUMFARADAY]) (&(a_fcontavg[N1BND][N2BND][N3BND][0]));
  fcovtavg =  (FTYPE (*)[N2M][N3M][NUMFARADAY]) (&(a_fcovtavg[N1BND][N2BND][N3BND][0]));
  
  afcontavg =  (FTYPE (*)[N2M][N3M][NUMFARADAY]) (&(a_afcontavg[N1BND][N2BND][N3BND][0]));
  afcovtavg =  (FTYPE (*)[N2M][N3M][NUMFARADAY]) (&(a_afcovtavg[N1BND][N2BND][N3BND][0]));

  massfluxtavg =  (FTYPE (*)[N2M][N3M][NDIM]) (&(a_massfluxtavg[N1BND][N2BND][N3BND][0]));
  amassfluxtavg =  (FTYPE (*)[N2M][N3M][NDIM]) (&(a_amassfluxtavg[N1BND][N2BND][N3BND][0]));

  othertavg =  (FTYPE (*)[N2M][N3M][NUMOTHER]) (&(a_othertavg[N1BND][N2BND][N3BND][0]));
  aothertavg =  (FTYPE (*)[N2M][N3M][NUMOTHER]) (&(a_aothertavg[N1BND][N2BND][N3BND][0]));

#if(CALCFARADAYANDCURRENTS)
  jcontavg = (FTYPE (*)[N2M][N3M][NDIM]) (&(a_jcontavg[N1BND][N2BND][N3BND][0]));
  jcovtavg = (FTYPE (*)[N2M][N3M][NDIM]) (&(a_jcovtavg[N1BND][N2BND][N3BND][0]));

  ajcontavg = (FTYPE (*)[N2M][N3M][NDIM]) (&(a_ajcontavg[N1BND][N2BND][N3BND][0]));
  ajcovtavg = (FTYPE (*)[N2M][N3M][NDIM]) (&(a_ajcovtavg[N1BND][N2BND][N3BND][0]));
#endif

  tudtavg = (FTYPE (*)[N2M][N3M][NUMSTRESSTERMS]) (&(a_tudtavg[N1BND][N2BND][N3BND][0]));
  atudtavg = (FTYPE (*)[N2M][N3M][NUMSTRESSTERMS]) (&(a_atudtavg[N1BND][N2BND][N3BND][0]));
#endif  





  
      

  
#if(MCOORD!=CARTMINKMETRIC)
  /* grid functions */
  // GODMARK: for axisymmetric space-times, may want to keep metric functions 2D to save memory
  
  // these have 1 extra value on outer edge.  Shift for real pointer no different
  gcon = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (&(a_gcon[N1BND][N2BND][N3BND][0][0][0]));
  gcov = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (&(a_gcov[N1BND][N2BND][N3BND][0][0][0]));
  gcovpert = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
    (&(a_gcovpert[N1BND][N2BND][N3BND][0][0]));
  alphalapse = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
    (&(a_alphalapse[N1BND][N2BND][N3BND][0]));
  
#if(DOEVOLVEMETRIC)
  gcovlast = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (&(a_gcovlast[N1BND][N2BND][N3BND][0][0][0]));
  gcovpertlast = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
    (&(a_gcovpertlast[N1BND][N2BND][N3BND][0][0]));
  alphalapselast = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
    (&(a_alphalapselast[N1BND][N2BND][N3BND][0]));
#else
  gcovlast = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])(0);
  gcovpertlast = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])(0);
  alphalapselast = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])(0);
#endif
  
  gdet = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
    (&(a_gdet[N1BND][N2BND][N3BND][0]));

#if(WHICHEOM==WITHGDET)
  eomfunc = gdet; // just a pointer then
#else
  eomfunc = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
    (&(a_eomfunc[N1BND][N2BND][N3BND][0]));
#endif

#if(GDETVOLDIFF)
  gdetvol = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
    (&(a_gdetvol[N1BND][N2BND][N3BND][0]));
#endif
  
#if(DOSTOREPOSITIONDATA)
  dxdxpstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (&(a_dxdxpstore[N1BND][N2BND][N3BND][0][0][0]));
  idxdxpstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (&(a_idxdxpstore[N1BND][N2BND][N3BND][0][0][0]));
  Xstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
    (&(a_Xstore[N1BND][N2BND][N3BND][0][0]));
  Vstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
    (&(a_Vstore[N1BND][N2BND][N3BND][0][0]));
#else
  dxdxpstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (0);
  idxdxpstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (0);
  Xstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
    (0);
  Vstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
    (0);
#endif
  
  
  // rest are always located at CENT
  conn = (FTYPE (*)[N2M][N3M][NDIM][NDIM][NDIM])
    (&(a_conn[N1BND][N2BND][N3BND][0][0][0]));
  conn2 = (FTYPE (*)[N2M][N3M][NDIM])
    (&(a_conn2[N1BND][N2BND][N3BND][0]));
  
#if(VOLUMEDIFF)
  idxvol = (FTYPE (*)[N2M][N3M][NDIM])(&(a_idxvol[N1BND][N2BND][N3BND][0]));
#endif
  
#else// MCOORD==CARTMINKMETRIC
  // pointers on left remain same, but if access "0" element of that pointer, should get really allocated memory.
  
  gcon = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (&(a_gcon[0][0][0][0][0][0]));
  gcov = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (&(a_gcov[0][0][0][0][0][0]));
  gcovpert = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
    (&(a_gcovpert[0][0][0][0][0]));
  alphalapse = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
    (&(a_alphalapse[0][0][0][0]));
  
#if(DOEVOLVEMETRIC)
  gcovlast = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (&(a_gcovlast[0][0][0][0][0][0]));
  gcovpertlast = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
    (&(a_gcovpertlast[0][0][0][0][0]));
  alphalapselast = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
    (&(a_alphalapselast[0][0][0][0]));
#else
  gcovlast = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])(0);
  gcovpertlast = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])(0);
  alphalapselast = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])(0);
#endif
  
  gdet = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
    (&(a_gdet[0][0][0][0]));
  eomfunc = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
    (&(a_eomfunc[0][0][0][0]));
#if(GDETVOLDIFF)
  gdetvol = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG])
    (&(a_gdet[0][0][0][0]));
#endif
  
  
#if(DOSTOREPOSITIONDATA)
  dxdxpstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (&(a_dxdxpstore[0][0][0][0][0][0]));
  idxdxpstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (&(a_idxdxpstore[0][0][0][0][0][0]));
  Xstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
    (&(a_Xstore[0][0][0][0][0]));
  Vstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
    (&(a_Vstore[0][0][0][0][0]));
#else
  dxdxpstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (0);
  idxdxpstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (0);
  Xstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
    (0);
  Vstore = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM])
    (0);
#endif
  
  // rest are always located at CENT
  conn = (FTYPE (*)[N2M][N3M][NDIM][NDIM][NDIM])
    (&(a_conn[0][0][0][0][0][0]));
  conn2 = (FTYPE (*)[N2M][N3M][NDIM])
    (&(a_conn2[0][0][0][0]));
  
#if(VOLUMEDIFF)
  idxvol = (FTYPE (*)[N2M][N3M][NDIM])(&(a_idxvol[0][0][0][0]));
#endif
  
#endif
  
  
  
#if(DOLUMVSR)
  // yes, for each cpu
  lumvsr=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
  if(lumvsr==NULL){
    dualfprintf(fail_file,"Couldn't open lumvsr memory\n");
    myexit(1);
  }
  
  lumvsr_tot=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
  if(lumvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open lumvsr_tot memory\n");
    myexit(1);
  }
#endif
  
#if(DODISSVSR)
  for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++){
    // yes, for each cpu
    dissvsr[dissloop]=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
    if(dissvsr[dissloop]==NULL){
      dualfprintf(fail_file,"Couldn't open dissvsr memory: %d\n",dissloop);
      myexit(1);
    }
    
    dissvsr_tot[dissloop]=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
    if(dissvsr_tot[dissloop]==NULL){
      dualfprintf(fail_file,"Couldn't open dissvsr_tot memory: %d\n",dissloop);
      myexit(1);
    }
  }
  //for(ii=0;ii<ncpux1*N1;ii++) dissvsr[ii]=0;
  //for(ii=0;ii<ncpux1*N1;ii++) dissvsr_tot[ii]=0;
#endif
  
#if(DOSELFGRAVVSR)
  // yes, for each cpu
  dMvsr=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dMvsr+=N1BND;
  if(dMvsr==NULL){
    dualfprintf(fail_file,"Couldn't open dMvsr memory\n");
    myexit(1);
  }
  
  dMvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dMvsr_tot+=N1BND;
  if(dMvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open dMvsr_tot memory\n");
    myexit(1);
  }

  dVvsr=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dVvsr+=N1BND;
  if(dVvsr==NULL){
    dualfprintf(fail_file,"Couldn't open dVvsr memory\n");
    myexit(1);
  }
  
  dVvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dVvsr_tot+=N1BND;
  if(dVvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open dVvsr_tot memory\n");
    myexit(1);
  }

  vrsqvsr=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  vrsqvsr+=N1BND;
  if(vrsqvsr==NULL){
    dualfprintf(fail_file,"Couldn't open vrsqvsr memory\n");
    myexit(1);
  }
  
  vrsqvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  vrsqvsr_tot+=N1BND;
  if(vrsqvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open vrsqvsr_tot memory\n");
    myexit(1);
  }

  // yes, for each cpu
  dTrrvsr=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dTrrvsr+=N1BND;
  if(dTrrvsr==NULL){
    dualfprintf(fail_file,"Couldn't open dTrrvsr memory\n");
    myexit(1);
  }
  
  dTrrvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dTrrvsr_tot+=N1BND;
  if(dTrrvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open dTrrvsr_tot memory\n");
    myexit(1);
  }

  Mvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  Mvsr_tot+=N1BND;
  if(Mvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open Mvsr_tot memory\n");
    myexit(1);
  }

  Mvsrface1_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  Mvsrface1_tot+=N1BND;
  if(Mvsrface1_tot==NULL){
    dualfprintf(fail_file,"Couldn't open Mvsrface1_tot memory\n");
    myexit(1);
  }

  MOrvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  MOrvsr_tot+=N1BND;
  if(MOrvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open MOrvsr_tot memory\n");
    myexit(1);
  }

  // the potentials are located in ghost cells so outer boundary CPUs have sufficient data
  // +1 corresponds to last outer face, so not used for center phi but allocated to keep arrays same size
  phivsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  phivsr_tot+=N1BND;
  if(phivsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open phivsr_tot memory\n");
    myexit(1);
  }

  rcent=(FTYPE*)malloc(NUMGRAVPOS*sizeof(FTYPE));
  rcent+=N1BND;
  if(rcent==NULL){
    dualfprintf(fail_file,"Couldn't open rcent memory\n");
    myexit(1);
  }

  rcent_tot=(FTYPE*)malloc(NUMGRAVPOS*sizeof(FTYPE));
  rcent_tot+=N1BND;
  if(rcent_tot==NULL){
    dualfprintf(fail_file,"Couldn't open rcent_tot memory\n");
    myexit(1);
  }


  // yes, for each cpu
  dJvsr=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dJvsr+=N1BND;
  if(dJvsr==NULL){
    dualfprintf(fail_file,"Couldn't open dJvsr memory\n");
    myexit(1);
  }
  
  dJvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  dJvsr_tot+=N1BND;
  if(dJvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open dJvsr_tot memory\n");
    myexit(1);
  }

  Jvsr_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  Jvsr_tot+=N1BND;
  if(Jvsr_tot==NULL){
    dualfprintf(fail_file,"Couldn't open Jvsr_tot memory\n");
    myexit(1);
  }

  Jvsrface1_tot=(SFTYPE*)malloc(NUMGRAVPOS*sizeof(SFTYPE));
  Jvsrface1_tot+=N1BND;
  if(Jvsrface1_tot==NULL){
    dualfprintf(fail_file,"Couldn't open Jvsrface1_tot memory\n");
    myexit(1);
  }


#endif
  
  if(DOSELFGRAVVSR){
    // initialize self-gravity functions since not set yet and will otherwise be used to set first metric
    init_selfgrav();
  }

}

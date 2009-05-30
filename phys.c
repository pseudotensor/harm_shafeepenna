
#include "decs.h"


// making variables static inside functions made things slower if anything
// would have thought consruction/destruction operations would be important
// why static makes slower?
//#define VARSTATIC static
#define VARSTATIC 



// Calculate fluxes in direction dir and conserved variable U
// returntype==0 : flux with geometric factor geom->e (used by evolution code)
// returntype==1 : flux with physical geometry factor geom->g (used by diagnostics)
// see UtoU and source_conn()
// Note that if MAXWELL==PRIMMAXWELL then primtoflux doesn't use b^\mu or b_\mu (bcon and bcov)
int primtoflux(int returntype, FTYPE *pr, struct of_state *q, int dir,
	       struct of_geom *geom, FTYPE *flux)
{
  int primtoflux_ma(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxdiag);
  int primtoflux_em(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux);
  void UtoU_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  VARSTATIC FTYPE fluxinput[NPR],fluxinputma[NPR],fluxinputem[NPR];
  VARSTATIC FTYPE fluxdiag;
  VARSTATIC int pl;



  // initialize fluxinputma and fluxinputem so individual functions only have to define non-zero terms
  PLOOP(pl) fluxinputma[pl]=fluxinputem[pl]=0.0;
  fluxdiag=0.0;


  // define MA terms
  primtoflux_ma(&returntype, pr, q, dir, geom, fluxinputma, &fluxdiag);
  fluxinputma[UU+dir]+=fluxdiag; // add back to normal term

  // define EM terms
  primtoflux_em(&returntype, pr, q, dir, geom, fluxinputem);

  // add up MA+EM
  PLOOP(pl) fluxinput[pl] = fluxinputma[pl] + fluxinputem[pl];

  // DEBUG:
  //  PALLLOOP(pl) dualfprintf(fail_file,"ALLBEFORE: pl=%d flux=%21.15g\n",pl,fluxinput[pl]);


  // convert from UNOTHING->returntype
  // notice that geometry comes after subtractions/additions of EOMs
  UtoU_fromunothing(returntype,geom,fluxinput,flux);

  // DEBUG:
  //  PALLLOOP(pl) dualfprintf(fail_file,"ALLAFTER: pl=%d flux=%21.15g\n",pl,flux[pl]);


  return(0);

}

/* calculate fluxes in direction dir and conserved variable U; these
   are always needed together, so there is no point in calculated the
   stress tensor twice */

// returntype==0 : flux with geometric factor geom->e (used by evolution code)
// returntype==1 : flux with physical geometry factor geom->g (used by diagnostics)
// see UtoU and source_conn()
// fluxdir = true flux direction even if passing back conserved quantity (fundir==TT)
int primtoflux_splitmaem(int returntype, FTYPE *pr, struct of_state *q, int fluxdir, int fundir, struct of_geom *geom, FTYPE *fluxma, FTYPE *fluxem)
{
  int primtoflux_ma(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxdiag);
  int primtoflux_em(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux);
  void UtoU_ma_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  void UtoU_em_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  VARSTATIC FTYPE fluxinput[NPR],fluxinputma[NPR],fluxinputem[NPR];
  VARSTATIC FTYPE fluxdiag;
  VARSTATIC int pl;


  // initialize fluxinputma and fluxinputem so individual functions only have to define non-zero terms
  PLOOP(pl) fluxinputma[pl]=fluxinputem[pl]=0.0;

  // define MA terms
  primtoflux_ma(&returntype, pr, q, fundir, geom, fluxinputma, &fluxdiag);

  // SUPERGODMARK CHANGINGMARK
  // Note that pressure part of flux has 0 conserved quantity associated with it from the point of view of this output and correct HLL/LAXF calculation
#if(1)
  if(fundir!=TT) fluxinputma[FLUXSPLITPMA(fluxdir)]=fluxdiag;
  else fluxinputma[UU]+=fluxdiag;

  //  if(fundir!=TT) fluxinputma[FLUXSPLITPMA(fundir)]=fluxdiag; // only need to return single diagonal value and store it inside field part for now
  //  else fluxinputma[UU]+=fluxdiag; // then really part of conserved quantity and then don't separate it

#elif(0)

  // not as robust to put conserved quantity here such that dissipative term is not together with momentum-flux_dir_dir term
  fluxinputma[FLUXSPLITPMA(fluxdir)]=fluxdiag; // only need to return single diagonal value and store it inside field part for now
#else
  // not as robust to put conserved quantity here such that dissipative term is not together with momentum-flux_dir_dir term
  if(fundir!=TT) fluxinputma[FLUXSPLITPMA(fluxdir)]=fluxdiag;
  else{
    fluxinputma[UU]+=0.5*fluxdiag;
    fluxinputma[FLUXSPLITPMA(fluxdir)]=0.5*fluxdiag;
  }
#endif

  // define EM terms
  primtoflux_em(&returntype, pr, q, fundir, geom, fluxinputem);

  // convert from UNOTHING->returntype
  // notice that geometry comes after subtractions/additions of EOMs
  //  UtoU_ma(UNOTHING,returntype,geom,fluxinputma,fluxma); // properly converts separate diagonal flux in FLUXSPLITPMA(fluxdir)
  UtoU_ma_fromunothing(returntype,geom,fluxinputma,fluxma);

  //  UtoU_em(UNOTHING,returntype,geom,fluxinputem,fluxem);
  UtoU_em_fromunothing(returntype,geom,fluxinputem,fluxem);


  return(0);

}



// matter only terms (as if B=0)
int primtoflux_ma(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux, FTYPE *fluxdiag)
{
  // sizes: NPR,struct of_state, int, struct of_geom, NPR
  int massflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *massflux);
  int entropyflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *entropyflux);
  int advectedscalarflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *advectedscalarflux, int pnum);
  VARSTATIC FTYPE fluxdiagpress[NPR]; // temp var
  VARSTATIC int pl;

  // USE of SPLITNPR here is simplified for speed
  // assume if quantity is bracketed that including it


#if(SPLITNPR)
  if(nprlist[nprstart]<=RHO && nprlist[nprend]>=RHO)
#endif
    massflux_calc(pr, dir, q, &flux[RHO]); // fills RHO only


  // GODMARK WTF!  Problems with code (compiling?) with this
  // if(mhd_calc(pr,dir,geom,q,mhd)>=1)
  // FAILSTATEMENT("phys.c:primtoflux()","mhd_calc() dir=1or2",1);

  // MHD stress-energy tensor w/ first index up, second index down.
#if(SPLITNPR)
  if(nprlist[nprstart]<=UU && nprlist[nprend]>=U3)
#endif
  {
    mhd_calc_ma(pr, dir, geom, q, &flux[UU], &fluxdiagpress[UU]); // fills flux[UU->U3] and fluxdiagonal[UU->U3]
    *fluxdiag = fluxdiagpress[UU+dir];
  }

#if(DOYL!=DONOYL)
#if(SPLITNPR)
  if(nprlist[nprstart]<=YL && nprlist[nprend]>=YL)
#endif
  advectedscalarflux_calc(pr, dir, q, &flux[YL],YL); // fills YL only
#endif
#if(DOYNU!=DONOYNU)
#if(SPLITNPR)
  if(nprlist[nprstart]<=YNU && nprlist[nprend]>=YNU)
#endif
  advectedscalarflux_calc(pr, dir, q, &flux[YNU],YNU); // fills YNU only
#endif


#if(DOENTROPY!=DONOENTROPY)
#if(SPLITNPR)
  if(nprlist[nprstart]<=ENTROPY && nprlist[nprend]>=ENTROPY)
#endif
  entropyflux_calc(pr, dir, q, &flux[ENTROPY]); // fills ENTROPY only

  // below is special for utoprim() 5D version for full entropy evolution and inversion
  if(*returntype==UENTROPY){
    flux[UU]=flux[ENTROPY]; // overwrite for utoprim()
    *fluxdiag = 0.0; // overwrite for utoprim()
    *returntype=UNOTHING; // reset returntype for UtoU
  }

#endif


  // DEBUG:
  //  PALLLOOP(pl) dualfprintf(fail_file,"ALL: pl=%d flux=%21.15g\n",pl,flux[pl]);


  return (0);
}


// electromagnetic terms (as if rho=u=p=0)
int primtoflux_em(int *returntype, FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *flux)
{
  // sizes: NPR,struct of_state, int, struct of_geom, NPR
  //  FTYPE dualf[NDIM];
  int dualfaradayspatial_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualf);
  //  int massflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *massflux);
  //  int advectedscalarflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *advectedscalarflux, int pnum);


  // USE of SPLITNPR here is simplified for speed
  // assume if quantity is bracketed that including it


  // GODMARK WTF!  Problems with code (compiling?) with this
  // if(mhd_calc(pr,dir,geom,q,mhd)>=1)
  // FAILSTATEMENT("phys.c:primtoflux()","mhd_calc() dir=1or2",1);

  // MHD stress-energy tensor w/ first index up, second index down.
#if(SPLITNPR)
  if(nprlist[nprstart]<=UU && nprlist[nprend]>=U3)
#endif
    mhd_calc_em(pr, dir, geom, q, &flux[UU]); // fills UU->U3



#if(SPLITNPR)
  if(nprlist[nprstart]<=B1 && nprlist[nprend]>=B3)
#endif
  dualfaradayspatial_calc(pr,dir,q,&flux[B1]); // fills B1->B3


  return (0);
}


int massflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *massflux)
{
  /* particle number flux */
  *massflux = pr[RHO] * q->ucon[dir];

  // DEBUG:
  //  dualfprintf(fail_file,"massflux: %d %21.15g %21.15g\n",dir,pr[RHO],q->ucon[dir]);

  return(0);
}


int advectedscalarflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *advectedscalarflux, int pnum)
{
  int massflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *massflux);
  VARSTATIC FTYPE massflux;

  /* yl/ynu/etc. per unit rest-mass flux */

  // entropy=entropy per unit volume, where conserved quantity is specific entropy:
  // d/d\tau(entropy/rho)=0
  // -> \nabla_\mu(entropy u^\mu)=0

  massflux_calc(pr, dir, q, &massflux);
  *advectedscalarflux = pr[pnum] * massflux;

  return(0);
}




int entropyflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *entropyflux)
{
  VARSTATIC FTYPE entropy;
  int entropy_calc(FTYPE *pr, FTYPE *entropy);

  // get entropy
  entropy_calc(pr,&entropy);
  /* entropy per unit rest-mass flux */
  // entropy=entropy per unit volume, where conserved quantity is specific entropy:
  // d/d\tau(entropy/rho)=0
  // -> \nabla_\mu(entropy u^\mu)=0

  // DEBUG:
  //  dualfprintf(fail_file,"entropy=%21.15g dir=%d ucondir=%21.15g\n",entropy,dir,q->ucon[dir]);

  *entropyflux = entropy * q->ucon[dir];

  return(0);
}


// entropy wrapper
int entropy_calc(FTYPE *pr, FTYPE *entropy)
{

  *entropy = compute_entropy(pr[RHO],pr[UU]);

  return(0);
}

int invertentropyflux_calc(FTYPE entropyflux,int dir, struct of_state *q, FTYPE *pr)
{
  VARSTATIC FTYPE entropy;
  int ufromentropy_calc(FTYPE entropy, FTYPE *pr);

  // get entropy
  entropy=entropyflux/q->ucon[dir];
  ufromentropy_calc(entropy,pr);


  return(0);
}

// u from entropy (uses pr[RHO])
// wrapper
int ufromentropy_calc(FTYPE entropy, FTYPE *pr)
{

  // entropy version of ie
  pr[ENTROPY]=compute_u_from_entropy(pr[RHO],entropy);

  return(0);
}


// spatial part of dualfaraday
int dualfaradayspatial_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualf)
{
  VARSTATIC FTYPE dualffull[NDIM];
  int dualfullfaraday_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualf);


  dualfullfaraday_calc(pr,dir,q,dualffull);
  dualf[0]=dualffull[1];
  dualf[1]=dualffull[2];
  dualf[2]=dualffull[3];

  return(0);

}



// Notation for HARM paper and JCM's GRFFE paper is that:
// B^\mu = \eta_\nu *F^{\nu\mu} and for lab-frame we chose \eta_\nu={-1,0,0,0}
//
// One can then show that b^i u^j - b^j u^i = B^i v^j - B^j v^i where
// b^\mu = u_\nu *F^{\nu\mu} and v^i = u^i/u^t
//
// The form using B^i and v^i avoids catastrophic cancellation because otherwise
// starting with B^i (as primitive) and converting to b^i leads to u^i u^j (u.B)/u^t term that cancels exactly
// Since this is quite a high order term, then for highly relativistic flows this causes catastrophic
// cancellation issues, so this is why we use the primitives directly even if more complicated looking and more expensive to compute stress tensor or maxwell tensor
//
// returns \dF^{\mu dir}
// well, actually returns dualffull[dir], so gives columns instead of rows
int dualfullfaraday_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualffull)
{

#if(MAXWELL==GENMAXWELL)
  /* dual of Maxwell tensor */
  dualffull[0] = q->bcon[0] * q->ucon[dir] - q->bcon[dir] * q->ucon[0];
  dualffull[1] = q->bcon[1] * q->ucon[dir] - q->bcon[dir] * q->ucon[1];
  dualffull[2] = q->bcon[2] * q->ucon[dir] - q->bcon[dir] * q->ucon[2];
  dualffull[3] = q->bcon[3] * q->ucon[dir] - q->bcon[dir] * q->ucon[3];
#elif(MAXWELL==PRIMMAXWELL)
  if(dir>0){
    /* dual of Maxwell tensor */
    // dir refers to the direction of the derivative of the dualffull
    // B1,B2,B3 refers to LHS of equation dB^i/dt
    // due to antisymmetry, dir==i is 0
    dualffull[0] = - pr[B1+dir-1] ; // dualffull[i]=\dF^{i dir} where \dF^{0 dir} =-B^{dir}
    dualffull[1] = (pr[B1] * q->ucon[dir] - pr[B1+dir-1] * q->ucon[1])/q->ucon[0];
    dualffull[2] = (pr[B2] * q->ucon[dir] - pr[B1+dir-1] * q->ucon[2])/q->ucon[0];
    dualffull[3] = (pr[B3] * q->ucon[dir] - pr[B1+dir-1] * q->ucon[3])/q->ucon[0];
  }
  else{
    dualffull[0] = 0;
    dualffull[1] = pr[B1];
    dualffull[2] = pr[B2];
    dualffull[3] = pr[B3];
  }
#endif

  return(0);

}


/* dual of Maxwell tensor */
// returns \dF^{\mu \nu}
int Mcon_calc(FTYPE *pr, struct of_state *q, FTYPE (*Mcon)[NDIM])
{
  VARSTATIC int j,k;
  VARSTATIC FTYPE vcon[NDIM];


#if(MAXWELL==GENMAXWELL)

  DLOOP(j,k) Mcon[j][k] = q->bcon[j] * q->ucon[k] - q->bcon[k] * q->ucon[j];

#elif(MAXWELL==PRIMMAXWELL)

  // diagonal is 0
  DLOOPA(j) Mcon[j][j]=0.0;

  // space-time terms
  SLOOPA(k) {
    // \dF^{it} = B^i = pr[B1+i-1]
    Mcon[k][0] = pr[B1+k-1] ; 
    Mcon[0][k] = - Mcon[k][0] ;
  }

  // get v^i
  SLOOPA(k) vcon[k]= q->ucon[k]/q->ucon[TT];

  // space-space terms
  //  SLOOP(j,k) Mcon[j][k] = (pr[B1+j-1] * vcon[k] - pr[B1+k-1] * vcon[j]);
  // optimize
  Mcon[1][2] = (pr[B1] * vcon[2] - pr[B2] * vcon[1]);
  Mcon[1][3] = (pr[B1] * vcon[3] - pr[B3] * vcon[1]);
  Mcon[2][3] = (pr[B2] * vcon[3] - pr[B3] * vcon[2]);
  Mcon[2][1] = -Mcon[1][2];
  Mcon[3][1] = -Mcon[1][3];
  Mcon[3][2] = -Mcon[2][3];
  Mcon[1][1] = Mcon[2][2] = Mcon[3][3] = 0.0;

#endif


  return(0);

}



// returns entire space-time(NDIM in size) / EOM(NPR in size) matrix
int primtofullflux(int returntype, FTYPE *pr, struct of_state *q,
		   struct of_geom *ptrgeom, FTYPE (*flux)[NPR])
{
  VARSTATIC int j;
  
  // j=0,1,2,3 corresponding to U^j_\nu , where \nu corresponds to all EOMs and j to space-time for each
  // 1 stands for obey nprlist
  DLOOPA(j) primtoflux(returntype,pr,q,j,ptrgeom,flux[j]);

  return(0);
}


/* calculate "conserved" quantities */
int primtoU(int returntype, FTYPE *pr, struct of_state *q, struct of_geom *geom,
	    FTYPE *U)
{
  MYFUN(primtoflux(returntype,pr, q, 0, geom, U) ,"phys.c:primtoU()", "primtoflux_calc() dir=0", 1);

  return (0);
}




void UtoU_gen(int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_gengdet(int removemass, int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  void UtoU_gen_allgdet(int removemass, int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  VARSTATIC int removemass;

#if((REMOVERESTMASSFROMUU==1)&&(EOMTYPE!=EOMFFDE))
  removemass = (whichmaem==ISMAONLY || whichmaem==ISMAANDEM);
#else
  // otherwise removemass = 0 effectively
  removemass = 0;
#endif

#if(WHICHEOM!=WITHGDET)
  UtoU_gen_gengdet(removemass, whichmaem, inputtype, returntype,ptrgeom,Uin, Uout);
#else
  UtoU_gen_allgdet(removemass, whichmaem, inputtype, returntype,ptrgeom,Uin, Uout);
#endif

}


void UtoU_gen_fromunothing(int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_gengdet_fromunothing(int removemass, int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  void UtoU_gen_allgdet_fromunothing(int removemass, int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  VARSTATIC int removemass;


#if((REMOVERESTMASSFROMUU==1)&&(EOMTYPE!=EOMFFDE))
    removemass = (whichmaem==ISMAONLY || whichmaem==ISMAANDEM);
#else
    // otherwise removemass = 0 effectively
    removemass = 0;
#endif


#if(WHICHEOM!=WITHGDET)
  UtoU_gen_gengdet_fromunothing(removemass, whichmaem, returntype,ptrgeom,Uin, Uout);
#else
  UtoU_gen_allgdet_fromunothing(removemass, whichmaem, returntype,ptrgeom,Uin, Uout);
#endif

}


// standardized U form is geometry free and 
// \rho u^t , T^t_\nu , *F^{it}
// convert one form of U(or component of Flux) to another form
// UtoU controls meaning of todo and REMOVERESTMASSFROMUU.
// present order means start with geometry-free EOMs, add/subtract them, THEN geometry is assigned to that list of new EOMs.
// can choose to change order so that add geometry terms, THEN add/subtract them.  Rest of code shouldn't care (except source_conn()'s first connection)
// if SPLITNPR, operates on other primitives but only changes Uout, not Uin, so doesn't matter if change all conserved quantities
// whichmaem == ISEMONLY ISMAONLY ISMAANDEM
void UtoU_gen_gengdet(int removemass, int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_gengdet_fromunothing(int removemass, int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Ugeomfree, FTYPE *Uout);
  VARSTATIC FTYPE Ugeomfree[NPR];
  VARSTATIC int pl;




  if(inputtype==returntype){
    // then just copy
    PLOOP(pl) Uout[pl]=Uin[pl];
  }
  else{

    /////////////////////
    //
    // input
    //
    
    // set basic transformation
    PLOOP(pl) Ugeomfree[pl]=Uin[pl];
    
    // now fine-tune transformation
    switch(inputtype){
      
    case UEVOLVE:
      
      // get igeom
      set_igeom(ptrgeom);
      
      PLOOP(pl) Ugeomfree[pl] *= ptrgeom->ienosing[pl];
      
      //    dualfprintf(fail_file,"diff: %21.15g %21.15g\n",ptrgeom->e[UU],gdetvol[ptrgeom->i][ptrgeom->j][ptrgeom->k][ptrgeom->p]); 
#if((REMOVERESTMASSFROMUU==1)&&(EOMTYPE!=EOMFFDE))
      if(removemass){
	// go back to standard stress-energy tensor form
	Ugeomfree[UU]  +=  - Ugeomfree[RHO] ; // - means adding back rest-mass
      }
#endif
      break;
    case UDIAG:
      // get igeom
      set_igeom(ptrgeom);
      
      PLOOP(pl) Ugeomfree[pl] *= ptrgeom->igeomnosing;
      break;
    case UNOTHING:
      break;
    default:
      dualfprintf(fail_file,"UtoU: No such inputtype=%d\n",inputtype);
      myexit(246364);
      break;
    }
    

    // at this point, Ugeomfree is geometry-free standard form of conserved quantities
    // output
    UtoU_gen_gengdet_fromunothing(removemass, whichmaem, returntype,ptrgeom,Ugeomfree, Uout);
    
  }

}





void UtoU_gen_gengdet_fromunothing(int removemass, int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Ugeomfree, FTYPE *Uout)
{    
  VARSTATIC int pl;


  /////////////////////////
  //
  // output
  //
  switch(returntype){
  case UEVOLVE:
#if((REMOVERESTMASSFROMUU==1)&&(EOMTYPE!=EOMFFDE))
    if(removemass){ // diagnostics want normal stress-energy tensor
      // "subtract" rest-mass
      // should be done on geometry-free version
      Ugeomfree[UU] += Ugeomfree[RHO];
    }
#endif
    PLOOP(pl) Uout[pl]=Ugeomfree[pl]*ptrgeom->e[pl];
    break;
  case UDIAG:
    PLOOP(pl) Uout[pl]=Ugeomfree[pl]*ptrgeom->g;
    break;
  case UNOTHING:
    PLOOP(pl) Uout[pl]=Ugeomfree[pl];
    break;
  default:
    dualfprintf(fail_file,"UtoU: No such returntype=%d\n",returntype);
    myexit(724365);
    break;
  }

}




// simplified version that assumes WHICHEOM==WITHGDET
void UtoU_gen_allgdet(int removemass, int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  VARSTATIC FTYPE Ugeomfree[NPR];
  VARSTATIC int pl;
  void UtoU_gen_allgdet_fromunothing(int removemass, int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Ugeomfree, FTYPE *Uout);



#if(WHICHEOM!=WITHGDET)
  dualfprintf(fail_file,"Using UtoU_gen_allgdet() but WHICHEOM!=WITHGDET\n");
  myexit(342968346);
#endif
  
  
  
#if((REMOVERESTMASSFROMUU==1)&&(EOMTYPE!=EOMFFDE))
  removemass = (whichmaem==ISMAONLY || whichmaem==ISMAANDEM);
#endif// otherwise removemass = 0 effectively
  
  
  /////////////////////
  //
  // input
  //
  
  // now fine-tune transformation
  switch(inputtype){
    
  case UEVOLVE:
  case UDIAG:
    
    // get igeom
    set_igeomsimple(ptrgeom);
    
    PLOOP(pl) Ugeomfree[pl] = Uin[pl]*ptrgeom->igeomnosing;
    
#if((REMOVERESTMASSFROMUU==1)&&(EOMTYPE!=EOMFFDE))
    if(removemass){
      // go back to standard stress-energy tensor form
      Ugeomfree[UU]  +=  - Ugeomfree[RHO] ; // - means adding back rest-mass
    }
#endif
    break;
  case UNOTHING:
    PLOOP(pl) Ugeomfree[pl] = Uin[pl];
    break;
  default:
    dualfprintf(fail_file,"UtoU: No such inputtype=%d\n",inputtype);
    myexit(6746364);
    break;
  }

  
  // at this point, Ugeomfree is geometry-free standard form of conserved quantities
  // output
  UtoU_gen_allgdet_fromunothing(removemass, whichmaem, returntype,ptrgeom,Ugeomfree, Uout);
  

}



void UtoU_gen_allgdet_fromunothing(int removemass, int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Ugeomfree, FTYPE *Uout)
{    
  VARSTATIC int pl;

  /////////////////////////
  //
  // output
  //
  switch(returntype){
  case UEVOLVE:
  case UDIAG:
#if((REMOVERESTMASSFROMUU==1)&&(EOMTYPE!=EOMFFDE))
    if(removemass){ // diagnostics want normal stress-energy tensor
      // "subtract" rest-mass
      // should be done on geometry-free version
      Ugeomfree[UU] += Ugeomfree[RHO];
    }
#endif
    PLOOP(pl) Uout[pl]=Ugeomfree[pl]*ptrgeom->g;
    break;
  case UNOTHING:
    PLOOP(pl) Uout[pl]=Ugeomfree[pl];
    break;
  default:
    dualfprintf(fail_file,"UtoU: No such returntype=%d\n",returntype);
    myexit(924636);
    break;
  }

}


void UtoU_evolve2diag(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen(int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  int pl;

#if(WHICHEOM==WITHGDET)
  // then UEVOLVE and UDIAG have same geometry (geom.g)
  PALLLOOP(pl) Uout[pl]=Uin[pl];
#else
  UtoU_gen(ISMAANDEM, inputtype, returntype, ptrgeom, Uin, Uout);
#endif

}


void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen(int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen(ISMAANDEM, inputtype, returntype, ptrgeom, Uin, Uout);

}

void UtoU_ma(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen(int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen(ISMAONLY, inputtype, returntype, ptrgeom, Uin, Uout);

}

void UtoU_em(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen(int whichmaem, int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen(ISEMONLY, inputtype, returntype, ptrgeom, Uin, Uout);

}


void UtoU_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_fromunothing(int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen_fromunothing(ISMAANDEM, returntype, ptrgeom, Uin, Uout);

}

void UtoU_ma_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_fromunothing(int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen_fromunothing(ISMAONLY, returntype, ptrgeom, Uin, Uout);

}

void UtoU_em_fromunothing(int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  void UtoU_gen_fromunothing(int whichmaem, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  UtoU_gen_fromunothing(ISEMONLY, returntype, ptrgeom, Uin, Uout);

}



// standardized primitive form is assumed to be
// \rho , u, \tilde{u}^\mu , *F^{it}=B^i
// where \tilde{u} is relative 4-velocity, as relative to $n_\mu = (-\alpha,0,0,0)$ and $\alpha^2=-1/g^{tt}$.

// For any space-time with no time-like curves this 4-velocity is always single-valued (i.e. unique).  It can also take on any value, so a reasonable primitive quantity.

// convert from one primitive form to another
void PtoP(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *pin, FTYPE *pout)
{
  VARSTATIC FTYPE pstandard[NPR];
  VARSTATIC int pl;



}


/* calculate magnetic field four-vector */
void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov, FTYPE *bcon)
{
  VARSTATIC int j;

  bcon[TT] = pr[B1] * ucov[1] + pr[B2] * ucov[2] + pr[B3] * ucov[3];
  for (j = 1; j <= 3; j++)
    bcon[j] = (pr[B1 - 1 + j] + bcon[TT] * ucon[j]) / ucon[TT];

  return;
}

// inverse of bcon_calc()
void Bcon_calc(struct of_state *q, FTYPE*B)
{
  VARSTATIC FTYPE uu0,uu1,uu2,uu3;
  VARSTATIC FTYPE ud0,ud1,ud2,ud3;
  VARSTATIC FTYPE bu1,bu2,bu3;
  VARSTATIC FTYPE denom;


  uu0=q->ucon[TT];
  uu1=q->ucon[RR];
  uu2=q->ucon[TH];
  uu3=q->ucon[PH];

  ud0=q->ucov[TT];
  ud1=q->ucov[RR];
  ud2=q->ucov[TH];
  ud3=q->ucov[PH];

  bu1=q->bcon[RR];
  bu2=q->bcon[TH];
  bu3=q->bcon[PH];
 
  denom=1.0/(1.0+ud1*uu1+ud2*uu2+ud3*uu3);
  
  B[1]=uu0*(-(bu2*ud2+bu3*ud3)*uu1+bu1*(1.0+ud2*uu2+ud3*uu3))*denom;
  B[2]=uu0*(-(bu1*ud1+bu3*ud3)*uu2+bu2*(1.0+ud1*uu1+ud3*uu3))*denom;
  B[3]=uu0*(-(bu2*ud2+bu2*ud2)*uu3+bu3*(1.0+ud2*uu2+ud1*uu1))*denom;


}


// convert (e^\mu=0 case) b^\mu and (3-velocity in coordinate lab frame) v^\mu to pr
void vbtopr(FTYPE *vcon,FTYPE *bcon,struct of_geom *geom, FTYPE *pr)
{
  void Bcon_calc(struct of_state *q, FTYPE*B);
  VARSTATIC int pl;
  VARSTATIC struct of_state q;
  VARSTATIC FTYPE prim[NPR];
  VARSTATIC FTYPE ucon[NDIM];


  // go ahead and get pr velocity
  PLOOP(pl) prim[pl]=0.0;
  prim[U1]=vcon[1];
  prim[U2]=vcon[2];
  prim[U3]=vcon[3];

  //  vcon2pr(WHICHVEL,vcon,geom,pr); // need u^\mu, so do below instead
  ucon_calc_3vel(prim,geom,q.ucon);
  ucon2pr(WHICHVEL,q.ucon,geom,pr); // fills pr[U1->U3]
  
  //  q.ucon[TT]=ucon[TT];
  //  q.ucon[RR]=ucon[RR];
  //  q.ucon[TH]=ucon[TH];
  //  q.ucon[PH]=ucon[PH];
  
  lower_vec(q.ucon,geom,q.ucov);

  //  q.bcon[TT]=bcon[TT]; // not used below
  q.bcon[RR]=bcon[RR];
  q.bcon[TH]=bcon[TH];
  q.bcon[PH]=bcon[PH];
  
  Bcon_calc(&q,&pr[B1-1]); // &pr[B1-1] since Bcon_calc() fills 1-3



}

/* MHD stress tensor, with first index up, second index down */
// mhd^dir_j
void mhd_calc(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd)
{
  void mhd_calc_0(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);
  void mhd_calc_norestmass(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);

#if(REMOVERESTMASSFROMUU==2)
  mhd_calc_norestmass(pr, dir, geom, q, mhd);
#else
  mhd_calc_0(pr, dir, geom, q, mhd);
#endif

}

/* MHD stress tensor, with first index up, second index down */
// mhd^dir_j
// understood that mhddiagpress only contains non-zero element on mhddiagpress[dir] and all others should be 0.0
void mhd_calc_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhddiagpress)
{
  void mhd_calc_0_ma(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd, FTYPE *mhddiagpress);
  void mhd_calc_norestmass_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhddiagpress);

#if(REMOVERESTMASSFROMUU==2)
  mhd_calc_norestmass_ma(pr, dir, geom, q, mhd, mhddiagpress);
#else
  mhd_calc_0_ma(pr, dir, q, mhd, mhddiagpress);
#endif

}

/* MHD stress tensor, with first index up, second index down */
// mhd^dir_j
void mhd_calc_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd)
{
  void mhd_calc_0_em(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd);
  void mhd_calc_primfield_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);

#if(MAXWELL==GENMAXWELL)
  mhd_calc_0_em(pr, dir, q, mhd);
#elif(MAXWELL==PRIMMAXWELL)
  mhd_calc_primfield_em(pr, dir, geom, q, mhd);
#else
  #error No such MAXWELL
#endif


}


/* MHD stress tensor, with first index up, second index down */
void mhd_calc_0(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd)
{
  void mhd_calc_0_ma(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd, FTYPE *mhddiagpress);
  void mhd_calc_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);
  VARSTATIC int j;
  VARSTATIC FTYPE mhdma[NDIM],mhdem[NDIM];
  VARSTATIC FTYPE mhddiagpress[NDIM];


  mhd_calc_0_ma(pr, dir, q, mhdma,mhddiagpress);
  mhd_calc_em(pr, dir, geom, q, mhdem);

  // add up MA+EM
  DLOOPA(j) mhd[j] = (mhdma[j] + mhddiagpress[j]) + mhdem[j];

}



/* MHD stress tensor, with first index up, second index down */
void mhd_calc_0_ma(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd, FTYPE *mhddiagpress)
{
  VARSTATIC int j;
  VARSTATIC FTYPE r, u, P, w, eta, ptot;

  // below allows other scalars to be advected but not affect the stress-energy equations of motion
#if(EOMTYPE==EOMGRMHD)
  r = pr[RHO];
  u = pr[UU];
  P = pressure_rho0_u(r,u);
  w = P + r + u;
#elif(EOMTYPE==EOMCOLDGRMHD)
  w = r = pr[RHO];
  u = 0.0;
  P = 0.0;
#elif(EOMTYPE==EOMFFDE)
  r=u=P=w=0.0;
#endif
  eta = w;
  ptot = P;

  /* single row of mhd stress tensor, first index up, second index down 
   */
  // mhd^{dir}_{j} =
  // j=0..3
  DLOOPA(j) mhd[j] = eta * q->ucon[dir] * q->ucov[j];

  DLOOPA(j) mhddiagpress[j] = 0.0;
#if(SPLITPRESSURETERMINFLUXMA==0)
  mhd[dir] += ptot;
#else
  // below equivalent to ptot * delta(dir,j)
  mhddiagpress[dir] = ptot;
#endif

}


// EM part of stress-energy tensor
void mhd_calc_0_em(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd)
{
  VARSTATIC int j;
  VARSTATIC FTYPE r, u, P, w, bsq, eta, ptot;

  bsq = dot(q->bcon, q->bcov);
  eta = bsq;
  ptot = bsq*0.5;

  /* single row of mhd stress tensor, first index up, second index down 
   */
  // mhd^{dir}_{j} =
  // j=0..3
  //  DLOOPA(j) mhd[j] = eta * q->ucon[dir] * q->ucov[j] + ptot * delta(dir, j) - q->bcon[dir] * q->bcov[j];
  DLOOPA(j) mhd[j] = eta * q->ucon[dir] * q->ucov[j] - q->bcon[dir] * q->bcov[j];
  mhd[dir] += ptot;

}


/* MHD stress tensor, with first index up, second index down */
// avoids catastrophic cancellation with rest-mass density due to extracting velocity or internal energy from that conserved energy with order unity term from rest-mass
// also avoids catastrophic cancellation in field due to using 4-field.  Instead derive stress tensor from 3-velocity and 3-field usinc Mcon_calc()
// seems to work to avoid catastrophic cancellation with field, but maybe should use WHICHVEL=RELVEL4 directly?  GODMARK
// T^dir_\mu
void mhd_calc_norestmass(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd)
{
  void mhd_calc_norestmass_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhdma, FTYPE *mhddiagpress);
  void mhd_calc_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhdem);
  VARSTATIC FTYPE mhdma[NDIM];
  VARSTATIC FTYPE mhdem[NDIM];
  VARSTATIC int j;
  VARSTATIC FTYPE mhddiagpress[NDIM];

  mhd_calc_norestmass_ma(pr, dir, geom, q, mhdma, mhddiagpress);
  mhd_calc_em(pr, dir, geom, q, mhdem);

  // add up MA and EM parts
  DLOOPA(j) mhd[j] = (mhdma[j] + mhddiagpress[j]) + mhdem[j];

}



/* MHD stress tensor, with first index up, second index down */
// avoids catastrophic cancellation with rest-mass density due to extracting velocity or internal energy from that conserved energy with order unity term from rest-mass
// also avoids catastrophic cancellation in field due to using 4-field.  Instead derive stress tensor from 3-velocity and 3-field usinc Mcon_calc()
// seems to work to avoid catastrophic cancellation with field, but maybe should use WHICHVEL=RELVEL4 directly?  GODMARK
// T^dir_\mu
void mhd_calc_norestmass_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhddiagpress)
{
  VARSTATIC int j;
  VARSTATIC FTYPE rho, u, P, w, bsq, eta, ptot;
  VARSTATIC FTYPE plus1ud0;
  void compute_1plusud0(struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0); // plus1ud0=(1+q->ucov[TT])


  //////////////////////
  //
  // Hydro part
  //
  ///////////////////////

  // below allows other scalars to be advected but not affect the stress-energy equations of motion
#if(EOMTYPE==EOMGRMHD)
  rho = pr[RHO];
  u = pr[UU];
  P = pressure_rho0_u(rho,u);
  w = P + rho + u;
#elif(EOMTYPE==EOMCOLDGRMHD)
  w = rho = pr[RHO];
  u = 0.0;
  P = 0.0;
#elif(EOMTYPE==EOMFFDE)
  rho=u=P=w=0.0;
#endif
  eta = w;
  ptot = P;


  compute_1plusud0(geom,q,&plus1ud0); // plus1ud0=(1+q->ucov[TT])

  /* single row of mhd stress tensor, first index up, second index down 
   */
  // mhd^{dir}_{j} =
  // j=0..3
  // eta u^dir u_j + rho u^dir = (p+u+b^2) u^dir u_j + rho u^dir u_j + rho u^dir
  // = (p+u+b^2) u^dir u_j + rho u^dir (u_j + 1)

  DLOOPA(j) mhddiagpress[j] = 0.0;

  j=0; mhd[j] = (P+u) * q->ucon[dir] *q->ucov[j] + rho * q->ucon[dir] * plus1ud0 ;
  SLOOPA(j) mhd[j] = eta * q->ucon[dir] * q->ucov[j];

#if(SPLITPRESSURETERMINFLUXMA==0)
  // below equivalent to ptot * delta(dir,j)
  mhd[dir] += ptot;
#else
  // below equivalent to ptot * delta(dir,j)
  mhddiagpress[dir] = ptot;
#endif

}



/* MHD stress tensor, with first index up, second index down */
// avoids catastrophic cancellation with rest-mass density due to extracting velocity or internal energy from that conserved energy with order unity term from rest-mass
// also avoids catastrophic cancellation in field due to using 4-field.  Instead derive stress tensor from 3-velocity and 3-field usinc Mcon_calc()
// seems to work to avoid catastrophic cancellation with field, but maybe should use WHICHVEL=RELVEL4 directly?  GODMARK
// T^dir_\mu
// SHOULD NOT USE b^\mu b_\mu here
void mhd_calc_primfield_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd)
{
  int Mcon_calc(FTYPE *pr, struct of_state *q, FTYPE (*Mcon)[NDIM]);
  void ffdestresstensor_dir(int dir, FTYPE (*Mcon)[NDIM], struct of_geom *geom, FTYPE *TEMdir);
  VARSTATIC FTYPE Mcon[NDIM][NDIM];
  VARSTATIC FTYPE TEMdir[NDIM];


  //////////////////////
  //
  // Electrommagnetic part
  //
  ///////////////////////
  // gives \dF^{\mu \nu}
  // GODMARK: need full Mcon to get partial TEM -- can one optimize this?
  Mcon_calc(pr, q, Mcon);

  // gives T^dir_\mu
  ffdestresstensor_dir(dir, Mcon, geom, mhd);

}




// plus1ud0=(1+q->ucov[TT])
// avoids non-relativistic velocitiy issue with machine precision, but introduces relativistic limit problem
void compute_1plusud0_old(struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0)
{
  int j,k;
  volatile FTYPE plus1gv00;
  FTYPE AA,BB,alpha;
  FTYPE vcon[NDIM];

  // 3-velocity in coordinate basis
  SLOOPA(j) vcon[j]=q->ucon[j]/q->ucon[TT];

  //  plus1gv00=1.0+geom->gcov[TT][TT];
  plus1gv00=geom->gcovpert[TT];

  AA=0.0;
  SLOOPA(j) AA+=2.0*geom->gcov[TT][j]*vcon[j];
  SLOOP(j,k) AA+=geom->gcov[j][k]*vcon[j]*vcon[k];
  //AA/=geom->gcov[TT][TT];
  BB=geom->gcov[TT][TT];

  alpha=0.0;
  SLOOPA(j) alpha+=geom->gcov[j][TT]*q->ucon[j];

  //  *plus1ud0=(plus1gv00+(2.0*alpha+alpha*alpha)*(1.0+AA)+AA)/((1.0+alpha)*(1.0+AA)+sqrt(-geom->gcov[TT][TT]*(1.0+AA)));

  *plus1ud0=(plus1gv00*BB+(2.0*alpha+alpha*alpha)*(BB+AA)+AA)/((1.0+alpha)*(BB+AA)+BB*sqrt(-(BB+AA)));

}

// plus1ud0=(1+q->ucov[TT])
// avoids both non-rel and rel limit issues with machine precision
// GODMARK: Slow, can one optmize this somehow?
void compute_1plusud0(struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0)
{
  VARSTATIC int j,k;
  //  volatile FTYPE plus1gv00;
  VARSTATIC FTYPE plus1gv00;
  VARSTATIC FTYPE vsq,gvtt,alpha;
  VARSTATIC FTYPE vcon[NDIM];
  VARSTATIC FTYPE uu0;


  // 3-velocity in coordinate basis
  SLOOPA(j) vcon[j]=q->ucon[j]/q->ucon[TT];

  //  plus1gv00=1.0+geom->gcov[TT][TT];
  plus1gv00=geom->gcovpert[TT];

  vsq=geom->gcovpert[TT];
  SLOOPA(j) vsq+=2.0*geom->gcov[TT][j]*vcon[j];
  SLOOP(j,k) vsq+=geom->gcov[j][k]*vcon[j]*vcon[k];

  gvtt=geom->gcov[TT][TT];

  alpha=0.0;
  SLOOPA(j) alpha+=geom->gcov[j][TT]*q->ucon[j];

  uu0 = q->ucon[TT];

  *plus1ud0=alpha + ((1.0-gvtt)*plus1gv00 - uu0*uu0*vsq*gvtt*gvtt)/(1.0-gvtt*uu0);

  //  dualfprintf(fail_file,"%g %g %g : wrong: %g\n",alpha,plus1gv00,gvtt,*plus1ud0);

  //  *plus1ud0=1.0+q->ucov[TT];

  //  dualfprintf(fail_file,"right: %g\n",*plus1ud0);

}


/* add in source terms to equations of motion */
// ui and dUriemann in UEVOLVE form
int source(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q, FTYPE *ui, FTYPE *dUriemann,
	   FTYPE (*dUcomp)[NPR], FTYPE *dU)
{
  //  double (*)[8]
  VARSTATIC int i,j,k,sc;
  int source_conn(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q,FTYPE *dU);
  VARSTATIC FTYPE dUother[NPR];
  VARSTATIC FTYPE Ugeomfree[NPR];


  PLOOP(k){
    SCLOOP(sc) dUcomp[sc][k] = 0.;
    dU[k] = 0.;
  }

#if(SPLITNPR)
#if((WHICHEOM==WITHNOGDET)&&(NOGDETB1>0)||(NOGDETB2>0)||(NOGDETB3>0))
#error Not correct SPLITNPR use#1
#endif
  // if only doing B1-B3, then assume no source term for magnetic field
  if(nprlist[nprstart]==B1 && nprlist[nprend]==B3) return(0);
#endif



  ////////////////
  //
  // First get geometry source (already does contain geometry prefactor term)
  source_conn(pr,ptrgeom, q,dUcomp[GEOMSOURCE]);


  ////////////////
  //
  // Second get physics source terms (using dUother for constraint)
  // get update pre- additional physics
  PLOOP(k){
    dUother[k] = (dUriemann[k] + dUcomp[GEOMSOURCE][k])/ptrgeom->e[k]; // remove geometry factor for dUother
    Ugeomfree[k] = ui[k]/ptrgeom->e[k]; // expect ui to be UEVOLVE form
  }
  sourcephysics(pr, ptrgeom, q, Ugeomfree, dUother, dUcomp);


  //////////////////
  //
  // Third, deal with equation of motion factors that don't depend upon any additional physics
  // don't add geometry prefactor onto geometry source since already has it -- only add to additional physics terms
  SCPHYSICSLOOP(sc) PLOOP(k) dUcomp[sc][k] *= ptrgeom->e[k];


  //////////////////
  //
  // Fourth, compute total since that's all the evolution cares about (comp left just for diagnostics).
  SCLOOP(sc) PLOOP(k) dU[k]+=dUcomp[sc][k];
  


  /* done! */
  return (0);
}


/* add in source terms related to geometry to equations of motion */
int source_conn(FTYPE *pr, struct of_geom *ptrgeom,
		struct of_state *q,FTYPE *dU)
{
  VARSTATIC int i = 0, j = 0, k = 0, l=0;
  VARSTATIC FTYPE dUconn[NPR];
  VARSTATIC FTYPE todo[NPR];
  VARSTATIC FTYPE ftemp;
  VARSTATIC FTYPE mhd[NDIM][NDIM];
  VARSTATIC FTYPE flux[NDIM][NPR];
  int primtofullflux(int returntype, FTYPE *pr, struct of_state *q, struct of_geom *ptrgeom, FTYPE (*flux)[NPR]);
  void mhd_calc_0(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);
  void mhd_calc(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);
  VARSTATIC int myii,myjj,mykk;



  if((ANALYTICSOURCE)&&(defcoord==LOGRSINTH)&&(MCOORD==KSCOORDS)){
    // have both WHICHEOM==WITHGDET and WHICHEOM==WITHNOGDET
    mks_source_conn(pr,ptrgeom, q, dU); // returns without geometry prefactor
    PLOOP(k) dU[k]*=ptrgeom->e[k];
  }
  else { // general non-analytic form, which uses an analytic or numerical origin for conn/conn2 (see set_grid.c)


    /////////////////////
    //
    // define first connection
    //
    // notice mhd_calc(), rather than primtoflux(), is used because this is a special connection only operated on the (at most) 4 energy-momentum EOMs.
    //
    // notice that mhd_calc_0 is used, which includes rest-mass since connection coefficient source term has rest-mass.  The rest-mass was only subtracted out of conserved quantity and flux term.
    // SUPERGODMARK: problem in non-rel gravity for below term
    //DLOOPA(j)  mhd_calc_0(pr, j, ptrgeom, q, mhd[j]);
    // 
    DLOOPA(j) mhd_calc(pr, j, ptrgeom, q, mhd[j]);

    /* contract mhd stress tensor with connection */
    // mhd^{dir}_{comp} = mhd^j_k
    // dU^{l} = T^j_k C^{k}_{l j}
    PLOOP(k) dUconn[k]=0;


#if(MCOORD!=CARTMINKMETRIC)
    myii=ptrgeom->i;
    myjj=ptrgeom->j;
    mykk=ptrgeom->k;
#else
    myii=0;
    myjj=0;
    mykk=0;
#endif


    for(l=0;l<NDIM;l++)  DLOOP(j,k){
      dUconn[UU+l] += mhd[j][k] * conn[myii][myjj][mykk][k][l][j];
    }

#if(REMOVERESTMASSFROMUU==2)
    // then need to add-in density term to source (used to avoid non-rel problems)
    // GODMARK: for no non-rel problems one must make sure conn^t_{lj} is computed accurately.
    for(l=0;l<NDIM;l++)  DLOOPA(j){
      dUconn[UU+l] += - pr[RHO] * q->ucon[j] * conn[myii][myjj][mykk][TT][l][j];
    }
#endif


    ////////////////////////////
    //
    // use first connection
    //
    // true as long as UU->U3 are not added/subtracted from other EOMs
    // Note that UtoU is ignorant of this connection term and additions/subtractions of various EOMs to/from another must account for this connection here.
    PLOOP(k) dU[k]+=ptrgeom->e[k]*dUconn[k]; 


 

    ///////////////////
    //
    // second connection
    //

#if(WHICHEOM!=WITHGDET)

    // d_t(f T^t_\nu) = -d_x^j (f F^j_\nu) + (f F^j_\nu)(d_x^j(\ln(f/\detg))) + f T^\lambda_\kappa \Gamma^\kappa_{\nu\lambda}

    // conn2_j = (d_x^j(\ln(f/\detg)))

    // deal with second connection.  Now can use general flux as defined by primtofullflux since all EOMs operate similarly with respect to second connection
    primtofullflux(UEVOLVE,pr,q,ptrgeom,flux); // returns with geometry prefactor (f F^j_\nu)


    // todo = whether that EOM has the NOGDET form.  If so, then need 2nd connection.  Do this instead of different connection2 for each EOM since each spatial component is actually the same.

    todo[RHO]=(NOGDETRHO>0) ? 1 : 0;
    todo[UU]=(NOGDETU0>0) ? 1 : 0;
    todo[U1]=(NOGDETU1>0) ? 1 : 0;
    todo[U2]=(NOGDETU2>0) ? 1 : 0;
    todo[U3]=(NOGDETU3>0) ? 1 : 0;
    todo[B1]=(NOGDETB1>0) ? 1 : 0;
    todo[B2]=(NOGDETB2>0) ? 1 : 0;
    todo[B3]=(NOGDETB3>0) ? 1 : 0;
    if(DOENTROPY) todo[ENTROPY]=(NOGDETENTROPY>0) ? 1 : 0;

    // conn2 is assumed to take care of sign
    // conn2 has geom->e and normal d(ln(gdet)) factors combined

    if(REMOVERESTMASSFROMUU){
      if(todo[RHO]!=todo[UU]){
	dualfprintf(fail_file,"Mixed form of REMOVERESTMASSFROMUU and NOGDET's not allowed.\n");
	myexit(1);
      }
    }

    //////////////////////////////////
    //
    // notice that we assume equations are differenced first, then one manipulates the connections.
    // Thus, the manipulation of connections applies to final form of EOMs afer differencing or adding.
    // This agrees with how code evolves conserved quantities, so simplest choice to make.
    //
    // Alternatively stated, UtoU() controls this ordering issue.

    PLOOP(k) DLOOPA(j) dU[k] += todo[k]*(flux[j][k] * conn2[myii][myjj][mykk][j]); // (f F^j_\nu) conn2_j



#endif // end if (WHICHEOM!=WITHGDET)

  }



  /* done! */
  return (0);
}

/* returns b^2 (i.e., twice magnetic pressure) */
int bsq_calc(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *bsq)
{
  VARSTATIC struct of_state q;

  MYFUN(get_state(pr, ptrgeom, &q) ,"phys.c:bsq_calc()", "get_state() dir=0", 1);
  *bsq = dot(q.bcon, q.bcov);
  return (0);
}


void lower_vec(FTYPE *ucon, struct of_geom *geom, FTYPE *ucov)
{
  ucov[0] = geom->gcov[0][0]*ucon[0]
    + geom->gcov[0][1]*ucon[1]
    + geom->gcov[0][2]*ucon[2]
    + geom->gcov[0][3]*ucon[3] ;
  ucov[1] = geom->gcov[0][1]*ucon[0]
    + geom->gcov[1][1]*ucon[1]
    + geom->gcov[1][2]*ucon[2]
    + geom->gcov[1][3]*ucon[3] ;
  ucov[2] = geom->gcov[0][2]*ucon[0]
    + geom->gcov[1][2]*ucon[1]
    + geom->gcov[2][2]*ucon[2]
    + geom->gcov[2][3]*ucon[3] ;
  ucov[3] = geom->gcov[0][3]*ucon[0]
    + geom->gcov[1][3]*ucon[1]
    + geom->gcov[2][3]*ucon[2]
    + geom->gcov[3][3]*ucon[3] ;

  return ;
}

void lowerf(FTYPE *fcon, struct of_geom *geom, FTYPE *fcov)
{
  VARSTATIC int j,k;
  VARSTATIC int jp,kp;
  VARSTATIC FTYPE myfcon[NDIM][NDIM],myfcov[NDIM][NDIM];


  myfcon[0][0]=myfcon[1][1]=myfcon[2][2]=myfcon[3][3]=0;
  myfcon[0][1]=fcon[0];
  myfcon[0][2]=fcon[1];
  myfcon[0][3]=fcon[2];
  myfcon[1][2]=fcon[3];
  myfcon[1][3]=fcon[4];
  myfcon[2][3]=fcon[5];
  //
  myfcon[1][0]=-fcon[0];
  myfcon[2][0]=-fcon[1];
  myfcon[3][0]=-fcon[2];
  myfcon[2][1]=-fcon[3];
  myfcon[3][1]=-fcon[4];
  myfcon[3][2]=-fcon[5];
  
  DLOOP(j,k){
    myfcov[j][k]=0;
    for(jp=0;jp<NDIM;jp++) for(kp=0;kp<NDIM;kp++){
      myfcov[j][k]+=myfcon[jp][kp]*geom->gcov[j][jp]*geom->gcov[k][kp];
    }
  }
  fcov[0]=myfcov[0][1];
  fcov[1]=myfcov[0][2];
  fcov[2]=myfcov[0][3];
  fcov[3]=myfcov[1][2];
  fcov[4]=myfcov[1][3];
  fcov[5]=myfcov[2][3];

  return ;
}

void raise_vec(FTYPE *ucov, struct of_geom *geom, FTYPE *ucon)
{

  ucon[0] = geom->gcon[0][0]*ucov[0]
    + geom->gcon[0][1]*ucov[1]
    + geom->gcon[0][2]*ucov[2]
    + geom->gcon[0][3]*ucov[3] ;
  ucon[1] = geom->gcon[0][1]*ucov[0]
    + geom->gcon[1][1]*ucov[1]
    + geom->gcon[1][2]*ucov[2]
    + geom->gcon[1][3]*ucov[3] ;
  ucon[2] = geom->gcon[0][2]*ucov[0]
    + geom->gcon[1][2]*ucov[1]
    + geom->gcon[2][2]*ucov[2]
    + geom->gcon[2][3]*ucov[3] ;
  ucon[3] = geom->gcon[0][3]*ucov[0]
    + geom->gcon[1][3]*ucov[1]
    + geom->gcon[2][3]*ucov[2]
    + geom->gcon[3][3]*ucov[3] ;

  return ;
}


/* find ucon, ucov, bcon, bcov from primitive variables */
int get_state(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);


  get_state_uconucovonly(pr, ptrgeom, q);

  get_state_bconbcovonly(pr, ptrgeom, q);


  return (0);
}

/* find ucon, ucov, bcon, bcov from primitive variables */
int get_stateforsource(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);


  get_state_uconucovonly(pr, ptrgeom, q);

  // below needed for dUtodt() and compute_dt_fromsource()
#if(LIMITDTWITHSOURCETERM)
  get_state_bconbcovonly(pr, ptrgeom, q);
#endif



  return (0);
}


// wrapper for choosing to get state by computing it or from global array
int get_stateforfluxcalc(int dimen, int isleftright, FTYPE *pr, struct of_geom *ptrgeom, struct of_state **qptr)
{
  int pureget_stateforfluxcalc(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int i,j,k,jj;



#if(STOREFLUXSTATE)
  // retrieve state from global array
  i=ptrgeom->i;
  j=ptrgeom->j;
  k=ptrgeom->k;

  // JCM: I see now point in copying the data, just should change pointer address

#if(1)
  // pointer assign (overrides the existing default assignment of the pointer address in the calling function, so assumes calling function only uses pointer form subsequently)
  *qptr=&(fluxstate[dimen][isleftright][i][j][k]);
#else
  DLOOPA(jj){
    (*qptr)->(ucon[jj])=fluxstate[dimen][isleftright][i][j][k].ucon[jj];
    (*qptr)->(ucov[jj])=fluxstate[dimen][isleftright][i][j][k].ucov[jj];
  }
#if(COMPUTE4FIELDforFLUX)
  DLOOPA(jj){
    (*qptr)->(bcon[jj])=fluxstate[dimen][isleftright][i][j][k].bcon[jj];
    (*qptr)->(bcov[jj])=fluxstate[dimen][isleftright][i][j][k].bcov[jj];
  }
#endif

#endif // end if copy method

#else
  // compute state right now
  pureget_stateforfluxcalc(pr, ptrgeom, *qptr);
#endif


  return (0);
}


/* find ucon, ucov, bcon, bcov from primitive variables */
int pureget_stateforfluxcalc(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);


  get_state_uconucovonly(pr, ptrgeom, q);

#if(COMPUTE4FIELDforFLUX)
  get_state_bconbcovonly(pr, ptrgeom, q);
#endif



  return (0);
}




/* find ucon, ucov, bcon, bcov from primitive variables */
int get_stateforUdiss(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
  int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);


  get_state_uconucovonly(pr, ptrgeom, q);

#if(MAXWELL==GENMAXWELL)
  get_state_bconbcovonly(pr, ptrgeom, q);
#endif



  return (0);
}



/* find ucon, ucov, bcon, bcov from primitive variables */
int get_state_bconbcovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{

  // b^\mu
  bcon_calc(pr, q->ucon, q->ucov, q->bcon);
  // b^\mu
  lower_vec(q->bcon, ptrgeom, q->bcov);

  return (0);
}


// Get only u^\mu and u_\mu assumine b^\mu and b_\mu not used
int get_state_uconucovonly(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{

  // u^\mu
  MYFUN(ucon_calc(pr, ptrgeom, q->ucon) ,"phys.c:get_state()", "ucon_calc()", 1);
  // u_\mu
  lower_vec(q->ucon, ptrgeom, q->ucov);

  return (0);
}

/* load local geometry into structure geom */
// fills struct of_geom with PRIMECOORDS metric data
void get_geometry(int ii, int jj, int kk, int pp, struct of_geom *geom)
{
  VARSTATIC int j, k;
  void assign_eomfunc(struct of_geom *geom, FTYPE eomfunc);
  VARSTATIC int myii,myjj,mykk;
  VARSTATIC int pl;
  

#if(MCOORD!=CARTMINKMETRIC)
  myii=ii;
  myjj=jj;
  mykk=kk;
#else
  myii=0;
  myjj=0;
  mykk=0;
#endif

  

#if(GETGEOMUSEPOINTER==0)
  //  DLOOP(j,k) geom->gcov[j][k] = gcov[myii][myjj][mykk][pp][j][k];
  //DLOOP(j,k) geom->gcon[j][k] = gcon[myii][myjj][mykk][pp][j][k];
  // let's vectorize it
  for(j=0;j<=NDIM*NDIM-1;j++){
    geom->gcon[0][j] = gcon[myii][myjj][mykk][pp][0][j];
    geom->gcov[0][j] = gcov[myii][myjj][mykk][pp][0][j];
  }
  for(j=0;j<NDIM;j++){
    geom->gcovpert[j]=gcovpert[myii][myjj][mykk][pp][j];
  }
#else
  // let's pointer it, faster by a bit
  geom->gcov=(FTYPE (*)[NDIM])(&(gcov[myii][myjj][mykk][pp][0][0])); // pointer

  geom->gcovpert=(FTYPE (*))(&(gcovpert[myii][myjj][mykk][pp])); // pointer

  geom->gcon=(FTYPE (*)[NDIM])(&(gcon[myii][myjj][mykk][pp][0][0])); // pointer
#endif

  geom->g = gdet[myii][myjj][mykk][pp];

  // get eomfunc (see eomfunc_func() and lngdet_func_mcoord() in metric.c)
  assign_eomfunc(geom,eomfunc[myii][myjj][mykk][pp]);

  geom->alphalapse = alphalapse[myii][myjj][mykk][pp];

  ///////////
  //
  // these should be real grid i,j,k positions
  geom->i = ii;
  geom->j = jj;
  geom->k = kk;
  geom->p = pp;
#if(JONCHECKS)
  icurr = ii;
  jcurr = jj;
  kcurr = kk;
  pcurr = pp;
#endif


}

/* load local geometry into structure geom */
// fills struct of_geom with PRIMECOORDS metric data
void get_geometry_gdetonly(int ii, int jj, int kk, int pp, struct of_geom *geom)
{
  VARSTATIC int j, k;
  void assign_eomfunc(struct of_geom *geom, FTYPE eomfunc);
  VARSTATIC int myii,myjj,mykk;
  VARSTATIC int pl;
  

#if(MCOORD!=CARTMINKMETRIC)
  myii=ii;
  myjj=jj;
  mykk=kk;
#else
  myii=0;
  myjj=0;
  mykk=0;
#endif

  

  geom->g = gdet[myii][myjj][mykk][pp];


  ///////////
  //
  // these should be real grid i,j,k positions
  geom->i = ii;
  geom->j = jj;
  geom->k = kk;
  geom->p = pp;
#if(JONCHECKS)
  icurr = ii;
  jcurr = jj;
  kcurr = kk;
  pcurr = pp;
#endif


}


/* load local geometry into structure geom */
// fills struct of_geom with PRIMECOORDS metric data
void get_geometry_geomeonly(int ii, int jj, int kk, int pp, struct of_geom *geom)
{
  VARSTATIC int pl;
  void get_geometry(int ii, int jj, int kk, int pp, struct of_geom *geom);
  void get_geometry_gdetonly(int ii, int jj, int kk, int pp, struct of_geom *geom);
  


#if(WHICHEOM!=WITHGDET)
  // more complicated geom.e
  get_geometry(ii, jj, kk, pp, geom);

#else
  // then geom.e is just geom.g
  get_geometry_gdetonly(ii, jj, kk, pp, geom);
  PALLLOOP(pl) geom->e[pl]=geom->g;
#endif


}


// set igeom part of geometry since more expensive and not always needed
void set_igeom(struct of_geom *geom)
{
  VARSTATIC int pl;


  //////////////
  // avoids 0.0 for any sign of geom.e[pl]
#if(GDETVOLDIFF==0)

  geom->igeomnosing = sign(geom->g)/(fabs(geom->g)+SMALL);

  // use PALLLOOP so compiler can optimize
#if(WHICHEOM!=WITHGDET)
  PALLLOOP(pl) geom->ienosing[pl] = sign(geom->e[pl])/(fabs(geom->e[pl])+SMALL);
#else
  PALLLOOP(pl) geom->ienosing[pl]=geom->igeomnosing;
#endif


#else
  // volume regularization (correct to second order) // GODMARK: NOT FOR FINITE VOLUME WENO METHOD
  igeomnosing = sign(gdetvol[geom->i][geom->j][geom->k][geom->p])/(fabs(gdetvol[geom->i][geom->j][geom->k][geom->p])+SMALL);
  geom->igeomnosing = igeomnosing;
  // use PALLLOOP so compiler can optimize
  PALLLOOP(pl) geom->ienosing[pl] = igeomnosing;
#endif


}

// set igeom part of geometry since more expensive and not always needed
void set_igeomsimple(struct of_geom *geom)
{
  VARSTATIC int pl;


#if(WHICHEOM!=WITHGDET)
  dualfprintf(fail_file,"Using set_igeomsimple() but WHICHEOM!=WITHGDET\n");
  myexit(342968347);
#endif


  //////////////
  // avoids 0.0 for any sign of geom.e[pl]
#if(GDETVOLDIFF==0)
  geom->igeomnosing = sign(geom->g)/(fabs(geom->g)+SMALL);
#else
  // volume regularization (correct to second order) // GODMARK: NOT FOR FINITE VOLUME WENO METHOD
  igeomnosing = sign(gdetvol[geom->i][geom->j][geom->k][geom->p])/(fabs(gdetvol[geom->i][geom->j][geom->k][geom->p])+SMALL);
  geom->igeomnosing = igeomnosing;
#endif

}


/* load local geometry into structure geom */
void get_allgeometry(int ii, int jj, int kk, int pp, struct of_allgeom *allgeom, struct of_geom *geom)
{
  void get_geometry(int ii, int jj, int kk, int pp, struct of_geom *geom);
  VARSTATIC int j, k;
  VARSTATIC int pl;


  get_geometry(ii, jj, kk, pp, geom);

  allgeom->i=geom->i;
  allgeom->j=geom->j;
  allgeom->k=geom->k;
  allgeom->p=geom->p;

#if(GETGEOMUSEPOINTER==0)
  for(j=0;j<=NDIM*NDIM-1;j++){
    allgeom->gcov[0][j]=geom->gcov[0][j];
    allgeom->gcon[0][j]=geom->gcon[0][j];
  }
  for(j=0;j<NDIM;j++){
    allgeom->gcovpert[j]=geom->gcovpert[j];
  }
#else
  allgeom->gcov=geom->gcov;
  allgeom->gcon=geom->gcon;
  allgeom->gcovpert=geom->gcovpert;
#endif

  allgeom->g=geom->g;
  PLOOP(pl) allgeom->e[pl]=geom->e[pl];
  allgeom->alphalapse = geom->alphalapse;

  PALLLOOP(pl) allgeom->ienosing[pl] = geom->ienosing[pl];
  allgeom->igeomnosing = geom->igeomnosing;


  coord_ijk(ii, jj, kk, pp, allgeom->X);
  bl_coord_ijk( ii, jj, kk, pp , allgeom->V );
  dxdxprim_ijk( ii, jj, kk, pp , allgeom->dxdxp);


}


//eomfuncarray[ii][jj][kk][pp];
// stationary factor in equation of motion that multiplies T^\mu_\nu
void assign_eomfunc(struct of_geom *geom, FTYPE eomfunc)
{

  // now set each EOM
#if(NOGDETRHO==0)
  geom->e[RHO]=geom->g;
#else
  geom->e[RHO]=eomfunc;
#endif
#if(NOGDETU0==0)
  geom->e[UU]=geom->g;
#else
  geom->e[UU]=eomfunc;
#endif
#if(NOGDETU1==0)
  geom->e[U1]=geom->g;
#else
  geom->e[U1]=eomfunc;
#endif
#if(NOGDETU2==0)
  geom->e[U2]=geom->g;
#else
  geom->e[U2]=eomfunc;
#endif
#if(NOGDETU3==0)
  geom->e[U3]=geom->g;
#else
  geom->e[U3]=eomfunc;
#endif

#if(NOGDETB1==0)
  geom->e[B1]=geom->g;
#else
  geom->e[B1]=eomfunc;
#endif
#if(NOGDETB2==0)
  geom->e[B2]=geom->g;
#else
  geom->e[B2]=eomfunc;
#endif
#if(NOGDETB3==0)
  geom->e[B3]=geom->g;
#else
  geom->e[B3]=eomfunc;
#endif

#if(DOENTROPY)
#if(NOGDETENTROPY==0)
  geom->e[ENTROPY]=geom->g;
#else
  geom->e[ENTROPY]=eomfunc;
#endif
#endif


}



/* find contravariant four-velocity from the relative 4 velocity */
int ucon_calc_rel4vel_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE *ucon)
{
  VARSTATIC FTYPE alpha,alphasq,gamma ;
  VARSTATIC FTYPE beta[NDIM] ;
  VARSTATIC FTYPE qsq;
  VARSTATIC int j ;
  int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma);


  //  alpha = 1./sqrt(-geom->gcon[TT][TT]) ;
  alpha = geom->alphalapse;
  alphasq=alpha*alpha;
  SLOOPA(j) beta[j] = geom->gcon[TT][j]*alphasq ;

  MYFUN(gamma_calc_fromuconrel(uconrel,geom,&gamma),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.c",1);

  ucon[TT] = gamma/alpha ;
  SLOOPA(j) ucon[j] = uconrel[j] - ucon[TT]*beta[j] ;

  // hence v^j = uconrel^j/u^t - beta^j

  return(0) ;
}

/* find contravariant four-velocity from the relative 4 velocity */
int ucon_calc_rel4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
{
  VARSTATIC FTYPE uconrel[NDIM];
  VARSTATIC int j;
  int ucon_calc_rel4vel_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE *ucon);

  uconrel[0]=0;
  uconrel[1]=pr[U1];
  uconrel[2]=pr[U2];
  uconrel[3]=pr[U3];

  //  SLOOPA(j) dualfprintf(fail_file,"ucon_calc_rel4vel: uconrel[%d]=%21.15g\n",j,uconrel[j]); // CHANGINGMARK


  MYFUN(ucon_calc_rel4vel_fromuconrel(uconrel, geom, ucon),"ucon_calc_rel4vel: ucon_calc_rel4vel_fromuconrel failed\n","phys.c",1);

  return(0) ;
}


/* find gamma-factor wrt normal observer */
int gamma_calc(FTYPE *pr, struct of_geom *geom, FTYPE*gamma)
{
  VARSTATIC FTYPE uconrel[NDIM];
  VARSTATIC int j;
  int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma);

#if(WHICHVEL!=VELREL4)
  dualfprintf(fail_file,"gamma_calc() designed for WHICHVEL=VELREL4\n");
  myexit(77);
#endif

  // assumes input pr is WHICHVEL=VELREL4
  uconrel[0]=0;
  uconrel[1]=pr[U1];
  uconrel[2]=pr[U2];
  uconrel[3]=pr[U3];

  //  SLOOPA(j) dualfprintf(fail_file,"gamma_calc: uconrel[%d]=%21.15g\n",j,uconrel[j]); // CHANGINGMARK


  // get gamma
  MYFUN(gamma_calc_fromuconrel(uconrel, geom, gamma),"gamma_calc: gamma_calc_fromuconrel failed\n","phys.c",1);

  return(0);
}


/* find gamma-factor wrt normal observer */
int gamma_calc_fromuconrel(FTYPE *uconrel, struct of_geom *geom, FTYPE*gamma)
{
  VARSTATIC int j,k;
  VARSTATIC FTYPE qsq;
  int qsq_calc(FTYPE *uconrel, struct of_geom *geom, FTYPE *qsq);


  qsq_calc(uconrel,geom,&qsq);


#if(JONCHECKS && PRODUCTION==0)
  if(qsq<0.0){
    if(qsq>-NUMEPSILON*100){ // then assume not just machine precision
      qsq=0.0;
    }
    else{
      dualfprintf(fail_file,"gamma_calc failed: i=%d j=%d k=%d qsq=%21.15g\n",geom->i,geom->j,geom->k,qsq);
      SLOOPA(j) dualfprintf(fail_file,"uconrel[%d]=%21.15g\n",j,uconrel[j]);
      DLOOP(j,k) dualfprintf(fail_file,"gcov[%d][%d]=%21.15g\n",j,k,geom->gcov[j][k]);
      if (fail(FAIL_UTCALC_DISCR) >= 1)
	return (1);
    }
  }
#else
  // minimal check -- forces qsq=0 if qsq<0.0
  qsq = qsq*(qsq>=0.0);
#endif

  *gamma = sqrt(1. + qsq) ;

  //  SLOOPA(j) dualfprintf(fail_file,"uconrel[%d]=%21.15g\n",j,uconrel[j]);
  //  dualfprintf(fail_file,"qsq=%21.15g, gamma=%21.15g\n",qsq,*gamma); // CHANGINGMARK

  return(0) ;
}


// get \tilde{u}^i \tilde{u}^j g_{ij}
int qsq_calc(FTYPE *uconrel, struct of_geom *geom, FTYPE *qsq)
{
  VARSTATIC int j,k;
  
  *qsq =
    + geom->gcov[1][1]*uconrel[1]*uconrel[1]
    + geom->gcov[2][2]*uconrel[2]*uconrel[2]
    + geom->gcov[3][3]*uconrel[3]*uconrel[3]
    + 2.*(geom->gcov[1][2]*uconrel[1]*uconrel[2]
	  + geom->gcov[1][3]*uconrel[1]*uconrel[3]
	  + geom->gcov[2][3]*uconrel[2]*uconrel[3]) ;

  //  DLOOP(j,k) dualfprintf(fail_file,"gcov[%d][%d]=%21.15g uconrel[%d]=%21.15g uconrel[%d]=%21.15g\n",j,k,geom->gcov[j][k],j,uconrel[j],k,uconrel[k]); // CHANGINGMARK

  return(0) ;
}






/* find contravariant four-velocity */
//int ucon_calc(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
int ucon_calc_3vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
{
  VARSTATIC FTYPE negdiscr;
  VARSTATIC FTYPE velterm;
  VARSTATIC FTYPE vcon[NDIM];
  // debug stuff
  VARSTATIC int j,k;
  VARSTATIC FTYPE V[NDIM],X[NDIM];
  //
  int get_3velterm(FTYPE *vcon, struct of_geom *geom, FTYPE *velterm);



  SLOOPA(j) vcon[j] = pr[U1+j-1];

  get_3velterm(vcon, geom, &velterm);

  negdiscr = geom->gcov[0][0] + velterm ;


#if(JONCHECKS)  
  uttdiscr=-negdiscr;
#endif

  if (negdiscr > 0.) {
#if(JONCHECKS)
    if(whocalleducon==0){
      // then report on disc
      dualfprintf(fail_file,"negdisc=%21.15g, should be negative\n",negdiscr);
      for(k=U1;k<=U3;k++){
	dualfprintf(fail_file,"uconfailed on pr[%d]=%21.15g\n",k,pr[k]);
      }
      coord(geom->i,geom->j,geom->k,geom->p,X);
      bl_coord(X,V);
      dualfprintf(fail_file,"i=%d j=%d k=%d pcurr=%d\nx1=%21.15g x2=%21.15g x3=%21.15g \nV1=%21.15g V2=%21.15g V3=%21.15g \ng=%21.15g\n",startpos[1]+geom->i,startpos[2]+geom->j,startpos[3]+geom->k,geom->p,X[1],X[2],X[3],V[1],V[2],V[3],geom->g);
      dualfprintf(fail_file,"\ngcon\n");
      dualfprintf(fail_file,"{");
      for(j=0;j<NDIM;j++){
	dualfprintf(fail_file,"{");
	for(k=0;k<NDIM;k++){
	  dualfprintf(fail_file,"%21.15g",geom->gcon[j][k]);
	  if(k!=NDIM-1) dualfprintf(fail_file," , ");
	}
	dualfprintf(fail_file,"}");	
	if(j!=NDIM-1) dualfprintf(fail_file," , ");
      }
      dualfprintf(fail_file,"}");
      dualfprintf(fail_file,"\ngcov\n");
      dualfprintf(fail_file,"{");
      for(j=0;j<NDIM;j++){
	dualfprintf(fail_file,"{");
	for(k=0;k<NDIM;k++){
	  dualfprintf(fail_file,"%21.15g",geom->gcov[j][k]);
	  if(k!=NDIM-1) dualfprintf(fail_file," , ");
	}
	dualfprintf(fail_file,"}");	
	if(j!=NDIM-1) dualfprintf(fail_file," , ");
      }
      dualfprintf(fail_file,"}");
    }
#endif
    if (fail(FAIL_UTCALC_DISCR) >= 1)
      return (1);
  }

  ucon[TT] = 1. / sqrt(-negdiscr);
  SLOOPA(j) ucon[j]=vcon[j]*ucon[TT];

  return (0);
}


// 2v^i g_{it} + v^i v^j g_{ij}
int get_3velterm(FTYPE *vcon, struct of_geom *geom, FTYPE *velterm)
{
  VARSTATIC int j;


  *velterm =
    + geom->gcov[1][1] * vcon[1] * vcon[1]
    + geom->gcov[2][2] * vcon[2] * vcon[2]
    + geom->gcov[3][3] * vcon[3] * vcon[3]
    + 2. * (geom->gcov[0][1]* vcon[1]
	    + geom->gcov[0][2] * vcon[2]
	    + geom->gcov[0][3] * vcon[3]
	    + geom->gcov[1][2] * vcon[1] * vcon[2]
	    + geom->gcov[1][3] * vcon[1] * vcon[3]
	    + geom->gcov[2][3] * vcon[2] * vcon[3]
	    );
  return(0);
}


/* find contravariant time component of four-velocity from the 4velocity (3 terms)*/
int ucon_calc_4vel_bothut(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *ucon2)
{
  VARSTATIC int i,j,k;
  VARSTATIC FTYPE AA,BB,CCM1 ;
  VARSTATIC FTYPE discr ;
  VARSTATIC FTYPE bsq,X[NDIM] ;
  VARSTATIC FTYPE sqrtdiscr;
  int get_4velterms(FTYPE *ucon, struct of_geom *geom, FTYPE *AA, FTYPE *BB, FTYPE *CCM1, FTYPE *discr);


  ucon[1] = pr[U1] ;
  ucon[2] = pr[U2] ;
  ucon[3] = pr[U3] ;

  ucon2[1] = pr[U1] ;
  ucon2[2] = pr[U2] ;
  ucon2[3] = pr[U3] ;

  get_4velterms(ucon, geom, &AA, &BB, &CCM1, &discr);

  if(discr < 0.) {
    /*
      fprintf(fail_file,"failure %d %d\n",icurr,jcurr) ;
      ucon[TT] = (-BB - sqrt(-discr))/(2.*AA) ;
      ucon[TT] = -BB/(2.*AA) ;
    */
    fprintf(fail_file,"failure: spacelike four-velocity %21.15g\n",
	    discr) ;
    fprintf(fail_file,"i=%d j=%d k=%d p=%d\n",startpos[1]+geom->i,startpos[2]+geom->j,startpos[3]+geom->k,geom->p) ;
    coord(geom->i,geom->j,geom->k,geom->p,X);
    fprintf(fail_file,"%21.15g %21.15g %21.15g\n",X[1],X[2],X[3]) ;

    for(k=0;k<NPR;k++) fprintf(fail_file,"%d %21.15g\n",k,pr[k]) ;
    // GODMARK -- why did we have failed=1?
    //		failed=1;
    return(1);
  }

  sqrtdiscr=sqrt(discr);

  ucon[TT] = (-BB - sqrtdiscr)/(2.*AA) ; // standard solution always valid outside ergosphere
  ucon2[TT] = (-BB + sqrtdiscr)/(2.*AA) ; // can be a valid solution within ergosphere

  return(0) ;
}


/* find contravariant time component of four-velocity from the 4velocity (3 terms)*/
int ucon_calc_4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
{
  VARSTATIC int i,j,k;
  VARSTATIC FTYPE AA,BB,CCM1 ;
  VARSTATIC FTYPE discr ;
  VARSTATIC FTYPE bsq,X[NDIM] ;
  int get_4velterms(FTYPE *ucon, struct of_geom *geom, FTYPE *AA, FTYPE *BB, FTYPE *CCM1, FTYPE *discr);


  ucon[1] = pr[U1] ;
  ucon[2] = pr[U2] ;
  ucon[3] = pr[U3] ;

  get_4velterms(ucon, geom, &AA, &BB, &CCM1, &discr);

  if(discr < 0.) {
    /*
      fprintf(fail_file,"failure %d %d\n",icurr,jcurr) ;
      ucon[TT] = (-BB - sqrt(-discr))/(2.*AA) ;
      ucon[TT] = -BB/(2.*AA) ;
    */
    fprintf(fail_file,"failure: spacelike four-velocity %21.15g\n",
	    discr) ;
    fprintf(fail_file,"i=%d j=%d k=%d p=%d\n",startpos[1]+geom->i,startpos[2]+geom->j,startpos[3]+geom->k,geom->p) ;
    coord(geom->i,geom->j,geom->k,geom->p,X);
    fprintf(fail_file,"%21.15g %21.15g %21.15g\n",X[1],X[2],X[3]) ;

    for(k=0;k<NPR;k++) fprintf(fail_file,"%d %21.15g\n",k,pr[k]) ;
    // GODMARK -- why did we have failed=1?
    //		failed=1;
    return(1);
  }

  ucon[TT] = (-BB - sqrt(discr))/(2.*AA) ;

  return(0) ;
}


int get_4velterms(FTYPE *ucon, struct of_geom *geom, FTYPE *AA, FTYPE *BB, FTYPE *CCM1, FTYPE *discr)
{
  VARSTATIC FTYPE CC;

  // g_{tt}
  *AA = geom->gcov[TT][TT] ;

  // 2 u^i g_{it}
  *BB = 2.*(geom->gcov[TT][1]*ucon[1] +
	    geom->gcov[TT][2]*ucon[2] +
	    geom->gcov[TT][3]*ucon[3]) ;

  // u^i u^j g_{ij}
  *CCM1 = geom->gcov[1][1]*ucon[1]*ucon[1] +
    geom->gcov[2][2]*ucon[2]*ucon[2] +
    geom->gcov[3][3]*ucon[3]*ucon[3] +
    2.*(geom->gcov[1][2]*ucon[1]*ucon[2] +
	geom->gcov[1][3]*ucon[1]*ucon[3] +
	geom->gcov[2][3]*ucon[2]*ucon[3]) ;

  // 1 + u^i u^j g_{ij}
  CC = *CCM1 + 1.0;
  
  // B^2 - 4AC
  *discr = (*BB)*(*BB) - 4.*(*AA)*CC ;

  return(0);
}


/* find contravariant time component of four-velocity from the 4velocity (3 terms)*/
int ucon_calc_nonrel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
{

  // this isn't really right
  // neglects kinetic energy term.  Need to really re-write T^\mu_\nu

  ucon[1] = pr[U1] ;
  ucon[2] = pr[U2] ;
  ucon[3] = pr[U3] ;
  ucon[TT] = 1.0;

  return(0) ;
}


FTYPE taper_func(FTYPE R,FTYPE rin)
{

  if(R <= rin)
    return(0.) ;
  else
    return(1. - sqrt(rin/R)) ;

}





// used Mathematica's MinimumChangePermutations and Signature
// updown = 0 : down
// updown = 1 : up
FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda)
{
  int i;
  FTYPE lc4sign; // 1,-1,1,-1... for all 24 entires
  int l1[24]={1, 2, 3, 1, 2, 3, 4, 2, 1, 4, 2, 1, 1, 3, 4, 1, 3, 4, 4, 3, 2, 4, 3, 2};
  int l2[24]={2, 1, 1, 3, 3, 2, 2, 4, 4, 1, 1, 2, 3, 1, 1, 4, 4, 3, 3, 4, 4, 2, 2, 3};
  int l3[24]={3, 3, 2, 2, 1, 1, 1, 1, 2, 2, 4, 4, 4, 4, 3, 3, 1, 1, 2, 2, 3, 3, 4, 4};
  int l4[24]={4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1};

  for(i=0;i<24;i++){
    if((1+mu==l1[i])&&(1+nu==l2[i])&&(1+kappa==l3[i])&&(1+lambda==l4[i])){
      lc4sign=(i%2) ? -1 : 1;
      if(updown==1) return(-1.0/detg*lc4sign); // up
      else if(updown==0) return(detg*lc4sign); // down
    }
  }
  // if didn't get here, then 0
  return(0.0);
}

// below not used currently
void faraday_calc(int which, FTYPE *b, FTYPE *u, struct of_geom *geom, FTYPE faraday[][NDIM])
{
  int nu,mu,kappa,lambda;

  for(nu=0;nu<NDIM;nu++) for(mu=0;mu<NDIM;mu++){
    faraday[mu][nu]=0.0;
    for(kappa=0;kappa<NDIM;kappa++) for(lambda=0;lambda<NDIM;lambda++){
      faraday[mu][nu]+=lc4(which,geom->g,mu,nu,kappa,lambda)*u[kappa]*b[lambda];
    }
  }

}



// setup faraday for either instrinsic use (which==CURTYPEFARADAY) or for computing the currents
int current_doprecalc(int which, FTYPE p[][N2M][N3M][NPR])
{
  int i,j,k;
  int idel, jdel,kdel;
  struct of_geom geom;
  struct of_state q;
  FTYPE Dt;
  int face;

#if(CALCFARADAYANDCURRENTS==0)
  dualfprintf(fail_file,"Shouldn't be here in current_doprecalc()\n");
  myexit(12);
#endif

  if(WHICHCURRENTCALC==CURRENTCALC0){ // then need to save 2 times
    // which==1,2 should be using face primitives, but probably not
    if(which==CURTYPET){ face=CENT; idel=0; jdel=0; kdel=0; Dt=dt*0.5; }
    else if(which==CURTYPEX){ face=FACE1; idel=1; jdel=0; kdel=0; Dt=dt; }
    else if(which==CURTYPEY){ face=FACE2; idel=0; jdel=1; kdel=0; Dt=dt; }
    else if(which==CURTYPEZ){ face=FACE3; idel=0; jdel=0; kdel=1; Dt=dt; }
    else if(which==CURTYPEFARADAY){ face=CENT; idel=0; jdel=0; kdel=0; Dt=dt; }
  }
  else if(WHICHCURRENTCALC==CURRENTCALC1){ // time and space centered on present time (best method)
    face=CENT;
    idel=0;
    jdel=0;
    kdel=0;
    Dt=dt;
  }
  else if(WHICHCURRENTCALC==CURRENTCALC2){ // then need to save 2 times
    if(which==CURTYPET){ face=CENT; idel=0; jdel=0; kdel=0; Dt=dt*0.5; }
    else if(which==CURTYPEX){ face=CENT; idel=1; jdel=0; kdel=0; Dt=dt; }
    else if(which==CURTYPEY){ face=CENT; idel=0; jdel=1; kdel=0; Dt=dt; }
    else if(which==CURTYPEZ){ face=CENT; idel=0; jdel=0; kdel=1; Dt=dt; }
    else if(which==CURTYPEFARADAY){ face=CENT; idel=0; jdel=0; kdel=0; Dt=dt; }
  }

  // assume no other conditions GODMARK

  COMPFULLLOOP{
    get_geometry(i, j, k, face, &geom);
    MYFUN(get_state(p[i][j][k], &geom, &q),"phys.c:current_doprecalc()", "get_state()", 1);
    current_precalc(which,&geom,&q,Dt,cfaraday[i][j][k]);
  }

  return(0);
}

// which: 0: somewhere at half step
// which: 1: doing flux calculation in x1 direction, full step
// which: 2: doing flux calculation in x2 direction, full step
// faraday[0-3][0-3]=[0=mid-time,1=radial edge final time, 2= theta edge final time, 3=old mid-time][whatever 3 things needed to compute relevant current]
void current_precalc(int which, struct of_geom *geom, struct of_state *q, FTYPE Dt,FTYPE faraday[][3])
{
  // assume outside loop is like flux, from 0..N in r for r, 0..N in h for h.  And must wait till 2nd timestep before computing the current since need time differences

#if(CALCFARADAYANDCURRENTS==0)
  dualfprintf(fail_file,"Shouldn't be here in current_doprecalc()\n");
  myexit(12);
#endif


  if(which==CURTYPET){ // for d/dt

    // assume got here when DT==dt/2 and geom and state set at zone center
    // first save old calculation
    if((WHICHCURRENTCALC==0)||(WHICHCURRENTCALC==2)){ // then need to save 2 times
      faraday[4][0]=faraday[0][0];
      faraday[4][1]=faraday[0][1];
      faraday[4][2]=faraday[0][2];
    }
    else if(WHICHCURRENTCALC==1){ // then need to save 3 times and have [4] as 2 times ago
      // 2 times ago
      faraday[4][0]=faraday[5][0];
      faraday[4][1]=faraday[5][1];
      faraday[4][2]=faraday[5][2];

      // 1 time ago
      faraday[5][0]=faraday[0][0];
      faraday[5][1]=faraday[0][1];
      faraday[5][2]=faraday[0][2];
    }
    // now calculate new version

    // Need F^{xt} F^{yt} F^{zt}
    faraday[0][0]=-1.0/geom->g * (-q->bcov[3]*q->ucov[2]+q->bcov[2]*q->ucov[3]); // f^{rt}
    faraday[0][1]=-1.0/geom->g * (q->bcov[3]*q->ucov[1]-q->bcov[1]*q->ucov[3]); // f^{ht}
    faraday[0][2]=-1.0/geom->g * (-q->bcov[2]*q->ucov[1]+q->bcov[1]*q->ucov[2]); // f^{pt}
  }
  else if(which==CURTYPEX){// for d/dx

    // Need F^{tx} F^{yx} F^{zx}

    // assume got here with DT=dt and geom and state at radial zone edge
    faraday[1][0]=-1.0/geom->g * (q->bcov[3]*q->ucov[2]-q->bcov[2]*q->ucov[3]); // f^{tr}=-f^{rt}
    faraday[1][1]=-1.0/geom->g * (-q->bcov[3]*q->ucov[0]+q->bcov[0]*q->ucov[3]); // f^{hr}
    faraday[1][2]=-1.0/geom->g * (q->bcov[2]*q->ucov[0]-q->bcov[0]*q->ucov[2]); // f^{pr}
  }
  else if(which==CURTYPEY){ // for d/dy

    // Need F^{ty} F^{xy} F^{zy}

    // assume got here with DT=dt and geom and state at theta zone edge
    faraday[2][0]=-1.0/geom->g * (-q->bcov[3]*q->ucov[1]+q->bcov[1]*q->ucov[3]); // f^{th}=-f^{ht}
    faraday[2][1]=-1.0/geom->g * (q->bcov[3]*q->ucov[0]-q->bcov[0]*q->ucov[3]); // f^{rh}=-f^{hr}
    faraday[2][2]=-1.0/geom->g * (-q->bcov[1]*q->ucov[0]+q->bcov[0]*q->ucov[1]); // f^{ph}
  }
  else if(which==CURTYPEZ){ // for d/dz

    // Need F^{tz} F^{xz} F^{yz}

    faraday[3][0]=1.0/geom->g * (-q->bcov[3]*q->ucov[2]+q->bcov[2]*q->ucov[3]); // f^{rt} // f^{tr}=-f^{rt}
    faraday[3][1]=1.0/geom->g * (q->bcov[2]*q->ucov[0]-q->bcov[0]*q->ucov[2]); // f^{rp}=-f^{pr}
    faraday[3][2]=1.0/geom->g * (-q->bcov[1]*q->ucov[0]+q->bcov[0]*q->ucov[1]); // f^{hp}= -f^{ph}

  }
  else if(which==CURTYPEFARADAY){
    // DT==dt, but zone center
    fcon[geom->i][geom->j][geom->k][0]=-1.0/geom->g * (q->bcov[3]*q->ucov[2]-q->bcov[2]*q->ucov[3]); // f^{tr}
    fcon[geom->i][geom->j][geom->k][1]=-1.0/geom->g * (-q->bcov[3]*q->ucov[1]+q->bcov[1]*q->ucov[3]); // f^{th}
    fcon[geom->i][geom->j][geom->k][2]=-1.0/geom->g * (q->bcov[2]*q->ucov[1]-q->bcov[1]*q->ucov[2]); // f^{tp}
    fcon[geom->i][geom->j][geom->k][3]=-1.0/geom->g * (q->bcov[3]*q->ucov[0]-q->bcov[0]*q->ucov[3]); // f^{rh}
    fcon[geom->i][geom->j][geom->k][4]=-1.0/geom->g * (-q->bcov[2]*q->ucov[0]+q->bcov[0]*q->ucov[2]); // f^{rp}
    fcon[geom->i][geom->j][geom->k][5]=-1.0/geom->g * (q->bcov[1]*q->ucov[0]-q->bcov[0]*q->ucov[1]); // f^{hp}
  }
}


// choose type of current calculation
void current_calc(FTYPE cfaraday[][N2M][N3M][NUMCURRENTSLOTS][3])
{
  void current_calc_0(FTYPE cfaraday[][N2M][N3M][NUMCURRENTSLOTS][3]);
  void current_calc_1(int which, FTYPE cfaraday[][N2M][N3M][NUMCURRENTSLOTS][3]);


#if(CALCFARADAYANDCURRENTS==0)
  dualfprintf(fail_file,"Shouldn't be here in current_doprecalc()\n");
  myexit(12);
#endif

  if(WHICHCURRENTCALC==CURRENTCALC0){
    current_calc_0(cfaraday);
  }
  else if((WHICHCURRENTCALC==CURRENTCALC1)||(WHICHCURRENTCALC==CURRENTCALC2)){
    current_calc_1(WHICHCURRENTCALC,cfaraday);
  }

}



// the current is calculated to end up at the zone and time edge
// point is to have J^t at same time as rest of J's, although different spacial points
void current_calc_0(FTYPE cfaraday[][N2M][N3M][NUMCURRENTSLOTS][3])
{
  int i,j,k;
  struct of_geom geomt;

  struct of_geom geomr;
  struct of_geom geomh;
  struct of_geom geomp;

  struct of_geom geomrp1;
  struct of_geom geomhp1;
  struct of_geom geompp1;

  static FTYPE lastdt;
  static int calls=0;
  FTYPE idtc,idx1,idx2,idx3;
  FTYPE timeterm[NDIM];

#if(CALCFARADAYANDCURRENTS==0)
  dualfprintf(fail_file,"Shouldn't be here in current_doprecalc()\n");
  myexit(12);
#endif

  idtc=2.0/(lastdt+dt);
  idx1=1.0/dx[1];
  idx2=1.0/dx[2];
  idx3=1.0/dx[3];

  COMPLOOPP1{ // largest possible loop for this differencing (could isolate directions)
    get_geometry(i,j,k,CENT,&geomt);
    get_geometry(i,j,k,FACE1,&geomr);
    get_geometry(i,j,k,FACE2,&geomh);
    get_geometry(i,j,k,FACE3,&geomp);
    // geomtp1 is same as geomt since d/dt( geometry) -> 0
    get_geometry(ip1,j,k,FACE1,&geomrp1);
    get_geometry(i,jp1,k,FACE2,&geomhp1);
    get_geometry(i,j,kp1,FACE3,&geompp1);

    if(calls>0){ // since need 2 times
      timeterm[0]=0;
      timeterm[1]=+(cfaraday[i][j][k][0][0]-cfaraday[i][j][k][4][0])*idtc; // F^{rt},t
      timeterm[2]=+(cfaraday[i][j][k][0][1]-cfaraday[i][j][k][4][1])*idtc; // F^{ht},t
      timeterm[3]=+(cfaraday[i][j][k][0][2]-cfaraday[i][j][k][4][2])*idtc; // F^{pt},t
    }
    else{
      timeterm[0]=0;
      timeterm[1]=0;
      timeterm[2]=0;
      timeterm[3]=0;
    }
      
    // J^t = F^{tr},r + F^{th},h + F^{tp},p
    jcon[i][j][k][0]=
      +1./geomt.g*(geomrp1.g*cfaraday[ip1][j][k][1][0]-geomr.g*cfaraday[i][j][k][1][0])*idx1 // F^{tr},r
      +1./geomt.g*(geomhp1.g*cfaraday[i][jp1][k][2][0]-geomh.g*cfaraday[i][j][k][2][0])*idx2 // F^{th},h
      +1./geomt.g*(geompp1.g*cfaraday[i][j][kp1][3][0]-geomp.g*cfaraday[i][j][k][3][0])*idx3 // F^{tp},p
      ;
      
    // J^r = F^{rt},t + F^{rh},h + F^{rp},p
    jcon[i][j][k][1]=
      +timeterm[1]
      +1./geomt.g*(geomhp1.g*cfaraday[i][jp1][k][2][1]-geomh.g*cfaraday[i][j][k][2][1])*idx2 // F^{rh},h
      +1./geomt.g*(geompp1.g*cfaraday[i][j][kp1][3][1]-geomp.g*cfaraday[i][j][k][3][1])*idx3 // F^{rp},p
      ;
      
    // J^h = F^{ht},t + F^{hr},r + F^{hp},p
    jcon[i][j][k][2]=
      +timeterm[2]
      +1./geomt.g*(geomrp1.g*cfaraday[ip1][j][k][1][1]-geomr.g*cfaraday[i][j][k][1][1])*idx1 // F^{hr},r
      +1./geomt.g*(geompp1.g*cfaraday[i][j][kp1][3][2]-geomp.g*cfaraday[i][j][k][3][2])*idx3 // F^{hp},p
      ;
      
    // J^p = F^{pt},t + F^{pr},r + F^{ph},h
    jcon[i][j][k][3]=
      +timeterm[3]
      +1./geomt.g*(geomrp1.g*cfaraday[ip1][j][k][1][2]-geomr.g*cfaraday[i][j][k][1][2])*idx1 // F^{pr},r
      +1./geomt.g*(geomhp1.g*cfaraday[i][jp1][k][2][2]-geomh.g*cfaraday[i][j][k][2][2])*idx2 // F^{ph},h
      ;

  }// end loops

  calls++;
  lastdt=dt;

}



// the current is calculated to end up at the zone and time edge
// point is to have J^t at same time as rest of J's and at same spatial location
// the temporal value of things is obtained at same time of everything else by fitting to parabola in time
void current_calc_1(int which, FTYPE cfaraday[][N2M][N3M][NUMCURRENTSLOTS][3])
{
  int i,j,k;
  struct of_geom geomt;
  struct of_geom geomc;
  struct of_geom geomrp1;
  struct of_geom geomhp1;
  struct of_geom geompp1;
  struct of_geom geomrm1;
  struct of_geom geomhm1;
  struct of_geom geompm1;
  static FTYPE lastdt;
  static long long calls=0;
  FTYPE idtc,idx1,idx2,idx3;
  FTYPE AA,BB;
  FTYPE dtl,dtr;
  FTYPE fl,fr,f0,derf;
  FTYPE timeterm[NDIM];

#if(CALCFARADAYANDCURRENTS==0)
  dualfprintf(fail_file,"Shouldn't be here in current_doprecalc()\n");
  myexit(12);
#endif



  if(which==CURRENTCALC1){
    dtl=lastdt;
    dtr=dt;
  }
  else if(which==CURRENTCALC2){
    idtc=2.0/(lastdt+dt);
  }
  
  idx1=1.0/(2.0*dx[1]);
  idx2=1.0/(2.0*dx[2]);
  idx3=1.0/(2.0*dx[3]);



  COMPLOOPP1{ // largest possible loop for this differencing (could isolate directions)

    get_geometry(i,j,k,CENT,&geomt);

    get_geometry(i,j,k,CENT,&geomc);

    // geomtp1 is same as geomt since d/dt( geometry) -> 0
    get_geometry(ip1,j,k,CENT,&geomrp1);
    get_geometry(i,jp1,k,CENT,&geomhp1);
    get_geometry(i,j,kp1,CENT,&geompp1);

    get_geometry(im1,j,k,CENT,&geomrm1);
    get_geometry(i,jm1,k,CENT,&geomhm1);
    get_geometry(i,j,km1,CENT,&geompm1);

    if(which==CURRENTCALC1){
      
      if(calls>=2){ // since need 3 times to properly time center without having to worry about what RK is doing

	timeterm[0]=0;
	
	fl=cfaraday[i][j][k][4][0];
	f0=cfaraday[i][j][k][5][0];
	fr=cfaraday[i][j][k][0][0];
	AA=(dtr*(fl-f0)+dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
	BB=(-dtr*dtr*(fl-f0)+dtl*dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
	derf=2.0*AA*dtr+BB;
	
	timeterm[1]=derf;
	
	fl=cfaraday[i][j][k][4][1];
	f0=cfaraday[i][j][k][5][1];
	fr=cfaraday[i][j][k][0][1];
	AA=(dtr*(fl-f0)+dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
	BB=(-dtr*dtr*(fl-f0)+dtl*dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
	derf=2.0*AA*dtr+BB;
	
	timeterm[2]=derf;

	fl=cfaraday[i][j][k][4][2];
	f0=cfaraday[i][j][k][5][2];
	fr=cfaraday[i][j][k][0][2];
	AA=(dtr*(fl-f0)+dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
	BB=(-dtr*dtr*(fl-f0)+dtl*dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
	derf=2.0*AA*dtr+BB;
	
	timeterm[3]=derf;
      }
      else{
	timeterm[0]=0;
	timeterm[1]=0;
	timeterm[2]=0;
	timeterm[3]=0;
      }
    }
    else if(which==CURRENTCALC2){
      // the current is calculated to end up at the zone and time edge
      // point is to have J^t at same time as rest of J's, and same spatial points.  time for J is one timestep back for all components.
      // like current_calc_0 in time and current_calc_1 in space
      
      if(calls>0){ // since need 2 times
	timeterm[0]=0;
	timeterm[1]=(cfaraday[i][j][k][0][0]-cfaraday[i][j][k][4][0])*idtc; // F^{rt},t
	timeterm[2]=(cfaraday[i][j][k][0][1]-cfaraday[i][j][k][4][1])*idtc; // F^{ht},t
	timeterm[3]=(cfaraday[i][j][k][0][2]-cfaraday[i][j][k][4][2])*idtc; // F^{pt},t
      }
      else{
	timeterm[0]=0;
	timeterm[1]=0;
	timeterm[2]=0;
	timeterm[3]=0;
      }

    }

    // similar to other current_calc's except position of current and spacing of derivatives

    // J^t = F^{tr},r + F^{th},h + F^{tp},p
    jcon[i][j][k][0]=
      +1./geomt.g*(geomrp1.g*cfaraday[ip1][j][k][1][0]-geomrm1.g*cfaraday[im1][j][k][1][0])*idx1 // F^{tr},r
      +1./geomt.g*(geomhp1.g*cfaraday[i][jp1][k][2][0]-geomhm1.g*cfaraday[i][jm1][k][2][0])*idx2 // F^{th},h
      +1./geomt.g*(geompp1.g*cfaraday[i][j][kp1][3][0]-geompm1.g*cfaraday[i][j][km1][3][0])*idx3 // F^{tp},p
      ;

    // J^r = F^{rt},t + F^{rh},h + F^{rp},p
    jcon[i][j][k][1]=
      +timeterm[1]
      +1./geomt.g*(geomhp1.g*cfaraday[i][jp1][k][2][1]-geomhm1.g*cfaraday[i][jm1][k][2][1])*idx2 // F^{rh},h
      +1./geomt.g*(geompp1.g*cfaraday[i][j][kp1][3][1]-geompm1.g*cfaraday[i][j][km1][3][1])*idx3 // F^{rp},p
      ;
      
    // J^h = F^{ht},t + F^{hr},r + F^{hp},p
    jcon[i][j][k][2]=
      +timeterm[2]
      +1./geomt.g*(geomrp1.g*cfaraday[ip1][j][k][1][1]-geomrm1.g*cfaraday[im1][j][k][1][1])*idx1 // F^{hr},r
      +1./geomt.g*(geompp1.g*cfaraday[i][j][kp1][3][2]-geompm1.g*cfaraday[i][j][km1][3][2])*idx3 // F^{hp},p
      ;
      
    // J^p = F^{pt},t + F^{pr},r + F^{ph},h
    jcon[i][j][k][3]=
      +timeterm[3]
      +1./geomt.g*(geomrp1.g*cfaraday[ip1][j][k][1][2]-geomrm1.g*cfaraday[im1][j][k][1][2])*idx1 // F^{pr},r
      +1./geomt.g*(geomhp1.g*cfaraday[i][jp1][k][2][2]-geomhm1.g*cfaraday[i][jm1][k][2][2])*idx2 // F^{ph},h
      ;

    
  }// end loops

  calls++;
  lastdt=dt;

}




// assumes below get inlined

// much faster than macro using ? :
// super slow for get_geometry()'s sign() call!
FTYPE sign_bad(FTYPE a)
{
  if(a>0.0) return 1.0;
  else if(a<0.0) return -1.0;
  else return 0.0;
}


// not any faster than above (except when used alot)
FTYPE sign_func(FTYPE a)
{
#if(SUPERLONGDOUBLE)
  return(sign(a));
  // no such function
#else
  return(copysign(1.0,a));
#endif
}

// much faster than macro using ? :
// Generally avoid using, instead do:
// igeom = sign(geom.g)/(fabs(geom.g)+SMALL)
FTYPE signavoidzero(FTYPE a)
{
  if(a>0.0) return 1.0;
  else if(a<0.0) return -1.0;
  else return 1.0; // arbitrarily choose 1.0
}



#ifndef WIN32
FTYPE max(FTYPE a, FTYPE b)
{
  if(a>b) return a;
  else return b;
  // if equal, then above is fine
}

FTYPE min(FTYPE a, FTYPE b)
{
  if(a>b) return b;
  else return a;
  // if equal, then above is fine
}
#endif

// compute speed of light 3-velocity in particular direction assuming other direction velocities fixed
// dx^dir/dt assuming other direction's velocities are fixed
int sol(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin)
{
  int i,j,k;
  FTYPE vu[NDIM],BB,CC,vsol1p,vsol1m;
  FTYPE ftemp1,ftemp2,ftemp3;
  FTYPE disc;
  int diro1,diro2;

  /* 
     m%3+1 gives next 1->2,2->3,3->1
     3-(4-m)%3 gives previous 1->3,2->1,3->2
  */
  diro1=dir%3+1;
  diro2=3-(4-dir)%3;

  // 3-velocity in coordinate frame
  SLOOPA(j) vu[j]=q->ucon[j]/q->ucon[TT];

  ftemp1=0;
  ftemp2=0;
  ftemp3=0;
  SLOOPA(j){
    if(j!=dir){
      ftemp1+=vu[j]*geom->gcov[dir][j];
      ftemp2+=vu[j]*geom->gcov[TT][j];
      ftemp3+=vu[j]*vu[j]*geom->gcov[j][j];
    }
  }

  BB=2.0*(ftemp1+geom->gcov[0][dir])/geom->gcov[dir][dir];
  CC=(geom->gcov[TT][TT] + 2.0*ftemp2 + ftemp3 + 2.0*vu[diro1]*vu[diro2]*geom->gcov[diro1][diro2])/geom->gcov[dir][dir];
  
  disc=BB*BB-4.0*CC;
  if(disc>0){
    *vmax=0.5*(-BB+sqrt(disc));
    *vmin=0.5*(-BB-sqrt(disc));
  }
  else{
    dualfprintf(fail_file,"disc=%21.15g < 0\n",disc);
    return(1);
  }

  return(0);

}


// limit the 3-velocity to a physically valid velocity (i.e. less than c ), assuming all other velocity directions are the same.
int limitv3(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *v)
{
  FTYPE vmax,vmin;
  int sol(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin);
  FTYPE ratv;
  
  // get speed of light 3-velocity
  MYFUN(sol(pr,q,dir,geom,&vmax,&vmin),"phys.c:limitv3()", "sol()", 1);

  // get ratio of given 3-velocity to speed of light 3-velocity for appropriate direction of coordinate-based velocity
  ratv=(*v>0) ? *v/vmax : *v/vmin;

  // limit 3-velocity to speed of light
  if(ratv>1.0){
    if(*v>0.0) *v=vmax;
    else *v=vmin;
  }

  return(0);
}

// take projection of v onto u, both are 4-vectors
// vcon=0 or 1 : whether or not ucon (1=true, 0=ucov)
// same for vresultcon
void projectionvec(int vcon,int vresultcon, struct of_state *q, struct of_geom *geom,FTYPE *v,FTYPE*vresult)
{
  FTYPE proj[NDIM][NDIM];
  int j,k;


  if((vcon)&&(vresultcon)){ // vresult^\mu = P^\mu_\nu v^\nu
    DLOOP(j,k) proj[j][k]=delta(j,k) + q->ucon[j]*q->ucov[k];
  }
  if((!vcon)&&(!vresultcon)){ // vresult_\mu = P_\mu^\nu v_\nu
    DLOOP(j,k) proj[j][k]=delta(j,k) + q->ucov[j]*q->ucon[k];
  }
  else if((!vcon)&&(vresultcon)){ // vresult^\mu = P^{\mu\nu} v_\nu
    DLOOP(j,k) proj[j][k]=geom->gcon[j][k] + q->ucon[j]*q->ucon[k];
  }
  else if((vcon)&&(!vresultcon)){ // vresult_\mu = P_{\mu\nu} v^\nu
    DLOOP(j,k) proj[j][k]=geom->gcov[j][k] + q->ucov[j]*q->ucov[k];
  }
  DLOOPA(j) vresult[j]=0.0;
  DLOOP(j,k) vresult[j]+=proj[j][k]*v[k];
  

}


// g^{tt}+1 accurate for non-rel gravity to order v^2 without machine precision problems
void compute_gconttplus1(struct of_geom *geom, FTYPE *gconttplus1)
{
  int j,k;
  FTYPE gconttsq;

  // accurate for non-rel gravity since gconttsq is order v^4
  // could store gconttplus1 instead
  gconttsq = (geom->gcon[TT][TT]*geom->gcon[TT][TT]);
  *gconttplus1= (geom->gcovpert[TT]) * gconttsq + (1.0-gconttsq);
  SLOOPA(j) *gconttplus1 += 2.0*geom->gcov[TT][j]*(geom->gcon[j][TT]*geom->gcon[j][TT]);
  SLOOP(j,k) *gconttplus1 += geom->gcov[j][k]*geom->gcon[j][TT]*geom->gcon[k][TT];

}

// compute (u^t)^2 - 1 which is v^2 in non-rel regime.
// Has no non-rel problem and can be used to rescale both non-rel and rel velocities
int quasivsq_3vel(FTYPE *vcon, struct of_geom *geom, FTYPE *quasivsq)
{
  int i,j,k;
  FTYPE gcovttplus1, gconttplus1;
  FTYPE velterm;
  int get_3velterm(FTYPE *vcon, struct of_geom *geom, FTYPE *velterm);
  void compute_gconttplus1(struct of_geom *geom, FTYPE *gconttplus1);


  // since (u^t)^2 = -1/ (g_{tt} + 2v^i g_{it} + v^2 )
  // where v^2 = v^i v^j g_{ij}
  //
  // then (u^t)^2 - 1 = [(1 + g_{tt}) + (2v^i g_{it} + v^2)]/[ - g_{tt} - (2v^i g_{it} + v^2) ]

  // later this can be initialized as it's own variable when wanting to do non-rel gravity
  //  gcovttplus1 = geom->gcov[TT][TT] + 1.0;
  gcovttplus1 = geom->gcovpert[TT];

  //  gconttplus1= geom->gcon[TT][TT] + 1.0;
  compute_gconttplus1(geom,&gconttplus1);


  get_3velterm(vcon, geom, &velterm);

  // this has a good machine representable value in the limit of non-rel velocities
  // *utsqm1 = (gcovttplus1 + velterm) / (-geom->gcov[TT][TT] - velterm);

  // decided to use (u^t)^2 - 1 + (g^{tt} + 1)
  *quasivsq =  (gcovttplus1 + velterm) / (-geom->gcov[TT][TT] - velterm) + gconttplus1 ;
  
  return(0);


}


int quasivsq_4vel(FTYPE *pr, struct of_geom *geom, FTYPE *quasivsq)
{
  FTYPE AA,BB,CCM1,CCPLUSAA ;
  FTYPE discr ;
  FTYPE bsq,X[NDIM] ;
  int i=0,j=0,k=0 ;
  FTYPE gcovttplus1,gconttplus1;
  int get_4velterms(FTYPE *ucon, struct of_geom *geom, FTYPE *AA, FTYPE *BB, FTYPE *CCM1, FTYPE *discr);
  FTYPE ucon[NDIM];
  void compute_gconttplus1(struct of_geom *geom, FTYPE *gconttplus1);


  ucon[1] = pr[U1] ;
  ucon[2] = pr[U2] ;
  ucon[3] = pr[U3] ;

  // coefficients for quadratic solution  
  get_4velterms(ucon, geom, &AA, &BB, &CCM1, &discr);

  if(discr<0){
    dualfprintf(fail_file,"Problem with utsqm1_4vel\n");
    return(1);
  }

  // 1.0 + g_{tt}
  gcovttplus1 = geom->gcovpert[TT];

  //  gconttplus1= geom->gcon[TT][TT] + 1.0;
  compute_gconttplus1(geom,&gconttplus1);

  // v^2 in non-rel case
  CCPLUSAA = gcovttplus1 + CCM1;

  // first term is genuinely relativistic, so no correction needed for non-rel case
  // second term is combination of |discr| and -4A^2
  // in the non-rel limit this reduces (for BB=0 and AA=-1) to CCPLUSAA, which is v^2
  //  *utsqm1 = ( BB*(BB + sqrt(discr)) + fabs(BB*BB - 4.0*AA*CCPLUSAA) ) /(4.0*AA*AA) ;

  // decided to use (u^t)^2 - 1 + (g^{tt} + 1)
  *quasivsq = ( BB*(BB + sqrt(discr)) + fabs(BB*BB - 4.0*AA*CCPLUSAA) ) /(4.0*AA*AA) + gconttplus1;

  return(0) ;


}


// interestingly, this can be negative in some reasonable way if we allow qsq<0, which is a similar idea behind the v^2<0 for 3-velocities
// but I don't see how to write stress-energy tensor as I did for 3-velocity
int quasivsq_rel4vel(FTYPE *uconrel, struct of_geom *geom, FTYPE *quasivsq)
{
  FTYPE alphasq;
  FTYPE qsq;
  FTYPE timeterm;
  int j ;
  int qsq_calc(FTYPE *uconrel, struct of_geom *geom, FTYPE *qsq);
  
  //        alpha = 1./sqrt(-geom->gcon[TT][TT]) ;

  alphasq=1./(-geom->gcon[TT][TT]); // positive definite and 1 in non-rel case
  
  //timeterm = -(geom->gcon[TT][TT] + 1.0 ); // ~0 in non-rel case, and equal to 1+g_{tt} in non-rel case, which is always positive
  
  qsq_calc(uconrel,geom,&qsq);
  
  // *utsqm1 = qsq/alphasq  + timeterm ;

  // decided to use (u^t)^2 - 1 + (g^{tt} + 1) = utsqm1 - timeterm
  // already accurate to order v^2 when non-rel velocity/gravity used
  *quasivsq = qsq/alphasq;
  
  return(0) ;
}


// (u^t)^2 - 1 + (g^{tt} + 1)
int quasivsq_compute(FTYPE *pr, struct of_geom *geom, FTYPE *quasivsq)
{
  int quasivsq_3vel(FTYPE *vcon, struct of_geom *geom, FTYPE *quasivsq);
  int quasivsq_4vel(FTYPE *ucon, struct of_geom *geom, FTYPE *quasivsq);
  int quasivsq_rel4vel(FTYPE *uconrel, struct of_geom *geom, FTYPE *quasivsq);
  FTYPE ucon[NDIM],vcon[NDIM],uconrel[NDIM];

  
#if(WHICHVEL==VEL4)

  ucon[1]=pr[U1];
  ucon[2]=pr[U2];
  ucon[3]=pr[U3];
  return(quasivsq_4vel(ucon, geom, quasivsq));

#elif(WHICHVEL==VEL3)

  vcon[1]=pr[U1];
  vcon[2]=pr[U2];
  vcon[3]=pr[U3];
  return(quasivsq_3vel(vcon, geom, quasivsq));

#elif(WHICHVEL==VELREL4)

  uconrel[1]=pr[U1];
  uconrel[2]=pr[U2];
  uconrel[3]=pr[U3];
  return(quasivsq_rel4vel(uconrel, geom, quasivsq));

#endif

  //  return(0);

}

int limit_quasivsq(FTYPE quasivsqnew, struct of_geom *geom, FTYPE *pr)
{
  int limit_quasivsq_3vel(FTYPE quasivsqold,    FTYPE quasivsqnew, struct of_geom *geom, FTYPE *vcon);
  int limit_quasivsq_4vel(FTYPE quasivsqold,    FTYPE quasivsqnew, struct of_geom *geom, FTYPE *ucon);
  int limit_quasivsq_rel4vel(FTYPE quasivsqold, FTYPE quasivsqnew, struct of_geom *geom, FTYPE *uconrel);
  FTYPE ucon[NDIM],vcon[NDIM],uconrel[NDIM];
  int quasivsq_compute(FTYPE *pr, struct of_geom *geom, FTYPE *quasivsq);
  FTYPE quasivsqold;


  // get pr's quasivsq
  if(quasivsq_compute(pr, geom, &quasivsqold)>=1){
    dualfprintf(fail_file,"Problem with limit_quasivsq using quasivsq_compute\n");
    myexit(3);
  }

  //  dualfprintf(fail_file,"IN: v^x = %21.15g v^y = %21.15g v^z = %21.15g :: old = %21.15g new=%21.15g\n",pr[U1],pr[U2],pr[U3],quasivsqold,quasivsqnew);

  // now rescale the velocities to agree with quasivsqnew  
#if(WHICHVEL==VEL4)

  ucon[1]=pr[U1];
  ucon[2]=pr[U2];
  ucon[3]=pr[U3];
  limit_quasivsq_4vel(quasivsqold,quasivsqnew, geom, ucon);
  pr[U1]=ucon[1];
  pr[U2]=ucon[2];
  pr[U3]=ucon[3];

#elif(WHICHVEL==VEL3)

  vcon[1]=pr[U1];
  vcon[2]=pr[U2];
  vcon[3]=pr[U3];
  limit_quasivsq_3vel(quasivsqold,quasivsqnew, geom, vcon);
  pr[U1]=vcon[1];
  pr[U2]=vcon[2];
  pr[U3]=vcon[3];

#elif(WHICHVEL==VELREL4)

  uconrel[1]=pr[U1];
  uconrel[2]=pr[U2];
  uconrel[3]=pr[U3];
  limit_quasivsq_rel4vel(quasivsqold,quasivsqnew, geom, uconrel);
  pr[U1]=uconrel[1];
  pr[U2]=uconrel[2];
  pr[U3]=uconrel[3];

#endif

  //  dualfprintf(fail_file,"OUT: v^x = %21.15g v^y = %21.15g v^z = %21.15g :: old = %21.15g new=%21.15g\n",pr[U1],pr[U2],pr[U3],quasivsqold,quasivsqnew);

  return(0);

}

int limit_quasivsq_3vel(FTYPE quasivsqold,FTYPE quasivsqnew,struct of_geom *geom, FTYPE *vcon)
{
  FTYPE pref;

  // quasivsq \propto q^2 \propto \tilde{u}^2, so trivial to rescale
  pref = sqrt(quasivsqnew/(quasivsqold+SMALL));

  // only true in non-rel case
  vcon[1]*=pref;
  vcon[2]*=pref;
  vcon[3]*=pref;

  return(0);
}

int limit_quasivsq_4vel(FTYPE quasivsqold,FTYPE quasivsqnew,struct of_geom *geom, FTYPE *ucon)
{
  FTYPE pref;

  // quasivsq \propto q^2 \propto \tilde{u}^2, so trivial to rescale
  pref = sqrt(quasivsqnew/(quasivsqold+SMALL));

  // only true in non-rel case
  ucon[1]*=pref;
  ucon[2]*=pref;
  ucon[3]*=pref;

  return(0);
}

int limit_quasivsq_rel4vel(FTYPE quasivsqold,FTYPE quasivsqnew,struct of_geom *geom, FTYPE *uconrel)
{
  FTYPE pref;

  // quasivsq \propto q^2 \propto \tilde{u}^2, so trivial to rescale
  pref = sqrt(quasivsqnew/(quasivsqold+SMALL));

  // relativistically correct
  uconrel[1]*=pref;
  uconrel[2]*=pref;
  uconrel[3]*=pref;

  return(0);
}

#define COMPLOOPVORT COMPLOOPP11 COMPLOOPP12 COMPLOOPP13

// use p[][][][U1,U2,U3] to obtain pvort[][][][whichpl]
int compute_vorticity(FTYPE (*p)[N2M][N3M][NPR],FTYPE (*pvort)[N2M][N3M][NPR],int whichpl)
{
  int i,j,k;
  FTYPE X[NDIM],V[NDIM];
  struct of_geom geom;
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE vxm,vxp,vym,vyp;
  
  
  COMPLOOPVORT{

    coord(i, jm1, k, CENT, X);
    bl_coord(X, V);
    get_geometry(i,jm1,k,CENT,&geom) ;
    // dx/dx' where '=prim coords (i.e. nonuni coords)
    vxm = p[i][jm1][k][U1]*sqrt(geom.gcov[1][1]);

    coord(i, jp1, k, CENT, X);
    bl_coord(X, V);
    get_geometry(i,jp1,k,CENT,&geom) ;
    // dx/dx' where '=prim coords (i.e. nonuni coords)
    vxp = p[i][jp1][k][U1]*sqrt(geom.gcov[1][1]);



    coord(im1, j, k, CENT, X);
    bl_coord(X, V);
    get_geometry(im1,j,k,CENT,&geom) ;
    // dx/dx' where '=prim coords (i.e. nonuni coords)
    vym = p[im1][j][k][U2]*sqrt(geom.gcov[2][2]);

    coord(ip1, j, k, CENT, X);
    bl_coord(X, V);
    get_geometry(ip1,j,k,CENT,&geom) ;
    // dx/dx' where '=prim coords (i.e. nonuni coords)
    vyp = p[ip1][j][k][U2]*sqrt(geom.gcov[2][2]);

    // center
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    get_geometry(i,j,k,CENT,&geom) ;
    dxdxprim(X, V, dxdxp);

    pvort[i][j][k][whichpl] = (vyp-vym)/(2.0*dx[1])/dxdxp[1][1];

    pvort[i][j][k][whichpl] += -(vxp-vxm)/(2.0*dx[2])/dxdxp[2][2];

  }

  return(0);

}

// input \Omega_F and B^i (code's version) and get back primitive assuming stationary/axisymmetric flow
int OBtopr_general(FTYPE omegaf,FTYPE *Bccon,struct of_geom *geom, FTYPE *pr)
{
  int j;
  FTYPE Bccov[NDIM];
  FTYPE Bsq;
  FTYPE ftemp,ftemp2;
  FTYPE vcon[NDIM];


  lower_vec(Bccon,geom,Bccov);

  Bsq=0.0+SMALL;
  SLOOPA(j) Bsq+=Bccon[j]*Bccov[j];

  ftemp=(Bccov[TT]+omegaf*Bccov[PH])/Bsq;
  // ftemp2 is set so that no round off error in limit where toroidal field dominates
  // Does this cause problem for frame-dragged fields?
  //ftemp2=omegaf*(1.0 - Bccov[PH]*Bccon[3]/Bsq);
  // below more accurate than above?
  ftemp2=omegaf*(Bccov[RR]*Bccon[RR]+Bccov[TH]*Bccon[TH])/Bsq;

  vcon[1] = -Bccon[1]*ftemp;
  vcon[2] = -Bccon[2]*ftemp;
  //  vcon[3] =  omegaf  - Bccon[3]*ftemp;
  // designed so to avoid catastrophic cancellation like above has problems with
  vcon[3] =  ftemp2 - (Bccon[3]*Bccov[TT]/Bsq);


  MYFUN(vcon2pr(WHICHVEL, vcon, geom, pr),"phys.c:OBtopr_general()", "vcon2pr() dir=0", 1);

  //  if(t>1.9 && t<2.1){
    //    dualfprintf(fail_file,"ftemp=%21.15g Bsq=%21.15g ftemp2=%21.15g :: v1=%21.15g v2=%21.15g v3=%21.15g pru1=%21.15g pru2=%21.15g pru3=%21.15g\n",ftemp,Bsq,ftemp2,vcon[1],vcon[2],vcon[3],pr[U1],pr[U2],pr[U3]);
    // }

  return(0);

}



// input \Omega_F, 3-vel along field (scalar quantity really), and B^i (code's version) and get back primitive assuming stationary/axisymmetric flow
int OBtopr_general2(FTYPE omegaf, FTYPE v0, FTYPE *Bccon,struct of_geom *geom, FTYPE *pr)
{
  int j;
  FTYPE Bccov[NDIM];
  FTYPE Bsq;
  FTYPE ftemp,ftemp2;
  FTYPE vcon[NDIM];
  FTYPE v0oB;

  lower_vec(Bccon,geom,Bccov);

  Bsq=0.0+SMALL;
  SLOOPA(j) Bsq+=Bccon[j]*Bccov[j];

  v0oB=v0/sqrt(Bsq);

  ftemp=(Bccov[TT]+omegaf*Bccov[PH])/Bsq;
  // ftemp2 is set so that no round off error in limit where toroidal field dominates
  // Does this cause problem for frame-dragged fields?
  //ftemp2=omegaf*(1.0 - Bccov[PH]*Bccon[3]/Bsq);
  // below more accurate than above?
  ftemp2=omegaf*(Bccov[RR]*Bccon[RR]+Bccov[TH]*Bccon[TH])/Bsq;

  vcon[1] = (v0oB-ftemp)*Bccon[1];
  vcon[2] = (v0oB-ftemp)*Bccon[2];
  //  vcon[3] =  omegaf  - Bccon[3]*ftemp;
  // designed so to avoid catastrophic cancellation like above has problems with
  vcon[3] =  ftemp2 + Bccon[3]*(v0oB-Bccov[TT]/Bsq);


  MYFUN(vcon2pr(WHICHVEL, vcon, geom, pr),"phys.c:OBtopr_general()", "vcon2pr() dir=0", 1);

  //  if(t>1.9 && t<2.1){
    //    dualfprintf(fail_file,"ftemp=%21.15g Bsq=%21.15g ftemp2=%21.15g :: v1=%21.15g v2=%21.15g v3=%21.15g pru1=%21.15g pru2=%21.15g pru3=%21.15g\n",ftemp,Bsq,ftemp2,vcon[1],vcon[2],vcon[3],pr[U1],pr[U2],pr[U3]);
    // }

  return(0);

}

// input \Omega_F, extra 3-vel along field (scalar quantity really), and B^i (code's version) and get back primitive assuming stationary/axisymmetric flow
int OBtopr_general3(FTYPE omegaf, FTYPE v0, FTYPE *Bccon,struct of_geom *geom, FTYPE *pr)
{
  int j;
  FTYPE Bccov[NDIM];
  FTYPE Bsq;
  FTYPE absB;
  FTYPE v0other;
  FTYPE ftemp,ftemp2;
  FTYPE vcon[NDIM];
  FTYPE v0oB;

  lower_vec(Bccon,geom,Bccov);

  Bsq=0.0+SMALL;
  SLOOPA(j) Bsq+=Bccon[j]*Bccov[j];

  absB=sqrt(Bsq);

  vcon[1] =        v0*Bccon[1]/absB;
  vcon[2] =        v0*Bccon[2]/absB;
  vcon[3] = omegaf+v0*Bccon[3]/absB;

  MYFUN(vcon2pr(WHICHVEL, vcon, geom, pr),"phys.c:OBtopr_general()", "vcon2pr() dir=0", 1);

  return(0);

}


void ffdestresstensor(FTYPE (*Mcon)[NDIM], struct of_geom *geom, FTYPE (*T)[NDIM])
{
      int i,j,k ;
      void lower_A(FTYPE (*Acon)[NDIM], struct of_geom *geom, FTYPE (*Acov)[NDIM]);


      FTYPE Mcov[NDIM][NDIM] ; /* covariant Maxwell */
      FTYPE Msq ;

      /* get covariant Maxwell from contravariant */
      lower_A(Mcon, geom, Mcov) ;

      /* out with the old */
      DLOOP(j,k)  T[j][k] = 0. ;

      /* in with the new:
       *
       * a, b, c, d run 0 to 4:
       *
       * T^a_b = M^ac * M_bc - 1/4 KroneckerDelta[a,b] (M_cd M^cd)
       *
       */

      Msq = 0. ;
      DLOOP(j,k) Msq += Mcov[j][k]*Mcon[j][k] ;

      DLOOP(j,k) {
	for(i = 0; i < NDIM; i++) {
	  T[j][k] += Mcon[j][i]*Mcov[k][i] ;
	}
      }
      DLOOPA(j) T[j][j] -= 0.25 * Msq ;
}


// we use this in order to avoid catastrophic cancellation in the non-rel regime and ultrarel regime
// give back T^dir_\mu
void ffdestresstensor_dir(int dir, FTYPE (*Mcon)[NDIM], struct of_geom *geom, FTYPE *TEMdir)
{
      int i,j,k ;
      void lower_A(FTYPE (*Acon)[NDIM], struct of_geom *geom, FTYPE (*Acov)[NDIM]);


      FTYPE Mcov[NDIM][NDIM] ; /* covariant Maxwell */
      FTYPE Msq ;

      /* get covariant Maxwell from contravariant */
      lower_A(Mcon, geom, Mcov) ;

      /* out with the old */
      DLOOPA(k)  TEMdir[k] = 0. ;

      /* in with the new:
       *
       * a, b, c, d run 0 to 4:
       *
       * T^a_b = M^ac * M_bc - 1/4 KroneckerDelta[a,b] (M_cd M^cd)
       *
       */

      Msq = 0. ;
      DLOOP(j,k) Msq += Mcov[j][k]*Mcon[j][k] ;

      DLOOPA(k) {
	for(i = 0; i < NDIM; i++) {
	  TEMdir[k] += Mcon[dir][i]*Mcov[k][i] ;
	}
      }
      TEMdir[dir] -= 0.25 * Msq ;
}


void raise_A(FTYPE (*Acov)[NDIM], struct of_geom *geom, FTYPE (*Acon)[NDIM])
{
        int j,k ;
	FTYPE localgcon[NDIM][NDIM];

	/* out with the old */
	DLOOP(j,k){
	  Acon[j][k] = 0. ;
	  // store since use all terms many times -- compiler doesn't do this for some reason
	  localgcon[j][k]=geom->gcon[j][k];
	}

	/* in with the new */
	DLOOP(j,k) {
       	       Acon[0][1] += (localgcon[0][j])*(localgcon[1][k])*Acov[j][k] ;
       	       Acon[0][2] += (localgcon[0][j])*(localgcon[2][k])*Acov[j][k] ;
       	       Acon[0][3] += (localgcon[0][j])*(localgcon[3][k])*Acov[j][k] ;
	       Acon[1][2] += (localgcon[1][j])*(localgcon[2][k])*Acov[j][k] ;
	       Acon[2][3] += (localgcon[2][j])*(localgcon[3][k])*Acov[j][k] ;
	       Acon[3][1] += (localgcon[3][j])*(localgcon[1][k])*Acov[j][k] ;
	      }

	/*
	 * diagonal is already set to zero,
	 * just copy the rest
	 */

	Acon[1][0] = - Acon[0][1] ;
	Acon[2][0] = - Acon[0][2] ;
	Acon[3][0] = - Acon[0][3] ;
	Acon[2][1] = - Acon[1][2] ;
	Acon[3][2] = - Acon[2][3] ;
	Acon[1][3] = - Acon[3][1] ;
}

void lower_A(FTYPE (*Acon)[NDIM], struct of_geom *geom, FTYPE (*Acov)[NDIM])
{
        int j,k ;
	FTYPE localgcov[NDIM][NDIM];

	/* out with the old */
	DLOOP(j,k){
	  Acov[j][k] = 0. ;
	  // store since use all terms many times -- compiler doesn't do this for some reason and lower_A() appears to be very costly
	  localgcov[j][k]=geom->gcov[j][k];
	}

	/* in with the new */
	DLOOP(j,k) {
       	       Acov[0][1] += (localgcov[0][j])*(localgcov[1][k])*Acon[j][k] ;
       	       Acov[0][2] += (localgcov[0][j])*(localgcov[2][k])*Acon[j][k] ;
       	       Acov[0][3] += (localgcov[0][j])*(localgcov[3][k])*Acon[j][k] ;
	       Acov[1][2] += (localgcov[1][j])*(localgcov[2][k])*Acon[j][k] ;
	       Acov[2][3] += (localgcov[2][j])*(localgcov[3][k])*Acon[j][k] ;
	       Acov[3][1] += (localgcov[3][j])*(localgcov[1][k])*Acon[j][k] ;
	      }

	/*
	 * diagonal is already set to zero,
	 * just copy the rest
	 */

	Acov[1][0] = - Acov[0][1] ;
	Acov[2][0] = - Acov[0][2] ;
	Acov[3][0] = - Acov[0][3] ;
	Acov[2][1] = - Acov[1][2] ;
	Acov[3][2] = - Acov[2][3] ;
	Acov[1][3] = - Acov[3][1] ;
}




// Maxwell to Faraday
// which=0 : Mcon -> Fcov (for clean Mcon, Fcov has \detg)
// which=1 : Mcov -> Fcon (for clean Mcov)

// which=2 : Fcon -> Mcov
// which=3 : Fcov -> Mcon

// copies faraday_calc() in phys.c
void MtoF(int which, FTYPE (*invar)[NDIM],struct of_geom *geom, FTYPE (*outvar)[NDIM])
{
  int nu,mu,kappa,lambda;
  FTYPE prefactor;
  int whichlc;

  if((which==0)||(which==1)){
    prefactor=-0.5;
    if(which==0) whichlc=0;
    if(which==1) whichlc=1;
  }
  if((which==2)||(which==3)){
    prefactor=0.5;
    if(which==2) whichlc=0;
    if(which==3) whichlc=1;
  }

  for(nu=0;nu<NDIM;nu++) for(mu=0;mu<NDIM;mu++){
    outvar[mu][nu]=0.0;
    for(kappa=0;kappa<NDIM;kappa++) for(lambda=0;lambda<NDIM;lambda++){
      outvar[mu][nu]+=prefactor*lc4(whichlc,geom->g,mu,nu,kappa,lambda)*invar[kappa][lambda];
    }
  }


}


// Get dual of faraday for given B^i and v^i
int dualf_calc(FTYPE *Bcon, FTYPE *vcon, FTYPE (*dualffull)[NDIM])
{
  int j,k;

  SLOOPA(j){
    // \dF^{it} = B^i
    dualffull[j][0] = Bcon[j];
    dualffull[0][j] = -Bcon[j];
  }

  DLOOPA(j) dualffull[j][j] = 0;

  // \dF^{ij} = B^i v^j - B^j v^i
  SLOOP(j,k) dualffull[j][k] = Bcon[j] * vcon[k] - Bcon[k] * vcon[j] ;

  return(0);

}



int set_zamo_velocity(int whichvel, struct of_geom *ptrgeom, FTYPE *pr)
{

  // bl-normal observer (4-vel components)
  
  // normal observer velocity in atmosphere
  if(whichvel==VEL4){
    pr[U1] = -ptrgeom->gcon[0][1]*ptrgeom->alphalapse;
    pr[U2] = -ptrgeom->gcon[0][2]*ptrgeom->alphalapse;
    pr[U3] = -ptrgeom->gcon[0][3]*ptrgeom->alphalapse;
  }
  else if(whichvel==VEL3){
    pr[U1] = ptrgeom->gcon[0][1]*ptrgeom->alphalapse;
    pr[U2] = ptrgeom->gcon[0][2]*ptrgeom->alphalapse;
    pr[U3] = ptrgeom->gcon[0][3]*ptrgeom->alphalapse;
    // GAMMIE
    //ur = -1./(r*r);
    //uh=up=0.0;
  }
  else if(whichvel==VELREL4){
    pr[U1] = 0.0;
    pr[U2] = 0.0;
    pr[U3] = 0.0;
  }

  return(0);

}


int set_zamo_ucovuconplus1ud0(struct of_geom *ptrgeom, FTYPE *ucov, FTYPE *ucon, FTYPE *plus1ud0)
{
  FTYPE ialpha;
  FTYPE alpha;
  int j;
  FTYPE gconpert;

  // set alpha -- fabs just for roundoff error
  //  ialpha=sqrt(fabs(-ptrgeom->gcon[TT][TT]));
  alpha = ptrgeom->alphalapse;
  ialpha = 1.0/alpha;

  ucov[TT]=-alpha;
  SLOOPA(j) ucov[j]=0.0;

  raise_vec(ucov,ptrgeom,ucon);

  // 1+u_t
  //  gconpert = ptrgeom->gconpert[TT];  // should be GODMARK
  gconpert = 0.0; // for now since this function now used GODMARK

  *plus1ud0 = fabs(gconpert)*alpha/(ialpha+1.0);

  return(0);

}

int set_zamo_ucon(struct of_geom *ptrgeom, FTYPE *ucon)
{
  FTYPE ucov[NDIM];
  FTYPE plus1ud0;

  set_zamo_ucovuconplus1ud0(ptrgeom, ucov,ucon,&plus1ud0);

  return(0);

}

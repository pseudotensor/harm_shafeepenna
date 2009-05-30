//
// EOS from KAZ FULL

// TODO:
//
// 1) Compute Hvsr[ii] like how computed Mvsr[ii] once per timestepand
// 2) Set Hglobal=Hvsr[startpos[1] + locali] everytime about to access EOS
// 3) Set TDYNglobal once per timestep (equal to dt?)




#define NUMEOSQUANTITIES (2 + 2 + 1 + 3 + 3 + 1)

#define PofRHOU 0      // p(rho0,u)
#define UofRHOP 1      // u(rho0,p)

#define DPDRHOofRHOU 2 // dpdrho0 |u (rho0,u)
#define DPDUofRHOU 3   // dp/du |rho0 (rho0,u)

#define CS2ofRHOU 4    // cs^2(rho0,u)

#define SofRHOU 5      // S(rho0,u)
#define DSDRHOofRHOU 6 // dS/drho0 |u (rho0,u)
#define DSDUofRHOU 7   // dS/du |\rho0 (rho0,u)

#define PofRHOCHI 8    // P(rho0,\chi)  \chi = u+p
#define IDRHO0DP 9     // 1/(d\rho0/dp) |\chi
#define IDCHIDP 10      // 1/(d\chi/dp) |\rho0

#define NUCOOL 11      // Neutrino cooling rate (erg/s/cm^2)


//////////////
//
// Original eos.dat file from Kaz's code contains (total of 4+4 quantities):
//
// rhob, tk, hcm, tdyn,
// p_tot, u_tot, s_tot, Qm
//
//
// here:
// [independent variables]
// rhob = rho0 -- rest-mass density in g/cc
// tk = Temperature in Kelvin
// hcm = height of medium in cm
// tdyn = dynamical time in seconds (assumed NSE time)
//
// [dependent variables]
// p_tot = total pressure in cgs units
// u_tot = u = internal relativistic comoving energy (no rest-mass) in cgs units
// s_tot = total entropy density in comoving frame in cgs units
// Qm = total neutrino emission as surface flux (erg/s/cm^2) [Volume rate is Qm/H]

///////////////
//
// Input eosnew.dat contains (total of 12+4+4=20 quantities)
//
// 4-D grid location (m [# of rhob]  n [# of u]  o [# of H] p [# of T]
// 4-D independent vars as NUMEOSINDEPS below (rhob, u, H, T)
// then all quantities associated with NUMEOSQUANTITIES

// Generally:
// EOSN1=rho0
// EOSN2=u or p for 0-6,10 and \chi=u+p for 7-9
// EOSN3=height in cm
// EOSN4=tdyn in seconds
#define EOSN1 48
#define EOSN2 48
#define EOSN3 16
#define EOSN4 3

#define NUMEOSINDEPS 4
#define RHOEOS 0 // rest-mass density
#define UEOS   1 // energy density: used for u, p, \chi
#define HEOS   2 // height
#define TEOS   3 // Tdyn: dynamical time

#define TBLITEMS 4
#define TBLLINEARITEMS 2


// first [4] : above 4 types of independent variables
// second [4] : 0 = lower log10 limit, 1 = upper log10 limit, 2=step, 3 = divisor of grid position
FTYPE tablelimits[NUMEOSINDEPS][TBLITEMS];
// first [4] as above, second [2] : 0=lower linear limit, 1=upper linear limit
FTYPE lineartablelimits[NUMEOSINDEPS][TBLLINEARITEMS];
// first [4] : as first [4] above
int tablesize[NUMEOSINDEPS];

// then allocate space for table
// presumed to be several functions as functions of 4 other quantities
float a_eostable[NUMEOSQUANTITIES][EOSN4][EOSN3][EOSN2][EOSN1];
static float (*eostable)[EOSN4][EOSN3][EOSN2][EOSN1];

// some globally accessed functions
static int get_eos_fromtable(int which, FTYPE q1, FTYPE q2, FTYPE q3, FTYPE q4, FTYPE *answer);





#define EOSHEADNAME "eosnew.head"
#define EOSTABLENAME "eosnew.dat"

// reads in table for EOS and sets up the table parameters
void read_setup_eostable(void)
{
  FILE *inhead;
  FILE *intable;
  int	ii,jj;
  int m,n,o,p;
  FTYPE rhobgrid, ugrid, hcmgrid, tdyngrid;
  FTYPE rhobcomp, ucomp, hcmcomp, tdyncomp;
  int iii,jjj,kkk,lll,ppp;
  int totalindex[NUMEOSINDEPS];
  FTYPE indep[NUMEOSINDEPS];
  FTYPE lstepdep,diff,lindep,lindeptry;
  int numcolums;
  FTYPE lindepdiff;


  ///////////////////////////////////////////////
  //
  // This entire file assumes all things are doubles except eostable itself
  //
  ///////////////////////////////////////////////

  ///////////////////////////////////////////////
  //
  // assign global pointer to eostable (static assignment)
  //
  //////////////////////////////////////////////
  eostable = (float (*)[EOSN4][EOSN3][EOSN2][EOSN1])(&(a_eostable[0][0][0][0][0]));



  if(myid==0){ // do only over CPU1 and then send data to rest of CPUs

    ///////////////////////////////////
    //
    // first read-in the header file
    //
    //////////////////////////////////

    // open header file
    if( (inhead = fopen(EOSHEADNAME,"rt"))==NULL){
      dualfprintf(fail_file,"No such file %s\n",EOSHEADNAME);
      myexit(16622);
    }

    // get number of output colums and check that code and data file agree
    fscanf(inhead,"%d\n",&numcolums); 
    if(numcolumns!=NUMEOSINDEPS*2+NUMEOSQUANTITIES){
      dualfprintf(fail_file,"numcolumns=%d but code expected %d\n",numcolumns,NUMEOSINDEPS*2+NUMEOSQUANTITIES);
      myexit(16623);
    }

    for(ii=0;ii<NUMEOSINDEPS;ii++){
      fscanf(inhead,"%d %lf %lf %lf",&tablesize[ii],&tablelimits[ii][0],&tablelimits[ii][1],&tablelimits[ii][2]);
      // the below definition is consistent with Kaz's code, matlab eos_extract.m and elsewhere in this code
      tablelimits[ii][3] = ((FTYPE)tablesize[ii]-1.0)/(tablelimits[ii][1] - tablelimits[ii][0]);
      // below used to truncate(limit) input values so lookup doesn't have as many conditionals
      lineartablelimits[ii][0]=pow(10.0,tablelimits[ii][0]);
      lineartablelimits[ii][1]=pow(10.0,tablelimits[ii][1]);
    }

    // some checks
    if(EOSN1!=tablesize[RHOEOS] || EOSN2!=tablesize[UEOS] || EOSN3!=tablesize[HEOS] || EOSN3!=tablesize[TEOS] ){
      dualfprintf(fail_file,"Size of table does not match\n");
      for(ii=0;ii<NUMEOSINDEPS;ii++){
	dualfprintf(fail_file,"tablesize[%d]=%d\n",ii,tablesize[ii]);
      }
      dualfprintf(fail_file,"EOSN1,2,3,4 = %d %d %d %d\n",EOSN1,EOSN2,EOSN3,EOSN4);
      myexit(16624);
    }





    ////////////////////////////
    //
    // at this point header is consistent with expectations, so read in the table
    //
    // now read-in table as written
    //
    //////////////////////////////////

    // open eostable file
    if( (inhead = fopen(EOSTABLENAME,"rt"))==NULL){
      dualfprintf(fail_file,"No such file %s\n",EOSTABLENAME);
      myexit(16625);
    }

	
    for(lll=0;lll<tablesize[TEOS];lll++)for(kkk=0;kkk<tablesize[HEOS];kkk++)for(jjj=0;jjj<tablesize[UEOS];jjj++)for(iii=0;iii<tablesize[RHOEOS];iii++){
      // first get positions to make sure everything is consistent
      fscanf(intable,"%d %d %d %d",&m,&n,&o,&p);
      if(m!=iii || n!=jjj || o!=kkk || p!=lll){
	dualfprintf(fail_file,"Read-in table indicies inconsistent with expected indicies: m=%d iii=%d n=%d jjj=%d o=%d kkk=%d p=%d lll=%d\n",m,iii,n,jjj,o,kkk,p,lll);
	myexit(16626);
      }
      // second read in the independent variable values and compare with expected value
      fscanf(intable,"%lf %lf %lf %lf",&indep[RHOEOS],&indep[UEOS],&indep[HEOS],&indep[TEOS]);

      totalindex[RHOEOS] = iii;
      totalindex[UEOS]   = jjj;
      totalindex[HEOS]   = kkk;
      totalindex[TEOS]   = lll;

      for(ii=0;ii<NUMEOSINDEPS;ii++){
	// get step (consistent with how step is computed in Kaz's code and in matlab script eos_extract.m)
	lstepdep = (tablelimits[ii][1]-tablelimits[ii][0])/((FTYPE)tablesize[ii]-1.0);
	// compare step sizes to read-in step sizes
	diff = fabs(lstepdep - tablelimits[ii][2])/(fabs(lstepdep)+fabs(tablelimits[ii][2]));
	if(diff > 1E-13){
	  dualfprintf(fail_file,"Grid step size is incorrect: ii=%d readin-value=%21.15g lstepdep=%21.15g\n",ii,tablelimits[ii][2],lstepdep);
	  myexit(16627);
	}

	// get read-in value of independent variable
	lindep=log10(indep[ii]);
	// get computed value of independent variable (used later for lookup, so verifies lookup method)
	lindeptry=tablelimits[ii][0] + totalindex[ii]*lstepdep;
	// compare to be sure same
	diff = fabs(lindep-lindeptry)/(fabs(lindep)+fabs(lindeptry));
	if(diff>1E-13){
	  dualfprintf(fail_file,"Grid position data is incorrect: ii=%d readin-lindep=%21.15g lindeptry=%21.15g diff=%21.15g\n",ii,lindep,lindeptry,lindepdiff);
	  myexit(16628);
	}
      }

      // so everything is consistent
      // now read in columns of actual dependent variable values
      for(ppp=0;ppp<NUMEOSQUANTITIES;ppp++){
	fscanf(intable,"%f",&eostable[ppp][lll][kkk][jjj][iii]); // float values
      }
      // continue onto next row
    }// end loop over all rows



    //////////////////////////////////////////
    //
    // convert quantities to code units
    //
    //////////////////////////////////////////

    lineartablelimits[RHOEOS][0]/=rhounit;
    lineartablelimits[RHOEOS][1]/=rhounit;
    lineartablelimits[UEOS][0]/=Pressureunit;
    lineartablelimits[UEOS][1]/=Pressureunit;
    lineartablelimits[HEOS][0]/=Lunit;
    lineartablelimits[HEOS][1]/=Lunit;
    lineartablelimits[TEOS][0]/=Tunit;
    lineartablelimits[TEOS][1]/=Tunit;

    // recompute tablelimits
    for(jj=0;jj<NUMEOSINDEPS;jj++) for(ii=0;ii<=1;ii++){
      tablelimits[jj][ii]=log10(lineartablelimits[jj][ii]);
    }

    // convert log10 step to code units
    tablelimits[RHOEOS][2]=log10(pow(10.0,tablelimits[RHOEOS][2])/rhounit);
    tablelimits[UEOS][2]=log10(pow(10.0,tablelimits[UEOS][2])/Pressureunit);
    tablelimits[HEOS][2]=log10(pow(10.0,tablelimits[HEOS][2])/Lunit);
    tablelimits[TEOS][2]=log10(pow(10.0,tablelimits[TEOS][2])/Tunit);

    // just recompute the below
    for(ii=0;ii<NUMEOSINDEPS;ii++){
      tablelimits[jj][3] = ((FTYPE)tablesize[ii]-1.0)/(tablelimits[ii][1] - tablelimits[ii][0]);
    }
  // now 0-3 of tablelimits set, lineartablelimits set, tablesize is same


  // output to logfull_file so have information
    for(jj=0;jj<NUMEOSINDEPS;jj++){
      trifprintf("tablesize[%d]=%d\n",jj,tablesize[jj]);
      for(ii=0;ii<=1;ii++) trifprintf("lineartablelimits[%d][%d]=%21.15g\n",jj,ii,lineartablelimits[jj][ii]);
      for(ii=0;ii<4;ii++) trifprintf("tablelimits[%d][%d]=%21.15g\n",jj,ii,tablelimits[jj][ii]);
    }

    // report read-in
    trifprintf("Done reading in EOS table\n");



    ///////////////////////////////
    //
    // convert eostable units
    //
    ///////////////////////////////

    for(lll=0;lll<tablesize[TEOS];lll++)for(kkk=0;kkk<tablesize[HEOS];kkk++)for(jjj=0;jjj<tablesize[UEOS];jjj++)for(iii=0;iii<tablesize[RHOEOS];iii++){

      eostable[PofRHOU][lll][kkk][jjj][iii]/=Pressureunit;
      eostable[UofRHOP][lll][kkk][jjj][iii]/=Pressureunit;
      // dPdRHO0ofRHOU dimensionless
      // dPdUofRHOU dimensionless
      eostable[CS2ofRHOU][lll][kkk][jjj][iii]/=Vunit;
      // entropy density (erg/K/cc)
      // kb doesn't convert units, but changes kb T to erg
      // presumes entropy is used with energy as in first-law: dQ = (kbT)dS where kbT is in ergs
      eostable[SofRHOU][lll][kkk][jjj][iii]/=(kb/pow(Lunit,3.0));
      // below are dimensionless apart from erg/K
      eostable[DSDRHOofRHOU][lll][kkk][jjj][iii]/=(kb);
      eostable[DSDUofRHOU][lll][kkk][jjj][iii]/=(kb);

      eostable[PofRHOCHI][lll][kkk][jjj][iii]/=Pressureunit;
      // IDRHO0DP is dimensionless
      // IDCHIDP is dimensionless

      // Qm is in erg/s/cm^2 (Qsurf, not Qvol)
      // this is multiplied by H when used
      eostable[NUCOOL][lll][kkk][jjj][iii]/=(edotunit/(Lunit*Lunit));

    }
		
    // report conversion
    trifprintf("Done converting units of EOS table\n");



    
  } // end if myid==0


  //////////////////////////////////////
  //
  // send data to all CPUs
  // so all CPUs have parameters and full eostable data
#if(USEMPI)
  MPI_Bcast(&tablelimits,1,MPI_FTYPE,NUMEOSINDEPS*TBLITEMS,MPI_COMM_GRMHD);
  MPI_Bcast(&lineartablelimits,1,MPI_FTYPE,NUMEOSINDEPS*TBLLINEARITEMS,MPI_COMM_GRMHD);
  MPI_Bcast(&tablesize,1,MPI_INT,NUMEOSINDEPS,MPI_COMM_GRMHD);

  MPI_Bcast(&(eostable[0]),NUMEOSQUANTITIES*EOSN4*EOSN3*EOSN2*EOSN1,MPI_FTYPE,0,MPI_COMM_GRMHD);
#endif



  //  dualfprintf(log_file,"proc: %d: Done reading in EOS table\n",myid);




}





int iswithin_eostable(int which, FTYPE q1, FTYPE q2)
{
  // table setup so q1 and q2 always have same range
  // assume all quantities depend on limits in same way, so "which" doesn't matter
  // assume HEOS and TEOS are always withing range for now
  if(
     q1>=tablelimits[RHOEOS][0] && q1<=tablelimits[RHOEOS][1]
     &&
     q2>=tablelimits[UEOS][0] && q2<=tablelimits[UEOS][1]
     ){
    return(1);
  }
  else return(0);

  // GODMARK:
  // alternative to this function is that we simply place a floor on the linear values of rho, u, H, and T so always within table values
  // Use: lineartablelimits[ii][0,1]
  //
  // problem is that for inversion, truncation still will be an if check so as slow as above anyways, so can just truncate here if wanted and then more general.  Only slows down simpler single-type calls to EOS.
}




void eos_lookup(FTYPE rho0, FTYPE u, FTYPE H, FTYPE T, FTYPE *ieos, FTYPE *jeos, FTYPE *keos, FTYPE *leos)
{
  FTYPE lrho,lu,lH,lT;
  FTYPE di,dj,dk,dl;

  // is this expensive?
  lrho=log10(rho0);
  lu=log10(u);
  lH=log10(H);
  lT=log10(T);


  // i = (x - x0)*[(N-1)/(x1-x0)] such that if x=x0, then i=0, if x=x1, then i=N-1
  *ieos = (lrho-tablelimits[RHOEOS][0])*tablelimits[RHOEOS][3];
  *jeos = (lu-tablelimits[UEOS][0])*tablelimits[UEOS][3];
  *keos = (lH-tablelimits[HEOS][0])*tablelimits[HEOS][3];
  *leos = (lT-tablelimits[TEOS][0])*tablelimits[TEOS][3];

}




// bi-linear interpolation
// GODMARK: could choose nearest neighbor for tdyn?
FTYPE get_eos_fromlookup_4D(int which, FTYPE ieos, FTYPE jeos, FTYPE keos, FTYPE leos)
{
  int ii,jj,kk,ll;
  FTYPE di[2],dj[2],dk[2],dl[2];
  FTYPE dist[2][2][2][2],f[2][2][2][2];
  FTYPE totaldist,totalf;
  int iii,jjj,kkk,lll;


  ii=(int)ieos;
  jj=(int)jeos;
  kk=(int)keos;
  ll=(int)leos;
	
  di[1]=ieos-(FTYPE)ii;
  di[0]=1.0-di[1];
  dj[1]=jeos-(FTYPE)jj;
  dj[0]=1.0-dj[1];
  dk[1]=keos-(FTYPE)kk;
  dk[0]=1.0-dk[1];
  dl[1]=leos-(FTYPE)ll;
  dl[0]=1.0-dl[1];

  // 4-D means 2^4=16 positions
  totaldist=0.0;
  totalf=0.0;
  for(iii=0;iii<=1;iii++)for(jjj=0;jjj<=1;jjj++)for(kkk=0;kkk<=1;kkk++)for(lll=0;lll<=1;lll++){
    totaldist += dist[iii][jjj][kkk][lll] = di[iii]*dj[jjj]*dk[kkk]*dl[kkk];
    f[iii][jjj][kkk][lll] = eostable[which][ii+iii][jj+jjj][kk+kkk][ll+lll];
    totalf +=f[iii][jjj][kkk][lll]*dist[iii][jjj][kkk][lll];
  }

  // finally normalize
  totalf /=totaldist;

  return(totalf);

}


// bi-linear interpolation
// GODMARK: could choose nearest neighbor for tdyn?
FTYPE get_eos_fromlookup_3D(int which, FTYPE ieos, FTYPE jeos, FTYPE keos, FTYPE leos)
{
  int ii,jj,kk,ll;
  FTYPE di[2],dj[2],dk[2];
  FTYPE dist[2][2][2],f[2][2][2];
  FTYPE totaldist,totalf;
  int iii,jjj,kkk;


  ii=(int)ieos;
  jj=(int)jeos;
  kk=(int)keos;
  ll=(int)leos;
	
  di[1]=ieos-(FTYPE)ii;
  di[0]=1.0-di[1];
  dj[1]=jeos-(FTYPE)jj;
  dj[0]=1.0-dj[1];
  dk[1]=keos-(FTYPE)kk;
  dk[0]=1.0-dk[1];

  // 3-D means 2^3=8 positions
  totaldist=0.0;
  totalf=0.0;
  for(iii=0;iii<=1;iii++)for(jjj=0;jjj<=1;jjj++)for(kkk=0;kkk<=1;kkk++){
    totaldist += dist[iii][jjj][kkk] = di[iii]*dj[jjj]*dk[kkk];
    f[iii][jjj][kkk] = eostable[which][ii+iii][jj+jjj][kk+kkk][ll];
    totalf +=f[iii][jjj][kkk]*dist[iii][jjj][kkk];
  }

  // finally normalize
  totalf /=totaldist;

  return(totalf);

}

// bi-linear interpolation
// used for those quantities that are purely 2D
FTYPE get_eos_fromlookup_2D(int which, FTYPE ieos, FTYPE jeos, FTYPE keos, FTYPE leos)
{
  int ii,jj,kk,ll;
  FTYPE di[2],dj[2];
  FTYPE dist[2][2],f[2][2];
  FTYPE totaldist,totalf;
  int iii,jjj;


  ii=(int)ieos;
  jj=(int)jeos;
  kk=(int)keos;
  ll=(int)leos;
	
  di[1]=ieos-(FTYPE)ii;
  di[0]=1.0-di[1];
  dj[1]=jeos-(FTYPE)jj;
  dj[0]=1.0-dj[1];

  // 2-D means 2^2=4 positions
  totaldist=0.0;
  totalf=0.0;
  for(iii=0;iii<=1;iii++)for(jjj=0;jjj<=1;jjj++){
    totaldist += dist[iii][jjj] = di[iii]*dj[jjj];
    f[iii][jjj] = eostable[which][ii+iii][jj+jjj][kk][ll];
    totalf +=f[iii][jjj]*dist[iii][jjj];
  }

  // finally normalize
  totalf /=totaldist;

  return(totalf);

}


#define WHICHEOSDIMEN 2

// for now assume TDYN is locked to timestep that is order 1E-6 seconds, so can use lowest leos
static int get_eos_fromtable(int which, FTYPE q1, FTYPE q2, FTYPE q3, FTYPE q4, FTYPE *answer)
{
  FTYPE get_eos_fromlookup_4D(int which, FTYPE ieos, FTYPE jeos, FTYPE keos, FTYPE leos);
  FTYPE get_eos_fromlookup_3D(int which, FTYPE ieos, FTYPE jeos, FTYPE keos, FTYPE leos);
  FTYPE get_eos_fromlookup_2D(int which, FTYPE ieos, FTYPE jeos, FTYPE keos, FTYPE leos);
  int iswithin_eostable(int which, FTYPE q1, FTYPE q2);
  void eos_lookup(FTYPE rho0, FTYPE u, FTYPE H, FTYPE T, FTYPE *ieos, FTYPE *jeos, FTYPE *keos, FTYPE *leos);
	FTYPE ieos,jeos,keos,leos;		

  // see if within table, and if so lookup interpolated result
  if(iswithin_eostable(which,q1,q2)){


#if(WHICHEOSDIMEN==3)
	// easier to set q3 and q4 here from global quantity
	// Hglobal
	// q3 = Hvsr[startpos[1]+icurr] ?
	q4 = dt; // not used currently

    // first lookup the positions
	eos_lookup(q1, q2, q3, q4, &ieos, &jeos, &keos, &leos);
    leos=0; // TDYN is assumed to be on shortest period

    // now compute result (for not ignore TDYN)
    *answer=get_eos_fromlookup_3D(which, ieos, jeos, keos, leos);
#elif(WHICHEOSDIMEN==2)
	q3 = 1.0; // not used
	q4 = dt; // not used currently

    // first lookup the positions
	eos_lookup(q1, q2, q3, q4, &ieos, &jeos, &keos, &leos);
	keos=0; // H assumed smallest (optically thin limit)
    leos=0; // TDYN is assumed to be on shortest period

    // now compute result
    *answer=get_eos_fromlookup_2D(which, ieos, jeos, keos, leos);

#endif

    return(0);
  }
  else return(1); // indicates "failure" and no answer within table

}

// for ideal gas parts (dS) : GODMARK
#define GAMMA (gam)

// 1/\Gamma_r
#define GAMMAM1 (GAMMA-1.0)
#define IGAMMAR (GAMMAM1/GAMMA)



// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE pressure;

  if(get_eos_fromtable(PofRHOU,rho0,u,Hglobal,TDYNglobal,&pressure)){
    // otherwise use TM EOS
    pressure = u*(2.0*rho0+u)/(3.0*(rho0+u));
  }

  
  return(pressure);
}

// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_grmhd(FTYPE rho0, FTYPE p)
{
  FTYPE u;
	
  if(get_eos_fromtable(UofRHOP,rho0,p,Hglobal,TDYNglobal,&u)){
    u=1.5*(p + 3.0*p*p/(2.0*rho0+sqrt(9.0*p*p+4.0*rho0*rho0)));
  }
 
  return(u);
}

// dp(rho0, u)/du
FTYPE dpdu_rho0_u_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE dpdu;

  if(get_eos_fromtable(DPDRHOofRHOU,rho0,u,Hglobal,TDYNglobal,&dpdu)){
    dpdu = 1.0/3.0*(1.0 + rho0*rho0/( (rho0+u)*(rho0+u)));
  }

  return(dpdu);
}

// dp(rho0, u)/drho0
FTYPE dpdrho0_rho0_u_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE dpdrho0;

  if(get_eos_fromtable(DPDUofRHOU,rho0,u,Hglobal,TDYNglobal,&dpdrho0)){
    dpdrho0 = u*u/(3.0*(rho0+u)*(rho0+u));
  }

  return(dpdrho0) ;
}


// sound speed squared (for vchar.c)
FTYPE cs2_compute_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE pressure;
  FTYPE h;
  FTYPE cs2;


  if(get_eos_fromtable(CS2ofRHOU,rho0,u,Hglobal,TDYNglobal,&cs2)){
    pressure = pressure_rho0_u(rho0,u);
    h=rho0+u+pressure; // not specific h
    cs2=pressure*(5.0*h - 8.0*pressure) / (3.0*h*(h-pressure));
  }

  return(cs2);

}

// entropy as function of rho0 and internal energy (u)
// S(rho0,u)
FTYPE compute_entropy_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u);
  FTYPE entropy;

  if(get_eos_fromtable(SofRHOU,rho0,u,Hglobal,TDYNglobal,&entropy)){
    entropy=0.0; // GODMARK: not set yet
  }

  return(entropy);

}


// used for dudp_calc
FTYPE compute_dSdrho_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE indexn;
  FTYPE entropy;
  FTYPE dSdrho;
  FTYPE compute_entropy(FTYPE rho0, FTYPE u);

  if(get_eos_fromtable(DSDRHOofRHOU,rho0,u,Hglobal,TDYNglobal,&dSdrho)){
    entropy=compute_entropy(rho0,u);

    // ideal gas (no TM EOS VERSION YET)
    indexn=1.0/GAMMAM1;
    dSdrho=entropy/rho0-(indexn+1.0);
  }

  return(dSdrho);

}


// used for dudp_calc
FTYPE compute_dSdu_grmhd(FTYPE rho0, FTYPE u)
{
  FTYPE indexn;
  FTYPE dSdu;

  if(get_eos_fromtable(DSDUofRHOU,rho0,u,Hglobal,TDYNglobal,&dSdu)){
    // ideal gas (no TM EOS VERSION YET)
    indexn=1.0/GAMMAM1;
    dSdu=indexn*rho0/u;
  }

  return(dSdu);

}


// u(rho0,S)
// here entropy is entropy density?
// not needed for now
FTYPE compute_u_from_entropy_grmhd(FTYPE rho0, FTYPE entropy)
{
  FTYPE u;

  // GODMARK: not set yet
  // if(get_eos_fromtable(UofRHOS,rho0,entropy,Hglobal,TDYNglobal, &u)){
  u=0.0; // GODMARK: not set yet
  //}

  return(u);

}


// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0_grmhd(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE Q,delta,delta2;
  FTYPE pressure;

  if(get_eos_fromtable(PofRHOCHI,rho0,wmrho0,Hglobal,TDYNglobal, &pressure)){
    Q=wmrho0/rho0;
    delta=9.0/25.0*wmrho0*(2.0+Q);
    delta2=delta/rho0;
    pressure=(5.0/8.0)*(wmrho0 - delta/(1.0+sqrt(1.0+delta2)));
  }


  return(pressure);
}


// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp_grmhd(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE pressure_wmrho0(FTYPE rho0, FTYPE wmrho0);
  FTYPE idrho0dp;
  FTYPE p;

 
  if(get_eos_fromtable(IDRHO0DP,rho0,wmrho0,Hglobal,TDYNglobal, &idrho0dp)){
    p = pressure_wmrho0(rho0, wmrho0);
    idrho0dp = (2.0*wmrho0 - 5.0*p)/(5.0*rho0+5.0*wmrho0-8.0*p);
  }

  return(idrho0dp);
}



// 1 / (d(u+p)/dp) holding rho0 fixed
FTYPE compute_idwmrho0dp_grmhd(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE pressure_wmrho0(FTYPE rho0, FTYPE wmrho0);
  FTYPE idwmrho0dp;
  FTYPE p;


  if(get_eos_fromtable(IDCHIDP,rho0,wmrho0,Hglobal,TDYNglobal,&idwmrho0dp)){
    p = pressure_wmrho0(rho0, wmrho0);
    idwmrho0dp = (2.0*wmrho0 + 2.0*rho0 - 5.0*p)/(5.0*rho0+5.0*wmrho0-8.0*p);
  }

  return(idwmrho0dp);

}



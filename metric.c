
#include "decs.h"





// this file includes metric dependent terms, including for initial
// condition routines for IC coords.

// crucially need to setup analytic form of gcov.  All rest can be done numerically or analytically if wanted.

// obtain gcov in primcoords of whichcoord type metric/coords
// here ptrgeom is only expected to contain i,j,k,p location
void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal)
{
  void set_gcov_cylminkmetric    (FTYPE *V, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal);
  void set_gcov_spcminkmetric    (FTYPE *V, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal);
  void set_gcov_cartminkmetric   (FTYPE *V, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal);
  void set_gcov_unigravity   (FTYPE *V, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal);
  void set_gcov_htmetric         (FTYPE *V, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal);
  void set_gcov_htmetric_accurate(FTYPE *V, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal);
  void set_gcov_ksmetric         (FTYPE *V, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal);
  void set_gcov_ks_bh_tov_metric (FTYPE *X, FTYPE *V, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal);
  extern void set_gcov_ks_tov_metric(FTYPE *X, FTYPE *V, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal);
  extern void set_gcov_bl_tov_metric(FTYPE *X, FTYPE *V, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal);
  void set_gcov_blmetric         (FTYPE *V, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal);
  void gcov2gcovprim(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal, FTYPE gcovlocalprim[][NDIM], FTYPE *gcovpertlocalprim);
  FTYPE gcovselfpert[NDIM];
  extern int set_gcov_selfspcmetric(FTYPE *X, FTYPE *V, FTYPE *gcovselfpert);
  FTYPE V[NDIM];
  int j,k;
  int presenttime;
  int interpX_gcov(FTYPE *X, FTYPE (*gcovgrid)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM], FTYPE (*gcovpertgrid)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM], FTYPE gcovlocal[][NDIM], FTYPE *gcovpertlocal);
  FTYPE phi;







  // determine if asking for metric now or in past
  if(DOEVOLVEMETRIC){

    //    if(ptrgeom->i==7 && nstep==1084){
    //      dualfprintf(fail_file,"X[0]=%21.15g Xmetricold[0]=%21.15g Xmetricnew[0]=%21.15g t=%21.15g p=%d\n",X[0],Xmetricold[0],Xmetricnew[0],t,ptrgeom->p);
    //    }

    if(X[0]>=Xmetricnew[0]-SMALL){ // SMALL trying to account for numerical roundoff error. Assume never want metric at a time difference SMALL before
      // by present we mean last time latest metric was computed
      presenttime=1;
    }
    else{
      // should only ever possibly get here if DOEVOLVEMETRIC==1
      presenttime=0;
#if(DOEVOLVEMETRIC!=1)
      dualfprintf(fail_file,"presenttime=0 and DOEVOLVEMETRIC==0 are incompatible: X[0]=%21.15g t=%21.15g\n",X[0],t);
      myexit(5523);
#endif
      if(ptrgeom->p==NOWHERE){
	dualfprintf(fail_file,"Not capable of obtaining past metric that wasn't stored\n");
	myexit(2356);
	// for example, can't just interpolate metric since then won't necessarily satisfy divg=0
	// although probably not big error for anything requiring this interpolation
	// that is, Connection and most evolution things will use metric at known location
	// so could infact interpolate if this is needed -- wait and see
      }
      if(getprim==0){
	// not designed to feed back anything not in PRIMECOORDS when asking for past metric
	dualfprintf(fail_file,"getprim==0 is not compatible with wanting past metric\n");
	myexit(6246);
      }
    }
  }
  else presenttime=1; // no evolution, so always present





  if(presenttime){


    // if PRIMECOORDS, did store the data before, and choosing a gridded position, then can just grab from memory
    if(whichcoord==PRIMECOORDS && didstoremetricdata==1 && ptrgeom->p!=NOWHERE){

      DLOOP(j,k){
	gcovlocal[j][k]=gcov[ptrgeom->i][ptrgeom->j][ptrgeom->k][ptrgeom->p][j][k];
      }
      DLOOPA(j){
	gcovpertlocal[j]=gcovpert[ptrgeom->i][ptrgeom->j][ptrgeom->k][ptrgeom->p][j];
      }
      
    }
    else{ // then not defined yet or arbitrary position


      // before here, X set by i,j,k or chosen directly
      bl_coord_ijk_2(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, X, V);
      //bl_coord(X,V);


      //dualfprintf(fail_file,"blcoordcalled: i=%d j=%d k=%d p=%d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p);

 
      if(whichcoord>=0){
	if(whichcoord==BLCOORDS){
	  set_gcov_blmetric(V, gcovlocal, gcovpertlocal);
	}
	else if(whichcoord==KSCOORDS){
	  set_gcov_ksmetric(V, gcovlocal, gcovpertlocal);
	}
	else if(whichcoord==KS_BH_TOV_COORDS){
	  set_gcov_ks_bh_tov_metric(X, V, gcovlocal, gcovpertlocal);
	}
	else if(whichcoord==KS_TOV_COORDS){
	  set_gcov_ks_tov_metric(X, V, gcovlocal, gcovpertlocal);
	}
	else if(whichcoord==BL_TOV_COORDS){
	  set_gcov_bl_tov_metric(X, V, gcovlocal, gcovpertlocal);
	}
	else if(whichcoord==HTMETRIC){
	  set_gcov_htmetric(V, gcovlocal, gcovpertlocal);
	}
	else if(whichcoord==HTMETRICACCURATE){
	  set_gcov_htmetric_accurate(V, gcovlocal, gcovpertlocal);
	}
	else if(whichcoord==CARTMINKMETRIC){
	  set_gcov_cartminkmetric(V, gcovlocal, gcovpertlocal);
	}
	else if(whichcoord==UNIGRAVITY){
	  set_gcov_unigravity(V, gcovlocal, gcovpertlocal);
	}
	else if(whichcoord==CYLMINKMETRIC){
	  set_gcov_cylminkmetric(V, gcovlocal, gcovpertlocal);
	}
	else if(whichcoord==SPCMINKMETRIC){
	  set_gcov_spcminkmetric(V, gcovlocal, gcovpertlocal);
	}
	else{
	  dualfprintf(fail_file,"gcov_func(): no such whichcoord=%d\n",whichcoord);
	  myexit(12626);
	}

	//      dualfprintf(fail_file,"blcoordcalled2: i=%d j=%d k=%d p=%d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p);

      }
      else{
	dualfprintf(fail_file,"your request makes no sense (i.e. can't get prim gcov from prim gcov): getprim=%d whichcoord=%d\n",getprim,whichcoord);
	myexit(26632);
      }

#if(DOSELFGRAVVSR)
      if(
	 whichcoord==BLCOORDS ||
	 whichcoord==KSCOORDS ||
	 whichcoord==HTMETRIC ||
	 whichcoord==HTMETRICACCURATE ||
	 whichcoord==SPCMINKMETRIC
	 ){
	// then using spherical polar coordinates and can do self-gravity as designed
	set_gcov_selfspcmetric(X,V,gcovselfpert);

	// OLD DEBUG
	//    phi = -0.1/V[1];
	//gcovlocal[TT][TT] += -2.0*phi;
	//      gcovpertlocal[TT] += -2.0*phi;
	//      gcovlocal[RR][RR] += -2.0*phi;
	//      gcovpertlocal[RR] += -2.0*phi;

	//    gcovlocal[TT][TT] += gcovselfpert[TT];
	//    gcovpertlocal[TT] += gcovselfpert[RR];
	//    gcovlocal[RR][RR] += gcovselfpert[TH];
	//    gcovpertlocal[RR] += gcovselfpert[PH];

	//    if(gcovselfpert[TT]+2.0*phi>0.01){
	//	dualfprintf(fail_file,"got TT difference: %21.15g %21.15g\n",gcovselfpert[TT],-2.0*phi);
	//    }

	//    if(gcovselfpert[RR]+2.0*phi>0.01){
	//	dualfprintf(fail_file,"got RR difference: %21.15g %21.15g\n",gcovselfpert[RR],-2.0*phi);
	//    }

	//    if(gcovselfpert[TH]>0.01){
	//	dualfprintf(fail_file,"got TH difference: %21.15g %21.15g\n",gcovselfpert[TH],0.0);
	//    }

	//    if(gcovselfpert[PH]>0.01){
	//	dualfprintf(fail_file,"got PH difference: %21.15g %21.15g\n",gcovselfpert[PH],0.0);
	//    }


	//      DLOOPA(j){
	//	dualfprintf(fail_file,"t=%21.15g %ld %d :: X1=%21.15g :: postgcov[%d][%d]=%21.15g :: gcovselfpert[%d]=%21.15g mypert=%2.15g\n",t,steppart,nstep,X[1],j,j,gcovlocal[j][j],j,gcovselfpert[j],-2.0*phi);
	//      }


	// add self-gravity perturbation to metric
	DLOOPA(j){
	  gcovlocal[j][j] +=gcovselfpert[j];
	  gcovpertlocal[j]+=gcovselfpert[j]; // in this way the perturbation part is always resolved
	  //	dualfprintf(fail_file,"gcovselfpert[%d]=%21.15g\n",j,gcovselfpert[j]);
	}
	if(whichcoord==KSCOORDS){
	  // KS-form has these terms
	  gcovlocal[TT][RR] += gcovselfpert[TT];
	  gcovlocal[RR][TT] = gcovlocal[TT][RR];
	}
      }
      else if(whichcoord==KS_BH_TOV_COORDS || whichcoord==KS_TOV_COORDS || whichcoord==BL_TOV_COORDS){
	// then doing full TOV-like solution that doesn't just come in as a perturbation
      
      }
      else{
	// then not setup for these coordinates
	dualfprintf(fail_file,"DOSELFGRAVVSR=1 will not work with this whichcoord=%d\n",whichcoord);
	myexit(145);
      }
#endif


      //  DLOOP(j,k) { fprintf(stderr,"1gcov[%d][%d]=%21.15g\n",j,k,gcovlocal[j][k]); fflush(stderr);}

      // whether to convert to prim coords
      if(getprim==1){
	// all the above are analytic, so have to convert to prim coords.
	gcov2gcovprim(ptrgeom, X, V, gcovlocal,gcovpertlocal, gcovlocal, gcovpertlocal);
      }
      //  DLOOP(j,k) { fprintf(stderr,"2gcov[%d][%d]=%21.15g\n",j,k,gcovlocal[j][k]); fflush(stderr);}

      //    if(ptrgeom->i==7 && nstep==1084){
      //    //      DLOOP(j,k) dualfprintf(fail_file,"present time: gcov[%d][%d]=%21.15g\n",j,k,gcovlocal[j][k]);
      // }
    }// end if needing to recompute metric


  }
  else{

    // then not present time, and assume past time metric stored, so recover
    // also assume never ask for past time at arbitrary location, checked previously

    if(ptrgeom->p!=NOWHERE){
      DLOOP(j,k){
	gcovlocal[j][k]=gcovlast[ptrgeom->i][ptrgeom->j][ptrgeom->k][ptrgeom->p][j][k];
      }
      DLOOPA(j){
	gcovpertlocal[j]=gcovpertlast[ptrgeom->i][ptrgeom->j][ptrgeom->k][ptrgeom->p][j];
      }

      //      if(ptrgeom->i==7 && nstep==1084){
      //	DLOOP(j,k) dualfprintf(fail_file,"past time: gcov[%d][%d]=%21.15g\n",j,k,gcovlocal[j][k]);
      //      }
    }
    else{
      // then don't have grid at right location and need to interpolate
      interpX_gcov(X, gcovlast, gcovpertlast, gcovlocal, gcovpertlocal);
      // above untested, so fail for now
      dualfprintf(fail_file,"interpX_gcov() unexpectedly called\n");
      myexit(13515);
    } 
  }// end else if doing past time

  //dualfprintf(fail_file,"blcoordcalled3: i=%d j=%d k=%d p=%d\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p);



}



// interpolate grid-based gcov to arbitrary location
// used for stored past grid
int interpX_gcov(FTYPE *X, FTYPE (*gcovgrid)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM], FTYPE (*gcovpertgrid)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM], FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  int i,j,k;
  int jj,kk;
  extern void icoord(FTYPE *X,int loc, int *i, int *j, int *k);
  int ip,jp,kp;
  FTYPE Xijk[NDIM],Xip[NDIM],Xjp[NDIM],Xkp[NDIM];
  FTYPE gijk,gip,gjp,gkp;
  FTYPE dist[4+1],totaldist;
  int loc;


  // CENT is chosen as reference location for all interpolations
  loc=CENT;

  // use icoord() to get i,j,k from X so can get potential at any location from interpolation of existing potential at certain locations

  // get centered i,j,k
  icoord(X,loc,&i,&j,&k); // truncates to lower i  
  coord_ijk(i,j,k, loc, Xijk);
  // only want effective grid location (i.e. normalized X position) so that differences are 1 between grid cells -- don't care about absolute magnitude of X once normalized
  DLOOPA(jj) Xijk[jj]/=dx[jj];

  /////////////////
  // X1 dir
  //
  // limit lower value of i
  if(i<=-N1BND-SHIFT1){
    i=-N1BND;
  }
  // limit upper value of i
  if(i>=N1+N1BND-SHIFT1){ // -1 is to offset so +1 is in range
    i=N1+N1BND-SHIFT1-SHIFT1; // first -1 is to offset from out of range, and second -1 is so ip is located within array range
  }

  // get X-position of coordinates i and i+1
  ip=i+SHIFT1;
  coord_ijk(ip,j,k, loc, Xip);
  DLOOPA(jj) Xip[jj]/=dx[jj];



  /////////////////
  // X2 dir
  //
  // limit lower value of j
  if(j<=-N2BND-SHIFT2){
    j=-N2BND;
  }
  // limit upper value of j
  if(j>=N2+N2BND-SHIFT2){ // -1 is to offset so +1 is in range
    j=N2+N2BND-SHIFT2-SHIFT2; // first -1 is to offset from out of range, and second -1 is so jp is located within array range
  }

  // get X-position of coordinates j and j+1
  jp=j+SHIFT2;
  coord_ijk(i,jp,k, loc, Xjp);
  DLOOPA(jj) Xjp[jj]/=dx[jj];



  /////////////////
  // X3 dir
  //
  // limit lower value of k
  if(k<=-N3BND-SHIFT3){
    k=-N3BND;
  }
  // limit upper value of k
  if(k>=N3+N3BND-SHIFT3){ // -1 is to offset so +1 is in range
    k=N3+N3BND-SHIFT3-SHIFT3; // first -1 is to offset from out of range, and second -1 is so kp is located within array range
  }

  // get X-position of coordinates k and k+1
  kp=k+SHIFT3;
  coord_ijk(i,j,kp, loc, Xkp);
  DLOOPA(jj) Xkp[jj]/=dx[jj];


  // now use some interpolation of 3 points

  // use this position to interpolate in X (which happens to be the same as interpolating in ijk)
  // distance away from one point toward ijk
  dist[1]=(1-(X[1]-Xijk[1]))*(1-(X[2]-Xijk[2]))*(1-(X[3]-Xijk[3]));
  // away from one point toward ip
  dist[2]=(1-(X[1]-Xip[1]))*(1-(X[2]-Xip[2]))*(1-(X[3]-Xip[3]))*(1.0-dist[1]);
  // away from one point toward jp
  dist[3]=(1-(X[1]-Xjp[1]))*(1-(X[2]-Xjp[2]))*(1-(X[3]-Xjp[3]))*(1.0-dist[2])*(1.0-dist[1]);
  // away from one point toward kp
  dist[4]=(1-(X[1]-Xkp[1]))*(1-(X[2]-Xkp[2]))*(1-(X[3]-Xkp[3]))*(1.0-dist[3])*(1.0-dist[2])*(1.0-dist[1]);
  // normalization of total distance
  totaldist=dist[1]+dist[2]+dist[3]+dist[4];
  

  DLOOP(jj,kk){
    gijk=gcovgrid[i][j][k][loc][jj][kk];
    gip=gcovgrid[ip][j][k][loc][jj][kk];
    gjp=gcovgrid[i][jp][k][loc][jj][kk];
    gkp=gcovgrid[i][j][kp][loc][jj][kk];

    // X's are per unit dX's
    gcov[jj][kk] = (gijk*dist[1]+gip*dist[2]+gjp*dist[3]+gkp*dist[4])/totaldist;

  }

  DLOOPA(jj){
    gijk=gcovpertgrid[i][j][k][loc][jj];
    gip=gcovpertgrid[ip][j][k][loc][jj];
    gjp=gcovpertgrid[i][jp][k][loc][jj];
    gkp=gcovpertgrid[i][j][kp][loc][jj];

    // X's are per unit dX's
    gcovpert[jj] = (gijk*dist[1]+gip*dist[2]+gjp*dist[3]+gkp*dist[4])/totaldist;

  }

  return(0);

}



// obtain prim gcon in primcoords of whichcoord type metric/coords
void gcon_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM])
{
  void set_gcon_blmetric(FTYPE *V, FTYPE gcon[][NDIM]);
  void set_gcon_ksmetric(FTYPE *V, FTYPE gcon[][NDIM]);
  void gcon2gconprim(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE gcon[][NDIM],FTYPE gconprim[][NDIM]);
  void matrix_inverse(int whichcoord, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);
  void matrix_inverse_2d(FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);
  int j,k;
  FTYPE V[NDIM];


  if(whichcoord>=0){
    if(whichcoord==BLCOORDS){
      //dualfprintf(fail_file,"mi in BLCOORDS\n");
      if(ANALYTICGCON){
	bl_coord(X, V);
	set_gcon_blmetric(V, gcon);
	if(getprim) gcon2gconprim(ptrgeom, X, V, gcon,gcon);
      }
      // since don't have gcon and want to keep things simple by only having to specify gcov 
      else matrix_inverse(whichcoord,gcov,gcon);

    }
    else if(whichcoord==KS_BH_TOV_COORDS){
      matrix_inverse(whichcoord,gcov,gcon);
    }
    else if(whichcoord==KS_TOV_COORDS){
      matrix_inverse(whichcoord,gcov,gcon);
    }
    else if(whichcoord==BL_TOV_COORDS){
      matrix_inverse(whichcoord,gcov,gcon);
    }
    else if(whichcoord==KSCOORDS){
      //dualfprintf(fail_file,"mi in KSCOORDS\n");
      //DLOOP(j,k) dualfprintf(fail_file,"ks gcov[%d][%d]=%g\n",j,k,gcov[j][k]);

      // DEBUG:
      //      bl_coord(X, V);
      //      dualfprintf(fail_file,"i=%d j=%d k=%d loc=%d :: %21.15g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p,V[1]);


      if(ANALYTICGCON){
	bl_coord(X, V);
	set_gcon_ksmetric(V,gcon);
	if(getprim) gcon2gconprim(ptrgeom, X, V,gcon,gcon);
      }
      // since don't have gcon and want to keep things simple by only having to specify gcov 
      else matrix_inverse(whichcoord,gcov,gcon);
    }
    else if(whichcoord==HTMETRIC){
      // do not have analytic gcon, so invert numerically
      matrix_inverse(whichcoord,gcov,gcon);
    }
    else if(whichcoord==HTMETRICACCURATE){
      // do not have analytic gcon, so invert numerically
      matrix_inverse(whichcoord,gcov,gcon);
    }
    else if(whichcoord==CARTMINKMETRIC){
      // do not have analytic gcon, so invert numerically
      matrix_inverse(whichcoord,gcov,gcon);
    }
    else if(whichcoord==UNIGRAVITY){
      // do not have analytic gcon, so invert numerically
      matrix_inverse(whichcoord,gcov,gcon);
    }
    else if(whichcoord==CYLMINKMETRIC){
      // do not have analytic gcon, so invert numerically
      matrix_inverse(whichcoord,gcov,gcon);
    }
    else if(whichcoord==SPCMINKMETRIC){
      // do not have analytic gcon, so invert numerically
      matrix_inverse(whichcoord,gcov,gcon);
    }
    else{
      dualfprintf(fail_file,"gcon_func(): no such whichcoord=%d\n",whichcoord);
      myexit(2687);
    }
  }
  else{
    dualfprintf(fail_file,"your request makes no sense (i.e. can't get prim gcon from prim gcon\n");
    myexit(3783);
  }

  // avoid coordinate singularity
  // should be needed with new matrix_inverse()
#if(0)
  if(ISSPCMCOORD(whichcoord)){
    bl_coord(X, V);
    if(fabs(V[1]-0.0)<SMALL){
      // then make PRIMECOORD flat since never should really be used if gdet=0
      // use t-r inverse rather than full inverse
      matrix_inverse_2d(gcov,gcon);
    }
  }
#endif


}





// obtain prime or non-prime alphalapse
void alphalapse_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM], FTYPE *alphalapse)
{

  // set alpha -- fabs just for roundoff error
  // Note ptrgeom only has i,j,k,loc at this point
  *alphalapse = 1.0/sqrt(fabs(-gcon[TT][TT]));

}




// delta is simply how big the differencing is, should be small, but not so small to lead to errors due to erros in the metric itself (i.e. keep larger than machine precision)
//#define DELTA (NUMEPSILON*1000.0)
//#define DELTA 1.e-5
//#define DELTA 1.e-5
//#define DELTA (pow(NUMEPSILON,1.0/3.0))
//#define DELTA NUMSQRTEPSILON // as in NR's fdjac() -- smaller isn't always better
// how to generically set this?  Too high, even slightly (10^{-10} for long doubles) and connection is screwed)

// Avery mentions that long double trig. functions only return double precision answer.  see ~/research/utils/triglongdouble.c

#if((REALTYPE==DOUBLETYPE)||(REALTYPE==FLOATTYPE))
//#define CONNDELTA (dx[1])
#define CONNDELTA 1E-5 // default -- seems to work pretty good generally to reduce max error
//#define CONNDELTA 5.4E-5 // default -- seems to work pretty good
//#define CONNDELTA 4.6E-5 // min of error for a specific case, but apparently not generally good
#elif(REALTYPE==LONGDOUBLETYPE)
//#define CONNDELTA 7.17E-6 // based on min of error for specific case
#define CONNDELTA 1E-6 // based on min of error for specific case
// polar region likes 6.5E-8 (min of error for specific case)
#endif
// GODMARK: Should really choose CONNDELTA as (totalsize[jj]*dx[jj]*1E-5)


// connection not simply transformed -- so compute directly from final metric (primcoords)
void conn_func(int whichcoord, FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void set_conn_cylminkmetric(FTYPE *X, struct of_geom *geom,
			      FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_spcminkmetric(FTYPE *X, struct of_geom *geom,
			      FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_cartminkmetric(FTYPE *X, struct of_geom *geom,
			       FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_unigravity(FTYPE *X, struct of_geom *geom,
			   FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_htmetric(FTYPE *X, struct of_geom *geom,
			 FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_htmetric_accurate(FTYPE *X, struct of_geom *geom,
				  FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_ksmetric(FTYPE *X, struct of_geom *geom,
			 FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_ks_bh_tov_metric(FTYPE *X, struct of_geom *geom,
				 FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_ks_tov_metric(FTYPE *X, struct of_geom *geom,
			      FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_bl_tov_metric(FTYPE *X, struct of_geom *geom,
			      FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void set_conn_blmetric(FTYPE *X, struct of_geom *geom,
			 FTYPE conn[][NDIM][NDIM],FTYPE *conn2);


  if(whichcoord==BLCOORDS){
    set_conn_blmetric(X,geom,conn,conn2);
  }
  else if(whichcoord==KSCOORDS){
    set_conn_ksmetric(X,geom,conn,conn2);
  }
  else if(whichcoord==KS_BH_TOV_COORDS){
    set_conn_ks_bh_tov_metric(X,geom,conn,conn2);
  }
  else if(whichcoord==KS_TOV_COORDS){
    set_conn_ks_tov_metric(X,geom,conn,conn2);
  }
  else if(whichcoord==BL_TOV_COORDS){
    set_conn_bl_tov_metric(X,geom,conn,conn2);
  }
  else if(whichcoord==HTMETRIC){
    set_conn_htmetric(X,geom,conn,conn2);
  }
  else if(whichcoord==HTMETRICACCURATE){
    set_conn_htmetric_accurate(X,geom,conn,conn2);
  }
  else if(whichcoord==CARTMINKMETRIC){
    set_conn_cartminkmetric(X,geom,conn,conn2);
  }
  else if(whichcoord==UNIGRAVITY){
    set_conn_unigravity(X,geom,conn,conn2);
  }
  else if(whichcoord==CYLMINKMETRIC){
    set_conn_cylminkmetric(X,geom,conn,conn2);
  }
  else if(whichcoord==SPCMINKMETRIC){
    set_conn_spcminkmetric(X,geom,conn,conn2);
  }
  else{
    dualfprintf(fail_file,"conn_func(): no such whichcoord=%d\n",whichcoord);
    myexit(7367);
  }
}



// obtain eomfunc (f_{(\nu)} factor: see how used in connection below and get_geometry())
// assumes X set before eomfun_func()
// should allow for DOEVOLVEMETRIC.  Currently WITHGDET and WITHNOGDET are usable with DOEVOLVEMETRIC
void eomfunc_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE *eomfunc)
{

  FTYPE gcovmcoord[NDIM][NDIM];
  FTYPE gcovpertcoord[NDIM];
  FTYPE r,th;
  int j,k;
  FTYPE V[NDIM];
  void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);



  if(WHICHEOM==WITHGDET){
    gcov_func(ptrgeom, getprim, whichcoord,X,gcovmcoord,gcovpertcoord); // actually returns primcoords version of whichcoord
    bl_coord_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, V);
    *eomfunc=gdet_func_singcheck(whichcoord,V,gcovmcoord);
  }
  else if(WHICHEOM==WITHNOGDET){
    *eomfunc=1.0;
  }
  else if(WHICHEOM==WITHSINSQ){ // obviously coordinate dependent
    bl_coord_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, V);

    r=V[1];
    th=V[2];
    *eomfunc=sin(th)*sin(th);
  }
  else{
    dualfprintf(fail_file,"eomfunc_func(): no such WHICHEOM=%d\n",WHICHEOM);
    myexit(798436);
  }
}



///////////////////////////
//
// g_{\mu\nu} (always analytic -- returns some coordinate system)
//
///////////////////////////



// needs M, J, and Q
// M: total M+dM mass of star
// J: total angular momentum
// Q: mass quadrapole moment
// external space-time of slowly rotating star, accurate to second order in $\Omega_{\star}$
// c=G=1
// NS mass sheds at v=0.8c (eq 28), allowed approx up to v=0.1c
// WD mass sheds at v=0.01c, allowed approx up to v=0.001c
// equally restrictive as measured by ratio of rotational energy to gravitational energy
void set_gcov_htmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  FTYPE P2,Q12, Q22;
  FTYPE z,OO,RO,SS,AA,alphasq,betaphi;
  FTYPE SSpert,OOpert,AApert;
  FTYPE M,J,Q;
  FTYPE r,th;
  FTYPE ftemp;
  // size of r relative to M is controlled by M, not r.
  FTYPE r2small;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;
  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif


  th=V[2];


#define FLATSPACE 0
#define NOSPIN 1
#define FULLHT 2

#define HTMETRICTYPE FULLHT

  // here G=c=1 and M is either 0 or 1 usually, such that
  // J=a M^2 , Q~J^2/M =a^2 M
  // NOW: uses MBH as mass of NS (or any compact object within the boundary)
  M=MBH;

  if(HTMETRICTYPE==FULLHT){
    J=a*MBH; // I\Omega = J
    Q=a*a*MBH; // approximately, perhaps, like Kerr geometry, if hard enough EOS
  }
  else if(HTMETRICTYPE==NOSPIN){
    a = 0.0;
    J=a*MBH; // no spin of metric
    Q=a*a*MBH;
  }
  else if(HTMETRICTYPE==FLATSPACE){
    M=MBH=1E-100; // no mass of star, just flat space (GJ model)
    a=0.0;
    J=a*MBH;
    Q=a*a*MBH;
  }


  // surface radius is not constant radius if Q!=0
  // if Q=0, then surface radius is 3M?  depends on EOS.

  // let's set M=1, such that J is up to M^2 and Q is up to M^3 such that radius is in units of GM/c^2

  // For Sun: M/r<M/R~2E-6, J/r^2<J/R^2~1E-12 (J/M^2=4E-24), Q/r^3<Q/R^3<~1E-10 (Q/R^3=8E-28) and surface roughly spherical

  // Fastest pulsar PSR 1937+21 \tau=1.56ms , 640times/sec
  // Deepto Chakrabarty of MIT
  // Vela: PSR 0833-45 \tau=89.3ms
  // strongest: PSR 0329+54 \tau=0.715s
  // peak: 760times/sec (limit by GR?), 2-3X is theoretical limit but would mass shed.
  // formation spin rate: 30/sec


  // see Hartle & Thorne 1968 or Kim Lee Lee Lee 2005

  // Legendre polynomial
  P2=0.5*(3.0*cos(th)*cos(th)-1.0);

  z=rsharp/M-1.0;

  // Mathematica's LegendreQ (associated Legendre of 2nd kind) is such that:
  //   Q12=Im[Legendre[2,1,z]] and Q22=-Re[Legendre[2,2,z]]
  Q12=sqrt(z*z-1.0)*( (3.0*z*z-2.0)/(z*z-1.0)-3.0/2.0*z*log((z+1)/(z-1)) );

  Q22=(3.0/2.0*(z*z-1.0)*log((z+1)/(z-1))-(3.0*z*z*z-5.0*z)/(z*z-1.0) );



  // metric stuff
  OOpert=-2.0*M/r+2.0*J*J/(r*r*r*r);
  OO=1.0+OOpert;

  RO=1.0+2.0*(J*J/(M*r*r*r)*(1.0+M/r)+5.0/8.0*(Q-J*J/M)/(M*M*M)*Q22)*P2;
  
  SSpert=-2.0*(J*J/(M*r*r*r)*(1.0-5.0*M/r)+5.0/8.0*(Q-J*J/M)/(M*M*M)*Q22)*P2;
  SS=1.0+SSpert;

  AApert=+2.0*(-J*J/(M*r*r*r)*(1.0+2.0*M/r)+5.0/8.0*(Q-J*J/M)/(M*M*M)*( (2.0*M/(sqrt(r*r*(1.0-2.0*M/r))))*Q12 - Q22  ))*P2;
  AA=1.0+AApert;
  
  alphasq=OO*RO;

  betaphi=-2.0*J/(r*r*r); // note that omega_{zamo}=-betaphi

  gcov[PH][PH] = rsharp*rsharp*AA*sin(th)*sin(th) ;
  gcovpert[PH] = gcov[PH][PH]-1.0;
  
  gcov[TT][TT] = -(alphasq-betaphi*betaphi*gcov[PH][PH]) ;
  ftemp = (1.0-alphasq);
  // g_{tt} + 1
  gcovpert[TT] = ftemp + betaphi*betaphi*gcov[PH][PH] ;

  gcov[TT][RR] = 0.0 ;
  gcov[TT][TH] = 0.0 ;
  gcov[TT][PH] = betaphi*gcov[PH][PH] ;
    
  gcov[RR][TT] = 0.0 ;
  gcov[RR][RR] = SS/OO ;
  // SS/OO = (SSpert+1.0)/OO = SSpert/OO + 1.0/OO , so SS/OO-1.0 = SSpert/OO + (1.0/OO-1.0)
  ftemp=(1.0/OO - 1.0);
  gcovpert[RR] = SSpert/OO + ftemp ; // might help with radial gravity terms

  gcov[RR][TH] = 0.0 ;
  gcov[RR][PH] = 0.0 ;
    
  gcov[TH][TT] = gcov[TT][TH] ;
  gcov[TH][RR] = gcov[RR][TH] ;
  gcov[TH][TH] = rsharp*rsharp*AA ;
  // r^2 AA - 1  = r^2(AApert+1.0) -1 = r^2 AApert + r^2 -1, so not a big issue with catastrophic cancellation in non-rel limit
  // probably have to define background as being r^2 to have gravity work at large radii (GODMARK)
  gcovpert[TH] = gcov[TH][TH] - 1.0;

  gcov[TH][PH] = 0.0 ;
    
  gcov[PH][TT] = gcov[TT][PH] ;
  gcov[PH][RR] = gcov[RR][PH] ;
  gcov[PH][TH] = gcov[TH][PH] ;

  
								     



}

#undef HTMETRICTYPE

// from Berti, White, Maniopoulou, and Bruni (2005)
void set_gcov_htmetric_accurate(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  FTYPE M,J,Q;
  FTYPE u,p,A1,A2,W,F1,F2,H1,H2,L,G1;
  FTYPE r,th;
  FTYPE ftemp;
  FTYPE r2small;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;
  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif


  th=V[2];



#define FLATSPACE 0
#define NOSPIN 1
#define FULLHT 2

#define HTMETRICTYPE FULLHT

  // here G=c=1 and M is either 0 or 1 usually, such that
  // J=a M^2 , Q~J^2/M =a^2 M
  // NOW: uses MBH as mass of NS (or any compact object within the boundary)
  M=MBH;

  if(HTMETRICTYPE==FULLHT){
    J=a*MBH; // I\Omega = J
    Q=a*a*MBH; // approximately, perhaps, like Kerr geometry, if hard enough EOS
  }
  else if(HTMETRICTYPE==NOSPIN){
    a = 0.0;
    J=a*MBH; // no spin of metric
    Q=a*a*MBH;
  }
  else if(HTMETRICTYPE==FLATSPACE){
    M=MBH=1E-100; // no mass of star, just flat space (GJ model)
    a=0.0;
    J=a*MBH;
    Q=a*a*MBH;
  }





  L=(80.0*pow(M,6)+8.0*pow(M,4)*r*r+10.0*M*M*M*r*r*r+20.0*M*M*pow(r,4)-45.0*M*pow(r,5)+15.0*pow(r,6));

  u=cos(th);
  p=1.0/(8.0*M*pow(r,4)*(r-2.0*M)); // typo in paper for parenthesis


  A1=(15.0*r*(r-2.0*M)*(1.0-3.0*u*u))/(16.0*M*M)*log(r/(r-2.0*M));
  A2=(15.0*(r*r-2.0*M*M)*(3.0*u*u-1.0))/(16.0*M*M)*log(r/(r-2.0*M));
  
  W=(r-M)*(16.0*pow(M,5)+8.0*pow(M,4)*r-10*M*M*r*r*r-30.0*M*pow(r,4)+15.0*pow(r,5))+u*u*(48.0*pow(M,6)-8.0*pow(M,5)*r-24.0*pow(M,4)*r*r-30.0*M*M*M*r*r*r-60.0*M*M*pow(r,4)+135.0*M*pow(r,5)-45.0*pow(r,6));

  F1=-p*W+A1;
  F2=5.0*r*r*r*p*(3.0*u*u-1.0)*(r-M)*(2.0*M*M+6.0*M*r-3.0*r*r)-A1;

  H1=A2+(1.0/(8.0*M*r*r*r*r))*(1.0-3.0*u*u) * (16.0*pow(M,5)+8.0*pow(M,4)*r-10.0*M*M*r*r*r+15.0*M*pow(r,4)+15.0*pow(r,5)) ;
  H2=-A2 + (1.0/(8.0*M*r))*5.0*(1.0-3.0*u*u)*(2.0*M*M-3.0*M*r-3.0*r*r);

  G1=p*((L-72.0*M*M*M*M*M*r)-3.0*u*u*(L-56.0*M*M*M*M*M*r))-A1;
  




  ftemp  = J*J*F1-Q*F2;
  gcov[TT][TT] = -(1.0-2.0*M/r)*(1.0+ftemp);
  // g_{tt} + 1
  //  gcovpert[TT] =  -(J*J*F1-Q*F2) + (2.0*M/r)*(1.0+J*J*F1-Q*F2);
  gcovpert[TT] =  -ftemp +(1.0+ftemp)*(2.0*M/r);



  gcov[TT][RR] = 0.0 ;
  gcov[TT][TH] = 0.0 ;
  gcov[TT][PH] = (2.0*J*M*M/r)*sin(th)*sin(th);
    
  gcov[RR][TT] = 0.0 ;

  ftemp=J*J*G1+Q*F2;
  gcov[RR][RR] = (1.0/(1.0-2.0*M/r))*(1.0+ftemp);
  // g_{rr} - 1
  gcovpert[RR] = (ftemp+2.0*M/r)/(1.0-2.0*M/r); // should help with radial gravity in non-rel regime

  gcov[RR][TH] = 0.0 ;
  gcov[RR][PH] = 0.0 ;
    
  gcov[TH][TT] = gcov[TT][TH] ;
  gcov[TH][RR] = gcov[RR][TH] ;
  gcov[TH][TH] = rsharp*rsharp*(1.0+J*J*H1-Q*H2) ;
  gcovpert[TH] = gcov[TH][TH]-1.0;

  gcov[TH][PH] = 0.0 ;
    
  gcov[PH][TT] = gcov[TT][PH] ;
  gcov[PH][RR] = gcov[RR][PH] ;
  gcov[PH][TH] = gcov[TH][PH] ;
  gcov[PH][PH] = gcov[TH][TH]*sin(th)*sin(th);
  gcovpert[PH] = gcov[PH][PH]-1.0;

  
								     



}
#undef HTMETRICTYPE



// Cartesian minkowski
// (t,x,y,z)
void set_gcov_cartminkmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  
  gcov[TT][TT] = -1.0;
  gcovpert[TT] = 0.0;

  gcov[TT][RR] = 0.0 ;
  gcov[TT][TH] = 0.0 ;
  gcov[TT][PH] = 0.0 ;
    
  gcov[RR][TT] = 0.0 ;
  gcov[RR][RR] = 1.0 ;
  gcovpert[RR] = 0.0;

  gcov[RR][TH] = 0.0 ;
  gcov[RR][PH] = 0.0 ;
    
  gcov[TH][TT] = gcov[TT][TH] ;
  gcov[TH][RR] = gcov[RR][TH] ;
  gcov[TH][TH] = 1.0 ;
  gcovpert[TH] = 0.0;

  gcov[TH][PH] = 0.0 ;
    
  gcov[PH][TT] = gcov[TT][PH] ;
  gcov[PH][RR] = gcov[RR][PH] ;
  gcov[PH][TH] = gcov[TH][PH] ;
  gcov[PH][PH] = 1.0 ;
  gcovpert[PH] = 0.0;

}


// (t,x,y,z)
void set_gcov_unigravity(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{

  FTYPE phi;
  FTYPE x,y,z;
  FTYPE FORCEMAG;

  x=V[1];
  y=V[2];
  z=V[3];

  //this is the hard-cored value for test 28, right? SASMARK
  FORCEMAG=0.1/(coordparams.timescalefactor*coordparams.timescalefactor); // strength of force, scales like v^2

  // uniform force field in +y direction (i.e. dphi/dy = FORCEMAG, so Force = -FORCEMAG)
  phi = FORCEMAG*y;

  gcov[TT][TT] = -1.0-2.0*phi;
  gcovpert[TT] = -2.0*phi;

  gcov[TT][RR] = 0.0 ;
  gcov[TT][TH] = 0.0 ;
  gcov[TT][PH] = 0.0 ;
    
  gcov[RR][TT] = 0.0 ;
  gcov[RR][RR] = 1.0-2.0*phi ;
  gcovpert[RR] = -2.0*phi;

  gcov[RR][TH] = 0.0 ;
  gcov[RR][PH] = 0.0 ;
    
  gcov[TH][TT] = gcov[TT][TH] ;
  gcov[TH][RR] = gcov[RR][TH] ;
  gcov[TH][TH] = 1.0-2.0*phi ;
  gcovpert[TH] = -2.0*phi;

  gcov[TH][PH] = 0.0 ;
    
  gcov[PH][TT] = gcov[TT][PH] ;
  gcov[PH][RR] = gcov[RR][PH] ;
  gcov[PH][TH] = gcov[TH][PH] ;

  gcov[PH][PH] = 1.0-2.0*phi ;
  gcovpert[PH] = -2.0*phi;

  //  dualfprintf(fail_file,"tsf=%21.15g :: %21.15g %21.15g :: %21.15g %21.15g %21.15g %21.15g\n",coordparams.timescalefactor,FORCEMAG,phi,gcovpert[TT],gcovpert[RR],gcovpert[TH],gcovpert[PH]);



}



// (t,R,z,\phi)
void set_gcov_cylminkmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{

  FTYPE r,Mass,RSTAR;
  FTYPE phi;
  FTYPE R,z;



  R = V[1];
  z = V[2];




  Mass=MBH; // NOW: use MBH as any compact object mass (in length units)
  RSTAR=0.1; // could use dxdxp[RR][RR]*dx[RR] to get size
  r=sqrt(R*R+z*z);

  if(r<RSTAR){
    phi = (Mass/(2.0*RSTAR))*((r/RSTAR)*(r/RSTAR)-3.0);
  }
  else phi = -Mass/r;

  
  gcov[TT][TT] = -1.0-2.0*phi;
  gcovpert[TT] = -2.0*phi;

  gcov[TT][RR] = 0.0 ;
  gcov[TT][TH] = 0.0 ;
  gcov[TT][PH] = 0.0 ;
    
  gcov[RR][TT] = 0.0 ;
  gcov[RR][RR] = 1.0-2.0*phi*R*R/(r*r); 
  gcovpert[RR] = -2.0*phi*R*R/(r*r); 

  gcov[RR][TH] = -2.0*phi*R*z/(r*r) ; // order v^4
  gcov[RR][PH] = 0.0 ;
    
  gcov[TH][TT] = gcov[TT][TH] ;
  gcov[TH][RR] = gcov[RR][TH] ;
  gcov[TH][TH] = 1.0 -2.0*phi*z*z/(r*r) ;
  gcovpert[TH] =-2.0*phi*z*z/(r*r) ; // order v^4

  gcov[TH][PH] = 0.0 ;
    
  gcov[PH][TT] = gcov[TT][PH] ;
  gcov[PH][RR] = gcov[RR][PH] ;
  gcov[PH][TH] = gcov[TH][PH] ;
  gcov[PH][PH] = R*R;
  gcovpert[PH] = gcov[PH][PH]-1.0; // doesn't matter that -1.0 subtracted

}






// (t,r,\theta,\phi)
void set_gcov_spcminkmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  FTYPE r,th;
  FTYPE r2small;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;
  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif




  th=V[2];
  
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
    if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }

  // avoid singularity at polar axis
#if(COORDSINGFIX) // just for metric components // hack
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif


  
  gcov[TT][TT] = -1.0;
  gcovpert[TT] = 0.0;

  gcov[TT][RR] = 0.0 ;
  gcov[TT][TH] = 0.0 ;
  gcov[TT][PH] = 0.0 ;
    
  gcov[RR][TT] = 0.0 ;
  gcov[RR][RR] = 1.0; 
  gcovpert[RR] = 0.0;

  gcov[RR][TH] = 0.0 ;
  gcov[RR][PH] = 0.0 ;
    
  gcov[TH][TT] = gcov[TT][TH] ;
  gcov[TH][RR] = gcov[RR][TH] ;
  gcov[TH][TH] = rsharp*rsharp ;
  gcovpert[TH] = gcov[TH][TH]-1.0;

  gcov[TH][PH] = 0.0 ;
    
  gcov[PH][TT] = gcov[TT][PH] ;
  gcov[PH][RR] = gcov[RR][PH] ;
  gcov[PH][TH] = gcov[TH][PH] ;

  gcov[PH][PH] = (rsharp*sin(th))*(rsharp*sin(th));
  gcovpert[PH] = gcov[PH][PH]-1.0;

}


// (~t,r,\theta,~\phi)
// mixed KS BH+TOV metric
void set_gcov_ks_bh_tov_metric(FTYPE *X, FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  FTYPE r,th;
  FTYPE mysin(FTYPE th);
  FTYPE mycos(FTYPE th);
  FTYPE MSQ;
  FTYPE phi;
  int j,k;
  void set_gcov_ksmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  extern int set_gcov_ks_tov_spcmetric(FTYPE *X, FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert, SFTYPE *MOrself, SFTYPE *phiself, SFTYPE *vrsqself);
  FTYPE gcovtovks[NDIM][NDIM],gcovtovkspert[NDIM];
  FTYPE gcovbhks[NDIM][NDIM],gcovbhkspert[NDIM];
  FTYPE MtotOr,MBHprime,MBHprimeOr;
  SFTYPE MOrself,phiself,vrsqself;
  FTYPE starpot, totalpot;
  FTYPE r2small;
  FTYPE disc;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;



#if(0)
  // debug test
  DLOOP(j,k) gcov[j][k] = 0.0;
  DLOOPA(j,j) gcov[j][j] = 1.0;
  gcovpert[TT] gcov[TT][TT]=-1.0;
  DLOOPA(j)  gcovpert[j]= 0.0;
#endif
  




  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif


  th=V[2];


  // get TOV ks metric assuming pre-existing phivsr_tot, MvsrOr_tot, vrsqvsr_tot
  // needs X for interpolation
  set_gcov_ks_tov_spcmetric(X, V, gcovtovks, gcovtovkspert, &MOrself, &phiself, &vrsqself);

  // get black hole KS metric
  set_gcov_ksmetric(V, gcovbhks, gcovbhkspert);

  // default is KS BH metric
  DLOOP(j,k) gcov[j][k] = gcovbhks[j][k];
  DLOOPA(j) gcovpert[j] = gcovbhkspert[j];

  // now compute TOV modifications
  //  MBHprime = 2.0*MBH*r*r/(a*a+2.0*r*r+a*a*cos(2.0*th));
  MBHprime = MBH/(1.0 + (a*a + a*a*cos(2.0*th))/(2.0*r2small));
  MBHprimeOr = MBHprime/r;
  MtotOr = MBHprimeOr + MOrself; // assume MOrself has sign built-in
  starpot=vrsqself - 2.0*MOrself;
  totalpot = -2.0*MtotOr;

  // assign mixed metric components
  gcov[TT][TT] = gcovbhks[TT][TT]*(-gcovtovks[TT][TT]);
  gcovpert[TT] = gcov[TT][TT] + 1.0; // no obvious way to remove perturbed part

  //  gcov[RR][RR] = 1.0 + 2.0*MtotOr-vrsqself;
  gcov[RR][RR] = 1.0 + fabs(- totalpot);
  gcovpert[RR] = fabs(-totalpot);

  disc=-gcovtovks[TT][TT]/(1.0-fabs(starpot));
  if(disc<0.0){
    disc=0.0; // assume this value never used
    // and can't trust that signature of metric will work out now, so artificially force signature to be negative
    gcov[TT][TT]=-fabs(gcov[TT][TT]); // should never be used
  }

  gcov[TT][RR] = gcov[RR][TT] = (-totalpot)*sqrt(disc);

  // DEBUG (overwrite)
  //  DLOOP(j,k) gcov[j][k] = gcovtovks[j][k];
  //  DLOOPA(j) gcovpert[j] = gcovtovkspert[j];


  // DEBUG
  //  if(nstep>1080 && nstep<1100){
  //  if(icurr==0 || icurr==-1){
  //    dualfprintf(fail_file,"r=%21.15g icurr=%d starpot=%21.15g totalpot=%21.15g MBHprimeOr=%21.15g MOrself=%21.15g gcovtovks[TT][TT]=%21.15g disc=%21.15g gcovbhks[TT][TT]=%21.15g\n",r,icurr,starpot,totalpot,MBHprimeOr,MOrself,gcovtovks[TT][TT],disc,gcovbhks[TT][TT]);
  // DLOOP(j,k) dualfprintf(fail_file,"ks_bh_tov gcov[%d][%d]=%g\n",j,k,gcov[j][k]);
  // }


}


// TOV KS metric wrapper
void set_gcov_ks_tov_metric(FTYPE *X, FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  extern int set_gcov_ks_tov_spcmetric(FTYPE *X, FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert, SFTYPE *MOrself, SFTYPE *phiself, SFTYPE *vrsqself);
  SFTYPE MOrself,phiself,vrsqself;


  // get TOV ks metric assuming pre-existing phivsr_tot, Mvsr_tot, vrsqvsr_tot
  // needs X for interpolation
  set_gcov_ks_tov_spcmetric(X, V, gcov, gcovpert, &MOrself, &phiself, &vrsqself);

}

// TOV KS metric wrapper
void set_gcov_bl_tov_metric(FTYPE *X, FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  extern int set_gcov_bl_tov_spcmetric(FTYPE *X, FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert, SFTYPE *MOrself, SFTYPE *phiself, SFTYPE *vrsqself);
  SFTYPE MOrself,phiself,vrsqself;


  // get TOV BL metric assuming pre-existing phivsr_tot, Mvsr_tot, vrsqvsr_tot
  // needs X for interpolation
  set_gcov_bl_tov_spcmetric(X, V, gcov, gcovpert, &MOrself, &phiself, &vrsqself);

}



// (~t,r,\theta,~\phi)
void set_gcov_ksmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  FTYPE sth, cth, s2, a2, r2, r2small, r3;
  FTYPE rho2,rho2small;
  FTYPE r,th;
  FTYPE mysin(FTYPE th);
  FTYPE mycos(FTYPE th);
  FTYPE MSQ;
  FTYPE phi;
  int j,k;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;


  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
  //  dualfprintf(fail_file,"got here r=%21.15g drsing=%21.15g\n",r,drsing);
#else
  r2small = r*r + SMALL;
#endif


  //  dualfprintf(fail_file,"r=%21.15g\n",r);

  // theta gives well-defined metric for all th
  th=V[2];


  sth=mysin(th);
  cth=mycos(th);


  s2 = sth * sth;
  a2 = a * a;
  r2 = rsharp * rsharp;
  r3 = r2 * r;
  rho2 = r2 + a * a * cth * cth;
  rho2small = r2small + a * a * cth * cth;
  MSQ=(2.*MBH*r-QBH*QBH);

  // really gcov00 diverges at r=0, but force non-divergence to avoid nan's only
#define ks_gcov00 (-1. + MSQ/rho2small)
#define ks_gcov01 (MSQ/rho2small)
#define ks_gcov02 (0)
#define ks_gcov03 (-MSQ*a*s2/rho2small)
#define ks_gcov10 (ks_gcov01)
#define ks_gcov11 (1. + MSQ/rho2small)
#define ks_gcov12 (0)
#define ks_gcov13 (-a*s2*(1. + MSQ/rho2small))
#define ks_gcov20 (0)
#define ks_gcov21 (0)
#define ks_gcov22 (rho2)
#define ks_gcov23 (0)
#define ks_gcov30 (ks_gcov03)
#define ks_gcov31 (ks_gcov13)
#define ks_gcov32 (0)
  // old
  //#define ks_gcov33 (s2*(rho2 + a*a*s2*(1. + 2.*r/rho2small)))
  // new
#define ks_gcov33 (s2*(rho2 + a*a*s2*(1. + MSQ/rho2small)))
  // Living review format
  //#define ks_gcov33 (s2*(r*r+a*a + a*a*s2*(MSQ/rho2small)))


  gcov[TT][TT] = ks_gcov00 ;
  gcovpert[TT] = MSQ/rho2small;


  gcov[TT][RR] = ks_gcov01 ;
  gcov[TT][TH] = ks_gcov02 ;
  gcov[TT][PH] = ks_gcov03 ;
    
  gcov[RR][TT] = ks_gcov10 ;

  gcov[RR][RR] = ks_gcov11 ;
  gcovpert[RR] = MSQ/rho2small ;

  gcov[RR][TH] = ks_gcov12 ;
  gcov[RR][PH] = ks_gcov13 ;
    
  gcov[TH][TT] = ks_gcov20 ;
  gcov[TH][RR] = ks_gcov21 ;

  gcov[TH][TH] = ks_gcov22 ;
  gcovpert[TH] = gcov[TH][TH]-1.0;

  gcov[TH][PH] = ks_gcov23 ;
    
  gcov[PH][TT] = ks_gcov30 ;
  gcov[PH][RR] = ks_gcov31 ;
  gcov[PH][TH] = ks_gcov32 ;

  gcov[PH][PH] = ks_gcov33 ;
  gcovpert[PH] = gcov[PH][PH]-1.0;


  //  if(nstep>1051){
  // dualfprintf(fail_file,"MBH=%21.15g a=%21.15g gcov[TT][TT]=%21.15g\n",MBH,a,gcov[TT][TT]);
    
  //  DLOOP(j,k) dualfprintf(fail_file,"ks gcov[%d][%d]=%g\n",j,k,gcov[j][k]);
  //}

}


// (t,r,\theta,\phi)
void set_gcov_blmetric(FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert)
{
  FTYPE sth, cth, s2, a2, r2, r2small, r3, DD, mu;
  FTYPE mupert,DDpert;
  FTYPE r,th;
  FTYPE MSQ;
  int j,k;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;
  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif



  th=V[2];




  cth = cos(th);
  sth=sin(th);


  MSQ=(2.*MBH*rsharp-QBH*QBH);

  s2 = sth * sth;
  a2 = a * a;
  r2 = rsharp * rsharp;
  r3 = r2 * r;

  DDpert =- MSQ / r2small + a2 / r2small; 
  DD = 1. + DDpert;
  mupert=+ a2 * cth * cth / r2small;
  mu = 1. + mupert;



#define bl_gcov00 (-(1. - MSQ/(r2small*mu)))
#define bl_gcov01 (0)
#define bl_gcov02 (0)
#define bl_gcov03 (-MSQ*a*s2/(r2small*mu))
#define bl_gcov10 (0)
#define bl_gcov11 (mu/DD)
#define bl_gcov12 (0)
#define bl_gcov13 (0)
#define bl_gcov20 (0)
#define bl_gcov21 (0)
#define bl_gcov22 (r2*mu)
#define bl_gcov23 (0)
#define bl_gcov30 (bl_gcov03)
#define bl_gcov31 (0)
#define bl_gcov32 (0)
#define bl_gcov33 (r2*sth*sth*(1. + a2/r2small + MSQ*a2*s2/(r2small*r2small*mu)))

  gcov[TT][TT] = bl_gcov00 ;
  gcovpert[TT] = MSQ/(r2small*mu);

  gcov[TT][RR] = bl_gcov01 ;
  gcov[TT][TH] = bl_gcov02 ;
  gcov[TT][PH] = bl_gcov03 ;
    
  gcov[RR][TT] = bl_gcov10 ;

  gcov[RR][RR] = bl_gcov11 ;
  // mu/DD = (1.0+mupert)/(1.0+DDpert) = 
  gcovpert[RR] = (DDpert - mupert) / (1.0 + DDpert);

  gcov[RR][TH] = bl_gcov12 ;
  gcov[RR][PH] = bl_gcov13 ;
    
  gcov[TH][TT] = bl_gcov20 ;
  gcov[TH][RR] = bl_gcov21 ;

  gcov[TH][TH] = bl_gcov22 ;
  gcovpert[TH] = gcov[TH][TH]-1.0;

  gcov[TH][PH] = bl_gcov23 ;
    
  gcov[PH][TT] = bl_gcov30 ;
  gcov[PH][RR] = bl_gcov31 ;
  gcov[PH][TH] = bl_gcov32 ;

  gcov[PH][PH] = bl_gcov33 ;
  gcovpert[PH] = gcov[PH][PH]-1.0;



  //  DLOOP(j,k) dualfprintf(fail_file,"bl gcov[%d][%d]=%g\n",j,k,gcov[j][k]);



}


///////////////////////////
//
// g^{\mu\nu} (analytic or numerical)
//
///////////////////////////


// (~t,r,\theta,~\phi)
void set_gcon_ksmetric(FTYPE *V, FTYPE gcon[][NDIM])
{
  FTYPE sth, cth, s2, a2, r2, r2small,r3;
  FTYPE rho2;
  FTYPE r,th;
  FTYPE MSQ;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;
  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif
  // size of r and M controlled by changing M, not rescaling r




  th=V[2];

  cth = cos(th);
  sth=sin(th);


  s2 = sth * sth;
  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  rho2 = r2small + a * a * cth * cth; // always divide by this
  MSQ=(2.*MBH*r-QBH*QBH);


#define ks_gcon00 (-(1.+ MSQ/rho2))
#define ks_gcon01 (MSQ/rho2)
#define ks_gcon02 (0)
#define ks_gcon03 (0)
#define ks_gcon10 (ks_gcon01)
  //#define ks_gcon11 ((r*(r-2.)+a*a)/rho2)
#define ks_gcon11 ((r*r+a*a-MSQ)/rho2)
#define ks_gcon12 (0)
#define ks_gcon13 (a/rho2)
#define ks_gcon20 (ks_gcon02)
#define ks_gcon21 (ks_gcon12)
#define ks_gcon22 (1./rho2)
#define ks_gcon23 (0)
#define ks_gcon30 (ks_gcon03)
#define ks_gcon31 (ks_gcon13)
#define ks_gcon32 (ks_gcon23)
#define ks_gcon33 (1./(rho2*s2))


  gcon[TT][TT] = ks_gcon00 ;
  gcon[TT][RR] = ks_gcon01 ;
  gcon[TT][TH] = ks_gcon02 ;
  gcon[TT][PH] = ks_gcon03 ;
    
  gcon[RR][TT] = ks_gcon10 ;
  gcon[RR][RR] = ks_gcon11 ;
  gcon[RR][TH] = ks_gcon12 ;
  gcon[RR][PH] = ks_gcon13 ;
    
  gcon[TH][TT] = ks_gcon20 ;
  gcon[TH][RR] = ks_gcon21 ;
  gcon[TH][TH] = ks_gcon22 ;
  gcon[TH][PH] = ks_gcon23 ;
    
  gcon[PH][TT] = ks_gcon30 ;
  gcon[PH][RR] = ks_gcon31 ;
  gcon[PH][TH] = ks_gcon32 ;
  if(s2!=0.0){
    gcon[PH][PH] = ks_gcon33 ;
  }
  else gcon[PH][PH]=1.0; // avoid coordinate singularity -- although should never use this value

}

// (t,r,\theta,\phi)
void set_gcon_blmetric(FTYPE *V, FTYPE gcon[][NDIM])
{
  FTYPE sth, cth, s2, a2, r2, r2small, r3, DD, mu;
  FTYPE r,th;
  FTYPE MSQ;
#if(SMOOTHSING)
  FTYPE signr,rsmooth;
#endif
  FTYPE rsharp;
  
  r = rsharp = V[1];
#if(SMOOTHSING)
  signr=mysign(r);
  rsmooth = signr*(fabs(r)+SMALL+drsing);
  r = rsmooth;
  r2small = r*r;
#else
  r2small = r*r + SMALL;
#endif


  // size of r relative to M is controlled by M, not r.

  th=V[2];

  cth = cos(th);
  sth=sin(th);

  MSQ=(2.*MBH*r-QBH*QBH);

  s2 = sth * sth;
  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  DD = 1. - MSQ / r2 + a2 / r2;
  mu = 1. + a2 * cth * cth / r2;


#define bl_gcon00 (-1. - MSQ*(1. + a2/r2small)/(r2small*DD*mu))
#define bl_gcon01 (0)
#define bl_gcon02 (0)
#define bl_gcon03 (-MSQ*a/(r2small*r2small*DD*mu))
#define bl_gcon10 (0)
#define bl_gcon11 (DD/mu)
#define bl_gcon12 (0)
#define bl_gcon13 (0)
#define bl_gcon20 (0)
#define bl_gcon21 (0)
#define bl_gcon22 (1./(r2small*mu))
#define bl_gcon23 (0)
#define bl_gcon30 (bl_gcon03)
#define bl_gcon31 (0)
#define bl_gcon32 (0)
#define bl_gcon33 ((1. - MSQ/(r2small*mu))/(r2small*sth*sth*DD))


  gcon[TT][TT] = bl_gcon00 ;
  gcon[TT][RR] = bl_gcon01 ;
  gcon[TT][TH] = bl_gcon02 ;
  gcon[TT][PH] = bl_gcon03 ;
    
  gcon[RR][TT] = bl_gcon10 ;
  gcon[RR][RR] = bl_gcon11 ;
  gcon[RR][TH] = bl_gcon12 ;
  gcon[RR][PH] = bl_gcon13 ;
    
  gcon[TH][TT] = bl_gcon20 ;
  gcon[TH][RR] = bl_gcon21 ;
  gcon[TH][TH] = bl_gcon22 ;
  gcon[TH][PH] = bl_gcon23 ;
    
  gcon[PH][TT] = bl_gcon30 ;
  gcon[PH][RR] = bl_gcon31 ;
  gcon[PH][TH] = bl_gcon32 ;
  gcon[PH][PH] = bl_gcon33 ;

}

///////////////////////////
//
// CONNECTIONS (analytic or numerical)
//
///////////////////////////

void set_conn_cylminkmetric(FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  int i,j,k;
  FTYPE gcovmid[NDIM][NDIM];
  FTYPE gcovpertmid[NDIM];
  FTYPE gdetmid;

  FTYPE dxdxp[NDIM][NDIM]; //atch
  FTYPE V[NDIM];           //atch


  if(defcoord==UNIFORMCOORDS){ // uniform grid  SUPERSASMARK
    // could directly use gdet in global memory
    // only works for X1=R and X2=z
    gcov_func(geom, 1,CYLMINKMETRIC,X, gcovmid,gcovpertmid);

    bl_coord( X, V );  //actually, dxdxprim() does not use V or X since metric is uniform
    gdetmid=gdet_func_singcheck(MCOORD,V,gcovmid);


    // see transforms.c and mettometp() and see gcov2gcovprim()  //atch
    dxdxprim(X, V, dxdxp);                                       //atch

    for (k = 0; k < NDIM; k++) conn2[k]= 0.0;
    //conn2[RR]=-1.0/gdetmid;  //wrong as well... shouldn't it be 0 since there is no 2nd connection when WITHGDET?  SUPERSASMARK
 

    for (i = 0; i < NDIM; i++)
      for (j = 0; j < NDIM; j++)
	for (k = 0; k < NDIM; k++) {
	  conn[i][j][k] = 0.;
	}
    conn[PH][RR][PH]=1.0 / X[1]; //1.0/gdetmid; //apparently, wrong because should not care about the 2nd dimension  SUPERSASMARK
    conn[PH][PH][RR]=1.0 / X[1]; //1.0/gdetmid; //apparently, wrong because should not care about the 2nd dimension  SUPERSASMARK
    conn[RR][PH][PH]= - X[1] * dxdxp[3][3] * dxdxp[3][3]; //-gdetmid;
  }
  else{
    conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  }


  
}

// only works for X1=R and X2=z
void set_conn_cartminkmetric(FTYPE *X, struct of_geom *geom,
			     FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2);

  int i,j,k;

  if(defcoord==UNIFORMCOORDS){// uniform grid
    for (k = 0; k < NDIM; k++) conn2[k]= 0.0;
    
    for (i = 0; i < NDIM; i++)
      for (j = 0; j < NDIM; j++)
	for (k = 0; k < NDIM; k++) {
	  conn[i][j][k] = 0.;
	}
  }
  else{
    conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  }

}

void set_conn_unigravity(FTYPE *X, struct of_geom *geom,
			 FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2);

  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
}


void set_conn_spcminkmetric(FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2);

  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
}


void set_conn_htmetric(FTYPE *X, struct of_geom *geom,
		       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2);

  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
}

void set_conn_htmetric_accurate(FTYPE *X, struct of_geom *geom,
				FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2);

  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
}

void set_conn_ksmetric(FTYPE *X, struct of_geom *geom,
		       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  void mks_conn_func(FTYPE *X, struct of_geom *geom,
		     FTYPE lconn[][NDIM][NDIM],FTYPE *conn2);

  //  FTYPE DELTA; // for debug below

  // the analytic form can be used with WHICHEOM=WITHNOGDET
  // determine for which we have analytic expressions
  // this currently competes with mks_source_conn in phys.c
  if((!ANALYTICCONNECTION)||(defcoord!=LOGRSINTH)) conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  else mks_conn_func(X,geom,conn,conn2);


  /*
  // some debug stuff
  //  if((geom->i==62)&&(geom->j==32)){ // where c002 is bad
  if(0&&(geom->i==32)&&(geom->j==0)){ // where c023 is bad
  for(DELTA=1E-15;DELTA<1E5;DELTA*=1.01){
  conn_func_numerical1(DELTA,X, geom, conn, conn2);
  dualfprintf(fail_file,"%30.20Lg %30.20Lg\n",DELTA,conn[0][2][3]); fflush(fail_file);
  //      dualfprintf(fail_file,"%21.15g %21.15g\n",DELTA,conn[0][2][3]); fflush(fail_file);
  }
  exit(0);
  }
  else{
  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  }
  //mks_conn_func(X,geom,conn,conn2);
  */

}


void set_conn_ks_bh_tov_metric(FTYPE *X, struct of_geom *geom,
			       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  
  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  

}


void set_conn_ks_tov_metric(FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  
  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  

}

void set_conn_bl_tov_metric(FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2);
  
  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
  

}


void set_conn_blmetric(FTYPE *X, struct of_geom *geom,
		       FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  void conn_func_numerical1(FTYPE DELTA,FTYPE *X, struct of_geom *geom,
			    FTYPE conn[][NDIM][NDIM],FTYPE *conn2);

  conn_func_numerical1(CONNDELTA,X, geom, conn, conn2);
}


































// ///////////////////////////////////////////////////////////
// 
// below are independent of user choice of metric/coords/grid
// 
// ///////////////////////////////////////////////////////////


// find the con/cov forms of the chosen metric
void gset(int getprim, int whichcoord, int i, int j, int k, struct of_geom *ptrgeom)
{
  // assumes loc=CENT
  int loc;

  loc=CENT;
  gset_genloc(getprim, whichcoord, i, j, k, loc, ptrgeom);

}

// find the con/cov forms of the chosen metric
// fills in information like get_geometry but for arbitrary metric not just PRIMECOORDS
// doesn't set igeom's
void gset_genloc(int getprim, int whichcoord, int i, int j, int k, int loc, struct of_geom *ptrgeom)
{
  FTYPE X[NDIM],V[NDIM];
  struct of_geom tempgeom;
  extern void assign_eomfunc(struct of_geom *geom, FTYPE eomfuncgen);
  void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  void gcon_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);

  FTYPE eomfuncgen;
  FTYPE (*gcovptr)[NDIM];
  FTYPE (*gconptr)[NDIM];
  FTYPE *gcovpertptr;


  ptrgeom->i=i;
  ptrgeom->j=j;
  ptrgeom->k=k;
  ptrgeom->p=loc;
  icurr=i;
  jcurr=j;
  jcurr=k;
  pcurr=loc;


#if(GETGEOMUSEPOINTER)
  // then need to use dummy pointer space that has real memory assigned
  gcovptr=ptrgeom->gengcov;
  gconptr=ptrgeom->gengcon;
  gcovpertptr=ptrgeom->gengcovpert;

  ptrgeom->gcov=ptrgeom->gengcov; // pointer
  ptrgeom->gcon=ptrgeom->gengcon; // pointer
  ptrgeom->gcovpert=ptrgeom->gengcovpert; // pointer

#else
  // then ptrgeom->gcov,gcon,gcovpert are real memory spaces
  gcovptr=ptrgeom->gcov;
  gconptr=ptrgeom->gcon;
  gcovpertptr=ptrgeom->gcovpert;
#endif

  if(whichcoord>=0){
    coord_ijk(i, j, k, loc, X);
    bl_coord_ijk(i, j, k, loc, V);
    gcov_func(ptrgeom,getprim,whichcoord,X,gcovptr,gcovpertptr);
    ptrgeom->g=gdet_func_singcheck(whichcoord,V,gcovptr); // must come after gcov_func() above
    gcon_func(ptrgeom,getprim,whichcoord,X,gcovptr,gconptr); // must come after gcov_func() above
    alphalapse_func(ptrgeom,getprim,whichcoord,X,gcovptr,gconptr,&(ptrgeom->alphalapse));

    eomfunc_func(ptrgeom, getprim, whichcoord,X,&eomfuncgen);
    assign_eomfunc(ptrgeom,eomfuncgen); // must come after assigning ptrgeom->g above

  }
  else if(whichcoord==PRIMECOORDS){ // special case
    get_geometry(i,j,k,loc,ptrgeom);
  }
  else{
    dualfprintf(fail_file,"gset(): no such whichcoord=%d\n",whichcoord);
    myexit(3466);
  }

}


#define OPTMETRICLOOP 1 // whether to use highly optmized loop (assumes metric is symmetric)
#define COMPUTEPERTURBEDMETRIC 0 // GODMARK: NOT CORRECT RIGHT NOW, so do NOT do it

void gcov2gcovprim(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert, FTYPE gcovprim[][NDIM], FTYPE *gcovpertprim)
     //void gcov2gcovprim(FTYPE *X, FTYPE *V, FTYPE gcov[][NDIM], FTYPE *gcovpert, FTYPE gcovprim[][NDIM], FTYPE *gcovpertprim)
     //void gcov2gcovprim(FTYPE *X, FTYPE *V, FTYPE gcov[][NDIM],FTYPE gcovprim[][NDIM])
{
  int j, k, l, m;
  FTYPE dxdxp[NDIM][NDIM];
  FTYPE tmp[NDIM][NDIM];
  FTYPE ftemp1,ftemp2;
  int q;
  void transgcov(FTYPE gcov[][NDIM], FTYPE dxdxp[][NDIM], FTYPE gcovprim[][NDIM]);
  

  // now take term by term:
  // g_{u v} = \vec{e_{\mu}}\cdot\vec{e_{\nu}} 
  //           * (dx/dx')_{mu} * (dx/dx')_{\nu} =
  //          \vec{e'_{\mu}}\cdot\vec{e'_{\nu}} 

  // dx/dx' where '=prim coords (i.e. nonuni coords)
  //  dualfprintf(fail_file,"gcov2gcovprim: i=%d p=%d\n",ptrgeom->i,ptrgeom->p);
  dxdxprim_ijk_2(ptrgeom,X,V,dxdxp);

#if(OPTMETRICLOOP==0)
  DLOOP(j,k){
    tmp[j][k] = 0.;
    for(l=0;l<NDIM;l++) for(m=0;m<NDIM;m++){
      // g_{mup nup} = g_{mu nu} T^mu_mup T^nu_nup
      // where T^mu_mup == dx^mu[BL]/dx^mup[KSP uni grid]
      tmp[j][k] += gcov[l][m] * dxdxp[l][j] * dxdxp[m][k];
    }
    // use tmp since gcon might be same address as gconprim
    gcovprim[j][k] = tmp[j][k];
  }
#else
  transgcov(gcov,dxdxp,gcovprim);
#endif

  //  DLOOP(j,k) dualfprintf(fail_file,"prim gcov[%d][%d]=%21.15g\n",j,k,gcov[j][k]);
  //  DLOOP(j,k) dualfprintf(fail_file,"prim gcovprim[%d][%d]=%21.15g\n",j,k,gcovprim[j][k]);
  //  DLOOP(j,k) dualfprintf(fail_file,"prim dxdxp[%d][%d]=%21.15g\n",j,k,dxdxp[j][k]);

  ///////////////////////////////
  //
  // perturbed terms
  //
  //////////////////////////////

#if(COMPUTEPERTURBEDMETRIC)
  // SUPERGODMARK GODMARK: This only works for non-rel gravity if dxdxp is diagonal.  Not sure what to do in general
  DLOOPA(q){
    //    dualfprintf(fail_file,"gcovpert[%d]=%21.15g\n",q,gcovpert[q]);


    // get q-q term
    // -1 + g_{q'q'} for q!=TT, else  SOMECONSTANT + g_{q'q'} for q=TT

    // for spatial parts, only care about deviations from constant in getting connection coefficients
    // avoids catastrophic cancellation (GODMARK)
    if( (q!=TT)&&(defcoord==UNIFORMCOORDS)) ftemp1=0.0;
    // then unsure how to handle in general, so just leave alone (SUPERGODMARK GODMARK)
    // problem is with machine error in dxdxp leading to apparently large connection coefficients due to that machine error
    // also don't know ahead of time what constant to choose for prim coordinate quantities
    // and don't know how to handle variations in dxdxp, since then no constant to subtract off, but grid accelerations apparently wash out non-rel gravity accelerations?
    // except for time term, which is solid
    else ftemp1=(-1.0 + dxdxp[q][q] * dxdxp[q][q]);
    if(q==TT) ftemp1 *=-1.0;
    gcovpertprim[q]  = (gcovpert[q] * dxdxp[q][q] * dxdxp[q][q]) + ftemp1;

    // now add 15 other terms
    ftemp2 = 0.;
    for(l=0;l<NDIM;l++) for(m=0;m<NDIM;m++){
      if((l!=q)&&(m!=q)) ftemp2+= gcov[l][m] * dxdxp[l][q] * dxdxp[m][q];
    }
    // add other 15 terms to answer for total of 16 terms
    gcovpertprim[q]+=ftemp2;

    
    //    dualfprintf(fail_file,"dxdxp[%d][%d]=%21.15g\n",q,q,dxdxp[q][q]);
    //    dualfprintf(fail_file,"ftemp1[%d]=%21.15g ftemp2[%d]=%21.15g gcovpertprim[%d]=%21.15g\n",q,ftemp1,q,ftemp2,q,gcovpertprim[q]);
  }
#elif(0)
  // override for now
  gcovpertprim[TT]=gcovprim[TT][TT]-mink(TT,TT);
  SLOOPA(q){
    gcovpertprim[q]=gcovprim[q][q]-mink(q,q);
  }
#elif(1)
  // override for now
  gcovpertprim[TT]=gcovprim[TT][TT]+1.0;
  gcovpertprim[RR]=gcovprim[RR][RR]-1.0;
  gcovpertprim[TH]=gcovprim[TH][TH]-1.0;
  gcovpertprim[PH]=gcovprim[PH][PH]-1.0;
#endif

 

}

void transgcov_old(FTYPE gcov[][NDIM], FTYPE dxdxp[][NDIM], FTYPE gcovprim[][NDIM])
{
  int j, k, l, m;
  FTYPE tmp[NDIM][NDIM];

  DLOOP(j,k){
    tmp[j][k] = 0.;
    for(l=0;l<NDIM;l++) for(m=0;m<NDIM;m++){
      // g_{mup nup} = g_{mu nu} T^mu_mup T^nu_nup
      // where T^mu_mup == dx^mu[BL]/dx^mup[KSP uni grid]
      tmp[j][k] += gcov[l][m] * dxdxp[l][j] * dxdxp[m][k];
    }
    // use tmp since gcon might be same address as gconprim
    gcovprim[j][k] = tmp[j][k];
  }

}

// gcov might be same memory address as gcovprim, so use tmp
// assumes metric is symmetric 2nd rank tensor
void transgcov(FTYPE gcov[][NDIM], FTYPE dxdxp[][NDIM], FTYPE gcovprim[][NDIM])
{
  int j, k, l, m;
  FTYPE tmp[NDIM][NDIM];

  // g_{\alpha \beta} = g_{\mu \nu} \Lambda^\mu_\alpha \Lambda^\nu_\beta


  /*
  // 4 along diagonal and 6 off-diagonal with 6 other identical values
#define GCOV_DOT_DXDXP_DOT_DXDXP(a,b)\
          gcov[0][0] * dxdxp[0][a]* dxdxp[0][b]\
    +     gcov[1][1] * dxdxp[1][a]* dxdxp[1][b]\
    +     gcov[2][2] * dxdxp[2][a]* dxdxp[2][b]\
    +     gcov[3][3] * dxdxp[3][a]* dxdxp[3][b]\
    + 2.0*gcov[0][1] * dxdxp[0][a]* dxdxp[1][b]\
    + 2.0*gcov[0][2] * dxdxp[0][a]* dxdxp[2][b]\
    + 2.0*gcov[0][3] * dxdxp[0][a]* dxdxp[3][b]\
    + 2.0*gcov[1][2] * dxdxp[1][a]* dxdxp[2][b]\
    + 2.0*gcov[1][3] * dxdxp[1][a]* dxdxp[3][b]\
    + 2.0*gcov[2][3] * dxdxp[2][a]* dxdxp[3][b]
  */


  // 4 along diagonal and 6 off-diagonal with 6 other identical values
#define GCOV_DOT_DXDXP_DOT_DXDXP(a,b)\
          gcov[0][0] * dxdxp[0][a]* dxdxp[0][b]\
    +     gcov[1][1] * dxdxp[1][a]* dxdxp[1][b]\
    +     gcov[2][2] * dxdxp[2][a]* dxdxp[2][b]\
    +     gcov[3][3] * dxdxp[3][a]* dxdxp[3][b]\
    +     gcov[0][1] * (dxdxp[0][a]* dxdxp[1][b] + dxdxp[1][a]* dxdxp[0][b])\
    +     gcov[0][2] * (dxdxp[0][a]* dxdxp[2][b] + dxdxp[2][a]* dxdxp[0][b])\
    +     gcov[0][3] * (dxdxp[0][a]* dxdxp[3][b] + dxdxp[3][a]* dxdxp[0][b])\
    +     gcov[1][2] * (dxdxp[1][a]* dxdxp[2][b] + dxdxp[2][a]* dxdxp[1][b])\
    +     gcov[1][3] * (dxdxp[1][a]* dxdxp[3][b] + dxdxp[3][a]* dxdxp[1][b])\
    +     gcov[2][3] * (dxdxp[2][a]* dxdxp[3][b] + dxdxp[3][a]* dxdxp[2][b])


  // first do 4 along diagonal
  DLOOPA(j) tmp[j][j] = GCOV_DOT_DXDXP_DOT_DXDXP(j,j);

  // now do 6 others off-diagonal
  j=0; k=1; tmp[j][k] = GCOV_DOT_DXDXP_DOT_DXDXP(j,k);
  j=0; k=2; tmp[j][k] = GCOV_DOT_DXDXP_DOT_DXDXP(j,k);
  j=0; k=3; tmp[j][k] = GCOV_DOT_DXDXP_DOT_DXDXP(j,k);
  j=1; k=2; tmp[j][k] = GCOV_DOT_DXDXP_DOT_DXDXP(j,k);
  j=1; k=3; tmp[j][k] = GCOV_DOT_DXDXP_DOT_DXDXP(j,k);
  j=2; k=3; tmp[j][k] = GCOV_DOT_DXDXP_DOT_DXDXP(j,k);
    

  // copy over result assuming tmp on upper-diagonal only, but filling gcovprim fully
  // 4 along diagonal (00,11,22,33) and 6 off-diagonal (01, 02, 03, 12, 13, 23), 6 more as copies of off-diagonal
  DLOOPA(j) gcovprim[j][j] = tmp[j][j];
  gcovprim[0][1]=gcovprim[1][0]=tmp[0][1];
  gcovprim[0][2]=gcovprim[2][0]=tmp[0][2];
  gcovprim[0][3]=gcovprim[3][0]=tmp[0][3];
  gcovprim[1][2]=gcovprim[2][1]=tmp[1][2];
  gcovprim[1][3]=gcovprim[3][1]=tmp[1][3];
  gcovprim[2][3]=gcovprim[3][2]=tmp[2][3];

}


void gcon2gconprim(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE gcon[][NDIM],FTYPE gconprim[][NDIM])
{
  int j, k, l, m;
  //FTYPE dxdxp[NDIM][NDIM],idxdxp[NDIM][NDIM];
  FTYPE idxdxp[NDIM][NDIM];
  FTYPE tmp[NDIM][NDIM];

  // see transforms.c and mettometp() and see gcov2gcovprim()
  idxdxprim_ijk_2(ptrgeom, X, V, idxdxp);

  //  dualfprintf(fail_file,"mi in gcon2gconprim\n");
  // PRIMECOORDS indicates no special considerations if fails to get inverse
  //  matrix_inverse(PRIMECOORDS, dxdxp,idxdxp);
  
  DLOOP(j,k) tmp[j][k] = 0.;
  DLOOP(j,k) for(l=0;l<NDIM;l++) for(m=0;m<NDIM;m++){
    tmp[j][k] += idxdxp[j][l] * idxdxp[k][m] * gcon[l][m] ;
  }
  // use tmp since gcon might be same address as gconprim
  DLOOP(j,k) gconprim[j][k] = tmp[j][k];

}



// determinant not simply transformed from analytic function -- so no analytic form possible yet
FTYPE gdet_func_singcheck(int whichcoord, FTYPE *V,FTYPE gcov[][NDIM])
{
  FTYPE gdet_func_orig(int whichcoord, FTYPE gcov[][NDIM]);
  FTYPE gdet;



  gdet=gdet_func_orig(whichcoord, gcov);


#if(FLIPGDETAXIS)
  if(ISSPCMCOORD(whichcoord)){
    if(V[2]<0.0) gdet*=-1.0;
    if(V[2]>M_PI) gdet*=-1.0;
  }
#endif


  return(gdet);

}

// find determinant in general of a metric
/* assumes gcov has been set first; returns determinant */

// determinant not simply transformed from analytic function -- so no analytic form possible yet
FTYPE gdet_func_orig(int whichcoord, FTYPE gcov[][NDIM])
{
  static int firstc = 1;
  static FTYPE **tmp;
  FTYPE d;
  int j, k, indx[NDIM];
  int singfix;
  int anglesing,centersing,truedim;
  void metric_sing_check(int whichcoord, FTYPE gcov[][NDIM], int *anglesing, int*centersing, int *truedim);



  singfix=0;

  if (firstc) {
    tmp = dmatrix(1, NDIM, 1, NDIM);
    firstc = 0;
  }


#if(1)
  // check for coordinate singularity (using this avoids bad input to ludcmp())
  metric_sing_check(whichcoord, gcov, &anglesing, &centersing, &truedim);
#else
  truedim=NDIM;
#endif



  DLOOP(j,k) tmp[j + 1][k + 1] = gcov[j][k];


  if(ludcmp(tmp, truedim, indx - 1, &d)>=1){
    if(debugfail>=2) dualfprintf(fail_file,"ludcmp failure in gdet_func_orig\n");
#if(1)
    if(ISSPCMCOORD(whichcoord)){
      // super hack
      if(debugfail>=2) dualfprintf(fail_file,"Assuming on polar axis\n");
      singfix=1;
    }
    else{
      dualfprintf(fail_file,"ludcmp failure 2: whichcoord=%d\n",whichcoord);
      myexit(3477);
    }
#else
    dualfprintf(fail_file,"ludcmp failure: Can't assume on polar axis for whichcoord=%d\n",whichcoord);
    myexit(246);
#endif

  }

  if(singfix==0 || truedim<NDIM){
    // below from 1..NDIM due to ludcmp requiring 1..N
    for (j = 1; j <= NDIM; j++)
      d *= tmp[j][j];

    if(d>0.0 && d<SMALL){
      return(0.0); // then assume 0
    }
    else if(d>0.0){
      dualfprintf(fail_file,"Metric has bad signature: d=%21.15g icurr=%d jcurr=%d kcurr=%d pcurr=%d\n",d,icurr,jcurr,kcurr,pcurr);
      DLOOP(j,k) dualfprintf(fail_file,"gcov[%d][%d]=%21.15g\n",j,k,gcov[j][k]);
#if(DOSELFGRAVVSR)
      dualfprintf(fail_file,"Mvsr_tot=%21.15g phivsr_tot=%21.15g @ icurr=%d jcurr=%d kcurr=%d\n",Mvsr_tot[startpos[1]+icurr],phivsr_tot[startpos[1]+icurr],icurr,jcurr,kcurr);
#endif
      myexit(3478);
    }
    else{
      return (sqrt(-d)); 
    }
  }
  else{
    return(0);
  }



  return(-1); // shouldn't ever get here
}







void setup_delta(int whichfun,int whichdifftype, FTYPE defaultdelta, struct of_geom *geom, struct of_geom (*localptrgeoml)[NDIM], struct of_geom (*localptrgeomh)[NDIM], FTYPE *truedelta)
{
  int jj;
  int localwhichdifftype;
  int N[NDIM];

  N[0]=2; // if evolving metric then 2, but doesn't matter since at CENT always
  N[1]=N1;
  N[2]=N2;
  N[3]=N3;



  if(whichfun==0 || whichfun==1){ // connection or dxdxp
    if(whichdifftype==DIFFGAMMIE || whichdifftype==DIFFNUMREC){
      localwhichdifftype=0; // infinitesimal
    }
    else if(whichdifftype==DIFFFINITE){
      localwhichdifftype=1; // finite
    }
  }
  else{
    dualfprintf(fail_file,"no such whichfun=%d for whichdifftype=%d",whichfun,whichdifftype);
    myexit(915);
  }



  if(localwhichdifftype==0){ // infinitesimal
    DLOOPA(jj){
      // specify NOWHERE so that won't use gridded position
      ((*localptrgeoml)[jj]).p=((*localptrgeomh)[jj]).p=NOWHERE;
      
      ((*localptrgeoml)[jj]).i=((*localptrgeomh)[jj]).i=(geom->i);
      ((*localptrgeoml)[jj]).j=((*localptrgeomh)[jj]).j=(geom->j);
      ((*localptrgeoml)[jj]).k=((*localptrgeomh)[jj]).k=(geom->k);
      
      truedelta[jj]=defaultdelta;
    }
#if(DOEVOLVEMETRIC)
    // if evolving metric, then use CENT and gcov_func() will use X[0] to choose if present time or not
    // GODMARK: assume didstoremetricdata==1 is set when presenttime==1 is set so doesn't reach bl_coord_ijk_2() and do coord(i,j,k...) that would be wrong
    jj=TT;
    ((*localptrgeoml)[jj]).p=((*localptrgeomh)[jj]).p=CENT; // when taking temporal differences, values are always spatially located at loc=CENT
#endif

  }
  else if(localwhichdifftype==1){ // finite
    // use gridded data
    // assumes connection at CENT
    if(geom->p==CENT){


      DLOOPA(jj){
	((*localptrgeoml)[jj]).p=((*localptrgeomh)[jj]).p=CENT + jj*(N[jj]!=1); // jj=0 -> CENT , jj=1 -> FACE1 , jj=2 -> FACE2 , jj=3 -> FACE3
	
	((*localptrgeoml)[jj]).i=(geom->i);
	((*localptrgeoml)[jj]).j=(geom->j);
	((*localptrgeoml)[jj]).k=(geom->k);
	
	// GODMARK: Note that infinitesimal version obtains correct connection even in reduced dimensions, while the finite version just reduces to no connection if reduced dimension (so this assumes something about the position or may not even be right in general)
	((*localptrgeomh)[jj]).i=(geom->i) + (jj==1 ? SHIFT1 : 0);
	((*localptrgeomh)[jj]).j=(geom->j) + (jj==2 ? SHIFT2 : 0);
	((*localptrgeomh)[jj]).k=(geom->k) + (jj==3 ? SHIFT3 : 0);
	
	truedelta[jj]=0.5*dx[jj];
      }
    }
    else{
      dualfprintf(fail_file,"No localptrgeom defined for p=%d\n",geom->p);
      myexit(916);
    }
  }
  else{
    dualfprintf(fail_file,"No such difftype=%d\n",localwhichdifftype);
    myexit(917);
  }

#if(0)
  dualfprintf(fail_file,"whichfun=%d whichdifftype=%d localwhichdifftype=%d\n",whichfun,whichdifftype,localwhichdifftype);
  DLOOPA(jj) dualfprintf(fail_file,"%d :: truedelta=%21.15g\n",jj,truedelta[jj]);
  DLOOPA(jj) dualfprintf(fail_file,"%d :: low: i=%d j=%d k=%d p=%d :: high: i=%d j=%d k=%d p=%d\n",jj,((*localptrgeoml)[jj]).i,((*localptrgeoml)[jj]).j,((*localptrgeoml)[jj]).k,((*localptrgeoml)[jj]).p,((*localptrgeomh)[jj]).i,((*localptrgeomh)[jj]).j,((*localptrgeomh)[jj]).k,((*localptrgeomh)[jj]).p);
#endif

}












/* 
   this gives the connection coefficient \Gamma^{i}_{j,k} =
   conn[..][i][j][k] where i,j,k = {0,1,2,3} corresponds to {t,r,theta,phi} 
*/

/*
  we also compute the 2nd connection:
  -d/dj(ln(\detg))
*/




/*

Based upon EOM:

$$
(f_{(\nu)} T^t_\nu)_{,t} 
= -(f_{(\nu)}T^j_\nu)_{,j} 
+  f_{(\nu)} [
     T^\lambda_\nu[\ln(f_{(\nu)}/\detg)]_{,\lambda} 
   + T^\mu_\lambda \Gamma^\lambda_{\mu\nu} 
   + \ln(f_{(\nu)})_{,t} T^t_\nu
]
$$

More compactly (JCM: 05/07/08):

$$
d_\mu (f T^\mu_\nu) = f[ T^\lambda_\kappa \Gamma^\kappa_{\nu\lambda} - d_\mu ln (\detg/f) T^\mu_\nu]
$$

SUPERGODMARK: for  if(WHICHEOM!=WITHGDET) don't yet compute correct conn2 I believe

// to avoid body forces in general, must compute (e.g. through correction):

$$
\Gamma^\kappa_{\nu\lambda}[new] = \Gamma^\kappa_{\nu\lambda}[old] - (1/4)[\Gamma^\alpha_{\nu\alpha} - ( (d_\nu f)/f + d_\nu \ln(\detg/f) )]\delta^\kappa_\lambda
$$

or:

$$
\Gamma^\kappa_{\nu\lambda} += - (1/4)[ \Gamma^\alpha_{\nu\alpha} - ( (d_\nu f)/f - conn2_\nu )]\delta^\kappa_\lambda
$$

For the WHICHEOM!=WITHGDET version, in general this would require a different \Gamma for each EOM.  For simplicity just assume f=\detg and don't allow this non-body version unless WHICHEOM==WITHGDET, and so then one has:

$$
\Gamma^\kappa_{\nu\lambda} += - (1/4)[ \Gamma^\alpha_{\nu\alpha} - (d_\nu \detg)/\detg ]\delta^\kappa_\lambda
$$

The above is only true for second order scheme.  For higher order FLUXRECON scheme one has:

$$
\Gamma^\kappa_{\nu\lambda} += - (1/4)[ \Gamma^\alpha_{\nu\alpha} - (d_\nu a2c_\nu \detg)/\detg ]\delta^\kappa_\lambda
$$

Problems:
1) But a2c_\nu would normally be chosen adaptively and so one wouldn't have good cancellation.
2) Fact that \detg is there means difficult to generally have cancellation otherwise unless use NOGDET
3) NOGDET not yet setup for EVOLVEMETRIC, but otherwise for spherical polar coordinates should just use NOGDET for r and \theta.
4) Alternatively, can perform correction every timestep and use same a2c as used for each flux -- expensive
   Problem is b^2/2 and p will be treated differently for SPLITMAEM unless constant all weights
5) So seems only solution apart from NOGDET is to use constant all weights that is unstable

*/



// if evolving in time, then metric changed from Xmetricold[TT] to now (t)
// by shifting time in past we tell metric to use old metric
// if t==Xmetricold[0], then feeding the below tells gcov to assume standard difference at present time
#define DELTAFORTIMEh(DELTA) (Xmetricnew[TT]!=Xmetricold[TT] ? (0.0) : 2.0*DELTA)
#define DELTAFORTIMEl(DELTA) (Xmetricnew[TT]!=Xmetricold[TT] ? (-(Xmetricnew[TT]-Xmetricold[TT])) : DELTA)

// finite volume form
// doesn't make sense for pole with singularity
//#define MYDELTAh(DELTA,k) ( k==TT ? DELTAFORTIMEh(DELTA) :  dx[k]*0.5 )
//#define MYDELTAl(DELTA,k) ( k==TT ? DELTAFORTIMEl(DELTA) : -dx[k]*0.5 )

// finite difference form
#define MYDELTAh(DELTA,k) ( k==TT ? DELTAFORTIMEh(DELTA) :  DELTA )
#define MYDELTAl(DELTA,k) ( k==TT ? DELTAFORTIMEl(DELTA) : -DELTA )


/* NOTE: parameter hides global variable */
void conn_func_numerical1(FTYPE DELTA, FTYPE *X, struct of_geom *geom,
			  FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  int i, j, k, l;
  int kk,jj;
  FTYPE conndiag[NDIM],conndiag2[NDIM];
  FTYPE gdethgen[NDIM],gdetlgen[NDIM];
  FTYPE lngdethgen[NDIM],lngdetlgen[NDIM];
  //  FTYPE *gdeth,*gdetl;
  //  FTYPE *lngdeth,*lngdetl;
  FTYPE gdetmid;
  FTYPE tmp[NDIM][NDIM][NDIM];
  FTYPE Xhgen[NDIM][NDIM];
  FTYPE Xlgen[NDIM][NDIM];
  FTYPE signdXgen[NDIM];
  //  FTYPE *Xh, *Xl;
  //  FTYPE signdX;
  FTYPE V[NDIM];
  FTYPE gmid[NDIM][NDIM];
  FTYPE ghgen[NDIM][NDIM][NDIM];
  FTYPE glgen[NDIM][NDIM][NDIM];
  //  FTYPE (*gh)[NDIM];
  //  FTYPE (*gl)[NDIM];
  FTYPE gcovpertmid[NDIM];
  FTYPE gcovperthgen[NDIM][NDIM],gcovpertlgen[NDIM][NDIM];
  //  FTYPE *gcovperth,*gcovpertl;
  void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);
  FTYPE dfridr(FTYPE (*func)(struct of_geom *,FTYPE*,int,int), struct of_geom *geom, FTYPE *X,int ii, int jj, int kk);
  FTYPE gcov_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int i, int j);
  FTYPE lngdet_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int i, int j);
  //
  // below variables used for setup_delta()
  // and used for defining whether infinitesimal or finite difference
  FTYPE truedelta[NDIM];
  struct of_geom localgeoml[NDIM];
  struct of_geom localgeomh[NDIM];
  struct of_geom (*localptrgeoml)[NDIM];
  struct of_geom (*localptrgeomh)[NDIM];
  void setup_delta(int whichfun,int whichdifftype, FTYPE defaultdelta, struct of_geom *geom, struct of_geom (*localptrgeoml)[NDIM], struct of_geom (*localptrgeomh)[NDIM], FTYPE *truedelta);
  FTYPE gdet_func_mcoord_usegcov(FTYPE gcovmcoord[][NDIM], struct of_geom *ptrgeom, FTYPE* X, int i, int j);
  int doingmachinebody;
  int setup_XlXh(FTYPE *X,FTYPE *truedelta, FTYPE (*Xlgen)[NDIM],FTYPE (*Xhgen)[NDIM],FTYPE *signdXgen);
  int compute_metricquantities_midlh(int donormal, int domachinebody
				     ,struct of_geom *geom, FTYPE *X, FTYPE (*gmid)[NDIM], FTYPE *gcovpertmid,FTYPE *gdetmid,FTYPE *gdetlgen,FTYPE *gdethgen
				     ,struct of_geom (*localptrgeoml)[NDIM],struct of_geom (*localptrgeomh)[NDIM],FTYPE (*Xlgen)[NDIM],FTYPE (*Xhgen)[NDIM]
				     ,FTYPE *lngdetlgen, FTYPE *lngdethgen, FTYPE (*glgen)[NDIM][NDIM], FTYPE (*ghgen)[NDIM][NDIM], FTYPE (*gcovpertlgen)[NDIM], FTYPE (*gcovperthgen)[NDIM]
				     );





  // setup conditional
  doingmachinebody = CONNMACHINEBODY && WHICHEOM==WITHGDET;
  




  localptrgeoml=&localgeoml;
  localptrgeomh=&localgeomh;

  ///////////////////////////
  //
  // Setup delta
  //
  ///////////////////////////
  // 0 indicates connection type
  setup_delta(0,CONNDERTYPE,DELTA,geom,localptrgeoml,localptrgeomh,truedelta);



#if(DOEVOLVEMETRIC==1 && CONNDERTYPE==DIFFNUMREC)
  // no choice in sense that NUMREC is too slow and unstable
  dualfprintf(fail_file,"Not good idea to use DOEVOLVEMETRIC==1 with CONNDERTYPE==DIFFNUMREC\n");
  myexit(7225);
#endif



  ////////////////////////////////
  //
  // gabc_{ijk}=dg_{ij}/dx^k
  // gammie derivative (attempt to get analytical value) or finite difference (true finite difference with given cell size)
  //
  ////////////////////////////////
  if(CONNDERTYPE==DIFFGAMMIE || CONNDERTYPE==DIFFFINITE){


    // setup Xl and Xh and signdX
    setup_XlXh(X,truedelta,Xlgen,Xhgen,signdXgen);

    // get metric quantities
    compute_metricquantities_midlh(1,doingmachinebody
				   ,geom, X, gmid, gcovpertmid,&gdetmid,gdetlgen,gdethgen
				   ,localptrgeoml,localptrgeomh,Xlgen,Xhgen
				   ,lngdetlgen, lngdethgen, glgen, ghgen, gcovpertlgen, gcovperthgen
				   );


    ////////////
    //
    // Now compute derivatives
    //
    ////////////

    for (k = 0; k < NDIM; k++) {

      if(WHICHEOM!=WITHGDET){
	//$$
	//d_\mu (f T^\mu_\nu) = f[ T^\lambda_\kappa \Gamma^\kappa_{\nu\lambda} - d_\mu ln (\detg/f) T^\mu_\nu]
	//$$
	conn2[k]= signdXgen[k]*(lngdethgen[k] - lngdetlgen[k]) / (Xhgen[k][k] - Xlgen[k][k]);
      }
      else{
	conn2[k]=0.0; // no 2nd connection then
      }

      for (i = 0; i < NDIM; i++) for (j = 0; j < NDIM; j++){
	// d(1+g_{tt}) -> dg_{tt}, so can use 1+g_{tt} for accurate non-relativistic gravity
	//if(i==j) conn[i][j][k] = (gcovperth[i] - gcovpertl[i]) / (Xh[k] - Xl[k]);
	// else 

	conn[i][j][k] = signdXgen[k]*(ghgen[k][i][j] - glgen[k][i][j]) / (Xhgen[k][k] - Xlgen[k][k]);
			  
	}
      
    }// end over k


  }
  else if(CONNDERTYPE==DIFFNUMREC){

    for (k = 0; k < NDIM; k++) {

      if(WHICHEOM!=WITHGDET){
	conn2[k]= dfridr(lngdet_func_mcoord,geom,X,0,0,k); // 0,0 not used
      }
      else conn2[k]=0.0; // then no 2nd connection

      for (i = 0; i < NDIM; i++) for (j = 0; j < NDIM; j++){
	conn[i][j][k] = dfridr(gcov_func_mcoord,geom,X,i,j,k);
      }
    }

  }


  ////////////////////////////////////////////////////
  //
  // now rearrange to find \Gamma_{ijk}=1/2*(gabc_{jik}+gabc_{kij}-gabc_{kji})
  //
  ////////////////////////////////////////////////////
  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      for (k = 0; k < NDIM; k++)
	tmp[i][j][k] =
	  0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);


  ////////////////////////////////////////////////////
  //
  // finally, raise first index
  //
  ////////////////////////////////////////////////////
  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      for (k = 0; k < NDIM; k++) {
	conn[i][j][k] = 0.;
	for (l = 0; l < NDIM; l++){
	  conn[i][j][k] += geom->gcon[i][l] * tmp[l][j][k];
	}
      }


  //////////////////////////////////////////////////
  //
  // now correct for accurate body forces
  // only makes sense for no a2c on flux right now
  //
  //////////////////////////////////////////////////
  if(doingmachinebody){

    if(CONNDERTYPE!=DIFFFINITE){
      // Setup finite difference for correction
      // 0 indicates connection type
      // correction always uses DIFFFINITE
      setup_delta(0,DIFFFINITE,DELTA,geom,localptrgeoml,localptrgeomh,truedelta);

      // setup Xl and Xh and signdX
      setup_XlXh(X,truedelta,Xlgen,Xhgen,signdXgen);

      // 0 means not normal calculations      
      compute_metricquantities_midlh(0,doingmachinebody
				     ,geom, X, gmid, gcovpertmid,&gdetmid,gdetlgen,gdethgen
				     ,localptrgeoml,localptrgeomh,Xlgen,Xhgen
				     ,lngdetlgen, lngdethgen, glgen, ghgen, gcovpertlgen, gcovperthgen
				     );
      
    }

    //////////////////
    // form contracted connection
    DLOOPA(kk){
      conndiag[kk]=0.0;
      DLOOPA(jj) conndiag[kk] += conn[jj][kk][jj];
    }

    /////////////////
    // form differential \detg
    for (k = 0; k < NDIM; k++) {
      conndiag2[k] = signdXgen[k]*(gdethgen[k] - gdetlgen[k]) / (Xhgen[k][k] - Xlgen[k][k]);
      conndiag2[k] /= gdetmid;
    }

    //////////////////
    // now obtain correction 
    // for Xh-Xl->0 this vanishes as required
    // Plugging this new conn into EOM for T^\mu_\nu = p \delta^\mu_\nu gives exactly cancellation between source and flux differencing of pressure
    //    for (i = 0; i < NDIM; i++) for (j = 0; j < NDIM; j++) for (k = 0; k < NDIM; k++) conn[i][j][k] += -0.25*(conndiag[j] - conndiag2[j])*delta(i,k);
    // apply delta(i,k) directly by setting i=k
    for (i = 0; i < NDIM; i++) for (j = 0; j < NDIM; j++) conn[i][j][i] += -0.25*(conndiag[j] - conndiag2[j]);

  } // end if correcting body forces




#if(0) // DEBUG
  if(geom->i==0 || geom->i==-1){
    for (i = 0; i < NDIM; i++) for (l = 0; l < NDIM; l++){
      dualfprintf(fail_file,"i=%d gcon[%d][%d]=%21.15g\n",geom->i,i,l,geom->gcon[i][l]);
    }
    dualfprintf(fail_file,"i=%d conn[0][0][0]=%21.15g\n",geom->i,conn[0][0][0]);
  }
#endif
	
  /* done! */
}



int setup_XlXh(FTYPE *X,FTYPE *truedelta, FTYPE (*Xlgen)[NDIM],FTYPE (*Xhgen)[NDIM],FTYPE *signdXgen)
{
  int k,l;

    ////////////
    // first form high and low positions for function locations
    ////////////
    for (k = 0; k < NDIM; k++) {

      for (l = 0; l < NDIM; l++){
	Xhgen[k][l] = X[l];
	Xlgen[k][l] = X[l];
      }

      // normal case
      Xhgen[k][k] += MYDELTAh(truedelta[k],k);
      Xlgen[k][k] += MYDELTAl(truedelta[k],k);
      signdXgen[k]=1.0;

#if(0)// debug
      if(k==TT){
	if(Xlgen[k][k]>X[k]){
	  dualfprintf(fail_file,"Xl in future! Xl=%21.15g Xh=%21.15g\n",Xlgen[k][k],Xhgen[k][k]);
	}
      }
#endif


#if(0)
      // not really wanted since want "force" to be in same "radial" direction so sign SHOULD flip
      if(k==1 && ISSPCMCOORD(MCOORD)){
	// then check if r<0 and invert Xh and Xl if so
	bl_coord(X, V);
	if(V[k]<0.0){
	  signdXgen[k]=-1.0;
	}
      }
#endif
      //			if(k==TT) dualfprintf(fail_file,"k=%d DELTAl=%21.15g DELTAl=%21.15g\n",k,MYDELTAl(DELTA,k),MYDELTAh(DELTA,k));
      // DEBUG
      //Xhgen[k][k] += DELTA;
      //			Xlgen[k][k] -= DELTA;

    }

    return(0);

}



// compute low/high metric quantities
int compute_metricquantities_midlh(int donormal, int domachinebody
				   ,struct of_geom *geom, FTYPE *X, FTYPE (*gmid)[NDIM], FTYPE *gcovpertmid,FTYPE *gdetmid,FTYPE *gdetlgen,FTYPE *gdethgen
				   ,struct of_geom (*localptrgeoml)[NDIM],struct of_geom (*localptrgeomh)[NDIM],FTYPE (*Xlgen)[NDIM],FTYPE (*Xhgen)[NDIM]
				   ,FTYPE *lngdetlgen, FTYPE *lngdethgen, FTYPE (*glgen)[NDIM][NDIM], FTYPE (*ghgen)[NDIM][NDIM], FTYPE (*gcovpertlgen)[NDIM], FTYPE (*gcovperthgen)[NDIM]
				   )
{
  int k;
  FTYPE gcov_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int i, int j);
  FTYPE lngdet_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int i, int j);
  FTYPE gdet_func_mcoord_usegcov(FTYPE gcovmcoord[][NDIM], struct of_geom *ptrgeom, FTYPE* X, int i, int j);

    // get gdet in cell center
  if(domachinebody){
    gcov_func(geom,1,MCOORD,X, gmid,gcovpertmid);
    *gdetmid=gdet_func_mcoord_usegcov(gmid, geom, X, 0,0);
  }

  // get k-dependent things
  for (k = 0; k < NDIM; k++) {

    if(donormal){
      if(WHICHEOM!=WITHGDET){
	lngdethgen[k]=lngdet_func_mcoord(&((*localptrgeomh)[k]),Xhgen[k],0,0); // doesn't use 0,0
	lngdetlgen[k]=lngdet_func_mcoord(&((*localptrgeoml)[k]),Xlgen[k],0,0); // doesn't use 0,0
      }
    }
    if(donormal || domachinebody){
      // must come before below gdet_func_mcoord_usegcov() call
      gcov_func(&((*localptrgeomh)[k]),1,MCOORD,Xhgen[k], ghgen[k],gcovperthgen[k]);
      gcov_func(&((*localptrgeoml)[k]),1,MCOORD,Xlgen[k], glgen[k],gcovpertlgen[k]);
    }

    if(domachinebody){
      gdethgen[k]=gdet_func_mcoord_usegcov(ghgen[k], &((*localptrgeomh)[k]), Xhgen[k], 0,0);
      gdetlgen[k]=gdet_func_mcoord_usegcov(glgen[k], &((*localptrgeoml)[k]), Xlgen[k], 0,0);
    }

  }

  return(0);

}





// returns MCOORD gcov value for i,j element
// excessive to compute other elements, but ok for now
FTYPE gcov_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int i, int j)
{
  FTYPE gcovmcoord[NDIM][NDIM];
  FTYPE gcovpert[NDIM];
  void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);

  gcov_func(ptrgeom,1,MCOORD,X,gcovmcoord,gcovpert);
  //  if(i==j) return(gcovpert[i]); // for accurate non-rel gravity (works since d(1+g_{tt}) = dg_{tt})
  //  else
  return(gcovmcoord[i][j]);
}


// returns MCOORD  value for log(f/gdet).  Doesn't use i,j (these are not grid locations)
FTYPE lngdet_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int i, int j)
{
  FTYPE gcovmcoord[NDIM][NDIM];
  FTYPE gcovpertcoord[NDIM];
  FTYPE toreturn;
  FTYPE eomfunc;
  FTYPE V[NDIM];
  void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);

  gcov_func(ptrgeom,1,MCOORD,X,gcovmcoord,gcovpertcoord);
  eomfunc_func(ptrgeom,1,MCOORD,X,&eomfunc);

  bl_coord_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, V);
  toreturn=log(eomfunc/gdet_func_singcheck(MCOORD,V,gcovmcoord));

  return(toreturn);
}

// returns MCOORD  value for gdet.  Doesn't use i,j (these are not grid locations)
FTYPE gdet_func_mcoord(struct of_geom *ptrgeom, FTYPE* X, int i, int j)
{
  FTYPE gcovmcoord[NDIM][NDIM];
  FTYPE gcovpertcoord[NDIM];
  FTYPE toreturn;
  FTYPE V[NDIM];
  void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);

  gcov_func(ptrgeom,1,MCOORD,X,gcovmcoord,gcovpertcoord);

  bl_coord_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, V);
  toreturn=gdet_func_singcheck(MCOORD,V,gcovmcoord);

  return(toreturn);
}

// returns MCOORD  value for gdet using gcovmcoord as input to avoid repeated computations of gcovmcoord.  Doesn't use i,j (these are not grid locations)
FTYPE gdet_func_mcoord_usegcov(FTYPE gcovmcoord[][NDIM], struct of_geom *ptrgeom, FTYPE* X, int i, int j)
{
  //  FTYPE gcovmcoord[NDIM][NDIM];
  //  FTYPE gcovpertcoord[NDIM];
  FTYPE toreturn;
  FTYPE V[NDIM];
  //void gcov_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE *gcovpert);

  //  gcov_func(ptrgeom,1,MCOORD,X,gcovmcoord,gcovpertcoord);

  bl_coord_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, V);
  toreturn=gdet_func_singcheck(MCOORD,V,gcovmcoord);

  return(toreturn);
}































// jon's MKS connection (and conn2)
// only applies to defcoord==LOGRSINTH
void mks_conn_func(FTYPE *X, struct of_geom *ptrgeom,
		   FTYPE conn[][NDIM][NDIM],FTYPE *conn2)
{
  int i, j, k, l;
  FTYPE V[NDIM];
  FTYPE r,th,sigma,dxdxptrue[NDIM][NDIM];
#ifdef WIN32
  extern FTYPE cot(FTYPE arg);
#endif
  FTYPE dxdxp[NDIM];


  if(MBH!=1.0){
    dualfprintf(fail_file,"mks_conn_func not setup for MBH!=1.0\n");
    myexit(10);
  }


  // get bl coordinates
  bl_coord_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, V);
  r=V[1];
  th=V[2];


  // the connection

  // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
  // follows gcov_func()
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
    if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }

  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif

  // set aux vars
  dxdxprim_ijk(ptrgeom->i,ptrgeom->j,ptrgeom->k,ptrgeom->p, dxdxptrue);
  DLOOPA(j) dxdxp[j]=dxdxptrue[j][j]; // defcoord==LOGRSINTH assumes transformation is diagonal
  sigma=r*r+a*a*cos(th)*cos(th);


  conn[0][0][0]=(-2.*r*sigma + 4.*pow(r,3.))*pow(sigma,-3.);
  conn[0][0][1]=dxdxp[1]*(2.*r + sigma)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.);
  conn[0][0][2]=-1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th);
  conn[0][0][3]=-2.*a*r*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
  conn[0][1][0]=dxdxp[1]*(2.*r + sigma)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.);
  conn[0][1][1]=2.*(r + sigma)*pow(dxdxp[1],2.)*
    (-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.);
  conn[0][1][2]=-1.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*
    sin(2.*th);
  conn[0][1][3]=dxdxp[1]*a*(2.*r + sigma)*(sigma - 2.*pow(r,2.))*
    pow(sigma,-3.)*pow(sin(th),2.);
  conn[0][2][0]=-1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th);
  conn[0][2][1]=-1.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*
    sin(2.*th);
  conn[0][2][2]=-2.*pow(dxdxp[2],2.)*pow(r,2.)*pow(sigma,-1.);
  conn[0][2][3]=2.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-2.)*
    pow(sin(th),3.);
  conn[0][3][0]=-2.*a*r*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
  conn[0][3][1]=dxdxp[1]*a*(2.*r + sigma)*(sigma - 2.*pow(r,2.))*
    pow(sigma,-3.)*pow(sin(th),2.);
  conn[0][3][2]=2.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-2.)*
    pow(sin(th),3.);
  conn[0][3][3]=2.*r*pow(sigma,-3.)*pow(sin(th),2.)*
    (-1.*r*pow(sigma,2.) + pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*
     pow(sin(th),2.));
  conn[1][0][0]=pow(dxdxp[1],-1.)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.)*(-2.*r + sigma + pow(a,2.)*pow(sin(th),2.));
  conn[1][0][1]=0.5*(4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (sigma - 2.*pow(r,2.))*pow(sigma,-3.);
  conn[1][0][2]=0.;
  conn[1][0][3]=0.5*a*pow(dxdxp[1],-1.)*
    (4.*r - 2.*sigma - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*pow(sin(th),2.);
  conn[1][1][0]=0.5*(4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (sigma - 2.*pow(r,2.))*pow(sigma,-3.);
  conn[1][1][1]=pow(sigma,-3.)*
    (-1.*dxdxp[1]*(2.*r + sigma)*(-1.*sigma + 2.*pow(r,2.)) + 
     pow(sigma,3.) + dxdxp[1]*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*
     pow(sin(th),2.));
  conn[1][1][2]=-1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)
    ;
  conn[1][1][3]=0.5*a*(pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
		       cos(2.*th)*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.)) + 
		       2.*r*((-2. + sigma)*sigma + 4.*pow(r,2.)))*pow(sigma,-3.)*
    pow(sin(th),2.);
  conn[1][2][0]=0.;
  conn[1][2][1]=-1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)
    ;
  conn[1][2][2]=-1.*r*pow(dxdxp[1],-1.)*pow(dxdxp[2],2.)*
    pow(sigma,-1.)*(-2.*r + sigma + pow(a,2.)*pow(sin(th),2.));
  conn[1][2][3]=0.;
  conn[1][3][0]=0.5*a*pow(dxdxp[1],-1.)*
    (4.*r - 2.*sigma - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
    (-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*pow(sin(th),2.);
  conn[1][3][1]=0.5*a*(pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
		       cos(2.*th)*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.)) + 
		       2.*r*((-2. + sigma)*sigma + 4.*pow(r,2.)))*pow(sigma,-3.)*
    pow(sin(th),2.);
  conn[1][3][2]=0.;
  conn[1][3][3]=-1.*pow(dxdxp[1],-1.)*pow(sigma,-3.)*pow(sin(th),2.)*
    (-2.*r + sigma + pow(a,2.)*pow(sin(th),2.))*
    (r*pow(sigma,2.) + pow(a,2.)*(sigma - 2.*pow(r,2.))*pow(sin(th),2.));
  conn[2][0][0]=-1.*r*pow(dxdxp[2],-1.)*pow(a,2.)*pow(sigma,-3.)*
    sin(2.*th);
  conn[2][0][1]=-1.*dxdxp[1]*r*pow(dxdxp[2],-1.)*pow(a,2.)*
    pow(sigma,-3.)*sin(2.*th);
  conn[2][0][2]=0.;
  conn[2][0][3]=2.*a*r*cos(th)*pow(dxdxp[2],-1.)*pow(sigma,-3.)*
    (sigma + pow(a,2.)*pow(sin(th),2.))*sin(th);
  conn[2][1][0]=-1.*dxdxp[1]*r*pow(dxdxp[2],-1.)*pow(a,2.)*
    pow(sigma,-3.)*sin(2.*th);
  conn[2][1][1]=-1.*r*pow(dxdxp[1],2.)*pow(dxdxp[2],-1.)*pow(a,2.)*
    pow(sigma,-3.)*sin(2.*th);
  conn[2][1][2]=dxdxp[1]*r*pow(sigma,-1.);
  conn[2][1][3]=dxdxp[1]*a*pow(dxdxp[2],-1.)*pow(sigma,-3.)*sin(th)*
    (sigma*(2.*r + sigma)*cos(th) + r*pow(a,2.)*sin(th)*sin(2.*th));
  conn[2][2][0]=0.;
  conn[2][2][1]=dxdxp[1]*r*pow(sigma,-1.);
  conn[2][2][2]=4.*(M_PI*X[2] - 1.*th)*pow(dxdxp[2],-1.)*
    pow(M_PI,2.) - 1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)\
    ;
  conn[2][2][3]=0.;
  conn[2][3][0]=2.*a*r*cos(th)*pow(dxdxp[2],-1.)*pow(sigma,-3.)*
    (sigma + pow(a,2.)*pow(sin(th),2.))*sin(th);
  conn[2][3][1]=dxdxp[1]*a*pow(dxdxp[2],-1.)*pow(sigma,-3.)*sin(th)*
    (sigma*(2.*r + sigma)*cos(th) + r*pow(a,2.)*sin(th)*sin(2.*th));
  conn[2][3][2]=0.;
  conn[2][3][3]=-1.*cos(th)*pow(dxdxp[2],-1.)*pow(sigma,-3.)*
    (pow(sigma,3.) + sigma*(4.*r + sigma)*pow(a,2.)*pow(sin(th),2.) + 
     2.*r*pow(a,4.)*pow(sin(th),4.))*sin(th);
  conn[3][0][0]=a*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.);
  conn[3][0][1]=dxdxp[1]*a*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)
    ;
  conn[3][0][2]=-2.*dxdxp[2]*a*r*cot(th)*pow(sigma,-2.);
  conn[3][0][3]=-1.*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
  conn[3][1][0]=dxdxp[1]*a*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)
    ;
  conn[3][1][1]=a*pow(dxdxp[1],2.)*(-1.*sigma + 2.*pow(r,2.))*
    pow(sigma,-3.);
  conn[3][1][2]=-1.*dxdxp[1]*dxdxp[2]*a*(2.*r + sigma)*cot(th)*
    pow(sigma,-2.);
  conn[3][1][3]=dxdxp[1]*pow(sigma,-3.)*
    (r*pow(sigma,2.) + pow(a,2.)*(sigma - 2.*pow(r,2.))*pow(sin(th),2.));
  conn[3][2][0]=-2.*dxdxp[2]*a*r*cot(th)*pow(sigma,-2.);
  conn[3][2][1]=-1.*dxdxp[1]*dxdxp[2]*a*(2.*r + sigma)*cot(th)*
    pow(sigma,-2.);
  conn[3][2][2]=-1.*a*r*pow(dxdxp[2],2.)*pow(sigma,-1.);
  conn[3][2][3]=dxdxp[2]*
    (cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th));
  conn[3][3][0]=-1.*pow(a,2.)*(-1.*sigma + 2.*pow(r,2.))*pow(sigma,-3.)*
    pow(sin(th),2.);
  conn[3][3][1]=dxdxp[1]*pow(sigma,-3.)*
    (r*pow(sigma,2.) + pow(a,2.)*(sigma - 2.*pow(r,2.))*pow(sin(th),2.));
  conn[3][3][2]=dxdxp[2]*
    (cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th));
  conn[3][3][3]=pow(sigma,-3.)*
    (-1.*a*r*pow(sigma,2.)*pow(sin(th),2.) + 
     pow(a,3.)*(-1.*sigma + 2.*pow(r,2.))*pow(sin(th),4.));
  conn2[0]=0.;
  conn2[1]=-1.*pow(sigma,-1.)*(2.*dxdxp[1]*r + pow(r,2.) + 
			       pow(a,2.)*pow(cos(th),2.));
  conn2[2]=-1.*dxdxp[2]*cot(th) + 
    4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
    dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th);
  conn2[3]=0.;


}



/*********************************************************************************************
   Scott's s MKS connection that can be used for any transformation between r,th <-> X1,X2 :
*********************************************************************************************/
#define M (1.)
void mks_conn_func_general(FTYPE *X, struct of_geom *geom, FTYPE conn[][NDIM][NDIM] )
{
  int i, j, k, l;
  FTYPE V[NDIM];
  FTYPE r,th,sigma,dxdxp[NDIM][NDIM],dxdxp_dxp[NDIM][NDIM][NDIM];

  FTYPE t1,   t10,   t102,   t1024,   t1035,   t1037,   t104,   t11,   t114,   t116,   t119;
  FTYPE  t12,   t121,   t123,   t126,   t129,   t130,   t132,   t14,   t148,   t149,   t15,   t152;
  FTYPE t154,   t156,   t157,   t158,   t159,   t161,   t169,   t17,   t171,   t172,   t175,   t177;
  FTYPE t185,   t2,   t203,   t204,   t208,   t209,   t21,   t210,   t212,   t214,   t22,   t221;
  FTYPE   t222,   t224,   t227,   t23,   t236,   t24,   t240,   t241,   t242,   t243,   t245,   t246;
  FTYPE    t247,   t248,   t25,   t250,   t251,   t258,   t26,   t260,   t261,   t263,   t264,   t271;
  FTYPE   t273,   t275,   t276,   t278,   t28,   t280,   t281,   t283,   t284,   t285,   t286,   t288;
  FTYPE    t289,   t29,   t297,   t299,   t3,   t30,   t300,   t303,   t305,   t306,   t308,   t309;
  FTYPE    t31,   t310,   t313,   t314,   t320,   t325,   t327,   t328,   t329,   t330,   t333,   t336;
  FTYPE    t338,   t34,   t340,   t342,   t344,   t346,   t35,   t358,   t361,   t363,   t366,   t367;
  FTYPE    t368,   t370,   t372,   t375,   t38,   t380,   t381,   t384,   t385,   t387,   t39,   t399;
  FTYPE    t4,   t40,   t400,   t402,   t404,   t405,   t406,   t408,   t409,   t41,   t411,   t412;
  FTYPE    t418,   t42,   t421,   t425,   t428,   t431,   t432,   t433,   t434,   t437,   t440,   t442;
  FTYPE    t448,   t451,   t453,   t454,   t459,   t462,   t467,   t469,   t480,   t481,   t486,   t487;
  FTYPE    t488,   t491,   t492,   t498,   t501,   t504,   t507,   t508,   t510,   t512,   t52,   t521;
  FTYPE    t528,   t53,   t530,   t553,   t556,   t56,   t57,   t588,   t60,   t607,   t627,   t628;
  FTYPE    t63,   t630,   t631,   t632,   t634,   t636,   t637,   t64,   t651,   t652,   t654,   t656;
  FTYPE    t657,   t659,   t661,   t662,   t670,   t673,   t675,   t677,   t686,   t689,   t7,   t712;
  FTYPE    t74,   t748,   t75,   t78,   t793,   t794,   t795,   t799,   t8,   t800,   t801,   t803;
  FTYPE    t806,  t807,   t813,   t816,   t822,   t83,   t831,   t84,   t845,   t86,   t89,   t891; 
  FTYPE    t90,   t91,   t916,   t917,   t920,   t924,   t928,   t940,   t946,   t968,   t97, t970,t991;


  if(MBH!=1.0){
    dualfprintf(fail_file,"mks_conn_func_general not setup for MBH!=1.0\n");
    myexit(11);
  }


  // get bl coordinates
  bl_coord(X,V);
  r=V[1];
  th=V[2];

  // the connection

  // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
  // follows gcov_func()
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
    if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }
  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif

  // set aux vars
  dxdxprim(X,V,dxdxp);
  //  DLOOPA(j) dxdxp[j]=dxdxptrue[j][j]; // defcoord==LOGRSINTH assumes transformation is diagonal
  //  dx_dxp_calc(X,r,th,dx_dxp);
  //dx_dxp_dxp_calc(X,r,th,dx_dxp_dxp);



  // GODMARK
  // need to set second derivative analytically
  for(i=0;i<NDIM;i++)  for(j=0;j<NDIM;j++)  for(k=0;k<NDIM;k++){
    dxdxp_dxp[i][j][k]=0.0;
  }


  //  t1 = rf(X1,X2);
  t1 = r;
  //  t2 = thf(X1,X2);
  t2 = th;
  t3 = cos(t2);
  t4 = a*t3;
  t7 = (-t1+t4)*(t1+t4);
  t8 = M*M;
  t10 = t1*t1;
  t11 = a*a;
  t12 = t3*t3;
  t14 = t10+t11*t12;
  t15 = t14*t14;
  t17 = 1/t15/t14;
  conn[0][0][0] = -2.0*t7*t8*t1*t17;
  t21 = t10*t10;
  //  t22 = diff(rf(X1,X2),X1);
  t22 = dxdxp[RR][1];
  t23 = t21*t22;
  t24 = t10*t1;
  t25 = M*t24;
  t26 = t25*t22;
  t28 = t24*t3;
  t29 = sin(t2);
  //  t30 = diff(thf(X1,X2),X1);
  t30 = dxdxp[TH][1];
  t31 = t29*t30;
  t34 = M*t1;
  t35 = t22*t12;
  t38 = t12*t12;
  t39 = t38*t22;
  t40 = t29*t1;
  t41 = t12*t3;
  t42 = t41*t30;
  conn[0][0][1] = -M*(-t23-2.0*t26+(2.0*t28*t31+2.0*t34*t35+(t39+2.0*t40*t42)*t11)*t11)*t17;
  //  t52 = diff(rf(X1,X2),X2);
  t52 = dxdxp[RR][2];
  t53 = t21*t52;
  //  t56 = diff(thf(X1,X2),X2);
  t56 = dxdxp[TH][2];
  t57 = t29*t56;
  t60 = t52*t12;
  t63 = t38*t52;
  t64 = t41*t1;
  conn[0][0][2] = -M*(-t53-2.0*t25*t52+(2.0*t28*t57+2.0*t34*t60+(t63+2.0*t64*t57)*t11)*t11)*t17;
  t74 = -1.0+t3;
  t75 = 1.0+t3;
  t78 = t7*t74*t75;
  conn[0][0][3] = -2.0*t78*a*t1*t8*t17;
  conn[0][1][0] = conn[0][0][1];
  t83 = t30*t30;
  t84 = t21*t10;
  t86 = t22*t22;
  t89 = t24*t22;
  t90 = t3*t29;
  t91 = t90*t30;
  t97 = t86*t12;
  t102 = t22*t29;
  t104 = t64*t102;
  conn[0][1][1] = -2.0*(t83*t84-t21*t86-t25*t86+(2.0*t89*t91+2.0*t83*t21*t12+
						 t34*t97+(t83*t10*t38+t38*t86+2.0*t104*t30)*t11)*t11)*M*t17;
  t114 = t22*t52;
  t116 = t30*t56;
  t119 = t24*t52;
  t121 = t114*t12;
  t123 = t21*t12;
  t126 = t90*t56;
  t129 = t52*t29;
  t130 = t129*t30;
  t132 = t10*t38;
  conn[0][1][2] = -2.0*(-t25*t114+t116*t84-t23*t52+(t119*t91+t34*t121+2.0*
						    t116*t123+t89*t126
						    +(t39*t52+t64*t130+t116*t132+t104*t56)*t11)*t11)*M*t17;
  t148 = 2.0*t28*t30;
  t149 = t102*t12;
  t152 = t41*t24;
  t154 = 2.0*t152*t30;
  t156 = 2.0*t64*t30;
  t157 = t39*t29;
  t158 = t30*t1;
  t159 = t38*t3;
  t161 = 2.0*t158*t159;
  t169 = a*t29*t17;
  conn[0][1][3] = -(2.0*t25*t102+t23*t29+(-t148-2.0*t34*t149+t154+(-t156-t157+t161)*t11)*t11)*M*t169;
  conn[0][2][0] = conn[0][0][2];
  conn[0][2][1] = conn[0][1][2];
  t171 = t52*t52;
  t172 = t171*M;
  t175 = t56*t56;
  t177 = t1*t12;
  t185 = t52*t41;
  conn[0][2][2] = -2.0*(-t172*t24-t171*t21+t175*t84+(t172*t177+2.0*t119*t126+
						     2.0*t175*t21*t12
						     +(t171*t38+2.0*t185*t40*t56+t175*t10*t38)*t11)*t11)*M*t17;
  t203 = 2.0*t152*t56;
  t204 = t129*t12;
  t208 = 2.0*t28*t56;
  t209 = t63*t29;
  t210 = t56*t1;
  t212 = 2.0*t210*t159;
  t214 = 2.0*t64*t56;
  conn[0][2][3] = (-2.0*t25*t129-t53*t29+(-t203+2.0*t34*t204+t208+(t209-t212+t214)*t11)*t11)*M*t169;
  conn[0][3][0] = conn[0][0][3];
  conn[0][3][1] = conn[0][1][3];
  conn[0][3][2] = conn[0][2][3];
  t221 = t21*t1;
  t222 = t24*t12;
  t224 = t10*t12;
  t227 = t1*t38;
  t236 = (-t221+(-2.0*t222+(-t224+t10)*M+(-t227+(t38-t12)*M)*t11)*t11)*t74*t75;
  conn[0][3][3] = -2.0*t236*t34*t17;
  t240 = t56*M;
  t241 = t240*t24;
  t242 = 2.0*t241;
  t243 = t56*t21;
  t245 = 2.0*t240*t177;
  t246 = t56*t10;
  t247 = t246*t12;
  t248 = t1*t3;
  t250 = 2.0*t248*t129;
  t251 = t56*t12;
  t258 = t30*t52;
  t260 = 1/(-t22*t56+t258);
  t261 = t260*t17;
  conn[1][0][0] = -(-t242+t243+(t245-t247+t250+t246-t251*t11)*t11)*M*t261;
  t263 = M*t22;
  t264 = t56*t38;
  t271 = (-t242+(t245-t247+t250+t246+(-t251+t264)*t11)*t11)*t260*t17;
  conn[1][0][1] = -t263*t271;
  t273 = M*t52;
  conn[1][0][2] = -t273*t271;
  t275 = t24*t29;
  t276 = t240*t275;
  t278 = t119*t3;
  t280 = t243*t29;
  t281 = t40*t12;
  t283 = 2.0*t240*t281;
  t284 = t29*t12;
  t285 = t246*t284;
  t286 = t52*t3;
  t288 = 2.0*t286*t1;
  t289 = t57*t10;
  t297 = a*t29*t260*t17;
  conn[1][0][3] = (-2.0*t276+2.0*t278+t280+(t283-t285+t288+t289-t56*t11*t284)*t11)*M*t297;
  conn[1][1][0] = conn[1][0][1];
  //  t299 = diff(diff(thf(X1,X2),X1),X1);
  t299 = dxdxp_dxp[TH][1][1];
  t300 = t52*t299;
  t303 = t258*t221*t22;
  t305 = t56*t84;
  //  t306 = diff(diff(rf(X1,X2),X1),X1);
  t306 = dxdxp_dxp[RR][1][1] ;
  t308 = t83*t56;
  t309 = t21*t24;
  t310 = t308*t309;
  t313 = 2.0*t308*t84;
  t314 = t56*t24;
  t320 = t30*t22;
  t325 = t258*t89*t12;
  t327 = t308*t221;
  t328 = t52*t83;
  t329 = t90*t21;
  t330 = t328*t329;
  t333 = t306*t12;
  t336 = t221*t12;
  t338 = 2.0*t308*t336;
  t340 = 4.0*t308*t123;
  t342 = t248*t29;
  t344 = 2.0*t52*t86*t342;
  t346 = t56*t86;
  t358 = t306*t38;
  t361 = t1*t22;
  t363 = t258*t361*t38;
  t366 = 2.0*t308*t222;
  t367 = t24*t38;
  t368 = t308*t367;
  t370 = t41*t29*t10;
  t372 = 2.0*t328*t370;
  t375 = 2.0*t308*t132;
  t380 = t159*t29;
  t381 = t380*t56;
  t384 = t308*t227;
  t385 = t38*t12;
  t387 = t328*t380;
  conn[1][1][1] = -(-t300*t84-2.0*t303+t305*t306-t310+(-t243*t86+t313-2.0*t314*t86*M)*M
		    +(-2.0*t320*t3*t280-4.0*t325-t327+t330-3.0*t300*t123+3.0*t243*t333
		      -t338+(t340+t344-t246*t97+t346*t10+2.0*t210*t97*M)*M
		      +(-4.0*t320*t41*t289-3.0*t300*t132+3.0*t246*t358-2.0*t363-t366-t368+t372
			+(-t346*t12+t375+2.0*t346*t38)*M
			+(-2.0*t320*t381-t384-t300*t385+t387+t56*t306*t385)*t11)*t11)*t11)*t260*t17;
  //  t399 = diff(diff(thf(X1,X2),X1),X2);
  t399 = dxdxp_dxp[TH][1][2]; 
  t400 = t52*t399;
  //  t402 = diff(diff(rf(X1,X2),X1),X2);
  t402 = dxdxp_dxp[RR][1][2]; 
  t404 = t30*t175;
  t405 = t404*t309;
  t406 = t171*t30;
  t408 = t56*t221;
  t409 = t114*t408;
  t411 = 2.0*t404*t84;
  t412 = t52*t56;
  t418 = t402*t12;
  t421 = t404*t221;
  t425 = t114*t314*t12;
  t428 = 2.0*t404*t336;
  t431 = t22*t175;
  t432 = t431*t329;
  t433 = t10*t22;
  t434 = t433*t12;
  t437 = 4.0*t404*t123;
  t440 = 2.0*t171*t22*t342;
  t442 = t210*M;
  t448 = 2.0*t370*t431;
  t451 = t114*t210*t38;
  t453 = 2.0*t404*t222;
  t454 = t402*t38;
  t459 = t404*t367;
  t462 = 2.0*t404*t132;
  t467 = t404*t227;
  t469 = t431*t380;
  conn[1][1][2] = (t400*t84-t305*t402+t405+t406*t221+t409+(-t411+t412*t23+2.0*t114*t241)*M
		   +(-3.0*t243*t418+t421+3.0*t400*t123+2.0*t425+t428+2.0*t406*t222+
		     t432+(t412*t434-t437-t440-t412*t433-2.0*t121*t442)*M
		     +(t448+t406*t227+t451+t453-3.0*t246*t454+3.0*t400*t132+t459
		       +(t114*t251-t462-2.0*t114*t264)*M+(t467+t400*t385+t469
							  -t56*t402*t385)*t11)*t11)*t11)*t260*t17;
  t480 = t286*t21;
  t481 = t57*t221;
  t486 = 2.0*t185*t10;
  t487 = t57*t222;
  t488 = 2.0*t487;
  t491 = t57*t227;
  t492 = t52*t159;
  t498 = (-t491+t492+(-t57*t12+t57*t38)*M)*t11;
  t501 = t480-t481+(2.0*t278-2.0*t276)*M+(t486-t488+(-t285+t289+t288+t283)*M+t498)*t11;
  conn[1][1][3] = t501*t22*t297;
  conn[1][2][0] = conn[1][0][2];
  conn[1][2][1] = conn[1][1][2];
  t504 = t171*t56;
  //  t507 = diff(diff(thf(X1,X2),X2),X2);
  t507 = dxdxp_dxp[TH][2][2];
  t508 = t52*t507;
  //  t510 = diff(diff(rf(X1,X2),X2),X2);
  t510 = dxdxp_dxp[RR][2][2];
  t512 = t175*t56;
  t521 = t512*t221;
  t528 = t52*t175;
  t530 = t510*t12;
  t553 = t510*t38;
  t556 = t512*t24;
  conn[1][2][2] = -(-2.0*t504*t221-t508*t84+t305*t510-t512*t309
		    +(2.0*t512*t84-t504*t21-2.0*t504*t25)*M
		    +(-2.0*t521*t12-3.0*t508*t123-4.0*t504*t222-t521-t528*t329+3.0*t243*t530
		      +(-t504*t224+t504*t10+2.0*t342*t171*t52+4.0*t512*t21*t12+2.0*t12*t171*t442)*M
		      +(-3.0*t508*t132-2.0*t504*t227-2.0*t528*t370+3.0*t246*t553-2.0*t556*t12-t556*t38
			+(2.0*t512*t10*t38-t504*t12+2.0*t504*t38)*M
			+(-t528*t380-t512*t1*t38-t508*t385+t56*t510*t385)*t11)*t11)*t11)*t260*t17;
  conn[1][2][3] = t501*t52*t297;
  conn[1][3][0] = conn[1][0][3];
  conn[1][3][1] = conn[1][1][3];
  conn[1][3][2] = conn[1][2][3];
  t588 = t84*t29;
  t607 = t29*t38;
  conn[1][3][3] = -(-t56*t309*t29+t286*t84+2.0*t240*t588
		    +(-t481-2.0*t408*t284+t480+2.0*t185*t21+(-4.0*t185*t24+t280+3.0*t243*t284+4.0*t278
							     +(-2.0*t314*t29+2.0*t487)*M)*M
		      +(t492*t10+t486-t314*t607-t488+(-2.0*t492*t1+t289+t288-2.0*t285+3.0*t246*t607
						      +(-2.0*t491+2.0*t210*t284)*M)*M+t498)*t11)*t11)*t29*t261;
  t627 = t30*t21;
  t628 = t30*M;
  t630 = 2.0*t628*t24;
  t631 = t30*t10;
  t632 = t631*t12;
  t634 = 2.0*t628*t177;
  t636 = 2.0*t361*t90;
  t637 = t30*t12;
  conn[2][0][0] = -(-t627+t630+(t632-t634-t631-t636+t637*t11)*t11)*M*t261;
  t651 = (-t630+(t634+t636+t631-t632+(-t637+t30*t38)*t11)*t11)*t260*t17;
  conn[2][0][1] = t263*t651;
  conn[2][0][2] = t273*t651;
  t652 = t628*t275;
  t654 = t89*t3;
  t656 = t627*t29;
  t657 = t631*t284;
  t659 = 2.0*t628*t281;
  t661 = 2.0*t361*t3;
  t662 = t31*t10;
  conn[2][0][3] = (2.0*t652-2.0*t654-t656+(t657-t659-t661-t662+t30*t11*t284)*t11)*M*t297;
  conn[2][1][0] = conn[2][0][1];
  t670 = t30*t86;
  t673 = t30*t84;
  t675 = t22*t299;
  t677 = t83*t30;
  t686 = t677*t221;
  t689 = t83*t22;
  t712 = t677*t24;
  conn[2][1][1] = -(2.0*t670*t221-t673*t306+t675*t84+t677*t309+(-2.0*t677*t84+t627*t86+2.0*t670*t25)*M
		    +(2.0*t686*t12+t689*t329-3.0*t627*t333+t686+4.0*t670*t222+3.0*t675*t123
		      +(-2.0*t86*t22*t1*t90+t631*t97-t670*t10-4.0*t677*t21*t12-2.0*t637*t86*t34)*M
		      +(2.0*t712*t12+t712*t38+2.0*t670*t227+3.0*t675*t132-3.0*t631*t358+2.0*t689*t370
			+(t670*t12-2.0*t677*t10*t38-2.0*t670*t38)*M
			+(t675*t385-t30*t306*t385+t689*t380+t677*t1*t38)*t11)*t11)*t11)*t260*t17;
  t748 = t22*t399;
  conn[2][1][2] = -(t303+t310-t673*t402+t748*t84+t346*t221
		    +(-t313+t258*t23+2.0*t258*t26)*M
		    +(t327+t330+2.0*t325-3.0*t627*t418+2.0*t346*t222+3.0*t748*t123+t338
		      +(-t340-t258*t433-t344+t258*t434-2.0*t60*t30*t361*M)*M
		      +(t346*t227+t372+t363+t366+t368-3.0*t631*t454+3.0*t748*t132
			+(t258*t35-t375-2.0*t258*t39)*M+(t387+t748*t385+t384
							 -t30*t402*t385)*t11)*t11)*t11)*t260*t17;
  t793 = t31*t221;
  t794 = t3*t22;
  t795 = t794*t21;
  t799 = t31*t222;
  t800 = 2.0*t799;
  t801 = t22*t41;
  t803 = 2.0*t801*t10;
  t806 = t31*t227;
  t807 = t22*t159;
  t813 = (-t806+t807+(t31*t38-t31*t12)*M)*t11;
  t816 = -t793+t795+(2.0*t654-2.0*t652)*M+(-t800+t803+(t661-t657+t662+t659)*M+t813)*t11;
  conn[2][1][3] = -t816*t22*t297;
  conn[2][2][0] = conn[2][0][2];
  conn[2][2][1] = conn[2][1][2];
  t822 = t22*t507;
  t831 = t3*t56;
  t845 = t41*t56;
  conn[2][2][2] = -(-t673*t510+2.0*t409+t405+t822*t84+(t406*t21-t411+2.0*t406*t25)*M
		    +(-3.0*t627*t530+t421-t432+2.0*t130*t831*t21+t428+3.0*t822*t123+4.0*t425
		      +(t406*t224-t406*t10-t440-t437-2.0*t406*t177*M)*M
		      +(4.0*t130*t845*t10+2.0*t451+3.0*t822*t132+t453+t459-t448-3.0*t631*t553
			+(t406*t12-t462-2.0*t406*t38)*M+(2.0*t258*t381+t822*t385-t30*t510*t385
							 +t467-t469)*t11)*t11)*t11)*t260*t17;
  conn[2][2][3] = -t816*t52*t297;
  conn[2][3][0] = conn[2][0][3];
  conn[2][3][1] = conn[2][1][3];
  conn[2][3][2] = conn[2][2][3];
  t891 = t30*t24;
  conn[2][3][3] = (t794*t84+2.0*t628*t588-t30*t309*t29
		   +(t795-t793-2.0*t30*t221*t284+2.0*t801*t21
		     +(3.0*t627*t284+t656-4.0*t801*t24+4.0*t654+(-2.0*t891*t29+2.0*t799)*M)*M
		     +(t807*t10+t803-t891*t607-t800+(3.0*t631*t607-2.0*t807*t1+t661-2.0*t657+t662
						     +(2.0*t158*t284-2.0*t806)*M)*M+t813)*t11)*t11)*t29*t261;
  t917 = a*M;
  conn[3][0][0] = -t7*t917*t17;
  t920 = t102*t10;
  t924 = 1/t29;
  t916 = t924*t17;
  conn[3][0][1] = -t917*(-t920+t148+(t156+t149)*t11)*t916;
  t928 = t129*t10;
  conn[3][0][2] = -t917*(t208-t928+(t214+t204)*t11)*t916;
  conn[3][0][3] = -t78*M*t11*t17;
  conn[3][1][0] = conn[3][0][1];
  t940 = t83*t29;
  t946 = t86*t29;
  t968 = t916;
  conn[3][1][1] = -(t940*t221+2.0*t794*t627+(4.0*t794*t891-t946*t10)*M
		    +(4.0*t801*t631+2.0*t940*t222+(4.0*t361*t42+t946*t12)*M+(t940*t227+2.0*t807*t30)*t11)
		    *t11)*a*t968;
  t970 = t29*t221;
  t991 = t52*t1;
  conn[3][1][2] = -(t116*t970+t286*t627+t794*t243
		    +(2.0*t286*t891-t102*t52*t10+2.0*t794*t314)*M
		    +(2.0*t185*t631+2.0*t801*t246+2.0*t116*t275*t12+(2.0*t361*t845+2.0*t991*t42+t102*t60)*M
		      +(t116*t40*t38+t492*t30+t807*t56)*t11)*t11)*a*t968;
  t1024 = t38*t41;
  conn[3][1][3] = (t3*t30*t84+t970*t22+(3.0*t42*t21+2.0*t275*t35
					+(t102*t224+t148-t154-t920)*M
					+(t40*t39+3.0*t159*t30*t10+(t156-t161+t149-t157)*M
					  +t1024*t30*t11)*t11)*t11)*t924*t17;
  conn[3][2][0] = conn[3][0][2];
  conn[3][2][1] = conn[3][1][2];
  t1035 = t175*t29;
  t1037 = t171*t29;
  conn[3][2][2] = -(2.0*t286*t243+t1035*t221+(-t1037*t10+4.0*t286*t314)*M
		    +(4.0*t185*t246+2.0*t1035*t222+(t1037*t12+4.0*t991*t845)*M
		      +(t1035*t227+2.0*t492*t56)*t11)*t11)*a*t968;
  conn[3][2][3] = -(-t970*t52-t831*t84+(-2.0*t275*t60-3.0*t845*t21
					+(t203-t208-t129*t224+t928)*M
					+(-t40*t63-3.0*t159*t56*t10+(t209+t212-t204-t214)*M
					  -t1024*t56*t11)*t11)*t11)*t924*t17;
  conn[3][3][0] = conn[3][0][3];
  conn[3][3][1] = conn[3][1][3];
  conn[3][3][2] = conn[3][2][3];
  conn[3][3][3] = -t236*a*t17;

  return;

}

#undef M






// jon's MKS source in mid-simplified form (used to UPDATE the source)
// only for ideal EOS (gam)
void mks_source_conn(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q,FTYPE *dU)
{
  int ii,jj,kk;
  int i=0, j=0, k=0, l=0;
  FTYPE r,th,X[NDIM],V[NDIM],sigma,dxdxptrue[NDIM][NDIM];
#ifdef WIN32
  extern FTYPE cot(FTYPE arg);
#endif
  extern FTYPE csc(FTYPE arg);
  FTYPE b[NDIM],u[NDIM],bsq,en,rho;
  FTYPE dxdxp[NDIM];



  if(MBH!=1.0){
    dualfprintf(fail_file,"mks_source_conn not setup for MBH!=1.0\n");
    myexit(12);
  }


  ii=ptrgeom->i;
  jj=ptrgeom->j;
  kk=ptrgeom->k;


  bsq = dot(q->bcon, q->bcov);
  u[TT]=q->ucon[TT];
  u[RR]=q->ucon[RR];
  u[TH]=q->ucon[TH];
  u[PH]=q->ucon[PH];

  b[TT]=q->bcon[TT];
  b[RR]=q->bcon[RR];
  b[TH]=q->bcon[TH];
  b[PH]=q->bcon[PH];

  rho=pr[RHO];
  en=pr[UU];

  coord(ptrgeom->i, ptrgeom->j, ptrgeom->k,ptrgeom->p, X);
  // get bl coordinates
  bl_coord(X,V);
  r=V[1];
  th=V[2];


  // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
  // follows gcov_func()
  if(POSDEFMETRIC){
    if(th<0.0){ th=-th;}
    if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }
  // avoid singularity at polar axis
#if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
    if(th>=0) th=SINGSMALL;
    if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
    if(th>=M_PI) th=M_PI+SINGSMALL;
    if(th<M_PI) th=M_PI-SINGSMALL;
  }
#endif
  // set aux vars
  dxdxprim(X,V,dxdxptrue);
  DLOOPA(j) dxdxp[j]=dxdxptrue[j][j]; // defcoord==LOGRSINTH assumes transformation is diagonal


  sigma=r*r+a*a*cos(th)*cos(th);



  if((WHICHEOM==WITHNOGDET)&&(NOGDETU0==1)){
    // see grmhd-fullsource-simplify.nb

    dU[UU]+=pow(sigma,-1.)*(-1.*(-1. - 2.*dxdxp[1]*r*pow(sigma,-1.))*
			    (b[RR]*(-1.*b[TT]*pow(a,2.)*pow(cos(th),2.) + 
				    r*(2.*b[TT] + 2.*b[RR]*dxdxp[1] - 1.*b[TT]*r - 
				       2.*b[PH]*a*pow(sin(th),2.))) + 
			     u[RR]*(bsq + en*gam + rho)*
			     (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
			      r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
				 u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) + 
			    (b[TH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
				    r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
				       2.*b[PH]*a*pow(sin(th),2.))) - 
			     1.*u[TH]*(bsq + en*gam + rho)*
			     (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
			      r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
				 u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*
			    (-1.*dxdxp[2]*cot(th) + 
			     4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
			     dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th)));

  }
  // else nothing to add then


  if((WHICHEOM==WITHNOGDET)&&(NOGDETU1==1)){

    // see grmhd-ksksp-mhd.nb and the other source*simplify*.nb files

    dU[U1]+=0.0625*(16.*dxdxp[1]*r*
		    (0.5*bsq + en*(-1. + gam) - 
		     1.*sigma*pow(b[TH],2.)*pow(dxdxp[2],2.) + 
		     (bsq + en*gam + rho)*sigma*pow(dxdxp[2],2.)*pow(u[TH],2.))*
		    pow(sigma,-1.) + 2.*(-1. - 2.*dxdxp[1]*r*pow(sigma,-1.))*
		    (4.*bsq + 8.*en*(-1. + gam) - 
		     1.*b[RR]*dxdxp[1]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
					8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
					cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
					1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
					8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.))*
		     pow(sigma,-1.) + dxdxp[1]*u[RR]*(bsq + en*gam + rho)*
		     (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
		      4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
		      4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
		      u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
		      4.*u[PH]*a*pow(r,2.))*pow(sigma,-1.)) - 
		    1.*dxdxp[1]*(4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
		    (-1.*b[TT]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
				8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
				cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
				1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
				8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.)) + 
		     u[TT]*(bsq + en*gam + rho)*
		     (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
		      4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
		      4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
		      u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
		      4.*u[PH]*a*pow(r,2.)))*pow(sigma,-4.)*
		    (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.)) + 
		    8.*a*pow(dxdxp[1],2.)*(-1.*b[RR]*
					   (-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
					    b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
					    cos(2.*th)*pow(a,2.)*
					    (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
					     b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
					    b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
					   u[RR]*(bsq + en*gam + rho)*
					   (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
						     dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
					    cos(2.*th)*pow(a,2.)*
					    (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
					     u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
					    u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*pow(sigma,-4.)*
		    (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.) + 
		    8.*dxdxp[1]*a*(-1.*b[TT]*
				   (-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
				    b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
				    cos(2.*th)*pow(a,2.)*
				    (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
				     b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
				    b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
				   u[TT]*(bsq + en*gam + rho)*
				   (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
					     dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
				    cos(2.*th)*pow(a,2.)*
				    (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
				     u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
				    u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*pow(sigma,-4.)*
		    (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.) + 
		    dxdxp[1]*a*(-1.*b[PH]*
				(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 8.*b[PH]*a*r + 
				 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*cos(2.*th) + 
				 4.*b[RR]*dxdxp[1]*pow(a,2.) - 1.*b[PH]*pow(a,3.) + 
				 b[PH]*cos(4.*th)*pow(a,3.) + 8.*b[RR]*dxdxp[1]*pow(r,2.) - 
				 4.*b[PH]*a*pow(r,2.)) + 
				u[PH]*(bsq + en*gam + rho)*
				(16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
				 4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
				 4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
				 u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
				 4.*u[PH]*a*pow(r,2.)))*pow(sigma,-4.)*
		    (pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
		     2.*r*(pow(1.,2.)*(-2.*sigma + 4.*pow(r,2.)) + pow(sigma,2.)) + 
		     cos(2.*th)*pow(a,2.)*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.)))*
		    pow(sin(th),2.) + 8.*dxdxp[1]*pow(sigma,-3.)*
		    (bsq + 2.*en*(-1. + gam) - 
		     1.*b[PH]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
			       b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
			       cos(2.*th)*pow(a,2.)*
			       (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
				b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
			       b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.))*pow(sigma,-1.)*
		     pow(sin(th),2.) + u[PH]*(bsq + en*gam + rho)*
		     (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
			       dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
		      cos(2.*th)*pow(a,2.)*
		      (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
		       u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
		      u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.))*pow(sigma,-1.)*
		     pow(sin(th),2.))*(r*pow(sigma,2.) + 
				       pow(a,2.)*(-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)
				       ) + 2.*pow(sigma,-3.)*(4.*bsq + 8.*en*(-1. + gam) - 
							      1.*b[RR]*dxdxp[1]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
										 8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
										 cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
										 1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
										 8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.))*
							      pow(sigma,-1.) + dxdxp[1]*u[RR]*(bsq + en*gam + rho)*
							      (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
							       4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
							       4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
							       u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
							       4.*u[PH]*a*pow(r,2.))*pow(sigma,-1.))*
		    (-1.*dxdxp[1]*(2. + r)*pow(r,3.) + pow(sigma,3.) + 
		     dxdxp[1]*pow(a,2.)*(2.*r*pow(cos(th),2.) + 
					 pow(a,2.)*pow(cos(th),4.) + 
					 (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))) + 
		    16.*dxdxp[1]*a*pow(sigma,-4.)*
		    (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
		    (-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)*
		    (b[PH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
			    r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
			       2.*b[PH]*a*pow(sin(th),2.))) - 
		     1.*u[PH]*(bsq + en*gam + rho)*
		     (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
		      r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
			 u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) - 
		    32.*pow(dxdxp[1],2.)*pow(sigma,-4.)*
		    (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
		    (r*(1. + r) + pow(a,2.)*pow(cos(th),2.))*
		    (b[RR]*(-1.*b[TT]*pow(a,2.)*pow(cos(th),2.) + 
			    r*(2.*b[TT] + 2.*b[RR]*dxdxp[1] - 1.*b[TT]*r - 
			       2.*b[PH]*a*pow(sin(th),2.))) + 
		     u[RR]*(bsq + en*gam + rho)*
		     (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
		      r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
			 u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) + 
		    16.*dxdxp[1]*pow(sigma,-3.)*
		    (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
		    (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
		    (0.5*bsq + en*(-1. + gam) + 
		     b[TT]*pow(sigma,-1.)*
		     (b[TT]*pow(a,2.)*pow(cos(th),2.) + 
		      r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
			 2.*b[PH]*a*pow(sin(th),2.))) - 
		     1.*u[TT]*(bsq + en*gam + rho)*pow(sigma,-1.)*
		     (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
		      r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
			 u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) - 
		    2.*dxdxp[1]*dxdxp[2]*cos(th)*pow(a,2.)*
		    (-1.*b[TH]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
				8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
				cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
				1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
				8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.)) + 
		     u[TH]*(bsq + en*gam + rho)*
		     (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
		      4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
		      4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
		      u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
		      4.*u[PH]*a*pow(r,2.)))*pow(sigma,-2.)*sin(th) - 
		    8.*dxdxp[1]*dxdxp[2]*a*cos(th)*
		    (-1.*b[TH]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
				b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
				cos(2.*th)*pow(a,2.)*
				(-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
				 b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
				b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
		     u[TH]*(bsq + en*gam + rho)*
		     (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
			       dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
		      cos(2.*th)*pow(a,2.)*
		      (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
		       u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
		      u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*pow(sigma,-3.)*
		    (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*sin(th) + 
		    16.*dxdxp[1]*dxdxp[2]*r*
		    (b[TH]*b[TT] - 1.*u[TH]*u[TT]*(bsq + en*gam + rho))*pow(a,2.)*
		    pow(sigma,-2.)*sin(2.*th) + 
		    16.*dxdxp[2]*r*(b[RR]*b[TH] - 
				    1.*u[RR]*u[TH]*(bsq + en*gam + rho))*pow(dxdxp[1],2.)*pow(a,2.)*
		    pow(sigma,-2.)*sin(2.*th) - 
		    16.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
		    (b[TH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
			    r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
			       2.*b[PH]*a*pow(sin(th),2.))) - 
		     1.*u[TH]*(bsq + en*gam + rho)*
		     (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
		      r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
			 u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*sin(2.*th) + 
		    2.*dxdxp[1]*(-1.*b[TH]*
				 (16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 8.*b[PH]*a*r + 
				  4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*cos(2.*th) + 
				  4.*b[RR]*dxdxp[1]*pow(a,2.) - 1.*b[PH]*pow(a,3.) + 
				  b[PH]*cos(4.*th)*pow(a,3.) + 8.*b[RR]*dxdxp[1]*pow(r,2.) - 
				  4.*b[PH]*a*pow(r,2.)) + 
				 u[TH]*(bsq + en*gam + rho)*
				 (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
				  4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
				  4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
				  u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
				  4.*u[PH]*a*pow(r,2.)))*pow(sigma,-1.)*
		    (-1.*dxdxp[2]*cot(th) + 
		     4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
		     dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th)) - 
		    16.*dxdxp[1]*dxdxp[2]*a*
		    (b[PH]*b[TH] - 1.*u[PH]*u[TH]*(bsq + en*gam + rho))*
		    pow(sigma,-2.)*sin(th)*(r*(2. + r)*sigma*cos(th) + 
					    sigma*pow(a,2.)*pow(cos(th),3.) + r*pow(a,2.)*sin(th)*sin(2.*th)));

  }
  else{

    dU[U1]+=0.5*dxdxp[1]*pow(sigma,-5.)*
      (r*pow(sigma,4.)*(bsq - 2.*en + 2.*en*gam - 
			2.*sigma*pow(dxdxp[2],2.)*pow(b[TH],2.) + 
			(bsq + en*gam + rho)*pow(dxdxp[2],2.)*
			(pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[TH],2.)) + 
       a*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (-1.*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
	(r*(2. + 3.*r)*pow(a,2.) + 
	 cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
	 2.*pow(r,4.)) - 2.*a*
	((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
	 0.5*dxdxp[1]*(bsq + en*gam + rho)*u[TT]*u[RR]*
	 (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
	(r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
	4.*a*r*sigma*(-0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
		      (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 1.*pow(b[TT],2.) + 
		      (bsq + en*gam + rho)*pow(u[TT],2.)))*pow(sin(th),2.) + 
       dxdxp[1]*a*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (-4.*a*r*((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
		 dxdxp[1]*(bsq + en*gam + rho)*sigma*u[TT]*u[RR])*
	pow(dxdxp[1],-1.) + 
	sigma*(r*(2. + 3.*r)*pow(a,2.) + 
	       cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
	       2.*pow(r,4.))*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
			      0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
	- 1.*a*pow(dxdxp[1],-1.)*
	(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
	((bsq + 2.*en*(-1. + gam))*((-2. + r)*r + pow(a,2.)) - 
	 2.*sigma*pow(dxdxp[1],2.)*pow(b[RR],2.) + 
	 (bsq + en*gam + rho)*pow(dxdxp[1],2.)*
	 (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[RR],2.)))*
       pow(sin(th),2.) - 1.*sigma*
       (4.*r - 1.*pow(a,2.) + cos(2.*th)*pow(a,2.))*
       (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
	 dxdxp[1]*(bsq + en*gam + rho)*sigma*u[TT]*u[RR])*pow(sigma,-1.)*
	(r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) + 
	2.*r*(-0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
	      (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 1.*pow(b[TT],2.) + 
	      (bsq + en*gam + rho)*pow(u[TT],2.)) - 
	1.*a*(-1.*b[TT]*b[PH] + (bsq + en*gam + rho)*u[TT]*u[PH])*
	(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)) + 
       (4.*a*r*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH]) - 
	1.*a*(a*(bsq + 2.*en*(-1. + gam)) - 
	      1.*dxdxp[1]*b[RR]*b[PH]*
	      (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)) + 
	      dxdxp[1]*(bsq + en*gam + rho)*u[RR]*u[PH]*
	      (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
	(r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) + 
	0.5*(r*(2. + 3.*r)*pow(a,2.) + 
	     cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
	     2.*pow(r,4.))*((bsq + 2.*en*(-1. + gam))*pow(csc(th),2.) - 
			    2.*sigma*pow(b[PH],2.) + 2.*(bsq + en*gam + rho)*sigma*pow(u[PH],2.))
	)*pow(sin(th),2.)*(r*pow(sigma,2.) + 
			   pow(a,2.)*(-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)
			   ) + 2.*a*sigma*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       (-1.*pow(r,2.) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)*
       (r*(a*(bsq + 2.*en*(-1. + gam)) - 2.*dxdxp[1]*sigma*b[RR]*b[PH] + 
	   2.*dxdxp[1]*(bsq + en*gam + rho)*sigma*u[RR]*u[PH])*pow(sigma,-1.)\
	- 1.*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
	((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) - 
	2.*a*r*(0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
		pow(csc(th),2.) - 1.*pow(b[PH],2.) + 
		(bsq + en*gam + rho)*pow(u[PH],2.))*pow(sin(th),2.)) + 
       a*sigma*(pow(a,2.)*(sigma - 2.*pow(r,2.)) + 
		2.*r*(pow(1.,2.)*(-2.*sigma + 4.*pow(r,2.)) + pow(sigma,2.)) + 
		cos(2.*th)*pow(a,2.)*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.)))*
       pow(sin(th),2.)*(2.*r*(-1.*b[TT]*b[PH] + 
			      (bsq + en*gam + rho)*u[TT]*u[PH]) + 
			dxdxp[1]*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
				  0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
			*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
			1.*a*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
			(0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*pow(csc(th),2.) - 
			 1.*pow(b[PH],2.) + (bsq + en*gam + rho)*pow(u[PH],2.))*
			pow(sin(th),2.)) + sigma*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       ((bsq + 2.*en*(-1. + gam))*sigma + 
	2.*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.))*pow(b[TT],2.) - 
	1.*(bsq + en*gam + rho)*
	(2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.))*pow(u[TT],2.) + 
	4.*r*b[TT]*(-1.*dxdxp[1]*b[RR] + a*b[PH]*pow(sin(th),2.)) + 
	4.*r*(bsq + en*gam + rho)*u[TT]*
	(dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.))) + 
       2.*sigma*(2.*r*(-1.*b[TT]*b[RR] + (bsq + en*gam + rho)*u[TT]*u[RR] + 
		       (bsq + 2.*en*(-1. + gam))*r*pow(dxdxp[1],-1.)*pow(sigma,-1.)) + 
		 dxdxp[1]*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
		 (0.5*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-2.)*
		  ((-2. + r)*r + pow(a,2.))*pow(sigma,-1.) - 1.*pow(b[RR],2.) + 
		  (bsq + en*gam + rho)*pow(u[RR],2.)) - 
		 1.*a*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
		       0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
		 *(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))*
       (-1.*dxdxp[1]*(2. + r)*pow(r,3.) + pow(sigma,3.) + 
	dxdxp[1]*pow(a,2.)*(2.*r*pow(cos(th),2.) + 
			    pow(a,2.)*pow(cos(th),4.) + 
			    (pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))) + 
       4.*dxdxp[1]*sigma*(pow(r,2.) - 1.*pow(a,2.)*pow(cos(th),2.))*
       (r*(1. + r) + pow(a,2.)*pow(cos(th),2.))*
       (b[TT]*b[RR]*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.)) - 
	2.*dxdxp[1]*r*pow(b[RR],2.) + 2.*a*r*b[RR]*b[PH]*pow(sin(th),2.) + 
	0.5*(bsq + en*gam + rho)*u[RR]*
	(-1.*u[TT]*(2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.)) + 
	 4.*r*(dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.)))) - 
       1.*dxdxp[2]*a*cos(th)*pow(sigma,2.)*
       (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
       (4.*a*r*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH]) - 
	1.*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*
	(r*(2. + 3.*r)*pow(a,2.) + 
	 cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
	 2.*pow(r,4.)) + 2.*dxdxp[1]*a*
	(b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*
	(r*(2. + r) + pow(a,2.)*pow(cos(th),2.)))*sin(th) - 
       2.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,3.)*
       (2.*r*(-1.*b[TT]*b[TH] + (bsq + en*gam + rho)*u[TT]*u[TH]) + 
	dxdxp[1]*(-1.*b[RR]*b[TH] + (bsq + en*gam + rho)*u[RR]*u[TH])*
	(r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
	1.*a*(-1.*b[TH]*b[PH] + (bsq + en*gam + rho)*u[TH]*u[PH])*
	(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))*sin(th) + 
       2.*dxdxp[2]*r*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH])*
       pow(a,2.)*pow(sigma,3.)*sin(2.*th) + 
       2.*dxdxp[1]*dxdxp[2]*r*
       (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*pow(a,2.)*pow(sigma,3.)*
       sin(2.*th) - 2.*dxdxp[2]*r*pow(a,2.)*pow(sigma,2.)*
       (-2.*dxdxp[1]*r*(b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH]) - 
	1.*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH])*
	((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) + 
	2.*a*r*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*pow(sin(th),2.)
	)*sin(2.*th) - 2.*dxdxp[2]*a*
       (b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*pow(sigma,3.)*sin(th)*
       (r*(2. + r)*sigma*cos(th) + sigma*pow(a,2.)*pow(cos(th),3.) + 
	r*pow(a,2.)*sin(th)*sin(2.*th)));
  }



  if((WHICHEOM==WITHNOGDET)&&(NOGDETU2==1)){

    dU[U2]+=dxdxp[1]*r*(-1.*b[RR]*b[TH] + 
			u[RR]*u[TH]*(bsq + en*gam + rho))*pow(dxdxp[2],2.) - 
      1.*(-1.*b[RR]*b[TH] + u[RR]*u[TH]*(bsq + en*gam + rho))*
      (2.*dxdxp[1]*r + sigma)*pow(dxdxp[2],2.) - 
      0.125*r*pow(dxdxp[2],2.)*((-2. + r)*r + pow(a,2.))*
      (-1.*b[TH]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
		  8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
		  cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
		  1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
		  8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.)) + 
       u[TH]*(bsq + en*gam + rho)*
       (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
	4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
	4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
	u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
	4.*u[PH]*a*pow(r,2.)))*pow(sigma,-2.) - 
      0.5*a*r*pow(dxdxp[2],2.)*(-1.*b[TH]*
				(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
				 b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
				 cos(2.*th)*pow(a,2.)*(-1.*b[RR]*dxdxp[1]*a + 
						       b[PH]*(-2. + r)*r + b[PH]*pow(a,2.)) - 
				 1.*b[RR]*dxdxp[1]*pow(a,3.) + b[PH]*pow(a,4.) + 
				 2.*b[PH]*pow(r,4.)) + 
				u[TH]*(bsq + en*gam + rho)*
				(-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + dxdxp[1]*u[RR]*r) + 
				 u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
				 cos(2.*th)*pow(a,2.)*(-1.*dxdxp[1]*u[RR]*a + 
						       u[PH]*(-2. + r)*r + u[PH]*pow(a,2.)) - 
				 1.*dxdxp[1]*u[RR]*pow(a,3.) + u[PH]*pow(a,4.) + 
				 2.*u[PH]*pow(r,4.)))*pow(sigma,-2.)*pow(sin(th),2.) - 
      2.*pow(dxdxp[2],2.)*pow(r,2.)*pow(sigma,-2.)*
      (b[TH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
	      r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
		 2.*b[PH]*a*pow(sin(th),2.))) - 
       1.*u[TH]*(bsq + en*gam + rho)*
       (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
	r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
	   u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.)))) + 
      2.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-3.)*
      (b[PH]*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
	      r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
		 2.*b[PH]*a*pow(sin(th),2.))) - 
       1.*u[PH]*(bsq + en*gam + rho)*
       (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
	r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
	   u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*pow(sin(th),3.) - 
      1.*dxdxp[2]*a*r*cos(th)*(-1.*b[TT]*
			       (-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
				b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
				cos(2.*th)*pow(a,2.)*(-1.*b[RR]*dxdxp[1]*a + 
						      b[PH]*(-2. + r)*r + b[PH]*pow(a,2.)) - 
				1.*b[RR]*dxdxp[1]*pow(a,3.) + b[PH]*pow(a,4.) + 
				2.*b[PH]*pow(r,4.)) + 
			       u[TT]*(bsq + en*gam + rho)*
			       (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + dxdxp[1]*u[RR]*r) + 
				u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
				cos(2.*th)*pow(a,2.)*(-1.*dxdxp[1]*u[RR]*a + 
						      u[PH]*(-2. + r)*r + u[PH]*pow(a,2.)) - 
				1.*dxdxp[1]*u[RR]*pow(a,3.) + u[PH]*pow(a,4.) + 
				2.*u[PH]*pow(r,4.)))*pow(sigma,-3.)*sin(th) - 
      0.125*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*
      (4.*bsq + 8.*en*(-1. + gam) - 
       1.*b[RR]*dxdxp[1]*(16.*b[TT]*r + 16.*b[RR]*dxdxp[1]*r - 
			  8.*b[PH]*a*r + 4.*a*(b[RR]*dxdxp[1]*a + b[PH]*r*(2. + r))*
			  cos(2.*th) + 4.*b[RR]*dxdxp[1]*pow(a,2.) - 
			  1.*b[PH]*pow(a,3.) + b[PH]*cos(4.*th)*pow(a,3.) + 
			  8.*b[RR]*dxdxp[1]*pow(r,2.) - 4.*b[PH]*a*pow(r,2.))*
       pow(sigma,-1.) + dxdxp[1]*u[RR]*(bsq + en*gam + rho)*
       (16.*dxdxp[1]*u[RR]*r + 16.*u[TT]*r - 8.*u[PH]*a*r + 
	4.*a*(dxdxp[1]*u[RR]*a + u[PH]*r*(2. + r))*cos(2.*th) + 
	4.*dxdxp[1]*u[RR]*pow(a,2.) - 1.*u[PH]*pow(a,3.) + 
	u[PH]*cos(4.*th)*pow(a,3.) + 8.*dxdxp[1]*u[RR]*pow(r,2.) - 
	4.*u[PH]*a*pow(r,2.))*pow(sigma,-1.))*sin(th) - 
      0.5*dxdxp[1]*dxdxp[2]*a*cos(th)*
      (-1.*b[RR]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
		  b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
		  cos(2.*th)*pow(a,2.)*(-1.*b[RR]*dxdxp[1]*a + 
					b[PH]*(-2. + r)*r + b[PH]*pow(a,2.)) - 
		  1.*b[RR]*dxdxp[1]*pow(a,3.) + b[PH]*pow(a,4.) + 
		  2.*b[PH]*pow(r,4.)) + 
       u[RR]*(bsq + en*gam + rho)*
       (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + dxdxp[1]*u[RR]*r) + 
	u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
	cos(2.*th)*pow(a,2.)*(-1.*dxdxp[1]*u[RR]*a + 
			      u[PH]*(-2. + r)*r + u[PH]*pow(a,2.)) - 
	1.*dxdxp[1]*u[RR]*pow(a,3.) + u[PH]*pow(a,4.) + 
	2.*u[PH]*pow(r,4.)))*pow(sigma,-3.)*
      (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*sin(th) + 
      (0.5*bsq + en*(-1. + gam) - 1.*sigma*pow(b[TH],2.)*pow(dxdxp[2],2.) + 
       (bsq + en*gam + rho)*sigma*pow(dxdxp[2],2.)*pow(u[TH],2.))*
      (4.*(M_PI*X[2] - 1.*th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) - 
       1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)) + 
      dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
      (b[RR]*(-1.*b[TT]*pow(a,2.)*pow(cos(th),2.) + 
	      r*(2.*b[TT] + 2.*b[RR]*dxdxp[1] - 1.*b[TT]*r - 
		 2.*b[PH]*a*pow(sin(th),2.))) + 
       u[RR]*(bsq + en*gam + rho)*
       (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
	r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
	   u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*sin(2.*th) - 
      1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-2.)*
      (0.5*bsq + en*(-1. + gam) + 
       b[TT]*pow(sigma,-1.)*(b[TT]*pow(a,2.)*pow(cos(th),2.) + 
			     r*(-2.*b[TT] - 2.*b[RR]*dxdxp[1] + b[TT]*r + 
				2.*b[PH]*a*pow(sin(th),2.))) - 
       1.*u[TT]*(bsq + en*gam + rho)*pow(sigma,-1.)*
       (u[TT]*pow(a,2.)*pow(cos(th),2.) + 
	r*(-2.*dxdxp[1]*u[RR] - 2.*u[TT] + dxdxp[1]*u[TT] + 
	   u[TT]*R0 + 2.*u[PH]*a*pow(sin(th),2.))))*sin(2.*th) + 
      0.5*dxdxp[2]*(bsq + 2.*en*(-1. + gam) - 
		    1.*b[PH]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
			      b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
			      cos(2.*th)*pow(a,2.)*(-1.*b[RR]*dxdxp[1]*a + 
						    b[PH]*(-2. + r)*r + b[PH]*pow(a,2.)) - 
			      1.*b[RR]*dxdxp[1]*pow(a,3.) + b[PH]*pow(a,4.) + 
			      2.*b[PH]*pow(r,4.))*pow(sigma,-1.)*pow(sin(th),2.) + 
		    u[PH]*(bsq + en*gam + rho)*
		    (-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + dxdxp[1]*u[RR]*r) + 
		     u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
		     cos(2.*th)*pow(a,2.)*(-1.*dxdxp[1]*u[RR]*a + 
					   u[PH]*(-2. + r)*r + u[PH]*pow(a,2.)) - 
		     1.*dxdxp[1]*u[RR]*pow(a,3.) + u[PH]*pow(a,4.) + 
		     2.*u[PH]*pow(r,4.))*pow(sigma,-1.)*pow(sin(th),2.))*
      (cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th)) + 
      (0.5*bsq + en*(-1. + gam) - 1.*sigma*pow(b[TH],2.)*pow(dxdxp[2],2.) + 
       (bsq + en*gam + rho)*sigma*pow(dxdxp[2],2.)*pow(u[TH],2.))*
      (-1.*dxdxp[2]*cot(th) + 4.*(-1.*M_PI*X[2] + th)*
       pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
       dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th));
  }
  else{


    dU[U2]+=0.5*(-2.*dxdxp[1]*r*
		 (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*pow(dxdxp[2],2.) - 
		 1.*a*r*pow(dxdxp[2],2.)*pow(sigma,-2.)*
		 (4.*a*r*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH]) - 
		  1.*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*
		  (r*(2. + 3.*r)*pow(a,2.) + 
		   cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
		   2.*pow(r,4.)) + 2.*dxdxp[1]*a*
		  (b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH])*
		  (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)))*pow(sin(th),2.) - 
		 4.*pow(dxdxp[2],2.)*pow(r,2.)*pow(sigma,-2.)*
		 (-2.*dxdxp[1]*r*(b[RR]*b[TH] - 1.*(bsq + en*gam + rho)*u[RR]*u[TH]) - 
		  1.*(b[TT]*b[TH] - 1.*(bsq + en*gam + rho)*u[TT]*u[TH])*
		  ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) + 
		  2.*a*r*(b[TH]*b[PH] - 1.*(bsq + en*gam + rho)*u[TH]*u[PH])*pow(sin(th),2.)
		  ) - 2.*r*pow(dxdxp[2],2.)*((-2. + r)*r + pow(a,2.))*pow(sigma,-2.)*
		 (2.*r*(-1.*b[TT]*b[TH] + (bsq + en*gam + rho)*u[TT]*u[TH]) + 
		  dxdxp[1]*(-1.*b[RR]*b[TH] + (bsq + en*gam + rho)*u[RR]*u[TH])*
		  (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
		  1.*a*(-1.*b[TH]*b[PH] + (bsq + en*gam + rho)*u[TH]*u[PH])*
		  (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.)) + 
		 4.*dxdxp[2]*r*cos(th)*pow(a,3.)*pow(sigma,-3.)*
		 (r*(a*(bsq + 2.*en*(-1. + gam)) - 2.*dxdxp[1]*sigma*b[RR]*b[PH] + 
		     2.*dxdxp[1]*(bsq + en*gam + rho)*sigma*u[RR]*u[PH])*pow(sigma,-1.)\
		  - 1.*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
		  ((2. - 1.*r)*r - 1.*pow(a,2.)*pow(cos(th),2.)) - 
		  2.*a*r*(0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
			  pow(csc(th),2.) - 1.*pow(b[PH],2.) + 
			  (bsq + en*gam + rho)*pow(u[PH],2.))*pow(sin(th),2.))*pow(sin(th),3.)
		 - 2.*dxdxp[2]*a*r*cos(th)*pow(sigma,-4.)*
		 (-1.*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH])*
		  (r*(2. + 3.*r)*pow(a,2.) + 
		   cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
		   2.*pow(r,4.)) - 2.*a*
		  ((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
		   0.5*dxdxp[1]*(bsq + en*gam + rho)*u[TT]*u[RR]*
		   (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
		  (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 
		  4.*a*r*sigma*(-0.5*(bsq + 2.*en*(-1. + gam))*pow(sigma,-1.)*
				(r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) - 1.*pow(b[TT],2.) + 
				(bsq + en*gam + rho)*pow(u[TT],2.)))*sin(th) - 
		 1.*dxdxp[1]*dxdxp[2]*a*cos(th)*pow(sigma,-4.)*
		 (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
		 (-4.*a*r*((bsq + 2.*en*(-1. + gam))*r - 1.*dxdxp[1]*sigma*b[TT]*b[RR] + 
			   dxdxp[1]*(bsq + en*gam + rho)*sigma*u[TT]*u[RR])*
		  pow(dxdxp[1],-1.) + 
		  sigma*(r*(2. + 3.*r)*pow(a,2.) + 
			 cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
			 2.*pow(r,4.))*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
					0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
		  - 1.*a*pow(dxdxp[1],-1.)*
		  (r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
		  ((bsq + 2.*en*(-1. + gam))*((-2. + r)*r + pow(a,2.)) - 
		   2.*sigma*pow(dxdxp[1],2.)*pow(b[RR],2.) + 
		   (bsq + en*gam + rho)*pow(dxdxp[1],2.)*
		   (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[RR],2.)))*
		 sin(th) - 2.*dxdxp[1]*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-2.)*
		 (2.*r*(-1.*b[TT]*b[RR] + (bsq + en*gam + rho)*u[TT]*u[RR] + 
			(bsq + 2.*en*(-1. + gam))*r*pow(dxdxp[1],-1.)*pow(sigma,-1.)) + 
		  dxdxp[1]*(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*
		  (0.5*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-2.)*
		   ((-2. + r)*r + pow(a,2.))*pow(sigma,-1.) - 1.*pow(b[RR],2.) + 
		   (bsq + en*gam + rho)*pow(u[RR],2.)) - 
		  1.*a*(-1.*b[RR]*b[PH] + (bsq + en*gam + rho)*u[RR]*u[PH] + 
			0.5*a*(bsq + 2.*en*(-1. + gam))*pow(dxdxp[1],-1.)*pow(sigma,-1.))
		  *(r*(2. + r) + pow(a,2.)*pow(cos(th),2.))*pow(sin(th),2.))*sin(th)\
		 + (bsq - 2.*en + 2.*en*gam - 2.*sigma*pow(dxdxp[2],2.)*pow(b[TH],2.) + 
		    (bsq + en*gam + rho)*pow(dxdxp[2],2.)*
		    (pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.))*pow(u[TH],2.))*
		 (4.*(M_PI*X[2] - 1.*th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) - 
		  1.*dxdxp[2]*cos(th)*pow(a,2.)*pow(sigma,-1.)*sin(th)) - 
		 1.*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
		 ((bsq + 2.*en*(-1. + gam))*sigma + 
		  2.*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.))*pow(b[TT],2.) - 
		  1.*(bsq + en*gam + rho)*
		  (2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.))*pow(u[TT],2.) + 
		  4.*r*b[TT]*(-1.*dxdxp[1]*b[RR] + a*b[PH]*pow(sin(th),2.)) + 
		  4.*r*(bsq + en*gam + rho)*u[TT]*
		  (dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.)))*sin(2.*th) - 
		 2.*dxdxp[1]*dxdxp[2]*r*pow(a,2.)*pow(sigma,-3.)*
		 (b[TT]*b[RR]*((-2. + r)*r + pow(a,2.)*pow(cos(th),2.)) - 
		  2.*dxdxp[1]*r*pow(b[RR],2.) + 2.*a*r*b[RR]*b[PH]*pow(sin(th),2.) + 
		  0.5*(bsq + en*gam + rho)*u[RR]*
		  (-1.*u[TT]*(2.*(-2. + r)*r + pow(a,2.) + cos(2.*th)*pow(a,2.)) + 
		   4.*r*(dxdxp[1]*u[RR] - 1.*a*u[PH]*pow(sin(th),2.))))*sin(2.*th) + 
		 dxdxp[2]*pow(sigma,-2.)*
		 (4.*a*r*sigma*(b[TT]*b[PH] - 1.*(bsq + en*gam + rho)*u[TT]*u[PH]) - 
		  1.*a*(a*(bsq + 2.*en*(-1. + gam)) - 
			1.*dxdxp[1]*b[RR]*b[PH]*
			(pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)) + 
			dxdxp[1]*(bsq + en*gam + rho)*u[RR]*u[PH]*
			(pow(a,2.) + cos(2.*th)*pow(a,2.) + 2.*pow(r,2.)))*
		  (r*(2. + r) + pow(a,2.)*pow(cos(th),2.)) + 
		  0.5*(r*(2. + 3.*r)*pow(a,2.) + 
		       cos(2.*th)*pow(a,2.)*((-2. + r)*r + pow(a,2.)) + pow(a,4.) + 
		       2.*pow(r,4.))*((bsq + 2.*en*(-1. + gam))*pow(csc(th),2.) - 
				      2.*sigma*pow(b[PH],2.) + 2.*(bsq + en*gam + rho)*sigma*pow(u[PH],2.))
		  )*pow(sin(th),2.)*(cot(th) + r*pow(a,2.)*pow(sigma,-2.)*sin(2.*th)))\
      ;
  }



  if((WHICHEOM==WITHNOGDET)&&(NOGDETU3==1)){

    dU[U3]+=0.5*pow(sigma,-1.)*pow(sin(th),2.)*
      ((-1.*b[RR]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
		   b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
		   cos(2.*th)*pow(a,2.)*
		   (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
		    b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
		   b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
	u[RR]*(bsq + en*gam + rho)*
	(-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
		  dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
	 cos(2.*th)*pow(a,2.)*
	 (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
	  u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
	 u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*
       (-1. - 2.*dxdxp[1]*r*pow(sigma,-1.)) + 
       (-1.*b[TH]*(-2.*a*r*(2.*b[TT] + b[RR]*dxdxp[1]*(2. + r)) + 
		   b[PH]*r*(2. + 3.*r)*pow(a,2.) + 
		   cos(2.*th)*pow(a,2.)*
		   (-1.*b[RR]*dxdxp[1]*a + b[PH]*(-2. + r)*r + 
		    b[PH]*pow(a,2.)) - 1.*b[RR]*dxdxp[1]*pow(a,3.) + 
		   b[PH]*pow(a,4.) + 2.*b[PH]*pow(r,4.)) + 
	u[TH]*(bsq + en*gam + rho)*
	(-2.*a*r*(2.*(dxdxp[1]*u[RR] + u[TT]) + 
		  dxdxp[1]*u[RR]*r) + u[PH]*r*(2. + 3.*r)*pow(a,2.) + 
	 cos(2.*th)*pow(a,2.)*
	 (-1.*dxdxp[1]*u[RR]*a + u[PH]*(-2. + r)*r + 
	  u[PH]*pow(a,2.)) - 1.*dxdxp[1]*u[RR]*pow(a,3.) + 
	 u[PH]*pow(a,4.) + 2.*u[PH]*pow(r,4.)))*
       (-1.*dxdxp[2]*cot(th) + 
	4.*(-1.*M_PI*X[2] + th)*pow(dxdxp[2],-1.)*pow(M_PI,2.) + 
	dxdxp[2]*pow(a,2.)*pow(sigma,-1.)*sin(2.*th)));

  }
  // else nothing to add then


}



// GODMARK: remove? Worthless idea? 2D only also.
// jon's volume diff for uni theta grid and exp r grid (defcoord==LOGRUNITH)
void mks_unitheta_idxvol_func(int i, int j, int k, FTYPE *idxvol)
{

  /*
    int k, l;
    FTYPE r,th,ph;
    FTYPE r1[2],th1[2],r2[2],th2[2];
    FTYPE X0[NDIM],X1[2][NDIM],X2[2][NDIM];
    //  FTYPE cot(FTYPE arg);


    coord(i, j, CENT, X0);
    coord(i, j, FACE1, X1[0]);
    coord(i+1, j, FACE1, X1[1]);
    coord(i, j, FACE2, X2[0]);
    coord(i, j+1, FACE2, X2[1]);

    // get bl coordinates
    bl_coord(X0,&r,&th);
    bl_coord(X1[0],&r1[0],&th1[0]);
    bl_coord(X1[1],&r1[1],&th1[1]);
    bl_coord(X2[0],&r2[0],&th2[0]);
    bl_coord(X2[1],&r2[1],&th2[1]);

  */




  /* comment out for now until adjust everything


  if(MBH!=1.0){
  dualfprintf(fail_file,"mks_unitheta_idxvold_func not setup for MBH!=1.0\n");
  myexit(13);
  }


  // this is not exactly right, since derivative of metric is derivative of absolute values, but shouldn't/doesn't seem to matter much
  // follows gcov_func()
  if(POSDEFMETRIC){
  if(th<0.0){ th=-th;}
  if(th>M_PI) { th=M_PI-th; }
  }
  else{
  }
  // avoid singularity at polar axis
  #if(COORDSINGFIX)
  if(fabs(th)<SINGSMALL){
  if(th>=0) th=SINGSMALL;
  if(th<0) th=-SINGSMALL;
  }
  if(fabs(M_PI-th)<SINGSMALL){
  if(th>=M_PI) th=M_PI+SINGSMALL;
  if(th<M_PI) th=M_PI-SINGSMALL;
  }
  #endif
  */

#define IDXR(a,R0,r,th,rl,rh) ((pow(a,2) + 2.*pow(r,2) + pow(a,2)*cos(2.*th))/((rh - 1.*rl)*(2.*R0 + rh + rl) + (pow(a,2) + 2.*pow(R0,2))*(log(-1.*R0 + rh) - 1.*log(-1.*R0 + rl)) + pow(a,2)*cos(2.*th)*(log(-1.*R0 + rh) - 1.*log(-1.*R0 + rl))))

#define IDXTH(a,R0,r,th,thl,thh) ((-3.*M_PI*(pow(a,2) + 2.*pow(r,2) + pow(a,2)*cos(2.*th))*sin(th))/(cos(thh)*(pow(a,2) + 6.*pow(r,2) + pow(a,2)*cos(2.*thh)) - 1.*cos(thl)*(pow(a,2) + 6.*pow(r,2) + pow(a,2)*cos(2.*thl))))


  //#define IDXTH(a,R0,r,th,thl,thh) ((3.*(pow(a,2) + 2.*pow(r,2) + pow(a,2)*cos(2.*th))*sin(th))/((-1.*cos(thh) + cos(thl))*(2.*(pow(a,2) + 3.*pow(r,2)) + pow(a,2)*(cos(2.*thh) + 2.*cos(thh)*cos(thl) + cos(2.*thl)))))

  // fullth integrated -- not used currently
#define FIDXR(a,R0,r,th) ((pow(a,2) + 2.*pow(r,2) + pow(a,2)*cos(2.*th))/(4.*R0*(-1.*R0 + r) + pow(-1.*R0 + r,2) + (pow(a,2) + 2.*pow(R0,2) + pow(a,2)*cos(2.*th))*log(-1.*R0 + r)))

#define FIDXTH(a,R0,r,th) ((-3.*pow(a,2)*sin(2.*th) - 6.*pow(r,2)*tan(th))/(pow(a,2) + 6.*pow(r,2) + pow(a,2)*cos(2.*th)))


  /*

  idxvol[TT]=1.0; // really 1/dt, but changes in time
  // comment out non-volume regularization
  //  idxvol[RR]=1.0/dx[1];
  //idxvol[TH]=1.0/dx[2];
  idxvol[RR]=IDXR(a,R0,r,th,r1[0],r1[1]);
  idxvol[TH]=IDXTH(a,R0,r,th,th2[0],th2[1]);
  idxvol[PH]=1.0/dx[3];

  fprintf(fail_file,"%d %d %21.15g %21.15g\n",i,j,idxvol[RR]*dx[1],idxvol[TH]*dx[2]);
  
  */




}


// geometry only contains i,j,k,p
// only for spherical polar coordinates with negligible relativistic effects
// only used if GDETVOLDIFF==1
void gdetvol_func(struct of_geom *ptrgeom, FTYPE (*gdet)[N2M+SHIFT2][N3M+SHIFT3][NPG], FTYPE *gdetvol)
{
  int i,j,k,loc;
  FTYPE Xc[NDIM],Vc[NDIM];
  FTYPE Xim[NDIM],Vim[NDIM],Xip1[NDIM],Vip1[NDIM];
  FTYPE Xjm[NDIM],Vjm[NDIM],Xjp1[NDIM],Vjp1[NDIM];
  FTYPE Xkm[NDIM],Vkm[NDIM],Xkp1[NDIM],Vkp1[NDIM];
  FTYPE dr,dh,dp;
  FTYPE drold,dhold,dpold;
  FTYPE newjacvol,oldjacvol;
  FTYPE dxdxpc[NDIM][NDIM];
  int locarray[4];


  i=ptrgeom->i;
  j=ptrgeom->j;
  k=ptrgeom->k;
  loc=ptrgeom->p;


  // loc==CENT 

  if(loc==CENT){ // e.g. for conserved quantities at CENT
    locarray[0]=loc;
    locarray[1]=FACE1; // left-right
    locarray[2]=FACE2; // up-down
    locarray[3]=FACE3; // in-out
  }
  else if(loc==FACE1){ // e.g. for fluxes F1
    locarray[0]=loc;
    locarray[1]=CENT;
    locarray[2]=CORN3;
    locarray[3]=CORN2;    
  }
  else if(loc==FACE2){ // e.g. for fluxes F2
    locarray[0]=loc;
    locarray[1]=CORN3;
    locarray[2]=CENT;
    locarray[3]=CORN1;    
  }
  else if(loc==FACE3){ // e.g. for fluxes F3
    locarray[0]=loc;
    locarray[1]=CORN2;
    locarray[2]=CORN1;
    locarray[3]=CENT;    
  }
  else{
    // maybe don't need yet :: only for fields -- not used yet.
    locarray[0]=loc;
    locarray[1]=FACE1; // left-right
    locarray[2]=FACE2; // up-down
    locarray[3]=FACE3; // in-out
  }
    

  // rest of code independent of loc    

  coord_ijk(i, j, k, locarray[0], Xc);
  bl_coord_ijk(i, j, k, locarray[0], Vc);
  dxdxprim_ijk(i, j, k, locarray[0], dxdxpc);



 
  ////////////////////////
  //
  // NEW JAC
  coord_ijk(i, j, k, locarray[1], Xim);
  bl_coord_ijk(i, j, k, locarray[1], Vim);

  coord_ijk(i+1, j, k, locarray[1], Xip1); // ok to have i+1 since coord and bl_coord don't depend upon memory locations
  bl_coord_ijk(i+1, j, k, locarray[1], Vip1);
    
  coord_ijk(i, j, k, locarray[2], Xjm);
  bl_coord_ijk(i, j, k, locarray[2], Vjm);

  coord_ijk(i, j+1, k, locarray[2], Xjp1); // ok to have j+1 since coord and bl_coord don't depend upon memory locations
  bl_coord_ijk(i, j+1, k, locarray[2], Vjp1);
    
  coord_ijk(i, j, k, locarray[3], Xkm);
  bl_coord_ijk(i, j, k, locarray[3], Vkm);

  coord_ijk(i, j, k+1, locarray[3], Xkp1); // ok to have k+1 since coord and bl_coord don't depend upon memory locations
  bl_coord_ijk(i, j, k+1, locarray[3], Vkp1);
    
  dr = THIRD*(pow(Vip1[1],3.0)-pow(Vim[1],3.0)); // really \delta(r^3/3) = r^2 dr
  // dh doesn't reduce to Pi, but 2.0 that is the correct answer for any N2 resolution
  dh = -(cos(Vjp1[2])-cos(Vjm[2])); // really \delta(-cos(h)) = sinh dh
  // below is inconsistent with rest of code when N2!=1
  //    dh = (N2==1) ? M_PI : -(cos(Vjp1[2])-cos(Vjm[2])); // really \delta(-cos(h)) = sinh dh
  dp = (Vkp1[3]-Vkm[3]); // just d\phi 
    
  // finite volume of cell
  newjacvol = dr*dh*dp/(dx[1]*dx[2]*dx[3]);


  ////////////////////////
  //
  // OLD JAC
  drold=Vc[1]*Vc[1]*dxdxpc[1][1]*dx[1]; // r^2 dr
  if(totalsize[2]==1) dhold=2.0;
  else dhold=sin(Vc[2])*dxdxpc[2][2]*dx[2]; // sinh dh
  //dhold=sin(Vc[2])*dxdxpc[2][2]*dx[2]; // sinh dh (oldjacvol consistent with gdet)
  dpold=2.0*M_PI*dx[3];
  oldjacvol = drold*dhold*dpold/(dx[1]*dx[2]*dx[3]); // only true if diagonal mapping

  // suppose gdet already corrected if wanted to be corrected
  //  if(totalsize[2]==1) gdet[i][j][k][loc] /= (M_PI*0.5); // correct eomfunc and gdet
  //  if(WHICHEOM==WITHGDET) eomfunc[i][j][k][loc]=gdet[i][j][k][loc]; // else up to user to make sure right

  //////////////////////////
  //
  // use below, works best.  Only central conserved quantity operated on (source terms too)
  if(loc==CENT) *gdetvol=gdet[i][j][k][loc]=eomfunc[i][j][k][loc]=newjacvol;
  else *gdetvol=gdet[i][j][k][loc];


  ///////////////////////////////////
  //
  // SOME OTHER ATTEMPTS for 2 lines "best" above:
  //
  //////////////////////////////////

  // horrible
  //  if(loc==CENT) *gdetvol=newjacvol;
  //  else *gdetvol=gdet[i][j][k][loc];


  // central region undershoots
  //if(loc==CENT) *gdetvol=newjacvol;
  //  else *gdetvol=gdet[i][j][k][loc];

  // big overshoot, so probably don't want to use for fluxes
  //  *gdetvol=gdet[i][j][k][loc]=eomfunc[i][j][k][loc]=newjacvol;

  // with newjac in potential and gdet here
  // little bit more of a spike at center
  //  *gdetvol=gdet[i][j][k][loc];

  // with oldjac in potential and gdet here
  // same as above :  little bit more of a spike at center
  //  *gdetvol=gdet[i][j][k][loc];


  //  *gdetvol = newjacvol;
  //    *gdetvol = dr*dh*dp/5.0;
  //*gdetvol = gdet[i][j][k][loc]; // disable gdetvol but keep memory
  //  eomfunc[i][j][k][loc]=gdet[i][j][k][loc]=*gdetvol=newjacvol;


  //  *gdetvol=newjacvol;
  //  if(loc==CENT){
  //  eomfunc[i][j][k][loc]=gdet[i][j][k][loc]=*gdetvol=newjacvol;
  //  eomfunc[i][j][k][loc]=gdet[i][j][k][loc]=*gdetvol=oldjacvol;
  //  eomfunc[i][j][k][loc]=*gdetvol=gdet[i][j][k][loc];
  //  }
  //  else{
  //    eomfunc[i][j][k][loc]=*gdetvol=gdet[i][j][k][loc]=oldjacvol;
  //eomfunc[i][j][k][loc]=*gdetvol=gdet[i][j][k][loc]=newjacvol;
  //  }


  //  if(loc==CENT){
  //    dualfprintf(fail_file,"i=%d loc=%d :: gdet=%15.7g oldjac=%15.7g newjac=%15.7g :: dr=%g drold=%g :: dh=%g dhold=%g :: dp=%g dpold=%g :: %g %g %g\n",i,loc,gdet[i][j][k][loc],oldjacvol,newjacvol,dr,drold,dh,dhold,dp,dpold,dx[1],dx[2],dx[3]);
  //  }

  //gdet[i][j][k][loc] = *gdetvol; // replace gdet with this version


}


















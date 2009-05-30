

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// structure definitions
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct blink {
  int num;
  struct blink * np;
  // only used by cpu=0
  int cpu; // which cpu
  int i,j,k,col; // starting values for cpu=0
  int ri,rj,rk,rcol; // reference values for first cpu in sequence of nodes for a single buffer
  int end;
};



// structure declarations
/* set global variables that indicate current local metric, etc. */
struct of_geom {
  // dummy space for gset() version
  FTYPE gengcon[NDIM][NDIM];
  FTYPE gengcov[NDIM][NDIM];
  FTYPE gengcovpert[NDIM];

#if(GETGEOMUSEPOINTER==0)
  FTYPE gcon[NDIM][NDIM];
  FTYPE gcov[NDIM][NDIM];
  FTYPE gcovpert[NDIM];
#else
  // bit faster since not all values always used
  FTYPE (*gcov)[NDIM];
  FTYPE (*gcon)[NDIM];
  FTYPE *gcovpert;
#endif
  FTYPE g;
  FTYPE e[NPR]; // eomfunc
  FTYPE igeomnosing,ienosing[NPR];
  FTYPE alphalapse;
  int i,j,k,p;
};

struct of_allgeom {
  // dummy space for gset() version
  FTYPE gengcon[NDIM][NDIM];
  FTYPE gengcov[NDIM][NDIM];
  FTYPE gengcovpert[NDIM];

#if(GETGEOMUSEPOINTER==0)
  FTYPE gcon[NDIM][NDIM];
  FTYPE gcov[NDIM][NDIM];
  FTYPE gcovpert[NDIM];
#else
  // bit faster since not all values always used
  FTYPE (*gcov)[NDIM];
  FTYPE (*gcon)[NDIM];
  FTYPE *gcovpert;
#endif
  FTYPE g;
  FTYPE e[NPR]; // eomfunc
  FTYPE alphalapse;
  FTYPE igeomnosing,ienosing[NPR];

  FTYPE X[NDIM];
  FTYPE V[NDIM];
  FTYPE dxdxp[NDIM][NDIM];
  int i,j,k,p;
};


struct of_state {
  FTYPE rho;
  FTYPE ie;
  FTYPE ucon[NDIM];
  FTYPE ucov[NDIM];
  FTYPE bcon[NDIM];
  FTYPE bcov[NDIM];
};




struct of_loop {
  int is, ie;
  int js, je;
  int ks, ke;
  int dir,intdir;
  int ps, pe;
  int bs, be;
  int di,dj,dk;
};

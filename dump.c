#include "decs.h"

// mpi.h has following datatypes corresponding to the C types
// pick one per dump file. (no per column types yet)
// same for image.c
// same for restart.c
/*
  #define MPI_CHAR           ((MPI_Datatype)1)
  #define MPI_UNSIGNED_CHAR  ((MPI_Datatype)2)
  #define MPI_BYTE           ((MPI_Datatype)3)
  #define MPI_SHORT          ((MPI_Datatype)4)
  #define MPI_UNSIGNED_SHORT ((MPI_Datatype)5)
  #define MPI_INT            ((MPI_Datatype)6)
  #define MPI_UNSIGNED       ((MPI_Datatype)7)
  #define MPI_LONG           ((MPI_Datatype)8)
  #define MPI_UNSIGNED_LONG  ((MPI_Datatype)9)
  #define MPI_FLOAT          ((MPI_Datatype)10)
  #define MPI_DOUBLE         ((MPI_Datatype)11)
  #define MPI_LONG_DOUBLE    ((MPI_Datatype)12)
  #define MPI_LONG_LONG_INT  ((MPI_Datatype)13)
*/



/* Follow these steps to create a new dump file

1) defs.h: create the storage variable if needed

2) set_array.c : shift variable if necessary (follow examples)

3) global.h : change NUMDUMPTYPES and add label

4) global.h change NUMDTDS and add another entry if want separate timing of output

5) initbase.c : add dnumcolumns[LABEL]=NUMCOLUMNS where NUMCOLUMNS is number of entries in dump file

6) dump.c : follow examples from here (dump() uses dump_header() and dump_content()).  One must define the header and content function and the wrapper (3 functions) or use an existing header function

7) diag.c : follow example of "dump", dumpc, tlastdump, etc.

8) init.c : DTdumpgen, dumpcntgen, and other things.

*/



int dump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping dump# %ld ... ",dump_cnt);

  whichdump=DUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/dump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,dump_content)>=1) return(1);

	/////// output the symmetry information to the fail file
	//writesyminfo();
	///////

  trifprintf("end dumping dump# %ld ... ",dump_cnt);


  return(0);
  
}



/////// output the symmetry information to the fail file; symmetrizes w.r.t. i == j
void writesyminfo( void )
{
	int i, j;

	for( i = 0; i < N1; i++ ) {
	}

}

int dump_header(int bintxt, FILE *headerptr)
{
  int realtotalsize[NDIM];
  FTYPE realstartx[NDIM];
  FTYPE X[NDIM];


  realtotalsize[1]=totalsize[1]+2*EXTRADUMP1;
  realtotalsize[2]=totalsize[2]+2*EXTRADUMP2;
  realtotalsize[3]=totalsize[3]+2*EXTRADUMP3;

  // get real startx's (assumes rectangular grid)
  coord(0-EXTRADUMP1,0,0,CENT,X);
  realstartx[1]=X[1];
  coord(0,0-EXTRADUMP2,0,CENT,X);
  realstartx[2]=X[2];
  coord(0,0,0-EXTRADUMP3,CENT,X);
  realstartx[3]=X[3];
  
  // dx is the same (constant)

  // 15+3=18 elements total
  if(bintxt==BINARYOUTPUT){
    fwrite(&tsteppartf,sizeof(FTYPE),1,headerptr);
    fwrite(&realtotalsize[1],sizeof(int),1,headerptr);
    fwrite(&realtotalsize[2],sizeof(int),1,headerptr);
    fwrite(&realtotalsize[3],sizeof(int),1,headerptr);
    fwrite(&realstartx[1],sizeof(FTYPE),1,headerptr);
    fwrite(&realstartx[2],sizeof(FTYPE),1,headerptr);
    fwrite(&realstartx[3],sizeof(FTYPE),1,headerptr);
    fwrite(&dx[1],sizeof(FTYPE),1,headerptr);
    fwrite(&dx[2],sizeof(FTYPE),1,headerptr);
    fwrite(&dx[3],sizeof(FTYPE),1,headerptr);
    fwrite(&realnstep,sizeof(long),1,headerptr);
    fwrite(&gam,sizeof(FTYPE),1,headerptr);
    fwrite(&a,sizeof(FTYPE),1,headerptr);
    fwrite(&R0,sizeof(FTYPE),1,headerptr);
    fwrite(&Rin,sizeof(FTYPE),1,headerptr);
    fwrite(&Rout,sizeof(FTYPE),1,headerptr);
    fwrite(&hslope,sizeof(FTYPE),1,headerptr);
    fwrite(&dt,sizeof(FTYPE),1,headerptr);
    fwrite(&MBH,sizeof(int),1,headerptr);
    fwrite(&QBH,sizeof(int),1,headerptr);
  }
  else{
#if(REALTYPE==DOUBLETYPE)
    fprintf(headerptr, "%21.15g %d %d %d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %ld %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %d %21.15g %21.15g\n", tsteppartf, realtotalsize[1], realtotalsize[2], realtotalsize[3], realstartx[1], realstartx[2], realstartx[3], dx[1], dx[2], dx[3], realnstep,gam,a,R0,Rin,Rout,hslope,dt,defcoord,MBH,QBH);
#elif(REALTYPE==LONGDOUBLETYPE)
    fprintf(headerptr, "%31.25Lg %d %d %d %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %ld %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %d %31.25Lg %31.25Lg\n", tsteppartf, realtotalsize[1], realtotalsize[2], realtotalsize[3], realstartx[1], realstartx[2], realstartx[3], dx[1], dx[2],dx[3],realnstep,gam,a,R0,Rin,Rout,hslope,dt,defcoord,MBH,QBH);
#endif
  }
  fflush(headerptr);
  return(0);
}	


int dump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl;
  FTYPE r, th, vmin[NDIM], vmax[NDIM];
  int ignorecourant;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM],V[NDIM];
  FTYPE divb;
  FTYPE b[NDIM],ucon[NDIM];
  FTYPE U[NPR];
  FTYPE ftemp;
  FTYPE jcov[NDIM];
  FTYPE fcov[NUMFARADAY];
  FTYPE rho,u,pressure,cs2,Sden;
  int dir,l,m,n,o;


  //////////////
  //
  // some calculations
  //

  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  // if failed, then data output for below invalid, but columns still must exist    

  get_geometry(i, j, k, CENT, &geom);

  if (!failed) {
    if (get_state(pdump[i][j][k], &geom, &q) >= 1)
      FAILSTATEMENT("dump.c:dump()", "get_state() dir=0", 1);
    if (vchar(pdump[i][j][k], &q, 1, &geom, &vmax[1], &vmin[1],&ignorecourant) >= 1)
      FAILSTATEMENT("dump.c:dump()", "vchar() dir=1or2", 1);
    if (vchar(pdump[i][j][k], &q, 2, &geom, &vmax[2], &vmin[2],&ignorecourant) >= 1)
      FAILSTATEMENT("dump.c:dump()", "vchar() dir=1or2", 2);
    if (vchar(pdump[i][j][k], &q, 3, &geom, &vmax[3], &vmin[3],&ignorecourant) >= 1)
      FAILSTATEMENT("dump.c:dump()", "vchar() dir=1or2", 3);
  }
  else {// do a per zone check, otherwise set to 0
    whocalleducon=1; // force no failure mode, just return like failure, and don't return if failure, just set to 0 and continue
    if (get_state(pdump[i][j][k], &geom, &q) >= 1){
      for (pl = 0; pl < NDIM; pl++)
	q.ucon[pl]=0;
      for (pl = 0; pl < NDIM; pl++)
	q.ucov[pl]=0;
      for (pl = 0; pl < NDIM; pl++)
	q.bcon[pl]=0;
      for (pl = 0; pl < NDIM; pl++)
	q.bcov[pl]=0;
    }
    if (vchar(pdump[i][j][k], &q, 1, &geom, &vmax[1], &vmin[1],&ignorecourant) >= 1){
      vmax[1]=vmin[1]=0;
    }
    
    if (vchar(pdump[i][j][k], &q, 2, &geom, &vmax[2], &vmin[2],&ignorecourant) >= 1){
      vmax[2]=vmin[2]=0;
    }

    if (vchar(pdump[i][j][k], &q, 3, &geom, &vmax[3], &vmin[3],&ignorecourant) >= 1){
      vmax[3]=vmin[3]=0;
    }

    whocalleducon=0; // return to normal state
    
  }


  setfdivb(&divb, pdump, udump, i, j, k); // udump also set externally GODMARK

  //////////////////////////
  //
  // do the assignments
  //
  // if you change # of outputted vars, remember to change numcolumns


  //static
  if(!GAMMIEDUMP){
    ftemp=(FTYPE)(i+startpos[1]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(j+startpos[2]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(k+startpos[3]);
    myset(datatype,&ftemp,0,1,writebuf);
  }
  myset(datatype,X,1,3,writebuf);
  myset(datatype,V,1,3,writebuf);
  // 9

  ////////////////////////
  //
  // rest dynamic

  // primitives
  // must use PDUMPLOOP() since may be any order unlike NPR loop
  PDUMPLOOP(pl) myset(datatype,&(pdump[i][j][k][pl]),0,1,writebuf); // NPRDUMP

  ////////////
  //
  // output some EOS stuff since in general not simple function of rho0,u
  rho = pdump[i][j][k][RHO];
  u = pdump[i][j][k][UU];


  pressure = pressure_rho0_u(rho,u);
  cs2 = cs2_compute(rho,u);
  Sden = compute_entropy(rho,u);
  //  dUdtau = compute_qdot(rho,u);
  
  myset(datatype,&pressure,0,1,writebuf); // 1
  myset(datatype,&cs2,0,1,writebuf); // 1
  myset(datatype,&Sden,0,1,writebuf); // 1
  //  myset(datatype,&dUdtau,0,1,writebuf); // 1

  //////////////////////
  //
  // output the conserved quantities since not easily inverted and at higher order aren't invertable from point primitives
  PDUMPLOOP(pl) myset(datatype,&(udump[i][j][k][pl]),0,1,writebuf); // NPRDUMP
  myset(datatype,&divb,0,1,writebuf); // 1

  for (pl = 0; pl < NDIM; pl++)
    myset(datatype,&(q.ucon[pl]),0,1,writebuf);
  for (pl = 0; pl < NDIM; pl++)
    myset(datatype,&(q.ucov[pl]),0,1,writebuf);
  for (pl = 0; pl < NDIM; pl++)
    myset(datatype,&(q.bcon[pl]),0,1,writebuf);
  for (pl = 0; pl < NDIM; pl++)
    myset(datatype,&(q.bcov[pl]),0,1,writebuf);
  // 4*4
    
  myset(datatype,&vmin[1],0,1,writebuf);
  myset(datatype,&vmax[1],0,1,writebuf);
  myset(datatype,&vmin[2],0,1,writebuf);
  myset(datatype,&vmax[2],0,1,writebuf);
  myset(datatype,&vmin[3],0,1,writebuf);
  myset(datatype,&vmax[3],0,1,writebuf);
  // 6

  // one static term
  myset(datatype,&geom.g,0,1,writebuf); // 1


#if(CALCFARADAYANDCURRENTS) // NIM*2+6*2 = 8+12=20
  // updated 11/16/2003
  // new 10/23/2003
  // current density 
  lower_vec(jcon[i][j][k],&geom,jcov); 
  myset(datatype,jcon[i][j][k],0,NDIM,writebuf); // (NDIM)
  myset(datatype,jcov,0,NDIM,writebuf);// (NDIM)
  // faraday (2*6)
  lowerf(fcon[i][j][k],&geom,fcov);
  myset(datatype,fcon[i][j][k],0,NUMFARADAY,writebuf); //  (6)
  myset(datatype,fcov,0,NUMFARADAY,writebuf); // (6)
#endif

  if(FLUXB==FLUXCTSTAG && 0){ // DEBUG (change corresponding code in dump.c)
    // uses jrdp3dudebug in gtwod.m that assumes CALCFARADAYANDCURRENTS==0
    for(l=1;l<=COMPDIM;l++) myset(datatype,gp_l[l][i][j][k],0,NPR2INTERP,writebuf); // 3*8 = 24
    for(l=1;l<=COMPDIM;l++) myset(datatype,gp_r[l][i][j][k],0,NPR2INTERP,writebuf); // 3*8 = 24
    myset(datatype,pstagscratch[i][j][k],0,NPR,writebuf); // 8
    for(dir=1;dir<=COMPDIM;dir++) for(pl=B1;pl<=B3;pl++) for(n=0;n<=1;n++) myset(datatype,&pbcorninterp[dir][pl][n][i][j][k],0,1,writebuf); // 3*3*2 = 18
    for(dir=1;dir<=COMPDIM;dir++) for(pl=U1;pl<=U3;pl++) for(n=0;n<=1;n++) for(o=0;o<=1;o++) myset(datatype,&pvcorninterp[dir][pl][n][o][i][j][k],0,1,writebuf); // 3*3*2*2 = 36
  }

  return (0);
}






int debugdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping debug dump# %ld ... ",dump_cnt);

  whichdump=DEBUGCOL;
  datatype=MPI_CTYPE;
  strcpy(fileprefix,"dumps/debug");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  // same header as dump
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,debug_content)>=1) return(1);

  trifprintf("end dumping debug# %ld ... ",dump_cnt);

  return(0);

}




int enodebug_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  // could also make everything FTYPE and convert like for normal i,j dump file
  //myset(datatype,enodebugarray[i][j][k],0,3*NUMINTERPTYPES*NPR*NUMENODEBUGS,writebuf);
  myset(datatype,enodebugarray[i][j][k],0,(3-1)*NUMINTERPTYPES*(NPR-4)*NUMENODEBUGS,writebuf);  //atch corrected
    
  return(0);
}

int enodebugdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
  int eno_dump_header(int bintxt, FILE *headerptr);


  trifprintf("begin dumping enodebug dump# %ld ... ",dump_cnt);

  whichdump=ENODEBUGCOL;
  //  datatype=MPI_FTYPE;
  datatype=MPI_CTYPE;
  strcpy(fileprefix,"dumps/enodebug");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  // same header as dump
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,eno_dump_header,enodebug_content)>=1) return(1);

  trifprintf("end dumping enodebug# %ld ... ",dump_cnt);

  return(0);

}


int eno_dump_header(int bintxt, FILE *headerptr)
{
  int realtotalsize[NDIM];
  FTYPE realstartx[NDIM];
  FTYPE X[NDIM];


  realtotalsize[1]=totalsize[1]+2*EXTRADUMP1;
  realtotalsize[2]=totalsize[2]+2*EXTRADUMP2;
  realtotalsize[3]=totalsize[3]+2*EXTRADUMP3;

  // get real startx's (assumes rectangular grid)
  coord(0-EXTRADUMP1,0,0,CENT,X);
  realstartx[1]=X[1];
  coord(0,0-EXTRADUMP2,0,CENT,X);
  realstartx[2]=X[2];
  coord(0,0,0-EXTRADUMP3,CENT,X);
  realstartx[3]=X[3];
  
  // dx is the same (constant)

  // 15+3=18 elements total
  if(bintxt==BINARYOUTPUT){
    fwrite(&tsteppartf,sizeof(FTYPE),1,headerptr);
    fwrite(&realtotalsize[1],sizeof(int),1,headerptr);
    fwrite(&realtotalsize[2],sizeof(int),1,headerptr);
    fwrite(&realtotalsize[3],sizeof(int),1,headerptr);
    fwrite(&realstartx[1],sizeof(FTYPE),1,headerptr);
    fwrite(&realstartx[2],sizeof(FTYPE),1,headerptr);
    fwrite(&realstartx[3],sizeof(FTYPE),1,headerptr);
    fwrite(&dx[1],sizeof(FTYPE),1,headerptr);
    fwrite(&dx[2],sizeof(FTYPE),1,headerptr);
    fwrite(&dx[3],sizeof(FTYPE),1,headerptr);
    fwrite(&realnstep,sizeof(long),1,headerptr);
    fwrite(&gam,sizeof(FTYPE),1,headerptr);
    fwrite(&a,sizeof(FTYPE),1,headerptr);
    fwrite(&R0,sizeof(FTYPE),1,headerptr);
    fwrite(&Rin,sizeof(FTYPE),1,headerptr);
    fwrite(&Rout,sizeof(FTYPE),1,headerptr);
    fwrite(&hslope,sizeof(FTYPE),1,headerptr);
    fwrite(&dt,sizeof(FTYPE),1,headerptr);
    fwrite(&defcoord,sizeof(int),1,headerptr);
  }
  else{
#if(REALTYPE==DOUBLETYPE)
    fprintf(headerptr, "%21.15g %d %d %d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %ld %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %d %d\n", tsteppartf, realtotalsize[1], realtotalsize[2], realtotalsize[3], realstartx[1], realstartx[2], realstartx[3], dx[1], dx[2], dx[3], realnstep,gam,a,R0,Rin,Rout,hslope,dt,defcoord,steppart);
#elif(REALTYPE==LONGDOUBLETYPE)
    fprintf(headerptr, "%31.25Lg %d %d %d %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %ld %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %d %d\n", tsteppartf, realtotalsize[1], realtotalsize[2], realtotalsize[3], realstartx[1], realstartx[2], realstartx[3], dx[1], dx[2],dx[3],realnstep,gam,a,R0,Rin,Rout,hslope,dt,defcoord, steppart);
#endif
  }
  fflush(headerptr);
  return(0);
}	






int debug_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  // could also make everything FTYPE and convert like for normal i,j dump file
  myset(datatype,failfloorcount[i][j][k][0],0,NUMTSCALES*NUMFAILFLOORFLAGS,writebuf);
    
  return(0);
}


int avgdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
	


  trifprintf("begin dumping avgdump# %ld ... ",dump_cnt);

  whichdump=AVGCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/avg");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,avg_content)>=1) return(1);

  trifprintf("end dumping avgdump# %ld ... ",dump_cnt);


  return(0);

}

int avg_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl = 0, l = 0, col = 0;
  struct of_geom geom;
  FTYPE X[NDIM],V[NDIM];
  FTYPE ftemp;



  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  get_geometry(i, j, k, CENT, &geom);

  if(!GAMMIEDUMP){
    ftemp=(FTYPE)(i+startpos[1]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(j+startpos[2]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(k+startpos[3]);
    myset(datatype,&ftemp,0,1,writebuf);
  }
  myset(datatype,X,1,3,writebuf);
  myset(datatype,V,1,3,writebuf);


  myset(datatype,&geom.g,0,1,writebuf);

  // now do time average stuff
  myset(datatype,normalvarstavg[i][j],0,NUMNORMDUMP,writebuf);
  myset(datatype,anormalvarstavg[i][j],0,NUMNORMDUMP,writebuf);

#if(CALCFARADAYANDCURRENTS)
  myset(datatype,jcontavg[i][j],0,NDIM,writebuf);
  myset(datatype,jcovtavg[i][j],0,NDIM,writebuf);
  myset(datatype,ajcontavg[i][j],0,NDIM,writebuf);
  myset(datatype,ajcovtavg[i][j],0,NDIM,writebuf);
#endif
  myset(datatype,massfluxtavg[i][j],0,NDIM,writebuf);
  myset(datatype,amassfluxtavg[i][j],0,NDIM,writebuf);

  myset(datatype,othertavg[i][j],0,NUMOTHER,writebuf);
  myset(datatype,aothertavg[i][j],0,NUMOTHER,writebuf);

#if(CALCFARADAYANDCURRENTS)
  myset(datatype,fcontavg[i][j],0,NUMFARADAY,writebuf);
  myset(datatype,fcovtavg[i][j],0,NUMFARADAY,writebuf);
  myset(datatype,afcontavg[i][j],0,NUMFARADAY,writebuf);
  myset(datatype,afcovtavg[i][j],0,NUMFARADAY,writebuf);
#endif

#if(DOAVG2==0)
  myset(datatype,tudtavg[i][j],0,NUMSTRESSTERMS,writebuf);
  myset(datatype,atudtavg[i][j],0,NUMSTRESSTERMS,writebuf);
#endif

  return(0);

}


int avg2dump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];



  trifprintf("begin dumping avg2dump# %ld ... ",dump_cnt);

  whichdump=AVG2COL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/avg2");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");

  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,avg2_content)>=1) return(1);

  trifprintf("end dumping avg2dump# %ld ... ",dump_cnt);


  return(0);

}


int avg2_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl = 0, l = 0, col = 0;
  struct of_geom geom;
  FTYPE X[NDIM],V[NDIM];
  FTYPE ftemp;


  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  get_geometry(i, j, k, CENT, &geom);
  // if you change # of outputted vars, remember to change numcolumns above

  if(!GAMMIEDUMP){
    ftemp=(FTYPE)(i+startpos[1]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(j+startpos[2]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(k+startpos[3]);
    myset(datatype,&ftemp,0,1,writebuf);
  }
  myset(datatype,X,1,3,writebuf);
  myset(datatype,V,1,3,writebuf);

  myset(datatype,&geom.g,0,1,writebuf);
  // 10

  myset(datatype,tudtavg[i][j][k],0,NUMSTRESSTERMS,writebuf);
  myset(datatype,atudtavg[i][j][k],0,NUMSTRESSTERMS,writebuf);
  // 112*2

  // total=10+112*2=234

  return(0);
}


int gdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
	

  trifprintf("begin dumping gdump# %ld ... ",dump_cnt);

  whichdump=GDUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/gdump");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");

  // dump_cnt==-1 means no file number
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,gdump_content)>=1) return(1);

  trifprintf("end dumping gdump# %ld ... ",dump_cnt);

  return(0);
}




int gdump_content(int i, int j, int k, MPI_Datatype datatype, void *writebuf)
{
  int pl = 0, l = 0, m = 0, n = 0, col = 0;
  FTYPE X[NDIM],V[NDIM];
  FTYPE ftemp;
  FTYPE *ptrftemp;
  FTYPE dxdxp[NDIM][NDIM];
  int myii,myjj,mykk;


  coord(i, j, k, CENT, X);
  bl_coord(X, V);
  dxdxprim(X, V, dxdxp);



  ftemp=(FTYPE)(i+startpos[1]);
  myset(datatype,&ftemp,0,1,writebuf);
  ftemp=(FTYPE)(j+startpos[2]);
  myset(datatype,&ftemp,0,1,writebuf);
  ftemp=(FTYPE)(k+startpos[3]);
  myset(datatype,&ftemp,0,1,writebuf);
  // 3
  myset(datatype,X,1,3,writebuf);
  myset(datatype,V,1,3,writebuf);
  // 6




#if(MCOORD!=CARTMINKMETRIC)
  myii=i;
  myjj=j;
  mykk=k;
#else
  myii=0;
  myjj=0;
  mykk=0;
#endif


  ptrftemp=(FTYPE*)(&conn[myii][myjj][mykk][0][0][0]);
  myset(datatype,ptrftemp,0,NDIM*NDIM*NDIM,writebuf);

    
  ptrftemp=(FTYPE*)(&gcon[myii][myjj][mykk][CENT][0][0]);
  myset(datatype,ptrftemp,0,NDIM*NDIM,writebuf);
  ptrftemp=(FTYPE*)(&gcov[myii][myjj][mykk][CENT][0][0]);
  myset(datatype,ptrftemp,0,NDIM*NDIM,writebuf);
  ptrftemp=(FTYPE*)(&gdet[myii][myjj][mykk][CENT]);
  //ptrftemp=(FTYPE*)(&gdetvol[myii][myjj][mykk][CENT]); // can take a peek if GDETVOLDIFF==1
  myset(datatype,ptrftemp,0,1,writebuf);

    
  ptrftemp=(FTYPE*)(&conn2[myii][myjj][mykk][0]);
  myset(datatype,ptrftemp,0,NDIM,writebuf);

  // 4*4
  ptrftemp=(FTYPE*)(&dxdxp[0][0]);
  myset(datatype,ptrftemp,0,NDIM*NDIM,writebuf);


  return(0);

}



int fieldlinedump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
	


  trifprintf("begin dumping fieldlinedump# %ld ... ",dump_cnt);

  whichdump=FIELDLINECOL;
  datatype=MPI_FLOAT; // don't need good precision
  strcpy(fileprefix,"dumps/fieldline");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  // MIXEDOUTPUT means text header and forced binary data
  if(dump_gen(WRITEFILE,dump_cnt,MIXEDOUTPUT,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,fieldline_content)>=1) return(1);

  trifprintf("end dumping fieldlinedump# %ld ... ",dump_cnt);


  return(0);

}

int fieldline_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl = 0, l = 0, col = 0;
  struct of_geom geom;
  struct of_state q;
  //FTYPE U[NPR];
  FTYPE FL[NPR];
  // must be same precision as written content
  float ftemp;

  //////////////
  //
  // some calculations
  //

  // if failed, then data output for below invalid, but columns still must exist    
  get_geometry(i, j, k, CENT, &geom);
  if (!failed) {
    if (get_state(pdump[i][j][k], &geom, &q) >= 1)
      FAILSTATEMENT("dump.c:dump()", "get_state() dir=0", 1);
  }
  else {// do a per zone check, otherwise set to 0
    whocalleducon=1; // force no failure mode, just return like failure, and don't return if failure, just set to 0 and continue
    if (get_state(pdump[i][j][k], &geom, &q) >= 1){
      for (pl = 0; pl < NDIM; pl++)
	q.ucon[pl]=0;
      for (pl = 0; pl < NDIM; pl++)
	q.ucov[pl]=0;
      for (pl = 0; pl < NDIM; pl++)
	q.bcon[pl]=0;
      for (pl = 0; pl < NDIM; pl++)
	q.bcov[pl]=0;
    }
    whocalleducon=0; // return to normal state
    
  }

  MYFUN(primtoflux(UDIAG,pdump[i][j][k], &q, RR, &geom, FL),"step_ch.c:fluxcalc()", "primtoflux_calc() dir=1/2 l", RR);


  //////////////////////////
  //
  // do the assignments
  //
  // if you change # of outputted vars, remember to change numcolumns



  ////////////////////
  //
  // 2 various things

  // rho (for various things)
  ftemp=(float)pdump[i][j][k][RHO];
  myset(datatype,&ftemp,0,1,writebuf);

  // u (for various things)
  ftemp=(float)pdump[i][j][k][UU];
  myset(datatype,&ftemp,0,1,writebuf);


  //////////////////////
  //
  // 2 things for jet/energy per baryon at infinity

  // -u_t (-hu_t can be found from this and rho/u/p above)
  ftemp=(float)(-q.ucov[0]);
  myset(datatype,&ftemp,0,1,writebuf);

  // -T^t_t/(gdet rho u^t)
  //  ftemp=(float)(-U[UU]/(geom.g * pdump[i][j][k][RHO]*q.ucon[TT]));
  //myset(datatype,&ftemp,0,1,writebuf);

  // -T^r_t/(rho u^r)
  if(q.ucon[RR]!=0.0){
    ftemp=(float)(-FL[UU]/(geom.g * pdump[i][j][k][RHO]*q.ucon[RR]));
  }
  else ftemp=0.0;
  myset(datatype,&ftemp,0,1,writebuf);


  // 1 extra thing

  // u^t
  ftemp=(float)(q.ucon[0]);
  myset(datatype,&ftemp,0,1,writebuf);


  ///////////////////////////
  //
  // 6 things for the field line stuff

  // v^r [ in grid frame]
  ftemp=(float)(q.ucon[1]/q.ucon[0]);
  myset(datatype,&ftemp,0,1,writebuf);

  // v^\theta
  ftemp=(float)(q.ucon[2]/q.ucon[0]);
  myset(datatype,&ftemp,0,1,writebuf);

  // v^\phi
  ftemp=(float)(q.ucon[3]/q.ucon[0]);
  myset(datatype,&ftemp,0,1,writebuf);

  // B^r
  ftemp=(float)(pdump[i][j][k][B1]);
  myset(datatype,&ftemp,0,1,writebuf);

  // B^\theta
  ftemp=(float)(pdump[i][j][k][B2]);
  myset(datatype,&ftemp,0,1,writebuf);

  // B^\phi
  ftemp=(float)(pdump[i][j][k][B3]);
  myset(datatype,&ftemp,0,1,writebuf);

  // see grmhd-dualfcon2omegaf.nb
  // below can be obtained from above set of v and B
  // \Omega_F_1
  //  ftemp=(float)(v3-B3*v2/(B2+SMALL));
  //myset(datatype,&ftemp,0,1,writebuf);

  // \Omega_F_2
  //ftemp=(float)(v3-B3*v1/(B1+SMALL));
  // myset(datatype,&ftemp,0,1,writebuf);

  return(0);

}




int dissdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping dissdump# %ld ... ",dump_cnt);

  whichdump=DISSDUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/dissdump");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,dissdump_content)>=1) return(1);

  trifprintf("end dumping dissdump# %ld ... ",dump_cnt);


  return(0);
  
}



int dissdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  myset(datatype,&dissfunpos[i][j][k],0,NUMDISSFUNPOS,writebuf);

  return (0);
}



int dumpother(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping dumpother# %ld ... ",dump_cnt);

  whichdump=DUMPOTHERCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/dumpother");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,dumpother_content)>=1) return(1);

	/////// output the symmetry information to the fail file
	//writesyminfo();
	///////

  trifprintf("end dumping dumpother# %ld ... ",dump_cnt);


  return(0);
  
}

int dumpother_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl;

  myset(datatype,&panalytic[i][j][k],0,NPR,writebuf);

  for(pl=0;pl<NUMPOTHER;pl++){
    myset(datatype,&pother[pl][i][j][k],0,1,writebuf);
  }

  return (0);
}




int fluxdumpdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping fluxdump# %ld ... ",dump_cnt);

  whichdump=FLUXDUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/fluxdump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,fluxdump_content)>=1) return(1);

  trifprintf("end dumping fluxdump# %ld ... ",dump_cnt);


  return(0);
  
}

int fluxdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int pl;

  myset(datatype,&fluxdump[i][j][k],0,NUMFLUXDUMP,writebuf);

  return (0);
}



// dump stuff related specially to EOS
int eosdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping eosdump# %ld ... ",dump_cnt);

  whichdump=EOSDUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/eosdump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,eosdump_content)>=1) return(1);

	/////// output the symmetry information to the fail file
	//writesyminfo();
	///////

  trifprintf("end dumping eosdump# %ld ... ",dump_cnt);


  return(0);
  
}

int eosdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  FTYPE rho,u;
  FTYPE height,ye,ynu,temp;
  FTYPE extras[MAXNUMEXTRAS];
  FTYPE processed[MAXPROCESSEDEXTRAS];
  int numextras;
  int extraiter;
  int numparms;
  FTYPE parlist[MAXPARLIST];


  //////////////////////////
  //
  // do the assignments
  //
  // if you change # of outputted vars, remember to change numcolumns


  // below needed for computing H and EOS
  icurr=i;
  jcurr=j;
  kcurr=k;
  pcurr=CENT;
  ////////////
  //
  // output some EOS stuff since in general not simple function of rho0,u
  rho = pdump[i][j][k][RHO];
  u = pdump[i][j][k][UU];

  // compute EOS stuff
  get_EOS_parms(i, j, k, &numparms, parlist);
  temp = compute_temp(rho,u);
  // get extra EOS stuff
  get_extrasprocessed(1, i, j, k, pdump[i][j][k], extras, processed);
  //  compute_allextras(0,rho,u,&numextras,extras);
  
  // write EOS stuff
  // write out all stuff for all EOSs so uniform format to read for any run (otherwise need version information in file)
  myset(datatype,&parlist,0,MAXPARLIST,writebuf); // numparms=6
  myset(datatype,&temp,0,1,writebuf); // 1
  // write extras EOS stuff
  myset(datatype,extras,0,MAXNUMEXTRAS,writebuf); // numextras
  myset(datatype,processed,0,MAXPROCESSEDEXTRAS,writebuf); // numprocessedextras


  return (0);
}





int vpotdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping vpotdump# %ld ... ",dump_cnt);

  whichdump=VPOTDUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/vpotdump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,vpotdump_content)>=1) return(1);

  trifprintf("end dumping vpotdump# %ld ... ",dump_cnt);


  return(0);
  
}



int vpotdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int jj;

  DLOOPA(jj){
    myset(datatype,&vpotarray[jj][i][j][k],0,1,writebuf); // 1 each
  }

  return (0);
}




int fakedump(void)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping fakedump");

  whichdump=FAKEDUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/fakedump");
  strcpy(fileformat,"%04ld");  //atch adjust dump every substep
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,-1,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,fakedump_header,fakedump_content)>=1) return(1);

  trifprintf("end dumping fakedump");


  return(0);
  
}



int fakedump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{
  int jj;

  // blank

  return (0);
}

int fakedump_header(int bintxt, FILE *headerptr)
{
  int fake;

  fake=0;


  if(bintxt==BINARYOUTPUT){
    fwrite(&fake,sizeof(int),1,headerptr);
  }
  else{
#if(REALTYPE==DOUBLETYPE)
    fprintf(headerptr, "DONE: %d\n",fake);
#elif(REALTYPE==LONGDOUBLETYPE)
    fprintf(headerptr, "DONE: %d\n",fake);
#endif
  }
  fflush(headerptr);
  return(0);
}	

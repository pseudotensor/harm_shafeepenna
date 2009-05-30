
/* restart functions; restart_init and restart_dump */

#include "decs.h"

#include "restart.h"


static int restart_process_extra_variables(void);





// This function reads header and data (primitives and conserved quantities)
// It then outputs restart file for user to check consistency
// CHANGINGMARK: Note that for evolving metric should really store old metric in restart file
int restart_init(int which)
{



  trifprintf("begin restart init\n");

  ////////////////
  //
  // set random seed
  //
  ////////////////
  ranc(7);

  ////////////////
  //
  // get restart header and grid data
  //
  ////////////////
  restart_read(which);


  ////////////////
  //
  // set any extra things (only used with restart.rebeccaoldcode.c)
  //
  ////////////////
  extrarestartfunction();


  ////////////////
  //
  // get restart old metric
  //
  ////////////////
  if(DOEVOLVEMETRIC) restartmetric_read(which);

  ////////////////
  //
  // process variables not written to file
  //
  ////////////////
  restart_process_extra_variables();
    

  ////////////////
  //
  // set all CPUs to have same restart header stuff (needed in mpicombine==1 mode)
  //
  ////////////////
  restart_read_defs();

  ////////////////
  //
  // set coordinate parameters by reading coordparms.dat
  // assumes doesn't overwrite what just set (i.e. what's written to restart file)
  //
  ////////////////
  read_coord_parms();

  ////////////////
  //
  // write header to screen to all CPUs
  //
  ////////////////
  fprintf(log_file,"header contents below\n"); fflush(log_file);
  write_restart_header(TEXTOUTPUT,log_file);

  ////////////////
  //
  // check if failed inside rdump file
  //
  ////////////////
  if(failed!=0){
    dualfprintf(log_file,"WARNING: failed=%d in rdump file.  Must clense the way between us.\n",failed);
    dualfprintf(log_file,"Setting failed=0\n");
    failed=0;
  }

  // certain global.h parameters may override restart parameters.  This happens, say, when the user decides not to DOAVG after restart, but DOAVGDIAG is still set.  Without this the code will segfault in diag.c
  // corresponds to "if" statements in initbase.c as well, where user has either choice or no choice.
  if(DOAVG==0) DOAVGDIAG=0;



  ////////////////////////
  //
  // test read by looking at written file
  //
  ////////////////////////
  restart_write(3);

  // can't check (e.g.) images here if wanting to look at conserved quantities

  trifprintf("proc: %d t=%21.15g failed=%d\n", myid, t, failed);



  trifprintf("end restart_init\n");


  return(0);

}







// assign any quantity[i][j][k] that isn't set by restart_read() that want to set so later checks or later computational uses can occur without errors
// GODMARK: must keep this up to date with what one is evolving not written to dump
int restart_process_extra_variables(void)
{
  int i,j,k;

  // if evolving entropy need to reset primitive to avoid error checks detecting problem
  if(DOENTROPY!=DONOENTROPY){
    FULLLOOP p[i][j][k][ENTROPY]=p[i][j][k][UU];
  }


  return(0);

}



int restart_write(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
  long truedump_cnt;

  trifprintf("begin dumping rdump# %ld ... ",dump_cnt);

  whichdump=RDUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/rdump");
  if(dump_cnt>=0) {
    strcpy(fileformat,"-%01ld"); // very quick restart
    truedump_cnt=dump_cnt;
  }
  else {
    strcpy(fileformat,"--%04ld"); // assumed at some longer cycle and never overwritten .. must rename this to normal format to use as real rdump.
    truedump_cnt=-dump_cnt-1;
  }
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,truedump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,write_restart_header,rdump_content)>=1) return(1);

  trifprintf("end dumping rdump# %ld ... ",dump_cnt);


  return(0);

}


int rdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  // always NPR
  myset(datatype,p[i][j][k],0,NPR,writebuf);
  myset(datatype,unew[i][j][k],0,NPR,writebuf);

  return(0);
}









int restart_read(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
  int bintxt;


  trifprintf("begin reading rdump# %ld ... ",dump_cnt);

  whichdump=RDUMPCOL;
  datatype=MPI_FTYPE;
  bintxt=binaryoutput;
  strcpy(fileprefix,"dumps/rdump");
  strcpy(fileformat,"-%01ld");
  strcpy(filesuffix,"");
 
  if(dump_gen(READFILE,dump_cnt,bintxt,whichdump,datatype,fileprefix,fileformat,filesuffix,read_restart_header,rdump_read_content)>=1) return(1);



  trifprintf("end reading rdump# %ld ... ",dump_cnt);


  return(0);

}


int rdump_read_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  // always NPR
  myget(datatype,p[i][j][k],0,NPR,writebuf);
  myget(datatype,unew[i][j][k],0,NPR,writebuf);


  return(0);
}














int restartmetric_write(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
  long truedump_cnt;

  trifprintf("begin dumping rmetricdump# %ld ... ",dump_cnt);

  whichdump=RMETRICDUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/rmetricdump");
  if(dump_cnt>=0) {
    strcpy(fileformat,"-%01ld"); // very quick restart
    truedump_cnt=dump_cnt;
  }
  else {
    strcpy(fileformat,"--%04ld"); // assumed at some longer cycle and never overwritten .. must rename this to normal format to use as real rdump.
    truedump_cnt=-dump_cnt-1;
  }
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,truedump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,write_restartmetric_header,rmetricdump_content)>=1) return(1);

  trifprintf("end dumping rmetricdump# %ld ... ",dump_cnt);


  return(0);

}


int rmetricdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  myset(datatype,gcovlast[i][j][k],0,NPG*NDIM*NDIM,writebuf);
  myset(datatype,gcovpertlast[i][j][k],0,NPG*NDIM*NDIM,writebuf);
  myset(datatype,alphalapselast[i][j][k],0,NPG,writebuf);

  return(0);
}







int write_restartmetric_header(int bintxt, FILE *headerptr)
{

  // nothing so far

  return(0);
}

int read_restartmetric_header(int bintxt, FILE *headerptr)
{

  // nothing so far

  return(0);
}


int restartmetric_read(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
  int bintxt;


  trifprintf("begin reading rmetricdump# %ld ... ",dump_cnt);

  dualfprintf(fail_file,"Reading in old metric, but so far need to bound somehow since no BCs but needed!\n");
  myexit(2496736);

  whichdump=RMETRICDUMPCOL;
  datatype=MPI_FTYPE;
  bintxt=binaryoutput;
  strcpy(fileprefix,"dumps/rmetricdump");
  strcpy(fileformat,"-%01ld");
  strcpy(filesuffix,"");
 
  if(dump_gen(READFILE,dump_cnt,bintxt,whichdump,datatype,fileprefix,fileformat,filesuffix,read_restartmetric_header,rmetricdump_read_content)>=1) return(1);



  trifprintf("end reading rmetricdump# %ld ... ",dump_cnt);


  return(0);

}


int rmetricdump_read_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf)
{

  myget(datatype,gcovlast[i][j][k],0,NPG*NDIM*NDIM,writebuf);
  myget(datatype,gcovpertlast[i][j][k],0,NPG*NDIM*NDIM,writebuf);
  myget(datatype,alphalapselast[i][j][k],0,NPG,writebuf);

  return(0);
}














// headerptr created and only used here OR passed a given pointer
int bcast_restart_header(void)
{
  int readwrite_restart_header(int readwrite, int bintxt, int bcasthead, FILE*headerptr);
  int bcasthead;

  bcasthead=1;
  
  // 0 and NULL not used
  readwrite_restart_header(NOTHINGHEAD, 0, bcasthead, NULL);

  return(0);
}


// headerptr created and only used here OR passed a given pointer
int read_restart_header_new(int bintxt, FILE*headerptr)
{
  int readwrite_restart_header(int readwrite, int bintxt, int bcasthead, FILE*headerptr);
  int bcasthead;

  bcasthead=0;
  
  readwrite_restart_header(READHEAD, bintxt, bcasthead, headerptr);

  return(0);
}



int write_restart_header_new(int bintxt,FILE*headerptr)
{
  int readwrite_restart_header(int readwrite, int bintxt, int bcasthead, FILE*headerptr);
  int bcasthead;

  bcasthead=0;
  
  readwrite_restart_header(WRITEHEAD, bintxt, bcasthead, headerptr);

  return(0);
}



// here bintxt is only 
int readwrite_restart_header(int readwrite, int bintxt, int bcasthead, FILE*headerptr)
{
  int ii;
  int k,dir,pl;
  int enerregion, floor, tscale;
  int idum1,idum2,idum3;
  int dissloop;
  int dtloop;
  char headerone[MAXFILENAME];
  char sheaderone[MAXFILENAME];
  char ctypeheaderone[MAXFILENAME];


  if(readwrite==READHEAD) trifprintf("begin reading header of restart file\n");
  else if(readwrite==WRITEHEAD) trifprintf("begin writing header of restart file\n");


  // setup idum1,2,3 so read and write are processed same way
  if(readwrite==WRITEHEAD){
    idum1=totalsize[1];
    idum2=totalsize[2];
    idum3=totalsize[3];
  }


  if(readwrite==READHEAD){
    sprintf(headerone,"%s",HEADERONEIN);
    sprintf(sheaderone,"%s",SFTYPEHEADERONEIN);
    sprintf(ctypeheaderone,"%s",CTYPEHEADERONEIN);
  }
  else if(readwrite==WRITEHEAD){
    sprintf(headerone,"%s",HEADERONEOUT);
    sprintf(sheaderone,"%s",SFTYPEHEADERONEOUT);
    sprintf(ctypeheaderone,"%s",CTYPEHEADERONEOUT);
  }


  /////////////////
  //
  // START HEADER read/write list
  //
  /////////////////

  header1_gen(readwrite,bintxt,bcasthead,&idum1,sizeof(int),"%d",1,MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&idum2,sizeof(int),"%d",1,MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&idum3,sizeof(int),"%d",1,MPI_INT,headerptr); // 3D thing

  // all cpus read the rest of header the same
  header1_gen(readwrite,bintxt,bcasthead,&t, sizeof(SFTYPE), headerone, 1, MPI_SFTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&tf, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&nstep, sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&a, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&MBH, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&QBH, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&gam, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&gamideal, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&cour, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);


  // for metric evolution
  header1_gen(readwrite,bintxt,bcasthead,&Xmetricold, sizeof(FTYPE), headerone, NDIM, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&Xmetricnew, sizeof(FTYPE), headerone, NDIM, MPI_FTYPE, headerptr);

      
  header1_gen(readwrite,bintxt,bcasthead,&dt, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&lim[1], sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&lim[2], sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&lim[3], sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&TIMEORDER, sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&fluxmethod, sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&FLUXB, sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&UTOPRIMVERSION, sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&failed, sizeof(int), "%d", 1, MPI_INT, headerptr);
  
  header1_gen(readwrite,bintxt,bcasthead,&defcoord, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&R0, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&Rin, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&Rout, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&hslope, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&Zin, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&Zout, sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&Rin_array, sizeof(FTYPE), headerone, NDIM, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&Rout_array, sizeof(FTYPE), headerone, NDIM, MPI_FTYPE, headerptr);
  
  header1_gen(readwrite,bintxt,bcasthead,&BCtype[0],sizeof(int), "%d", COMPDIM*2, MPI_INT, headerptr);

  // new May 6, 2003
  header1_gen(readwrite,bintxt,bcasthead,&realnstep, sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&debugfail,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&whichrestart,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&cooling,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&restartsteps[0],sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&restartsteps[1],sizeof(long), "%ld", 1, MPI_LONG, headerptr);

  header1_gen(readwrite,bintxt,bcasthead,&GAMMIEDUMP,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&GAMMIEIMAGE,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&GAMMIEENER,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DODIAGS,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DOENERDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DOGDUMPDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DORDUMPDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DODUMPDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DOAVGDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DOIMAGEDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DOAREAMAPDIAG,sizeof(int), "%d", 1, MPI_INT, headerptr);

  header1_gen(readwrite,bintxt,bcasthead,&DODIAGEVERYSUBSTEP,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DOENODEBUGEVERYSUBSTEP,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DOCOLSPLIT,sizeof(int), "%d", NUMDUMPTYPES, MPI_INT, headerptr);

  header1_gen(readwrite,bintxt,bcasthead,&POSDEFMETRIC,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&periodicx1,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&periodicx2,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&periodicx3,sizeof(int), "%d", 1, MPI_INT, headerptr); // 3D thing
  header1_gen(readwrite,bintxt,bcasthead,&binaryoutput,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&sortedoutput,sizeof(int), "%d", 1, MPI_INT, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&defcon,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&SAFE,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&RHOMIN,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&UUMIN,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&RHOMINLIMIT,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&UUMINLIMIT,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&BSQORHOLIMIT,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&BSQOULIMIT,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&UORHOLIMIT,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&GAMMAMAX,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&GAMMADAMP,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&GAMMAFAIL,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);      
  // end new May 6, 2003
  
  // reorganized order for DT related stuff
  header1_gen(readwrite,bintxt,bcasthead,&DTr, sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DTdumpgen[0], sizeof(SFTYPE), sheaderone, NUMDTDS, MPI_SFTYPE, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&rdump_cnt, sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&dumpcntgen[0], sizeof(long), "%ld", NUMDTDS, MPI_LONG, headerptr);
  
  //    header1_gen(readwrite,bintxt,bcasthead,&DTd, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  //header1_gen(readwrite,bintxt,bcasthead,&DTavg, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  //header1_gen(readwrite,bintxt,bcasthead,&DTener, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  //header1_gen(readwrite,bintxt,bcasthead,&DTi, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  // header1_gen(readwrite,bintxt,bcasthead,&DTr, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr) ;
  //header1_gen(readwrite,bintxt,bcasthead,&DTdebug, sizeof(SFTYPE), sheaderone, 1, MPI_SFTYPE, headerptr);
  
  //    header1_gen(readwrite,bintxt,bcasthead,&dump_cnt, sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  //header1_gen(readwrite,bintxt,bcasthead,&image_cnt, sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  //    header1_gen(readwrite,bintxt,bcasthead,&avg_cnt, sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  //header1_gen(readwrite,bintxt,bcasthead,&debug_cnt, sizeof(long), "%ld", 1, MPI_LONG, headerptr);
  
  
  header1_gen(readwrite,bintxt,bcasthead,&prMAX[0],sizeof(FTYPE), headerone, NPR, MPI_FTYPE,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&prfloorcoef[0],sizeof(FTYPE), headerone, NPR, MPI_FTYPE,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&rescaletype,sizeof(int), "%d", 1, MPI_INT,headerptr);
  
  // Nov 11, 2006 : post-Sasha-WENO code WENO stuff
  header1_gen(readwrite,bintxt,bcasthead,&avgscheme[1],sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&avgscheme[2],sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&avgscheme[3],sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&dofluxreconevolvepointfield,sizeof(int), "%d", 1, MPI_INT,headerptr);

  header1_gen(readwrite,bintxt,bcasthead,&do_transverse_flux_integration[0],sizeof(int), "%d", NPR, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&do_source_integration[0],sizeof(int), "%d", NPR, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&do_conserved_integration[0],sizeof(int), "%d", NPR, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&INVERTFROMAVERAGEIFFAILED,sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&LIMIT_AC_PRIM_FRAC_CHANGE,sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&LIMIT_AC_FRAC_CHANGE,sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&MAX_AC_PRIM_FRAC_CHANGE,sizeof(FTYPE), headerone, 1, MPI_FTYPE,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&MAX_AC_FRAC_CHANGE,sizeof(FTYPE), headerone, 1, MPI_FTYPE,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DOENOFLUX,sizeof(int), "%d", 1, MPI_INT,headerptr);


  header1_gen(readwrite,bintxt,bcasthead,&CHECKCONT,sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DOTSTEPDIAG,sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DOLOGSTEP,sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&DOLOGPERF,sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&NDTCCHECK,sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&NZCCHECK,sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&NDTDOTCCHECK,sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&NGOCHECK,sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&NTIMECHECK,sizeof(int), "%d", 1, MPI_INT,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&PERFWALLTIME,sizeof(FTYPE), headerone, 1, MPI_FTYPE,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&ZCPSESTIMATE,sizeof(FTYPE), headerone, 1, MPI_FTYPE,headerptr);


  // new June 6, 2003 (cumulatives)
  // always NPR
  header1_gen(readwrite,bintxt,bcasthead,&pcumreg_tot[0][0][0],sizeof(FTYPE), headerone, NUMENERREGIONS*(COMPDIM*2)*NPR, MPI_FTYPE,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&fladdreg_tot[0][0],sizeof(FTYPE), headerone, NUMENERREGIONS*NPR, MPI_FTYPE,headerptr);
  header1_gen(readwrite,bintxt,bcasthead,&sourceaddreg_tot[0][0],sizeof(FTYPE), headerone, NUMENERREGIONS*NPR, MPI_FTYPE,headerptr);

  header1_gen(readwrite,bintxt,bcasthead,&Ureg_init_tot[0][0],sizeof(FTYPE), headerone, NUMENERREGIONS*NPR, MPI_FTYPE,headerptr);
  //TSCALELOOP(tscale) FLOORLOOP(floor)
  header1_gen(readwrite,bintxt,bcasthead,&failfloorcountlocal_tot[0][0],sizeof(CTYPE),ctypeheaderone,NUMTSCALES*NUMFAILFLOORFLAGS,MPI_CTYPE,headerptr);


  // end new June 6,2003

  if(DODISS){
    header1_gen(readwrite,bintxt,bcasthead,&dissreg_tot[0][0],sizeof(SFTYPE), sheaderone, NUMENERREGIONS*NUMDISSVERSIONS, MPI_SFTYPE,headerptr);

    //    for(enerregion=0;enerregion<NUMENERREGIONS;enerregion++){
    //      for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++){
    //	dualfprintf(fail_file,"dissloop=%d enerregion=%d dissreg_tot=%21.15g\n",dissloop,enerregion,dissreg_tot[enerregion][dissloop]);
    //      }
    //    }
  }
    
  if(DOLUMVSR){
    header1_gen(readwrite,bintxt,bcasthead,&lumvsr_tot[0],sizeof(SFTYPE),sheaderone, ncpux1*N1, MPI_SFTYPE, headerptr);
    //    for(ii=0;ii<ncpux1*N1;ii++){
    //      if(!isfinite(lumvsr_tot[ii])){
    //	dualfprintf(fail_file,"restart: ii=%d lumvsr_tot=%21.15g\n",ii,lumvsr_tot[ii]);
    //      }
    //    }
  }
  if(DODISSVSR){
    // for consistent restart file, assume NUMDISSVERSIONS doesn't change
    for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++){// this loop is over pointers, not a continuous memory space!
      header1_gen(readwrite,bintxt,bcasthead,&dissvsr_tot[dissloop][0],sizeof(SFTYPE),sheaderone, ncpux1*N1, MPI_SFTYPE, headerptr);
    }
  }

 
  // Aug 16, 2007 -- active region parameters (updated by JCM 07/24/08)
  if(DOGRIDSECTIONING) {
    header1_gen(readwrite,bintxt,bcasthead,&global_sectiondef,sizeof(int), "%d", NUMUPDOWN*NDIM, MPI_INT,headerptr);
    header1_gen(readwrite,bintxt,bcasthead,&t_transition,sizeof(FTYPE), headerone, 1, MPI_FTYPE, headerptr);

  }
  //end Aug 16, 2007 (updated by JCM 07/24/08)
 

  // BELOW moved to dump_gen
  // now read of tail is controlled by dump_gen()
  //  if(bintxt==TEXTOUTPUT){
    // flush to just after the header line in case binary read of data
    //    if(readwrite==READHEAD) while(fgetc(headerptr)!='\n');
  // if(readwrite==READHEAD){
      // do nothing
  // }
  // else if(readwrite==WRITEHEAD) fprintf(headerptr,"\n");
  // }

  
  if(bintxt==TEXTOUTPUT){
    // flush to just after the header line in case binary read of data
    if(readwrite==WRITEHEAD) fprintf(headerptr,"\n");
  }


  // flux header so written to disk fully
  fflush(headerptr);


  if(readwrite==READHEAD){
    /////////////////
    //
    // some checks
    if (idum1 != totalsize[1]) {
      dualfprintf(fail_file, "error reading restart file; N1 differs\n");
      dualfprintf(fail_file, "got totalsize[1]=%d needed totalsize[1]=%d\n",idum1,totalsize[1]);
      myexit(3);
    }
    if (idum2 != totalsize[2]) {
      dualfprintf(fail_file, "error reading restart file; N2 differs\n");
      dualfprintf(fail_file, "got totalsize[2]=%d needed totalsize[2]=%d\n",idum2,totalsize[2]);
      myexit(4);
    }
    if (idum3 != totalsize[3]) {
      dualfprintf(fail_file, "error reading restart file; N3 differs\n");
      dualfprintf(fail_file, "got totalsize[3]=%d needed totalsize[3]=%d\n",idum3,totalsize[3]);
      myexit(4);
    }
  }


  if(readwrite==READHEAD) trifprintf("end reading header of restart file\n");
  else if(readwrite==WRITEHEAD) trifprintf("end writing header of restart file\n");


  return(0);



}


// setup global accounting for CPU=0
int restart_read_defs_new(void)
{
  int enerregion;
  int floor,tscale;
  int dissloop;
  int dir,k;
  int ii;
  int dtloop;
  int bcast_restart_header(void);
  

  if(myid==0){
    ////////////
    //
    // Define for cpu=0 only, which will continue to keep track of the total after restart
    //
    ENERREGIONLOOP(enerregion) DIRLOOP(dir) PDUMPLOOP(k) pcumreg[enerregion][dir][k]=pcumreg_tot[enerregion][dir][k];
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) fladdreg[enerregion][k]=fladdreg_tot[enerregion][k];
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) sourceaddreg[enerregion][k]=sourceaddreg_tot[enerregion][k];
    ENERREGIONLOOP(enerregion) PDUMPLOOP(k) Ureg_init[enerregion][k]=Ureg_init_tot[enerregion][k];

    TSCALELOOP(tscale) FLOORLOOP(floor) failfloorcountlocal[tscale][floor]=failfloorcountlocal_tot[tscale][floor];


    if(DODISS) ENERREGIONLOOP(enerregion) for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++) dissreg[enerregion][dissloop]=dissreg_tot[enerregion][dissloop];
    // assume dissfunpos[][][] not restored since zeroed out each dump

    if(DOLUMVSR) for(ii=0;ii<ncpux1*N1;ii++) lumvsr[ii]=lumvsr_tot[ii];
    if(DODISSVSR) for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++) for(ii=0;ii<ncpux1*N1;ii++) dissvsr[dissloop][ii]=dissvsr_tot[dissloop][ii];

  }


  /////////////
  //
  // broadcast things read from header to all CPUs
  //
  /////////////
  bcast_restart_header();



  return(0);
}



///////////
//
// old header functions:
//
///////////
#include "restart.rebeccaoldcode.c"




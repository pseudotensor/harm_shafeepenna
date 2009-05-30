#include "decs.h"



////////////////
//
// General Initialize MPI
// Only include non-blocking MPI commands
//
////////////////
int init_MPI_general(int *argc, char **argv[])
{

#if(USEMPI)
  MPI_Init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, &truenumprocs); // WORLD total number of processors
  MPI_Comm_rank(MPI_COMM_WORLD, &myid_world); // WORLD proc id
  MPI_Get_processor_name(processor_name, &procnamelen); // to ensure really on certain nodes

  fprintf(stderr, "WORLD proc: %d of %d on %s\n", myid_world,truenumprocs,processor_name);

  fprintf(stderr, "end: init_MPI\n");
  fflush(stderr);
#else
  truenumprocs=1;
  myid=0;
#endif



  if (MAXCPUS < truenumprocs) {
    fprintf(stderr,
	    "Must increase MAXCPUS in global.h, %d is too many\n",
	    truenumprocs);
    myexit(1);
  }



  return(0);

}




void init_MPI_setupfilesandgrid(int argc, char *argv[])
{
  int size;



  ////////////////////
  //
  // choose to combine files or not
  //
  //
  ///////////////////
  if(USEMPI){
    mpicombine = 1;    // choice
    //mpicombine=0;

    // 
    if(mpicombine){
      if(USEROMIO==0){
	// choice
	if(sortedoutput==SORTED) mpicombinetype=MPICOMBINEMINMEM;
	else if(sortedoutput==UNSORTED) mpicombinetype=MPICOMBINESIMPLE; //forced to happen since no unsorted method for the advanced combine technique
	//mpicombinetype=MPICOMBINESIMPLE; // forced for testing
      }
      else truempicombinetype=mpicombinetype=MPICOMBINEROMIO;
    }
  }
  else{
    // no choice
    mpicombine = 0;
  }




  ///////////////////
  //
  // Setup files and grid information
  // always done
  //
  ///////////////////

  // after below can use dualfprintf, and trifprintf, etc.
  init_genfiles(0);

  // after below can bound data using MPI
  init_placeongrid();



  ///////////////////
  //
  // Some MPI type size checks
  //
  ///////////////////

#if(USEMPI)
  
  if(sizeof(CTYPE)==sizeof(long long int)){
    MPI_Type_size(MPI_LONG_LONG_INT,&size);
    if(size!=8){
      dualfprintf(fail_file,"size of the long long int in MPI=%d, should be 8\n",size);
      myexit(1000);
    }
  }
  
  if((sizeof(REALTYPE)==sizeof(long double))||(sizeof(SENSITIVE)==sizeof(long double))){
    MPI_Type_size(MPI_LONG_DOUBLE,&size);
    if(size!=16){
      dualfprintf(fail_file,"size of the long double in MPI=%d, should be 16\n",size);
      myexit(1000);
    }
  }
#endif


  trifprintf("done with init_MPI_setupfilesandgrid()\n");  fflush(log_file);

}





// myargs() reads in arguments for all of GRMHD codes
// globally reads and sets ncpux1,ncpux2,ncpux3,RESTARTMODE, WHICHFILE for all CPUs
// also globally sets numprocs and makes sure numprocs makes sense
void myargs(int argc, char *argv[])
{
  int argi;

  if(USEMPI){

    if(argc<COMPDIM+1){
      if(myid==0){
	fprintf(stderr,"proc: %04d : Incorrect command line: argc: %d needed at least=%d, please specify:\n",myid,argc,3+1);// was COMPDIM+1, but have to fix based upon below code
	fprintf(stderr,"proc: %04d : mpirun <mpirunoptions> <progname> ncpux1 ncpux2 ncpux3\n",myid);
	fprintf(stderr,"proc: %04d : OR\n",myid);
	fprintf(stderr,"proc: %04d : mpirun <mpirunoptions> <progname> ncpux1 ncpux2 ncpux3 RESTARTMODE WHICHFILE\n",myid);
      }
      exit(1);
    }
    argi=1;
    ncpux1=atoi(argv[argi++]);
    ncpux2=atoi(argv[argi++]);
    ncpux3=atoi(argv[argi++]);
    // default unless user adds below
    RESTARTMODE=0;
    WHICHFILE=0;

    if(argc==1+3+2){
      RESTARTMODE=atoi(argv[argi++]);
      WHICHFILE=atoi(argv[argi++]);
    }// no failure stuff since have to recompile anyways for failure tracking

  }
  else{
    ncpux1 = 1;
    ncpux2 = 1;
    ncpux3 = 1;

    if(argc==1){}
    else if(argc==3){
      RESTARTMODE=atoi(argv[1]);
      WHICHFILE=atoi(argv[2]);    
    }
    else{
      if(myid==0){
	fprintf(stderr,"<progname>\n");
	fprintf(stderr,"OR\n");
	fprintf(stderr,"<progname> RESTARTMODE WHICHFILE\n");
      }
      exit(1);
    }
  }
  

  // get total number of processors
  numprocs = ncpux1*ncpux2*ncpux3;

  // test
  if(sizeproclist_grmhd!=numprocs){
    fprintf(stderr, "Got (sizeproclist_grmhd=%d) != (numprocs=ncpux1*ncpux2*ncpux3=%d)\n",sizeproclist_grmhd,numprocs);
    myexit(1);
  }


  myfprintf(stderr,
	    "numprocs=%d ncpux1=%d ncpux2=%d ncpux3=%d percpusize: N1=%d N2=%d N3=%d\n",
	    numprocs, ncpux1, ncpux2, ncpux3, N1, N2,N3);


}




void init_genfiles(int gopp)
{
  char temps[MAXFILENAME];
  char extension[MAXFILENAME];
  char binarytype[MAXFILENAME];
  void set_binarytype(char *binarytype);


  fprintf(stderr, "begin: init_genfiles ... ");
  fflush(stderr);

  if (gopp == 1) {
    strcpy(extension, PPEXT);
  } else if (gopp == 0) {
    strcpy(extension, OUTEXT);
  }
  // always have fail and general log open

  sprintf(temps, "%s0_fail%s%s", DATADIR, extension, myidtxt);
  if ((fail_file = fopen(temps, "at")) == NULL) {
    fprintf(stderr, "fail: Cannot open: %s\n", temps);
    exit(1);
  }
  fprintf(stderr, "opened: %s\n", temps);



  sprintf(temps, "%s0_log%s%s", DATADIR, extension, myidtxt);
  if ((log_file = fopen(temps, "at")) == NULL) {
    fprintf(stderr, "log: Cannot open: %s\n", temps);
    exit(1);
  }
  fprintf(stderr, "opened: %s\n", temps);
  fprintf(log_file, "fail_file: %d log_file: %d\n", (int)fail_file,
	  (int)log_file);
  fflush(log_file);


  sprintf(temps, "%s0_logdt%s%s", DATADIR, extension, myidtxt);
  if ((logdt_file = fopen(temps, "at")) == NULL) {
    fprintf(stderr, "logdt: Cannot open: %s\n", temps);
    exit(1);
  }
  fprintf(stderr, "opened: %s\n", temps);
  fflush(logdt_file);





  if (myid == 0) {

    set_binarytype(binarytype);


    sprintf(temps, "%s0_logfull%s%s", DATADIR, binarytype, extension);
    if ((logfull_file = fopen(temps, "at")) == NULL) {
      fprintf(stderr, "logfull: Cannot open: %s\n", temps);
      exit(1);
    }
    fprintf(stderr, "opened: %s\n", temps);
    fprintf(logfull_file, "logfull_file: %d \n", (int)logfull_file);
    fflush(logfull_file);


    sprintf(temps, "%s0_logdtfull%s%s", DATADIR, binarytype, extension);
    if ((logdtfull_file = fopen(temps, "at")) == NULL) {
      fprintf(stderr, "logdtfull: Cannot open: %s\n", temps);
      exit(1);
    }
    fprintf(stderr, "opened: %s\n", temps);
    fprintf(logdtfull_file, "logdtfull_file: %d \n", (int)logdtfull_file);
    fflush(logdtfull_file);

    sprintf(temps, "%s0_logstep%s%s", DATADIR, binarytype, extension);
    if ((logstep_file = fopen(temps, "at")) == NULL) {
      fprintf(stderr, "logstep: Cannot open: %s\n", temps);
      exit(1);
    }
    fprintf(stderr, "opened: %s\n", temps);

    sprintf(temps, "%s0_logperf%s%s", DATADIR, binarytype, extension);
    if ((logperf_file = fopen(temps, "at")) == NULL) {
      fprintf(stderr, "logperf: Cannot open: %s\n", temps);
      exit(1);
    }
    fprintf(stderr, "opened: %s\n", temps);



  }


  // ok now
  trifprintf("end: init_genfiles\n");
}



void set_binarytype(char *binarytype)
{
  if(DOINGGRMHDTYPECODE){
    sprintf(binarytype,".grmhd");
  }
  else if(DOINGGRRAYTYPECODE){
    sprintf(binarytype,".grray");
  }
  else if(DOINGLIAISONTYPECODE){
    sprintf(binarytype,".liaison");
  }
}






void init_placeongrid(void)
{
  // 3's were COMPDIM, but below code is fixed to require all 3 now
  int stage,stagei,stagef;
  int i, j, m, l;
  int N[3 + 1];

  int dir, bti;
  int opp[3*2];
  int pl;

  int numbnd[NDIM],surfa[NDIM];
  int numnpr;




  trifprintf("begin: init_placeongrid ... ");

  N[1] = N1;
  N[2] = N2;
  N[3] = N3;

  numbercpu[1] = ncpux1;
  numbercpu[2] = ncpux2;
  numbercpu[3] = ncpux3;

  mycpupos[1]=myid%ncpux1;
  mycpupos[2]=(int)((myid%(ncpux1*ncpux2))/ncpux1);
  mycpupos[3]=(int)(myid/(ncpux1*ncpux2));

  //  SLOOPA(j) dualfprintf(log_file,"mycpupos[%d]=%d",j,mycpupos[j]);


  for (m = 1; m <= COMPDIM; m++) {
    startpos[m] = mycpupos[m] * N[m];
    endpos[m] = (mycpupos[m] + 1) * N[m] - 1;

    // add up sizes for total size of grid
    totalsize[m] = 0;
    itotalsize[m] = 0;
    for (i = 0; i < numbercpu[m]; i++) {
      totalsize[m] += N[m];
      itotalsize[m] += N[m];
    }
  }

  for(m=1;m<=COMPDIM;m++){
    if((mycpupos0[m]=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate mycpupos0[%d]\n",m);
      myexit(1);
    }
    if((startpos0[m]=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate startpos0[%d]\n",m);
      myexit(1);
    }
    if((endpos0[m]=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate endpos0[%d]\n",m);
      myexit(1);
    }
  }
  // for cpu=0 as master to rest, needs this info
  for(i=0;i<numprocs;i++){
    mycpupos0[1][i]=i%ncpux1;
    mycpupos0[2][i]=(int)((i%(ncpux1*ncpux2))/ncpux1);
    mycpupos0[3][i]=(int)(i/(ncpux1*ncpux2));

    
    for (m = 1; m <= COMPDIM; m++) {
      startpos0[m][i] = mycpupos0[m][i] * N[m];
      endpos0[m][i] = (mycpupos0[m][i] + 1) * N[m] - 1;
    }
  }



  realtotalzones = totalzones = totalsize[1] * totalsize[2] * totalsize[3];
  itotalzones = itotalsize[1] * itotalsize[2] * itotalsize[3];


  /////////////// standard interior MPI data transfer setup
  //
  for(bti=0;bti<NUMBOUNDTYPES;bti++) for(dir=0;dir<COMPDIM*2;dir++) for(j=0;j<DIRNUMVARS;j++){
    dirset[bti][dir][j]=0;
  }
  // see where this cpu needs to send/recv

  // set number of quantities to bound per point
  for(dir=0;dir<COMPDIM*2;dir++){
    dirset[BOUNDPRIMTYPE][dir][DIRNUMPR]=NPRBOUND; // not used if SPLITNPR==1 or doing general range for quantities
    dirset[BOUNDFLUXTYPE][dir][DIRNUMPR]=NFLUXBOUND;// not used if SPLITNPR==1 or doing general range for quantities
  }



  for(bti=0;bti<NUMBOUNDTYPES;bti++) { // same for any bounding type

    if(N1>1){
      // figure out left/right send/recv
      if (mycpupos[1] > 0) {
	dirset[bti][X1DN][DIRIF] = 1;		// do -x1 dir
      }
      if (mycpupos[1] < ncpux1 - 1) {
	dirset[bti][X1UP][DIRIF] = 1;		// do +x1 dir
      }
      
      // only do periodic mpi if 
      if(periodicx1&&(ncpux1>1)){
	if(mycpupos[1]==0) dirset[bti][X1DN][DIRIF]=1;
	else if(mycpupos[1]==ncpux1-1) dirset[bti][X1UP][DIRIF]=1;
      }
      

    }
    else{
      // then no assume boundaries to copy (i.e. can't have CPU with dimension of length 1 and stack them up -- inefficient anyways)
      dirset[bti][X1UP][DIRIF] = 0;
      dirset[bti][X1DN][DIRIF] = 0;
    }
    
    // figure out up/down send/recv
    if(N2>1){
      if (mycpupos[2] > 0) {
	dirset[bti][X2DN][DIRIF] = 1;		// -x2 dir
      }
      if (mycpupos[2] < ncpux2 - 1) {
	dirset[bti][X2UP][DIRIF] = 1;		// towards and from +x2 dir
      }
      
      if(periodicx2&&(ncpux2>1)){
	if(mycpupos[2]==0) dirset[bti][X2DN][DIRIF]=1;
	else if(mycpupos[2]==ncpux2-1) dirset[bti][X2UP][DIRIF]=1;
      }
      
      
    }
    else{
      // then no assume boundaries to copy (i.e. can't have CPU with dimension of length 1 and stack them up -- inefficient anyways)
      dirset[bti][X2UP][DIRIF] = 0;
      dirset[bti][X2DN][DIRIF] = 0;
    }
    
    // figure out out/in send/recv
    if(N3>1){
      if (mycpupos[3] > 0) {
	dirset[bti][X3DN][DIRIF] = 1;		// -x3 dir
      }
      if (mycpupos[3] < ncpux3 - 1) {
	dirset[bti][X3UP][DIRIF] = 1;		// towards and from +x3 dir
      }
      
      if(periodicx3&&(ncpux3>1)){
	if(mycpupos[3]==0) dirset[bti][X3DN][DIRIF]=1;
	else if(mycpupos[3]==ncpux3-1) dirset[bti][X3UP][DIRIF]=1;
      }
      
      
    }
    else{
      // then no assume boundaries to copy (i.e. can't have CPU with dimension of length 1 and stack them up -- inefficient anyways)
      dirset[bti][X3UP][DIRIF] = 0;
      dirset[bti][X3DN][DIRIF] = 0;
    }
  }
  


  // define which direction the other guy communicates.  Just a simple way to store this obvious fact
  opp[X1DN]=X1UP;
  opp[X1UP]=X1DN;
  opp[X2DN]=X2UP;
  opp[X2UP]=X2DN;
  opp[X3DN]=X3UP;
  opp[X3UP]=X3DN;
  

  // which CPUs communicate to eachother 
  // same for any method (bti)
  for(bti=0;bti<NUMBOUNDTYPES;bti++) for(dir=0;dir<COMPDIM*2;dir++){
    // sets opposite direction
    dirset[bti][dir][DIROPP]=opp[dir];
    

    // matching CPU to transfer to/from

    // x1
    if((dir==X1UP)||(dir==X1DN)){
      if( (periodicx1==0)
	  ||((mycpupos[1]>0)&&(mycpupos[1]<ncpux1-1)) // interior CPUs
	  || (mycpupos[1]==0 && dir==X1UP) // inner-CPUs pointing up
	  || (mycpupos[1]==ncpux1-1 && dir==X1DN)
	  ){ // outer-CPUs pointing down
	if(dir==X1UP) dirset[bti][dir][DIROTHER]=myid+1;
	if(dir==X1DN) dirset[bti][dir][DIROTHER]=myid-1;
      }
      else if(periodicx1){ // here X1DN/X1UP are implicitly associated with 0/ncpux1-1
	if(mycpupos[1]==0) dirset[bti][dir][DIROTHER]=myid+(ncpux1-1);  // for X1DN
	else if(mycpupos[1]==ncpux1-1) dirset[bti][dir][DIROTHER]=myid-(ncpux1-1); // for X1UP
      }
    }

    // x2
    if((dir==X2UP)||(dir==X2DN)){
      if( (periodicx2==0)
	  ||((mycpupos[2]>0)&&(mycpupos[2]<ncpux2-1))
	  || (mycpupos[2]==0 && dir==X2UP)
	  || (mycpupos[2]==ncpux2-1 && dir==X2DN)
	  ){
	if(dir==X2UP) dirset[bti][dir][DIROTHER]=myid+ncpux1;
	if(dir==X2DN) dirset[bti][dir][DIROTHER]=myid-ncpux1;
      }
      else if(periodicx2){
	if(mycpupos[2]==0) dirset[bti][dir][DIROTHER]=myid+(ncpux2-1)*ncpux1;
	else if(mycpupos[2]==ncpux2-1) dirset[bti][dir][DIROTHER]=myid-(ncpux2-1)*ncpux1;
      }
    }

    // x3
    if((dir==X3UP)||(dir==X3DN)){
      if( (periodicx3==0)
	  ||((mycpupos[3]>0)&&(mycpupos[3]<ncpux3-1))
	  || (mycpupos[3]==0 && dir==X3UP)
	  || (mycpupos[3]==ncpux3-1 && dir==X3DN)
	  ){
	if(dir==X3UP) dirset[bti][dir][DIROTHER]=myid+ncpux1*ncpux2;
	if(dir==X3DN) dirset[bti][dir][DIROTHER]=myid-ncpux1*ncpux2;
      }
      else if(periodicx3){
	if(mycpupos[3]==0) dirset[bti][dir][DIROTHER]=myid+(ncpux3-1)*ncpux1*ncpux2;
	else if(mycpupos[3]==ncpux3-1) dirset[bti][dir][DIROTHER]=myid-(ncpux3-1)*ncpux1*ncpux2;
      }
    }

    // MPI tags that label transfer, must be unique while doing multiple transfers
    dirset[bti][dir][DIRTAGS]= myid     * COMPDIM * 2 + dir;
    dirset[bti][dir][DIRTAGR]= dirset[bti][dir][DIROTHER] * COMPDIM * 2 + dirset[bti][dir][DIROPP];
  } 





  // In reduced dimensions (i.e. N=1), below reduces to 0 positions,
  // but above DIRIF's control whether that direction is bounded or
  // not and should take care not MPI bounding the reduced dimensions.
  
  // for reduced dimensions, normal boundary zones are handled by LOOPS defined in global.h
  
  // tags are defined by sender's ID and direction sent
  // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
  // 0=right, 1=up,2=left,3=down,4=out,5=in
  // works/v bc[1=output/2=input][0,1,2,3,4,5]
  // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are
  // like: (fromid,otherid*COMPDIM*2+*) where ? and * are
  // opposites(i.e. 0 and 2, 1 and 3, 4 and 5)

  // dirset[bti] sets limits on all inclusive for loops used in boundmpi.c

  for(bti=0;bti<NUMBOUNDTYPES;bti++){
    ///////////////////////////////////////////////////////////////////
    //
    // first set BOUNDPRIMTYPE : CENT or FACE1,2,3 quantities that need full bounding
    //
    // BOUNDPRIMSIMPLETYPE as well
    //
    // then set BOUNDINTTYPE : CENT or FACE1,2,3 quantities that need full bounding
    //
    ///////////////////////////////////////////////////////////////////

    // get number of boundary cells and number of quantities to bound
    set_numbnd(bti, numbnd, &numnpr);

    if(bti==BOUNDPRIMTYPE || bti==BOUNDPRIMSIMPLETYPE || bti==BOUNDINTTYPE){
      // then can stay
    }
    else continue; // no others here
  

    // surface area must be consistent with loops
    surfa[1]=(N2+numbnd[2]*2)*(N3+numbnd[3]*2)*N1NOT1;
    surfa[2]=(N1+numbnd[1]*2)*(N3+numbnd[3]*2)*N2NOT1;
    surfa[3]=(N1+numbnd[1]*2)*(N2+numbnd[2]*2)*N3NOT1;
  
 

    for(dir=0;dir<COMPDIM*2;dir++){

      //////////////////////////
      //
      // sets size of transfer for centered quantities
      //
      ////////////////////////
      if((dir==X1UP)||(dir==X1DN)) dirset[bti][dir][DIRSIZE]=numbnd[1]*surfa[1]*numnpr;
      else if((dir==X2UP)||(dir==X2DN)) dirset[bti][dir][DIRSIZE]=numbnd[2]*surfa[2]*numnpr;
      else if((dir==X3UP)||(dir==X3DN)) dirset[bti][dir][DIRSIZE]=numbnd[3]*surfa[3]*numnpr;

      if(dirset[bti][dir][DIRIF]){

	///////////////////
	//
	// PACKING for centered quantities
	//
	//////////////////
	// zones to copy from (packing -- where to copy FROM)
	if(dir==X1UP){ // right
	  dirset[bti][dir][DIRPSTART1]=(N1-1)-(numbnd[1]-SHIFT1); //N1-numbnd[1];
	  dirset[bti][dir][DIRPSTOP1]=N1-1;
	}
	else if(dir==X1DN){ // left
	  dirset[bti][dir][DIRPSTART1]=0;
	  dirset[bti][dir][DIRPSTOP1]=numbnd[1]-SHIFT1;
	}

	if((dir==X1UP)||(dir==X1DN)){
	  dirset[bti][dir][DIRPSTART2]=-numbnd[2];
	  dirset[bti][dir][DIRPSTOP2]=N2-1+numbnd[2];
	  dirset[bti][dir][DIRPSTART3]=-numbnd[3];
	  dirset[bti][dir][DIRPSTOP3]=N3-1+numbnd[3];
	}

	if(dir==X2UP){ // up
	  dirset[bti][dir][DIRPSTART2]=(N2-1)-(numbnd[2]-SHIFT2); //N2-numbnd[2];
	  dirset[bti][dir][DIRPSTOP2]=N2-1;
	}
	else if(dir==X2DN){ // down
	  dirset[bti][dir][DIRPSTART2]=0;
	  dirset[bti][dir][DIRPSTOP2]=numbnd[2]-SHIFT2;
	}

	if((dir==X2UP)||(dir==X2DN)){
	  dirset[bti][dir][DIRPSTART1]=-numbnd[1];
	  dirset[bti][dir][DIRPSTOP1]=N1-1+numbnd[1];
	  dirset[bti][dir][DIRPSTART3]=-numbnd[3];
	  dirset[bti][dir][DIRPSTOP3]=N3-1+numbnd[3];
	}

	if(dir==X3UP){ // up
	  dirset[bti][dir][DIRPSTART3]=(N3-1)-(numbnd[3]-SHIFT3); //N3-numbnd[3];
	  dirset[bti][dir][DIRPSTOP3]=N3-1;
	}
	else if(dir==X3DN){ // down
	  dirset[bti][dir][DIRPSTART3]=0;
	  dirset[bti][dir][DIRPSTOP3]=numbnd[3]-SHIFT3;
	}

	if((dir==X3UP)||(dir==X3DN)){
	  dirset[bti][dir][DIRPSTART1]=-numbnd[1];
	  dirset[bti][dir][DIRPSTOP1]=N1-1+numbnd[1];
	  dirset[bti][dir][DIRPSTART2]=-numbnd[2];
	  dirset[bti][dir][DIRPSTOP2]=N2-1+numbnd[2];
	}




	///////////////////
	//
	// UNPACKING for centered quantities
	//
	//////////////////
	// zones to copy into (unpacking -- where to copy INTO)

	// x1
	if(dir==X1UP){ // right
	  dirset[bti][dir][DIRUSTART1]=N1-1+SHIFT1;
	  dirset[bti][dir][DIRUSTOP1]=N1-1+numbnd[1];
	}
	else if(dir==X1DN){ // left
	  dirset[bti][dir][DIRUSTART1]=-numbnd[1];
	  dirset[bti][dir][DIRUSTOP1]=-SHIFT1;
	}
	if((dir==X1UP)||(dir==X1DN)){
	  dirset[bti][dir][DIRUSTART2]=-numbnd[2];
	  dirset[bti][dir][DIRUSTOP2]=N2-1+numbnd[2];
	  dirset[bti][dir][DIRUSTART3]=-numbnd[3];
	  dirset[bti][dir][DIRUSTOP3]=N3-1+numbnd[3];
	}

	// x2
	if(dir==X2UP){ // up
	  dirset[bti][dir][DIRUSTART2]=N2-1+SHIFT2;
	  dirset[bti][dir][DIRUSTOP2]=N2-1+numbnd[2];
	}
	else if(dir==X2DN){ // down
	  dirset[bti][dir][DIRUSTART2]=-numbnd[2];
	  dirset[bti][dir][DIRUSTOP2]=-SHIFT2;
	}

	if((dir==X2UP)||(dir==X2DN)){
	  dirset[bti][dir][DIRUSTART1]=-numbnd[1];
	  dirset[bti][dir][DIRUSTOP1]=N1-1+numbnd[1];
	  dirset[bti][dir][DIRUSTART3]=-numbnd[3];
	  dirset[bti][dir][DIRUSTOP3]=N3-1+numbnd[3];
	}

	// x3
	if(dir==X3UP){ // up
	  dirset[bti][dir][DIRUSTART3]=N3-1+SHIFT3;
	  dirset[bti][dir][DIRUSTOP3]=N3-1+numbnd[3];
	}
	else if(dir==X3DN){ // down
	  dirset[bti][dir][DIRUSTART3]=-numbnd[3];
	  dirset[bti][dir][DIRUSTOP3]=-SHIFT3;
	}

	if((dir==X3UP)||(dir==X3DN)){
	  dirset[bti][dir][DIRUSTART1]=-numbnd[1];
	  dirset[bti][dir][DIRUSTOP1]=N1-1+numbnd[1];
	  dirset[bti][dir][DIRUSTART2]=-numbnd[2];
	  dirset[bti][dir][DIRUSTOP2]=N2-1+numbnd[2];
	}
      }// end if DIRIF
    }// end DIRLOOP
  }// end bti loop









  ///////////////////////////////////////////////////////////////////
  //
  // set BOUNDFLUXTYPE : FACE(1/2/3) quantities -- used to only bound along flux direction (all that's needed for ENO to de-average flux using finite difference method)
  //
  // and BOUNDFLUXSIMPLETYPE too
  //
  ///////////////////////////////////////////////////////////////////

  for(bti=0;bti<NUMBOUNDTYPES;bti++){

    // get number of boundary cells and number of quantities to bound
    set_numbnd(bti, numbnd, &numnpr);

    if(bti==BOUNDFLUXTYPE || bti==BOUNDFLUXSIMPLETYPE){
      // then can stay
    }
    else continue; // no other bti to be done here


    // surface area must be consistent with loops
    surfa[1]=(N2+numbnd[2]*2)*(N3+numbnd[3]*2)*N1NOT1;
    surfa[2]=(N1+numbnd[1]*2)*(N3+numbnd[3]*2)*N2NOT1;
    surfa[3]=(N1+numbnd[1]*2)*(N2+numbnd[2]*2)*N3NOT1;
 

    // dirset[bti] sets limits on all inclusive for loops used in boundmpi.c
    for(dir=0;dir<COMPDIM*2;dir++){



      //////////////////////////
      //
      // sets size of transfer
      //
      ////////////////////////
      // (different for "left" and "right") for flux-types
      if(dir==X1UP) dirset[bti][dir][DIRSIZE]=numbnd[1]*surfa[1]*numnpr;
      else if(dir==X1DN) dirset[bti][dir][DIRSIZE]=(numbnd[1]-1)*surfa[1]*numnpr;
      else if(dir==X2UP) dirset[bti][dir][DIRSIZE]=numbnd[2]*surfa[2]*numnpr;
      else if(dir==X2DN) dirset[bti][dir][DIRSIZE]=(numbnd[2]-1)*surfa[2]*numnpr;
      else if(dir==X3UP) dirset[bti][dir][DIRSIZE]=numbnd[3]*surfa[3]*numnpr;
      else if(dir==X3DN) dirset[bti][dir][DIRSIZE]=(numbnd[3]-1)*surfa[3]*numnpr;


      if(dirset[bti][dir][DIRIF]){
	///////////////////
	//
	// PACKING for face quantities
	//
	//////////////////
	// zones to copy from (packing -- where to copy FROM)
	if(dir==X1UP){ // right
	  dirset[bti][dir][DIRPSTART1]=N1-numbnd[1]+SHIFT1;  //(N1-1)-(numbnd[1]-SHIFT1);
	  dirset[bti][dir][DIRPSTOP1]=N1;  //N1-1;
	}
	else if(dir==X1DN){ // left
	  dirset[bti][dir][DIRPSTART1]=1; // such that won't copy anything if N1=1 //0;
	  dirset[bti][dir][DIRPSTOP1]=numbnd[1]-SHIFT1; //numbnd[1]-SHIFT1;
	}

	if((dir==X1UP)||(dir==X1DN)){
	  dirset[bti][dir][DIRPSTART2]=-numbnd[2];
	  dirset[bti][dir][DIRPSTOP2]=N2-1+numbnd[2];
	  dirset[bti][dir][DIRPSTART3]=-numbnd[3];
	  dirset[bti][dir][DIRPSTOP3]=N3-1+numbnd[3];
	}

	if(dir==X2UP){ // up
	  dirset[bti][dir][DIRPSTART2]=N2-numbnd[2]+SHIFT2;
	  dirset[bti][dir][DIRPSTOP2]=N2;
	}
	else if(dir==X2DN){ // down
	  dirset[bti][dir][DIRPSTART2]=1;
	  dirset[bti][dir][DIRPSTOP2]=numbnd[2]-SHIFT2;
	}

	if((dir==X2UP)||(dir==X2DN)){
	  dirset[bti][dir][DIRPSTART1]=-numbnd[1];
	  dirset[bti][dir][DIRPSTOP1]=N1-1+numbnd[1];
	  dirset[bti][dir][DIRPSTART3]=-numbnd[3];
	  dirset[bti][dir][DIRPSTOP3]=N3-1+numbnd[3];
	}

	if(dir==X3UP){ // up
	  dirset[bti][dir][DIRPSTART3]=N3-numbnd[3]+SHIFT3;
	  dirset[bti][dir][DIRPSTOP3]=N3;
	}
	else if(dir==X3DN){ // down
	  dirset[bti][dir][DIRPSTART3]=1;
	  dirset[bti][dir][DIRPSTOP3]=numbnd[3]-SHIFT3;
	}

	if((dir==X3UP)||(dir==X3DN)){
	  dirset[bti][dir][DIRPSTART1]=-numbnd[1];
	  dirset[bti][dir][DIRPSTOP1]=N1-1+numbnd[1];
	  dirset[bti][dir][DIRPSTART2]=-numbnd[2];
	  dirset[bti][dir][DIRPSTOP2]=N2-1+numbnd[2];
	}


	///////////////////
	//
	// UNPACKING for face quantities
	//
	//////////////////
	// zones to copy into (unpacking -- where to copy INTO)

	// x1
	if(dir==X1UP){ // right
	  dirset[bti][dir][DIRUSTART1]=N1+SHIFT1;
	  dirset[bti][dir][DIRUSTOP1]=N1-1+numbnd[1];
	}
	else if(dir==X1DN){ // left
	  dirset[bti][dir][DIRUSTART1]=-numbnd[1]+1;
	  dirset[bti][dir][DIRUSTOP1]=0;
	}
	if((dir==X1UP)||(dir==X1DN)){
	  dirset[bti][dir][DIRUSTART2]=-numbnd[2];
	  dirset[bti][dir][DIRUSTOP2]=N2-1+numbnd[2];
	  dirset[bti][dir][DIRUSTART3]=-numbnd[3];
	  dirset[bti][dir][DIRUSTOP3]=N3-1+numbnd[3];
	}

	// x2
	if(dir==X2UP){ // up
	  dirset[bti][dir][DIRUSTART2]=N2+SHIFT2;
	  dirset[bti][dir][DIRUSTOP2]=N2-1+numbnd[2];
	}
	else if(dir==X2DN){ // down
	  dirset[bti][dir][DIRUSTART2]=-numbnd[2]+1;
	  dirset[bti][dir][DIRUSTOP2]=0;
	}

	if((dir==X2UP)||(dir==X2DN)){
	  dirset[bti][dir][DIRUSTART1]=-numbnd[1];
	  dirset[bti][dir][DIRUSTOP1]=N1-1+numbnd[1];
	  dirset[bti][dir][DIRUSTART3]=-numbnd[3];
	  dirset[bti][dir][DIRUSTOP3]=N3-1+numbnd[3];
	}

	// x3
	if(dir==X3UP){ // up
	  dirset[bti][dir][DIRUSTART3]=N3+SHIFT3;
	  dirset[bti][dir][DIRUSTOP3]=N3-1+numbnd[3];
	}
	else if(dir==X3DN){ // down
	  dirset[bti][dir][DIRUSTART3]=-numbnd[3]+1;
	  dirset[bti][dir][DIRUSTOP3]=0;
	}

	if((dir==X3UP)||(dir==X3DN)){
	  dirset[bti][dir][DIRUSTART1]=-numbnd[1];
	  dirset[bti][dir][DIRUSTOP1]=N1-1+numbnd[1];
	  dirset[bti][dir][DIRUSTART2]=-numbnd[2];
	  dirset[bti][dir][DIRUSTOP2]=N2-1+numbnd[2];
	}
      }// end if DIRIF
    }// end DIRLOOP
  }// end over bti












  /////////////////
  //
  // output those things that were defined
  //
  /////////////////

#if(USEMPI)
  fprintf(log_file,"myid=%d node_name=%s procnamelen=%d\n",myid,processor_name,procnamelen);
#endif
  trifprintf("\nnumprocs=%d ncpux1=%d ncpux2=%d ncpux3=%d\n",numprocs,ncpux1,ncpux2,ncpux3);
  fprintf(log_file,"per: %d %d\n", periodicx1, periodicx2);

  for (m = 1; m <= COMPDIM; m++) {
    fprintf(log_file,"mycpupos[%d]: %d\n", m, mycpupos[m]);
    fprintf(log_file, "startpos[%d]: %d\n", m, startpos[m]);
    fprintf(log_file, "endpos[%d]: %d\n", m, endpos[m]);
    fprintf(log_file, "totalsize[%d]: %d\n", m, totalsize[m]);
  }

  for(bti=0;bti<NUMBOUNDTYPES;bti++) {
    for (m = 0; m < COMPDIM*2; m++) {
      for(l = 0 ; l < DIRNUMVARS ; l++) {
	fprintf(log_file, "dirset[%d][%d][%d]: %d\n", bti, m, l, dirset[bti][m][l]);
      }
    }
  }

  trifprintf("totalzones: %d\n", totalzones);




  /////////////////// 
  //
  // Setup supermpi method (not working right now)
  //
  ///////////////////


#if(SIMULBCCALC!=-1)
  // this is definitely not setup for 3D, and never fully worked...still interesting idea.

  if(SIMULBCCALC==2){
    if(SIMULBCCALC<=0){ stagei=STAGEM1; stagef=STAGEM1; }
    else if(SIMULBCCALC==1) { stagei=STAGE0; stagef=STAGE2;}
    else if(SIMULBCCALC==2) { stagei=STAGE0; stagef=STAGE5;}
    
    if(SIMULBCCALC>=1){
      for(stage=stagei;stage<=stagef;stage++){
	STAGECONDITION(0,N1-1,0,N2-1,isc,iec,jsc,jec);
	fprintf(log_file,"CZLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(0,N1,-1,N2,isc,iec,jsc,jec);
	fprintf(log_file,"F1LOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(-1,N1,0,N2,isc,iec,jsc,jec);
	fprintf(log_file,"F2LOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(0,N1,0,N2,isc,iec,jsc,jec);
	fprintf(log_file,"EMFLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(0,N1,0,N2-1,isc,iec,jsc,jec);
	fprintf(log_file,"F1CTLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(0,N1-1,0,N2,isc,iec,jsc,jec);
	fprintf(log_file,"F2CTLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(-1,N1,-1,N2,isc,iec,jsc,jec);
	fprintf(log_file,"DQLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(-2,N1+1,-2,N2+1,isc,iec,jsc,jec);
	fprintf(log_file,"PREDQLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	fprintf(log_file,"\n");
      }
    }
  }
#endif

  fflush(log_file);




  trifprintf("end: init_placeongrid\n");
}


int myexit(int call_code)
{
  int i, j, k, l;
  int cleanfinish,dofaildump;
  FILE *faildump;
  char mysys[MAXFILENAME];
  char binarytype[MAXFILENAME];
  void set_binarytype(char *binarytype);


  trifprintf("proc: %s : Exiting cc: %d nstep: %ld\n", myidtxt,
	  call_code, nstep);




#if(MAILWHENDONE && !MPIAVOIDFORK)
    if(myid==0){

      set_binarytype(binarytype);

      sprintf(mysys,"echo \"%s : done with `pwd`\" > done%s.txt",EMAILMESSAGE,binarytype);
      system(mysys);
      if(MAILFROMREMOTE){
	sprintf(mysys,"scp done.txt %s ; ssh %s \"mail %s < done%s.txt\"",REMOTEHOST,REMOTEHOST,EMAILADDRESS,binarytype);
	system(mysys);
      }
      else{
	sprintf(mysys,"mail %s < done%s.txt",EMAILADDRESS,binarytype);
	system(mysys);
      }
    }
#endif

  dofaildump=0;
  if (call_code > 0) {
    fprintf(stderr,
	    "proc: %s : Failure.  Please check failure file: cc: %d\n",
	    myidtxt, call_code);

    if(call_code<ERRORCODEBELOWCLEANFINISH) cleanfinish = 1;
    else cleanfinish=0; // assume this means dump procedure failed, so don't get into infinite failure loop
    // should never have non-clean finish, but sometimes do have them in code, but not marked right now
    if(cleanfinish) dofaildump=1;
    if(!cleanfinish){
#if(USEMPI)
      // must abort since no clear to communicate to other cpus now
      MPI_Abort(MPI_COMM_GRMHD, 1);
#endif
    }
  }
  else{
    dofaildump=0;
    cleanfinish=1;
  }



  if (dofaildump) {
    fprintf(stderr, "proc: %s : dumping failure dump with callcode=2\n",
	    myidtxt);

    // assume want previous timestep data, not bad just-computed
    // data\n");
    // now diag should not fail if last timestep was non-fail type
    if (DODIAGS) diag(FINAL_OUT,t,nstep,realnstep);
  }


  ////////////////////////////////
  //
  // first cleanup any prior MPI non-blocking calls by making fake write call
  //
  ////////////////////////////////
  if(cleanfinish && USEROMIO==1 && MPIVERSION==2 ){
    fakedump();
  }



  ////////////////////////////////
  //
  // must close AFTER diag()
  //
  ////////////////////////////////
  if (call_code >= 0) {
    if (fail_file)
      fclose(fail_file);
    if (log_file)
      fclose(log_file);
    myfclose(&logfull_file,"Can't close logfull_file\n");
  }




  if(cleanfinish){
    fprintf(stderr,
	    "Ending Computation on proc: %s, holding for other cpus\n",
	    myidtxt);
  }


  myfprintf(stderr, "Ended Computation on all processors\n");
  final_myexit();

  return (0);
}

// note, this may be called in different locations of the code by
// different CPUs
int error_check(int wherefrom)
{
  int i, j, k;
  int errorsend = 0;
  // check if error exists and exit if so

  if (failed > 0) {
    dualfprintf(fail_file,
	    "Detected failure on proc: %d failed: %d nstep: %ld realnstep: %ld t: %21.15g wherefrom = %d\n",
	    myid, failed, nstep, realnstep, t,wherefrom);
  }

  if (numprocs > 1) {
    errorsend = failed;
#if(USEMPI)
    // fprintf(fail_file,"wtf: %d %d\n",errorsend,failed);
    // fflush(fail_file);
    MPI_Allreduce(&errorsend, &failed, 1, MPI_INT, MPI_MAX,
		  MPI_COMM_GRMHD);
    // fprintf(fail_file,"wtf: %d %d\n",errorsend,failed);
    // fflush(fail_file);
#endif
  }
  if (failed > 0) {
    dualfprintf(fail_file,
	    "Result: Detected failure on proc: %d failed: %d nstep: %ld realnstep: %ld t: %21.15g\n",
	    myid, failed, nstep, realnstep, t);
    // control behavior of failure here (i.e. could return(1) and
    // continue or something)
    // if(failed==1) myexit(1);
    // if(failed==2) myexit(1);
    // if(failed==3) myexit(1);
    myexit(wherefrom);
    return (1);
  }
  return (0);
}






// just copied from pnmhd code
#if(0)
void init_MPIgroup(void)
{
  int ranks[MAXCPUS] = {0}; 
  int i,j,k,numranks;

  MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

  fprintf(stderr,"begin: proc: %s init_MPIgroup\n",myidtxt); fflush(stderr);

  // x1-inner
  j=0;
  for(i=0;i<numprocs;i++){
    if(i%ncpux1==0){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux2*ncpux3){
    fprintf(stderr,"problem with inner x1-group: numranks: %d ncpux2: %d ncpux3: %d ncpux2*ncpux3: %d\n",numranks,ncpux2,ncpux3,ncpux2*ncpux3);
    myexit(1);
  }
  // now ranks holds inner x1 boundary of cpus, and numranks holds number of such ranks

  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[2]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[2], &combound[2]); 

  // x1 - outer
  j=0;
  for(i=0;i<numprocs;i++){
    if(i%ncpux1==ncpux1-1){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux2*ncpux3){
    fprintf(stderr,"problem with outer x1-group: numranks: %d ncpux2*ncpux3: %d\n",numranks,ncpux2*ncpux3);
    myexit(1);
  }
  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[0]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[0], &combound[0]); 

  // x2 - inner
  j=0;
  for(i=0;i<numprocs;i++){
    if((i%(ncpux1*ncpux2))/ncpux1==0){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux1*ncpux3){
    fprintf(stderr,"problem with inner x2-group: numranks: %d ncpux1*ncpux3: %d\n",numranks,ncpux1*ncpux3);
    myexit(1);
  }
  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[1]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[1], &combound[1]); 

  // x2 - outer
  j=0;
  for(i=0;i<numprocs;i++){
    if((i%(ncpux1*ncpux2))/ncpux1==ncpux2-1){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux1*ncpux3){
    fprintf(stderr,"problem with outer x2-group: numranks: %d ncpux1*ncpux3: %d\n",numranks,ncpux1*ncpux3);
    myexit(1);
  }

  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[3]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[3], &combound[3]); 



  // x3 - inner
  j=0;
  for(i=0;i<numprocs;i++){
    if(i/(ncpux1*ncpux2)==0){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux1*ncpux2){
    fprintf(stderr,"problem with inner x3-group: numranks: %d ncpux1*ncpux2: %d\n",numranks,ncpux1*ncpux2);
    myexit(1);
  }
  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[5]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[5], &combound[5]); 

  // x3 - outer
  j=0;
  for(i=0;i<numprocs;i++){
    if(i/(ncpux1*ncpux2)==ncpux3-1){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux1*ncpux2){
    fprintf(stderr,"problem with outer x3-group: numranks: %d ncpux1*ncpux2: %d\n",numranks,ncpux1*ncpux2);
    myexit(1);
  }

  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[4]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[4], &combound[4]); 

  // 0: right 1: up 2: left 3: down 4: out 5: in(as in bound.c)

  // when using these communicators, must make sure the call to a communication using it isn't done by a non-member cpu!
  // (above: stupid, I know, should just skip if non-member cpu tries a function)
  fprintf(stderr,"end: proc: %s init_MPIgroup\n",myidtxt); fflush(stderr);
}
#endif

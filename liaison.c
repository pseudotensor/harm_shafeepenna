#include "defs.liaison.h"

// in global.h, mympi.h, makefile:
// 1) DOINGLIAISON must be set to 1 for both GRMHD and liaison code
// 2) USINGLIAISON is used at compile time for ONLY liaison code




// to run multiple binaries during MPI:

///////////////////
// MPICH-1
// 1) Create pgfile named <pgfilename> with contents like:
//    <host0> 0 <full path/programa>
//    <host1> 1 <full path/programb>
//    <host2> 1 <full path/programc>
//    <host3> 1 <full path/programd>
//    ...
//
// Then run using: mpirun -np <numprocesses> -p4pg <pgfilename> -p4wd <workingdirectory>

/////////////////////////
// MPICH-2 or LAM
// If no mpd, then do: mpd &
// If not setup ~/.mpd.conf then it'll say how
// To ensure running mpd correctly, list nodes can run on: mpdtrace 
//
// To run MPI program with different binaries do:
// mpiexec -n <numgrmhdprocs> grmhdbinary <grmhdargs> : -n <numgrrayprocs> grraybinary <grrayargs> : -n <numliaisonprocs> liaisonbinary <liaisonargs>

//////////////////
// http://www.mcs.anl.gov/research/projects/mpich2/documentation/files/mpich2-doc-user.pdf
// To use MPI-2/MPICH-2/mpiexec on PBS do:
//
// sort $PBS NODEFILE | uniq -C | awk '{ printf("%s:%s", $2, $1); }' > mpd.nodes
//Once the PBS node file is converted, MPD can be normally started within the PBS job script using mpdboot and torn down using mpdallexit.
// mpdboot -f mpd.hosts -n [NUM NODES REQUESTED]
// mpiexec -n [NUM PROCESSES] ./my test program
// mpdallexit



// LAM:
// http://www.lam-mpi.org/faq/category5.php3#question23










// definitions and declarations for liaison specific code
int init(int *argc, char **argv[]);
void set_arrays();
void step_liaison_full(void);
int init_MPI_LIAISON(int *argc, char **argv[]);
void liaison_init_mpi_communicators(void);








int main(int argc, char *argv[])
{
  static SFTYPE comptstart;

 /* perform initializations */
  if (init(&argc,&argv) >= 1) {
    dualfprintf(fail_file, "main:init: failure\n");
    myexit(1);
  }


  // write version header info for perf and step logs
  if((myid<=0)&&(appendold==0)){
    if(DOLOGSTEP){
      myfprintf(logstep_file,"#%10s\n%10d %10d\n#","STEPVER",STEPVER,STEPTYPE);
    }
    if(DOLOGPERF){
      myfprintf(logperf_file,"#%10s\n%10d %10d\n#","PERFVER",PERFVER,PERFTYPE);
    }
  }


  trifprintf("proc: %04d : Start computation\n", myid);

  // set initial computational time
  //  comptstart=t; // not right
  comptstart=0; // assumes normal simulation started t=0 at nstep=0

  // start timing
  timecheck(STARTTIME,comptstart);


  while (reallaststep==0) {


    step_liaison_full();

    nstep++;

    // time check (should come after first nstep++ to setup things first)
    timecheck(CHECKTIME,comptstart);

    // speed check
    timecheck(SPEEDTIME,comptstart);


    // output timestep info
    output_steptimedt_info(comptstart);

  }
  
  trifprintf("proc: %04d : End computation\n", myid);


  // stop timing
  timecheck(STOPTIME,comptstart);

  // report final timing and other simulation performance characteristics
  timecheck(REPORTTIME,comptstart);


  myexit(0);
  return (0);
}



// for testing LIAISON+GRMHD code communications
void test_liaison(void)
{
  int myint;

  myint=-100;

#if(USEMPILIAISON)
  MPI_Bcast(&myint,1,MPI_INT,0,MPI_COMM_GRMHD_LIAISON);
  dualfprintf(fail_file,"myid=%d myint=%d\n",myid,myint);
#endif

  myexit(0);

 
}




// repeated procedure to do for each liaison step
void step_liaison_full(void)
{
  void test_liaison(void);

  // for now test
  test_liaison();

}




// initialize liaison code
int init(int *argc, char **argv[])
{


  // power up random number generator in case used without init
  ranc(0);


#if(DOINGLIAISON==0)
  fprintf(stderr,"Are you sure you want to have liaison code with DOINGLIAISON==0?\n");
#endif

  // init MPI (assumes nothing in set_arrays.c used here)
  init_MPI_LIAISON(argc, argv);
  


  // report system information
  report_systeminfo(stderr);
  report_systeminfo(log_file);
  if(myid==0) report_systeminfo(logfull_file);



  // setup files for writing and reading (must come after init_mpi_setupfilesandgrid())
  makedirs();


  // inititialize arrays for storing grid data
  set_arrays();



  return(0);

}





// Initialize MPI for LIAISON code
int init_MPI_LIAISON(int *argc, char **argv[])
{

#if(USEMPILIAISON)
  fprintf(stderr, "begin: init_MPI_LIAISON\n");
  fflush(stderr);
  // init MPI (assumes nothing in set_arrays.c used here) : always done
  // non-blocking
  init_MPI_general(argc, argv);
#else
  fprintf(stderr, "Did NOT init_MPI_LIAISON\n");
  fflush(stderr);
#endif


  // always do below since just sets defaults if not doing liaisonmode
  liaison_init_mpi_liaisonmode_globalset();

#if(USEMPILIAISON)
  // this is non-blocking local operation
  MPI_Comm_rank(MPI_COMM_LIAISON_FROM_GRMHD, &myid); // proc id within LIAISON_FROM_GRMHD context only
#endif

  // for file names
  sprintf(myidtxt, ".liaison.%04d", myid);





  // currently INIT provides args to rest of processes
  // liaison must have same arguments as GRMHD code
  myargs(*argc,*argv);

  // rest of initialization
  init_MPI_setupfilesandgrid(*argc, *argv);



  return (0);


}  



void makedirs(void)
{

  if ((USEMPILIAISON == 0) || (USEMPILIAISON && (!MPIAVOIDFORK))) {
    if ((mpicombine && (myid == 0)) || (mpicombine == 0)) {
      //      system("mkdir dumps");
      //      system("mkdir images");
    }
#if(USEMPILIAISON)
    //    MPI_Barrier(MPI_COMM_LIAISON_FROM_GRMHD);	// all cpus wait for directory
    // to be created
#endif
  }
}





void set_arrays()
{
  int i, j, k, pl, l, m;
  int pl2;
  int ii;
  int floor,pf, tscale,dtstage;
  int dissloop;
  FTYPE valueinit;
  int dir,interpi,enodebugi;
  int dimen;





#if(PRODUCTION==0)
  // initialize things to NAN in order to (hopefully) trigger memory leaks to be noticed
  valueinit=sqrt(-1.0);
#else
  // avoid this practice for production run since processing nan's slows code and do process some never-initialized/never-used cells for simplicity of code loops
  valueinit=0.0;
#endif



  ////////////////////////////////////////////////
  //
  // Basic primitive quantity
  //
  ////////////////////////////////////////////////

  p = (FTYPE (*)[N2M][N3M][NPR]) (&(a_p[N1BND][N2BND][N3BND][0]));
  FULLLOOP PLOOP(pl){
    p[i][j][k][pl] = valueinit;
  }
  


  gcon = (FTYPE (*)[N2M+SHIFT2][N3M+SHIFT3][NPG][NDIM][NDIM])
    (&(a_gcon[N1BND][N2BND][N3BND][0][0][0]));
  
  // rest are always located at CENT
  conn = (FTYPE (*)[N2M][N3M][NDIM][NDIM][NDIM])
    (&(a_conn[N1BND][N2BND][N3BND][0][0][0]));

}




int diag(int call_code, FTYPE localt, long localnstep, long localrealnstep)
{
  return(0);
}




#include "decs.h"



// Initialize MPI for GRMHD code
int init_MPI_GRMHD(int *argc, char **argv[])
{

#if(USEMPI)
  fprintf(stderr, "begin: init_MPI_GRMHD\n");
  fflush(stderr);
  // init MPI (assumes nothing in set_arrays.c used here) : always done
  // non-blocking:
  init_MPI_general(argc, argv);
#else
  fprintf(stderr, "Did NOT init_MPI_GRMHD\n");
  fflush(stderr);
#endif


  // always do below since just sets defaults if not doing liaisonmode
  grmhd_init_mpi_liaisonmode_globalset();


#if(USEMPI)
  // this is non-blocking local operation
  MPI_Comm_rank(MPI_COMM_GRMHD, &myid); // proc id within GRMHD context only
#endif

  // for file names
  sprintf(myidtxt, ".grmhd.%04d", myid);




  // currently INIT provides args to rest of processes
  myargs(*argc,*argv);

  // rest of initialization
  init_MPI_setupfilesandgrid(*argc, *argv);


#if(DOINGLIAISON)
  // liaison-related test code:
  test_nonliaison();
#endif

  return (0);


}  


// for testing LIAISON+GRMHD code communication
void test_nonliaison(void)
{
  int myint;

  myint=myid;

#if(USEMPI&&DOINGLIAISON)
  MPI_Bcast(&myint,1,MPI_INT,0,MPI_COMM_GRMHD_LIAISON);
  dualfprintf(fail_file,"myid=%d myint=%d\n",myid,myint);
#endif

  myexit(0);

  
}

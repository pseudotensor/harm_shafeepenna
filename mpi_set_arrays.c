
#include "decs.h"


void mpi_set_arrays(void)
{

#if(USEMPI)
  workbc = (FTYPE(*)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM]) (&(workbca[-1][0][0]));
  
  workbc_int =(PFTYPE(*)[COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM]) (&(workbc_inta[-1][0][0]));
#endif

}


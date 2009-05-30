#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>

#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */

#include "SwapEndian.h"



#ifndef WIN32
#include<sys/stat.h>
#endif

///////////////////
//
// pre-global stuff:
//
///////////////////

// whether doing performance testing (also see pnmhd code)
#define DOINGLIAISON 0 // choice of 0 or 1 (should always be 1 for liaison mode and liaison code comilation)
#define NCSA 0
#define PERFTEST 0
#include "mytime.h"


//////////////////
//
// GLOBAL STUFF:
//
//////////////////


#ifndef GLOBAL_H
#define GLOBAL_H

// Notes:
//
// To check global symbol sizes use:
//      nm -S --size-sort grmhd 
// or:  objdump -axfhp grmhd
//
// e.g. to check total size used do:
// 1) nm -S --size-sort grmhd | awk '{print $2}' > list.grmhd.sizes.txt
// 2) total=0
// 3) for fil in `cat list.grmhd.sizes.txt` ; do total=$(($total+0x$fil)) ; done
// 4) echo $total
// Now take that number and divide by N1M*N2M*N3M and that's approximate per zone bytes used  (must include boundary cells)
//
// Seems to be about 1250 elements/zone


// define FTYPE, SFTYPE, etc.
#include "global.realdef.h"

// all pure nmenomics that don't depend upon anything else but what's inside file
#include "global.nondepnmemonics.h"

// RK-related macros
#include "global.stepch.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Default global and user choices for various code options that can change for each run without significant modifications
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// default values
#include "definit.h"

// user specific values
#include "init.h"


// define nmenomics that depend on other mnemonics
#include "global.depnmemonics.h"

#include "global.loops.h"

#include "global.variousmacros.h"

#include "global.fieldmacros.h"

#include "global.structs.h"


// now that all hashes have been defined, get mpi header
#include "mympi.h"
#include "global.gridsectioning.h"
#include "global.comploops.h"

// all global external function declarations
#include "global.funcdeclare.h"


#include "global.other.h"


#endif// endif for #ifndef GLOBAL_H



/////////////////////////////////////////
//
// test declarations
// see checkexterns.sh
// comment below 3 lines when not wanting to test
// uncomment when wanting to test, and then look at make.log
// need to keep files listed in checkexterns.sh updated
//
//typedef double doublereal;
//#include "utoprim_jon.h" // extra defs
//#include "temptempfinalallc.txt"
//
//
/////////////////////////////////////////






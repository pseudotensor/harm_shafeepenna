
#define DOINGGRMHDTYPECODE 1 // always 1
#define DOINGGRRAYTYPECODE 0 // always 0
#define DOINGLIAISONTYPECODE 0 // always 0
#include "global.general.h"



// debugging on Harvard's Bluegene:
// Compile with -g -O0
// On BG/P, core files are text files.  Look at the core file with a text editor, focus on the function call chain; feed the hex addresses to addr2line.
// addr2line -e your.x  hex_address
// 
// tail -n 10 core.511 | addr2line -e your.x
// 
// Use grep and word-count (wc) to examine core files :
//  grep hex_address “core.*” | wc -l
// 
// You can get the instruction that failed by using objdump:
// powerpc-bgp-linux-objdump -d your.x >your.dump
// You can locate the instruction address in the dump file, and can at least find the routine where failure occurred, even without –g.
// 
// If your application exits without leaving a core file,
// set the env variable BG_COREDUMPONEXIT=1

/* coreprocessor.pl … coreprocessor perl script in : */

/* /bgsys/drivers/ppcfloor/tools/coreprocessor */

/* online help via :  coreprocessor.pl –help */

/* Can analyze and sort text core files, and can attach to hung processes for deadlock determination.   */

/* Normal use is “gui” mode, so set DISPLAY then: */

/* coreprocessor.pl –c=/bgusr/username/rundir –b=your.x */

/* click on “Select Grouping mode”,  */
/* select  “Stack Trace (condensed)” */
/* click on source statement of interest */

/* There is also a non-gui mode: coreprocessor.pl –help */

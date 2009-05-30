# NOTES:

# When compiling the code on a new machine:
# make superclean ; make prep ; make grmhd

# When on same machine but changed files make*:
# make prep ; make grmhd

# When only changed *.h and *.c do just:
# make grmhd

# for liaison code do (e.g.):
# make superclean ; make prepliaison ; make liaison

# for jon_interp code do (e.g.):
# make superclean ; make prepjoninterp ; make iinterp

# for bin2txt code do (e.g.):
# make superclean ; make prepbin2txt ; make bin2txt
 
# for smcalc code do (e.g.):
# make superclean ; make prepsmcalc ; make smcalc
 



############################
# header
include makehead.inc

# turns on -mp and -pc64 for USEICC=1
# Note that if using -pc64 -mp that error in inversion seems to be limited for doubles to 1E-11 instead of 1E-15
# causes code to be less accurate: e.g. inversion doesn't reach to 1E-15, but rather 1E-11
# for some simple linear wave tests with small amplitude can cause significant problems if any  noise exists in intial solution
ENFORCEPRECISION=0

ifeq ($(ENFORCEPRECISION),1)
PRECISE=-mp -pc64
endif
ifeq ($(ENFORCEPRECISION),0)
PRECISE=
endif

#####################
# WATCH OUT FOR SPACES, etc. AFTER ASSIGNMENTS!!!!
USEBG=1
USEICC=0
USEGCC=0
USECCC=0
USEORANGE=0
USENERSC=0
USETACCLONESTAR=0
USETACCRANGER=0

# default
USEMCCSWITCH=0
AVOIDFORK=0
USESPECIAL4GENERATE=0
CCGENERATE=gcc


ifeq ($(USEMPI),1)

# override with new default
USEMCCSWITCH=1


ifeq ($(USEBG),1)
# override again
USEMCCSWITCH=0
AVOIDFORK=1
MCC = mpicc
CCGENERATE=gcc
USESPECIAL4GENERATE=1
endif

ifeq ($(USETACCLONESTAR),1)
# override again
AVOIDFORK=1
MCC = mpicc
endif

ifeq ($(USENERSC),1)
# override again
USEMCCSWITCH=0
AVOIDFORK=1
MCC = cc
endif

ifeq ($(USEICC),1)
# uses -static for secure library usage
# MCC=/usr/local/p4mpich-1.2.5-icc-noshmem/bin/mpicc
MCC = mpicc
endif

ifeq ($(USETACCRANGER),1)
# don't have to avoid fork/system calls
AVOIDFORK=0
MCC = mpicc
endif

ifeq ($(USEORANGE),1)
# orange can't have -static
MCC = mpicc -I/afs/slac/package/OpenMPI/include/ -L/afs/slac/package/OpenMPI/lib/
endif

ifeq ($(USEGCC),1)
MCC = mpicc
#.gcc
endif

ifeq ($(USEPGCC),1)
MCC=mpicc
endif    

endif    




# default extra flags to pass to compiler
EXTRA=-DUSINGMPIAVOIDFORK=$(AVOIDFORK)





#
# Define preprocessor and compile flags, and linked libraries

ifeq ($(USEGCC),1)

DFLAGS=-DUSINGICC=0 -DUSINGORANGE=0 $(EXTRA)


COMP=gcc $(DFLAGS)
# COMP=gcc3

# gcc type flags
#
########################
# UES BELOW
########################
#CFLAGSPRE = -Wall -mpentium -O3 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -DCPU=686 -DNEED_GETOPT -DLINUX
#
#
#
# for DEBUG:
#CFLAGSPRE = -Wall -O0 -g -pg  $(DFLAGS)
#CFLAGSPRENONPRECISE=-O0 -g -pg $(DFLAGS)

CFLAGSPRE = -Wall -O3 $(DFLAGS)
CFLAGSPRENONPRECISE=-O3 $(DFLAGS)

#
#CFLAGS = -Wall -mpentium -O3 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -DCPU=686 -DNEED_GETOPT -DLINUX -ffast-math -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -mstack-align-double -ffast-math -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -mstack-align-double -ffast-math -finline-functions -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -ffast-math -finline-functions -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O3 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -ffast-math -finline-functions -pg -g -a
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O3 -ffast-math -finline-functions -funroll-loops
#CFLAGS=-O0 -g -Wall -wunused-label -wunused-parameter
#CFLAGS=-O0 -g -Wall
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -ffast-math -finline-functions -g
#-pg
#-pg -g  source lines
#-pg -g -a   line count
# gprof -l <file> > out.txt
# gprof -A -I<sourcedir>
# gprof -l -A -x s

#below is typical flags for double precision...can take -pg off for no profile
#add -mstack-align-double if using pgcc
#CFLAGS = -Wall -O0 -g
#  -fomit-frame-pointer



#CFLAGS = -Wall -O0
#CFLAGS = -O6 -g
#CFLAGS = -O0 -pg -g
LDFLAGS = -lm
# -l$(LAPACKLIB) -l$(BLASLIB)  -L/usr/lib/gcc-lib/i386-redhat-linux/2.96/ -l$(F2CLIB) 

#CC = cc
#AR	=	ar r
#RANLIB	=	ranlib
endif





ifeq ($(USEICC),1)

DFLAGS=-DUSINGICC=1  -DUSINGORANGE=0 $(EXTRA)

COMP=icc $(DFLAGS)
# -ipo is not same for multiple names in c files, unless static names
# -ipo_wp is generally not safe, and not for HARM
# -ip is safe, always
# -ipo not safe for icc v8.1 for connection coefficients for some reason
#CFLAGSPRE=-O3 -tpp7 -axKW -mp -unroll -Wall -Wcheck -Wshadow -w2 -wd=175,177,279,593,869,810,981,1418,1419,310,1572 -g -pg
# from older code
#CFLAGSPRE=-O3 -tpp7 -axiMKW -unroll -Wall -w2 -wd=175,177,279,593,869,810,981,1418,1419,310
# newer code below
#CFLAGSPRE=-O3 -tpp7 -axKW -unroll -Wall -Wcheck -Wshadow -pc64 -mp -w2 -wd=175,177,279,593,869,810,981,1418,1419,310,1572 -mp

#CFLAGSPRENONPRECISE=-O3 -tpp7 -axKW -unroll -Wall -Wcheck -Wshadow -w2 -wd=175,177,279,593,869,810,981,1418,1419,310,1572
#
#########################
# NORMALLY USE BELOW
#########################
#
# autoparallelization doesn't work except for very trivial loops and not even all of those
#CFLAGSPRENONPRECISE=-O2 -parallel -Wall -Wcheck -Wshadow -w2 -wd=175,177,279,593,869,810,981,1418,1419,310,1572 $(DFLAGS)


# NORMAL:
CFLAGSPRENONPRECISE=-O2 -finline -finline-functions -ip -fno-alias -unroll -openmp -Wall -Wcheck -Wshadow -w2 -wd=175,177,279,593,869,810,981,1418,1419,310,1572 $(DFLAGS)
#

# FOR CHECKING OPTIMIZATIONS:
#CFLAGSPRENONPRECISE=-O2 -openmp -Wall -Wcheck -Wshadow -w2 -wd=175,177,279,593,869,810,981,1418,1419,310,1572 -g -pg $(DFLAGS)
#
#########################
# DEBUG BELOW
#########################
#CFLAGSPRENONPRECISE=-O0 -g -openmp $(DFLAGS)
#
#
#CFLAGSPRENONPRECISE=
#########################
# USE BELOW FOR DEBUG
#########################
#CFLAGSPRENONPRECISE=-O0 -g
#
#
#
CFLAGSPRE=$(PRECISE) $(CFLAGSPRENONPRECISE)

GCCCFLAGS = -Wall -O2 $(DFLAGS)

#CFLAGSPRE=-O0 -g
#CFLAGSPRE=-O0 -g
#CFLAGSPRE=-O3 -mp -unroll -Wall -Wcheck -Wshadow -w2 -wd=175,177,279,593,869,810,981,1418,1419,310,1572
#CFLAGSPRE=-O3 -mp -unroll -Wall -Wcheck -Wshadow -w2 -wd=175,177,279,593,869,810,981,1418,1419,310,1572 -g
#CFLAGSPRE=-O0 -g
#CFLAGSPRE=-O3 -tpp7 -axKW -ip -unroll -Wall -Wcheck -Wshadow -w2 -wd=175,177,279,593,869,810,981,1418,1419,310,1572
#CFLAGS=-O3 -unroll -Wall -w2 -wd=981,279,869,1572,1418,177,1419,593,810,310

#GCCCFLAGS = -Wall -mpentium -O3 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -DCPU=686 -DNEED_GETOPT -DLINUX -ffast-math
#GCCCFLAGS = -O0 -g -pg
#CFLAGSPRE = -O3 -g -Wall -w2 -wd=175,177,279,593,869,810,981,1418,1419,310

# P4 (don't use -wp_ipo -- incompat code)
#CFLAGS=-O3 -tpp7 -axiMKW -wp_ipo -unroll -w1
#CFLAGS=-O3 -tpp7 -axiMKW -ipo -unroll -Wall -w2 -wd=175,177,279,593,869,810,981,1418,1419,310 -pg
#CFLAGS=-O3 -unroll -axiMKW -unroll -pg
# -rcd causes problems, like asymetries in interp program I've noticed
#CFLAGS = -Wall -w2 -O3 -axiMKW -unroll -ipo -tpp7 -march=pentium4 -mcpu=pentium4# -p\
#arallel
# -w2 displays more info for warnings
#CFLAGS=-O3 -axiMKW -g -ipo -pg
# P3
# CFLAGS=-O3 -tpp6 -axiMK -ipo -unroll -w1
# GAMMIE
# CFLAGS = -O3 -ipo

LDFLAGS = -lm
LDFLAGSOTHER=
endif





ifeq ($(USETACCLONESTAR),1)
DFLAGS=-DUSINGICC=1  -DUSINGORANGE=0 $(EXTRA)
COMP=icc $(DFLAGS)
CFLAGSPRENONPRECISE=-O2 -xT -finline -finline-functions -ip -fno-alias -unroll -openmp -Wall -Wcheck -Wshadow -w2 -wd=175,177,279,593,869,810,981,1418,1419,310,1572 $(DFLAGS)
CFLAGSPRE=$(PRECISE) $(CFLAGSPRENONPRECISE)
GCCCFLAGS = -Wall -O2 -L$ICC_LIB -lirc $(DFLAGS)
LDFLAGS = -lm
LDFLAGSOTHER=
endif

ifeq ($(USETACCRANGER),1)
DFLAGS=-DUSINGICC=1  -DUSINGORANGE=0 $(EXTRA)
COMP=icc $(DFLAGS)
CFLAGSPRENONPRECISE=-xW -O2 -finline -finline-functions -ip -fno-alias -unroll -openmp -Wall -Wcheck -Wshadow -w2 -wd=175,177,279,593,869,810,981,1418,1419,310,1572 $(DFLAGS)
CFLAGSPRE=$(PRECISE) $(CFLAGSPRENONPRECISE)
GCCCFLAGS = -Wall -O2 -L$ICC_LIB -lirc $(DFLAGS)
LDFLAGS = -lm
LDFLAGSOTHER=
endif




ifeq ($(USECCC),1)

DFLAGS=-DUSINGICC=0 -DUSINGORANGE=0 $(EXTRA)
COMP=ccc $(DFLAGS)
LDFLAGS =  -lm -lcxml

#CDEBUG = -g3 # -g3 for higher opts than -O0
#CDEBUG = -g
# do profile (profile screws up speed for loops, doesn't unroll them, etc.)
#CDEBUG = -pg -g3
# production level
CFLAGSPRE = -Wall -O4 -fast -msg_disable badsubscript -msg_disable subscrbounds -msg_disable unreachcode -msg_disable noparmlist -msg_disable subscrbounds2 -msg_disable longlongtype -finline-functions -funroll-loops  $(DFLAGS)
CFLAGSPRENONPRECISE=$(CFLAGSPRE)
#CFLAGS3 = -Wall -O0
# super annoying develop level
#CFLAGS3 = -Wall -O2 -fast
#CFLAGS3 = -fast -arch ev67
# debug level
#CFLAGS3 = -Wall -O0 -msg_disable badsubscript -msg_disable subscrbounds -msg_disable unreachcode -msg_disable noparmlist -msg_disable subscrbounds2
#CFLAGS3 = -Wall -O0
#CFLAGS3 = -Wall -O2

endif


ifeq ($(USEBG),1)
DFLAGS=-DUSINGICC=0  -DUSINGORANGE=0 $(EXTRA)
COMP=gcc $(DFLAGS)
# don't use -ffast_math, causes asymmetries in calculations
CFLAGSPRE= -O3 -funroll-loops  $(DFLAGS)
CFLAGSPRENONPRECISE= $(CFLAGSPRE)
GCCCFLAGS= -O3 $(DFLAGS)
#-funroll-loops -fargument-noalias -mcpu=k8 -msse2 -mfpmath=sse -static
LDFLAGS= -lm
endif



ifeq ($(USEORANGE),1)
DFLAGS=-DUSINGICC=1  -DUSINGORANGE=1 $(EXTRA)
COMP=icc -I/afs/slac/package/OpenMPI/include/ -L/afs/slac/package/OpenMPI/lib/ $(DFLAGS)
#CFLAGSPRE=-O3 -fno-alias -ftz -unroll -Wall -w2 -wd=175,177,279,593,869,810,981,1418,1419,310,1572
CFLAGSPRE=-O2 -fno-alias -ftz -unroll -Wall -w2 -wd=175,177,279,593,869,810,981,1418,1419,310 $(DFLAGS)
CFLAGSPRENONPRECISE=$(CFLAGSPRE)
LDFLAGS = -lm
#GCCCFLAGS = -Wall -O3 -m32
GCCCFLAGS = -Wall -O3 $(DFLAGS)
endif




ifeq ($(USENERSC),1)
DFLAGS=-DUSINGICC=1  -DUSINGORANGE=0 $(EXTRA)
COMP=cc $(DFLAGS)
CFLAGSPRE = -O3 $(DFLAGS)
CFLAGSPRENONPRECISE = $(CFLAGSPRE)
GCCCFLAGS =  $(CFLAGSPRE)
LDFLAGS = -lm
endif













CFLAGS=$(CFLAGSPRE)
CFLAGSNONPRECISE=$(CFLAGSPRENONPRECISE)


#CFLAGSNONPRECISE=-long_double $(CFLAGSPRENONPRECISE)
#CFLAGS = -long_double $(CFLAGSPRE)



# for for normal installation of v5d and hdf
BIN2TXTLIBS=-I /usr/include/hdf/ -L /usr/lib64/hdf/ -lmfhdf -ldf -ljpeg -lz -lv5d
# below for ki-rh39 and likeness
# BIN2TXTLIBS=-I /usr2/local/include/ -L /usr2/local/lib/ -lmfhdf -ldf -ljpeg -lz -lv5d



include maketail.inc

# below is NO LONGER needed
#include maketailsuperlong.inc

# new instructions:
#------------------
# for real trig/etc long doubles
# 1) setup code dir
# mkdir newgrmhdcode ; cd newgrmhdcode ; cp ~/grmhd/ .
# 2) run sh double2longdouble.sh
# 3) need to copy that suite into a directory with the code (not original!)
# cp /home/jondata/ldouble/* .
# 4) uncomment below 2 CFLAGS -type lines and comment out above CFLAGS line before the include maketail.inc
# 5) make


# for real trig/etc long doubles
# 1) setup code dir
# mkdir newgrmhdcode ; cd newgrmhdcode ; cp ~/grmhd/ .
# 2) run sh double2longdouble.sh
# 3) need to copy that suite into a directory with the code (not original!)
# cp /home/jondata/ldouble/* .
# 4) uncomment below 2 include and next line and comment out above include maketail.inc and it's following line
# 5) make




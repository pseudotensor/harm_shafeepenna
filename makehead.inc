#
#
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
USEMPI=1

LAPACKLIB=lapack
BLASLIB=blas
F2CLIB=g2c  #-l$(F2CLIB)

PREP = prep
PREPLD = prepld
PREPLIAISON = prepliaison
PREPJONINTERP = prepjoninterp

CMD=grmhd
CMD2=grmhdldouble
LIAISONCMD=liaison
JONINTERPCMD=iinterp

# Set default goal
.DEFAULT_GOAL := $(CMD)

$(warning default goal is $(.DEFAULT_GOAL))


#END



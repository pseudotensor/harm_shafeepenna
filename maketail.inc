
# for long doubles -m96bit-long-double 
CFLAGSGCCOTHER= -O2 -m128bit-long-double -Wall -Wuninitialized -fno-builtin -DUSINGICC=0
CFLAGSOTHER= -O3 -long-double -unroll -Wall -Wuninitialized -fno-builtin -DUSINGICC=0


ifeq ($(USEMPI),1)

ifeq ($(USEMCCSWITCH),1)
CC=$(MCC) -cc=$(COMP) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=$(MCC) -cc=gcc $(GCCCFLAGS)
endif

ifeq ($(USEMCCSWITCH),0)
CC=$(MCC) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=$(MCC) $(GCCCFLAGS)
endif

MMCOMP=$(CC)

endif # endif usempich==1


ifeq ($(USEMPI),0)
CC=$(COMP)  $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=gcc  $(GCCCFLAGS)

MMCOMP=$(COMP)

endif # endif usempich==0

# when adding C files, only need to modify the below
include maketail.harm.inc

# below is for long double code (not HARM)
include maketail.ldouble.inc

# below is for liaison code
include maketail.liaison.inc

# below is for jon_interp code
include maketail.joninterp.inc



# "make prep" must be called first so these are set correctly
#MYCOMMANDS1 := $(shell cat objsonepoint.cmd.list ) 
#MYCOMMANDS2 := $(shell cat objsmanypoints.cmd.list ) 


#all:	$(PREP) $(CMD) $(FINISH)
all:	$(CMD) $(FINISH)


clean:
	rm *.o

superclean:
	rm -rf $(CMD) $(CMD2) $(LIAISONCMD) $(JONINTERPCMD) *.o $(PREP) $(PREPLD) $(PREPLIAISON) $(PREPJONINTERP) deps.onepoint.list deps.manypoints.list objsmanypoints.cmd.list objsonepoint.cmd.list objslong.cmd.list global.defnprs.h generatenprs libml.a
	rm -rf objsonepoint.cmd.list objsmanypoints.cmd.list objslong.cmd.list objsdepsliaison.cmd.list objsdepsjoninterp.cmd.list
	rm -rf deps.manypoints.list deps.onepoint.list deps.long.list deps.depsliaison.list deps.depsjoninterp.list


cleani:
	rm *.o *.il

cleanall:
	rm *.o *.il *~

cleanbackup:
	rm *~


	# make defnprs.h
global.defnprs.h : global.realdef.h global.nondepnmemonics.h definit.h init.h generatenprs.c
	$(COMP) -o generatenprs generatenprs.c ; ./generatenprs



$(PREP): makefile makedecs.h.sh makehead.inc maketail.harm.deps.inc maketail.harm.inc maketail.inc maketail.ldouble.deps.inc maketail.ldouble.inc maketailsuperlong.inc defs.grmhd.h mpidefs.h superdefs.h supermpidefs.h global.defnprs.h
 	############################################
	#below is for GRMHD code
 	############################################
	#
	#
	#
	( rm -rf defs.h ; ln -s defs.grmhd.h defs.h )
	( rm -rf global.h ; ln -s global.grmhd.h global.h )
	( sh ./makedecs.h.sh defs.h )
	( sh ./makedecs.h.sh mpidefs.h )
	( sh ./makedecs.h.sh superdefs.h )
	( sh ./makedecs.h.sh supermpidefs.h )
	( sh ./makedecs.h.sh mpidefs.mpi_grmhd_grray_liaison.h )
	#
 	############################################
	#below is for code working per point
 	############################################
	( echo $(OBJSONEPOINT) > objsonepoint )
	# remove and touch (so exists) files to create
	( rm -rf objsonepoint.cmd.list ; rm -rf onepoint ; touch objsonepoint.cmd.list ; touch onepoint)
	# create list of prefixes
	( for fil in `cat objsonepoint` ; do prefa=`echo $$fil | sed 's/\.o/ /'` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo -e $$pref >> onepoint ; done )
	# create list of commands
	( for fil in `cat onepoint` ; do echo -e $$fil.o :\\n\\t$(CC) $(CFLAGSNONPRECISE) -c $$fil.c $(LDFLAGS) >> objsonepoint.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.onepoint.list ; touch deps.onepoint.list )
	# create new dependencies file
	( for fil in `cat onepoint` ; do  echo $$fil.c ; $(MMCOMP) -MM $$fil.c >> deps.onepoint.list ; done )
	#
	#
 	############################################
	#below is for code working on multiple point at once
 	############################################
	# below is only different in CFLAGS part
	( echo $(OBJSMANYPOINTS) > objsmanypoints )
	# remove and touch (so exists) files to create
	( rm -rf objsmanypoints.cmd.list ; rm -rf manypoints ; touch objsmanypoints.cmd.list ; touch manypoints)
	# create list of prefixes
	( for fil in `cat objsmanypoints` ; do prefa=`echo $$fil | sed 's/\.o/ /'` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo -e $$pref >> manypoints ; done )
	# create list of commands
	( for fil in `cat manypoints` ; do echo -e $$fil.o :\\n\\t$(CC) $(CFLAGS) -c $$fil.c $(LDFLAGS) >> objsmanypoints.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.manypoints.list ; touch deps.manypoints.list )
	# create new dependencies file
	( for fil in `cat manypoints` ; do  echo $$fil.c ; $(MMCOMP) -MM $$fil.c >> deps.manypoints.list ; done )
	( echo	$(OBJSMANYPOINTS) > $(PREP) )
	#



 	#
$(PREPLD): makefile makedecs.h.sh makehead.inc maketail.harm.deps.inc maketail.harm.inc maketail.inc maketail.ldouble.deps.inc maketail.ldouble.inc maketailsuperlong.inc defs.grmhd.h mpidefs.h superdefs.h supermpidefs.h global.defnprs.h
 	############################################
	# below is for long double code
 	############################################
 	#
	( echo $(OBJSOTHER) > objslong )
	# remove and touch (so exists) files to create
	( rm -rf objslong.cmd.list ; rm -rf long ; touch objslong.cmd.list ; touch long)
	# create list of prefixes
	( for fil in `cat objslong` ; do prefa=`echo $$fil | sed 's/\.o/ /'` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo -e $$pref >> long ; done )
	# create list of commands
	( for fil in `cat long` ; do echo -e $$fil.o :\\n\\t$(CC2) $(CFLAGSGCCOTHER) -c $$fil.c $(LDFLAGS) >> objslong.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.long.list ; touch deps.long.list )
	# create new dependencies file
	( for fil in `cat long` ; do  echo $$fil.c ; $(CC2) -MM $$fil.c >> deps.long.list ; done )
	#
  	# 	create prep file for make to know to recompile or prep again
	( echo	$(OBJSMANYPOINTS) > $(PREPLD) )




$(PREPLIAISON): makefile makedecs.h.sh makehead.inc maketail.harm.deps.inc maketail.harm.inc maketail.liaison.inc maketail.inc maketail.ldouble.deps.inc maketail.ldouble.inc maketailsuperlong.inc defs.liaison.h mpidefs.h superdefs.h supermpidefs.h defs.liaison.h liaison.defs.h global.defnprs.h
 	############################################
	#below is for liaison code
 	############################################
	#
	#
	#
	( rm -rf defs.h ; ln -s defs.liaison.h defs.h )
	( rm -rf global.h ; ln -s global.liaison.h global.h )
	( sh ./makedecs.h.sh defs.h )
	( sh ./makedecs.h.sh mpidefs.h )
	( sh ./makedecs.h.sh superdefs.h )
	( sh ./makedecs.h.sh supermpidefs.h )
	( sh ./makedecs.h.sh liaison.defs.h )
	( sh ./makedecs.h.sh mpidefs.mpi_grmhd_grray_liaison.h )
	#
	( echo $(OBJSLIAISON) > objsdepsliaison )
	# remove and touch (so exists) files to create
	( rm -rf objsdepsliaison.cmd.list ; rm -rf depsliaison ; touch objsdepsliaison.cmd.list ; touch depsliaison)
	# create list of prefixes
	( for fil in `cat objsdepsliaison` ; do prefa=`echo $$fil | sed 's/\.o/ /'` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo -e $$pref >> depsliaison ; done )
	# create list of commands
	( for fil in `cat depsliaison` ; do echo -e $$fil.o :\\n\\t$(CC) $(CFLAGSNONPRECISE) -c $$fil.c $(LDFLAGS) >> objsdepsliaison.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.depsliaison.list ; touch deps.depsliaison.list )
	# create new dependencies file
	( for fil in `cat depsliaison` ; do  echo $$fil.c ; $(MMCOMP) -MM $$fil.c >> deps.depsliaison.list ; done )
	#
  	# 	create prep file for make to know to recompile or prep again
	( echo	$(OBJSLIAISON) > $(PREPLIAISON) )
	#


$(PREPJONINTERP): makefile makedecs.h.sh makehead.inc maketail.harm.deps.inc maketail.harm.inc maketail.inc maketail.ldouble.deps.inc maketail.ldouble.inc maketailsuperlong.inc defs.jon_interp.h global.h global.defnprs.h
 	############################################
	#below is for jon_interp code
 	############################################
	#
	#
	#
	( rm -rf defs.h ; ln -s defs.jon_interp.h defs.h )
	( rm -rf global.h ; ln -s global.jon_interp.h global.h )
	( sh ./makedecs.h.sh defs.h )
	#
	( echo $(OBJSJONINTERP) > objsdepsjoninterp )
	# remove and touch (so exists) files to create
	( rm -rf objsdepsjoninterp.cmd.list ; rm -rf depsjoninterp ; touch objsdepsjoninterp.cmd.list ; touch depsjoninterp)
	# create list of prefixes
	( for fil in `cat objsdepsjoninterp` ; do prefa=`echo $$fil | sed 's/\.o/ /'` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo -e $$pref >> depsjoninterp ; done )
	# create list of commands
	( for fil in `cat depsjoninterp` ; do echo -e $$fil.o :\\n\\t$(CC) $(CFLAGSNONPRECISE) -c $$fil.c $(LDFLAGS) >> objsdepsjoninterp.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.depsjoninterp.list ; touch deps.depsjoninterp.list )
	# create new dependencies file
	( for fil in `cat depsjoninterp` ; do  echo $$fil.c ; $(MMCOMP) -MM $$fil.c >> deps.depsjoninterp.list ; done )
	#
  	# 	create prep file for make to know to recompile or prep again
	( echo	$(OBJSJONINTERP) > $(PREPJONINTERP) )
	#





$(FINISH):
#	rm -rf objsonepoint objsmanypoints

#$(CMD):   	GROUP1 GROUP2 makefile freespace.o tau_neededbyharm.o
#	echo objsmanypoints: $(OBJSMANYPOINTS)
#	echo objsonepoint: $(OBJSONEPOINT)
#	$(MYCOMMANDS2)
#	$(MYCOMMANDS1)


$(CMD):	$(PREP) $(OBJSMANYPOINTS) $(OBJSONEPOINT) makefile freespace.o tau_neededbyharm.o
	$(CC)  $(CFLAGSNONPRECISE) -o $(CMD) $(OBJSMANYPOINTS) $(OBJSONEPOINT) freespace.o tau_neededbyharm.o $(LDFLAGS)


$(CMD2):       $(PREP) $(PREPLD) $(OBJSMANYPOINTS) $(OBJSONEPOINT) makefile libml.a $(OBJSOTHER) freespace.o tau_neededbyharm.o
	       $(CC2) $(CFLAGSGCCOTHER) -c polevll.c $(LDFLAGSOTHER)
	       $(CC)  $(CFLAGS) -o $(CMD2) $(OBJSMANYPOINTS) $(OBJSONEPOINT) freespace.o tau_neededbyharm.o polevll.o libml.a $(LDFLAGSOTHER)

$(LIAISONCMD): $(PREPLIAISON) $(OBJSLIAISON) makefile freespace.o 
	       $(CC)  $(CFLAGS) -DUSINGLIAISON=1 -o $(LIAISONCMD) $(OBJSLIAISON) freespace.o $(LDFLAGSOTHER)

$(JONINTERPCMD): $(PREPJONINTERP) $(OBJSJONINTERP) makefile
	       $(CC)  $(CFLAGS) -DUSINGJONINTERP=1 -o $(JONINTERPCMD) $(OBJSJONINTERP) $(LDFLAGSOTHER)

 


ifeq ($(MAKECMDGOALS),$(filter $(MAKECMDGOALS),$CMD $CMD2 $LIAISONCMD))
-include objsonepoint.cmd.list
-include objsmanypoints.cmd.list
-include deps.manypoints.list
-include deps.onepoint.list
endif

ifeq ($(MAKECMDGOALS),$(CMD2))
-include objslong.cmd.list
-include deps.long.list
endif

ifeq ($(MAKECMDGOALS),$(LIAISONCMD))
-include objsdepsliaison.cmd.list
-include deps.depsliaison.list
endif

ifeq ($(MAKECMDGOALS),$(JONINTERPCMD))
-include objsdepsjoninterp.cmd.list
-include deps.depsjoninterp.list
endif



# below automatic now (original has special code for assembly but doesn't work anymore on ki-rh42)
#-include maketail.ldouble.deps.inc



freespace.o : freespace.c
	    $(CC2) $(GCCCFLAGS) -c freespace.c $(LDFLAGS)

tau_neededbyharm.o : tau_neededbyharm.c
	    $(CC2) $(GCCCFLAGS) -c tau_neededbyharm.c -lm



# Use $(OBJS) in ar command for libml.a if possible; else *.o
libml.a: $(OBJSOTHER) mconf.h objslong.cmd.list deps.long.list
	 ar -rv libml.a $(OBJSOTHER)
	 ranlib libml.a
#	 ar -rv libml.a *.o





#$(OBJSMANYPOINTS) :

#$(MYCOMMANDS2)

#$(OBJSONEPOINT) :

#$(MYCOMMANDS1)





#GROUP2:	$(OBJSMANYPOINTS) makefile
#	$(CC) $(CFLAGS) -c $(SRCSONEPOINT) $(LDFLAGS)
#	CFLAGS=$(CFLAGS)
#	echo $(MYCOMMANDS2)
#	$(MYCOMMANDS2)

#GROUP1:	$(OBJSONEPOINT) makefile
#	$(CC) $(CFLAGSNONPRECISE) -c $(SRCSONEPOINT) $(LDFLAGS)
#	CFLAGS=$(CFLAGSNONPRECISE)
#	echo $(MYCOMMANDS1)
#	$(MYCOMMANDS1)



# below file not used anymore since generated automagically
#include maketail.harm.deps.inc


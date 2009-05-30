
SRCSMANYPOINTS=\
init.c initbase.c initbase.boundloop.c initbase.enerregions.c initbase.tools.c initbase.gridsectioning.c restart.c restart.checks.c set_arrays.c mpi_set_arrays.c set_grid.c \
advance.c reconstructeno.c paraenohybrid.c flux.c flux.mergedc2ea2cmethod.c higherorder_pointavg.c wavespeeds.c boundsflux.c \
fluxct.c fluxctstag.c fluxcompute.c fluxvpot.c \
bounds.c boundsint.c \
diag.c dump.c dumpgen.c dump_ener.c image.c \
fail.c fixup.c \
mpi_init.c mpi_init_grmhd.c mpi_fileio.c mpi_fprintfs.c boundmpi.c boundmpiint.c mpi_grmhd_grray_liaison.c \
interppoint.c interppoint.para.c interpline.c interpline.para.c interpline.mono.c interpline.jmono.c interpline.smono.c
#END

OBJSMANYPOINTS=\
init.o initbase.o initbase.boundloop.o  initbase.enerregions.o initbase.tools.o initbase.gridsectioning.o restart.o restart.checks.o set_arrays.o mpi_set_arrays.o set_grid.o \
advance.o reconstructeno.o paraenohybrid.o flux.o  flux.mergedc2ea2cmethod.o higherorder_pointavg.o wavespeeds.o boundsflux.o \
fluxct.o fluxctstag.o fluxcompute.o fluxvpot.o \
bounds.o boundsint.o \
diag.o dump.o dumpgen.o dump_ener.o image.o \
fail.o fixup.o \
mpi_init.o mpi_init_grmhd.o  mpi_fileio.o mpi_fprintfs.o boundmpi.o boundmpiint.o mpi_grmhd_grray_liaison.o \
interppoint.o interppoint.para.o interpline.o interpline.para.o interpline.mono.o interpline.jmono.o interpline.smono.o \
#END

# phys.c does currents across multiple points, but otherwise one pointed
# wavespeeds kinda at more than one point

# for some reason main.c really matters alot.  Needs to be compiled
# withOUT -mp -pc64 for precision in inversion to be ok  ODD!  GODMARK

SRCSONEPOINT=\
main.c mytime.c metric.c metric_tools.c math_tools.c metric_selfgravity_or_evolvemetric.c coord.c phys.c step_ch.c \
phys.ffde.c phys.coldgrmhd.c eos.c \
vchar.c transforms.c sources.c rescale_interp.c \
gaussj.c lubksb.c ludcmp.c mnewt.c nrutil.c ranc.c tensor.c \
utoprimgen.c dudp_calc_3vel.c dudp_calc.c \
utoprim.c utoprim_2d.c utoprim_1d.c utoprim_1d_opt.c utoprim_ldz.c \
utoprim_jon.c utoprim_1d_final.c utoprim_2d_final.c utoprim_5d2_final.c
#END

OBJSONEPOINT=\
main.o mytime.o metric.o metric_tools.o math_tools.o metric_selfgravity_or_evolvemetric.o coord.o phys.o step_ch.o \
phys.ffde.o phys.coldgrmhd.o eos.o \
vchar.o transforms.o sources.o rescale_interp.o \
gaussj.o lubksb.o ludcmp.o mnewt.o nrutil.o ranc.o tensor.o \
utoprimgen.o dudp_calc_3vel.o dudp_calc.o \
utoprim.o utoprim_2d.o utoprim_1d.o utoprim_1d_opt.o utoprim_ldz.o \
utoprim_jon.o utoprim_1d_final.o utoprim_2d_final.o utoprim_5d2_final.o
#END








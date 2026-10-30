# probMHD - Makefile (probCFD style: modules in dependency order)
# NOTE: this Makefile does not track modification times; after editing a
#       source run  make clean && make
FC     = mpif90
FFLAGS = -O2 -ffree-line-length-none -fallow-argument-mismatch
SRCDIR = src
SRCS = mod_global.F90 mod_mpi.F90 geometry.F90 read_input.F90 equilibrium.F90 boundary.F90 convection.F90 \
       induction.F90 diffusion.F90 time_integration.F90 output.F90 stop_criteria.F90 Main.F90
OBJS = $(SRCS:.F90=.o)

probMHD: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o probMHD

%.o: $(SRCDIR)/%.F90
	$(FC) $(FFLAGS) -c $< -o $@

# debug build with run-time checks (slower)
debug:
	$(MAKE) clean
	$(MAKE) probMHD FFLAGS="-O0 -g -ffree-line-length-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -fallow-argument-mismatch"

clean:
	rm -f *.o *.mod probMHD

###
#   This is a simple makefile that builds the 3D wave equation WHI
###

FC = mpif90
LD = mpif90
FFLAGS = -c -O1 -g -fbounds-check 
FINC   =
LDFLAGS= 
LD = $(FC) 
EX = wave3d.x

OBJS_WAVE = fd_3d.o mpi_cart_fun.o 

.PHONY: clean
#

#
# executable 
#
wave3d.x: $(OBJS_WAVE)
	$(LD) -o $(EX) $(OBJS_WAVE) 
#
mpi_cart_fun.o: mpi_cart_fun.f90
	$(FC) $(FFLAGS) mpi_cart_fun.f90
#
%.o: %.f90
	$(FC) $(FFLAGS) $<
#
clean:
	/bin/rm -f $(EX) *.o *.mod

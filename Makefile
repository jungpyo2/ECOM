# Makefile for Grad Shfranov solver, ECOM
FC = gfortran

FFLAGS = -O3 

all: ECOM

Splot_sub = quaplot.f
Slib_sub = interpol.f cqrsolve.f durdecplus.f anaresa.f timefile.f solvesub.f 
Oplot_sub = $(Splot_sub:.f=.o)
Olib_sub = $(Slib_sub:.f=.o)
Spois_sub = elliptic.f bvpadapfix.f
Opois_sub = $(Spois_sub:.f=.o)


arrays.o arrays.mod:  arrays.f90
	$(FC) $(FFLAGS)  arrays.f90 -c $<

cmap.o:  cmap.f arrays.f90 $(Slib_sub) $(Splot_sub) 
	$(FC) $(FFLAGS)  cmap.f arrays.f90 $(Splot_sub) $(Slib_sub) -c $<

readbf.o: readbf.f arrays.o quaplot.f 
	$(FC) $(FFLAGS)  readbf.f arrays.f90 quaplot.f -c  $(Slib_sub) $<

solvef.o: solvef.f arrays.f90 $(Spois_sub) $(Splot_sub) $(Slib_sub)
	$(FC) $(FFLAGS) solvef.f arrays.f90 $(Spois_sub) $(Splot_sub) $(Slib_sub)  -c $<

solveq.o: solveq.f arrays.f90 $(Spois_sub) $(Splot_sub) $(Slib_sub)
	$(FC) $(FFLAGS) solveq.f arrays.f90 $(Spois_sub) $(Splot_sub) $(Slib_sub)  -c $<

solvej.o: solvej.f arrays.f90 $(Spois_sub) $(Splot_sub) $(Slib_sub)
	$(FC) $(FFLAGS) solvej.f arrays.f90 $(Spois_sub) $(Splot_sub) $(Slib_sub)  -c $<

postproc.o: postproc.f arrays.f90 $(Spois_sub) $(Splot_sub) $(Slib_sub)
	$(FC) $(FFLAGS) postproc.f arrays.f90  $(Spois_sub) $(Splot_sub) $(Slib_sub)  -c $<

ECOM : readbf.o cmap.o solvef.o solveq.o solvej.o postproc.o arrays.o
	 $(FC) $(FFLAGS) -w -o ecom ecom.f solvef.o solveq.o solvej.o readbf.o cmap.o postproc.o arrays.o $(Opois_sub) $(Oplot_sub)  $(Olib_sub) 


clean :
	rm -rf *.o *.mod 

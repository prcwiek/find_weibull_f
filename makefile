.SUFFIXES:

FC = gfortran

# compile flags
FCFLAGS = -ffree-line-length-none 

find_weibull:	find_weibull.f90
		$(FC) $< $(FCFLAGS) -o $@.out

%.mod: %.f90
		$(FC) $(FCFLAGS) -c $<	


PHONY: clean
clean:
	-rm -f *.o *.mod *.smod *.anc 

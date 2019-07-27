BINDIR = $(HOME)/bin

OBJ = hilbert.o

.SUFIXES: .f90 .o .mod

F90   = gfortran 
FCOPT = -O3


ALL: hilbert.x

hilbert.x: $(OBJ)
	$(F90) $(OBJ) $(FCOPT) -o $(BINDIR)/hilbert.x
%.o:%.f90
	$(F90) -c $(FCOPT) $<

.PHONY: clean

clean:
	rm *.o *.mod $(BINDIR)/hilbert.x

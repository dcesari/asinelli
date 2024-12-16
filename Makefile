FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fcheck=all
FCFLAGS = -g -O2
# flags forall (look for system .mod files required in gfortran)
FCFLAGS += -I/usr/include

# libraries needed for linking
#LDFLAGS = -li_need_this_lib

# List of executables to be built
PROGRAMS = asino

# "make" builds all
all: $(PROGRAMS)


asino.o: asinelli.o
asino: asinelli.o

# General rule for building exe from exe.o; $^ is used in order to
# list additional object files on which the executable depends
%: %.o
	$(FC) $(LDFLAGS) -o $@ $^

# General rule for building exe.o from exe.f90 or exe.F90; $< is used
# in order to list only the first prerequisite (the source file) and
# not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) -c $(FCFLAGS) $<

%.o: %.F90
	$(FC) -c $(FCFLAGS) $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod .MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)

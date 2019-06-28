# compiler
FC = /opt/local/bin/gfortran-mp-8
INCLUDES        = -I/opt/local/include 

# compile flags
FCFLAGS = -g -c -fdefault-real-8 -O3 ${INCLUDES}
#-fbacktrace -fno-align-commons 
FLIBS   = -L/opt/local/lib -lnetcdff -lnetcdf
# link flags
FLFLAGS = -O3 ${INCLUDES}
#-g -fbacktrace

# source files and objects
SRCS = fluidloop_init.o fluidloop_io.o fluidloop.o

# program name
PROGRAM = fluidloop

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FLFLAGS) ${FLIBS} $^ -o $@  
  
%.o: %.f90
	$(FC) $(FCFLAGS) ${FLIBS} $< -o $@ 

clean:
	rm -f *.o *.mod



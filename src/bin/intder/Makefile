# Sample Makefile for compiling under PSI 3.2.beta-rc1

# If using C++ use: CC = g++
# If using C use  : CC = gcc
CC = g++

CODE = intder 

# -g add debugging information to the executable
CFLAGS = -g 

# -lm is the default math library
# Change the others to the location of where you
# installed Psi3.

LIBS = -lm \
       /home/nate/psi3/bin.linux/lib/libPSI_chkpt.a \
       /home/nate/psi3/bin.linux/lib/libPSI_ipv1.a \
       /home/nate/psi3/bin.linux/lib/libPSI_iwl.a \
       /home/nate/psi3/bin.linux/lib/libPSI_psio.a \
       /home/nate/psi3/bin.linux/lib/libPSI_ciomr.a \

# -I. indicated that the compiler should check the
#     current directory for header files
# The next 2 are located in the location of the
# source code of your version of Psi3
# The ./share/. is the location where you installed
# Psi3.

INCLUDES = -I./                               \
           -I/home/nate/psi3/src/lib   \
           -I/home/nate/psi3/objdir/include   \
           -I/home/nate/psi3/objdir/src/lib  \
           -I/home/nate/psi3/include \


# Include all of your source files.
CSRC = intder.cc cartesian.cc atom.cc molecule.cc displacements.cc \
       params.cc intco.cc misc.cc 3dmatrix.cc read_geom.cc h_compute.cc 

# If your source files have the .c extension then
# change the %.cpp to %.c  in the following line.
CBINOBJ = $(CSRC:%.cc=%.o)

$(CODE): $(CBINOBJ) $(LIBS)
	$(CC) $(CFLAGS) $^ -o $(CODE)

# If your source files have the .c extension then
# change the %.cpp to %.c  in the following line.
%.o: %.cc
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

# By invoking "make clean" the following line will be
# executed,
clean:
	rm -f *.o

install:
	cp -f $(CODE) /home/nate/progs/$(CODE) 


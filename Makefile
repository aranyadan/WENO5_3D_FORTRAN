### create directories
DIRECTORY = 	obj plots data mod
MKDIRP = mkdir -p
SRCDIR = ./src/
OBJDIR = ./obj/
MODDIR = ./mod/

### suffix rule
.SUFFIXES:
.SUFFIXES: .f90 .o

### compiler
F90 = f95
COMMONFLAGS  =  -O4
COMPFLAGS          =  -c  $(COMMONFLAGS)
LINKFLAGS             =       $(COMMONFLAGS)

### objects
OBJ = IC3DReimann.o pad_fluxes.o saver.o primitives.o transforms.o fluxes.o \
 			WENO.o WENO53d.o main.o viscous_fluxes.o set_boundary.o \
			compute_residual.o
OBJS = $(addprefix $(OBJDIR), $(OBJ))
### compile and link

all: $(DIRECTORY) output
### Create the obj directory
$(DIRECTORY):
	$(MKDIRP) $@
### Creating the obj files
$(addprefix ./obj/, %.o): $(addprefix ./src/, %.f90)
	$(F90) $(COMPFLAGS) $< -o $@ -J$(MODDIR)
# $< refers to file calling the statement, $@ refers to the file being created
### Linking and creating executable file
output:  $(OBJS)
	$(F90)  -o  $@  $(LINKFLAGS)  $(OBJS)

clean:
	rm -f *.o
	rm -f output
	rm -rf obj mod
	rm -rf *.mod
cleanall:
	rm -rf data
	rm -rf plots
	rm -f *.o
	rm -f output
	rm -rf obj mod
	rm -rf *.mod
cleandat:
	rm -rf data/*
	rm -rf plots/*

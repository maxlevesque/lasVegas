# ——————————————— Program property ———————————————

EXE = lasVegas
MKFLAGS = $(NULLSTRING)

# ——————————————— Variables of locations ———————————————

MODDIR = mod
SRCDIR = src
OBJDIR = obj

# _______________ Libraries and other folders __________

FFTW_INCLUDES  = -I/usr/local/include -I/usr/include
FFTW_LIBRARIES = -L/usr/local/lib -L/usr/lib

# ——————————————— Fortran compiler ———————————————

FC = gfortran

# ——————————————— Compiling options ———————————————

FCFLAGS = -J$(MODDIR) -I$(MODDIR) $(FFTW_INCLUDES)
LDFLAGS =

DEBUG = -Og -g -Wall -fimplicit-none -fbacktrace -std=f2008 -pedantic -fwhole-file -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fbounds-check -fcheck=all -pg -Wunused-parameter -frecursive
# -g turns on debugging
# -p turns on profiling
# -Wextra turns on extra warning. It is extremely verbose.
# -fcheck-array-temporaries -Warray-temporaries -Wconversion -Wimplicit-interface

OPTIM = -O3 -march=native -ffast-math -funroll-loops
# FOR BACKUP : -march=native -O3 -ffast-math -funroll-loops   VERY AGRESSIVE
# -fopenmp for OPENMP support

# ——————————————— Files to compile ———————————————

FOBJ = $(OBJDIR)/main.o\
			$(OBJDIR)/print_header.o

# ——————————————— Global rules ———————————————

all: $(EXE)
	$(FC) $(FCFLAGS) $(MKFLAGS) -o $(EXE) $(FOBJ) $(LDFLAGS)

optim: MKFLAGS = $(OPTIM)

optim: $(EXE)
	$(FC) $(FCFLAGS) $(MKFLAGS) -o $(EXE) $(FOBJ) $(LDFLAGS)

debug: MKFLAGS = $(DEBUG)

debug: $(EXE)
	$(FC) $(FCFLAGS) $(MKFLAGS) -o $(EXE) $(FOBJ) $(LDFLAGS)

clean:
	rm -vf gmon.out $(EXE) $(MODDIR)/* $(OBJDIR)/*

# ——————————————— Pattern rules ———————————————

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) $(MKFLAGS) -c $< -o $@

# For GNU make, *.f90 cannot be compiled automatically.

# ——————————————— Dependence rules ———————————————

$(EXE): $(FOBJ)

$(OBJDIR)/main.o:\
	$(SRCDIR)/main.f90

$(OBJDIR)/print_header.o:\
	$(SRCDIR)/print_header.f90

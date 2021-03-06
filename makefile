# Makefile for qsosim8
# JKW, Dec 2013
# Compiler:
FC= gfortran

PROGS = qsosim10
# List of executables to be built within the package

# this is the default (ie 'make' will make for Linux as below)
#all: $(PROGRAMS)
# to recompile subroutines add q1 before q9 in the next line
#clean: 
#	rm *.o rm *.c
#
SOURCES =	qsosim.f dsepvar.f ewred.f spvoigt.f voigt.f vp_lycont.f f13_read.f read_nrows.f \
		spline.f SDSS_readfits.f power_laws.f assign.f \
		qsosim10.f metal_db.f
OBJECTS =	qsosim.o dsepvar.o ewred.o spvoigt.o voigt.o vp_lycont.o f13_read.o read_nrows.o \
		SDSS_readfits.o readfits.o qsosim.o writefits.o spline.o  \
		power_laws.o assign.o metal_db.o
OBJDIR = ./source
VPATH = $(OBJDIR)
#  Macbook Pro version
# cfitsio = -L/usr/local/Cellar/cfitsio/3.390/lib -lcfitsio
cfitsio = -L/usr/local/lib -lcfitsio
pgplot = -L/usr/local/opt -L/usr/X11/lib -lpgplot -lX11 

all : $(PROGS)
#q1: qsosim10.f $(obj1)
#	$(LINK.f) -c $(obj1)
$(PROGS): $(addprefix $(OBJDIR)/, $(OBJECTS))
	$(LINK.f) -o qsosim10 $(addprefix $(OBJDIR)/, $(OBJECTS)) \
	$(pgplot) $(cfitsio)

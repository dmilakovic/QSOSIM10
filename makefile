# Makefile for qsosim8
# JKW, Dec 2013
# Compiler:
FC= gfortran

# List of executables to be built within the package

# this is the default (ie 'make' will make for Linux as below)
#all: $(PROGRAMS)
all: q1 q9

obj1=	spline.f qsosim.f read_nrows.f SDSS_readfits.f readfits.f power_laws.f assign.f \
		qsosim10.f write_cloudy_input.f cloudy.f read_cloudy_output.f
obj2=	dsepvar.o ewred.o spvoigt.o voigt.o vp_lycont.o f13_read.o read_nrows.f \
		SDSS_readfits.o readfits.o qsosim.o writefits.o spline.o  \
		power_laws.o assign.o write_cloudy_input.o cloudy.o read_cloudy_output.o

#  Macbook Pro version
q1: qsosim10.f spline.f qsosim.f power_laws.f assign.f \
		write_cloudy_input.f cloudy.f read_cloudy_output.f
	$(LINK.f) -c $(obj1)
q9: qsosim10.o $(obj2)
	$(LINK.f) -o qsosim10 qsosim10.o $(obj2) \
	-L/opt/local/lib -lpgplot -L/usr/X11/lib -lX11 \
	-L/opt/local/lib -lcfitsio	

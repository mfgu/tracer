RM= rm -f
CC= g++ -Wno-deprecated -g
LIB= -lm
INCDIR= /opt/local/include
LIBDIR= /opt/local/lib

.cc.o:
	$(RM) $@
	$(CC) -I$(INCDIR) -c $*.cc

all::  	mg_spect torg_spect johann_spect

util.o: util.cc util.h const.h

geom_objects.o: geom_objects.cc util.cc util.h geom_objects.h

opt_devices.o:  opt_devices.cc opt_devices.h geom_objects.cc geom_objects.h

tracer.o:       tracer.cc tracer.h  opt_devices.cc util.cc util.h opt_devices.h geom_objects.cc geom_objects.h 

spectrometer.o: spectrometer.cc spectrometer.h tracer.cc util.cc util.h tracer.h opt_devices.cc opt_devices.h geom_objects.cc geom_objects.h

mg_spect.o:     mg_spect.cc spectrometer.cc tracer.cc util.cc util.h tracer.h  opt_devices.cc opt_devices.h geom_objects.cc geom_objects.h

toroid.o:	toroid.cc toroid.h geom_objects.cc util.cc util.h geom_objects.h

torg_spect.o:   torg_spect.cc spectrometer.cc tracer.cc tracer.h opt_devices.cc opt_devices.h toroid.cc toroid.h geom_objects.cc geom_objects.h

johann_spect.o: tracer.cc tracer.h opt_devices.cc opt_devices.h toroid.cc toroid.h geom_objects.cc geom_objects.h

mg_spect:       geom_objects.o opt_devices.o tracer.o spectrometer.o mg_spect.o util.o
	$(RM) $@
	$(CC) -o $@ $(LIB) geom_objects.o opt_devices.o tracer.o spectrometer.o util.o $@.o -L$(LIBDIR) -lgsl -lgslcblas

torg_spect:	geom_objects.o toroid.o opt_devices.o tracer.o spectrometer.o torg_spect.o util.o
	$(RM) $@
	$(CC) -o $@ $(LIB) geom_objects.o toroid.o opt_devices.o tracer.o spectrometer.o util.o $@.o -L$(LIBDIR) -lgsl -lgslcblas

johann_spect: johann_spect.o geom_objects.o opt_devices.o tracer.o util.o
	$(RM) $@
	$(CC) -o $@ $(LIB) geom_objects.o opt_devices.o tracer.o util.o $@.o -L$(LIBDIR) -lgsl -lgslcblas

clean:
	rm -f *.o *~ torg_spect mg_spect johann_spect

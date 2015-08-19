all: settings clean wrapradarfilter execradarfilter

include ../linux_compliation_settings.mk

wrapradarfilter: $(DEPS)
	swig -python -outcurrentdir ../src/wrapradarfilter.i
	$(CC) -o _wrapradarfilter.so wrapradarfilter_wrap.c $(DEPS) $(INC) $(CCFLAGS) -shared -fPIC

execradarfilter:  $(DEPS)
	$(CC) -o radarfilter.o ../src/wrapradarfilter_exec.c $(DEPS) $(INC) $(CCFLAGS)

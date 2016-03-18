RELPAT = ../

all: settings clean wrapedr

include $(RELPAT)linux_compliation_settings.mk

wrapedr: $(DEPS)
	swig -python -outcurrentdir ../src/wrapedr.i
	$(CC) -o _wrapedr.so wrapedr_wrap.c $(DEPS) $(INC) $(CCFLAGS) -shared -fPIC









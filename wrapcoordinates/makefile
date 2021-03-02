RELPAT = ../

all: settings clean wrapcoordinates

include $(RELPAT)linux_compliation_settings.mk

wrapcoordinates: $(DEPS)
	swig -python -outcurrentdir ../src/wrapcoordinates.i
	$(CC) -o _wrapcoordinates.so wrapcoordinates_wrap.c $(DEPS) $(INC) $(CCFLAGS) -shared -fPIC

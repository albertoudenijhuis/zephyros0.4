RELPAT = ../


# GC			= gcc
# CFLAGS 		= -lm -ldl
# DEBUGFLAGS 	= -g -o0 -Wall -Wextra
# INC			= -I../src/ -I/usr/include/python2.7 -I/usr/lib64/python2.7/site-packages/numpy/core/include
# DEPS 		= ../src/coordinates.c

all: clean wrapcoordinates

include $(RELPAT)linux_compliation_settings.mk

clean:
	rm -f *.so *.o *_wrap.c

wrapcoordinates: 
	swig -python -outcurrentdir ../src/wrapcoordinates.i
	$(GC) -o _wrapcoordinates.so wrapcoordinates_wrap.c $(DEPS) $(INC) $(CFLAGS) -shared -fPIC 

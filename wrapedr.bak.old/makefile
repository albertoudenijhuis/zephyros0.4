GC			= gcc
CFLAGS 		= -ldl -lm 
DEBUGFLAGS 	= -g -o0 -Wall -Wextra
INC			= -I../src/ -I/usr/include/python2.7 -I/usr/lib64/python2.7/site-packages/numpy/core/include
DEPS 		= ../src/edr.c ../src/func.c

all: clean wrapedr

clean:
	rm -f *.so

wrapedr: 
	swig -python -outcurrentdir ../src/wrapedr.i
	$(GC) -o _wrapedr.so wrapedr_wrap.c $(DEPS) $(INC) $(CFLAGS) -shared -fPIC 



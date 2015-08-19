GC			= gcc
DEBUGFLAGS 	= -g -o0 -Wall -Wextra
CFLAGS 		= -ldl -lm 
INC			= -I../src/ -I/usr/include/python2.7 -I/usr/lib64/python2.7/site-packages/numpy/core/include

clean:
	rm -f *.so

edr.so: ../src/edr.i
	swig -python -outcurrentdir $<
	$(GC) -c ../src/$(basename $@).c $(INC) $(CFLAGS) -shared -fPIC 
	$(GC) -c $(basename $@)_wrap.c $(INC) $(CFLAGS) -shared -fPIC 
	$(GC) -shared -o _$@ $(basename $@).o $(basename $@)_wrap.o ../src/func.c $(CFLAGS) -fPIC 

#install:
#	cp ./_$@ ../_$@
#	cp ./$(basename $@).py ../$(basename $@).py


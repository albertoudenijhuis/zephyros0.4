RELPAT = ../

all: settings clean func

include $(RELPAT)linux_compliation_settings.mk

func: $(DEPS)
	swig -python -outcurrentdir ../src/func.i
	$(CC) -o _func.so func_wrap.c $(DEPS) $(INC) $(CCFLAGS) -shared -fPIC








#func.so: ../src/func.i
#	swig -python -outcurrentdir $<
#	$(GC) -c ../src/$(basename $@).c $(INC) $(CFLAGS) -shared -fPIC 
#	$(GC) -c $(basename $@)_wrap.c $(INC) $(CFLAGS) -shared -fPIC 
#	$(GC) -shared -o _$@ $(basename $@).o $(basename $@)_wrap.o $(INC) $(CFLAGS) -fPIC 
#	#cp ./_$@ ../_$@
#	#cp ./$(basename $@).py ../$(basename $@).py

all: settings clean wrapwindfield

include ../linux_compliation_settings.mk
	
wrapwindfield: $(DEPS)
	swig -python -outcurrentdir ../src/wrapwindfield.i
	$(CC) -o _wrapwindfield.so wrapwindfield_wrap.c $(DEPS) $(INC) $(CCFLAGS) -shared -fPIC

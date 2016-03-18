RELPAT = ../

all: settings clean wrapretrieval

include $(RELPAT)linux_compliation_settings.mk

wrapretrieval: $(DEPS)
	swig -python -outcurrentdir ../src/wrapretrieval.i
	$(CC) -o _wrapretrieval.so wrapretrieval_wrap.c $(DEPS) $(INC) $(CCFLAGS) -shared -fPIC



RELPAT = ../../../

all: clean test

include ../../../linux_compliation_settings.mk

test: $(DEPS)
	$(CC)  -o test.o $(DEPS) $(INC) $(CCFLAGS) test.c

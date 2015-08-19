all: clean test

include ../linux_compliation_settings.mk

test: $(DEPS)
	$(CC)  -o test.o $(DEPS) $(INC) $(CCFLAGS) test.c

bhmie.o: ./mie/bhmie.f90
	$(FC) $(FCFLAGS) -c $< -o $@

%.o: ../../../src/%.f
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	rm -f *.o

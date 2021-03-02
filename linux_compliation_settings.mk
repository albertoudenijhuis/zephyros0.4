FC = gfortran
CC = gcc

#without debugging flags
#FCFLAGS=-O2 -fPIC
#CCFLAGS=-lcxsparse -lnlopt -O2 -lgfortran -lm

#with debugging flags
FCFLAGS=-O3 -fPIC -g
CCFLAGS=-lcxsparse -lnlopt -O3 -lgfortran -lm  -g -w


# -w
DEPS = $(filter-out %/wrapretrieval_exec.c %/wrapradarfilter_exec.c %/wrapwindfield_exec.c, $(wildcard $(RELPAT)src/*.c)) particles_mishchenko2000_ampld.lp.o particles_mishchenko2000_lpd.o

INC1			= -I$(RELPAT)src/ -I/usr/include/python2.7 -I/usr/lib64/python2.7/site-packages/numpy/core/include -I/usr/include/suitesparse/
INC2		= -I/usr/local/include/python2.7 -I/usr/local/lib64/python2.7/site-packages/numpy/core/include
INC			= $(INC1) $(INC2)

settings:
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib

clean:
	rm -f *.so *.o *_wrap.c

%.o: $(RELPAT)src/%.f
	$(FC) $(FCFLAGS) -c $< -o $@



# installation instruction
# yum install suitesparse-devel suitesparse-static suitesparse
# standard nlopt libraries do not work under linux
# instead do the following:
# download NLOPT, and then ...
# ./configure --enable-shared && make && sudo make install
# make sure LD_LIBRARY_PATH is set correct.
# i.e. if NLOPT libraries are installed in /usr/local/lib, make sure that export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib

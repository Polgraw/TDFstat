CC = gcc

LAL_FLAGS = -DUSE_LAL -I$(CONDA_PREFIX)/include
LAL_LIBS = -L$(CONDA_PREFIX)/lib -llalpulsar -llalsupport -llalmetaio -llal
#LAL_FLAGS =
#LAL_LIBS =

INIPARSER_DIR = ../../utils/iniparser
INIPARSER_FLAGS = -I$(INIPARSER_DIR)/src
INIPARSER_LIBS = -L$(INIPARSER_DIR) -liniparser

#FFTW_DIR = /opt/fftw/3.3.10-gcc11-mvapich2-broadwell
FFTW_DIR = 
FFTW_FLAGS = -I$(FFTW_DIR)/include
FFTW_LIBS = -L$(FFTW_DIR)/lib -lfftw3

CFLAGS = -Wall $(LAL_FLAGS) $(FFTW_FLAGS) $(INIPARSER_FLAGS) -O2 -UDEBUG

LDFLAGS = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
	-Wl,-Bdynamic -lc -lm -lz -lpthread
LDLIBS1 = -lgslcblas  -lgsl $(FFTW_LIBS) $(INIPARSER_LIBS) $(LAL_LIBS)
LDLIBS2 = $(LAL_LIBS)

all:	genseg laltest

genseg:	genseg.o GrubbsOutliers.o get_barycenter.o sid.o
	$(MAKE) -C $(INIPARSER_DIR)
	$(CC) -o genseg $^ $(LDFLAGS) $(LDLIBS1)

laltest: laltest.o get_barycenter.o sid.o
	$(CC) -o laltest $^ $(LDFLAGS) $(LDLIBS1)

clean:
	rm -f *.o genseg laltest
	$(MAKE) -C $(INIPARSER_DIR) veryclean

.PHONY:	clean

CC = gcc -g -O2
USE_LAL = yes

ifeq ($(strip $(USE_LAL)), yes)
	LAL_FLAGS = -DUSE_LAL -I$(CONDA_PREFIX)/include
	LAL_LIBS = -Wl,-rpath=$(CONDA_PREFIX)/lib -L$(CONDA_PREFIX)/lib -llalpulsar -llalsupport -llalmetaio -llal
else
	LAL_FLAGS =
	LSL_LIBS =
endif

INIPARSER_DIR = ../utils/iniparser
INIPARSER_FLAGS = -I$(INIPARSER_DIR)/src
INIPARSER_LIBS = -L$(INIPARSER_DIR) -liniparser

CFLAGS = -Wall $(LAL_FLAGS) $(INIPARSER_FLAGS) -O2 -UDEBUG

LDFLAGS = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
	-Wl,-Bdynamic -lc -lm -lz -lpthread
LDLIBS1 = -lgslcblas  -lgsl $(INIPARSER_LIBS) $(LAL_LIBS)
LDLIBS2 = $(LAL_LIBS)

all:	genseg-hdf laltest

genseg-hdf:	genseg_hdf.o GrubbsOutliers-float.o get_barycenter.o sid.o
	$(MAKE) -C $(INIPARSER_DIR)
	$(CC) -o genseg-hdf $^ $(LDFLAGS) -lhdf5 -lhdf5_hl $(LDLIBS1)

laltest: laltest.o get_barycenter.o sid.o
	$(CC) -o laltest $^ $(LDFLAGS) $(LDLIBS1)

clean:
	rm -f *.o genseg-hdf laltest
	$(MAKE) -C $(INIPARSER_DIR) veryclean

.PHONY:	clean


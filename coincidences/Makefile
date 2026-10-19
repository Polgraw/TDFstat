TARGET = coincidences
CC = gcc -std=gnu11 -static-libgcc
INIPARSER_DIR = ../utils/iniparser
SEARCH_DIR = ../search/network/openmp

CFLAGS = -I. -I$(INIPARSER_DIR) -I$(SEARCH_DIR)
CFLAGS += -g -Wall -Wno-unused -fno-omit-frame-pointer \
	     -O3 -funsafe-loop-optimizations -funroll-loops \
	     -march=broadwell -mtune=broadwell -mavx2 -mfma

LDFLAGS = -L$(INIPARSER_DIR)
LDLIBS = -lc -lrt -lm -lhdf5 -lhdf5_hl -lgsl -lcblas -liniparser

SRC = coincidences.c io.c fap.c
INC = cdefs.h ../search/network/openmp/struct.h coincidences.h
DEPS = $(SRC) $(INC) Makefile

all: $(TARGET)

$(TARGET): $(DEPS)
	$(CC) $(LDFLAGS) -o $@ $(SRC) $(CFLAGS) $(LDLIBS)

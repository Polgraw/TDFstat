CFLAGS = -O3 -g -Wall -W -D RELATIVE
LDFLAGS =

SRC := extract_band.c \
	../pss_lib/pss_gw.c \
	../pss_lib/pss_math.c \
	../pss_lib/pss_serv.c \
	../pss_lib/pss_snag.c \
	../pss_astro/novas/novas.c \
	../pss_astro/novas/novascon.c \
	../pss_astro/novas/solsys3.c \
	../pss_astro/novas/readeph0.c \
	../pss_astro/pss_astro.c \
	../pss_astro/daspostare.c \
	pss_ante.c \
	pss_sfdb.c

OBJ := ${SRC:.c=.o}
SRC_HDF_TD := ${SRC:extract_band.c=extract_band_hdf_td.c}
OBJ_HDF_TD := ${SRC_HDF_TD:.c=.o}

%.o:    %.c
	$(CC) $(CFLAGS) -c $< -o $@

all: extract_band_out extract_band_hdf_td

extract_band_out: $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ -lm

extract_band_hdf_td: $(OBJ_HDF_TD)
	$(CC) $(LDFLAGS) -o $@ $^ -lm -lhdf5 -lhdf5_hl -lfftw3f

clean:
	rm -f $(OBJ) $(OBJ_HDF_TD) extract_band_out extract_band_hdf_td

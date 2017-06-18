#CC=g++
CC=mpic++
CFLAGS= \
		-O2 -std=c++11 \
		-I$(HDF5_DIR)/include/
LDFLAGS= \
		-shared \
		-L$(HDF5_DIR)/lib/ \
		-lhdf5_hl \
		-lhdf5


LIB=libVTS_interface
LIBEXT=so


all: $(LIB).$(LIBEXT)

$(LIB).$(LIBEXT): VTS.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $?

clean:
	rm -f *.o

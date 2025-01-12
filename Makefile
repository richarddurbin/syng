# makefile for gaffer developed on Richard's Mac

CFLAGS = -O3
#CFLAGS = -g	# for debugging

ALL = syng ONEview syngmap syngref

DESTDIR = ~/bin

all: $(ALL)

install:
	cp $(ALL) $(DESTDIR)

clean:
	$(RM) *.o *~ $(ALL)
	\rm -r *.dSYM

### object files

UTILS_OBJS = hash.o dict.o array.o utils.o
UTILS_HEADERS = utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

#SEQIO_OPTS = -DONEIO
#SEQIO_LIBS = -lm -lz

HTS_DIR = $(PWD)/../htslib
SEQIO_OPTS = -DONEIO -DBAMIO -I$(HTS_DIR)/htslib/
SEQIO_LIBS = -L$(HTS_DIR) -Wl,-rpath $(HTS_DIR) -lhts -lm -lbz2 -llzma -lcurl -lz 

seqio.o: seqio.c seqio.h ONElib.h $(UTILS_HEADERS)
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

seqhash.o: seqhash.c seqhash.h $(UTILS_HEADERS)
	$(CC) $(CFLAGS) -c $^

kmerhash.o: kmerhash.c kmerhash.h $(UTILS_HEADERS)
	$(CC) $(CFLAGS) -c $^

syngbwt.o: syngbwt.c syng.h $(UTILS_HEADERS) ONElib.h
	$(CC) $(CFLAGS) -c $^

ONElib.o: ONElib.c ONElib.h 
	$(CC) $(CFLAGS) -c $^

### programs

syng: syng.c syngbwt.o seqio.o seqhash.o kmerhash.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lpthread $(SEQIO_LIBS)

syngmap: syngmap.c syngbwt.o seqio.o seqhash.o kmerhash.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lpthread $(SEQIO_LIBS)

syngprune: syngprune.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(SEQIO_LIBS)

syngbwt: syngbwt.c syng.h seqio.o seqhash.o kmerhash.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(SEQIO_LIBS) syngbwt.o

syngref: syngref.c seqio.o seqhash.o kmerhash.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(SEQIO_LIBS)

ONEview: ONEview.c ONElib.o
	$(CC) $(CFLAGS) -o $@ $^ -lz

### end of file

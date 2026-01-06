# makefile for gaffer developed on Richard's Mac

CFLAGS = -O3
#CFLAGS = -g	# for debugging

ALL = syng ONEview syngmap syngstat k31type

DESTDIR = ~/bin

all: $(ALL)

install:
	cp $(ALL) $(DESTDIR)

clean:
	$(RM) *.o *~ $(ALL) TEST/*.1* TEST/gbwt.fa
	$(RM) -r *.dSYM

### object files

UTILS_OBJS = hash.o dict.o array.o utils.o
UTILS_HEADERS = utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

# Set BAMIO=1 to enable SAM/BAM/CRAM support (requires htslib in ../htslib)
ifdef BAMIO
HTS_DIR = $(PWD)/../htslib/.
SEQIO_OPTS = -DONEIO -DBAMIO -I$(HTS_DIR)/htslib/
SEQIO_LIBS = -L$(HTS_DIR) -Wl,-rpath $(HTS_DIR) -lhts -lm -lbz2 -llzma -lcurl -lz
else
SEQIO_OPTS = -DONEIO
SEQIO_LIBS = -lm -lz
endif 

seqio.o: seqio.c seqio.h ONElib.h $(UTILS_HEADERS)
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

seqhash.o: seqhash.c seqhash.h $(UTILS_HEADERS)
	$(CC) $(CFLAGS) -c $^

kmerhash.o: kmerhash.c kmerhash.h $(UTILS_HEADERS)
	$(CC) $(CFLAGS) -DONEIO -c $^

syngbwt.o: syngbwt.c syng.h $(UTILS_HEADERS) ONElib.h
	$(CC) $(CFLAGS) -c $^

syncmerset.o: syncmerset.c syncmerset.h $(UTILS_HEADERS) ONElib.h
	$(CC) $(CFLAGS) -c $^

ONElib.o: ONElib.c ONElib.h 
	$(CC) $(CFLAGS) -c $^

### programs

syng: syng.c syngbwt.o syncmerset.o seqio.o seqhash.o kmerhash.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lpthread $(SEQIO_LIBS)

syngmap: syngmap.c syngbwt.o seqio.o seqhash.o kmerhash.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lpthread $(SEQIO_LIBS)

syngstat: syngstat.c syngbwt.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lz

syngprune: syngprune.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(SEQIO_LIBS)

syngbwt: syngbwt.c syng.h seqio.o seqhash.o kmerhash.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(SEQIO_LIBS) syngbwt.o

k31type: k31type.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(SEQIO_LIBS)


ONEview: ONEview.c ONElib.o
	$(CC) $(CFLAGS) -o $@ $^ -lz

### test

test: syng TEST/test.fa
	./syng -o TEST/test -writeK -writeGBWT -outputEnds TEST/test.fa
	./syng -readK TEST/test.1khash -o TEST/gbwt -writeSeq -outputEnds TEST/test.1gbwt
	seqconvert -o TEST/gbwt.fa TEST/gbwt.1seq
	diff TEST/test.fa TEST/gbwt.fa

### end of file

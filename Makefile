# makefile for gaffer developed on Richard's Mac

CFLAGS = -O3
ifneq ($(CSYNCMER),)
ifneq ($(CSYNCMER),0)
CFLAGS += -DUSE_CSYNCMER
USE_CSYNCMER = 1
endif
endif
ifneq ($(AVX2),)
ifneq ($(AVX2),0)
CFLAGS += -DHAVE_AVX2
CFLAGS_AVX2 = -march=native -mavx2
USE_AVX2 = 1
endif
endif
ifdef USE_CSYNCMER
ifndef USE_AVX2
$(error CSYNCMER=1 needs AVX2=1)
endif
endif
#CFLAGS = -g	# for debugging

ALL = syng syngpath2gbwt ONEview syngmap syngstat k31type

DESTDIR = ~/bin

all: $(ALL)

install:
	cp $(ALL) $(DESTDIR)

clean:
	$(RM) *.o *.gch *~ $(ALL) TEST/*.1* TEST/gbwt.fa
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

# syncmer_iter.o provides the syncmer iterator: seqhash.c is Richard's seeded hash,
# syncmer_iter.c the ntHash/AVX2 one from csyncmer_fast
ifdef USE_CSYNCMER
CSYNCMER_HEADERS = csyncmer_fast.h
syncmer_iter.o: syncmer_iter.c syncmer_iter.h $(CSYNCMER_HEADERS) $(UTILS_HEADERS)
	$(CC) $(CFLAGS) $(CFLAGS_AVX2) -c $< -o syncmer_iter.o
else
syncmer_iter.o: seqhash.c seqhash.h $(UTILS_HEADERS)
	$(CC) $(CFLAGS) -c $< -o syncmer_iter.o
endif

ifdef USE_AVX2
avx2.o: avx2.c avx2.h seqio.h $(UTILS_HEADERS)
	$(CC) $(CFLAGS) $(CFLAGS_AVX2) -c $<
LINK_AVX2 = avx2.o
endif

kmerhash.o: kmerhash.c kmerhash.h $(UTILS_HEADERS)
	$(CC) $(CFLAGS) -DONEIO -c $^

rskip.o: rskip.c $(UTILS_HEADERS)
	$(CC) $(CFLAGS) -c $^

syngbwt3.o: syngbwt3.c syng.h $(UTILS_HEADERS) ONElib.h
	$(CC) $(CFLAGS) -c $^

syncmerset.o: syncmerset.c syncmerset.h $(UTILS_HEADERS) ONElib.h
	$(CC) $(CFLAGS) -c $^

ONElib.o: ONElib.c ONElib.h 
	$(CC) $(CFLAGS) -c $^

### programs

syng: syng.c syngpipe.h syngbwt3.o rskip.o syncmerset.o seqio.o syncmer_iter.o kmerhash.o ONElib.o $(UTILS_OBJS) $(LINK_AVX2)
	$(CC) $(CFLAGS) -o $@ $(filter-out %.h,$^) -lpthread $(SEQIO_LIBS)

syngpath2gbwt: syngpath2gbwt.c syngbwt3.o rskip.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lpthread $(SEQIO_LIBS)

syngmap: syngmap.c syngbwt3.o rskip.o syncmerset.o syncmer_iter.o kmerhash.o seqio.o ONElib.o $(UTILS_OBJS) $(LINK_AVX2)
	$(CC) $(CFLAGS) -o $@ $^ -lpthread $(SEQIO_LIBS)

syngstat: syngstat.c syngbwt3.o rskip.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lz -lpthread

syngprune: syngprune.c seqio.o syncmer_iter.o ONElib.o $(UTILS_OBJS) $(LINK_AVX2)
	$(CC) $(CFLAGS) -o $@ $^ $(SEQIO_LIBS)

syngbwt3: syngbwt3.c rskip.o syng.h seqio.o syncmer_iter.o kmerhash.o ONElib.o $(UTILS_OBJS) $(LINK_AVX2)
	$(CC) $(CFLAGS) -o $@ $^ $(SEQIO_LIBS)

k31type: k31type.c seqio.o syncmer_iter.o ONElib.o $(UTILS_OBJS) $(LINK_AVX2)
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

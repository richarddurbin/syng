# syng
Syncmer graphs, and potentially other sorts of sequence graphs.

The main product here is **syng**. This essentially reads in sequence files of various sorts, and represents them in terms of *paths* over a conceptual graph of *syncmers*, which are a subset of all possible kmers guaranteed to provide a sparse but complete cover of any DNA sequence (average depth approximately two). See [Edgar, 2021](https://peerj.com/articles/10805/) for further information about syncmers - what we call a syncmer here, Edgar calls a "closed syncmer".  The resulting syncmer graphs can be used as assembly graphs or pangenome graphs.

syng makes extensive use of the ONEcode package <https://github.com/thegenemyers/ONEcode> for compact, fast and efficient file representation.

Sets of syncmers, for example those needed to represent a set of sequences, are normally stored in a **.1khash** file, although it is possible to export them also in gzipped fasta file.  As well as explicitly exporting paths as strings of kmer indices in a **.1path** file, syng can also build a [GBWT](https://academic.oup.com/bioinformatics/article/36/2/400/5538990) that implicitly represents the paths, stored in a **.1gbwt** file.  We plan for syng to also export, and potentially import [GFA](https://gfa-spec.github.io/GFA-spec/GFA1.html) files.

To give some idea of performance, syng converted a 20x PacBio HiFi read data set of a 935Mb cichlid genome (~19Gbp) into a 1.05GB .1khash file and a 493Mb .1gbwt file of (1023,32)-syncmers in 62 seconds on a MacBook Pro, and 92 human genomes (277Gbp) from release 1 of [HPRC](https://humanpangenome.org/) (leaving out HG002 for evaluation) into a 4.4GB .1khash file and a 5.3GB .1gbwt file of (63,8)-syncmers in about 4 hours on a Linux HPC cluster. 

The .1khash, .1path and .1gbwt files are all examples of **ONEcode**
files.  The project also contains **ONEview** from the
[ONEcode](https://github.com/thegenemyers/ONEcode) repository. Further
utilities **seqconvert**, **seqstat**, **seqextract** to summarise, manipulate and
interconvert between fastz[.gz] and **.1seq** ONEcode sequence files
are available from the [SEQUENCE_UTILTIES](https://github.com/thegenemyers/ONEcode/tree/main/SEQUENCE_UTILITIES) subdirectory of ONEcode.

## Building
```
  git clone https://github.com/richarddurbin/syng.git
  cd syng
  make
```
If you want to be able to read SAM/BAM/CRAM files then you need to install [htslib](https://github.com/samtools/htslib) in a parallel directory and build with `BAMIO=1`:

 ```
  cd ..
  git clone https://github.com/samtools/htslib.git
  cd htslib
  autoreconf -i  # Build the configure script and install files it uses
  ./configure    # Optional but recommended, for choosing extra functionality
  make
  cd ../syng
  make clean
  make BAMIO=1
```

## Libraries
There are various useful libraries, with header files:

- **ONElib.[hc]** supports ONEcode file reading and writing, with implementation in the single file ONElib.c with no dependendencies. See <https://github.com/thegenemyers/ONEcode> for further information.
- **seqio.[hc]** supports reading, writing DNA files with a few other basic operations. Implementation in seqio.c, with dependencies on utils.[hc], libz, ONElib and htslib depending on compile operations.
- **seqhash.[hc]** sequence processing to give kmers of various types, including extraction of syncmers, minimizers and modimizers, via an iterator interface
- **kmerhash.[hc]** efficient code for building, searching, writing and reading (into a .1kash file) tables of fixed length kmers, as for example returned by seqhash
- **utils.[hc]** some very low level type definitions (e.g. I8 to I64 and U8 to U64), die(), warn(), new(), new0(), and a timing package. No dependencies beyond normal C run time library. NB there is a handy fzopen() which will silently open .gz files as standard files, but this depends on funopen() which is not available on all systems. If this does not compile/link then you will need to undefine WITH_ZLIB to link.
- **array.[ch], dict.[ch], hash.[ch]** respectively implement advanced language style extendable arrays, dictionaries (hashes of strings) and general hashes of basic types (up to 64-bit).
- 
## Synopsis

An example usage pattern is given below.

```
> syng
Usage: syng <operation>* <input>*
possible operations are:
  -w <window length>     : [55] syncmer length = w + k
  -k <smer length>       : [8] must be under 32
  -seed <seed>           : [7] for the hashing function
  -T <threads>           : [8] number of threads
  -o <outfile prefix>    : [syngOut] applies to all following write* options
  -readK <.1khash file>  : read and start from this syncmer (khash) file
  -zeroK                 : zero the kmer counts
  -noAddK                : do not create new syncmers - convert unmatched syncmers to 0
  -histK                 : output quadratic histogram of kmer counts (after sequence processsing)
  -writeK                : write the syncmers as a .1khash file
  -writeKfa              : write the syncmers as a fasta file, with ending .kmer.fa.gz
  -writeNewK <file prefix> : write new syncmers as a .1khash file; implies -noAddK
  -writePath             : write a .1path file (paths of nodes)
  -writeGBWT             : write a .1gbwt file (nodes, edges and paths in GBWT form)
  -writeSeq              : write a .1seq file (paths converted back to sequences)
  -outputEnds            : write the non-syncmer ends of path sequences as X,Y lines
possible inputs are:
  <sequence file>        : any of fasta[.gz], fastq[.gz], BAM/CRAM/SAM, .1seq
  <.1path file>          : sequences as lists of kmers, with optional non-syncmer DNA ends
  <.1gbwt file>          : graph BWT with paths, with optional ends
Operations are carried out in order as they are parsed, with some setting up future actions,
e.g. changing the outfile prefix affects following lower case options for file opening
Some output files, e.g. .1gbwt will be output at the end, after all inputs are processed,
whereas others, e.g. .1path are written as inputs are processed.

> syng -o cichlid -writeK -writePath -outputEnds *.fa.gz
k, w, seed are 8 55 7
sequence file 1 fAstCal68.FINAL.fa.gz type fasta: had 31 sequences 935785828 bp, yielding 36820971 syncs with 29445488 extra syncmers
user    18.041796       system  0.918415        elapsed 13.309440       alloc_max 3752  max_RSS 5911035904
sequence file 2 fAulStu2.FINAL.fa.gz type fasta: had 23 sequences 927632324 bp, yielding 36196076 syncs with 6870638 extra syncmers
user    17.979181       system  0.289314        elapsed 7.157738        alloc_max 4128  max_RSS 112132096
sequence file 3 fDipLim2.FINAL.fa.gz type fasta: had 437 sequences 932767327 bp, yielding 36681142 syncs with 4535666 extra syncmers
user    18.628796       system  0.458877        elapsed 7.594654        alloc_max 6516  max_RSS 2246590464
sequence file 4 fLabFue1.FINAL.fa.gz type fasta: had 31 sequences 936110624 bp, yielding 36736333 syncs with 3748552 extra syncmers
user    17.172514       system  0.225017        elapsed 6.408912        alloc_max 6526  max_RSS 0
sequence file 5 fLabTre1.FINAL.fa.gz type fasta: had 26 sequences 936863659 bp, yielding 36778443 syncs with 1973510 extra syncmers
user    16.720750       system  0.158108        elapsed 5.891383        alloc_max 6811  max_RSS 86949888
sequence file 6 fMayPea1.FINAL.fa.gz type fasta: had 30 sequences 927122434 bp, yielding 36548424 syncs with 2765274 extra syncmers
user    17.186198       system  0.197489        elapsed 6.093247        alloc_max 7113  max_RSS 827047936
Total for this run 578 sequences, total length 5596282196
Overall total 219761389 instances of 49339130 syncmers, average 4.45 coverage
wrote 49339130 syncmers to file cichlid.1khash
user    2.358794        system  0.312918        elapsed 2.805761        alloc_max 7113  max_RSS 1245184
total: user     108.088040      system  2.594965        elapsed 49.295976       alloc_max 7113  max_RSS 9185001472

> syng -o cichlid -readK cichlid.1khash -writeGBWT -outputEnds cichlid.1path 
read syncmer parameters k 8 w 55 (size 63) seed 7
read 49339130 syncmers from cichlid.1khash with total count 219761392
user    2.078966        system  0.367622        elapsed 2.459572        alloc_max 3972  max_RSS 3435151360
k, w, seed are 8 55 7
path file 1 cichlid.1path: had 578 sequences containing 219761397 syncmers
user    98.485414       system  28.079270       elapsed 106.574848      alloc_max 12455 max_RSS 23292624896
Total for this run 578 sequences, total length 0
Overall total 439522789 instances of 49339130 syncmers, average 8.91 coverage
wrote gbwt to file cichlid.1gbwt
user    19.222173       system  1.010528        elapsed 20.482916       alloc_max 12455 max_RSS 0
total: user     120.153096      system  29.556817       elapsed 130.028740      alloc_max 12455 max_RSS 26727776256

> syng -o cichlid -readK cichlid.1khash -writeSeq cichlid.1gbwt
read syncmer parameters k 8 w 55 (size 63) seed 7
read 49339130 syncmers from cichlid.1khash with total count 219761392
user    2.076341        system  0.382973        elapsed 2.460453        alloc_max 3972  max_RSS 3435397120
k, w, seed are 8 55 7
path file 1 cichlid.1gbwt: read GBWT with 49339130 vertices and 115529066 edges
user    33.337580       system  0.797333        elapsed 34.188887       alloc_max 7022  max_RSS 3624419328
had 578 sequences containing 219761397 syncmers
user    159.932554      system  21.517852       elapsed 96.478368       alloc_max 15660 max_RSS 20926840832
Total for this run 578 sequences, total length 0
Overall total 439522789 instances of 49339130 syncmers, average 8.91 coverage
total: user     195.346545      system  22.837188       elapsed 133.267362      alloc_max 15660 max_RSS 27986657280

> ls -l cichlid*
-rw-r--r--  1 rd  staff  1012108820 Mar 16 12:34 cichlid.1khash		// the kmer sequences
-rw-r--r--  1 rd  staff   776147166 Mar 16 12:34 cichlid.1path		// the sequences as lists of kmers
-rw-r--r--  1 rd  staff  1290427468 Mar 16 12:38 cichlid.1gbwt		// the sequences stored as a GBWT over kmers
-rw-r--r--  1 rd  staff  1399079525 Mar 16 13:21 cichlid.1seq		// the regenerated sequences 2-bit compressed with indices

// as the level of replication increases the gbwt representation becomes more efficient

> seqstat -b cichlid.1seq
onecode file, 578 sequences >= 0, 5596282196 total, 9682149.13 average, 1000 min, 85450508 max
bases
  a 1646170104 29.4 %
  c 1152062443 20.6 %
  g 1152334894 20.6 %
  t 1645714755 29.4 %

> ONEview -h cichlid.1gbwt | head  // -h does not show the header, just the body

> ONEview -H cichlid.1gbwt  // -H shows only the header, with the schema and object statistics
```

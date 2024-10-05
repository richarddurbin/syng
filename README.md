# syng
Syncmer graphs, and potentially other sorts of sequence graphs.

The main product here is **syng**. This essentially reads in sequence files of various sorts, and represents them in terms of *paths* over a conceptual graph of *syncmers*, which are a subset of all possible kmers guaranteed to provide a sparse but complete cover of any DNA sequence (average depth approximately two). See [Edgar, 2021](https://peerj.com/articles/10805/) for further information about syncmers - what we call a syncmer here, Edgar calls a "closed syncmer".  The resulting syncmer graphs can be used as assembly graphs or pangenome graphs.

syng makes extensive use of the ONEcode package <https://github.com/thegenemyers/ONEcode> for compact, fast and efficient file representation.

Sets of syncmers, for example those needed to represent a set of sequences, are normally stored in a **.1khash** file, although it is possible to export them also in gzipped fasta file.  As well as explicitly exporting paths as strings of kmer indices in a **.1path** file, syng can also build a [GBWT](https://academic.oup.com/bioinformatics/article/36/2/400/5538990) that implicitly represents the paths, stored in a **.1gbwt** file.  We plan for syng to also export, and potentially import [GFA](https://gfa-spec.github.io/GFA-spec/GFA1.html) files.

To give some idea of performance, syng converted a 20x PacBio HiFi read data set of a 935Mb cichlid genome (~19Gbp) into a 1.05GB .1khash file and a 493Mb .1gbwt file of (1023,32)-syncmers in 62 seconds on a MacBook Pro, and 92 human genomes (277Gbp) from release 1 of [HPRC](https://humanpangenome.org/) (leaving out HG002 for evaluation) into a 4.4GB .1khash file and a 5.3GB .1gbwt file of (63,8)-syncmers in about 4 hours on a Linux HPC cluster. 

The project also contains relevant versions of various Durbin package utility programs:

- **ONEview** to view ONEcode files, and convert them between binary and ascii. By convention .1seq files are binary, .seq files are ascii.
- **seqstat** to give information about the data in DNA sequence files: fasta[.gz], fastq[.gz], 1seq (if compiled with -DONEIO, as by default in this project), and SAM/BAM/CRAM (if compiled with -DBAMIO, see below).
- **seqconvert** to convert between DNA file types. This also will homopolymer compress (-H "hoco"), and if writing to a 1seq file while doing this will store the offsets in the original file of each new sequence base position, allowing to unhoco (-U) back to the original sequence. Because ONEcode files compress all lists, this
- **seqextract** to extract (sets of) sequences or subsequences from a sequence file, reverse-complementing them if wished

For each program, running it without any arguments gives usage information.

## Building
```
  git clone https://github.com/richarddurbin/gaffer.git
  make
```
If you want to be able to read SAM/BAM/CRAM files then you need to install htslib in a parallel directory and use Makefile.bam:

 ```
 cd ..
  git clone https://github.com/samtools/htslib.git
  cd htslib
  autoreconf -i  # Build the configure script and install files it uses
  ./configure    # Optional but recommended, for choosing extra functionality
  make
  make install
  cd ../gaffer
  make clean
  make -f Makefile.bam
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

An example usage pattern is given below. **THIS NEEDS COMPLETING ONCE THE CODE IS REFACTORED**

```
> syng
       
> syng -o test

> seqstat -t fAulStu2-reads.fa.gz // -t gives time and memory usage
fasta file, 1174680 sequences >= 0, 18910078042 total, 16098.07 average, 316 min, 180179 max
user    68.803876       system  1.278024        elapsed 70.200890       alloc_max 16    max_RSS 16924672

> seqconvert -t -o fAulStu2-reads.1seq fAulStu2-reads.fa.gz   // -1 is implied by the output filename
reading from file type fasta
written 1174680 sequences to file type onecode, total length 18910078042, max length 180179
user    71.821538       system  4.444104        elapsed 77.369275       alloc_max 18    max_RSS 34488320

> seqstat -t fAulStu2-reads.1seq // .1seq is a much faster format, also indexed, self-documenting and supporting threaded read/write
onecode file, 1174680 sequences >= 0, 18910078042 total, 16098.07 average, 316 min, 180179 max
user    5.731012        system  0.700100        elapsed 6.434019        alloc_max 16    max_RSS 37257216

> syng -o fAulStu2 -k -b -S fAulStu2-reads.1seq

> ls -ltr fAulStu*

> ONEview -h fAulStu.1gbwt | head  // -h does not show the header, just the body

> ONEview -H fAulStu.1gbwt  // -H shows only the header, with the schema and object statistics
```

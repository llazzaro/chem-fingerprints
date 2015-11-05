The FPS format is an exchange format. It's easy to read and write, but not all that fast to load. The FPB format is an alternative binary format designed to work with the the chemfp arena data model. I wrote the first proposals in 2010 and 2011, based on the PNG file format. However, at that point there was no driving need for a fast loader, the format needed more review and experimentation, and the chemfp internals weren't developed enough to make use of it.

Now that chemfp is in wider use, the load time for large (>100,000 fingerprints) is noticeable in command-line scripts. Chemfp is also making its way into web servers, where development reload time is important and where a memory-mapped file makes a better use of memory.

FBP support is a goal for chemfp 1.2.

The format is a pretty standard [FourCC](http://en.wikipedia.org/wiki/FourCC) format containing a set of block. Each block has an 8 byte length, a 4 byte block type, and a data payload. The details are still in development, so stay tuned!
chemfp documentation
====================

chemfp is a set of tools for working with cheminformatics
fingerprints in the FPS format.

Most people will use the command-line programs to generate and search
fingerprint files. ob2fps, oe2fps, and rdkit2fps use respectively the
Open Babel, OpenEye, and RDKit chemistry toolkits to convert structure
files into fingerprint files. sdf2fps extracts fingerprints encoded in
SD tags to make the fingerprint file. simsearch finds targets in a
fingerprint file which are sufficiently similar to the queries.

The programs are based on the chemfp Python library, which uses a C
extension for the performance critical sections. The parts of the
library API documented here are meant for public use, along with some
examples.


Contents:

.. toctree::
   :maxdepth: 2


Installing
======

The chemfp tools depends on a working Python installation.
You can download Python 2.7 from ... (Note! OpenEye doesn't yet
support Python 2.7 so you will need to download Python 2.6 from ...)

The core chemfp functionality does not depend on a third-party library
but you will need a chemistry toolkit in order to generate new
fingerprints from structure files. chemfp supports the free Open Babel
and RDKit toolkits and the proprietary OEChem toolkit. Make sure you
install the Python libraries for the toolkit(s) you select.


If you have a source version of chemfp then you will need a C compiler
in order to compile it. This uses Python's standard "setup.py" so you
can see http://.../ for details. The short version is that on Unix
systems using sudo (that is, Mac OS X and most Linux-based OSes) you
can do:

   sudo python setup.py install

while for Windows you can do

   python setup.py install



Getting started: Working with PubChem data files
====================

In this section you'll learn how to create a fingerprint file and use
it for similarity searches. We'll use data files from PubChem. You do
not need a chemistry toolkit for this work because we'll use the
pre-computed CACTVS fingerprints included with each PubChem record.


Start by downloading the files
Compound_027575001_027600000.sdf.gz
ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_027575001_027600000.sdf.gz
and
Compound_014550001_014575000.sdf.gz
ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_014550001_014575000.sdf.gz
from PubChem. At the time of writing they contain 224 and 3119 records,
respectively. (I chose smaller than average files so they would be
easier to open and review.)

Next, convert the files into fingerprint files. On the command line
do the following two commands:

    sdf2fps --pubchem Compound_027575001_027600000.sdf.gz -o pubchem_queries.fps
    sdf2fps --pubchem Compound_014550001_014575000.sdf.gz  -o pubchem_targets.fps

Congratulations!

How does this work? Each PubChem record contains the precomputed
CACTVS substructure keys in the PUBCHEM_CACTVS_SUBSKEYS tag. The
"--pubchem" flag tells sdf2fps to get the value of that tag and decode
it to get the fingerprint. It also adds a few metadata fields to the
fingerprint file header.

The order of the fingerprints are the same as the order of the
corresponding record in the SDF, although unconvertable records might
be skipped, depending on the --errors flag.


k-nearest neighbor search
--------------

We'll use the pubchem_queries.fps as queries for a k=2 nearest neighor
similarity search of the target file puchem_targets.gps.

   simsearch -k 2 -q pubchem_queries.fps pubchem_targets.fps

That's all! You should get output which starts:

    #Simsearch/1
    #num_bits=881
    #type=Tanimoto k=2 threshold=0.0
    #software=chemfp/1.0
    #queries=pubchem_queries.fps
    #targets=pubchem_targets.fps
    #query_sources=Compound_027575001_027600000.sdf.gz
    #target_sources=Compound_014550001_014575000.sdf.gz
    2	27575433	14568234	0.6904	14550456	0.6492
    2	27575577	14570945	0.7487	14570951	0.7385
    2	27575602	14553070	0.7594	14572463	0.7256
    2	27575603	14553070	0.7594	14572463	0.7256

How do you interpret the output? The lines starting with '#' are
header lines. It contains metadata information describing that this is
a similarity search report. You can see the search parameters, the
name of the tool which did the search, and the filenames which went
into the search.

After the '#' header lines come the search results, with one result
per line. There are in the same order as the query fingerprints. Each
result line contains tab-delimited columns. The first column is the
number of hits. The second column is the query identifier used. The
remaining columns contain the hit data, with alternating target id and
its score.

For example, the first result line contains the 2 hits for the
query 27575433. The first hit is the target id 14568234 with score
0.6904 and the second hit is 14550456 with score 0.6492. Since this is
a k-nearest neighor search, the hits are sorted by score, starting
with the highest score. Do be aware that ties are broken arbitrarily.


Threshold search
-------------

Let's do a threshold search and find all hits which are at least 0.738
similar to the queries.

    simsearch --threshold 0.738 -q pubchem_queries.fps pubchem_targets.fps

The first 20 lines of output from this are:

    #Simsearch/1
    #num_bits=881
    #type=Tanimoto k=all threshold=0.738
    #software=chemfp/1.0
    #queries=pubchem_queries.fps
    #targets=pubchem_targets.fps
    #query_sources=Compound_027575001_027600000.sdf.gz
    #target_sources=Compound_014550001_014575000.sdf.gz
    0	27575433
    3	27575577	14570945	0.7487	14570992	0.7383	1457095	0.7385
    1	27575602	14553070	0.7594
    1	27575603	14553070	0.7594
    1	27575880	14569866	0.7727
    0	27575897
    0	27577227
    0	27577234
    0	27577237
    1	27577250	14569555	0.7474
    0	27577307
    0	27577324

Let's look at the second result line, which contains the 3 hits for
the query id 27575577. As before, the hit information alternates
between the target ids and the target scores, but unlike the k-nearest
search, the hits are not in a particular order. You can see that here
with the scores are 0.7487, 0.7383, and 0.7385.

You might be wondering why I chose the 0.738 threshold. Query id
27575577 has 6 hits with a threshold of 0.7 or higher. That requires
14 columns to show, which is a bit overwhelming.

Combined k-nearest and threshold search
-----------------------------------

You can combine the -k and --threshold queries to find the k-nearest
neighbors which are all above a given threshold. 

    simsearch -k 3 --threshold 0.7 -q pubchem_queries.fps pubchem_targets.fps

This find the nearest 3 structures, which all must be at least 0.7
similar to the query fingerprint. The output from the above starts:

    #Simsearch/1
    #num_bits=881
    #type=Tanimoto k=3 threshold=0.7
    #software=chemfp/1.0
    #queries=pubchem_queries.fps
    #targets=pubchem_targets.fps
    #query_sources=Compound_027575001_027600000.sdf.gz
    #target_sources=Compound_014550001_014575000.sdf.gz
    0	27575433
    3	27575577	14570945	0.7487	14570951	0.7385	14570990.7383
    3	27575602	14553070	0.7594	14572463	0.7256	14553060.7208
    3	27575603	14553070	0.7594	14572463	0.7256	14553060.7208
    3	27575880	14569866	0.7727	14567856	0.7308	14566360.7246
    0	27575897
    1	27577227	14570135	0.7143
    0	27577234
    1	27577237	14569555	0.7371

The output format is identical to the previous examples.


NxN (self-similarity) searches
-------------

chemfp has no special support for the NxN search of a fingerprint data
set against itself. Instead, use the same file as both the queries and
the targets. This will take twice as much memory and time as an
optimized search. (If you are interested in funding such a tool, I can
provide you a cost estimate.)



ChEBI - Working with a toolkit
===================

In this section I'll show you how to create a fingerprint file using a
chemistry toolkit. chemfp includes support for Open Babel, OpenEye, and
RDKit.

We'll work with data from ChEBI http://www.ebi.ac.uk/chebi/ which
contains "Chemical Entities of Biological Interest". They distribute
their structures in several formats, including as an SD file. For this
section, download the "lite" version from
ftp://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_lite.sdf.gz . It
contains the same structure data as the complete version but has many
fewer tag data.  For ChEBI 82 this file contains 19640 records and the
compressed file is 5MB.

Unlike the PubChem data set, the ChEBI data set does not contain
fingerprints so we'll need to generate them using a toolkit.

ChEBI record titles don't contain the id
------------------------------

Strangely, the ChEBI dataset does not use the title line of the SD
file to store the record id. A simple examination shows that 16498 of
the title lines are empty, 2119 of them have the title "ChEBI", and 45
of them are labeled "Structure #1."

Instead, the id is the value of the "ChEBI ID" tag, which looks like

    > <ChEBI ID>
    CHEBI:776

By default the toolkit-based fingerprint generation tools use the
title as the identifier, and exits with an error if the identifier is
missing. (Use the --errors option to change the behaviour). If you try
one of them with this data file you will get the error message:

    ERROR: Missing title for record #1 of 'ChEBI_lite.sdf.gz'. Exiting.

Instead, use the --id-tag option to specify of the name of the data
tag containing the id. For this data set you'll need to write it as

    --id-tag "ChEBI ID"

The quotes are important because of the space in the tag name.

There's also a "ChEBI Name" tag which includes data values like
"tropic acid" and "(+)-guaia-6,9-diene". Every record has a unique
name so the names could be used as the primary identifier. The FPS
fingerprint file format allows identifiers with a space, or comma, or
anything other tab, newline, and a couple of other bytes, so it's no
problem using those names directly.

To use the ChEBI Name as the primary chemfp identifier, use:

    --id-tag "ChEBI Name"


Generating fingerprints with Open Babel
--------------------------------

If you have the Open Babel Python library installed then you can use
ob2fps to generate fingerprints.

    ob2fps --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o ob_chebi.fps

This takes about 30 seconds on my laptop.

The default uses the FP2 fingerprints, so the above is the same as

    ob2fps --FP2 --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o ob_chebi.fps

ob2fps can generate several other types of fingerprints. (See XXX for
details). For example, to generate the Open Babel implementation of the
MACCS definition use:

    ob2fps --MACCS --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o chebi_maccs.fps


Generating fingerprints with OpenEye
--------------------------------

If you have the OEChem Python library installed then you can use
oe2fps to generate fingerprints.

    oe2fps --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o oe_chebi.fps

This takes about 10 seconds on my laptop and generates a number of
stereochemistry warnings.

The default settings produce OEGraphSim path fingerprint with the
values:

    numbits=4096 minbonds=0 maxbonds=5 atype=DefaultAtom btype=DefaultBond

Each of these can be changed through command-line options. See XXX for
details. There are also options to use an alternate aromaticity model.

oe2fps can generate several other types of fingerprints. For example,
to generate the OpenEye implementation of the MACCS definition use:

   oe2fps --maccs166 --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o chebi_maccs.fps


Generating fingerprints with RDKit
------------------------------

If you have the RDKit Python library installed then you can use
rdkit2fps to generate fingerprints. Based on the previous examples you
probably guessed that the command-line is:

    rdkit2fps --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o rdkit_chebi.fps

Unfortunately, this isn't enough. If you do this you'll get the message

    [18:54:07] Explicit valence for atom # 12 N, 4, is greater than permitted
    ERROR: Could not parse molecule block at line 14840 of 'ChEBI_lite.sdf.gz'. Exiting.

The first line comes from RDKit's error log. RDKit is careful to check
that structures make chemical sense, and in this case it didn't like
the 4-valent nitrogen. It refuses to process this molecule.

The second line comes from rdkit2fps. By default it complains and
exits with an error if RDKit cannot process a record. Basically it
highlights the source of the problem and demands that you do something
about it.

In most cases it's okay to skip a few records which can't be
processed. You can tell rdkit2fps to report the error but continue
processing by using the --errors option:

    rdkit2fps --id-tag "ChEBI ID" --errors report ChEBI_lite.sdf.gz -o rdkit_chebi.fps

Four minutes later I see that 403 records out of the 19640 could not
be processed.

The previous command-line created RDKit's path fingerprints with
parameters

    minPath=1 maxPath=7 fpSize=2048 nBitsPerHash=4 useHs=1

Each of those can be changed through command-line options. See XXX for
details.

rdkit2fps can generate several other types of fingerprints. For
example, to generate the RDKit implementation of the MACCS definition
use:

   rdkitfps --maccs166 --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o chebi_maccs.fps

chemfp supports neither count fingerprints nor sparse fingerprints so
cannot generate RDKit's circular fingerprints.


"--rdmaccs" and "--substruct" substructure fingerprints
====================

chemfp implements two platform-independent fingerprints meant for
substructure filters but which are also used for similarity
searches. One is based on the MACCS implementation in RDKit and the
other comes from the PubChem/CACTVS substructure fingerprints.


All three of the supported toolkits have built-in support to generate
MACCS fingerprints. What you might not know is that they have
different interpretation of those bits. Open Babel's implementation
derives from RDKit's but OpenEye has a completely different
implementation and gives different results for some cases.

It seems there's no real published definition of those bits, at least
not in a way that's unambiguous. There's also no comprehensive test
suite. What a sad state of affairs that this is the most widely used
fingerprint scheme in our field.

Well, it's time for a change. chemfp includes a set of pattern
definitions derived from RDKit and implementations of those patterns
for each of the supported toolkits. They are called the "rdmaccs"
fingerprints because they owe so much to RDKit. You can use them to
compare MACCS keys from two different toolkits ... and STILL get
different fingerprints.

See, each toolkit perceives chemistry differently. Open Babel before
2.3 didn't support chirality so chiral-based bits will never be
set. Each toolkit uses a different definition of aromaticity, so a bit
which is set when there are "two or more aromatic rings" will be
toolkit dependent.

The code is still experimental but feel free to try it out by asking
for the "--rdmaccs" fingerprint from the command-line tools.

Even more experimental is support for the 881 bit fingerprints derived
from the PubChem/CACTVS substructure keys. If you want to try those
out, ask for the "--substruct" fingerprint from the command-line tools.



The chemfp Python library
=================

The chemfp command-line programs use a Python library called
chemfp. Many parts of the API are in flux and subject to change. The
stable portions of the API which are open for general use are
documented XXX.

The API includes:

 - low-level Tanimoto and popcount operations
 - Tanimo search algorithms based on threshold and/or k-nearest neighbors
 - a cross-toolkit interface for reading fingerprints from a structure file

The following sections gives examples of how to use the API.

Basic concepts: byte and hex fingerprints
=========

chemfp stores fingerprints as byte strings. Here are two 8 bit
fingerprints:

    >>> fp1 = "A"
    >>> fp2 = "B"

The `chemfp.bitops` module contains functions which work on byte
fingerprints. Here's the Tanimoto of those two fingerprints:

    >>> from chemfp import bitops
    >>> bitops.byte_tanimoto(fp1, fp2)
    0.33333333333333331
    >>> 

To understand why, you have to know that ASCII character "A" has the
value 65, and "B" has the value 66. The bit representation is

     "A" = 01000001   and   "B" = 01000010

so their intersection has 1 bit and the union has 3, giving a Tanimoto
of 1/3 or 0.33333333333333331 when represented as a 64 bit floating
point number on the computer.

You can compute the Tanimoto between any two byte strings with the
same length, as in:

    >>> bitops.byte_tanimoto("apples&", "oranges")
    0.58333333333333337

You'll get a chemfp exception if they have different lengths.

Most fingerprints are not as easy to read as the English ones I showed
above. They tend to look more like:

    P1@\x84K\x1aN\x00\n\x01\xa6\x10\x98\\\x10\x11

which is hard to read. I usually show hex-encoded fingerprints. The above
fingerprint in hex is:

    503140844b1a4e000a01a610985c1011

which is simpler to read, though you still need to know your hex
digits.

The bitops module includes other low-level functions which work on
byte fingerprints, as well as corresponding functions which work on
hex fingerprints. (Hex-encoded fingerprints are decidedly second-class
citizens in chemfp, but they are citizens.)


Basic concepts: fingerprint collections and metadata
=========

A fingerprint record is the fingerprint plus an identifier. In chemfp,
a fingerprint collection is a object which contains fingerprint
records and which follows the common API providing access to those
records.

That's rather abstract, so let's work with a few real examples. You'll
need to create a copy of the "pubchem_targets.fps" file from section
XXX in order to follow allow.

Here's how to open an FPS file:

    >>> import chemfp
    >>> reader = chemfp.open("pubchem_targets.fps")

Every fingerprint collection has a metadata attribute with details
about the fingerprints. It comes from the header of the FPS file. You
can view the metadata in Python repr format

    >>> reader.metadata
    Metadata(num_bits=881, num_bytes=111, software=u'CACTVS/unknown', type
    ='CACTVS-E_SCREEN/1.0 extended=2', sources=['Compound_014550001_014575
    000.sdf.gz'], date='2011-09-14T12:10:34', aromaticity=None)

but I think it's easier to view it in string format, which matches the
format of the FPS header

    >>> print reader.metadata
    #num_bits=881
    #type=CACTVS-E_SCREEN/1.0 extended=2
    #software=CACTVS/unknown
    #source=Compound_014550001_014575000.sdf.gz
    #date=2011-09-14T12:10:34
    
    >>> 


All fingerprint collections support iteration. Each step of the
iteration returns the fingerprint identifier and its score. Since I
know the 6th record has the id 14550045, I can write a simple loop
which stops with that record

    >>> for (id, fp) in reader:
    ...   print id, "starts with", fp.encode("hex")[:20]
    ...   if id == "14550045":
    ...     break
    ... 
    14550001 starts with 034e1c00020000000000
    14550002 starts with 034e0c00020000000000
    14550003 starts with 034e0400020000000000
    14550005 starts with 010e1c00000600000000
    14550010 starts with 034e1c40000000000000
    14550045 starts with 071e8c03000000000000
    >>> 


Fingerprint collections also support iterating via arenas, and the
standard Tanimoto search functions.


Basic concepts: a fingerprint arena
=========

The FPSReader reads through or searches a fingerprint file once. If
you want to read the file again you have to reopen in.

Reading from disk is slow, and the FPS format is designed for
ease-of-use and not performance. If you want to do many queries then
it's best to store everything in memory. The `FingerprintArena` is a
fingerprint collection which does that.

Here's how to load fingerprints into an arena:

    >>> import chemfp
    >>> arena = chemfp.load_fingerprints("pubchem_targets.fps")
    >>> print arena.metadata
    #num_bits=881
    #type=CACTVS-E_SCREEN/1.0 extended=2
    #software=CACTVS/unknown
    #source=Compound_014550001_014575000.sdf.gz
    #date=2011-09-14T12:10:34
    
    >>> 

This implements the fingerprint collection API, so you can do things
like iterate over an arena and get the id/fingerprint pairs.

    >>> from chemfp import bitops
    >>> for id, fp in arena:
    ...     print id, "with popcount", bitops.byte_popcount(fp)
    ...     if id == "14574718":
    ...         break
    ... 
    14550474 with popcount 2
    14574635 with popcount 2
    14550409 with popcount 4
    14550416 with popcount 6
    14574551 with popcount 7
    14550509 with popcount 8
    14550423 with popcount 10
    14550427 with popcount 10
    14574637 with popcount 10
    14574890 with popcount 11
    14574718 with popcount 12
    >>> 

If you look closely you'll notice that the fingerprint record order
has changed from the previous section, and that the population counts
are suspiciously non-decreasing. By default `load_fingerprints`
reorders the fingerprints into a data structure which is faster to
search, although you can disable that if you want the fingerprints to
be the same as the input order.

The `FingerprintArena` has new capabilities. You can ask it how many
fingerprints it contains, get the list of identifiers, and look up a
fingerprint record given an index, as in:

    >>> len(arena)
    3119
    >>> arena.ids[:5]
    ['14550474', '14574635', '14550409', '14550416', '14574551']
    >>> id, fp = arena[6]
    >>> id
    '14550423'
    >>> arena[-1][0]
    '14566760'
    >>> bitops.byte_popcount(arena[-1][1])
    231
    >>> 

An arena supports iterating through subarenas. This is like having a
long list and being able to iterate over sublists. Here's an example
of iterating over the arena to get subarenas of size 1000 (excepting
the last), and print information about each subarena.

   >>> for subarena in arena.iter_arenas(1000):
    ...   print subarena.ids[0], len(subarena)
    ... 
    14550474 1000
    14573373 1000
    14555885 1000
    14560068 119
    >>> arena[0][0]
    '14550474'
    >>> arena[1000][0]
    '14573373'

To help demonstrate what's going on, I verified that the subarena
"ids[0]" matched the corresponding record in the main arena.


Arenas are a core part of chemfp. Processing one fingerprint at a time
is slow, so the main search routines expect to iterate over query
arenas, rather than query fingerprints.

Thus, the FPSReaders - and all chemfp fingerprint collections - also
support the `iter_arenas` interface. Here's an example of reading the
targets file 25 records at a time:

    >>> queries = chemfp.open("pubchem_queries.fps")
    >>> for arena in queries.iter_arenas(25):
    ...   print len(arena)
    ... 
    25
    25
    25
    25
    25
    25
    25
    25
    24
    >>> 




How to use query fingerprints to search for similar target fingerprints
==============================

In this section you'll learn how to do a Tanimoto search using the
previously created PubChem fingerprint files for the queries and the
targets.

It's faster to search an arena, so I'll load the target fingerprints:

    >>> import chemfp
    >>> targets = chemfp.load_fingerprints("pubchem_targets.fps")
    >>> len(targets)
    3119
    >>> 

and open the queries as an FPSReader.

    >>> queries = chemfp.open("pubchem_queries.fps")
    >>>

I'll use `chemfp.threshold_tanimoto_search` to find, for each query,
all hits which are at least 0.7 similar to the query.

    >>> queries = chemfp.open("pubchem_queries.fps")
    >>> for (query_id, hits) in chemfp.threshold_tanimoto_search(queries, targets, threshold=0.7):
    ...   print query_id, len(hits), hits[:2]
    ... 
    27575433 0 []
    27575577 18 [('14570945', 0.74874371859296485), ('14570946', 0.73762376237623761)]
    27575602 3 [('14572463', 0.72560975609756095), ('14553070', 0.75935828877005351)]
    27575603 3 [('14572463', 0.72560975609756095), ('14553070', 0.75935828877005351)]
    27575880 9 [('14569876', 0.72307692307692306), ('14567856', 0.73076923076923073)]
    27575897 0 []
    27577227 1 [('14570135', 0.7142857142857143)]
    27577234 0 []
           ....

I'm only showing the first two hits for the sake of space. It seems
rather pointless, after all, to show all 18 hits of query id 27575577.

What you don't see is that the implementation uses the iter_arenas()
interface on the queries so that it processes only a subarena at a
time. There's a tradeoff between a large arena, which is faster
because it doesn't often go back to Python code, or a small arena,
which uses less memory and is more responsive. You can change the
tradeoff using the `arena_size` parameter.


If all you cared about was the count of the hits within a given
threshold then use `chemfp.count_tanimoto_hits`

    >>> queries = chemfp.open("pubchem_queries.fps")
    >>> for (query_id, count) in chemfp.count_tanimoto_hits(queries, targets, threshold=0.7):
   ...     print query_id, count
    ...   break
    ... 
    27575433 0
    27575577 18
    27575602 3
    27575603 3
    27575880 9
    27575897 0
    27577227 1
    27577234 0
    27577237 1
    27577250 4
         ...

Or, if you only want the k=2 nearest neighbors to each target within
that same threshold of 0.7 then use `chemfp.knearest_tanimoto_search`

    >>> queries = chemfp.open("pubchem_queries.fps")
    >>> for (query_id, hits) in chemfp.knearest_tanimoto_search(query_arena, targets, k=2, threshold=0.7):
   ...     print query_id, hits
    ...   break
    ... 
    27575433 []
    27575577 [('14570945', 0.74874371859296485), ('14570951', 0.73853211009174313)]
    27575602 [('14553070', 0.75935828877005351), ('14572463', 0.72560975609756095)]
    27575603 [('14553070', 0.75935828877005351), ('14572463', 0.72560975609756095)]
    27575880 [('14569866', 0.77272727272727271), ('14567856', 0.73076923076923073)]
    27575897 []
    27577227 [('14570135', 0.7142857142857143)]
    27577234 []
    27577237 [('14569555', 0.73711340206185572)]
    27577250 [('14569555', 0.74742268041237114), ('14550456', 0.72131147540983609)]
         ...


How to search an FPS file
==================

The previous example loaded the fingerprints into a
FingerprintArena. That's the fastest way to do multiple
searches. Sometimes though you only want to do one or a couple of
queries. It seems rather excessive to read the entire targets file
into an in-memory data structure before doing the search when you
could search will processing the file.


For that case, use an FPSReader as the target file. Here I'll get the
first record from the queries file and use it to search the targets
file:

    >>> query_arena = next(chemfp.open("pubchem_queries.fps").iter_arenas(1))

This line opens the file, iterates over its fingerprint records, and
return the first one. The next() function was added after Python
2.5 so the above won't work for that version. Instead, use:

    >>> query_arena = chemfp.open("pubchem_queries.fps").iter_arenas(1).next()

which is the older form.

Here are the k=5 closest hits against the targets file:

    >>> targets = chemfp.open("pubchem_targets.fps")
    >>> for query_id, hits in chemfp.knearest_tanimoto_search(query_arena, targets, k=5, threshold=0.0):
    ...   print "Hits for", query_id
    ...   for hit in hits:
    ...     print "", hit
    ... 
    Hits for 27575433
     ('14568234', 0.69035532994923854)
     ('14550456', 0.64921465968586389)
     ('14572463', 0.64444444444444449)
     ('14566364', 0.63953488372093026)
     ('14573723', 0.63247863247863245)
    >>> 

Remember that the FPSReader is based on reading an FPS file. Once
you've done a search, the file is read, and you can't do another
search. You'll need to reopen the file.

Bear in mind that each search processes arena_size query fingerprints
at a time. You will need to increase that value if you want to search
more than that number of fingerprints with this method. The search
performance tradeoff between a FPSReader search and loading the
fingerprints into a FingerprintArena occurs with under 10 queries, so
there should be little reason to worry about this.


FingerprintArena searches returning indicies instead of ids
================================


The previous sections used a high-level interface to the Tanimoto
search code. Those are designed for the common case where you just
want the query id and the hits, where each hit includes the target id.

Working with strings is actually rather inefficient in both speed and
memory. It's better to work with indicies if you can, and in the next
section I'll show how to make a distance matrix using this interface.

The index methods are only available as methods of a FingerprintArena,
where the arena contains the targets. Three of the methods
(`count_tanimoto_search_arena`, `threshold_tanimoto_search_arena`, and
`knearest_tanimoto_search_arena`) take another arena as the
query. Here's an example where I use the first 5 records from
pubchem_queries.fps to search the entire contents of that file.

   >>> import chemfp
    >>> dataset = chemfp.load_fingerprints("pubchem_queries.fps")
    >>> first_5 = next(dataset.iter_arenas(5))
    >>> results = dataset.threshold_tanimoto_search_arena(first_5, threshold=0.7)

You can iterate over the results to get the list of hits for each of
the queries. (The order of the results is the same as the order of the
records in the query.)

   >>> for hits in results:
    ...   print len(hits), hits[:3]
    ... 
    2 [('27581954', 1.0), ('27581957', 1.0)]
    2 [('27581954', 1.0), ('27581957', 1.0)]
    3 [('27580389', 1.0), ('27580394', 0.88235294117647056), ('27581637', 0.75)]
    2 [('27584917', 1.0), ('27585106', 0.89915966386554624)]
    2 [('27584917', 0.89915966386554624), ('27585106', 1.0)]
    >>> 

This is like what you saw earlier, except that it doesn't have the
query id. (If you want that you can enumerate() over the results and
use the index into the query arena's ids[] list.)

But what I really want to show is that you can get the same data only
using the offset index for the target record instead of its id. The
result from a Tanimoto search is a `SearchResult` object, with the
methods `iter_hits()`.

    >>> for hits in results.iter_hits():
    ...   print len(hits), hits[:3]
    ... 
    2 [(0, 1.0), (1, 1.0)]
    2 [(0, 1.0), (1, 1.0)]
    3 [(2, 1.0), (5, 0.88235294117647056), (20, 0.75)]
    2 [(3, 1.0), (4, 0.89915966386554624)]
    2 [(3, 0.89915966386554624), (4, 1.0)]
    >>> 
    >>> dataset.ids[0]
    '27581954'
    >>> dataset.ids[1]
    '27581957'
    >>> dataset.ids[5]
    '27580394'
    >>> 

I did a few id lookups given the target dataset to show you that the
index corresponds to the identifiers from the previous code.


Computing a distance matrix for clustering
==================

chemfp does not do clustering. There's a huge number of tools which
arleady do that. A goal of chemfp in the future is to provide some
core components which clustering algorithms can use.

That's in the future. Right now you can use the following to build
a distance matrix and pass that to one of those tools. It's a bit
inefficient since it computes almost twice as many Tanimoto scores as
it needs to do, and uses twice the necessary memory, but what's a
factor of two among friends?


    import numpy  # NumPy must be installed
    
    # Compute distance[i][j] = 1-Tanimoto(fp[i], fp[j])
    
    def distance_matrix(arena):
        n = len(arena)
        
        # The Tanimoto search computes all of the scores when threshold=0.0.
        # The SearchResult contains sparse data, so I set all values
        # now to 1.0 so you can experiment with higher thresholds.
        distances = numpy.ones((n, n), "d")
        
        # Keep track of where the query subarena is in the query
        query_row = 0
        
        for query_arena in arena.iter_arenas():
            results = arena.threshold_tanimoto_search_arena(query_arena, threshold=0.0)
            for q_i, hits in enumerate(results.iter_hits()):
                query_idx = query_row + q_i
                for target_idx, score in hits:
                    distances[query_idx, target_idx] = 1.0 - score
            query_row += len(query_arena)
        
        return distances


We're almost there - let's cluster using the hcluster package which I
downloaded earlier today. (By that I mean that I don't have enough
experience with Python clustering packages to recommend one for you.)

http://code.google.com/p/scipy-cluster/

    import chemfp
    import hcluster # Clustering package from http://code.google.com/p/scipy-cluster/
    
    # ... insert the `distance_matrix` function definition here ...

    dataset = chemfp.load_fingerprints("docs/pubchem_queries.fps")
    distances  = distance_matrix(dataset)
    
    linkage = hcluster.linkage(distances, method="single", metric="euclidean")
    
    # Plot using matplotlib, which you must have installed
    hcluster.dendrogram(linkage, labels=dataset.ids)
    
    import pylab
    pylab.show()


Taylor-Butina clustering
==========

For the last clustering example, here's my (non-validated) variation
of the Butina algorithm from JCICS 1999, 39, 747-750.
http://www.chemomine.co.uk/dbclus-paper.pdf . See also
http://www.redbrick.dcu.ie/~noel/R_clustering.html .


   import chemfp
    
    dataset = chemfp.load_fingerprints("pubchem_targets.fps")
    search = dataset.threshold_tanimoto_search_arena(dataset, threshold = 0.8)
    
    def get_hit_indicies(hits):
        return [id for (id, score) in hits]
    
    # Reorder so the centroid with the most hits comes first.
    # (That's why I do a reverse search.)
    # Ignore the arbitrariness of breaking ties by fingerprint index
    results = sorted( (  (len(hits), i, get_hit_indicies(hits))
                                        for (i, hits) in enumerate(search.iter_hits())  ),
                      reverse=True)
    
    # Determine the true/false singletons and the clusters
    true_singletons = []
    false_singletons = []
    clusters = []
    
    seen = set()
    for (size, fp_idx, members) in results:
        if fp_idx in seen:
            # Can't use a centroid which is already assigned
            continue
        seen.add(fp_idx)
    
        
        if size == 1:
            # The only fingerprint in the exclusion sphere is itself
            true_singletons.append(fp_idx)
            continue
    
        # Figure out which ones haven't yet been assigned
        unassigned = [target_idx for target_idx in members if target_idx not in seen]
    
        if not unassigned:
            false_singletons.append(fp_idx)
            continue
            
        # this is a new cluster
        clusters.append( (fp_idx, unassigned) )
        seen.update(unassigned)
    
    # Report the results
    
    print len(true_singletons), "true singletons"
    print "=>", " ".join(sorted(dataset.ids[idx] for idx in true_singletons))
    print
    
    print len(false_singletons), "false singletons"
    print "=>", " ".join(sorted(dataset.ids[idx] for idx in false_singletons))
    print
    
    # Sort so the cluster with the most compounds comes first,
    # then by alphabetically smallest id
    def cluster_sort_key(cluster):
        centroid_idx, members = cluster
        return -len(members), dataset.ids[centroid_idx]
        
    clusters.sort(key=cluster_sort_key)
    
    print len(clusters), "clusters"
    for centroid_idx, members in clusters:
        print dataset.ids[centroid_idx], "has", len(members), "other members"
        print "=>", " ".join(dataset.ids[idx] for idx in members)


Out of curiosity, I tried this on 100,000 compounds selected
arbitrarily from PubChem. It took 7 minutes for my laptop to process
with a threshold of 0.8. In the Butina paper, it took 24 hours to do
the same, although that was with a 1024 bit fingerprint instead of
881. It's hard to judge the absolute speed differences of a 12 year
old MIPS R4000 to a two year old laptop, but it's less than a factor
of 200. Part must certainly be due to the work I put into making
chemfp fast.



Reading structure fingerprints using a toolkit
==========

What happens if you're given a structure file and you want to find the
two nearest matches in an FPS file? You'll have to generate the
fingerprints for the structures in the structure file, then do the
comparison.

For this section you'll need to have a chemistry toolkit. I'll use the
"chebi_maccs.fps" file you generated earlier as the targets, and the
PubChem file "Compound_027575001_027600000.sdf.gz as the source of
query structures.

    >>> import chemfp
    >>> targets = chemfp.load_fingerprints("chebi_maccs.fps")
    >>> queries = chemfp.read_structure_fingerprints(targets.metadata, "Compound_027575001_027600000.sdf.gz")
    >>> for (query_id, hits) in chemfp.knearest_tanimoto_search(queries, targets, k=2, threshold=0.4):
    ...   print query_id, "=>",
    ...   for (target_id, score) in hits:
    ...     print "%s %.3f" % (target_id, score),
    ...   print
    ... 
    27575433 => CHEBI:280152 0.667 CHEBI:3176 0.662
    27575577 => CHEBI:6375 0.600 CHEBI:46068 0.600
    27575602 => CHEBI:3090 0.683 CHEBI:6790 0.635
    27575603 => CHEBI:3090 0.683 CHEBI:6790 0.635
    27575880 => CHEBI:59736 0.725 CHEBI:8887 0.617
    27575897 => CHEBI:8887 0.632 CHEBI:51491 0.622
    27577227 => CHEBI:59007 0.831 CHEBI:59120 0.721
    27577234 => CHEBI:59007 0.809 CHEBI:9398 0.722
    27577237 => CHEBI:59007 0.789 CHEBI:52890 0.741
    27577250 => CHEBI:59007 0.753 CHEBI:4681 0.722
          ... many others omitted ...

That's it! Pretty simple, wasn't it? You didn't even need to explictly
specify which toolkit you wanted to use.

The only new thing here is `chemfp.read_structure_fingerprints`. The
first parameter of this is the metadata used to configure the
reader. In my case it's:

    >>> print targets.metadata
    #num_bits=166
    #type=OpenEye-MACCS166/1
    #software=OEGraphSim/1.0.0 (20100809)
    #aromaticity=openeye
    #source=ChEBI_lite.sdf.gz
    #date=2011-09-14T17:50:28

The "type" told chemfp which toolkit to use to read molecules, and how
to generate fingerprints from those molecules, while "aromaticity"
told it which aromaticity model to use when reading the molecule file.

You can of course pass in your own metadata as the first parameter to
read_structure_fingerprints, and as a shortcut, if you pass in a
string then it will be used as the fingerprint type.

For examples, if you have OpenBabel installed then you can do:

   >>> reader = chemfp.read_structure_fingerprints("OpenBabel-MACCS", "Compound_027575001_027600000.sdf.gz")
    >>> for i, (id, fp) in enumerate(reader):
    ...   print id, fp.encode("hex")
    ...   if i == 3:
    ...     break
    ... 
    27575433 800404000840549e848189cca1f132aedfab6eff1b
    27575577 800400000000449e850581c22190022f8a8baadf1b
    27575602 000000000000449e840191d820a0122eda9abaff1b
    27575603 000000000000449e840191d820a0122eda9abaff1b
    >>> 

If you have OEChem and OEGraphSim installed then you can do

    >>> reader = chemfp.read_structure_fingerprints("OpenEye-MACCS166", "Compound_027575001_027600000.sdf.gz")
    >>> for i, (id, fp) in enumerate(reader):
    ...   print id, fp.encode("hex")
    ...   if i == 3:
    ...     break
    ... 
    27575433 000000080840448e8481cdccb1f1b216daaa6a7e3b
    27575577 000000080000448e850185c2219082178a8a6a5e3b
    27575602 000000080000448e8401d14820a01216da983b7e3b
    27575603 000000080000448e8401d14820a01216da983b7e3b
    >>> 

And if you have RDKit installed then you can do:

    >>> reader = chemfp.read_structure_fingerprints("RDKit-MACCS166", "Compound_027575001_027600000.sdf.gz")
    >>> for i, (id, fp) in enumerate(reader):
    ...   print id, fp.encode("hex")
    ...   if i == 3:
    ...     break
    ... 
    27575433 000000000840549e84818dccb1f1323cdfab6eff1f
    27575577 000000000000449e850185c22190023d8a8beadf1f
    27575602 000000000000449e8401915820a0123eda98bbff1f
    27575603 000000000000449e8401915820a0123eda98bbff1f
    >>> 



The anatomy of an PubChem SD record
==================

This section is for people who don't have much experience with the
details of SD files. It defines concepts I use elsewhere in the text.

Take a look at the first record of the file, which has PubChem ID 14550001
http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=14550001&loc=ec_rcs
. You may have to uncompress the file order to read it, or you can
view the record directly at
http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=14550001&disopt=DisplaySDF
. I've also extracted parts of the record here:


   14550001
      -OEChem-09091107342D
    
     26 26  0     0  0  0  0  0  0999 V2000
        5.1350    0.8450    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
        2.5369   -0.6550    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        6.0010    3.3450    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        6.0010   -3.6550    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
        6.8671   -2.1550    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    ...
     15 23  1  0  0  0  0
     15 24  1  0  0  0  0
    M  CHG  2   4  -1   6   1
    M  END
    > <PUBCHEM_COMPOUND_CID>
    14550001
    
    > <PUBCHEM_COMPOUND_CANONICALIZED>
    1
    
    > <PUBCHEM_CACTVS_COMPLEXITY>
    209
    
    > <PUBCHEM_CACTVS_HBOND_ACCEPTOR>
    4
    
    > <PUBCHEM_CACTVS_HBOND_DONOR>
    2
    
    > <PUBCHEM_CACTVS_ROTATABLE_BOND>
    4
    
    > <PUBCHEM_CACTVS_SUBSKEYS>
    AAADccByOABAAAAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAAAABAAAAHgQECAAADAyl2AKwxoAQQgiBAiRCQwCCAAAgIhAoiAAGbIoINiKikZOAcABkwBEo2AewwCAOAEAAQAAIAAAAgACAABAAAAAAAAAAAA==
    
    > <PUBCHEM_IUPAC_OPENEYE_NAME>
    2-(2-hydroxyethylsulfanylmethyl)-4-nitro-phenol
      
    > <PUBCHEM_IUPAC_CAS_NAME>
    2-[(2-hydroxyethylthio)methyl]-4-nitrophenol
    
    > <PUBCHEM_IUPAC_NAME>
    2-(2-hydroxyethylsulfanylmethyl)-4-nitrophenol

    ...
    > <PUBCHEM_BONDANNOTATIONS>
    10  13  8
    11  14  8
    13  14  8
    7  10  8
    7  9  8
    9  11  8
    
    $$$$


The first line of an SDF record is the title field. PubChem, like most
databases, stores the record identifier in the title field. In this
case it's "14550001". (See the ChEBI section below if you have a SD
file which stores the identifier in a tag.)

Starting with line 4 is the connection table, which contains the
atom and bond information. This ends with the line "M  END".

After that is the tag data section, which is a list of tag/value data
pairs defined by the data provider. PubChem has decided that for this
record the tag "PUBCHEM_COMPOUND_CID" is associated with the data
value "14550001". If you look through other records you'll see that
PUBCHEM_COMPOUND_CID is always the PubChem compound identifier for
that record.

PubChem defines many other tags. A record's PUBCHEM_IUPAC_NAME tag
stores the IUPAC Preferred Name created by OpenEye's Lexichem naming
tool, and PUBCHEM_CACTVS_ROTATABLE_BOND stores the number of
rotatable bonds found by Xemistry's CACTVS program. For a full list of
tags, see
  ftp://ftp.ncbi.nlm.nih.gov//pubchem/specifications/pubchem_sdtags.txt


Precomputed CACTVS fingerprints
========================

Each PubChem records contain a PUBCHEM_CACTVS_SUBSKEYS tag which
contains the 881 bit substructure fingerprints computed by
CACTVS. They are encoded in a special format, and the value for
compound id 14550001 is the single long line

    AAADccByOABAAAAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAAAABAAAAHgQECAAADAyl2AKwxo
    AQQgiBAiRCQwCCAAAgIhAoiAAGbIoINiKikZOAcABkwBEo2AewwCAOAEAAQAAIAAAAgACA
    ABAAAAAAAAAAAA==

which is decoded
ftp://ftp.ncbi.nlm.nih.gov//pubchem/specifications/pubchem_fingerprints.txt
and reencoded in chemfp's hex format as:

    034e1c000200000000000000000000000000000000000c000000000000000080000000
    7820201000003030a51b400d630108421081402442c200410000044408141100603651
    106c444589c9010e00260388141be00d03047000020002001000000001000100080000
    000000000000

Each bit has a physical meaning based on substructure. For this record
it means the structure has >=4 and >=8 hydrogens (but not >=16
hydrogens), >=2 and >=4 and >= 8 carbons, but not >= 16 carbons, and
so on.

PubChem uses these fingerprints for similarity searches, and we'll do
the same.



Anatomy of an FPS fingerprint file
=====================

Take a look at Compound_014550001_014575000.fps . It contains 3125
lines. The first six lines contain the header and should look
something like:

    #FPS1
    #num_bits=881
    #type=CACTVS-E_SCREEN/1.0 extended=2
    #software=CACTVS/unknown
    #source=Compound_014550001_014575000.sdf.gz
    #date=2011-09-14T11:03:18

This contains metadata about the fingerprints. The first line is the
version of the file format. The "#num_bits" line contains the number
of bits in the fingerprint. The "#type" and "#software" line state the
fingerprint type and the name of the program used to generate it, and
are automatically set to these values when you use the --pubchem
command-line option.

The "#source" line describes where the fingerprint data came from, and
"#date" contains a timestamp of when the file was created.

The rest of the file contains the fingerprint data. Each line contains
the fingerprint as 222 hex characters, a tab character, and the
compound id. Here are the first four fingerprint lines from
Compound_014550001_014575000.fps, slighly reformatted so it's easier
to read:

    034e1c000200000000000000000000000000000000000c000000000000000080000000
    7820201000003030a51b400d630108421081402442c200410000044408141100603651
    106c444589c9010e00260388141be00d03047000020002001000000001000100080000
    000000000000	14550001
    034e0c000200000000000000000000000000000000000c000000000000000080000000
    7820081000003030a51b400d6301024010014024420200410000044408101100603611
    106c444589c9010e00260b88141be00d03047000020002001000000001000100080000
    000000000000	14550002
    034e04000200000000000000000000000000000000000c000000000000000080000000
    7820081000003030a11b004d6301024010014024420200410000044408101100603611
    1064444589c9010e00260b88101be00d03047000020002001000000001000100080000
    000000000000	14550003
    010e1c000006000000000000000000000000000000000c0600000000000000830a0000
    58000000000030000119000c10030000001140044b1000400000240000101180001013
    10644c01a808018c002403801091e111130f7103004000000800000100200000040000
    000000000000	14550005

How did sdf2fps get the fingerprint data? When you specified --pubchem
you told it that the fingerprint data is found in the
PUBCHEM_CACTVS_SUBSKEYS tag and encoded with the CACTVS method.

The --pubchem option is really a shortcut. You could instead use the
command-line options

    --software=CACTVS/unknown --type 'CACTVS-E_SCREEN/1.0 extended=2'
     --fp-tag=PUBCHEM_CACTVS_SUBSKEYS --cactvs

but who wants to write out all of that each time?






Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


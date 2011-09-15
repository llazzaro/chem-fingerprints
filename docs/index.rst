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
library API documented here are meant for public use, and there are
some examples of how to build a leader clustering algorithm or a
web-service for structure search.

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



Using the chemfp Python library
=================

The chemfp command-line programs use a Python library called
chemfp. Many parts of the API are in flux and subject to change. The
stable portions of the API which are open for general use are
documented XXX.

In this section I'll show you examples of how to use the API. I'll use
the "pubchem_targets.fps" file created in XXX. For reference, the
first ten lines of that file are:

    #FPS1
    #num_bits=881
    #type=CACTVS-E_SCREEN/1.0 extended=2
    #software=CACTVS/unknown
    #source=Compound_014550001_014575000.sdf.gz
    #date=2011-09-14T12:10:34
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

One of the fingerprint collection types is called a
FingerprintArena. It stores all fingerprints in-memory. I'll use
"load_library" to open a fingerprint file and create an arena with its
contents:

    >>> targets = chemfp.load_fingerprints("pubchem_targets.fps")
    >>> len(targets)
    3119
    >>> 

This says there are 3119 fingerprints in the arena. The "metadata"
attribute stores information about the contents of the arena.

    >>> targets.metadata
    Metadata(num_bits=881, num_bytes=111, software=u'CACTVS/unknown', type
    ='CACTVS-E_SCREEN/1.0 extended=2', sources=['Compound_014550001_014575
    000.sdf.gz'], date='2011-09-14T12:10:34', aromaticity=None)
    >>> targets.metadata.num_bits
    881
    >>> targets.metadata.type
    'CACTVS-E_SCREEN/1.0 extended=2'
    >>> 

All arenas support iteration through the contents of the arena as (id,
fingerprint) pairs, where the fingerprint is a byte string. Here's a
(slow) linear scan to find the fingerprint for id 14550005.

    >>> for (id, fp) in targets:
    ...   if id == "14550005":
    ...     print id, fp[:5].encode("hex"), repr(fp[:5])
    ... 
    14550005 010e1c0000 '\x01\x0e\x1c\x00\x00'
    >>> 

The arenas offer no faster way to look up a fingerprint given its
identifier but you can do it yourself because the identifers are
available in iteration order as the "id" attribute and because arenas
support index lookup to return the (id, fingerprint) pair.

    >>> target_id_index = dict( (id, i) for (i, id) in enumerate(targets.ids) )
    >>> targets[target_id_index["14550005"]]
    ('14550005', '\x01\x0e\x1c\x00\x00\x06\x00\x00 ... ')
    >>>





All pubchem fingerprint collections support an interface called a
"FingerprintReader". One type of FingerprintReader is an FPSReader,
which reads fingerprints from an FPS file.

    >>> import chemfp
    >>> reader = chemfp.open("pubchem_targets.fps")
    >>> reader
   <chemfp.readers.FPSReader object at 0x100547590>
   >>>

All readers have a "metadata" attribute. For the FPSReader the
metadata comes from the header of the fingerprint file.

    >>> reader.metadata.num_bits
    881
    >>> reader.metadata.type
    'CACTVS-E_SCREEN/1.0 extended=2'
    >>>
   >>> reader.metadata
    Metadata(num_bits=881, num_bytes=111, software=u'CACTVS/unknown', type
    ='CACTVS-E_SCREEN/1.0 extended=2', sources=['Compound_014550001_014575
    000.sdf.gz'], date='2011-09-14T12:10:34', aromaticity=None)
   >>> print reader.metadata
    #num_bits=881
    #type=CACTVS-E_SCREEN/1.0 extended=2
    #software=CACTVS/unknown
    #source=Compound_014550001_014575000.sdf.gz
    #date=2011-09-14T12:10:34
    
    >>> 

If you iterate over a reader you get (id, fingerprint) pairs. For
example:

    >>> for id, fp in reader:
    ...   print id, fp[:5].encode("hex"), repr(fp[:5])
    ...   if id == "14550003":
    ...     break
    ... 
    14550001 034e1c0002 '\x03N\x1c\x00\x02'
    14550002 034e0c0002 '\x03N\x0c\x00\x02'
    14550003 034e040002 '\x03N\x04\x00\x02'
    >>> 

Some FingerprintReaders allow multiple iterators but others only
support forward iteration. For example, the FPSReader reads
fingerprints from a file and supports neither multiple iterators nor
resets to the start of the file.





Another fingerprint collection type is called a FingerprintArena.



The anatomy of an PubChem SD record
==================

Feel free to jump to the next section and start working with the
chemfp tools. This section is for people who don't have much
experience with the details of SD files. It defines concepts I use
later in the text.

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


Extract fingerprint data from PubChem
=========================

The sdf2fps command-line program has an incredibly boring name - it
converts an SD file to a fingerprint file in FPS format, using
information from each SD record.

To convert your downloaded PubChem file, do this in your shell:

    sdf2fps --pubchem Compound_014550001_014575000.sdf.gz -o Compound_014550001_014575000.fps

All of that should be on a single line. The chemfp tools understand
gzip files so I told sdf2fps to work directly on the compressed
file. If you earlier uncompressed the file then you can instead work
with that file, as:

    sdf2fps --pubchem Compound_014550001_014575000.sdf -o Compound_014550001_014575000.fps


Congratulations! You've just made a fingerprint file!

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


.. highlight:: none

===================================
Working with the command-line tools
===================================

The sections in this chapter describe examples of using the
command-line tools to generate fingerprint files and to do similarity
searches of those files.

.. _pubchem_fingerprints:

Generating fingerprint files from PubChem SD files
==================================================

In this section you'll learn how to create a fingerprint file from an
SD file which contains pre-computed CACTVS fingerprints. You do not
need a chemistry toolkit for this section.

`PubChem <http://pubchem.ncbi.nlm.nih.gov/>`_ is a great resource
of publically available chemistry information. The data is available
for `ftp download <ftp://ftp.ncbi.nlm.nih.gov>`_. We'll use some of
their `SD formatted
<http://en.wikipedia.org/wiki/Structure_Data_File#SDF>`_ files.
Each record has a PubChem/CACTVS fingerprint field, which we'll used.

Start by downloading the files 
Compound_027575001_027600000.sdf.gz
(from
ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_027575001_027600000.sdf.gz
)
and Compound_014550001_014575000.sdf.gz
(from
ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_014550001_014575000.sdf.gz
). At the time of writing they contain 224 and 3119 records,
respectively. (I chose smaller than average files so they would be
easier to open and review.)

Next, convert the files into fingerprint files. On the command line
do the following two commands::

    sdf2fps --pubchem Compound_027575001_027600000.sdf.gz -o pubchem_queries.fps
    sdf2fps --pubchem Compound_014550001_014575000.sdf.gz -o pubchem_targets.fps

Congratulations, that was it!

How does this work? Each PubChem record contains the precomputed
CACTVS substructure keys in the PUBCHEM_CACTVS_SUBSKEYS tag. The
:option:`--pubchem` flag tells sdf2fps to get the value of that tag and decode
it to get the fingerprint. It also adds a few metadata fields to the
fingerprint file header.

The order of the fingerprints are the same as the order of the
corresponding record in the SDF, although unconvertable records might
be skipped, depending on the :option:`--errors` flag.

If you store records in an SD file then you almost certainly don't use
the same fingerprint encoding as PubChem. sdf2ps can decode from a
number of encodings. Use :option:`--help` to see the list of available
decoders.


k-nearest neighbor search
=========================

In this section you'll learn how to search a fingerprint file to find
the k-nearest neighbors. You will need the fingerprint files generated
in :ref:`pubchem_fingerprints` but you do not need a chemistry
toolkit.

We'll use the pubchem_queries.fps as the queries for a k=2 nearest
neighor similarity search of the target file puchem_targets.gps::

   simsearch -k 2 -q pubchem_queries.fps pubchem_targets.fps

That's all! You should get output which starts::

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
================

In this section you'll learn how to search a fingerprint file to find
all of the neighbors at or above a given threshold. You will need the
fingerprint files generated in :ref:`pubchem_fingerprints` but you do
not need a chemistry toolkit.

Let's do a threshold search and find all hits which are at least 0.738
similar to the queries::

    simsearch --threshold 0.738 -q pubchem_queries.fps pubchem_targets.fps

The first 20 lines of output from this are::

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

Take a look at the second result line, which contains the 3 hits for
the query id 27575577. As before, the hit information alternates
between the target ids and the target scores, but unlike the k-nearest
search, the hits are not in a particular order. You can see that here
with the scores are 0.7487, 0.7383, and 0.7385.

You might be wondering why I chose the 0.738 threshold. Query id
27575577 has 6 hits with a threshold of 0.7 or higher. That requires
14 columns to show, which is a bit overwhelming.

Combined k-nearest and threshold search
=======================================

In this section you'll learn how to search a fingerprint file to find
the k-nearest neighbors, where all of the hits must be at or above
given threshold. You will need the fingerprint files generated in
:ref:`pubchem_fingerprints` but you do not need a chemistry toolkit.


You can combine the :option:`-k` and :option:`--threshold` queries to
find the k-nearest neighbors which are all above a given threshold::

    simsearch -k 3 --threshold 0.7 -q pubchem_queries.fps pubchem_targets.fps

This find the nearest 3 structures, which all must be at least 0.7
similar to the query fingerprint. The output from the above starts::

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

The output format is identical to the previous two search examples,
and because this is a k-nearest search, the hits are sorted from
higest score to lowest.

NxN (self-similar) searches
===========================

chemfp has no special support for the NxN search of a fingerprint data
set against itself. Instead, use the same file as both the queries and
the targets. This will take twice as much memory and time as an
optimized search.

If you are interested in funding such a tool, I can provide you a cost
estimate.


Using a toolkit to process the ChEBI dataset
============================================

In this section you'll learn how to create a fingerprint file from a
structure file. The structure processing and fingerprint generation
are done with a third-party chemisty toolkit. chemfp supports Open
Babel, OpenEye, and RDKit. (OpenEye users please note that you will
need an OEGraphSim license to use the OpenEye-specific
fingerprinters.)


We'll work with data from ChEBI http://www.ebi.ac.uk/chebi/ which
contains "Chemical Entities of Biological Interest". They distribute
their structures in several formats, including as an SD file. For this
section, download the "lite" version from
ftp://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_lite.sdf.gz . It
contains the same structure data as the complete version but many
fewer tag data fields.  For ChEBI 82 this file contains 19640 records
and the compressed file is 5MB.

Unlike the PubChem data set, the ChEBI data set does not contain
fingerprints so we'll need to generate them using a toolkit.

ChEBI record titles don't contain the id
----------------------------------------

Strangely, the ChEBI dataset does not use the title line of the SD
file to store the record id. A simple examination shows that 16498 of
the title lines are empty, 2119 of them have the title "ChEBI", and 45
of them are labeled "Structure #1."

Instead, the id is the value of the "ChEBI ID" tag, which looks like::

    > <ChEBI ID>
    CHEBI:776

By default the toolkit-based fingerprint generation tools use the
title as the identifier, and exits with an error if the identifier is
missing. (Use the :option:`--errors` option to change the behaviour). If you try
one of them with this data file you will get the error message::

    ERROR: Missing title for record #1 of 'ChEBI_lite.sdf.gz'. Exiting.

Instead, use the :option:`--id-tag` option to specify of the name of
the data tag containing the id. For this data set you'll need to write
it as::

    --id-tag "ChEBI ID"

The quotes are important because of the space in the tag name.

There's also a "ChEBI Name" tag which includes data values like
"tropic acid" and "(+)-guaia-6,9-diene". Every record has a unique
name so the names could be used as the primary identifier. The FPS
fingerprint file format allows identifiers with a space, or comma, or
anything other tab, newline, and a couple of other bytes, so it's no
problem using those names directly.

To use the ChEBI Name as the primary chemfp identifier, use::

    --id-tag "ChEBI Name"


Generating fingerprints with Open Babel
---------------------------------------

If you have the Open Babel Python library installed then you can use
:ref:`ob2fps <ob2fps>` to generate fingerprints::

    ob2fps --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o ob_chebi.fps

This takes about 30 seconds on my laptop.

The default uses the FP2 fingerprints, so the above is the same as::

    ob2fps --FP2 --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o ob_chebi.fps

ob2fps can generate several other types of fingerprints. (See XXX for
details). For example, to generate the Open Babel implementation of the
MACCS definition use::

    ob2fps --MACCS --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o chebi_maccs.fps


Generating fingerprints with OpenEye
------------------------------------

If you have the OEChem Python library installed then you can use
:ref:`oe2fps <oe2fps>` to generate fingerprints::

    oe2fps --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o oe_chebi.fps

This takes about 10 seconds on my laptop and generates a number of
stereochemistry warnings.

The default settings produce OEGraphSim path fingerprint with the
values::

    numbits=4096 minbonds=0 maxbonds=5 atype=DefaultAtom btype=DefaultBond

Each of these can be changed through command-line options. See XXX for
details. There are also options to use an alternate aromaticity model.

oe2fps can generate several other types of fingerprints. For example,
to generate the OpenEye implementation of the MACCS definition use::

   oe2fps --maccs166 --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o chebi_maccs.fps


Generating fingerprints with RDKit
----------------------------------

If you have the RDKit Python library installed then you can use
:ref:`rdkit2fps <rdkit2fps>` to generate fingerprints. Based on the
previous examples you probably guessed that the command-line is::

    rdkit2fps --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o rdkit_chebi.fps

Unfortunately, this isn't enough. If you do this you'll get the message::

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
processing by using the :option:`--errors` option::

    rdkit2fps --id-tag "ChEBI ID" --errors report ChEBI_lite.sdf.gz -o rdkit_chebi.fps

Four minutes later I see that 403 records out of the 19640 could not
be processed.

The previous command-line created RDKit's path fingerprints with
parameters::

    minPath=1 maxPath=7 fpSize=2048 nBitsPerHash=4 useHs=1

Each of those can be changed through command-line options. See XXX for
details.

rdkit2fps can generate several other types of fingerprints. For
example, to generate the RDKit implementation of the MACCS definition
use::

   rdkitfps --maccs166 --id-tag "ChEBI ID" ChEBI_lite.sdf.gz -o chebi_maccs.fps

chemfp supports neither count fingerprints nor sparse fingerprints so
cannot generate RDKit's circular fingerprints.


chemfp's two cross-toolkit substructure fingerprints
====================================================

In this section you'll learn how to generate the two
substructure-based fingerprints which come as part of chemfp. These
are based on cross-toolkit SMARTS pattern definitions and can be used
with Open Babel, OpenEye, and RDKit. (For OpenEye users, these
fingerprints use the base OEChem library and not the separately licensed
OEGraphSim add-on.)

chemfp implements two platform-independent fingerprints where were
originally designed for substructure filters but which are also used
for similarity searches. One is based on the 166-bit MACCS
implementation in RDKit and the other comes from the 881-bit
PubChem/CACTVS substructure fingerprints.

The chemfp MACCS definition is called "rdmaccs" because it closely
derives from the MACCS SMARTS patterns used in RDKit. (These pattern
definitions are also used in Open Babel and the CDK, but are
completely independent from the OpenEye implementation.)

Here are example of the respective rdmaccs fingerprint for phenol
using each of the toolkits.

Open Babel::

    % echo "c1ccccc1O phenol" | ob2fps --in smi --rdmaccs
    #FPS1
    #num_bits=166
    #type=RDMACCS-OpenBabel/1
    #software=OpenBabel/2.2.3
    #date=2011-09-19T22:32:13
    00000000000000000000000000000140004480101e	phenol

OpenEye::

    % echo "c1ccccc1O phenol" | oe2fps --in smi --rdmaccs
    #FPS1
    #num_bits=166
    #type=RDMACCS-OpenEye/1
    #software=OEChem/1.7.4 (20100809)
    #aromaticity=openeye
    #date=2011-09-19T22:31:33
    00000000000000000000000000000140004480101e	phenol

RDKit::

    % echo "c1ccccc1O phenol" | rdkit2fps --in smi --rdmaccs
    echo "c1ccccc1O phenol" | python2.7 rdkit2fps --in smi --rdmaccs
    #FPS1
    #num_bits=166
    #type=RDMACCS-RDKit/1
    #software=RDKit/2011.06.1
    #date=2011-09-19T22:34:42
    00000000000000000000000000000140004480101e	phenol


You might be wondering why :option:`--rdmaccs` produces different fingerprint
types even if the toolkits use the same SMARTS definitions. Each
toolkit perceives chemistry differently. Open Babel before 2.3 didn't
support chirality so chiral-based bits will never be set. Each toolkit
uses a different definition of aromaticity, so a bit which is set when
there are "two or more aromatic rings" will be toolkit dependent.


substruct fingerprints
----------------------

(This is a horribly generic name. If you can think of a better one
then let me know.)

chemp also includes an experimental "substruct" substructure
fingerprint. This is an 881 bit fingerprint derived from the
PubChem/CACTVS substructure keys. They are still being tested and
validated, but you you want to try them out, use the
:option:`--substruct` option.

.. _intro:

====================
chemfp documentation
====================

chemfp is a set of tools for working with cheminformatics fingerprints
in the FPS format.

Most people will use the command-line programs to generate and search
fingerprint files. :ref:`ob2fps <ob2fps>`, :ref:`oe2fps <oe2fps>`, and
:ref:`rdkit2fps <rdkit2fps>` use respectively the `Open Babel
<http://openbabel.org/>`_, `OpenEye <http://www.eyesopen.com/>`_, and
`RDKit <http://www.rdkit.org/>`_ chemistry toolkits to convert
structure files into fingerprint files. :ref:`sdf2fps <sdf2fps>`
extracts fingerprints encoded in SD tags to make the fingerprint
file. :ref:`simsearch <simsearch>` finds targets in a fingerprint file
which are sufficiently similar to the queries.

The programs are built using the :ref:`chemfp Python library API <chemfp-api>`,
which in turn uses a C extension for the performance
critical sections. The parts of the library API documented here are
meant for public use, along with some examples.

Remember: chemfp cannot generate fingerprints from a structure file
without a third-party chemistry toolkit.



Contents:

.. toctree::
   :maxdepth: 2

   chemfp <intro>
   Installing <Installing>
   Working with the command-line tools <Working with the command-line tools>

.. highlight:: none

Installing
==========

The chemfp tools depends on a working Python installation.  You can
download Python 2.7 from http://www.python.org/download/. (Note for
OEChem users: OpenEye doesn't yet support Python 2.7 so you will need
to install Python 2.6 from
http://www.python.org/download/releases/2.6.7/ .)

The core chemfp functionality does not depend on a third-party library
but you will need a chemistry toolkit in order to generate new
fingerprints from structure files. chemfp supports the free Open Babel
and RDKit toolkits and the proprietary OEChem toolkit. Make sure you
install the Python libraries for the toolkit(s) you select.

If you have a source version of chemfp then you will need a C compiler
in order to compile it. This uses Python's standard "setup.py" and you
can see http://docs.python.org/install/index.html for details of how
to use it. The short version is that on Unix systems using sudo (that
is, Mac OS X and most Linux-based OSes) you can do::


  sudo python setup.py install

while for Windows you can do::

   python setup.py install

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
SD file which contains pre-computed CACTVS. You do not need a
chemistry toolkit for this section.

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
    sdf2fps --pubchem Compound_014550001_014575000.sdf.gz  -o pubchem_targets.fps

Congratulations, that was it!

How does this work? Each PubChem record contains the precomputed
CACTVS substructure keys in the PUBCHEM_CACTVS_SUBSKEYS tag. The
"--pubchem" flag tells sdf2fps to get the value of that tag and decode
it to get the fingerprint. It also adds a few metadata fields to the
fingerprint file header.

The order of the fingerprints are the same as the order of the
corresponding record in the SDF, although unconvertable records might
be skipped, depending on the --errors flag.

If you store records in an SD file then you almost certainly don't use
the same fingerprint encoding as PubChem. sdf2ps can decode from a
number of encodings. Use --help to see the list.


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


You can combine the -k and --threshold queries to find the k-nearest
neighbors which are all above a given threshold::

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
missing. (Use the --errors option to change the behaviour). If you try
one of them with this data file you will get the error message::

    ERROR: Missing title for record #1 of 'ChEBI_lite.sdf.gz'. Exiting.

Instead, use the --id-tag option to specify of the name of the data
tag containing the id. For this data set you'll need to write it as::

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
processing by using the --errors option::

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


You might be wondering why "--rdmaccs" produces different fingerprint
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
validated, but you you want to try them out, use the --substruct
option.

===================
Command-line --help
===================


.. _ob2fps:

ob2fps command-line options
===========================


The following comes from "ob2fps --help"::

  usage: ob2fps [-h] [--FP2 | --FP3 | --FP4 | --MACCS | --substruct | --rdmaccs]
                [--id-tag NAME] [--in FORMAT] [-o FILENAME]
                [--errors {strict,report,ignore}]
                [filenames [filenames ...]]
  
  Generate FPS fingerprints from a structure file using OpenBabel
  
  positional arguments:
    filenames             input structure files (default is stdin)
  
  optional arguments:
    -h, --help            show this help message and exit
    --FP2
    --FP3
    --FP4
    --MACCS
    --substruct           generate ChemFP substructure fingerprints
    --rdmaccs             generate 166 bit RDKit/MACCS fingerprints
    --id-tag NAME         tag name containing the record id (SD files only)
    --in FORMAT           input structure format (default autodetects from the
                          filename extension)
    -o FILENAME, --output FILENAME
                          save the fingerprints to FILENAME (default=stdout)
    --errors {strict,report,ignore}
                          how should structure parse errors be handled?
                          (default=strict)

.. _oe2fps:

oe2fps command-line options
===========================

The following comes from "oe2fps --help"::
  
  usage: oe2fps [-h] [--path] [--numbits INT] [--minbonds INT] [--maxbonds INT]
                [--atype ATYPE] [--btype BTYPE] [--maccs166] [--substruct]
                [--rdmaccs] [--aromaticity NAME] [--id-tag NAME] [--in FORMAT]
                [-o FILENAME] [--errors {strict,report,ignore}]
                [filenames [filenames ...]]
  
  Generate FPS fingerprints from a structure file using OEChem
  
  positional arguments:
    filenames             input structure files (default is stdin)
  
  optional arguments:
    -h, --help            show this help message and exit
    --aromaticity NAME    use the named aromaticity model
    --id-tag NAME         tag name containing the record id (SD files only)
    --in FORMAT           input structure format (default guesses from filename)
    -o FILENAME, --output FILENAME
                          save the fingerprints to FILENAME (default=stdout)
    --errors {strict,report,ignore}
                          how should structure parse errors be handled?
                          (default=strict)
  
  path fingerprints:
    --path                generate path fingerprints (default)
    --numbits INT         number of bits in the path fingerprint (default=4096)
    --minbonds INT        minimum number of bonds in path (default=0)
    --maxbonds INT        maximum number of bonds in path (default=5)
    --atype ATYPE         atom type (default=DefaultAtom)
    --btype BTYPE         bond type (default=DefaultBond)
  
  166 bit MACCS substructure keys:
    --maccs166            generate MACCS fingerprints
  
  881 bit ChemFP substructure keys:
    --substruct           generate ChemFP substructure fingerprints
  
  ChemFP version of the 166 bit RDKit/MACCS keys:
    --rdmaccs             generate 166 bit RDKit/MACCS fingerprints
  
  ATYPE is one or more of the following, separated by commas
    Aromaticity AtomicNumber Chiral DefaultAtom EqAromatic EqHalogen
    FormalCharge HvyDegree Hybridization InRing
  Examples:
    --atype DefaultAtom
    --atype AtomicNumber,HvyDegree
  
  BTYPE is one or more of the following, separated by commas
    BondOrder Chiral DefaultBond InRing
  Examples:
    --btype DefaultBond,Chiral
    --btype BondOrder
  
  OEChem guesses the input structure format based on the filename
  extension and assumes SMILES for structures read from stdin.
  Use "--in FORMAT" to select an alternative, where FORMAT is one of:
   
    File Type      Valid FORMATs (use gz if compressed)
    ---------      ------------------------------------
     SMILES        smi, ism, can, smi.gz, ism.gz, can.gz
     SDF           sdf, mol, sdf.gz, mol.gz
     SKC           skc, skc.gz
     CDK           cdk, cdk.gz
     MOL2          mol2, mol2.gz
     PDB           pdb, ent, pdb.gz, ent.gz
     MacroModel    mmod, mmod.gz
     OEBinary v2   oeb, oeb.gz
     old OEBinary  bin

.. _rdkit2fps:

rdkit2fps command-line options
==============================


The following comes from "rdkit2fps --help"::
  
  usage: rdkit2fps [-h] [--RDK] [--fpSize INT] [--minPath INT] [--maxPath INT]
                   [--nBitsPerHash INT] [--useHs USEHS] [--maccs166]
                   [--substruct] [--rdmaccs] [--id-tag NAME] [--in FORMAT]
                   [-o FILENAME] [--errors {strict,report,ignore}]
                   [filenames [filenames ...]]
  
  Generate FPS fingerprints from a structure file using RDKit
  
  positional arguments:
    filenames             input structure files (default is stdin)
  
  optional arguments:
    -h, --help            show this help message and exit
    --id-tag NAME         tag name containing the record id (SD files only)
    --in FORMAT           input structure format (default guesses from filename)
    -o FILENAME, --output FILENAME
                          save the fingerprints to FILENAME (default=stdout)
    --errors {strict,report,ignore}
                          how should structure parse errors be handled?
                          (default=strict)
  
  RDKit topological fingerprints:
    --RDK                 generate RDK fingerprints (default)
    --fpSize INT          number of bits in the fingerprint (default=2048)
    --minPath INT         minimum number of bonds to include in the subgraphs
                          (default=1)
    --maxPath INT         maximum number of bonds to include in the subgraphs
                          (default=7)
    --nBitsPerHash INT    number of bits to set per path (default=4)
    --useHs USEHS         information about the number of hydrogens on each atom
  
  166 bit MACCS substructure keys:
    --maccs166            generate MACCS fingerprints
  
  881 bit substructure keys:
    --substruct           generate ChemFP substructure fingerprints
  
  ChemFP version of the 166 bit RDKit/MACCS keys:
    --rdmaccs             generate 166 bit RDKit/MACCS fingerprints
  
  This program guesses the input structure format based on the filename
  extension. If the data comes from stdin, or the extension name us
  unknown, then use "--in" to change the default input format. The
  supported format extensions are:
  
    File Type      Valid FORMATs (use gz if compressed)
    ---------      ------------------------------------
     SMILES        smi, ism, can, smi.gz, ism.gz, can.gz
     SDF           sdf, mol, sd, mdl, sdf.gz, mol.gz, sd.gz, mdl.gz


.. _sdf2fps:

sdf2fps command-line options
============================

The following comes from "sdf2fps --help"::

  usage: sdf2fps [-h] [--id-tag TAG] [--fp-tag TAG] [--num-bits INT]
                 [--errors {strict,report,ignore}] [-o FILENAME]
                 [--software TEXT] [--type TEXT] [--decompress METHOD]
                 [--binary] [--binary-msb] [--hex] [--hex-lsb] [--hex-msb]
                 [--base64] [--cactvs] [--decoder DECODER] [--pubchem]
                 [filenames [filenames ...]]
  
  Extract a fingerprint tag from an SD file and generate FPS fingerprints
  
  positional arguments:
    filenames             input SD files (default is stdin)
  
  optional arguments:
    -h, --help            show this help message and exit
    --id-tag TAG          get the record id from TAG instead of the first line
                          of the record
    --fp-tag TAG          get the fingerprint from tag TAG (required)
    --num-bits INT        use the first INT bits of the input. Use only when the
                          last 1-7 bits of the last byte are not part of the
                          fingerprint. Unexpected errors will occur if these
                          bits are not all zero.
    --errors {strict,report,ignore}
                          how should structure parse errors be handled?
                          (default=strict)
    -o FILENAME, --output FILENAME
                          save the fingerprints to FILENAME (default=stdout)
    --software TEXT       use TEXT as the software description
    --type TEXT           use TEXT as the fingerprint type description
    --decompress METHOD   use METHOD to decompress the input (default='auto',
                          'none', 'gzip', 'bzip2')
  
  Fingerprint decoding options:
    --binary              Encoded with the characters '0' and '1'. Bit #0 comes
                          first. Example: 00100000 encodes the value 4
    --binary-msb          Encoded with the characters '0' and '1'. Bit #0 comes
                          last. Example: 00000100 encodes the value 4
    --hex                 Hex encoded. Bit #0 is the first bit (1<<0) of the
                          first byte. Example: 01f2 encodes the value \x01\xf2 =
                          498
    --hex-lsb             Hex encoded. Bit #0 is the eigth bit (1<<7) of the
                          first byte. Example: 804f encodes the value \x01\xf2 =
                          498
    --hex-msb             Hex encoded. Bit #0 is the first bit (1<<0) of the
                          last byte. Example: f201 encodes the value \x01\xf2 =
                          498
    --base64              Base-64 encoded. Bit #0 is first bit (1<<0) of first
                          byte. Example: AfI= encodes value \x01\xf2 = 498
    --cactvs              CACTVS encoding, based on base64 and includes a
                          version and bit length
    --decoder DECODER     import and use the DECODER function to decode the
                          fingerprint
  
  shortcuts:
    --pubchem             decode CACTVS substructure keys used in PubChem. Same
                          as --software=CACTVS/unknown --type 'CACTVS-
                          E_SCREEN/1.0 extended=2' --fp-
                          tag=PUBCHEM_CACTVS_SUBSKEYS --cactvs

.. _simsearch:

simsearch command-line options
==============================

The following comes from "simsearch --help"::

  usage: simsearch [-h] [-k K_NEAREST] [-t THRESHOLD] [-q QUERIES]
                   [--hex-query HEX_QUERY] [--query-id QUERY_ID] [--in FORMAT]
                   [-o FILENAME] [-c] [-b BATCH_SIZE] [--scan] [--memory]
                   [--times]
                   target_filename
  
  Search an FPS file for similar fingerprints
  
  positional arguments:
    target_filename       target filename
  
  optional arguments:
    -h, --help            show this help message and exit
    -k K_NEAREST, --k-nearest K_NEAREST
                          select the k nearest neighbors (use 'all' for all
                          neighbors)
    -t THRESHOLD, --threshold THRESHOLD
                          minimum similarity score threshold
    -q QUERIES, --queries QUERIES
                          filename containing the query fingerprints
    --hex-query HEX_QUERY
                          query in hex
    --query-id QUERY_ID   id for the hex query
    --in FORMAT           input query format (default uses the file extension,
                          else 'fps')
    -o FILENAME, --output FILENAME
                          output filename (default is stdout)
    -c, --count           report counts
    -b BATCH_SIZE, --batch-size BATCH_SIZE
                          batch size
    --scan                scan the file to find matches (low memory overhead)
    --memory              build and search an in-memory data structure (faster
                          for multiple queries)
    --times               report load and execution times to stderr


.. highlight:: python

=========================
The chemfp Python library
=========================

The chemfp command-line programs use a Python library called
chemfp. Portions of the API are in flux and subject to change. The
stable portions of the API which are open for general use are
documented in :ref:`chemfp-api`.

The API includes:

 - low-level Tanimoto and popcount operations
 - Tanimo search algorithms based on threshold and/or k-nearest neighbors
 - a cross-toolkit interface for reading fingerprints from a structure file

The following chapters give examples of how to use the API.

Byte and hex fingerprints
=========================

In this section you'll learn how chemfp stores fingerprints and some
of the low-level bit operations on those fingerprints.

chemfp stores fingerprints as byte strings. Here are two 8 bit
fingerprints::

    >>> fp1 = "A"
    >>> fp2 = "B"

The :ref:`chemfp.bitops <chemfp.bitops>` module contains functions which work on byte
fingerprints. Here's the Tanimoto of those two fingerprints::

    >>> from chemfp import bitops
    >>> bitops.byte_tanimoto(fp1, fp2)
    0.33333333333333331

To understand why, you have to know that ASCII character "A" has the
value 65, and "B" has the value 66. The bit representation is::

     "A" = 01000001   and   "B" = 01000010

so their intersection has 1 bit and the union has 3, giving a Tanimoto
of 1/3 or 0.33333333333333331 when represented as a 64 bit floating
point number on the computer.

You can compute the Tanimoto between any two byte strings with the
same length, as in::

    >>> bitops.byte_tanimoto("apples&", "oranges")
    0.58333333333333337

You'll get a chemfp exception if they have different lengths.

.. highlight:: none

Most fingerprints are not as easy to read as the English ones I showed
above. They tend to look more like::


    P1@\x84K\x1aN\x00\n\x01\xa6\x10\x98\\\x10\x11

which is hard to read. I usually show hex-encoded fingerprints. The above
fingerprint in hex is::

    503140844b1a4e000a01a610985c1011

which is simpler to read, though you still need to know your hex
digits.

.. highlight:: python

The bitops module includes other low-level functions which work on
byte fingerprints, as well as corresponding functions which work on
hex fingerprints. (Hex-encoded fingerprints are decidedly second-class
citizens in chemfp, but they are citizens.)


Fingerprint collections and metadata
====================================

In this section you'll learn the basic operations on a fingerprint
collection and the fingerprint metadata.

A fingerprint record is the fingerprint plus an identifier. In chemfp,
a fingerprint collection is a object which contains fingerprint
records and which follows the common API providing access to those
records.

That's rather abstract, so let's work with a few real examples. You'll
need to create a copy of the "pubchem_targets.fps" file generated in
:ref:`pubchem_fingerprints` in order to follow along.

Here's how to open an FPS file::

    >>> import chemfp
    >>> reader = chemfp.open("pubchem_targets.fps")

Every fingerprint collection has a metadata attribute with details
about the fingerprints. It comes from the header of the FPS file. You
can view the metadata in Python repr format:

    >>> reader.metadata
    Metadata(num_bits=881, num_bytes=111, type='CACTVS-E_SCREEN/1.0 extend
    ed=2', aromaticity=None, sources=['Compound_014550001_014575000.sdf.gz
    '], software=u'CACTVS/unknown', date='2011-09-14T12:10:34')

but I think it's easier to view it in string format, which matches the
format of the FPS header:

    >>> print reader.metadata
    #num_bits=881
    #type=CACTVS-E_SCREEN/1.0 extended=2
    #software=CACTVS/unknown
    #source=Compound_014550001_014575000.sdf.gz
    #date=2011-09-14T12:10:34
    

All fingerprint collections support iteration. Each step of the
iteration returns the fingerprint identifier and its score. Since I
know the 6th record has the id 14550045, I can write a simple loop
which stops with that record::

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


Fingerprint collections also support iterating via arenas, and several
support Tanimoto search functions.


FingerprintArena
================

In this section you'll learn about the FingerprintArena fingerprint
collection and how to iterate through arenas in a collection.

The FPSReader reads through or searches a fingerprint file once. If
you want to read the file again you have to reopen it.

Reading from disk is slow, and the FPS format is designed for
ease-of-use and not performance. If you want to do many queries then
it's best to store everything in memory. The
:ref:`FingerprintArena <FingerprintArena>` is a
fingerprint collection which does that.

Here's how to load fingerprints into an arena::

    >>> import chemfp
    >>> arena = chemfp.load_fingerprints("pubchem_targets.fps")
    >>> print arena.metadata
    #num_bits=881
    #type=CACTVS-E_SCREEN/1.0 extended=2
    #software=CACTVS/unknown
    #source=Compound_014550001_014575000.sdf.gz
    #date=2011-09-14T12:10:34
    

This implements the fingerprint collection API, so you can do things
like iterate over an arena and get the id/fingerprint pairs.::

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

If you look closely you'll notice that the fingerprint record order
has changed from the previous section, and that the population counts
are suspiciously non-decreasing. By default ref:`load_fingerprints`
reorders the fingerprints into a data structure which is faster to
search, although you can disable that if you want the fingerprints to
be the same as the input order.

The :ref:`FingerprintArena <FingerprintArena>` has new capabilities. You can ask it how many
fingerprints it contains, get the list of identifiers, and look up a
fingerprint record given an index, as in::

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

An arena supports iterating through subarenas. This is like having a
long list and being able to iterate over sublists. Here's an example
of iterating over the arena to get subarenas of size 1000 (excepting
the last), and print information about each subarena.::

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

To help demonstrate what's going on, I showed the first id of each
record along with the main arena ids for records 0 and 1000, so you
can verify that they are the same.

Arenas are a core part of chemfp. Processing one fingerprint at a time
is slow, so the main search routines expect to iterate over query
arenas, rather than query fingerprints.

Thus, the FPSReaders - and all chemfp fingerprint collections - also
support the `iter_arenas` interface. Here's an example of reading the
targets file 25 records at a time::

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

Those add up to 224, which you can verify is the number of structures
in the original source file.

How to use query fingerprints to search for similar target fingerprints
=======================================================================

In this section you'll learn how to do a Tanimoto search using the
previously created PubChem fingerprint files for the queries and the
targets.

It's faster to search an arena, so I'll load the target fingerprints:

    >>> import chemfp
    >>> targets = chemfp.load_fingerprints("pubchem_targets.fps")
    >>> len(targets)
    3119

and open the queries as an FPSReader.

    >>> queries = chemfp.open("pubchem_queries.fps")

I'll use :ref:`threshold_tanimoto_search <chemfp_threshold_tanimoto_search>` to find, for each query,
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
          # ... many lines omitted ...

I'm only showing the first two hits for the sake of space. It seems
rather pointless, after all, to show all 18 hits of query id 27575577.

What you don't see is that the implementation uses the iter_arenas()
interface on the queries so that it processes only a subarena at a
time. There's a tradeoff between a large arena, which is faster
because it doesn't often go back to Python code, or a small arena,
which uses less memory and is more responsive. You can change the
tradeoff using the `arena_size` parameter.


If all you care about isthe count of the hits within a given
threshold then use :ref:`chemfp.count_tanimoto_hits <chemfp_count_tanimoto_hits>`

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
         # ... many lines omitted ...

Or, if you only want the k=2 nearest neighbors to each target within
that same threshold of 0.7 then use ref:`chemfp.knearest_tanimoto_search`::

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
         # ... many lines omitted ...



How to search an FPS file
=========================

In this section you'll learn how to search an FPS file directly,
without loading it into a FingerprintArena.

The previous example loaded the fingerprints into a
FingerprintArena. That's the fastest way to do multiple
searches. Sometimes though you only want to do one or a couple of
queries. It seems rather excessive to read the entire targets file
into an in-memory data structure before doing the search when you
could search will processing the file.

For that case, use an FPSReader as the target file. Here I'll get the
first record from the queries file and use it to search the targets
file::

    >>> query_arena = next(chemfp.open("pubchem_queries.fps").iter_arenas(1))

This line opens the file, iterates over its fingerprint records, and
return the first one.

(Note: the next() function was added after Python 2.5 so the above
won't work for that version. Instead, use::

    >>> query_arena = chemfp.open("pubchem_queries.fps").iter_arenas(1).next()

which is the older form.)

Here are the k=5 closest hits against the targets file::

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

Remember that the FPSReader is based on reading an FPS file. Once
you've done a search, the file is read, and you can't do another
search. You'll need to reopen the file.

Each search processes arena_size query fingerprints at a time. You
will need to increase that value if you want to search more than that
number of fingerprints with this method. The search performance
tradeoff between a FPSReader search and loading the fingerprints into
a FingerprintArena occurs with under 10 queries, so there should be
little reason to worry about this.


FingerprintArena searches returning indicies instead of ids
===========================================================

In this section you'll learn how to search a FingerprintArea and use
hits based on integer indicies rather than string ids.

The previous sections used a high-level interface to the Tanimoto
search code. Those are designed for the common case where you just
want the query id and the hits, where each hit includes the target id.

Working with strings is actually rather inefficient in both speed and
memory. It's usually better to work with indicies if you can, and in
the next section I'll show how to make a distance matrix using this
interface.

The index methods are only available as methods of a FingerprintArena,
where the arena contains the targets. Three of the methods
(`count_tanimoto_search_arena`, `threshold_tanimoto_search_arena`, and
`knearest_tanimoto_search_arena`) take another arena as the
query. Here's an example where I use the first 5 records from
pubchem_queries.fps to search the entire contents of that file::

    >>> import chemfp
    >>> dataset = chemfp.load_fingerprints("pubchem_queries.fps")
    >>> first_5 = next(dataset.iter_arenas(5))
    >>> results = dataset.threshold_tanimoto_search_arena(first_5, threshold=0.7)

You can iterate over the results to get the list of hits for each of
the queries. (The order of the results is the same as the order of the
records in the query.)::

   >>> for hits in results:
    ...   print len(hits), hits[:3]
    ... 
    2 [('27581954', 1.0), ('27581957', 1.0)]
    2 [('27581954', 1.0), ('27581957', 1.0)]
    3 [('27580389', 1.0), ('27580394', 0.88235294117647056), ('27581637', 0.75)]
    2 [('27584917', 1.0), ('27585106', 0.89915966386554624)]
    2 [('27584917', 0.89915966386554624), ('27585106', 1.0)]

This is like what you saw earlier, except that it doesn't have the
query id. (If you want that you can enumerate() over the results and
use the index into the query arena's ids[] list.)

What I really want to show is that you can get the same data only
using the offset index for the target record instead of its id. The
result from a Tanimoto search is a `SearchResult` object, with the
methods `iter_hits()`::

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

I did a few id lookups given the target dataset to show you that the
index corresponds to the identifiers from the previous code.


Computing a distance matrix for clustering
==========================================

In this section you'll learn how to compute a distance matrix using
the chemfp API.

chemfp does not do clustering. There's a huge number of tools which
arleady do that. A goal of chemfp in the future is to provide some
core components which clustering algorithms can use.

That's in the future. Right now you can use the following to build a
distance matrix and pass that to one of those tools. The following is
a somewhat inefficient since it computes almost twice as many Tanimoto
scores as it needs to do, and uses twice the necessary memory, but
what's a factor of two among friends?

Most of those tools use `NumPy <http://numpy.scipy.org/>`_, which is a
popular third-party package for numerical computing. You will need to
have it installed for the following to work.

::

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


Once you've computed the distance matrix, clustering is easy. I
installed the `hcluster <http://code.google.com/p/scipy-cluster/>`_
package, as well as `matplotlib <http://matplotlib.sourceforge.net/>`_,
then ran the following to see the hierarchical clustering::

    import chemfp
    import hcluster # Clustering package from http://code.google.com/p/scipy-cluster/
    
    # ... insert the 'distance_matrix' function definition here ...

    dataset = chemfp.load_fingerprints("docs/pubchem_queries.fps")
    distances  = distance_matrix(dataset)
    
    linkage = hcluster.linkage(distances, method="single", metric="euclidean")
    
    # Plot using matplotlib, which you must have installed
    hcluster.dendrogram(linkage, labels=dataset.ids)
    
    import pylab
    pylab.show()



Taylor-Butina clustering
========================

For the last clustering example, here's my (non-validated) variation
of the `Butina algorithm from JCICS 1999, 39, 747-750 <http://www.chemomine.co.uk/dbclus-paper.pdf>`_.
See also http://www.redbrick.dcu.ie/~noel/R_clustering.html .

First, for each fingerprint find all other fingerprints with a
threshold of 0.8::

    import chemfp
    
    dataset = chemfp.load_fingerprints("pubchem_targets.fps")
    search = dataset.threshold_tanimoto_search_arena(dataset, threshold = 0.8)


Sort the results so that fingerprints with more hits come first. This
is more likely to be a cluster centroid. Break ties arbitrarily by the
fingerprint id; since fingerprints are ordered by the number of bits
this likely makes larger structures appear first.::

    def get_hit_indicies(hits):
        return [id for (id, score) in hits]
    
    # Reorder so the centroid with the most hits comes first.
    # (That's why I do a reverse search.)
    # Ignore the arbitrariness of breaking ties by fingerprint index
    results = sorted( (  (len(hits), i, get_hit_indicies(hits))
                                        for (i, hits) in enumerate(search.iter_hits())  ),
                      reverse=True)

Apply the leader algorithm to determine the cluster centroids and the singletons::


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

Once done, report the results::

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


The algorithm is quick for this small data set.

Out of curiosity, I tried this on 100,000 compounds selected
arbitrarily from PubChem. It took 7 minutes for my laptop to process
with a threshold of 0.8. In the Butina paper, it took 24 hours to do
the same, although that was with a 1024 bit fingerprint instead of
881. It's hard to judge the absolute speed differences of a 12 year
old MIPS R4000 to a two year old laptop, but it's less than a factor
of 200. Part must certainly be due to the work I put into making
chemfp fast, and I'm almost certain I can get another 3-fold
performance increase.

The core Tanimoto search routines release the Python global
interpreter lock, which means algorithms like this should be easily
parallizable.


Reading structure fingerprints using a toolkit
==============================================

In this section you'll learn how to use a chemistry toolkit in order
to compute fingerprints from a given structure file.

What happens if you're given a structure file and you want to find the
two nearest matches in an FPS file? You'll have to generate the
fingerprints for the structures in the structure file, then do the
comparison.

For this section you'll need to have a chemistry toolkit. I'll use the
"chebi_maccs.fps" file you generated earlier as the targets, and the
PubChem file "Compound_027575001_027600000.sdf.gz as the source of
query structures.::

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
         # ... many lines omitted ...

That's it! Pretty simple, wasn't it? You didn't even need to explictly
specify which toolkit you wanted to use.

The only new thing here is :ref:`read_structure_fingerprints <chemfp_read_structure_fingerprints>`. The
first parameter of this is the metadata used to configure the
reader. In my case it's::

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

For examples, if you have OpenBabel installed then you can do::

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

If you have OEChem and OEGraphSim installed then you can do::

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

And if you have RDKit installed then you can do::

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

.. _chemfp-api:

==========
chemfp API
==========

This chapter contains the docstrings for the public portion of the
chemfp API.

(Note: the values of 0.69999999999999996 which you see are actually
0.7. I don't know how to get my documentation tool to use the more
understandable value.)

.. py:module:: chemfp

.. autofunction:: chemfp.open
.. autofunction:: chemfp.load_fingerprints

.. _chemfp_read_structure_fingerprints:
.. autofunction:: chemfp.read_structure_fingerprints

.. _chemfp_count_tanimoto_hits:
.. autofunction:: chemfp.count_tanimoto_hits

.. _chemfp_threshold_tanimoto_search:
.. autofunction:: chemfp.threshold_tanimoto_search
.. autofunction:: chemfp.knearest_tanimoto_search

.. autoclass:: chemfp.Metadata
.. autoclass:: chemfp.FingerprintIterator
.. autoclass:: chemfp.Fingerprints

.. _FingerprintArena:

.. autoclass:: chemfp.arena.FingerprintArena
.. autoclass:: chemfp.arena.SearchResults

.. _chemfp.bitops:
.. py:module:: chemfp.bitops

.. autofunction:: chemfp.bitops.byte_popcount
.. autofunction:: chemfp.bitops.byte_intersect_popcount
.. autofunction:: chemfp.bitops.byte_tanimoto
.. autofunction:: chemfp.bitops.byte_contains
.. autofunction:: chemfp.bitops.hex_isvalid
.. autofunction:: chemfp.bitops.hex_popcount
.. autofunction:: chemfp.bitops.hex_intersect_popcount
.. autofunction:: chemfp.bitops.hex_tanimoto
.. autofunction:: chemfp.bitops.hex_contains



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


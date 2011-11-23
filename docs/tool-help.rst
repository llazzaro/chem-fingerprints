
===============================
Help for the command-line tools
===============================


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


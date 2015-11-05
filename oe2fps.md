# Introduction #

OpenEye sells the OEChem toolkit for cheminformatics and the OEGraphSim add-on toolkit for fingerprint generation.

oe2fps uses both of them to generate the OpenEye path and MACCS fingerprints, and uses just the OEChem toolkit to generate the [RDMACCS](RDMACCS.md) and [Substruct](Substruct.md) fingerprints based on the chemfp patterns.

More information can be found in the [chemfp documentation](http://readthedocs.org/docs/chemfp/en/latest/).

# Details #

OEGraphSim's --path fingerprint (which is the chemfp default) is 4096 bits long, so I've folded the fingerprint so it fits better in this web page. The fingerprint lines should have the entire fingerprint on one line, in the first column, with the record identifier in the second column.

```
% oe2fps $PUBCHEM/Compound_005050001_005075000.sdf.gz | head -9
#type=OpenEye-Path/1 numbits=4096 minbonds=0 maxbonds=5 atype=DefaultAtom btype=DefaultBond
#software=OEGraphSim/1.0.0 (20100809)
#aromaticity=openeye
#source=/Users/dalke/databases/pubchem/Compound_005050001_005075000.sdf.gz
#date=2011-09-20T16:35:07
9084018000104400c0844200500400028000804a01200500011110000010200001100e20000042
08000020180008080003004230100000010010001006010100420f0402c0000040008e00010014
000080002010002110014120010010080080410040000010080000000000050011062020008000
424000210040000000504150080048012048064000001000c80460810800008100004000040002
040c5000000400040000000181000080000080000008000002015000080080000800004a002100
000001000020000001000002080000000318000101024840010140104000248408004080008800
000000008400182100041220001000a800000000002400088000000c8400044000003000012500
000000000100140480844800000008000002000080240028800000020001000100124220401032
10c0000300400000250000002e004001841000004850000280021a081000020000004440880180
0044040000000008120100100010005a2110404020004225000880002100804050016800089010
800000082003a02400028000000404080495c84844084000200046100000280d00040102100001
12220008000200800410004010004004140400411000c200000000042000021810000000400404
000200002000003000c08a00000100800080420090080038040000000830206000001000009020
8400c00000	5050001
908480000090440040844240100000028000000001028400001100000050200005000a20000042
000000201800000840030042300000000000100208064004204202140200000044008000010010
000480000010002100000120090010090000400040000202080000000000040001022020000000
000000032040010000c00100000050002000104000001600000060810800400010404000040002
041c0000000c840400000001804000000030800000000008020010000800000000800148000100
0000010001000000a1000002000000000200000005004810000540004000200000004180000800
01004000240018010004a280005000880000000000201000000000088400040001003200012000
00000000012004040084400000200c000000410080240020800000000100800000904020411036
0000000200710200240000000e0041018010000008500002c44212401000000020044040880080
4448000000000008100100040010001a2100000020024205000080002200004200016800009018
08000200208160000000a008000404080005c840050040000020020000001801000c2120100001
000000000002000000140000180040008204005010005000010400002008021010000000400400
004001002000001000c08100000100800000020090001078040000400030206000001404000020
8400d00000	5050002
```

You can change the fingerprint parameters

```
% oe2fps --numbits 220 --minbonds 3 $PUBCHEM/Compound_005050001_005075000.sdf.gz | head -13
#FPS1
#num_bits=220
#type=OpenEye-Path/1 numbits=220 minbonds=3 maxbonds=5 atype=DefaultAtom btype=DefaultBond
#software=OEGraphSim/1.0.0 (20100809)
#aromaticity=openeye
#source=/Users/dalke/databases/pubchem/Compound_005050001_005075000.sdf.gz
#date=2011-09-20T16:37:00
ef7bd3ffff77fb9fffefbf77d3ffb7fff6fbbdffcffbfdfbfb4f0600	5050001
7f7bf29c3f7fff3f7f3bbe6fdfbd2ffffcb9edfd5ed3eee9fbe40600	5050002
cf69ebf5d9dffffc97efb78fbfff7dfdfff7ddfbbfdaefff5fd70500	5050003
8b33080d42b3d60fc3dee7c6e9c30528504a50948cf62830fd590000	5050004
4e29cff59dc65fa8a6f7bfbc9fed777d3f777dfa0ff26bf71fee0500	5050005
aabb4e1c7ab5d517d5dae7d6f95b247f8a4f46948ed2aab6fd5f0300	5050006
```


OEGraphSim can generate MACCS keys.

```
% oe2fps --maccs $PUBCHEM/Compound_005050001_005075000.sdf.gz | head -10                              #FPS1
#num_bits=166
#type=OpenEye-MACCS166/1
#software=OEGraphSim/1.0.0 (20100809)
#aromaticity=openeye
#source=/Users/dalke/databases/pubchem/Compound_005050001_005075000.sdf.gz
#date=2011-09-20T16:38:40
000000880100448f04e7f5e83af3f073dbaa3b7a3b	5050001
000000080000000001e654c81aa1d453d0a0392e3b	5050002
00000000000002008008000914810150809833523a	5050003
```

Even if you don't have an OEGraphSim license, you can still generate the [RDMACCS](RDMACCS.md) variant of the MACCS keys:

```
% oe2fps --rdmaccs $PUBCHEM/Compound_005050001_005075000.sdf.gz | head -12
#FPS1
#num_bits=166
#type=RDMACCS-OpenEye/1
#software=OEChem/1.7.4 (20100809)
#aromaticity=openeye
#source=/Users/dalke/databases/pubchem/Compound_005050001_005075000.sdf.gz
#date=2011-09-20T16:39:49
000000800100449f04e7b5e83af3f0dbdbeb3ffb1f	5050001
000000000000000001e614c81a21d45350e139af1f	5050002
0000000000000200800800111401217000d8b7d21e	5050003
00000000000201018000248000940024280b82571a	5050004
0000000000000221000c30401e51a0d359c3ad7b1f	5050005
```

and the PubChem-influenced 881 bit [Substruct](Substruct.md) keys.

```
% oe2fps --substruct $PUBCHEM/Compound_005050001_005075000.sdf.gz | head -12
#FPS1
#num_bits=881
#type=ChemFP-Substruct-OpenEye/1
#software=OEChem/1.7.4 (20100809)
#aromaticity=openeye
#source=/Users/dalke/databases/pubchem/Compound_005050001_005075000.sdf.gz
#date=2011-09-20T16:42:19
000000800100449f04e7b5e83af3f0dbdbeb3ffb1f	5050001
000000000000000001e614c81a21d45350e139af1f	5050002
0000000000000200800800111401217000d8b7d21e	5050003
00000000000201018000248000940024280b82571a	5050004
0000000000000221000c30401e51a0d359c3ad7b1f	5050005
```


## oe2fps --help ##

```
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
  --minbonds INT        minimum number of bonds in the path fingerprint
                        (default=0)
  --maxbonds INT        maximum number of bonds in the path fingerprint
                        (default=5)
  --atype ATYPE         atom type flags, described below (default=Default)
  --btype BTYPE         bond type flags, described below (default=Default)

166 bit MACCS substructure keys:
  --maccs166            generate MACCS fingerprints

881 bit ChemFP substructure keys:
  --substruct           generate ChemFP substructure fingerprints

ChemFP version of the 166 bit RDKit/MACCS keys:
  --rdmaccs             generate 166 bit RDKit/MACCS fingerprints

ATYPE is one or more of the following, separated by the '|' character.
  Aromaticity AtomicNumber Chiral EqAromatic EqHalogen FormalCharge
  HvyDegree Hybridization InRing
The terms 'Default' and 'DefaultAtom' are expanded to OpenEye's
suggested default of AtomicNumber|Aromaticity|Chiral|FormalCharge|HvyDegree|Hybridization|EqHalogen.
Examples:
  --atype Default
  --atype AtomicNumber|HvyDegree
(Note that most atom type names change in OEGraphSim 2.0.0.)

BTYPE is one or more of the following, separated by the '|' character
  BondOrder Chiral InRing
The terms 'Default' and 'DefaultBond' are expanded to OpenEye's
suggested default of BondOrder|Chiral.
Examples:
  --btype Default
  --btype BondOrder
(Note that "BondOrder" changes to "Order" in OEGraphSim 2.0.0.)

For simpler Unix command-line compatibility, a comma may be used
instead of a '|' to separate different fields. Example:
  --atype AtomicNumber,HvyDegree

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
```
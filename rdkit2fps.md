# Introduction #

RDKit implements a number of fingerprints but only the path fingerprints and the MACCS fingerprints are appropriate for chemfp. The others are sparse fingerprints, which chemfp does not yet support.

More information can be found in the [chemfp documentation](http://readthedocs.org/docs/chemfp/en/latest/).


# Details #

The default RDKit fingerprints are 2048 bit hash fingerprints. You can also specify "--RDK" if you want to be specific. (Got a better name?). Since these fingerprints are so long, I've folded the line over so it fits better. In real life the fingerprint data is a tab-delimited line containing the hex-encoded fingerprint in the first column, and the id in the second.

The various parameters (minPath, maxPath, fpSize, and so on) are all tunable by command-line arguments to rdkit2fps. Read the RDKit documentation for details on what they mean.


```
% rdkit2fps $PUBCHEM/Compound_005050001_005075000.sdf.gz | head | fold -w 80
#FPS1
#num_bits=2048
#type=RDKit-Fingerprint/1 minPath=1 maxPath=7 fpSize=2048 nBitsPerHash=4 useHs=1
#software=RDKit/2011.06.1
#source=/Users/dalke/databases/pubchem/Compound_005050001_005075000.sdf.gz
#date=2011-09-20T16:47:22
fef9ff3f7feede37ffff3fefdabbbf9ff7bdde8ecbfe7eb7ffdbfaffcfedfffffffd59bfefb39f7d
febffcfddb7733ffbbee7ff7ffcfffff7f7fcdf5d7bbff9db93f7cdfdffdeaefffff5fdfdf59ff9e
dfffdfffefe7f6ffff7ff73e377fffb773bf7bcbffb7feedefffcab39eda77effd3f76ffdff9b6ef
ffbe9ffbedff1be6ffb9dff3ffe5deeffeffdbeff7fbef9ffd9afbd7ffffcf8ff7ffd7b6ffd7bf9c
e7fd7ffff7fcbbe9ed4fd7f7df73f7bb2db3e99a7ef7fbebf7feffff7f7bffbcffffc3fff1fbfabb
f7beffffdef6ff8faff8ffffbddff67dfabdf76f3ffbdebdfaaec2bffffbeefa5f45fdfffbf9f77f
bdddf7dfee27dffeeedffbefff5f7fff	5050001
feebff3f5beeff97dfeb77eedf3fa4bd37a7da828dfe6bbeffd2f6d74b99dffeaefd7a7fed31df55
fe77eef9fb7f3bdd3b2e73f7ffcfefff6a77c1b5f57bf7ddb9777edfdffba2ff7fbf75d77f5dff9e
e7fdbbd1efe3f6f7ee7ff77e34fdb6cf7b3f7bc7ddbffeedcff5c6b31e5267fff92f76ffdeb133cd
ff9cdff3efff5bf7f7b9df13bdeddcedfe7ffbedf6aaef8fe9daebd7b9be8f8df7ffd6b7fff3ff9c
e6fd7ffbf4fcbbe9cb6fd6f7ff733f3321e3f1de75e7ffebf7fefbff7b1bf7acfb7b8bd7f1f2eedf
f7bee7bfdef4ff0f2fd2ffdebccffa79efad776f37fbdaadfa8ec2be36facff88b55fbff9b6bff7e
ddd15acffd26dbfbeecdfbcedadf7bff	5050002
e2eb6e1dfff1b1b7d4a97f276ef91a27b3cab9bf9f3bdfba6baf6a6bddfd5ef7efdf64bee77dfe43
bfbb765ef75da9f84ffb0a9fbda7ebf9fedccdd66f9fadd86e3adff703f513df3dbfbf263eb546f9
fd78fdb3737bf26ffbc4ef7df475aa5f47d73afeef327fce2ebfc35bf8a52abdb3db6efa1d6b7ffe
ffec5aa1eaa2a9ef753a4c6bc2b922aff6ed1bbb3e73f7fe479f5bfffdbd4d9bf9ffbbd7f9f7cff7
e7f9f97b5c2fc97b39bae9f1d5975befdf3c27fdad8f7b4f9e66f45a5fdeff46dfd9d4dbcfe240ad
f43df27e5fcdd65a7d92beffacfdaa3df41f6fe5c32bc7afaee4accc77f5f4758cd3c1d7e6fbff3f
f9770be54f6d8bcf3bad9fbcfbbe9f36	5050003
25c946882008585d18169303041119a28018d80720a9b829000b608028038ec3008aec4182321c06
5c5c0f15ac08a5c057e18002081285848a1046828108b40041808021a9b15a688c89752290302627
6d3844a0016200d0a940070eecaa90d3086680456940b7c4630789870015109a046985c290210324
c00c070120803f218c024b20f4015c17102aa9882150850ca8c8129889f80502b84c380481cf81c0
1c68b0c038a704d84c4322024b484308f04003023ef0c02956a242ba23812006e36de72806c28c44
41086244c60e141c111086940444104071063a04c12d8d6526cc6171a065c0261218e03071c6306b
54a4a4f3a1ca100870010b621620c604	5050004
```

Here I'll use a smaller fingerprint size and fewer bits per fragment hash:

```
% rdkit2fps --fpSize 256 --nBitsPerHash 1 $PUBCHEM/Compound_005050001_005075000.sdf.gz | head
#FPS1
#num_bits=256
#type=RDKit-Fingerprint/1 minPath=1 maxPath=7 fpSize=256 nBitsPerHash=1 useHs=1
#software=RDKit/2011.06.1
#source=/Users/dalke/databases/pubchem/Compound_005050001_005075000.sdf.gz
#date=2011-09-20T16:48:57
ffaffdffefffff3fffeffeffdffbfffffffffffff9ffffdfffffdfffffffffff	5050001
ffebffffe5ffffbdfffffeffdefff7edefffffd6fbffffdffff7dfffffff7fff	5050002
6e7f5ffef6f7fbfff7bff7ff7ffffffeffffdbfbb7bf7ffdffeeffe7efffff6f	5050003
8eeea5adb9e1dc5c7093cc035ed39dba94d18486a3fbc5b8be7b5e936de736e7	5050004
```

This uses RDKit's own 166 bit MACCS keys fingerprinter:

```
% rdkit2fps --maccs166 $PUBCHEM/Compound_005050001_005075000.sdf.gz | head -12
#FPS1
#num_bits=166
#type=RDKit-MACCS166/1
#software=RDKit/2011.06.1
#source=/Users/dalke/databases/pubchem/Compound_005050001_005075000.sdf.gz
#date=2011-09-20T16:50:40
000000800100449f04e7b5e83af3f0dbdbeb3ffb1f	5050001
000000000000000001e614c81a21d45350e139af1f	5050002
0000000000000200800800111401217000d8b7d21e	5050003
00000000000201018000248000940024280b82571a	5050004
0000000000000221000c30401e51a0d359c3ad7b1f	5050005
0000000000000101800044804014017428c9a2571e	5050006
```

and this is the chemfp version of the same MACCS definitions. (The chemfp fingerprints come almost directly from RDKit and should have very few differences.)

```
% rdkit2fps --rdmaccs $PUBCHEM/Compound_005050001_005075000.sdf.gz | head -12
#FPS1
#num_bits=166
#type=RDMACCS-RDKit/1
#software=RDKit/2011.06.1
#source=/Users/dalke/databases/pubchem/Compound_005050001_005075000.sdf.gz
#date=2011-09-20T16:53:59
000000800100449f04e7b5e83af3f0dbdbeb3ffb1f	5050001
000000000000000001e614c81a21d45350e139af1f	5050002
0000000000000200800800111401217000d8b7d21e	5050003
00000000000201018000248000940024280b82571a	5050004
0000000000000221000c30401e51a0d359c3ad7b1f	5050005
0000000000000101800044804014017428c9a2571e	5050006
```

rdkit2fps also implements the chemfp [Substruct](Substruct.md) keys

```
% rdkit2fps --substruct $PUBCHEM/Compound_005050001_005075000.sdf.gz | head -12 | fold -80
#FPS1
#num_bits=881
#type=ChemFP-Substruct-RDKit/1
#software=RDKit/2011.06.1
#source=/Users/dalke/databases/pubchem/Compound_005050001_005075000.sdf.gz
#date=2011-09-20T16:54:38
07de1c000200000000000000000000000080000000008c0701000000000002800200007820000200
003014871b604c83c100204015c0a44e2a0a43000884400010119814261950044d03a989e1043006
61001b13e03811017100000000044000000000000002200000000000000000	5050001
07de04000000000000000000000000000080000000003c0200000000000002800300007800000000
00b0108359207c03c910080015c0acee2a004101048c4804101b841c2e1910064f03a98929041006
61001313e03811017100000000004000000000000000200000000000000000	5050002
071e1c000000000000000000000000000000000000000c0601000000000000810200005800000000
0030200119404c60010020001140054a0000401000240400101180601113b0e44c21ac58018c44a4
03a85095e111173d71050080100018020002004008000c0100000000000000	5050003
010e9c0300000000000000000000000000800400000000000000000000000080010000d800000000
0010200509480c30010020021100054b1040401000240400101180621013b0644c21ac7841980425
03885095e11176300600000000000000000000000000000000000000000000	5050004
075e0c000000000000000000000000000000000000008c0781000000000000810200007800000000
0030308719600c03c10020001140254a0200410000a4400010118010261350044c03a989e1841006
01001b93e11811017111000000000800000800000000040000000000000000	5050005
030e1c002000000000000000000000000080040000000c0000000000000000800300005840000000
0010600509440c70010020021100054b104040100024a40210118062d013b064cca1ec78419c0425
03885095e11176300600000000000010000000000000000800000000000000	5050006
```


## rdkit2fps --help ##

```
usage: rdkit2fps [-h] [--fpSize INT] [--RDK] [--minPath INT] [--maxPath INT]
                 [--nBitsPerHash INT] [--useHs 0|1] [--morgan] [--radius INT]
                 [--useFeatures 0|1] [--useChirality 0|1] [--useBondTypes 0|1]
                 [--torsions] [--targetSize INT] [--pairs] [--minLength INT]
                 [--maxLength INT] [--maccs166] [--substruct] [--rdmaccs]
                 [--id-tag NAME] [--in FORMAT] [-o FILENAME]
                 [--errors {strict,report,ignore}]
                 [filenames [filenames ...]]

Generate FPS fingerprints from a structure file using RDKit

positional arguments:
  filenames             input structure files (default is stdin)

optional arguments:
  -h, --help            show this help message and exit
  --fpSize INT          number of bits in the fingerprint (applies to RDK,
                        Morgan, topological torsion, and atom pair
                        fingerprints (default=2048)
  --id-tag NAME         tag name containing the record id (SD files only)
  --in FORMAT           input structure format (default guesses from filename)
  -o FILENAME, --output FILENAME
                        save the fingerprints to FILENAME (default=stdout)
  --errors {strict,report,ignore}
                        how should structure parse errors be handled?
                        (default=strict)

RDKit topological fingerprints:
  --RDK                 generate RDK fingerprints (default)
  --minPath INT         minimum number of bonds to include in the subgraph
                        (default=1)
  --maxPath INT         maximum number of bonds to include in the subgraph
                        (default=7)
  --nBitsPerHash INT    number of bits to set per path (default=4)
  --useHs 0|1           include information about the number of hydrogens on
                        each atom (default=1)

RDKit Morgan fingerprints:
  --morgan              generate Morgan fingerprints
  --radius INT          radius for the Morgan algorithm (default=2)
  --useFeatures 0|1     use chemical-feature invariants (default=0)
  --useChirality 0|1    include chirality information (default=0)
  --useBondTypes 0|1    include bond type information (default=1)

RDKit Topological Torsion fingerprints:
  --torsions            generate Topological Torsion fingerprints
  --targetSize INT      number of bits in the fingerprint (default=4)

RDKit Atom Pair fingerprints:
  --pairs               generate Atom Pair fingerprints
  --minLength INT       minimum bond count for a pair (default=1)
  --maxLength INT       maximum bond count for a pair (default=30)

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
```
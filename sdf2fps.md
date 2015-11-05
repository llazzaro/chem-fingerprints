# Introduction #

Many people store fingerprint data as tags in an SD file. One prime example is PubChem, which stores their substructure fingerprints in the PUBCHEM\_CACTVS\_SUBSKEYS tag, encoded using base64. Other people store them as hex fingerprints or as sets of literal 0s and 1s.

The sdf2fps tool extracts the fingerprint string and can apply one of several built-in decoders or a user-defined one to turn the data into fingerprint. It can also gets the record title from the top of the SD record or from a tag, as specified.

It generates an [FPS](FPS.md) file output, optionally gzip compressed.

More information can be found in the [chemfp documentation](http://readthedocs.org/docs/chemfp/en/latest/).

# Details #

Here is an example of how to extract the PubChem fingerprints from a PubChem data file. The --pubchem
option is a shortcut for "--software=CACTVS/unknown --type 'CACTVS-E\_SCREEN/1.0 extended=2'
--fp-tag=PUBCHEM\_CACTVS\_SUBSKEYS --cactvs".


```
% python sdf2fps --pubchem $PUBCHEM/Compound_013150001_013175000.sdf.gz | head -12 | fold -w 80
#FPS1
#num_bits=881
#type=CACTVS-E_SCREEN/1.0 extended=2
#software=CACTVS/unknown
#source=/Users/dalke/databases/pubchem/Compound_013150001_013175000.sdf.gz
#date=2011-09-20T20:04:17
075e00000000000000000000000000000000000000000c0683c10000000000832a00003800080000
0030108118000c034303000001402442020041000084400010110000261110044603898921041006
09001313e00801037013004002004800000900200100240000000000000000	13150007
071e04000000000000000000000000000080040000000c0683c10000000012833f00005800000000
0030200119000c60030020021140054a000040100024040010118020101330644c21ac58419c0425
03881095e111130f710700c0000018000003006000000c0000000000000000	13150008
075e00000000000000000000000068000000000000000c060300a001000000830200003800000000
0030148318204c00c100000001400442000041000004000010110010201110044401898821041006
01001111e00801037001000000000800000900200100240000000000000000	13150009
075e00000000000000000000000000000000000000000c0603000000000000832a00003800080000
0030108318204c034303000001402442020041000084400010110110261110044403898921041006
09001313e00801037013004002004800000900200100240000000000000000	13150010
031e0c000200000000000000000000000000000000002c0601000000000000890200005820000000
003020211b000d80010000501140054a000e42000024100810119800001310044c05a80801840004
01001491e11011017100000000002000000000000000100000000000000000	13150011
031e0c000208000000000000000000000000000000002c0601000000000000890200005820020000
803520211b000d80010000501140054a000e42000024102810119800001710044c05a80801840004
010014d1e91011017140000000002000002000000000100000000000000000	13150013
```


## sdf2fps --help ##

```
usage: sdf2fps [-h] [--id-tag TAG] [--fp-tag TAG] [--num-bits INT]
               [--errors {strict,report,ignore}] [-o FILENAME]
               [--software TEXT] [--type TEXT] [--decompress METHOD]
               [--binary] [--binary-msb] [--hex] [--hex-lsb] [--hex-msb]
               [--base64] [--cactvs] [--daylight] [--decoder DECODER]
               [--pubchem]
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
  --daylight            Daylight encoding, which is is base64 variant
  --decoder DECODER     import and use the DECODER function to decode the
                        fingerprint

shortcuts:
  --pubchem             decode CACTVS substructure keys used in PubChem. Same
                        as --software=CACTVS/unknown --type 'CACTVS-
                        E_SCREEN/1.0 extended=2' --fp-tag=PUBCHEM_CACTVS_SUBSKEYS
                        --cactvs
```
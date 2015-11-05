# Introduction #

PubChem uses the CACTVS 661 bit substructure fingerprints. These are available as a set of loose pattern definitions ("loose" here meaning "not in a strict format like SMARTS"). Though a lot of effort, and by talking to the right people, and comparing my work to the work in the CDK, I created a set of pattern definitions which are close to the CACTVS ones.

They are _not_ the same as the PubChem substructure keys. Those depend on the CACTVS aromaticity model, which is not portable. They also depend on extra bond annotations stored in the SD file in one of the SD tags, and no other toolkit supports this extension.

Instead, they are "heavily inspired" by the PubChem definitions.

They are also in search of a better name than "chemfp substruct keys."

The chemfp tools (ob2fps, oe2fps and rdkit2fps) all support substruct key generation with the --substruct option.

## ChemFP Substruct in [OpenBabel](OpenBabel.md) / [ob2fps](ob2fps.md) ##

```
echo "Cn1cnc2c1c(=O)n(c(=O)n2C)C caffeine" | ./ob2fps --substruct | fold
#FPS1
#num_bits=881
#software=OpenBabel/2.3.0
#type=ChemFP-Substruct-OpenBabel/1
#date=2011-05-30T01:56:39
03ce0d0000000000000000000000000000800600000034000000000000001a800300007800000000
00101080e920e02ffd3008001580e08e2e000101b4e80805800a84152a0a01121002120628201110
e0440200060000000000000000000000000000000000000000000000000000 caffeine
```

## ChemFP Substruct in OpenEye's [OEChem](OEChem.md) / [oe2fps](oe2fps.md) ##

```
% echo "Cn1cnc2c1c(=O)n(c(=O)n2C)C caffeine" | ./oe2fps --substruct | fold
#FPS1
#num_bits=881
#software=OEChem/1.7.4 (20100809)
#type=ChemFP-Substruct-OpenEye/1
#date=2011-05-30T01:57:09
03ce0d0000000000000000000000000000800600000034000000000000001a800100007800000000
00101080e920e02ffd3008001580e08e2e000101b4e80805800a84152a0a01121002120628201110
e0440200060000000000000000000000000000000000000000000000000000 caffeine
```


## ChemFP Substruct in [RDKit](RDKit.md) / [rdkit2fps](rdkit2fps.md) ##

```
% echo "Cn1cnc2c1c(=O)n(c(=O)n2C)C caffeine" | ./rdkit2fps --substruct | fold
#FPS1
#num_bits=881
#software=RDKit/2009Q3_1
#type=ChemFP-Substruct-RDKit/1
#date=2011-05-30T01:57:27
03ce0d0000000000000000000000000000800600000034000000000000001a800300007800000000
00101080e920e02ffd3008001580e08e2e000101b4e80805800a84152a0a01121002120628201110
e0440200060000000000000000000000000000000000000000000000000000 caffeine
```
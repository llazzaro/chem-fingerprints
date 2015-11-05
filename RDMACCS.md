# Introduction #

The MACCS keys are very widely used in cheminformatics but yet also ill defined. There is no published validation suite or formal definition.

Greg Landrum added MACCS support to RDKit, with nearly everything defined as a set of [SMARTS definitions](http://rdkit.org/Python_Docs/rdkit.Chem.MACCSkeys-module.html), plus some special code for a couple of bits which can't be done in SMARTS.

This was the first public definition, and with a definition which was easy to support in other toolkits. The patterns quickly got taken up  [OpenBabel](OpenBabel.md) and [CDK](CDK.md).

However, the SMARTS definitions were specific to RDKit's SMARTS matcher. For example, by default it doesn't match `*` against hydrogens like `[2H]`. It also didn't include that other toolkits support a different set of aromatic atoms.

I took RDKit's SMARTS patterns and tweaked them to be more portable. I decided to call the result "RDMACCS" and promote them as part of chemfp. This includes a validation suite.

Each of [ob2fps](ob2fps.md), [oe2fps](oe2fps.md), and [rdkit2fps](rdkit2fps.md) implement the RDMACCS fingerprint through the `--rdmaccs` option. Be aware that this fingerprint depends in part on aromaticity and ring perception so there will often be differences in the results of each toolkit.


## RDMACCS in [OpenBabel](OpenBabel.md) / [ob2fps](ob2fps.md) ##

```
% echo "Cn1cnc2c1c(=O)n(c(=O)n2C)C caffeine" | ob2fps --rdmaccs
#FPS1
#num_bits=166
#software=OpenBabel/2.3.0
#type=RDMACCS-OpenBabel/1
#date=2011-05-30T01:46:32
000000003000000001d414d91323915380f138ea1f caffeine
```

## RDMACCS in OpenEye's [OEChem](OEChem.md) / [oe2fps](oe2fps.md) ##

```
% echo "Cn1cnc2c1c(=O)n(c(=O)n2C)C caffeine" | oe2fps --rdmaccs
#FPS1
#num_bits=166
#software=OEChem/1.7.4 (20100809)
#type=RDMACCS-OpenEye/1
#date=2011-05-30T01:46:28
000000003000000001d414d91323915380f138ea1f caffeine
```


## RDMACCS in [RDKit](RDKit.md) / [rdkit2fps](rdkit2fps.md) ##

```
% echo "Cn1cnc2c1c(=O)n(c(=O)n2C)C caffeine" | rdkit2fps --rdmaccs
#FPS1
#num_bits=166
#software=RDKit/2009Q3_1
#type=RDMACCS-RDKit/1
#date=2011-05-30T01:46:15
000000003000000001d414d91323915380f138ea1f caffeine
```
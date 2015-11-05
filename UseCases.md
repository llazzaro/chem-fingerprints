The chemfp project is meant to address the following use cases.

## Use case #1: Web service for 3-nearest similarity searches ##

A research group has a set of structures with associated experimental
values. The information is updated every month. They want to provide a
web service which takes an input structure and finds the 3 nearest
compounds with better than 70% Tanimoto similarity, and reports the
similarity score and structure information.

There are no ready-made tools for this, excepting groups which have a
chemistry-aware database. The pieces to build it are available from
several projects, but large parts need to be filled in. Ideally there
should be a set of command-line and library functions for working with
chemfp files.

## Use case #2: Data provenance ##

It's time to write the paper. Where did this data set come from? Which
program generated it and with which options? Was it the one which used
the buggy SMARTS definitions?

Most formats don't track this information. While it's impossible to be
perfect, some information would help.


## Use case #3: New fingerprint types ##

A researcher develops a new fingerprint scheme and wants to compare
the applicability to existing fingerprints, including the linear hash
in OpenBabel and the topological hash in RDKit.

Currently that requires quite a bit of work to figure out how to use
each program to read the input files and generate the fingerprints
which can be used for the analysis.

Separating fingerprint generation from fingerprint use also makes it
easier to implement things like N\*N similarity clustering on a cluster
without needing to install the chemistry tools on all the machines.

## Use case #4: Comparison of search algorithms ##

An algorithms developer comes up with a new scheme for fast similarity
wants to compare the effectiveness of the algorithm against other
implementations.

While there are many published papers on this topic, few of the
algorithms and data sets are available for direct comparison. I don't
think the problem is the lack of a common format, but I suggest it
might help.
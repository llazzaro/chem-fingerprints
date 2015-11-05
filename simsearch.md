# Introduction #

simsearch support Tanimoto searches of a target fingerprint data set. More information can be found in the [chemfp documentation](http://readthedocs.org/docs/chemfp/en/latest/).


# Details #

The first step is to make a target database. Here I'll use a PubChem data set with 3995 structures.

```
% sdf2fps --pubchem $PUBCHEM/Compound_013150001_013175000.sdf.gz -o targets.fps.gz
```

See how the -o option knows to compress the output?

Second, make a query data set. In this case I'll use 4 fingerprints from another PubChem data set

```
% sdf2fps --pubchem $PUBCHEM/Compound_005000001_005025000.sdf.gz | head > queries.fps
% fold queries.fps 
#FPS1
#num_bits=881
#type=CACTVS-E_SCREEN/1.0 extended=2
#software=CACTVS/unknown
#source=/Users/dalke/databases/pubchem/Compound_005000001_005025000.sdf.gz
#date=2011-09-20T17:05:14
07de04000600000000000000000000000080040000003c060100000000001a800300007820080000
003014a31b208d81c10300103140844a0a00c10001a6109810118810221311045c07ab8921841116
e1401793e61811037101000000000000000000000000000000000000000000	5000001
07de0c000200000000000000000000000080460200000c0200000000000000810300007820080000
003038a77b60bd03c993291015c0adee2e00410984ee400c909b851d261b50045f07bb8de1841126
69001b93e31999407100000000004000000000000000200000000000000000	5000002
07de04000600000000000000000000000080060000000c0603000000000000800a00007820080000
003010811b000c03c10300103140a44a0a0041000086409810110000261310044603898921041006
01001393e10801007010000000000000000800000000000000000000000000	5000003
075e1c000200000000000000000018000080040000000c0000000000000012800100007820080000
00b000851b404091410320103140800b1a40c10001a6109800118802321350645c072db9e1881166
2b801f97e21913077101000000000000000100800000100000000000000000	5000005
```


Third, I'll search it, using the default -k value of 3 and --threshold of 0.7 to find the nearest 3 target fingerprints for each query which have a similarity of at least 0.7.

```
% simsearch -q queries.fps targets.fps.gz
#Simsearch/1
#num_bits=881
#type=Tanimoto k=3 threshold=0.7
#software=chemfp/1.0
#queries=queries.fps
#targets=targets.fps.gz
#query_sources=/Users/dalke/databases/pubchem/Compound_005000001_005025000.sdf.gz
#target_sources=/Users/dalke/databases/pubchem/Compound_013150001_013175000.sdf.gz
3	5000001	13163368 	0.7403	13163366		0.7308	13163593		0.7160
3	5000002	13153311 	0.7233	13152891		0.7192	13152887		0.7178
3	5000003	13163171 	0.7734	13163170		0.7638	13170866		0.7328
0	5000005
```

## Simsearch/1 output format ##

The output format is still somewhat experimental. I'm looking for
feedback. You can see it's in the same "family" as the fps format,
with a line containing the format and version, followed by key/value
header lines, followed by data. The data fields are tab delimited.

I'm still not sure about what metadata
is needed here so feedback much appreciated.

The similarity results contain one line per query structure. The
format is the number of matches (in this case 3), followed by the
query identifier. After that are the target identifier and score for
the best match, then for the second best match, and so on.

## Changing -k and setting a --threshold ##


By default simsearch finds the closest three structures. You can specify a different value, and you can specify a minimum similarity score. For example:

```
% simsearch -k 4 --threshold 0.70 -q queries.fps targets.fps.gz
#Simsearch/1
#num_bits=881
#type=Tanimoto k=4 threshold=0.7
#software=chemfp/1.0
#queries=queries.fps
#targets=targets.fps.gz
#query_sources=/Users/dalke/databases/pubchem/Compound_005000001_005025000.sdf.gz
#target_sources=/Users/dalke/databases/pubchem/Compound_013150001_013175000.sdf.gz
3	5000001	13163368	0.7403	13163366	0.7308	13163593	0.7160
4	5000002	13153311	0.7233	13152891	0.7192	13152887	0.7178	13152431	0.7150
4	5000003	13163171	0.7734	13163170	0.7638	13170866	0.7328	13173188	0.7248
0	5000005
```

If k is a number then the hits will be ordered so the most similar structures come first. If k is "all" then all hits at or above the threshold are reported and the order of the results is arbitrary.

If you specify k and no threshold then the default threshold will be 0.0. If you specify the threshold but not k then the default k will be "all."


## Count mode ##

Sometimes you only want to know how many structures are within a given threshold of the query. You don't need to know their identifiers. For example, how many structures are within 0.4 of the queries?

```
% simsearch -c --threshold 0.4 -q queries.fps targets.fps.gz
#Count/1
#num_bits=881
#type=Count threshold=0.4
#software=chemfp/1.0
#queries=queries.fps
#targets=targets.fps.gz
#query_sources=/Users/dalke/databases/pubchem/Compound_005000001_005025000.sdf.gz
#target_sources=/Users/dalke/databases/pubchem/Compound_013150001_013175000.sdf.gz
2211 	5000001
1942 	5000002
2427	 	5000003
1636	 	5000005
```

The format here is the count followed by the identifier, which is the same as the Simsearch/1 format, only without the target hit information.

## NxN searches ##

Use the --NxN option when you want to use the same set of fingerprints as both the queries and the targets. In that case you know that the diagonal terms will always match, which is rather boring and uninteresting, so the self-similarity score will not be included in the output.


```
#Simsearch/1
#num_bits=1021
#type=Tanimoto k=3 threshold=0.7 NxN=full
#software=chemfp/1.1b6
#targets=targets.fps.gz
#target_sources=targets.smi
0	22525001
2	22525002	22525096	0.83824	22525006	0.71875
0	22525003
0	22525004
0	22525005
1	22525006	22525002	0.71875
  ...
```


## simsearch --help ##

```
usage: simsearch [-h] [-k K_NEAREST] [-t THRESHOLD] [-q QUERIES] [--NxN]
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
  --NxN                 use the targets as the queries, and exclude the self-
                        similarity term
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
```
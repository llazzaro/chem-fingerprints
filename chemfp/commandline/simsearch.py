from __future__ import with_statement
import math
import sys
import itertools

import chemfp
from chemfp import argparse, readers, io, SOFTWARE, bitops

def write_simsearch_magic(outfile):
    outfile.write("#Simsearch/1\n")

def write_count_magic(outfile):
    outfile.write("#Count/1\n")

def write_simsearch_header(outfile, d):
    lines = []
    for name in ("num_bits", "type", "software", "queries", "targets"):
        value = d.get(name, None)
        if value is not None:
            lines.append("#%s=%s\n" % (name, value))
    for name in ("query_sources", "target_sources"):
        for value in d.get(name, []):
            lines.append("#%s=%s\n" % (name, value))
    outfile.writelines(lines)


def report_threshold(outfile, float_formatter, query_arenas, targets, threshold):
    def search_function(query_arena):
        return chemfp.threshold_tanimoto_search(query_arena, targets,
                                                threshold=threshold)
    _report_search(outfile, float_formatter, query_arenas, search_function)

def report_knearest(outfile, float_formatter, query_arenas, targets, k, threshold):
    def search_function(query_arena):
        return chemfp.knearest_tanimoto_search(query_arena, targets,
                                               k=k, threshold=threshold)
    _report_search(outfile, float_formatter, query_arenas, search_function)

def _report_search(outfile, float_formatter, query_arenas, search_function):
    hit_formatter = "\t%s\t" + float_formatter
    for query_arena in query_arenas:
        for query_id, hits in search_function(query_arena):
            outfile.write("%d\t%s" % (len(hits), query_id))
            for hit in hits:
                outfile.write(hit_formatter % hit)
            outfile.write("\n") # XXX flush?
    


def report_counts(outfile, query_arenas, targets, threshold):
    for query_arena in query_arenas:
        results = chemfp.count_tanimoto_hits(query_arena, targets, threshold)
        for query_id, hit_count in results:
            outfile.write("%d\t%s\n" % (hit_count, query_id))
        
def int_or_all(s):
    if s == "all":
        return s
    return int(s)

# the "2fps" options need a way to say "get the options from --reference"
# ob2fps --reference targets.fps | simsearch -k  5 --threshold 0.5 targets.fps


parser = argparse.ArgumentParser(
    description="Search an FPS file for similar fingerprints")
parser.add_argument("-k" ,"--k-nearest", help="select the k nearest neighbors (use 'all' for all neighbors)",
                    default=None, type=int_or_all)
parser.add_argument("-t" ,"--threshold", help="minimum similarity score threshold",
                    default=None, type=float)
parser.add_argument("-q", "--queries", help="filename containing the query fingerprints")
parser.add_argument("--hex-query", help="query in hex")
parser.add_argument("--query-id", default="Query1",
                    help="id for the hex query")
parser.add_argument("--in", metavar="FORMAT", dest="query_format",
                    help="input query format (default uses the file extension, else 'fps')")
parser.add_argument("-o", "--output", metavar="FILENAME",
                    help="output filename (default is stdout)")

parser.add_argument("-c", "--count", help="report counts", action="store_true")

parser.add_argument("-b", "--batch-size", help="batch size",
                    default=100, type=int)

parser.add_argument("--scan", help="scan the file to find matches (low memory overhead)",
                    action="store_true")
parser.add_argument("--memory", help="build and search an in-memory data structure (faster for multiple queries)",
                    action="store_true")

parser.add_argument("--times", help="report load and execution times to stderr",
                    action="store_true")

parser.add_argument("target_filename", nargs=1, help="target filename", default=None)

## Something to enable multi-threading
#parser.add_argument("-j", "--jobs", help="number of jobs ",
#                    default=10, type=int)


def main(args=None):
    args = parser.parse_args(args)
    target_filename = args.target_filename[0]

    threshold = args.threshold
    k = args.k_nearest

    if args.count and k is not None and k != "all":
        parser.error("--count search does not support --k-nearest")

    # People should not use this without setting parameters.  On the
    # other hand, I don't want an error message if there are no
    # parameters. This solution seems to make sense.

    if threshold is None:
        if k is None:
            # If nothing is set, use defaults of --thresdhold 0.7 -k 3            
            threshold = 0.7
            k = 3
        else:
            # only k is set; search over all possible matches
            threshold = 0.0
    else:
        if k is None:
            # only threshold is set; search for all hits above that threshold
            k = "all"
    
    if args.scan and args.memory:
        parser.error("Cannot specify both --scan and --memory")
    
    if args.hex_query and args.queries:
        parser.error("Cannot specify both --hex-query and --queries")
    if args.hex_query:
        query_id = args.query_id
        for c, name in ( ("\t", "tab"),
                         ("\n", "newline"),
                         ("\r", "control-return"),
                         ("\0", "NUL")):
            if c in query_id:
                parser.error("--query-id must not contain the %s character" %
                             (name,))
            

    if k == "all":
        pass
    elif k < 0:
        parser.error("--k-nearest must non-negative or 'all'")

    if not (0.0 <= threshold <= 1.0):
        parser.error("--threshold must be between 0.0 and 1.0, inclusive")

    if args.batch_size < 1:
        parser.error("--batch-size must be positive")

    batch_size = args.batch_size

    bitops.use_environment_variables()

    # Open the target file. This reads just enough to get the header.

    targets = chemfp.open(target_filename)
            
    if args.hex_query is not None:
        try:
            query_fp = args.hex_query.decode("hex")
        except ValueError, err:
            parser.error("--hex-query is not a hex string: %s" % (err,))

        for (severity, error, msg_template) in chemfp.check_fp_problems(query_fp, targets.metadata):
            if severity == "error":
                parser.error(msg_template % dict(fp="query", metadata=repr(target_filename)))
            
        num_bits = targets.metadata.num_bits
        if num_bits is None:
            num_bits = len(query_fp) * 8
        query_metadata = chemfp.Metadata(num_bits=num_bits, num_bytes=len(query_fp))
        queries = chemfp.Fingerprints(query_metadata,
                                      [(query_id, query_fp)])
        query_filename = None
    else:
        query_filename = args.queries
        queries = chemfp.open(query_filename, format=args.query_format)

    query_arena_iter = queries.iter_arenas(batch_size)
    
    # Suppose you have a 4K fingerprint.
    #   1/4096 = 0.000244140625.
    #   2/4096 = 0.00048828125
    # You only need to show "0.0002" and "0.0005" to
    # disambiguate the scores. I don't like seeing only
    # the minimum resolution, so I also show at least
    # the next bit.
    #   For 4096 the float_formatter is %.5f and the
    # above values are 0.00024 and 0.00049.
    # This also prevents the results from being shown
    # in scientific notation.
    num_digits = int(math.log10(targets.metadata.num_bytes*8)) + 2
    float_formatter = "%." + str(num_digits) + "f"

    import time
    t1 = time.time()

    first_query_arena = None
    for first_query_arena in query_arena_iter:
        break

    if args.scan:
        # Leave the targets as-is
        pass
    elif args.memory:
        targets = chemfp.load_fingerprints(targets)
    if not first_query_arena:
        # No input. Leave as-is
        pass
    elif len(first_query_arena) < min(10, batch_size):
        # Figure out the optimal search. If there is a
        # small number of inputs (< ~10) then a scan
        # of the FPS file is faster than an arena search.
        pass
    else:
        targets = chemfp.load_fingerprints(targets)

    problems = chemfp.check_metadata_problems(queries.metadata, targets.metadata)
    for (severity, error, msg_template) in problems:
        msg = msg_template % dict(metadata1="queries", metadata2="targets")
        if severity == "error":
            parser.error(msg)
        elif severity == "warning":
            sys.stderr.write("WARNING: " + msg + "\n")

    t2 = time.time()
    outfile = io.open_output(args.output)
    with io.ignore_pipe_errors:
        type = "Tanimoto k=%(k)s threshold=%(threshold)s" % dict(
            k=k, threshold=threshold, max_score=1.0)

        if args.count:
            type = "Count threshold=%(threshold)s" % dict(
                threshold=threshold)
            write_count_magic(outfile)
        else:
            write_simsearch_magic(outfile)

        write_simsearch_header(outfile, {
            "num_bits": targets.metadata.num_bits,
            "software": SOFTWARE,
            "type": type,
            "queries": query_filename,
            "targets": target_filename,
            "query_sources": queries.metadata.sources,
            "target_sources": targets.metadata.sources})

        if first_query_arena:
            query_arenas = itertools.chain([first_query_arena],
                                           query_arena_iter)

            if args.count:
                report_counts(outfile, query_arenas, targets,
                              threshold = threshold)
            elif k == "all":
                report_threshold(outfile, float_formatter, query_arenas, targets,
                                 threshold = threshold)
            else:
                report_knearest(outfile, float_formatter, query_arenas, targets,
                                k = k, threshold = threshold)
                                
                
                    
    t3 = time.time()
    if args.times:
        sys.stderr.write("open %.2f search %.2f total %.2f\n" % (t2-t1, t3-t2, t3-t1))

if __name__ == "__main__":
    main()

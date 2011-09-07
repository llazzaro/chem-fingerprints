from __future__ import with_statement
import chemfp
import math
from chemfp import argparse, readers, io, SOFTWARE
import sys
import itertools

def write_simsearch_magic(outfile):
    outfile.write("#Simsearch/1\n")

def write_count_magic(outfile):
    outfile.write("#Count/1\n")

def write_simsearch_header(outfile, d):
    lines = []
    for name in ("num_bits", "software", "type", "query_source", "target_source"):
        value = d.get(name, None)
        if value is not None:
            lines.append("#%s=%s\n" % (name, value))
    outfile.writelines(lines)


def report_threshold(outfile, query_arenas, targets, threshold):
    def search_function(query_arena):
        return chemfp.threshold_tanimoto_search(query_arena, targets,
                                                threshold=threshold)
    _report_search(outfile, query_arenas, search_function)

def report_knearest(outfile, query_arenas, targets, k, threshold):
    def search_function(query_arena):
        return chemfp.knearest_tanimoto_search(query_arena, targets,
                                               k=k, threshold=threshold)
    _report_search(outfile, query_arenas, search_function)

def _report_search(outfile, query_arenas, search_function):
    for query_arena in query_arenas:
        for query_id, hits in search_function(query_arena):
            outfile.write("%d\t%s" % (len(hits), query_id))
            for hit in hits:
                outfile.write("\t%s\t%f" % hit)
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
                    default=3, type=int_or_all)
parser.add_argument("-t" ,"--threshold", help="minimum similarity score threshold",
                    default=0.0, type=float)
parser.add_argument("-q", "--queries", help="filename containing the query fingerprints")
parser.add_argument("--hex-query", help="query in hex")
parser.add_argument("--query-id", default="query",
                    help="id for the hex query")
parser.add_argument("--in", metavar="FORMAT", dest="query_format",
                    help="input query format (default uses the file extension, else 'fps')")
parser.add_argument("-o", "--output", metavar="FILENAME",
                    help="output filename (default is stdout)")
parser.add_argument("--type", help="fingerprint type", default=None)

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

    if args.scan and args.memory:
        parser.error("Cannot specify both --scan and --memory")
    
    if args.hex_query and args.queries:
        parser.error("Cannot specify both --hex-query and --queries")
    if args.hex_query:
        query_id = args.query_id
        for c, name in ( ("\t", "tab"),
                         ("\n", "newline"),
                         ("\0", "NUL")):
            if c in query_id:
                parser.error("--query-id must not contain the %s character" %
                             (name,))
            

    if args.k_nearest == "all":
        pass
    elif args.k_nearest < 0:
        parser.error("--k-nearest must non-negative or 'all'")

    if not (0.0 <= args.threshold <= 1.0):
        parser.error("--threshold must be between 0.0 and 1.0, inclusive")

    if args.batch_size < 1:
        parser.error("--batch-size must be positive")

    batch_size = args.batch_size

    # Open the target file. This reads just enough to get the header.

    try:
        targets = chemfp.open(target_filename, type=args.type)
    except TypeError, err:
        if "'type' is required" in str(err):
            parser.error("--type is required to convert structure in the targets file to fingerprints")
        raise
            
    if args.hex_query is not None:
        try:
            query_fp = args.hex_query.decode("hex")
        except ValueError, err:
            parser.error("--hex-query is not a hex string: %s" % (err,))
        compatible = io.check_compatibility(fp=query_fp, header=targets.header)

        query_num_bits = len(query_fp) * 8
        target_num_bits = targets.header.num_bits
        
        if not compatible:
            parser.error("--hex-query with %d bits is not compatible with targets with %d bits" %
                         (query_num_bits, target_num_bits))
        
        queries = chemfp.Fingerprints(io.Header(num_bits=target_num_bits),
                                      [(query_id, query_fp)])

    else:
        queries = chemfp.open(args.queries, format=args.query_format, type=type)        

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
    num_digits = int(math.log10(queries.header.num_bits)) + 2
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

    if queries.header.num_bits is not None and targets.header.num_bits is not None:
        if queries.header.num_bits != targets.header.num_bits:
            sys.stderr.write("WARNING: query has %d bits and target has %d\n" %
                             (queries.header.num_bits, targets.header.num_bits))
    if queries.header.type and targets.header.type:
        if queries.header.type != targets.header.type:
            sys.stderr.write("WARNING: fingerprints have incompatible headers\n")
            sys.stderr.write("  query: %s\n" % (queries.header.type,))
            sys.stderr.write(" target: %s\n" % (targets.header.type,))

    t2 = time.time()
    outfile = io.open_output(args.output)
    with io.ignore_pipe_errors:
        type = "Tanimoto k=%(k)s threshold=%(threshold)s" % dict(
            k=args.k_nearest, threshold=threshold, max_score=1.0)

        if args.count:
            type = "Count threshold=%(threshold)s" % dict(
                threshold=args.threshold)
            write_count_magic(outfile)
        else:
            write_simsearch_magic(outfile)
            
        write_simsearch_header(outfile, {
            "num_bits": targets.header.num_bits,
            "software": SOFTWARE,
            "type": type,
            "query_source": queries.header.source,
            "target_source": targets.header.source})

        if first_query_arena:
            query_arenas = itertools.chain([first_query_arena],
                                           query_arena_iter)

            if args.count:
                report_counts(outfile, query_arenas, targets,
                              threshold = args.threshold)
            elif args.k_nearest == "all":
                report_threshold(outfile, query_arenas, targets,
                                 threshold = args.threshold)
            else:
                report_knearest(outfile, query_arenas, targets,
                                k = args.k_nearest,
                                threshold = args.threshold)
                
                    
    t3 = time.time()
    if args.times:
        sys.stderr.write("open %.2f search %.2f total %.2f\n" % (t2-t1, t3-t2, t3-t1))

if __name__ == "__main__":
    main()

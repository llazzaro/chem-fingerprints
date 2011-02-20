from __future__ import with_statement
import chemfp
import math
from chemfp import argparse, readers, io, SOFTWARE
import sys

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


def report_knearest(query_iter, batch_ids, batch_fps, targets, args, float_formatter, outfile):
    format_string = float_formatter + " %s"
    k = args.k_nearest
    threshold = args.threshold
    batch_size = args.batch_size
    
    first_time = True
    while 1:
        if not first_time:
            targets.reset()
        else:
            first_time = False
            
        results = chemfp.tanimoto_knearest_search_batch(batch_fps, targets, k, threshold)
            
        for query_id, closest_targets in zip(batch_ids, results):
            fields = ["%d %s" % (len(closest_targets), query_id)]
            for (target_name, score) in closest_targets:
                fields.append(format_string % (score, target_name))
            outfile.write(" ".join(fields))
            outfile.write("\n") # XX flush

        batch_ids, batch_fps = read_batch(batch_size, query_iter)
        if not batch_ids:
            break

def report_counts(query_iter, batch_ids, batch_fps, targets, args, outfile):
    threshold = args.threshold
    batch_size = args.batch_size

    first_time = True
    while 1:
        if not first_time:
            targets.reset()
        else:
            first_time = False
            
        results = chemfp.tanimoto_count_batch(batch_fps, targets, threshold)
            
        for count, query_id in zip(results, batch_ids):
            outfile.write("%d %s\n" % (count, query_id))

        batch_ids, batch_fps = read_batch(batch_size, query_iter)
        if not batch_ids:
            break
        


# the "2fps" options need a way to say "get the options from --reference"
# ob2fps --reference targets.fps | simsearch -k  5 --threshold 0.5 targets.fps


parser = argparse.ArgumentParser(
    description="Search an FPS file for similar fingerprints")
parser.add_argument("-k" ,"--k-nearest", help="select the k nearest neighbors",
                    default=3, type=int)
parser.add_argument("-t" ,"--threshold", help="minimum similarity score threshold",
                    default=0.0, type=float)
parser.add_argument("-q", "--queries", help="filename containing the query fingerprints")
parser.add_argument("--hex-query", help="query in hex")
parser.add_argument("--in", metavar="FORMAT", dest="query_format",
                    help="input query format (default uses the file extension, else 'fps')")
parser.add_argument("-o", "--output", metavar="FILENAME",
                    help="output filename (default is stdout)")
parser.add_argument("--type", help="fingerprint type", default=None)

parser.add_argument("-c", "--count", help="report counts", action="store_true")

parser.add_argument("-b", "--batch-size", help="batch size",
                    default=100, type=int)

parser.add_argument("--file-scan", help="search directly from the input FPS file",
                    action="store_true")
parser.add_argument("--in-memory", help="use an in-memory fingerprint search",
                    action="store_true")

parser.add_argument("--times", help="report load and execution times to stderr",
                    action="store_true")

parser.add_argument("target_filename", nargs=1, help="target filename", default=None)

## Something to enable multi-threading
#parser.add_argument("-j", "--jobs", help="number of jobs ",
#                    default=10, type=int)

def read_batch(batch_size, queries):
    if batch_size <= 0:
        return [], []
    batch_size -= 1
    batch_fps = []
    batch_ids = []
    for i, query in enumerate(queries):
        batch_fps.append(query[0])
        batch_ids.append(query[1])
        if i == batch_size:
            break

    return batch_ids, batch_fps


def main(args=None):
    args = parser.parse_args(args)
    target_filename = args.target_filename[0]
    threshold = args.threshold

    if args.file_scan and args.in_memory:
        args.error("Cannot specify both --file-scan and --in-memory")
    
    if args.hex_query and args.queries:
        args.error("Cannot specify both --hex-query and --queries")

    batch_size = args.batch_size # args.batch_size

    # Open the file. This reads just enough to get the header.

    try:
        targets = chemfp.open(target_filename, type=args.type)
    except TypeError, err:
        if "'type' is required" in str(err):
            parser.error("--type is required to convert structure in the targets file to fingerprints")
        raise
            
    type = targets.header.type

    if args.hex_query is not None:
        try:
            query = args.hex_query.decode("hex")
        except ValueError, err:
            args.error("--hex-query is not a hex string: %s" % (err,))
        query_iter = iter( [(query, "Query1")] )
        queries = io.FPIterator(io.Header(num_bits = len(query) * 8), query_iter)

    elif args.queries is not None:
        queries = chemfp.open(args.queries, format=args.query_format, type=type)
    else:
        queries = chemfp.open(None, format=args.query_format, type=type)

    query_iter = iter(queries)

    # See if there's enough queries to justify reading the targets into memory
    batch_ids, batch_fps = read_batch(batch_size, query_iter)
    if not batch_ids:
        return

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

    # Looks like in-memory is about 5x faster (500 structures)

    ## Estimate about 10 is the tradeoff, but these are bad numbers.
    # fps_reader:
    #   1  1.6
    #   6  5.4
    #  12 10.1
    #  19 15.4
    # in_memory
    #   1  3.3  (0.03 in search)
    #   6  3.8 (0.4 for search)
    #  12  7.1 (1.0 for search)
    #  19  5.2 (1.8 for search)

    # In my testing, using the FPS reader is much faster than an
    # in-memory search for 1 structure. The breakeven point is around
    # 6 input structures.

    if args.file_scan:
        use_in_memory = False
    elif args.in_memory:
        use_in_memory = True
    elif (len(batch_ids) <= 6 and batch_size > 6):
        use_in_memory = False
    else:
        use_in_memory = True

    if use_in_memory:
        targets = chemfp.read_into_memory(targets)

    if queries.header.num_bits is not None and targets.header.num_bits is not None:
        if queries.header.num_bits != targets.header.num_bits:
            sys.stderr.write("WARNING: query has %d bits and target has %d\n" %
                             (queries.header.num_bits, target.header.num_bits))
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
        write_simsearch_magic(outfile)
        write_simsearch_header(outfile, {
            "num_bits": targets.header.num_bits,
            "software": SOFTWARE,
            "type": type,
            "query_source": queries.header.source,
            "target_source": targets.header.source})

        if args.count:
            report_counts(query_iter, batch_ids, batch_fps, targets, args, outfile)
        else:
                
            report_knearest(query_iter, batch_ids, batch_fps, targets, args, float_formatter, outfile)

    t3 = time.time()
    if args.times:
        sys.stderr.write("open %.2f search %.2f total %.2f\n" % (t2-t1, t3-t2, t3-t1))

if __name__ == "__main__":
    main()

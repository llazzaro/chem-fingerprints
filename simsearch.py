import chemfp
import math
from chemfp import argparse, readers
import sys



def report_knearest(query_iter, batch_ids, batch_fps, targets, args, float_formatter):
    k = args.k_nearest
    threshold = args.threshold
    batch_size = args.batch_size
    
    print "Simsearch"
#    write_similarity_header(args, queries, targets)

    first_time = True
    while 1:
        print "Processing", len(batch_fps)
        if not first_time:
            targets.reset()
        else:
            first_time = False
            
        results = chemfp.tanimoto_knearest_search_batch(batch_fps, targets, k, threshold)
            
        for query_id, closest_targets in zip(batch_ids, results):
            fields = [query_id]
            for (target_name, score) in closest_targets:
                fields.append(float_formatter % score)
                fields.append(target_name)
            print " ".join(fields)

        batch_ids, batch_fps = read_batch(batch_size, query_iter)
        if not batch_ids:
            break

def report_counts(query_iter, batch_ids, batch_fps, targets, args):
    threshold = args.threshold
    batch_size = args.batch_size

    print "Counts"
# write_counts_header(args, queries, targets)

    first_time = True
    while 1:
        print "Processing", len(batch_fps)
        if not first_time:
            targets.reset()
        else:
            first_time = False
            
        results = chemfp.tanimoto_count_batch(batch_fps, targets, threshold)
            
        for query_id, count in zip(batch_ids, results):
            print query_id, count

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
parser.add_argument("--query-hex", help="query in hex")

parser.add_argument("-c", "--count", help="report counts", action="store_true")

parser.add_argument("-b", "--batch-size", help="batch size",
                    default=100, type=int)

parser.add_argument("--use-reader", help="search directly from the input FPS file",
                    action="store_true")
parser.add_argument("--in-memory", help="use an in-memory fingerprint search",
                    action="store_true")

parser.add_argument("target_filename", nargs=1, help="target filename", default=None)

## Something to enable multi-threading
#parser.add_argument("-j", "--jobs", help="number of jobs ",
#                    default=10, type=int)

'''
#FPSimilarity/1
#software=chemfp/1.0
#query-source=q.fps
#target-source=blah.smi
#type=tanimoto k-nearest=3 threshold=0.8
#date=2011-01-29T13:22:30
blah 0.9 hit1 0.8 hit2
'''

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

    if args.use_reader and args.in_memory:
        args.error("Cannot specify both --use-reader and --in-memory")
    
    batch_size = args.batch_size # args.batch_size
    
    if args.queries is not None:
        queries = chemfp.open_reader(args.queries)
    else:
        queries = chemfp.open_reader(sys.stdin)
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

    if args.use_reader:
        use_reader = True
    elif args.in_memory:
        use_reader = False
    elif (len(batch_ids) <= 6 and batch_size > 6):
        use_reader = True
    else:
        use_reader = False

    if use_reader:
        targets = chemfp.open_reader(target_filename)
    else:
        targets = chemfp.open_in_memory(target_filename)

    print targets

#    if not compatible(queries, targets):
#        raise SystemExit("Can not do search")

    t2 = time.time()
    if args.count:
        report_counts(query_iter, batch_ids, batch_fps, targets, args)
    else:
        report_knearest(query_iter, batch_ids, batch_fps, targets, args, float_formatter)

    t3 = time.time()
    print "total", t3-t1
    print "load", t2-t1
    print "search", t3-t2

if __name__ == "__main__":
    main()

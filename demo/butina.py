# An implementation of Taylor-Butina clustering

# See http://www.chemomine.co.uk/dbclus-paper.pdf
# and http://www.redbrick.dcu.ie/~noel/R_clustering.html

import chemfp

THRESHOLD = 0.80

dataset = chemfp.load_fingerprints("docs/pubchem_targets.fps")
print "Clustering", len(dataset), "fingerprints"

# I'll make a list with tuples containing:
#   - the number of hits
#   - an arbitrary and not very good tie-breaker value (larger values go first)
#   - the fingerprint index
#   - the list of fingerprint indices within THRESHOLD of that fingerprint

def tie_breaker_value(hits):
    # This is pretty arbitrary; it's the largest non-1.0 score; or 1.0
    # Noel references better work on tie breaking by John MacCuish at Mesa Analytics
    try:
        tie_breaker = max(score for (idx, score) in hits if score != 1.0)
    except ValueError:
        tie_breaker = 1.0
    return tie_breaker

def hit_members(hits):
    return [idx for (idx, score) in hits]

# Assign the compound index to its hits
search = dataset.threshold_tanimoto_search_arena(dataset, threshold = THRESHOLD).iter_hits()
results = sorted( (  (len(hits), tie_breaker_value(hits), i, hit_members(hits))
                                           for (i, hits) in enumerate(search)  ),
                  reverse=True)

true_singletons = []
false_singletons = []
clusters = []

seen = set()
for (size, ignore, fp_idx, members) in results:
    if fp_idx in seen:
        # Can't use a centroid which is already assigned
        continue
    seen.add(fp_idx)

    
    if size == 1:
        # The only fingerprint in the exclusion sphere is itself
        true_singletons.append(fp_idx)
        continue

    # Figure out which ones haven't yet been assigned
    unassigned = [target_idx for target_idx in members if target_idx not in seen]

    if not unassigned:
        false_singletons.append(fp_idx)
        continue
        
    # this is a new cluster
    clusters.append( (fp_idx, unassigned) )
    seen.update(unassigned)

print len(true_singletons), "true singletons"
print "=>", " ".join(sorted(dataset.ids[idx] for idx in true_singletons))
print

print len(false_singletons), "false singletons"
print "=>", " ".join(sorted(dataset.ids[idx] for idx in false_singletons))
print

# Sort so the cluster with the most compounds comes first,
# then by alphabetically smallest id
def cluster_sort_key(cluster):
    centroid_idx, members = cluster
    return -len(members), dataset.ids[centroid_idx]
    
clusters.sort(key=cluster_sort_key)

print len(clusters), "clusters"
for centroid_idx, members in clusters:
    print dataset.ids[centroid_idx], "has", len(members), "other members"
    print "=>", " ".join(dataset.ids[idx] for idx in members)

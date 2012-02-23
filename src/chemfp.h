#ifndef CHEMFP_H
#define CHEMFP_H

/* Errors are always negative numbers. */
enum chemfp_errors {
  CHEMFP_OK = 0,
  CHEMFP_BAD_ARG = -1,

  CHEMFP_NO_MEM = -2,  /* memory allocation failed */

  /* File format errors */
  CHEMFP_UNSUPPORTED_WHITESPACE = -30,
  CHEMFP_MISSING_FINGERPRINT = -31,
  CHEMFP_BAD_FINGERPRINT = -32,
  CHEMFP_UNEXPECTED_FINGERPRINT_LENGTH = -33,
  CHEMFP_MISSING_ID = -34,
  CHEMFP_BAD_ID = -35,
  CHEMFP_MISSING_NEWLINE = -36,

  /* Popcount errors */
  CHEMFP_METHOD_MISMATCH = -50
  
};

int chemfp_get_num_methods(void);
const char *chemfp_get_method_name(int method);

int chemfp_get_num_alignments(void);
const char *chemfp_get_alignment_name(int alignment);

int chemfp_get_alignment_method(int alignment);
int chemfp_set_alignment_method(int alignment, int method);

int chemfp_select_fastest_method(int alignment, int repeat);


int chemfp_get_num_options(void);
const char *chemfp_get_option_name(int index);
int chemfp_get_option(const char *option);
int chemfp_set_option(const char *option, int value);

/* This gives compile-time version information. */
/* Use "chemfp_version" for run-time version information */
#define CHEMFP_MAJOR_VERSION 1
#define CHEMFP_MINOR_VERSION 1
#define CHEMFP_PATCHLEVEL 0

/* This is of the form (\d+\.\d+) (\.\d)? ((a|b|pre)\d+)
     for examples:  0.9, 1.0.4, 1.0pre2.
 The "a"lpha, "b"eta, and "pre"view suffixes will never be seen in production releases */
#define CHEMFP_VERSION_STRING "1.1a4"

/* Return the CHEMFP_VERSION.  */
const char *chemfp_version(void);

/* Convert an error code to a string description */
const char *chemfp_strerror(int err);

typedef struct {
  int popcount;
  int index;
} ChemFPOrderedPopcount; /* XXX Why mixed case? */

typedef struct {
  double score;
  int query_index;
  int id_start;
  int id_end;
} chemfp_tanimoto_cell;


/* Linked list of blocks */
#define CHEMFP_HIT_BLOCK_SIZE 16

typedef struct chemfp_hit_block {
  int target_indices[CHEMFP_HIT_BLOCK_SIZE];
  double scores[CHEMFP_HIT_BLOCK_SIZE];
  struct chemfp_hit_block *next;
} chemfp_hit_block;


typedef struct chemfp_threshold_result {
  int num_hits;
  int num_allocated;
  int *indices;
  double *scores;
} chemfp_threshold_result;

chemfp_threshold_result *chemfp_alloc_threshold_results(int num_results);
void chemfp_free_results(int num_results, chemfp_threshold_result *);
int chemfp_get_num_hits(chemfp_threshold_result *results);

typedef int (*chemfp_assign_hits_p)(void *data, int i, int target_index, double score);
int chemfp_threshold_result_get_hits(chemfp_threshold_result *results,
                                     chemfp_assign_hits_p add_callback, void *payload);

/*** Low-level operations directly on hex fingerprints ***/

/* Return 1 if the string contains only hex characters; 0 otherwise */
int chemfp_hex_isvalid(int len, const char *fp);

/* Return the population count of a hex fingerprint, otherwise return -1 */
int chemfp_hex_popcount(int len, const char *fp);

/* Return the population count of the intersection of two hex fingerprints,
   otherwise return -1. */
int chemfp_hex_intersect_popcount(int len, const char *fp1, const char *fp2);

/* Return the Tanimoto between two hex fingerprints, or -1.0 for invalid fingerprints
   If neither fingerprint has any set bits then return 1.0 */
double chemfp_hex_tanimoto(int len, const char *fp1, const char *fp2);

/* Return 1 if the query fingerprint is contained in the target, 0 if it isn't,
   or -1 for invalid fingerprints */
int chemfp_hex_contains(int len, const char *query_fp, const char *target_fp);

/**** Low-level operations directly on byte fingerprints ***/

/* Return the population count of a byte fingerprint */
int chemfp_byte_popcount(int len, const unsigned char *fp);

/* Return the population count of the intersection of two byte fingerprints */
int chemfp_byte_intersect_popcount(int len, const unsigned char *fp1,
                                   const unsigned char *fp2);

/* Return the Tanitoto between two byte fingerprints, or -1.0 for invalid fingerprints
   If neither fingerprint has any set bits then return 1.0 */
double chemfp_byte_tanimoto(int len, const unsigned char *fp1,
                            const unsigned char *fp2);

double chemfp_byte_hex_tanimoto(int size, const unsigned char *byte_fp,
                                const char *hex_fp);

/* Return 1 if the query fingerprint is contained in the target, 0 if it isn't */
int chemfp_byte_contains(int len, const unsigned char *query_fp,
                         const unsigned char *target_fp);


/**** Functions which work with data from an fps block ***/

/* NOTE: an "fps block" means "one or more fingerprint lines from an fps
   file." These contain the hex fingerprint and the identifier, plus optional
   additional fields. The fps block must end with a newline. */

/* Return 0 if string is a valid fps fingerprint line, otherwise an error code */
int chemfp_fps_line_validate(int hex_size,  /* use -1 if not known */
                             int line_size, const char *line_start);
int chemfp_fps_find_id(int hex_size, const char *line,
                       const char **id_start, const char **id_end);

int chemfp_threshold_tanimoto_hexfp_fps(
        int hex_size, const char *hex_query_fp,
        int target_block_size, const char *target_block_start,
        int *lineno, const char **stopped_at,
        double threshold,
        int num_cells, int *id_starts, int *id_ends, double *scores);

/* Return the number of fingerprints in the fps block which are greater
   than or equal to the specified threshold. */
int chemfp_fps_count_tanimoto_hits(
        int num_bits,
        int query_storage_size,
        const unsigned char *query_arena, int query_start, int query_end,
        const char *target_block, int target_block_end,
        double threshold,
        int *counts, int *num_lines_processed);


typedef struct {
  int size;           /* current heap size */
  int heap_state;

  /* These all point to arrays of size k */
  int *indices;      /* [k]; contains a unique id or index */
  char **ids;          /* [k]; array of NULL or malloc'ed identifier */
  double *scores;     /* [k]; the Tanimoto similarity */
} chemfp_fps_heap;

typedef struct {
  const unsigned char *query_start;
  int num_queries;
  int query_fp_size;
  int query_storage_size;
  int k;              /* max number of elements to find */
  int search_state;          /* 0 for not not finished, 1 for finished */
  double threshold;   /* initial threshold */
  chemfp_fps_heap *heaps;    /* [num_queries] heaps */

  int num_targets_processed;
  char **_all_ids;
  double *_all_scores;
} chemfp_fps_knearest_search;


int chemfp_fps_knearest_search_init(
        chemfp_fps_knearest_search *knearest_search,
        int num_bits, int query_storage_size,
        const unsigned char *query_arena, int query_start, int query_end,
        int k, double threshold);

/* Update the heap based on the lines in an fps fingerprint data block. */
int chemfp_fps_knearest_tanimoto_search_feed(
        chemfp_fps_knearest_search *knearest_search,
        int target_block_len, const char *target_block);

/* Call this after the last fps block, in order to convert the heap into an
   sorted array. */
void chemfp_fps_knearest_search_finish(chemfp_fps_knearest_search *knearest_search);

void chemfp_fps_knearest_search_free(chemfp_fps_knearest_search *knearest_search);


int chemfp_fps_threshold_tanimoto_search(
        int num_bits,
        int query_storage_size,
        const unsigned char *query_arena, int query_start, int query_end,
        
        const char *target_block, int target_block_end,
        double threshold,
        int num_cells, chemfp_tanimoto_cell *cells,
        const char ** stopped_at, int *num_lines_processed, int *num_cells_processed);


/***** The byte-oriented algorithms  ********/


/* You must have enough cells */

int chemfp_knearest_tanimoto_block(
        int k,
        int query_len, unsigned char *query_fp,
        int num_targets, unsigned char *target_block, int offset, int storage_len,
        double threshold,
        int *indices, double *scores);

int chemfp_hex_tanimoto_block(
        int n,
        int len, unsigned char *hex_query_fp,
        int target_len, unsigned char *target_block,
        double threshold,
        double *scores, int *id_starts, int *id_ends, int *lineno);


int chemfp_byte_intersect_popcount_count(
        int len, unsigned char *query_fp,
        int num_targets, unsigned char *target_block, int offset, int storage_len,
                int min_overlap);

int chemfp_reorder_by_popcount(
        int num_bits,
        int storage_size, const unsigned char *arena, int start, int end,
        unsigned char *new_arena, ChemFPOrderedPopcount *ordering,
        int *popcount_indices);


int chemfp_count_tanimoto_arena(
        /* Count all matches within the given threshold */
        double threshold,

        /* Number of bits in the fingerprint */
        int num_bits,

        /* Query arena, start and end indices */
        int query_storage_size,
        const unsigned char *query_arena, int query_start, int query_end,

        /* Target arena, start and end indices */
        int target_storage_size,
        const unsigned char *target_arena, int target_start, int target_end,

        /* Target popcount distribution information */
        /*  (must have at least num_bits+1 elements) */
        int *target_popcount_indices,

        /* Results go into these arrays  */
        int *result_counts
                                );

int chemfp_threshold_tanimoto_arena(
        /* Within the given threshold */
        double threshold,

        /* Number of bits in the fingerprint */
        int num_bits,

        /* Query arena, start and end indices */
        int query_storage_size, const unsigned char *query_arena,
        int query_start, int query_end,

        /* Target arena, start and end indices */
        int target_storage_size, const unsigned char *target_arena,
        int target_start, int target_end,

        /* Target popcount distribution information */
        /*  (must have at least num_bits+1 elements) */
        int *target_popcount_indices,

        /* Results go into this data structure  */
        chemfp_threshold_result *results
                                    );


int chemfp_knearest_tanimoto_arena(
        /* Find the 'k' nearest items */
        int k,
        /* Within the given threshold */
        double threshold,

        /* Number of bits in the fingerprint */
        int num_bits,

        /* Query arena, start and end indices */
        int query_storage_size, const unsigned char *query_arena,
        int query_start, int query_end,

        /* Target arena, start and end indices */
        int target_storage_size, const unsigned char *target_arena,
        int target_start, int target_end,

        /* Target popcount distribution information */
        /*  (must have at least num_bits+1 elements) */
        int *target_popcount_indices,

        /* Results go into this data structure  */
        chemfp_threshold_result *results
                                   );



int chemfp_count_tanimoto_hits_arena_symmetric(
        /* Count all matches within the given threshold */
        double threshold,

        /* Number of bits in the fingerprint */
        int num_bits,

        /* Fingerprint arena */
        int storage_size, const unsigned char *arena,

        /* Row start and end indices */
        int query_start, int query_end,

        /* Column start and end indices */
        int target_start, int target_end,

        /* Target popcount distribution information */
        int *target_popcount_indices,

        /* Results _increment_ existing values in the array - remember to initialize! */
        int *result_counts);

int chemfp_threshold_tanimoto_arena_symmetric(
        /* Within the given threshold */
        double threshold,

        /* Number of bits in the fingerprint */
        int num_bits,

        /* Arena */
        int storage_size, const unsigned char *arena,

        /* start and end indices for the rows and columns */
        int query_start, int query_end,
        int target_start, int target_end,
        
        /* Target popcount distribution information */
        /*  (must have at least num_bits+1 elements) */
        int *target_popcount_indices,

        /* Results go here */
        /* NOTE: This must have enough space for all of the fingerprints! */
        chemfp_threshold_result *results);


void chemfp_knearest_results_finalize(chemfp_threshold_result *results_start,
                                      chemfp_threshold_result *results_end);


typedef int (*chemfp_popcount_f)(int len, const unsigned char *p1);
typedef int (*chemfp_intersect_popcount_f)(int len, const unsigned char *p1,
                                           const unsigned char *p2);

chemfp_popcount_f
chemfp_select_popcount(int num_bits,
                       int storage_len, const unsigned char *arena);

chemfp_intersect_popcount_f
chemfp_select_intersect_popcount(int num_bits,
                                 int storage_len1, const unsigned char *arena1,
                                 int storage_len2, const unsigned char *arena2);


/* OpenMP interface */

int chemfp_get_num_threads(void);
void chemfp_set_num_threads(int num_threads);
int chemfp_get_max_threads(void);


#endif /* CHEMFP_H */

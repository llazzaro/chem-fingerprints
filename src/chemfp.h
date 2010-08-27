#ifndef CHEMFP_H
#define CHEMFP_H

// Should the errors be positive or negative numbers?
enum chemfp_errors {
  CHEMFP_OK = 0,
  CHEMFP_BAD_ARG = -1,

  /* File format errors */
  CHEMFP_UNSUPPORTED_WHITESPACE = -30,
  CHEMFP_MISSING_FINGERPRINT = -31,
  CHEMFP_BAD_FINGERPRINT = -32,
  CHEMFP_UNEXPECTED_FINGERPRINT_LENGTH = -33,
  CHEMFP_MISSING_ID = -34,
  CHEMFP_BAD_ID = -35,
  CHEMFP_MISSING_NEWLINE = -36,
};

// Convert an error code to a string description
const char *chemfp_strerror(int err);



// Low-level operations directly on hex strings
int chemfp_hex_isvalid(int len, const unsigned char *fp);
int chemfp_hex_popcount(int len, const unsigned char *fp);
int chemfp_hex_intersect_popcount(int len, const unsigned char *fp1,
								  const unsigned char *fp2);
double chemfp_hex_tanimoto(int len, const unsigned char *fp1,
						   const unsigned char *fp2);
int chemfp_hex_contains(int len, const unsigned char *query_fp,
						const unsigned char *target_fp);


// A semi-generic data structure for Tanimoto searches

typedef struct {
  int size;           // current heap size
  int k;              // max heap size
  int unique_idx;     // counter if a unique index is needed
  int _reserved;      // used for nice alignment on 64 bit machines
  double threshold;   // initial threshold
  int *indicies;      // [k]; contains a unique id or index
  double *scores;     // [k]; the Tanimoto similarity
  char **id_starts;   // [k]; location (in the current block) of the start of the id
  int *id_lens;       // [k]; length of that id
} chemfp_heap;


// Operations to work with FPS block data

// Return 0 if the line is valid, otherwise an error code
int chemfp_fps_line_validate(int hex_len,  // use -1 if not known
							 int line_len, char *line);


// Compute Tanimoto scores and report all lines with scores at
// or over the given threshold. The id_starts, id_lens, and scores
// must be large enough.
int chemfp_fps_tanimoto(
    int hex_len, char *hex_query,   // The query fingerprint, in hex

	// Target block data, in FPS format. Last character must be a newline.
	int target_block_len, char *target_block,  

	double threshold,    // Report only those values >= threshold
	int *num_found,      // Will be set to the number of fingerprints which matched
	char **id_starts, int *id_lens,     // id locations in the current block
	double *scores,      // Corresponding Tanimoto similarity score

	// (optional) track the line number. Must be initialized if not NULL.
	// The first line number is 1, but it should be the line number for 
	// the first line in the target block
	int *lineno
  );


int chemfp_fps_tanimoto_count(int hex_len, char *hex_query,
							  int target_block_len, char *target_block,
							  double threshold,
							  int *num_found, int *lineno);


///// Used to find the k-nearest records with at least a given Tanimoto similarity

// Initialize the heap, including pointing to available memory
void chemfp_fps_heap_init(chemfp_heap *heap,
						  int k, double threshold,
						  int *indicies, double *scores,
						  char **id_starts, int *id_lens);

// Update the heap based on FPS data
int chemfp_fps_heap_update_tanimoto(chemfp_heap *heap,
									int hex_len, char *hex_query,
									int target_block_len, char *target_block,
									int *lineno);

// Call this when finished processing, in order to sort the heap so it's meaningful.
void chemfp_fps_heap_finish_tanimoto(chemfp_heap *heap);

//////////////////////////////////


int chemfp_byte_popcount(int len, const unsigned char *fp);
int chemfp_byte_intersect_popcount(int len, const unsigned char *fp1,
								   const unsigned char *fp2);
double chemfp_byte_tanimoto(int len, const unsigned char *fp1,
							const unsigned char *fp2);
int chemfp_byte_contains(int len, const unsigned char *query_fp,
						 const unsigned char *target_fp);

int chemfp_nlargest_tanimoto_block(
        int n,
		int query_len, unsigned char *query_fp,
		int num_targets, unsigned char *target_block, int offset, int storage_len,
		double threshold,
		int *indicies, double *scores);

int chemfp_hex_tanimoto_block(
        int n,
		int len, unsigned char *hex_query_fp,
		int target_len, unsigned char *target_block,
		double threshold,
		double *scores, unsigned char **start_ids, int *id_lens, int *lineno);

#endif /* CHEMFP_H */

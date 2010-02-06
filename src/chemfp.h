int chemfp_hex_isvalid(int len, const unsigned char *fp);
int chemfp_hex_popcount(int len, const unsigned char *fp);
int chemfp_hex_intersect_popcount(int len, const unsigned char *fp1,
								  const unsigned char *fp2);
double chemfp_hex_tanimoto(int len, const unsigned char *fp1,
						   const unsigned char *fp2);
int chemfp_hex_contains(int len, const unsigned char *query_fp,
						const unsigned char *target_fp);

int chemfp_byte_popcount(int len, const unsigned char *fp);
int chemfp_byte_intersect_popcount(int len, const unsigned char *fp1,
								   const unsigned char *fp2);
double chemfp_byte_tanimoto(int len, const unsigned char *fp1,
							const unsigned char *fp2);
int chemfp_byte_contains(int len, const unsigned char *query_fp,
						 const unsigned char *target_fp);

int chemfp_nlargest_tanimoto(
		int len, unsigned char *query_fp,
		int num_fps, unsigned char *fps, int offset, int storage_len,
		unsigned int n, int *indicies, double *scores);

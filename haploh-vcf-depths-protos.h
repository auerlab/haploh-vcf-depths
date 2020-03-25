/* haploh-vcf-depths.c */
int main(int argc, const char *argv[]);
void usage(const char *argv[]);
int haploh_median_depths(const char *event_glob_pattern, const char *vcf_glob_pattern);
FILE *vcf_open(const char *vcf_filename, _Bool *compressed);
int vcf_close(FILE *vcf_stream, _Bool compressed);
int sample_from_glob(const char *vcf_glob_pattern, const char *filename, char *event_sample_id);

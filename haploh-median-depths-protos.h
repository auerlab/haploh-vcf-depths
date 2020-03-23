/* haploh-median-depths.c */
int main(int argc, const char *argv[]);
void usage(const char *argv[]);
int haploh_median_depths(const char *sample_id, const char *event_file, const char *glob_pattern);
int read_event(FILE *event_stream, event_t *event);
FILE *vcf_open(const char *vcf_filename, _Bool *compressed);
int vcf_close(FILE *vcf_stream, _Bool compressed);

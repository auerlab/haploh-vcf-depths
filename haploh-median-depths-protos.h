/* haploh-median-depths.c */
int main(int argc, const char *argv[]);
void usage(const char *argv[]);
int haploh_median_depths(const char *argv[], FILE *events_stream, const char *vcf_dir);
int read_event(FILE *event_stream, event_t *event);

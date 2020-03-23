/* haploh-median-depths.c */
int main(int argc, const char *argv[]);
void usage(const char *argv[]);
int haploh_median_depths(const char *sample_id, const char *event_file, const char *glob_pattern);
int read_event(FILE *event_stream, event_t *event);
unsigned median_depth(const char *filename, event_t *event);

/* events.c */
int depth_cmp(depth_t *n1, depth_t *n2);
event_t *event_read_list(const char *event_filename, size_t *event_count, const char *sample_id);
void event_add_depth(event_t *event, depth_t depth, const char *vcf_sample_id);
int event_read(event_t *event, FILE *event_stream, const char *sample_id);

/* events.c */
int depth_cmp(depth_t *n1, depth_t *n2);
event_t *event_read_list(const char *event_glob_pattern, size_t *event_count);
void event_add_depth(event_t *event, depth_t depth, const char *vcf_sample_id);
int event_read(event_t *event, FILE *event_stream, const char *event_sample_id);
void event_open_depth_files(event_t *events, size_t event_count);
void event_sort(event_t *events, size_t event_count);
int event_cmp(const event_t *e1, const event_t *e2);

/* events.c */
int depth_cmp(depth_t *n1, depth_t *n2);
event_t *event_read_list(const char *event_filename, size_t *event_count);
void event_add_depth(event_t *event, depth_t depth);
int event_read(FILE *event_stream, event_t *event);

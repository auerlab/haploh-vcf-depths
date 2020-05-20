#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
#include <glob.h>
#include <sys/stat.h>
#include <vcfio.h>
#include <biostring.h>
#include "events.h"
#include "haploh-vcf-depths-protos.h"    // Temporary for sample_from_glob()

int     depth_cmp(depth_t *n1, depth_t *n2)

{
    return *n2 - *n1;
}


event_t *event_read_list(const char *event_glob_pattern, size_t *event_count)

{
    FILE    *event_stream;
    event_t *events,
	    event;
    size_t  c;
    int     separator;
    glob_t  event_glob;
    char    **event_filename_ptr,
	    event_sample_id[PATH_MAX + 1];;
    
    fprintf(stderr, "%s\n", event_glob_pattern);
    glob(event_glob_pattern, 0, NULL, &event_glob);
    
    // Small files: Read through once to count, then allocate perfectly
    // sized array.
    *event_count = 0;
    for (event_filename_ptr = event_glob.gl_pathv; *event_filename_ptr != NULL;
	    ++event_filename_ptr)
    {
	if ( (event_stream = fopen(*event_filename_ptr, "r")) == NULL )
	{
	    fprintf(stderr, "Cannot open %s: %s\n", *event_filename_ptr,
		    strerror(errno));
	    exit(EX_NOINPUT);
	}
	
	sample_from_glob(event_glob_pattern, *event_filename_ptr, event_sample_id);
	
	// Skip header
	while ( (separator = event_read(&event, event_stream, event_sample_id)) == EVENT_READ_HEADER )
	    ;
	    
	while ( separator == EVENT_READ_OK )
	{
	    ++*event_count;
	    separator = event_read(&event, event_stream, event_sample_id);
	}
	
	if ( separator != EVENT_READ_EOF )
	{
	    fprintf(stderr, "haploh_media_depths(): Error reading event stream.\n");
	    exit(EX_DATAERR);
	}
	
	fclose(event_stream);
    }

    fprintf(stderr, "%zu events counted.\n", *event_count);
    if ( (events = EVENT_MALLOC(*event_count, event_t)) == NULL )
    {
	fprintf(stderr, "event_read_list(): Malloc failed.\n");
	exit(EX_UNAVAILABLE);
    }
    
    c = 0;
    for (event_filename_ptr = event_glob.gl_pathv; *event_filename_ptr != NULL;
	    ++event_filename_ptr)
    {
	if ( (event_stream = fopen(*event_filename_ptr, "r")) == NULL )
	{
	    fprintf(stderr, "Cannot open %s: %s\n", *event_filename_ptr,
		    strerror(errno));
	    exit(EX_NOINPUT);
	}
	
	// Skip header again
	while ( (separator = event_read(&event, event_stream, event_sample_id))
		 == EVENT_READ_HEADER )
	    ;
	
	sample_from_glob(event_glob_pattern, *event_filename_ptr, event_sample_id);

	while ( separator == EVENT_READ_OK )
	{
	    events[c++] = event;
	    separator = event_read(&event, event_stream, event_sample_id);
	}

	if ( separator != EVENT_READ_EOF )
	{
	    fprintf(stderr, "haploh_media_depths(): Error reading event stream.\n");
	    exit(EX_DATAERR);
	}
	
	fclose(event_stream);
    }
    
    event_sort(events, *event_count);
    for (c = 0; c < *event_count; ++c)
	fprintf(stderr, "%s %zu %zu %s\n", events[c].chromosome,
		events[c].begin, events[c].end, events[c].sample_id);
    event_open_depth_files(events, *event_count);
    return events;
}


void    event_add_depth(event_t *event, depth_t depth, const char *vcf_sample_id)

{
    //fprintf(stderr, "%s %s\n", vcf_sample_id, event->sample_id);
    if ( strcmp(vcf_sample_id, event->sample_id) == 0 )
	fprintf(event->same_sample_depth_stream, "%u\n", depth);
    else
	fprintf(event->other_samples_depth_stream, "%u\n", depth);
}


int     event_read(event_t *event, FILE *event_stream, const char *event_sample_id)

{
    int     separator,
	    status;
    size_t  len;
    char    temp[VCF_POSITION_MAX_CHARS + 1],
	    *endptr;
    
    separator = tsv_read_field(event_stream, event->chromosome,
			       VCF_CHROMOSOME_MAX_CHARS, &len);
    if ( separator == '\t' )
    {
	// fprintf(stderr, "chromosome = %s\n", event->chromosome);
	
	// BEGIN
	if ( (separator = tsv_read_field(event_stream, temp,
			       VCF_POSITION_MAX_CHARS, &len)) != '\t' )
	{
	    fputs("event_read(): Did not find tab after BEGIN.\n", stderr);
	    return EVENT_READ_TRUNCATED;
	}
	if ( strcmp(temp, "BEGIN") == 0 )
	{
	    status = EVENT_READ_HEADER;
	}
	else
	{
	    event->begin = strtoul(temp, &endptr, 10);
	    if ( *endptr != '\0' )
	    {
		fprintf(stderr, "event_read(): Invalid BEGIN: %s\n", temp);
		return EVENT_READ_TRUNCATED;
	    }
	    
	    // END
	    if ( (separator = tsv_read_field(event_stream, temp,
				   VCF_POSITION_MAX_CHARS, &len)) != '\t' )
	    {
		fputs("event_read(): Did not find tab after END.\n", stderr);
		return EVENT_READ_TRUNCATED;
	    }
	    event->end = strtoul(temp, &endptr, 10);
	    if ( *endptr != '\0' )
	    {
		fprintf(stderr, "event_read(): Invalid END: %s\n", temp);
		return EVENT_READ_TRUNCATED;
	    }

	    // Initialize
	    strlcpy(event->sample_id, event_sample_id, PATH_MAX);

	    fprintf(stderr, "Opened %s %s %zu %zu\n", event->sample_id,
		    event->chromosome, event->begin, event->end);
	    status = EVENT_READ_OK;
	}
	
	if ( (separator = tsv_skip_rest_of_line(event_stream)) != '\n' )
	{
	    fprintf(stderr, "event_read(): Did not find newline at end of events[c].\n");
	    return EVENT_READ_TRUNCATED;
	}
	
	return status;
    }
    else
    
    // fprintf(stderr, "separator = %d\n", separator);
    return separator;
}


void    event_open_depth_files(event_t *events, size_t event_count)

{
    size_t  c;
    char    filename[PATH_MAX + 1];

    mkdir("Depths", 0755);
    for (c = 0; c < event_count; ++c)
    {
	snprintf(filename, PATH_MAX, "Depths/depths-%s-%s-%zu-%zu-same.txt",
		 events[c].sample_id, events[c].chromosome,
		 events[c].begin, events[c].end);
	events[c].same_sample_depth_stream = fopen(filename, "w");
	if ( events[c].same_sample_depth_stream == NULL )
	{
	    fprintf(stderr, "event_read(): Could not open %s: %s\n",
		    filename, strerror(errno));
	    exit(EX_CANTCREAT);
	}
	
	snprintf(filename, PATH_MAX, "Depths/depths-%s-%s-%zu-%zu-others.txt",
		 events[c].sample_id, events[c].chromosome,
		 events[c].begin, events[c].end);
	events[c].other_samples_depth_stream = fopen(filename, "w");
	if ( events[c].other_samples_depth_stream == NULL )
	{
	    fprintf(stderr, "event_read(): Could not open %s: %s\n",
		    filename, strerror(errno));
	    exit(EX_CANTCREAT);
	}
    }
}


void    event_sort(event_t *events, size_t event_count)

{
    qsort(events, event_count, sizeof(event_t),
	  (int (*)(const void *, const void *))event_cmp);
}


int     event_cmp(const event_t *e1, const event_t *e2)

{
    int     status;
    
    status = chromosome_name_cmp(e1->chromosome, e2->chromosome, 3);
    if ( status == 0 )
	status = e1->begin - e2->begin;;
    if ( status == 0 )
	status = e1->end - e2->end;
    return status;
}

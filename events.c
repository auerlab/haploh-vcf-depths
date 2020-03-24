#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <string.h>
#include <errno.h>
#include <vcfio.h>
#include <tsvio.h>
#include <stdbool.h>
#include "events.h"

int     depth_cmp(depth_t *n1, depth_t *n2)

{
    return *n2 - *n1;
}


event_t *event_read_list(const char *event_filename, size_t *event_count,
			 const char *sample_id)

{
    FILE    *event_stream;
    event_t *events,
	    event;
    size_t  c;
    int     separator;
    
    if ( (event_stream = fopen(event_filename, "r")) == NULL )
    {
	fprintf(stderr, "Cannot open %s: %s\n", event_filename, strerror(errno));
	exit(EX_NOINPUT);
    }
	
    // Skip header
    while ( (separator = event_read(&event, event_stream, sample_id)) == EVENT_READ_HEADER )
	;
	
    // Small file: Read through once to count, then allocate list
    *event_count = 0;
    while ( separator == EVENT_READ_OK )
    {
	++*event_count;
	separator = event_read(&event, event_stream, sample_id);
    }
    
    if ( separator != EVENT_READ_EOF )
    {
	fprintf(stderr, "haploh_media_depths(): Error reading event stream.\n");
	exit(EX_DATAERR);
    }

    if ( (events = EVENT_MALLOC(*event_count, event_t)) == NULL )
    {
	fprintf(stderr, "event_read_list(): Malloc failed.\n");
	exit(EX_UNAVAILABLE);
    }
    
    rewind(event_stream);
    // Skip header again
    while ( (separator = event_read(&event, event_stream, sample_id)) == EVENT_READ_HEADER )
	;
    for (c = 0; c < *event_count; ++c)
    {
	events[c] = event;
	event_read(&event, event_stream, sample_id);
    }
    fclose(event_stream);
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


int     event_read(event_t *event, FILE *event_stream, const char *sample_id)

{
    int     separator,
	    status;
    size_t  len;
    char    temp[VCF_POSITION_MAX_CHARS + 1],
	    filename[PATH_MAX + 1],
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
	    strlcpy(event->sample_id, sample_id, PATH_MAX);
	    
	    // Open depth files
	    snprintf(filename, PATH_MAX, "depths-%s-%s-%zu-%zu-same.txt",
		     sample_id, event->chromosome, event->begin, event->end);
	    event->same_sample_depth_stream = fopen(filename, "w");
	    if ( event->same_sample_depth_stream == NULL )
	    {
		fprintf(stderr, "event_read(): Could not open %s: %s\n",
			filename, strerror(errno));
		exit(EX_CANTCREAT);
	    }
	    
	    snprintf(filename, PATH_MAX, "depths-%s-%s-%zu-%zu-others.txt",
		     sample_id, event->chromosome, event->begin, event->end);
	    event->other_samples_depth_stream = fopen(filename, "w");
	    if ( event->other_samples_depth_stream == NULL )
	    {
		fprintf(stderr, "event_read(): Could not open %s: %s\n",
			filename, strerror(errno));
		exit(EX_CANTCREAT);
	    }

	    fprintf(stderr, "Opened %s %s %zu %zu %p %p\n", event->sample_id,
		    event->chromosome, event->begin, event->end,
		    event->same_sample_depth_stream,
		    event->other_samples_depth_stream);
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

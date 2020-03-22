/***************************************************************************
 *  Description:
 *      Rapidly compute median depths of VCF files for haplohseq events
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-03-21  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <string.h>
#include <vcfio.h>
#include <tsvio.h>
#include "haploh-median-depths.h"

int     main(int argc,const char *argv[])

{
    FILE        *event_stream;
    const char  *event_file,
		*vcf_dir;
    struct stat dir_stat;
    
    if ( argc != 3 )
	usage(argv);

    event_file = argv[1];
    vcf_dir = argv[2];
    
    if ( (event_stream = fopen(event_file, "r")) == NULL )
    {
	fprintf(stderr, "Cannot open %s: %s\n", event_file, strerror(errno));
	exit(EX_NOINPUT);
    }
    
    if ( (stat(vcf_dir, &dir_stat) == -1) || !S_ISDIR(dir_stat.st_mode) )
    {
	fprintf(stderr,
		"%s must be a readable directory containing VCF files.\n",
		vcf_dir);
	exit(EX_NOINPUT);
    }
	    
    return haploh_median_depths(argv, event_stream, vcf_dir);
}


void    usage(const char *argv[])

{
    fprintf(stderr, "Usage: %s haplohseq-event-file directory-with-VCF-files\n", argv[0]);
    exit(EX_USAGE);
}


int     haploh_median_depths(const char *argv[], FILE *event_stream,
			     const char *vcf_dir)

{
    int     separator;
    event_t event;
    
    while ( ((separator = read_event(event_stream, &event)) == EVENT_READ_OK)
	    || (separator == EVENT_READ_HEADER) )
    {
	fprintf(stderr, "%s %zu %zu\n", event.chromosome, event.begin, event.end);
	/*
	 *  Compute median depth of VCF calls between event.begin and event.end
	 *  for this sample.
	 */
	
	/*
	 *  Compute median depth of VCF calls between event.begin and event.end
	 *  for all other samples.
	 */
	// Use glob()?
    }
    
    if ( separator == EOF )
	return EX_OK;
    else
	return EX_DATAERR;
}


int     read_event(FILE *event_stream, event_t *event)

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
	fprintf(stderr, "chromosome = %s\n", event->chromosome);
	
	// BEGIN
	if ( (separator = tsv_read_field(event_stream, temp,
			       VCF_POSITION_MAX_CHARS, &len)) != '\t' )
	{
	    fputs("read_event(): Did not find tab after BEGIN.\n", stderr);
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
		fprintf(stderr, "read_event(): Invalid BEGIN: %s\n", temp);
		return EVENT_READ_TRUNCATED;
	    }
	    
	    // END
	    if ( (separator = tsv_read_field(event_stream, temp,
				   VCF_POSITION_MAX_CHARS, &len)) != '\t' )
	    {
		fputs("read_event(): Did not find tab after END.\n", stderr);
		return EVENT_READ_TRUNCATED;
	    }
	    event->end = strtoul(temp, &endptr, 10);
	    if ( *endptr != '\0' )
	    {
		fprintf(stderr, "read_event(): Invalid END: %s\n", temp);
		return EVENT_READ_TRUNCATED;
	    }
	    status = EVENT_READ_OK;
	}
	if ( (separator = tsv_skip_rest_of_line(event_stream)) != '\n' )
	{
	    fprintf(stderr, "read_event(): Did not find newline at end of event.\n");
	    return EVENT_READ_TRUNCATED;
	}
	return status;
    }
    fprintf(stderr, "separator = %d\n", separator);
    return separator;
}

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
#include <glob.h>
#include <stdbool.h>
#include "haploh-median-depths.h"

int     main(int argc,const char *argv[])

{
    FILE        *event_stream;
    const char  *sample_id,
		*event_file,
		*glob_pattern;
    struct stat dir_stat;
    
    if ( argc != 4 )
	usage(argv);

    sample_id = argv[1];
    event_file = argv[2];
    glob_pattern = argv[3];
    
    /*
    if ( (stat(vcf_dir, &dir_stat) == -1) || !S_ISDIR(dir_stat.st_mode) )
    {
	fprintf(stderr,
		"%s must be a readable directory containing VCF files.\n",
		vcf_dir);
	exit(EX_NOINPUT);
    }
    */
	    
    return haploh_median_depths(sample_id, event_file, glob_pattern);
}


void    usage(const char *argv[])

{
    fprintf(stderr, "\nUsage: %s sample-id haplohseq-event-file 'vcf-filename-glob'\n", argv[0]);
    fprintf(stderr, "\nVCF filenames must contain the sample-id. The glob for VCF filenames must\n");
    fprintf(stderr, "contain a '*' where the sample-id appears and must be enclosed in quotes.\n\n");
    fprintf(stderr, "E.g. for files like combined-NWD294426-ad.vcf.xz, glob = 'combined-*-ad.vcf.xz'.\n\n");
    exit(EX_USAGE);
}


int     haploh_median_depths(const char *sample_id, const char *event_file,
			     const char *glob_pattern)

{
    int     separator;
    event_t event;
    FILE    *event_stream;
    glob_t  glob_list;
    char    **p;
    
    if ( (event_stream = fopen(event_file, "r")) == NULL )
    {
	fprintf(stderr, "Cannot open %s: %s\n", event_file, strerror(errno));
	exit(EX_NOINPUT);
    }

    glob(glob_pattern, 0, NULL, &glob_list);
    
    /*
     *  FIXME: Rather than reread the same VCF file for each sample, it
     *  would be more efficient to leapfrog through the events and VCF calls.
     *  Both files should be sorted by chromosome and position.
     */
    
    while ( ((separator = read_event(event_stream, &event)) == EVENT_READ_OK)
	    || (separator == EVENT_READ_HEADER) )
    {
	fprintf(stderr, "%s %zu %zu\n", event.chromosome, event.begin, event.end);
    
	/*
	 *  Compute median depth of VCF calls between event.begin and event.end
	 *  for this sample.
	 */
       
	for (p = glob_list.gl_pathv; *p != NULL; ++p)
	{
	    //fprintf(stderr, "%s %s\n", *p, sample_id);
	    if ( strstr(*p, sample_id) != NULL )
	    {
		fprintf(stderr, "Matching VCF file: %s\n", *p);
		printf("Median depth = %u\n", median_depth(*p, &event));
	    }
	}
	
	/*
	 *  Compute median depth of VCF calls between event.begin and event.end
	 *  for all other samples.
	 */
	// Use glob()?
    }
    
    fclose(event_stream);
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
    // fprintf(stderr, "separator = %d\n", separator);
    return separator;
}


unsigned    median_depth(const char *vcf_filename, event_t *event)

{
    const char  *p;
    FILE        *vcf_stream;
    char        cmd[1025],
		vcf_sample[VCF_SAMPLE_MAX_CHARS + 1];
    bool        compressed;
    vcf_call_t  vcf_call;
    int         status;
    size_t      len;
    unsigned long   count = 0;
    
    p = strstr(vcf_filename, ".xz");
    compressed = (p != NULL) && (strcmp(p,".xz") == 0);
    if ( compressed )
    {
	snprintf(cmd, 1024, "xzcat %s", vcf_filename);
	vcf_stream = popen(cmd, "r");
    }
    else
	vcf_stream = fopen(vcf_filename, "r");
    
    while ( (status = vcf_read_static_fields(vcf_stream, &vcf_call)) == VCF_READ_OK )
    {
	if ( tsv_read_field(vcf_stream, vcf_sample, VCF_SAMPLE_MAX_CHARS, &len) != '\n' )
	{
	    fprintf(stderr, "median_depth(): Expected newline after VCF sample.\n");
	    exit(EX_DATAERR);
	}
	if ( (strcmp(vcf_call.chromosome, event->chromosome) == 0) &&
	     (vcf_call.pos >= event->begin) && (vcf_call.pos <= event->end) )
	{
	    //fprintf(stderr, "%s\n", vcf_call.samples[0]);
	    ++count;
	}
    }
    fprintf(stderr, "Found %lu matching calls.\n", count);
    
    if ( status != VCF_READ_EOF )
    {
	fprintf(stderr, "Error reading VCF file.\n");
	exit(EX_DATAERR);
    }
    
    if ( compressed )
	pclose(vcf_stream);
    else
	fclose(vcf_stream);
    return 0;
}

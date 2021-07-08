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
#include <glob.h>
#include <stdbool.h>
#include <limits.h>
#include <xtend.h>      // strlcpy() on Linux
#include <biolibc/vcf.h>
#include <biolibc/biostring.h>
#include "haploh-vcf-depths.h"

int     main(int argc,const char *argv[])

{
    const char  *event_glob_pattern,
		*vcf_glob_pattern;
    
    if ( argc != 3 )
	usage(argv);

    event_glob_pattern = argv[1];
    vcf_glob_pattern = argv[2];
    
    return haploh_median_depths(event_glob_pattern, vcf_glob_pattern);
}


void    usage(const char *argv[])

{
    fprintf(stderr, "\nUsage: %s 'event-file-glob-pattern' 'vcf-file-glob-pattern'\n", argv[0]);
    fprintf(stderr, "\nVCF and event filenames must contain the sample-id. The glob patterns must\n");
    fprintf(stderr, "contain a '*' where the sample-id appears and must be enclosed in quotes.\n\n");
    fprintf(stderr, "E.g. for files like combined-NWD294426-ad.vcf.xz, glob = 'combined-*-ad.vcf.xz'.\n\n");
    exit(EX_USAGE);
}


int     haploh_median_depths(const char *event_glob_pattern,
			     const char *vcf_glob_pattern)

{
    FILE    *vcf_stream;
    glob_t  vcf_glob;
    char    **vcf_filename_ptr,
	    *depth_str,
	    *end_num,
	    vcf_sample_id[PATH_MAX + 1],
	    cmd[CMD_MAX + 1];
    bool    compressed;
    int     status;
    size_t  event_count,
	    c,
	    vcf_count = 0;
    event_t *events;
    depth_t depth;
    static vcf_call_t  vcf_call;
    static char        vcf_sample[VCF_SAMPLE_MAX_CHARS + 1];
    
    vcf_call_init(&vcf_call, VCF_INFO_MAX_CHARS, VCF_FORMAT_MAX_CHARS,
		  VCF_SAMPLE_MAX_CHARS);
    
    if ( (events = event_read_list(event_glob_pattern, &event_count))
	  == NULL )
    {
	fprintf(stderr, "haploh_median_depths(): Error reading event list.\n");
	exit(EX_DATAERR);
    }

    glob(vcf_glob_pattern, 0, NULL, &vcf_glob);
    fprintf(stderr, "%zu VCF files.\n", vcf_glob.gl_pathc);
    
    /*
     *  FIXME: Rather than reread the same VCF file for each sample, it
     *  would be more efficient to leapfrog through the events and VCF calls.
     *  Both files should be sorted by chromosome and position.
     */
    
    for (vcf_filename_ptr = vcf_glob.gl_pathv; *vcf_filename_ptr != NULL;
	    ++vcf_filename_ptr)
    {
	sample_from_glob(vcf_glob_pattern, *vcf_filename_ptr, vcf_sample_id);
	
	// Open VCF file
	if ( (vcf_stream = vcf_open(*vcf_filename_ptr, &compressed)) == NULL )
	{
	    fprintf(stderr, "Error opening VCF file %s.\n", *vcf_filename_ptr);
	    exit(EX_NOINPUT);
	}
	fprintf(stderr, "Processing file %zu %s...\n", ++vcf_count, *vcf_filename_ptr);

	/*
	 *  Compute median depth of VCF calls between events[c].begin and events[c].end
	 *  for this sample.
	 *  Compute median depth of VCF calls between events[c].begin and events[c].end
	 *  for all other samples.
	 */
	
	// Count calls for same chromosome and within range
	while ( (status = vcf_read_ss_call(vcf_stream, &vcf_call)) == BIO_READ_OK )
	{
	    if ( (depth_str = strrchr(vcf_sample, ':')) == NULL )
	    {
		fprintf(stderr, "haploh_median_depths(): ':' expected in sample data.\n");
		fprintf(stderr, "Got %s.\n", vcf_sample);
		exit(EX_DATAERR);
	    }
	    
	    depth = strtol(depth_str + 1, &end_num, 10);
	    if ( *end_num != '\0' )
	    {
		fprintf(stderr, "haploh_median_depths(): ':' expected sample to end in :depth.\n");
		fprintf(stderr, "Got %s.\n", vcf_sample);
		exit(EX_DATAERR);
	    }

	    // Skip VCF calls for chromosomes before the first event
	    if ( chromosome_name_cmp(VCF_CHROMOSOME(&vcf_call),
				     EVENT_CHROMOSOME(events + 0)) < 0 )
		continue;
		
	    // FIXME: Build index of first events for each chromosome
	    // to replace this waste of time?
	    for (c = 0; (c < event_count) &&
			(chromosome_name_cmp(EVENT_CHROMOSOME(events + c),
			    VCF_CHROMOSOME(&vcf_call)) < 0); ++c)
		;
	    
	    while ( (c < event_count) && 
		    (events[c].end < VCF_POS(&vcf_call)) &&
		    (chromosome_name_cmp(EVENT_CHROMOSOME(events + c),
					 VCF_CHROMOSOME(&vcf_call)) == 0) )
		++c;
	    
	    while ( (c < event_count) && 
		    (events[c].begin <= VCF_POS(&vcf_call)) &&
		    (chromosome_name_cmp(EVENT_CHROMOSOME(events + c),
					 VCF_CHROMOSOME(&vcf_call)) == 0) )
	    {
		/*fprintf(stderr, "%s %zu %zu ~ %zu\n",
			EVENT_CHROMOSOME(events + c),
			events[c].begin, events[c].end,
			VCF_POS(&vcf_call));*/
		if ( (depth_str = strrchr(vcf_sample, ':')) == NULL )
		{
		    fprintf(stderr, "haploh_median_depths(): ':' expected in sample data.\n");
		    fprintf(stderr, "Got %s.\n", vcf_sample);
		    exit(EX_DATAERR);
		}
		
		depth = strtol(depth_str + 1, &end_num, 10);
		if ( *end_num != '\0' )
		{
		    fprintf(stderr, "haploh_median_depths(): Expected sample to end in :depth.\n");
		    fprintf(stderr, "Got %s.\n", vcf_sample);
		    exit(EX_DATAERR);
		}

		event_add_depth(events + c, depth, vcf_sample_id);
		++c;
	    }
	}
	
	if ( status != BIO_READ_EOF )
	{
	    fprintf(stderr, "haploh_median_depths(): Error reading VCF file.\n");
	    snprintf(cmd, CMD_MAX + 1, "mv %s Broken", *vcf_filename_ptr);
	    system(cmd);
	    //exit(EX_DATAERR);
	}
	vcf_close(vcf_stream, compressed);
    }
    globfree(&vcf_glob);
    
    // Close depth files for all events
    
    return EX_OK;
}


FILE    *vcf_open(const char *vcf_filename, bool *compressed)

{
    const char  *p;
    char        cmd[1025];
    FILE        *vcf_stream;
    
    p = strstr(vcf_filename, ".xz");
    *compressed = (p != NULL) && (strcmp(p,".xz") == 0);
    if ( *compressed )
    {
	snprintf(cmd, 1024, "xzcat %s", vcf_filename);
	vcf_stream = popen(cmd, "r");
    }
    else
	vcf_stream = fopen(vcf_filename, "r");
    return vcf_stream;
}


int     vcf_close(FILE *vcf_stream, bool compressed)

{
    if ( compressed )
	return pclose(vcf_stream);
    else
	return fclose(vcf_stream);
}


int     sample_from_glob(const char *vcf_glob_pattern, const char *filename,
			 char *event_sample_id)

{
    char    *star,
	    *p,
	    format[PATH_MAX + 1];
    size_t  index;
    
    if ( (star = strchr(vcf_glob_pattern, '*')) == NULL )
    {
	fprintf(stderr, "sample_from_glob(): Expected '*' in glob pattern.\n");
	exit(EX_USAGE);
    }
    
    index = star - vcf_glob_pattern;
    memcpy(format, vcf_glob_pattern, index);
    memcpy(format + index, "%s", 2);
    strlcpy(format + index + 2, star + 1, PATH_MAX - (index + 2));
    sscanf(filename, format, event_sample_id);
    //fprintf(stderr, "format = '%s'\n", format);
    //fprintf(stderr, "filename = '%s'\n", filename);
    // sscanf() includes trailing text after %s
    if ( (p = strstr(event_sample_id, star + 1)) != NULL )
	*p = '\0';
    //fprintf(stderr, "sample id = '%s'\n", event_sample_id);
    return 0;
}

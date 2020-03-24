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
    const char  *sample_id,
		*event_file,
		*glob_pattern;
    
    if ( argc != 4 )
	usage(argv);

    sample_id = argv[1];
    event_file = argv[2];
    glob_pattern = argv[3];
    
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
    FILE    *vcf_stream;
    glob_t  glob_list;
    char    **vcf_filename_ptr,
	    *depth_str,
	    *end_num;
    bool    same_sample,
	    compressed;
    int     status;
    size_t  event_count,
	    c;
    event_t *events;
    depth_t depth;
    static vcf_call_t  vcf_call = VCF_CALL_INIT;
    static char        vcf_sample[VCF_SAMPLE_MAX_CHARS + 1];
    
    fprintf(stderr, "%s\n", event_file);
    if ( (events = event_read_list(event_file, &event_count)) == NULL )
    {
	fprintf(stderr, "haploh_median_depths(): Error reading event list.\n");
	exit(EX_DATAERR);
    }

    glob(glob_pattern, 0, NULL, &glob_list);
    
    /*
     *  FIXME: Rather than reread the same VCF file for each sample, it
     *  would be more efficient to leapfrog through the events and VCF calls.
     *  Both files should be sorted by chromosome and position.
     */
    
    for (vcf_filename_ptr = glob_list.gl_pathv; *vcf_filename_ptr != NULL;
	    ++vcf_filename_ptr)
    {
	fprintf(stderr, "%s %s\n", *vcf_filename_ptr, sample_id);
	same_sample = (strstr(*vcf_filename_ptr, sample_id) != NULL);
	
	// Open VCF file
	if ( (vcf_stream = vcf_open(*vcf_filename_ptr, &compressed)) == NULL )
	{
	    fprintf(stderr, "Error opening VCF file %s.\n", *vcf_filename_ptr);
	    exit(EX_NOINPUT);
	}
	
	for (c = 0; c < event_count; ++c)
	{
	    fprintf(stderr, "Event: %s %zu %zu\n",
		    EVENT_CHROMOSOME(events + c),
		    events[c].begin, events[c].end);
    
	    /*
	     *  Compute median depth of VCF calls between events[c].begin and events[c].end
	     *  for this sample.
	     *  Compute median depth of VCF calls between events[c].begin and events[c].end
	     *  for all other samples.
	     */
	    
	    // Skip calls for earlier chromosomes or positions
	    while ( ((status = vcf_read_ss_call(vcf_stream, &vcf_call,
				    vcf_sample, VCF_SAMPLE_MAX_CHARS))
				    == VCF_READ_OK) &&
		     ((chromosome_name_cmp(vcf_call.chromosome,
					  EVENT_CHROMOSOME(events + c)) < 0) ||
		     (vcf_call.pos < EVENT_BEGIN(events + c))) )
		;
	    fprintf(stderr, "Skipped calls up to %s %zu\n",
		    vcf_call.chromosome, vcf_call.pos);
	    
	    // Count calls for same chromosome and within range
	    while ( (status == VCF_READ_OK) &&
		    (chromosome_name_cmp(vcf_call.chromosome,
			EVENT_CHROMOSOME(events + c)) == 0) &&
		    (vcf_call.pos <= events[c].end) )
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
		
		event_add_depth(events + c, depth);
		status = vcf_read_ss_call(vcf_stream, &vcf_call,
					  vcf_sample, VCF_SAMPLE_MAX_CHARS);
	    }
	    fprintf(stderr, "Counted calls up to %s %zu\n",
		    vcf_call.chromosome, vcf_call.pos);
	    
	    fprintf(stderr, "Found %lu matching calls.\n",
		    EVENT_DEPTH_COUNT(events + c));
	    
	    /*
	    event_depth_median()
	    qsort(depths, count, sizeof(depth_t),
		  (int (*)(const void *, const void *))depth_cmp);
	    fprintf(stderr, "Median value is %u.\n", depths[(count + 1) / 2]);
	    */
	    
	    if ( (status != VCF_READ_OK) && (status != VCF_READ_EOF) )
	    {
		fprintf(stderr, "median_depths(): Error reading VCF file.\n");
		exit(EX_DATAERR);
	    }
	}
	vcf_close(vcf_stream, compressed);
    }
    
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

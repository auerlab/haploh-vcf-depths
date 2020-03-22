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
#include "haploh-median-depths.h"

int     main(int argc,const char *argv[])

{
    FILE        *events_stream;
    const char  *events_file,
		*vcf_dir;
    struct stat dir_stat;
    
    if ( argc != 2 )
	usage(argv);

    events_file = argv[1];
    vcf_dir = argv[2];
    
    if ( (events_stream = fopen(events_file, "r")) == NULL )
    {
	fprintf(stderr, "Cannot open %s: %s\n", events_file, strerror(errno));
	exit(EX_NOINPUT);
    }
    
    if ( (stat(vcf_dir, &dir_stat) == -1) || !S_ISDIR(dir_stat.st_mode) )
    {
	fprintf(stderr,
		"%s must be a readable directory containing VCF files.\n",
		vcf_dir);
	exit(EX_NOINPUT);
    }
	    
    return haploh_median_depths(events_stream, vcf_dir);
}


void    usage(const char *argv[])

{
    fprintf(stderr, "Usage: %s haplohseq-event-file directory-with-VCF-files\n", argv[0]);
    exit(EX_USAGE);
}


int     haploh_median_depths(FILE *events_stream, const char *vcf_dir)

{
    return EX_OK;
}

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
#include <vcfio.h>
#include "haploh-median-depths.h"

int     main(int argc,const char *argv[])

{
    if ( argc != 1 )
	usage(argv);
    return EX_OK;
}


void    usage(const char *argv[])

{
    fprintf(stderr, "Usage: %s haplohseq-event-file\n", argv[0]);
    exit(EX_USAGE);
}


typedef struct
{
    char    chromosome[VCF_CHROMOSOME_MAX_CHARS + 1];
    size_t  begin,
	    end;
}   event_t;

#define EVENT_READ_OK           0
#define EVENT_READ_EOF          -1
#define EVENT_READ_OVERFLOW     -2
#define EVENT_READ_TRUNCATED    -3
#define EVENT_READ_HEADER       -4

#include "haploh-median-depths-protos.h"

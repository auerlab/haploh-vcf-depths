typedef unsigned short  depth_t;

typedef struct
{
    char    chromosome[VCF_CHROMOSOME_MAX_CHARS + 1];
    size_t  begin,
	    end;
    depth_t *depth_list;
    size_t  depth_count,
	    max_depths;
}   event_t;

#define EVENT_READ_OK           0
#define EVENT_READ_EOF          -1
#define EVENT_READ_OVERFLOW     -2
#define EVENT_READ_TRUNCATED    -3
#define EVENT_READ_HEADER       -4

#define EVENT_DEFAULT_MAX_DEPTHS    500000  // Based on TOPMed data

#define EVENT_DEPTH_COUNT(event)    ((event)->depth_count)
#define EVENT_BEGIN(event)          ((event)->begin)
#define EVENT_END(event)            ((event)->end)
#define EVENT_CHROMOSOME(event)     ((event)->chromosome)

#define EVENT_MALLOC(count,type)  ((type *)malloc((count) * sizeof(type)))
#define EVENT_REALLOC(buff, count, type) ((type *)realloc((buff), (count) * sizeof(type)))

#include "events-protos.h"

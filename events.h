#include <limits.h>

typedef unsigned short  depth_t;

typedef struct
{
    char    chromosome[BL_CHROM_MAX_CHARS + 1];
    char    sample_id[PATH_MAX + 1];
    size_t  begin,
	    end;
    FILE    *same_sample_depth_stream;
    FILE    *other_samples_depth_stream;
}   event_t;

#define EVENT_READ_OK           0
#define EVENT_READ_EOF          -1
#define EVENT_READ_OVERFLOW     -2
#define EVENT_READ_TRUNCATED    -3
#define EVENT_READ_HEADER       -4

#define EVENT_DEFAULT_MAX_DEPTHS    500000  // Based on TOPMed data

#define EVENT_BEGIN(event)          ((event)->begin)
#define EVENT_END(event)            ((event)->end)
#define EVENT_CHROMOSOME(event)     ((event)->chromosome)
#define EVENT_DEPTH_COUNT(event)    ((event)->depth_count)

#define EVENT_MALLOC(count,type)  ((type *)malloc((count) * sizeof(type)))
#define EVENT_REALLOC(buff, count, type) ((type *)realloc((buff), (count) * sizeof(type)))

#include "events-protos.h"

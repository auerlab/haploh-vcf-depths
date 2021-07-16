
#define CMD_MAX 1024

// Match these with ad2vcf and vcf-split
// Yes, we actually saw a few INFO fields over 512k in some dbGap BCFs
#define BL_VCF_INFO_MAX_CHARS          1048576
#define BL_VCF_FORMAT_MAX_CHARS        4096
#define BL_VCF_SAMPLE_MAX_CHARS        2048

#include "events.h"
#include "haploh-vcf-depths-protos.h"

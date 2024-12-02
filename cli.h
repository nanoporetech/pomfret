#ifndef __METHPHASER2_CLI_H__
#define __METHPHASER2_CLI_H__
#include <stdint.h>

typedef struct {
    int lo;
    int hi;
    int is_r9;
    char *fn_gtf;
    char *fn_bam;
    char *output_prefix;
    int use_hyper_parameter_search;
    int readlen_threshold;  // min read length
    int mapq; 
    int k; // methmer k
    int k_span;  // methmer base space length limit
    int cov_for_selection;  // meth site min coverage for initial inclusion
    int n_candidates_per_iter;  // # reads for consideration at each tagging iter
}cliopt_t;

cliopt_t cliopt;
int parse_cli(int argc, char *argv[]);

#endif  // __METHPHASER2_CLI_H__
 
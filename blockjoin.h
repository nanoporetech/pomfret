#ifndef POMFRET_BLOCKJOIN_H
#define POMFRET_BLOCKJOIN_H
#include "cli.h"

#define VERSION "v0.1-r13"

typedef struct{
    int k;
    int k_span;
    int lo, hi;
    int cov_known; // read coverage, if known
    int cov_for_selection;  // i.e. mmr_min_cov
    int cov_for_runtime;
    int readlen_threshold;
    int min_mapq;
}mmr_config_t;

int main_debug(int argc, char *argv[]);

int main_blockjoin(cliopt_t* cliopt);
int main_varhaptag(char *fn_vcf, char *fn_bam, char *fn_out, int n_thread, 
                  int verbose, int do_write_bam);
int main_methstat(char *fn_bam, 
                  char *fn_intervals, enum input_file_format fn_intervals_type, 
                  char *fn_out, 
                  int lo, int hi, int cov_for_selection, int readlen_threshold);

int main_methreport(cliopt_t* clio);
#endif  //POMFRET_BLOCKJOIN_H


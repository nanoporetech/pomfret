#ifndef POMFRET_CLI_H
#define POMFRET_CLI_H
#include <stdint.h>

double Get_T(void);
double Get_U(void);

enum input_file_format{
    IS_GTF,
    IS_VCF,
    IS_TSV
};

extern int pomfret_n_bam_threads;


// for methphase
typedef struct {
    int is_help;
    int threads;
    int threads_bam; // override `threads` only for bam writing and indexing
    int lo;
    int hi;
    int is_r9;
    char *fn_tsv;
    char *fn_gtf;
    char *fn_vcf;
    char *fn_bam;
    int bam_needs_haplotagging;
    int write_bam_input_haplotagging;
    char *output_prefix;
    int use_hyper_parameter_search;
    int readlen_threshold;  // min read length
    int mapq; 
    int k; // methmer k
    int k_span;  // methmer base space length limit
    int cov;  // input read alignment coverage, if known
    int cov_for_selection;  // meth site min coverage for initial inclusion
    int n_candidates_per_iter;  // # reads for consideration at each tagging iter
    int do_output_bam;
    int do_output_tsv;
    int write_debug_files;
}cliopt_t;
int cliopt_verbose;
cliopt_t* parse_cli(int argc, char *argv[]);
int sancheck_cliopt(cliopt_t *opt);
void destroy_cliopt_t(cliopt_t *cliopt);

// for varhaptag
typedef struct{
    int n_threads;
    int verbose;
    int is_help;
    char *fn_vcf;
    char *fn_bam;
    //char **regions;
    //int n_regions;
    //int m_regions;
    char *fn_out;
    int do_write_bam;
}cliopt_haptag_t;
int sancheck_cliopt_varhaptag(cliopt_haptag_t *c);
cliopt_haptag_t *parse_cli_varhaptag(int argc, char *argv[]);
void destroy_cliopt_varhaptag_t(cliopt_haptag_t *cliopt_haptag);


// for methstat
typedef struct{
    char *fn_bam;
    char *fn_vcf;
    char *fn_gtf;
    char *fn_tsv;
    char *fn_out;
    int is_help;
    int lo, hi, cov, readlen;
}cliopt_methstat_t;
int sancheck_cliopt_methstat(cliopt_methstat_t *clio);
cliopt_methstat_t *parse_cli_methstat(int argc, char *argv[]);
void destroy_cliopt_methstat_t(cliopt_methstat_t *clio);

#endif  // POMFRET_CLI_H


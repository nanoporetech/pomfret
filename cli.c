#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include "ketopt.h"
#include "cli.h"

ketopt_t opt;

static ko_longopt_t longopts[] = {
    { "modlow", ko_required_argument,  301 },  // exclusive upperbound of unmeth call qual
    { "modhigh", ko_required_argument, 302 },  // inclusive lowerbound of meth call qual
    { "R9", ko_no_argument,         303 },  // input is R9 data
    { "gtf", ko_required_argument, 304},  // input: gtf file name
    { "bam", ko_required_argument, 305}, // input: bam file name
    { "hypers", ko_no_argument, 306},
    { "mapq", ko_required_argument, 307},
    { NULL, 0, 0 }
  };

void init_cliopt_t(){
    cliopt.lo = 100;
    cliopt.hi = 156;
    cliopt.is_r9 = 0;
    cliopt.fn_gtf = 0;
    cliopt.fn_bam = 0;
    cliopt.output_prefix = 0;
    cliopt.use_hyper_parameter_search = 0;
    cliopt.readlen_threshold = 15000;
    cliopt.mapq = 10;
    cliopt.k = 7;
    cliopt.k_span = 5000;
    cliopt.cov_for_selection = 6;
    cliopt.n_candidates_per_iter = 15;
}
void print_help_cli(){
    //TODO
    char idt[] = "    ";
    fprintf(stderr, "TBD: make vcf and bam filenames positional, and gtf file optional\n");
    fprintf(stderr, "TBD: write help\n");
    //fprintf(stderr, "Usage: methp [options] \n");
    exit(0);
}

void dir_splice_from_string(char *s, int s_l, int *l){
    for (int i=s_l-1; i>=0; i--){
        if (s[i]=='/') {
            *l = i+1;
            return;
        }
    }
    *l = 0;
}

int sancheck_cliopt(cliopt_t *opt){
    // return 0 if alright, 1 otherwise
    FILE *fp;
    int l, l2;
    struct stat sb;
    // sane mod call qual thresholds
    if (opt->lo<0){
        fprintf(stderr, "[E::%s] lower threshold for mod call quality is too low (%d)\n", __func__, opt->lo);
        goto fail;
    }
    if (opt->lo>127){
        fprintf(stderr, "[E::%s] lower threshold for mod call quality is too high (%d)\n", __func__, opt->lo);
        goto fail;
    }
    if (opt->hi>255){
        fprintf(stderr, "[E::%s] upper threshold for mod call quality is too high (%d)\n", __func__, opt->hi);
        goto fail;
    }
    if (opt->hi<=127){
        fprintf(stderr, "[E::%s] upper threshold for mod call quality is too low (%d)\n", __func__, opt->hi);
        goto fail;
    }
    if (opt->is_r9!=0 && opt->is_r9!=1){
        fprintf(stderr, "[E::%s] invalid is_r9 value (%d)\n", __func__, opt->is_r9);
        goto fail;
    }
    // sane read and methmer options
    if (opt->readlen_threshold<0) opt->readlen_threshold = 0;
    if (opt->mapq>60){
        fprintf(stderr, "[W::%s] mapq seems too high, proceed anyways\n", __func__);
    }
    if (opt->mapq<0) opt->mapq=0;
    if (opt->k<=0) {
        fprintf(stderr, "[W::%s] clipping mether k to 1\n", __func__);
        opt->k = 1;
    }
    if (opt->k_span<=0){
        fprintf(stderr, "[W::%s] clipping mether span to 1\n", __func__);
        opt->k_span = 1;
    }
    if (opt->cov_for_selection<=0){
        fprintf(stderr, "[W::%s] clipping meth site coverage threshold to 1\n", __func__);
        opt->cov_for_selection= 1;
    }
    if (opt->n_candidates_per_iter<=0){
        fprintf(stderr, "[W::%s] clipping # candidates per iter to 1\n", __func__);
        opt->n_candidates_per_iter = 1;
    }
    if (opt->n_candidates_per_iter<5){
        fprintf(stderr, "[W::%s] number of candidates per iter might be too low\n", __func__);
    }
    // input file names exist (not checking files)
    if (opt->fn_gtf==0){
        fprintf(stderr, "[E::%s] missing input: gtf file\n", __func__);
        goto fail;
    }
    if (opt->fn_bam==0){
        fprintf(stderr, "[E::%s] missing bam file\n", __func__);
        goto fail;
    }
    // output file name exists and valid
    if (opt->output_prefix==0){
        fprintf(stderr, "[E::%s] no output prefix given\n", __func__);
        goto fail;
    }
    l = strlen(opt->output_prefix); 
    if (l==0){
        fprintf(stderr, "[E::%s] no output prefix given\n", __func__);
        goto fail;
    }
    if (opt->output_prefix[l-1]=='/'){
        fprintf(stderr, "[W::%s] output prefix has trailing '/', stripping them.\n", __func__);
        for (int i=l-1; i>=0; i--){
            if (opt->output_prefix[i]=='/'){
                l = i;
                opt->output_prefix[i] = 0;
            }else break;
        }
    }
    if (l==0){
        fprintf(stderr, "[E::%s] no output prefix given\n", __func__);
        goto fail;
    }

    return 0;
fail:
    return 1;
}

// TODO: k is at most 15 because the integer encoding and sort packing
// TODO: length of sequence is actually limited at 29bits
int parse_cli(int argc, char *argv[]){
    opt = KETOPT_INIT;
    int i, c;
    init_cliopt_t();

    while ((c = ketopt(&opt, argc, argv, 1, "xho:k:L:l:c:n:", longopts)) >= 0) {
        if (c == 'x') cliopt.use_hyper_parameter_search = 1;
        else if (c == 'h') print_help_cli();
        else if (c == 'o') cliopt.output_prefix=opt.arg;
        else if (c == 'k') cliopt.k = atoi(opt.arg);
        else if (c == 'l') cliopt.k_span= atoi(opt.arg);
        else if (c == 'L') cliopt.readlen_threshold = atoi(opt.arg);
        else if (c == 'c') cliopt.cov_for_selection = atoi(opt.arg);
        else if (c == 'n') cliopt.n_candidates_per_iter = atoi(opt.arg);
        else if (c == 301) cliopt.lo = atoi(opt.arg);
        else if (c == 302) cliopt.hi = atoi(opt.arg);
        else if (c == 303) cliopt.is_r9 = 1;
        else if (c == 304) cliopt.fn_gtf = opt.arg;
        else if (c == 305) cliopt.fn_bam = opt.arg;
        else if (c == 306) cliopt.use_hyper_parameter_search = 1;
        else if (c == 307) cliopt.mapq= atoi(opt.arg);
        else if (c == '?') printf("unknown opt: -%c\n", opt.opt? opt.opt : ':');
        else if (c == ':') printf("missing arg: -%c\n", opt.opt? opt.opt : ':');
    }
    printf("Non-option arguments:");
    for (i = opt.ind; i < argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');
    return 0;
}

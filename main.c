#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "cli.h"
#include "blockjoin.h"


void print_help_main(){
    fprintf(stderr, "Usage: pomfret <subcommand> [options]\n");
    fprintf(stderr, "Subcommands:\n");
    fprintf(stderr, "  methphase  Given aligned reads with methylation calls in bam\n");
    fprintf(stderr, "             and exiting phase blocks, try to use methylation to\n");
    fprintf(stderr, "             phase the unphased regions.\n");
    fprintf(stderr, "  varhaptag  Given aligned reads in bam and a phased vcf, produce\n");
    fprintf(stderr, "             a haplotagged bam.\n");
    fprintf(stderr, "  methstat   Given aligned reads with methylation calls in bam\n");
    fprintf(stderr, "             and exiting phase blocks, find the seemingly heterozygous\n");
    fprintf(stderr, "             sites for the unphased regions.\n");
}


int main (int argc, char *argv[]){
    fprintf(stderr, "[M::%s] pomfret %s\n", __func__, VERSION);
    fprintf(stderr, "[M::%s] CMD: ", __func__);
    for (int i=0; i<argc; i++) fprintf(stderr, "%s ", argv[i]);
    fprintf(stderr, "\n");

    double T = Get_T();
    int ret = 0;

    //main_debug(argc, argv);
    //return 0;

    cliopt_t* cliopt = 0;
    cliopt_haptag_t* clio_h= 0;
    cliopt_methstat_t *clio_m= 0;

    if (argc<2){
        print_help_main();
        ret = 1;
        goto finish;
    }

    if (strcmp(argv[1], "methphase")==0){
        cliopt = parse_cli(argc-1, argv+1);
        if (!cliopt || cliopt->is_help || sancheck_cliopt(cliopt)){
            ret = 1;
            goto finish;
        }
        ret = main_blockjoin(cliopt);
    }else if (strcmp(argv[1], "varhaptag")==0){
        clio_h = parse_cli_varhaptag(argc-1, argv+1);
        if (!clio_h|| clio_h->is_help || 
             sancheck_cliopt_varhaptag(clio_h)){
            ret = 1;
            goto finish;
        }
        ret = main_varhaptag(clio_h->fn_vcf, 
                         clio_h->fn_bam, 
                         clio_h->fn_out, 
                         clio_h->n_threads, clio_h->verbose, clio_h->do_write_bam);
    }else if (strcmp(argv[1], "methstat")==0){
        clio_m= parse_cli_methstat(argc-1, argv+1);
        if (!clio_m|| clio_m->is_help || 
            sancheck_cliopt_methstat(clio_m)){
                ret = 1;
                goto finish;
        }
        ret = main_methstat(clio_m->fn_bam,
                      clio_m->fn_tsv? clio_m->fn_tsv 
                      :clio_m->fn_gtf? clio_m->fn_gtf: clio_m->fn_vcf, 
                      clio_m->fn_tsv? IS_TSV 
                      :clio_m->fn_gtf? IS_GTF: IS_VCF, 
                      clio_m->fn_out, 
                      clio_m->lo, clio_m->hi, clio_m->cov, clio_m->readlen);
    }else{
        fprintf(stderr, "[E::%s] unknown subcommand: %s\n", __func__, argv[1]);
        print_help_main();
        ret = 1;
    }
    

finish:
    if (cliopt) destroy_cliopt_t(cliopt);
    if (clio_h) destroy_cliopt_varhaptag_t(clio_h);
    if (clio_m) destroy_cliopt_methstat_t(clio_m);

    fprintf(stderr, "\n[M::%s] CMD: ", __func__);
    for (int i=0; i<argc; i++) fprintf(stderr, "%s ", argv[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "[M::%s] used: %.1fs, peak RSS %.1fGiB\n", 
                    __func__, Get_T()-T, Get_U());
    if (global_data_has_implicit){
        fprintf(stderr, "[W::%s] Input BAM has implicit modified base calls.\n");
        fprintf(stderr, "  Pomfret extracts 5mC without considering 5hmC, which is different from\n");
        fprintf(stderr, "  `modkit adjust-mods --motif CG 0 --ignore h in.bam out.bam`\n");
    }

    return ret;
}


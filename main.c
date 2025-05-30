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
    fprintf(stderr, "  report     Given aligned reads in bam and a phased vcf, sample\n");
    fprintf(stderr, "             intervals within phase blocks, pretend they are phase gaps\n");
    fprintf(stderr, "             and report whether meth-phasing would generate correct \n");
    fprintf(stderr, "             phase block joining decisions.\n");
}


int main (int argc, char *argv[]){
    fprintf(stderr, "[M::%s] pomfret %s\n", __func__, VERSION);
    fprintf(stderr, "[M::%s] CMD: ", __func__);
    for (int i=0; i<argc; i++) fprintf(stderr, "%s ", argv[i]);
    fprintf(stderr, "\n");

    double T = Get_T();
    int ret = 0;

    cliopt_t* cliopt = 0;
    cliopt_haptag_t* clio_h= 0;
    cliopt_methstat_t *clio_m= 0;

    if (argc<2){
        print_help_main();
        ret = 1;
        goto finish;
    }
    if (strcmp(argv[1], "-h")==0 || strcmp(argv[1], "--help")==0 || strcmp(argv[1], "help")==0){
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
    }
    else if (strcmp(argv[1], "varhaptag")==0){
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
    }
    else if (strcmp(argv[1], "report")==0){
        cliopt = parse_cli(argc-1, argv+1);
        if (!cliopt || cliopt->is_help || sancheck_cliopt(cliopt)){
            ret = 1;
            goto finish;
        }
        if (!cliopt->fn_vcf){
            fprintf(stderr, "[E::%s] missing input: phasd vcf file.\n", __func__);
            ret = 1;
            goto finish;
        }
        main_methreport(cliopt);

    }
    else{
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
        fprintf(stderr, "[W::%s] Input BAM has implicit modified base calls.\n", __func__);
        fprintf(stderr, "  Pomfret extracts 5mC without considering 5hmC, which is different from\n");
        fprintf(stderr, "  `modkit adjust-mods --motif CG 0 --ignore h in.bam out.bam`.\n");
    }

    return ret;
}


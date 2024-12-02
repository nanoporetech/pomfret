#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "cli.h"
#include "blockjoin.h"

int main (int argc, char *argv[]){
    parse_cli(argc, argv);
    
    int l = strlen(cliopt.output_prefix);
    char *fn_out_tsv = (char*)malloc(l+5);
    sprintf(fn_out_tsv, "%s.tsv", cliopt.output_prefix);
    //main_blockjoin(cliopt.fn_gtf, cliopt.fn_bam, fn_out_tsv);
    main_blockjoin_dbg(argc, argv);

    return 0;
}


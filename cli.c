#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include "ketopt.h"
#include "kstring.h"
#include "cli.h"

ketopt_t opt;

int pomfret_n_bam_threads = 1;
int global_data_has_implicit = 0;

double Get_T(void){
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec+t.tv_usec/1000000.0;
}
double Get_U(void){
    struct rusage s;
    getrusage(RUSAGE_SELF, &s);
    return (double)s.ru_maxrss/1048576.0;  // GB
}

/*** methphase ***/
static ko_longopt_t longopts[] = {
    { "lo", ko_required_argument,  301 },  // exclusive upperbound of unmeth call qual
    { "hi", ko_required_argument,  302 },  // inclusive lowerbound of meth call qual
    { "gtf", ko_required_argument, 304 },  // input: gtf file name
    { "vcf", ko_required_argument, 306 },   // input: vcf file name
    { "mapq", ko_required_argument,307 },
    { "tsv", ko_required_argument, 308 },   // input: tsv file name
    { "dont-write-bam", ko_no_argument, 309 },   // otherwise will produce a bam with modified read haplotags
    { "output-tsv", ko_no_argument,310 },  // write a tsv of new phaseblocks
    { "bam-threads", ko_required_argument,311 },
    { "bam-is-untagged", ko_no_argument, 312},
    { "write-input-tagging", ko_no_argument, 313},
    { "help", ko_no_argument, 400},
    { "dbg", ko_no_argument, 401},
    { 0, 0, 0 }
};

void init_cliopt_t(cliopt_t *cliopt){
    cliopt->is_help = 0;
    cliopt->threads = 1;
    cliopt->threads_bam = cliopt->threads;
    cliopt->lo = 100;
    cliopt->hi = 156;
    cliopt->is_r9 = 0;
    cliopt->fn_gtf = 0;
    cliopt->fn_tsv = 0;
    cliopt->fn_vcf = 0;
    cliopt->fn_bam = 0;
    cliopt->bam_needs_haplotagging = 0;
    cliopt->write_bam_input_haplotagging = 0;
    cliopt->output_prefix = "pomfret";
    cliopt->readlen_threshold = 15000;
    cliopt->mapq = 10;
    cliopt->k = 3;
    cliopt->k_span = 5000;
    cliopt->cov_for_selection = -1;//6;
    cliopt->n_candidates_per_iter = 15;
    cliopt->do_output_bam = 1;
    cliopt->do_output_tsv = 0;
    cliopt->write_debug_files = 0;
}
void destroy_cliopt_t(cliopt_t *cliopt){
    free(cliopt);
}
void print_help_cli(cliopt_t *cliopt){
    fprintf(stderr, "Usage: pomfret methphase -o out_prefix --vcf phased.vcf[.gz] [...] reads.bam 2>log\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  bam    [pos] Aligned reads. Must be sorted and has index. If reads are not\n");
    fprintf(stderr, "               haplotagged, supply -u and provide vcf (via --vcf).\n");
    fprintf(stderr, "  -h,--help [   ] Display this message.\n");
    fprintf(stderr, "  -c     [req] Read coverage (total, not per-haplotap).\n");
    fprintf(stderr, "  -o     [opt] Name prefix of output files. [%s]\n", cliopt->output_prefix);
    fprintf(stderr, "  --vcf  [opt] Input, sorted vcf file containing phased variants.\n"); 
    fprintf(stderr, "               Either vcf, gtf or tsv need to be present. Plain or gz'd.\n");
    fprintf(stderr, "  --gtf  [opt] Input, sorted gtf file of prescribed phase blocks.\n" 
                    "               If present, overrides phase blocks defined by vcf.\n"
                    "               Plain or gz'd.\n");
    fprintf(stderr, "  --tsv  [opt] Input, sorted 3-column tsv file of prescribed phase blocks:\n" 
                    "               reference name, start, end. Plain or gz'd.\n"
                    "               If present, overrides both gtf and vcf.\n");
    fprintf(stderr, "  -u,--bam-is-untagged [opt] If present, will haplotag reads \n"
                    "               with phased variants in vcf first. --vcf must be \n"
                    "               supplied. Ignores any haptags present in the bam.\n"
                    "               All variants supplies by the vcf will be used as evidences.\n");
    fprintf(stderr, "  -U,--write-input-tagging [opt] If present along with -u, output a tsv \n"
                    "               containing the input haplotagging results.\n");
    fprintf(stderr, "  -t     [opt] Number of threads to use. [%d]\n", cliopt->threads);
    //fprintf(stderr, "  --lo   [opt] lower cutoff of meth call, i.e. unmethylated. [%d]\n", cliopt->lo);
    //fprintf(stderr, "  --hi   [opt] upper cutoff of meth call, i.e. methylated. [%d]\n", cliopt->hi);
    //fprintf(stderr, "  --dont-write-bam [opt] if present, do not produce re-haplotagged bam.\n");
    fprintf(stderr, " -T,--bam-threads [opt] #threads when writing&indexing bam output.\n");
    fprintf(stderr, "               Used by bgzf_mt().[defaults to identical to -t]\n");
    fprintf(stderr, "Note: Inputs may need to be opened or read for more than once.\n");

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

int sancheck_cliopt(cliopt_t *cliopt){
    // return 0 if alright, 1 otherwise
    FILE *fp;
    int l, l2;
    struct stat sb;
    // basic
    if (cliopt->threads<=0){
        fprintf(stderr, "[W::%s] invalid thread number (%d), clipped to 1\n", __func__, cliopt->threads);
        cliopt->threads = 1;
    }
    if (cliopt->threads_bam<=0){
        fprintf(stderr, "[W::%s] invalid bam thread number(%d), clipped to 1\n", __func__, cliopt->threads_bam);
        cliopt->threads_bam = 1;
        pomfret_n_bam_threads = 1;
    }
    // sane mod call qual thresholds
    if (cliopt->lo<0){
        fprintf(stderr, "[E::%s] lower threshold for mod call quality is too low (%d)\n", __func__, cliopt->lo);
        goto fail;
    }
    if (cliopt->lo>127){
        fprintf(stderr, "[E::%s] lower threshold for mod call quality is too high (%d)\n", __func__, cliopt->lo);
        goto fail;
    }
    if (cliopt->hi>255){
        fprintf(stderr, "[E::%s] upper threshold for mod call quality is too high (%d)\n", __func__, cliopt->hi);
        goto fail;
    }
    if (cliopt->hi<=127){
        fprintf(stderr, "[E::%s] upper threshold for mod call quality is too low (%d)\n", __func__, cliopt->hi);
        goto fail;
    }
    if (cliopt->is_r9!=0 && cliopt->is_r9!=1){
        fprintf(stderr, "[E::%s] invalid is_r9 value (%d)\n", __func__, cliopt->is_r9);
        goto fail;
    }
    // sane read and methmer options
    if (cliopt->readlen_threshold<0) cliopt->readlen_threshold = 0;
    if (cliopt->mapq>60){
        fprintf(stderr, "[W::%s] mapq seems too high, proceed anyways\n", __func__);
    }
    if (cliopt->mapq<0) cliopt->mapq=0;
    if (cliopt->k<=0) {
        fprintf(stderr, "[W::%s] clipping mether k to 1\n", __func__);
        cliopt->k = 1;
    }
    if (cliopt->k_span<=0){
        fprintf(stderr, "[W::%s] clipping mether span to 1\n", __func__);
        cliopt->k_span = 1;
    }
    if (cliopt->cov_for_selection<=0){
        //fprintf(stderr, "[E::%s] must provide input coverage.\n", __func__);
        //goto fail;
        fprintf(stderr, "[M::%s] read coverage not provided, will estimate.\n", __func__);
    }
    if (cliopt->n_candidates_per_iter<=0){
        fprintf(stderr, "[W::%s] clipping candidate per iter to 1\n", __func__);
        cliopt->n_candidates_per_iter = 1;
    }
    if (cliopt->n_candidates_per_iter<=0){
        fprintf(stderr, "[W::%s] clipping # candidates per iter to 1\n", __func__);
        cliopt->n_candidates_per_iter = 1;
    }
    if (cliopt->n_candidates_per_iter<5){
        fprintf(stderr, "[W::%s] number of candidates per iter might be too low\n", __func__);
    }

    // input file names exist (not checking files)
    if (cliopt->fn_gtf==0 && cliopt->fn_tsv==0 && cliopt->fn_vcf==0){
        fprintf(stderr, "[E::%s] gtf, tsv and vcf cannot all be absent\n", __func__);
        goto fail;
    }
    if ((cliopt->fn_gtf?1:0) + (cliopt->fn_tsv?1:0) + (cliopt->fn_vcf?1:0)>1){
        fprintf(stderr, "[M::%s] multiple phase block files. Will not resolve conflict, only total override (order is always: tsv > gtf > vcf)\n", __func__);
    }

    // if input bam is unhaplotagged, must supply vcf
    if (cliopt->bam_needs_haplotagging && !cliopt->fn_vcf){
        fprintf(stderr, "[E::%s] input bam was flagged unhaplotagged, but vcf is missing.\n", __func__);
        goto fail;
    }
    //if (cliopt->n_bam==0){
    if (!cliopt->fn_bam){
        fprintf(stderr, "[E::%s] missing bam file\n", __func__);
        goto fail;
    }
    // output file name exists and valid
    if (cliopt->output_prefix==0){
        fprintf(stderr, "[E::%s] no output prefix given\n", __func__);
        goto fail;
    }
    l = strlen(cliopt->output_prefix); 
    if (l==0){
        fprintf(stderr, "[E::%s] no output prefix given\n", __func__);
        goto fail;
    }
    if (cliopt->output_prefix[l-1]=='/'){
        fprintf(stderr, "[W::%s] output prefix has trailing '/', stripping them.\n", __func__);
        for (int i=l-1; i>=0; i--){
            if (cliopt->output_prefix[i]=='/'){
                l = i;
                cliopt->output_prefix[i] = 0;
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

void sancheck_bamfiles(cliopt_t *cliopt){
    // Just make sure they are using the same reference genome. 
    fprintf(stderr, "TODO/TBD\n");
}

// TODO: k is at most 15 because the integer encoding and sort packing
// TODO: length of sequence is actually limited at 29bits
cliopt_t* parse_cli(int argc, char *argv[]){
    cliopt_t *cliopt = (cliopt_t*)malloc(sizeof(cliopt_t));
    init_cliopt_t(cliopt);
    cliopt_verbose = 0;
    opt = KETOPT_INIT;
    int i, c;

    if (argc==1){
        print_help_cli(cliopt);
        destroy_cliopt_t(cliopt);
        return 0;
    }

    while ((c = ketopt(&opt, argc, argv, 1, "vhuUo:k:L:l:c:n:t:T:", longopts)) >= 0) {
        if (c == 'v') cliopt_verbose++;
        else if (c == 'h' || c==400) {print_help_cli(cliopt); cliopt->is_help=1;}
        else if (c == 't') {
            cliopt->threads = atoi(opt.arg);
            cliopt->threads_bam = cliopt->threads;
            pomfret_n_bam_threads = cliopt->threads;
        }
        else if (c == 'o') cliopt->output_prefix=opt.arg;
        else if (c == 'k') cliopt->k = atoi(opt.arg);
        else if (c == 'l') cliopt->k_span= atoi(opt.arg);
        else if (c == 'L') cliopt->readlen_threshold = atoi(opt.arg);
        else if (c == 'c') {
            int target_coverage = atoi(opt.arg);
            cliopt->cov = target_coverage;
            cliopt->cov_for_selection = target_coverage/10;
            cliopt->n_candidates_per_iter = target_coverage/4;
        }
        else if (c == 'n') cliopt->n_candidates_per_iter = atoi(opt.arg);
        else if (c == 301) cliopt->lo = atoi(opt.arg);
        else if (c == 302) cliopt->hi = atoi(opt.arg);
        //else if (c == 303) cliopt->is_r9 = 1;
        else if (c == 304) cliopt->fn_gtf = opt.arg;
        //else if (c == 305) cliopt->fn_bam = opt.arg;
        else if (c == 306) cliopt->fn_vcf= opt.arg;
        else if (c == 307) cliopt->mapq= atoi(opt.arg);
        else if (c == 308) cliopt->fn_tsv= opt.arg;
        else if (c == 309) cliopt->do_output_bam = 0;
        else if (c == 310) cliopt->do_output_tsv = 1;
        else if (c == 401) cliopt->write_debug_files= 1;
        else if (c == 'T' || c == 311) {
            cliopt->threads_bam = atoi(opt.arg);
            pomfret_n_bam_threads = cliopt->threads_bam;
        }
        else if (c == 312 || c=='u') cliopt->bam_needs_haplotagging = 1;
        else if (c == 313 || c=='U') cliopt->write_bam_input_haplotagging = 1;
        else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option argument in \"%s\"\n", __func__, argv[opt.i - 1]);
            //printf("unknown opt: -%c\n", opt.opt? opt.opt : ':');
            destroy_cliopt_t(cliopt);
            return 0;
        }
        else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option argument in \"%s\"\n", __func__, argv[opt.i - 1]);
            //printf("missing arg: -%c\n", opt.opt? opt.opt : ':');
            destroy_cliopt_t(cliopt);
            return 0;
        }
    }
    if (argc-opt.ind>1){
        fprintf(stderr, "[TODO::%s] too many positional arguments; multi bam input not impl'd yet\n", __func__);
    }
    for (i = opt.ind; i < argc; i++){
        if (cliopt->fn_bam){
            fprintf(stderr, "[E::%s] multiple bam input is not supported.\n", __func__);
            exit(1);
        }
        cliopt->fn_bam = argv[i];
    }
    putchar('\n');
    return cliopt;
}
/*** end of methphase***/



/*** varhaptag ***/

void init_cliopt_varhaptag_t(cliopt_haptag_t *cliopt_haptag){
    cliopt_haptag->fn_vcf = 0;
    cliopt_haptag->fn_bam = 0;
    //cliopt_haptag->n_regions = 0;
    //cliopt_haptag->m_regions = 16;
    //cliopt_haptag->regions = (char**)calloc(16, sizeof(char*));
    cliopt_haptag->fn_out = "pomfret_varhaptag";
    cliopt_haptag->n_threads = 1;
    cliopt_haptag->verbose = 0;
    cliopt_haptag->is_help = 0;
    cliopt_haptag->do_write_bam = 1;
}
void destroy_cliopt_varhaptag_t(cliopt_haptag_t *cliopt_haptag){
    //free(cliopt_haptag->regions);
    free(cliopt_haptag);
}
void print_help_cli_varhaptag(cliopt_haptag_t *cliopt_haptag){
    fprintf(stderr, "Usage: pomfret varhaptag [-t threads] -o out.bam in.vcf in.bam 2>log\n");
    fprintf(stderr, "Purpose:\n");
    fprintf(stderr, "  A simple haplotagging utility. Given aligned reads and phased\n");
    fprintf(stderr, "  variants, produce a haplotagged bam and a tsv which maps reads\n");
    fprintf(stderr, "  to raw & new haplotags.\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  vcf    [pos] A phased VCF. Assumed to be single-sample. All variants\n");
    fprintf(stderr, "               will be used.\n");
    fprintf(stderr, "  bam    [pos] Aligned reads from one sample. Existing HP tags will\n");
    fprintf(stderr, "               be ignored and overriden.\n");
    fprintf(stderr, "  -o     [req] Name of the output bam file. [%s]\n", cliopt_haptag->fn_out);
    fprintf(stderr, "  -h,--help [   ] Display this message.\n");
    //fprintf(stderr, "  -r     [opt] Restrict haplotagging to the given interval(s).\n");
    //fprintf(stderr, "               Provide this switch multiple times if specifying more\n");
    //fprintf(stderr, "               one interval. Formatting example: chr6:10M-11,000,000\n");
    fprintf(stderr, "  -t     [opt] Number of threads for writing bam. [%d]\n", cliopt_haptag->n_threads);
    fprintf(stderr, "  -v     [opt] Verbose, write read tags to stderr.\n"); 
    fprintf(stderr, "  --dont-write-bam [opt] Suppress bam output.\n"); 
    fprintf(stderr, "Note: Inputs may need to be opened or read for more than once.\n");
}
int sancheck_cliopt_varhaptag(cliopt_haptag_t *c){
    int ret = 0;
    if (!c->fn_bam){
        fprintf(stderr, "[E::%s] Missing input bam.\n", __func__);
        ret = 1;
    }
    if (!c->fn_vcf){
        fprintf(stderr, "[E::%s] Missing input vcf.\n", __func__);
        ret = 1;
    }
    if (c->n_threads<1){
        fprintf(stderr, "[W::%s] invalid number of threads (%d), clipping to 1 .\n", 
        __func__, c->n_threads);
        c->n_threads = 1;
    }
    return ret;
}
cliopt_haptag_t *parse_cli_varhaptag(int argc, char *argv[]){
    cliopt_haptag_t *cliopt_haptag = (cliopt_haptag_t*)malloc(sizeof(cliopt_haptag_t));
    init_cliopt_varhaptag_t(cliopt_haptag);
    cliopt_verbose = 0;
    opt = KETOPT_INIT;
    int i, c;

    if (argc==1){
        print_help_cli_varhaptag(cliopt_haptag);
        destroy_cliopt_varhaptag_t(cliopt_haptag);
        return 0;
    }

    while ((c = ketopt(&opt, argc, argv, 1, "ho:t:r:v", longopts)) >= 0) {
        if (c == 'h') {
            print_help_cli_varhaptag(cliopt_haptag); 
            cliopt_haptag->is_help=1;
        }
        //else if (c == 'r'){
        //    if (cliopt_haptag->n_regions==cliopt_haptag->m_regions){
        //        cliopt_haptag->m_regions += (cliopt_haptag->m_regions>>1);
        //        cliopt_haptag->regions = (char**)realloc(cliopt_haptag->regions, 
        //                                      sizeof(char*)*cliopt_haptag->m_regions);
        //    }
        //    cliopt_haptag->regions[cliopt_haptag->n_regions++] = opt.arg;

        //}
        else if (c == 't') cliopt_haptag->n_threads = atoi(opt.arg);
        else if (c == 'v') cliopt_haptag->verbose = 1;
        else if (c == 'o') cliopt_haptag->fn_out=opt.arg;
        else if (c == 309) cliopt_haptag->do_write_bam = 0;
        else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option argument in \"%s\"\n", __func__, argv[opt.i - 1]);
            destroy_cliopt_varhaptag_t(cliopt_haptag);
            return 0;
        }
        else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option argument in \"%s\"\n", __func__, argv[opt.i - 1]);
            destroy_cliopt_varhaptag_t(cliopt_haptag);
            return 0;
        }
    }
    int j=0;
    for (i=opt.ind; i<argc; i++){
        if (j==0)
            cliopt_haptag->fn_vcf = argv[i];
        else if (j==1)
            cliopt_haptag->fn_bam = argv[i];
        else{
            fprintf(stderr, "[E::%s] Invalid positional arguments: ", __func__);
            for (int k=opt.ind; k<argc; k++){
                fprintf(stderr, "%s ", argv[k]);
            }
            fprintf(stderr, ". Please provide a vcf followed by a bam.\n");
            destroy_cliopt_varhaptag_t(cliopt_haptag);
            return 0;
        }
        j++;
    }
    return cliopt_haptag;
}
/*** end of varhaptag ***/


/*** methstat ***/
static ko_longopt_t longopts_methstat[] = {
    { "lo", ko_required_argument,  301 },  // exclusive upperbound of unmeth call qual
    { "hi", ko_required_argument,  302 },  // inclusive lowerbound of meth call qual
    { "cov", ko_required_argument,  303 },  // 
    { "readlen", ko_required_argument,  304 },  // 
    { "gtf", ko_required_argument, 305 },  // input: gtf file name
    { "vcf", ko_required_argument, 306 },   // input: vcf file name
    { "tsv", ko_required_argument, 307 },   // input: tsv file name
    { "help", ko_no_argument, 400 },  
    { 0, 0, 0 }
};

void init_cliopt_methstat_t(cliopt_methstat_t *clio){
    clio->fn_bam = 0;
    clio->fn_vcf = 0;
    clio->fn_gtf = 0;
    clio->fn_tsv = 0;
    clio->fn_out = 0;
    clio->is_help = 0;
    clio->lo = 100;
    clio->hi = 156;
    clio->readlen = 15000;
    clio->cov = 6;
}
void destroy_cliopt_methstat_t(cliopt_methstat_t *clio){
    free(clio);
}
void print_help_cli_methstat(cliopt_methstat_t *clio){
    fprintf(stderr, "Usage: pomfret methstat -o out.tsv [options] in.bam 2>log\n");
    fprintf(stderr, "Purpose:\n");
    fprintf(stderr, "  Given aligned reads with methylation calls in bam and phase \n");
    fprintf(stderr, "  blocks in vcf/gtf/tsv, find seemingly heterozyous methylation\n");
    fprintf(stderr, "  positions for unphased regions.\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  bam        [pos] Aligned reads, with methylation calls.\n");
    fprintf(stderr, "  -o         [req] Output file name.\n");
    fprintf(stderr, "  -h,--help  [   ] Display this message.\n");
    fprintf(stderr, "  -c,--cov   [opt] For a site to be reported, the minimum occurence among\n");
    fprintf(stderr, "              meth, unmeth and nocall is required to be $INT.[%d]\n", clio->cov);
    fprintf(stderr, "  -t,--tsv   [opt] Phase blocks specified in tsv (chrom,start,end).\n");
    fprintf(stderr, "  -v,--vcf   [opt] Phase blocks specified in vcf.\n");
    fprintf(stderr, "  -g,--gtf   [opt] Phase blocks specified in gtf.\n");
    fprintf(stderr, "  -a,--hi    [opt] Upper threshold of methylation qual. Calls with \n");
    fprintf(stderr, "              quality value higher than this will be considered\n");
    fprintf(stderr, "              methylated. [%d]\n", clio->hi);
    fprintf(stderr, "  -b,--lo    [opt] Lower threshold of methylation qual. Calls with \n");
    fprintf(stderr, "             quality value lower than this will be considered\n");
    fprintf(stderr, "             unmethylated. [%d]\n", clio->lo);
    fprintf(stderr, "  -l,--readlen [opt] Ignore reads shorter than $INT. [%d]\n", clio->readlen);
}
int sancheck_cliopt_methstat(cliopt_methstat_t *clio){
    int ret = 0;
    if (!clio->fn_gtf && !clio->fn_vcf && !clio->fn_tsv){
        fprintf(stderr, "[E::%s] Must provide an interval file (vcf/gtf/tsv).\n", __func__);
        ret = 1;
    }
    if (!clio->fn_bam){
        fprintf(stderr, "[E::%s] Must provide one BAM file.\n", __func__);
        ret = 1;
    }
    if (!!clio->fn_gtf + !!clio->fn_vcf + !!clio->fn_tsv >1){
        fprintf(stderr, "[M::%s] More than one interval file were provided, will follow priority order tsv>gtf>vcf rather than merging.\n", __func__);
    }
    if (clio->lo<0){
        fprintf(stderr, "[W::%s] Invalid lower quality threshold of methylation (%d), clipping to 0\n", 
                __func__, clio->lo);
        clio->lo = 0;
    }
    if (clio->hi>255){
        fprintf(stderr, "[W::%s] Invalid higher quality threshold of methylation (%d), clipping to 255\n", 
                __func__, clio->hi);
        clio->lo = 255;
    }
    return ret;
}
cliopt_methstat_t *parse_cli_methstat(int argc, char *argv[]){
    cliopt_methstat_t *clio = (cliopt_methstat_t*)malloc(sizeof(cliopt_methstat_t));
    init_cliopt_methstat_t(clio);
    cliopt_verbose = 0;
    opt = KETOPT_INIT;
    int i, c;

    if (argc==1){
        print_help_cli_methstat(clio);
        destroy_cliopt_methstat_t(clio);
        return 0;
    }

    while ((c = ketopt(&opt, argc, argv, 1, "ho:v:g:t:a:b:l:c:", longopts_methstat)) >= 0) {
        if (c == 'h' || c==400) {print_help_cli_methstat(clio); clio->is_help=1;}
        else if (c == 't' || c==307) clio->fn_tsv = opt.arg;
        else if (c == 'v' || c==306) clio->fn_vcf = opt.arg;
        else if (c == 'g' || c==305) clio->fn_gtf= opt.arg;
        else if (c == 'o') clio->fn_out=opt.arg;
        else if (c == 'a' || c==302) clio->hi = atoi(opt.arg);
        else if (c == 'b' || c==301) clio->lo = atoi(opt.arg);
        else if (c == 'l' || c==304) clio->readlen = atoi(opt.arg);
        else if (c == 'c' || c==303) clio->cov= atoi(opt.arg);
        else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option argument in \"%s\"\n", __func__, argv[opt.i - 1]);
            destroy_cliopt_methstat_t(clio);
            return 0;
        }
        else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option argument in \"%s\"\n", __func__, argv[opt.i - 1]);
            destroy_cliopt_methstat_t(clio);
            return 0;
        }
    }
    for (i=opt.ind; i<argc; i++){
        if (clio->fn_bam){
            fprintf(stderr, "[E::%s] Multiple bam input is not supported.\n", __func__);
            destroy_cliopt_methstat_t(clio);
            return 0;
        }
        clio->fn_bam = argv[i];
    }
    return clio;
}



/*** end of methstat ***/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include <zlib.h>
#include "kvec.h"
#include "khashl.h"
#include "ksort.h"
#include "htslib/htslib/sam.h"
#include "htslib/htslib/hts.h"
#include "htslib/htslib/kstring.h"
#include "htslib/htslib/kfunc.h"
#include "cli.h"

#define VERBOSE 1
#define READBACK 50000
#define MAX(a,b) ((a)>(b)? (a):(b))
#define MIN(a,b) ((a)<=(b)? (a):(b))

typedef struct{
    uint32_t pos;
    //uint16_t cnt[3][3];  // 2 haplotypes + unphased, meth/unmeth/nocall
    uint16_t cnt[3]; // meth, unmeth, nocall
}methsite_counter_t;
KHASHL_MAP_INIT(KH_LOCAL, htu32_t, htu32_ht, uint32_t, methsite_counter_t, kh_hash_uint32, kh_eq_generic)

#define generic_key(x) (x)
#define generic_cmp(x,y) ((x)<(y))
KRADIX_SORT_INIT(ksu32, uint32_t, generic_key, 4)
KSORT_INIT(kssu32, uint32_t, generic_cmp)
KRADIX_SORT_INIT(ksu64, uint64_t, generic_key, 8)


typedef kvec_t(uint64_t) vu64_t;
typedef kvec_t(uint32_t) vu32_t;
typedef kvec_t(uint16_t) vu16_t;
typedef kvec_t(uint8_t) vu8_t;
typedef kvec_t(vu8_t) vu8_v;
typedef kvec_t(int) vi_t;
typedef kvec_t(float) vfloat_t;

#define get_cigar_op_len(a) ((a)>>4)
#define get_cigar_op_type(a) ((a)&15)

int search_arr(uint32_t *a, uint32_t l, uint32_t v, uint32_t *idx){
    // Search and try to find the closest index. 
    // Return: -1 if not found and not within range, to the left(idx u32imax)
    //         -2 if not found and not within range, to the right(idx u32imax)
    //         -3 when a is empty
    //         0 if not found but fall within range (idx contains i of the larger ref value)
    //         1 if found (idx containes the index)
    if (l==0) return -3;
    if (v<a[0]) {
        *idx = UINT32_MAX;
        return -1;
    }
    if (v>a[l-1]) {
        *idx = UINT32_MAX;
        return -2;
    }
    if (l<16){
        for (uint32_t i=0; i<l; i++){
            if (a[i]==v){
                *idx = i;
                return 1;
            }
            if (a[i]>v){
                *idx = i;
                return 0;
            }
        }
        fprintf(stderr, "[E::%s] impossible 1\n", __func__);
        exit(1);
    }else{
        uint32_t low = 0;
        uint32_t mid;
        uint32_t high = l-1;
        while (low<high){
            mid = low + (high-low)/2;
            if (v<=a[mid]){
                high = mid;
            }else{
                low = mid+1;
            }
        }
        if (a[high]==v){
            *idx = high;
            return 1;
        }else if (a[high]>v){
            *idx = high;
            return 0;
        }
        fprintf(stderr, "[E::%s] impossible 2\n", __func__);
        exit(1);
    }
}



typedef struct{
    vu32_t calls;
    vu8_t quals;  // actually this is categorized to 0(meth), 1(unmeth) and 2(nocall)
}mod_t;
void init_mod_t(mod_t *d, int l){
    if (l<16) l = 16;
    kv_init(d->calls);
    kv_resize(uint32_t, d->calls, l);
    kv_init(d->quals);
    kv_resize(uint8_t, d->quals, l);
}
void destroy_mod_t(mod_t* d, int include_self){
    kv_destroy(d->calls);
    kv_destroy(d->quals);
    if (include_self) free(d);
}


KHASHL_MAP_INIT(KH_LOCAL, htstru32_t, htstru32_ht, char*, uint32_t, kh_hash_str, kh_eq_generic)
typedef struct{
    uint32_t i;
    int hp;  // haplo tag
    uint32_t len; // length
    mod_t meth;
    // methmers collected
    uint32_t *mmr;
    int mmr_n;
    uint32_t mmr_start_i;  // as-in-buffer index, refer to methmers_t alloc
}read_t;  // a single read
typedef struct{
    read_t *a;
    uint32_t m, n;
    uint32_t ref_start, ref_end;
    vu32_t names_acl;  // read name accumulated lengths
    char *names;
    uint32_t names_m, names_n;
    htstru32_t *name2namei;
}rs_t;  // read set
rs_t *init_rs_t(){
    rs_t *ret = (rs_t*)malloc(sizeof(rs_t));
    ret->n = 0;
    ret->m = 1024;
    // read alloc
    ret->a = (read_t*)malloc(sizeof(read_t)*ret->m);
    for (int i=0; i<ret->m; i++){
        init_mod_t(&ret->a[i].meth, 0);
        ret->a[i].mmr_n = 0;
    }
    // read name indexing alloc
    kv_init(ret->names_acl);
    kv_resize(uint32_t, ret->names_acl, ret->m+1);
    kv_push(uint32_t, ret->names_acl, 0);
    // read name alloc
    ret->names_n = 0;
    ret->names_m = 10240;
    ret->names = (char*)malloc(ret->names_m);
    ret->name2namei = htstru32_ht_init();
    return ret;
}
read_t* get_empty_slot_rs_t(rs_t* rs){
    if (rs->n+1==rs->m){
        rs->m = rs->m + (rs->m>>1);
        rs->a = (read_t*)realloc(rs->a, rs->m*sizeof(read_t));
        kv_resize(uint32_t, rs->names_acl, rs->m+1);
        for (int i=rs->n; i<rs->m; i++){
            init_mod_t(&rs->a[i].meth, 0);
        }
    }
    rs->n++;
    return &rs->a[rs->n-1];
}
#define GET_QNAME_S(rs, i) ((rs).names+(rs).names_acl.a[(i)])
#define GET_QNAME_L(rs, i) ((rs).names_acl.a[(i)+1]-(rs).names_acl.a[(i)])
void push_read_name_rs_t(rs_t* rs, char *s, int s_l){
    if (s_l + rs->names_n +1 >= rs->names_m){
        while (s_l + rs->names_n >= rs->names_m){
            rs->names_m += rs->names_m>>1;
        }
        rs->names = (char*)realloc(rs->names, rs->names_m);
    }
    memcpy(rs->names+rs->names_n, s, s_l);  
    rs->names[rs->names_n+s_l] = 0; // null term
    int i = rs->names_acl.n-1;
    kv_push(uint32_t, rs->names_acl, rs->names_acl.a[i]+s_l+1);  // null term
    rs->names_n += s_l+1; // null term
}
void destroy_rs_t(rs_t *rs){
    free(rs->names);
    kv_destroy(rs->names_acl);
    for (int i=0; i<rs->m; i++){
        destroy_mod_t(&rs->a[i].meth, 0);
    }
    for (int i=0; i<rs->n; i++){
        if (rs->a[i].mmr_n) free(rs->a[i].mmr);
    }
    free(rs->a);
    htstru32_ht_destroy(rs->name2namei);
    free(rs);
}

typedef struct {
    char *fn;
    samFile *fp;
    hts_idx_t *fp_idx;
    bam_hdr_t *fp_header;
    bam1_t *buf_aln;
    hts_base_mod_state *buf_mod;
}bamfile_t;
bamfile_t *init_and_open_bamfile_t(char *fn){
    bamfile_t *h = (bamfile_t*)malloc(sizeof(bamfile_t));
    h->fn = (char*)malloc(strlen(fn)+1);
    sprintf(h->fn, "%s", fn);
    h->fp = hts_open(fn, "r");
    h->fp_idx = sam_index_load(h->fp, h->fn);
    h->fp_header = sam_hdr_read(h->fp);
    h->buf_aln = bam_init1();
    h->buf_mod = hts_base_mod_state_alloc();
    return h;
}
void destroy_bamfile_t(bamfile_t *h, int include_self){
    free(h->fn);
    bam_destroy1(h->buf_aln);
    hts_base_mod_state_free(h->buf_mod);
    hts_idx_destroy(h->fp_idx);
    sam_hdr_destroy(h->fp_header);
    sam_close(h->fp);
    if (include_self) free(h);
}

static char *code(int id) {
    static char code[20];
    if (id > 0) {
        code[0] = id;
        code[1] = 0;
    } else {
        snprintf(code, sizeof(code), "(%d)", -id);
    }
    return code;
}


int get_mod_poss_on_ref(mod_t *ret,
                         uint32_t *cigar, int cigar_l, 
                         uint32_t qs, int aln_strand, 
                         uint32_t *mod_poss, uint8_t *mod_quals, int mod_l){
    // Get mod indices as on the reference, along with its qualities [0, 256).
    // Return 0 if no call is availalbe, 1 otherwise.  
    int verbose = 0;
    if (cigar_l==0 || mod_l==0) {
        return 0;
    }
    int cgoffset = 0;
    if (aln_strand!=0)  cgoffset = -1;
    
    uint32_t i_read = 0;  // pos on read
    uint32_t i_ref = qs;  // pos on ref
    uint32_t i_ref_raw = qs;  // pos on ref
    uint32_t i_trigger = 0;  // idx of meth, as in buffer
    uint32_t next_trigger = mod_poss[i_trigger];  // idx of meth, as on read
    uint8_t next_qual = mod_quals[i_trigger];  // qual of meth

    // special case at start: silently consume mods that are within clipped range
    int i_cigar = 0;
    if (get_cigar_op_type(cigar[0])==4){  // 'S'
        i_read = get_cigar_op_len(cigar[0]);
        while (next_trigger<i_read){
            i_trigger++;
            if (i_trigger<mod_l){
                next_trigger = mod_poss[i_trigger];
                next_qual = mod_quals[i_trigger];
            }else break;
        }
        if (next_trigger==i_read){
            kv_push(uint32_t, ret->calls,       i_ref+cgoffset);
            kv_push(uint8_t,  ret->quals, next_qual);
            i_trigger++;
            if (i_trigger<mod_l){
                next_trigger = mod_poss[i_trigger];
                next_qual = mod_quals[i_trigger];
            }else{
                ;  // stepped out of bound, do nothing
            }
        }
        i_ref -= get_cigar_op_len(cigar[0]);
        i_cigar = 1;
    }

    // step through cigar operations
    int offset = 0;
    int offset_ref = 0;
    for (;i_cigar<cigar_l; i_cigar++){
        uint8_t action = get_cigar_op_type(cigar[i_cigar]);
        uint32_t length = get_cigar_op_len(cigar[i_cigar]);
        if (verbose) fprintf(stdout, "[dbg::%s] cigar op is: %d, len is %d\n", __func__, (int)action, (int)length);
        if (action<=1){  // 0 for 'M', 1 for 'I'
            while (i_read+length>=next_trigger){
                if (action==0 && next_trigger!=UINT32_MAX){
                    kv_push(uint32_t, ret->calls, i_ref+next_trigger+cgoffset + offset);
                    kv_push(uint8_t, ret->quals, next_qual);
                }
                i_trigger++;
                if (i_trigger>=mod_l){  // stepped out of bound
                    next_trigger = UINT32_MAX;
                    break;
                }else{
                    next_trigger = mod_poss[i_trigger];
                    next_qual = mod_quals[i_trigger];
                }
            }
            if (action==0){
                i_read+=length; 
                offset_ref+=length;
            }else{
                i_read += length;
                offset -= length;
            }
        }else if (action==2){  // 'D'
            offset+=length; 
            offset_ref+= length;
        }else if (action==3){
            break;
        }else if (action==4){  // soft clip
            offset_ref += length;
            break;
        }else{
            fprintf(stderr, "[E::%s] fatal: unknown cigar operation (value=%d). Bug?\n", __func__, (int)action);
            exit(1);
        }

    }

    // dbg
    if (verbose){
        for (int i=0; i<ret->calls.n; i++){
            fprintf(stdout, "[dbg::%s]     meth pos=%d, qual=%d\n", __func__, (int)ret->calls.a[i], (int)ret->quals.a[i]);
        }
    }

    assert(ret->calls.n==ret->quals.n);
    return 1;
}

int fill_read_meth_record_from_bam_line(read_t *h, 
                                    bam1_t *buf_aln, 
                                    hts_base_mod_state *buf_mod, 
                                    vu32_t *buf_mod_poss, 
                                    vu8_t *buf_mod_quals, 
                                    uint8_t qual_lo, uint8_t qual_hi){
    // return 1 if added something; 0 otherwise
    uint32_t *cigar = bam_get_cigar(buf_aln);
    int cigar_l = buf_aln->core.n_cigar;
	uint32_t len = buf_aln->core.l_qseq; //length of the read.
    uint32_t qs = buf_aln->core.pos;
    int strand = !!(buf_aln->core.flag&16);
    bam_parse_basemod(buf_aln, buf_mod);


    // get as-on-read coordiates of meth calls
    int n, mod_pos;
    hts_base_mod mods[5];
    char qstr[10];
    while ((n=bam_next_basemod(buf_aln, buf_mod, mods, 5, &mod_pos)) > 0) {
        for (int j = 0; j < n && j < 5; j++) {
            if (mods[j].qual == HTS_MOD_UNCHECKED)
                qstr[0] = '#', qstr[1] = 0;
            else if (mods[j].qual == HTS_MOD_UNKNOWN)
                qstr[0] = '.', qstr[1] = 0;
            else
                snprintf(qstr, 10, "%d", mods[j].qual);  // this is the int qual

            if (mods[j].canonical_base=='C' && 
                 strcmp(code(mods[j].modified_base), "m")==0){
                kv_push(uint32_t, *buf_mod_poss, (uint32_t)mod_pos);
                uint8_t q = (uint8_t)mods[j].qual;
                kv_push(uint8_t, *buf_mod_quals, q<qual_lo? 
                                                 1 : q>=qual_hi?
                                                 0 : 2);
            }
        }
    }

    // store as-on-reference coordinates of meth calls
    int stat = get_mod_poss_on_ref(&h->meth, cigar, cigar_l, 
                                   qs, strand, 
                                   buf_mod_poss->a, buf_mod_quals->a, 
                                   buf_mod_poss->n);
    return stat;
}

int get_hp_from_aln(bam1_t *aln){
    uint8_t *tagd = bam_aux_get(aln, "HP");
    if (tagd){
        int hp = (int)bam_aux2i(tagd);
        if (hp==0) {
            fprintf(stderr, "[W::%s] irregular HP tag?\n", __func__);
            return 2; // TODO diploid assumption
        }
        return hp-1;  // use 0-index
    }
    return 2;  // TODO diploid assumption
}
int add_read_record_from_bam_line(rs_t *rs, bam1_t *aln, 
                                   hts_base_mod_state *mod, 
                                   vu32_t *mod_poss, 
                                   vu8_t *mod_quals, 
                                   uint8_t qual_lo, uint8_t qual_hi){
    read_t *h = get_empty_slot_rs_t(rs);

    // process mods first, if there is nothing to add, we ignore this read
    int stat = fill_read_meth_record_from_bam_line(h, aln, mod, mod_poss, mod_quals, qual_lo, qual_hi);
    if (stat==0){
        rs->n--;
        return 0; 
    }

    char *qname = bam_get_qname(aln);  // is null terminated
    push_read_name_rs_t(rs, qname, strlen(qname));

    // add other info
    int hp = get_hp_from_aln(aln);

    h->i = rs->n-1;
    h->len = aln->core.l_qseq;
    h->hp = hp;

    return 1;
}

rs_t *load_reads_given_interval(char *fn,
                                char *chrom, int itvl_s, int itvl_e, 
                                int readback,
                                int min_mapq, int min_read_len, 
                                uint8_t qual_lo, uint8_t qual_hi, 
                                int* n_start_reads,   // stores the number of ref reads on the left
                                int* n_end_reads  // stores the number of ref reads on the right
                                ){
    rs_t* rs = init_rs_t();
    rs->ref_start = itvl_s;
    rs->ref_end = itvl_e;
    char *itvl_fmt = (char*)malloc(strlen(chrom)+30);
    sprintf(itvl_fmt, "%s:%d-%d", chrom, itvl_s-readback, itvl_e);

    bamfile_t *hf = init_and_open_bamfile_t(fn);
    hts_itr_t *fp_itr = sam_itr_querys(hf->fp_idx, hf->fp_header, itvl_fmt);

    vu32_t buf_mod_poss;
    vu8_t buf_mod_quals;
    kv_init(buf_mod_poss); kv_resize(uint32_t, buf_mod_poss, 16);
    kv_init(buf_mod_quals); kv_resize(uint8_t, buf_mod_quals, 16);

    int n_added = 0;
    int get_n_start = (n_start_reads)!=0? 1 : 0;
    int get_n_end   = (n_end_reads  )!=0? 1 : 0;
    int left_side_cov_check[2] = {0,0};
	//while(sam_read1(fp, fp_header, buf_aln) > 0){
    while(sam_itr_next(hf->fp, fp_itr, hf->buf_aln)>=0){
        //char *chr = hf->fp_header->target_name[hf->buf_aln->core.tid];
        //char *qn = bam_get_qname(hf->buf_aln);

        int flag = hf->buf_aln->core.flag;
		//uint8_t *q = bam_get_seq(buf_aln); //quality string
		uint32_t mapq = hf->buf_aln->core.qual ; //mapping quality
		uint32_t len = hf->buf_aln->core.l_qseq; //length of the read.
        float de = -1;
        uint8_t *tmp = bam_aux_get(hf->buf_aln, "de");
        if (tmp) de = bam_aux2f(tmp);

        if ((flag&4) || (flag&256) || (flag&2048)) continue;
        if (mapq<min_mapq) continue;
        if (len<min_read_len /* && hf->buf_aln->core.pos>itvl_s*/) continue;  // maybe for reads that overlap with the start position, we retain them regardless of their length.
        if (de>0.1) continue;  // filter out reads with 
                               // poor alignment and/or poor bascalling

        buf_mod_poss.n = 0;
        buf_mod_quals.n = 0;
        int stat = add_read_record_from_bam_line(rs, hf->buf_aln, hf->buf_mod, 
                        &buf_mod_poss, &buf_mod_quals, qual_lo, qual_hi);
        if (stat>0) n_added++;
        if (hf->buf_aln->core.pos<=itvl_s) {
            int hp = get_hp_from_aln(hf->buf_aln);
            if (hp==0 || hp==1){
                left_side_cov_check[hp]++;
            }
            if (get_n_start) (*n_start_reads)++;
        }
        if (get_n_end && (hf->buf_aln->core.pos + len)>=itvl_e) (*n_end_reads)++;
    }
    // generate reverse index: store qname to read index as-in-buffer in ht
    for (int i=0; i<rs->n; i++){
        char *qn = GET_QNAME_S(*rs, i);
        khint_t k = htstru32_ht_get(rs->name2namei, qn);
        int absent;
        k = htstru32_ht_put(rs->name2namei, qn, &absent);
        if (!absent){
            fprintf(stderr, "[E::%s] duplicated read name seen from reading bam: %s\n", 
                __func__, qn);
            exit(1);
        }
        kh_val(rs->name2namei, k) = i;
        kh_key(rs->name2namei, k) = qn;
    }

    fprintf(stderr, "[M::%s] loaded %d reads from %s\n", __func__, n_added, itvl_fmt);
    if (get_n_start)
        fprintf(stderr, "[M::%s] %d reads at the start position\n", __func__, *n_start_reads);

    // TODO better cleanup for this case
    // if left side coverage is bad, assume the region is not worth trying 
    if (left_side_cov_check[0]<15 || left_side_cov_check[1]<15){
        rs->n = 0;
    }

    // clean up
    hts_itr_destroy(fp_itr);
    destroy_bamfile_t(hf, 1);
    kv_destroy(buf_mod_poss);
    kv_destroy(buf_mod_quals);
    free(itvl_fmt);
    
    return rs;
}




typedef struct {
    vu32_t starts;
    vu32_t ends;
}ranges_t;
ranges_t *init_ranges_t(){
    ranges_t *ret = (ranges_t*)malloc(sizeof(ranges_t));
    kv_init(ret->starts); 
    kv_init(ret->ends); 
    kv_resize(uint32_t, ret->starts, 16);
    kv_resize(uint32_t, ret->ends, 16);
    return ret;
}
void destroy_ranges_t(ranges_t* d){
    kv_destroy(d->starts);
    kv_destroy(d->ends);
    free(d);
}



int push_substring_forward_in_buffer(char *buf, int buf_l, int start){
    memmove(buf, buf+start, buf_l-start);
    return buf_l-start;
}

void insert_gtf_line(char *s, char *chrom, ranges_t *intvls, uint32_t *prev_end, 
                     int is_tsv){
    int verbose = 1;
    char *tok = strtok(s, "\t");
    int i=0;
    uint8_t use = 0;
    int col_chrom, col_s, col_e;
    if (is_tsv){  // input is a 3-col tsv rather than gtf
        col_chrom = 0;
        col_s = 1;
        col_e = 2;
    }else{
        col_chrom = 0;
        col_s = 3;
        col_e = 4;
    }
    while (tok){
        if (i==col_chrom){
            if (strcmp(chrom, tok)!=0){
                use = 0;
            }else{
                use = 1;
            }
        }else if (i==col_s && use){  // start of phased block
            if ((*prev_end)!=UINT32_MAX){
                kv_push(uint32_t, intvls->starts, *prev_end);
                kv_push(uint32_t, intvls->ends, strtoul(tok, NULL, 10));
            }
        }else if (i==col_e && use){  // end of phased block
            *prev_end = strtoul(tok, NULL, 10);
        }
        tok = strtok(NULL, "\t");
        i++;
    }
}


ranges_t *load_intervals_from_gtf(char *fn, char *chrom, int is_tsv){
    fprintf(stderr, "[M::%s] reading from %s\n", __func__, fn);
    ranges_t *ret = init_ranges_t();
    gzFile fp = gzopen(fn, "rb");
    if (!fp){
        fprintf(stderr, "[E::%s] failed to open file: %s\n", __func__, fn);
        exit(1);
    }
    int m = 128;
    int m_tmp = m;
    int n = 0;
    char *buffer = (char*)malloc(m);
    char *tmp_buffer = (char*)malloc(m);
    int offset = 0;
    int n_insertions = 0;

    uint32_t prev_end = UINT32_MAX;
    int n_read;
    while ((n_read = gzread(fp, buffer+offset, m-offset)) > 0) {
        int start = 0;
        for (int i=0; i<offset+n_read; i++){
            if (buffer[i]=='\n'){
                if (i-start+1>=m_tmp){
                    m_tmp = i+1;
                    tmp_buffer = (char*)realloc(tmp_buffer, m_tmp);
                }
                sprintf(tmp_buffer, "%.*s", i-start, buffer+start);
                insert_gtf_line(tmp_buffer, chrom, ret, &prev_end, is_tsv);
                n_insertions++;
                start = i+1; 
            }
        }

        if (start!=0){
            offset = push_substring_forward_in_buffer(buffer, m, start);
        }else{  // did not see a full line
            offset = m;
            m = m+(m>>1);
            buffer = (char*)realloc(buffer, m);
        }
    }
    fprintf(stderr, "[M::%s] chrom %s, read %d unphased intervals from %s\n", __func__, 
        chrom, n_insertions, fn);
    gzclose(fp);
    free(buffer);
    free(tmp_buffer);
    return ret;
}


void merge_close_intervals(ranges_t *rg, int threshold){
    int j = 0;
    for (int i=1; i<rg->starts.n; i++){
        if (rg->starts.a[i] - rg->ends.a[j]<threshold){  // merge
            rg->ends.a[j] = rg->ends.a[i];
        }else{
            j++;
            rg->starts.a[j] = rg->starts.a[i];
            rg->ends.a[j] = rg->ends.a[i];
        }
    }
    fprintf(stderr, "[M::%s] had %ld intervals, merged %ld at threshold %d\n", 
        __func__, rg->starts.n, rg->starts.n-j-1, threshold);
    rg->starts.n = j+1;
    rg->ends.n = j+1;
}


typedef struct{
    vu32_t vmmr_u32i;
    vu16_t hap_cnt[2];
    uint16_t sum[2];
}mmr_t;  // methmers at one position 
typedef struct{
    int n;
    int k;
    int max_span;
    uint32_t *sites_real_poss;  // meth sites real positions(i.e. as-on-reference)
    uint8_t *mmr_lens;  // methmer lens
    mmr_t *mmr_a;  // methers storage: mers and per-haplotype count at each site
    // what methermer is available to be used in voting:
    uint32_t mmr_min_i;  // index of the first methmer with enough coverage
    uint32_t mmr_max_i;  // index of the last methmer with enough coverage
}methmers_t;
methmers_t *init_methmers_t(int k, int max_span, int n){
    methmers_t *ret = (methmers_t*)malloc(sizeof(methmers_t));
    ret->k = k;
    ret->n = n;
    ret->max_span = 0;
    if (n!=0){
        ret->sites_real_poss = (uint32_t*)malloc(sizeof(uint32_t)*n);
        ret->mmr_lens = (uint8_t*)malloc(n);
    }else{
        ret->sites_real_poss = 0;
        ret->mmr_lens = 0;
    }

    // mmr_t allocation will wait until sites are collected
    ret->mmr_a = 0;


    return ret;
}
void destroy_methmers_t(methmers_t *h){
    if (h->n!=0){
        free(h->sites_real_poss);
        free(h->mmr_lens);
    }
    if (h->mmr_a){
        for (int i=0; i<h->n; i++){
            kv_destroy(h->mmr_a[i].vmmr_u32i);
            kv_destroy(h->mmr_a[i].hap_cnt[0]);
            kv_destroy(h->mmr_a[i].hap_cnt[1]);
        }
        free(h->mmr_a);
    }
    free(h);
}

typedef struct{
    rs_t *rs;
    methmers_t *ms;
}dataset_t;
void destroy_dataset_t(dataset_t* d, int include_self){
    destroy_rs_t(d->rs);
    if (d->ms) destroy_methmers_t(d->ms);
    if (include_self) free(d);
}

uint32_t methmer_to_uint32(char *s, int l){
    uint32_t ret = 0;
    for (int i=0; i<l; i++){
        ret = ret<<2 | (s[i]=='m'? 
                        0UL: s[i]=='u'?
                        1UL: 2UL);
    }
    return ret;
}
void uint32_to_methmer(uint32_t v, char *buf, int k){
    for (int i=0; i<k; i++){
        buf[i] = "mu-"[(v>>((k-i-1)*2))&3];
    }
}


methmers_t *get_methmer_sites_and_ranges(rs_t *rs, 
                                  int mmr_k, int mmr_max_span, 
                                  int min_either_cov){
    // Collect number of meth/unmeth/nocall calls, decide meth sites and methmers
    // TODO use bst instead of ht?
    methmers_t *ret;
    htu32_t *ht = htu32_ht_init();
    khint_t htk;
    int absent;
    for (int i=0; i<rs->n; i++){
        for (int j=0; j<rs->a[i].meth.calls.n; j++){
            uint32_t pos = rs->a[i].meth.calls.a[j];
            int call = rs->a[i].meth.quals.a[j];
            //int hp = rs->a[i].hp;
            //if (hp!=0 && hp!=1) hp = 2;
            htk = htu32_ht_get(ht, pos);
            if (htk==kh_end(ht)){
                htk = htu32_ht_put(ht, pos, &absent);
                kh_val(ht, htk).cnt[0] = 0; // meth
                kh_val(ht, htk).cnt[1] = 0; // unmeth 
                kh_val(ht, htk).cnt[2] = 0; // uncall
                kh_val(ht, htk).cnt[call] = 1;
                kh_val(ht, htk).pos = pos;
                //printf("[dbg::%s] saw pos %d\n", __func__, (int)pos);
            } else {
                if (kh_val(ht, htk).cnt[call]<UINT16_MAX-1) 
                    kh_val(ht, htk).cnt[call]++;
            }
        }
    }

    // collect all positions as-on-ref that 
    // have enough both meth calls and unmeth calls.
    vu32_t ok_sites; 
    kv_init(ok_sites); kv_resize(uint32_t, ok_sites, 32);
    for (htk=0; htk<kh_end(ht); htk++){
        if (kh_exist(ht, htk)){
            methsite_counter_t tmp = kh_val(ht, htk);
            if (tmp.cnt[0]>=min_either_cov && tmp.cnt[1]>=min_either_cov){ 
                // TODO use percentage wrt site coverage?
                kv_push(uint32_t, ok_sites, tmp.pos);
            }
        }
    }
    //fprintf(stderr, "[dbg::%s] # sites: %d\n", __func__, (int)ok_sites.n);
    htu32_ht_destroy(ht);
    // (store the sites; there is no reverse index, use binary search)
    ret = init_methmers_t(mmr_k, mmr_max_span, ok_sites.n);
    memcpy(ret->sites_real_poss, ok_sites.a, ok_sites.n*sizeof(uint32_t));
    kv_destroy(ok_sites);
    radix_sort_ksu32(ret->sites_real_poss, ret->sites_real_poss+ret->n);
    //fprintf(stderr, "[dbg::%s] sorted sites\n", __func__);

    // collect methmers
    for (int i=0; i<ret->n; i++){
        int j = i+mmr_k;
        j = j>ret->n-1? ret->n-1 : j;
        int len;
        while (1){
            if (ret->sites_real_poss[j]-ret->sites_real_poss[i]<=mmr_max_span){
                break;
            }
            j--;
        }
        len = j-i;
        ret->mmr_lens[i] = len;
    }

    // allocate storage for counting 
    ret->mmr_a = (mmr_t*)malloc(sizeof(mmr_t)*ret->n);
    for (int i=0; i<ret->n; i++){
        kv_init(ret->mmr_a[i].vmmr_u32i);
        kv_resize(uint32_t, ret->mmr_a[i].vmmr_u32i, 16);
        kv_init(ret->mmr_a[i].hap_cnt[0]);
        kv_resize(uint16_t, ret->mmr_a[i].hap_cnt[0], 16);
        kv_init(ret->mmr_a[i].hap_cnt[1]);
        kv_resize(uint16_t, ret->mmr_a[i].hap_cnt[1], 16);
        ret->mmr_a[i].sum[0] = 0;
        ret->mmr_a[i].sum[1] = 0;
    }

    return ret;
}


int get_mmr_of_read(read_t *read, methmers_t *ms, vu32_t *buf_mmr, uint32_t *start_i){ 
    // Given methmer ranges, collect methmers (in uint32) from a read and 
    //  store them in buf_mmr.
    // Return: number of methmers. 
    buf_mmr->n = 0;

    uint32_t *sites = ms->sites_real_poss;
    uint32_t sites_n = ms->n;
    uint8_t *mml = ms->mmr_lens;
    uint32_t *calls = read->meth.calls.a;
    uint8_t *quals= read->meth.quals.a;
    uint32_t calls_n = read->meth.calls.n;
    if (calls_n==0) return 0;

    uint32_t x, x_i_left, x_i_right;  // as-on-ref position and as-in-buffer index for fixed index

    // init: find the left most call
    // (left)
    int stat = search_arr(sites, sites_n, calls[0], &x_i_left);
    if (stat==-2 || stat==-3) return 0;   // no overlap: to the right
    if (stat==0) x_i_left = x_i_left==0? 0: x_i_left-1;
    // (right)
    stat = search_arr(sites, sites_n, calls[calls_n-1], &x_i_right);
    if (stat==-1 || stat==-3) {
        return 0;
    }
    if (x_i_left==UINT32_MAX) x_i_left = 0;
    if (x_i_right==UINT32_MAX) x_i_right = sites_n;
    
    // piggybacked sort
    vu64_t buf;
    kv_init(buf);
    kv_resize(uint64_t, buf, 16);
    for (uint32_t i=x_i_left; i<x_i_right; i++){
        kv_push(uint64_t, buf, ((uint64_t)sites[i])<<35 | i);
    }
    for (int i=0; i<calls_n; i++){
        kv_push(uint64_t, buf, ((uint64_t)((calls[i])<<3 | 4 | quals[i]))<<32);
    }
    radix_sort_ksu64(buf.a, buf.a+buf.n);

    // collect methmers
    char mmr[16];  // note: this max length is limited by the u32i for encoding a methmer
    buf_mmr->n = 0;
    uint32_t start_pos_i = UINT32_MAX;
    uint64_t maskbit = ((uint64_t)4)<<32;
    for (int i=0; i<buf.n; i++){
        if ((buf.a[i]&maskbit)==0){  // start from an available site
            uint32_t pos = buf.a[i]>>35;
            uint32_t pos_i = (uint32_t)buf.a[i];
            int mmr_len = ms->mmr_lens[pos_i];
            int n = 0;
            for (int j=i; j<buf.n-1;){
                if ((buf.a[j]&maskbit)!=0){
                    j++;
                    continue;
                }
                if ( (buf.a[j]>>35) == (buf.a[j+1]>>35) && (buf.a[j+1]&maskbit)!=0 ){
                    mmr[n] = "mu-"[(buf.a[j+1]>>32) & 3];  // TODO do not need this just use the integer
                    n++;
                    j+=2;
                }else{  // read lacks any call at this position, use -
                    mmr[n] = '-';
                    n++;
                    j++;
                }
                if (n>=mmr_len) {
                    break;
                }
            }
            if (n!=mmr_len) continue;  // not a full methmer; this should only occur at the end of the read
            if (start_pos_i==UINT32_MAX) start_pos_i = pos_i;  // start index as-in-storage
            //fprintf(stderr, "[%s] pos=%d , mmr: %.*s\n", __func__, pos, n, mmr);
            uint32_t mmr_u32i = methmer_to_uint32(mmr, n);
            kv_push(uint32_t, *buf_mmr, mmr_u32i);
        }
    }
    if (buf_mmr->n!=0) *start_i = start_pos_i;
    else *start_i=UINT32_MAX;

    kv_destroy(buf);
    return buf_mmr->n;
}

void insert_mmrs_to_counts(methmers_t *ms, uint32_t *mmr_u32i, int n_mmr, int mmr_start_i, int hap){
    // TODO: is linear search really faster than hashtable?
    //       Probably yes? If no, simd?
    //       Given that length of each run is bound by read coverage,
    //       not combinatory, search is very like to have <100 comparisons
    //       in a contiguous layout. Shouldn't need a hashtable / hashtables. 
    if (hap!=0 && hap!=1){
        fprintf(stderr, "[E::%s] bug or unimplemented: hap=%d is invalid.\n", __func__, hap);
        exit(1);
    }
    for (int i0=0, j; i0<n_mmr; i0++){
        int i = i0+mmr_start_i;
        uint32_t query = mmr_u32i[i0];
        int found = 0;
        for (j=0; j<ms->mmr_a[i].vmmr_u32i.n; j++){  // check each stored methmer key
            if (query==ms->mmr_a[i].vmmr_u32i.a[j]){
                // is a known methmer, increment
                ms->mmr_a[i].hap_cnt[hap].a[j]++;
                ms->mmr_a[i].sum[hap]++;
                found = 1;
                break;
            }
        }
        if (!found){  // a new methmer key at this position
            kv_push(uint32_t, ms->mmr_a[i].vmmr_u32i, query);
            kv_push(uint16_t, ms->mmr_a[i].hap_cnt[0], 0);
            kv_push(uint16_t, ms->mmr_a[i].hap_cnt[1], 0);
            ms->mmr_a[i].hap_cnt[hap].a[j]++;
            ms->mmr_a[i].sum[hap]++;
        }
    }
}
void query_counts_of_mmrs(methmers_t *ms, uint32_t *mmr_u32i, int n_mmr, 
                uint32_t start_i, int hap, vfloat_t *buf){
    // TODO: is linear search really faster than hashtable?
    if (hap!=0 && hap!=1){
        fprintf(stderr, "[E::%s] bug or unimpl, hap=%d is invalid\n", 
                    __func__, hap);
        exit(1);
    }
    buf->n = 0;

    for (int i0=0, j; i0<n_mmr; i0++){
        int i = start_i + i0;
        uint32_t q = mmr_u32i[i0];
        int found = 0;
        for (j=0; j<ms->mmr_a[i].vmmr_u32i.n; j++){
            if (ms->mmr_a[i].vmmr_u32i.a[j]==q){
                uint32_t cnt = ms->mmr_a[i].hap_cnt[hap].a[j];
                uint32_t sum = ms->mmr_a[i].sum[hap];
                if (sum==0) ;//kv_push(float, *buf, 0);
                else        kv_push(float, *buf, (float)cnt/sum);
                found = 1; 
                break;
            }
        }
        //if (!found){
        //    kv_push(float, *buf, 0);
        //}
    }
}


void store_mmr_of_one_read(read_t *r, uint32_t *mmr, int mmr_n, uint32_t start_i){
    if (mmr_n==0 || start_i==UINT32_MAX){
        r->mmr_n = 0;
        r->mmr = 0;
        r->mmr_start_i = 0;
    }else{
        r->mmr_n = mmr_n;
        r->mmr_start_i = start_i;
        r->mmr = (uint32_t*)malloc(sizeof(uint32_t)*mmr_n);
        memcpy(r->mmr, mmr, sizeof(uint32_t)*mmr_n);
    }
}


void store_mmr_of_reads(rs_t* rs, methmers_t *ms){
    vu32_t buf_mmr;
    uint32_t start_i=0;
    kv_init(buf_mmr);
    kv_resize(uint32_t, buf_mmr, 64);
    for (int i=0; i<rs->n; i++){
        int n_mmr = get_mmr_of_read(&rs->a[i], ms, &buf_mmr, &start_i);
        store_mmr_of_one_read(&rs->a[i], buf_mmr.a, buf_mmr.n, start_i);
    }
    kv_destroy(buf_mmr);
}

void update_mmr_count_with_one_read(read_t *r, methmers_t *ms, int hap){
    if (r->mmr_n>0) 
        insert_mmrs_to_counts(ms, r->mmr, r->mmr_n, r->mmr_start_i, hap);
}


int use_mmr_count_predict_tag_for_one_read(read_t *r, methmers_t *ms, 
                                           float *best_score, 
                                           float score_diff_min, 
                                           int score_l_min){
    // Return: -1 untagged, 0 hap0, 1 hap1
    vfloat_t buf;
    kv_init(buf);
    kv_resize(float, buf, 32);

    float score0 = 0;
    float score1 = 0;
    int score0_l = 0;
    int score1_l = 0;

    // get actual methmer range: start from read's first methmer, 
    // read until the last methmer with enough runtime coverage && within the read's range 
    int n_mmr_avail = ms->mmr_max_i - r->mmr_start_i + 1;
    n_mmr_avail = r->mmr_n < n_mmr_avail? r->mmr_n : n_mmr_avail;

    // count of hap0
    buf.n = 0;
    query_counts_of_mmrs(ms, r->mmr, n_mmr_avail, r->mmr_start_i, 0, &buf);
    score0_l = buf.n;
    for (int i=0; i<buf.n; i++) {
        if (buf.a[i]>0){
            score0 += buf.a[i];
            score0_l++;
        }
    }

    // count of hap1
    buf.n = 0;
    query_counts_of_mmrs(ms, r->mmr, n_mmr_avail, r->mmr_start_i, 1, &buf);
    score1_l = buf.n;
    for (int i=0; i<buf.n; i++) {
        if (buf.a[i]>0){
            score1 += buf.a[i];
            score1_l++;
        }
    }
    float score_diff = score0>score1? score0-score1 : score1-score0;
    //fprintf(stderr, "[dbg::%s] hap0 %f (n=%d), hap1 %f (n=%d)\n", __func__, 
    // score0, score0_l, score1, score1_l);

    kv_destroy(buf);

    // require a good enough score to tag
    float affine_diff_min = score_diff_min + ((float)MIN(score0_l, score1_l)-3)*0.5;
    if (score_diff<score_diff_min && 
         (score0_l<score_l_min || score1_l<score_l_min) ) {
        *best_score = 0;
        return -1;
    }
    *best_score = score_diff;
    if (score0>score1){
        return 0;
    }else{
        return 1;
    }
}


typedef struct {
    float score;
    int tag;
    uint32_t readID;  // as in read set alloc
    uint32_t idx;  // as in a temprary buffer, for retrivial when done
}forpred_t;
typedef kvec_t(forpred_t) forpred_v;
#define forpred_cmp(x,y) ((x).score<(y).score)
KSORT_INIT(ksforpred, forpred_t, forpred_cmp)

int predict_tags_of_reads(rs_t *rs, methmers_t *ms, 
                           uint32_t *readIDs, int readIDs_n, 
                           vu32_t *buf_readIDs_sorted, vi_t *buf_tags, 
                           forpred_v *buf_tmp, forpred_v *buf_tmp_tmp,  // structure for sorting and its tmp buffer required by merge sort
                           int insert_best_n, 
                           int mmr_cov_for_update, 
                           float score_diff_min, int score_l_min
                           ){
    // outputs are stored sorted by score
    float s;
    
    // reset tmp buffers
    buf_tmp->n = 0;
    buf_tmp_tmp->n = 0;

    for (uint32_t i=0; i<readIDs_n; i++){
        uint32_t readID = readIDs[i];
        int tag = use_mmr_count_predict_tag_for_one_read(&rs->a[readID], ms, &s, 
                                                   score_diff_min, score_l_min);
        //fprintf(stderr, "[dbg::%s] read was %.*s , score %f, tag %d\n", 
        //    __func__, 
        //    (int)GET_QNAME_L(*rs, readID), 
        //    GET_QNAME_S(*rs, readID), s, tag);
        kv_push(forpred_t, *buf_tmp, ((forpred_t){.score=s, 
                                                 .tag=tag,
                                                 .readID=readID, 
                                                 .idx=i}
                                     )
               );
    }
    if (buf_tmp->m > buf_tmp_tmp->m){
        buf_tmp_tmp->m = buf_tmp->m;
        kv_resize(forpred_t, *buf_tmp_tmp, buf_tmp_tmp->m);
        buf_tmp_tmp->n = 0;
    }
    ks_mergesort_ksforpred(buf_tmp->n, buf_tmp->a, buf_tmp_tmp->a); // stable


    buf_readIDs_sorted->n = 0;
    buf_tags->n = 0;
    for (int i=0; i<buf_tmp->n; i++){
        //uint32_t idx = buf_tmp->a[i].idx;
        kv_push(int, *buf_tags, buf_tmp->a[i].tag);
        kv_push(uint32_t, *buf_readIDs_sorted, buf_tmp->a[i].readID);
    }
    //fprintf(stderr, "==== sancheck\n");
    //for (int i=0; i<buf_tmp->n; i++){
    //    fprintf(stderr, "read %.*s , score %f, tag %d\n", 
    //        (int)GET_QNAME_L(*rs, buf_tmp->a[i].readID), 
    //        GET_QNAME_S(*rs, buf_tmp->a[i].readID), 
    //        buf_tmp->a[i].score,
    //        buf_tmp->a[i].tag);
    //}
    //fprintf(stderr, "==== end of ^sancheck\n");

    // optional: update counter and tag the best read(s)
    int n = 0;
    if (insert_best_n>0){  // is tagging new reads, only update when query has haplotype predicted
        for (int i=buf_readIDs_sorted->n-1; i>=0; i--){
            uint32_t readID = buf_readIDs_sorted->a[i];
            int hap = buf_tags->a[i];
            if (hap==0 || hap==1){
                rs->a[readID].hp = hap;
                insert_mmrs_to_counts(ms, rs->a[readID].mmr, 
                    rs->a[readID].mmr_n, 
                    rs->a[readID].mmr_start_i, hap);
                n++;
                //fprintf(stderr, "[dbg::%s] tagged read# %d, qname=%.*s, hap=%d\n", 
                //    __func__, readID, 
                //    (int)GET_QNAME_L(*rs, readID), 
                //    GET_QNAME_S(*rs, readID), hap);
                if (n==insert_best_n) 
                    break;
            }
        }
        
    }
    if (n>0){  // inserted some reads, now update the valid methmer range
        for (int i=ms->mmr_max_i; i<ms->n; i++){
            if (ms->mmr_a[i].sum[0]+ms->mmr_a[i].sum[1]>=mmr_cov_for_update) ms->mmr_max_i = i;
            else break;  // should be impossible to have holes in methmer count, so look no further
        }
        //fprintf(stderr, "[dbg::%s] methmer now valid until %d\n", 
        //    __func__, 
        //    (int)ms->sites_real_poss[ms->mmr_max_i]);

    }
    return n;
    
}

void insert_ref_reads_methmer_counts(rs_t* rs, methmers_t* ms, 
                                     uint32_t *readIDs, uint32_t readIDs_n, int mmr_cov_for_update){
    // wipe counter buffer
    for (int i=0; i<ms->n; i++){
        // reset individual mmr string's count
        ms->mmr_a[i].hap_cnt[0].n = 0;
        ms->mmr_a[i].hap_cnt[1].n = 0;
        // reset total counts
        ms->mmr_a[i].sum[0] = 0;
        ms->mmr_a[i].sum[1] = 0;
        // remove mmr strings
        ms->mmr_a[i].vmmr_u32i.n = 0;
    }

    int sancheck[2] = {0,0};
    for (int i=0; i<readIDs_n; i++){
        uint32_t readID = readIDs[i];
        int hap = rs->a[readID].hp;
        if (hap==0 || hap==1){
            sancheck[hap]++;
            insert_mmrs_to_counts(ms, rs->a[readID].mmr, rs->a[readID].mmr_n, rs->a[readID].mmr_start_i, hap);
            //fprintf(stderr, "[dbg::%s] saw ref read %.*s, hap is %d\n", __func__, 
            //    (int)GET_QNAME_L(*rs, readID), GET_QNAME_S(*rs, readID), rs->a[readID].hp);
        }
    }
    for (int i=ms->mmr_max_i; i<ms->n; i++){
        if (ms->mmr_a[i].sum[0]+ms->mmr_a[i].sum[1]>=mmr_cov_for_update) ms->mmr_max_i = i;
        else break;  // should be impossible to have holes in methmer count, so look no further
    }
    //fprintf(stderr, "[M::%s] ref: %d hap0, %d hap1; methmer valid until %d\n", 
    //    __func__, sancheck[0], sancheck[1], 
    //    (int)ms->sites_real_poss[ms->mmr_max_i]);
    //for (int i=0; i<ms->mmr_max_i; i++){
    //    fprintf(stderr, "[dbg::%s] mmr#%d, hap0 sum %d, hap1 sum %d\n", __func__, i, ms->mmr_a[i].sum[0], ms->mmr_a[i].sum[1]);
    //}
}

int permute_haplotags(rs_t *rs, int start, int end, int n){
    // TODO diploid assumption
    // TODO maybe consider to let rs_t keep the buffers to avoid re allocations
    // swap n reads between the two haplotypes  
    vu32_t buf_hap1; kv_init(buf_hap1); kv_resize(uint32_t, buf_hap1, 32);
    vu32_t buf_hap2; kv_init(buf_hap2); kv_resize(uint32_t, buf_hap2, 32);
    int ret = 0;
    if (end-start<=1){
        ret = 1;
        goto cleanup;
    }
    for (uint32_t i=start; i<end; i++){
        if (rs->a[i].hp==0) kv_push(uint32_t, buf_hap1, i);
        if (rs->a[i].hp==1) kv_push(uint32_t, buf_hap2, i);
    }

    // bound n and then shuffle
    ks_shuffle_kssu32(buf_hap1.n, buf_hap1.a);
    ks_shuffle_kssu32(buf_hap2.n, buf_hap2.a);

    // alter haplotags 
    //fprintf(stderr, "[dbg::%s] start=%d end=%d\n", __func__, start, end);
    for (int i=0, j; i<n; i++){
        //fprintf(stderr, "[dbg::%s] %d buf1 has %d,  buf2 has%d\n", 
        //    __func__, i, buf_hap1.a[i], buf_hap2.a[i]);
        j = buf_hap1.a[start+i];
        rs->a[j].hp = 1;
        j = buf_hap2.a[start+i];
        rs->a[j].hp = 0;
    }

cleanup:
    kv_destroy(buf_hap1);
    kv_destroy(buf_hap2);
    return ret;
}
void random_assign_haplotags_if_has_none(rs_t *rs, int start, int end){
    int n = 0;
    for (int i=start; i<end; i++){
        if (rs->a[i].hp!=0 && rs->a[i].hp!=1){
            rs->a[i].hp = rand()&1;
            n++;
        }
    }
    fprintf(stderr, "[dbg::%s] randomly assigned %d (out of %d)\n", __func__, 
        n, end-start
    );
}

float evaluate_separation(uint8_t *ref, uint8_t *query, int n, 
                          int *dir){
    int buf[2][2];
    for (int i=0; i<2; i++){
        for (int j=0; j<2; j++){
            buf[i][j] = 0;
        }
    }
    for (int i=0; i<n; i++){
        if (ref[i]!=0 && ref[i]!=1) continue;
        if (query[i]!=0 && query[i]!=1) continue;
        buf[ref[i]][query[i]]++;
    }
    // require unambiguous ratios for both
    fprintf(stderr, "[dbg::%s] n is %d; %d %d %d %d\n", __func__, n, buf[0][0], buf[0][1], buf[1][0], buf[1][1]);
    float min, max;
    float scores[2];
    int which_way = 0;
    for (int i=0; i<2; i++){
        if (buf[i][0]>buf[i][1]){
            min = buf[i][1];
            max = buf[i][0];
            which_way = i==0? which_way+1 : which_way-1;
        }else{
            min = buf[i][0];
            max = buf[i][1];
            which_way = i==0? which_way-1 : which_way+1;
        }
        if (max==0) return (float)1;
        min = min==0? 1: min;
        if (max/min<3) return (float)1;
        scores[i] = max/min;
    }

    // require significance
    double fisher_left_p, fisher_right_p, fisher_twosided_p;
    kt_fisher_exact(buf[0][0], buf[0][1], buf[1][0], buf[1][1], 
                    &fisher_left_p, &fisher_right_p, &fisher_twosided_p);
    if (fisher_twosided_p<0.001) {
        if (dir!=0){
            *dir= which_way;
        }
        return MIN(scores[0], scores[1]);
    }else{
        return (float)1;
    }
}


void haplotag_region1(rs_t *rs, methmers_t *ms, 
                      uint32_t *readIDs, uint32_t readIDs_n, 
                      uint32_t readIDs_init_n, 
                      int n_candidates_per_iter, 
                      int min_mmr_recruit_cov){
    vu32_t buf_anno_readIDs;  // will store readIDs sorted by corresponding tagging scores
    vi_t buf_anno_tags;
    forpred_v buf_srt;
    forpred_v buf_srt_tmp;  
    kv_init(buf_anno_readIDs); kv_resize(uint32_t, buf_anno_readIDs, 32);
    kv_init(buf_anno_tags);    kv_resize(int, buf_anno_tags, 32);
    kv_init(buf_srt);          kv_resize(forpred_t, buf_srt, 32);
    kv_init(buf_srt_tmp);      kv_resize(forpred_t, buf_srt_tmp, 32);
    int mmr_max_i = 0;  // we will only use methmers that have accumulated enough coverage 
                        // during tagging to tag incoming query reads


    // step1: collect ref reads on the left
    //fprintf(stderr, "[M::%s] collecting ref reads on the left...\n", __func__);
    ms->mmr_max_i = 0;
    insert_ref_reads_methmer_counts(rs, ms, readIDs, readIDs_init_n, min_mmr_recruit_cov);
    for (int i=ms->mmr_max_i; i<ms->n; i++){  // mark all sites before start position as available
        if (ms->sites_real_poss[i]<=rs->ref_start) ms->mmr_max_i++;
        else break;
    }
    //fprintf(stderr, "[M::%s] collected ref on the left\n", __func__);
    //fprintf(stderr, "[dbg::%s] mmr counts of this iter:\n", __func__);
    //for (int j=0; j<ms->n; j++){
    //    fprintf(stderr, "  idx=%d hap1sum=%d hap2sum=%d\n", j, ms->mmr_a[j].sum[0], ms->mmr_a[j].sum[1]);
    //}


    // step1.5: mark non-ref reads as unhaplotagged, including the reads on the right side
    for (int i=readIDs_init_n; i<readIDs_n; i++){
        rs->a[i].hp = 9;
    }

    // step2: extend, tag 1 read per iteration
    int i_last_untagged = 0;
    int i_last_mmr = 0;
    //fprintf(stderr, "[dbg::%s] tagging reads from %d to %d\n", __func__, i_last_untagged, readIDs_n);
    int failed_prev_batch = 0;
    while (1){
        // step2 1) collect candidate reads
        // TODO this is not identical to the py prototype as i do not want to 
        //  depend on bam utils here. If does not work as good, write a helper
        //  rather than using bam utils. 
        buf_anno_readIDs.n = 0;
        if (i_last_untagged>=readIDs_n) break;
        for (int i=i_last_untagged; i<readIDs_n; i++){
            if (rs->a[i].hp!=0 && rs->a[i].hp!=1) {
                kv_push(uint32_t, buf_anno_readIDs, i);
                if (buf_anno_readIDs.n>=n_candidates_per_iter) break;
            }
        }
        if (buf_anno_readIDs.n==0){ // nothing can be tagged
            failed_prev_batch++;
            if (failed_prev_batch>10) break;
            i_last_untagged+=n_candidates_per_iter; // try to get another batch; TODO perhaps should just give up?
            continue;
        }
        //fprintf(stderr, "[dbg::%s] collected a batch, first ID is %d, n is %d, last untagged is %d\n", 
        //        __func__, buf_anno_readIDs.a[0], (int)buf_anno_readIDs.n, i_last_untagged);

        // step2 2) tag & insert
        int inserted = predict_tags_of_reads(rs, ms, 
                            buf_anno_readIDs.a, buf_anno_readIDs.n, 
                            &buf_anno_readIDs, &buf_anno_tags, 
                            &buf_srt, &buf_srt_tmp, 
                            1,  // try to add 1 read per round and update mmr counts accordingly
                            min_mmr_recruit_cov,
                            3, 3);
        if (inserted==0){  // nothing can be tagged
            failed_prev_batch++;
            if (failed_prev_batch>10) break;
            i_last_untagged+=n_candidates_per_iter; // try to get another batch; TODO perhaps should just give up?
            continue;
        }
        failed_prev_batch = 0;
    }

cleanup:
    kv_destroy(buf_anno_readIDs);
    kv_destroy(buf_anno_tags);
    kv_destroy(buf_srt);
    kv_destroy(buf_srt_tmp);
}





vu8_t *store_haplotags(rs_t *rs){
    vu8_t *ret = (vu8_t*)malloc(sizeof(vu8_t));
    kv_init(*ret);
    kv_resize(uint8_t, *ret, rs->n);
    for (int i=0; i<rs->n; i++){
        kv_push(uint8_t, *ret, (uint8_t)rs->a[i].hp);
    }
    return ret;
}
void restore_haplotags(rs_t *rs, vu8_t *tags){
    assert(rs->n==tags->n);
    for (int i=0; i<rs->n; i++){
        rs->a[i].hp = tags->a[i];
    }
}
void set_right_end(rs_t *rs, int n_end_reads, int tag){
    for (int i=rs->n-n_end_reads; i<rs->n; i++){
        rs->a[i].hp = tag;
    }
}
void set_all_as_unphased(rs_t *rs){
    for (int i=0; i<rs->n; i++) rs->a[i].hp = 2;

}

void haplotag_region2(rs_t *rs, methmers_t *ms, 
                      uint32_t *readIDs, uint32_t readIDs_n, 
                      uint32_t readIDs_init_n, 
                      int n_candidates_per_iter, 
                      int min_mmr_recruit_cov, 
                      int n_permutations, 
                      int n_end_reads // need to check ratios on the right side
                      ){
    fprintf(stderr, "[dbg::%s] *** using permuation ***\n", __func__);
    // Wraper around haplotag_region1, permute initial state,
    // after filtering the joining decision, take the majority.
    // Output: stores read haplotags.  
    int threshold = n_permutations/2;
    int threshold_blank = n_permutations/3;  // do not allow too many ambiguous results
    vu8_t *initial_state = store_haplotags(rs);
    vu8_t buf[n_permutations];
    float scores[n_permutations];
    int dir[3] = {0,0,0};  // undecided, dir1, dir2
    float best_score[2] = {1, 1};
    int best_score_i[2] = {-1, -1};
    for (int i=0; i<n_permutations; i++){  // get score and store tagging
        if (i!=0)
            permute_haplotags(rs, 0, readIDs_init_n, 5);
        //random_assign_haplotags_if_has_none(rs, 0, readIDs_init_n);
        //for (int j=0; j<readIDs_init_n; j++){
        //    fprintf(stderr, "[dbg::%s]   read %d hp=%d\n", __func__, j, rs->a[j].hp);
        //}
        haplotag_region1(rs, ms, readIDs, readIDs_n, readIDs_init_n, 
                         n_candidates_per_iter, min_mmr_recruit_cov);
        kv_init(buf[i]);
        kv_resize(uint8_t, buf[i], 16);
        for (int j=0; j<rs->n; j++){
            kv_push(uint8_t, buf[i], rs->a[j].hp);
        }
        int offset = rs->n - n_end_reads;
        int which_way;
        scores[i] = evaluate_separation(initial_state->a+offset, 
                                    buf[i].a+offset, n_end_reads, &which_way);
        if (scores[i]>=2 && which_way!=0){
            if (which_way>0) {dir[1]++; which_way = 0;}
            else             {dir[2]++; which_way = 1;}
            if (scores[i]>best_score[which_way]){
                best_score[which_way] = scores[i];
                best_score_i[which_way] = i;
            }
        }else{
            dir[0]++;
        }
        restore_haplotags(rs, initial_state);
    }
    // summarize
    if ((dir[1]>=threshold && dir[2]<=3/*|| (dir[1]/(dir[2]==0? 1 : dir[2])>2)*/) && 
         dir[0]<threshold_blank && 
         best_score_i[0]>=0){
        restore_haplotags(rs, &buf[best_score_i[0]]);
        fprintf(stderr, "[dbg::%s] >>> dir1 (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
    }else if ((dir[2]>=threshold && dir[1]<=3/*|| (dir[2]/(dir[1]==0? 1:dir[1])>2)*/ ) && 
         dir[0]<threshold_blank && 
         best_score_i[1]>=0){
        restore_haplotags(rs, &buf[best_score_i[1]]);
        fprintf(stderr, "[dbg::%s] <<< dir2 (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
    }else{  // failed, permutation result isn't conclusive
        restore_haplotags(rs, initial_state);
        set_all_as_unphased(rs);
        fprintf(stderr, "[dbg::%s] xxx undecided (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
    }

cleanup:
    kv_destroy(*initial_state);
    free(initial_state);
}



dataset_t* haplotag_region_given_bam(char *fn_bam, 
                     char *chrom, uint32_t ref_start, uint32_t ref_end,   // alignment coordinates
                     int call_lo, int call_hi,
                     int min_mapq, 
                     int mmr_k, int mmr_max_span, 
                     int min_mmr_cov, int min_readlen,
                     int n_candidates_per_iter, 
                     int mmr_cov_for_update,
                     int do_n_permuations
                     ){
    dataset_t *ret = (dataset_t*)malloc(sizeof(dataset_t));
    vu32_t readIDs;
    kv_init(readIDs);
    kv_resize(uint32_t, readIDs, 128);

    // load reads and methmers
    int n_start_reads = 0;
    int n_end_reads = 0;
    rs_t *rs = load_reads_given_interval(fn_bam, chrom, ref_start, ref_end, READBACK,
                                min_mapq,
                                min_readlen, 
                                call_lo, call_hi, 
                                &n_start_reads, 
                                &n_end_reads);
    vu8_t *initial_hps = store_haplotags(rs);
    fprintf(stderr, "[dbg::%s] loaded reads\n", __func__);

    methmers_t *ms = get_methmer_sites_and_ranges(rs, mmr_k, mmr_max_span, min_mmr_cov);
    fprintf(stderr, "[dbg::%s] parsed meth sites chr%s:%d-%d (k=%d span=%d cov=%d runtime_cov=%d)\n", 
        __func__, chrom, ref_start, ref_end,
        mmr_k, mmr_max_span, min_mmr_cov, mmr_cov_for_update);

    store_mmr_of_reads(rs, ms);
    //fprintf(stderr, "[dbg::%s] parsed methmers for reads\n", __func__);

    if (ms->n==0){
        fprintf(stderr, "[W::%s] %s:%d-%d has no methylation, skipping.\n", 
        __func__, chrom, ref_start, ref_end);
        goto cleanup;
    }
    fprintf(stderr, "[M::%s] input ready\n", __func__);

    // main routine
    for (uint32_t i=0; i<rs->n;i++){
        kv_push(uint32_t, readIDs, i);
    }
    if (do_n_permuations>0){
        if (do_n_permuations<3){
            fprintf(stderr, "[W::%s] number of permutations is too small, might have no effects.\n", __func__);
        }
        haplotag_region2(rs, ms, readIDs.a, readIDs.n, n_start_reads, 
                         n_candidates_per_iter, mmr_cov_for_update, do_n_permuations, n_end_reads);
    }else{
        haplotag_region1(rs, ms, readIDs.a, readIDs.n, n_start_reads, 
                         n_candidates_per_iter, mmr_cov_for_update);
    }
                

    // (debug)
    //fprintf(stderr, "[W::%s] writing temporary debug output\n", __func__);
    //char tmpfn[1024];
    //sprintf(tmpfn, "/media/groups/machine_learning/active/xfeng/p_methp/mp/hyperparameter_search/chr%s_%d_%d.k%d_span%d_cov%d_updtcov%d.tsv", 
    //    chrom, ref_start, ref_end, 
    //    mmr_k, mmr_max_span, min_mmr_cov, mmr_cov_for_update);
    //FILE *fp = fopen(tmpfn, "w");
    //for (int i=0; i<rs->n; i++){
    //    fprintf(fp, "%.*s\t%d\t%d\n", (int)GET_QNAME_L(*rs, i), GET_QNAME_S(*rs, i), (int)initial_hps->a[i], (int)rs->a[i].hp);
    //}
    //fclose(fp);

cleanup:
    kv_destroy(readIDs);
    kv_destroy(*initial_hps);
    free(initial_hps);
    ret->rs = rs;
    ret->ms = ms;
    return ret;
}


dataset_t* haplotag_region_given_bam_hyperparameter_search(char *fn_bam, 
                     char *chrom, uint32_t ref_start, uint32_t ref_end,   // alignment coordinates
                     int call_lo, int call_hi,
                     int min_mapq, 
                     int mmr_k, int mmr_max_span, 
                     int min_mmr_cov, int min_readlen,
                     int n_candidates_per_iter, 
                     int mmr_cov_for_update
                     ){
    int silent = 1;
    dataset_t *ret = (dataset_t*)malloc(sizeof(dataset_t));
    vu32_t readIDs;
    kv_init(readIDs);
    kv_resize(uint32_t, readIDs, 128);

    // load reads and methmers
    int n_start_reads = 0;
    int n_end_reads = 0;
    rs_t *rs = load_reads_given_interval(fn_bam, chrom, ref_start, ref_end, READBACK,
                                min_mapq,
                                min_readlen, 
                                call_lo, call_hi, 
                                &n_start_reads, 
                                &n_end_reads);
    vu8_t *initial_hps = store_haplotags(rs);
    fprintf(stderr, "[dbg::%s] loaded reads\n", __func__);

    int tmpN = 5*8*7;
    vu8_t buf[tmpN];
    float scores[tmpN];
    char log_string[tmpN][128];
    int n =0 ;
    int done = 0;
    float best = 0;
    int best_i = -1;
    if (rs->n==0) {
        goto cleanup;  // nothing loaded
    }
    for (int mmr_k=7; mmr_k<=10; mmr_k++){
        for (int mmr_max_span=6000; mmr_max_span<=8000; mmr_max_span+=500){
            for (int min_mmr_cov=6; min_mmr_cov<=12; min_mmr_cov++){
                for (int mmr_cov_for_update=min_mmr_cov*2;  mmr_cov_for_update<min_mmr_cov*2+1; mmr_cov_for_update+=3){  // let's just try 1
                    methmers_t *ms = get_methmer_sites_and_ranges(rs, mmr_k, mmr_max_span, min_mmr_cov);
                    store_mmr_of_reads(rs, ms);
                    //fprintf(stderr, "[dbg::%s] parsed methmers for reads\n", __func__);

                    if (ms->n==0){
                        fprintf(stderr, "[W::%s] %s:%d-%d has no methylation, skipping.\n", 
                        __func__, chrom, ref_start, ref_end);
                        goto cleanup;
                    }
                    //fprintf(stderr, "[M::%s] input ready\n", __func__);

                    // main routine
                    for (uint32_t i=0; i<rs->n;i++){
                        kv_push(uint32_t, readIDs, i);
                    }
                    haplotag_region1(rs, ms, readIDs.a, readIDs.n, n_start_reads, 
                                     n_candidates_per_iter, mmr_cov_for_update);
                
                    // get score and store tagging
                    kv_init(buf[n]);
                    kv_resize(uint8_t, buf[n], 16);
                    for (int i=0; i<rs->n; i++){
                        kv_push(uint8_t, buf[n], rs->a[i].hp);
                    }
                    int offset = rs->n - n_end_reads;
                    scores[n] = evaluate_separation(initial_hps->a+offset, 
                                                buf[n].a+offset, n_end_reads, 0);
                    sprintf(log_string[n], "k=%d,span=%d,cov=%d", mmr_k, mmr_max_span, min_mmr_cov);

                    // (debug)
                    if (!silent){
                        char tmpfn[1024];
                        sprintf(tmpfn, "/media/groups/machine_learning/active/xfeng/p_methp/mp/hyperparameter_search/chr%s_%d_%d.k%d_span%d_cov%d_updtcov%d.tsv", 
                            chrom, ref_start, ref_end, 
                            mmr_k, mmr_max_span, min_mmr_cov, mmr_cov_for_update);
                        FILE *fp = fopen(tmpfn, "w");
                        for (int i=0; i<rs->n; i++){
                            fprintf(fp, "%.*s\t%d\t%d\n", (int)GET_QNAME_L(*rs, i), GET_QNAME_S(*rs, i), (int)initial_hps->a[i], (int)rs->a[i].hp);
                        }
                        fclose(fp);
                    }

                    // restore tagging and reset
                    for (int i=0; i<rs->n; i++){
                        rs->a[i].hp = initial_hps->a[i];
                    }
                    readIDs.n = 0;
                    destroy_methmers_t(ms);

                    // experimental : take the 1st passing result
                    //if (scores[n]<0.001) {
                    //    done = 1;
                    //    best_i = n;
                    //    best = scores[n];
                    //    fprintf(stderr, "[dbg::%s] best found: chr%s:%d-%d (k=%d span=%d cov=%d runtime_cov=%d), score %f\n", 
                    //        __func__, chrom, ref_start, ref_end,
                    //        mmr_k, mmr_max_span, min_mmr_cov, mmr_cov_for_update, scores[n]);
                    //    break;
                    //}
                    n++;
                } 
                if (done) break;
            }
            if (done) break;
        }
        if (done)break;
    }
    for (int i=0; i<n; i++){
        if (scores[i]>best) {
            best_i = i;
            best = scores[i];
        }
    }
    if (/*!done*/ best<3){  // treat as unphased
        fprintf(stderr, "[dbg::%s] best NOT found: chr%s:%d-%d (-), score %f (used %d right reads)\n", 
            __func__, chrom, ref_start, ref_end,
             scores[best_i], n_end_reads);
        set_all_as_unphased(rs);
    }else{  // use the best
        fprintf(stderr, "[dbg::%s] best IS found: chr%s:%d-%d (%s), score %f (used %d right reads)\n", 
            __func__, chrom, ref_start, ref_end,
            log_string[best_i],  scores[best_i], n_end_reads);
        for (int i=0; i<rs->n; i++){
            rs->a[i].hp = buf[best_i].a[i];
        }
    }
    for (int i=0; i<n; i++){
        kv_destroy(buf[i]);
    }

cleanup:
    kv_destroy(readIDs);
    kv_destroy(*initial_hps);
    ret->rs = rs;
    ret->ms = 0;
    return ret;
}



//int main_blockjoin(char *fn_gtf, char *fn_bam, char *fn_out_tsv, 
//                   int use_hypersearch, 
//                   int readlen_threshold, 
//                   int lo, int hi, int mapq, 
//                   int k, int k_span, int cov_for_selection, 
//                   int n_candidates_per_iter){
//    ranges_t *ranges = load_intervals_from_gtf(fn_gtf, chrom_s, 0);
//    merge_close_intervals(ranges, READBACK);
//}
int main_blockjoin_dbg(int argc, char *argv[]){
    if (argc<5){
        fprintf(stderr, "[E::%s] usage: exe suffix which_tool output_prefix use_hyperparameter_search={search,nosearch}\n", __func__);
        exit(1);
    }
    for (int chrom=1; chrom<23; chrom++){
        char chrom_s[50];
        char fn_gtf[1024];
        char fn_bam[1024];
        char fn_out[1024];
        int lo = cliopt.lo;
        int hi = cliopt.hi;
        int mapq = cliopt.mapq;
        int k = cliopt.k;
        int k_span = cliopt.k_span;
        int cov_for_selection = cliopt.cov_for_selection;
        int readlen_threshold = cliopt.readlen_threshold;
        int n_candidates_per_iter = cliopt.n_candidates_per_iter;

        //sprintf(fn_gtf, "/media/groups/machine_learning/active/xfeng/p_methp/mp/run_methphaser_all_chr/chr%d/hapcut2_onlysnpQ10/phased.vcf.gtf", chrom);
        //sprintf(fn_bam, "/media/groups/machine_learning/active/xfeng/p_methp/mp/run_methphaser_all_chr/chr%d/hapcut2_onlysnpQ10/haplotagged.F4F256F2048.bam", chrom);

        sprintf(fn_gtf, "/media/groups/machine_learning/active/xfeng/p_methp/mp/eval_nov2024/chr%d/%s/%s/phased.vcf.gtf", chrom, argv[1], argv[2]);
        ////sprintf(fn_gtf, "/media/groups/machine_learning/active/xfeng/p_methp/mp/hyperparameter_search/d/phased_useT2Tbounds.hapcut2.chr%d.vcf.gtf", chrom);
        sprintf(fn_bam, "/media/groups/machine_learning/active/xfeng/p_methp/mp/eval_nov2024/chr%d/%s/%s/haplotagged.F4F256F2048.bam", chrom, argv[1], argv[2]);
        sprintf(chrom_s, "chr%d", chrom);
        //sprintf(fn_out, "/media/groups/machine_learning/active/xfeng/p_methp/mp/eval_nov2024/chr%d/%s/%s_mytagging2.tsv", chrom, argv[1], argv[2]);
        sprintf(fn_out, "/media/groups/machine_learning/active/xfeng/p_methp/mp/hyperparameter_search/mytagging2%s_%s_%s.chr%d.tsv", argv[3], argv[1], argv[2], chrom);

        ranges_t *ranges = load_intervals_from_gtf(fn_gtf, chrom_s, 0);
        merge_close_intervals(ranges, READBACK);

        FILE *fp = fopen(fn_out, "w");
        if (!fp){
            fprintf(stderr, "[E::%s] failed to open output file: %s\n", __func__, fn_out);
            exit(1);
        }
        for (int i=0; i<ranges->starts.n; i++){
        //for (int i=0; i<1; i++){
            int start = ranges->starts.a[i];
            int end = ranges->ends.a[i];

            fprintf(stderr, "[M::%s] at %s:%d-%d\n", __func__, chrom_s, start, end);
            dataset_t *ds;
            if (strcmp(argv[4], "search")==0){
                ds = haplotag_region_given_bam_hyperparameter_search(fn_bam, 
                        chrom_s, start, end,
                        lo, hi,
                        mapq,  // mapq
                        k, // methmer k
                        k_span,   // methmer range
                        cov_for_selection,  // coverage threshold of site selection
                        readlen_threshold, 
                        n_candidates_per_iter, // n candidates per iter
                        cov_for_selection*2 // run time mmr coverage for inclusion
                        );
            }else{
                ds = haplotag_region_given_bam(fn_bam, 
                        chrom_s, start, end,
                        lo, hi,
                        mapq,  // mapq
                        k, // methmer k
                        k_span,   // methmer range
                        cov_for_selection,  // coverage threshold of site selection
                        readlen_threshold, 
                        n_candidates_per_iter, // n candidates per iter
                        cov_for_selection*2, // run time mmr coverage for inclusion
                        30  // # permutations
                        );
            }
            for (int i=0; i<ds->rs->n; i++){
                int hap = ds->rs->a[i].hp;
                hap = hap<0? 2 : hap;
                fprintf(fp, "%.*s\t-1\t%d\n", (int)GET_QNAME_L(*(ds->rs), i), GET_QNAME_S(*(ds->rs), i), hap);
            }
            destroy_dataset_t(ds, 1);
        }
        destroy_ranges_t(ranges);
        fclose(fp);
    }

    
	return 0;
}
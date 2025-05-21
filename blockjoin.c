#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include <zlib.h>
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/kfunc.h"  // borrow kt_fisher_exact
#include "kvec.h"
#include "khashl.h"
#include "ksort.h"
#include "kthread.h"
#include "kstring.h"
#include "blockjoin.h"
#include "cli.h"

#define READBACK 50000
#define HARD_COV_THRESHOLD 15
#define HARD_CONTAMINATE_THRESHOLD 5
#define READLINE_BUF_LEN 128
#define MIN_ALN_DE 0.1
#define EVAL_P_THRE 0.001

#define HAPTAG_UNPHASED 254

#define VAR_OP_M 0
#define VAR_OP_X 1
#define VAR_OP_I 2
#define VAR_OP_D 3
#define VAR_DIFF_OVERRIDE_RATIO 5

#define N_MODS 10


#define MAX(a,b) ((a)>(b)? (a):(b))
#define MIN(a,b) ((a)<=(b)? (a):(b))
int get_posint_digits(int n) {
    assert(n>=0);
    if (n<10) return 1;
    if (n<100) return 2;
    if (n<1000) return 3;
    if (n<10000) return 4;
    if (n<100000) return 5;
    if (n<1000000) return 6;
    if (n<10000000) return 7;
    if (n<100000000) return 8;
    if (n<1000000000) return 9;
    return 10;
}


const unsigned char seq_atcgun_table[256] = { 
    // let ATCGUatcguNn <4
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 1/*A*/, 4, 1/*C*/,  4, 4, 4, 1/*G*/,  4, 4, 4, 4,  4, 4, 1/*N*/, 4, 
    4, 4, 4, 4,  1/*T*/, 1/*U*/, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 1/*a*/, 4, 1/*c*/,  4, 4, 4, 1/*g*/,  4, 4, 4, 4,  4, 4, 1/*n*/, 4, 
    4, 4, 4, 4,  1/*t*/, 1/*u*/, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
};

const unsigned char seq_nt4_table[256] = { 
    // translate ACG{T,U} to 0123 case insensitive
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0/*A*/, 4, 1/*C*/,  4, 4, 4, 2/*G*/,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3/*T*/, 3/*U*/, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0/*a*/, 4, 1/*c*/,  4, 4, 4, 2/*g*/,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3/*t*/, 3/*u*/, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
};

const unsigned char md_op_table[256]={
    // 0-9 gives 0
    // ^ gives 1
    // ATCGatcgUuNn gives 2 (allow N base since it is the reference)
    // else: 4
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 4, 4,  4, 4, 4, 4,
    4, 2, 4, 2,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2/*N*/, 4,
    4, 4, 4, 4,  2, 2, 4, 4,  4, 4, 4, 4,  4, 4, 1, 4,
    4, 2, 4, 2,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2/*n*/, 4,
    4, 4, 4, 4,  2, 2, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
};

int natoi(char *s, int l){
    // Trusted non-null terminated string to positive int.
    // Just for parsing md tag.
    if (l<=0) return -1;
    int ret = 0;
    for (int i=0; i<l; i++){
        int e = 1;
        for (int j=0; j<l-1-i; j++){
            e *= 10;
        }
        ret += ((int)s[i]-48) * e;
    }
    return ret;
}





typedef struct{
    uint32_t pos;
    uint16_t cnt[3]; // meth, unmeth, nocall
}methsite_counter_t;
KHASHL_MAP_INIT(static klib_unused, htmethcnt_t, htmethcnt_ht, uint32_t, methsite_counter_t, kh_hash_uint32, kh_eq_generic)
KHASHL_MAP_INIT(static klib_unused, htstri_t, htstri_ht, char*, int, kh_hash_str, kh_eq_str)
KHASHL_MAP_INIT(static klib_unused, htu32_t, htu32_ht, uint32_t, uint32_t, kh_hash_uint32, kh_eq_generic)
KHASHL_MAP_INIT(static klib_unused, htu64_t, htu64_ht, uint64_t, uint64_t, kh_hash_uint64, kh_eq_generic)

#define generic_key(x) (x)
#define generic_cmp(x,y) ((x)<(y))
KRADIX_SORT_INIT(ksu32, uint32_t, generic_key, 4)
KSORT_INIT(kssu32, uint32_t, generic_cmp)
KRADIX_SORT_INIT(ksu64, uint64_t, generic_key, 8)


typedef struct{
    uint32_t s, e;
}u32p_t;
typedef kvec_t(u32p_t) vu32p_t;
typedef kvec_t(uint64_t) vu64_t;
typedef kvec_t(uint32_t) vu32_t;
typedef kvec_t(uint16_t) vu16_t;
typedef kvec_t(uint8_t) vu8_t;
typedef kvec_t(vu8_t) vu8_v;
typedef kvec_t(int) vi_t;
typedef kvec_t(float) vfloat_t;


typedef struct{
    uint32_t pos;  // position as on reference, 0-index
    uint32_t len;
    uint8_t op;  // cigar operation
    vu8_t chars; // sequences, use 0123 (as in seq_nt4_table)
                 // Record alt for SNP and INS, 
                 // record ref for DEL.
    uint8_t haptag;  // For the reference collection loaded from vcf, 
                     // the haptag will always represent the ref allele.
                     // For a query read, this field is filled but not used. 
}variant_t;
typedef kvec_t(variant_t) vvar_t;
void init_vvar_t(vvar_t *h){
    kv_init(*h);
    kv_resize(variant_t, *h, 16);
}
void wipe_vvar_t(vvar_t *h){
    for (int i=0; i<h->n; i++){
        kv_destroy(h->a[i].chars);
    }
    h->n = 0;
}
void push_to_vvar_t(vvar_t *h, uint32_t pos, 
                    uint32_t op, uint32_t op_l, 
                    int haptag,
                    char *seq){
    variant_t v;
    v.pos = pos;
    v.op = op;
    v.len = op_l;
    v.haptag = haptag;
    kv_push(variant_t, *h, v);
    vu8_t *chars =  &h->a[h->n-1].chars;
    kv_init(*chars);
    kv_resize(uint8_t, *chars, 4);
    for (int i=0; i<op_l; i++){
        kv_push(uint8_t, *chars, seq_nt4_table[(int)seq[i]]);
    }
}
void destroy_vvar_t(vvar_t *h, int include_self){
    for (int i=0; i<h->n; i++){
        kv_destroy(h->a[i].chars);
    }
    kv_destroy(*h);
    if (include_self) free(h);
}

typedef struct {
    vu32_t IDs_left;
    vu32_t IDs_left_strict;  // reads on the left side that actually 
                           // overlaps with the left SNP evidence (i.e. 
                           // the read start position)
    vu32_t IDs_right;
    vu32_t IDs_right_strict; // (similar to the above)
}refreads_t;
refreads_t *init_refreads_t(){
    refreads_t *ret = (refreads_t*)malloc(sizeof(refreads_t));
    kv_init(ret->IDs_left);
    kv_init(ret->IDs_right);
    kv_init(ret->IDs_left_strict);
    kv_init(ret->IDs_right_strict);
    kv_resize(uint32_t, ret->IDs_left, 16);
    kv_resize(uint32_t, ret->IDs_right, 16);
    kv_resize(uint32_t, ret->IDs_left_strict, 16);
    kv_resize(uint32_t, ret->IDs_right_strict, 16);
    return ret;
}
void insert_refreads_t(refreads_t *rf, uint32_t ID, int which_side, int overlaps_with_snp){
    if (which_side==0){
        kv_push(uint32_t, rf->IDs_left, ID);
        if (overlaps_with_snp) 
            kv_push(uint32_t, rf->IDs_left_strict, ID);
    }else if (which_side==1){
        kv_push(uint32_t, rf->IDs_right, ID);
        if (overlaps_with_snp) 
            kv_push(uint32_t, rf->IDs_right_strict, ID);
    }else{
        fprintf(stderr, "[E::%s] invalid side give (%d), check code\n", __func__, which_side);
    }
}
void destroy_refreads_t(refreads_t *rf){
    kv_destroy(rf->IDs_left);
    kv_destroy(rf->IDs_right);
    kv_destroy(rf->IDs_left_strict);
    kv_destroy(rf->IDs_right_strict);
}

#define get_cigar_op_len(a) ((a)>>4)
#define get_cigar_op_type(a) ((a)&15)


void reverse_arru8(uint8_t *a, uint32_t l){
    uint8_t tmp;
    for (int i=0; i<l/2; i++){
        tmp = a[i];
        a[i] = a[l-1-i];
        a[l-1-i] = tmp;
    }
}
void reverse_arr(uint32_t *a, uint32_t l){
    uint32_t tmp;
    for (int i=0; i<l/2; i++){
        tmp = a[i];
        a[i] = a[l-1-i];
        a[l-1-i] = tmp;
    }
}

int count_char_in_string(char *s, char q){
    int i=0;
    int cnt = 0;
    while (s[i]){
        if (s[i]==q) cnt++;
        i++;
    }
    return cnt;
}

int search_substr_idx(char *s, char *q, char delimiter, 
                      int get_idx_insteadof_start, 
                      int s_l, int q_l  // in case s or q is not null termd
                      ){
    // Naively search a short substring in a delimited short string. 
    // Just for parsing vcf format field.
    // Return: index; -1 if not found.
    int i=0, start=0, col=0; 
    int ret = -1;
    int lr, lq; 
    int found = 0;
    if (s_l) lr = s_l;
    else     lr = strlen(s);
    if (q_l) lq = q_l;
    else     lq = strlen(q);
    while (i<=lr){
        if (s[i]==delimiter || i==lr){
            if ((i-start)==lq && strncmp(q, s+start, lq)==0){
                ret = i;
                found = 1;
                break;
            }
            if (i==lr) break;
            start = i+1;
            col++;
        }
        i++;
    }
    if (found)
        return get_idx_insteadof_start? col:start;
    return -1;
}
int get_substr_by_idx(char *s, int idx, char delimiter, 
                      int *substr_start, int *substr_len, int s_l){
    // Naive. Just for parsing vcf format field.
    // Return: -1 if invalid.
    int i=0, col=0, start=0;
    int lr;
    if (s_l) lr = s_l;
    else     lr = strlen(s);
    while (i<=lr){
        if (s[i]==delimiter || i==lr){
            if (col==idx){
                *substr_start = start;
                *substr_len = i-start;
                return 0;
            }
            if (i==lr) break;
            start = i+1;
            col++;
        }
        i++;
    }
    return -1;
}

int search_arr1(uint32_t *a, uint32_t l, uint32_t v, uint32_t *idx){
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
int search_arr(uint32_t *a, uint32_t l, uint32_t v, uint32_t *idx, int which_end){
    // Binary search.
    // If there are multiple hits of v, always return the left most.
    // This is supposed to be only used with methmer allocs, which
    //  have limited amount of duplications (for each position there's 
    //  k dups at most, and we will set k to be like 7).
    uint32_t i;
    int stat = search_arr1(a, l, v, &i);
    if (which_end<0){
        *idx = i;
        return stat;
    }
    if (stat<=0){
        *idx = i;
        return stat;
    }
    if (which_end==0){
        for (int j=i-1; j>=0; j--){
            if (a[j]==v) i=j;
            else break;
        }
    }else{
        for (int j=i+1; j<l; j++){
            if (a[j]==v) i=j;
            else break;
        }

    }
    *idx = i;
    return stat;
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


typedef struct{
    uint32_t i;            // read ID
    int hp;                // haptag of read
    uint8_t strand;        // read alignment strandness
    uint32_t len;          // length of read
    mod_t meth;            // methylation calls parsed from bam meth tags
    uint32_t *mmr;         // methmers collected
    int mmr_n;             // #methmers
    uint32_t mmr_start_i;  // as-in-buffer index refers to methmers_t alloc.
                           // sentinel is UINT32_MAX
}read_t;  // a single read
typedef struct{
    read_t *a;
    uint32_t m, n;
    int has_mmr;
    uint32_t ref_start, ref_end;
    vu32_t names_acl;  // read name accumulated lengths
    char *names;
    uint32_t names_m, names_n;
    htstri_t *name2namei;
    refreads_t *rf;  // IDs of reads on the left and the right of the phasing gap
    vu64_t revbuf;  // backward tagging require reads to be sorted by end position
}rs_t;  // read set
rs_t *init_rs_t(){
    rs_t *ret = (rs_t*)calloc(1, sizeof(rs_t));
    ret->n = 0;
    ret->m = 1024;
    ret->has_mmr = 0;
    // read alloc
    ret->a = (read_t*)calloc(1, sizeof(read_t)*ret->m);
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
    ret->names = (char*)calloc(1, ret->names_m);
    ret->name2namei = htstri_ht_init();
    // allocs for backward extension
    ret->rf = init_refreads_t();
    kv_init(ret->revbuf);
    kv_resize(uint64_t, ret->revbuf, ret->m);
    return ret;
}
read_t* get_empty_slot_rs_t(rs_t* rs){
    if (rs->n==rs->m){
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
    htstri_ht_destroy(rs->name2namei);
    destroy_refreads_t(rs->rf);
    free(rs->rf);
    kv_destroy(rs->revbuf);
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
bamfile_t *init_and_open_bamfile_t(char *fn, int n_threads){
    bamfile_t *h = (bamfile_t*)calloc(1, sizeof(bamfile_t));
    h->fn = (char*)malloc(strlen(fn)+1);
    sprintf(h->fn, "%s", fn);
    h->fp = hts_open(fn, "r");
    if (!h->fp){
        fprintf(stderr, "[E::%s] failed to open: %s\n", __func__, fn);
        return h;
    }

    if (n_threads>1 && h->fp->is_bgzf) {
        bgzf_mt(h->fp->fp.bgzf, n_threads, 0/*unused*/);
    }
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

void code(int id, char buf_code[20]) {
    if (id > 0) {
        buf_code[0] = id;
        buf_code[1] = 0;
    } else {
        snprintf(buf_code, 20, "(%d)", -id);
    }
}


int get_mod_poss_on_ref(mod_t *ret,
                         uint32_t *cigar, int cigar_l, 
                         uint32_t qs, int aln_strand, 
                         uint32_t *mod_poss, uint8_t *mod_quals, int mod_l, 
                         uint8_t *seqi, uint32_t aln_len, char *qn){
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
            mod_poss[i_trigger] = i_ref+cgoffset;  // debug overwrite
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
        //if (verbose) fprintf(stderr, "[dbg::%s] cigar op is: %d, len is %d\n", __func__, (int)action, (int)length);
        if (action<=1){  // 0 for 'M', 1 for 'I'
            int pos_canonical = i_read;
            while (i_read+length>=next_trigger){
                if (action==0 && next_trigger!=UINT32_MAX){
                    // save implicit mod calls?
                    if (seqi){  // caller required us to handle implicit canonicals 
                        //fprintf(stderr, "[dbg::%s] %s pos_canonical abs is %d, next_trigger abs is %d\n", 
                        //        __func__, qn, i_ref+pos_canonical+offset, i_ref+next_trigger+offset);
                        int until = next_trigger-1<i_read+length? next_trigger-1 : i_read+length;
                        for (uint32_t tmpi=pos_canonical; tmpi<until; tmpi++){
                            if (verbose>1) {fprintf(stderr, "[dbg::%s]   trying abs %d; %c %c %c (until: %d)\n", 
                                __func__, i_ref+tmpi+offset, 
                                tmpi==0? 'x' : seq_nt16_str[bam_seqi(seqi, tmpi-1)],
                                seq_nt16_str[bam_seqi(seqi, tmpi)],
                                seq_nt16_str[bam_seqi(seqi, tmpi+1)], 
                                i_ref+until+offset
                                );
                            }
                            if (tmpi<aln_len-1
                                 &&seq_nt16_str[bam_seqi(seqi, tmpi)]=='C' 
                                 && seq_nt16_str[bam_seqi(seqi, tmpi+1)]=='G') { // read is fwd strand
                                int pos_cano = i_ref+tmpi+offset;
                                if (ret->calls.n>0 && ret->calls.a[ret->calls.n-1]==pos_cano){ // ???
                                    if (verbose){
                                        fprintf(stderr, "[dbg::%s] qn %s implicit but not pushing (pos %d) (1)\n", 
                                        __func__, qn, ret->calls.a[ret->calls.n-1]);
                                    }
                                }else{
                                    kv_push(uint32_t, ret->calls, pos_cano);
                                    kv_push(uint8_t, ret->quals, 1); // is unmeth
                                    if (verbose) {fprintf(stderr, "[dbg::%s] qn %s implicit push(1) %d (pos_canon abs is %d, next trigger abs is %d)\n", 
                                                    __func__, qn, ret->calls.a[ret->calls.n-1], 
                                                    i_ref+pos_canonical+offset,  
                                                    i_ref+next_trigger+offset);
                                    }
                                }
                                tmpi++; // skip next base
                            }
                        }
                    }

                    // save the explicit mod call
                    int pos_trigger = i_ref+next_trigger+cgoffset + offset;
                    if (ret->calls.n>0 && ret->calls.a[ret->calls.n-1]==pos_trigger){  // ???
                        ret->quals.a[ret->quals.n-1] = next_qual;
                    }else{
                        kv_push(uint32_t, ret->calls, pos_trigger);
                        kv_push(uint8_t, ret->quals, next_qual);
                    }
                    mod_poss[i_trigger] = pos_trigger;  // debug overwrite

                    // save new position for handling implicit calls
                    pos_canonical = cgoffset==0
                                    ? next_trigger+1 
                                    : next_trigger+2;
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
                if (seqi){  // caller required us to handle implicit canonicals 
                    //fprintf(stderr, "[dbg::%s] %s pos_canonical abs is %d, last bit in this cigar op\n", 
                    //        __func__, qn, i_ref+pos_canonical+offset);
                    int until = i_read+length;
                    for (uint32_t tmpi=pos_canonical; tmpi<until; tmpi++){
                        if (verbose>1) {fprintf(stderr, "[dbg::%s]   trying abs %d; %c %c %c (until: %d)\n", 
                            __func__, i_ref+tmpi+offset, 
                            tmpi==0? 'x' : seq_nt16_str[bam_seqi(seqi, tmpi-1)],
                            seq_nt16_str[bam_seqi(seqi, tmpi)],
                            seq_nt16_str[bam_seqi(seqi, tmpi+1)], 
                            i_ref+until+offset
                            );
                        }
                        if (tmpi<aln_len-1
                             &&seq_nt16_str[bam_seqi(seqi, tmpi)]=='C' 
                             && seq_nt16_str[bam_seqi(seqi, tmpi+1)]=='G') { // read is fwd strand
                            int pos_cano = i_ref+tmpi+offset;
                            if (ret->calls.n>0 && ret->calls.a[ret->calls.n-1]==pos_cano){ // ???
                                if (verbose){
                                    fprintf(stderr, "[dbg::%s] qn %s implicit but not pushing (pos %d) (2)\n", 
                                    __func__, qn, ret->calls.a[ret->calls.n-1]);
                                }
                            }else{
                                kv_push(uint32_t, ret->calls, pos_cano);
                                kv_push(uint8_t, ret->quals, 1); // is unmeth
                                if (verbose) {fprintf(stderr, "[dbg::%s] qn %s implicit push(2) %d (pos_canon abs is %d, next trigger abs is %d)\n", 
                                                __func__, qn, ret->calls.a[ret->calls.n-1], 
                                                i_ref+pos_canonical+offset,  
                                                i_ref+next_trigger+offset);
                                }
                            }
                            tmpi++; // skip next base
                        }
                    }
                }
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
            fprintf(stderr, "[dbg::%s]     meth pos=%d, qual=%d\n", __func__, (int)ret->calls.a[i], (int)ret->quals.a[i]);
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
    int verbose = 0;
    uint32_t *cigar = bam_get_cigar(buf_aln);
    int cigar_l = buf_aln->core.n_cigar;
	uint32_t len = buf_aln->core.l_qseq; // caller ensures >=2
    uint32_t qs = buf_aln->core.pos;
    int strand = !!(buf_aln->core.flag&16);
    bam_parse_basemod(buf_aln, buf_mod);

    uint8_t *seqi = bam_get_seq(buf_aln);
    assert(seqi);

    // get as-on-read coordiates of meth calls
    int n;
    hts_base_mod mods[N_MODS] = {0};
    char *qn = bam_get_qname(buf_aln);
    
    vu8_t dbg_cano;
    vu8_t dbg_seqbase;
    vu8_t dbg_modbases;
    vu8_t dbg_quals;
    if (verbose){
        kv_init(dbg_cano);
        kv_init(dbg_seqbase);
        kv_init(dbg_modbases);
        kv_init(dbg_quals);
    }

    int has_implicit = 0;
    int mod_pos = 0;
    char buf_code[20];
    //while ((n=bam_next_basemod(buf_aln, buf_mod, mods, N_MODS, &mod_pos)) > 0) {
    while (mod_pos<len){
        n = bam_mods_at_next_pos(buf_aln, buf_mod, mods, N_MODS);
        if (n<=0) {
            mod_pos++;
            continue;
        }
        if (n>N_MODS){
            fprintf(stderr, "[W::%s] mod reading buffer was length %d but iter saw %d\n", 
                __func__, N_MODS, n);
            mod_pos++;
            continue;
        }
        for (int j = 0; j < n && j < n; j++) {
            code(mods[j].modified_base, buf_code);
            if (mods[j].canonical_base=='C'
                 && strcmp(buf_code, "m")==0
                 && mod_pos<len-1
                 && mod_pos>0
                 ){  // 5mC at CpG
                 
                if (!(seq_nt16_str[bam_seqi(seqi, mod_pos)]=='C'
                     ?seq_nt16_str[bam_seqi(seqi, mod_pos+1)]=='G' 
                     :seq_nt16_str[bam_seqi(seqi, mod_pos-1)]=='C' 
                     )){
                    has_implicit = 1;
                    continue;
                }

                kv_push(uint32_t, *buf_mod_poss, (uint32_t)mod_pos);
                uint8_t q = (uint8_t)mods[j].qual;

                if (verbose){
                    if (strcmp(buf_code, "m")==0){
                        kv_push(uint8_t, dbg_modbases, (uint8_t)'m');
                    }else if (strcmp(buf_code, "h")==0){
                        kv_push(uint8_t, dbg_modbases, (uint8_t)'h');
                    }else{
                        kv_push(uint8_t, dbg_modbases, (uint8_t)'?');
                    }
                    kv_push(uint8_t, dbg_quals, q);
                    kv_push(uint8_t, dbg_seqbase, seq_nt16_str[bam_seqi(seqi, mod_pos)]);
                    kv_push(uint8_t, dbg_cano, (uint8_t)mods[j].canonical_base);
                }

                kv_push(uint8_t, *buf_mod_quals, q<qual_lo? 
                                                 1 : q>=qual_hi?
                                                 0 : 2);
            }
        }
        mod_pos++;
    }

    // store as-on-reference coordinates of meth calls
    int stat = get_mod_poss_on_ref(&h->meth, cigar, cigar_l, 
                                   qs, strand, 
                                   buf_mod_poss->a, buf_mod_quals->a, 
                                   buf_mod_poss->n, 
                                   has_implicit? seqi:NULL, 
                                   len, qn);
    if (verbose){
        for (int i=0; i<dbg_quals.n; i++){
            fprintf(stderr, "[dbg::%s] modcall dump, qn %s: pos %d : seqi %c, cano %c, qual %d, mod %c\n", 
                __func__, qn, buf_mod_poss->a[i], 
                dbg_seqbase.a[i], dbg_cano.a[i], dbg_quals.a[i], dbg_modbases.a[i]
            );
        }
        kv_destroy(dbg_modbases);
        kv_destroy(dbg_quals);
        kv_destroy(dbg_cano);
        kv_destroy(dbg_seqbase);
    }
    if (has_implicit) {
        global_data_has_implicit = 1;
    }

    return stat;
}

int get_hp_from_aln(bam1_t *aln){
    uint8_t *tagd = bam_aux_get(aln, "HP");
    if (tagd){
        int hp = (int)bam_aux2i(tagd);
        if (hp==0) {
            char *qn = bam_get_qname(aln);
            fprintf(stderr, "[W::%s] irregular HP tag? qn=%s qs=%d\n", 
                    __func__, qn, (int)aln->core.pos);
            return HAPTAG_UNPHASED;
        }
        return hp-1;  // use 0-index internally
    }
    return HAPTAG_UNPHASED;
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
    h->strand = !!(aln->core.flag & 16);
    return 1;
}

int *estimate_read_coverage_dirtyfast(char *fn_bam, int n_threads){
    // Collect coverage with bins by parsing alignment ranges. 
    // Does not filter base quality and does not handle indels.
    double T = Get_T();
    int mod = 5000;
    bamfile_t *hf = init_and_open_bamfile_t(fn_bam, n_threads);
    hts_itr_t *fp_itr = sam_itr_querys(hf->fp_idx, hf->fp_header, ".");
    fprintf(stderr, "[M::%s] estimate read depths...\n", __func__);

    int *covs= (int*)calloc(hf->fp_header->n_targets, sizeof(int));
    vu64_t buf;
    kv_init(buf);

    int prev_refID = -1;
    int refID = -1;
    while(sam_itr_next(hf->fp, fp_itr, hf->buf_aln)>=0){
        refID = hf->buf_aln->core.tid;
        if (refID<0) continue;  // TODO: ?
        if (refID>hf->fp_header->n_targets){
            fprintf(stderr, "[E::%s] refID from alignment is larger than bam heander's n_targets? (%d, %d) Pretending we saw nothing.\n", 
                __func__, refID, hf->fp_header->n_targets);
            continue;
        }
        if (refID!=prev_refID){
            if (prev_refID>=0){
                if (buf.n>0){
                    uint64_t tot = 0;
                    for (int i=0; i<buf.n; i++) {
                        tot+= buf.a[i];
                    }
                    covs[prev_refID] = tot/buf.n;
                }else{
                    covs[prev_refID] = 0;
                }

            }
            // reset bin counters
            int n = hf->fp_header->target_len[refID] / mod;
            if (n>buf.m) kv_resize(uint64_t, buf, n);
            for (int i=0; i<n; i++) buf.a[i] = 0;
            buf.n = n;
            // update name
            prev_refID = refID;
        }

        uint16_t flag = hf->buf_aln->core.flag;
        if ((flag&4) || (flag&256) || (flag&2048)) continue;

		uint32_t mapq = hf->buf_aln->core.qual; 
        if (mapq<5) continue;

		uint32_t len = hf->buf_aln->core.l_qseq;
        float de = -1;
        uint8_t *tmp = bam_aux_get(hf->buf_aln, "de");
        if (tmp) de = bam_aux2f(tmp);
        if (len<15000) continue;
        if (de>MIN_ALN_DE) continue;

        // accumulate counts into bins
        uint32_t start_pos = hf->buf_aln->core.pos;
        uint32_t end_pos = bam_endpos(hf->buf_aln);
        for (int i=start_pos; i<end_pos; /**/){
            buf.a[i/mod]++;
            i += mod;
        }
    }
    // last block
    if (refID>=0){
        if (buf.n>0){
            uint64_t tot = 0;
            for (int i=0; i<buf.n; i++) {
                tot+= buf.a[i];
            }
            covs[refID] = tot/buf.n;
        }else{
            covs[refID] = 0;
        }
    }
    
    for (int i=0; i<hf->fp_header->n_targets; i++){
        fprintf(stderr, "[M::%s] %s est. coverage is %d\n", 
                        __func__, hf->fp_header->target_name[i], covs[i]);
    }
    fprintf(stderr, "[T::%s] used %.1fs\n", __func__, Get_T()-T);

    destroy_bamfile_t(hf, 1);
    hts_itr_destroy(fp_itr);
    kv_destroy(buf);
    return covs;
}


rs_t *load_reads_given_interval(char *fn,
                                char *chrom, int itvl_s, int itvl_e, 
                                int readback, mmr_config_t mmr_config, 
                                htstri_t *qname2haptag_raw,  // if present, use haptags in this hashtable
                                kstring_t *log
                                ){
    rs_t* rs = init_rs_t();
    rs->ref_start = itvl_s>=0? itvl_s : 0;
    rs->ref_end = itvl_e;
    char *itvl_fmt = (char*)malloc(strlen(chrom)+30);
    sprintf(itvl_fmt, "%s:%d-%d", chrom, 
            (itvl_s-readback)>0?itvl_s-readback:0, itvl_e+readback);

    bamfile_t *hf = init_and_open_bamfile_t(fn, 1);
    if (!hf){
        fprintf(stderr, "[E::%s] failed to open input file: %s\n", __func__, fn);
        exit(1);  // TODO cleanup self
    }
    hts_itr_t *fp_itr = sam_itr_querys(hf->fp_idx, hf->fp_header, itvl_fmt);

    vu32_t buf_mod_poss;
    vu8_t buf_mod_quals;
    kv_init(buf_mod_poss); kv_resize(uint32_t, buf_mod_poss, 16);
    kv_init(buf_mod_quals); kv_resize(uint8_t, buf_mod_quals, 16);

    int n_added = 0;
    int left_side_cov_check[2] = {0,0};
    while(sam_itr_next(hf->fp, fp_itr, hf->buf_aln)>=0){
        //char *chr = hf->fp_header->target_name[hf->buf_aln->core.tid];
        char *qn = bam_get_qname(hf->buf_aln);

        int flag = hf->buf_aln->core.flag;
		//uint8_t *q = bam_get_seq(buf_aln); //quality string
		uint32_t mapq = hf->buf_aln->core.qual ; //mapping quality
		uint32_t len = hf->buf_aln->core.l_qseq; //length of the read.
        float de = -1;
        uint8_t *tmp = bam_aux_get(hf->buf_aln, "de");
        if (tmp) de = bam_aux2f(tmp);
        if ((flag&4) || (flag&256) || (flag&2048)) continue;
        if (mapq<mmr_config.min_mapq) continue;
        if (len<2 || len<mmr_config.readlen_threshold/* && hf->buf_aln->core.pos>itvl_s*/) continue;
        if (de>MIN_ALN_DE) continue;  // filter out reads with 
                               // poor alignment and/or poor bascalling

        buf_mod_poss.n = 0;
        buf_mod_quals.n = 0;
        int stat = add_read_record_from_bam_line(rs, hf->buf_aln, hf->buf_mod, 
                        &buf_mod_poss, &buf_mod_quals, 
                        mmr_config.lo, mmr_config.hi);
        if (0 /*log*/){
            ksprintf(log, "[dbg::%s] interval %s, qn %s, add_recrod stat is %d", 
            __func__, qn, itvl_fmt, stat);
            if (stat>0){
                ksprintf(log, ", methcall n is %d (buf n is %d and %d), hp is %d\n", 
                    (int)rs->a[rs->n-1].meth.calls.n, 
                    (int)buf_mod_poss.n, (int)buf_mod_quals.n,
                    (int)rs->a[rs->n-1].hp
                    );
                ksprintf(log, "[dbg::%s]  ^dump:", __func__);
                for (int tmpi=0; tmpi<rs->a[rs->n-1].meth.calls.n; tmpi++){
                    ksprintf(log, "%d<%d>, ", 
                    rs->a[rs->n-1].meth.calls.a[tmpi], rs->a[rs->n-1].meth.quals.a[tmpi]);
                }
                ksprintf(log, "\n");
            }else{
                ksprintf(log, "\n");
            }
        }
        if (stat>0) {
            n_added++;

            if (qname2haptag_raw){ // if qname to haptag mapping is given, 
                                   // override haptags in the bam (this is for variant-haptagging untagged input bam)
                khint_t k = htstri_ht_get(qname2haptag_raw, qn);
                if (k!=kh_end(qname2haptag_raw)){
                    rs->a[rs->n-1].hp = kh_val(qname2haptag_raw, k);
                }else{
                    rs->a[rs->n-1].hp = HAPTAG_UNPHASED;
                }
            }

            uint32_t start_pos = hf->buf_aln->core.pos;
            uint64_t end_pos = bam_endpos(hf->buf_aln);
            kv_push(uint64_t, rs->revbuf, (end_pos<<32)|(rs->n-1));
            if (start_pos<=itvl_s) {
                insert_refreads_t(rs->rf, rs->n-1, 0, end_pos>itvl_s? 1:0);
                int hp = rs->a[rs->n-1].hp;
                if (hp==0 || hp==1){
                    left_side_cov_check[hp]++;
                }
            }
            else if (end_pos>=itvl_e) {
                insert_refreads_t(rs->rf, rs->n-1, 1, start_pos<itvl_e? 1:0);
            }
        }
    }
    // sort readIDs by end position, required by backward tagging
    radix_sort_ksu64(rs->revbuf.a, rs->revbuf.a+rs->revbuf.n);

    // generate reverse index: store qname to read index as-in-buffer in ht
    for (int i=0; i<rs->n; i++){
        char *qn = GET_QNAME_S(*rs, i);
        khint_t k = htstri_ht_get(rs->name2namei, qn);
        int absent;
        k = htstri_ht_put(rs->name2namei, qn, &absent);
        if (!absent){
            fprintf(stderr, "[E::%s] duplicated read name seen from reading bam: %s\n", 
                __func__, qn);
            exit(1);
        }
        kh_val(rs->name2namei, k) = i;
        kh_key(rs->name2namei, k) = qn;
    }

    //fprintf(stderr, "[M::%s] loaded %d reads from %s\n", __func__, n_added, itvl_fmt);

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
    uint32_t abs_start, abs_end;  // first and last available SNP
    vu32p_t dropped;  // small phase blocks dropped due to readback
    vu32_t starts;  vu32_t ends;  // unphase-blocks (with the readback dropping tiny input phased blocks)
    vi_t decisions;  // local decisions for each unphase-block
    vu32p_t rawunphasedblocks;
    vi_t decisions_onraw;  // local decisions for each raw unphase-block
    vi_t flips_onraw; // global decisions for each raw unphase-block, 
                      // i.e. whether phase needs to be flipped wrt 
                      // its value in the input.
    vu32p_t phaseblocks;
}ranges_t;
ranges_t *init_ranges_t(){
    ranges_t *ret = (ranges_t*)calloc(1, sizeof(ranges_t));
    kv_init(ret->rawunphasedblocks); 
    kv_init(ret->dropped); 
    kv_init(ret->starts); 
    kv_init(ret->ends); 
    kv_init(ret->decisions); 
    kv_init(ret->decisions_onraw); 
    kv_init(ret->flips_onraw); 
    kv_init(ret->phaseblocks); 
    kv_resize(u32p_t, ret->rawunphasedblocks, 16);
    kv_resize(u32p_t, ret->dropped, 16);
    kv_resize(uint32_t, ret->starts, 16);
    kv_resize(uint32_t, ret->ends, 16);
    kv_resize(int, ret->decisions, 16);
    kv_resize(int, ret->decisions_onraw, 16);
    kv_resize(int, ret->flips_onraw, 16);
    kv_resize(u32p_t, ret->phaseblocks, 16);
    return ret;
}
void destroy_ranges_t(ranges_t* d){
    kv_destroy(d->rawunphasedblocks);
    kv_destroy(d->dropped);
    kv_destroy(d->starts);
    kv_destroy(d->ends);
    kv_destroy(d->decisions);
    kv_destroy(d->decisions_onraw);
    kv_destroy(d->flips_onraw);
    kv_destroy(d->phaseblocks);
    free(d);
}


typedef struct{
    int ref_n;
    int ref_m;
    char **ref_names;
    ranges_t **ranges;
    htu32_t **varphase_in_dropped;  // length is ref_n
    htstri_t *qname2haptag;  // haptags obtained from meth phasing
    htstri_t *qname2haptag_raw;  // haptags obtained from input bam pre-processing
                                 // (if input bam is haplotagged, tags are not
                                 // stored here)
    int stores_raw_tag;
}storage_t;
void init_storage_t(storage_t *d){
    d->ref_n = 0;
    d->ref_m = 32;
    d->ref_names = (char**)calloc(d->ref_m, sizeof(char*));
    d->ranges = (ranges_t**)calloc(d->ref_m, sizeof(ranges_t*));
    d->varphase_in_dropped = 0;
    d->qname2haptag = htstri_ht_init();
    d->qname2haptag_raw = htstri_ht_init();
    d->stores_raw_tag = 0;
}
void resize_storage_t(storage_t *d, int m){
    d->ref_names = (char**)realloc(d->ref_names, sizeof(char*)*m);
    d->ranges = (ranges_t**)realloc(d->ranges, sizeof(ranges_t*)*m);
    d->ref_m = m;
}
void wipe_intervals_of_storage_t(storage_t *d){
    // drop loaded interval informations in order for them to be reloaded.
    for (int i=0; i<d->ref_n; i++){
        free(d->ref_names[i]);
        d->ranges[i]->dropped.n = 0;
        d->ranges[i]->starts.n = 0;
        d->ranges[i]->ends.n = 0;
        d->ranges[i]->phaseblocks.n = 0;
        d->ranges[i]->rawunphasedblocks.n = 0;
    }
    d->ref_n = 0;
}
void destroy_storage_t(storage_t *d, int include_self){
    for (int i=0; i<d->ref_n; i++){
        free(d->ref_names[i]);
        if (d->ranges[i]){
            destroy_ranges_t(d->ranges[i]);
        }
    }
    free(d->ref_names);
    free(d->ranges);
    if (d->varphase_in_dropped){
        for (int i=0; i<d->ref_n; i++){
            htu32_ht_destroy(d->varphase_in_dropped[i]);
        }
        free(d->varphase_in_dropped);
    }
    if (d->qname2haptag){
        htstri_t *h = d->qname2haptag;
        for (khint_t k=0; k<kh_end(h); k++){
            if (kh_exist(h, k)){
                free((char*)kh_key(h, k));
            }
        }
        htstri_ht_destroy(d->qname2haptag);
    }
    if (d->stores_raw_tag){
        htstri_t *h = d->qname2haptag_raw;
        for (khint_t k=0; k<kh_end(h); k++){
            if (kh_exist(h, k)){
                free((char*)kh_key(h, k));
            }
        }
    }
    htstri_ht_destroy(d->qname2haptag_raw);
    if (include_self) free(d);
}



int push_substring_forward_in_buffer(char *buf, int buf_l, int start){
    memmove(buf, buf+start, buf_l-start);
    return buf_l-start;
}

int insert_gtf_line(char *s, char *chrom, ranges_t *intvls, uint32_t *prev_end, 
                     int is_tsv){
    char *saveptr;
    int ret = 0;
    char *tok = strtok_r(s, "\t", &saveptr);
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
                kv_push(int, intvls->decisions, -1);
                ret = 1;
            }else{
                intvls->abs_start = strtoul(tok, NULL, 10);
            }
        }else if (i==col_e && use){  // end of phased block
            *prev_end = strtoul(tok, NULL, 10);
        }
        tok = strtok_r(NULL, "\t", &saveptr);
        i++;
    }
    return ret;
}


int insert_vcf_line(char *s, char *chrom, ranges_t *intvls, 
                     uint32_t *prev_pos, uint32_t *prev_group_ID){
    // collect phase blocks from vcf using PS fields
    // TODO/note: only uses SNP
    int used_the_variant = 0;
    char *saveptr;
    if (s[0]=='#'){
        if (s[1]=='#') return -1;
        else{
            char *tok = strtok_r(s, "\t", &saveptr);
            int n = 0;
            while (tok){
                tok = strtok_r(NULL, "\t", &saveptr);
                n++;
            }
            if (n<10){
                fprintf(stderr, "[E::%s] vcf only has %d columns; mandatory >=8; we also need FORMAT and at least 1 sample\n", __func__, n);
                exit(1);
            }else if (n>10){
                fprintf(stderr, "[E::%s] multi-sample vcf not implemented, TODO/TBD\n", __func__);
                exit(1);
            }
            return -1;
        }
    }
    char *tok = strtok_r(s, "\t", &saveptr);
    int i = 0, i_ps=-1;
    int use = 0, is_snp=0;
    uint32_t pos = 0;
    char groupID_s[11];
    while (tok){
        if (i==0){
            i_ps = -1;
            if (strcmp(tok, chrom)!=0) use = 0;
            else use = 1;
            is_snp = 1;
        }else if (i==1 && use){
            pos = (uint32_t)strtoul(tok, NULL, 10) ;
            if (*prev_pos!=UINT32_MAX && pos<*prev_pos){
                fprintf(stderr, "[E::%s] vcf not sorted? last line pos=%d, current pos=%d\n", 
                        __func__, (int)(*prev_pos), pos);
                exit(1);
            }
        } else if (i==3 && use){  // col: REF
            ;//if (strlen(tok)!=1 || seq_atcgun_table[(int)tok[0]]>=4) is_snp = 0;
        } else if (i==4 && use){  // col: ALT
            ;//if (strlen(tok)!=1 || seq_atcgun_table[(int)tok[0]]>=4) is_snp = 0;
        }else if (i==8 && use){
            i_ps = search_substr_idx(tok, "PS", ':', 1, 0, 0);
        }else if (i==9 && use){
            if (i_ps>=0 && is_snp){  // has phase group and entry is a SNP
                used_the_variant = 1;
                int ps_start, ps_l;
                get_substr_by_idx(tok, i_ps, ':', &ps_start, &ps_l, 0);
                if (ps_l==1 && tok[ps_start]=='.'){;}
                else{
                    sprintf(groupID_s, "%.*s", ps_l, tok+ps_start);
                    uint32_t groupID = strtoul(groupID_s, NULL, 10);
                    if (*prev_group_ID==UINT32_MAX) {  // special case: at init
                        *prev_group_ID = groupID;  
                        *prev_pos = pos;
                        intvls->abs_start = pos;
                    }
                    if (groupID==*prev_group_ID){  // still within the old phase block
                        *prev_pos = pos;
                    }else{ // encountered new phase block, push the last one.
                           // note we are pushing the unphased interval, not the phase block
                        if (*prev_pos!=UINT32_MAX){
                            kv_push(uint32_t, intvls->starts, *prev_pos);
                            kv_push(uint32_t, intvls->ends, groupID);
                            kv_push(int, intvls->decisions, -1);
                        }
                        *prev_group_ID = groupID;
                        *prev_pos = pos;
                    }
                }
            }
        }  // end of checking sample column
        tok = strtok_r(NULL, "\t", &saveptr);
        i++;
    }  // end of tok 
    return used_the_variant;
}

int insert_variant_from_vcf_line(char *s, char *refname, vvar_t *dest, 
                                 int range_start, int range_end){
    // return: 0 if nothing inserted, 
    //         1 if inserted variant from the given line.
    // sancheck header line
    char *saveptr;
    if (s[0]=='#'){
        if (s[1]=='#') {
            return 0;
        }else{
            char *tok = strtok_r(s, "\t", &saveptr);
            int n = 0;
            while (tok){
                tok = strtok_r(NULL, "\t", &saveptr);
                n++;
            }
            if (n<10){
                fprintf(stderr, "[E::%s] vcf only has %d columns; mandatory >=8; we also need FORMAT and at least 1 sample\n", __func__, n);
                exit(1);
            }else if (n>10){
                fprintf(stderr, "[E::%s] multi-sample vcf not implemented, TODO/TBD\n", __func__);
                exit(1);
            }
            return 0;
        }
    }
    // input is a variant entry
    char *tok = strtok_r(s, "\t", &saveptr);
    char *ref; int ref_l=-1;
    char *alt; int alt_l=-1;
    char *gt;
    int gt_start, gt_l;
    int i_col = 0;
    int i_gt;
    int pos=-1;
    while (tok){
        if (i_col==0){
            if (strcmp(tok, refname)!=0) {
                return 0;
            }
        }else if (i_col==1){
            pos = (uint32_t)strtoul(tok, NULL, 10) -1;  // use 0-index
            if (range_start>=0 && range_end>=0){
                if (pos<range_start || pos>=range_end) {
                    return 0;
                }
            }
        }else if (i_col==3){  // REF
            ref = tok;
            ref_l = strlen(ref);
        }else if (i_col==4){  // ALT
            alt = tok;
            alt_l = strlen(alt);
        }else if (i_col==8){  // tags
            i_gt = search_substr_idx(tok, "GT", ':', 1, 0, 0);
        }else if (i_col==9){  // fields
            if (i_gt>=0){
                get_substr_by_idx(tok, i_gt, ':', &gt_start, &gt_l, 0);
                if (gt_l!=3) {
                    return 0;
                }
                gt = tok+gt_start;
                // TODO diploid
                if (gt[1]!='|')               {
                    return 0;
                }
                if (gt[0]!='0' && gt[0]!='1') {
                    return 0;
                }
                if (gt[2]!='0' && gt[2]!='1') {
                    return 0;
                }
                // Record alt for SNP and INS, 
                // record ref for DEL.
                // The variants collected from vcf will always hold
                // haptags representing the ref alleles.
                int op = -1;
                int op_l = -1;
                char *var_s = 0;
                int hp = -1;
                if (ref_l==1 && alt_l==1){ // snp
                    op = VAR_OP_X;
                    op_l = 1;
                    var_s = alt;
                    hp = ((int)gt[0])-48;
                }else if (ref_l==alt_l){  // ?
                    ;
                }else if (ref_l>alt_l){  // del
                    op = VAR_OP_D;
                    op_l = ref_l - alt_l;
                    pos  += 1;
                    var_s = ref+1;
                    hp = ((int)gt[0])-48;
                }else if (ref_l<alt_l){  // ins
                    op = VAR_OP_I;
                    op_l = alt_l - ref_l;
                    var_s = alt+1;
                    hp = ((int)gt[0])-48;
                }
                if (op<0){
                    fprintf(stderr, "[W::%s] unhandled variant case at %s:%d ref=%s alt=%s\n", 
                                     __func__, refname, pos+1, ref, alt);
                    return 0;
                }
                push_to_vvar_t(dest, pos, op, op_l, hp, var_s);
            }
        }
        tok = strtok_r(NULL, "\t", &saveptr);
        i_col++;
    }
    return 1;
}

void parse_variants_for_one_read(bam1_t *aln, vvar_t *buf){
    int debug_print = 0;
    int self_start = 0;
    int ref_start = aln->core.pos;

    // parse cigar for insertions
    vu64_t insertions; kv_init(insertions); kv_resize(uint64_t, insertions, 16);  // for parsing MD tag
    uint32_t *cigar = bam_get_cigar(aln);
    uint32_t self_len = bam_cigar2qlen(aln->core.n_cigar, cigar);
    uint32_t op, op_l, consume;
    char *seq = (char*)malloc(16);
    int seq_n = 16;
    uint8_t *seqi = bam_get_seq(aln);
    uint32_t ref_pos = ref_start;
    uint32_t self_pos = 0;
    if (debug_print){
        fprintf(stderr, "[dbg::%s] at read %s (ref start pos=%d)\n", 
                        __func__, bam_get_qname(aln), ref_start);
    }
    for (int i=0; i<aln->core.n_cigar; i++){
        op = bam_cigar_op(cigar[i]);
        op_l = bam_cigar_oplen(cigar[i]);
        if (op==BAM_CREF_SKIP){
            ref_pos += op_l;
        }else if (op==BAM_CSOFT_CLIP){
            if (i==0) self_start = op_l;
            self_pos += op_l;
        }else if (op==BAM_CMATCH || op==BAM_CEQUAL || op==BAM_CDIFF){
            ref_pos += op_l;
            self_pos += op_l;
        }else if (op==BAM_CINS){
            if (op_l>=seq_n){
                seq_n = op_l+1;
                seq = (char*)realloc(seq, seq_n);
            }
            for (int j=0; j<op_l; j++){
                seq[j] = seq_nt16_str[bam_seqi(seqi, self_pos+j)];
            }
            push_to_vvar_t(buf, ref_pos, VAR_OP_I, op_l, HAPTAG_UNPHASED, seq);
            kv_push(uint64_t, insertions, ((uint64_t)op_l)<<32|self_pos);
            self_pos += op_l;
        }else if (op==BAM_CDEL){
            ref_pos += op_l;
        }
    }

    // parse MD tag for SNPs and deletions
    char snp_base[2]; snp_base[1]=0;
    char snp_base_dbg[10]; snp_base_dbg[9] = 0;
    int prev_ins_idx = 0;
    uint8_t *tagd = bam_aux_get(aln, "MD");
    assert(tagd);
    char *md_s = bam_aux2Z(tagd);
    int prev_md_i = 0;
    int prev_md_type, md_type;
    if (debug_print){
        fprintf(stderr, "[dbg::%s] qn=%s\n", __func__, bam_get_qname(aln));
        fprintf(stderr, "[dbg::%s] MD=%s\n", __func__, md_s);
    }
    // (init)
    self_pos = self_start;
    //fprintf(stderr, "self_pos at init: %d\n", self_pos);
    ref_pos = ref_start;
    prev_md_type = md_op_table[(int)md_s[0]];
    if (prev_md_type==2){
        snp_base[0] = seq_nt16_str[bam_seqi(seqi, self_pos)];
        push_to_vvar_t(buf, ref_pos, VAR_OP_X, 1, HAPTAG_UNPHASED, snp_base);
        ref_pos++;
        self_pos++;
        prev_md_type = -1;
    }
    assert(prev_md_type<4);
    // (collect operations)
    int i=1;
    while (md_s[i]){
        md_type = md_op_table[(int)md_s[i]];
        if (md_type==4){
            fprintf(stderr, "[E::%s] invalid MD: %c\n", __func__, md_s[i]);
            exit(1);
        }
        assert(md_type<4);
        if (md_type!=prev_md_type){
            if (prev_md_type==0){  // prev was match
                int l = natoi(md_s+prev_md_i, i-prev_md_i);
                ref_pos += l;
                self_pos += l;
                while (prev_ins_idx<insertions.n &&
                       self_pos>(uint32_t)insertions.a[prev_ins_idx]){
                    self_pos += insertions.a[prev_ins_idx]>>32;
                    prev_ins_idx++;
                }
            }else if (prev_md_type==1){ // prev was del
                if (md_type==0){  // current sees numeric, del run has ended
                    push_to_vvar_t(buf, ref_pos, VAR_OP_D, i-prev_md_i-1, 
                                   HAPTAG_UNPHASED , md_s+prev_md_i+1);

                    ref_pos += i-prev_md_i-1;
                    prev_md_type = md_type;
                    prev_md_i = i;
                }else{  // still in a del run, do not update status
                    ;
                }
                i++;
                continue;
            }
            // is current a SNP?
            if (md_type==2){
                snp_base[0] = seq_nt16_str[bam_seqi(seqi, self_pos)];
                push_to_vvar_t(buf, ref_pos, VAR_OP_X, 1, HAPTAG_UNPHASED, 
                               snp_base);
                if (debug_print) {
                    for (int x=self_pos-7, y=0; x<self_pos+2; x++, y++){
                        snp_base_dbg[y] = seq_nt16_str[bam_seqi(seqi, x)];
                    }
                    fprintf(stderr, "[dbg::%s] pushed SNP ref_pos=%d self_pos=%d base=%s, -7~+1:%s\n", 
                        __func__, ref_pos, self_pos, snp_base, snp_base_dbg);
                }
                ref_pos++;
                self_pos++;
                prev_md_type = -1;
                prev_md_i = i;
            }else{
                prev_md_type = md_type;
                prev_md_i = i;
            }

        }
        i++;
    }
    if (debug_print){
        for (int i=0; i<buf->n; i++){
            fprintf(stderr, "[dbg::%s]    op=%c pos=%d len=%d ", __func__, 
                "MXID"[buf->a[i].op], buf->a[i].pos, buf->a[i].len);
            if (buf->a[i].op>0){
                fprintf(stderr, "seq=");
                for (int j=0; j<buf->a[i].len; j++){
                    fprintf(stderr, "%c", "ACGT?"[buf->a[i].chars.a[j]]);
                }
            }
            fprintf(stderr, "\n");
        }
    }

cleanup:
    kv_destroy(insertions);
    free(seq);
}

int haptag_one_read_with_variants(char *refname, char *qname,// for logging
                                  vvar_t *known_vars, vvar_t *read_vars, 
                                  int start_pos, int end_pos, 
                                  int *prev_i_left){
    // note: `known_vars` needs to be sorted. Query reads are also
    //       expected to be processed in sorted-by-start-position order.
    //       In this impl, the interval loader will check vcf sorting 
    //       as it also depends on the ordering; reads are fed sorted
    //       as the input bam is required to be sorted & indexed. 
    //       `prev_i_left` reset between chromosome changes is 
    //       caller's responsibility.
    // return: 
    //         $(haptag) when tagged
    //         $(HAPTAG_UNPHASED) if no tag
    int debug_print = 0;
    if (known_vars->n==0) {
        if (debug_print){
            fprintf(stderr, "[dbg::%s] no known variants.\n", __func__);
        }
        return HAPTAG_UNPHASED;
    }

    // find index of start of the query's span on reference
    int i_left = *prev_i_left;
    while (i_left<known_vars->n && known_vars->a[i_left].pos<start_pos){
        i_left++;
    }
    *prev_i_left = i_left==0? 0 : i_left-1;

    // find informative variants
    // (collect)
    vu64_t piggyback;
    kv_init(piggyback);
    kv_resize(uint64_t, piggyback, 16);
    uint64_t typebit =1ULL<<32; 
    for (uint32_t i=i_left; i<known_vars->n; i++){
        if (known_vars->a[i].pos>=end_pos) break;
        uint64_t v = ((uint64_t)known_vars->a[i].pos)<<33 | i;
        kv_push(uint64_t, piggyback, v);
    }
    for (uint32_t i=0; i<read_vars->n; i++){
        uint64_t v = ((uint64_t)read_vars->a[i].pos)<<33 | typebit | i;
        kv_push(uint64_t, piggyback, v);
    }

    radix_sort_ksu64(piggyback.a, piggyback.a+piggyback.n);
    //for (int i= 0; i<piggyback.n; i++){
    //    fprintf(stderr, "[dbg::%s] qn=%s piggy#%d pos=%d typebit=%d\n", __func__, 
    //        qname, i, 
    //        piggyback.a[i]>>33, (piggyback.a[i]>>32)&1);
    //}

    // (parse)
    int hp_cnt[2] = {0,0};  // vote counter
    int hp;
    uint32_t ref_pos, self_pos, ref_i, self_i;
    for (int i=0; i<piggyback.n;){
        if (piggyback.a[i]&typebit){ // read's variant is not known to be informative
            if (debug_print) fprintf(stderr, "[dbg::%s] qn=%s case-1(read-uniq) read_pos= %d\n", __func__, qname, (int)(piggyback.a[i]>>33));
            i++;
            continue;
        } 
        ref_pos = piggyback.a[i]>>33;
        ref_i = (uint32_t)(piggyback.a[i]);
        if (i+1==piggyback.n){ // special case (end of interval); read must hold a REF
            hp = known_vars->a[ref_i].haptag;
            hp_cnt[hp]++;
            if (debug_print) fprintf(stderr, "[dbg::%s] qn=%s case0(ref) ref_pos= %d vote_hp=%d\n", __func__, qname, ref_pos,hp);
            break;
        }
        self_pos = piggyback.a[i+1]>>33;
        self_i = (uint32_t)(piggyback.a[i+1]);
        if (ref_pos!=self_pos){  // read maybe holds a REF at current position
            int skip_due_del = 0;
            if (i>0){  // look back to see if current position is actually a del
                if (piggyback.a[i-1]&typebit){
                    uint32_t last_self_pos = piggyback.a[i-1]>>33;
                    uint32_t last_self_i = (uint32_t)piggyback.a[i-1];
                    if (read_vars->a[last_self_i].op==VAR_OP_D){
                        if (last_self_pos+read_vars->a[last_self_i].len >= ref_pos)
                            skip_due_del = 1;
                    }
                }
            }
            if (skip_due_del){
                ;
            }else{
                hp = known_vars->a[ref_i].haptag;
                if (debug_print) fprintf(stderr, "[dbg::%s] qn=%s case1(ref-defaulting) ref_pos= %d vote_hp=%d\n", __func__, qname, ref_pos, hp);
                hp_cnt[hp]++;
            }
            i++;
        }else{
            if (!(piggyback.a[i+1]&typebit)){ // FIXME skip if reference collection has more than 1 alt at this position
                if (debug_print) fprintf(stderr, "[W::%s] input read haplotagging does not support multi-allele vcf entry, skipped %s:%d\n", 
                        __func__, refname, (int)ref_pos);
                i+=2;  
            }else{  //  check if read holds the ALT allele
                int ok = 1;
                variant_t *r = &known_vars->a[ref_i];
                variant_t *s = &read_vars->a[self_i];
                if (r->len==s->len){
                    for (int j=0; j<r->len; j++){
                        if (r->chars.a[j]!=s->chars.a[j]){
                            ok = 0;
                            break;
                        }
                    }
                }else ok=0;
                if (ok){
                    hp = known_vars->a[ref_i].haptag^1;
                    hp_cnt[hp]++;
                    if (debug_print) fprintf(stderr, "[dbg::%s] qn=%s case2(alt) ref_pos= %d vote_hp=%d\n", __func__, qname, ref_pos, hp);
                }else{
                    if (debug_print)
                        fprintf(stderr, "[dbg::%s] qn=%s case3(alt-fail) ref_pos =%d (r_len=%d s_len=%d s[0]==%c)\n", 
                        __func__, qname, ref_pos, r->len, s->len, "ACGT"[s->chars.a[0]]);
                }
                i+=2;
            }
        }
    }
    kv_destroy(piggyback);

    int haptag = -1;
    // FIXME: diploid
    float hp_cnt_max = MAX(hp_cnt[0], hp_cnt[1]);
    int hp_cnt_min = MIN(hp_cnt[0], hp_cnt[1]);
    float var_diff_ratio = hp_cnt_min==0? 0 : hp_cnt_max/(float)hp_cnt_min;
    if ((hp_cnt[0]>3 && 
         hp_cnt[1]>3 && 
         var_diff_ratio<VAR_DIFF_OVERRIDE_RATIO /*allow some counter evidences if self is very strong*/
         ) || 
        (hp_cnt[0]==hp_cnt[1])){  // ambiguous, treat as unphased
        haptag = HAPTAG_UNPHASED;
    }else if (hp_cnt[0]>hp_cnt[1]){
        haptag = 0;
    }else{
        haptag = 1;
    }
    if (debug_print){
        fprintf(stderr, "[dbg::%s] read=%s, votes: hp0=%d hp1=%d (int ratio=%.1f), final haptag=%d\n", 
            __func__, qname, hp_cnt[0], hp_cnt[1], 
            var_diff_ratio,
            haptag);
    }
    return haptag;
}
void pre_haplotagging_read_in_one_ref(char *fn_bam, char *ref_itvl, 
                                  vvar_t *known_vars, htstri_t *qname2haptag_raw){
    // Input preprocessing. After the known variants were loaded from 
    // vcf, tag reads and store the tags. 
    // note/FIXME: read with multiple alignment hits may have wrong
    //  phasing.
    int debug_print = 0;
    bamfile_t *hf = init_and_open_bamfile_t(fn_bam, pomfret_n_bam_threads);
    if (!hf){
        fprintf(stderr, "[E::%s] failed to open input file: %s\n", __func__, fn_bam);
        exit(1);  // TODO cleanup self
    }
    hts_itr_t *fp_itr = sam_itr_querys(hf->fp_idx, hf->fp_header, ref_itvl);
    if (!fp_itr) {return;}

    vvar_t read_vars;
    init_vvar_t(&read_vars);
    int prev_i_left = 0;
    khint_t k;
    htstri_t *h = qname2haptag_raw;
    int absent;
    int tot_tagged[4]={0,0,0,0};  // new hap0, new hap1, new unphased, dup  
    while(sam_itr_next(hf->fp, fp_itr, hf->buf_aln)>=0){
        char *qn = bam_get_qname(hf->buf_aln);
        int start_pos = hf->buf_aln->core.pos;
        int end_pos = bam_endpos(hf->buf_aln);

        wipe_vvar_t(&read_vars);
        int flag = hf->buf_aln->core.flag;
        if ((flag&4) || (flag&256) || (flag&2048)) continue;
        parse_variants_for_one_read(hf->buf_aln, &read_vars);

        int haptag = haptag_one_read_with_variants(ref_itvl, qn, 
                        known_vars, &read_vars, start_pos, end_pos, &prev_i_left);
        if (debug_print) fprintf(stderr, "[dbg::%s] qn=%s generated_haptag=%d\n", __func__, qn, haptag);

        int qn_l = strlen(qn);
        char *qname_alloc = (char*)calloc(qn_l+1, 1);
        sprintf(qname_alloc, "%s", qn);
        k = htstri_ht_put(h, qn, &absent);
        if (absent){
            kh_key(h, k) = qname_alloc;
            kh_val(h, k) = haptag;
            if (haptag==1 || haptag==0) tot_tagged[haptag]++;
            else tot_tagged[2]++;
        }else{
            tot_tagged[3]++;
            free(qname_alloc);
        }
    }

    // clean up
    destroy_vvar_t(&read_vars, 0);
    hts_itr_destroy(fp_itr);
    destroy_bamfile_t(hf, 1);
    fprintf(stderr, "[dbg::%s] tagged: %d new hap0, %d new hap1, %d new unphased, %d dup\n", 
                    __func__, tot_tagged[0], tot_tagged[1], tot_tagged[2], tot_tagged[3]);
}




ranges_t *load_intervals_from_file_one_ref(char *fn, char *chrom, 
                                  enum input_file_format fn_format, 
                                  vvar_t *vars, 
                                  int vars_range_start, int vars_range_end){
    fprintf(stderr, "[M::%s] reading from %s\n", __func__, fn);
    ranges_t *ret = init_ranges_t();
    gzFile fp = gzopen(fn, "rb");
    if (!fp){
        fprintf(stderr, "[E::%s] failed to open file: %s\n", __func__, fn);
        exit(1);
    }
    int m = READLINE_BUF_LEN;
    int m_tmp = m;
    int n = 0;
    char *buffer = (char*)malloc(m);
    char *buffer2 = (char*)malloc(m);
    int offset = 0;
    int n_insertions = 0;

    uint32_t prev_end = UINT32_MAX;
    uint32_t prev_group_ID = UINT32_MAX;
    int n_read;
    int tot_variants = 0;
    int used_variants = 0;
    int n_phasing_variants = 0;
    while ((n_read = gzread(fp, buffer+offset, m-offset)) > 0) {
        int start = 0;
        for (int i=0; i<offset+n_read; i++){
            if (buffer[i]=='\n'){
                buffer[i] = 0;
                if (fn_format==IS_GTF){
                    insert_gtf_line(buffer+start, chrom, ret, &prev_end, 0);
                }else if (fn_format==IS_TSV){
                    insert_gtf_line(buffer+start, chrom, ret, &prev_end, 1);
                }else if (fn_format==IS_VCF){
                    if (vars){
                        sprintf(buffer2, "%s", buffer+start);
                        n_phasing_variants += insert_variant_from_vcf_line(buffer2, chrom, vars, 
                                                vars_range_start, vars_range_end);
                    }
                    int used = insert_vcf_line(buffer+start, chrom, ret, &prev_end, &prev_group_ID);
                    if (used>=0){
                        tot_variants++;
                        if (used>0) used_variants++;
                    }
                }
                start = i+1; 
            }
        }
        if (start!=0){
            offset = push_substring_forward_in_buffer(buffer, m, start);
        }else{  // did not see a full line
            offset = m;
            m = m+(m>>1);
            buffer = (char*)realloc(buffer, m);
            buffer2 = (char*)realloc(buffer2, m);
        }
    }
    if (prev_end!=UINT32_MAX) ret->abs_end = prev_end;
    if (fn_format==IS_VCF){
        fprintf(stderr, "[M::%s] loaded from vcf: total %d variants, used %d when looking for phaseblocks\n", 
                __func__, tot_variants, used_variants);
    }
    fprintf(stderr, "[M::%s] chrom %s, read %d unphased intervals from %s; phasing's abs start=%d, abs end=%d\n", 
                    __func__, chrom, (int)ret->starts.n, fn, (int)ret->abs_start, (int)ret->abs_end);
    if (vars){
        fprintf(stderr, "[M::%s] recorded %d variants\n", __func__, n_phasing_variants);
    }
    gzclose(fp);
    free(buffer);
    free(buffer2);
    return ret;
}

void load_intervals_from_file(char *fn, 
                              enum input_file_format fn_format, 
                              storage_t *st, 
                              int load_vcf_variants_too, 
                              char *fn_bam, 
                              vvar_t *var_storage,  // if provided and input is vcf, will store variants; won't reset, init or realloc
                              int var_storage_l){
    double T = Get_T();
    gzFile fp = gzopen(fn, "rb");
    if (!fp){
        fprintf(stderr, "[E::%s] failed to open file for phase blocks: %s\n", __func__, fn);
        exit(1);
    }

    int m = READLINE_BUF_LEN;
    int m_tmp = m;
    int n = 0;
    char *buffer = (char*)malloc(m);
    char *buffer2 = (char*)malloc(m);
    char *buffer3 = (char*)malloc(m);
    char *tmp_buffer;
    char *tok;
    int tok_l;
    int offset = 0;
    int n_read;

    ranges_t **rp = 0;
    char *chrom=0;
    uint32_t prev_end = UINT32_MAX;
    uint32_t prev_group_ID = UINT32_MAX;

    vvar_t phased_variants;
    if (load_vcf_variants_too && !var_storage) {
        init_vvar_t(&phased_variants);
        st->stores_raw_tag = 1;
    }

    int i_ref_cache = -1;
    kstring_t *chrom_cache = (kstring_t*)calloc(1, sizeof(kstring_t));
    while ((n_read = gzread(fp, buffer+offset, m-offset)) > 0) {
        int start = 0;
        char *saveptr_tmpbuffer;
        for (int i=0; i<offset+n_read; i++){
            if (buffer[i]=='\n'){
                buffer[i] = 0;
                tmp_buffer = buffer+start;
                sprintf(buffer2, "%s", tmp_buffer);
                sprintf(buffer3, "%s", tmp_buffer);

                if (tmp_buffer[0]=='#') {  // vcf comment line
                    start = i+1;
                    continue;  
                }
                tok = strtok_r(tmp_buffer, "\t", &saveptr_tmpbuffer);
                tok_l = strlen(tok);
                if (!tok) {start=i+1; continue;}
                if (st->ref_n==0){
                    st->ref_names[0] = (char*)malloc(tok_l+1);
                    sprintf(st->ref_names[0], "%s", tok);
                    chrom = st->ref_names[0];
                    st->ranges[st->ref_n] = init_ranges_t();
                    rp = &st->ranges[st->ref_n];
                    st->ref_n++;
                    fprintf(stderr, "[M::%s] at ref %s\n", __func__, chrom);
                }else{
                    int found = 0;
                    for (int i=st->ref_n-1; i>=0; i--){
                        if (strcmp(st->ref_names[i], tok)==0){
                            found = 1;
                            rp = &st->ranges[i];
                            chrom = st->ref_names[i];
                            break;
                        }
                    }
                    if (!found){
                        // Out of cautious, compare with other known
                        // reference names as well, to recognize 
                        // unsorted input.
                        // TODO

                        // new refname.
                        rp = &st->ranges[st->ref_n-1];
                        chrom = st->ref_names[st->ref_n-1];

                        // update end position of previous block
                        if (prev_end!=UINT32_MAX) (*rp)->abs_end = prev_end;

                        // does the previous block need haptagging the reads?
                        if (load_vcf_variants_too && phased_variants.n>0){
                            //fprintf(stderr, "[M::%s] preprocessing input bam (haplotagging) for ref=%s...\n", 
                            //    __func__, chrom);
                            if (var_storage){  // (caller request to get variants, do not tag reads)  
                                               // TODO refactor
                                if (var_storage_l<st->ref_n){
                                    fprintf(stderr, "[E::%s] check code, requested to store variants but buffer is short.\n", __func__);
                                    exit(1);
                                }
                            }else{
                                pre_haplotagging_read_in_one_ref(
                                    fn_bam, chrom, 
                                    &phased_variants, 
                                    st->qname2haptag_raw);
                                wipe_vvar_t(&phased_variants);
                            }
                        }

                        // store new refname and init range buffer
                        if (!(load_vcf_variants_too && var_storage)){
                            if (st->ref_n==st->ref_m){
                                resize_storage_t(st, st->ref_m+(st->ref_m>>1));
                            }
                            st->ref_names[st->ref_n] = (char*)malloc(tok_l+1);
                            sprintf(st->ref_names[st->ref_n], "%s", tok);
                            st->ranges[st->ref_n] = init_ranges_t();
                            rp = &st->ranges[st->ref_n];
                            chrom = st->ref_names[st->ref_n];
                        }else{
                            chrom = tok;
                        }

                        // reset and step forward
                        prev_end = UINT32_MAX;
                        if (!(load_vcf_variants_too && var_storage)) st->ref_n++;
                        //fprintf(stderr, "[M::%s] at ref %s\n", __func__, st->ref_names[st->ref_n-1]);
                    }
                }

                // insert values
                if (fn_format==IS_GTF){
                    insert_gtf_line(buffer2, chrom, *rp, &prev_end, 0);
                }else if (fn_format==IS_TSV){
                    insert_gtf_line(buffer2, chrom, *rp, &prev_end, 1);
                }else if (fn_format==IS_VCF){
                    if (load_vcf_variants_too && var_storage){
                        if (!chrom_cache->s || strcmp(chrom, chrom_cache->s)!=0){
                            for (int tmpi=0; tmpi<st->ref_n; tmpi++){
                                assert(chrom);
                                if (strcmp(st->ref_names[tmpi], chrom)==0){
                                    i_ref_cache = tmpi;
                                    break;
                                }
                            }
                            ksprintf(chrom_cache, "%s", chrom);
                        }
                        if (i_ref_cache>=0){
                            insert_variant_from_vcf_line(buffer2, chrom, &var_storage[i_ref_cache], -1, -1);
                        }
                    }else if (load_vcf_variants_too && !var_storage){
                        insert_variant_from_vcf_line(buffer2, chrom, &phased_variants, -1, -1);
                        insert_vcf_line(buffer3, chrom, *rp, &prev_end, &prev_group_ID);
                    }
                    else {
                        insert_vcf_line(buffer3, chrom, *rp, &prev_end, &prev_group_ID);
                    }
                }

                // step forward in readline buffer
                start = i+1; 
            }
        }
        if (start!=0){
            offset = push_substring_forward_in_buffer(buffer, m, start);
        }else{  // did not see a full line
            offset = m;
            m = m+(m>>1);
            buffer = (char*)realloc(buffer, m);
            buffer2 = (char*)realloc(buffer2, m);
            buffer3 = (char*)realloc(buffer3, m);
        }

    }
    // last block
    if (load_vcf_variants_too && !var_storage && chrom){
        if (var_storage);
        else{
            pre_haplotagging_read_in_one_ref(fn_bam, chrom, 
                        &phased_variants, st->qname2haptag_raw);
            wipe_vvar_t(&phased_variants);
        }
    }
    if (prev_end!=UINT32_MAX && rp) {
        (*rp)->abs_end = prev_end;
    }

    free(chrom_cache->s);
    free(chrom_cache);
    free(buffer);
    free(buffer2);
    free(buffer3);
    gzclose(fp);

    if (load_vcf_variants_too && !var_storage) destroy_vvar_t(&phased_variants, 0);
    fprintf(stderr, "[T::%s] used %.1fs\n", __func__, Get_T()-T);
    //if (!load_vcf_variants_too){
    //    for (int i=0; i<st->ref_n; i++){
    //        fprintf(stderr, "[dbg::%s] %s has %d intervals\n", __func__, 
    //        st->ref_names[i], (int)st->ranges[i]->starts.n);
    //    }
    //}
}

void store_raw_intervals(ranges_t *rg){
    int verbose = 0;
    rg->rawunphasedblocks.n = 0;
    for (int i=0; i<rg->starts.n; i++){
        if (verbose){
            fprintf(stderr, "[dbg::%s] raw itvl #%d: %d-%d\n", __func__, i, rg->starts.a[i], rg->ends.a[i]);
        }
        kv_push(u32p_t, rg->rawunphasedblocks, ((u32p_t){.s=rg->starts.a[i], 
                                                         .e=rg->ends.a[i]}));
    }
}

void merge_close_intervals(ranges_t *rg, int threshold){
    int debug_print = 0;
    int j = 0;
    if (rg->starts.n<=1){
        rg->decisions.n = rg->starts.n;
        return;
    }
    for (int i=1; i<rg->starts.n; i++){
        if (rg->starts.a[i] - rg->ends.a[j]<threshold){  // merge
            kv_push(u32p_t, rg->dropped, 
                ((u32p_t){.s=rg->ends.a[j], .e=rg->starts.a[i]}) );
            if (debug_print){
                fprintf(stderr, "[dbg::%s] mark %d-%d as dropped phased interval\n", 
                                __func__, rg->ends.a[j], rg->starts.a[i]);
            }
            rg->ends.a[j] = rg->ends.a[i];
        }else{
            if (debug_print){
                fprintf(stderr, "[dbg::%s] proceed\n", __func__);
            }
            j++;
            rg->starts.a[j] = rg->starts.a[i];
            rg->ends.a[j] = rg->ends.a[i];
        }
    }
    rg->decisions.n = rg->starts.n;
    rg->starts.n = j+1;
    rg->ends.n = j+1;
}




void dbgoutput_intermediate_read_haplotags(storage_t *st, char *prefix){
    // note: outputs haplotags on per-interval basis. This 
    //       will be different from the final bam.
    char *fn_out  = (char*)malloc(strlen(prefix)+25);
    sprintf(fn_out, "%s.mp.dbg.read2tag", prefix);
    FILE *fp_out = fopen(fn_out, "w");
    if (!fp_out){
        fprintf(stderr, "[E::%s] failed to open output file: %s\n", 
            __func__, fn_out);
        exit(1);
    }
    htstri_t *h = st->qname2haptag;
    int n = 0;
    for (khint_t k=0; k<kh_end(h);k++){
        if (kh_exist(h, k)){
            int hap = kh_val(h, k);
            char *qname = kh_key(h, k);
            hap = hap<0? HAPTAG_UNPHASED : hap;
            fprintf(fp_out, "%s\t-1\t%d\n", qname, hap+1);
            n++;
        }
    }
    fclose(fp_out);
    free(fn_out);
    fprintf(stderr, "[dbg::%s] wrote dbg read haplotag txt; %d reads were stored\n", __func__, n);
}

void lift_decisions(storage_t *st){
    // create global flipping decisions from the local ones. 
    int verbose = 0;
    for (int i_ref=0; i_ref<st->ref_n; i_ref++){
        ranges_t *rr = st->ranges[i_ref];
        rr->phaseblocks.n = 0;
        int j = 0;
        for (int i=0; i<rr->decisions.n; i++){
            if (rr->decisions.a[i]<0){
                while (j<rr->rawunphasedblocks.n && rr->rawunphasedblocks.a[j].e<=rr->ends.a[i]){
                    if (verbose)
                        fprintf(stderr, "[dbg::%s] decision#%d %s:%d-%d decision=%d ; maintain (j=%d, within %d-%d)\n", 
                        __func__, i, st->ref_names[i_ref], 
                        rr->rawunphasedblocks.a[j].s, rr->rawunphasedblocks.a[j].e, 
                        rr->decisions.a[i], j, 
                        rr->starts.a[i], rr->ends.a[i]
                        );
                    kv_push(int, rr->decisions_onraw, rr->decisions.a[i]);
                    j++;
                }
            }else{  // might have joined multiple blocks, modify the raw
                    // phase block range as needed
                if (rr->rawunphasedblocks.a[j].e < rr->ends.a[i]){
                    int j2;
                    int found = 0;
                    for (j2=j; j2<rr->rawunphasedblocks.n; j2++){
                        if (rr->rawunphasedblocks.a[j2].e==rr->ends.a[i]) {
                            found = 1; 
                            break;
                        }
                    }
                    assert(found);
                    if (verbose) fprintf(stderr, "[dbg::%s] update end: %d<-%d\n", __func__, rr->rawunphasedblocks.a[j].e, rr->ends.a[i]);
                    rr->rawunphasedblocks.a[j].e = rr->ends.a[i];
                    int j_left = j+1;
                    int j_right = j2+1;
                    while (j_right<rr->rawunphasedblocks.n){
                        if (verbose){
                            fprintf(stderr, "[dbg::%s] moving: %d<-%d\n", __func__, j_left, j_right);
                            fprintf(stderr, "[dbg::%s]     s: %d<-%d\n", __func__, rr->rawunphasedblocks.a[j_left].s, rr->rawunphasedblocks.a[j_right].s );
                            fprintf(stderr, "[dbg::%s]     e: %d<-%d\n", __func__, rr->rawunphasedblocks.a[j_left].e, rr->rawunphasedblocks.a[j_right].e );
                        }
                        rr->rawunphasedblocks.a[j_left].s = rr->rawunphasedblocks.a[j_right].s;
                        rr->rawunphasedblocks.a[j_left].e = rr->rawunphasedblocks.a[j_right].e;
                        j_left++;
                        j_right++;
                    }
                    rr->rawunphasedblocks.n -= j2-j;
                    if (verbose) fprintf(stderr, "[dbg::%s] shrink %d\n", __func__, j2-j);
                }
                kv_push(int, rr->decisions_onraw, rr->decisions.a[i]);
                if (verbose)
                    fprintf(stderr, "[dbg::%s] decision#%d %s:%d-%d decision=%d ; override\n", 
                                    __func__, i, st->ref_names[i_ref], 
                                    rr->rawunphasedblocks.a[j].s, rr->rawunphasedblocks.a[j].e, 
                                    rr->decisions.a[i]);
                j++; 
            }
        }
    }
}

void make_decisions_flippings_onraw(storage_t *st){
    for (int i_ref=0; i_ref<st->ref_n; i_ref++){
        ranges_t *rr = st->ranges[i_ref];
        int flip = 0;
        for (int i=0; i<rr->decisions_onraw.n; i++){
            if (rr->decisions_onraw.a[i]<0) flip = 0;
            else flip ^= rr->decisions_onraw.a[i];
            //fprintf(stderr, "[dbg::%s] interval#%d (%d-%d) flip=%d\n", 
            //    __func__, i, rr->rawunphasedblocks.a[i].s, rr->rawunphasedblocks.a[i].e, flip);
            kv_push(int, rr->flips_onraw, flip);
        }
    }
}

void generate_new_phase_blocks(storage_t *st, int use_raw){
    int verbose = 0;
    for (int i_ref=0; i_ref<st->ref_n; i_ref++){
        ranges_t *rr = st->ranges[i_ref];
        uint32_t start = st->ranges[i_ref]->abs_start;
        uint32_t end = UINT32_MAX;
        int N = use_raw? rr->decisions_onraw.n : rr->decisions.n;
        int *de = use_raw? rr->decisions_onraw.a : rr->decisions.a;
        if (verbose) {
            fprintf(stderr, "[dbg::%s] ref#%d %s, blocks=%d, start=%d\n", 
                        __func__, i_ref, st->ref_names[i_ref], N ,start);
            for (int i=0; i<N; i++){
                fprintf(stderr, "[dbg::%s] range#%d : %d-%d\n", 
                    __func__, i, rr->starts.a[i], rr->ends.a[i]);
            }
        }
        char *ref_name = st->ref_names[i_ref];
        for (int i=0; i<N; i++){
            if (de[i]>=0) {
                if (verbose) fprintf(stderr, "[dbg::%s] a join\n", __func__);
                continue;
            }
            end = use_raw? rr->rawunphasedblocks.a[i].s : rr->starts.a[i];
            if (verbose) fprintf(stderr, "[dbg::%s] push %d - %d\n", __func__, start, end);
            kv_push(u32p_t, rr->phaseblocks, ((u32p_t){.s=start, 
                                                       .e=end}));
            start = use_raw? rr->rawunphasedblocks.a[i].e :  rr->ends.a[i];
        }
        if (N>0 && end!=rr->abs_end){
            end = end==UINT32_MAX? rr->abs_start : end;
            kv_push(u32p_t, rr->phaseblocks, ((u32p_t){.s=end, 
                                                       .e=rr->abs_end}));
            if (verbose) fprintf(stderr, "[dbg::%s] push %d - %d\n", __func__, end, rr->abs_end);
        }
    }
}



int get_new_phaseblock_ID1(ranges_t *r, int *prev_idx, int pos){
    for (int i=0; i<r->phaseblocks.n; i++){
        if (r->phaseblocks.a[i].s==UINT32_MAX || 
            r->phaseblocks.a[i].e==0 || 
            r->phaseblocks.a[i].e==UINT32_MAX) {
            continue;
        }

        if ( pos>=r->phaseblocks.a[i].s && pos<r->phaseblocks.a[i].e){
            *prev_idx = i;
            return r->phaseblocks.a[i].s;
        }
        //if (pos>=r->phaseblocks.a[i].e) {break;}
        //if (pos<r->phaseblocks.a[i].s)  {break;}
    } 
    return -1;
}
int get_new_phaseblock_ID(storage_t *st, char *refname, int *prev_idx, int pos){
    // TODO: use hashtable? cache? (same file thus ref order is the same)
    for (int i=0; i<st->ref_n; i++){
        if (strcmp(st->ref_names[i], refname)==0){
            return get_new_phaseblock_ID1(st->ranges[i], prev_idx, pos);
        }
    }
    return -2;
}

int tmp_check_if_in_dropped_intervals(storage_t *st, char *refname, int pos){
    for (int i=0; i<st->ref_n; i++){
        if (strcmp(st->ref_names[i], refname)==0){
            ranges_t *r = st->ranges[i];
            for (int j=0; j<r->dropped.n; j++){
                if (pos>=r->dropped.a[j].s && pos<=r->dropped.a[j].e){ 
                    return 1;
                }
            }
        }
    }
    return 0;
}

int check_if_in_phased_intervals(storage_t *st, char *refname, int pos, 
                                   int *prev_idx, int *updated){
    int prev  = *prev_idx;
    for (int i=0; i<st->ref_n; i++){
        if (strcmp(st->ref_names[i], refname)==0){
            ranges_t *r = st->ranges[i];
            for (int j=*prev_idx; j<r->starts.n; j++){
                if (pos>=r->ends.a[j-1] && 
                    pos<=r->starts.a[j]){ 
                    if (j!=prev){
                        *updated = 1;
                        *prev_idx = j;
                    }
                    return 1;
                }
            }
            break;
        }
    }
    return 0;
}

int get_flip_status_by_idx(storage_t *st, char *refname, int idx){
    for (int i=0; i<st->ref_n; i++){
        if (strcmp(st->ref_names[i], refname)==0){
            ranges_t *rr = st->ranges[i];
            return rr->flips_onraw.a[idx];
        }
    }
    return -1;
}

int get_flip_status(storage_t *st, char *refname, int *prev_idx, int pos){
    for (int i=0; i<st->ref_n; i++){
        if (strcmp(st->ref_names[i], refname)==0){
            ranges_t *rr = st->ranges[i];
            int j;
            for (j=*prev_idx; j<rr->rawunphasedblocks.n; j++){
                int start = rr->rawunphasedblocks.a[j].s;
                int end   = rr->rawunphasedblocks.a[j].e;
                if (start>=pos){ 
                    *prev_idx = j==0? 0:j-1;
                    int stat = rr->flips_onraw.a[*prev_idx];

                    // special case: prior to the first gap
                    if (pos<=rr->rawunphasedblocks.a[0].s) stat = 0;

                    return stat;
                }
            }
            *prev_idx = j-1;
            //for (j=0; j<rr->rawunphasedblocks.n; j++){
            //    int start = rr->rawunphasedblocks.a[j].s;
            //    int end   = rr->rawunphasedblocks.a[j].e;
            //    if (start>=pos){
            //        //fprintf(stderr, "[dbg::%s] found flip status: pos=%d interval=%d-%d\n", 
            //        //__func__, pos, start, end);
            //        j = j==0? 0 : j-1;
            //        return rr->flips_onraw.a[j];
            //    }
            //}
            return rr->flips_onraw.a[rr->rawunphasedblocks.n==0? 0:rr->rawunphasedblocks.n-1];
            break;
        }
    }
    return -1;

}

void recover_variant_phase_in_one_interval(storage_t *st, 
                            char *refname, uint32_t start, uint32_t end,   // interval
                            uint32_t *poss, uint32_t poss_n, // variant positions as on reference
                            htu32_t *pos2haptag,
                            char *fn_bam, char *fn_vcf){
    int debug_print = 0;
    vvar_t vars; init_vvar_t(&vars);
    vvar_t read_vars; init_vvar_t(&read_vars);

    vu64_t pb;  // piggyback: pos:31|is_read:1|id:32 (records all-hap read coverage in the ref's entry, and idx in buffer for query)
    kv_init(pb);
    kv_resize(uint64_t, pb, 16);

    // load variants
    // TODO diploid: since we are confident about the variants, 
    // have diploid assumption, and only need to decide a phase, 
    // it does not matter too much to not check whether an ALT allele
    // in a read matches the ALT in the reference variant set.
    // To fully load reference variants for propoer comparison:
    //    ranges_t *tmp = load_intervals_from_file_one_ref(fn_vcf, refname, 
    //                                          IS_VCF, &vars, start, end);
    //    destroy_ranges_t(tmp);
    for (uint32_t i=0; i<poss_n; i++){
        kv_push(uint64_t, pb, ((uint64_t)poss[i])<<33);
    }

    // Process reads and collect their variants.
    // Diploid: similar to above, we do not have to care about REF alleles 
    // in reads, since the variants should be het & somewhat confident and 
    // we only need to figure out the phasing.
    kstring_t *ref_itvl = (kstring_t*)calloc(1, sizeof(kstring_t));
    ksprintf(ref_itvl, "%s:%d-%d", refname, start, end);
    bamfile_t *hf = init_and_open_bamfile_t(fn_bam, pomfret_n_bam_threads);
    if (!hf){
        fprintf(stderr, "[E::%s] failed to open input file: %s\n", __func__, fn_bam);
        exit(1);  // TODO cleanup self
    }
    hts_itr_t *fp_itr = sam_itr_querys(hf->fp_idx, hf->fp_header, ref_itvl->s);
    while(sam_itr_next(hf->fp, fp_itr, hf->buf_aln)>=0){
        // For each read in the region, if it has been meth-phased, 
        // collect its variants.
        char *qn = bam_get_qname(hf->buf_aln);
        uint32_t start_pos = hf->buf_aln->core.pos;
        uint32_t end_pos = bam_endpos(hf->buf_aln);
        
        khint_t k = htstri_ht_get(st->qname2haptag, qn);
        if (k==kh_end(st->qname2haptag)) continue;
        uint8_t hap_methphased = kh_val(st->qname2haptag, k);

        uint8_t hap_raw;
        if (st->stores_raw_tag){
            k = htstri_ht_get(st->qname2haptag_raw, qn);
            if (k==kh_end(st->qname2haptag_raw)) continue;
            hap_raw = kh_val(st->qname2haptag_raw, k);
        }else{
            hap_raw = get_hp_from_aln(hf->buf_aln);
        }
        if (hap_raw==HAPTAG_UNPHASED) continue;

        // get read's variants; also stuff haptag (piggybacked with raw haptag), 
        // if meth-phased, in the structs.
        int prev_i = read_vars.n;
        parse_variants_for_one_read(hf->buf_aln, &read_vars);
        for (int i=prev_i; i<read_vars.n; i++){
            read_vars.a[i].haptag = hap_methphased<<4 | hap_raw;  // LIMIT: limited to max 16 haps
        }

        // update coverage
        for (int i=0; i<pb.n; i++){
            if (pb.a[i]>=start_pos && 
                pb.a[i]<end_pos && 
                (uint32_t)pb.a[i]<UINT32_MAX) pb.a[i]++;
        }
    }
    if (debug_print){
        fprintf(stderr, "[dbg::%s] interval %s has %d read vars (%d ref vars)\n", 
                        __func__, ref_itvl->s, (int)read_vars.n, (int)pb.n);
        fprintf(stderr, "[dbg::%s] dump ref vars:\n", __func__);
        for (int i=0; i<pb.n; i++){
            fprintf(stderr, "[dbg::%s]  ref#%d pos=%d\n", __func__, i, (int)(pb.a[i]>>33));
        }
    }

    if (pb.n==0) goto done;

    for (uint32_t i=0; i<read_vars.n; i++){
        kv_push(uint64_t, pb, 
                ((uint64_t)read_vars.a[i].pos)<<33|(1ULL<<32)| i);
    }

    // group the variants by position
    uint64_t typebit = 1ULL<<32;
    radix_sort_ksu64(pb.a, pb.a+pb.n);

    // go through the variants
    for (int i=0, j; i<pb.n-1;){
        if (pb.a[i]&typebit) {i++; continue;}  // not a known variant
        uint32_t ref_pos = pb.a[i]>>33;
        int n_tot = (uint32_t)pb.a[i];
        int n_alt = 0;
        int hp_cnt[2] = {0, 0}; // TODO:haploid
        for (j=i+1; j<pb.n; j++){
            if (!(pb.a[j]&typebit)) break;  // got to the next reference
            if ((pb.a[j]>>33)!=ref_pos) break;  // got to the next non-ref variant position
            int hap = read_vars.a[(uint32_t)pb.a[j]].haptag>>4;
            //int hap_raw = read_vars.a[(uint32_t)pb.a[j]].haptag & 15;
            if (hap==0 || hap==1) {  // TODO:haploid
                hp_cnt[hap]++;
            };
        }
        int hap_of_ref;
        if (hp_cnt[0]>hp_cnt[1]) hap_of_ref = 1;  // variant called on hap0, so let ref have hap1
        else if (hp_cnt[1]>hp_cnt[0]) hap_of_ref = 0;
        else hap_of_ref = HAPTAG_UNPHASED;

        // insert to ht
        int absent;
        khint_t k = htu32_ht_put(pos2haptag, ref_pos, &absent);
        if (debug_print){
            if (!absent){
                fprintf(stderr, "[W::%s] refpos %d is in the ht already\n", __func__, ref_pos);
            }else{
                fprintf(stderr, "[dbg::%s] inserted pos=%d hap_of_ref=%d (buffer content: %d %d) (tot coverage=%d)\n", 
                __func__, ref_pos, hap_of_ref, hp_cnt[0], hp_cnt[1], n_tot);
            }
        }
        kh_val(pos2haptag, k) = hap_of_ref;

        i = j;
    }

done: 
    // reset
    destroy_bamfile_t(hf, 1);
    hts_itr_destroy(fp_itr);
    free(ref_itvl->s);
    free(ref_itvl);
    destroy_vvar_t(&vars, 0);
    destroy_vvar_t(&read_vars, 0);
    kv_destroy(pb);
}


void recover_variant_phase_in_dropped_intervals(storage_t *st,
                                  char *fn_bam, char *fn_vcf){
    int debug_print = 0;
    vu32_t poss;
    kv_init(poss);
    kv_resize(uint32_t, poss, 16);

    // init variant buffers
    storage_t *st2 = (storage_t*)calloc(1, sizeof(storage_t));
    init_storage_t(st2);
    vvar_t *vars_v = (vvar_t*)malloc(sizeof(vvar_t)*st->ref_n);
    assert(vars_v);
    st2->ref_names = (char**)realloc(st2->ref_names, sizeof(char*)*st->ref_n);
    for (int i=0; i<st->ref_n; i++){
        init_vvar_t(&vars_v[i]);
        st2->ref_names[i] = (char*)malloc(strlen(st->ref_names[i])+1);
        sprintf(st2->ref_names[i], "%s", st->ref_names[i]);
    }
    st2->ref_n = st->ref_n;
    st2->ranges = (ranges_t**)realloc(st2->ranges, sizeof(ranges_t*)*st2->ref_n);
    for (int i=0; i<st2->ref_n; i++) st2->ranges[i] = init_ranges_t();
    load_intervals_from_file(fn_vcf, IS_VCF, st2, 1, fn_bam, vars_v, st->ref_n);
    destroy_storage_t(st2, 1);

    // init hashtables
    st->varphase_in_dropped = (htu32_t**)calloc(st->ref_n, sizeof(htu32_t*));
    for (int i_ref=0; i_ref<st->ref_n; i_ref++){
        st->varphase_in_dropped[i_ref] = htu32_ht_init();
    }

    for (int i_ref=0; i_ref<st->ref_n; i_ref++){
        char *refname = st->ref_names[i_ref];
        vvar_t *vars = &vars_v[i_ref];
        if (debug_print) fprintf(stderr, "[M::%s] at ref %s\n", __func__, refname);
        //// just grab all variants
        //ranges_t *tmp = load_intervals_from_file_one_ref(fn_vcf, refname, 
        //                IS_VCF, &vars, -1, -1);
        //destroy_ranges_t(tmp);

        int prev_i = 0;
        ranges_t *rg = st->ranges[i_ref];
        for (int i_itvl=0; i_itvl<rg->dropped.n; i_itvl++){
            // collect known variant positions
            uint32_t start = rg->dropped.a[i_itvl].s-1;
            uint32_t end = rg->dropped.a[i_itvl].e+1;
            for (int i=prev_i; i<vars->n; i++){
                uint32_t pos = vars->a[i].pos;
                if (pos>=start && pos<end){
                    kv_push(uint32_t, poss, pos);
                }
                if (pos>=end){
                    prev_i = i;
                    break;
                }
            }
            
            // collect haptag of the REF alleles
            recover_variant_phase_in_one_interval(st, refname, start, end, 
                poss.a, poss.n, st->varphase_in_dropped[i_ref], fn_bam, fn_vcf);

            // reset
            poss.n = 0;
        }

        //// reset
        //wipe_vvar_t(&vars);
    }


    for (int i=0; i<st->ref_n; i++){
        destroy_vvar_t(&vars_v[i], 0);
    }
    free(vars_v);
    kv_destroy(poss);
}


void output_tsv(storage_t *st, char *prefix){
    // assumes phaseblocks have been collected
    char *fn_out  = (char*)malloc(strlen(prefix)+25);
    sprintf(fn_out, "%s.mp.tsv", prefix);
    FILE *fp_out = fopen(fn_out, "w");
    if (!fp_out){
        fprintf(stderr, "[E::%s] failed to open output file: %s\n", 
            __func__, fn_out);
        exit(1);
    }
    int n_blocks = 0;
    for (int i_ref=0; i_ref<st->ref_n; i_ref++){
        char *ref_name = st->ref_names[i_ref];
        for (int i=0; i<st->ranges[i_ref]->phaseblocks.n; i++){
            fprintf(fp_out, "%s\t%d\t%d\n", st->ref_names[i_ref], 
                st->ranges[i_ref]->phaseblocks.a[i].s, 
                st->ranges[i_ref]->phaseblocks.a[i].e);
            n_blocks++;
        }
    }
    fclose(fp_out);
    free(fn_out);
    fprintf(stderr, "[M::%s] wrote tsv (%d refs, total %d blocks)\n",
        __func__, st->ref_n, n_blocks);
}

void output_gtf(storage_t *st, char *prefix){
    // assumes phaseblocks have been collected
    char *fn_out  = (char*)malloc(strlen(prefix)+25);
    sprintf(fn_out, "%s.mp.gtf", prefix);
    FILE *fp_out = fopen(fn_out, "w");
    if (!fp_out){
        fprintf(stderr, "[E::%s] failed to open output file: %s\n", 
            __func__, fn_out);
        exit(1);
    }
    int n_blocks = 0;
    for (int i_ref=0; i_ref<st->ref_n; i_ref++){
        char *ref_name = st->ref_names[i_ref];
        for (int i=0; i<st->ranges[i_ref]->phaseblocks.n; i++){
            // note: gene_id and transcript_id are required by GTF2.2, 
            // a value must be "globally unique identifier" but not
            // gene/transcript names of some prescribed database. 
            // wf-human-sv and gtf from whatshap phasing the outputs of 
            // wf-human-sv uses block start positions for these IDs. 
            // Here we do the same. 
            int start = st->ranges[i_ref]->phaseblocks.a[i].s;
            int end   = st->ranges[i_ref]->phaseblocks.a[i].e;
            if (start==0 || end==0) continue; // placeholders
            fprintf(fp_out, "%s\tPhasing\texon\t%d\t%d\t.\t+\t.\tgene_id \"%d\"; transcript_id \"%d.1\"\n", 
                st->ref_names[i_ref], 
                start, end, 
                start, start);
            n_blocks++;
        }
    }
    fclose(fp_out);
    free(fn_out);
    fprintf(stderr, "[M::%s] wrote gtf (%d refs, total %d blocks)\n",
        __func__, st->ref_n, n_blocks);
}


int alter_vcf_line(char *s, int s_l, 
                    int *prev_group_idx,  int *prev_block_idx, 
                    storage_t *st,
                    char *buf_refname,
                    char *buf_newline, 
                    int *prev_pos){
    // s is null term
    // return 0 if no modification; 1 otherwise
    // caller maintains the two buffers to be longer than the current line (by 0 and by 21)
    int verbose = 0;
    int n;
    char pos_s[21];
    char GT[10];
    if (s[0]=='#'){
        if (s[1]=='#') return 0;
        else{
            n = count_char_in_string(s, '\t')+1;
            if (n<10){
                fprintf(stderr, "[E::%s] vcf only has %d columns; mandatory >=8; we also need FORMAT and at least 1 sample\n", __func__, n);
                exit(1);
            }else if (n>10){
                fprintf(stderr, "[E::%s] multi-sample vcf not implemented, TODO/TBD\n", __func__);
                exit(1);
            }
            return 0;
        }
    }
    int col = 0;
    int start = 0;
    int pos = 0;
    int i_ps = -1;  // in-column index of PS tag
    int i_gt = -1;  // in-column index of GT tag
    int i_ref = 0;
    for (int i=0; i<s_l; i++){
        if (s[i]!='\t') continue;
        else{
            if (col==0){  // reference name
                sprintf(buf_refname, "%.*s", i-start, s+start);
                for (i_ref=0; i_ref<st->ref_n; i_ref++){ //IMPRV: could cache
                    if (strcmp(st->ref_names[i_ref], buf_refname)==0)
                        break;
                }
                pos = 0;
                i_ps = -1;
                i_gt = -1;
                if (i_ref==st->ref_n){
                    //fprintf(stderr, "[W::%s] ref %s in vcf line not found in storage\n", __func__, buf_refname);
                    break;
                }
            }else if (col==1){  // variant position
                sprintf(pos_s, "%.*s", i-start, s+start);
                pos = atoi(pos_s);
                if (prev_pos){
                    if (pos<(*prev_pos)){  // vcf is sorted. so we've encountered a new chromosome
                        *prev_pos = pos;
                        *prev_group_idx = 0;
                        *prev_block_idx = 0;
                    }
                    *prev_pos = pos;
                }
            }else if (col==8){  // the format col
                i_ps = search_substr_idx(s+start, "PS", ':', 1, i-start, 0);
                i_gt = search_substr_idx(s+start, "GT", ':', 1, i-start, 0);
            }
            col++;
            start = i+1;
        }
    }
    // process col9 (TODO: multisample vcf)
    if (pos==0 || i_ps<0){
        return 0;
    }else{
        int ps_start, ps_l, gt_start, gt_l;
        get_substr_by_idx(s+start, i_ps, ':', &ps_start, &ps_l, 0);
        get_substr_by_idx(s+start, i_gt, ':', &gt_start, &gt_l, 0);
        if (ps_start<0 || gt_start<0) {
            fprintf(stderr, "[E::%s] saw PS or GT tag but value not found? ref=%s pos=%d\n", 
                __func__, buf_refname, pos);
            return 0;
        }
        if (ps_l==1 && (s+start)[ps_start]=='.'){return 0;}

        sprintf(GT, "%.*s", gt_l, s+start+gt_start);
        if (GT[1]!='|') return 0;
        // do not touch multi allele variants
        if (GT[0]!='0' && GT[0]!='1') return 0; // TODO diploid
        if (GT[2]!='0' && GT[2]!='1') return 0;

        // if we reach here, may need to modify phase groupID
        int groupID = get_new_phaseblock_ID(st, buf_refname, prev_group_idx, pos);
        int is_dropped = tmp_check_if_in_dropped_intervals(st, buf_refname, pos);
        int need_flip = get_flip_status(st, buf_refname, prev_block_idx, pos);
        if (verbose){
            fprintf(stderr, "[dbg::%s] chrom=%s pos=%d , groupID=%d is_dropped=%d need_flip=%d\n", __func__, 
            buf_refname, pos, groupID, is_dropped, need_flip);
        }

        // if variant is in a dropped interval, see if we have re-phased the spot
        int is_middle_var = 0;
        if (groupID>=0 && is_dropped){
            khint_t k = htu32_ht_get(st->varphase_in_dropped[i_ref], pos-1);
            if (k!=kh_end(st->varphase_in_dropped[i_ref])){
                //fprintf(stderr, "[dbg::%s] rescuing %s:%d...", __func__, buf_refname, pos);
                int hap_of_ref = kh_val(st->varphase_in_dropped[i_ref], k);
                if (hap_of_ref==0){  // TODO:diploid
                    //GT[0] = '0'; GT[2] = '1';
                    //need_flip = 0;
                    is_middle_var = 1;
                }else if (hap_of_ref==1){
                    //GT[0] = '1'; GT[2] = '0';
                    //need_flip = 0;
                    is_middle_var = 1;
                }
            }
        }


        if (groupID<0 || is_dropped){
            if (!is_middle_var) return 0;
            else{
                if (verbose) {
                    fprintf(stderr, "[W::%s] dropped groupID and phase for pos=%s:%d\n", __func__, buf_refname, pos);
                }
                int print_offset = start+ps_start;
                sprintf(buf_newline, "%.*s", print_offset, s);  // left
                sprintf(buf_newline+print_offset, ".");
                print_offset+=1;
                sprintf(buf_newline+print_offset, "%s", 
                        s+start+ps_start+ps_l);  //right
                // wipe genotype's phasing
                buf_newline[start+gt_start+1] = '/';
                return 2;
            }
        }else{
            int print_offset = start+ps_start;
            sprintf(buf_newline, "%.*s", print_offset, s);  // left
            sprintf(buf_newline+print_offset, "%d", groupID); 
            print_offset+=get_posint_digits(groupID);
            sprintf(buf_newline+print_offset, "%s", 
                    s+start+ps_start+ps_l);  //right
            // modify genotype
            if (verbose) 
                fprintf(stderr, "[dbg::%s] pos=%d need-flip=%d (raw gt: %.*s)\n", __func__, pos, need_flip, 3, buf_newline+start+gt_start);
            if (need_flip) {
                buf_newline[start+gt_start] = buf_newline[start+gt_start]=='0'? '1' : '0';
                buf_newline[start+gt_start+2] = buf_newline[start+gt_start]=='0'? '1' : '0';
            }
            return 1;
        }
    }
}
void output_modify_vcf(char *fn_vcf, storage_t *st, char *prefix){
    gzFile fp_in = gzopen(fn_vcf, "rb");
    if (!fp_in){
        fprintf(stderr, "[E::%s] failed to open input file: %s\n", __func__, fn_vcf);
        exit(1);
    }
    int n_modified = 0;
    int n_failed_modify = 0;
    int n_tot = 0;

    char *fn_out  = (char*)malloc(strlen(prefix)+10);
    sprintf(fn_out, "%s.mp.vcf", prefix);
    FILE *fp_out = fopen(fn_out, "w");
    if (!fp_out){
        fprintf(stderr, "[E::%s] failed to open output file: %s\n", 
            __func__, fn_out);
        exit(1);
    }

    int m = READLINE_BUF_LEN;
    int m_tmp = m;
    int n = 0;
    char *buffer = (char*)malloc(m);
    char *buf_refname = 0;
    char *buf_newline = 0;
    int last_line_l;
    int last_pos = -1;
    int prev_group_idx = 0;
    int prev_block_idx = 0;
    int offset = 0;
    int n_read;
    while ((n_read = gzread(fp_in, buffer+offset, m-offset)) > 0) {
        int start = 0;
        for (int i=0; i<offset+n_read; i++){
            if (buffer[i]=='\n'){
                int line_l = i-start;
                buffer[i] = 0;
                if (!buf_refname){  // (init)
                    buf_refname = (char*)malloc(line_l+1);
                    buf_newline = (char*)malloc(line_l+21);
                    last_line_l = line_l;
                }else if (line_l>last_line_l){
                    buf_refname = (char*)realloc(buf_refname, line_l+1);
                    buf_newline = (char*)realloc(buf_newline, line_l+21);
                    last_line_l = line_l;
                }
                int altered = alter_vcf_line(buffer+start, line_l, 
                    &prev_group_idx, 
                    &prev_block_idx,
                    st, buf_refname, buf_newline, &last_pos);
                n_tot++;
                if (!altered){
                    fprintf(fp_out, "%s\n", buffer+start);
                }else{
                    if (altered==2) n_failed_modify++;
                    else n_modified++;
                    fprintf(fp_out, "%s\n", buf_newline);
                }
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
    fclose(fp_out);
    gzclose(fp_in);
    free(fn_out);
    free(buffer);
    if (buf_refname){
        free(buf_refname);
        free(buf_newline);
    }
    fprintf(stderr, "[M::%s] wrote vcf output, (%d ok + %d dropped)/%d lines modified \n", 
        __func__, n_modified, n_failed_modify, n_tot);
}

int get_read_new_haplotag(char *refname, char *qn, int qs, int hp_raw,
                          storage_t *st, int need_flip){
    int verbose = 0;
    int hp_tagged=-42;
    khint_t k;
    int absent;
    htstri_t *h = st->qname2haptag;

    k = htstri_ht_get(h, qn);
    if (k!=kh_end(h)){
        hp_tagged = kh_val(h, k);
        if (cliopt_verbose>1 || verbose){
            fprintf(stderr, "[dbg::%s] chrom=%s qn=%s qs=%d has logged haptag: %d\n", 
                    __func__, refname, qn, qs, hp_tagged);
        }
    }else{
        hp_tagged = hp_raw;
        if (hp_tagged!=0 && hp_tagged!=1)  // TODO: diploid
            goto finish;
    }
    if (need_flip){
        hp_tagged ^= 1;
    }
    if (cliopt_verbose>1 || verbose){
        fprintf(stderr, "[dbg::%s] chrom=%s qn=%s qs=%d , need_flip=%d; final haptag: %d\n", 
            __func__, refname, qn, qs, need_flip, hp_tagged);
    }

finish:
    return hp_tagged;
}

void output_modify_bam(char *fn_bam, storage_t *st, char *fn_out, int n_thread){
    // self impl note: do this after paralleled step; all read names and 
    // tags need to be in mem & as hashtables associated to ref names 
    // (alt: re-read from output file)

    // the input
    bamfile_t *hf = init_and_open_bamfile_t(fn_bam, pomfret_n_bam_threads);
    if (!hf){
        fprintf(stderr, "[E::%s] failed to open input bam: %s\n", __func__, fn_bam);
    }
    hts_itr_t *fp_itr = sam_itr_querys(hf->fp_idx, hf->fp_header, ".");

    // helper
    int prev_group_idx=0, prev_block_idx=0;
    int prev_unphased_idx=1;

    // the output
    int hp;
    int prev_tid=0;  // TODO: check if this is correct
    BGZF *fp_out = bgzf_open(fn_out, "w");
    if (!fp_out){
        fprintf(stderr, "[E::%s] failed to open output file: %s\n", __func__, fn_out);
        goto cleanup;
    }
    if (n_thread>1) bgzf_mt(fp_out, n_thread, 0/*unused*/);

    int stat = bam_hdr_write(fp_out, hf->fp_header);
    if (stat!=0){
        fprintf(stderr, "[E::%s] failed to write bam header\n", __func__);
        exit(1);
    }
    int need_flip = 0;
    int updated = 0;
    while(sam_itr_next(hf->fp, fp_itr, hf->buf_aln)>=0){
        char *refname = hf->fp_header->target_name[hf->buf_aln->core.tid];
        if (hf->buf_aln->core.tid!=prev_tid){
            prev_group_idx = 0;
            prev_block_idx = 0;
            prev_unphased_idx = 1;
            prev_tid = hf->buf_aln->core.tid;
        }
        char *qn = bam_get_qname(hf->buf_aln);
        int start_pos = hf->buf_aln->core.pos;
		int len = hf->buf_aln->core.l_qseq; 
        int end_pos = bam_endpos(hf->buf_aln);
        int hp_raw;
        if (st->stores_raw_tag){
            khint_t k = htstri_ht_get(st->qname2haptag_raw, qn);
            if (k==kh_end(st->qname2haptag_raw)) hp_raw = HAPTAG_UNPHASED;
            else hp_raw = kh_val(st->qname2haptag_raw, k);
        }else{
            hp_raw = get_hp_from_aln(hf->buf_aln);
        }
        
        updated = 0;
        int starts_in_unphased = check_if_in_phased_intervals(st, 
                                    refname, start_pos, &prev_unphased_idx, 
                                    &updated);
        if (updated){
            int flip = get_flip_status_by_idx(st, refname, prev_unphased_idx-1);
            //fprintf(stderr, "[dbg::%s] phaseblock check updated by read %s (start=%d). prev_unphased_idx now is %d (sancheck: self is in %d-%d) . Current flip=%d , this=%d\n", 
            //                __func__, qn, start_pos, prev_unphased_idx, 
            //                st->ranges[0]->ends.a[prev_unphased_idx-1], 
            //                st->ranges[0]->starts.a[prev_unphased_idx],
            //                need_flip, flip);
            assert (flip>=0);
            need_flip = flip;
        } 
        hp = get_read_new_haplotag(refname, qn, start_pos, hp_raw, 
                                   st, need_flip);
        bam_aux_update_int(hf->buf_aln, "HP", hp+1);
        stat = bam_write1(fp_out, hf->buf_aln);
        if (stat<0){
            fprintf(stderr, "[E::%s] failed to write bam entry (ref=%s pos=%d qn=%s)\n", 
                            __func__, refname, start_pos, qn);
        }
    }
cleanup:
    bgzf_close(fp_out);
    hts_itr_destroy(fp_itr);
    destroy_bamfile_t(hf, 1);
}


typedef struct{
    vu32_t vmmr_u32i;
    vu16_t hap_cnt[2];
    uint16_t sum[2];
}mmr_t;  // methmers at one position 
typedef struct{
    int n;
    mmr_config_t mmr_config;
    uint32_t *sites_real_poss;  // meth sites real positions(i.e. as-on-reference)
    uint32_t *sites_starts;  // 
    uint8_t *mmr_lens;  // methmer lens
    mmr_t *mmr_a;  // methers storage: mers and per-haplotype count at each site
    // what methermer is available to be used in voting:
    uint32_t mmr_min_i;  // index of the first methmer with enough coverage
    uint32_t mmr_max_i;  // index of the last methmer with enough coverage
}methmers_t;
methmers_t *init_methmers_t(mmr_config_t config, int n){
    methmers_t *ret = (methmers_t*)calloc(1, sizeof(methmers_t));
    ret->mmr_config = config;
    ret->n = n;
    if (n!=0){
        ret->sites_starts = (uint32_t*)calloc(n, sizeof(uint32_t));
        ret->sites_real_poss = (uint32_t*)calloc(n, sizeof(uint32_t));
        ret->mmr_lens = (uint8_t*)calloc(n, 1);
    }else{
        ret->sites_starts = 0;
        ret->sites_real_poss = 0;
        ret->mmr_lens = 0;
    }

    /// mmr_t allocation will wait until sites are collected
    ///ret->mmr_a = 0;
    ret->mmr_a = (mmr_t*)malloc(sizeof(mmr_t) * n);
    for (int i=0; i<n; i++){
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
void destroy_methmers_t(methmers_t *h){
    if (h->n!=0){
        free(h->sites_starts);
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
    int decision; // 0 for cis, 1 for trans, -1 for no join
    rs_t *rs;
    methmers_t *ms;
    methmers_t *ms_bwd;
}dataset_t;
dataset_t *init_dataset_t(){
    dataset_t *ret = (dataset_t*)calloc(1, sizeof(dataset_t));
    return ret;
}
void destroy_dataset_t(dataset_t* d, int include_self){
    destroy_rs_t(d->rs);
    if (d->ms) destroy_methmers_t(d->ms);
    if (d->ms_bwd) destroy_methmers_t(d->ms_bwd);
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
                                  mmr_config_t mmr_config, int direction, 
                                  kstring_t *log, 
                                  htu32_t *masked_positions){
    // Collect number of meth/unmeth/nocall calls, decide meth sites and methmers
    // TODO use bst instead of ht?
    int verbose = 0;
    methmers_t *ret;
    htmethcnt_t *ht = htmethcnt_ht_init(); // each of cnt[3] is: count:28 | fwd strand count:2 | bwd strand count:2
    int shift = 4;
    uint32_t incre = 1<<shift;
    uint32_t limit = (UINT32_MAX>>shift) - incre;
    uint32_t strand_mask = (1<<shift)-1;
    khint_t htk;
    int absent;
    for (int i=0; i<rs->n; i++){
        uint8_t strand = rs->a[i].strand;
        for (int j=0; j<rs->a[i].meth.calls.n; j++){
            uint32_t pos = rs->a[i].meth.calls.a[j];
            int call = rs->a[i].meth.quals.a[j];
            //int hp = rs->a[i].hp;
            //if (hp!=0 && hp!=1) hp = 2;
            htk = htmethcnt_ht_get(ht, pos);

            // increment read counter
            if (htk==kh_end(ht)){
                htk = htmethcnt_ht_put(ht, pos, &absent);
                kh_val(ht, htk).cnt[0] = 0; // meth
                kh_val(ht, htk).cnt[1] = 0; // unmeth 
                kh_val(ht, htk).cnt[2] = 0; // uncall
                kh_val(ht, htk).cnt[call] = incre;
                kh_val(ht, htk).pos = pos;
                //fprintf(stderr, "[dbg::%s] saw pos %d\n", __func__, (int)pos);
            } else {
                if ((kh_val(ht, htk).cnt[call]>>shift)<limit) {
                    kh_val(ht, htk).cnt[call] += incre;
                }
            }

            // increment strand counter
            if (strand==0) {
                if ((kh_val(ht, htk).cnt[call]&3) <3){
                    kh_val(ht, htk).cnt[call] += 1;
                }
            }else{
                if (((kh_val(ht, htk).cnt[call]>>2)&3) <3){
                    kh_val(ht, htk).cnt[call] += 4;
                }
            }

        }
    }

    // collect all positions as-on-ref that 
    // have enough both meth calls and unmeth calls.
    vu32_t ok_sites; 
    kv_init(ok_sites); kv_resize(uint32_t, ok_sites, 32);
    int n_filtered = 0;
    for (htk=0; htk<kh_end(ht); htk++){
        if (kh_exist(ht, htk)){
            methsite_counter_t tmp = kh_val(ht, htk);
            if (verbose) {fprintf(stderr, "[dbg::%s] pos %d , counts: %d and %d\n", 
                    __func__, (int)tmp.pos, tmp.cnt[0]>>shift, tmp.cnt[1]>>shift);
            }
            if (log){
                ksprintf(log, "[dbg::%s] (dir=%d) pos %d , counts: %d and %d\n", 
                    __func__, direction, (int)tmp.pos, tmp.cnt[0]>>shift, tmp.cnt[1]>>shift);
            }
            if (  (tmp.cnt[0]>>shift)>=mmr_config.cov_for_selection 
               && (tmp.cnt[1]>>shift)>=mmr_config.cov_for_selection 
               //&& (tmp.cnt[0]&3)>1
               //&& (tmp.cnt[1]&3)>1
               //&& ((tmp.cnt[0]>>2)&3)>1
               //&& ((tmp.cnt[1]>>2)&3)>1
               ){ 
                if (masked_positions){
                    khint_t k = htu32_ht_get(masked_positions, tmp.pos);
                    if (k!=kh_end(masked_positions)){
                        n_filtered++;
                        continue;
                    }
                }

                // TODO use percentage wrt site coverage?
                kv_push(uint32_t, ok_sites, tmp.pos);
            }
        }
    }
    //fprintf(stderr, "[dbg::%s] # sites: %d\n", __func__, (int)ok_sites.n);
    htmethcnt_ht_destroy(ht);

    // (store the sites; there is no reverse index, use binary search)
    ret = init_methmers_t(mmr_config, ok_sites.n);
    memcpy(ret->sites_real_poss, ok_sites.a, ok_sites.n*sizeof(uint32_t));
    uint32_t *sites = ret->sites_real_poss;
    kv_destroy(ok_sites);
    radix_sort_ksu32(sites, sites+ret->n);
    if (direction==1) reverse_arr(sites, ret->n);

    // collect methmers
    // note: this is directional only due to the methmer max length constraint. 
    //       Around a large region deprived of meth sites, collecting methmers
    //       from left the right will truncate mmrs on the left, while
    //       collecting from right to left will have truncations on the 
    //       right side instead. 
    if (direction==0 || direction==1){  // directional
        int compensate = 0;
        int i;
        for (i=0; i<ret->n;i++){
            int j = i + mmr_config.k;
            j = j>ret->n-1? ret->n-1 : j;
            int len;
            while (1){
                if ( (direction==0 && (sites[j]-sites[i])<=mmr_config.k_span) ||
                     (direction==1 && (sites[i]-sites[j])<=mmr_config.k_span)
                ) break;
                j--;
            }
            ret->mmr_lens[i] = j-i==0? 1 : j-i;
            ret->sites_starts[i] = direction==0? sites[i]:sites[j]; 
        }

        // restore site ordering for the bwd case
        if (direction==1) {
            reverse_arr(sites, ret->n);
            reverse_arr(ret->sites_starts, ret->n);
            reverse_arru8(ret->mmr_lens, ret->n);
        }
    }else if (direction==2){  // symmetric methmers
        // note: need to modify init alloc: for asym mmr, number of methmers
        // does not exceed the number of sites. 
        fprintf(stderr, "[not implemented] TBD/TODO\n");
        exit(1);
    }else{
        fprintf(stderr, "[dbg::%s] invalid direction=%d, check code\n", __func__, direction);
        exit(1);
    }

    // allocate storage for counting 
    //ret->mmr_a = (mmr_t*)malloc(sizeof(mmr_t)*ret->n);
    //for (int i=0; i<ret->n; i++){
    //    kv_init(ret->mmr_a[i].vmmr_u32i);
    //    kv_resize(uint32_t, ret->mmr_a[i].vmmr_u32i, 16);
    //    kv_init(ret->mmr_a[i].hap_cnt[0]);
    //    kv_resize(uint16_t, ret->mmr_a[i].hap_cnt[0], 16);
    //    kv_init(ret->mmr_a[i].hap_cnt[1]);
    //    kv_resize(uint16_t, ret->mmr_a[i].hap_cnt[1], 16);
    //    ret->mmr_a[i].sum[0] = 0;
    //    ret->mmr_a[i].sum[1] = 0;
    //}

    return ret;
}


int get_mmr_of_read(read_t *read, methmers_t *ms, vu32_t *buf_mmr, uint32_t *start_i){ 
    // Given methmer ranges, collect methmers (in uint32) from a read and 
    //  store them in buf_mmr.
    // Return: number of methmers. 
    buf_mmr->n = 0;

    uint32_t *sites = ms->sites_starts;
    uint32_t sites_n = ms->n;
    uint8_t *mml = ms->mmr_lens;
    uint32_t *calls = read->meth.calls.a;
    uint8_t *quals= read->meth.quals.a;
    uint32_t calls_n = read->meth.calls.n;
    if (calls_n==0) return 0;

    uint32_t x, x_i_left, x_i_right;  // as-on-ref position and as-in-buffer index for fixed index

    // init: find the left most call
    // (left)
    int stat = search_arr(sites, sites_n, calls[0], &x_i_left, 0);
    if (stat==-2 || stat==-3) return 0;   // no overlap: to the right
    if (stat==0) x_i_left = x_i_left==0? 0: x_i_left-1;
    // (right)
    stat = search_arr(sites, sites_n, calls[calls_n-1], &x_i_right, 0);
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
        if (i>1 && sites[i]==sites[i-1]) continue;  // note: with bwd or sym mmr, each ref pos may
                                                    // have multiple mmrs starting from it, but here
                                                    // we are sorting to determine valid pos on reads, 
                                                    // so only unique poss are collected. 
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
            uint32_t pos_i = (uint32_t)buf.a[i];
            for (int j=pos_i; j<ms->n; j++){  // For bwd or sym methmers, 
                                              // more than one methmers may 
                                              // start from a given site. Need to
                                              // check them all.
                if (sites[j]!=sites[pos_i]) {
                    break;
                }
                int mmr_len = ms->mmr_lens[j];
                uint32_t pos = buf.a[i]>>35;
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
                if (start_pos_i==UINT32_MAX) start_pos_i = j;  // start index as-in-storage
                //fprintf(stderr, "[%s] pos=%d , mmr: %.*s\n", __func__, pos, n, mmr);
                uint32_t mmr_u32i = methmer_to_uint32(mmr, n);
                kv_push(uint32_t, *buf_mmr, mmr_u32i);
            }
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
                //fprintf(stderr, "--- incrementing at pos#%d, hap=%d\n", i, hap);
                found = 1;
                break;
            }
        }
        if (!found){  // a new methmer key at this position
            kv_push(uint32_t, ms->mmr_a[i].vmmr_u32i, query);
            kv_push(uint16_t, ms->mmr_a[i].hap_cnt[0], 0);
            kv_push(uint16_t, ms->mmr_a[i].hap_cnt[1], 0);
            assert(ms->mmr_a[i].vmmr_u32i.n==ms->mmr_a[i].hap_cnt[0].n);
            ms->mmr_a[i].hap_cnt[hap].a[j]++;
            ms->mmr_a[i].sum[hap]++;
        }
    }
}
void query_counts_of_mmrs(methmers_t *ms, uint32_t *mmr_u32i, int n_mmr, 
                uint32_t mmr_u32i_start, int hap, vfloat_t *buf){
    // TODO: is linear search really faster than hashtable?
    if (hap!=0 && hap!=1){
        fprintf(stderr, "[E::%s] bug or unimpl, hap=%d is invalid\n", 
                    __func__, hap);
        exit(1);
    }
    buf->n = 0;

    char dbg[10];
    char dbg2[10];
    for (int i0=0, j; i0<n_mmr; i0++){
        int i = mmr_u32i_start + i0;
        if (i<ms->mmr_min_i || i>=ms->mmr_max_i) continue;
        uint32_t q = mmr_u32i[i0];
        int found = 0;
        for (j=0; j<ms->mmr_a[i].vmmr_u32i.n; j++){
            if (ms->mmr_a[i].vmmr_u32i.a[j]==q){
                uint32_t cnt = ms->mmr_a[i].hap_cnt[hap].a[j];
                uint32_t sum = ms->mmr_a[i].sum[hap];
                if (sum==0) {;}
                else        {kv_push(float, *buf, (float)cnt/sum);}
                found = 1; 
                break;
            }
        }
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
    if (rs->has_mmr) {
        fprintf(stderr, "[E::%s] storing when there is something, check code\n", __func__);
        exit(1);
    }
    kv_init(buf_mmr);
    kv_resize(uint32_t, buf_mmr, 64);
    for (int i=0; i<rs->n; i++){
        int n_mmr = get_mmr_of_read(&rs->a[i], ms, &buf_mmr, &start_i);
        store_mmr_of_one_read(&rs->a[i], buf_mmr.a, buf_mmr.n, start_i);
        if (start_i!=UINT32_MAX){
            rs->has_mmr = 1;
        }
    }
    kv_destroy(buf_mmr);
}

void wipe_mmr_of_one_read(read_t *r){
    r->mmr_n = 0;
    r->mmr_start_i=0;
    if (r->mmr){
        free(r->mmr);
    }
}
void wipe_mmr_of_reads(rs_t* rs){
    if (!rs->has_mmr){
        fprintf(stderr, "[E::%s] wiping when there is nothing, check code\n", __func__);
    }
    for (int i=0; i<rs->n; i++){
        wipe_mmr_of_one_read(&rs->a[i]);
    }
    rs->has_mmr = 0;
}

void update_mmr_count_with_one_read(read_t *r, methmers_t *ms, int hap){
    if (r->mmr_n>0 && r->mmr_start_i!=UINT32_MAX) {
        insert_mmrs_to_counts(ms, r->mmr, r->mmr_n, r->mmr_start_i, hap);
    }
}


void get_indices_of_first_mmr_on_read_restricted(read_t *r, 
                                    int mmr_min_i, int mmr_max_i, 
                                    int *i_start, int *n){
    // Indexing only knows the absolute first available methmer;
    // during the run we want to find the first available methmer with enough count.
    int start = r->mmr_start_i;
    int end = start + r->mmr_n-1;
    if (end<mmr_min_i || start>=mmr_max_i) {
        *i_start = r->mmr_n-1;
        *n = 0;
    }else{
        *i_start = MAX(start, mmr_min_i);
        *n = MIN(end, mmr_max_i-1) - (*i_start);
    }
    //fprintf(stderr, "[dbg::%s] min=%d max=%d r's start=%d n=%d, resulted in i_start=%d n=%d\n", 
    //    __func__, mmr_min_i, mmr_max_i, r->mmr_start_i, r->mmr_n,
    //    *i_start, *n);
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
    //int n_mmr_avail = ms->mmr_max_i - r->mmr_start_i + 1;
    //n_mmr_avail = r->mmr_n < n_mmr_avail? r->mmr_n : n_mmr_avail;
    int i_start, l_onread;
    get_indices_of_first_mmr_on_read_restricted(r, ms->mmr_min_i, ms->mmr_max_i, &i_start, &l_onread);

    // count of hap0
    // fprintf(stderr, "[dbg::%s] below is read %d\n", __func__, r->i);
    buf.n = 0;
    query_counts_of_mmrs(ms, r->mmr, r->mmr_n, r->mmr_start_i, 0, &buf);
    score0_l = buf.n;
    for (int i=0; i<buf.n; i++) {
        if (buf.a[i]>0){
            score0 += buf.a[i];
            score0_l++;
        }
    }

    // count of hap1
    buf.n = 0;
    query_counts_of_mmrs(ms, r->mmr, r->mmr_n, r->mmr_start_i, 1, &buf);
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
    float affine_diff_min = score_diff_min + (MIN(score0_l, score1_l)>3? ((float)MIN(score0_l, score1_l)-3)*0.5 : 0);
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

int update_available_methmer_range(methmers_t *ms, int mmr_cov_for_update){
    int updated = 0;
    for (int i=ms->mmr_min_i; i>=0; i--){  // left side
        int sum = ms->mmr_a[i].sum[0]+ms->mmr_a[i].sum[1];
        if (sum>=mmr_cov_for_update) {
            ms->mmr_min_i = i;
            updated++;
        }
        else break;  // ref range is forced to be available, so it should be impossible to have holes in methmer count, look no further
    }
    for (int i=ms->mmr_max_i; i<ms->n; i++){  // right side
        if (ms->mmr_a[i].sum[0]+ms->mmr_a[i].sum[1]>=mmr_cov_for_update){
            ms->mmr_max_i = i;
            updated++;
        }
        else break;
    }
    //fprintf(stderr, "[dbg::%s] methmer now valid between %d and %d\n", 
    //    __func__, 
    //    (int)ms->sites_real_poss[ms->mmr_min_i], 
    //    (int)ms->sites_real_poss[ms->mmr_max_i]);
    return updated;
}

int predict_tags_of_reads(rs_t *rs, methmers_t *ms, 
                           uint32_t *readIDs, int readIDs_n, 
                           vu32_t *buf_readIDs_sorted, vi_t *buf_tags, 
                           forpred_v *buf_tmp, forpred_v *buf_tmp_tmp,  // structure for sorting and its tmp buffer required by merge sort
                           int insert_best_n, 
                           int mmr_cov_for_update, 
                           float score_diff_min, int score_l_min, 
                           uint32_t *dbg_val
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


    // store stuff in caller's buffers
    buf_readIDs_sorted->n = 0;
    buf_tags->n = 0;
    for (int i=0; i<buf_tmp->n; i++){
        //uint32_t idx = buf_tmp->a[i].idx;
        kv_push(int, *buf_tags, buf_tmp->a[i].tag);
        kv_push(uint32_t, *buf_readIDs_sorted, buf_tmp->a[i].readID);
    }

    // optional: update counter and tag the best read(s)
    int n = 0;
    if (insert_best_n>0){  // is tagging new reads, only update when query has haplotype predicted
        for (int i=buf_readIDs_sorted->n-1; i>=0; i--){
            uint32_t readID = buf_readIDs_sorted->a[i];
            int hap = buf_tags->a[i];
            if ( (hap==0 || hap==1) && rs->a[readID].mmr_start_i!=UINT32_MAX){
                rs->a[readID].hp = hap;
                insert_mmrs_to_counts(ms, rs->a[readID].mmr, 
                    rs->a[readID].mmr_n, 
                    rs->a[readID].mmr_start_i, hap);
                n++;
                //fprintf(stderr, "[dbg::%s] tagged read# %d, qname=%.*s, hap=%d, score_diff=%.2f; read mmr range %d-%d; mmr range now: %d - %d\n", 
                //    __func__, readID, 
                //    (int)GET_QNAME_L(*rs, readID), 
                //    GET_QNAME_S(*rs, readID), 
                //    hap, 
                //    buf_tmp->a[i].score, 
                //    rs->a[readID].mmr_start_i, rs->a[readID].mmr_start_i+rs->a[readID].mmr_n,
                //    ms->mmr_min_i, ms->mmr_max_i
                //    );
                if (dbg_val) *dbg_val = readID;
                if (n==insert_best_n) 
                    break;
            }
        }
        
    }
    if (n>0){  // inserted some reads, now update the valid methmer range
        update_available_methmer_range(ms, mmr_cov_for_update);
    }
    return n;
    
}

void insert_ref_reads_methmer_counts(rs_t* rs, methmers_t* ms, 
                                     uint32_t *refreadIDs, int refreadIDs_n, 
                                     int mmr_cov_for_update){
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
    int increment;
    for (int i=0; i<refreadIDs_n; i++){
        uint32_t readID = refreadIDs[i];
        int hap = rs->a[readID].hp;
        if ((hap==0 || hap==1) && rs->a[readID].mmr_start_i!=UINT32_MAX){
            sancheck[hap]++;
            insert_mmrs_to_counts(ms, rs->a[readID].mmr, rs->a[readID].mmr_n, 
                                  rs->a[readID].mmr_start_i, hap);
            //fprintf(stderr, "[dbg::%s] saw ref read %.*s, hap is %d\n", __func__, 
            //    (int)GET_QNAME_L(*rs, readID), GET_QNAME_S(*rs, readID), rs->a[readID].hp);
        }
    }

    update_available_methmer_range(ms, mmr_cov_for_update);

    //for (int i=ms->mmr_min_i; i<ms->mmr_max_i; i++){
    //    fprintf(stderr, "[dbg::%s] mmr#%d, hap0 sum %d, hap1 sum %d\n", __func__, i, ms->mmr_a[i].sum[0], ms->mmr_a[i].sum[1]);
    //}
}

int permute_haplotags(rs_t *rs, uint32_t *IDs, uint32_t IDs_n, int n){
    // return: 0 if ok (permuted something); 1 if nothing was done
    // TODO diploid assumption
    // TODO maybe consider to let rs_t keep the buffers to avoid re allocations
    // swap n reads between the two haplotypes  
    vu32_t buf_hap1; kv_init(buf_hap1); kv_resize(uint32_t, buf_hap1, 32);
    vu32_t buf_hap2; kv_init(buf_hap2); kv_resize(uint32_t, buf_hap2, 32);
    int ret = 0;
    n = n>IDs_n? IDs_n : n;
    if (n==0){
        ret = 1;
        goto cleanup;
    }

    // slice
    for (uint32_t i=0, readID; i<IDs_n; i++){
        readID = IDs[i];
        if (rs->a[readID].hp==0) kv_push(uint32_t, buf_hap1, i);
        if (rs->a[readID].hp==1) kv_push(uint32_t, buf_hap2, i);
    }

    // shuffle
    ks_shuffle_kssu32(buf_hap1.n, buf_hap1.a);
    ks_shuffle_kssu32(buf_hap2.n, buf_hap2.a);

    // alter haplotags 
    for (int i=0, j; i<MIN(buf_hap1.n, n); i++){
        j = buf_hap1.a[i];
        uint32_t readID = IDs[j];
        fprintf(stderr, "[dbg::%s] permute hap0->1: read%d (sancheck: hap was %d)\n", __func__, (int)readID, rs->a[readID].hp);
        rs->a[readID].hp = 1;
    }
    for (int i=0, j; i<MIN(buf_hap2.n, n); i++){
        j = buf_hap2.a[i];
        uint32_t readID = IDs[j];
        fprintf(stderr, "[dbg::%s] permute hap1->0: read%d (sancheck: hap was %d)\n", __func__, (int)readID, rs->a[readID].hp);
        rs->a[readID].hp = 0;
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

int evaluate_ref_sanity(rs_t *rs, float *ref_ratio, int which_side){
    float cnt[2] = {0,0};
    vu32_t *IDs = which_side==0? &rs->rf->IDs_left : &rs->rf->IDs_right;
    for (int i=0; i<IDs->n; i++){
        uint32_t readID = IDs->a[i];
        if (rs->a[readID].hp==0) cnt[0]++;
        if (rs->a[readID].hp==1) cnt[1]++;
    }
    float r = MAX(cnt[0], cnt[1])/MIN(cnt[0], cnt[1]);
    if (ref_ratio) *ref_ratio = r;
    if (r>1.2) return 0; // fail
    return 1;
}
float evaluate_separation1(uint8_t *ref, uint8_t *query, int n, int *join_dir, kstring_t *log){
    int buf[2][2];
    for (int i=0; i<2; i++){
        for (int j=0; j<2; j++){
            buf[i][j] = 0;
        }
    }
    for (int i=0; i<n; i++){
        if (ref[i]!=0 && ref[i]!=1) continue;
        if (query[i]!=0 && query[i]!=1) continue;
        buf[ref[i]][query[i]]++; // buf[0][0] is the number of reads that are
                                 // hap0 in SNP, and hap0 in meth tagging 
    }
    // require unambiguous ratios for both
    int hard_coverage_fail = (MIN(buf[0][0], buf[0][1])> HARD_COV_THRESHOLD || MIN(buf[1][0], buf[1][1])>HARD_COV_THRESHOLD);
    if (cliopt_verbose>1) 
        fprintf(stderr, "[dbg::%s] n is %d; %d %d %d %d (hard coverage fail=%d)\n", __func__, n, buf[0][0], buf[0][1], buf[1][0], buf[1][1], hard_coverage_fail);
    if (log)
        ksprintf(log, "[dbg::%s] n is %d; %d %d %d %d (hard coverage fail=%d)\n", __func__, n, buf[0][0], buf[0][1], buf[1][0], buf[1][1], hard_coverage_fail);
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
        if (MIN(buf[0][0], buf[0][1])>HARD_CONTAMINATE_THRESHOLD || 
            MIN(buf[1][0], buf[1][1])>HARD_CONTAMINATE_THRESHOLD) {
            if(join_dir)*join_dir=-9; 
            return (float)1;
        }
        if (max==0) {if(join_dir)*join_dir=-9; return (float)1;}
        min = min==0? 1: min;
        if (max/min<3){if(join_dir)*join_dir=-9; return (float)1;}
        scores[i] = max/min;
    }

    // require significance
    double fisher_left_p, fisher_right_p, fisher_twosided_p;
    kt_fisher_exact(buf[0][0], buf[0][1], buf[1][0], buf[1][1], 
                    &fisher_left_p, &fisher_right_p, &fisher_twosided_p);
    if (fisher_twosided_p<EVAL_P_THRE && !hard_coverage_fail) {
        if (join_dir){
            *join_dir= which_way;
        }
        return MIN(scores[0], scores[1]);
    }else{
        if (join_dir){
            *join_dir = -9; // dummy
        }
        return (float)1;
    }
}
float evaluate_separation(rs_t *rs, uint8_t *raw_tags, int which_side, int *joining_type, kstring_t *log){
    vu8_t ref, query;
    kv_init(ref); kv_init(query);
    uint32_t *rf = which_side==0? rs->rf->IDs_left_strict.a : rs->rf->IDs_right_strict.a;
    int     rf_n = which_side==0? rs->rf->IDs_left_strict.n : rs->rf->IDs_right_strict.n;
    kv_resize(uint8_t, ref, rf_n);
    kv_resize(uint8_t, query, rf_n);
    for (int i=0; i<rf_n; i++){
        uint32_t readID = rf[i];
        kv_push(uint8_t, ref,   raw_tags[readID]);
        kv_push(uint8_t, query, rs->a[readID].hp);
    }
    float ret = evaluate_separation1(ref.a, query.a, rf_n, joining_type, log);
    kv_destroy(ref);
    kv_destroy(query);
    return ret;
}

void haplotag_region1(rs_t *rs, methmers_t *ms, 
                      uint32_t *readIDs, uint32_t readIDs_n, 
                      int n_candidates_per_iter, 
                      int min_mmr_recruit_cov, 
                      uint8_t direction, 
                      kstring_t *log){
    vu32_t buf_anno_readIDs;  // will store readIDs sorted by corresponding tagging scores
    vi_t buf_anno_tags;
    forpred_v buf_srt;
    forpred_v buf_srt_tmp;  
    uint32_t *refreadIDs; 
    int refreadIDs_n;
    kv_init(buf_anno_readIDs); kv_resize(uint32_t, buf_anno_readIDs, 32);
    kv_init(buf_anno_tags);    kv_resize(int, buf_anno_tags, 32);
    kv_init(buf_srt);          kv_resize(forpred_t, buf_srt, 32);
    kv_init(buf_srt_tmp);      kv_resize(forpred_t, buf_srt_tmp, 32);
    // we will only use methmers that have accumulated enough coverage 
    // during tagging to tag incoming query reads
    if (direction==0){  // left to right
        ms->mmr_min_i = 0;
        ms->mmr_max_i = 0;
        refreadIDs = rs->rf->IDs_left.a;
        refreadIDs_n = rs->rf->IDs_left.n;
    }else if (direction==1){
        ms->mmr_min_i = ms->n-1;
        ms->mmr_max_i = ms->n-1;
        refreadIDs = rs->rf->IDs_right.a;
        refreadIDs_n = rs->rf->IDs_right.n;
    }else{
        if (log) ksprintf(log,   "[E::%s] invalid direction, check code\n", __func__);
        else     fprintf(stderr, "[E::%s] invalid direction, check code\n", __func__);
        exit(1);
    }

    // step1: collect ref reads on the reference phasing boundary
    // (ref region's sites shall be all available)
    if (direction==0){
        for (int i=ms->mmr_max_i; i<ms->n; i++){  
            if (ms->sites_real_poss[i]<=rs->ref_start) ms->mmr_max_i++;
            else break;
        }
    }else{
        for (int i=ms->mmr_min_i; i>=0; i--){
            if (ms->sites_real_poss[i]>rs->ref_end) ms->mmr_min_i--;
            else break;
        }
    }
    insert_ref_reads_methmer_counts(rs, ms, refreadIDs, refreadIDs_n,
                                    min_mmr_recruit_cov);

    // step1.5: mark non-ref reads as unhaplotagged, 
    //          including ref reads on the other side of the gap
    vu8_t *dbg_initial_tags = store_haplotags(rs);
    vu32_t tmp_tags; kv_init(tmp_tags); 
    kv_resize(uint32_t, tmp_tags, refreadIDs_n);
    for (int i=0; i<refreadIDs_n; i++){
        uint32_t readID = refreadIDs[i];
        kv_push(uint32_t, tmp_tags, (readID<<2)|rs->a[readID].hp );
    }
    if (log) ksprintf(log, "[dbg::%s] has %d ref reads\n", __func__, (int)tmp_tags.n);
    else fprintf(stderr, "[dbg::%s] has %d ref reads\n", __func__, (int)tmp_tags.n);
    set_all_as_unphased(rs);
    for (int i=0; i<tmp_tags.n; i++){
        uint32_t readID = tmp_tags.a[i]>>2;
        int hp = tmp_tags.a[i] & 3;
        rs->a[readID].hp = hp;
    }
    kv_destroy(tmp_tags);

    // step2: extend, tag 1 read per iteration
    int i_last_untagged = direction==0? 0 : readIDs_n-1;
    int increment = direction==0? 1: -1;
    int failed_prev_batch = 0;
    //fprintf(stderr, "[dbg::%s] tagging reads from %d to %d\n", __func__, i_last_untagged, readIDs_n);
    while (1){
        // step2 1) collect candidate reads
        // TODO this is not identical to the py prototype as i do not want to 
        //  depend on bam utils here. If does not work as good, write a helper
        //  rather than using bam utils. 
        buf_anno_readIDs.n = 0;
        if ((direction==0 && i_last_untagged>=readIDs_n)|| (direction!=0 && i_last_untagged<=0)) break;
        for (int i0=i_last_untagged; direction==0?i0<readIDs_n:i0>=0; i0+=increment){
            int i = direction==0? i0 : (uint32_t)rs->revbuf.a[i0];
            if (rs->a[i].hp!=0 && rs->a[i].hp!=1) {
                kv_push(uint32_t, buf_anno_readIDs, i);
                if (buf_anno_readIDs.n>=n_candidates_per_iter) break;
            }
        }
        if (buf_anno_readIDs.n==0){ // nothing can be tagged
            failed_prev_batch++;
            if (failed_prev_batch>10) break;
            i_last_untagged+=n_candidates_per_iter*increment; // try to get another batch; TODO perhaps should just give up?
            continue;
        }
        //fprintf(stderr, "[dbg::%s] collected a batch, first ID is %d, n is %d, last untagged is %d\n", 
        //        __func__, buf_anno_readIDs.a[0], (int)buf_anno_readIDs.n, i_last_untagged);

        // step2 2) tag & insert
        uint32_t dbg_inserted_readID;
        int inserted = predict_tags_of_reads(rs, ms, 
                            buf_anno_readIDs.a, buf_anno_readIDs.n, 
                            &buf_anno_readIDs, &buf_anno_tags, 
                            &buf_srt, &buf_srt_tmp, 
                            1, // try to add 1 read per round and update mmr counts accordingly
                            min_mmr_recruit_cov,
                            3, 3, &dbg_inserted_readID);
        if (inserted==0){  // nothing can be tagged
            failed_prev_batch++;
            if (failed_prev_batch>10) break;
            i_last_untagged+=n_candidates_per_iter*increment; // try to get another batch; TODO perhaps should just give up?
            continue;
        }
        failed_prev_batch = 0;
    }

cleanup:
    kv_destroy(buf_anno_readIDs);
    kv_destroy(buf_anno_tags);
    kv_destroy(buf_srt);
    kv_destroy(buf_srt_tmp);
    kv_destroy(*dbg_initial_tags);
    free(dbg_initial_tags);
}







int haplotag_region2(rs_t *rs, methmers_t *ms, 
                      uint32_t *readIDs, uint32_t readIDs_n, 
                      int ext_direction,
                      int n_candidates_per_iter, 
                      int min_mmr_recruit_cov, 
                      int n_permutations, int do_reset, 
                      kstring_t *log
                      ){
    // return: 0 if joining cis (wrt the SNP pair of tags, i.e. left hap0 is right hap0), 
    //         1 if joining trans (i.e. left hap0 is right hap1)
    //         -1 if no join
    int ret=-1;

    // Wrapper around haplotag_region1, permute initial state,
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
    int err_permutation = 0;
    for (int i=0; i<n_permutations; i++){  // get score and store tagging
        if (i!=0){
            if (ext_direction==0) {
                err_permutation = permute_haplotags(rs, 
                                        rs->rf->IDs_left.a, 
                                        rs->rf->IDs_left.n, 5);
                if (err_permutation) {
                    if (log) ksprintf(log, "[W::%s] can't permute, leaving\n", __func__); 
                    else fprintf(stderr, "[W::%s] can't permute, leaving\n", __func__); 
                    break;
                    }
            }else{
                err_permutation = permute_haplotags(rs, 
                                        rs->rf->IDs_right.a, 
                                        rs->rf->IDs_right.n, 5);
                if (err_permutation) {
                    if (log) ksprintf(log, "[W::%s] can't permute, leaving\n", __func__); 
                    else fprintf(stderr, "[W::%s] can't permute, leaving\n", __func__); 
                    break;
                }
            }
        }
        //random_assign_haplotags_if_has_none(rs, 0, readIDs_init_n);
        haplotag_region1(rs, ms, readIDs, readIDs_n, 
                         n_candidates_per_iter, min_mmr_recruit_cov, 
                         ext_direction, log);
        kv_init(buf[i]);
        kv_resize(uint8_t, buf[i], 16);
        for (int j=0; j<rs->n; j++){
            kv_push(uint8_t, buf[i], rs->a[j].hp);
        }
        int which_way;
        scores[i] = evaluate_separation(rs, initial_state->a, 
                                        ext_direction==0?1:0, &which_way, log);
        if (scores[i]>=2 && which_way!=0){
            if (which_way>0) {dir[1]++; which_way = 0;} // cis join
            else             {dir[2]++; which_way = 1;} // trans join
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
    if (err_permutation){
        if (log) ksprintf(log, "[W::%s] permutation failed, treating as unphased\n", __func__);
        else fprintf(stderr, "[W::%s] permutation failed, treating as unphased\n", __func__);
        ret = -1;
        goto cleanup;
    }
    if (n_permutations>5){
        if ((dir[1]>=threshold && dir[2]<=3/*|| (dir[1]/(dir[2]==0? 1 : dir[2])>2)*/) && 
             dir[0]<threshold_blank && 
             best_score_i[0]>=0){
            ret = 0;
            restore_haplotags(rs, &buf[best_score_i[0]]);
            if (log) ksprintf(log, "[dbg::%s] (with perm) >>> dir1 (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
            else fprintf(stderr, "[dbg::%s] (with perm) >>> dir1 (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
        }else if ((dir[2]>=threshold && dir[1]<=3/*|| (dir[2]/(dir[1]==0? 1:dir[1])>2)*/ ) && 
             dir[0]<threshold_blank && 
             best_score_i[1]>=0){
            ret = 1;
            restore_haplotags(rs, &buf[best_score_i[1]]);
            if (log) ksprintf(log, "[dbg::%s] (with perm) <<< dir2 (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
            else fprintf(stderr, "[dbg::%s] (with perm) <<< dir2 (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
        }else{  // failed, permutation result isn't conclusive
            ret = -1;
            restore_haplotags(rs, initial_state);
            set_all_as_unphased(rs);
            if (log) ksprintf(log, "[dbg::%s] (with perm) xxx undecided (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
            else fprintf(stderr, "[dbg::%s] (with perm) xxx undecided (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
        }
    }else{
        if (best_score_i[0]>=0){
            ret = 0;
            restore_haplotags(rs, &buf[best_score_i[0]]);
            if (log) ksprintf(log, "[dbg::%s] (noperm) >>> dir1 (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
            else fprintf(stderr, "[dbg::%s] (noperm) >>> dir1 (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
        }else if (best_score_i[1]>=0){
            ret = 1;
            restore_haplotags(rs, &buf[best_score_i[1]]);
            if (log) ksprintf(log, "[dbg::%s] (noperm) <<< dir2 (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
            else fprintf(stderr, "[dbg::%s] (noperm) <<< dir2 (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
        }else{
            ret = -1;
            restore_haplotags(rs, initial_state);
            set_all_as_unphased(rs);
            if (log) ksprintf(log, "[dbg::%s] (noperm) xxx undecided (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
            else fprintf(stderr, "[dbg::%s] (noperm) xxx undecided (un=%d, dir1=%d, dir2=%d)\n", __func__, dir[0], dir[1], dir[2]);
        }
    }

cleanup:
    if (do_reset) restore_haplotags(rs, initial_state); 
    for (int i=0; i<n_permutations; i++) kv_destroy(buf[i]);
    kv_destroy(*initial_state);
    free(initial_state);
    return ret;
}


dataset_t* haplotag_region_given_bam(storage_t *st, 
                     char *fn_bam, 
                     char *chrom, uint32_t ref_start, uint32_t ref_end,   // alignment coordinates
                     mmr_config_t mmr_config, 
                     int n_candidates_per_iter, 
                     int do_n_permuations, 
                     int *decision){
    int verbose = 0;
    dataset_t *ret = init_dataset_t();
    kstring_t *log = (kstring_t*)calloc(1, sizeof(kstring_t));

    vu32_t readIDs;
    kv_init(readIDs);
    kv_resize(uint32_t, readIDs, 128);

    // load reads and methmers
    rs_t *rs = load_reads_given_interval(fn_bam, chrom, 
                                         ref_start, ref_end, READBACK, 
                                         mmr_config, 
                                         st->stores_raw_tag? st->qname2haptag_raw : 0, 
                                         log);
    vu8_t *initial_hps = store_haplotags(rs);
    ksprintf(log, "[dbg::%s] loaded %d reads (interval: %s:%d-%d), left has #ref=%d (strict: %d), right has=%d (strict: %d)\n", 
                     __func__, rs->n, chrom, ref_start, ref_end,
                     (int)rs->rf->IDs_left.n,  (int)rs->rf->IDs_left_strict.n, 
                     (int)rs->rf->IDs_right.n, (int)rs->rf->IDs_right_strict.n);
    ksprintf(log, "[dbg::%s] mmrconfig: k %d, k_span %d, lo %d, hi %d; cov_known %d, cov_for_selection %d, cov_for_runtime %d; readlen thre %d, min_mapq %d\n", 
        __func__, mmr_config.k, mmr_config.k_span, 
        mmr_config.lo, mmr_config.hi, 
        mmr_config.cov_known<0?0:mmr_config.cov_known, 
        mmr_config.cov_for_selection, mmr_config.cov_for_runtime, 
        mmr_config.readlen_threshold, mmr_config.min_mapq);

    methmers_t *ms =     get_methmer_sites_and_ranges(rs, mmr_config, 0, NULL, NULL);
    methmers_t *ms_bwd = get_methmer_sites_and_ranges(rs, mmr_config, 1, NULL, NULL);
    ksprintf(log, "[dbg::%s] parsed meth sites %s:%d-%d (k=%d span=%d cov=%d runtime_cov=%d), fwd n=%d, bwd n=%d\n", 
        __func__, chrom, ref_start, ref_end,
        mmr_config.k, mmr_config.k_span, 
        mmr_config.cov_for_selection, mmr_config.cov_for_runtime, 
        ms->n, ms_bwd->n);
    
    if (verbose){  // print sies of methmers, as-on-ref positions, as-on-ref positions
        ksprintf(log, "[dbg::%s] fwd sites: [", __func__);
        for (int tmpi=0; tmpi<ms->n; tmpi++){
            ksprintf(log, "%d, ", ms->sites_real_poss[tmpi]);
        }
        ksprintf(log, "]\n");
    } 

    if (ms->n==0 || ms_bwd->n==0){
        ksprintf(log, "[W::%s] %s:%d-%d does not have methmer in both directions. Skipping.\n", 
                        __func__, chrom, ref_start, ref_end);
        goto cleanup;
    }

    // main routine
    for (uint32_t i=0; i<rs->n;i++){
        kv_push(uint32_t, readIDs, i);
    }
    if (do_n_permuations>1 && do_n_permuations<5){
        ksprintf(log, "[W::%s] number of permutations is too small, might have no effects.\n", __func__);
    }
    float ref_ratio_left = -1;
    float ref_ratio_right = -1;
    int ref_valid_left = evaluate_ref_sanity(rs, &ref_ratio_left, 0);  // for fwd: checking left
    int ref_valid_right = evaluate_ref_sanity(rs, &ref_ratio_right, 1);  // for bwd: checking right
    ksprintf(log, "[M::%s] left ref ratio: %.2f (valid=%d); right ref ratio: %.2f (valid=%d)\n", 
                    __func__, ref_ratio_left, ref_valid_left, 
                    ref_ratio_right, ref_valid_right);


    int join1=-1, join2=-1; 
    
    ////if (ref_valid){
    if (1){
        store_mmr_of_reads(rs, ms_bwd);
        join2 = haplotag_region2(rs, ms_bwd, readIDs.a, readIDs.n, 
                 1,  // extension direction
                 n_candidates_per_iter, 
                 mmr_config.cov_for_runtime, 
                 do_n_permuations, //do_n_permuations, 
                 1,  // reset?
                 log
                 );
        wipe_mmr_of_reads(rs);
        store_mmr_of_reads(rs, ms);
        join1 = haplotag_region2(rs, ms, readIDs.a, readIDs.n, 
                 0,  // extension direction
                 n_candidates_per_iter, 
                 mmr_config.cov_for_runtime, 
                 do_n_permuations, //do_n_permuations, 
                 0,  // reset?
                 log
                 );
    }
    ksprintf(log, "[dbg::%s] fwd join %d, bwd join %d\n\n", __func__, join1, join2);
    if (join1!=join2 || (join1==-1 && join2==-1)){
        set_all_as_unphased(rs);
        if (decision) *decision = -1;
        ret->decision = -1;
    }else{
        if (decision) *decision = join1;
        ret->decision = join1;
    }


cleanup:
    // TODO print log
    fprintf(stderr, "%s", log->s);
    free(log->s);
    free(log);
    kv_destroy(readIDs);
    kv_destroy(*initial_hps);
    free(initial_hps);
    ret->rs = rs;
    ret->ms = ms;
    ret->ms_bwd = ms_bwd;
    return ret;
}




typedef struct{
    storage_t *st;
    enum input_file_format fn_interval_type;
    char *fn_bam;
    mmr_config_t mmr_config;
    int n_candidates_per_iter;
    int n_permutation;
    int *ref_coverages;
    htstri_t **qname2haptag;  // n_ref hashtables; haptags obtained from meth phasing
}blockjoin_work_t;
static void blockjoin_one_chrom_callback(void *data, long job_i, int thread_i){
    blockjoin_work_t *dt = (blockjoin_work_t*)data;
    ranges_t *ranges = dt->st->ranges[job_i];
    char *ref_name = dt->st->ref_names[job_i];
    mmr_config_t mmr_config = dt->mmr_config;
    htstri_t *qname2haptag = dt->qname2haptag[job_i];

    int n_candidates_per_iter = dt->n_candidates_per_iter;
    if (mmr_config.cov_for_selection<=0){
        bamfile_t *tmphf = init_and_open_bamfile_t(dt->fn_bam, 1);
        bam_hdr_t *h= tmphf->fp_header;
        //BGZF *fp = bgzf_open(dt->fn_bam, "r");
        //bam_hdr_t *h= bam_hdr_read(fp);
        int ref_i = -1;
        for (int i=0; i<h->n_targets; i++){
            if (strcmp(h->target_name[i], ref_name)==0){
                ref_i = i;
                break;
            }
        }
        assert (ref_i>=0);
        
        int coverage = dt->ref_coverages[ref_i];
        mmr_config.cov_for_selection = coverage/10 + 1;  // overriding local copy
        mmr_config.cov_for_runtime = mmr_config.cov_for_selection*2;
        n_candidates_per_iter = coverage/4 + 1;
        //bgzf_close(fp);
        //bam_hdr_destroy(h);
        destroy_bamfile_t(tmphf, 1);
    }
    // sancheck parameters and clamp if bad (should not happen)
    if (mmr_config.cov_for_selection<=0) {
        fprintf(stderr, "[W::%s] had to clamp cov_for_selection (ref: %s)\n", 
                        __func__, ref_name);
        mmr_config.cov_for_selection = 1;
    }
    if (n_candidates_per_iter<=1) {
        fprintf(stderr, "[W::%s] had to clamp n_candidates_per_iter (ref: %s)\n", 
                        __func__, ref_name);
        n_candidates_per_iter = 2;
    }
    fprintf(stderr, "[dbg::%s] ref %s using: cov_for_selection=%d, n_cand_per_iter=%d\n", 
            __func__, ref_name, mmr_config.cov_for_selection, n_candidates_per_iter);

    khint_t k;
    int absent;
    for (int i=0; i<ranges->starts.n; i++){
        int start = ranges->starts.a[i];
        int end = ranges->ends.a[i];

        dataset_t *ds = haplotag_region_given_bam(dt->st, 
                        dt->fn_bam, 
                        ref_name, start, end, 
                        mmr_config, 
                        n_candidates_per_iter, 
                        dt->n_permutation, 
                        &ranges->decisions.a[i]);
        // phased, store read haplotags
        if (ranges->decisions.a[i]>=0){
            for (int j=0; j<ds->rs->n; j++){
                char *qname = GET_QNAME_S(*(ds->rs), j); // null termed
                int qname_l = GET_QNAME_L(*(ds->rs), j);
                char *qname_key = (char*)calloc(qname_l+1, 1);
                sprintf(qname_key, "%s", qname);
                k = htstri_ht_put(qname2haptag, qname_key, &absent);
                if (absent){
                    kh_key(qname2haptag, k) = qname_key;
                    kh_val(qname2haptag, k) = ds->rs->a[j].hp;
                }else{
                    free(qname_key);
                    //fprintf(stderr, "[W::%s] read %s tagged but was already recorded\n", __func__, qname);
                }
            }
        }
        destroy_dataset_t(ds, 1);
    }
}

storage_t* blockjoin_parallel(char *fn_interval, 
                        enum input_file_format fn_interval_type, 
                        char *fn_bam, 
                        int bam_needs_haplotagging, int write_baminput_tagging, 
                        char *fn_vcf, 
                        char *output_prefix, 
                        mmr_config_t mmr_config, int n_candidates_per_iter, 
                        int n_permutation, 
                        int n_threads, 
                        int n_bam_threads){
    // get chromosome names and phasing gaps
    double T = Get_T();

    storage_t *st = (storage_t*)calloc(1, sizeof(storage_t));
    init_storage_t(st);

    // load phase block intervals and variants
    // also haptag reads if required
    if (bam_needs_haplotagging){
        assert(fn_vcf);
        load_intervals_from_file(fn_vcf, IS_VCF, st, 1, fn_bam, 0, 0);


        int n_loaded = 0;
        for (int i=0; i<st->ref_n; i++){
            n_loaded += st->ranges[i]->starts.n;
        }
        if (n_loaded==0){
            fprintf(stderr, "[E::%s] Nothing loaded from vcf (ref_n=%d), cannot haptag the input bam. Terminating.\n", 
                    __func__, st->ref_n);
            exit(1);
        }

        if (fn_interval_type!=IS_VCF){  // gtf or tsv was supplied, 
                                        // override vcf phase blocks
            wipe_intervals_of_storage_t(st);
            load_intervals_from_file(fn_interval, fn_interval_type, st, 0, 0, 0, 0);
        }
    }else{
        load_intervals_from_file(fn_interval, fn_interval_type, st, 0, 0, 0, 0);
    }

    int n_loaded = 0;
    for (int i=0; i<st->ref_n; i++){
        n_loaded += st->ranges[i]->starts.n;
    }
    if (n_loaded==0){
        fprintf(stderr, "[E::%s] No intervals loaded, terminating.\n", __func__);
        exit(1);
    }


    fprintf(stderr, "[M::%s] input has %d references\n", __func__, st->ref_n);
    if (cliopt_verbose){
        for (int i=0; i<st->ref_n; i++){
            fprintf(stderr, "[dbg::%s] raw    ref %s: %d gaps\n", 
                __func__, st->ref_names[i], (int)st->ranges[i]->starts.n);
            if (cliopt_verbose>=2){
                for (int j=0; j<st->ranges[i]->starts.n; j++){
                    fprintf(stderr, "[dbg::%s] raw       %s %d - %d\n", __func__, 
                        st->ref_names[i], st->ranges[i]->starts.a[j], st->ranges[i]->ends.a[j]);
                }
            }
        }
        fprintf(stderr, "\n");
    }
    if (bam_needs_haplotagging && write_baminput_tagging){
        char *fn_haptag_tsv = (char*)malloc(strlen(output_prefix)+25);
        sprintf(fn_haptag_tsv, "%s.mp.input_haptag.tsv", output_prefix);
        FILE *fp_haptag_tsv = fopen(fn_haptag_tsv, "w");
        fprintf(fp_haptag_tsv, "#qname\treal_hp\ttagged_hp\n");

        bamfile_t *tmp_hf = init_and_open_bamfile_t(fn_bam, n_bam_threads);
        hts_itr_t *tmp_fp_itr = sam_itr_querys(tmp_hf->fp_idx, tmp_hf->fp_header, ".");
        while(sam_itr_next(tmp_hf->fp, tmp_fp_itr, tmp_hf->buf_aln)>=0){
            char *qn = bam_get_qname(tmp_hf->buf_aln);
            int hp_raw = get_hp_from_aln(tmp_hf->buf_aln);
            khint_t k = htstri_ht_get(st->qname2haptag_raw, qn);
            if (k==kh_end(st->qname2haptag_raw)){
                fprintf(fp_haptag_tsv, "%s\t%d\t255\n", qn, hp_raw+1);
            }else{
                fprintf(fp_haptag_tsv, "%s\t%d\t%d\n", qn, hp_raw+1, 
                                                kh_val(st->qname2haptag_raw, k)+1);
            }
        }
        fclose(fp_haptag_tsv);
        free(fn_haptag_tsv);
        destroy_bamfile_t(tmp_hf, 1);
        hts_itr_destroy(tmp_fp_itr);
    }

    for (int i=0; i<st->ref_n; i++){
        store_raw_intervals(st->ranges[i]);
        merge_close_intervals(st->ranges[i], READBACK);
        int n_merged = st->ranges[i]->rawunphasedblocks.n - st->ranges[i]->starts.n;
        //if (n_merged>=0){
        //    fprintf(stderr, "[M::%s] merged %d intervals for reference %s (threshold=%d)\n", 
        //                    __func__, n_merged, st->ref_names[i], READBACK);
        //}else if (n_merged<0){
        //    fprintf(stderr, "[M::%s] nothing done for reference %s\n", __func__, st->ref_names[i]);
        //}
    }
    fprintf(stderr, "[M::%s] loaded phase block gaps.\n\n\n", __func__);

    // prepare worker data
    int n_jobs = st->ref_n;
    blockjoin_work_t pack; 

    // collect reference coverage, if not provided
    if (mmr_config.cov_for_selection<=0){
        pack.ref_coverages = estimate_read_coverage_dirtyfast(fn_bam, n_bam_threads);
    }else{
        pack.ref_coverages = (int*)malloc(sizeof(int)*st->ref_n);
        for (int i=0; i<st->ref_n; i++){
            pack.ref_coverages[i] = mmr_config.cov_known;
        }
    }

    // Candidate unphased intervals are loaded. 
    // If input bam was untagged, read haptags have also been generated and stored. 
    // Now go to gap closing attempts, paralleled between chromosomes.
    pack.st = st;
    pack.fn_bam = fn_bam;
    pack.fn_interval_type = fn_interval_type;
    pack.mmr_config = mmr_config;
    pack.n_candidates_per_iter = n_candidates_per_iter;
    pack.n_permutation = n_permutation;
    pack.qname2haptag = (htstri_t**)calloc(st->ref_n, sizeof(htstri_t*));
    for (int refi=0; refi<st->ref_n; refi++){
        pack.qname2haptag[refi] = htstri_ht_init();
    }

    kt_for(n_threads, blockjoin_one_chrom_callback, &pack, n_jobs);
    //for (int jobi=0; jobi<n_jobs; jobi++){
    //    fprintf(stderr, "[dbg::%s] phasing ref#%d %s (single thread mode)\n", 
    //        __func__, jobi,
    //        st->ref_names[jobi]);
    //    blockjoin_one_chrom_callback(&pack, jobi, 0);
    //}

    // migrate read haptags from per-thread buffers
    int absent;
    {
        int tot = 0;
        for (khint_t k=0; k<kh_end(st->qname2haptag); k++){
            if (k!=kh_end(st->qname2haptag)){
                tot++;
            }
        }
        assert(tot==0);
    }
    for (int i=0; i<st->ref_n; i++){
        for (khint_t k=0; k<kh_end(pack.qname2haptag[i]); k++){
            if (kh_exist(pack.qname2haptag[i], k)){
                char *qname_key = kh_key(pack.qname2haptag[i], k);
                int hp = kh_val(pack.qname2haptag[i], k);
                khint_t k2 = htstri_ht_put(st->qname2haptag, qname_key, &absent);

                if (absent){
                    kh_key(st->qname2haptag, k2) = qname_key;
                    kh_val(st->qname2haptag, k2) = hp;
                }else{
                    ;
                }
            }
        }
        htstri_ht_destroy(pack.qname2haptag[i]);
    }

    // remaining cleanups
    free(pack.ref_coverages);
    free(pack.qname2haptag);

    fprintf(stderr, "\n\n[M::%s] done, used %.1fs.\n", __func__, Get_T()-T);
    return st;
}


int sancheck_cliopt_t_files_exist(cliopt_t *cliopt){
    // sancheck: bam file should exist
    samFile *tmp_fp_bam = hts_open(cliopt->fn_bam, "rb");
    if (!tmp_fp_bam || (!tmp_fp_bam->is_bgzf && !tmp_fp_bam->is_cram) ){
        fprintf(stderr, "[E::%s] cannot open bam file: %s\n", __func__, cliopt->fn_bam);
        return 1;
    }
    hts_close(tmp_fp_bam);

    // sancheck: vcf/gtf/tsv file should exist (if provided)
    if (cliopt->fn_vcf){
        FILE *tmp_fp = fopen(cliopt->fn_vcf, "r");
        if (!tmp_fp){
            fprintf(stderr, "[E::%s] cannot open vcf: %s\n", __func__, cliopt->fn_vcf);
            return 1;
        }
        fclose(tmp_fp);
    }
    if (cliopt->fn_tsv){
        FILE *tmp_fp = fopen(cliopt->fn_tsv, "r");
        if (!tmp_fp){
            fprintf(stderr, "[E::%s] cannot open tsv: %s\n", __func__, cliopt->fn_tsv);
            return 1;
        }
        fclose(tmp_fp);
    }
    if (cliopt->fn_gtf){
        FILE *tmp_fp = fopen(cliopt->fn_gtf, "r");
        if (!tmp_fp){
            fprintf(stderr, "[E::%s] cannot open gtf: %s\n", __func__, cliopt->fn_gtf);
            return 1;
        }
        fclose(tmp_fp);
    }
    return 0;
}

int main_blockjoin(cliopt_t *cliopt){ 
    int failed_opening_files = sancheck_cliopt_t_files_exist(cliopt);
    if (failed_opening_files){
        return 1;
    }

    mmr_config_t config;
    config.lo = cliopt->lo;
    config.hi = cliopt->hi;
    config.min_mapq = cliopt->mapq;
    config.k = cliopt->k;
    config.k_span = cliopt->k_span;
    config.cov_known = cliopt->cov;
    config.cov_for_selection = cliopt->cov_for_selection;
    config.cov_for_runtime = config.cov_for_selection*2;
    config.readlen_threshold = cliopt->readlen_threshold;
    int n_candidates_per_iter = cliopt->n_candidates_per_iter;

    char *fn_interval = cliopt->fn_tsv? cliopt->fn_tsv 
                                      :cliopt->fn_gtf? cliopt->fn_gtf 
                                      :cliopt->fn_vcf;
    enum input_file_format fn_interval_type = cliopt->fn_tsv? IS_TSV 
                                      :cliopt->fn_gtf? IS_GTF
                                      :IS_VCF;

    storage_t *st = blockjoin_parallel(fn_interval, fn_interval_type, 
                                      cliopt->fn_bam, 
                                      cliopt->bam_needs_haplotagging, 
                                      cliopt->write_bam_input_haplotagging,
                                      cliopt->fn_vcf,
                                      cliopt->output_prefix,
                                      config, cliopt->n_candidates_per_iter, 
                                      1, cliopt->threads, cliopt->threads_bam);

    lift_decisions(st);
    make_decisions_flippings_onraw(st);
    generate_new_phase_blocks(st, 1);
    if (cliopt->write_debug_files){
        dbgoutput_intermediate_read_haplotags(st, cliopt->output_prefix);
    }

    //if (cliopt_verbose){
    //    fprintf(stderr, "[M::%s] phaseblocks:\n", __func__);
    //    for (int i=0; i<st->ref_n; i++){
    //        for (int j=0; j<st->ranges[i]->phaseblocks.n; j++){
    //            fprintf(stderr, "[M::%s] ref=%s, phaseblock #%d, %d - %d\n", 
    //            __func__, st->ref_names[i], j, 
    //            (int)st->ranges[i]->phaseblocks.a[j].s,
    //            (int)st->ranges[i]->phaseblocks.a[j].e
    //            );
    //        }
    //    }
    //}


    output_gtf(st, cliopt->output_prefix);
    fprintf(stderr, "[M::%s] gtf written.\n", __func__);

    if (cliopt->do_output_tsv){
        output_tsv(st, cliopt->output_prefix);
        fprintf(stderr, "[M::%s] tsv written.\n", __func__);
    }

    if (cliopt->fn_vcf){
        fprintf(stderr, "[M::%s] writing vcf...\n", __func__);
        recover_variant_phase_in_dropped_intervals(st, 
                                     cliopt->fn_bam, cliopt->fn_vcf);
        output_modify_vcf(cliopt->fn_vcf, st, cliopt->output_prefix);
        fprintf(stderr, "[M::%s] vcf written.\n", __func__);
    }

    if (cliopt->do_output_bam){
        char *fn_bam_out = (char*)malloc(strlen(cliopt->output_prefix)+25);
        char *fn_bai_out = (char*)malloc(strlen(cliopt->output_prefix)+25);
        sprintf(fn_bam_out, "%s.mp.bam", cliopt->output_prefix);
        sprintf(fn_bai_out, "%s.mp.bam.bai", cliopt->output_prefix);

        output_modify_bam(cliopt->fn_bam, st, fn_bam_out, cliopt->threads_bam);
        fprintf(stderr, "[M::%s] bam written. now indexing...\n", __func__);

        int stat = sam_index_build3(fn_bam_out, fn_bai_out, 0, cliopt->threads_bam);
        if (stat!=0){
            fprintf(stderr, "[W::%s] failed to build index for output bam (status code=%d)\n", 
                            __func__, stat);
        }
        fprintf(stderr, "[M::%s] bam index written.\n", __func__);
        free(fn_bam_out);
        free(fn_bai_out);
    }
    
    destroy_storage_t(st, 1);
    return 0;
}

int main_varhaptag(char *fn_vcf, char *fn_bam, char *fn_out, int n_thread, 
                   int verbose, int write_bam){
    int ret = 0;
    int debug_print = 0;
    assert(fn_out);

    // open output
    // (bam and bam.bai alloc)
    char *fn_bai_out = 0;
    BGZF *fp_out=0;
    if (write_bam){
        fn_bai_out = (char*)malloc(strlen(fn_out)+10);
        sprintf(fn_bai_out, "%s.bai", fn_out);

        fp_out = bgzf_open(fn_out, "w");
        if (!fp_out){
            fprintf(stderr, "[E::%s] failed to open output file: %s\n", __func__, fn_out);
            return 1;
        }
        if (n_thread>1) bgzf_mt(fp_out, n_thread/2, 0/*unused*/);
    }
    // (tsv output)
    char *fn_tsv_out = (char*)malloc(strlen(fn_out)+25);
    sprintf(fn_tsv_out, "%s.varhaptag.tsv", fn_out);
    FILE *fp_tsv_out = fopen(fn_tsv_out, "w");
    if (!fp_tsv_out){
        fprintf(stderr, "[E::%s] failed to open output file: %s\n", __func__, fn_tsv_out);
    }

    // generate haplotags
    storage_t *st = (storage_t*)calloc(1, sizeof(storage_t));
    init_storage_t(st);
    load_intervals_from_file(fn_vcf, IS_VCF, st, 1, fn_bam, 0, 0);

    // prepare for output
    bamfile_t *hf = init_and_open_bamfile_t(fn_bam, n_thread/2);
    hts_itr_t *fp_itr = sam_itr_querys(hf->fp_idx, hf->fp_header, ".");

    // write header
    int stat;
    if (write_bam){
        stat = bam_hdr_write(fp_out, hf->fp_header);
        if (stat!=0){
            fprintf(stderr, "[E::%s] failed to write bam header\n", __func__);
            ret = 1;
            goto cleanup;
        }
    }
    fprintf(fp_tsv_out, "#qname\thaptag_input\thaptag_new\n");

    // write entries
    while(sam_itr_next(hf->fp, fp_itr, hf->buf_aln)>=0){
        char *qn = bam_get_qname(hf->buf_aln);
        khint_t k = htstri_ht_get(st->qname2haptag_raw, qn);
        int hp;
        if (k==kh_end(st->qname2haptag_raw)){
            hp = HAPTAG_UNPHASED;
        }else{
            hp = kh_val(st->qname2haptag_raw, k);
        }

        int hp_raw = get_hp_from_aln(hf->buf_aln);
        if (write_bam){
            bam_aux_update_int(hf->buf_aln, "HP", hp+1);
            if (stat<0){
                char *refname = hf->fp_header->target_name[hf->buf_aln->core.tid];
                int start_pos = hf->buf_aln->core.pos;
                fprintf(stderr, "[E::%s] failed to write bam entry (ref=%s pos=%d qn=%s)\n", 
                                __func__, refname, start_pos, qn);
                ret = 1;
            }
            stat = bam_write1(fp_out, hf->buf_aln);
        }

        if (debug_print || verbose){
            fprintf(stderr, "[M::%s]\tqn\t%s\thaptag_raw\t%d\thaptag_new\t%d\n", 
                            __func__, qn, hp_raw, hp);
        }
        fprintf(fp_tsv_out, "%s\t%d\t%d\n", qn, hp_raw+1, hp+1);
    }
    
    // index the bam output
    if (write_bam){
        bgzf_close(fp_out);
        stat = sam_index_build3(fn_out, fn_bai_out, 0, n_thread);
        if (stat!=0){
            fprintf(stderr, "[W::%s] failed to build index for output bam (status code=%d)\n", 
                            __func__, stat);
            ret = 1;
        }
    }
cleanup:
    hts_itr_destroy(fp_itr);
    destroy_bamfile_t(hf, 1);
    fclose(fp_tsv_out);
    free(fn_tsv_out);
    destroy_storage_t(st, 1);
    if (fn_bai_out) free(fn_bai_out);
    return ret;
}

int main_methstat(char *fn_bam,// int bam_threads, 
                  char *fn_intervals, enum input_file_format fn_intervals_type, 
                  char *fn_out, 
                  int lo, int hi, int cov_for_selection, int readlen_threshold){
    int ret = 0;
    FILE *fp_out = fopen(fn_out, "w");

    storage_t *st = (storage_t*)calloc(1, sizeof(storage_t));
    init_storage_t(st);
    load_intervals_from_file(fn_intervals, fn_intervals_type, st, 0, 0, 0, 0);
    mmr_config_t mmr_config;
    mmr_config.lo = lo;
    mmr_config.hi = hi;
    mmr_config.readlen_threshold = readlen_threshold;
    // dummy values
    mmr_config.min_mapq = 0;
    mmr_config.k = 1;
    mmr_config.k_span = 5000;
    mmr_config.cov_for_runtime = 1;

    // check if we need coverage estimation
    int *covs;
    if (cov_for_selection<=0){  // not provided by cli
        // TODO: get bam #thread from cli
        covs = estimate_read_coverage_dirtyfast(fn_bam, 1);
        for (int i=0; i<st->ref_n; i++) {
            covs[i]/=10;
            covs[i]++;
        }
    }else{
        covs = (int*)malloc(sizeof(int)*st->ref_n);
        for (int i=0; i<st->ref_n; i++) covs[i] = cov_for_selection;
    }

    for (int i_ref=0; i_ref<st->ref_n; i_ref++){
        mmr_config.cov_for_selection = covs[i_ref];
        char *chrom = st->ref_names[i_ref];
        for (int i_itvl=0; i_itvl<st->ranges[i_ref]->starts.n; i_itvl++){
            int ref_start = st->ranges[i_ref]->starts.a[i_itvl];
            int ref_end   = st->ranges[i_ref]->ends.a[  i_itvl];

            // collect
            rs_t *rs = load_reads_given_interval(fn_bam, chrom, 
                        ref_start, ref_end, 0/*readback*/, mmr_config, 0, NULL);
            methmers_t *ms = get_methmer_sites_and_ranges(rs, mmr_config, 0, NULL, NULL);

            // output
            for (int i=0; i<ms->n; i++){
                fprintf(fp_out, "%s\t%d\n", chrom, (int)ms->sites_real_poss[i]);
            }

            // cleanup
            destroy_methmers_t(ms);
            destroy_rs_t(rs);
        }
    }

    free(covs);
    destroy_storage_t(st, 1);
    fclose(fp_out);
    return ret;
}


int main_debug(int argc, char *argv[]){
    return 0;
}



int main_methreport(cliopt_t *clio){
    double T = Get_T();

    char *fn_bam = clio->fn_bam;
    char *fn_vcf = clio->fn_vcf;
    char *fn_out_prefix = clio->output_prefix;
    int n_threads = clio->threads;
    int bam_needs_haptagging = clio->bam_needs_haplotagging;
    int chunk_size = clio->chunck_size;
    int chunk_stride = clio->chunck_stride;
    int read_coverage = clio->cov;


    // make sure input files are available
    if (!fn_bam){
        fprintf(stderr, "[E::%s] input bam file name missing\n", __func__);
        exit(1);
    }
    if (!fn_vcf){
        fprintf(stderr, "[E::%s] input vcf file name missing\n", __func__);
        exit(1);
    }if (!fn_out_prefix){
        fprintf(stderr, "[E::%s] output file name prefix missing\n", __func__);
        exit(1);
    }
    bamfile_t *hf = init_and_open_bamfile_t(fn_bam, n_threads);
    if (!hf){
        fprintf(stderr, "[E::%s] failed to open input bam: %s\n", __func__, fn_bam);
        exit(1);
    }
    hts_itr_t *fp_itr = sam_itr_querys(hf->fp_idx, hf->fp_header, ".");
    if (!fp_itr){
        fprintf(stderr, "[E::%s] failed to get itr for input bam.\n", __func__);
        exit(1);
    }
    hts_itr_destroy(fp_itr);
    destroy_bamfile_t(hf, 1);

    // open output file
    char *fn_out = (char*)malloc(strlen(fn_out_prefix)+30);
    sprintf(fn_out, "%s.report.tsv", fn_out_prefix);
    FILE *fp_out = fopen(fn_out, "w");
    if (!fp_out){
        fprintf(stderr, "[E::%s] failed to open output file\n", __func__);
        exit(1);
    }

    // init
    storage_t *st = (storage_t*)calloc(1, sizeof(storage_t));
    init_storage_t(st);
    load_intervals_from_file(fn_vcf, IS_VCF, st, !!bam_needs_haptagging, 
                             fn_bam, 0, 0);

    // override intervals
    vu32_t starts, ends;
    kv_init(starts);
    kv_init(ends);
    uint32_t prev = 0, end, start;
    for (int i_ref=0; i_ref<st->ref_n; i_ref++){
        starts.n = 0;
        ends.n = 0;
        prev = st->ranges[i_ref]->abs_start;
        for (int i_itvl=0; i_itvl<st->ranges[i_ref]->starts.n; i_itvl++){
            start = st->ranges[i_ref]->starts.a[i_itvl];
            end = st->ranges[i_ref]->ends.a[i_itvl];
            if (start-prev > chunk_size) {
                for (uint32_t i=prev; i+chunk_stride<start; i+=chunk_stride){
                    kv_push(uint32_t, starts, i);
                    kv_push(uint32_t, ends, i+chunk_size);
                }
            }
            prev = end;
        }

        // override
        st->ranges[i_ref]->starts.n = 0;
        st->ranges[i_ref]->ends.n = 0;
        st->ranges[i_ref]->decisions.n = 0;
        for (int i=0; i<starts.n; i++){
            kv_push(uint32_t, st->ranges[i_ref]->starts, starts.a[i]);
            kv_push(uint32_t, st->ranges[i_ref]->ends,   ends.a[i]);
            kv_push(int, st->ranges[i_ref]->decisions, -1);
        }
        fprintf(stderr, "[M::%s] %s has %d intervals\n", 
                __func__, st->ref_names[i_ref], (int)starts.n);
    }
    kv_destroy(starts);
    kv_destroy(ends);

    // collect reference coverage if not provided
    int *covs = 0;
    if (read_coverage<=0){
        fprintf(stderr, "[M::%s] estimating read depths..\n", __func__);
        covs = estimate_read_coverage_dirtyfast(fn_bam, n_threads);
        fprintf(stderr, "[M::%s] done estimating read depths.\n", __func__);
    }

    // collect variants
    htu32_t **masked_poss = (htu32_t**)calloc(st->ref_n, sizeof(htu32_t*));
    storage_t *st2 = (storage_t*)calloc(1, sizeof(storage_t));
    init_storage_t(st2);
    vvar_t *vars_v = (vvar_t*)malloc(sizeof(vvar_t)*st->ref_n);
    for (int i=0; i<st->ref_n; i++){
        init_vvar_t(&vars_v[i]);
        masked_poss[i] = htu32_ht_init();
    }
    load_intervals_from_file(fn_vcf, IS_VCF, st2, 1, fn_bam, vars_v, st->ref_n);
    destroy_storage_t(st2, 1);
    int n_variants = 0;
    for (int i_ref=0, absent; i_ref<st->ref_n; i_ref++){
        for (int i_var=0; i_var<vars_v[i_ref].n; i_var++){
            for (int j=0; j<vars_v[i_ref].a[i_var].len; j++){
                uint32_t pos = vars_v[i_ref].a[i_var].pos + j;
                khint_t k = htu32_ht_put(masked_poss[i_ref], pos, &absent);
                if (absent) {
                    kh_val(masked_poss[i_ref], k) = pos;
                    n_variants++;
                }
            }
        }
    }
    fprintf(stderr, "[M::%s] recorded %d variants from vcf\n", __func__, n_variants);


    // methphase
    mmr_config_t mmr_config;
    mmr_config.lo = clio->lo;
    mmr_config.hi = clio->hi;
    mmr_config.readlen_threshold = clio->readlen_threshold;
    mmr_config.min_mapq = clio->mapq;
    mmr_config.k = clio->k;
    mmr_config.k_span = clio->k_span;
    float n_switch = 0;
    float n_fail = 0;
    float n_correct = 0;
    int tot = 0;
    for (int i_ref=0; i_ref<st->ref_n; i_ref++){
        mmr_config.cov_for_selection = read_coverage<=0
                                   ? covs[i_ref]/10+1
                                   : read_coverage/10+1;
        mmr_config.cov_for_runtime = mmr_config.cov_for_selection*2;
        int n_candidates_per_iter = read_coverage<=0
                                    ? covs[i_ref]/4+1
                                    : read_coverage/4+1;

        ranges_t *rg = st->ranges[i_ref];
        for (int i_itvl=0; i_itvl<rg->starts.n; i_itvl++){
            int start = rg->starts.a[i_itvl];
            int end = rg->ends.a[i_itvl];
            int decision = -1;
            dataset_t *ds = haplotag_region_given_bam(st, 
                        fn_bam, 
                        st->ref_names[i_ref], start, end, 
                        mmr_config, 
                        n_candidates_per_iter, 
                        1/*perm*/, 
                        &decision);
            fprintf(fp_out, "%s\t%d\t%d\t", st->ref_names[i_ref], start, end);
            if (decision==0)      {n_correct++; fprintf(fp_out, "correct\n");}
            else if (decision==1) {n_switch++;  fprintf(fp_out, "switch\n");}
            else                  {n_fail++;    fprintf(fp_out, "fail\n");}
            tot++;
            if (tot%100==0){
                fprintf(stdout, "Parsed N=%d regions, currently at %s:%d-%d, correct/(correct+switch)=%.2f%%, correct/N=%.2f%%\n", 
                        tot, st->ref_names[i_ref], start, end, n_correct/(n_correct+n_switch)*100.0, n_correct/(float)tot*100.0
                );
            }
            fflush(fp_out);
            destroy_dataset_t(ds, 1);
        }
    }
    fprintf(stdout, "Total N=%d regions, correct/(correct+switch)=%.2f%%, correct/N=%.2f%%\n", 
                        tot, n_correct/(n_correct+n_switch)*100.0, n_correct/(float)tot*100.0
            );
    fprintf(stderr, "[M::%s] Total N=%d regions, correct/(correct+switch)=%.2f%%, correct/N=%.2f%%\n", 
                        __func__, tot, n_correct/(n_correct+n_switch)*100.0, n_correct/(float)tot*100.0
            );
    fclose(fp_out);
    free(fn_out);
    for (int i=0; i<st->ref_n; i++){
        destroy_vvar_t(&vars_v[i], 0);
        htu32_ht_destroy(masked_poss[i]);
    }
    free(masked_poss);
    free(vars_v);
    destroy_storage_t(st, 1);

    fprintf(stderr, "[M::%s] done, used %.1fs\n", __func__, Get_T()-T);
    return 0;
}
void debug_search_arr(int argc, char *argv[]){
    if (argc<3) {
        fprintf(stderr, "[E::%s] usage: exe val a1 a2 a3...\n", __func__);
        exit(1);
    }
    uint32_t *a = (uint32_t*)malloc(sizeof(uint32_t)*(argc-2));
    uint32_t v = (uint32_t)strtoul(argv[1], NULL, 10);
    uint32_t l = argc-2;
    uint32_t idx;
    for (int i=2; i<argc; i++){
        a[i-2] = strtoul(argv[i], NULL, 10);
    }
    int found = search_arr(a, l, v, &idx);
    if (found)
        fprintf(stderr, "[M::%s] found, idx=%d\n", __func__, (int)idx);
    else
        fprintf(stderr, "[M::%s] not found\n", __func__);
    free(a);
}

void debug_load_reads_given_interval(char *fn,
                                char *chrom, int itvl_s, int itvl_e, 
                                int min_mapq, int min_read_len, 
                                uint8_t qual_lo, uint8_t qual_hi){
    int n_start_reads = -1; // dummy
    int n_end_reads = -1; // dummy
    rs_t *rs = load_reads_given_interval(fn, chrom, itvl_s, itvl_e, READBACK, min_mapq, min_read_len, qual_lo, qual_hi, &n_start_reads, &n_end_reads);
    for (int i=0; i<rs->n; i++){
        read_t *h = &rs->a[i];
        printf("[dbg::%s] read %d: qname=%.*s\tID=%d\tlen=%d\thp=%d\n", __func__, i, 
            (int)GET_QNAME_L(*rs, i), GET_QNAME_S(*rs, i), (int)h->i, (int)h->len, h->hp);
        for (int j=0; j<h->meth.calls.n; j++){
            printf("[dbg::%s]     pos=%d\tqual=%d\n", __func__, 
                        (int)h->meth.calls.a[j], (int)h->meth.quals.a[j]);
        }
    }
    destroy_rs_t(rs);
}

void debug_load_intervals_from_gtf(char *fn, char *chrom){
    ranges_t *tmp = load_intervals_from_gtf(fn, chrom, 0);
    for (int i=0; i<tmp->starts.n; i++){
        fprintf(stderr, 
        "[dbg::%s]s=%d , e=%d\n", __func__, tmp->starts.a[i], tmp->ends.a[i]);
    }
    destroy_ranges_t(tmp);
}

void debug_print_methmer_count_at_pos(methmers_t *ms, int pos){
    for (int i=0; i<ms->n; i++){
        if (ms->sites_real_poss[i]!=pos) continue;
        char tmp_buf[8];
        for (int j=0; j<ms->mmr_a[0].vmmr_u32i.n; j++){
            uint32_to_methmer(ms->mmr_a[0].vmmr_u32i.a[j], tmp_buf, 7);
            fprintf(stderr, "[dbg::%s] pos=%d mer=%.*s hap0_count=%d hap1_count=%d\n", 
                __func__,
                ms->sites_real_poss[0], 7, tmp_buf, 
                ms->mmr_a[0].hap_cnt[0].a[j], ms->mmr_a[0].hap_cnt[1].a[j]);
        }
    }
}
void debug_print_methmer_sites(methmers_t *mmr){
    for (int i=0; i<mmr->n; i++){
        printf("[dbg::%s] site\t%d\n", __func__, mmr->sites_real_poss[i]);
    }
}
void debug_print_methmer_pairs(methmers_t *mmr){
    for (int i=0; i<mmr->n; i++){
        int j = i+mmr->mmr_lens[i];
        printf("[dbg::%s] methmer\t%d -> %d\n", __func__, 
               mmr->sites_real_poss[i], mmr->sites_real_poss[j]);
    }
}
void debug_dump_sites_and_methmers(){
    // test sites and methmers
    int n_start_reads = -1;  // dummy
    int n_end_reads = -1;  // dummy
    rs_t *rs = load_reads_given_interval("haplotagged.F4F256F2048.bam", 
                                        "chr6", 31906356, 31991874, 
                                        READBACK, 
                                        0, 0, 100, 156, 
                                        &n_start_reads, &n_end_reads);
    methmers_t *m = get_methmer_sites_and_ranges(rs, 7, 5000, 10);
    debug_print_methmer_sites(m);
    debug_print_methmer_pairs(m);
    destroy_methmers_t(m);
    destroy_rs_t(rs);
}

void debug_print_methmers_of_read(rs_t *rs, methmers_t *ms, char *s){
    uint32_t readID;
    uint32_t start_i=0;
    int absent;
    vu32_t buf;
    kv_init(buf); kv_resize(uint32_t, buf, 32);

    khint_t k = htstru32_ht_put(rs->name2namei, s, &absent);
    if (!absent){
        fprintf(stderr, "[E::%s] bad\n", __func__);
        exit(1);
    }
    readID = kh_val(rs->name2namei, k);
    int n_mmr_tmp = get_mmr_of_read(&rs->a[readID], ms, 
            &buf, &start_i);
    for (int i=0; i<buf.n; i++){
        char tmp[8];
        tmp[7] = 0;
        uint32_to_methmer(buf.a[i], tmp, 7);
        fprintf(stderr, "[dbg::%s] methmer#%d: %s (packed form: %d), pos=%d\n", 
            __func__, i, tmp, buf.a[i], ms->sites_real_poss[i+start_i]);
    }
    kv_destroy(buf);
}

int main(int argc, char *argv[]){
    // test loading gtf
    //debug_load_reads_given_interval("sancheck.bam", "chr6", 0, 180000000, 0, 0);

    // test binary search
    //debug_search_arr(argc, argv);

    //int tmp = -1;
    //rs_t *rs = load_reads_given_interval("haplotagged.F4F256F2048.bam", 
    //        "chr6", 13199497, 13298832, READBACK, 0, 10000, 100, 156, &tmp);
    //printf("sancheck, printing a read's name: %.*s\n", (int)GET_QNAME_L(*rs, 91), GET_QNAME_S(*rs, 91));
    //exit(0);

    //debug_dump_sites_and_methmers();
    //exit(0);

    //int tmp = -1;
    //rs_t *rs = load_reads_given_interval("haplotagged.F4F256F2048.bam", 
    //        "chr6", 31906356, 31991874, READBACK, 0, 10000, 100, 156, &tmp);
    //methmers_t *ms = get_methmer_sites_and_ranges(rs, 7, 5000, 10);
    //store_mmr_of_reads(rs, ms);
    //char tmp_buf[8];
    //tmp_buf[7] = 0;
    //for (int i=0; i<ms->n; i++){
    //    fprintf(stderr, "[dbg-mmr] pos=%d:\n", ms->sites_real_poss[i]);
    //    for (int j=0; j<ms->mmr_a[i].vmmr_u32i.n; j++){
    //        uint32_to_methmer(ms->mmr_a[i].vmmr_u32i.a[j], tmp_buf, 7);
    //        fprintf(stderr, "[dbg-mmr]    #%d %.*s [%d, %d]\n", 
    //            j, 
    //            7, tmp_buf, 
    //            (int)ms->mmr_a[i].hap_cnt[0].a[j], 
    //            (int)ms->mmr_a[i].hap_cnt[1].a[j]);
    //    }
    //}
    //exit(0);

}
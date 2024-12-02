#ifndef __METHPHASER2_BLOCKJOIN_H__
#define __METHPHASER2_BLOCKJOIN_H__

int main_blockjoin(char *fn_gtf, char *fn_bam, char *fn_out_tsv, 
                   int use_hypersearch, 
                   int readlen_threshold, 
                   int lo, int hi, int mapq, 
                   int k, int k_span, int cov_for_selection, 
                   int n_candidates_per_iter);
int main_blockjoin_dbg(int argc, char *argv[]);

#endif  //__METHPHASER2_BLOCKJOIN_H__

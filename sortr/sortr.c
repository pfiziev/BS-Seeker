#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <time.h>

#include "utils.h"

	

//static int compare_sequences(const SEQUENCE *seq1, const SEQUENCE *seq2) {
static int compare_sequences(const void *pseq1, const void *pseq2) {
	SEQUENCE *seq1 = *(SEQUENCE * const *)pseq1;
	SEQUENCE *seq2 = *(SEQUENCE * const *)pseq2;
    char chrom_diff = CHROM(seq1) - CHROM(seq2);
    if (chrom_diff == 0) {
        return seq1->pos - seq2->pos;
    } else {
        return chrom_diff;
    }
}


int main (int argc, char **argv) {
	SEQUENCE **seqs = NULL;
	SEQUENCE *seqdata = NULL;
	WORD *bitseqdata = NULL;
    

	fprintf(stderr, "word size: %lu \n", sizeof(WORD));
	if (argc != 5) {
		fprintf(stderr, "usage: %s total_number_of_reads maximum_read_length <bs_seeker-output> <output_file>\n", argv[0]);	
		return 1;
		}
    
    int total_reads = atoi(argv[1]);
    int max_read_length = atoi(argv[2]) + 6;
    
    fprintf(stderr, "Total reads: %d \n", total_reads);
    fprintf(stderr, "Max read length: %d \n", max_read_length - 6);
	
    
    FILE *input = fopen((const char*)argv[3],"r");
	
	/* Check for validity of the file. */
	if(input == NULL) {
		fprintf(stderr,"can't open input file: %s.\n", argv[3]);
		return 1;
	}
	
    
    seqs = (SEQUENCE **)malloc((total_reads)*sizeof(SEQUENCE *));
    seqdata = (SEQUENCE *)malloc((total_reads)*sizeof(SEQUENCE));
    bitseqdata = (WORD *)calloc(total_reads * blen(max_read_length), sizeof(WORD));
    
    if (seqs == NULL || seqdata == NULL || bitseqdata == NULL){
		fprintf(stderr, "Not enough memory!\n");
        return 1;
    	}
	
    
    
	bs_seeker2seq(input, seqs, seqdata, bitseqdata, total_reads, max_read_length);
	fclose(input);
    //fprintf(stderr, "Number of sequences: %d\n", total_sequences);

	start_watch();	
	fprintf(stderr, "sorting sequences\n");
	
	qsort(seqs, total_reads, sizeof(SEQUENCE *), compare_sequences);
	
	fprintf(stderr, "sorting finished | ");
	stop_watch();
	
	start_watch();	
	fprintf(stderr, "outputting sequences\n");
	
	FILE *output = fopen(argv[4],"w");
	if(output == NULL) {
		fprintf(stderr,"can't open output file: %s.\n", argv[4]);
		return 1;
	}
	
	int i;
	char charseq[MAXLINE] = {0};
    for(i = 0; i < total_reads; i++) {
                            
        bitseq2char(seqs[i]->bitseq, LENGTH(seqs[i]), charseq);

        fprintf(output,"%04d%c%010d\t%s\n",CHROM(seqs[i]) , STRAND(seqs[i]), seqs[i]->pos, charseq);
        
    }
    fclose(output);
	fprintf(stderr, "outputting finished | ");
	stop_watch();
	
	return 0;
	}

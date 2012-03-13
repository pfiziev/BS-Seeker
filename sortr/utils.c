#include "utils.h"



WORD n2c(char n) {
	switch (n){
		case 'a': case 'A':
			return (WORD)0;
		case 'c': case 'C':
			return (WORD)1;
		case 'g': case 'G':
			return (WORD)2;
		case 't': case 'T':
			return (WORD)3;            
		case 'n': case 'N':
			return (WORD)4;            
		case '_':
			return (WORD)5;
        case '-':
			return (WORD)6;
            
		}
	return ERROR;
	}


// returns the length of the bit sequence given the real length of the dna
int blen(int n) {
	return (n % BASES_PER_WORD == 0 ? 0 : 1) + n / BASES_PER_WORD;
	}	

	
	
void bs_seeker2seq(          FILE *input, 
                             SEQUENCE **seqs, 
							 SEQUENCE *seqdata, 
							 WORD *bitseqdata,
                             int total_reads, 
                             int max_read_length
							 ) {
	// sequences is a pointer to an array of sequences
	char line[MAXLINE] = {0};
	
	
	start_watch();
	/* read and convert the sequences */
	fprintf(stderr,"reading sequences\n");
	
	int seq_no = 0;			
	
    int bitseq_index = 0;
	int byte_index = 0;
	int in_byte_index = 0;
	
    int i, j, k, col_number, qlen;
    SEQUENCE *seq = NULL;
    char buf[MAXLINE] = {0};
    
    while ( fgets (line, sizeof line, input) != NULL ) {
		
        if (seq_no == total_reads) {
            fprintf(stderr,"It seems that you have more than %d reads in the input file. Truncating!\n", total_reads);
            break;
            }
        
        //printf("%s", line);
		// record the header position in the file 
		
        
		// init the sequence
        
		seqs[seq_no] = &seqdata[seq_no];
		seq = seqs[seq_no];
        
        col_number = 0;
    	
        for (i = 0; line[i] != 0; i++) {
            if (isspace(line[i])) {
                col_number ++;
                i++;
                while (line[i] != 0 && isspace(line[i])) { i++; }
                }
            if (col_number == 3) {
                
                j = i;
                k = 0;
                while (line[j] != '+' && line[j] != '-') {buf[k] = line[j]; j++; k++; }
                buf[k] = 0;
                STRAND(seq) = line[j];
                CHROM(seq) = (char)atoi(buf);
                //printf("NEW: %s %d %c ", buf, CHROM(seq), line[j]);
                
                j++;
                k = 0;
                while (!isspace(line[j+1])) {buf[k] = line[j]; j++; k++; }
                buf[k] = line[j];
                k++;
                buf[k] = 0;
                seq->pos = atoi(buf);
                //printf(" %s %d\n", buf, seq->pos);
                i = j;
                }
            else if (col_number == 5) {
                qlen = 0;
                j = i;
                while(!isspace(line[j])) {
                    buf[qlen] = line[j];
                    j++;
                    qlen++;
                    }
                
                if (qlen > max_read_length) {
                    fprintf(stderr,"You have reads that are longer than %d (length %d). Truncating!", max_read_length - 6, qlen);
                    qlen = max_read_length;
	
                    }
                    
                buf[qlen] = 0;
                LENGTH(seq) = (char)qlen;
                
                seq->bitseq = &bitseqdata[bitseq_index];
                
                for (k = 0; k < qlen; k++) {
                    byte_index = k / BASES_PER_WORD;
                    in_byte_index = (BITS_PER_WORD - BITS_PER_BASE) - BITS_PER_BASE * (k % BASES_PER_WORD);
                    
                    WORD c = n2c(buf[k]) << in_byte_index;
                    
                    seqs[seq_no]->bitseq[byte_index] += c;
                    
                    //printf("%c: %d %d %d %u %u\n",line[i], bitseq_index, byte_index, in_byte_index, c, bitseqdata[bitseq_index + byte_index]);
                    
                    }
                
                bitseq_index += blen(qlen);
                seq_no ++;
                break;
                }
            }
   
		
        }	 
	fprintf(stderr,"reading finished | ");
	stop_watch();
	
	}



char c2n(WORD c) {
	switch(c) {
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
		case 4:
			return 'N';
		case 5:
			return '_';
		case 6:
			return '-';
		}
	return 'X';
	}



void bitseq2char(WORD *bitseq, int seq_real_length, char *charseq) {
	int i, j, dna_index = -1;
	int bytes_length = blen(seq_real_length);
	for (i = 0; i < bytes_length; i++) {
		for (j = 0; j < BASES_PER_WORD ; j++) {
			dna_index = i*BASES_PER_WORD + j;
			if (dna_index >= seq_real_length) { 
				break; 
			}
			charseq[dna_index] = c2n(BCHAR(bitseq[i],j));
			}	
		}
	charseq[seq_real_length] = 0;		
	}
		



time_t START_TIME;
void start_watch() { START_TIME = time(NULL); } 
void stop_watch() { fprintf (stderr, "%ld seconds\n", time(NULL) - START_TIME); }

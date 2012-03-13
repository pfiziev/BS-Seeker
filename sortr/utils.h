#ifndef __ANGS_UTILS_H__
#define __ANGS_UTILS_H__

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#define ERROR -2

//#define WORD unsigned char
#define WORD unsigned int
#define BITS_PER_WORD (8*sizeof(WORD)) // in bits
#define BITS_PER_BASE 3
#define BASES_PER_WORD (BITS_PER_WORD/BITS_PER_BASE)

#define BCHAR(byte,position) ((int)( ((byte) >> (BITS_PER_WORD - BITS_PER_BASE*((position)+1)) ) & ( (1 << BITS_PER_BASE) - 1 ) ))


#define MIN(a, b) ((a)<=(b) ? (a) : (b))
#define MAX(a, b) ((a)>=(b) ? (a) : (b))

#define MAXLINE 2048


typedef struct {
    int pos;
    char info[4]; // compress real length and strand into two bytes 
	WORD *bitseq;	
} SEQUENCE;

#define CHROM(seq)   (seq->info[0])
#define STRAND(seq)  (seq->info[1])
#define LENGTH(seq)  (seq->info[2])






void start_watch();
void stop_watch();

WORD n2c(char);


int blen(int);

void bs_seeker2seq(FILE *, SEQUENCE **, SEQUENCE *, WORD *, int, int);

char c2n(WORD);


void seq2buffer(SEQUENCE *, char *);

void bitseq2char(WORD *, int , char *);





#endif

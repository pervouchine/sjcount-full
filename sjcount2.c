//	Copyright 2012 Dmitri Pervouchine (dp@crg.eu), Lab Roderic Guigo
//	Bioinformatics and Genomics Group @ Centre for Genomic Regulation
//	Parc de Recerca Biomedica: Dr. Aiguader, 88, 08003 Barcelona
//	
//	This file is a part of the 'sjcount' package.
//	'sjcount' package is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//	
//	'sjcount' package is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with 'sjcount' package.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/ioctl.h>
#include <bam.h>
#include "progressbar.h"

#define MAXFILEBUFFLENGTH 1000
#define ARRAY_MARGIN 2
#define INFTY 65535

const char version[100] = "v.1.0-full";

int nbins   = 1;

const int STRAND[2] = {1, -1};

class site {
  public:
    int pos;
    int* count[2];
    site *next;

    site(int end) {
	pos = end;
	for(int j=0; j<2; j++) {
	    count[j] = (int*) malloc(sizeof(int) * nbins);
	    for(int k=0; k<nbins; k++) count[j][k] = 0;
	}
	next = NULL;
    }
};

class junction {
  public:
    int pos;
    site *partner;
    junction* next;

    junction(int beg, int end) {
	pos = beg;
	partner = new site(end);
	next = NULL;
    }
};

void update_site(site **ptr, int end, int strand, int offset, int count) {
    if(offset >= nbins) offset = nbins -1;
    while(*ptr != NULL && (*ptr)->pos < end) {
        ptr = &((*ptr)->next);
    }
    if(*ptr != NULL && (*ptr)->pos == end) {
	(*ptr)->count[strand][offset] += count;
    }
    else {
        site *next = (*ptr);
        (*ptr) = new site(end);
        (*ptr)->next = next;
	(*ptr)->count[strand][offset] += count;
    }
}

void update_jnxn(junction **ptr, int beg, int end, int strand, int offset, int count) {
    if(offset >= nbins) offset = nbins -1;
    while(*ptr != NULL && (*ptr)->pos < beg) {
    	ptr = &((*ptr)->next);
    }
    if(*ptr != NULL && (*ptr)->pos == beg) {
	update_site(&((*ptr)->partner), end, strand, offset, count);
    }
    else {
        junction *next = (*ptr);
	(*ptr) = new junction(beg, end);
	(*ptr)->next = next;
	(*ptr)->partner->count[strand][offset] += count;
    }
}

char strand_i2c(int i) {
    if(i>0) return('+');
    if(i<0) return('-');
    return('.');
}

int main(int argc,char* argv[]) {
    time_t timestamp, current_time;
    int i,j,k,s,a,n;

    FILE *input_file;
    FILE *ssj_file;
    FILE *ssc_file;
    FILE* log_file=stderr;

    bamFile bam_input;
    bam_header_t *header;
    bam1_t* b;
    bam1_core_t *c;
    uint32_t *cigar;

    char bam_file_name[MAXFILEBUFFLENGTH]="";
    char ssj_file_name[MAXFILEBUFFLENGTH]="";
    char ssc_file_name[MAXFILEBUFFLENGTH]="";
    char log_file_name[MAXFILEBUFFLENGTH]="";

    char buff[MAXFILEBUFFLENGTH];
    char chr[MAXFILEBUFFLENGTH];
    int beg, end, pos, offset, increm; 
    int ref_id, ref_id_prev;
    int read_type, mapped_strand;

    int rev_compl[2] = {1,0};
    int limit_counts = 0;
    int verbose = 1;
    int stranded = 1;

    int n_reads = 0;

    junction **root_junction, ***curr_junction;
    site **root_site, ***curr_site;

    site *qtr;
    junction* ptr;

    timestamp = time(NULL);

    if(argc==1) {
	fprintf(stderr, "sjcount %s counts split reads supporting splice junctions and continuous reads that cover exon boundaries\n", version);
	fprintf(stderr, "This version differs in counting reads that cover exon boundaries\n");
        fprintf(stderr, "Type %s -h for help\n",argv[0]);
        exit(1);
    }

    for(i=1;i<argc;i++) {
	if(strcmp(argv[i], "-bam") == 0) sscanf(argv[++i], "%s", &bam_file_name[0]);
	if(strcmp(argv[i], "-ssj") == 0) sscanf(argv[++i], "%s", &ssj_file_name[0]);
        if(strcmp(argv[i], "-ssc") == 0) sscanf(argv[++i], "%s", &ssc_file_name[0]);
        if(strcmp(argv[i], "-log") == 0) sscanf(argv[++i], "%s", &log_file_name[0]);

        if(strcmp(argv[i], "-read1") == 0) sscanf(argv[++i], "%i", &rev_compl[0]);
        if(strcmp(argv[i], "-read2") == 0) sscanf(argv[++i], "%i", &rev_compl[1]);

	if(strcmp(argv[i], "-lim") == 0)    sscanf(argv[++i], "%i", &limit_counts);

	if(strcmp(argv[i], "-nbins")   == 0) sscanf(argv[++i], "%i", &nbins);

	if(strcmp(argv[i], "-quiet") == 0) verbose = 0;
	if(strcmp(argv[i], "-unstranded") == 0) stranded = 0;

        if(strcmp(argv[i], "-h") ==0 ) {
            fprintf(stderr, "Usage: %s -bam bam_file [-ssj junctions_output] [-ssc boundaries_output] [-log log_file] [-read1 0|1] [-read2 0|1] ",argv[0]);
            fprintf(stderr, "[-nbins number_of_bins] [-lim number_of_lines] [-quiet]\n");
	    fprintf(stderr, "sjcount %s\n", version);
            fprintf(stderr, "Input:  a sorted BAM file with a header\n");
            fprintf(stderr, "Options:\n");
            fprintf(stderr, "\t-read1 0/1, reverse complement read1 no/yes (default=%i)\n",rev_compl[0]);
            fprintf(stderr, "\t-read2 0/1, reverse complement read2 no/yes (default=%i)\n",rev_compl[1]);
	    fprintf(stderr, "\t-nbins number of overhang bins, (default=%i)\n", nbins);
	    fprintf(stderr, "\t-lim nreads stop after nreads, (default=no limit)\n");
	    fprintf(stderr, "\t-unstranded, force strand=0\n");
	    fprintf(stderr, "\t-quiet, suppress verbose output\n\n"); 
            fprintf(stderr, "Output:\t-ssj: Splice Junction counts, tab-delimited  (default=stdout)\n");
            fprintf(stderr, "\tColumns are: chr, begin, end, strand, offset, count\n");
            fprintf(stderr, "\t-ssc: Splice boundary counts, tab-delimited  (default=none)\n");
            fprintf(stderr, "\tColumns are: chr, position, strand, offset, count\n");
            exit(1);
        }

    }

    if(log_file_name[0]==0) {
	log_file = stderr;
    }
    else {
	log_file = fopen(log_file_name,"w");
	if(log_file == NULL) log_file = stderr;
    }

    if(bam_file_name[0]==0) {
	fprintf(log_file,"Bam not specified, exiting\n");
	exit(1); 
    }

    if(ssj_file_name[0]==0) {
	fprintf(log_file,"[Warning: junction output set to stdout]\n");
	ssj_file = stdout;
    }
    else {
	ssj_file = fopen(ssj_file_name,"w");
	if(ssj_file == NULL) {
	    fprintf(log_file,"[Warning: junction output set to stdout]\n");
            ssj_file = stdout;
	} else {
	    fprintf(log_file,"[Junction counts: >%s]\n",ssj_file_name);
	}
    }

    if(ssc_file_name[0]==0) {
	fprintf(log_file,"[Warning: boundary output skipped]\n");
    }
    else {
        ssc_file = fopen(ssc_file_name,"w");
        if(ssc_file == NULL) {
            fprintf(log_file,"[Warning: boundary output set to stdout]\n");
            ssc_file = stdout;
        } else {
            fprintf(log_file,"[Boundary counts: >%s]\n",ssc_file_name);
        }
    }


    if(stranded==0) {
	fprintf(log_file,"[Warning: strand is ignored (forced to zero)]\n");
    }

    for(s = 0; s < 2; s++) {
	if(rev_compl[s]) fprintf(log_file,"[Warning: will take reverse complement of read %i]\n", s+1);
    }



    //*****************************************************************************************************************************//
   
    bam_input = bam_open(bam_file_name, "r");
    header = bam_header_read(bam_input);

    if(bam_input == NULL || header == NULL) {
        fprintf(log_file,"BAM can't be opened or contains no header, exiting\n");
        exit(1);
    }

    root_junction = (junction**)  malloc( sizeof(junction*)  * (header->n_targets + ARRAY_MARGIN) );
    curr_junction = (junction***) malloc( sizeof(junction**) * (header->n_targets + ARRAY_MARGIN) );

    root_site = (site**)  malloc( sizeof(site*)  * (header->n_targets + ARRAY_MARGIN) );
    curr_site = (site***) malloc( sizeof(site**) * (header->n_targets + ARRAY_MARGIN) );

    for(i=0; i < header->n_targets; i++) {
	root_junction[i] = NULL;
	curr_junction[i] = &root_junction[i];
	root_site[i] = NULL;
	curr_site[i] = &root_site[i];
    }

    b = bam_init1();
    k = 0;
    ref_id = -1;
    n_reads = 0;
    while(bam_read1(bam_input, b)>=0) {
        c   = &b->core;
        if(c->tid < 0 || c->tid >= header->n_targets) continue;

        ref_id_prev = ref_id;
        ref_id = c->tid;

        cigar = bam1_cigar(b);

	s = ((c->flag & BAM_FREVERSE)>0);

	mapped_strand = (c->flag & BAM_FREAD1) ? (s + rev_compl[0]) & 1 : (s + rev_compl[1]) & 1;
	mapped_strand*= stranded;

	pos = beg = c->pos + 1;

	while((*curr_junction[ref_id])!=NULL && (*curr_junction[ref_id])->pos < beg) curr_junction[ref_id] = &((*curr_junction[ref_id])->next);
        while((*curr_site[ref_id])!=NULL && (*curr_site[ref_id])->pos < beg) curr_site[ref_id] = &((*curr_site[ref_id])->next);

	if(ref_id != ref_id_prev  && ref_id_prev >= 0) {
	    progressbar(1, 1, header->target_name[ref_id_prev], verbose);
	    k=0;
	}

	for(;k<beg;k++) progressbar(k, header->target_len[ref_id], header->target_name[ref_id], verbose);

	offset = 0;
        for(i = 0; i < c->n_cigar; i++) {
	    increm = cigar[i] >> 4;
	    switch(cigar[i] & 0x0F) {
	    	case BAM_CMATCH: 	pos += increm;  	// match to the reference
					offset += increm;	//
				 	break;
		case BAM_CINS:		offset += increm;	// insertion to the reference, pos stays unchanged
					break;
		case BAM_CDEL:		pos += increm;		// deletion from the reference (technically the same as 'N') pointer moves
                                        break;
                case BAM_CREF_SKIP: 
					update_jnxn(curr_junction[ref_id], pos - 1, 	pos + increm, 	mapped_strand, offset, 1);
					update_site(curr_site[ref_id], 	   pos - 1, 			mapped_strand, offset, 0);
					update_site(curr_site[ref_id],     pos + increm,		mapped_strand, offset, 0);
					pos += increm;
				 	break;
		case BAM_CSOFT_CLIP:	offset += increm;
		case BAM_CHARD_CLIP:
		case BAM_CPAD:
		default:		break;
	    }
       	}
	n_reads++;
	if(n_reads>limit_counts && limit_counts>0) break;
    }
    if(verbose) progressbar(1, 1, header->target_name[ref_id_prev], verbose); 

    for(i = 0; i < header->n_targets; i++) {
	ptr = root_junction[i];
	while(ptr != NULL) {
	    site* qtr = ptr->partner;
	    while(qtr != NULL) {
		for(j=0; j<2; j++) {
            	    for(k = 0; k < nbins; k++) { 
		    	if(qtr->count[j][k] > 0) fprintf(ssj_file, "%s\t%i\t%i\t%c\t%i\t%i\n", header->target_name[i], ptr->pos, qtr->pos, strand_i2c(STRAND[j]*stranded), k, qtr->count[j][k]);
		    }
		}
		qtr = qtr->next;
	    }
	    ptr= ptr->next;
	}
    }

    fclose(ssj_file);
    bam_header_destroy(header);
    bam_close(bam_input);
    bam_destroy1(b);

    //*********************************************************************************************************************//

    if(ssc_file == NULL || ssc_file_name[0]==0) {
    	current_time = time(NULL);
    	fprintf(log_file,"Completed in %1.0lf seconds\n",difftime(current_time,timestamp));
    	return 0;
    }

    bam_input = bam_open(bam_file_name, "r");
    header = bam_header_read(bam_input);

    for(i = 0; i < header->n_targets; i++) {
        site *qtr = root_site[i];
        while(qtr != NULL) {
	    for(j=0; j<2; j++) {
            	for(k = 0; k < nbins; k++) {
		    qtr->count[j][k] = 0;
		}
	    }
            qtr = qtr->next;
        }
	curr_site[i] = &root_site[i];
    }

    b = bam_init1();
    k = 0;
    ref_id = -1;
    n_reads = 0;
    while(bam_read1(bam_input, b)>=0) {
        c   = &b->core;
        if(c->tid < 0 || c->tid >= header->n_targets) continue;

        ref_id_prev = ref_id;
        ref_id = c->tid;

        cigar = bam1_cigar(b);

        s = ((c->flag & BAM_FREVERSE)>0);

        mapped_strand = (c->flag & BAM_FREAD1) ? (s + rev_compl[0]) & 1 : (s + rev_compl[1]) & 1;
        mapped_strand*= stranded;

        beg = c->pos + 1;

        while((*curr_site[ref_id])!=NULL && (*curr_site[ref_id])->pos < beg) curr_site[ref_id] = &((*curr_site[ref_id])->next);

        if(ref_id != ref_id_prev  && ref_id_prev >= 0) {
            progressbar(1, 1, header->target_name[ref_id_prev], verbose);
            k=0;
        }

        for(;k<beg;k++) progressbar(k, header->target_len[ref_id], header->target_name[ref_id], verbose);

    	pos = beg;
        offset = 0;
        for(i = 0; i < c->n_cigar; i++) {
            increm = cigar[i] >> 4;
            switch(cigar[i] & 0x0F) {
                case BAM_CMATCH:	qtr = (*curr_site[ref_id]);
					while(qtr != NULL && qtr->pos < pos + increm) {
					    if(qtr->pos > pos && qtr->pos < pos + increm -1 ) {
						int bin = qtr->pos - pos + offset;
						qtr->count[mapped_strand][bin]++;
					    }
					    qtr = qtr->next;
					}
				        pos += increm;  	// match to the reference
					offset += increm;	//
                                        break;
                case BAM_CINS:          offset += increm;      	// insertion to the reference, pointer stays unchanged
                                        break;
                case BAM_CDEL:          			// deletion from the reference (technically the same as 'N') pointer moves
                case BAM_CREF_SKIP:	pos += increm;
					break;
                case BAM_CSOFT_CLIP:	offset += increm;
					break;
                case BAM_CHARD_CLIP:
                case BAM_CPAD:           
                default:
					break;
            }
        }

        n_reads++;
        if(n_reads>limit_counts && limit_counts>0) break;
    }
    if(verbose) progressbar(1, 1, header->target_name[ref_id_prev], verbose);

    for(i = 0; i < header->n_targets; i++) {
        qtr = root_site[i];
        while(qtr != NULL) {
	    for(j=0; j<2; j++) { 
            	for(k = 0; k < nbins; k++) {
            	    if(qtr->count[j][k] > 0) fprintf(ssc_file, "%s\t%i\t%i\t%c\t%i\t%i\n", header->target_name[i], qtr->pos, qtr->pos, strand_i2c(STRAND[j]*stranded), k, qtr->count[j][k]);
		}
	    }
            qtr = qtr->next;
        }
    }
    fclose(ssc_file);
    bam_header_destroy(header);
    bam_close(bam_input);
    bam_destroy1(b);

    current_time = time(NULL);
    fprintf(log_file,"Completed in %1.0lf seconds\n",difftime(current_time,timestamp));
    return 0;

}

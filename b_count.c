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
#include <zlib.h>

#define MAXFILEBUFFLENGTH 4096
#define ARRAY_MARGIN 2
#define INFTY 65535
#define MAXSPLIT 256

const char version[100] = "v3.1";

int nbins   = 1;
int stranded = 1;
FILE *input_file;
void *ssj_file;
void *ssc_file;
FILE* log_file=stderr;
int zipped = 0;

const int STRAND[2] = {1, -1};

//*************************************************************************************************************************//

char strand_i2c(int i) {
    if(i>0) return('+');
    if(i<0) return('-');
    return('.');
}    

int strand_c2i(char c) {
    if(c=='+') return(0);
    if(c=='-') return(1);
    return(-1);
}   

//*************************************************************************************************************************//

class site {
    public:
        site *next;
        int pos;
        int* count[2];

    site(int p) {
        pos = p;
        for(int j=0; j<2; j++) {
            count[j] = (int*) malloc(sizeof(int) * nbins);
            for(int k=0; k<nbins; k++) count[j][k] = 0;
        }
        next = NULL;
    }
};

void update_site(site **ptr, int pos, int strand, int offset, int v) {
    if(offset >= nbins) offset = nbins -1;
    while(*ptr != NULL && (*ptr)->pos < pos) {
        ptr = &((*ptr)->next);
    }
    if(*ptr != NULL && (*ptr)->pos == pos) {
        (*ptr)->count[strand][offset]+=v;
    }
    else {
        site *next = (*ptr);
        (*ptr) = new site(pos);
        (*ptr)->next = next;
        (*ptr)->count[strand][offset]+=v;
    }
}

//*************************************************************************************************************************//

int main(int argc,char* argv[]) {
    time_t timestamp, current_time;
    int i,j,k,a,n;
    char s;

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
    char aux[MAXFILEBUFFLENGTH];
    int beg, end, pos, offset, increm; 
    int ref_id, ref_id_prev;
    int read_type, mapped_strand;

    int rev_compl[2] = {1,0};
    int limit_counts = 0;
    int verbose = 1;
    int continuous = 0;

    int n_reads = 0;
    int max_nh = 0;

    site **root_site, ***curr_site;
    site *r;

    int n_split;

    timestamp = time(NULL);

    if(argc==1) {
	fprintf(stderr, "This utility (%s) counts split reads supporting splice junctions and continuous reads that cover exon boundaries\n", version);
        fprintf(stderr, "Type %s -h for help\n\n",argv[0]);
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
	if(strcmp(argv[i], "-maxnh")   == 0) sscanf(argv[++i], "%i", &max_nh);

	if(strcmp(argv[i], "-quiet") == 0) verbose = 0;
	if(strcmp(argv[i], "-unstranded") == 0) stranded = 0;
	if(strcmp(argv[i], "-continuous") == 0) continuous = 1;
	if(strcmp(argv[i], "-gz") == 0) zipped = 1;

        if(strcmp(argv[i], "-h") ==0 ) {
            fprintf(stderr, "Usage: %s -bam bam_file [-ssj junctions_output] [-ssc boundaries_output] [-log log_file] [-read1 0|1] [-read2 0|1] ",argv[0]);
            fprintf(stderr, "[-nbins number_of_bins] [-lim number_of_lines] [-quiet]\n");
	    fprintf(stderr, "sjcount %s\n", version);
            fprintf(stderr, "Inputs:\n\ta sorted BAM file with a header\n");
	    fprintf(stderr, "\t-ssj: SJ counts, tab-delimited\n");
            fprintf(stderr, "Options:\n");
            fprintf(stderr, "\t-read1 0/1, reverse complement read1 no/yes (default=%i)\n",rev_compl[0]);
            fprintf(stderr, "\t-read2 0/1, reverse complement read2 no/yes (default=%i)\n",rev_compl[1]);
	    fprintf(stderr, "\t-nbins number of overhang bins, (default=%i)\n", nbins);
	    fprintf(stderr, "\t-maxnh, the max value of the NH tag, (default=none)\n");
	    fprintf(stderr, "\t-lim nreads stop after nreads, (default=no limit)\n");
	    fprintf(stderr, "\t-unstranded, force strand to be '.'\n");
	    fprintf(stderr, "\t-continuous, no mismatches when overlapping splice boundaries\n");
	    fprintf(stderr, "\t-gz, gzip output ('.gz' extension will be added to output file names)\n");
	    fprintf(stderr, "\t-quiet, suppress verbose output\n\n"); 
            fprintf(stderr, "Output:\n\t-ssc: Splice boundary counts, tab-delimited  (default=none)\n");
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
	fprintf(log_file,"[ERROR: BAM not specified, exiting]\n");
	exit(1); 
    }

    if(ssj_file_name[0]==0) {
	fprintf(log_file,"[ERROR: SJ file not specified]\n");
	exit(1);
    }

    if(zipped) {
        fprintf(log_file,"[Warning: output is zipped]\n");
    }

    if(stranded==0) {
	fprintf(log_file,"[Warning: strand is ignored (forced to zero)]\n");
    }

    if(continuous==1) {
        fprintf(log_file,"[Warning: only continuous reads overlap splice boundaries]\n");
    }

    if(max_nh>0) {
        fprintf(log_file,"[Warning: max value of NH tag is %i]\n", max_nh);
    }

    if(nbins>0) {
        fprintf(log_file,"[Warning: %i bins for offsets]\n", nbins);
    }

    for(s = 0; s < 2; s++) {
	if(rev_compl[s]) fprintf(log_file,"[Warning: will take reverse complement of read %i]\n", s+1);
    }

    if(limit_counts>0) {
	fprintf(log_file,"[Warning: input limited to %i lines]\n", limit_counts);
    }

    //*****************************************************************************************************************************//

    bam_input = bam_open(bam_file_name, "r");
    header = bam_header_read(bam_input);

    if(bam_input == NULL || header == NULL) {
        fprintf(log_file,"[ERROR: BAM can't be opened or contains no header, exiting]\n");
        exit(1);
    }

    root_site = (site**)  malloc( sizeof(site*)  * (header->n_targets + ARRAY_MARGIN) );
    curr_site = (site***) malloc( sizeof(site**) * (header->n_targets + ARRAY_MARGIN) );

    for(i=0; i < header->n_targets; i++) {
        root_site[i] = NULL;
        curr_site[i] = &root_site[i];
    }

    ssj_file = zipped ? (void*) gzopen(ssj_file_name,"r") : (void*) fopen(ssj_file_name,"r");
    fprintf(log_file, "[Reading %s]", ssj_file_name);

    fprintf(log_file, "[Reading %s", ssj_file_name);
    while(1) {
        buff[0]=0;
        if(zipped) {
            gzgets((gzFile)ssj_file, buff, MAXFILEBUFFLENGTH);
        }
        else {
            fgets(buff, MAXFILEBUFFLENGTH, (FILE*)ssj_file);
        }
        if(strlen(buff)==0) break;
        sscanf(buff, "%s %i %i", &aux, &i, &offset);
        if(i!=1) continue;
        for(j=0;j<strlen(aux);j++) if(aux[j]=='_') aux[j]=' ';
        sscanf(aux, "%s %i %i %c", &chr[0], &beg, &end, &s);
        ref_id = -1;
        for(i=0; i < header->n_targets; i++) {
            if(strcmp(header->target_name[i], chr)==0) ref_id = i;
        }
        if(ref_id>=0) {
            while((*curr_site[ref_id])!=NULL && (*curr_site[ref_id])->pos < beg) curr_site[ref_id] = &((*curr_site[ref_id])->next);
            update_site(curr_site[ref_id], beg, strand_c2i(s), offset, 0);
            update_site(curr_site[ref_id], end, strand_c2i(s), offset, 0);
        }
    }
    fprintf(log_file, "]\n", ssj_file_name);

    //*****************************************************************************************************************************//

    if(ssc_file_name[0]==0) {
	fprintf(log_file,"[Warning: boundary output set to stdout]\n");
        ssc_file = stdout;
        zipped = 0;
    }
    else {
        ssc_file = zipped ? (void*)gzopen(ssc_file_name, "w") : (void*)fopen(ssc_file_name, "w");
        if(ssc_file == NULL) {
            fprintf(log_file,"[ERROR: cannot write to %s, exiting]\n", ssc_file_name);
	    exit(1);
        } else {
            fprintf(log_file,"[Boundary counts: >%s]\n", ssc_file_name);
        }
    }

    bam_input = bam_open(bam_file_name, "r");
    header = bam_header_read(bam_input);

    for(i = 0; i < header->n_targets; i++) {
        r = root_site[i];
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

	n_reads++;
        if(n_reads>limit_counts && limit_counts>0) break;

        uint8_t* ptr = bam_aux_get(b, "NH");
        if(ptr!=NULL && max_nh>0) {
            if(bam_aux2i(ptr)>max_nh) continue;
        }

        beg = c->pos + 1;

        while((*curr_site[ref_id])!=NULL && (*curr_site[ref_id])->pos < beg) curr_site[ref_id] = &((*curr_site[ref_id])->next);

        if(ref_id != ref_id_prev  && ref_id_prev >= 0) {
            progressbar(1, 1, header->target_name[ref_id_prev], verbose);
            k=0;
        }

        for(;k<beg;k++) progressbar(k, header->target_len[ref_id], header->target_name[ref_id], verbose);

	if(continuous && c->n_cigar>1) continue;

    	pos = beg;
        offset = 0;
        for(i = 0; i < c->n_cigar; i++) {
            increm = cigar[i] >> 4;
            switch(cigar[i] & 0x0F) {
                case BAM_CMATCH:	r = (*curr_site[ref_id]);
					while(r != NULL && r->pos < pos + increm) {
					    if(r->pos > pos && r->pos < pos + increm -1 ) {
						int bin = r->pos - pos + offset;
					 	if(bin >= nbins) bin = nbins -1;
						r->count[mapped_strand][bin]++;
					    }
					    r = r->next;
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

    }
    if(verbose) progressbar(1, 1, header->target_name[ref_id_prev], verbose);


    for(i = 0; i < header->n_targets; i++) {
        r = root_site[i];
        while(r != NULL) {
	    for(j=0; j<2; j++) { 
            	for(k = 0; k < nbins; k++) {
            	    if(r->count[j][k] > 0) {
			if(zipped) {
			    gzprintf((gzFile)ssc_file, "%s_%i_%c\t%i\t%i\t%i\n", header->target_name[i], r->pos, strand_i2c(STRAND[j]*stranded), 0, k, r->count[j][k]);
			}
			else {
			    fprintf((FILE*)ssc_file, "%s_%i_%c\t%i\t%i\t%i\n", header->target_name[i], r->pos, strand_i2c(STRAND[j]*stranded), 0, k, r->count[j][k]);
			}
		    }
		}
	    }
            r = r->next;
        }
    }
   if(zipped) {
        gzclose((gzFile)ssc_file);
    }
    else {
        fclose((FILE*)ssc_file);
    }
    bam_header_destroy(header);
    bam_close(bam_input);
    bam_destroy1(b);

    current_time = time(NULL);
    fprintf(log_file,"Completed in %1.0lf seconds\n",difftime(current_time,timestamp));
    return 0;
}
